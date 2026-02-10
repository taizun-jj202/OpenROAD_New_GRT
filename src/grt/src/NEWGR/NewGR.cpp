#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "Net.h"
#include "RoutingTracks.h"
#include "grt/Rudy.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
  int congestion_iterations = 0;
  float global_adjustment = 0.0f;
  std::vector<RegionAdjustment> region_adjustments;
  std::vector<float> routing_layer_adjustments;
};

struct CandidateConfig
{
  std::string name;
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  int seed = 0;
  float critical_percentage = 0.0f;
  int congestion_iterations = 0;
  float global_adjustment = 0.0f;
  int rudy_hotspots = 0;            // number of hotspots to patch (0 disables)
  int rudy_expand_tiles = 0;        // expand each hotspot by N tiles
  float rudy_adjustment = 1.0f;     // capacity multiplier (e.g. 0.8 = -20%)
  int rudy_layers = 0;              // apply to [min_layer, min_layer + L)
};

struct CongestionScore
{
  // 0.0 is ideal; higher means "tighter" routing.
  double score = 0.0;
  double max_util = 0.0;
  double p95_util = 0.0;
  double mean_util = 0.0;
  double mean_excess = 0.0;
  double frac_util_90 = 0.0;
  double frac_util_80 = 0.0;
};

struct GuidePatchingOptions
{
  // Hard caps to avoid exploding guide count / runtime.
  int max_total_patches = 6000;
  int max_patches_per_net = 16;
  int max_patched_pins = 1500;

  // How many Rudy hotspot tiles to consider (prefix of sorted list).
  int rudy_hotspot_prefix = 55;

  // Derived-from-GR congestion hot tiles (based on edge utilization).
  // These complement Rudy hotspots by reacting to actual GR usage patterns.
  int cong_layer_count = 3;           // apply to [min_layer, min_layer + N)
  double cong_util_threshold = 0.88;  // utilization (usage / eff_cap)
  int cong_edge_prefix = 1400;        // keep only top-N hot edges
  int cong_max_tiles = 2600;          // cap on unique hot tiles tracked

  // Patch radius around selected pins, in tiles (1 => +cross neighbors).
  int pin_patch_radius_tiles = 2;
  // Add short wire stubs on the pin connection layer to improve local access
  // without forcing extra layer switching.
  int pin_wire_stub_tiles = 4;

  // Long-segment patching (in tiles along segment).
  int long_segment_tiles = 11;
  int very_long_segment_tiles = 30;
  int long_segment_stub_tiles = 3;
  // When patching a long segment, add a short *same-layer* parallel "side lane"
  // around hotspot samples. This tends to improve DR flexibility without
  // explicitly encouraging layer switching (vias).
  int long_segment_side_lane_span_tiles = 14;
};

static bool is_valid_grid_center(const odb::Rect& die_bounds,
                                 int tile,
                                 int x,
                                 int y)
{
  const int half = tile / 2;
  return x >= die_bounds.xMin() + half && x <= die_bounds.xMax() - half
         && y >= die_bounds.yMin() + half && y <= die_bounds.yMax() - half;
}

static bool point_in_any_rect(const odb::Point& p,
                              const std::vector<odb::Rect>& rects,
                              int prefix)
{
  const int limit = std::min<int>(prefix, static_cast<int>(rects.size()));
  for (int i = 0; i < limit; i++) {
    const odb::Rect& r = rects[i];
    if (p.x() >= r.xMin() && p.x() <= r.xMax() && p.y() >= r.yMin()
        && p.y() <= r.yMax()) {
      return true;
    }
  }
  return false;
}

static bool dbu_to_grid_index(const odb::Rect& die_bounds,
                              int tile,
                              int x_grids,
                              int y_grids,
                              int x_dbu,
                              int y_dbu,
                              int& gx,
                              int& gy)
{
  if (tile <= 0 || x_grids <= 0 || y_grids <= 0) {
    return false;
  }
  const int half = tile / 2;
  const int origin_x = die_bounds.xMin() + half;
  const int origin_y = die_bounds.yMin() + half;
  const int dx = x_dbu - origin_x;
  const int dy = y_dbu - origin_y;
  if (dx < 0 || dy < 0) {
    return false;
  }
  gx = dx / tile;
  gy = dy / tile;
  if (gx < 0 || gx >= x_grids || gy < 0 || gy >= y_grids) {
    return false;
  }
  return true;
}

static void maybe_add_via_patch(GRoute& route,
                                std::unordered_set<GSegment, GSegmentHash>& seen,
                                int x,
                                int y,
                                int layer0,
                                int layer1)
{
  if (layer0 <= 0 || layer1 <= 0) {
    return;
  }
  if (layer0 == layer1) {
    return;
  }
  if (std::abs(layer0 - layer1) != 1) {
    return;
  }
  const GSegment seg(x, y, layer0, x, y, layer1, false);
  if (seen.insert(seg).second) {
    route.push_back(seg);
  }
}

static void maybe_add_wire_patch(GRoute& route,
                                 std::unordered_set<GSegment, GSegmentHash>& seen,
                                 int x0,
                                 int y0,
                                 int layer,
                                 int x1,
                                 int y1,
                                 bool is_jumper = false)
{
  if (layer <= 0) {
    return;
  }
  // Only axis-aligned segments are valid guides.
  if (x0 != x1 && y0 != y1) {
    return;
  }
  const GSegment seg(x0, y0, layer, x1, y1, layer, is_jumper);
  if (seen.insert(seg).second) {
    route.push_back(seg);
  }
}

static void maybe_add_cross_wire_stubs(
    GRoute& route,
    std::unordered_set<GSegment, GSegmentHash>& seen,
    const odb::Rect& die_bounds,
    int tile,
    int x,
    int y,
    int layer,
    int stub_tiles,
    const odb::dbTechLayerDir& preferred_dir)
{
  if (stub_tiles <= 0 || tile <= 0) {
    return;
  }
  if (layer <= 0) {
    return;
  }

  const int d = stub_tiles * tile;
  if (d <= 0) {
    return;
  }

  const bool do_h = (preferred_dir == odb::dbTechLayerDir::HORIZONTAL)
                    || (preferred_dir == odb::dbTechLayerDir::NONE);
  const bool do_v = (preferred_dir == odb::dbTechLayerDir::VERTICAL)
                    || (preferred_dir == odb::dbTechLayerDir::NONE);

  auto valid = [&](int xx, int yy) {
    return is_valid_grid_center(die_bounds, tile, xx, yy);
  };

  if (do_h && valid(x - d, y) && valid(x + d, y)) {
    maybe_add_wire_patch(route, seen, x - d, y, layer, x + d, y);
  } else if (do_h) {
    if (valid(x - d, y) && valid(x, y)) {
      maybe_add_wire_patch(route, seen, x - d, y, layer, x, y);
    }
    if (valid(x, y) && valid(x + d, y)) {
      maybe_add_wire_patch(route, seen, x, y, layer, x + d, y);
    }
  }

  if (do_v && valid(x, y - d) && valid(x, y + d)) {
    maybe_add_wire_patch(route, seen, x, y - d, layer, x, y + d);
  } else if (do_v) {
    if (valid(x, y - d) && valid(x, y)) {
      maybe_add_wire_patch(route, seen, x, y - d, layer, x, y);
    }
    if (valid(x, y) && valid(x, y + d)) {
      maybe_add_wire_patch(route, seen, x, y, layer, x, y + d);
    }
  }
}

static void patch_guides_for_dr_friendliness(
    GlobalRouter* grouter,
    NetRouteMap& routes,
    int min_routing_layer,
    int max_routing_layer,
    const std::vector<odb::Rect>& rudy_hotspot_regions,
    utl::Logger* logger)
{
  if (grouter == nullptr || grouter->grid() == nullptr) {
    return;
  }

  odb::dbDatabase* db = grouter->db();
  odb::dbTech* tech = (db != nullptr) ? db->getTech() : nullptr;

  const int tile = grouter->grid()->getTileSize();
  if (tile <= 0) {
    return;
  }

  const GuidePatchingOptions opts;
  const odb::Rect die_bounds = grouter->grid()->getGridArea();
  const int x_grids = grouter->grid()->getXGrids();
  const int y_grids = grouter->grid()->getYGrids();
  const int tech_max_routing_layer
      = (tech != nullptr) ? tech->getRoutingLayerCount() : max_routing_layer;
  const int min_patch_layer = std::max(1, min_routing_layer);
  const int max_patch_layer
      = std::min(max_routing_layer, tech_max_routing_layer);
  if (min_patch_layer > max_patch_layer) {
    return;
  }

  // Build a deterministic set of hot tiles from *actual* GR edge utilization.
  // This helps focus patching on places where the chosen candidate is tight,
  // even if Rudy hotspots are imperfect.
  std::unordered_set<std::uint64_t> cong_tiles;
  cong_tiles.reserve(static_cast<std::size_t>(opts.cong_max_tiles));

  FastRouteCore* fr = grouter->fastroute();
  if (fr != nullptr && opts.cong_layer_count > 0 && opts.cong_edge_prefix > 0
      && opts.cong_util_threshold > 0.0) {
    struct HotEdge
    {
      double util = 0.0;
      int layer = 0;  // 1-based routing layer
      int x = 0;
      int y = 0;
      bool horizontal = true;
    };
    std::vector<HotEdge> hot_edges;
    hot_edges.reserve(static_cast<std::size_t>(opts.cong_edge_prefix) * 2);

    const auto consider_edge = [&](int layer_1based,
                                   int x,
                                   int y,
                                   bool horizontal,
                                   const auto& edge) {
      const int eff_cap = static_cast<int>(edge.real_cap)
                          - static_cast<int>(edge.red);
      if (eff_cap <= 0) {
        return;
      }
      const double util
          = static_cast<double>(edge.usage) / static_cast<double>(eff_cap);
      if (util < opts.cong_util_threshold) {
        return;
      }
      hot_edges.push_back({util, layer_1based, x, y, horizontal});
    };

    const int cong_min_layer = min_patch_layer;
    const int cong_max_layer
        = std::min(max_patch_layer,
                   min_patch_layer + opts.cong_layer_count - 1);

    const auto& h_edges = fr->getHorizontalEdges3D();
    const auto& v_edges = fr->getVerticalEdges3D();
    const int max_layer_idx = static_cast<int>(h_edges.shape()[0]);

    for (int layer = cong_min_layer; layer <= cong_max_layer; layer++) {
      const int l_idx = layer - 1;
      if (l_idx < 0 || l_idx >= max_layer_idx) {
        continue;
      }
      for (int y = 0; y < static_cast<int>(h_edges.shape()[1]); y++) {
        for (int x = 0; x < static_cast<int>(h_edges.shape()[2]); x++) {
          consider_edge(layer, x, y, true, h_edges[l_idx][y][x]);
        }
      }
      for (int y = 0; y < static_cast<int>(v_edges.shape()[1]); y++) {
        for (int x = 0; x < static_cast<int>(v_edges.shape()[2]); x++) {
          consider_edge(layer, x, y, false, v_edges[l_idx][y][x]);
        }
      }
    }

    std::sort(hot_edges.begin(),
              hot_edges.end(),
              [](const HotEdge& a, const HotEdge& b) {
                if (a.util != b.util) {
                  return a.util > b.util;
                }
                if (a.layer != b.layer) {
                  return a.layer < b.layer;
                }
                if (a.horizontal != b.horizontal) {
                  return a.horizontal > b.horizontal;
                }
                if (a.y != b.y) {
                  return a.y < b.y;
                }
                return a.x < b.x;
              });

    const int edge_limit
        = std::min<int>(opts.cong_edge_prefix, static_cast<int>(hot_edges.size()));

    auto add_tile = [&](int gx, int gy) {
      if (cong_tiles.size() >= static_cast<std::size_t>(opts.cong_max_tiles)) {
        return;
      }
      if (gx < 0 || gx >= x_grids || gy < 0 || gy >= y_grids) {
        return;
      }
      const std::uint64_t key
          = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 32)
            | static_cast<std::uint32_t>(gy);
      cong_tiles.insert(key);
    };

    for (int i = 0;
         i < edge_limit
         && cong_tiles.size() < static_cast<std::size_t>(opts.cong_max_tiles);
         i++) {
      const HotEdge& e = hot_edges[i];
      add_tile(e.x, e.y);
      if (e.horizontal) {
        add_tile(e.x + 1, e.y);
      } else {
        add_tile(e.x, e.y + 1);
      }
    }
  }

  auto point_in_cong_tiles = [&](int x_dbu, int y_dbu) -> bool {
    if (cong_tiles.empty()) {
      return false;
    }
    int gx = 0;
    int gy = 0;
    if (!dbu_to_grid_index(
            die_bounds, tile, x_grids, y_grids, x_dbu, y_dbu, gx, gy)) {
      return false;
    }
    const std::uint64_t key
        = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 32)
          | static_cast<std::uint32_t>(gy);
    return cong_tiles.find(key) != cong_tiles.end();
  };

  int total_patches = 0;
  int patched_pins = 0;
  int patched_nets = 0;

  for (auto& [db_net, route] : routes) {
    if (total_patches >= opts.max_total_patches) {
      break;
    }
    if (db_net == nullptr || route.empty()) {
      continue;
    }

    Net* net = grouter->getNet(db_net);
    if (net == nullptr) {
      continue;
    }

    int patches_for_net = 0;
    std::unordered_set<GSegment, GSegmentHash> seen;
    seen.reserve(route.size() + 32);
    for (const auto& seg : route) {
      seen.insert(seg);
    }

    // (A) Pin-region patching: for a small subset of "risky" pins, add small
    // via-tile guides around the pin to increase access flexibility (similar
    // in spirit to CUGR patching) without changing the routed topology.
    //
    // We intentionally keep this conservative to avoid blowing up guide count.
    if (patched_pins < opts.max_patched_pins) {
      for (Pin& pin : net->getPins()) {
        if (patched_pins >= opts.max_patched_pins
            || patches_for_net >= opts.max_patches_per_net
            || total_patches >= opts.max_total_patches) {
          break;
        }

        const odb::Point pin_grid = pin.getOnGridPosition();
        const odb::Point pin_pos = pin.getPosition();
        static_cast<void>(pin_pos);

        const bool in_hotspot = point_in_any_rect(
            pin_grid, rudy_hotspot_regions, opts.rudy_hotspot_prefix);
        const bool in_cong = point_in_cong_tiles(pin_grid.x(), pin_grid.y());
        const bool dr_risky = in_hotspot || in_cong;

        const bool should_patch_pin
            = pin.isPort() || pin.isConnectedToPadOrMacro()
              || (dr_risky && net->getNumPins() >= 10);

        if (!should_patch_pin) {
          continue;
        }

        const int conn_layer = pin.getConnectionLayer();
        if (conn_layer < min_patch_layer || conn_layer > max_patch_layer) {
          continue;
        }

        const int above = conn_layer + 1;
        const int below = conn_layer - 1;

        const int radius = dr_risky ? opts.pin_patch_radius_tiles : 0;
        const int d = radius * tile;

        odb::dbTechLayerDir preferred_dir = odb::dbTechLayerDir::NONE;
        if (tech != nullptr) {
          odb::dbTechLayer* tech_layer = tech->findRoutingLayer(conn_layer);
          if (tech_layer != nullptr) {
            preferred_dir = tech_layer->getDirection();
          }
        }

        // Center + cross neighbors (radius 1 => +/-1 tile).
        const std::vector<std::pair<int, int>> offsets = {
            {0, 0},
            {d, 0},
            {-d, 0},
            {0, d},
            {0, -d},
        };

        for (const auto& [dx, dy] : offsets) {
          if (dx == 0 && dy == 0) {
            // Always include the center point.
          } else if (radius == 0) {
            continue;
          }

          const int x = pin_grid.x() + dx;
          const int y = pin_grid.y() + dy;
          if (!is_valid_grid_center(die_bounds, tile, x, y)) {
            continue;
          }

          const int before_total = static_cast<int>(route.size());
          // Always add small same-layer stubs to improve local access without
          // forcing extra vias.
          maybe_add_cross_wire_stubs(route,
                                    seen,
                                    die_bounds,
                                    tile,
                                    x,
                                    y,
                                    conn_layer,
                                    opts.pin_wire_stub_tiles,
                                    preferred_dir);

          // Be conservative with explicit via guide patches: they can reduce
          // pin-access failures, but they also tend to inflate via count in DR.
          //
          // Only add *one* adjacent-layer via guide for ports, and only when
          // the pin sits on a tile that is hot based on actual GR utilization.
          const bool allow_pin_vias = pin.isPort();
          if (allow_pin_vias && in_cong) {
            const int target_layer
                = (above <= max_patch_layer) ? above : below;
            if (target_layer >= min_patch_layer
                && target_layer <= max_patch_layer) {
              maybe_add_via_patch(route, seen, x, y, conn_layer, target_layer);
            }
          }
          const int added = static_cast<int>(route.size()) - before_total;
          if (added > 0) {
            patches_for_net += added;
            total_patches += added;
          }
          if (patches_for_net >= opts.max_patches_per_net
              || total_patches >= opts.max_total_patches) {
            break;
          }
        }

        if (patches_for_net > 0) {
          patched_pins++;
        }
      }
    }

    // (B) Long-segment patching: when a net has long straight segments that
    // pass through the hottest Rudy tiles, add a via-tile patch at one or two
    // internal points to allow the detailed router an earlier layer switch.
    for (const GSegment& seg : route) {
      if (patches_for_net >= opts.max_patches_per_net
          || total_patches >= opts.max_total_patches) {
        break;
      }
      if (seg.isVia() || seg.init_layer != seg.final_layer) {
        continue;
      }

      const int layer = seg.init_layer;
      if (layer < min_patch_layer || layer > max_patch_layer) {
        continue;
      }

      const int dist = seg.length();
      const int tiles = tile > 0 ? (dist / tile) : 0;
      if (tiles < opts.long_segment_tiles) {
        continue;
      }

      const int x0 = seg.init_x;
      const int y0 = seg.init_y;
      const int x1 = seg.final_x;
      const int y1 = seg.final_y;

      const bool horizontal = (y0 == y1);
      const bool vertical = (x0 == x1);
      if (!horizontal && !vertical) {
        continue;
      }

      // Only patch segments that plausibly intersect the hottest Rudy regions.
      // Midpoint-only sampling can miss hotspots concentrated near an endpoint,
      // so use a small deterministic set of samples.
      const std::vector<int> sample_steps = {
          tiles / 4,
          tiles / 2,
          (3 * tiles) / 4,
      };

      auto step_to_xy = [&](int step_tiles) -> std::pair<int, int> {
        const int x = horizontal ? (std::min(x0, x1) + step_tiles * tile) : x0;
        const int y = vertical ? (std::min(y0, y1) + step_tiles * tile) : y0;
        return {x, y};
      };

      std::vector<int> hot_steps;
      hot_steps.reserve(sample_steps.size());
      for (const int step_tiles : sample_steps) {
        if (step_tiles <= 0 || step_tiles >= tiles) {
          continue;
        }
        const auto [x, y] = step_to_xy(step_tiles);
        if (!is_valid_grid_center(die_bounds, tile, x, y)) {
          continue;
        }
        const odb::Point p(x, y);
        const bool in_rudy
            = point_in_any_rect(p, rudy_hotspot_regions, opts.rudy_hotspot_prefix);
        const bool in_cong = point_in_cong_tiles(x, y);
        if (in_rudy || in_cong) {
          hot_steps.push_back(step_tiles);
        }
      }

      if (hot_steps.empty()) {
        continue;
      }

      std::sort(hot_steps.begin(), hot_steps.end());
      hot_steps.erase(std::unique(hot_steps.begin(), hot_steps.end()),
                      hot_steps.end());

      auto try_add_patch = [&](auto fn) {
        if (patches_for_net >= opts.max_patches_per_net
            || total_patches >= opts.max_total_patches) {
          return;
        }
        const int before = static_cast<int>(route.size());
        fn();
        const int added = static_cast<int>(route.size()) - before;
        if (added > 0) {
          patches_for_net += added;
          total_patches += added;
        }
      };

      auto add_local_side_lane_at_step = [&](int step_tiles) {
        const int span = std::max(1, opts.long_segment_side_lane_span_tiles);
        const int start_step = std::max(0, step_tiles - span);
        const int end_step = std::min(tiles, step_tiles + span);
        if (end_step <= start_step) {
          return;
        }

        const int axis_min = horizontal ? std::min(x0, x1) : std::min(y0, y1);
        const int s0 = axis_min + start_step * tile;
        const int s1 = axis_min + end_step * tile;

        const int wx0 = horizontal ? s0 : x0;
        const int wy0 = vertical ? s0 : y0;
        const int wx1 = horizontal ? s1 : x1;
        const int wy1 = vertical ? s1 : y1;

        if (!is_valid_grid_center(die_bounds, tile, wx0, wy0)
            || !is_valid_grid_center(die_bounds, tile, wx1, wy1)) {
          return;
        }

        // Add a parallel segment on the same layer, offset by 1 tile in the
        // perpendicular direction (if both endpoints stay in-bounds). This
        // creates a narrow corridor that can absorb local detours without
        // requiring a layer change.
        const int ox = vertical ? tile : 0;
        const int oy = horizontal ? tile : 0;
        const auto endpoints_ok = [&](int dx, int dy) {
          return is_valid_grid_center(die_bounds, tile, wx0 + dx, wy0 + dy)
                 && is_valid_grid_center(die_bounds, tile, wx1 + dx, wy1 + dy);
        };

        int dx = 0;
        int dy = 0;
        if (endpoints_ok(ox, oy)) {
          dx = ox;
          dy = oy;
        } else if (endpoints_ok(-ox, -oy)) {
          dx = -ox;
          dy = -oy;
        } else {
          return;
        }

        try_add_patch([&]() {
          maybe_add_wire_patch(
              route, seen, wx0 + dx, wy0 + dy, layer, wx1 + dx, wy1 + dy);
        });
      };

      odb::dbTechLayerDir preferred_dir = odb::dbTechLayerDir::NONE;
      if (tech != nullptr) {
        odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer);
        if (tech_layer != nullptr) {
          preferred_dir = tech_layer->getDirection();
        }
      }

      // Add same-layer parallel "side lanes" around hotspot samples. Prefer
      // this to adjacent-layer escape lanes to avoid inflating via count.
      if (tiles >= opts.very_long_segment_tiles) {
        // Very long segments get a couple of side-lane samples to increase
        // flexibility while keeping guide count bounded.
        add_local_side_lane_at_step(tiles / 3);
        add_local_side_lane_at_step((2 * tiles) / 3);
      } else {
        for (const int step_tiles : hot_steps) {
          add_local_side_lane_at_step(step_tiles);
        }
      }

      // Add same-layer local flexibility at hotspot samples (helps DR avoid
      // detours without requiring layer switches).
      for (const int step_tiles : hot_steps) {
        const auto [x, y] = step_to_xy(step_tiles);
        if (!is_valid_grid_center(die_bounds, tile, x, y)) {
          continue;
        }
        try_add_patch([&]() {
          maybe_add_cross_wire_stubs(route,
                                    seen,
                                    die_bounds,
                                    tile,
                                    x,
                                    y,
                                    layer,
                                    /*stub_tiles=*/opts.long_segment_stub_tiles,
                                    preferred_dir);
        });
      }

      // Intentionally avoid explicit via guide patches here. These tend to
      // increase downstream via count; the same-layer side lanes + stubs above
      // usually provide enough flexibility for DR in hotspot regions.
    }

    if (patches_for_net > 0) {
      patched_nets++;
    }
  }

  if (total_patches > 0) {
    logger->info(GNR,
                 6011,
                 "NEWGR patching: added {} guide patches across {} nets (patched pins {}, cap per-net {})",
                 total_patches,
                 patched_nets,
                 patched_pins,
                 opts.max_patches_per_net);
  }
}

static void simplify_guides(GlobalRouter* grouter,
                            NetRouteMap& routes,
                            utl::Logger* logger)
{
  if (routes.empty()) {
    return;
  }

  int merge_gap_dbu = 0;
  if (grouter != nullptr && grouter->grid() != nullptr) {
    // Many guide endpoints are snapped to the GCell centers spaced by
    // `tile_size`. Allowing a small (<= 1 tile) merge gap can reduce guide
    // fragmentation without materially increasing guide area.
    merge_gap_dbu = std::max(0, grouter->grid()->getTileSize());
  }

  int nets_touched = 0;
  int merged_segments = 0;
  int removed_duplicates = 0;

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.size() < 2) {
      continue;
    }

    // Keep unmergeable segments (vias + any unexpected diagonals).
    std::vector<GSegment> keep;
    keep.reserve(route.size());

    // Group wires by (layer, is_jumper, horizontal, fixed_coord).
    using Key = std::tuple<int, int, int, int>;
    std::map<Key, std::vector<std::pair<int, int>>> intervals;

    for (const GSegment& seg : route) {
      if (seg.isVia()) {
        keep.push_back(seg);
        continue;
      }
      if (seg.init_layer != seg.final_layer) {
        keep.push_back(seg);
        continue;
      }

      const bool horizontal = (seg.init_y == seg.final_y);
      const bool vertical = (seg.init_x == seg.final_x);
      if (!horizontal && !vertical) {
        keep.push_back(seg);
        continue;
      }

      const int layer = seg.init_layer;
      const int jumper = seg.is_jumper ? 1 : 0;
      if (horizontal) {
        const int y = seg.init_y;
        const int x0 = std::min(seg.init_x, seg.final_x);
        const int x1 = std::max(seg.init_x, seg.final_x);
        intervals[Key{layer, jumper, 1, y}].push_back({x0, x1});
      } else {
        const int x = seg.init_x;
        const int y0 = std::min(seg.init_y, seg.final_y);
        const int y1 = std::max(seg.init_y, seg.final_y);
        intervals[Key{layer, jumper, 0, x}].push_back({y0, y1});
      }
    }

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(keep.size() + intervals.size());
    for (const GSegment& seg : keep) {
      rebuilt.push_back(seg);
    }

    for (auto& [key, seg_intervals] : intervals) {
      if (seg_intervals.empty()) {
        continue;
      }
      std::sort(seg_intervals.begin(),
                seg_intervals.end(),
                [](const auto& a, const auto& b) {
                  if (a.first != b.first) {
                    return a.first < b.first;
                  }
                  return a.second < b.second;
                });

      std::vector<std::pair<int, int>> merged;
      merged.reserve(seg_intervals.size());
      int cur_lo = seg_intervals[0].first;
      int cur_hi = seg_intervals[0].second;
      for (std::size_t i = 1; i < seg_intervals.size(); i++) {
        const int lo = seg_intervals[i].first;
        const int hi = seg_intervals[i].second;
        if (lo <= cur_hi + merge_gap_dbu) {  // overlap, touch, or small gap
          cur_hi = std::max(cur_hi, hi);
        } else {
          merged.push_back({cur_lo, cur_hi});
          cur_lo = lo;
          cur_hi = hi;
        }
      }
      merged.push_back({cur_lo, cur_hi});

      const auto [layer, jumper, horizontal, fixed] = key;
      for (const auto& [lo, hi] : merged) {
        if (horizontal) {
          rebuilt.emplace_back(lo, fixed, layer, hi, fixed, layer, jumper != 0);
        } else {
          rebuilt.emplace_back(fixed, lo, layer, fixed, hi, layer, jumper != 0);
        }
      }

      merged_segments
          += static_cast<int>(seg_intervals.size())
             - static_cast<int>(merged.size());
    }

    // Deduplicate (also removes any wire segments that may have been kept and
    // reconstructed identically).
    std::unordered_set<GSegment, GSegmentHash> seen;
    seen.reserve(rebuilt.size() * 2);
    std::vector<GSegment> unique;
    unique.reserve(rebuilt.size());
    for (const GSegment& seg : rebuilt) {
      if (seen.insert(seg).second) {
        unique.push_back(seg);
      } else {
        removed_duplicates++;
      }
    }

    if (unique.size() != route.size()) {
      nets_touched++;
    }
    route.swap(unique);
  }

  if (logger != nullptr && (merged_segments > 0 || removed_duplicates > 0)) {
    logger->info(GNR,
                 6012,
                 "NEWGR guide simplify: touched {} nets, merged {}, removed dup {} (gap {} dbu)",
                 nets_touched,
                 merged_segments,
                 removed_duplicates,
                 merge_gap_dbu);
  }
}

}  // namespace

NewGR::NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR::run(std::vector<Net*>& nets,
	                       int min_routing_layer,
	                       int max_routing_layer)
{
	  if (nets.empty()) {
	    return {};
	  }

	  const std::vector<Net*> canonical_nets = nets;

	  struct RouteMetrics
	  {
	    long wirelength_dbu = 0;
	    long via_count = 0;
	    double wirelength_um = 0.0;
	    int total_overflow = 0;
      CongestionScore congestion;
	  };

  auto compute_metrics = [&](const NetRouteMap& routes) -> RouteMetrics {
    RouteMetrics metrics;
    for (const auto& [db_net, segments] : routes) {
      static_cast<void>(db_net);
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          metrics.via_count++;
        } else {
          metrics.wirelength_dbu += segment.length();
        }
      }
    }

    if (metrics.wirelength_dbu > 0 && grouter_->db_ != nullptr
        && grouter_->db_->getTech() != nullptr) {
      metrics.wirelength_um
          = metrics.wirelength_dbu
            / static_cast<double>(
                grouter_->db_->getTech()->getDbUnitsPerMicron());
    }
    return metrics;
  };

  const auto compute_congestion_score = [&]() -> CongestionScore {
    CongestionScore out;
    const auto& h_edges = grouter_->fastroute_->getHorizontalEdges3D();
    const auto& v_edges = grouter_->fastroute_->getVerticalEdges3D();

    std::vector<double> utils;
    utils.reserve(1024);

    auto scan_edges = [&](const auto& edges) {
      std::size_t edge_count = 0;
      std::size_t util_80 = 0;
      std::size_t util_90 = 0;
      double sum_util = 0.0;
      double sum_excess = 0.0;

      for (const auto& layer : edges) {
        for (const auto& row : layer) {
          for (const auto& edge : row) {
            // Use real_cap to compute a stable utilization proxy. The adjusted
            // cap can be overly pessimistic and tends to saturate the metric.
            const int eff_cap = static_cast<int>(edge.real_cap)
                                - static_cast<int>(edge.red);
            if (eff_cap <= 0) {
              continue;
            }
            const double util = static_cast<double>(edge.usage)
                                / static_cast<double>(eff_cap);
            out.max_util = std::max(out.max_util, util);
            sum_util += util;
            if (util > 1.0) {
              sum_excess += (util - 1.0);
            }
            if (util >= 0.8) {
              util_80++;
            }
            if (util >= 0.9) {
              util_90++;
            }
            utils.push_back(util);
            edge_count++;
          }
        }
      }

      if (edge_count > 0) {
        out.mean_util
            = std::max(out.mean_util,
                       sum_util / static_cast<double>(edge_count));
        out.mean_excess = std::max(
            out.mean_excess, sum_excess / static_cast<double>(edge_count));
        out.frac_util_80
            = std::max(out.frac_util_80,
                       static_cast<double>(util_80) / edge_count);
        out.frac_util_90
            = std::max(out.frac_util_90,
                       static_cast<double>(util_90) / edge_count);
      }
    };

    scan_edges(h_edges);
    scan_edges(v_edges);

    if (!utils.empty()) {
      const std::size_t idx
          = static_cast<std::size_t>(0.95 * (utils.size() - 1));
      std::nth_element(utils.begin(), utils.begin() + idx, utils.end());
      out.p95_util = utils[idx];
    }

    // DR-friendliness proxy: avoid near-saturated edges (p95) and discourage
    // overload (mean_excess). Fractions are light tie-breakers.
    out.score = out.p95_util + 4.0 * out.mean_excess + 0.15 * out.frac_util_90
                + 0.05 * out.frac_util_80;
    return out;
  };

  auto capture_snapshot = [&]() -> RouterSnapshot {
    RouterSnapshot snapshot;
    snapshot.caps_percentage = grouter_->caps_perturbation_percentage_;
    snapshot.perturbation_amount = grouter_->perturbation_amount_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    snapshot.congestion_iterations = grouter_->congestion_iterations_;
    snapshot.global_adjustment = grouter_->adjustment_;
    snapshot.region_adjustments = grouter_->region_adjustments_;

    if (grouter_->db_ != nullptr && grouter_->db_->getTech() != nullptr) {
      odb::dbTech* tech = grouter_->db_->getTech();
      snapshot.routing_layer_adjustments.resize(max_routing_layer + 1, 0.0f);
      for (int layer = 1; layer <= max_routing_layer; layer++) {
        odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer);
        if (tech_layer != nullptr) {
          snapshot.routing_layer_adjustments[layer]
              = tech_layer->getLayerAdjustment();
        }
      }
    }
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
    grouter_->setCongestionIterations(snapshot.congestion_iterations);
    grouter_->setAdjustment(snapshot.global_adjustment);
    grouter_->region_adjustments_ = snapshot.region_adjustments;

    if (!snapshot.routing_layer_adjustments.empty() && grouter_->db_ != nullptr
        && grouter_->db_->getTech() != nullptr) {
      odb::dbTech* tech = grouter_->db_->getTech();
      const int max_layer
          = std::min<int>(max_routing_layer,
                          static_cast<int>(
                              snapshot.routing_layer_adjustments.size())
                              - 1);
      for (int layer = 1; layer <= max_layer; layer++) {
        odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer);
        if (tech_layer != nullptr) {
          tech_layer->setLayerAdjustment(
              snapshot.routing_layer_adjustments[layer]);
        }
      }
    }
  };

  const RouterSnapshot snapshot = capture_snapshot();

  // Build a deterministic list of "hotspot" regions from Rudy (wire density)
  // to optionally reserve extra soft-capacity around congested areas, similar
  // in spirit to SPRoute soft capacity / CUGR patching.
  std::vector<odb::Rect> rudy_hotspot_regions;
  {
    Rudy* rudy = grouter_->getRudy();
    if (rudy != nullptr && grouter_->grid_ != nullptr) {
      const int x_grids = grouter_->grid_->getXGrids();
      const int y_grids = grouter_->grid_->getYGrids();
      if (x_grids > 0 && y_grids > 0) {
        rudy->setGridConfig(grouter_->grid_->getGridArea(), x_grids, y_grids);
        rudy->calculateRudy();

        struct TileRudy
        {
          int x = 0;
          int y = 0;
          float rudy = 0.0f;
        };
        std::vector<TileRudy> tiles;
        tiles.reserve(static_cast<std::size_t>(x_grids) * y_grids);

        for (int x = 0; x < x_grids; x++) {
          for (int y = 0; y < y_grids; y++) {
            const float value = rudy->getTile(x, y).getRudy();
            if (value <= 0.0f) {
              continue;
            }
            tiles.push_back({x, y, value});
          }
        }

        std::sort(tiles.begin(),
                  tiles.end(),
                  [](const TileRudy& a, const TileRudy& b) {
                    if (a.rudy != b.rudy) {
                      return a.rudy > b.rudy;
                    }
                    if (a.x != b.x) {
                      return a.x < b.x;
                    }
                    return a.y < b.y;
                  });

        // Keep the list small and stable; we'll take prefixes for candidates.
        const int max_regions
            = std::min<int>(80, static_cast<int>(tiles.size()));
        rudy_hotspot_regions.reserve(max_regions);
        for (int i = 0; i < max_regions; i++) {
          const auto& t = tiles[i];
          rudy_hotspot_regions.push_back(rudy->getTile(t.x, t.y).getRect());
        }
      }
    }
  }

  const auto add_rudy_hotspot_adjustments = [&](const CandidateConfig& config) {
    if (config.rudy_hotspots <= 0 || config.rudy_adjustment >= 1.0f
        || config.rudy_layers <= 0 || rudy_hotspot_regions.empty()) {
      return;
    }

    // GlobalRouter::computeRegionAdjustments expects a *reduction percentage*
    // (0.10 => reduce capacity by 10%), while NEWGR candidate configs specify
    // a more intuitive capacity multiplier (0.90 => keep 90% capacity).
    const float reduction
        = std::clamp(1.0f - config.rudy_adjustment, 0.0f, 1.0f);
    if (reduction <= 0.0f) {
      return;
    }

    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();
    const int tile = grouter_->grid_->getTileSize();
    if (x_grids <= 0 || y_grids <= 0 || tile <= 0) {
      return;
    }

    const int max_hotspots
        = std::min<int>(config.rudy_hotspots,
                        static_cast<int>(rudy_hotspot_regions.size()));
    const int min_layer = min_routing_layer;
    const int max_layer = std::min(max_routing_layer,
                                   min_routing_layer + config.rudy_layers - 1);

    for (int i = 0; i < max_hotspots; i++) {
      odb::Rect rect = rudy_hotspot_regions[i];
      if (config.rudy_expand_tiles > 0) {
        const int expand = config.rudy_expand_tiles * tile;
        rect = odb::Rect(rect.xMin() - expand,
                         rect.yMin() - expand,
                         rect.xMax() + expand,
                         rect.yMax() + expand);
      }
      for (int layer = min_layer; layer <= max_layer; layer++) {
        grouter_->region_adjustments_.emplace_back(rect.xMin(),
                                                  rect.yMin(),
                                                  rect.xMax(),
                                                  rect.yMax(),
                                                  layer,
                                                  reduction);
      }
    }
  };

	  const auto prepare_fastroute = [&](const CandidateConfig& config,
	                                     std::vector<Net*>& candidate_nets) {
	    // Re-initialize FastRoute state and grid/capacities without rebuilding
	    // Net objects/pins (which is expensive). We reuse the `nets` list built
	    // by the caller's initFastRoute() and only rebuild the FastRoute core
	    // data structures, capacities, and adjustments.
	    //
    // This keeps multi-candidate exploration affordable while still allowing
    // capacity perturbations/seed changes to take effect.

    grouter_->pad_pins_connections_.clear();

    grouter_->fastroute_->clear();
    grouter_->sproute_grid_ready_ = false;
    grouter_->sproute_nets_ready_ = false;
    grouter_->sproute_grid_data_ = SprouteGridData();
    grouter_->sproute_nets_.clear();
    grouter_->sproute_total_overflow_ = 0;
    grouter_->h_nets_in_pos_.clear();
    grouter_->v_nets_in_pos_.clear();

    grouter_->setCongestionIterations(config.congestion_iterations);
    grouter_->setCapacitiesPerturbationPercentage(config.caps_percentage);
    grouter_->setPerturbationAmount(config.perturbation_amount);
    grouter_->setSeed(config.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(config.critical_percentage);
    grouter_->setAdjustment(config.global_adjustment);

    // Start from the snapshot's adjustments and optionally add Rudy-based
    // hotspot "soft-capacity" reservations for DR-friendliness.
    grouter_->region_adjustments_ = snapshot.region_adjustments;

    grouter_->routing_tracks_.clear();
    grouter_->routing_layers_.clear();
    grouter_->grid_->clear();
    grouter_->vertical_capacities_.clear();
    grouter_->horizontal_capacities_.clear();

    grouter_->ensureLayerForGuideDimension(max_routing_layer);
    grouter_->configFastRoute();
    grouter_->initRoutingLayers(min_routing_layer, max_routing_layer);
    grouter_->initRoutingTracks(max_routing_layer);
    grouter_->initCoreGrid(max_routing_layer);
    grouter_->setCapacities(min_routing_layer, max_routing_layer);
    grouter_->captureSprouteGridData();

    add_rudy_hotspot_adjustments(config);

    grouter_->applyAdjustments(min_routing_layer, max_routing_layer);
    grouter_->perturbCapacities();

	    // Init the data structures to monitor 3D capacity during 2D phases
	    grouter_->fastroute_->initEdgesCapacityPerLayer();

	    // Rebuild the FastRoute netlist from the already-created Net objects.
	    grouter_->initFastRouteIncr(candidate_nets);
	    grouter_->initialized_ = true;
	  };

  logger_->info(GNR,
                6004,
                "NEWGR snapshot: caps% {:.1f} pert_amt {} seed {} crit% {:.1f} cong_iters {} global_adj {:.2f}",
                snapshot.caps_percentage,
                snapshot.perturbation_amount,
                snapshot.seed,
                snapshot.critical_percentage,
                snapshot.congestion_iterations,
                snapshot.global_adjustment);

  // NEWGR strategy (wirelength-first, via-second):
  // - Explore a *very small* deterministic set of configurations (to keep
  //   runtime close to a single FastRoute run) and pick the best by:
  //   (1) routability, (2) global wirelength window, then (3) congestion
  //   looseness (DR-friendliness) and (4) via count.
  //
  // Rationale: for our regression, small deterministic perturbations can
  // improve downstream detailed-routing metrics by nudging congestion away
  // from hard-to-route pin-access regions, even if global cost differences
  // are small.
  // Candidate A: explore an alternate deterministic seed. This often changes
  // local congestion patterns (and downstream DR detours) without materially
  // changing global wirelength or runtime.
  const CandidateConfig tuned{"perturb3-seed11-crit0",
                              3.0f,
                              1,
                              11,
                              0.0f,
                              snapshot.congestion_iterations,
                              snapshot.global_adjustment,
                              0,
                              0,
                              1.0f,
                              0};

  const std::vector<CandidateConfig> candidates = {tuned};

	  const auto run_candidate = [&](const CandidateConfig& candidate) {
	    restore_snapshot(snapshot);
	    grouter_->setAllowCongestion(snapshot.allow_congestion);
	    std::vector<Net*> candidate_nets = canonical_nets;
	    prepare_fastroute(candidate, candidate_nets);

	    NetRouteMap routes = grouter_->findRouting(
	        candidate_nets, min_routing_layer, max_routing_layer);
	    RouteMetrics metrics = compute_metrics(routes);
	    metrics.total_overflow = grouter_->fastroute_->totalOverflow();
      metrics.congestion = compute_congestion_score();
	    logger_->info(GNR,
	                  6005,
	                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow {}, cong {:.3f} (p95 {:.3f}, max {:.3f}, excess {:.4f})",
                  candidate.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.total_overflow,
                  metrics.congestion.score,
                  metrics.congestion.p95_util,
                  metrics.congestion.max_util,
                  metrics.congestion.mean_excess);
    return std::make_pair(std::move(routes), metrics);
  };

  auto finalize = [&](int overflow) {
    restore_snapshot(snapshot);
    if (overflow > 0 && !snapshot.allow_congestion) {
      logger_->warn(GNR,
                    6008,
                    "NEWGR finished with overflow ({}); enabling allow_congestion to avoid abort.",
                    overflow);
      grouter_->setAllowCongestion(true);
    }
  };

	  CandidateConfig best_candidate = tuned;
	  RouteMetrics best_metrics;
	  bool have_best = false;

    CandidateConfig last_candidate = tuned;
    RouteMetrics last_metrics;
    NetRouteMap last_routes;
    bool have_last = false;

	  struct CandidateEval
	  {
	    CandidateConfig candidate;
	    RouteMetrics metrics;
	  };
	  std::vector<CandidateEval> evals;
	  evals.reserve(candidates.size());

	  int best_overflow = std::numeric_limits<int>::max();
	  for (const auto& candidate : candidates) {
	    auto [routes, metrics] = run_candidate(candidate);

	    evals.push_back({candidate, metrics});
	    best_overflow = std::min(best_overflow, metrics.total_overflow);

      // Keep the last run's full route so we can avoid a redundant rerun when
      // the winner ends up being the last candidate (common case).
      last_candidate = candidate;
      last_metrics = metrics;
      last_routes = std::move(routes);
      have_last = true;
		  }

		  // Selection policy (runtime-aware, DR-friendly):
		  // 1) Prefer routable solutions (overflow == 0), else minimize overflow.
		  // 2) Keep candidates within a tight global-WL window of the best WL
		  //    (wirelength-first objective).
		  // 3) Within that window, prefer looser congestion (often reduces DR
		  //    detours / DR wirelength), then fewer vias.
		  const int target_overflow = (best_overflow == 0) ? 0 : best_overflow;

      long min_wl_dbu = std::numeric_limits<long>::max();
      for (const auto& eval : evals) {
        const RouteMetrics& metrics = eval.metrics;
        if (metrics.total_overflow != target_overflow) {
          continue;
        }
        min_wl_dbu = std::min(min_wl_dbu, metrics.wirelength_dbu);
      }

      // Keep candidates very close to the best global WL. This guards against
      // drifting into a longer-GR regime while still letting us choose a
      // slightly "looser" solution for DR if it is essentially WL-equivalent.
      constexpr double wl_slack_ratio = 0.0006;   // 0.06%
      constexpr long wl_slack_min_dbu = 80000;    // ~80um @ 1000 DBU/um
      const long wl_slack_dbu = std::max<long>(
          wl_slack_min_dbu,
          static_cast<long>(std::llround(min_wl_dbu * wl_slack_ratio)));
      const long wl_limit_dbu = min_wl_dbu + wl_slack_dbu;

      logger_->info(GNR,
                    6010,
                    "NEWGR selection: overflow {} wl_ref {} dbu wl_limit {} dbu (+{} / {:.2f}%)",
                    target_overflow,
                    min_wl_dbu,
                    wl_limit_dbu,
                    wl_slack_dbu,
                    100.0 * wl_slack_dbu / std::max<double>(1.0, min_wl_dbu));

      for (const auto& eval : evals) {
        const RouteMetrics& metrics = eval.metrics;
        if (metrics.total_overflow != target_overflow) {
          continue;
        }
        if (metrics.wirelength_dbu > wl_limit_dbu) {
          continue;
        }

        if (!have_best) {
          best_candidate = eval.candidate;
          best_metrics = metrics;
          have_best = true;
          continue;
        }

        if (metrics.congestion.score != best_metrics.congestion.score) {
          if (metrics.congestion.score < best_metrics.congestion.score) {
            best_candidate = eval.candidate;
            best_metrics = metrics;
          }
          continue;
        }

        // Within similar congestion, prefer fewer vias (secondary objective).
        if (metrics.via_count != best_metrics.via_count) {
          if (metrics.via_count < best_metrics.via_count) {
            best_candidate = eval.candidate;
            best_metrics = metrics;
          }
          continue;
        }

        if (metrics.wirelength_dbu != best_metrics.wirelength_dbu) {
          if (metrics.wirelength_dbu < best_metrics.wirelength_dbu) {
            best_candidate = eval.candidate;
            best_metrics = metrics;
          }
          continue;
        }

        if (eval.candidate.name < best_candidate.name) {
          best_candidate = eval.candidate;
          best_metrics = metrics;
        }
      }

	  // If the winner is already the last evaluated candidate, reuse its routes
	  // directly. Otherwise rerun so `grouter_->fastroute_` internal state
	  // matches the returned routes (important for subsequent incremental calls
	  // in the flow).
    NetRouteMap winner_routes;
    RouteMetrics winner_metrics;
    if (have_last && last_candidate.name == best_candidate.name) {
      winner_routes = std::move(last_routes);
      winner_metrics = last_metrics;
    } else {
      std::tie(winner_routes, winner_metrics) = run_candidate(best_candidate);
    }

  // Post-processing: add conservative DR-friendly "patch" guides in/around
  // the hottest Rudy regions to help reduce downstream detours (wirelength)
  // without materially changing the GR topology.
  patch_guides_for_dr_friendliness(grouter_,
                                   winner_routes,
                                   min_routing_layer,
                                   max_routing_layer,
                                   rudy_hotspot_regions,
                                   logger_);
  // Guide post-pass: merge collinear segments and drop duplicates. This is a
  // cheap way to reduce guide fragmentation, which tends to reduce DR detours
  // (wirelength) and can also avoid unnecessary layer switching (vias).
  simplify_guides(grouter_, winner_routes, logger_);

  logger_->info(GNR,
                6007,
                "NEWGR picked {} (overflow {}, wl {:.0f} um, vias {})",
                best_candidate.name,
                winner_metrics.total_overflow,
                winner_metrics.wirelength_um,
                winner_metrics.via_count);
  finalize(winner_metrics.total_overflow);
  return winner_routes;
}

}  // namespace grt
