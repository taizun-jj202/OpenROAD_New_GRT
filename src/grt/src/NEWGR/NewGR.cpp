#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "boost/functional/hash.hpp"
#include "FastRoute.h"
#include "Grid.h"
#include "grt/GRoute.h"
#include "stt/SteinerTreeBuilder.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct GuidePatchStats
{
  int nets_patched = 0;
  int vias_added = 0;
  int wires_added = 0;
};

int pick_adjacent_layer(const int layer,
                        const int min_routing_layer,
                        const int max_routing_layer)
{
  if (layer + 1 <= max_routing_layer) {
    return layer + 1;
  }
  if (layer - 1 >= min_routing_layer) {
    return layer - 1;
  }
  return -1;
}

struct CongestionPatchStats
{
  int nets_patched = 0;
  int points_added = 0;
};

struct EndpointAccessPatchStats
{
  int nets_patched = 0;
  int endpoints_patched = 0;
  int vias_added = 0;
};

// Adds "standalone" 1-GCell guide patches as *degenerate* (zero-length) wire
// segments on adjacent layers at congested positions along long guides.
//
// This follows the spirit of the CUGR post-processing patching:
// provide extra escape points / track switching opportunities in places where
// edge spare resources are low, without forcing the global route topology to
// detour (and thus without directly increasing global wirelength).
//
// Important note: these patches are not required to be connected. Detailed
// routing can still exploit them because guides are treated as *allowed
// regions* per layer, and vias are legal at overlapping guide areas.
void add_congestion_patches(NetRouteMap& routes,
                            FastRouteCore* fastroute,
                            const Grid* grid,
                            odb::dbTech* tech,
                            utl::Logger* logger,
                            const int min_routing_layer,
                            const int max_routing_layer)
{
  if (fastroute == nullptr || grid == nullptr || tech == nullptr
      || logger == nullptr) {
    return;
  }

  // Congestion patching is enabled by default in NEWGR. The segments added
  // here are *degenerate points* (0-length wires) and add neither wirelength
  // nor vias by construction. The intent is to give the detailed router extra
  // layer-escape opportunities in congested/low-resource areas without
  // perturbing the global route topology.

  // Keep conservative defaults; can be tuned via debug if needed.
  constexpr int kMaxTotalPatches = 120000;  // guardrail
  constexpr int kMaxPatchesPerNet = 24;
  constexpr int kMinSegmentLenTiles = 8;
  constexpr int kEndMarginTiles = 1;
  constexpr int kSampleStepTiles = 2;
  constexpr int kMinPatchSpacingTiles = 4;

  // In FastRoute-style flows, pin-access issues usually dominate on the
  // bottom layers. Keep patching focused to avoid bloating guides.
  const int max_patch_layer
      = std::min(max_routing_layer, min_routing_layer + 2);

  // Threshold is expressed in "spare tracks" on the edge.
  constexpr int kWirePatchSpareThreshold = 2;
  constexpr int kAdjLayerMinSpare = 1;

  const int tile_size = grid->getTileSize();
  if (tile_size <= 0) {
    return;
  }

  const int grid_xmin = grid->getXMin();
  const int grid_ymin = grid->getYMin();
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();

  auto to_grid = [&](const int x_dbu, const int y_dbu) -> std::pair<int, int> {
    const int gx = (x_dbu - grid_xmin) / tile_size;
    const int gy = (y_dbu - grid_ymin) / tile_size;
    return {gx, gy};
  };

  auto to_dbu = [&](const int gx, const int gy) -> std::pair<int, int> {
    // FastRoute represents each GCell as a node at the tile center:
    //   x = x_corner + tile_size * (gx + 0.5)
    // Use the same convention to keep patches on-grid.
    const int x_dbu = grid_xmin + gx * tile_size + tile_size / 2;
    const int y_dbu = grid_ymin + gy * tile_size + tile_size / 2;
    return {x_dbu, y_dbu};
  };

  auto in_grid = [&](const int gx, const int gy) -> bool {
    return gx >= 0 && gx < x_grids && gy >= 0 && gy < y_grids;
  };

  auto min_spare_at_point_preferred_dir
      = [&](const int layer, const int gx, const int gy) -> int {
    if (layer < min_routing_layer || layer > max_routing_layer) {
      return 0;
    }
    if (!in_grid(gx, gy)) {
      return 0;
    }

    odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer);
    if (tech_layer == nullptr) {
      return 0;
    }

    int spare = std::numeric_limits<int>::max();
    const auto dir = tech_layer->getDirection();
    if (dir == odb::dbTechLayerDir::HORIZONTAL) {
      // Edge arrays have size (x_grids-1) in X for horizontal edges.
      if (gx + 1 < x_grids) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy, gx + 1, gy, layer));
      }
      if (gx - 1 >= 0) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx - 1, gy, gx, gy, layer));
      }
    } else if (dir == odb::dbTechLayerDir::VERTICAL) {
      if (gy + 1 < y_grids) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy, gx, gy + 1, layer));
      }
      if (gy - 1 >= 0) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy - 1, gx, gy, layer));
      }
    }

    if (spare == std::numeric_limits<int>::max()) {
      return 0;
    }
    return spare;
  };

  auto min_spare_on_segment_dir
      = [&](const int layer,
            const int gx,
            const int gy,
            const bool horizontal) -> int {
    if (layer < min_routing_layer || layer > max_routing_layer) {
      return 0;
    }
    if (!in_grid(gx, gy)) {
      return 0;
    }

    int spare = std::numeric_limits<int>::max();
    if (horizontal) {
      if (gx + 1 < x_grids) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy, gx + 1, gy, layer));
      }
      if (gx - 1 >= 0) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx - 1, gy, gx, gy, layer));
      }
    } else {
      if (gy + 1 < y_grids) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy, gx, gy + 1, layer));
      }
      if (gy - 1 >= 0) {
        spare = std::min(
            spare, fastroute->getAvailableResources(gx, gy - 1, gx, gy, layer));
      }
    }

    if (spare == std::numeric_limits<int>::max()) {
      return 0;
    }
    return spare;
  };

  CongestionPatchStats stats;
  int total_patches = 0;

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty() || total_patches >= kMaxTotalPatches) {
      continue;
    }

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(route.size() * 2 + 64);
    for (const auto& seg : route) {
      existing.insert(seg);
    }

    const std::vector<GSegment> original = route;

    int net_patches = 0;
    int net_last_patch_gx = std::numeric_limits<int>::min();
    int net_last_patch_gy = std::numeric_limits<int>::min();
    bool net_counted = false;

    for (const auto& seg : original) {
      if (total_patches >= kMaxTotalPatches
          || net_patches >= kMaxPatchesPerNet) {
        break;
      }

      if (seg.isVia() || seg.isJumper()) {
        continue;
      }
      if (seg.init_layer != seg.final_layer) {
        continue;
      }

      const int layer = seg.init_layer;
      if (layer < min_routing_layer || layer > max_patch_layer) {
        continue;
      }

      const int len_dbu = seg.length();
      if (len_dbu <= 0 || (len_dbu % tile_size) != 0) {
        continue;
      }

      const int len_tiles = len_dbu / tile_size;
      if (len_tiles < kMinSegmentLenTiles) {
        continue;
      }

      const bool horizontal
          = (seg.init_y == seg.final_y && seg.init_x != seg.final_x);
      const bool vertical
          = (seg.init_x == seg.final_x && seg.init_y != seg.final_y);
      if (!horizontal && !vertical) {
        continue;
      }

      const int margin = kEndMarginTiles * tile_size;
      const int x_start
          = std::min(seg.init_x, seg.final_x) + (horizontal ? margin : 0);
      const int x_end
          = std::max(seg.init_x, seg.final_x) - (horizontal ? margin : 0);
      const int y_start
          = std::min(seg.init_y, seg.final_y) + (vertical ? margin : 0);
      const int y_end
          = std::max(seg.init_y, seg.final_y) - (vertical ? margin : 0);

      if (horizontal && x_end <= x_start) {
        continue;
      }
      if (vertical && y_end <= y_start) {
        continue;
      }

      auto [gx0, gy0] = to_grid(horizontal ? x_start : seg.init_x,
                                vertical ? y_start : seg.init_y);
      auto [gx1, gy1] = to_grid(horizontal ? x_end : seg.final_x,
                                vertical ? y_end : seg.final_y);

      if (!in_grid(gx0, gy0) || !in_grid(gx1, gy1)) {
        continue;
      }

      if (horizontal && gy0 != gy1) {
        continue;
      }
      if (vertical && gx0 != gx1) {
        continue;
      }

      if (horizontal && gx1 < gx0) {
        std::swap(gx0, gx1);
      }
      if (vertical && gy1 < gy0) {
        std::swap(gy0, gy1);
      }

      int last_patch_along = std::numeric_limits<int>::min();
      const int begin = horizontal ? gx0 : gy0;
      const int end = horizontal ? gx1 : gy1;
      const int fixed = horizontal ? gy0 : gx0;

      for (int c = begin; c <= end; c += kSampleStepTiles) {
        if (total_patches >= kMaxTotalPatches
            || net_patches >= kMaxPatchesPerNet) {
          break;
        }

        if (c - last_patch_along < kMinPatchSpacingTiles) {
          continue;
        }

        const int gx = horizontal ? c : fixed;
        const int gy = horizontal ? fixed : c;

        const int spare = min_spare_on_segment_dir(layer, gx, gy, horizontal);
        if (spare >= kWirePatchSpareThreshold) {
          continue;
        }

        // Also apply a per-net spacing constraint in 2D to avoid clustering.
        if (std::abs(gx - net_last_patch_gx) + std::abs(gy - net_last_patch_gy)
            < kMinPatchSpacingTiles) {
          continue;
        }

        bool added = false;
        for (int adj = layer - 1; adj <= layer + 1; adj += 2) {
          if (adj < min_routing_layer || adj > max_routing_layer) {
            continue;
          }
          if (min_spare_at_point_preferred_dir(adj, gx, gy)
              < kAdjLayerMinSpare) {
            continue;
          }

          const auto [x_dbu, y_dbu] = to_dbu(gx, gy);
          const int via_lo = std::min(layer, adj);
          const int via_hi = std::max(layer, adj);
          const GSegment escape_via(x_dbu, y_dbu, via_lo, x_dbu, y_dbu, via_hi);
          if (existing.insert(escape_via).second) {
            route.push_back(escape_via);
            stats.points_added++;
            added = true;
          }
        }

        if (added) {
          last_patch_along = c;
          net_last_patch_gx = gx;
          net_last_patch_gy = gy;
          net_patches++;
          total_patches++;
          if (!net_counted) {
            stats.nets_patched++;
            net_counted = true;
          }
        }
      }
    }
  }

  if (stats.points_added > 0) {
    logger->info(GNR,
                 6009,
                 "NEWGR congestion patches: nets {}, added escape vias {}",
                 stats.nets_patched,
                 stats.points_added);
  }
}

// Add small "pin-access" escape vias at *existing* leaf endpoints when local
// spare resources are low.
//
// Motivation: detailed routing often introduces wirelength detours (and extra
// vias) near macros/pads due to restricted pin access. Adding a few extra guide
// points near those pins gives the detailed router more flexibility to connect
// without forcing a global-route detour.
//
// This is intentionally conservative to avoid bloating guides across the whole
// design: it only targets non-core pins and only triggers when local spare
// resources are below a threshold.
void add_endpoint_access_patches(NetRouteMap& routes,
                            FastRouteCore* fastroute,
                            const Grid* grid,
                            odb::dbTech* tech,
                            utl::Logger* logger,
                            const int min_routing_layer,
                            const int max_routing_layer)
{
  if (routes.empty() || fastroute == nullptr || grid == nullptr
      || tech == nullptr || logger == nullptr) {
    return;
  }

  const int tile_size = grid->getTileSize();
  if (tile_size <= 0) {
    return;
  }

  // Guardrails: keep patching bounded and focused.
  constexpr int kMaxTotalVias = 18000;
  constexpr int kMaxViasPerNet = 36;
  constexpr int kMaxEndpointsPerNet = 8;
  // IMPORTANT: keep patches connected to the existing guide graph.
  // Standalone guide islands can break incremental routing/parasitics updates
  // (net connectivity checks). Therefore we only add escape vias at the
  // existing endpoint coordinate (no 3x3 halo).
  constexpr int kPatchRadiusTiles = 0;

  // Threshold is expressed in "spare tracks" on the edge.
  constexpr int kPinSpareThreshold = 1;

  // Focus on the lowest layers where pin-access is usually the bottleneck.
  const int max_patch_layer
      = std::min(max_routing_layer, min_routing_layer + 1);

  const int grid_xmin = grid->getXMin();
  const int grid_ymin = grid->getYMin();
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();

  auto to_grid = [&](const int x_dbu, const int y_dbu) -> std::pair<int, int> {
    const int gx = (x_dbu - grid_xmin) / tile_size;
    const int gy = (y_dbu - grid_ymin) / tile_size;
    return {gx, gy};
  };

  auto to_dbu = [&](const int gx, const int gy) -> std::pair<int, int> {
    const int x_dbu = grid_xmin + gx * tile_size + tile_size / 2;
    const int y_dbu = grid_ymin + gy * tile_size + tile_size / 2;
    return {x_dbu, y_dbu};
  };

  auto in_grid = [&](const int gx, const int gy) -> bool {
    return gx >= 0 && gx < x_grids && gy >= 0 && gy < y_grids;
  };

  auto min_spare_around_cell = [&](const int layer,
                                   const int gx,
                                   const int gy) -> int {
    if (layer < min_routing_layer || layer > max_routing_layer) {
      return 0;
    }
    if (!in_grid(gx, gy)) {
      return 0;
    }

    int spare = std::numeric_limits<int>::max();
    if (gx + 1 < x_grids) {
      spare = std::min(
          spare, fastroute->getAvailableResources(gx, gy, gx + 1, gy, layer));
    }
    if (gx - 1 >= 0) {
      spare = std::min(
          spare, fastroute->getAvailableResources(gx - 1, gy, gx, gy, layer));
    }
    if (gy + 1 < y_grids) {
      spare = std::min(
          spare, fastroute->getAvailableResources(gx, gy, gx, gy + 1, layer));
    }
    if (gy - 1 >= 0) {
      spare = std::min(
          spare, fastroute->getAvailableResources(gx, gy - 1, gx, gy, layer));
    }

    if (spare == std::numeric_limits<int>::max()) {
      return 0;
    }
    return spare;
  };

  struct NodeKey
  {
    int x = 0;
    int y = 0;
    int layer = 0;

    bool operator==(const NodeKey& other) const
    {
      return x == other.x && y == other.y && layer == other.layer;
    }
  };

  struct NodeKeyHash
  {
    size_t operator()(const NodeKey& k) const noexcept
    {
      size_t h = 0;
      // Cheap combination; coordinates are DBU values and routing levels.
      h ^= std::hash<int>{}(k.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<int>{}(k.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<int>{}(k.layer) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  EndpointAccessPatchStats stats;
  int total_vias = 0;

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty() || total_vias >= kMaxTotalVias) {
      continue;
    }

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(route.size() * 2 + 64);
    for (const auto& seg : route) {
      existing.insert(seg);
    }

    int net_vias = 0;
    int net_endpoints = 0;
    bool net_counted = false;

    // Identify leaf endpoints (degree==1) in the guide graph. These are a good
    // proxy for pins without depending on Net/Pin object lifetimes.
    std::unordered_map<NodeKey, int, NodeKeyHash> degree;
    degree.reserve(route.size() * 2 + 64);

    for (const auto& seg : route) {
      if (seg.length() == 0) {
        continue;
      }

      const NodeKey a{seg.init_x, seg.init_y, seg.init_layer};
      const NodeKey b{seg.final_x, seg.final_y, seg.final_layer};
      degree[a]++;
      degree[b]++;
    }

    // Iterate over leaf nodes and add small patches near them on constrained
    // regions of the lowest layers.
    for (const auto& [node, deg] : degree) {
      if (total_vias >= kMaxTotalVias || net_vias >= kMaxViasPerNet
          || net_endpoints >= kMaxEndpointsPerNet) {
        break;
      }

      if (deg != 1) {
        continue;
      }
      if (node.layer < min_routing_layer || node.layer > max_patch_layer) {
        continue;
      }

      auto [gx, gy] = to_grid(node.x, node.y);
      if (!in_grid(gx, gy)) {
        continue;
      }

      // Only patch if this endpoint's neighborhood is locally constrained.
      if (min_spare_around_cell(node.layer, gx, gy) > kPinSpareThreshold) {
        continue;
      }

      stats.endpoints_patched++;
      net_endpoints++;

      for (int dl : {-1, 1}) {
        const int adj = node.layer + dl;
        if (adj < min_routing_layer || adj > max_routing_layer) {
          continue;
        }
        if (tech->findRoutingLayer(adj) == nullptr) {
          continue;
        }

        const int via_lo = std::min(node.layer, adj);
        const int via_hi = std::max(node.layer, adj);

        for (int dx = -kPatchRadiusTiles; dx <= kPatchRadiusTiles; dx++) {
          for (int dy = -kPatchRadiusTiles; dy <= kPatchRadiusTiles; dy++) {
            if (total_vias >= kMaxTotalVias || net_vias >= kMaxViasPerNet) {
              break;
            }

            const int px = gx + dx;
            const int py = gy + dy;
            if (!in_grid(px, py)) {
              continue;
            }

            const auto [x_dbu, y_dbu] = to_dbu(px, py);
            const GSegment escape_via(x_dbu, y_dbu, via_lo, x_dbu, y_dbu, via_hi);
            if (!existing.insert(escape_via).second) {
              continue;
            }

            route.push_back(escape_via);
            stats.vias_added++;
            total_vias++;
            net_vias++;

            if (!net_counted) {
              stats.nets_patched++;
              net_counted = true;
            }
          }
        }
      }
    }
  }

  if (stats.vias_added > 0) {
    logger->info(GNR,
                 6010,
                 "NEWGR endpoint-access patches: nets {}, endpoints {}, added escape vias {}",
                 stats.nets_patched,
                 stats.endpoints_patched,
                 stats.vias_added);
  }
}

// Add small, optional "lane" patches (via-up, short parallel segment, via-down)
// along long guides. This intentionally increases guide flexibility for the
// detailed router without forcing a longer global route topology.
//
// Inspired by CUGR "patching" concepts: provide extra escape points / track
// switching opportunities for long segments, especially on lower layers where
// pin-access constraints are tighter.
void add_lane_patches(NetRouteMap& routes,
                      GlobalRouter* grouter,
                      utl::Logger* logger,
                      const int min_routing_layer,
                      const int max_routing_layer)
{
  if (grouter == nullptr || logger == nullptr) {
    return;
  }

  // Keep patches opt-in to avoid perturbing baseline metrics. Enable with:
  //   `utl::set_debug_level(GNR, "newgrGuidePatches", 1)` (or equivalent)
  // in the calling environment.
  if (!logger->debugCheck(GNR, "newgrGuidePatches", 1)) {
    return;
  }

  const int tile_size = grouter->getTileSize();
  if (tile_size <= 0) {
    return;
  }

  constexpr int kMinSegmentLenTiles = 8;   // only long segments
  constexpr int kEndMarginTiles = 2;       // keep away from endpoints/pins
  constexpr int kMinLaneSpanTiles = 3;     // ensure the lane has meaningful span
  constexpr int kMaxPatchesPerNet = 6;     // runtime guardrail
  constexpr int kMaxTotalPatchedSegs = 25000;

  // Focus patches on the lower routing layers where DR often struggles most.
  const int max_patch_layer
      = std::min(max_routing_layer, min_routing_layer + 2);

  GuidePatchStats stats;
  int total_patched = 0;

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty()) {
      continue;
    }
    if (total_patched >= kMaxTotalPatchedSegs) {
      break;
    }

    // Iterate over a snapshot to avoid patching patches.
    const std::vector<GSegment> original = route;

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(route.size() * 2 + 32);
    for (const auto& seg : route) {
      existing.insert(seg);
    }

    int net_patches = 0;
    bool net_counted = false;

    for (const auto& seg : original) {
      if (total_patched >= kMaxTotalPatchedSegs
          || net_patches >= kMaxPatchesPerNet) {
        break;
      }

      if (seg.isVia() || seg.isJumper()) {
        continue;
      }
      if (seg.init_layer != seg.final_layer) {
        continue;
      }

      const int layer = seg.init_layer;
      if (layer < min_routing_layer || layer > max_patch_layer) {
        continue;
      }

      const int adj_layer
          = pick_adjacent_layer(layer, min_routing_layer, max_routing_layer);
      if (adj_layer < 0 || adj_layer == layer) {
        continue;
      }

      const int len_dbu = seg.length();
      if (len_dbu <= 0 || len_dbu % tile_size != 0) {
        continue;  // off-grid or degenerate
      }

      const int len_tiles = len_dbu / tile_size;
      if (len_tiles < kMinSegmentLenTiles) {
        continue;
      }

      const bool horizontal
          = (seg.init_y == seg.final_y && seg.init_x != seg.final_x);
      const bool vertical
          = (seg.init_x == seg.final_x && seg.init_y != seg.final_y);
      if (!horizontal && !vertical) {
        continue;
      }

      const int margin = kEndMarginTiles * tile_size;
      if (len_dbu <= (2 * margin + kMinLaneSpanTiles * tile_size)) {
        continue;
      }

      int x1 = seg.init_x;
      int y1 = seg.init_y;
      int x2 = seg.final_x;
      int y2 = seg.final_y;
      if (horizontal) {
        x1 = seg.init_x + margin;
        x2 = seg.final_x - margin;
      } else {
        y1 = seg.init_y + margin;
        y2 = seg.final_y - margin;
      }

      const int span_dbu = std::abs(x1 - x2) + std::abs(y1 - y2);
      if (span_dbu <= 0 || span_dbu % tile_size != 0) {
        continue;
      }
      const int span_tiles = span_dbu / tile_size;
      if (span_tiles < kMinLaneSpanTiles) {
        continue;
      }

      const int via_lo = std::min(layer, adj_layer);
      const int via_hi = std::max(layer, adj_layer);

      // Create a parallel lane on the adjacent layer spanning between two
      // internal points, with vias to allow DR to hop up/down.
      const GSegment via_a(x1, y1, via_lo, x1, y1, via_hi);
      const GSegment lane(x1, y1, adj_layer, x2, y2, adj_layer);
      const GSegment via_b(x2, y2, via_lo, x2, y2, via_hi);

      int added_here = 0;
      if (existing.insert(via_a).second) {
        route.push_back(via_a);
        stats.vias_added++;
        added_here++;
      }
      if (existing.insert(lane).second) {
        route.push_back(lane);
        stats.wires_added++;
        added_here++;
      }
      if (existing.insert(via_b).second) {
        route.push_back(via_b);
        stats.vias_added++;
        added_here++;
      }

      if (added_here > 0) {
        net_patches++;
        total_patched++;
        if (!net_counted) {
          stats.nets_patched++;
          net_counted = true;
        }
      }
    }
  }

  if (stats.vias_added > 0 || stats.wires_added > 0) {
    logger->info(GNR,
                 6008,
                 "NEWGR guide patches: nets {}, added vias {}, added segs {}",
                 stats.nets_patched,
                 stats.vias_added,
                 stats.wires_added);
  }
}

struct RouteMetrics
{
  long wirelength_dbu = 0;
  long via_count = 0;
  double wirelength_um = 0.0;
};

struct CandidateSettings
{
  std::string name;
  int seed = 0;
  float caps_perturbation_percentage = 0.0f;
  int perturbation_amount = 1;
  int congestion_iterations = 50;
  float critical_nets_percentage = 10.0f;
  bool allow_congestion = false;
};

struct CandidateResult
{
  CandidateSettings settings;
  RouteMetrics metrics;
  int overflow = std::numeric_limits<int>::max();
  NetRouteMap routes;
};

struct GuideUnionStats
{
  int nets_checked = 0;
  int nets_merged = 0;
  int nets_replaced = 0;
  int nets_reverted_disconnected = 0;
  int nets_skipped_small = 0;
  int nets_skipped_no_improvement = 0;
  int nets_skipped_over_budget = 0;
  int segments_added = 0;
  int wires_added = 0;
  int vias_added = 0;
};

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  int congestion_iterations = 50;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
};

struct GridPoint
{
  int gx = 0;
  int gy = 0;
  bool operator==(const GridPoint& other) const
  {
    return gx == other.gx && gy == other.gy;
  }
};

struct GridPointHash
{
  std::size_t operator()(const GridPoint& p) const
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, p.gx);
    boost::hash_combine(seed, p.gy);
    return seed;
  }
};

int pick_routing_layer_by_dir(odb::dbTech* tech,
                              const int min_routing_layer,
                              const int max_routing_layer,
                              const odb::dbTechLayerDir::Value dir)
{
  if (tech == nullptr) {
    return min_routing_layer;
  }
  for (int layer = min_routing_layer; layer <= max_routing_layer; layer++) {
    odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer);
    if (tech_layer == nullptr) {
      continue;
    }
    if (tech_layer->getDirection() == dir) {
      return layer;
    }
  }
  return min_routing_layer;
}

struct LShapeScore
{
  int blocked_edges = 0;
  int min_available = std::numeric_limits<int>::max();
};

LShapeScore score_lshape(const GridPoint& a,
                         const GridPoint& b,
                         const GridPoint& turn,
                         const int h_layer,
                         const int v_layer,
                         odb::dbNet* db_net,
                         GlobalRouter* grouter)
{
  LShapeScore score;
  FastRouteCore* fastroute = (grouter != nullptr) ? grouter->fastroute() : nullptr;
  if (fastroute == nullptr || db_net == nullptr) {
    return score;
  }

  auto eval_segment = [&](const GridPoint& p0,
                          const GridPoint& p1,
                          const int layer) {
    const int req = fastroute->getDbNetLayerEdgeCost(db_net, layer);
    if (p0.gx == p1.gx && p0.gy == p1.gy) {
      return;
    }

    if (p0.gy == p1.gy) {
      const int y = p0.gy;
      const int x0 = std::min(p0.gx, p1.gx);
      const int x1 = std::max(p0.gx, p1.gx);
      for (int x = x0; x < x1; x++) {
        const int avail
            = fastroute->getAvailableResources(x, y, x + 1, y, layer);
        score.min_available = std::min(score.min_available, avail);
        if (avail < req) {
          score.blocked_edges++;
        }
      }
      return;
    }

    if (p0.gx == p1.gx) {
      const int x = p0.gx;
      const int y0 = std::min(p0.gy, p1.gy);
      const int y1 = std::max(p0.gy, p1.gy);
      for (int y = y0; y < y1; y++) {
        const int avail
            = fastroute->getAvailableResources(x, y, x, y + 1, layer);
        score.min_available = std::min(score.min_available, avail);
        if (avail < req) {
          score.blocked_edges++;
        }
      }
      return;
    }
  };

  // Segment a->turn
  if (a.gy == turn.gy) {
    eval_segment(a, turn, h_layer);
  } else if (a.gx == turn.gx) {
    eval_segment(a, turn, v_layer);
  }

  // Segment turn->b
  if (turn.gy == b.gy) {
    eval_segment(turn, b, h_layer);
  } else if (turn.gx == b.gx) {
    eval_segment(turn, b, v_layer);
  }

  return score;
}

std::uint64_t encode_unit_edge_key(int gx0,
                                   int gy0,
                                   int gx1,
                                   int gy1,
                                   int layer)
{
  if (gx0 > gx1 || (gx0 == gx1 && gy0 > gy1)) {
    std::swap(gx0, gx1);
    std::swap(gy0, gy1);
  }

  constexpr std::uint64_t kCoordMask = (1u << 14) - 1u;
  const std::uint64_t l = static_cast<std::uint64_t>(layer & 0xff);
  const std::uint64_t x0 = static_cast<std::uint64_t>(gx0) & kCoordMask;
  const std::uint64_t y0 = static_cast<std::uint64_t>(gy0) & kCoordMask;
  const std::uint64_t x1 = static_cast<std::uint64_t>(gx1) & kCoordMask;
  const std::uint64_t y1 = static_cast<std::uint64_t>(gy1) & kCoordMask;

  return (l << 56) | (x0 << 42) | (y0 << 28) | (x1 << 14) | y1;
}

std::unordered_set<std::uint64_t> collect_unit_grid_edges(const GRoute& route,
                                                          const Grid* grid)
{
  std::unordered_set<std::uint64_t> edges;
  if (grid == nullptr) {
    return edges;
  }

  const int tile = grid->getTileSize();
  if (tile <= 0) {
    return edges;
  }

  const int x0_corner = grid->getXMin();
  const int y0_corner = grid->getYMin();
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return edges;
  }

  auto to_grid_idx = [&](const int x_dbu, const int y_dbu) -> std::pair<int, int> {
    int gx = (x_dbu - x0_corner) / tile;
    int gy = (y_dbu - y0_corner) / tile;
    gx = std::max(0, std::min(x_grids - 1, gx));
    gy = std::max(0, std::min(y_grids - 1, gy));
    return {gx, gy};
  };

  edges.reserve(route.size() * 6 + 64);
  for (const auto& seg : route) {
    if (seg.isVia()) {
      continue;
    }
    const int layer = seg.init_layer;

    auto [gx0, gy0] = to_grid_idx(seg.init_x, seg.init_y);
    auto [gx1, gy1] = to_grid_idx(seg.final_x, seg.final_y);
    if (gx0 == gx1 && gy0 == gy1) {
      continue;
    }

    if (gy0 == gy1) {
      const int y = gy0;
      const int x_min = std::min(gx0, gx1);
      const int x_max = std::max(gx0, gx1);
      for (int x = x_min; x < x_max; x++) {
        edges.insert(encode_unit_edge_key(x, y, x + 1, y, layer));
      }
    } else if (gx0 == gx1) {
      const int x = gx0;
      const int y_min = std::min(gy0, gy1);
      const int y_max = std::max(gy0, gy1);
      for (int y = y_min; y < y_max; y++) {
        edges.insert(encode_unit_edge_key(x, y, x, y + 1, layer));
      }
    }
  }

  return edges;
}

bool find_slack_shortest_path(const GridPoint& start,
                              const GridPoint& goal,
                              const int h_layer,
                              const int v_layer,
                              odb::dbNet* db_net,
                              GlobalRouter* grouter,
                              const std::unordered_set<std::uint64_t>& base_edges,
                              std::vector<GridPoint>& out_path)
{
  out_path.clear();
  if (grouter == nullptr || db_net == nullptr) {
    return false;
  }
  FastRouteCore* fastroute = grouter->fastroute();
  const Grid* grid = grouter->grid();
  if (fastroute == nullptr || grid == nullptr) {
    return false;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return false;
  }

  const int dx = std::abs(start.gx - goal.gx);
  const int dy = std::abs(start.gy - goal.gy);
  const int manhattan = dx + dy;

  const int margin = std::max(6, std::min(32, manhattan / 3));
  const int xmin = std::max(0, std::min(start.gx, goal.gx) - margin);
  const int xmax = std::min(x_grids - 1, std::max(start.gx, goal.gx) + margin);
  const int ymin = std::max(0, std::min(start.gy, goal.gy) - margin);
  const int ymax = std::min(y_grids - 1, std::max(start.gy, goal.gy) + margin);

  const int width = xmax - xmin + 1;
  const int height = ymax - ymin + 1;
  if (width <= 0 || height <= 0) {
    return false;
  }

  constexpr int kMaxStates = 160000;  // width*height*3
  if (width * height * 3 > kMaxStates) {
    return false;
  }

  enum Dir : uint8_t
  {
    kNone = 0,
    kH = 1,
    kV = 2
  };

  auto state_index = [&](const int gx, const int gy, const int dir) -> int {
    const int lx = gx - xmin;
    const int ly = gy - ymin;
    return ((ly * width + lx) * 3) + dir;
  };

  auto decode_state = [&](const int idx) -> std::tuple<int, int, int> {
    const int cell = idx / 3;
    const int dir = idx % 3;
    const int ly = cell / width;
    const int lx = cell % width;
    return {xmin + lx, ymin + ly, dir};
  };

  auto heuristic = [&](const int gx, const int gy) -> int {
    return std::abs(gx - goal.gx) + std::abs(gy - goal.gy);
  };

  const int req_h = fastroute->getDbNetLayerEdgeCost(db_net, h_layer);
  const int req_v = fastroute->getDbNetLayerEdgeCost(db_net, v_layer);

  // Radical WL strategy:
  // Allow the "directness-first" A* to traverse mildly-overused edges (avail < req)
  // with a steep penalty instead of hard-blocking them. This intentionally
  // injects shorter rectilinear corridors into the guide set even when the
  // global congestion DB is pessimistic; the detailed router can still choose
  // the original (detour) guides if needed.
  //
  // Keep true hard blocks (avail <= 0) forbidden to avoid routing through
  // blocked/absent resources (e.g., macro blockages).
  constexpr bool kAllowSoftOveruse = true;
  constexpr int kSoftOveruseBasePenalty = 40;
  constexpr int kSoftOverusePerTrackPenalty = 6;

  auto step_ok_and_penalty = [&](const int gx0,
                                 const int gy0,
                                 const int gx1,
                                 const int gy1,
                                 const int layer,
                                 const int req,
                                 int& out_penalty) -> bool {
    out_penalty = 0;
    if (req <= 0) {
      return true;
    }

    const int x0 = std::min(gx0, gx1);
    const int y0 = std::min(gy0, gy1);
    const int x1 = std::max(gx0, gx1);
    const int y1 = std::max(gy0, gy1);

    int avail = fastroute->getAvailableResources(x0, y0, x1, y1, layer);
    const std::uint64_t key = encode_unit_edge_key(gx0, gy0, gx1, gy1, layer);
    if (base_edges.find(key) != base_edges.end()) {
      // Credit this edge with the net's own demand, approximating a rip-up.
      avail += req;
    }

    if (avail <= 0) {
      return false;
    }

    if (avail < req) {
      if (!kAllowSoftOveruse) {
        return false;
      }
      const int deficit = req - avail;
      out_penalty += kSoftOveruseBasePenalty
                     + deficit * kSoftOverusePerTrackPenalty;
      return true;
    }

    const int slack = avail - req;
    constexpr int kSoftSlackThreshold = 1;
    if (slack < kSoftSlackThreshold) {
      out_penalty += (kSoftSlackThreshold - slack);
    }

    return true;
  };

  constexpr int kStepCost = 1;
  constexpr int kTurnPenalty = 2;

  struct Node
  {
    int f = 0;
    int g = 0;
    int idx = 0;
  };

  struct NodeCmp
  {
    bool operator()(const Node& a, const Node& b) const { return a.f > b.f; }
  };

  const int nstates = width * height * 3;
  constexpr int kInf = std::numeric_limits<int>::max() / 4;
  std::vector<int> gscore(nstates, kInf);
  std::vector<int> parent(nstates, -1);

  std::priority_queue<Node, std::vector<Node>, NodeCmp> pq;

  const int start_idx = state_index(start.gx, start.gy, kNone);
  gscore[start_idx] = 0;
  pq.push(Node{heuristic(start.gx, start.gy), 0, start_idx});

  int goal_idx = -1;

  while (!pq.empty()) {
    const Node cur = pq.top();
    pq.pop();

    if (cur.g != gscore[cur.idx]) {
      continue;
    }

    auto [gx, gy, dir] = decode_state(cur.idx);
    if (gx == goal.gx && gy == goal.gy) {
      goal_idx = cur.idx;
      break;
    }

    auto try_relax = [&](const int ngx, const int ngy, const int ndir) {
      if (ngx < xmin || ngx > xmax || ngy < ymin || ngy > ymax) {
        return;
      }

      const bool is_h = (ngy == gy) && (ngx != gx);
      const bool is_v = (ngx == gx) && (ngy != gy);
      if (!is_h && !is_v) {
        return;
      }

      const int layer = is_h ? h_layer : v_layer;
      const int req = is_h ? req_h : req_v;

      int slack_penalty = 0;
      if (!step_ok_and_penalty(gx, gy, ngx, ngy, layer, req, slack_penalty)) {
        return;
      }

      int step_cost = kStepCost + slack_penalty;
      if (dir != kNone && dir != ndir) {
        step_cost += kTurnPenalty;
      }

      const int next_idx = state_index(ngx, ngy, ndir);
      const int tentative_g = cur.g + step_cost;
      if (tentative_g < gscore[next_idx]) {
        gscore[next_idx] = tentative_g;
        parent[next_idx] = cur.idx;
        const int f = tentative_g + heuristic(ngx, ngy);
        pq.push(Node{f, tentative_g, next_idx});
      }
    };

    try_relax(gx + 1, gy, kH);
    try_relax(gx - 1, gy, kH);
    try_relax(gx, gy + 1, kV);
    try_relax(gx, gy - 1, kV);
  }

  if (goal_idx < 0) {
    return false;
  }

  std::vector<GridPoint> rev;
  rev.reserve(static_cast<std::size_t>(manhattan) + 8);
  int cur = goal_idx;
  while (cur >= 0) {
    auto [gx, gy, _dir] = decode_state(cur);
    rev.push_back(GridPoint{gx, gy});
    cur = parent[cur];
  }

  if (rev.empty()) {
    return false;
  }
  out_path.assign(rev.rbegin(), rev.rend());
  return out_path.size() >= 1;
}

void add_wire(GRoute& route,
              const odb::Point& p0,
              const odb::Point& p1,
              const int layer)
{
  if (p0 == p1) {
    return;
  }
  route.emplace_back(p0.x(), p0.y(), layer, p1.x(), p1.y(), layer);
}

void add_via(GRoute& route,
             const odb::Point& p,
             const int layer0,
             const int layer1)
{
  if (layer0 == layer1) {
    return;
  }
  route.emplace_back(p.x(), p.y(), layer0, p.x(), p.y(), layer1);
}

RouteMetrics compute_route_metrics(const GRoute& segments, odb::dbTech* tech);

// Add "rectilinear corridor" guides based on an RSMT topology (FLUTE when
// available), but *without* FastRoute's detours. These extra corridors are
// injected only for nets whose base guides appear significantly detoured vs
// their RSMT length.
//
// Rationale (radical WL strategy): detailed routing wirelength can increase
// when the guide set only contains congestion-safe detours. By injecting the
// shorter rectilinear corridors as *additional* guides, we give the detailed
// router a chance to realize a shorter path when tracks/vias allow it, while
// still keeping the overflow-free detour guides present as a fallback.
void add_rsmt_rectilinear_corridors(NetRouteMap& routes,
                                    GlobalRouter* grouter,
                                    stt::SteinerTreeBuilder* stt_builder,
                                    odb::dbTech* tech,
                                    utl::Logger* logger,
                                    const int min_routing_layer,
                                    const int max_routing_layer)
{
  if (routes.empty() || grouter == nullptr || tech == nullptr || logger == nullptr
      || stt_builder == nullptr) {
    return;
  }

  const Grid* grid = grouter->grid();
  if (grid == nullptr) {
    return;
  }

  const int tile_size = grid->getTileSize();
  if (tile_size <= 0) {
    return;
  }

  const int x_corner = grid->getXMin();
  const int y_corner = grid->getYMin();
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return;
  }

  const int h_layer = pick_routing_layer_by_dir(
      tech,
      min_routing_layer,
      max_routing_layer,
      odb::dbTechLayerDir::HORIZONTAL);
  const int v_layer = pick_routing_layer_by_dir(
      tech, min_routing_layer, max_routing_layer, odb::dbTechLayerDir::VERTICAL);

  // Guardrails.
  constexpr int kMaxTotalAddedSegments = 220000;
  constexpr int kMaxAddedSegmentsPerNet = 220;
  const long min_net_wl_dbu = static_cast<long>(tile_size) * 25;

  // Only target "detoured" nets: base_wl must exceed RSMT length by a ratio.
  // Lower threshold to escape WL local minima: add RSMT corridors for more nets
  // so DR has more opportunities to realize a shorter path.
  constexpr double kDetourRatioThreshold = 1.06;

  // Keep the corridors "rectilinear": never detour outside the Manhattan
  // rectangle for an RSMT edge; optionally add both L-shapes when both look
  // similarly routable (more options for DR).
  constexpr int kMaxBlockedEdgesForAlt = 4;

  auto to_grid = [&](const odb::Point& p) -> GridPoint {
    GridPoint gp;
    gp.gx = (p.x() - x_corner) / tile_size;
    gp.gy = (p.y() - y_corner) / tile_size;
    gp.gx = std::max(0, std::min(x_grids - 1, gp.gx));
    gp.gy = std::max(0, std::min(y_grids - 1, gp.gy));
    return gp;
  };

  auto to_dbu = [&](const GridPoint& gp) -> odb::Point {
    const int x_dbu = x_corner + gp.gx * tile_size + tile_size / 2;
    const int y_dbu = y_corner + gp.gy * tile_size + tile_size / 2;
    return odb::Point(x_dbu, y_dbu);
  };

  auto is_better_score = [](const LShapeScore& a, const LShapeScore& b) -> bool {
    if (a.blocked_edges != b.blocked_edges) {
      return a.blocked_edges < b.blocked_edges;
    }
    if (a.min_available != b.min_available) {
      return a.min_available > b.min_available;
    }
    return false;
  };

  int total_added = 0;
  int nets_augmented = 0;

  for (auto& [db_net, route] : routes) {
    if (total_added >= kMaxTotalAddedSegments) {
      break;
    }
    if (db_net == nullptr || route.empty()) {
      continue;
    }

    const RouteMetrics base_m = compute_route_metrics(route, tech);
    if (base_m.wirelength_dbu < min_net_wl_dbu) {
      continue;
    }

    std::vector<PinGridLocation> pin_locs = grouter->getPinGridPositions(db_net);
    if (pin_locs.size() < 2) {
      continue;
    }

    std::vector<GridPoint> pins_grid;
    pins_grid.reserve(pin_locs.size());
    std::unordered_set<GridPoint, GridPointHash> seen;
    seen.reserve(pin_locs.size() * 2);
    for (const auto& loc : pin_locs) {
      const GridPoint gp = to_grid(loc.grid_pt);
      if (seen.insert(gp).second) {
        pins_grid.push_back(gp);
      }
    }
    if (pins_grid.size() < 2) {
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(pins_grid.size());
    ys.reserve(pins_grid.size());
    for (const auto& gp : pins_grid) {
      xs.push_back(gp.gx);
      ys.push_back(gp.gy);
    }

    constexpr int kFluteAccuracy = 3;
    const stt::Tree rsmt = stt_builder->flute(xs, ys, kFluteAccuracy);
    if (rsmt.branchCount() <= 0) {
      continue;
    }

    // Estimate RSMT length (grid-manhttan units), then translate to DBU.
    long rsmt_len_tiles = 0;
    {
      std::unordered_set<std::uint64_t> seen_edges;
      seen_edges.reserve(static_cast<std::size_t>(rsmt.branchCount()) * 2);

      for (int i = 0; i < rsmt.branchCount(); i++) {
        const int j = rsmt.branch[i].n;
        if (j < 0 || j >= rsmt.branchCount() || j == i) {
          continue;
        }
        int a = i;
        int b = j;
        if (a > b) {
          std::swap(a, b);
        }
        const std::uint64_t key
            = (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
        if (!seen_edges.insert(key).second) {
          continue;
        }
        rsmt_len_tiles += std::abs(rsmt.branch[a].x - rsmt.branch[b].x)
                          + std::abs(rsmt.branch[a].y - rsmt.branch[b].y);
      }
    }
    const long rsmt_len_dbu = rsmt_len_tiles * static_cast<long>(tile_size);
    if (rsmt_len_dbu <= 0) {
      continue;
    }

    const double detour_ratio
        = static_cast<double>(base_m.wirelength_dbu)
          / static_cast<double>(std::max<long>(1, rsmt_len_dbu));
    if (detour_ratio < kDetourRatioThreshold) {
      continue;
    }

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(route.size() * 2 + 256);
    for (const auto& seg : route) {
      existing.insert(seg);
    }

    int added_for_net = 0;
    bool net_counted = false;

    auto add_wire_unique = [&](const odb::Point& p0,
                               const odb::Point& p1,
                               const int layer) {
      if (p0 == p1) {
        return;
      }
      const GSegment seg(p0.x(), p0.y(), layer, p1.x(), p1.y(), layer);
      if (existing.insert(seg).second) {
        route.push_back(seg);
        added_for_net++;
        total_added++;
      }
    };

    auto add_via_unique = [&](const odb::Point& p,
                              const int layer0,
                              const int layer1) {
      if (layer0 == layer1) {
        return;
      }
      const int via_lo = std::min(layer0, layer1);
      const int via_hi = std::max(layer0, layer1);
      const GSegment via(p.x(), p.y(), via_lo, p.x(), p.y(), via_hi);
      if (existing.insert(via).second) {
        route.push_back(via);
        added_for_net++;
        total_added++;
      }
    };

    std::unordered_set<std::uint64_t> emitted_edges;
    emitted_edges.reserve(static_cast<std::size_t>(rsmt.branchCount()) * 2);

    auto emit_lshape = [&](const GridPoint& a, const GridPoint& b, const GridPoint& t) {
      const odb::Point p_a = to_dbu(a);
      const odb::Point p_b = to_dbu(b);
      const odb::Point p_t = to_dbu(t);

      if (a.gx == t.gx) {
        add_wire_unique(p_a, p_t, v_layer);
      } else if (a.gy == t.gy) {
        add_wire_unique(p_a, p_t, h_layer);
      }
      if (t.gx == b.gx) {
        add_wire_unique(p_t, p_b, v_layer);
      } else if (t.gy == b.gy) {
        add_wire_unique(p_t, p_b, h_layer);
      }
      if (h_layer != v_layer && a.gx != b.gx && a.gy != b.gy) {
        add_via_unique(p_t, h_layer, v_layer);
      }
    };

    for (int i = 0; i < rsmt.branchCount(); i++) {
      if (total_added >= kMaxTotalAddedSegments
          || added_for_net >= kMaxAddedSegmentsPerNet) {
        break;
      }

      const int j = rsmt.branch[i].n;
      if (j < 0 || j >= rsmt.branchCount() || j == i) {
        continue;
      }
      int a = i;
      int b = j;
      if (a > b) {
        std::swap(a, b);
      }
      const std::uint64_t key
          = (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
      if (!emitted_edges.insert(key).second) {
        continue;
      }

      const GridPoint p0{std::max(0, std::min(x_grids - 1, rsmt.branch[a].x)),
                         std::max(0, std::min(y_grids - 1, rsmt.branch[a].y))};
      const GridPoint p1{std::max(0, std::min(x_grids - 1, rsmt.branch[b].x)),
                         std::max(0, std::min(y_grids - 1, rsmt.branch[b].y))};
      if (p0 == p1) {
        continue;
      }

      const GridPoint turn1{p0.gx, p1.gy};
      const GridPoint turn2{p1.gx, p0.gy};
      const LShapeScore s1
          = score_lshape(p0, p1, turn1, h_layer, v_layer, db_net, grouter);
      const LShapeScore s2
          = score_lshape(p0, p1, turn2, h_layer, v_layer, db_net, grouter);

      const GridPoint best = is_better_score(s1, s2) ? turn1 : turn2;
      emit_lshape(p0, p1, best);

      // Optional: also emit the alternative L-shape when both are "almost
      // unblocked", giving DR a second rectilinear corridor without full maze
      // detours. This is intentionally aggressive (radical exploration) but
      // bounded by strict budgets above.
      if (s1.blocked_edges <= kMaxBlockedEdgesForAlt
          && s2.blocked_edges <= kMaxBlockedEdgesForAlt && !(turn1 == turn2)
          && total_added < kMaxTotalAddedSegments
          && added_for_net + 3 < kMaxAddedSegmentsPerNet) {
        const GridPoint alt = (best == turn1) ? turn2 : turn1;
        emit_lshape(p0, p1, alt);
      }

      if (!net_counted && added_for_net > 0) {
        nets_augmented++;
        net_counted = true;
      }
    }
  }

  if (nets_augmented > 0) {
    logger->info(GNR,
                 6014,
                 "NEWGR RSMT rect-corridors: nets {}, added segs {}",
                 nets_augmented,
                 total_added);
  } else {
    logger->info(GNR, 6015, "NEWGR RSMT rect-corridors: no nets met detour criteria");
  }
}

GRoute build_pattern_direct_route(odb::dbNet* db_net,
                                  GlobalRouter* grouter,
                                  stt::SteinerTreeBuilder* stt_builder,
                                  odb::dbTech* tech,
                                  const int min_routing_layer,
                                  const int max_routing_layer,
                                  const GRoute* base_route_for_slack)
{
  GRoute route;
  if (db_net == nullptr || grouter == nullptr) {
    return route;
  }

  const Grid* grid = grouter->grid();
  if (grid == nullptr) {
    return route;
  }

  const int tile_size = grid->getTileSize();
  if (tile_size <= 0) {
    return route;
  }

  const int x_corner = grid->getXMin();
  const int y_corner = grid->getYMin();

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return route;
  }

  auto to_grid = [&](const odb::Point& p) -> GridPoint {
    GridPoint gp;
    gp.gx = (p.x() - x_corner) / tile_size;
    gp.gy = (p.y() - y_corner) / tile_size;
    gp.gx = std::max(0, std::min(x_grids - 1, gp.gx));
    gp.gy = std::max(0, std::min(y_grids - 1, gp.gy));
    return gp;
  };

  auto to_dbu = [&](const GridPoint& gp) -> odb::Point {
    const int x_dbu = x_corner + gp.gx * tile_size + tile_size / 2;
    const int y_dbu = y_corner + gp.gy * tile_size + tile_size / 2;
    return odb::Point(x_dbu, y_dbu);
  };

  const int h_layer = pick_routing_layer_by_dir(
      tech, min_routing_layer, max_routing_layer, odb::dbTechLayerDir::HORIZONTAL);
  const int v_layer = pick_routing_layer_by_dir(
      tech, min_routing_layer, max_routing_layer, odb::dbTechLayerDir::VERTICAL);

  std::vector<PinGridLocation> pin_locs = grouter->getPinGridPositions(db_net);
  if (pin_locs.size() < 2) {
    return route;
  }

  std::vector<GridPoint> pins_grid;
  pins_grid.reserve(pin_locs.size());
  std::unordered_set<GridPoint, GridPointHash> seen;
  seen.reserve(pin_locs.size() * 2);
  for (const auto& loc : pin_locs) {
    const GridPoint gp = to_grid(loc.grid_pt);
    if (seen.insert(gp).second) {
      pins_grid.push_back(gp);
    }
  }
  if (pins_grid.size() < 2) {
    return route;
  }

  // Build a *direct* per-net candidate with a different topology generator
  // than FastRoute's congestion-driven tree:
  // - Prefer a FLUTE RSMT (in grid units) when available to reduce raw
  //   topology wirelength for multi-pin nets.
  // - Fall back to the previous MST topology if FLUTE isn't available.
  //
  // This is intentionally "radical": it can produce a very different guide
  // graph than the base routes, giving the detailed router a stronger chance
  // to realize shorter wiring.
  stt::Tree rsmt;
  bool use_flute = false;
  if (stt_builder != nullptr) {
    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(pins_grid.size());
    ys.reserve(pins_grid.size());
    for (const auto& gp : pins_grid) {
      xs.push_back(gp.gx);
      ys.push_back(gp.gy);
    }

    // Accuracy 3 is the OpenROAD default and is usually a good WL/CPU tradeoff.
    constexpr int kFluteAccuracy = 3;
    rsmt = stt_builder->flute(xs, ys, kFluteAccuracy);
    use_flute = rsmt.branchCount() > 0;
  }

  std::unordered_map<long long, uint8_t> endpoint_layer_mask;
  endpoint_layer_mask.reserve(static_cast<size_t>(pins_grid.size()) * 8);

  const std::unordered_set<std::uint64_t> base_edges
      = (base_route_for_slack != nullptr)
            ? collect_unit_grid_edges(*base_route_for_slack, grid)
            : std::unordered_set<std::uint64_t>{};

  std::unordered_set<GSegment, GSegmentHash> existing;
  existing.reserve(256);

  auto add_wire_unique = [&](const odb::Point& p0,
                             const odb::Point& p1,
                             const int layer) {
    if (p0 == p1) {
      return;
    }
    const GSegment seg(p0.x(), p0.y(), layer, p1.x(), p1.y(), layer);
    if (existing.insert(seg).second) {
      route.push_back(seg);
    }
  };

  auto add_via_unique = [&](const odb::Point& p,
                            const int layer0,
                            const int layer1) {
    if (layer0 == layer1) {
      return;
    }
    const int via_lo = std::min(layer0, layer1);
    const int via_hi = std::max(layer0, layer1);
    const GSegment via(p.x(), p.y(), via_lo, p.x(), p.y(), via_hi);
    if (existing.insert(via).second) {
      route.push_back(via);
    }
  };

  auto add_endpoint_layer = [&](const odb::Point& p, const int layer) {
    const long long key
        = (static_cast<long long>(p.x()) << 32) ^ static_cast<unsigned>(p.y());
    uint8_t mask = endpoint_layer_mask[key];
    if (layer == h_layer) {
      mask |= 0x1;
    }
    if (layer == v_layer) {
      mask |= 0x2;
    }
    endpoint_layer_mask[key] = mask;
  };

  auto is_better_score = [](const LShapeScore& a, const LShapeScore& b) -> bool {
    if (a.blocked_edges != b.blocked_edges) {
      return a.blocked_edges < b.blocked_edges;
    }
    if (a.min_available != b.min_available) {
      return a.min_available > b.min_available;
    }
    return false;
  };

  auto emit_edge = [&](const GridPoint& a, const GridPoint& b) {
    if (a == b) {
      return;
    }

    // Radical directness-first alternative: route with an A* that optimizes
    // for shortest length under a *hard* capacity-feasibility constraint on
    // the current congestion DB (crediting the base route for this net).
    //
    // If it fails (too tight / too large search), fall back to the old
    // L-shape chooser.
    // Raise this threshold to give more large nets a chance to get a direct,
    // length-driven corridor (even if some edges are soft-overused). This is
    // intentionally more expensive and may increase GR runtime.
    constexpr std::size_t kMaxPinsForSlackMaze = 60;
    if (pins_grid.size() <= kMaxPinsForSlackMaze) {
      std::vector<GridPoint> path;
      if (find_slack_shortest_path(
              a, b, h_layer, v_layer, db_net, grouter, base_edges, path)
          && path.size() >= 2) {
        auto step_dir = [&](const GridPoint& p0, const GridPoint& p1) -> int {
          if (p0.gy == p1.gy && p0.gx != p1.gx) {
            return 1;  // H
          }
          if (p0.gx == p1.gx && p0.gy != p1.gy) {
            return 2;  // V
          }
          return 0;
        };

        int cur_dir = step_dir(path[0], path[1]);
        if (cur_dir != 0) {
          int cur_layer = (cur_dir == 1) ? h_layer : v_layer;
          GridPoint seg_start = path[0];
          for (std::size_t i = 1; i < path.size(); i++) {
            const GridPoint prev = path[i - 1];
            const GridPoint cur = path[i];
            const int d = step_dir(prev, cur);
            if (d == 0) {
              break;
            }
            if (d != cur_dir) {
              const odb::Point p0 = to_dbu(seg_start);
              const odb::Point p1 = to_dbu(prev);
              add_wire_unique(p0, p1, cur_layer);
              add_endpoint_layer(p0, cur_layer);
              add_endpoint_layer(p1, cur_layer);

              const odb::Point p_turn = to_dbu(prev);
              add_via_unique(p_turn, h_layer, v_layer);

              seg_start = prev;
              cur_dir = d;
              cur_layer = (cur_dir == 1) ? h_layer : v_layer;
            }
          }

          const odb::Point p0 = to_dbu(seg_start);
          const odb::Point p1 = to_dbu(path.back());
          add_wire_unique(p0, p1, cur_layer);
          add_endpoint_layer(p0, cur_layer);
          add_endpoint_layer(p1, cur_layer);
          return;
        }
      }
    }

    const GridPoint turn1{a.gx, b.gy};
    const GridPoint turn2{b.gx, a.gy};
    const LShapeScore s1
        = score_lshape(a, b, turn1, h_layer, v_layer, db_net, grouter);
    const LShapeScore s2
        = score_lshape(a, b, turn2, h_layer, v_layer, db_net, grouter);
    const GridPoint turn = is_better_score(s1, s2) ? turn1 : turn2;

    const odb::Point p_a = to_dbu(a);
    const odb::Point p_b = to_dbu(b);
    const odb::Point p_t = to_dbu(turn);

    if (a.gy == turn.gy) {
      add_wire_unique(p_a, p_t, h_layer);
      add_endpoint_layer(p_a, h_layer);
      add_endpoint_layer(p_t, h_layer);
    } else if (a.gx == turn.gx) {
      add_wire_unique(p_a, p_t, v_layer);
      add_endpoint_layer(p_a, v_layer);
      add_endpoint_layer(p_t, v_layer);
    }

    if (turn.gy == b.gy) {
      add_wire_unique(p_t, p_b, h_layer);
      add_endpoint_layer(p_t, h_layer);
      add_endpoint_layer(p_b, h_layer);
    } else if (turn.gx == b.gx) {
      add_wire_unique(p_t, p_b, v_layer);
      add_endpoint_layer(p_t, v_layer);
      add_endpoint_layer(p_b, v_layer);
    }

    if (h_layer != v_layer && a.gx != b.gx && a.gy != b.gy) {
      add_via_unique(p_t, h_layer, v_layer);
    }
  };

  if (use_flute) {
    std::unordered_set<std::uint64_t> seen_edges;
    seen_edges.reserve(static_cast<std::size_t>(rsmt.branchCount()) * 2);

    for (int i = 0; i < rsmt.branchCount(); i++) {
      const int j = rsmt.branch[i].n;
      if (j < 0 || j >= rsmt.branchCount() || j == i) {
        continue;
      }

      int a = i;
      int b = j;
      if (a > b) {
        std::swap(a, b);
      }
      const std::uint64_t key
          = (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
      if (!seen_edges.insert(key).second) {
        continue;
      }

      const GridPoint p0{std::max(0, std::min(x_grids - 1, rsmt.branch[a].x)),
                         std::max(0, std::min(y_grids - 1, rsmt.branch[a].y))};
      const GridPoint p1{std::max(0, std::min(x_grids - 1, rsmt.branch[b].x)),
                         std::max(0, std::min(y_grids - 1, rsmt.branch[b].y))};
      emit_edge(p0, p1);
    }
  } else {
    const int n = static_cast<int>(pins_grid.size());

    std::vector<int> parent(n, -1);
    std::vector<int> best_dist(n, std::numeric_limits<int>::max());
    std::vector<bool> in_tree(n, false);
    best_dist[0] = 0;

    for (int iter = 0; iter < n; iter++) {
      int u = -1;
      int u_best = std::numeric_limits<int>::max();
      for (int i = 0; i < n; i++) {
        if (!in_tree[i] && best_dist[i] < u_best) {
          u_best = best_dist[i];
          u = i;
        }
      }
      if (u == -1) {
        break;
      }
      in_tree[u] = true;

      const GridPoint pu = pins_grid[u];
      for (int v = 0; v < n; v++) {
        if (in_tree[v]) {
          continue;
        }
        const GridPoint pv = pins_grid[v];
        const int dist = std::abs(pu.gx - pv.gx) + std::abs(pu.gy - pv.gy);
        if (dist < best_dist[v]) {
          best_dist[v] = dist;
          parent[v] = u;
        }
      }
    }

    for (int v = 1; v < n; v++) {
      const int u = parent[v];
      if (u < 0 || u >= n) {
        continue;
      }
      emit_edge(pins_grid[u], pins_grid[v]);
    }
  }

  if (h_layer != v_layer) {
    for (const auto& [key, mask] : endpoint_layer_mask) {
      if ((mask & 0x3) != 0x3) {
        continue;
      }
      const int x = static_cast<int>(key >> 32);
      const int y = static_cast<int>(key & 0xffffffffLL);
      add_via_unique(odb::Point(x, y), h_layer, v_layer);
    }
  }

  return route;
}

RouteMetrics compute_route_metrics(const GRoute& segments,
                                  odb::dbTech* tech)
{
  RouteMetrics metrics;
  for (const GSegment& segment : segments) {
    if (segment.isVia()) {
      metrics.via_count++;
    } else {
      metrics.wirelength_dbu += std::abs(segment.final_x - segment.init_x)
                                + std::abs(segment.final_y - segment.init_y);
    }
  }

  if (metrics.wirelength_dbu > 0 && tech != nullptr) {
    metrics.wirelength_um
        = metrics.wirelength_dbu
          / static_cast<double>(tech->getDbUnitsPerMicron());
  }
  return metrics;
}

void union_shorter_guides(NetRouteMap& base_routes,
                          const NetRouteMap& direct_routes,
                          odb::dbTech* tech,
                          GlobalRouter* grouter,
                          utl::Logger* logger)
{
  if (grouter == nullptr || logger == nullptr) {
    return;
  }

  // Budget guardrails: guide union can balloon the DR search space.
  constexpr int kMaxTotalAddedSegments = 160000;
  constexpr int kMaxAddedSegmentsPerNet = 280;
  // Aggressive union: accept any non-worsening direct candidate. The intent is
  // to inject "short corridors" into the guide set and let the detailed router
  // pick shorter realizations even when the congestion-safe global topology
  // detoured.
  constexpr double kMinRelativeImprovement = 0.0;
  // Via count is secondary for this phase; allow some extra vias to escape
  // DR-wirelength local minima.
  constexpr int kMaxExtraViasPerNet = 40;

  // More radical than "union": if the direct candidate is significantly
  // shorter and appears capacity-feasible on the current congestion DB,
  // replace the base guide topology for that net.
  constexpr double kReplaceMinRelativeImprovement = 0.003;  // 0.3% shorter
  constexpr int kReplaceMaxBlockedEdges = 0;

  const int tile_size = grouter->getTileSize();
  const long min_net_wl_dbu = std::max<long>(0, static_cast<long>(tile_size) * 20);

  GuideUnionStats stats;
  int total_added = 0;

  auto segment_is_line = [](const GSegment& segment) -> bool {
    const int dimensionality = (segment.init_x != segment.final_x)
                               + (segment.init_y != segment.final_y)
                               + (segment.init_layer != segment.final_layer);
    return dimensionality == 1;
  };

  auto segments_connect = [](const GSegment& segment1,
                             const GSegment& segment2) -> bool {
    auto [s1_min_x, s1_max_x]
        = std::minmax(segment1.init_x, segment1.final_x);
    auto [s1_min_y, s1_max_y]
        = std::minmax(segment1.init_y, segment1.final_y);
    auto [s1_min_z, s1_max_z]
        = std::minmax(segment1.init_layer, segment1.final_layer);
    auto [s2_min_x, s2_max_x]
        = std::minmax(segment2.init_x, segment2.final_x);
    auto [s2_min_y, s2_max_y]
        = std::minmax(segment2.init_y, segment2.final_y);
    auto [s2_min_z, s2_max_z]
        = std::minmax(segment2.init_layer, segment2.final_layer);
    return (s1_max_x >= s2_min_x && s1_min_x <= s2_max_x)
           && (s1_max_y >= s2_min_y && s1_min_y <= s2_max_y)
           && (s1_max_z >= s2_min_z && s1_min_z <= s2_max_z);
  };

  auto route_is_connected = [&](const GRoute& route) -> bool {
    const int total_segments = static_cast<int>(route.size());
    if (total_segments <= 1) {
      return true;
    }
    if (!segment_is_line(route[0])) {
      return false;
    }

    std::vector<int> parent(total_segments);
    std::vector<int> rank(total_segments, 0);
    for (int i = 0; i < total_segments; i++) {
      parent[i] = i;
    }
    int groups = 1;

    std::function<int(int)> find = [&](int x) -> int {
      if (parent[x] != x) {
        parent[x] = find(parent[x]);
      }
      return parent[x];
    };

    std::function<void(int, int)> unite = [&](int u, int v) {
      int root_u = find(u);
      int root_v = find(v);
      if (root_u == root_v) {
        return;
      }
      if (rank[root_u] > rank[root_v]) {
        parent[root_v] = root_u;
      } else if (rank[root_u] < rank[root_v]) {
        parent[root_u] = root_v;
      } else {
        parent[root_v] = root_u;
        rank[root_u]++;
      }
      groups--;
    };

    for (int i = 1; i < total_segments; i++) {
      if (!segment_is_line(route[i])) {
        return false;
      }
      groups++;
      for (int j = i - 1; j >= 0 && groups > 1; --j) {
        if (segments_connect(route[i], route[j])) {
          unite(i, j);
          if (groups == 1) {
            break;
          }
        }
      }
    }
    return groups == 1;
  };

  auto count_blocked_edges = [&](odb::dbNet* db_net,
                                 const std::unordered_set<std::uint64_t>& credit_edges,
                                 const GRoute& route) -> int {
    FastRouteCore* fastroute = (grouter != nullptr) ? grouter->fastroute() : nullptr;
    const Grid* grid = (grouter != nullptr) ? grouter->grid() : nullptr;
    if (fastroute == nullptr || grid == nullptr || db_net == nullptr) {
      return std::numeric_limits<int>::max();
    }

    const int tile = grid->getTileSize();
    if (tile <= 0) {
      return std::numeric_limits<int>::max();
    }

    const int x0_corner = grid->getXMin();
    const int y0_corner = grid->getYMin();
    const int x_grids = grid->getXGrids();
    const int y_grids = grid->getYGrids();

    auto to_grid_idx = [&](const int x_dbu, const int y_dbu) -> std::pair<int, int> {
      int gx = (x_dbu - x0_corner) / tile;
      int gy = (y_dbu - y0_corner) / tile;
      gx = std::max(0, std::min(x_grids - 1, gx));
      gy = std::max(0, std::min(y_grids - 1, gy));
      return {gx, gy};
    };

    int blocked = 0;
    for (const auto& seg : route) {
      if (seg.isVia()) {
        continue;
      }

      const int layer = seg.init_layer;
      const int req = fastroute->getDbNetLayerEdgeCost(db_net, layer);

      auto [gx0, gy0] = to_grid_idx(seg.init_x, seg.init_y);
      auto [gx1, gy1] = to_grid_idx(seg.final_x, seg.final_y);
      if (gx0 == gx1 && gy0 == gy1) {
        continue;
      }

      if (gy0 == gy1) {
        const int y = gy0;
        const int x_min = std::min(gx0, gx1);
        const int x_max = std::max(gx0, gx1);
        for (int x = x_min; x < x_max; x++) {
          int avail = fastroute->getAvailableResources(x, y, x + 1, y, layer);
          if (credit_edges.find(encode_unit_edge_key(x, y, x + 1, y, layer))
              != credit_edges.end()) {
            avail += req;
          }
          if (avail < req) {
            blocked++;
          }
        }
      } else if (gx0 == gx1) {
        const int x = gx0;
        const int y_min = std::min(gy0, gy1);
        const int y_max = std::max(gy0, gy1);
        for (int y = y_min; y < y_max; y++) {
          int avail = fastroute->getAvailableResources(x, y, x, y + 1, layer);
          if (credit_edges.find(encode_unit_edge_key(x, y, x, y + 1, layer))
              != credit_edges.end()) {
            avail += req;
          }
          if (avail < req) {
            blocked++;
          }
        }
      }
    }
    return blocked;
  };

  for (auto& [db_net, base] : base_routes) {
    stats.nets_checked++;

    if (total_added >= kMaxTotalAddedSegments) {
      stats.nets_skipped_over_budget++;
      break;
    }

    const auto it = direct_routes.find(db_net);
    if (it == direct_routes.end()) {
      continue;
    }
    const GRoute& direct = it->second;

    const RouteMetrics base_m = compute_route_metrics(base, tech);
    const RouteMetrics direct_m = compute_route_metrics(direct, tech);

    if (base_m.wirelength_dbu < min_net_wl_dbu) {
      stats.nets_skipped_small++;
      continue;
    }

    const long base_wl = base_m.wirelength_dbu;
    const long direct_wl = direct_m.wirelength_dbu;
    if (direct_wl >= static_cast<long>(base_wl * (1.0 - kMinRelativeImprovement))) {
      stats.nets_skipped_no_improvement++;
      continue;
    }

    const int base_vias = static_cast<int>(base_m.via_count);
    const int direct_vias = static_cast<int>(direct_m.via_count);
    if (direct_vias > base_vias + kMaxExtraViasPerNet) {
      continue;
    }

    const long min_abs_gain_dbu = static_cast<long>(tile_size) * 4;
    const bool can_replace
        = (base_wl - direct_wl) >= min_abs_gain_dbu
          && direct_wl
                 < static_cast<long>(base_wl * (1.0 - kReplaceMinRelativeImprovement))
          && route_is_connected(direct);
    if (can_replace) {
      const std::unordered_set<std::uint64_t> credit_edges
          = collect_unit_grid_edges(base, (grouter != nullptr) ? grouter->grid()
                                                               : nullptr);
      const int blocked_edges = count_blocked_edges(db_net, credit_edges, direct);
      if (blocked_edges <= kReplaceMaxBlockedEdges) {
        base = direct;
        stats.nets_replaced++;
        continue;
      }
    }

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(base.size() * 2 + 64);
    for (const auto& seg : base) {
      existing.insert(seg);
    }

    const std::size_t original_size = base.size();

    // Ensure we never introduce disconnected components:
    // grow a connected component starting from existing base endpoints.
    struct RouteNodeKey
    {
      int x = 0;
      int y = 0;
      int layer = 0;

      bool operator==(const RouteNodeKey& other) const
      {
        return x == other.x && y == other.y && layer == other.layer;
      }
    };

    struct RouteNodeKeyHash
    {
      std::size_t operator()(const RouteNodeKey& k) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.x);
        boost::hash_combine(seed, k.y);
        boost::hash_combine(seed, k.layer);
        return seed;
      }
    };

    auto endpoint_a = [](const GSegment& s) -> RouteNodeKey {
      return RouteNodeKey{static_cast<int>(s.init_x),
                          static_cast<int>(s.init_y),
                          static_cast<int>(s.init_layer)};
    };
    auto endpoint_b = [](const GSegment& s) -> RouteNodeKey {
      return RouteNodeKey{static_cast<int>(s.final_x),
                          static_cast<int>(s.final_y),
                          static_cast<int>(s.final_layer)};
    };

    std::unordered_set<RouteNodeKey, RouteNodeKeyHash> connected_nodes;
    connected_nodes.reserve(base.size() * 2 + 64);
    for (const auto& seg : base) {
      connected_nodes.insert(endpoint_a(seg));
      connected_nodes.insert(endpoint_b(seg));
    }

    int added_for_net = 0;
    int wires_added_for_net = 0;
    int vias_added_for_net = 0;
    bool progressed = true;
    // Multiple passes to allow adding a chain of segments where only the first
    // one touches the base.
    while (progressed && added_for_net < kMaxAddedSegmentsPerNet
           && (total_added + added_for_net) < kMaxTotalAddedSegments) {
      progressed = false;
      for (const auto& seg : direct) {
        if (added_for_net >= kMaxAddedSegmentsPerNet
            || (total_added + added_for_net) >= kMaxTotalAddedSegments) {
          break;
        }
        if (existing.find(seg) != existing.end()) {
          continue;
        }

        const RouteNodeKey a = endpoint_a(seg);
        const RouteNodeKey b = endpoint_b(seg);
        if (connected_nodes.find(a) == connected_nodes.end()
            && connected_nodes.find(b) == connected_nodes.end()) {
          continue;  // would create a disconnected component
        }

        if (existing.insert(seg).second) {
          base.push_back(seg);
          connected_nodes.insert(a);
          connected_nodes.insert(b);
          added_for_net++;
          if (seg.isVia()) {
            vias_added_for_net++;
          } else {
            wires_added_for_net++;
          }
          progressed = true;
        }
      }
    }

    if (added_for_net > 0) {
      // GlobalRouter enforces connectedness for *every* net route. If our
      // additions break this invariant, drop them.
      if (!route_is_connected(base)) {
        base.resize(original_size);
        stats.nets_reverted_disconnected++;
        continue;
      }

      stats.nets_merged++;
      total_added += added_for_net;
      stats.segments_added += added_for_net;
      stats.wires_added += wires_added_for_net;
      stats.vias_added += vias_added_for_net;
    }
  }

  if (stats.segments_added > 0 || stats.nets_replaced > 0) {
    logger->info(
        GNR,
        6012,
        "NEWGR guide union/replace: unioned {}, replaced {}, reverted {}, added segs {} (wires {}, vias {}), budget {}",
        stats.nets_merged,
        stats.nets_replaced,
        stats.nets_reverted_disconnected,
        stats.segments_added,
        stats.wires_added,
        stats.vias_added,
        kMaxTotalAddedSegments);
  } else {
    logger->info(GNR, 6013, "NEWGR guide union: no nets met merge criteria");
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

  auto compute_metrics = [&](const NetRouteMap& routes) -> RouteMetrics {
    RouteMetrics metrics;
    for (const auto& [db_net, segments] : routes) {
      static_cast<void>(db_net);
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          metrics.via_count++;
        } else {
          metrics.wirelength_dbu
              += std::abs(segment.final_x - segment.init_x)
                 + std::abs(segment.final_y - segment.init_y);
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

  auto capture_snapshot = [&]() -> RouterSnapshot {
    RouterSnapshot snapshot;
    snapshot.caps_percentage = grouter_->caps_perturbation_percentage_;
    snapshot.perturbation_amount = grouter_->perturbation_amount_;
    snapshot.congestion_iterations = grouter_->congestion_iterations_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setCongestionIterations(snapshot.congestion_iterations);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
  };

  RouterSnapshot snapshot = capture_snapshot();

  auto run_candidate = [&](const CandidateSettings& settings,
                           bool keep_routes) -> CandidateResult {
    CandidateResult result;
    result.settings = settings;

    restore_snapshot(snapshot);
    grouter_->setAllowCongestion(settings.allow_congestion);
    grouter_->setCapacitiesPerturbationPercentage(
        settings.caps_perturbation_percentage);
    grouter_->setPerturbationAmount(settings.perturbation_amount);
    grouter_->setCongestionIterations(settings.congestion_iterations);
    grouter_->setSeed(settings.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        settings.critical_nets_percentage);

    std::vector<Net*> candidate_nets
        = grouter_->initFastRoute(min_routing_layer, max_routing_layer);
    NetRouteMap routes;
    if (!candidate_nets.empty()) {
      routes = grouter_->findRouting(
          candidate_nets, min_routing_layer, max_routing_layer);
    }

    result.metrics = compute_metrics(routes);
    result.overflow = (grouter_->fastroute_ != nullptr)
                          ? grouter_->fastroute_->totalOverflow()
                          : std::numeric_limits<int>::max();
    if (keep_routes) {
      result.routes = std::move(routes);
    }

    logger_->info(GNR,
                  6006,
                  "NEWGR candidate {}: wl {:.0f} um, vias {}, overflow {}",
                  result.settings.name,
                  result.metrics.wirelength_um,
                  result.metrics.via_count,
                  result.overflow);

    return result;
  };

  CandidateSettings baseline;
  baseline.name = "baseline";
  baseline.seed = snapshot.seed;
  baseline.caps_perturbation_percentage = snapshot.caps_percentage;
  baseline.perturbation_amount = snapshot.perturbation_amount;
  baseline.congestion_iterations = snapshot.congestion_iterations;
  baseline.critical_nets_percentage = snapshot.critical_percentage;
  baseline.allow_congestion = false;

  // Wirelength-focused candidate with runtime guardrails:
  // - Keep congestion iterations bounded (fewer detours + faster runtime).
  // - Disable critical-net ordering to avoid STA overhead and additional
  //   ripup/re-route churn.
  // - Use modest capacity perturbation to avoid pathological tie-breaking.
  CandidateSettings wl_lean = baseline;
  wl_lean.name = "wl-lean";
  // Seed affects tie-breaking in routing and can meaningfully change final WL.
  // Keep a deterministic seed. Empirically, small changes here can shift the
  // rip-up/reroute tie-breaking enough to reduce DR wirelength without adding
  // extra global-routing passes.
  wl_lean.seed = 37;
  // Keep perturbation modest: helps with tie-breaking without pushing the
  // router into longer detours.
  wl_lean.caps_perturbation_percentage
      = std::max(1.0f, snapshot.caps_percentage);
  wl_lean.perturbation_amount = std::max(1, snapshot.perturbation_amount);
  wl_lean.congestion_iterations
      = std::min(30, snapshot.congestion_iterations);
  wl_lean.critical_nets_percentage = 0.0f;
  wl_lean.allow_congestion = false;

  // Directness-first candidate used only to provide an alternative guide tree
  // to the detailed router (via guide union). This run is allowed to end with
  // overflow; the overflow-free base candidate is still selected and returned.
  CandidateSettings wl_short = wl_lean;
  wl_short.name = "wl-short";
  wl_short.allow_congestion = true;
  wl_short.congestion_iterations = std::min(10, wl_lean.congestion_iterations);
  wl_short.caps_perturbation_percentage = 0.0f;
  wl_short.critical_nets_percentage = 0.0f;
  wl_short.seed = 11;

  auto is_better = [&](const CandidateResult& current,
                       const CandidateResult& best) -> bool {
    const bool current_ok = current.overflow == 0;
    const bool best_ok = best.overflow == 0;

    if (current_ok != best_ok) {
      return current_ok;  // prefer overflow-free
    }
    if (!current_ok && current.overflow != best.overflow) {
      return current.overflow < best.overflow;
    }
    if (current.metrics.wirelength_dbu != best.metrics.wirelength_dbu) {
      return current.metrics.wirelength_dbu < best.metrics.wirelength_dbu;
    }
    return current.metrics.via_count < best.metrics.via_count;
  };

  // One-pass default: run the WL-lean configuration and only fall back if
  // overflow remains.
  CandidateResult direct_fastroute = run_candidate(wl_short, /*keep_routes=*/true);
  CandidateResult primary = run_candidate(wl_lean, /*keep_routes=*/true);

  odb::dbTech* tech = (grouter_->db_ != nullptr) ? grouter_->db_->getTech()
                                                 : nullptr;
  stt::SteinerTreeBuilder* stt_builder = grouter_->stt_builder_;
  auto build_direct_candidate = [&](const NetRouteMap& base_routes) -> NetRouteMap {
    NetRouteMap direct_routes;
    const int tile = (grouter_ != nullptr) ? grouter_->getTileSize() : 0;
    const long min_net_wl_dbu = std::max<long>(0, static_cast<long>(tile) * 20);

    for (const auto& [db_net, base_segments] : base_routes) {
      const auto it_fast = direct_fastroute.routes.find(db_net);
      const GRoute* fast_route
          = (it_fast != direct_fastroute.routes.end()) ? &it_fast->second
                                                       : nullptr;

      const RouteMetrics base_m = compute_route_metrics(base_segments, tech);
      if (base_m.wirelength_dbu < min_net_wl_dbu) {
        continue;
      }

      // Prefer the "slack-shortest" per-net candidate (capacity-feasible,
      // directness-first). Fall back to the overflow-tolerant FastRoute run if
      // needed (keeps prior behavior for pathological nets).
      GRoute pattern_route = build_pattern_direct_route(db_net,
                                                        grouter_,
                                                        stt_builder,
                                                        tech,
                                                        min_routing_layer,
                                                        max_routing_layer,
                                                        &base_segments);
      if (!pattern_route.empty()) {
        direct_routes[db_net] = std::move(pattern_route);
        continue;
      }

      if (fast_route != nullptr && !fast_route->empty()) {
        direct_routes[db_net] = *fast_route;
      }
    }
    return direct_routes;
  };

  if (primary.overflow == 0) {
    logger_->info(GNR,
                  6005,
                  "NEWGR selected {}: wl {:.0f} um, vias {}, overflow {}",
                  primary.settings.name,
                  primary.metrics.wirelength_um,
                  primary.metrics.via_count,
                  primary.overflow);
    NetRouteMap routes = std::move(primary.routes);
    NetRouteMap direct_routes = build_direct_candidate(routes);
    union_shorter_guides(
        routes, direct_routes, tech, grouter_, logger_);
    add_rsmt_rectilinear_corridors(routes,
                                   grouter_,
                                   stt_builder,
                                   tech,
                                   logger_,
                                   min_routing_layer,
                                   max_routing_layer);
    add_congestion_patches(routes,
                           grouter_->fastroute(),
                           grouter_->grid(),
                           tech,
                           logger_,
                           min_routing_layer,
                           max_routing_layer);
    add_endpoint_access_patches(routes,
                                grouter_->fastroute(),
                                grouter_->grid(),
                                tech,
                                logger_,
                                min_routing_layer,
                                max_routing_layer);
    add_lane_patches(
        routes, grouter_, logger_, min_routing_layer, max_routing_layer);
    return routes;
  }

  CandidateResult recovery = run_candidate(baseline, /*keep_routes=*/true);
  const bool primary_is_best = is_better(primary, recovery);
  CandidateResult best = primary_is_best ? std::move(primary)
                                        : std::move(recovery);

  // Ensure GlobalRouter/FastRoute internal state (congestion DB, etc.) matches
  // the returned routes. This only triggers when the "primary" candidate wins
  // despite running the recovery candidate after it.
  if (primary_is_best) {
    best = run_candidate(best.settings, /*keep_routes=*/true);
  }

  logger_->info(GNR,
                6007,
                "NEWGR selected {}: wl {:.0f} um, vias {}, overflow {}",
                best.settings.name,
                best.metrics.wirelength_um,
                best.metrics.via_count,
                best.overflow);

  NetRouteMap routes = std::move(best.routes);
  NetRouteMap direct_routes = build_direct_candidate(routes);
  union_shorter_guides(
      routes, direct_routes, tech, grouter_, logger_);
  add_rsmt_rectilinear_corridors(routes,
                                 grouter_,
                                 stt_builder,
                                 tech,
                                 logger_,
                                 min_routing_layer,
                                 max_routing_layer);
  add_congestion_patches(routes,
                         grouter_->fastroute(),
                         grouter_->grid(),
                         tech,
                         logger_,
                         min_routing_layer,
                         max_routing_layer);
  add_endpoint_access_patches(routes,
                              grouter_->fastroute(),
                              grouter_->grid(),
                              tech,
                              logger_,
                              min_routing_layer,
                              max_routing_layer);
  add_lane_patches(
      routes, grouter_, logger_, min_routing_layer, max_routing_layer);
  return routes;
}

}  // namespace grt
