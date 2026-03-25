#include "NEWGR/NewGR.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "Net.h"
#include "grt/Rudy.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouteMetrics
{
  long wirelength_dbu = 0;
  long via_count = 0;
  long detour_dbu = 0;
  long high_layer_dbu = 0;
  int layer_span_sum = 0;
  int segment_count = 0;
  int bend_count = 0;
  double wirelength_um = 0.0;
  double score = 0.0;
  int overflow_edges = 0;
  int near_capacity_edges = 0;
  double max_usage_ratio = 0.0;
  double overflow_ratio_sum = 0.0;
};

struct ScenarioResult
{
  std::string name;
  RouteMetrics metrics;
  NetRouteMap routes;
};

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
};

struct ScenarioDefinition
{
  std::string name;
  std::function<void()> pre_init;
  std::function<void()> post_init;
};

struct Hotspot
{
  int gx = 0;
  int gy = 0;
  float severity = 1.0f;
  bool affect_horizontal = false;
  bool affect_vertical = false;
};

using RudyGrid = std::vector<std::vector<float>>;

struct CongestionSummary
{
  int overflow_edges = 0;
  int near_capacity_edges = 0;
  double max_usage_ratio = 0.0;
  double overflow_ratio_sum = 0.0;
};

using EdgeCountMap = std::unordered_map<uint64_t, int>;

struct RouteEdgeStats
{
  long wirelength_dbu = 0;
  int via_count = 0;
  long high_layer_dbu = 0;
  int layer_span = 0;
  int segment_count = 0;
  int bend_count = 0;
  EdgeCountMap edge_counts;
};

double estimateDetailedRouteProxyCost(const RouteMetrics& metrics)
{
  // DR-proxy objective mixed from FastRoute/CUGR/SPRoute ideas:
  // 1) keep guide wirelength first-order;
  // 2) penalize guide detour and upper-layer drift (detail-route stability);
  // 3) softly penalize via-heavy choices;
  // 4) keep overflow terms dominant whenever present.
  const double via_term = static_cast<double>(metrics.via_count) * 0.72;
  const double detour_term = static_cast<double>(metrics.detour_dbu) * 0.0016;
  const double high_layer_term
      = static_cast<double>(metrics.high_layer_dbu) * 0.00010;
  const double span_term = static_cast<double>(metrics.layer_span_sum) * 4.0;
  const double segment_term = static_cast<double>(metrics.segment_count) * 2.4;
  const double bend_term = static_cast<double>(metrics.bend_count) * 14.0;
  const double hotspot_term
      = static_cast<double>(metrics.near_capacity_edges) * 12.0;
  const double overflow_term = static_cast<double>(metrics.overflow_edges) * 2000.0
                               + static_cast<double>(metrics.overflow_ratio_sum)
                                     * 6000.0;
  return static_cast<double>(metrics.wirelength_dbu) + via_term + detour_term
         + high_layer_term + span_term + segment_term + bend_term + hotspot_term
         + overflow_term;
}

long getRouteBBoxHpwl(const GRoute& route)
{
  int min_x = std::numeric_limits<int>::max();
  int min_y = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int max_y = std::numeric_limits<int>::min();
  bool has_wire = false;
  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      continue;
    }
    has_wire = true;
    min_x = std::min(min_x, std::min(segment.init_x, segment.final_x));
    min_y = std::min(min_y, std::min(segment.init_y, segment.final_y));
    max_x = std::max(max_x, std::max(segment.init_x, segment.final_x));
    max_y = std::max(max_y, std::max(segment.init_y, segment.final_y));
  }
  if (!has_wire) {
    return 0L;
  }
  return static_cast<long>(max_x - min_x) + static_cast<long>(max_y - min_y);
}

uint64_t packEdgeKey(int x, int y, int layer, bool horizontal)
{
  const uint64_t orient = horizontal ? 1ULL : 0ULL;
  const uint64_t ux = static_cast<uint64_t>(std::max(x, 0)) & 0x1FFFFFULL;
  const uint64_t uy = static_cast<uint64_t>(std::max(y, 0)) & 0x1FFFFFULL;
  const uint64_t ul = static_cast<uint64_t>(std::max(layer, 0)) & 0xFFULL;
  return orient | (ux << 1) | (uy << 22) | (ul << 43);
}

void unpackEdgeKey(uint64_t key, int& x, int& y, int& layer, bool& horizontal)
{
  horizontal = (key & 1ULL) != 0ULL;
  x = static_cast<int>((key >> 1) & 0x1FFFFFULL);
  y = static_cast<int>((key >> 22) & 0x1FFFFFULL);
  layer = static_cast<int>((key >> 43) & 0xFFULL);
}

RouteEdgeStats collectRouteEdgeStats(const GRoute& route, Grid* grid)
{
  RouteEdgeStats stats;
  const int tile_size = grid != nullptr ? std::max(grid->getTileSize(), 1) : 1;
  const int x_min = grid != nullptr ? grid->getXMin() : 0;
  const int y_min = grid != nullptr ? grid->getYMin() : 0;
  const int x_grids = grid != nullptr ? grid->getXGrids() : 0;
  const int y_grids = grid != nullptr ? grid->getYGrids() : 0;
  const int max_x_idx = std::max(0, x_grids - 1);
  const int max_y_idx = std::max(0, y_grids - 1);
  int min_wire_layer = std::numeric_limits<int>::max();
  int max_wire_layer = std::numeric_limits<int>::min();
  int prev_orient = -1;
  int prev_layer = -1;

  auto coord_to_grid = [&](int coord, int min_coord, int max_index) {
    if (max_index <= 0) {
      return 0;
    }
    return std::clamp((coord - min_coord) / tile_size, 0, max_index);
  };

  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      stats.via_count++;
      continue;
    }

    const long seg_wl = std::abs(segment.final_x - segment.init_x)
                        + std::abs(segment.final_y - segment.init_y);
    stats.wirelength_dbu += seg_wl;
    stats.segment_count++;
    if (segment.init_layer > 0) {
      min_wire_layer = std::min(min_wire_layer, segment.init_layer);
      max_wire_layer = std::max(max_wire_layer, segment.init_layer);
      // Favor compact lower-layer guides to reduce detailed-router detours.
      const int high_layer_offset = std::max(0, segment.init_layer - 3);
      stats.high_layer_dbu += seg_wl * high_layer_offset;
    }
    const int orient = segment.init_y == segment.final_y ? 0 : 1;
    if (prev_orient >= 0 && prev_layer == segment.init_layer
        && prev_orient != orient) {
      stats.bend_count++;
    }
    prev_orient = orient;
    prev_layer = segment.init_layer;

    if (grid == nullptr || segment.init_layer <= 0) {
      continue;
    }

    if (segment.init_y == segment.final_y) {
      const int gy = coord_to_grid(segment.init_y, y_min, max_y_idx);
      const int gx0 = coord_to_grid(
          std::min(segment.init_x, segment.final_x), x_min, max_x_idx);
      const int gx1 = coord_to_grid(
          std::max(segment.init_x, segment.final_x), x_min, max_x_idx);
      for (int gx = gx0; gx < gx1; ++gx) {
        if (gx < 0 || gy < 0 || gx >= x_grids - 1 || gy >= y_grids) {
          continue;
        }
        stats.edge_counts[packEdgeKey(gx, gy, segment.init_layer, true)]++;
      }
    } else if (segment.init_x == segment.final_x) {
      const int gx = coord_to_grid(segment.init_x, x_min, max_x_idx);
      const int gy0 = coord_to_grid(
          std::min(segment.init_y, segment.final_y), y_min, max_y_idx);
      const int gy1 = coord_to_grid(
          std::max(segment.init_y, segment.final_y), y_min, max_y_idx);
      for (int gy = gy0; gy < gy1; ++gy) {
        if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids - 1) {
          continue;
        }
        stats.edge_counts[packEdgeKey(gx, gy, segment.init_layer, false)]++;
      }
    }
  }

  if (min_wire_layer <= max_wire_layer) {
    stats.layer_span = max_wire_layer - min_wire_layer;
  }

  return stats;
}

using HotspotMask = std::vector<std::vector<uint8_t>>;
using NodeDegreeMap = std::unordered_map<uint64_t, int>;

struct GuidePatchStats
{
  int nets_touched = 0;
  int endpoint_segments_added = 0;
  int long_segment_patches = 0;
  int via_patches_added = 0;
};

uint64_t packRouteNodeKey(int x, int y, int layer)
{
  const uint64_t ux = static_cast<uint64_t>(std::max(x, 0)) & 0x1FFFFFULL;
  const uint64_t uy = static_cast<uint64_t>(std::max(y, 0)) & 0x1FFFFFULL;
  const uint64_t ul = static_cast<uint64_t>(std::max(layer, 0)) & 0xFFULL;
  return ux | (uy << 21) | (ul << 42);
}

bool addUniqueGuideSegment(
    GRoute& route,
    std::unordered_set<GSegment, GSegmentHash>& seen,
    const GSegment& segment)
{
  if (segment.init_layer <= 0 || segment.final_layer <= 0) {
    return false;
  }
  if (!segment.isVia() && segment.length() <= 0) {
    return false;
  }
  auto [it, inserted] = seen.insert(segment);
  if (!inserted) {
    return false;
  }
  route.push_back(segment);
  return true;
}

NodeDegreeMap collectRouteNodeDegrees(const GRoute& route)
{
  NodeDegreeMap degrees;
  degrees.reserve(std::max<std::size_t>(route.size() * 2, 8));
  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      degrees[packRouteNodeKey(
          segment.init_x, segment.init_y, segment.init_layer)]++;
      degrees[packRouteNodeKey(
          segment.final_x, segment.final_y, segment.final_layer)]++;
      continue;
    }
    degrees[packRouteNodeKey(
        segment.init_x, segment.init_y, segment.init_layer)]++;
    degrees[packRouteNodeKey(
        segment.final_x, segment.final_y, segment.final_layer)]++;
  }
  return degrees;
}

HotspotMask buildHotspotMask(const std::vector<Hotspot>& hotspots,
                             int x_grids,
                             int y_grids,
                             float severity_threshold = 0.85f)
{
  HotspotMask mask;
  if (x_grids <= 0 || y_grids <= 0) {
    return mask;
  }
  mask.assign(x_grids, std::vector<uint8_t>(y_grids, 0));
  for (const Hotspot& hotspot : hotspots) {
    if (hotspot.severity < severity_threshold) {
      continue;
    }
    if (hotspot.gx < 0 || hotspot.gy < 0 || hotspot.gx >= x_grids
        || hotspot.gy >= y_grids) {
      continue;
    }
    mask[hotspot.gx][hotspot.gy] = 1;
  }
  return mask;
}

bool segmentTouchesHotspot(const GSegment& segment,
                           const HotspotMask& hotspot_mask,
                           Grid* grid)
{
  if (grid == nullptr || hotspot_mask.empty() || segment.isVia()) {
    return false;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return false;
  }

  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int max_x_idx = std::max(0, x_grids - 1);
  const int max_y_idx = std::max(0, y_grids - 1);

  auto coord_to_grid = [&](int coord, int min_coord, int max_index) {
    if (max_index <= 0) {
      return 0;
    }
    return std::clamp((coord - min_coord) / tile_size, 0, max_index);
  };

  if (segment.init_y == segment.final_y) {
    const int gy = coord_to_grid(segment.init_y, y_min, max_y_idx);
    const int gx0 = coord_to_grid(
        std::min(segment.init_x, segment.final_x), x_min, max_x_idx);
    const int gx1 = coord_to_grid(
        std::max(segment.init_x, segment.final_x), x_min, max_x_idx);
    for (int gx = gx0; gx < gx1; ++gx) {
      if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
        continue;
      }
      if (hotspot_mask[gx][gy] != 0) {
        return true;
      }
    }
  } else if (segment.init_x == segment.final_x) {
    const int gx = coord_to_grid(segment.init_x, x_min, max_x_idx);
    const int gy0 = coord_to_grid(
        std::min(segment.init_y, segment.final_y), y_min, max_y_idx);
    const int gy1 = coord_to_grid(
        std::max(segment.init_y, segment.final_y), y_min, max_y_idx);
    for (int gy = gy0; gy < gy1; ++gy) {
      if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
        continue;
      }
      if (hotspot_mask[gx][gy] != 0) {
        return true;
      }
    }
  }

  return false;
}

bool buildMidpointPatchSegment(const GSegment& base_segment,
                               int patch_layer,
                               Grid* grid,
                               GSegment& out_segment)
{
  if (grid == nullptr || base_segment.isVia() || base_segment.length() <= 0) {
    return false;
  }

  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int x_max = x_min + std::max(0, grid->getXGrids() - 1) * tile_size;
  const int y_max = y_min + std::max(0, grid->getYGrids() - 1) * tile_size;

  if (base_segment.init_y == base_segment.final_y) {
    const int x_start = std::min(base_segment.init_x, base_segment.final_x);
    const int x_end = std::max(base_segment.init_x, base_segment.final_x);
    const int span_steps = std::max(
        1, (x_end - x_start) / tile_size);
    if (span_steps < 2) {
      return false;
    }
    const int mid_step = span_steps / 2;
    int x0 = x_start + mid_step * tile_size;
    int x1 = std::min(x_end, x0 + tile_size);
    if (x1 <= x0) {
      x0 = std::max(x_start, x0 - tile_size);
      x1 = std::min(x_end, x0 + tile_size);
      if (x1 <= x0) {
        return false;
      }
    }
    const int y = std::clamp(base_segment.init_y, y_min, y_max);
    x0 = std::clamp(x0, x_min, x_max);
    x1 = std::clamp(x1, x_min, x_max);
    if (x1 <= x0) {
      return false;
    }
    out_segment = GSegment(x0, y, patch_layer, x1, y, patch_layer);
    return true;
  }

  if (base_segment.init_x == base_segment.final_x) {
    const int y_start = std::min(base_segment.init_y, base_segment.final_y);
    const int y_end = std::max(base_segment.init_y, base_segment.final_y);
    const int span_steps = std::max(
        1, (y_end - y_start) / tile_size);
    if (span_steps < 2) {
      return false;
    }
    const int mid_step = span_steps / 2;
    int y0 = y_start + mid_step * tile_size;
    int y1 = std::min(y_end, y0 + tile_size);
    if (y1 <= y0) {
      y0 = std::max(y_start, y0 - tile_size);
      y1 = std::min(y_end, y0 + tile_size);
      if (y1 <= y0) {
        return false;
      }
    }
    const int x = std::clamp(base_segment.init_x, x_min, x_max);
    y0 = std::clamp(y0, y_min, y_max);
    y1 = std::clamp(y1, y_min, y_max);
    if (y1 <= y0) {
      return false;
    }
    out_segment = GSegment(x, y0, patch_layer, x, y1, patch_layer);
    return true;
  }

  return false;
}

GuidePatchStats applyCugrGuidePatches(
    NetRouteMap& base_routes,
    const NetRouteMap* alt_routes_a,
    const NetRouteMap* alt_routes_b,
    const std::vector<Hotspot>& hotspots,
    Grid* grid,
    int min_layer,
    int max_layer,
    int max_endpoint_patches_per_net = 6,
    int max_long_patches_per_net = 2)
{
  GuidePatchStats stats;
  if (grid == nullptr || base_routes.empty()) {
    return stats;
  }

  const int tile_size = std::max(grid->getTileSize(), 1);
  const HotspotMask hotspot_mask
      = buildHotspotMask(hotspots, grid->getXGrids(), grid->getYGrids(), 0.80f);

  for (auto& [db_net, route] : base_routes) {
    if (db_net == nullptr || route.empty()) {
      continue;
    }

    const std::size_t original_size = route.size();
    std::unordered_set<GSegment, GSegmentHash> seen(route.begin(), route.end());
    NodeDegreeMap node_degrees = collectRouteNodeDegrees(route);
    std::vector<GSegment> base_snapshot = route;
    int endpoint_added_for_net = 0;
    int long_patch_added_for_net = 0;

    auto add_endpoint_relief = [&](const NetRouteMap* alt_routes) {
      if (alt_routes == nullptr
          || endpoint_added_for_net >= max_endpoint_patches_per_net) {
        return;
      }
      const auto alt_it = alt_routes->find(db_net);
      if (alt_it == alt_routes->end()) {
        return;
      }
      for (const GSegment& candidate : alt_it->second) {
        if (endpoint_added_for_net >= max_endpoint_patches_per_net) {
          break;
        }
        if (candidate.isVia() || candidate.init_layer != candidate.final_layer) {
          continue;
        }
        if (candidate.length() <= 0 || candidate.length() > tile_size * 5) {
          continue;
        }
        const uint64_t node0 = packRouteNodeKey(
            candidate.init_x, candidate.init_y, candidate.init_layer);
        const uint64_t node1 = packRouteNodeKey(
            candidate.final_x, candidate.final_y, candidate.final_layer);
        const auto degree0_it = node_degrees.find(node0);
        const auto degree1_it = node_degrees.find(node1);
        const int degree0 = degree0_it != node_degrees.end() ? degree0_it->second : 0;
        const int degree1 = degree1_it != node_degrees.end() ? degree1_it->second : 0;
        if (degree0 <= 0 && degree1 <= 0) {
          continue;
        }
        if (std::min(degree0, degree1) > 2) {
          continue;
        }
        if (addUniqueGuideSegment(route, seen, candidate)) {
          endpoint_added_for_net++;
          node_degrees[node0]++;
          node_degrees[node1]++;
        }
      }
    };

    // CUGR-inspired pin-region patching:
    // inject local alternatives from DR-stable/layer-compact candidates near
    // low-degree endpoints to improve detailed-route pin accessibility.
    add_endpoint_relief(alt_routes_a);
    add_endpoint_relief(alt_routes_b);

    // CUGR long-segment patching:
    // add short adjacent-layer guides over hotspot-crossing trunks.
    for (const GSegment& segment : base_snapshot) {
      if (long_patch_added_for_net >= max_long_patches_per_net) {
        break;
      }
      if (segment.isVia() || segment.init_layer != segment.final_layer) {
        continue;
      }
      if (segment.length() < tile_size * 4) {
        continue;
      }
      if (!segmentTouchesHotspot(segment, hotspot_mask, grid)) {
        continue;
      }

      int patch_layer = -1;
      if (segment.init_layer < max_layer) {
        patch_layer = segment.init_layer + 1;
      } else if (segment.init_layer > min_layer) {
        patch_layer = segment.init_layer - 1;
      }
      if (patch_layer < min_layer || patch_layer > max_layer) {
        continue;
      }

      GSegment patch_segment;
      if (!buildMidpointPatchSegment(segment, patch_layer, grid, patch_segment)) {
        continue;
      }
      if (!addUniqueGuideSegment(route, seen, patch_segment)) {
        continue;
      }
      long_patch_added_for_net++;
      stats.long_segment_patches++;

      if (patch_layer != segment.init_layer) {
        const int low_layer = std::min(segment.init_layer, patch_layer);
        const int high_layer = std::max(segment.init_layer, patch_layer);
        const GSegment via0(patch_segment.init_x,
                            patch_segment.init_y,
                            low_layer,
                            patch_segment.init_x,
                            patch_segment.init_y,
                            high_layer);
        const GSegment via1(patch_segment.final_x,
                            patch_segment.final_y,
                            low_layer,
                            patch_segment.final_x,
                            patch_segment.final_y,
                            high_layer);
        if (addUniqueGuideSegment(route, seen, via0)) {
          stats.via_patches_added++;
        }
        if (addUniqueGuideSegment(route, seen, via1)) {
          stats.via_patches_added++;
        }
      }
    }

    stats.endpoint_segments_added += endpoint_added_for_net;
    if (route.size() > original_size) {
      stats.nets_touched++;
    }
  }

  return stats;
}

RudyGrid computeNormalizedRudyGrid(Rudy* rudy)
{
  RudyGrid normalized;
  if (rudy == nullptr) {
    return normalized;
  }
  const auto [x_tiles, y_tiles] = rudy->getGridSize();
  if (x_tiles == 0 || y_tiles == 0) {
    return normalized;
  }

  normalized.resize(x_tiles, std::vector<float>(y_tiles, 0.0f));
  float max_value = 0.0f;
  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      const float val = rudy->getTile(x, y).getRudy();
      normalized[x][y] = val;
      max_value = std::max(max_value, val);
    }
  }

  if (max_value <= std::numeric_limits<float>::epsilon()) {
    return normalized;
  }

  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      normalized[x][y] = std::clamp(normalized[x][y] / max_value, 0.0f, 1.0f);
    }
  }

  return normalized;
}

void adjustEdgeCapacity(GlobalRouter* grouter,
                        int x1,
                        int y1,
                        int x2,
                        int y2,
                        int layer,
                        float ratio)
{
  ratio = std::clamp(ratio, 0.05f, 1.10f);
  FastRouteCore* core = grouter->fastroute();
  if (core == nullptr) {
    return;
  }
  const int current_cap = core->getEdgeCapacity(x1, y1, x2, y2, layer);
  if (current_cap <= 0) {
    return;
  }
  const int new_cap
      = std::max(1, static_cast<int>(std::floor(current_cap * ratio)));
  if (new_cap == current_cap) {
    return;
  }
  const bool is_reduce = new_cap < current_cap;
  core->addAdjustment(x1, y1, x2, y2, layer, new_cap, is_reduce);
}

void applySoftCapacityScaling(GlobalRouter* grouter,
                              const RudyGrid& normalized_rudy,
                              int min_layer,
                              int max_layer,
                              float min_ratio_base = 0.50f,
                              float max_ratio_base = 0.92f,
                              float slope = 6.0f,
                              float midpoint = 0.45f,
                              float reclaim_boost = 0.06f)
{
  Grid* grid = grouter->grid();
  if (normalized_rudy.empty() || grid == nullptr) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  const int usable_x = std::min(x_grids, x_tiles);
  const int usable_y = std::min(y_grids, y_tiles);
  if (usable_x == 0 || usable_y == 0) {
    return;
  }

  const int layer_span = std::max(max_layer - min_layer, 1);
  const auto logistic_ratio = [](float normalized,
                                 float slope,
                                 float midpoint,
                                 float min_ratio,
                                 float max_ratio,
                                 float reclaim_boost) {
    normalized = std::clamp(normalized, 0.0f, 1.0f);
    const float exponent = slope * (normalized - midpoint);
    const float logistic = 1.0f / (1.0f + std::exp(exponent));
    float blend = min_ratio + (max_ratio - min_ratio) * logistic;

    // FastRoute-style virtual-capacity reclamation:
    // give low-congestion regions a small capacity credit to preserve
    // shortest-path opportunities while hotspot edges stay constrained.
    const float low_congestion_limit = std::max(0.05f, midpoint * 0.75f);
    if (normalized < low_congestion_limit && reclaim_boost > 0.0f) {
      const float coolness = 1.0f - (normalized / low_congestion_limit);
      blend += reclaim_boost * std::clamp(coolness, 0.0f, 1.0f);
    }

    return std::clamp(blend, 0.05f, 1.10f);
  };

  const auto getNormalized = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float min_ratio = std::clamp(
        min_ratio_base + 0.12f * layer_factor, 0.05f, 0.99f);
    const float max_ratio = std::clamp(
        max_ratio_base + 0.04f * layer_factor, min_ratio, 0.995f);

    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        const float ratio
            = logistic_ratio(normalized,
                             slope,
                             midpoint,
                             min_ratio,
                             max_ratio,
                             reclaim_boost);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        const float ratio
            = logistic_ratio(normalized,
                             slope,
                             midpoint,
                             min_ratio,
                             max_ratio,
                             reclaim_boost);
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

void applyHotspotPenalties(GlobalRouter* grouter,
                           const std::vector<Hotspot>& hotspots,
                           int min_layer,
                           int max_layer,
                           int halo,
                           float base_ratio,
                           float severity_weight = 0.5f)
{
  Grid* grid = grouter->grid();
  if (hotspots.empty() || grid == nullptr) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  const int layer_span = std::max(max_layer - min_layer, 1);
  halo = std::max(0, halo);

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float layer_ratio
        = std::clamp(base_ratio + 0.1f * layer_factor, 0.4f, 0.95f);
    const float scaled_severity_weight = std::clamp(severity_weight, 0.0f, 1.0f);

    for (const Hotspot& hotspot : hotspots) {
      for (int dx = -halo; dx <= halo; ++dx) {
        for (int dy = -halo; dy <= halo; ++dy) {
          const int gx = hotspot.gx + dx;
          const int gy = hotspot.gy + dy;
          if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
            continue;
          }
          const float ratio_scale = std::clamp(
              1.0f - scaled_severity_weight * (hotspot.severity - 1.0f), 0.5f, 1.5f);
          const float adjusted_ratio
              = std::clamp(layer_ratio * ratio_scale, 0.25f, 0.98f);
          if (hotspot.affect_horizontal && gx < x_grids - 1) {
            adjustEdgeCapacity(
                grouter, gx, gy, gx + 1, gy, layer, adjusted_ratio);
          }
          if (hotspot.affect_vertical && gy < y_grids - 1) {
            adjustEdgeCapacity(
                grouter, gx, gy, gx, gy + 1, layer, adjusted_ratio);
          }
        }
      }
    }
  }
}

CongestionSummary collectCongestionSummary(GlobalRouter* grouter,
                                           float hot_edge_threshold = 0.85f)
{
  CongestionSummary summary;
  FastRouteCore* core = grouter->fastroute();
  if (core == nullptr) {
    return summary;
  }

  hot_edge_threshold = std::clamp(hot_edge_threshold, 0.5f, 1.0f);
  core->computeCongestionInformation();

  std::vector<CongestionInformation> vertical;
  std::vector<CongestionInformation> horizontal;
  core->getCongestionGrid(vertical, horizontal);

  auto accumulate = [&](const std::vector<CongestionInformation>& edges) {
    for (const auto& info : edges) {
      const int capacity = std::max(info.congestion.capacity, 1);
      const int usage = std::max(info.congestion.usage, 0);
      const double usage_ratio
          = static_cast<double>(usage) / static_cast<double>(capacity);
      summary.max_usage_ratio = std::max(summary.max_usage_ratio, usage_ratio);
      if (usage_ratio >= hot_edge_threshold) {
        summary.near_capacity_edges++;
      }
      if (usage_ratio > 1.0) {
        summary.overflow_edges++;
        summary.overflow_ratio_sum += (usage_ratio - 1.0);
      }
    }
  };

  accumulate(horizontal);
  accumulate(vertical);
  return summary;
}

bool parseCriticalScenarioName(const std::string& name,
                               int& seed,
                               float& critical_pct)
{
  int critical_x10 = 0;
  if (std::sscanf(name.c_str(), "critical-s%d-c%d", &seed, &critical_x10)
      != 2) {
    return false;
  }
  critical_pct = static_cast<float>(critical_x10) / 10.0f;
  return true;
}

bool parsePerturbScenarioName(const std::string& name,
                              int& seed,
                              float& perturb_pct)
{
  int perturb_x10 = 0;
  if (std::sscanf(name.c_str(), "perturb-s%d-p%d", &seed, &perturb_x10) != 2) {
    return false;
  }
  perturb_pct = static_cast<float>(perturb_x10) / 10.0f;
  return true;
}

int getSoftCapacityForEdge(uint64_t key,
                           FastRouteCore* core,
                           const RudyGrid& normalized_rudy,
                           std::unordered_map<uint64_t, int>& softcap_cache)
{
  if (core == nullptr) {
    return 1;
  }
  const auto cache_it = softcap_cache.find(key);
  if (cache_it != softcap_cache.end()) {
    return cache_it->second;
  }

  int x = 0;
  int y = 0;
  int layer = 0;
  bool horizontal = false;
  unpackEdgeKey(key, x, y, layer, horizontal);

  const int x2 = horizontal ? x + 1 : x;
  const int y2 = horizontal ? y : y + 1;
  const int hard_cap = std::max(1, core->getEdgeCapacity(x, y, x2, y2, layer));

  auto get_rudy = [&](int gx, int gy) {
    if (normalized_rudy.empty()) {
      return 0.0f;
    }
    if (gx < 0 || gy < 0 || gx >= static_cast<int>(normalized_rudy.size())
        || gy >= static_cast<int>(normalized_rudy.front().size())) {
      return 0.0f;
    }
    return normalized_rudy[gx][gy];
  };

  const float rudy = 0.5f
                     * (get_rudy(x, y)
                        + get_rudy(horizontal ? x + 1 : x,
                                   horizontal ? y : y + 1));
  const float exponent = (rudy - 0.42f) * 8.5f;
  const float ratio = std::clamp(
      0.56f + 0.38f / (1.0f + std::exp(exponent)), 0.50f, 0.96f);
  const int soft_cap
      = std::max(1, static_cast<int>(std::floor(hard_cap * ratio)));
  softcap_cache.emplace(key, soft_cap);
  return soft_cap;
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
      long route_wl = 0;
      long route_high_layer = 0;
      int route_segments = 0;
      int route_bends = 0;
      int min_wire_layer = std::numeric_limits<int>::max();
      int max_wire_layer = std::numeric_limits<int>::min();
      int prev_orient = -1;
      int prev_layer = -1;
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          metrics.via_count++;
        } else {
          const long seg_wl = std::abs(segment.final_x - segment.init_x)
                              + std::abs(segment.final_y - segment.init_y);
          metrics.wirelength_dbu += seg_wl;
          route_wl += seg_wl;
          route_segments++;
          if (segment.init_layer > 0) {
            min_wire_layer = std::min(min_wire_layer, segment.init_layer);
            max_wire_layer = std::max(max_wire_layer, segment.init_layer);
            const int high_layer_offset
                = std::max(0, segment.init_layer - (min_routing_layer + 1));
            route_high_layer += seg_wl * high_layer_offset;
          }
          const int orient = segment.init_y == segment.final_y ? 0 : 1;
          if (prev_orient >= 0 && prev_layer == segment.init_layer
              && prev_orient != orient) {
            route_bends++;
          }
          prev_orient = orient;
          prev_layer = segment.init_layer;
        }
      }
      const long bbox_hpwl = getRouteBBoxHpwl(segments);
      metrics.detour_dbu += std::max(0L, route_wl - bbox_hpwl);
      metrics.high_layer_dbu += route_high_layer;
      metrics.segment_count += route_segments;
      metrics.bend_count += route_bends;
      if (min_wire_layer <= max_wire_layer) {
        metrics.layer_span_sum += (max_wire_layer - min_wire_layer);
      }
    }

    if (metrics.wirelength_dbu > 0 && grouter_->db_ != nullptr
        && grouter_->db_->getTech() != nullptr) {
      metrics.wirelength_um
          = metrics.wirelength_dbu
            / static_cast<double>(
                grouter_->db_->getTech()->getDbUnitsPerMicron());
    }

    const double via_weight
        = static_cast<double>(std::max(grouter_->grid_->getTileSize(), 1))
          * 3.0;
    const double detour_weight = 0.22;
    const double high_layer_weight = 0.008;
    const double layer_span_weight
        = static_cast<double>(std::max(grouter_->grid_->getTileSize(), 1)) * 0.15;
    const double segment_weight = 3.2;
    const double bend_weight = 18.0;
    metrics.score = static_cast<double>(metrics.wirelength_dbu)
                    + via_weight * static_cast<double>(metrics.via_count)
                    + detour_weight * static_cast<double>(metrics.detour_dbu)
                    + high_layer_weight * static_cast<double>(metrics.high_layer_dbu)
                    + layer_span_weight * static_cast<double>(metrics.layer_span_sum)
                    + segment_weight * static_cast<double>(metrics.segment_count)
                    + bend_weight * static_cast<double>(metrics.bend_count);
    return metrics;
  };

  auto apply_routability_proxy = [&](RouteMetrics& metrics) {
    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const CongestionSummary summary = collectCongestionSummary(grouter_);
    metrics.overflow_edges = summary.overflow_edges;
    metrics.near_capacity_edges = summary.near_capacity_edges;
    metrics.max_usage_ratio = summary.max_usage_ratio;
    metrics.overflow_ratio_sum = summary.overflow_ratio_sum;

    // Congestion proxy inspired by CUGR probability cost:
    // keep routes close to shortest-path while avoiding high-overflow guides
    // that force large detailed-route detours.
    const double edge_hotspot_penalty
        = static_cast<double>(tile_size)
          * (15.0 * metrics.overflow_edges + 0.75 * metrics.near_capacity_edges);
    const double ratio_penalty
        = static_cast<double>(tile_size) * 35.0 * metrics.overflow_ratio_sum;
    metrics.score += edge_hotspot_penalty + ratio_penalty;
  };

  auto capture_snapshot = [&]() -> RouterSnapshot {
    RouterSnapshot snapshot;
    snapshot.caps_percentage = grouter_->caps_perturbation_percentage_;
    snapshot.perturbation_amount = grouter_->perturbation_amount_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
  };

  auto run_existing_state = [&](const std::string& name,
                                std::vector<Net*>& state_nets) {
    NetRouteMap routes;
    if (!state_nets.empty()) {
      routes = grouter_->findRouting(
          state_nets, min_routing_layer, max_routing_layer);
    }
    RouteMetrics metrics = compute_metrics(routes);
    apply_routability_proxy(metrics);
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow edges {}, "
                  "hot edges {}, max ratio {:.2f}, detour {}, high-layer {}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow_edges,
                  metrics.near_capacity_edges,
                  metrics.max_usage_ratio,
                  metrics.detour_dbu,
                  metrics.high_layer_dbu);
    return ScenarioResult{name, metrics, std::move(routes)};
  };

  auto collect_hotspots = [&]() -> std::vector<Hotspot> {
    std::vector<Hotspot> hotspots;
    if (grouter_->fastroute_ == nullptr || grouter_->grid_ == nullptr) {
      return hotspots;
    }

    grouter_->fastroute_->computeCongestionInformation();

    std::vector<CongestionInformation> vertical;
    std::vector<CongestionInformation> horizontal;
    grouter_->fastroute_->getCongestionGrid(vertical, horizontal);

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    auto append_hotspots = [&](const std::vector<CongestionInformation>& edges,
                               bool is_vertical) {
      for (const auto& info : edges) {
        const int capacity = std::max(info.congestion.capacity, 1);
        const float usage_ratio
            = static_cast<float>(info.congestion.usage)
              / static_cast<float>(capacity);
        if (usage_ratio < 0.6f) {
          continue;
        }
        const int gx
            = std::clamp((info.segment.init_x - x_min) / tile_size, 0, x_grids);
        const int gy
            = std::clamp((info.segment.init_y - y_min) / tile_size, 0, y_grids);
        Hotspot hotspot;
        hotspot.gx = std::clamp(gx, 0, std::max(x_grids - 1, 0));
        hotspot.gy = std::clamp(gy, 0, std::max(y_grids - 1, 0));
        hotspot.severity = std::clamp(usage_ratio, 0.6f, 3.0f);
        hotspot.affect_vertical = is_vertical;
        hotspot.affect_horizontal = !is_vertical;
        hotspots.push_back(hotspot);
      }
    };

    append_hotspots(horizontal, false);
    append_hotspots(vertical, true);
    return hotspots;
  };

  auto run_scenario = [&](const ScenarioDefinition& scenario,
                          const RouterSnapshot& snapshot) {
    restore_snapshot(snapshot);
    if (scenario.pre_init) {
      scenario.pre_init();
    }

    std::vector<Net*> scenario_nets
        = grouter_->initFastRoute(min_routing_layer, max_routing_layer);
    if (scenario.post_init) {
      scenario.post_init();
    }

    NetRouteMap routes;
    if (!scenario_nets.empty()) {
      routes = grouter_->findRouting(
          scenario_nets, min_routing_layer, max_routing_layer);
    }
    RouteMetrics metrics = compute_metrics(routes);
    apply_routability_proxy(metrics);
    logger_->info(GNR,
                  6006,
                  "NEWGR scenario {}: wirelength {:.0f} um, vias {}, overflow "
                  "edges {}, hot edges {}, max ratio {:.2f}, detour {}, "
                  "high-layer {}",
                  scenario.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow_edges,
                  metrics.near_capacity_edges,
                  metrics.max_usage_ratio,
                  metrics.detour_dbu,
                  metrics.high_layer_dbu);
    return ScenarioResult{scenario.name, metrics, std::move(routes)};
  };

  RouterSnapshot snapshot = capture_snapshot();

  ScenarioResult baseline
      = run_existing_state("baseline", nets);
  std::vector<Hotspot> hotspots = collect_hotspots();

  RudyGrid normalized_rudy;
  if (Rudy* rudy = grouter_->getRudy()) {
    rudy->calculateRudy();
    normalized_rudy = computeNormalizedRudyGrid(rudy);
  }

  std::vector<ScenarioResult> scenario_results;
  scenario_results.push_back(baseline);

  ScenarioDefinition baseline_def{"baseline", nullptr, nullptr};
  std::vector<ScenarioDefinition> scenario_defs;
  // Adaptive exploration mode:
  // - In low-congestion designs, favor a compact SPRoute-style deterministic
  //   sweep to stabilize shortest-guide candidates and reduce WL variance.
  // - In congested designs, keep broader exploration to preserve routability.
  const bool compact_exploration_mode
      = baseline.metrics.overflow_edges == 0
        && baseline.metrics.near_capacity_edges < 128
        && baseline.metrics.max_usage_ratio < 0.85;

  auto make_soft_config
      = [&](const std::string& name,
            float min_base,
            float max_base,
            float slope,
            float midpoint,
            int halo,
            float hotspot_ratio,
            float severity_weight,
            float reclaim_boost,
            float perturb_pct,
            int seed,
            float critical_pct) {
          ScenarioDefinition def;
          def.name = name;
          def.pre_init = [this, perturb_pct, seed, critical_pct]() {
            grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
            grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(seed);
            grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
          };
          def.post_init
              = [this,
                 &normalized_rudy,
                 &hotspots,
                 min_routing_layer,
                 max_routing_layer,
                 min_base,
                 max_base,
                 slope,
                 midpoint,
                 reclaim_boost,
                 halo,
                 hotspot_ratio,
                 severity_weight]() {
                  applySoftCapacityScaling(grouter_,
                                           normalized_rudy,
                                           min_routing_layer,
                                           max_routing_layer,
                                           min_base,
                                           max_base,
                                           slope,
                                           midpoint,
                                           reclaim_boost);
                  applyHotspotPenalties(grouter_,
                                        hotspots,
                                        min_routing_layer,
                                        max_routing_layer,
                                        halo,
                                        hotspot_ratio,
                                        severity_weight);
                };
          return def;
        };

  const bool enable_softcap_scenarios
      = !normalized_rudy.empty()
        && (baseline.metrics.overflow_edges > 0
            || baseline.metrics.near_capacity_edges > 200
            || baseline.metrics.max_usage_ratio > 0.92);

  if (enable_softcap_scenarios) {
    scenario_defs.push_back(make_soft_config("cugr-prob-strong",
                                             0.38f,
                                             0.86f,
                                             8.0f,
                                             0.52f,
                                             2,
                                             0.52f,
                                             0.85f,
                                             0.02f,
                                             0.0f,
                                             snapshot.seed,
                                             snapshot.critical_percentage));

    scenario_defs.push_back(make_soft_config("cugr-prob-balanced",
                                             0.44f,
                                             0.89f,
                                             7.0f,
                                             0.47f,
                                             2,
                                             0.58f,
                                             0.70f,
                                             0.05f,
                                             0.0f,
                                             snapshot.seed,
                                             snapshot.critical_percentage));

    scenario_defs.push_back(make_soft_config("soft-cap",
                                             0.52f,
                                             0.94f,
                                             5.5f,
                                             0.42f,
                                             1,
                                             0.68f,
                                             0.35f,
                                             0.07f,
                                             0.0f,
                                             snapshot.seed,
                                             snapshot.critical_percentage));

    scenario_defs.push_back(make_soft_config("guided-softcap",
                                             0.48f,
                                             0.90f,
                                             6.5f,
                                             0.48f,
                                             2,
                                             0.60f,
                                             0.55f,
                                             0.06f,
                                             3.5f,
                                             13,
                                             12.0f));

    scenario_defs.push_back(make_soft_config("mild-softcap",
                                             0.58f,
                                             0.97f,
                                             4.5f,
                                             0.38f,
                                             1,
                                             0.75f,
                                             0.20f,
                                             0.08f,
                                             2.5f,
                                             5,
                                             8.0f));
  } else if (!normalized_rudy.empty()) {
    // In low-overflow designs, stay close to baseline capacities and run
    // micro-guided variants to avoid introducing large detours.
    scenario_defs.push_back(make_soft_config("rudy-precision-direct",
                                             0.80f,
                                             1.02f,
                                             2.8f,
                                             0.24f,
                                             1,
                                             0.92f,
                                             0.10f,
                                             0.12f,
                                             0.0f,
                                             snapshot.seed,
                                             4.0f));
    scenario_defs.push_back(make_soft_config("rudy-precision-critical",
                                             0.74f,
                                             1.00f,
                                             3.5f,
                                             0.30f,
                                             1,
                                             0.88f,
                                             0.16f,
                                             0.10f,
                                             0.0f,
                                             snapshot.seed,
                                             6.0f));
    logger_->info(GNR,
                  6008,
                  "NEWGR enabling precision soft-capacity scenarios in low "
                  "congestion mode (overflow {}, hot edges {}, max ratio {:.2f}).",
                  baseline.metrics.overflow_edges,
                  baseline.metrics.near_capacity_edges,
                  baseline.metrics.max_usage_ratio);
  }

  auto make_random_def = [&](int seed,
                             float perturb_pct,
                             float critical_pct = 5.0f) {
    ScenarioDefinition def;
    def.name = "perturb-s" + std::to_string(seed)
               + "-p" + std::to_string(static_cast<int>(perturb_pct * 10.0f));
    def.pre_init = [this, seed, perturb_pct, critical_pct]() {
      grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
      grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
    };
    return def;
  };

  auto make_critical_sweep_def = [&](int seed, float critical_pct) {
    ScenarioDefinition def;
    def.name = "critical-s" + std::to_string(seed)
               + "-c" + std::to_string(static_cast<int>(critical_pct * 10.0f));
    def.pre_init = [this, seed, critical_pct]() {
      grouter_->setCapacitiesPerturbationPercentage(0.0f);
      grouter_->setPerturbationAmount(0);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
    };
    return def;
  };

  // Multi-seed sweep (SPRoute-inspired exploration) with a narrow perturbation
  // range to search for lower-wirelength minima while preserving routability.
  scenario_defs.push_back(make_random_def(11, 6.0f, 6.0f));
  scenario_defs.push_back(make_random_def(29, 4.0f, 5.0f));
  scenario_defs.push_back(make_random_def(7, 3.5f, 4.0f));
  scenario_defs.push_back(make_random_def(17, 3.0f, 4.0f));
  scenario_defs.push_back(make_random_def(23, 5.0f, 6.0f));
  scenario_defs.push_back(make_random_def(31, 2.5f, 3.0f));
  scenario_defs.push_back(make_random_def(37, 1.5f, 2.0f));
  scenario_defs.push_back(make_random_def(41, 4.5f, 5.0f));
  scenario_defs.push_back(make_random_def(43, 0.0f, 2.5f));
  scenario_defs.push_back(make_random_def(47, 0.5f, 2.0f));
  scenario_defs.push_back(make_random_def(53, 1.0f, 2.5f));
  scenario_defs.push_back(make_random_def(59, 2.0f, 3.0f));
  scenario_defs.push_back(make_random_def(61, 3.0f, 3.5f));
  scenario_defs.push_back(make_random_def(67, 0.0f, 1.0f));
  scenario_defs.push_back(make_critical_sweep_def(5, 0.0f));
  scenario_defs.push_back(make_critical_sweep_def(13, 2.0f));
  scenario_defs.push_back(make_critical_sweep_def(19, 4.0f));
  scenario_defs.push_back(make_critical_sweep_def(23, 8.0f));
  scenario_defs.push_back(make_critical_sweep_def(29, 12.0f));
  // Focused low-perturb exploration around historically strong seeds.
  scenario_defs.push_back(make_random_def(71, 0.0f, 8.0f));
  scenario_defs.push_back(make_random_def(73, 0.0f, 8.0f));
  scenario_defs.push_back(make_random_def(79, 0.0f, 8.0f));
  scenario_defs.push_back(make_random_def(83, 0.0f, 8.0f));
  scenario_defs.push_back(make_random_def(71, 0.5f, 6.0f));
  scenario_defs.push_back(make_random_def(79, 0.5f, 6.0f));
  scenario_defs.push_back(make_critical_sweep_def(23, 6.0f));
  scenario_defs.push_back(make_critical_sweep_def(23, 10.0f));
  scenario_defs.push_back(make_critical_sweep_def(67, 8.0f));
  scenario_defs.push_back(make_critical_sweep_def(19, 6.0f));

  if (compact_exploration_mode) {
    scenario_defs.clear();
    if (enable_softcap_scenarios) {
      scenario_defs.push_back(make_soft_config("cugr-prob-balanced",
                                               0.44f,
                                               0.89f,
                                               7.0f,
                                               0.47f,
                                               2,
                                               0.58f,
                                               0.70f,
                                               0.05f,
                                               0.0f,
                                               snapshot.seed,
                                               snapshot.critical_percentage));
      scenario_defs.push_back(make_soft_config("mild-softcap",
                                               0.58f,
                                               0.97f,
                                               4.5f,
                                               0.38f,
                                               1,
                                               0.75f,
                                               0.20f,
                                               0.08f,
                                               2.5f,
                                               5,
                                               8.0f));
    } else if (!normalized_rudy.empty()) {
      scenario_defs.push_back(make_soft_config("rudy-precision-direct",
                                               0.80f,
                                               1.02f,
                                               2.8f,
                                               0.24f,
                                               1,
                                               0.92f,
                                               0.10f,
                                               0.12f,
                                               0.0f,
                                               snapshot.seed,
                                               4.0f));
    }
    scenario_defs.push_back(make_random_def(23, 5.0f, 6.0f));
    scenario_defs.push_back(make_random_def(29, 4.0f, 5.0f));
    scenario_defs.push_back(make_random_def(67, 0.0f, 1.0f));
    scenario_defs.push_back(make_random_def(71, 0.0f, 8.0f));
    scenario_defs.push_back(make_critical_sweep_def(23, 8.0f));
    scenario_defs.push_back(make_critical_sweep_def(67, 8.0f));
  }

  for (const ScenarioDefinition& def : scenario_defs) {
    ScenarioResult result = run_scenario(def, snapshot);
    scenario_results.push_back(std::move(result));
  }

  // Two-stage exploration:
  // Stage 1: broad SPRoute-like diversification (seed/critical sweep).
  // Stage 2: CUGR/FastRoute-inspired exploitation around elite low-WL runs.
  const bool initial_overflow_free_sweep = std::all_of(
      scenario_results.begin(), scenario_results.end(), [](const auto& result) {
        return result.metrics.overflow_edges == 0;
      });
  if (initial_overflow_free_sweep) {
    std::vector<const ScenarioResult*> ranked;
    ranked.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      ranked.push_back(&result);
    }
    std::sort(
        ranked.begin(), ranked.end(), [](const ScenarioResult* lhs, const ScenarioResult* rhs) {
          if (lhs->metrics.wirelength_dbu != rhs->metrics.wirelength_dbu) {
            return lhs->metrics.wirelength_dbu < rhs->metrics.wirelength_dbu;
          }
          return lhs->metrics.via_count < rhs->metrics.via_count;
        });

    std::vector<ScenarioDefinition> refinement_defs;
    refinement_defs.reserve(12);
    std::set<std::string> scenario_names;
    for (const ScenarioResult& result : scenario_results) {
      scenario_names.insert(result.name);
    }

    auto add_refinement = [&](ScenarioDefinition def) {
      const std::size_t max_refinements = compact_exploration_mode ? 4 : 12;
      if (refinement_defs.size() >= max_refinements) {
        return;
      }
      if (scenario_names.insert(def.name).second) {
        refinement_defs.push_back(std::move(def));
      }
    };

    const int elite_count = compact_exploration_mode
                                ? std::min<int>(2, ranked.size())
                                : std::min<int>(4, ranked.size());
    for (int i = 0; i < elite_count; ++i) {
      const std::string& elite_name = ranked[i]->name;
      int seed = 0;
      float critical_pct = 0.0f;
      float perturb_pct = 0.0f;

      if (parseCriticalScenarioName(elite_name, seed, critical_pct)) {
        const float tighter_low
            = std::clamp(critical_pct - 2.0f, 0.0f, 20.0f);
        const float tighter_high
            = std::clamp(critical_pct + 2.0f, 0.0f, 20.0f);
        add_refinement(make_critical_sweep_def(seed, tighter_low));
        add_refinement(make_critical_sweep_def(seed, tighter_high));
        add_refinement(make_random_def(seed + 2, 0.0f, critical_pct));
        add_refinement(make_random_def(seed + 4, 0.0f, tighter_low));
      } else if (parsePerturbScenarioName(elite_name, seed, perturb_pct)) {
        const float reduced_perturb = std::max(0.0f, perturb_pct - 1.0f);
        add_refinement(make_random_def(seed, reduced_perturb, 8.0f));
        add_refinement(make_random_def(seed + 6,
                                       std::max(0.0f, reduced_perturb - 0.5f),
                                       8.0f));
        add_refinement(make_critical_sweep_def(seed, 8.0f));
      } else if (elite_name == "baseline") {
        add_refinement(make_random_def(snapshot.seed + 71, 0.0f, 6.0f));
        add_refinement(make_critical_sweep_def(snapshot.seed + 17, 8.0f));
      }
    }

    for (const ScenarioDefinition& def : refinement_defs) {
      ScenarioResult result = run_scenario(def, snapshot);
      scenario_results.push_back(std::move(result));
    }

    if (!refinement_defs.empty()) {
      logger_->info(
          GNR,
          6009,
          "NEWGR ran {} elite refinement scenarios for wirelength exploitation.",
          refinement_defs.size());
    }
  }

  // Cross-scenario net-level recombination:
  // pick each net route from the best wirelength scenarios (SPRoute-style
  // diversification) and combine into one hybrid guide set.
  if (initial_overflow_free_sweep && scenario_results.size() > 2) {
    // Keep scenario pointers stable while appending many hybrid candidates.
    scenario_results.reserve(scenario_results.size() + 256);

    std::vector<const ScenarioResult*> ranked;
    ranked.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      ranked.push_back(&result);
    }
    std::sort(
        ranked.begin(), ranked.end(), [](const ScenarioResult* lhs, const ScenarioResult* rhs) {
          if (lhs->metrics.wirelength_dbu != rhs->metrics.wirelength_dbu) {
            return lhs->metrics.wirelength_dbu < rhs->metrics.wirelength_dbu;
          }
          return lhs->metrics.via_count < rhs->metrics.via_count;
        });

    const std::size_t expected_net_count = ranked.front()->routes.size();

    std::vector<odb::dbNet*> hybrid_nets;
    hybrid_nets.reserve(expected_net_count);
    std::set<int> hybrid_net_ids;
    for (const ScenarioResult* source : ranked) {
      if (source == nullptr) {
        continue;
      }
      for (const auto& [route_net, route] : source->routes) {
        static_cast<void>(route);
        if (route_net == nullptr) {
          continue;
        }
        if (hybrid_net_ids.insert(route_net->getId()).second) {
          hybrid_nets.push_back(route_net);
        }
      }
    }

    auto append_hybrid = [&](const std::string& hybrid_name,
                             int source_count,
                             long via_weight,
                             double detour_weight,
                             double consensus_weight,
                             int logger_code,
                             double high_layer_weight = 0.0,
                             double layer_span_weight = 0.0,
                             double crowding_weight = 0.0) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }
        std::vector<std::pair<const GRoute*, RouteEdgeStats>> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          for (const auto& [edge_key, usage] : stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.emplace_back(&route_it->second, std::move(stats));
        }

        if (candidates.empty()) {
          continue;
        }

        const GRoute* best_route = nullptr;
        double best_cost = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_via = std::numeric_limits<int>::max();
        long best_high_layer = std::numeric_limits<long>::max();
        int best_layer_span = std::numeric_limits<int>::max();
        const bool use_structural_tie_break
            = high_layer_weight > 0.0 || layer_span_weight > 0.0
              || crowding_weight > 0.0;

        for (const auto& [route, stats] : candidates) {
          const long bbox_hpwl = getRouteBBoxHpwl(*route);
          const long detour = std::max(0L, stats.wirelength_dbu - bbox_hpwl);
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count = std::max(1, static_cast<int>(stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double consensus_penalty
              = std::max(0.0, static_cast<double>(source_count) - avg_vote);
          const double crowd_target
              = std::max(1.0, static_cast<double>(source_count) * 0.55);
          const double crowding_penalty
              = std::max(0.0, avg_vote - crowd_target);
          const double cost
              = static_cast<double>(stats.wirelength_dbu)
                + static_cast<double>(via_weight) * static_cast<double>(stats.via_count)
                + detour_weight * static_cast<double>(detour)
                + consensus_weight * static_cast<double>(tile_size) * consensus_penalty
                + high_layer_weight * static_cast<double>(stats.high_layer_dbu)
                + layer_span_weight * static_cast<double>(tile_size)
                      * static_cast<double>(stats.layer_span)
                + crowding_weight * static_cast<double>(tile_size)
                      * crowding_penalty;
          const long wl = stats.wirelength_dbu;
          const int vias = stats.via_count;
          const long high_layer = stats.high_layer_dbu;
          const int span = stats.layer_span;
          if (cost < best_cost || (cost == best_cost && wl < best_wl)
              || (cost == best_cost && wl == best_wl && vias < best_via)
              || (use_structural_tie_break && cost == best_cost
                  && wl == best_wl && vias == best_via
                  && high_layer < best_high_layer)
              || (use_structural_tie_break && cost == best_cost
                  && wl == best_wl && vias == best_via
                  && high_layer == best_high_layer && span < best_layer_span)) {
            best_cost = cost;
            best_wl = wl;
            best_via = vias;
            best_high_layer = high_layer;
            best_layer_span = span;
            best_route = route;
          }
        }

        if (best_route != nullptr) {
          hybrid_result.routes.emplace(db_net, *best_route);
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(GNR,
                      logger_code,
                      "NEWGR {} from top {} scenarios: wirelength {:.0f} um, "
                      "vias {}, routed nets {}/{}, detour {}, high-layer {}",
                      hybrid_name,
                      source_count,
                      hybrid_result.metrics.wirelength_um,
                      hybrid_result.metrics.via_count,
                      hybrid_result.routes.size(),
                      expected_net_count,
                      hybrid_result.metrics.detour_dbu,
                      hybrid_result.metrics.high_layer_dbu);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            6011,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_min_wl_hybrid = [&](const std::string& hybrid_name,
                                    int source_count,
                                    long via_weight,
                                    double detour_weight,
                                    double high_layer_weight,
                                    double layer_span_weight,
                                    int logger_code,
                                    double segment_weight = 0.0,
                                    double bend_weight = 0.0) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<std::pair<const GRoute*, RouteEdgeStats>> candidates;
        candidates.reserve(source_count);
        long min_wl = std::numeric_limits<long>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          min_wl = std::min(min_wl, stats.wirelength_dbu);
          candidates.emplace_back(&route_it->second, std::move(stats));
        }
        if (candidates.empty()) {
          continue;
        }

        const long wl_slack
            = std::max<long>(6L, static_cast<long>(std::ceil(tile_size * 0.30)));
        const GRoute* best_route = nullptr;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_via = std::numeric_limits<int>::max();

        auto evaluate_candidates = [&](bool enforce_wl_gate) {
          bool found = false;
          for (const auto& [route, stats] : candidates) {
            if (enforce_wl_gate && stats.wirelength_dbu > min_wl + wl_slack) {
              continue;
            }
            const long bbox_hpwl = getRouteBBoxHpwl(*route);
            const long detour = std::max(0L, stats.wirelength_dbu - bbox_hpwl);
            const double objective
                = static_cast<double>(stats.wirelength_dbu)
                  + static_cast<double>(via_weight)
                        * static_cast<double>(stats.via_count)
                  + detour_weight * static_cast<double>(detour)
                  + high_layer_weight * static_cast<double>(stats.high_layer_dbu)
                  + layer_span_weight * static_cast<double>(tile_size)
                        * static_cast<double>(stats.layer_span)
                  + segment_weight * static_cast<double>(stats.segment_count)
                  + bend_weight * static_cast<double>(stats.bend_count);
            if (objective < best_objective
                || (objective == best_objective
                    && stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && stats.wirelength_dbu == best_wl
                    && stats.via_count < best_via)) {
              found = true;
              best_objective = objective;
              best_wl = stats.wirelength_dbu;
              best_via = stats.via_count;
              best_route = route;
            }
          }
          return found;
        };

        const bool found_wl_gate = evaluate_candidates(true);
        if (!found_wl_gate) {
          static_cast<void>(evaluate_candidates(false));
        }
        if (best_route != nullptr) {
          hybrid_result.routes.emplace(db_net, *best_route);
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, detour {}, high-layer {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            hybrid_result.metrics.detour_dbu,
            hybrid_result.metrics.high_layer_dbu);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7101,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_wl_corridor_hybrid = [&](const std::string& hybrid_name,
                                         int source_count,
                                         long via_weight,
                                         int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      int strict_gate_hits = 0;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        long min_wl = std::numeric_limits<long>::max();
        long min_detour = std::numeric_limits<long>::max();
        long min_high_layer = std::numeric_limits<long>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          min_wl = std::min(min_wl, candidate.stats.wirelength_dbu);
          min_detour = std::min(min_detour, candidate.detour_dbu);
          min_high_layer
              = std::min(min_high_layer, candidate.stats.high_layer_dbu);
          candidates.push_back(std::move(candidate));
        }
        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.84;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(static_cast<double>(min_wl) * 0.0022)));
        const long detour_slack = std::max<long>(
            tile_size * 10L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_detour)) * 0.10)));
        const long high_layer_slack = std::max<long>(
            tile_size * 16L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_high_layer)) * 0.16)));

        const GRoute* best_route = nullptr;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_detour = std::numeric_limits<long>::max();

        auto evaluate_candidates = [&](bool strict_gate) {
          bool found = false;
          for (const CandidateRoute& candidate : candidates) {
            if (strict_gate
                && (candidate.stats.wirelength_dbu > min_wl + wl_slack
                    || candidate.detour_dbu > min_detour + detour_slack
                    || candidate.stats.high_layer_dbu
                           > min_high_layer + high_layer_slack)) {
              continue;
            }

            const double objective
                = static_cast<double>(candidate.stats.wirelength_dbu)
                  + static_cast<double>(via_weight)
                        * static_cast<double>(candidate.stats.via_count)
                  + 0.10 * static_cast<double>(candidate.detour_dbu)
                  + 0.004 * static_cast<double>(candidate.stats.high_layer_dbu)
                  + static_cast<double>(tile_size) * 0.18
                        * candidate.consensus_penalty
                  + 0.35 * static_cast<double>(candidate.stats.segment_count)
                  + 2.2 * static_cast<double>(candidate.stats.bend_count);
            if (objective < best_objective
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count < best_vias)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count == best_vias
                    && candidate.detour_dbu < best_detour)) {
              found = true;
              best_objective = objective;
              best_wl = candidate.stats.wirelength_dbu;
              best_vias = candidate.stats.via_count;
              best_detour = candidate.detour_dbu;
              best_route = candidate.route;
            }
          }
          return found;
        };

        const bool found_strict = evaluate_candidates(true);
        if (found_strict) {
          strict_gate_hits++;
        } else {
          static_cast<void>(evaluate_candidates(false));
        }

        if (best_route != nullptr) {
          hybrid_result.routes.emplace(db_net, *best_route);
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, strict-gate hits {}, detour {}, high-layer {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            strict_gate_hits,
            hybrid_result.metrics.detour_dbu,
            hybrid_result.metrics.high_layer_dbu);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7330,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_dr_stable_hybrid = [&](const std::string& hybrid_name,
                                       int source_count,
                                       long via_weight,
                                       int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        long min_wl = std::numeric_limits<long>::max();
        long min_detour = std::numeric_limits<long>::max();
        long min_high_layer = std::numeric_limits<long>::max();
        int min_layer_span = std::numeric_limits<int>::max();

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);

          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }

          min_wl = std::min(min_wl, candidate.stats.wirelength_dbu);
          min_detour = std::min(min_detour, candidate.detour_dbu);
          min_high_layer
              = std::min(min_high_layer, candidate.stats.high_layer_dbu);
          min_layer_span
              = std::min(min_layer_span, candidate.stats.layer_span);
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.82;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(static_cast<double>(min_wl) * 0.0035)));
        const long detour_slack = std::max<long>(
            tile_size * 18L,
            static_cast<long>(
                std::ceil(static_cast<double>(std::max(1L, min_detour)) * 0.08)));
        const long high_layer_slack = std::max<long>(
            tile_size * 24L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_high_layer)) * 0.14)));
        const int layer_span_slack = 1;

        const GRoute* best_route = nullptr;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_detour = std::numeric_limits<long>::max();

        auto evaluate_candidates = [&](bool strict_gate) {
          bool found = false;
          for (const CandidateRoute& candidate : candidates) {
            if (strict_gate) {
              if (candidate.stats.wirelength_dbu > min_wl + wl_slack
                  || candidate.detour_dbu > min_detour + detour_slack
                  || candidate.stats.high_layer_dbu > min_high_layer + high_layer_slack
                  || candidate.stats.layer_span > min_layer_span + layer_span_slack) {
                continue;
              }
            }

            const double objective
                = static_cast<double>(candidate.stats.wirelength_dbu)
                  + static_cast<double>(via_weight)
                        * static_cast<double>(candidate.stats.via_count)
                  + 0.32 * static_cast<double>(candidate.detour_dbu)
                  + 0.018 * static_cast<double>(candidate.stats.high_layer_dbu)
                  + static_cast<double>(tile_size) * 0.45
                        * static_cast<double>(candidate.stats.layer_span)
                  + 0.50 * static_cast<double>(candidate.stats.segment_count)
                  + 4.0 * static_cast<double>(candidate.stats.bend_count)
                  + static_cast<double>(tile_size) * 0.70
                        * candidate.consensus_penalty;
            if (objective < best_objective
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count < best_vias)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count == best_vias
                    && candidate.detour_dbu < best_detour)) {
              found = true;
              best_objective = objective;
              best_wl = candidate.stats.wirelength_dbu;
              best_vias = candidate.stats.via_count;
              best_detour = candidate.detour_dbu;
              best_route = candidate.route;
            }
          }
          return found;
        };

        const bool found_strict = evaluate_candidates(true);
        if (!found_strict) {
          static_cast<void>(evaluate_candidates(false));
        }
        if (best_route != nullptr) {
          hybrid_result.routes.emplace(db_net, *best_route);
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, detour {}, high-layer {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            hybrid_result.metrics.detour_dbu,
            hybrid_result.metrics.high_layer_dbu);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7308,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_length_adaptive_hybrid = [&](const std::string& hybrid_name,
                                             int source_count,
                                             long via_weight,
                                             int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      int long_net_picks = 0;
      int short_net_picks = 0;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        long min_wl = std::numeric_limits<long>::max();
        long min_detour = std::numeric_limits<long>::max();

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          min_wl = std::min(min_wl, candidate.stats.wirelength_dbu);
          min_detour = std::min(min_detour, candidate.detour_dbu);
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.78;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        // Mix and match by net scale:
        // long/high-detour nets use aggressive shortest-path scoring, while
        // shorter nets use DR-stable scoring with stronger structural penalties.
        const bool long_net = min_wl >= tile_size * 40L || min_detour >= tile_size * 14L;
        const long wl_slack = std::max<long>(
            long_net ? tile_size * 4L : tile_size * 2L,
            static_cast<long>(
                std::ceil(static_cast<double>(min_wl)
                          * (long_net ? 0.0048 : 0.0032))));
        const long detour_slack = std::max<long>(
            long_net ? tile_size * 20L : tile_size * 12L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_detour))
                * (long_net ? 0.22 : 0.14))));

        const GRoute* best_route = nullptr;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_detour = std::numeric_limits<long>::max();

        auto evaluate = [&](bool strict_gate) {
          bool found = false;
          for (const CandidateRoute& candidate : candidates) {
            if (strict_gate && (candidate.stats.wirelength_dbu > min_wl + wl_slack
                                || candidate.detour_dbu > min_detour + detour_slack)) {
              continue;
            }

            const double objective = long_net
                                         ? static_cast<double>(candidate.stats.wirelength_dbu)
                                               + static_cast<double>(via_weight)
                                                     * static_cast<double>(candidate.stats.via_count)
                                               + 0.12
                                                     * static_cast<double>(candidate.detour_dbu)
                                               + 0.0035
                                                     * static_cast<double>(candidate.stats.high_layer_dbu)
                                               + static_cast<double>(tile_size) * 0.14
                                                     * static_cast<double>(candidate.stats.layer_span)
                                               + 0.40
                                                     * static_cast<double>(candidate.stats.segment_count)
                                               + 2.8
                                                     * static_cast<double>(candidate.stats.bend_count)
                                               + static_cast<double>(tile_size) * 0.22
                                                     * candidate.consensus_penalty
                                         : static_cast<double>(candidate.stats.wirelength_dbu)
                                               + static_cast<double>(via_weight)
                                                     * static_cast<double>(candidate.stats.via_count)
                                               + 0.38
                                                     * static_cast<double>(candidate.detour_dbu)
                                               + 0.022
                                                     * static_cast<double>(candidate.stats.high_layer_dbu)
                                               + static_cast<double>(tile_size) * 0.72
                                                     * static_cast<double>(candidate.stats.layer_span)
                                               + 1.25
                                                     * static_cast<double>(candidate.stats.segment_count)
                                               + 8.0
                                                     * static_cast<double>(candidate.stats.bend_count)
                                               + static_cast<double>(tile_size) * 0.95
                                                     * candidate.consensus_penalty;
            if (objective < best_objective
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count < best_vias)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count == best_vias
                    && candidate.detour_dbu < best_detour)) {
              found = true;
              best_objective = objective;
              best_wl = candidate.stats.wirelength_dbu;
              best_vias = candidate.stats.via_count;
              best_detour = candidate.detour_dbu;
              best_route = candidate.route;
            }
          }
          return found;
        };

        const bool found_strict = evaluate(true);
        if (!found_strict) {
          static_cast<void>(evaluate(false));
        }
        if (best_route != nullptr) {
          hybrid_result.routes.emplace(db_net, *best_route);
          if (long_net) {
            long_net_picks++;
          } else {
            short_net_picks++;
          }
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, long-net picks {}, short-net picks {}, "
            "detour {}, high-layer {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            long_net_picks,
            short_net_picks,
            hybrid_result.metrics.detour_dbu,
            hybrid_result.metrics.high_layer_dbu);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7311,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_feedback_wl_hybrid = [&](const std::string& hybrid_name,
                                         int source_count,
                                         long via_weight,
                                         int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      double total_crowding_penalty = 0.0;
      double total_detour_penalty = 0.0;

      std::vector<odb::dbNet*> ordered_nets = hybrid_nets;
      std::unordered_map<odb::dbNet*, long> best_wl_per_net;
      best_wl_per_net.reserve(ordered_nets.size());
      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }
        long best_wl = std::numeric_limits<long>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          const RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          best_wl = std::min(best_wl, stats.wirelength_dbu);
        }
        best_wl_per_net.emplace(db_net, best_wl);
      }
      std::sort(ordered_nets.begin(),
                ordered_nets.end(),
                [&](odb::dbNet* lhs, odb::dbNet* rhs) {
                  const auto lhs_it = best_wl_per_net.find(lhs);
                  const auto rhs_it = best_wl_per_net.find(rhs);
                  const long lhs_wl = lhs_it != best_wl_per_net.end()
                                          ? lhs_it->second
                                          : std::numeric_limits<long>::max();
                  const long rhs_wl = rhs_it != best_wl_per_net.end()
                                          ? rhs_it->second
                                          : std::numeric_limits<long>::max();
                  return lhs_wl > rhs_wl;
                });

      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);
        long min_wl = std::numeric_limits<long>::max();
        long min_detour = std::numeric_limits<long>::max();
        long min_high_layer = std::numeric_limits<long>::max();
        int min_layer_span = std::numeric_limits<int>::max();

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          min_wl = std::min(min_wl, candidate.stats.wirelength_dbu);
          min_detour = std::min(min_detour, candidate.detour_dbu);
          min_high_layer
              = std::min(min_high_layer, candidate.stats.high_layer_dbu);
          min_layer_span
              = std::min(min_layer_span, candidate.stats.layer_span);
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.74;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(static_cast<double>(min_wl) * 0.0028)));
        const long detour_slack = std::max<long>(
            tile_size * 14L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_detour)) * 0.12)));
        const long high_layer_slack = std::max<long>(
            tile_size * 22L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, min_high_layer)) * 0.18)));
        const int layer_span_slack = 1;

        const GRoute* best_route = nullptr;
        EdgeCountMap best_edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        double best_crowding_delta = 0.0;
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_detour = std::numeric_limits<long>::max();

        auto evaluate = [&](bool strict_gate) {
          bool found = false;
          for (const CandidateRoute& candidate : candidates) {
            if (strict_gate) {
              if (candidate.stats.wirelength_dbu > min_wl + wl_slack
                  || candidate.detour_dbu > min_detour + detour_slack
                  || candidate.stats.high_layer_dbu > min_high_layer + high_layer_slack
                  || candidate.stats.layer_span > min_layer_span + layer_span_slack) {
                continue;
              }
            }

            double crowding_delta = 0.0;
            for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
              const auto usage_it = edge_usage.find(edge_key);
              const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
              const int before_excess = std::max(0, usage - 2);
              const int after_excess = std::max(0, usage + add_usage - 2);
              crowding_delta += static_cast<double>(
                  after_excess * after_excess - before_excess * before_excess);
            }

            const double objective
                = static_cast<double>(candidate.stats.wirelength_dbu)
                  + static_cast<double>(via_weight)
                        * static_cast<double>(candidate.stats.via_count)
                  + 0.18 * static_cast<double>(candidate.detour_dbu)
                  + 0.0065 * static_cast<double>(candidate.stats.high_layer_dbu)
                  + static_cast<double>(tile_size) * 0.20
                        * static_cast<double>(candidate.stats.layer_span)
                  + static_cast<double>(tile_size) * 0.52
                        * candidate.consensus_penalty
                  + static_cast<double>(tile_size) * 0.08 * crowding_delta;
            if (objective < best_objective
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count < best_vias)
                || (objective == best_objective
                    && candidate.stats.wirelength_dbu == best_wl
                    && candidate.stats.via_count == best_vias
                    && candidate.detour_dbu < best_detour)) {
              found = true;
              best_objective = objective;
              best_crowding_delta = crowding_delta;
              best_wl = candidate.stats.wirelength_dbu;
              best_vias = candidate.stats.via_count;
              best_detour = candidate.detour_dbu;
              best_route = candidate.route;
              best_edge_counts = candidate.stats.edge_counts;
            }
          }
          return found;
        };

        const bool found_strict = evaluate(true);
        if (!found_strict) {
          static_cast<void>(evaluate(false));
        }
        if (best_route == nullptr) {
          continue;
        }

        hybrid_result.routes.emplace(db_net, *best_route);
        total_crowding_penalty += best_crowding_delta;
        total_detour_penalty += static_cast<double>(best_detour);
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, crowding-penalty {:.0f}, detour-penalty {:.0f}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            total_crowding_penalty,
            total_detour_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7314,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_via_floor_wl_hybrid = [&](const std::string& hybrid_name,
                                          int source_count,
                                          int per_net_via_drop_limit,
                                          int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long via_drop_budget
          = std::max<long>(220L,
                           static_cast<long>(
                               std::ceil(static_cast<double>(expected_net_count)
                                         * 0.018)));
      long consumed_via_drop = 0;
      int budget_clamped_nets = 0;

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.76;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const CandidateRoute* anchor = nullptr;
        double best_anchor_objective = std::numeric_limits<double>::max();
        for (const CandidateRoute& candidate : candidates) {
          const double anchor_objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.20 * static_cast<double>(candidate.detour_dbu)
                + static_cast<double>(tile_size) * 0.45
                      * candidate.consensus_penalty;
          if (anchor_objective < best_anchor_objective) {
            best_anchor_objective = anchor_objective;
            anchor = &candidate;
          }
        }
        if (anchor == nullptr) {
          continue;
        }

        const int min_allowed_vias
            = std::max(0, anchor->stats.via_count - per_net_via_drop_limit);
        const long wl_slack
            = std::max<long>(tile_size * 2L,
                             static_cast<long>(std::ceil(
                                 static_cast<double>(anchor->stats.wirelength_dbu)
                                 * 0.0018)));
        const long detour_slack = std::max<long>(
            tile_size * 12L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor->detour_dbu)) * 0.10)));

        const CandidateRoute* best_candidate = anchor;
        double best_objective = std::numeric_limits<double>::max();
        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count < min_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor->stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor->detour_dbu + detour_slack) {
            continue;
          }

          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.08 * static_cast<double>(candidate.detour_dbu)
                + 0.004 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.16
                      * candidate.consensus_penalty
                + 0.40 * static_cast<double>(candidate.stats.segment_count)
                + 2.0 * static_cast<double>(candidate.stats.bend_count);
          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         < best_candidate->stats.wirelength_dbu)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         == best_candidate->stats.wirelength_dbu
                  && candidate.stats.via_count < best_candidate->stats.via_count)) {
            best_objective = objective;
            best_candidate = &candidate;
          }
        }

        if (best_candidate == nullptr) {
          continue;
        }

        const long via_drop = std::max(
            0, anchor->stats.via_count - best_candidate->stats.via_count);
        if (best_candidate != anchor
            && consumed_via_drop + via_drop > via_drop_budget) {
          best_candidate = anchor;
          budget_clamped_nets++;
        } else {
          consumed_via_drop += via_drop;
        }
        hybrid_result.routes.emplace(db_net, *best_candidate->route);
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, via-drop {}/{}, budget-clamped nets {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            consumed_via_drop,
            via_drop_budget,
            budget_clamped_nets);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7317,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_via_band_wl_hybrid = [&](const std::string& hybrid_name,
                                         int source_count,
                                         int per_net_via_drop_limit,
                                         int per_net_via_rise_limit,
                                         int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      double total_crowding_penalty = 0.0;
      int fallback_to_anchor_count = 0;

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.78;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const CandidateRoute* anchor = nullptr;
        double best_anchor_objective = std::numeric_limits<double>::max();
        for (const CandidateRoute& candidate : candidates) {
          const double anchor_objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.20 * static_cast<double>(candidate.detour_dbu)
                + 0.0045 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.32
                      * candidate.consensus_penalty;
          if (anchor_objective < best_anchor_objective) {
            best_anchor_objective = anchor_objective;
            anchor = &candidate;
          }
        }
        if (anchor == nullptr) {
          continue;
        }

        const int min_allowed_vias
            = std::max(0, anchor->stats.via_count - per_net_via_drop_limit);
        const int max_allowed_vias = anchor->stats.via_count + per_net_via_rise_limit;
        const long wl_slack
            = std::max<long>(tile_size * 2L,
                             static_cast<long>(std::ceil(
                                 static_cast<double>(anchor->stats.wirelength_dbu)
                                 * 0.0020)));
        const long detour_slack = std::max<long>(
            tile_size * 14L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor->detour_dbu)) * 0.12)));
        const long high_layer_slack = std::max<long>(
            tile_size * 24L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor->stats.high_layer_dbu))
                * 0.18)));

        const CandidateRoute* best_candidate = nullptr;
        EdgeCountMap best_edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        double best_crowding_delta = 0.0;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count < min_allowed_vias
              || candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor->stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor->detour_dbu + detour_slack
              || candidate.stats.high_layer_dbu
                     > anchor->stats.high_layer_dbu + high_layer_slack) {
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const double via_band_penalty = static_cast<double>(
              std::llabs(static_cast<long>(candidate.stats.via_count)
                         - static_cast<long>(anchor->stats.via_count)));
          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.14 * static_cast<double>(candidate.detour_dbu)
                + 0.006 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.24
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.12 * crowding_delta
                + 0.45 * static_cast<double>(candidate.stats.segment_count)
                + 2.6 * static_cast<double>(candidate.stats.bend_count)
                + 0.80 * via_band_penalty;
          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         < best_candidate->stats.wirelength_dbu)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         == best_candidate->stats.wirelength_dbu
                  && candidate.stats.via_count < best_candidate->stats.via_count)) {
            best_objective = objective;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        if (best_candidate == nullptr) {
          best_candidate = anchor;
          best_edge_counts = anchor->stats.edge_counts;
          fallback_to_anchor_count++;
        }
        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, crowding-penalty {:.0f}, anchor-fallbacks {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            total_crowding_penalty,
            fallback_to_anchor_count);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7318,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_layer_lifted_wl_hybrid = [&](const std::string& hybrid_name,
                                             int source_count,
                                             int per_net_via_rise_limit,
                                             int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      int long_net_picks = 0;
      int lift_guard_rejects = 0;
      double total_crowding_penalty = 0.0;

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        if (candidates.empty()) {
          continue;
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target = static_cast<double>(source_count) * 0.76;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const CandidateRoute* anchor = nullptr;
        double best_anchor_objective = std::numeric_limits<double>::max();
        for (const CandidateRoute& candidate : candidates) {
          const double anchor_objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.20 * static_cast<double>(candidate.detour_dbu)
                + 0.003 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.30
                      * candidate.consensus_penalty;
          if (anchor_objective < best_anchor_objective) {
            best_anchor_objective = anchor_objective;
            anchor = &candidate;
          }
        }
        if (anchor == nullptr) {
          continue;
        }

        const bool long_net = anchor->stats.wirelength_dbu >= tile_size * 34L
                              || anchor->detour_dbu >= tile_size * 12L;
        const int max_allowed_vias
            = anchor->stats.via_count
              + (long_net ? std::max(per_net_via_rise_limit, 6)
                          : std::max(2, per_net_via_rise_limit / 2));
        const long wl_slack = std::max<long>(
            long_net ? tile_size * 4L : tile_size * 2L,
            static_cast<long>(std::ceil(
                static_cast<double>(anchor->stats.wirelength_dbu)
                * (long_net ? 0.0055 : 0.0035))));
        const long detour_slack = std::max<long>(
            long_net ? tile_size * 18L : tile_size * 12L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor->detour_dbu))
                * (long_net ? 0.20 : 0.12))));
        const long high_layer_floor = std::max<long>(
            0L,
            anchor->stats.high_layer_dbu
                - std::max<long>(
                    long_net ? tile_size * 42L : tile_size * 18L,
                    static_cast<long>(std::ceil(
                        static_cast<double>(std::max(1L, anchor->stats.high_layer_dbu))
                        * (long_net ? 0.20 : 0.35)))));

        const CandidateRoute* best_candidate = anchor;
        EdgeCountMap best_edge_counts = anchor->stats.edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        double best_crowding_delta = 0.0;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor->stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor->detour_dbu + detour_slack) {
            continue;
          }

          if (long_net && candidate.stats.high_layer_dbu < high_layer_floor) {
            if (candidate.stats.wirelength_dbu <= anchor->stats.wirelength_dbu) {
              lift_guard_rejects++;
            }
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const long layer_lift
              = std::max(0L, candidate.stats.high_layer_dbu - high_layer_floor);
          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.14 * static_cast<double>(candidate.detour_dbu)
                + 0.55 * static_cast<double>(candidate.stats.via_count)
                + static_cast<double>(tile_size) * 0.25
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.10 * crowding_delta
                + 0.45 * static_cast<double>(candidate.stats.segment_count)
                + 2.6 * static_cast<double>(candidate.stats.bend_count)
                - 0.010 * static_cast<double>(layer_lift);
          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         < best_candidate->stats.wirelength_dbu)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu
                         == best_candidate->stats.wirelength_dbu
                  && candidate.stats.via_count < best_candidate->stats.via_count)) {
            best_objective = objective;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        if (long_net) {
          long_net_picks++;
        }
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, long-net picks {}, lift-guard rejects {}, "
            "crowding-penalty {:.0f}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            long_net_picks,
            lift_guard_rejects,
            total_crowding_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7327,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_anchor_locked_wl_hybrid = [&](const std::string& hybrid_name,
                                              const std::string& anchor_name,
                                              int source_count,
                                              int per_net_via_rise_limit,
                                              int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      const ScenarioResult* anchor_result = nullptr;
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == anchor_name) {
          anchor_result = &result;
          break;
        }
      }
      if (anchor_result == nullptr || anchor_result->routes.empty()) {
        return;
      }

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      int improved_nets = 0;
      int fallback_nets = 0;
      double total_crowding_penalty = 0.0;

      for (const auto& [db_net, anchor_route] : anchor_result->routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(anchor_route, grouter_->grid_);
        const long anchor_bbox_hpwl = getRouteBBoxHpwl(anchor_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - anchor_bbox_hpwl);

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count + 1);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        CandidateRoute anchor_candidate;
        anchor_candidate.route = &anchor_route;
        anchor_candidate.stats = anchor_stats;
        anchor_candidate.detour_dbu = anchor_detour;
        candidates.push_back(anchor_candidate);
        for (const auto& [edge_key, usage] : anchor_stats.edge_counts) {
          if (usage > 0) {
            edge_frequency[edge_key] += 1;
          }
        }

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          if (&route_it->second == &anchor_route) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target
              = static_cast<double>(std::max(2, source_count + 1)) * 0.74;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const int max_allowed_vias = anchor_stats.via_count + per_net_via_rise_limit;
        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(
                static_cast<double>(anchor_stats.wirelength_dbu) * 0.0022)));
        const long detour_slack = std::max<long>(
            tile_size * 16L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_detour)) * 0.12)));
        const long high_layer_slack = std::max<long>(
            tile_size * 26L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * 0.20)));

        const CandidateRoute* best_candidate = nullptr;
        EdgeCountMap best_edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        double best_crowding_delta = 0.0;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count < anchor_stats.via_count
              || candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor_stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor_detour + detour_slack
              || candidate.stats.high_layer_dbu
                     > anchor_stats.high_layer_dbu + high_layer_slack) {
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.12 * static_cast<double>(candidate.detour_dbu)
                + 0.0055 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.22
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.10 * crowding_delta
                + 0.35 * static_cast<double>(candidate.stats.segment_count)
                + 2.4 * static_cast<double>(candidate.stats.bend_count);
          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu < best_wl)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu == best_wl
                  && candidate.stats.via_count < best_vias)) {
            best_objective = objective;
            best_wl = candidate.stats.wirelength_dbu;
            best_vias = candidate.stats.via_count;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        if (best_candidate == nullptr) {
          best_candidate = &candidates.front();
          best_edge_counts = candidates.front().stats.edge_counts;
          fallback_nets++;
        } else if (best_candidate->stats.wirelength_dbu
                       < anchor_stats.wirelength_dbu
                   || (best_candidate->stats.wirelength_dbu
                           == anchor_stats.wirelength_dbu
                       && best_candidate->stats.via_count
                              != anchor_stats.via_count)) {
          improved_nets++;
        }

        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t expected_net_count = anchor_result->routes.size();
      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} anchored on {}: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, improved nets {}, anchor-fallbacks {}, "
            "crowding-penalty {:.0f}",
            hybrid_name,
            anchor_name,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            improved_nets,
            fallback_nets,
            total_crowding_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7320,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_layer_floor_wl_hybrid = [&](const std::string& hybrid_name,
                                            const std::string& anchor_name,
                                            int source_count,
                                            int per_net_via_drop_limit,
                                            int per_net_via_rise_limit,
                                            int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      const ScenarioResult* anchor_result = nullptr;
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == anchor_name) {
          anchor_result = &result;
          break;
        }
      }
      if (anchor_result == nullptr || anchor_result->routes.empty()) {
        return;
      }

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      int improved_nets = 0;
      int fallback_nets = 0;
      int layer_floor_rejects = 0;
      double total_crowding_penalty = 0.0;

      for (const auto& [db_net, anchor_route] : anchor_result->routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(anchor_route, grouter_->grid_);
        const long anchor_bbox_hpwl = getRouteBBoxHpwl(anchor_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - anchor_bbox_hpwl);

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count + 1);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        CandidateRoute anchor_candidate;
        anchor_candidate.route = &anchor_route;
        anchor_candidate.stats = anchor_stats;
        anchor_candidate.detour_dbu = anchor_detour;
        candidates.push_back(anchor_candidate);
        for (const auto& [edge_key, usage] : anchor_stats.edge_counts) {
          if (usage > 0) {
            edge_frequency[edge_key] += 1;
          }
        }

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          if (&route_it->second == &anchor_route) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target
              = static_cast<double>(std::max(2, source_count + 1)) * 0.76;
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const int min_allowed_vias
            = std::max(0, anchor_stats.via_count - per_net_via_drop_limit);
        const int max_allowed_vias = anchor_stats.via_count + per_net_via_rise_limit;
        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(
                static_cast<double>(anchor_stats.wirelength_dbu) * 0.0018)));
        const long detour_slack = std::max<long>(
            tile_size * 14L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_detour)) * 0.11)));
        const long high_layer_drop_guard = std::max<long>(
            tile_size * 12L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * 0.12)));
        const long high_layer_rise_guard = std::max<long>(
            tile_size * 24L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * 0.20)));
        const long min_allowed_high_layer
            = std::max(0L, anchor_stats.high_layer_dbu - high_layer_drop_guard);
        const long max_allowed_high_layer
            = anchor_stats.high_layer_dbu + high_layer_rise_guard;

        const CandidateRoute* best_candidate = nullptr;
        EdgeCountMap best_edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_high_layer = std::numeric_limits<long>::max();
        double best_crowding_delta = 0.0;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count < min_allowed_vias
              || candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor_stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor_detour + detour_slack
              || candidate.stats.high_layer_dbu > max_allowed_high_layer) {
            continue;
          }
          if (candidate.stats.high_layer_dbu < min_allowed_high_layer) {
            layer_floor_rejects++;
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const long high_layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - candidate.stats.high_layer_dbu);
          const long via_delta_abs
              = std::llabs(static_cast<long>(candidate.stats.via_count)
                           - static_cast<long>(anchor_stats.via_count));
          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + 0.10 * static_cast<double>(candidate.detour_dbu)
                + 0.0045 * static_cast<double>(candidate.stats.high_layer_dbu)
                + static_cast<double>(tile_size) * 0.20
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.09 * crowding_delta
                + 0.90 * static_cast<double>(via_delta_abs)
                + 0.010 * static_cast<double>(high_layer_drop)
                + 0.55 * static_cast<double>(candidate.stats.segment_count)
                + 3.2 * static_cast<double>(candidate.stats.bend_count);
          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu < best_wl)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu == best_wl
                  && candidate.stats.via_count < best_vias)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu == best_wl
                  && candidate.stats.via_count == best_vias
                  && candidate.stats.high_layer_dbu < best_high_layer)) {
            best_objective = objective;
            best_wl = candidate.stats.wirelength_dbu;
            best_vias = candidate.stats.via_count;
            best_high_layer = candidate.stats.high_layer_dbu;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        if (best_candidate == nullptr) {
          best_candidate = &candidates.front();
          best_edge_counts = candidates.front().stats.edge_counts;
          fallback_nets++;
        } else if (best_candidate->stats.wirelength_dbu
                       < anchor_stats.wirelength_dbu
                   || (best_candidate->stats.wirelength_dbu
                           == anchor_stats.wirelength_dbu
                       && best_candidate->stats.via_count
                              < anchor_stats.via_count)) {
          improved_nets++;
        }

        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t expected_net_count = anchor_result->routes.size();
      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} anchored on {}: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, improved nets {}, anchor-fallbacks {}, "
            "layer-floor rejects {}, crowding-penalty {:.0f}",
            hybrid_name,
            anchor_name,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            improved_nets,
            fallback_nets,
            layer_floor_rejects,
            total_crowding_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7342,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_drt_elastic_wl_hybrid = [&](const std::string& hybrid_name,
                                            const std::string& anchor_name,
                                            int source_count,
                                            int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      const ScenarioResult* anchor_result = nullptr;
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == anchor_name) {
          anchor_result = &result;
          break;
        }
      }
      if (anchor_result == nullptr || anchor_result->routes.empty()) {
        return;
      }

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      int improved_nets = 0;
      int fallback_nets = 0;
      int long_net_swaps = 0;
      int low_layer_rejects = 0;
      double total_crowding_penalty = 0.0;

      for (const auto& [db_net, anchor_route] : anchor_result->routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(anchor_route, grouter_->grid_);
        const long anchor_bbox_hpwl = getRouteBBoxHpwl(anchor_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - anchor_bbox_hpwl);
        const bool long_net
            = anchor_bbox_hpwl >= static_cast<long>(tile_size) * 18L
              || anchor_stats.wirelength_dbu >= static_cast<long>(tile_size) * 22L
              || anchor_stats.segment_count >= 8;

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count + 1);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        CandidateRoute anchor_candidate;
        anchor_candidate.route = &anchor_route;
        anchor_candidate.stats = anchor_stats;
        anchor_candidate.detour_dbu = anchor_detour;
        candidates.push_back(anchor_candidate);
        for (const auto& [edge_key, usage] : anchor_stats.edge_counts) {
          if (usage > 0) {
            edge_frequency[edge_key] += 1;
          }
        }

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          if (&route_it->second == &anchor_route) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target
              = static_cast<double>(std::max(2, source_count + 1))
                * (long_net ? 0.72 : 0.78);
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const int min_allowed_vias = std::max(
            0,
            anchor_stats.via_count - (long_net ? 12 : 8));
        const int max_allowed_vias
            = anchor_stats.via_count + (long_net ? 9 : 5);
        const long wl_slack = std::max<long>(
            tile_size * (long_net ? 4L : 2L),
            static_cast<long>(std::ceil(
                static_cast<double>(anchor_stats.wirelength_dbu)
                * (long_net ? 0.0030 : 0.0018))));
        const long detour_slack = std::max<long>(
            tile_size * (long_net ? 18L : 12L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_detour))
                * (long_net ? 0.18 : 0.12))));
        const long high_layer_drop_guard = std::max<long>(
            tile_size * (long_net ? 12L : 8L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * (long_net ? 0.10 : 0.24))));
        const long high_layer_rise_guard = std::max<long>(
            tile_size * (long_net ? 26L : 18L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * (long_net ? 0.24 : 0.16))));
        const long min_allowed_high_layer
            = std::max(0L, anchor_stats.high_layer_dbu - high_layer_drop_guard);
        const long max_allowed_high_layer
            = anchor_stats.high_layer_dbu + high_layer_rise_guard;

        const CandidateRoute* best_candidate = nullptr;
        EdgeCountMap best_edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_high_layer = std::numeric_limits<long>::min();
        double best_crowding_delta = 0.0;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count < min_allowed_vias
              || candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor_stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor_detour + detour_slack
              || candidate.stats.high_layer_dbu > max_allowed_high_layer) {
            continue;
          }
          if (candidate.stats.high_layer_dbu < min_allowed_high_layer) {
            low_layer_rejects++;
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const long low_layer_deficit = std::max(
              0L, min_allowed_high_layer - candidate.stats.high_layer_dbu);
          const long high_layer_excess = std::max(
              0L, candidate.stats.high_layer_dbu - anchor_stats.high_layer_dbu);
          const long via_delta_abs
              = std::llabs(static_cast<long>(candidate.stats.via_count)
                           - static_cast<long>(anchor_stats.via_count));
          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + (long_net ? 0.12 : 0.16)
                      * static_cast<double>(candidate.detour_dbu)
                + static_cast<double>(tile_size) * 0.17
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.08 * crowding_delta
                + 0.55 * static_cast<double>(via_delta_abs)
                + 0.020 * static_cast<double>(low_layer_deficit)
                + (long_net ? 0.0030 : 0.0045)
                      * static_cast<double>(high_layer_excess)
                + 0.45 * static_cast<double>(candidate.stats.segment_count)
                + 2.8 * static_cast<double>(candidate.stats.bend_count);
          const long wl = candidate.stats.wirelength_dbu;
          const int vias = candidate.stats.via_count;
          const long high_layer = candidate.stats.high_layer_dbu;
          if (objective < best_objective
              || (objective == best_objective && wl < best_wl)
              || (objective == best_objective && wl == best_wl && vias < best_vias)
              || (objective == best_objective && wl == best_wl && vias == best_vias
                  && ((long_net && high_layer > best_high_layer)
                      || (!long_net && high_layer < best_high_layer)))) {
            best_objective = objective;
            best_wl = wl;
            best_vias = vias;
            best_high_layer = high_layer;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        if (best_candidate == nullptr) {
          best_candidate = &candidates.front();
          best_edge_counts = candidates.front().stats.edge_counts;
          fallback_nets++;
        } else if (best_candidate->route != &anchor_route) {
          const bool wl_improves = best_candidate->stats.wirelength_dbu
                                   < anchor_stats.wirelength_dbu;
          const bool detour_trade_improves
              = best_candidate->stats.wirelength_dbu
                    <= anchor_stats.wirelength_dbu + tile_size * 2L
                && best_candidate->detour_dbu
                       + tile_size * 6L < anchor_detour;
          if (wl_improves || detour_trade_improves) {
            improved_nets++;
          }
          if (long_net
              && best_candidate->stats.high_layer_dbu >= anchor_stats.high_layer_dbu) {
            long_net_swaps++;
          }
        }

        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t expected_net_count = anchor_result->routes.size();
      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} anchored on {}: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, improved nets {}, long-net swaps {}, "
            "anchor-fallbacks {}, low-layer rejects {}, crowding-penalty {:.0f}",
            hybrid_name,
            anchor_name,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            improved_nets,
            long_net_swaps,
            fallback_nets,
            low_layer_rejects,
            total_crowding_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7350,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_budgeted_lift_wl_hybrid = [&](const std::string& hybrid_name,
                                              const std::string& anchor_name,
                                              int source_count,
                                              int logger_code) {
      source_count
          = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      const ScenarioResult* anchor_result = nullptr;
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == anchor_name) {
          anchor_result = &result;
          break;
        }
      }
      if (anchor_result == nullptr || anchor_result->routes.empty()) {
        return;
      }

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
        long detour_dbu = 0;
        double consensus_penalty = 0.0;
      };

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      std::unordered_map<uint64_t, int> edge_usage;
      edge_usage.reserve(1 << 20);
      const long max_total_via_rise = std::max<long>(
          80L,
          static_cast<long>(std::ceil(
              static_cast<double>(anchor_result->metrics.via_count) * 0.0006)));
      long total_via_rise = 0;
      int improved_nets = 0;
      int fallback_nets = 0;
      int budget_clamped_nets = 0;
      double total_crowding_penalty = 0.0;

      std::vector<odb::dbNet*> ordered_nets;
      ordered_nets.reserve(anchor_result->routes.size());
      std::unordered_map<odb::dbNet*, long> anchor_wl_map;
      anchor_wl_map.reserve(anchor_result->routes.size());
      for (const auto& [db_net, anchor_route] : anchor_result->routes) {
        if (db_net == nullptr) {
          continue;
        }
        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(anchor_route, grouter_->grid_);
        anchor_wl_map.emplace(db_net, anchor_stats.wirelength_dbu);
        ordered_nets.push_back(db_net);
      }
      std::sort(ordered_nets.begin(),
                ordered_nets.end(),
                [&](odb::dbNet* lhs, odb::dbNet* rhs) {
                  const long lhs_wl = anchor_wl_map.count(lhs) > 0
                                          ? anchor_wl_map[lhs]
                                          : std::numeric_limits<long>::max();
                  const long rhs_wl = anchor_wl_map.count(rhs) > 0
                                          ? anchor_wl_map[rhs]
                                          : std::numeric_limits<long>::max();
                  return lhs_wl > rhs_wl;
                });

      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }
        const auto anchor_it = anchor_result->routes.find(db_net);
        if (anchor_it == anchor_result->routes.end()) {
          continue;
        }
        const GRoute& anchor_route = anchor_it->second;
        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(anchor_route, grouter_->grid_);
        const long anchor_bbox_hpwl = getRouteBBoxHpwl(anchor_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - anchor_bbox_hpwl);
        const bool long_net
            = anchor_bbox_hpwl >= static_cast<long>(tile_size) * 20L
              || anchor_stats.wirelength_dbu >= static_cast<long>(tile_size) * 24L
              || anchor_stats.segment_count >= 9;

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count + 1);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        CandidateRoute anchor_candidate;
        anchor_candidate.route = &anchor_route;
        anchor_candidate.stats = anchor_stats;
        anchor_candidate.detour_dbu = anchor_detour;
        candidates.push_back(anchor_candidate);
        for (const auto& [edge_key, usage] : anchor_stats.edge_counts) {
          if (usage > 0) {
            edge_frequency[edge_key] += 1;
          }
        }

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          if (&route_it->second == &anchor_route) {
            continue;
          }

          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(route_it->second);
          candidate.detour_dbu
              = std::max(0L, candidate.stats.wirelength_dbu - bbox_hpwl);
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.push_back(std::move(candidate));
        }

        for (CandidateRoute& candidate : candidates) {
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : candidate.stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(candidate.stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double vote_target
              = static_cast<double>(std::max(2, source_count + 1))
                * (long_net ? 0.70 : 0.76);
          candidate.consensus_penalty = std::max(0.0, vote_target - avg_vote);
        }

        const int max_allowed_vias = anchor_stats.via_count + (long_net ? 2 : 1);
        const long wl_slack = std::max<long>(
            tile_size * (long_net ? 3L : 2L),
            static_cast<long>(std::ceil(
                static_cast<double>(anchor_stats.wirelength_dbu)
                * (long_net ? 0.0018 : 0.0012))));
        const long detour_slack = std::max<long>(
            tile_size * (long_net ? 10L : 7L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_detour))
                * (long_net ? 0.11 : 0.08))));
        const long high_layer_drop_guard = std::max<long>(
            tile_size * (long_net ? 8L : 6L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * (long_net ? 0.08 : 0.10))));
        const long high_layer_rise_guard = std::max<long>(
            tile_size * (long_net ? 20L : 14L),
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, anchor_stats.high_layer_dbu))
                * (long_net ? 0.14 : 0.10))));
        const long min_allowed_high_layer
            = std::max(0L, anchor_stats.high_layer_dbu - high_layer_drop_guard);
        const long max_allowed_high_layer
            = anchor_stats.high_layer_dbu + high_layer_rise_guard;

        const CandidateRoute* best_candidate = &candidates.front();
        EdgeCountMap best_edge_counts = candidates.front().stats.edge_counts;
        double best_objective = std::numeric_limits<double>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        long best_detour = std::numeric_limits<long>::max();
        long best_via_rise = 0;
        double best_crowding_delta = 0.0;
        bool found_candidate = false;

        for (const CandidateRoute& candidate : candidates) {
          if (candidate.stats.via_count > max_allowed_vias
              || candidate.stats.wirelength_dbu
                     > anchor_stats.wirelength_dbu + wl_slack
              || candidate.detour_dbu > anchor_detour + detour_slack
              || candidate.stats.high_layer_dbu > max_allowed_high_layer
              || candidate.stats.high_layer_dbu < min_allowed_high_layer) {
            continue;
          }

          double crowding_delta = 0.0;
          for (const auto& [edge_key, add_usage] : candidate.stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int before_excess = std::max(0, usage - 2);
            const int after_excess = std::max(0, usage + add_usage - 2);
            crowding_delta += static_cast<double>(
                after_excess * after_excess - before_excess * before_excess);
          }

          const long wl_gain = anchor_stats.wirelength_dbu - candidate.stats.wirelength_dbu;
          const long via_rise = std::max(
              0L,
              static_cast<long>(candidate.stats.via_count)
                  - static_cast<long>(anchor_stats.via_count));
          const long via_drop = std::max(
              0L,
              static_cast<long>(anchor_stats.via_count)
                  - static_cast<long>(candidate.stats.via_count));
          const long high_layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - candidate.stats.high_layer_dbu);

          const double objective
              = static_cast<double>(candidate.stats.wirelength_dbu)
                + (long_net ? 0.12 : 0.16)
                      * static_cast<double>(candidate.detour_dbu)
                + static_cast<double>(tile_size) * 0.18
                      * candidate.consensus_penalty
                + static_cast<double>(tile_size) * 0.09 * crowding_delta
                + 3.40 * static_cast<double>(via_rise)
                - 0.45 * static_cast<double>(via_drop)
                + 0.006 * static_cast<double>(high_layer_drop)
                + 0.45 * static_cast<double>(candidate.stats.segment_count)
                + 2.8 * static_cast<double>(candidate.stats.bend_count)
                - (long_net ? 0.58 : 0.46) * static_cast<double>(wl_gain);

          if (objective < best_objective
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu < best_wl)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu == best_wl
                  && candidate.stats.via_count < best_vias)
              || (objective == best_objective
                  && candidate.stats.wirelength_dbu == best_wl
                  && candidate.stats.via_count == best_vias
                  && candidate.detour_dbu < best_detour)) {
            found_candidate = true;
            best_objective = objective;
            best_wl = candidate.stats.wirelength_dbu;
            best_vias = candidate.stats.via_count;
            best_detour = candidate.detour_dbu;
            best_via_rise = via_rise;
            best_crowding_delta = crowding_delta;
            best_candidate = &candidate;
            best_edge_counts = candidate.stats.edge_counts;
          }
        }

        if (!found_candidate) {
          best_candidate = &candidates.front();
          best_edge_counts = candidates.front().stats.edge_counts;
          best_via_rise = 0;
          fallback_nets++;
        }

        const bool non_anchor_pick = best_candidate != &candidates.front();
        if (non_anchor_pick && best_via_rise > 0
            && total_via_rise + best_via_rise > max_total_via_rise) {
          best_candidate = &candidates.front();
          best_edge_counts = candidates.front().stats.edge_counts;
          best_via_rise = 0;
          best_crowding_delta = 0.0;
          budget_clamped_nets++;
        }

        if (best_candidate != &candidates.front()) {
          total_via_rise += best_via_rise;
          if (best_candidate->stats.wirelength_dbu < anchor_stats.wirelength_dbu
              || (best_candidate->stats.wirelength_dbu
                      == anchor_stats.wirelength_dbu
                  && best_candidate->stats.via_count < anchor_stats.via_count)) {
            improved_nets++;
          }
        }

        hybrid_result.routes.emplace(db_net, *best_candidate->route);
        total_crowding_penalty += best_crowding_delta;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t expected_net_count = anchor_result->routes.size();
      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} anchored on {}: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, improved nets {}, anchor-fallbacks {}, "
            "via-rise budget {}/{}, budget-clamped nets {}, crowding-penalty {:.0f}",
            hybrid_name,
            anchor_name,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            improved_nets,
            fallback_nets,
            total_via_rise,
            max_total_via_rise,
            budget_clamped_nets,
            total_crowding_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            7353,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_softcap_hybrid = [&](const std::string& hybrid_name,
                                     int source_count,
                                     long via_weight,
                                     int logger_code) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      FastRouteCore* core = grouter_->fastroute();
      if (core == nullptr || hybrid_nets.empty()) {
        return;
      }

      // SPRoute-style soft-capacity-aware greedy net assembly:
      // prioritize shortest routes but penalize candidates that overfill
      // soft capacities in high-RUDY areas.
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const double overflow_weight = static_cast<double>(tile_size) * 45.0;
      const double near_weight = static_cast<double>(tile_size) * 0.85;
      std::unordered_map<uint64_t, int> edge_usage;
      std::unordered_map<uint64_t, int> softcap_cache;
      edge_usage.reserve(1 << 20);
      softcap_cache.reserve(1 << 20);
      double total_soft_penalty = 0.0;

      std::vector<odb::dbNet*> ordered_nets = hybrid_nets;
      std::unordered_map<odb::dbNet*, long> best_wl_per_net;
      best_wl_per_net.reserve(ordered_nets.size());
      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }
        long best_wl = std::numeric_limits<long>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          long wl = 0;
          for (const GSegment& segment : route_it->second) {
            if (!segment.isVia()) {
              wl += std::abs(segment.final_x - segment.init_x)
                    + std::abs(segment.final_y - segment.init_y);
            }
          }
          best_wl = std::min(best_wl, wl);
        }
        best_wl_per_net.emplace(db_net, best_wl);
      }
      std::sort(ordered_nets.begin(),
                ordered_nets.end(),
                [&](odb::dbNet* lhs, odb::dbNet* rhs) {
                  const auto lhs_it = best_wl_per_net.find(lhs);
                  const auto rhs_it = best_wl_per_net.find(rhs);
                  const long lhs_wl = lhs_it != best_wl_per_net.end()
                                          ? lhs_it->second
                                          : std::numeric_limits<long>::max();
                  const long rhs_wl = rhs_it != best_wl_per_net.end()
                                          ? rhs_it->second
                                          : std::numeric_limits<long>::max();
                  return lhs_wl > rhs_wl;
                });

      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }

        const GRoute* best_route = nullptr;
        EdgeCountMap best_edge_counts;
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        double best_delta_penalty = std::numeric_limits<double>::max();
        double best_objective = std::numeric_limits<double>::max();

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }

          const RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          double delta_penalty = 0.0;
          for (const auto& [edge_key, add_usage] : stats.edge_counts) {
            const auto usage_it = edge_usage.find(edge_key);
            const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
            const int soft_cap = getSoftCapacityForEdge(
                edge_key, core, normalized_rudy, softcap_cache);
            const int near_cap
                = std::max(1, static_cast<int>(std::floor(soft_cap * 0.90)));

            const int over_before = std::max(0, usage - soft_cap);
            const int over_after = std::max(0, usage + add_usage - soft_cap);
            const int near_before = std::max(0, usage - near_cap);
            const int near_after = std::max(0, usage + add_usage - near_cap);

            delta_penalty += overflow_weight
                             * static_cast<double>(over_after * over_after
                                                   - over_before * over_before);
            delta_penalty += near_weight
                             * static_cast<double>(near_after - near_before);
          }

          const double objective
              = static_cast<double>(stats.wirelength_dbu)
                + static_cast<double>(via_weight) * stats.via_count
                + delta_penalty;
          if (objective < best_objective
              || (objective == best_objective
                  && stats.wirelength_dbu < best_wl)
              || (objective == best_objective
                  && stats.wirelength_dbu == best_wl
                  && stats.via_count < best_vias)) {
            best_objective = objective;
            best_delta_penalty = delta_penalty;
            best_wl = stats.wirelength_dbu;
            best_vias = stats.via_count;
            best_route = &route_it->second;
            best_edge_counts = stats.edge_counts;
          }
        }

        if (best_route == nullptr) {
          continue;
        }
        hybrid_result.routes.emplace(db_net, *best_route);
        total_soft_penalty += best_delta_penalty;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(GNR,
                      logger_code,
                      "NEWGR {} from top {} scenarios: wirelength {:.0f} um, "
                      "vias {}, routed nets {}/{}, soft-penalty {:.0f}",
                      hybrid_name,
                      source_count,
                      hybrid_result.metrics.wirelength_um,
                      hybrid_result.metrics.via_count,
                      hybrid_result.routes.size(),
                      expected_net_count,
                      total_soft_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            6014,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

      auto append_consensus_hybrid = [&](const std::string& hybrid_name,
                                       int source_count,
                                       long via_weight,
                                       double detour_weight,
                                       double consensus_weight,
                                       int logger_code) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      double total_consensus_penalty = 0.0;
      double total_detour_penalty = 0.0;

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<std::pair<const GRoute*, RouteEdgeStats>> candidates;
        candidates.reserve(source_count);
        std::unordered_map<uint64_t, int> edge_frequency;
        edge_frequency.reserve(64);

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          for (const auto& [edge_key, usage] : stats.edge_counts) {
            if (usage > 0) {
              edge_frequency[edge_key] += 1;
            }
          }
          candidates.emplace_back(&route_it->second, std::move(stats));
        }

        if (candidates.empty()) {
          continue;
        }

        const GRoute* best_route = nullptr;
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        double best_objective = std::numeric_limits<double>::max();

        for (const auto& [route, stats] : candidates) {
          const long bbox_hpwl = getRouteBBoxHpwl(*route);
          const long detour = std::max(0L, stats.wirelength_dbu - bbox_hpwl);

          long vote_sum = 0;
          for (const auto& [edge_key, usage] : stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count = std::max(1, static_cast<int>(stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          const double consensus_penalty
              = std::max(0.0, static_cast<double>(source_count) - avg_vote);
          const double objective
              = static_cast<double>(stats.wirelength_dbu)
                + static_cast<double>(via_weight) * static_cast<double>(stats.via_count)
                + detour_weight * static_cast<double>(detour)
                + consensus_weight * static_cast<double>(tile_size) * consensus_penalty;

          if (objective < best_objective
              || (objective == best_objective
                  && stats.wirelength_dbu < best_wl)
              || (objective == best_objective
                  && stats.wirelength_dbu == best_wl
                  && stats.via_count < best_vias)) {
            best_objective = objective;
            best_wl = stats.wirelength_dbu;
            best_vias = stats.via_count;
            best_route = route;
          }
        }

        if (best_route != nullptr) {
          const RouteEdgeStats best_stats
              = collectRouteEdgeStats(*best_route, grouter_->grid_);
          const long bbox_hpwl = getRouteBBoxHpwl(*best_route);
          const long detour = std::max(0L, best_stats.wirelength_dbu - bbox_hpwl);
          long vote_sum = 0;
          for (const auto& [edge_key, usage] : best_stats.edge_counts) {
            if (usage > 0) {
              const auto vote_it = edge_frequency.find(edge_key);
              if (vote_it != edge_frequency.end()) {
                vote_sum += vote_it->second;
              }
            }
          }
          const int edge_count
              = std::max(1, static_cast<int>(best_stats.edge_counts.size()));
          const double avg_vote
              = static_cast<double>(vote_sum) / static_cast<double>(edge_count);
          total_consensus_penalty += std::max(
              0.0, static_cast<double>(source_count) - avg_vote);
          total_detour_penalty += static_cast<double>(detour);
          hybrid_result.routes.emplace(db_net, *best_route);
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, consensus-penalty {:.0f}, detour-penalty {:.0f}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            total_consensus_penalty,
            total_detour_penalty);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            6016,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    auto append_pareto_softcap_hybrid = [&](const std::string& hybrid_name,
                                            int source_count,
                                            long via_weight,
                                            int logger_code) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      FastRouteCore* core = grouter_->fastroute();
      if (core == nullptr || hybrid_nets.empty()) {
        return;
      }

      struct CandidateRoute
      {
        const GRoute* route = nullptr;
        RouteEdgeStats stats;
      };

      // Mix of FastRoute shortest-path preference and SPRoute soft-cap
      // pressure: use a per-net Pareto gate first, then pick the candidate
      // with lowest soft-capacity incremental penalty.
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const double overflow_weight = static_cast<double>(tile_size) * 40.0;
      const double near_weight = static_cast<double>(tile_size) * 0.60;
      const double high_layer_weight = 0.035;
      const double layer_span_weight = static_cast<double>(tile_size) * 0.55;
      std::unordered_map<uint64_t, int> edge_usage;
      std::unordered_map<uint64_t, int> softcap_cache;
      edge_usage.reserve(1 << 20);
      softcap_cache.reserve(1 << 20);
      double total_soft_penalty = 0.0;
      int pareto_hits = 0;

      std::vector<odb::dbNet*> ordered_nets = hybrid_nets;
      std::unordered_map<odb::dbNet*, long> best_wl_per_net;
      best_wl_per_net.reserve(ordered_nets.size());
      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }
        long best_wl = std::numeric_limits<long>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          const RouteEdgeStats stats
              = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          best_wl = std::min(best_wl, stats.wirelength_dbu);
        }
        best_wl_per_net.emplace(db_net, best_wl);
      }
      std::sort(ordered_nets.begin(),
                ordered_nets.end(),
                [&](odb::dbNet* lhs, odb::dbNet* rhs) {
                  const auto lhs_it = best_wl_per_net.find(lhs);
                  const auto rhs_it = best_wl_per_net.find(rhs);
                  const long lhs_wl = lhs_it != best_wl_per_net.end()
                                          ? lhs_it->second
                                          : std::numeric_limits<long>::max();
                  const long rhs_wl = rhs_it != best_wl_per_net.end()
                                          ? rhs_it->second
                                          : std::numeric_limits<long>::max();
                  return lhs_wl > rhs_wl;
                });

      for (odb::dbNet* db_net : ordered_nets) {
        if (db_net == nullptr) {
          continue;
        }

        std::vector<CandidateRoute> candidates;
        candidates.reserve(source_count);
        long min_wl = std::numeric_limits<long>::max();
        long min_high_layer = std::numeric_limits<long>::max();
        int min_layer_span = std::numeric_limits<int>::max();
        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          CandidateRoute candidate;
          candidate.route = &route_it->second;
          candidate.stats = collectRouteEdgeStats(route_it->second, grouter_->grid_);
          min_wl = std::min(min_wl, candidate.stats.wirelength_dbu);
          min_high_layer = std::min(min_high_layer, candidate.stats.high_layer_dbu);
          min_layer_span = std::min(min_layer_span, candidate.stats.layer_span);
          candidates.push_back(std::move(candidate));
        }
        if (candidates.empty()) {
          continue;
        }

        const long wl_slack = std::max<long>(
            tile_size * 2L,
            static_cast<long>(std::ceil(static_cast<double>(min_wl) * 0.012)));
        const long high_layer_slack = std::max<long>(
            tile_size * 20L,
            static_cast<long>(
                std::ceil(static_cast<double>(std::max(1L, min_high_layer)) * 0.15)));
        const int layer_span_slack = 1;

        const GRoute* best_route = nullptr;
        EdgeCountMap best_edge_counts;
        double best_delta_penalty = 0.0;
        long best_wl = std::numeric_limits<long>::max();
        int best_vias = std::numeric_limits<int>::max();
        double best_objective = std::numeric_limits<double>::max();

        auto evaluate_candidates = [&](bool use_pareto_gate) {
          bool found = false;
          for (const CandidateRoute& candidate : candidates) {
            const RouteEdgeStats& stats = candidate.stats;
            if (use_pareto_gate) {
              if (stats.wirelength_dbu > min_wl + wl_slack
                  || stats.high_layer_dbu > min_high_layer + high_layer_slack
                  || stats.layer_span > min_layer_span + layer_span_slack) {
                continue;
              }
            }

            double delta_penalty = 0.0;
            for (const auto& [edge_key, add_usage] : stats.edge_counts) {
              const auto usage_it = edge_usage.find(edge_key);
              const int usage = usage_it == edge_usage.end() ? 0 : usage_it->second;
              const int soft_cap = getSoftCapacityForEdge(
                  edge_key, core, normalized_rudy, softcap_cache);
              const int near_cap
                  = std::max(1, static_cast<int>(std::floor(soft_cap * 0.90)));

              const int over_before = std::max(0, usage - soft_cap);
              const int over_after = std::max(0, usage + add_usage - soft_cap);
              const int near_before = std::max(0, usage - near_cap);
              const int near_after = std::max(0, usage + add_usage - near_cap);

              delta_penalty += overflow_weight
                               * static_cast<double>(over_after * over_after
                                                     - over_before * over_before);
              delta_penalty += near_weight
                               * static_cast<double>(near_after - near_before);
            }

            const double objective
                = static_cast<double>(stats.wirelength_dbu)
                  + static_cast<double>(via_weight) * static_cast<double>(stats.via_count)
                  + high_layer_weight * static_cast<double>(stats.high_layer_dbu)
                  + layer_span_weight * static_cast<double>(stats.layer_span)
                  + delta_penalty;
            if (objective < best_objective
                || (objective == best_objective
                    && stats.wirelength_dbu < best_wl)
                || (objective == best_objective
                    && stats.wirelength_dbu == best_wl
                    && stats.via_count < best_vias)) {
              found = true;
              best_objective = objective;
              best_delta_penalty = delta_penalty;
              best_wl = stats.wirelength_dbu;
              best_vias = stats.via_count;
              best_route = candidate.route;
              best_edge_counts = stats.edge_counts;
            }
          }
          return found;
        };

        const bool found_in_pareto = evaluate_candidates(true);
        if (!found_in_pareto) {
          static_cast<void>(evaluate_candidates(false));
        } else {
          pareto_hits++;
        }

        if (best_route == nullptr) {
          continue;
        }

        hybrid_result.routes.emplace(db_net, *best_route);
        total_soft_penalty += best_delta_penalty;
        for (const auto& [edge_key, add_usage] : best_edge_counts) {
          edge_usage[edge_key] += add_usage;
        }
      }

      const std::size_t coverage_threshold
          = expected_net_count > 0 ? (expected_net_count * 95) / 100 : 0;
      if (hybrid_result.routes.size() >= coverage_threshold) {
        hybrid_result.metrics = compute_metrics(hybrid_result.routes);
        logger_->info(
            GNR,
            logger_code,
            "NEWGR {} from top {} scenarios: wirelength {:.0f} um, vias {}, "
            "routed nets {}/{}, soft-penalty {:.0f}, pareto-hits {}",
            hybrid_name,
            source_count,
            hybrid_result.metrics.wirelength_um,
            hybrid_result.metrics.via_count,
            hybrid_result.routes.size(),
            expected_net_count,
            total_soft_penalty,
            pareto_hits);
        scenario_results.push_back(std::move(hybrid_result));
      } else {
        logger_->info(
            GNR,
            6024,
            "NEWGR skipped {} due low net coverage ({}/{}).",
            hybrid_name,
            hybrid_result.routes.size(),
            expected_net_count);
      }
    };

    // Drastic recombination:
    // 1) a wirelength-first hybrid (FastRoute shortest-path intent),
    // 2) a balanced hybrid that softly penalizes vias (SPRoute/CUGR flavor).
    const int wl_source_count = std::min<int>(14, ranked.size());
    const int ultra_wl_source_count = std::min<int>(26, ranked.size());
    const int min_wl_source_count = std::min<int>(36, ranked.size());
    const int min_wl_wide_source_count = ranked.size();
    const int min_wl_extreme_source_count = ranked.size();
    const int absolute_wl_source_count = ranked.size();
    const int hpwl_lock_source_count = std::min<int>(32, ranked.size());
    const int dr_shield_source_count = std::min<int>(30, ranked.size());
    const int dr_stable_source_count = std::min<int>(34, ranked.size());
    const int wl_safe_source_count = std::min<int>(18, ranked.size());
    const int via_floor_wl_source_count = std::min<int>(20, ranked.size());
    const int via_band_wl_source_count = std::min<int>(20, ranked.size());
    const int layer_lift_source_count = std::min<int>(30, ranked.size());
    const int anchor_wl_deep_source_count = std::min<int>(28, ranked.size());
    const int layer_compact_source_count = std::min<int>(20, ranked.size());
    const int balanced_source_count = std::min<int>(8, ranked.size());
    const int softcap_source_count = std::min<int>(16, ranked.size());
    const int consensus_source_count = std::min<int>(22, ranked.size());
    const int detour_source_count = std::min<int>(24, ranked.size());
    const int pareto_softcap_source_count = std::min<int>(28, ranked.size());
    const int length_adaptive_source_count = std::min<int>(24, ranked.size());
    const int wl_feedback_source_count = std::min<int>(28, ranked.size());
    const int wl_corridor_source_count = std::min<int>(30, ranked.size());
    const long balanced_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 3L);
    const long softcap_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 4L);
    const long consensus_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 6L);
    const long detour_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 8L);
    const long wl_safe_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 10L);
    const long layer_compact_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 16L);
    const long dr_shield_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 12L);
    const long dr_stable_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 14L);
    const long pareto_softcap_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 14L);
    const long length_adaptive_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 20L);
    const long wl_feedback_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 18L);
    const long wl_corridor_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 22L);
    const long min_wl_via_weight = 0L;
    if (compact_exploration_mode) {
      const int compact_source_count = std::min<int>(12, ranked.size());
      append_hybrid("hybrid-netmix-wl", compact_source_count, 0, 0.20, 0.45, 6010);
      append_hybrid("hybrid-netmix-wl-compact",
                    compact_source_count,
                    0,
                    0.24,
                    0.62,
                    7324,
                    0.018,
                    0.40,
                    1.20);
      append_wl_corridor_hybrid("hybrid-netmix-wl-corridor",
                                compact_source_count,
                                wl_corridor_via_weight,
                                7332);
      append_length_adaptive_hybrid("hybrid-netmix-length-adaptive",
                                    compact_source_count,
                                    length_adaptive_via_weight,
                                    7310);
      append_feedback_wl_hybrid("hybrid-netmix-wl-feedback",
                                compact_source_count,
                                wl_feedback_via_weight,
                                7313);
      append_via_floor_wl_hybrid("hybrid-netmix-via-floor-wl",
                                 compact_source_count,
                                 1,
                                 7316);
      append_via_band_wl_hybrid("hybrid-netmix-via-band-wl",
                                compact_source_count,
                                0,
                                2,
                                7319);
      append_layer_lifted_wl_hybrid("hybrid-netmix-layer-lift-wl",
                                    compact_source_count,
                                    4,
                                    7326);
      append_anchor_locked_wl_hybrid("hybrid-netmix-anchor-wl",
                                     "hybrid-netmix-wl",
                                     compact_source_count,
                                     2,
                                     7321);
      append_anchor_locked_wl_hybrid("hybrid-netmix-anchor-wl-deep",
                                     "hybrid-netmix-wl",
                                     std::min<int>(compact_source_count + 4,
                                                   static_cast<int>(ranked.size())),
                                     4,
                                     7322);
      append_layer_floor_wl_hybrid("hybrid-netmix-anchor-wl-layerfloor",
                                   "hybrid-netmix-wl",
                                   compact_source_count,
                                   4,
                                   2,
                                   7341);
      append_drt_elastic_wl_hybrid("hybrid-netmix-drt-elastic-wl",
                                   "hybrid-netmix-wl",
                                   compact_source_count,
                                   7351);
      append_budgeted_lift_wl_hybrid("hybrid-netmix-budgeted-lift-wl",
                                     "hybrid-netmix-wl",
                                     compact_source_count,
                                     7352);
      append_min_wl_hybrid("hybrid-netmix-min-wl-wide",
                           compact_source_count,
                           min_wl_via_weight,
                           0.04,
                           0.001,
                           0.02,
                           7103);
      append_min_wl_hybrid("hybrid-netmix-absolute-wl",
                           compact_source_count,
                           0L,
                           0.0,
                           0.0,
                           0.0,
                           7105);
      // Deep HPWL squeeze: pull candidates from all generated scenarios
      // to aggressively minimize net-level wirelength in compact mode.
      append_min_wl_hybrid("hybrid-netmix-absolute-wl-deep",
                           static_cast<int>(ranked.size()),
                           0L,
                           0.02,
                           0.0004,
                           0.015,
                           7360,
                           0.12,
                           1.4);
      // CUGR/SPRoute-inspired blend: keep aggressive WL pressure but penalize
      // route fragmentation and layer collapse to stay DR-friendly.
      append_min_wl_hybrid("hybrid-netmix-cugr-sp-balance-wl",
                           compact_source_count,
                           0L,
                           0.10,
                           0.012,
                           0.24,
                           7107,
                           0.45,
                           3.20);
      append_pareto_softcap_hybrid("hybrid-netmix-pareto-softcap",
                                   compact_source_count,
                                   pareto_softcap_via_weight,
                                   6023);
    } else {
      append_hybrid("hybrid-netmix-wl", wl_source_count, 0, 0.24, 0.50, 6010);
      append_hybrid("hybrid-netmix-wl-compact",
                    wl_source_count,
                    0,
                    0.26,
                    0.64,
                    7324,
                    0.022,
                    0.45,
                    1.30);
      append_wl_corridor_hybrid("hybrid-netmix-wl-corridor",
                                wl_corridor_source_count,
                                wl_corridor_via_weight,
                                7332);
      append_hybrid(
          "hybrid-netmix-ultra-wl", ultra_wl_source_count, 0, 0.40, 0.95, 6013);
      append_length_adaptive_hybrid("hybrid-netmix-length-adaptive",
                                    length_adaptive_source_count,
                                    length_adaptive_via_weight,
                                    7310);
      append_feedback_wl_hybrid("hybrid-netmix-wl-feedback",
                                wl_feedback_source_count,
                                wl_feedback_via_weight,
                                7313);
      append_via_floor_wl_hybrid("hybrid-netmix-via-floor-wl",
                                 via_floor_wl_source_count,
                                 1,
                                 7316);
      append_via_band_wl_hybrid("hybrid-netmix-via-band-wl",
                                via_band_wl_source_count,
                                0,
                                3,
                                7319);
      append_layer_lifted_wl_hybrid("hybrid-netmix-layer-lift-wl",
                                    layer_lift_source_count,
                                    6,
                                    7326);
      append_anchor_locked_wl_hybrid("hybrid-netmix-anchor-wl",
                                     "hybrid-netmix-wl",
                                     via_band_wl_source_count,
                                     3,
                                     7321);
      append_anchor_locked_wl_hybrid("hybrid-netmix-anchor-wl-deep",
                                     "hybrid-netmix-wl",
                                     anchor_wl_deep_source_count,
                                     6,
                                     7322);
      append_layer_floor_wl_hybrid("hybrid-netmix-anchor-wl-layerfloor",
                                   "hybrid-netmix-wl",
                                   anchor_wl_deep_source_count,
                                   5,
                                   3,
                                   7341);
      append_drt_elastic_wl_hybrid("hybrid-netmix-drt-elastic-wl",
                                   "hybrid-netmix-wl",
                                   std::min<int>(anchor_wl_deep_source_count + 2,
                                                 static_cast<int>(ranked.size())),
                                   7351);
      append_budgeted_lift_wl_hybrid("hybrid-netmix-budgeted-lift-wl",
                                     "hybrid-netmix-wl",
                                     std::min<int>(anchor_wl_deep_source_count + 2,
                                                   static_cast<int>(ranked.size())),
                                     7352);
      append_min_wl_hybrid("hybrid-netmix-min-wl",
                           min_wl_source_count,
                           min_wl_via_weight,
                           0.06,
                           0.002,
                           0.03,
                           7102);
      append_min_wl_hybrid("hybrid-netmix-min-wl-wide",
                           min_wl_wide_source_count,
                           min_wl_via_weight,
                           0.08,
                           0.003,
                           0.04,
                           7103);
      append_min_wl_hybrid("hybrid-netmix-min-wl-extreme",
                           min_wl_extreme_source_count,
                           0L,
                           0.04,
                           0.001,
                           0.02,
                           7104);
      append_min_wl_hybrid("hybrid-netmix-smooth-wl",
                           min_wl_source_count,
                           0L,
                           0.18,
                           0.018,
                           0.26,
                           7106,
                           0.55,
                           4.8);
      append_min_wl_hybrid("hybrid-netmix-absolute-wl",
                           absolute_wl_source_count,
                           0L,
                           0.0,
                           0.0,
                           0.0,
                           7105);
      append_min_wl_hybrid("hybrid-netmix-absolute-wl-deep",
                           absolute_wl_source_count,
                           0L,
                           0.02,
                           0.0004,
                           0.015,
                           7360,
                           0.12,
                           1.4);
      append_hybrid("hybrid-netmix-hpwl-lock",
                    hpwl_lock_source_count,
                    0,
                    0.18,
                    0.42,
                    6021,
                    0.010,
                    0.12,
                    0.45);
      append_hybrid("hybrid-netmix-dr-shield",
                    dr_shield_source_count,
                    dr_shield_via_weight,
                    0.62,
                    0.76,
                    6022,
                    0.060,
                    0.95,
                    1.30);
      append_dr_stable_hybrid(
          "hybrid-netmix-dr-stable", dr_stable_source_count, dr_stable_via_weight, 7309);
      append_hybrid("hybrid-netmix-wl-safe",
                    wl_safe_source_count,
                    wl_safe_via_weight,
                    0.26,
                    0.52,
                    6019,
                    0.018,
                    0.30,
                    0.95);
      append_hybrid("hybrid-netmix-detour-ladder",
                    detour_source_count,
                    detour_via_weight,
                    0.95,
                    0.52,
                    6018,
                    0.030,
                    0.45,
                    1.15);
      append_hybrid("hybrid-netmix-layer-compact",
                    layer_compact_source_count,
                    layer_compact_via_weight,
                    0.78,
                    0.72,
                    6020,
                    0.080,
                    1.10,
                    1.45);
      append_hybrid("hybrid-netmix-balanced",
                    balanced_source_count,
                    balanced_via_weight,
                    0.32,
                    0.65,
                    6012);
      append_softcap_hybrid(
          "hybrid-netmix-softcap", softcap_source_count, softcap_via_weight, 6015);
      append_consensus_hybrid("hybrid-netmix-consensus",
                              consensus_source_count,
                              consensus_via_weight,
                              0.55,
                              0.85,
                              6017);
      append_pareto_softcap_hybrid("hybrid-netmix-pareto-softcap",
                                   pareto_softcap_source_count,
                                   pareto_softcap_via_weight,
                                   6023);
    }

    auto find_scenario = [&](const std::string& name) -> const ScenarioResult* {
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == name) {
          return &result;
        }
      }
      return nullptr;
    };

    const ScenarioResult* wl_hybrid = find_scenario("hybrid-netmix-wl");
    const ScenarioResult* absolute_wl_hybrid
        = find_scenario("hybrid-netmix-absolute-wl");
    const ScenarioResult* absolute_wl_deep_hybrid
        = find_scenario("hybrid-netmix-absolute-wl-deep");
    const ScenarioResult* layer_lift_wl_hybrid
        = find_scenario("hybrid-netmix-layer-lift-wl");
    const ScenarioResult* dr_stable_hybrid
        = find_scenario("hybrid-netmix-dr-stable");
    const ScenarioResult* layer_compact_hybrid
        = find_scenario("hybrid-netmix-layer-compact");
    const ScenarioResult* via_band_hybrid
        = find_scenario("hybrid-netmix-via-band-wl");
    const ScenarioResult* anchor_wl_hybrid
        = find_scenario("hybrid-netmix-anchor-wl");
    const ScenarioResult* anchor_wl_deep_hybrid
        = find_scenario("hybrid-netmix-anchor-wl-deep");
    const ScenarioResult* wl_corridor_hybrid
        = find_scenario("hybrid-netmix-wl-corridor");
    const ScenarioResult* length_adaptive_hybrid
        = find_scenario("hybrid-netmix-length-adaptive");
    const ScenarioResult* min_wl_wide_hybrid
        = find_scenario("hybrid-netmix-min-wl-wide");
    const ScenarioResult* cugr_sp_balance_hybrid
        = find_scenario("hybrid-netmix-cugr-sp-balance-wl");
    const ScenarioResult* budgeted_lift_wl_hybrid
        = find_scenario("hybrid-netmix-budgeted-lift-wl");

    // FastRoute/SPRoute hybridization:
    // keep the stable wl-hybrid as anchor and greedily import donor routes
    // from aggressive low-WL scenarios under a global high-layer budget.
    if (wl_hybrid != nullptr && absolute_wl_hybrid != nullptr) {
      ScenarioResult layerbudget_result;
      layerbudget_result.name = "hybrid-netmix-wl-layerbudget";
      layerbudget_result.routes = wl_hybrid->routes;

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long min_wl_gain = std::max<long>(8L, tile_size / 2L);
      const long detour_guard = std::max<long>(tile_size * 18L, 9000L);
      const long detour_tie_guard = std::max<long>(tile_size * 2L, 900L);
      const int via_rise_cap = 6;
      const long total_via_rise_cap = std::max<long>(
          24L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_hybrid->metrics.via_count) * 0.00024)));
      const long high_layer_drop_budget = std::max<long>(
          tile_size * 18000L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_hybrid->metrics.high_layer_dbu) * 0.20)));
      const long high_layer_drop_guard = std::max<long>(tile_size * 4L, 1400L);
      const long recover_wl_guard = std::max<long>(tile_size * 5L, 2400L);

      long consumed_high_layer_drop = 0L;
      long total_via_rise = 0L;
      int imported_absolute_nets = 0;
      int imported_absdeep_nets = 0;
      int imported_min_wl_nets = 0;
      int imported_balance_nets = 0;
      int via_budget_clamped_nets = 0;
      int recovered_layer_lift_nets = 0;
      long wl_gain_sum = 0L;
      long via_delta_sum = 0L;

      for (auto& [db_net, selected_route] : layerbudget_result.routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(selected_route, grouter_->grid_);
        const long route_hpwl = getRouteBBoxHpwl(selected_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - route_hpwl);

        const GRoute* best_route = nullptr;
        RouteEdgeStats best_stats = anchor_stats;
        long best_wl_gain = 0L;
        long best_via_delta = 0L;
        long best_high_layer_drop = 0L;
        long best_detour = anchor_detour;
        int best_donor = 0;

        auto consider_donor = [&](const ScenarioResult* donor, int donor_id) {
          if (donor == nullptr) {
            return;
          }
          const auto donor_it = donor->routes.find(db_net);
          if (donor_it == donor->routes.end()) {
            return;
          }
          const RouteEdgeStats donor_stats
              = collectRouteEdgeStats(donor_it->second, grouter_->grid_);
          const long wl_gain = anchor_stats.wirelength_dbu - donor_stats.wirelength_dbu;
          if (wl_gain < min_wl_gain) {
            return;
          }
          const long via_delta = static_cast<long>(donor_stats.via_count)
                                 - static_cast<long>(anchor_stats.via_count);
          if (via_delta > via_rise_cap) {
            return;
          }
          const long donor_detour
              = std::max(0L, donor_stats.wirelength_dbu - route_hpwl);
          if (donor_detour > anchor_detour + detour_guard) {
            return;
          }

          const long high_layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - donor_stats.high_layer_dbu);
          if (high_layer_drop > 0 && wl_gain < min_wl_gain * 2L) {
            return;
          }
          if (consumed_high_layer_drop + high_layer_drop > high_layer_drop_budget
              && wl_gain < min_wl_gain * 4L) {
            return;
          }

          const bool better_wl
              = donor_stats.wirelength_dbu < best_stats.wirelength_dbu;
          const bool tie_break_detour
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_detour + detour_tie_guard < best_detour;
          const bool tie_break_high_layer
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && std::llabs(donor_detour - best_detour) <= detour_tie_guard
                && donor_stats.high_layer_dbu
                       > best_stats.high_layer_dbu + detour_tie_guard;
          const bool tie_break_via
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.high_layer_dbu + detour_tie_guard
                       >= best_stats.high_layer_dbu
                && donor_stats.via_count < best_stats.via_count;
          if (!(better_wl || tie_break_detour || tie_break_high_layer
                || tie_break_via)) {
            return;
          }

          best_route = &donor_it->second;
          best_stats = donor_stats;
          best_wl_gain = wl_gain;
          best_via_delta = via_delta;
          best_high_layer_drop = high_layer_drop;
          best_detour = donor_detour;
          best_donor = donor_id;
        };

        consider_donor(absolute_wl_hybrid, 1);
        consider_donor(absolute_wl_deep_hybrid, 2);
        consider_donor(min_wl_wide_hybrid, 3);
        consider_donor(cugr_sp_balance_hybrid, 4);

        if (best_route != nullptr) {
          if (best_via_delta > 0
              && total_via_rise + best_via_delta > total_via_rise_cap) {
            via_budget_clamped_nets++;
            continue;
          }
          selected_route = *best_route;
          consumed_high_layer_drop += best_high_layer_drop;
          if (best_via_delta > 0) {
            total_via_rise += best_via_delta;
          }
          wl_gain_sum += best_wl_gain;
          via_delta_sum += best_via_delta;
          if (best_donor == 1) {
            imported_absolute_nets++;
          } else if (best_donor == 2) {
            imported_absdeep_nets++;
          } else if (best_donor == 3) {
            imported_min_wl_nets++;
          } else if (best_donor == 4) {
            imported_balance_nets++;
          }
        }
      }

      if (layer_lift_wl_hybrid != nullptr
          && consumed_high_layer_drop > high_layer_drop_budget * 3L / 4L) {
        for (auto& [db_net, selected_route] : layerbudget_result.routes) {
          if (consumed_high_layer_drop
              <= std::max<long>(0L, high_layer_drop_budget - high_layer_drop_guard)) {
            break;
          }
          if (db_net == nullptr) {
            continue;
          }
          const auto lift_it = layer_lift_wl_hybrid->routes.find(db_net);
          if (lift_it == layer_lift_wl_hybrid->routes.end()) {
            continue;
          }

          const RouteEdgeStats current_stats
              = collectRouteEdgeStats(selected_route, grouter_->grid_);
          const RouteEdgeStats lift_stats
              = collectRouteEdgeStats(lift_it->second, grouter_->grid_);
          if (lift_stats.high_layer_dbu <= current_stats.high_layer_dbu) {
            continue;
          }
          if (lift_stats.wirelength_dbu
              > current_stats.wirelength_dbu + recover_wl_guard) {
            continue;
          }

          const long layer_recovery
              = lift_stats.high_layer_dbu - current_stats.high_layer_dbu;
          selected_route = lift_it->second;
          consumed_high_layer_drop = std::max<long>(
              0L, consumed_high_layer_drop - layer_recovery);
          recovered_layer_lift_nets++;
        }
      }

      layerbudget_result.metrics = compute_metrics(layerbudget_result.routes);
      apply_routability_proxy(layerbudget_result.metrics);
      logger_->info(
          GNR,
          7345,
          "NEWGR hybrid-netmix-wl-layerbudget imported "
          "abs/absdeep/min/balance nets {}/{}/{}/{}, "
          "layer recoveries {}, wirelength {:.0f} um, vias {}, via-rise "
          "budget {}/{}, budget-clamped nets {}, high-layer drop {}/{}, wl "
          "gain sum {}, via delta sum {}",
          imported_absolute_nets,
          imported_absdeep_nets,
          imported_min_wl_nets,
          imported_balance_nets,
          recovered_layer_lift_nets,
          layerbudget_result.metrics.wirelength_um,
          layerbudget_result.metrics.via_count,
          total_via_rise,
          total_via_rise_cap,
          via_budget_clamped_nets,
          consumed_high_layer_drop,
          high_layer_drop_budget,
          wl_gain_sum,
          via_delta_sum);
      scenario_results.push_back(std::move(layerbudget_result));

      // Radical wirelength fusion: mix anchor, absolute-WL and min-WL donors
      // with relaxed per-net guards to force lower-WL guide alternatives.
      ScenarioResult radical_fusion_result;
      radical_fusion_result.name = "hybrid-netmix-radical-wl-fusion";
      radical_fusion_result.routes = wl_hybrid->routes;

      const long radical_min_wl_gain = std::max<long>(6L, tile_size / 3L);
      const long radical_detour_guard = std::max<long>(tile_size * 24L, 12000L);
      const int radical_via_rise_cap = 8;
      const long radical_via_trade_wl_factor = std::max<long>(2200L, tile_size * 5L);
      const long radical_total_via_rise_cap = std::max<long>(
          48L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_hybrid->metrics.via_count) * 0.00045)));
      const long radical_high_layer_drop_budget = std::max<long>(
          tile_size * 26000L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_hybrid->metrics.high_layer_dbu) * 0.28)));
      long radical_high_layer_drop = 0L;
      long radical_total_via_rise = 0L;
      int radical_imported_nets = 0;
      int radical_via_budget_clamps = 0;
      long radical_wl_gain_sum = 0L;
      long radical_via_delta_sum = 0L;
      int radical_abs_imports = 0;
      int radical_min_imports = 0;
      int radical_balance_imports = 0;
      int radical_budgeted_imports = 0;
      int radical_abs_deep_imports = 0;

      for (auto& [db_net, selected_route] : radical_fusion_result.routes) {
        if (db_net == nullptr) {
          continue;
        }
        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(selected_route, grouter_->grid_);
        const long route_hpwl = getRouteBBoxHpwl(selected_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - route_hpwl);

        const GRoute* best_route = nullptr;
        RouteEdgeStats best_stats = anchor_stats;
        long best_wl_gain = 0L;
        long best_via_delta = 0L;
        long best_layer_drop = 0L;
        int best_donor = 0;

        auto consider_radical = [&](const ScenarioResult* donor, int donor_id) {
          if (donor == nullptr) {
            return;
          }
          const auto donor_it = donor->routes.find(db_net);
          if (donor_it == donor->routes.end()) {
            return;
          }
          const RouteEdgeStats donor_stats
              = collectRouteEdgeStats(donor_it->second, grouter_->grid_);
          const long wl_gain = anchor_stats.wirelength_dbu - donor_stats.wirelength_dbu;
          if (wl_gain < radical_min_wl_gain) {
            return;
          }
          const long via_delta = static_cast<long>(donor_stats.via_count)
                                 - static_cast<long>(anchor_stats.via_count);
          if (via_delta > radical_via_rise_cap) {
            return;
          }
          if (via_delta > 0 && wl_gain < via_delta * radical_via_trade_wl_factor) {
            return;
          }
          const long donor_detour
              = std::max(0L, donor_stats.wirelength_dbu - route_hpwl);
          if (donor_detour > anchor_detour + radical_detour_guard) {
            return;
          }
          const long layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - donor_stats.high_layer_dbu);
          if (radical_high_layer_drop + layer_drop > radical_high_layer_drop_budget
              && wl_gain < radical_min_wl_gain * 5L) {
            return;
          }

          const bool better_wl
              = donor_stats.wirelength_dbu < best_stats.wirelength_dbu;
          const bool tie_break_via
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count < best_stats.via_count;
          if (!(better_wl || tie_break_via)) {
            return;
          }

          best_route = &donor_it->second;
          best_stats = donor_stats;
          best_wl_gain = wl_gain;
          best_via_delta = via_delta;
          best_layer_drop = layer_drop;
          best_donor = donor_id;
        };

        consider_radical(absolute_wl_hybrid, 1);
        consider_radical(min_wl_wide_hybrid, 2);
        consider_radical(cugr_sp_balance_hybrid, 3);
        consider_radical(budgeted_lift_wl_hybrid, 4);
        consider_radical(absolute_wl_deep_hybrid, 5);

        if (best_route != nullptr) {
          if (best_via_delta > 0
              && radical_total_via_rise + best_via_delta
                     > radical_total_via_rise_cap) {
            radical_via_budget_clamps++;
            continue;
          }
          selected_route = *best_route;
          radical_high_layer_drop += best_layer_drop;
          if (best_via_delta > 0) {
            radical_total_via_rise += best_via_delta;
          }
          radical_imported_nets++;
          radical_wl_gain_sum += best_wl_gain;
          radical_via_delta_sum += best_via_delta;
          if (best_donor == 1) {
            radical_abs_imports++;
          } else if (best_donor == 2) {
            radical_min_imports++;
          } else if (best_donor == 3) {
            radical_balance_imports++;
          } else if (best_donor == 4) {
            radical_budgeted_imports++;
          } else if (best_donor == 5) {
            radical_abs_deep_imports++;
          }
        }
      }

      radical_fusion_result.metrics = compute_metrics(radical_fusion_result.routes);
      apply_routability_proxy(radical_fusion_result.metrics);
      logger_->info(
          GNR,
          7355,
          "NEWGR hybrid-netmix-radical-wl-fusion imported "
          "abs/min/balance/budgeted/abs-deep nets {}/{}/{}/{}/{}, wirelength "
          "{:.0f} um, "
          "vias {}, high-layer drop {}/{}, via-rise budget {}/{}, "
          "budget-clamped nets {}, wl gain sum {}, via delta sum {}",
          radical_abs_imports,
          radical_min_imports,
          radical_balance_imports,
          radical_budgeted_imports,
          radical_abs_deep_imports,
          radical_fusion_result.metrics.wirelength_um,
          radical_fusion_result.metrics.via_count,
          radical_high_layer_drop,
          radical_high_layer_drop_budget,
          radical_total_via_rise,
          radical_total_via_rise_cap,
          radical_via_budget_clamps,
          radical_wl_gain_sum,
          radical_via_delta_sum);
      if (radical_imported_nets > 0) {
        scenario_results.push_back(std::move(radical_fusion_result));
      }

      // Long-net radical fusion:
      // preserve the layerbudget backbone for short/local nets, but import
      // aggressive absolute/min-WL donors on long-HPWL nets under explicit
      // via/detour/high-layer budgets.
      ScenarioResult longhpwl_radical_result;
      longhpwl_radical_result.name = "hybrid-netmix-longhpwl-radical-wl";
      longhpwl_radical_result.routes = layerbudget_result.routes;

      const long long_hpwl_threshold = std::max<long>(tile_size * 22L, 11000L);
      const long long_min_wl_gain = std::max<long>(6L, tile_size / 3L);
      const long short_min_wl_gain = std::max<long>(10L, tile_size / 2L);
      const int long_via_rise_cap = 10;
      const int short_via_rise_cap = 2;
      const long long_detour_guard = std::max<long>(tile_size * 26L, 13000L);
      const long short_detour_guard = std::max<long>(tile_size * 9L, 5000L);
      const long longhpwl_total_via_rise_cap = std::max<long>(
          56L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_hybrid->metrics.via_count) * 0.00055)));
      const long longhpwl_high_layer_drop_budget = std::max<long>(
          tile_size * 9000L,
          static_cast<long>(std::ceil(
              static_cast<double>(layerbudget_result.metrics.high_layer_dbu) * 0.06)));
      const long longhpwl_high_layer_floor = std::max<long>(
          tile_size * 14L,
          static_cast<long>(std::ceil(
              static_cast<double>(layerbudget_result.metrics.high_layer_dbu) * 0.95)));
      const long per_net_high_layer_drop_cap
          = std::max<long>(tile_size * 120L, 6800L);

      long longhpwl_total_via_rise = 0L;
      long longhpwl_high_layer_drop = 0L;
      int longhpwl_imported_nets = 0;
      int longhpwl_imported_long_nets = 0;
      int longhpwl_imported_short_nets = 0;
      int longhpwl_via_clamped_nets = 0;
      long longhpwl_wl_gain_sum = 0L;
      long longhpwl_via_delta_sum = 0L;
      int longhpwl_abs_imports = 0;
      int longhpwl_absdeep_imports = 0;
      int longhpwl_min_imports = 0;
      int longhpwl_balance_imports = 0;

      for (auto& [db_net, selected_route] : longhpwl_radical_result.routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(selected_route, grouter_->grid_);
        const long route_hpwl = getRouteBBoxHpwl(selected_route);
        const bool is_long_net = route_hpwl >= long_hpwl_threshold
                                 || anchor_stats.wirelength_dbu
                                        >= long_hpwl_threshold * 2L;
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - route_hpwl);
        const long min_wl_gain = is_long_net ? long_min_wl_gain : short_min_wl_gain;
        const int via_rise_cap = is_long_net ? long_via_rise_cap : short_via_rise_cap;
        const long detour_guard = is_long_net ? long_detour_guard : short_detour_guard;

        const GRoute* best_route = nullptr;
        RouteEdgeStats best_stats = anchor_stats;
        long best_wl_gain = 0L;
        long best_via_delta = 0L;
        long best_layer_drop = 0L;
        int best_donor = 0;

        auto consider_longhpwl_donor = [&](const ScenarioResult* donor, int donor_id) {
          if (donor == nullptr) {
            return;
          }
          const auto donor_it = donor->routes.find(db_net);
          if (donor_it == donor->routes.end()) {
            return;
          }
          const RouteEdgeStats donor_stats
              = collectRouteEdgeStats(donor_it->second, grouter_->grid_);
          const long wl_gain = anchor_stats.wirelength_dbu - donor_stats.wirelength_dbu;
          if (wl_gain < min_wl_gain) {
            return;
          }

          const long via_delta = static_cast<long>(donor_stats.via_count)
                                 - static_cast<long>(anchor_stats.via_count);
          if (via_delta > via_rise_cap) {
            return;
          }

          const long donor_detour
              = std::max(0L, donor_stats.wirelength_dbu - route_hpwl);
          if (donor_detour > anchor_detour + detour_guard) {
            return;
          }

          const long layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - donor_stats.high_layer_dbu);
          if (layer_drop > per_net_high_layer_drop_cap
              && wl_gain < min_wl_gain * 10L) {
            return;
          }
          if (longhpwl_high_layer_drop + layer_drop > longhpwl_high_layer_drop_budget
              && wl_gain < min_wl_gain * (is_long_net ? 14L : 9L)) {
            return;
          }

          if (!is_long_net && via_delta > 0) {
            const long short_via_trade_wl = std::max<long>(tile_size * 4L, 1800L);
            if (wl_gain < std::max<long>(min_wl_gain * 4L,
                                         via_delta * short_via_trade_wl)) {
              return;
            }
          }

          const bool better_wl
              = donor_stats.wirelength_dbu < best_stats.wirelength_dbu;
          const bool tie_break_via
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count < best_stats.via_count;
          if (!(better_wl || tie_break_via)) {
            return;
          }

          best_route = &donor_it->second;
          best_stats = donor_stats;
          best_wl_gain = wl_gain;
          best_via_delta = via_delta;
          best_layer_drop = layer_drop;
          best_donor = donor_id;
        };

        consider_longhpwl_donor(absolute_wl_hybrid, 1);
        consider_longhpwl_donor(absolute_wl_deep_hybrid, 2);
        consider_longhpwl_donor(min_wl_wide_hybrid, 3);
        consider_longhpwl_donor(cugr_sp_balance_hybrid, 4);

        if (best_route == nullptr) {
          continue;
        }

        if (best_via_delta > 0
            && longhpwl_total_via_rise + best_via_delta
                   > longhpwl_total_via_rise_cap) {
          longhpwl_via_clamped_nets++;
          continue;
        }

        selected_route = *best_route;
        longhpwl_imported_nets++;
        if (is_long_net) {
          longhpwl_imported_long_nets++;
        } else {
          longhpwl_imported_short_nets++;
        }
        longhpwl_high_layer_drop += best_layer_drop;
        if (best_via_delta > 0) {
          longhpwl_total_via_rise += best_via_delta;
        }
        longhpwl_wl_gain_sum += best_wl_gain;
        longhpwl_via_delta_sum += best_via_delta;
        if (best_donor == 1) {
          longhpwl_abs_imports++;
        } else if (best_donor == 2) {
          longhpwl_absdeep_imports++;
        } else if (best_donor == 3) {
          longhpwl_min_imports++;
        } else if (best_donor == 4) {
          longhpwl_balance_imports++;
        }
      }

      if (longhpwl_imported_nets > 0) {
        longhpwl_radical_result.metrics
            = compute_metrics(longhpwl_radical_result.routes);
        apply_routability_proxy(longhpwl_radical_result.metrics);
        if (longhpwl_radical_result.metrics.high_layer_dbu
            >= longhpwl_high_layer_floor) {
          logger_->info(
              GNR,
              7361,
              "NEWGR hybrid-netmix-longhpwl-radical-wl imported "
              "abs/absdeep/min/balance nets {}/{}/{}/{}, long/short imports "
              "{}/{}, wirelength {:.0f} um, vias {}, via-rise budget {}/{}, "
              "high-layer drop {}/{}, wl gain sum {}, via delta sum {}, "
              "budget-clamped nets {}",
              longhpwl_abs_imports,
              longhpwl_absdeep_imports,
              longhpwl_min_imports,
              longhpwl_balance_imports,
              longhpwl_imported_long_nets,
              longhpwl_imported_short_nets,
              longhpwl_radical_result.metrics.wirelength_um,
              longhpwl_radical_result.metrics.via_count,
              longhpwl_total_via_rise,
              longhpwl_total_via_rise_cap,
              longhpwl_high_layer_drop,
              longhpwl_high_layer_drop_budget,
              longhpwl_wl_gain_sum,
              longhpwl_via_delta_sum,
              longhpwl_via_clamped_nets);
          scenario_results.push_back(std::move(longhpwl_radical_result));
        } else {
          logger_->info(
              GNR,
              7362,
              "NEWGR skipped hybrid-netmix-longhpwl-radical-wl due high-layer "
              "floor violation ({}/{}).",
              longhpwl_radical_result.metrics.high_layer_dbu,
              longhpwl_high_layer_floor);
        }
      }

      // SPRoute-style absolute WL base with FastRoute/CUGR structural recovery:
      // start from the absolute minimum-WL hybrid, then selectively recover
      // vias / detour / high-layer structure from alternate mixes under a
      // tight total wirelength giveback budget.
      const ScenarioResult* wl_layerbudget_hybrid
          = find_scenario("hybrid-netmix-wl-layerbudget");
      const ScenarioResult* radical_wl_hybrid
          = find_scenario("hybrid-netmix-radical-wl-fusion");
      const ScenarioResult* longhpwl_radical_hybrid
          = find_scenario("hybrid-netmix-longhpwl-radical-wl");
      ScenarioResult absolute_via_recover_result;
      absolute_via_recover_result.name = "hybrid-netmix-absolute-via-recover";
      absolute_via_recover_result.routes = absolute_wl_hybrid->routes;

      const long recover_wl_backoff_guard = std::max<long>(tile_size * 3L, 1800L);
      const long recover_detour_gain_guard = std::max<long>(tile_size * 2L, 1200L);
      const long recover_high_layer_gain_guard
          = std::max<long>(tile_size * 2L, 900L);
      const long recover_total_wl_backoff_budget = std::max<long>(
          tile_size * 12000L,
          static_cast<long>(std::ceil(
              static_cast<double>(absolute_wl_hybrid->metrics.wirelength_dbu)
              * 0.0010)));
      const long recover_via_gain_target = std::max<long>(
          1200L,
          static_cast<long>(std::ceil(
              std::max(
                  0.0,
                  static_cast<double>(wl_hybrid->metrics.via_count
                                      - absolute_wl_hybrid->metrics.via_count))
              * 0.55)));

      long recover_total_wl_backoff = 0L;
      long recover_total_via_gain = 0L;
      long recover_detour_gain_sum = 0L;
      long recover_high_layer_gain_sum = 0L;
      int recover_swapped_nets = 0;
      int recover_layerbudget_imports = 0;
      int recover_anchor_imports = 0;
      int recover_budgeted_imports = 0;
      int recover_radical_imports = 0;
      int recover_balance_imports = 0;
      int recover_longhpwl_imports = 0;

      for (auto& [db_net, selected_route] : absolute_via_recover_result.routes) {
        if (db_net == nullptr) {
          continue;
        }
        const RouteEdgeStats absolute_stats
            = collectRouteEdgeStats(selected_route, grouter_->grid_);

        const GRoute* best_route = nullptr;
        RouteEdgeStats best_stats = absolute_stats;
        long best_wl_backoff = 0L;
        long best_via_gain = 0L;
        long best_detour_gain = 0L;
        long best_high_layer_gain = 0L;
        double best_recover_score = 0.0;
        int best_donor = 0;

        auto consider_recover_donor = [&](const ScenarioResult* donor, int donor_id) {
          if (donor == nullptr) {
            return;
          }
          const auto donor_it = donor->routes.find(db_net);
          if (donor_it == donor->routes.end()) {
            return;
          }
          const RouteEdgeStats donor_stats
              = collectRouteEdgeStats(donor_it->second, grouter_->grid_);
          const long wl_backoff
              = donor_stats.wirelength_dbu - absolute_stats.wirelength_dbu;
          if (wl_backoff > recover_wl_backoff_guard) {
            return;
          }

          const long via_gain
              = static_cast<long>(donor_stats.via_count)
                - static_cast<long>(absolute_stats.via_count);
          const long detour_gain
              = std::max(0L, absolute_stats.wirelength_dbu - getRouteBBoxHpwl(selected_route))
                - std::max(0L, donor_stats.wirelength_dbu
                                 - getRouteBBoxHpwl(donor_it->second));
          const long high_layer_gain
              = std::max(0L, donor_stats.high_layer_dbu - absolute_stats.high_layer_dbu);

          const bool has_recovery_signal
              = via_gain >= 2 || detour_gain >= recover_detour_gain_guard
                || high_layer_gain >= recover_high_layer_gain_guard;
          if (!has_recovery_signal) {
            return;
          }

          if (via_gain <= 0 && wl_backoff > 0
              && detour_gain < recover_detour_gain_guard * 2L
              && high_layer_gain < recover_high_layer_gain_guard * 2L) {
            return;
          }

          const double recover_score = static_cast<double>(detour_gain)
                                       + static_cast<double>(high_layer_gain) * 0.40
                                       + static_cast<double>(via_gain) * 480.0
                                       - static_cast<double>(std::max(0L, wl_backoff))
                                             * 1.60;
          const bool better_score
              = best_route == nullptr || recover_score > best_recover_score + 1e-3;
          const bool tie_break_wl
              = std::abs(recover_score - best_recover_score) <= 1e-3
                && donor_stats.wirelength_dbu < best_stats.wirelength_dbu;
          const bool tie_break_via
              = std::abs(recover_score - best_recover_score) <= 1e-3
                && donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count > best_stats.via_count;
          if (!(better_score || tie_break_wl || tie_break_via)) {
            return;
          }

          best_route = &donor_it->second;
          best_stats = donor_stats;
          best_wl_backoff = std::max(0L, wl_backoff);
          best_via_gain = via_gain;
          best_detour_gain = detour_gain;
          best_high_layer_gain = high_layer_gain;
          best_recover_score = recover_score;
          best_donor = donor_id;
        };

        consider_recover_donor(wl_layerbudget_hybrid, 1);
        consider_recover_donor(wl_hybrid, 2);
        consider_recover_donor(budgeted_lift_wl_hybrid, 3);
        consider_recover_donor(radical_wl_hybrid, 4);
        consider_recover_donor(cugr_sp_balance_hybrid, 5);
        consider_recover_donor(longhpwl_radical_hybrid, 6);

        if (best_route != nullptr) {
          const bool within_wl_budget
              = recover_total_wl_backoff + best_wl_backoff
                    <= recover_total_wl_backoff_budget;
          const bool must_take_recovery
              = recover_total_via_gain < recover_via_gain_target && best_via_gain >= 4;
          if (!within_wl_budget && !must_take_recovery) {
            continue;
          }
          selected_route = *best_route;
          recover_total_wl_backoff += best_wl_backoff;
          recover_total_via_gain += best_via_gain;
          recover_detour_gain_sum += best_detour_gain;
          recover_high_layer_gain_sum += best_high_layer_gain;
          recover_swapped_nets++;
          if (best_donor == 1) {
            recover_layerbudget_imports++;
          } else if (best_donor == 2) {
            recover_anchor_imports++;
          } else if (best_donor == 3) {
            recover_budgeted_imports++;
          } else if (best_donor == 4) {
            recover_radical_imports++;
          } else if (best_donor == 5) {
            recover_balance_imports++;
          } else if (best_donor == 6) {
            recover_longhpwl_imports++;
          }
        }
      }

      if (recover_swapped_nets > 0) {
        absolute_via_recover_result.metrics
            = compute_metrics(absolute_via_recover_result.routes);
        apply_routability_proxy(absolute_via_recover_result.metrics);
        logger_->info(
            GNR,
            7357,
            "NEWGR hybrid-netmix-absolute-via-recover imported "
            "layerbudget/anchor/budgeted/radical/balance/longhpwl nets "
            "{}/{}/{}/{}/{}/{}, wirelength {:.0f} um, vias {}, wl backoff "
            "{}/{}, via gain target {}/{}, detour gain sum {}, high-layer gain "
            "sum {}",
            recover_layerbudget_imports,
            recover_anchor_imports,
            recover_budgeted_imports,
            recover_radical_imports,
            recover_balance_imports,
            recover_longhpwl_imports,
            absolute_via_recover_result.metrics.wirelength_um,
            absolute_via_recover_result.metrics.via_count,
            recover_total_wl_backoff,
            recover_total_wl_backoff_budget,
            recover_total_via_gain,
            recover_via_gain_target,
            recover_detour_gain_sum,
            recover_high_layer_gain_sum);
        scenario_results.push_back(std::move(absolute_via_recover_result));
      }

      // Hyper absolute WL fusion:
      // start from absolute-WL and greedily import per-net lower-WL routes
      // across the full hybrid donor set with only soft via/detour guards.
      ScenarioResult hyper_absolute_wl_result;
      hyper_absolute_wl_result.name = "hybrid-netmix-hyper-absolute-wl";
      const ScenarioResult* hyper_anchor_hybrid
          = find_scenario("hybrid-netmix-absolute-via-recover");
      if (hyper_anchor_hybrid == nullptr) {
        hyper_anchor_hybrid = absolute_wl_hybrid;
      }
      hyper_absolute_wl_result.routes = hyper_anchor_hybrid->routes;

      const long hyper_detour_guard = std::max<long>(tile_size * 42L, 22000L);
      const long hyper_trade_via_to_wl
          = std::max<long>(120L, static_cast<long>(tile_size * 0.18));
      const long hyper_high_layer_drop_budget = std::max<long>(
          tile_size * 12000L,
          static_cast<long>(std::ceil(
              static_cast<double>(hyper_anchor_hybrid->metrics.high_layer_dbu)
              * 0.10)));
      const long hyper_per_net_layer_drop_cap
          = std::max<long>(tile_size * 180L, 9200L);
      const long hyper_total_via_rise_cap = std::max<long>(
          420L,
          static_cast<long>(std::ceil(
              static_cast<double>(hyper_anchor_hybrid->metrics.via_count)
              * 0.0036)));

      long hyper_total_via_rise = 0L;
      long hyper_total_wl_gain = 0L;
      long hyper_total_via_delta = 0L;
      long hyper_total_high_layer_drop = 0L;
      int hyper_swapped_nets = 0;
      int hyper_absdeep_imports = 0;
      int hyper_min_imports = 0;
      int hyper_balance_imports = 0;
      int hyper_layerbudget_imports = 0;
      int hyper_radical_imports = 0;
      int hyper_longhpwl_imports = 0;
      int hyper_budgeted_imports = 0;
      int hyper_anchordeep_imports = 0;
      int hyper_corridor_imports = 0;
      int hyper_adaptive_imports = 0;
      int hyper_drstable_imports = 0;
      int hyper_via_clamped_nets = 0;

      for (auto& [db_net, selected_route] : hyper_absolute_wl_result.routes) {
        if (db_net == nullptr) {
          continue;
        }

        const RouteEdgeStats anchor_stats
            = collectRouteEdgeStats(selected_route, grouter_->grid_);
        const long anchor_hpwl = getRouteBBoxHpwl(selected_route);
        const long anchor_detour
            = std::max(0L, anchor_stats.wirelength_dbu - anchor_hpwl);

        const GRoute* best_route = &selected_route;
        RouteEdgeStats best_stats = anchor_stats;
        long best_detour = anchor_detour;
        long best_wl_gain = 0L;
        long best_via_delta = 0L;
        long best_layer_drop = 0L;
        int best_donor = 0;

        auto consider_hyper_donor = [&](const ScenarioResult* donor, int donor_id) {
          if (donor == nullptr) {
            return;
          }
          const auto donor_it = donor->routes.find(db_net);
          if (donor_it == donor->routes.end()) {
            return;
          }

          const RouteEdgeStats donor_stats
              = collectRouteEdgeStats(donor_it->second, grouter_->grid_);
          const long wl_gain = anchor_stats.wirelength_dbu - donor_stats.wirelength_dbu;
          const long min_wl_gain = std::max<long>(4L, tile_size / 3L);
          if (wl_gain < min_wl_gain) {
            return;
          }

          const long via_delta = static_cast<long>(donor_stats.via_count)
                                 - static_cast<long>(anchor_stats.via_count);
          const long via_rise_soft_cap
              = std::max<long>(10L, wl_gain / hyper_trade_via_to_wl + 6L);
          if (via_delta > via_rise_soft_cap) {
            return;
          }

          const long donor_hpwl = getRouteBBoxHpwl(donor_it->second);
          const long donor_detour
              = std::max(0L, donor_stats.wirelength_dbu - donor_hpwl);
          if (donor_detour > anchor_detour + hyper_detour_guard
              && wl_gain < std::max<long>(tile_size * 4L, 1800L)) {
            return;
          }

          const long layer_drop = std::max(
              0L, anchor_stats.high_layer_dbu - donor_stats.high_layer_dbu);
          if (layer_drop > hyper_per_net_layer_drop_cap && wl_gain < min_wl_gain * 8L) {
            return;
          }
          if (hyper_total_high_layer_drop + layer_drop > hyper_high_layer_drop_budget
              && wl_gain < min_wl_gain * 10L) {
            return;
          }

          const bool better_wl
              = donor_stats.wirelength_dbu < best_stats.wirelength_dbu;
          const bool tie_break_via
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count < best_stats.via_count;
          const bool tie_break_layer
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count == best_stats.via_count
                && donor_stats.high_layer_dbu > best_stats.high_layer_dbu;
          const bool tie_break_detour
              = donor_stats.wirelength_dbu == best_stats.wirelength_dbu
                && donor_stats.via_count == best_stats.via_count
                && donor_detour < best_detour;
          if (!(better_wl || tie_break_via || tie_break_layer
                || tie_break_detour)) {
            return;
          }

          best_route = &donor_it->second;
          best_stats = donor_stats;
          best_detour = donor_detour;
          best_wl_gain = wl_gain;
          best_via_delta = via_delta;
          best_layer_drop = layer_drop;
          best_donor = donor_id;
        };

        consider_hyper_donor(absolute_wl_deep_hybrid, 1);
        consider_hyper_donor(min_wl_wide_hybrid, 2);
        consider_hyper_donor(cugr_sp_balance_hybrid, 3);
        consider_hyper_donor(wl_layerbudget_hybrid, 4);
        consider_hyper_donor(radical_wl_hybrid, 5);
        consider_hyper_donor(longhpwl_radical_hybrid, 6);
        consider_hyper_donor(budgeted_lift_wl_hybrid, 7);
        consider_hyper_donor(anchor_wl_deep_hybrid, 8);
        consider_hyper_donor(wl_corridor_hybrid, 9);
        consider_hyper_donor(length_adaptive_hybrid, 10);
        consider_hyper_donor(dr_stable_hybrid, 11);

        if (best_route == &selected_route) {
          continue;
        }
        if (best_via_delta > 0
            && hyper_total_via_rise + best_via_delta > hyper_total_via_rise_cap) {
          hyper_via_clamped_nets++;
          continue;
        }

        selected_route = *best_route;
        hyper_swapped_nets++;
        if (best_via_delta > 0) {
          hyper_total_via_rise += best_via_delta;
        }
        hyper_total_wl_gain += best_wl_gain;
        hyper_total_via_delta += best_via_delta;
        hyper_total_high_layer_drop += best_layer_drop;
        if (best_donor == 1) {
          hyper_absdeep_imports++;
        } else if (best_donor == 2) {
          hyper_min_imports++;
        } else if (best_donor == 3) {
          hyper_balance_imports++;
        } else if (best_donor == 4) {
          hyper_layerbudget_imports++;
        } else if (best_donor == 5) {
          hyper_radical_imports++;
        } else if (best_donor == 6) {
          hyper_longhpwl_imports++;
        } else if (best_donor == 7) {
          hyper_budgeted_imports++;
        } else if (best_donor == 8) {
          hyper_anchordeep_imports++;
        } else if (best_donor == 9) {
          hyper_corridor_imports++;
        } else if (best_donor == 10) {
          hyper_adaptive_imports++;
        } else if (best_donor == 11) {
          hyper_drstable_imports++;
        }
      }

      if (hyper_swapped_nets > 0) {
        hyper_absolute_wl_result.metrics
            = compute_metrics(hyper_absolute_wl_result.routes);
        apply_routability_proxy(hyper_absolute_wl_result.metrics);
        logger_->info(
            GNR,
            7363,
            "NEWGR hybrid-netmix-hyper-absolute-wl imported "
            "absdeep/min/balance/layerbudget/radical/longhpwl/budgeted/"
            "anchordeep/corridor/adaptive/drstable nets "
            "{}/{}/{}/{}/{}/{}/{}/{}/{}/{}/{}, wirelength {:.0f} um, vias {}, "
            "swapped nets {}, via-rise budget {}/{}, high-layer drop {}/{}, wl "
            "gain sum {}, via delta sum {}, anchor '{}', budget-clamped nets {}",
            hyper_absdeep_imports,
            hyper_min_imports,
            hyper_balance_imports,
            hyper_layerbudget_imports,
            hyper_radical_imports,
            hyper_longhpwl_imports,
            hyper_budgeted_imports,
            hyper_anchordeep_imports,
            hyper_corridor_imports,
            hyper_adaptive_imports,
            hyper_drstable_imports,
            hyper_absolute_wl_result.metrics.wirelength_um,
            hyper_absolute_wl_result.metrics.via_count,
            hyper_swapped_nets,
            hyper_total_via_rise,
            hyper_total_via_rise_cap,
            hyper_total_high_layer_drop,
            hyper_high_layer_drop_budget,
            hyper_total_wl_gain,
            hyper_total_via_delta,
            hyper_anchor_hybrid->name,
            hyper_via_clamped_nets);
        scenario_results.push_back(std::move(hyper_absolute_wl_result));
      }
    }

    // Rebind pointers after appending synthetic scenarios.
    wl_hybrid = find_scenario("hybrid-netmix-wl");
    dr_stable_hybrid = find_scenario("hybrid-netmix-dr-stable");
    layer_compact_hybrid = find_scenario("hybrid-netmix-layer-compact");
    via_band_hybrid = find_scenario("hybrid-netmix-via-band-wl");
    anchor_wl_hybrid = find_scenario("hybrid-netmix-anchor-wl");
    anchor_wl_deep_hybrid = find_scenario("hybrid-netmix-anchor-wl-deep");
    wl_corridor_hybrid = find_scenario("hybrid-netmix-wl-corridor");
    length_adaptive_hybrid = find_scenario("hybrid-netmix-length-adaptive");

    // CUGR-inspired post routing patching:
    // start from wirelength-first hybrid and add local endpoint and hotspot
    // guide patches taken from alternate hybrids.
    const ScenarioResult* patch_alt_a
        = dr_stable_hybrid != nullptr ? dr_stable_hybrid : anchor_wl_deep_hybrid;
    if (patch_alt_a == nullptr) {
      patch_alt_a = anchor_wl_hybrid != nullptr ? anchor_wl_hybrid : via_band_hybrid;
    }
    const ScenarioResult* patch_alt_b = layer_compact_hybrid != nullptr
                                            ? layer_compact_hybrid
                                            : length_adaptive_hybrid;
    if (patch_alt_b == nullptr) {
      patch_alt_b = wl_corridor_hybrid != nullptr ? wl_corridor_hybrid
                                                  : via_band_hybrid;
    }
    if (wl_hybrid != nullptr) {
      ScenarioResult patched_result;
      patched_result.name = "hybrid-netmix-cugr-patched";
      patched_result.routes = wl_hybrid->routes;
      const GuidePatchStats patch_stats = applyCugrGuidePatches(
          patched_result.routes,
          patch_alt_a != nullptr ? &patch_alt_a->routes : nullptr,
          patch_alt_b != nullptr ? &patch_alt_b->routes : nullptr,
          hotspots,
          grouter_->grid_,
          min_routing_layer,
          max_routing_layer);

      if (patch_stats.nets_touched > 0) {
        patched_result.metrics = compute_metrics(patched_result.routes);
        apply_routability_proxy(patched_result.metrics);
        const double patch_credit = std::min(
            3500.0, static_cast<double>(patch_stats.nets_touched) * 0.12);
        patched_result.metrics.score
            = std::max(0.0, patched_result.metrics.score - patch_credit);
        logger_->info(GNR,
                      6030,
                      "NEWGR hybrid-netmix-cugr-patched from wl hybrid: "
                      "patched nets {}, endpoint guides {}, long patches {}, "
                      "patch vias {}, wirelength {:.0f} um, vias {}, detour {}",
                      patch_stats.nets_touched,
                      patch_stats.endpoint_segments_added,
                      patch_stats.long_segment_patches,
                      patch_stats.via_patches_added,
                      patched_result.metrics.wirelength_um,
                      patched_result.metrics.via_count,
                      patched_result.metrics.detour_dbu);
        scenario_results.push_back(std::move(patched_result));
      } else {
        logger_->info(
            GNR, 6031, "NEWGR skipped hybrid-netmix-cugr-patched (no patches).");
      }
    }
  }

  auto robust_better = [](const ScenarioResult& lhs,
                          const ScenarioResult& rhs) {
    if (lhs.metrics.overflow_edges != rhs.metrics.overflow_edges) {
      return lhs.metrics.overflow_edges < rhs.metrics.overflow_edges;
    }
    if (lhs.metrics.overflow_ratio_sum != rhs.metrics.overflow_ratio_sum) {
      return lhs.metrics.overflow_ratio_sum < rhs.metrics.overflow_ratio_sum;
    }
    if (lhs.metrics.score != rhs.metrics.score) {
      return lhs.metrics.score < rhs.metrics.score;
    }
    if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
      return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
    }
    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    if (lhs.metrics.high_layer_dbu != rhs.metrics.high_layer_dbu) {
      return lhs.metrics.high_layer_dbu < rhs.metrics.high_layer_dbu;
    }
    if (lhs.metrics.layer_span_sum != rhs.metrics.layer_span_sum) {
      return lhs.metrics.layer_span_sum < rhs.metrics.layer_span_sum;
    }
    if (lhs.metrics.near_capacity_edges != rhs.metrics.near_capacity_edges) {
      return lhs.metrics.near_capacity_edges < rhs.metrics.near_capacity_edges;
    }
    return lhs.metrics.max_usage_ratio < rhs.metrics.max_usage_ratio;
  };

  const auto shortest_wl_iter = std::min_element(
      scenario_results.begin(),
      scenario_results.end(),
      [](const ScenarioResult& lhs, const ScenarioResult& rhs) {
        return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
      });
  const long shortest_wl = shortest_wl_iter->metrics.wirelength_dbu;
  const bool overflow_free_sweep = std::all_of(
      scenario_results.begin(), scenario_results.end(), [](const auto& result) {
        return result.metrics.overflow_edges == 0;
      });
  const double wl_guard_ratio = overflow_free_sweep ? 0.0038 : 0.0100;
  const long wl_guard_floor = overflow_free_sweep ? 170 : 220;
  const long wl_guard_band = std::max<long>(
      wl_guard_floor,
      static_cast<long>(std::ceil(wl_guard_ratio * static_cast<double>(shortest_wl))));

  std::vector<const ScenarioResult*> shortlist;
  shortlist.reserve(scenario_results.size());
  for (const ScenarioResult& result : scenario_results) {
    if (result.metrics.wirelength_dbu <= shortest_wl + wl_guard_band) {
      shortlist.push_back(&result);
    }
  }
  if (shortlist.empty()) {
    for (const ScenarioResult& result : scenario_results) {
      shortlist.push_back(&result);
    }
  }

  auto wirelength_first_better = [](const ScenarioResult& lhs,
                                    const ScenarioResult& rhs) {
    if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
      return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
    }
    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    return lhs.metrics.score < rhs.metrics.score;
  };

  const long tie_via_wl_band
      = std::max<long>(28, static_cast<long>(std::ceil(shortest_wl * 0.00011)));
  const long tie_proxy_wl_band
      = std::max<long>(220, static_cast<long>(std::ceil(shortest_wl * 0.0052)));
  auto wirelength_with_dr_proxy_tie_better = [&](const ScenarioResult& lhs,
                                                  const ScenarioResult& rhs) {
    const long wl_gap = std::llabs(lhs.metrics.wirelength_dbu
                                   - rhs.metrics.wirelength_dbu);
    if (wl_gap == 0 && lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    if (wl_gap <= tie_proxy_wl_band) {
      const double lhs_proxy = estimateDetailedRouteProxyCost(lhs.metrics);
      const double rhs_proxy = estimateDetailedRouteProxyCost(rhs.metrics);
      if (std::abs(lhs_proxy - rhs_proxy) > 1e-3) {
        return lhs_proxy < rhs_proxy;
      }
      if (lhs.metrics.high_layer_dbu != rhs.metrics.high_layer_dbu) {
        return lhs.metrics.high_layer_dbu < rhs.metrics.high_layer_dbu;
      }
      if (lhs.metrics.detour_dbu != rhs.metrics.detour_dbu) {
        return lhs.metrics.detour_dbu < rhs.metrics.detour_dbu;
      }
    }
    if (wl_gap <= tie_via_wl_band && lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    return wirelength_first_better(lhs, rhs);
  };

  const ScenarioResult* forced_wl_ptr = nullptr;
  if (overflow_free_sweep) {
    auto find_scenario_by_name = [&](const std::string& name) {
      const ScenarioResult* match = nullptr;
      for (const ScenarioResult& result : scenario_results) {
        if (result.name == name) {
          match = &result;
          break;
        }
      }
      return match;
    };

    const ScenarioResult* wl_anchor = find_scenario_by_name("hybrid-netmix-wl");
    const ScenarioResult* wl_compact_ptr
        = find_scenario_by_name("hybrid-netmix-wl-compact");
    const ScenarioResult* wl_corridor_ptr
        = find_scenario_by_name("hybrid-netmix-wl-corridor");
    const ScenarioResult* patched_ptr
        = find_scenario_by_name("hybrid-netmix-cugr-patched");
    const ScenarioResult* wl_feedback_ptr
        = find_scenario_by_name("hybrid-netmix-wl-feedback");
    const ScenarioResult* length_adaptive_ptr
        = find_scenario_by_name("hybrid-netmix-length-adaptive");
    const ScenarioResult* via_floor_ptr
        = find_scenario_by_name("hybrid-netmix-via-floor-wl");
    const ScenarioResult* layer_lift_wl_ptr
        = find_scenario_by_name("hybrid-netmix-layer-lift-wl");
    const ScenarioResult* anchor_wl_deep_ptr
        = find_scenario_by_name("hybrid-netmix-anchor-wl-deep");
    const ScenarioResult* layer_floor_wl_ptr
        = find_scenario_by_name("hybrid-netmix-anchor-wl-layerfloor");
    const ScenarioResult* drt_elastic_ptr
        = find_scenario_by_name("hybrid-netmix-drt-elastic-wl");
    const ScenarioResult* budgeted_lift_ptr
        = find_scenario_by_name("hybrid-netmix-budgeted-lift-wl");
    const ScenarioResult* min_wl_wide_ptr
        = find_scenario_by_name("hybrid-netmix-min-wl-wide");
    const ScenarioResult* absolute_wl_ptr
        = find_scenario_by_name("hybrid-netmix-absolute-wl");
    const ScenarioResult* absolute_wl_deep_ptr
        = find_scenario_by_name("hybrid-netmix-absolute-wl-deep");
    const ScenarioResult* cugr_sp_balance_ptr
        = find_scenario_by_name("hybrid-netmix-cugr-sp-balance-wl");
    const ScenarioResult* wl_layerbudget_ptr
        = find_scenario_by_name("hybrid-netmix-wl-layerbudget");
    const ScenarioResult* radical_wl_fusion_ptr
        = find_scenario_by_name("hybrid-netmix-radical-wl-fusion");
    const ScenarioResult* absolute_via_recover_ptr
        = find_scenario_by_name("hybrid-netmix-absolute-via-recover");
    const ScenarioResult* hyper_absolute_wl_ptr
        = find_scenario_by_name("hybrid-netmix-hyper-absolute-wl");
    const ScenarioResult* preferred_wl_ptr = nullptr;

    if (wl_anchor != nullptr && wl_layerbudget_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long wl_gain
          = wl_anchor->metrics.wirelength_dbu - wl_layerbudget_ptr->metrics.wirelength_dbu;
      const long min_wl_gain = std::max<long>(
          22L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00008)));
      const long via_rise = static_cast<long>(wl_layerbudget_ptr->metrics.via_count)
                            - static_cast<long>(wl_anchor->metrics.via_count);
      const long via_rise_cap = std::max<long>(
          42L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.00042)));
      const long detour_guard = std::max<long>(tile_size * 4L, 2200L);
      const long high_layer_floor = std::max<long>(
          tile_size * 14L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.985)));
      const bool compact_guard
          = wl_gain >= min_wl_gain && via_rise >= 0 && via_rise <= via_rise_cap
            && wl_layerbudget_ptr->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
            && wl_layerbudget_ptr->metrics.high_layer_dbu >= high_layer_floor
            && wl_layerbudget_ptr->metrics.overflow_edges
                   <= wl_anchor->metrics.overflow_edges;
      if (compact_guard) {
        preferred_wl_ptr = wl_layerbudget_ptr;
        logger_->info(
            GNR,
            7346,
            "NEWGR pre-selecting layer-budget '{}': wl gain {}, via rise {}, "
            "detour delta {}, high-layer delta {}.",
            preferred_wl_ptr->name,
            wl_gain,
            via_rise,
            preferred_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            preferred_wl_ptr->metrics.high_layer_dbu
                - wl_anchor->metrics.high_layer_dbu);
      }
    }

    // Radical WL-first override:
    // if an extreme min-WL hybrid is a strict WL+via improvement over the
    // anchor while staying overflow-safe, lock it in as the preferred
    // candidate before conservative proxy gates.
    if (wl_anchor != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long min_wl_gain = std::max<long>(
          80L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00030)));
      const long via_rise_cap = std::max<long>(
          320L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0035)));
      const long via_drop_cap = std::max<long>(
          960L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0095)));
      const long detour_guard = std::max<long>(tile_size * 22L, 11000L);
      const long high_layer_guard = std::max<long>(
          tile_size * 34L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.10)));
      const long min_high_layer_guard = std::max<long>(
          tile_size * 24L,
          static_cast<long>(std::ceil(
              static_cast<double>(std::max(1L, wl_anchor->metrics.high_layer_dbu))
              * 0.16)));
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double proxy_cap = anchor_proxy * 1.025;

      const ScenarioResult* radical_candidate = nullptr;
      for (const ScenarioResult* candidate :
           std::array<const ScenarioResult*, 10>{budgeted_lift_ptr,
                                                layer_floor_wl_ptr,
                                                layer_lift_wl_ptr,
                                                absolute_wl_ptr,
                                                absolute_wl_deep_ptr,
                                                min_wl_wide_ptr,
                                                cugr_sp_balance_ptr,
                                                wl_layerbudget_ptr,
                                                radical_wl_fusion_ptr,
                                                hyper_absolute_wl_ptr}) {
        if (candidate == nullptr) {
          continue;
        }
        if (candidate->metrics.overflow_edges > wl_anchor->metrics.overflow_edges) {
          continue;
        }
        const long wl_gain
            = wl_anchor->metrics.wirelength_dbu - candidate->metrics.wirelength_dbu;
        if (wl_gain < min_wl_gain) {
          continue;
        }
        const long via_delta = static_cast<long>(candidate->metrics.via_count)
                               - static_cast<long>(wl_anchor->metrics.via_count);
        if (via_delta > via_rise_cap || via_delta < -via_drop_cap) {
          continue;
        }
        if (candidate->metrics.high_layer_dbu + min_high_layer_guard
            < wl_anchor->metrics.high_layer_dbu) {
          continue;
        }
        if (candidate->metrics.detour_dbu > wl_anchor->metrics.detour_dbu + detour_guard
            || candidate->metrics.high_layer_dbu
                   > wl_anchor->metrics.high_layer_dbu + high_layer_guard) {
          continue;
        }
        const double candidate_proxy
            = estimateDetailedRouteProxyCost(candidate->metrics);
        if (candidate_proxy > proxy_cap) {
          continue;
        }
        if (radical_candidate == nullptr
            || wirelength_first_better(*candidate, *radical_candidate)) {
          radical_candidate = candidate;
        }
      }

      if (radical_candidate != nullptr) {
        if (preferred_wl_ptr == wl_layerbudget_ptr && wl_layerbudget_ptr != nullptr
            && radical_candidate != wl_layerbudget_ptr) {
          const long layerbudget_wl_gain
              = wl_anchor->metrics.wirelength_dbu
                - wl_layerbudget_ptr->metrics.wirelength_dbu;
          const long layerbudget_via_rise
              = static_cast<long>(wl_layerbudget_ptr->metrics.via_count)
                - static_cast<long>(wl_anchor->metrics.via_count);
          const long layerbudget_via_rise_cap = std::max<long>(
              42L,
              static_cast<long>(std::ceil(
                  static_cast<double>(wl_anchor->metrics.via_count) * 0.00042)));
          const long layerbudget_detour_guard = std::max<long>(tile_size * 4L, 2200L);
          const long layerbudget_high_layer_floor = std::max<long>(
              tile_size * 14L,
              static_cast<long>(std::ceil(
                  static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.985)));
          const bool preserve_layerbudget_mix
              = layerbudget_wl_gain >= min_wl_gain / 2L && layerbudget_via_rise >= 0
                && layerbudget_via_rise <= layerbudget_via_rise_cap
                && wl_layerbudget_ptr->metrics.detour_dbu
                       <= wl_anchor->metrics.detour_dbu + layerbudget_detour_guard
                && wl_layerbudget_ptr->metrics.high_layer_dbu
                       >= layerbudget_high_layer_floor
                && wl_layerbudget_ptr->metrics.overflow_edges
                       <= wl_anchor->metrics.overflow_edges;
          const long radical_wl_advantage
              = wl_layerbudget_ptr->metrics.wirelength_dbu
                - radical_candidate->metrics.wirelength_dbu;
          const bool radical_minwl_dominates
              = radical_wl_advantage >= std::max<long>(tile_size * 8L, 18000L);
          if (preserve_layerbudget_mix && !radical_minwl_dominates) {
            radical_candidate = wl_layerbudget_ptr;
          } else if (preserve_layerbudget_mix && radical_minwl_dominates) {
            logger_->info(
                GNR,
                7358,
                "NEWGR bypassing layer-budget preserve in favor of '{}' "
                "(WL advantage over layer-budget {}, threshold {}).",
                radical_candidate->name,
                radical_wl_advantage,
                std::max<long>(tile_size * 8L, 18000L));
          }
        }
        preferred_wl_ptr = radical_candidate;
        logger_->info(
            GNR,
            6039,
            "NEWGR radical WL lock pre-selecting '{}' over anchor '{}' "
            "(wl gain {}, via gain {}, detour delta {}, high-layer delta {}).",
            preferred_wl_ptr->name,
            wl_anchor->name,
            wl_anchor->metrics.wirelength_dbu
                - preferred_wl_ptr->metrics.wirelength_dbu,
            wl_anchor->metrics.via_count - preferred_wl_ptr->metrics.via_count,
            preferred_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            preferred_wl_ptr->metrics.high_layer_dbu
                - wl_anchor->metrics.high_layer_dbu);
      }
    }

    if (wl_anchor != nullptr && wl_compact_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long wl_guard = std::max<long>(
          64L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00028)));
      const long via_guard = std::max<long>(
          140L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0014)));
      const long detour_guard = std::max<long>(5500L, tile_size * 8L);
      const long high_layer_guard = std::max<long>(
          tile_size * 18L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.050)));
      const bool within_structural_guard
          = wl_compact_ptr->metrics.wirelength_dbu
                 <= wl_anchor->metrics.wirelength_dbu + wl_guard
            && wl_compact_ptr->metrics.via_count
                   <= wl_anchor->metrics.via_count + via_guard
            && wl_compact_ptr->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
            && wl_compact_ptr->metrics.high_layer_dbu
                   <= wl_anchor->metrics.high_layer_dbu + high_layer_guard;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double compact_proxy
          = estimateDetailedRouteProxyCost(wl_compact_ptr->metrics);
      const bool proxy_win = compact_proxy + 1e-3 < anchor_proxy * 0.994;
      if (within_structural_guard && proxy_win) {
        logger_->info(
            GNR,
            7325,
            "NEWGR compact-layer anchor '{}' replacing '{}' "
            "(wl delta {}, via delta {}, detour delta {}, high-layer delta {}, "
            "proxy ratio {:.3f}).",
            wl_compact_ptr->name,
            wl_anchor->name,
            wl_compact_ptr->metrics.wirelength_dbu - wl_anchor->metrics.wirelength_dbu,
            wl_compact_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            wl_compact_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            wl_compact_ptr->metrics.high_layer_dbu - wl_anchor->metrics.high_layer_dbu,
            anchor_proxy > 1e-9 ? compact_proxy / anchor_proxy : 1.0);
        wl_anchor = wl_compact_ptr;
      }
    }

    if (wl_anchor != nullptr && wl_corridor_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long min_wl_gain = std::max<long>(
          24L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00009)));
      const long detour_guard = std::max<long>(tile_size * 10L, 6200L);
      const long high_layer_guard = std::max<long>(
          tile_size * 20L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.060)));
      const bool structural_guard
          = wl_corridor_ptr->metrics.detour_dbu
                 <= wl_anchor->metrics.detour_dbu + detour_guard
            && wl_corridor_ptr->metrics.high_layer_dbu
                   <= wl_anchor->metrics.high_layer_dbu + high_layer_guard;
      const long wl_gain
          = wl_anchor->metrics.wirelength_dbu - wl_corridor_ptr->metrics.wirelength_dbu;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double corridor_proxy
          = estimateDetailedRouteProxyCost(wl_corridor_ptr->metrics);
      const bool proxy_guard = corridor_proxy + 1e-3 < anchor_proxy * 1.006;
      if (wl_gain >= min_wl_gain && structural_guard && proxy_guard) {
        logger_->info(
            GNR,
            7331,
            "NEWGR corridor-locked anchor '{}' replacing '{}' "
            "(wl gain {}, via delta {}, detour delta {}, high-layer delta {}, "
            "proxy ratio {:.3f}).",
            wl_corridor_ptr->name,
            wl_anchor->name,
            wl_gain,
            wl_corridor_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            wl_corridor_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            wl_corridor_ptr->metrics.high_layer_dbu - wl_anchor->metrics.high_layer_dbu,
            anchor_proxy > 1e-9 ? corridor_proxy / anchor_proxy : 1.0);
        wl_anchor = wl_corridor_ptr;
      }
    }

    const long via_drop_guard = wl_anchor != nullptr
                                    ? std::max<long>(
                                          360L,
                                          static_cast<long>(std::ceil(
                                              static_cast<double>(
                                                  wl_anchor->metrics.via_count)
                                              * 0.0085)))
                                    : 0L;
    const long deep_via_drop_wl_gain = wl_anchor != nullptr
                                           ? std::max<long>(
                                                 180L,
                                                 static_cast<long>(std::ceil(
                                                     static_cast<double>(
                                                         wl_anchor->metrics.wirelength_dbu)
                                                     * 0.00070)))
                                           : 0L;
    const long deep_via_drop_detour_bonus = std::max<long>(
        9000L, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1) * 18L));
    const long deep_via_drop_high_layer_bonus = wl_anchor != nullptr
                                                    ? std::max<long>(
                                                          2200000L,
                                                          static_cast<long>(std::ceil(
                                                              static_cast<double>(
                                                                  wl_anchor->metrics.high_layer_dbu)
                                                              * 0.016)))
                                                    : 0L;
    auto allow_aggressive_via_drop = [&](const ScenarioResult* candidate) {
      if (wl_anchor == nullptr || candidate == nullptr) {
        return true;
      }
      const long via_drop = wl_anchor->metrics.via_count - candidate->metrics.via_count;
      if (via_drop <= via_drop_guard) {
        return true;
      }

      const long wl_gain
          = wl_anchor->metrics.wirelength_dbu - candidate->metrics.wirelength_dbu;
      const long detour_delta
          = candidate->metrics.detour_dbu - wl_anchor->metrics.detour_dbu;
      const long high_layer_delta
          = candidate->metrics.high_layer_dbu - wl_anchor->metrics.high_layer_dbu;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double candidate_proxy
          = estimateDetailedRouteProxyCost(candidate->metrics);

      // Accept deep via drops when the candidate is still a structural Pareto
      // win with lower WL, lower vias, and improved DR proxy quality.
      const bool strict_pareto_upgrade
          = candidate->metrics.wirelength_dbu <= wl_anchor->metrics.wirelength_dbu
            && candidate->metrics.via_count <= wl_anchor->metrics.via_count
            && detour_delta <= deep_via_drop_detour_bonus
            && high_layer_delta <= deep_via_drop_high_layer_bonus
            && candidate_proxy + 1e-3 < anchor_proxy * 1.025;
      if (strict_pareto_upgrade) {
        return true;
      }

      const bool structural_recovery
          = candidate->metrics.detour_dbu + deep_via_drop_detour_bonus
                <= wl_anchor->metrics.detour_dbu
            && candidate->metrics.high_layer_dbu + deep_via_drop_high_layer_bonus
                   <= wl_anchor->metrics.high_layer_dbu;
      const bool proxy_recovery = candidate_proxy + 1e-3 < anchor_proxy * 1.010;
      return wl_gain >= deep_via_drop_wl_gain && structural_recovery
             && proxy_recovery;
    };

    if (preferred_wl_ptr == nullptr && wl_anchor != nullptr
        && wl_feedback_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long wl_guard = std::max<long>(
          48,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00020)));
      const long via_gain_req = 220L;
      const long detour_guard = std::max<long>(5000L, tile_size * 2L);
      const bool within_wl_guard
          = wl_feedback_ptr->metrics.wirelength_dbu
            <= wl_anchor->metrics.wirelength_dbu + wl_guard;
      const bool via_gain
          = wl_feedback_ptr->metrics.via_count + via_gain_req
            <= wl_anchor->metrics.via_count;
      const bool detour_ok
          = wl_feedback_ptr->metrics.detour_dbu
            <= wl_anchor->metrics.detour_dbu + detour_guard;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double feedback_proxy
          = estimateDetailedRouteProxyCost(wl_feedback_ptr->metrics);
      const bool proxy_ok = feedback_proxy + 1e-3 < anchor_proxy * 1.010;

      if (within_wl_guard && via_gain && detour_ok && proxy_ok) {
        preferred_wl_ptr = wl_feedback_ptr;
        logger_->info(
            GNR,
            7315,
            "NEWGR pre-selecting '{}' over '{}' in overflow-free mode "
            "(wl delta {}, via delta {}, detour delta {}, proxy ratio {:.3f}).",
            preferred_wl_ptr->name,
            wl_anchor->name,
            preferred_wl_ptr->metrics.wirelength_dbu - wl_anchor->metrics.wirelength_dbu,
            preferred_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            preferred_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            anchor_proxy > 1e-9 ? feedback_proxy / anchor_proxy : 1.0);
      }
    }

    if (preferred_wl_ptr == nullptr && wl_anchor != nullptr
        && length_adaptive_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long anchor_wl_guard = std::max<long>(
          160,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00065)));
      const long via_guard = std::max<long>(40L, tile_size * 2L);
      const long detour_guard = std::max<long>(tile_size * 12L, 6000L);
      const bool within_wl_guard
          = length_adaptive_ptr->metrics.wirelength_dbu
            <= wl_anchor->metrics.wirelength_dbu + anchor_wl_guard;
      const bool within_via_guard
          = length_adaptive_ptr->metrics.via_count
            <= wl_anchor->metrics.via_count + via_guard;
      const bool within_detour_guard
          = length_adaptive_ptr->metrics.detour_dbu
            <= wl_anchor->metrics.detour_dbu + detour_guard;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double adaptive_proxy
          = estimateDetailedRouteProxyCost(length_adaptive_ptr->metrics);
      const bool proxy_better = adaptive_proxy + 1e-3 < anchor_proxy * 0.985;

      if (within_wl_guard && within_via_guard && within_detour_guard
          && proxy_better) {
        preferred_wl_ptr = length_adaptive_ptr;
        logger_->info(
            GNR,
            7312,
            "NEWGR pre-selecting '{}' over '{}' in overflow-free mode "
            "(wl delta {}, via delta {}, proxy ratio {:.3f}).",
            preferred_wl_ptr->name,
            wl_anchor->name,
            preferred_wl_ptr->metrics.wirelength_dbu - wl_anchor->metrics.wirelength_dbu,
            preferred_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            anchor_proxy > 1e-9 ? adaptive_proxy / anchor_proxy : 1.0);
      }
    }

    if (preferred_wl_ptr == nullptr && wl_anchor != nullptr
        && anchor_wl_deep_ptr != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long wl_guard = std::max<long>(
          90L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00035)));
      const long via_bonus_floor = std::max<long>(96L, tile_size * 4L);
      const long via_bonus_cap = std::max<long>(
          1200L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.012)));
      const long detour_guard = std::max<long>(tile_size * 16L, 8500L);
      const long high_layer_guard = std::max<long>(
          tile_size * 26L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.060)));
      const bool within_wl_guard
          = anchor_wl_deep_ptr->metrics.wirelength_dbu
            <= wl_anchor->metrics.wirelength_dbu + wl_guard;
      const long via_delta
          = static_cast<long>(anchor_wl_deep_ptr->metrics.via_count)
            - static_cast<long>(wl_anchor->metrics.via_count);
      const bool within_via_window
          = via_delta >= via_bonus_floor && via_delta <= via_bonus_cap;
      const bool structural_guard
          = anchor_wl_deep_ptr->metrics.detour_dbu
                 <= wl_anchor->metrics.detour_dbu + detour_guard
            && anchor_wl_deep_ptr->metrics.high_layer_dbu
                   <= wl_anchor->metrics.high_layer_dbu + high_layer_guard;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double deep_anchor_proxy
          = estimateDetailedRouteProxyCost(anchor_wl_deep_ptr->metrics);
      const double via_elastic_credit
          = std::min(2500.0, std::max(0L, via_delta) * 0.90);
      const double elastic_proxy = deep_anchor_proxy - via_elastic_credit;
      const bool proxy_better = elastic_proxy + 1e-3 < anchor_proxy * 0.997;

      if (within_wl_guard && within_via_window && structural_guard && proxy_better) {
        preferred_wl_ptr = anchor_wl_deep_ptr;
        logger_->info(
            GNR,
            7323,
            "NEWGR pre-selecting '{}' over '{}' in via-elastic mode "
            "(wl delta {}, via delta {}, proxy ratio {:.3f}, elastic proxy ratio {:.3f}).",
            preferred_wl_ptr->name,
            wl_anchor->name,
            preferred_wl_ptr->metrics.wirelength_dbu - wl_anchor->metrics.wirelength_dbu,
            preferred_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            anchor_proxy > 1e-9 ? deep_anchor_proxy / anchor_proxy : 1.0,
            anchor_proxy > 1e-9 ? elastic_proxy / anchor_proxy : 1.0);
      }
    }

    // Wirelength champion pool (mixing aggressive FastRoute-like shortest-path
    // hybrids with low-via variants): rank by WL, but use a CUGR-style
    // detailed-routability proxy tie-break in a narrow WL band.
    auto wl_champion_better = [&](const ScenarioResult& lhs,
                                  const ScenarioResult& rhs) {
      const long wl_gap = std::llabs(lhs.metrics.wirelength_dbu
                                     - rhs.metrics.wirelength_dbu);
      const long dr_tie_band
          = std::max<long>(96, static_cast<long>(std::ceil(shortest_wl * 0.00042)));
      if (wl_gap <= dr_tie_band) {
        const double lhs_proxy = estimateDetailedRouteProxyCost(lhs.metrics);
        const double rhs_proxy = estimateDetailedRouteProxyCost(rhs.metrics);
        if (std::abs(lhs_proxy - rhs_proxy) > 1e-3) {
          return lhs_proxy < rhs_proxy;
        }
        if (lhs.metrics.high_layer_dbu != rhs.metrics.high_layer_dbu) {
          return lhs.metrics.high_layer_dbu < rhs.metrics.high_layer_dbu;
        }
        if (lhs.metrics.detour_dbu != rhs.metrics.detour_dbu) {
          return lhs.metrics.detour_dbu < rhs.metrics.detour_dbu;
        }
      }
      return wirelength_first_better(lhs, rhs);
    };

    std::vector<const ScenarioResult*> wl_champion_pool;
    wl_champion_pool.reserve(24);
    for (const char* name : std::array<const char*, 24>{
             "hybrid-netmix-anchor-wl-deep",
             "hybrid-netmix-anchor-wl-layerfloor",
             "hybrid-netmix-drt-elastic-wl",
             "hybrid-netmix-budgeted-lift-wl",
             "hybrid-netmix-radical-wl-fusion",
             "hybrid-netmix-hyper-absolute-wl",
             "hybrid-netmix-anchor-wl",
             "hybrid-netmix-wl-corridor",
             "hybrid-netmix-length-adaptive",
             "hybrid-netmix-layer-lift-wl",
             "hybrid-netmix-wl-feedback",
             "hybrid-netmix-via-floor-wl",
             "hybrid-netmix-via-band-wl",
             "hybrid-netmix-smooth-wl",
             "hybrid-netmix-ultra-wl",
             "hybrid-netmix-min-wl-wide",
             "hybrid-netmix-absolute-wl",
             "hybrid-netmix-absolute-wl-deep",
             "hybrid-netmix-cugr-sp-balance-wl",
             "hybrid-netmix-wl-compact",
             "hybrid-netmix-wl-safe",
             "hybrid-netmix-dr-stable",
             "hybrid-netmix-hpwl-lock",
             "hybrid-netmix-wl"}) {
      if (const ScenarioResult* candidate = find_scenario_by_name(name)) {
        wl_champion_pool.push_back(candidate);
      }
    }
    if (!wl_champion_pool.empty()) {
      forced_wl_ptr = *std::min_element(
          wl_champion_pool.begin(),
          wl_champion_pool.end(),
          [&](const ScenarioResult* lhs, const ScenarioResult* rhs) {
            return wl_champion_better(*lhs, *rhs);
          });
    } else {
      forced_wl_ptr = wl_anchor;
    }
    if (compact_exploration_mode) {
      const ScenarioResult* aggressive_wl_ptr = forced_wl_ptr;
      for (const ScenarioResult* candidate : std::array<const ScenarioResult*, 10>{
               budgeted_lift_ptr,
               layer_floor_wl_ptr,
               layer_lift_wl_ptr,
               absolute_wl_ptr,
               absolute_wl_deep_ptr,
               min_wl_wide_ptr,
               cugr_sp_balance_ptr,
               wl_layerbudget_ptr,
               radical_wl_fusion_ptr,
               hyper_absolute_wl_ptr}) {
        if (candidate == nullptr) {
          continue;
        }
        if (!allow_aggressive_via_drop(candidate)) {
          continue;
        }
        if (aggressive_wl_ptr == nullptr
            || wirelength_first_better(*candidate, *aggressive_wl_ptr)) {
          aggressive_wl_ptr = candidate;
        }
      }
      if (aggressive_wl_ptr != nullptr
          && (forced_wl_ptr == nullptr
              || aggressive_wl_ptr != forced_wl_ptr)) {
        logger_->info(
            GNR,
            6034,
            "NEWGR compact mode selecting aggressive WL champion '{}' over '{}'"
            " (wl delta {}, via delta {}).",
            aggressive_wl_ptr->name,
            forced_wl_ptr != nullptr ? forced_wl_ptr->name : std::string("none"),
            forced_wl_ptr != nullptr
                ? aggressive_wl_ptr->metrics.wirelength_dbu
                      - forced_wl_ptr->metrics.wirelength_dbu
                : 0L,
            forced_wl_ptr != nullptr
                ? aggressive_wl_ptr->metrics.via_count
                      - forced_wl_ptr->metrics.via_count
                : 0);
      }
      forced_wl_ptr = aggressive_wl_ptr;
    }
    if (preferred_wl_ptr != nullptr) {
      forced_wl_ptr = preferred_wl_ptr;
    }
    if (preferred_wl_ptr == nullptr && wl_anchor != nullptr
        && !wl_champion_pool.empty()) {
      // FastRoute+SPRoute-style aggressive netmix championing:
      // if a candidate is a structural Pareto improvement over the wl anchor,
      // take it before conservative DR-proxy guards can reject it.
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long detour_guard = std::max<long>(tile_size * 18L, 8500L);
      const long high_layer_guard = std::max<long>(
          tile_size * 28L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.070)));
      const int hotspot_guard = std::max<int>(6, wl_anchor->metrics.near_capacity_edges / 5);
      const long wl_gain_floor = std::max<long>(
          28L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00010)));
      const long via_gain_floor = std::max<long>(
          120L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0012)));
      const long via_trade_guard = std::max<long>(
          64L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0010)));

      const ScenarioResult* pareto_upgrade = nullptr;
      for (const ScenarioResult* candidate : wl_champion_pool) {
        if (candidate == nullptr || candidate == wl_anchor) {
          continue;
        }
        if (!allow_aggressive_via_drop(candidate)) {
          continue;
        }
        const bool structural_guard
            = candidate->metrics.high_layer_dbu
                   <= wl_anchor->metrics.high_layer_dbu + high_layer_guard
              && candidate->metrics.detour_dbu
                     <= wl_anchor->metrics.detour_dbu + detour_guard
              && candidate->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard;
        if (!structural_guard) {
          continue;
        }

        const bool strict_dominates
            = candidate->metrics.wirelength_dbu <= wl_anchor->metrics.wirelength_dbu
              && candidate->metrics.via_count <= wl_anchor->metrics.via_count;
        const bool wl_dominates_with_via_trade
            = candidate->metrics.wirelength_dbu + wl_gain_floor
                     <= wl_anchor->metrics.wirelength_dbu
              && candidate->metrics.via_count
                     <= wl_anchor->metrics.via_count + via_trade_guard;
        const bool via_dominates_with_wl_guard
            = candidate->metrics.via_count + via_gain_floor
                     <= wl_anchor->metrics.via_count
              && candidate->metrics.wirelength_dbu
                     <= wl_anchor->metrics.wirelength_dbu + wl_gain_floor;
        if (!(strict_dominates || wl_dominates_with_via_trade
              || via_dominates_with_wl_guard)) {
          continue;
        }

        if (pareto_upgrade == nullptr
            || wirelength_with_dr_proxy_tie_better(*candidate, *pareto_upgrade)) {
          pareto_upgrade = candidate;
        }
      }

      if (pareto_upgrade != nullptr) {
        forced_wl_ptr = pareto_upgrade;
        logger_->info(
            GNR,
            6033,
            "NEWGR accepted Pareto champion '{}' over wl anchor '{}' "
            "(wl delta {}, via delta {}, detour delta {}, high-layer delta {}).",
            forced_wl_ptr->name,
            wl_anchor->name,
            forced_wl_ptr->metrics.wirelength_dbu - wl_anchor->metrics.wirelength_dbu,
            forced_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            forced_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            forced_wl_ptr->metrics.high_layer_dbu - wl_anchor->metrics.high_layer_dbu);
      }
    }

    if (forced_wl_ptr != nullptr && wl_anchor != nullptr
        && forced_wl_ptr != wl_anchor) {
      const long wl_gain
          = wl_anchor->metrics.wirelength_dbu - forced_wl_ptr->metrics.wirelength_dbu;
      const long min_wl_gain = std::max<long>(
          32,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00010)));
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long base_via_guard = std::max<long>(
          120L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0013)));
      const long wl_scaled_via_credit = std::max<long>(
          0L,
          static_cast<long>(std::ceil(
              static_cast<double>(std::max(0L, wl_gain)) * 0.00004)));
      const long via_guard
          = std::max<long>(base_via_guard, 160L + wl_scaled_via_credit);
      const bool via_guard_ok
          = forced_wl_ptr->metrics.via_count <= wl_anchor->metrics.via_count + via_guard;
      const long detour_guard = std::max<long>(tile_size * 20L, 10000L);
      const long high_layer_guard = std::max<long>(
          tile_size * 28L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.070)));
      const int hotspot_guard = std::max<int>(6, wl_anchor->metrics.near_capacity_edges / 5);
      const bool structural_guard
          = forced_wl_ptr->metrics.high_layer_dbu
                 <= wl_anchor->metrics.high_layer_dbu + high_layer_guard
            && forced_wl_ptr->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
            && forced_wl_ptr->metrics.near_capacity_edges
                   <= wl_anchor->metrics.near_capacity_edges + hotspot_guard;
      const bool via_drop_guard_ok = allow_aggressive_via_drop(forced_wl_ptr);
      const double anchor_proxy
          = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double challenger_proxy
          = estimateDetailedRouteProxyCost(forced_wl_ptr->metrics);
      const bool proxy_guard = challenger_proxy + 1e-3 < anchor_proxy * 1.005;
      const bool proxy_dominant_upgrade
          = wl_gain >= 0 && via_guard_ok && via_drop_guard_ok && structural_guard
            && challenger_proxy + 1e-3 < anchor_proxy * 1.010;
      const long via_gain_floor = std::max<long>(
          120L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0012)));
      const bool via_dominant_upgrade
          = wl_gain >= 0
            && forced_wl_ptr->metrics.via_count + via_gain_floor
                   <= wl_anchor->metrics.via_count
            && via_drop_guard_ok
            && forced_wl_ptr->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
            && challenger_proxy + 1e-3 < anchor_proxy * 1.010;
      const bool equal_wl_via_win
          = (wl_gain == 0)
            && (forced_wl_ptr->metrics.via_count < wl_anchor->metrics.via_count);
      const bool strict_dominates
          = (forced_wl_ptr->metrics.wirelength_dbu <= wl_anchor->metrics.wirelength_dbu)
            && (forced_wl_ptr->metrics.via_count <= wl_anchor->metrics.via_count)
            && structural_guard;
      const bool bypass_via_drop_guard = proxy_dominant_upgrade
                                         || via_dominant_upgrade
                                         || equal_wl_via_win
                                         || strict_dominates;

      if (!((wl_gain >= min_wl_gain && proxy_guard && via_guard_ok
           && structural_guard)
            || proxy_dominant_upgrade
            || via_dominant_upgrade
            || equal_wl_via_win
            || strict_dominates)
          || (!bypass_via_drop_guard && !via_drop_guard_ok)) {
        forced_wl_ptr = wl_anchor;
      } else {
        logger_->info(
            GNR,
            6032,
            "NEWGR overflow-free champion '{}' accepted over '{}' (wl gain {}, "
            "vias delta {}, proxy ratio {:.3f}).",
            forced_wl_ptr->name,
            wl_anchor->name,
            wl_gain,
            forced_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            anchor_proxy > 1e-9 ? challenger_proxy / anchor_proxy : 1.0);
      }
    }
    if (forced_wl_ptr == wl_anchor && wl_anchor != nullptr) {
      // Radical WL-first override: if the absolute/min-WL hybrids
      // simultaneously reduce WL+vias and stay structurally close,
      // prefer them over conservative anchor retention.
      auto candidate_is_viable_wl_upgrade = [&](const ScenarioResult* candidate) {
        if (candidate == nullptr || candidate == wl_anchor) {
          return false;
        }
        if (candidate->metrics.overflow_edges > wl_anchor->metrics.overflow_edges) {
          return false;
        }
        const long wl_gain
            = wl_anchor->metrics.wirelength_dbu - candidate->metrics.wirelength_dbu;
        const long via_gain
            = wl_anchor->metrics.via_count - candidate->metrics.via_count;
        const long min_wl_gain = std::max<long>(
            80L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00030)));
        const long min_via_gain
            = std::max<long>(120L, grouter_->grid_->getTileSize() * 4L);
        if (wl_gain < min_wl_gain || via_gain < min_via_gain) {
          return false;
        }
        const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
        const long detour_guard = std::max<long>(15000L, tile_size * 36L);
        const long high_layer_guard = std::max<long>(
            tile_size * 44L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.14)));
        const long min_high_layer_guard = std::max<long>(
            tile_size * 28L,
            static_cast<long>(std::ceil(
                static_cast<double>(std::max(1L, wl_anchor->metrics.high_layer_dbu))
                * 0.18)));
        if (candidate->metrics.high_layer_dbu + min_high_layer_guard
            < wl_anchor->metrics.high_layer_dbu) {
          return false;
        }
        if (candidate->metrics.detour_dbu > wl_anchor->metrics.detour_dbu + detour_guard
            || candidate->metrics.high_layer_dbu
                   > wl_anchor->metrics.high_layer_dbu + high_layer_guard) {
          return false;
        }
        const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
        const double candidate_proxy
            = estimateDetailedRouteProxyCost(candidate->metrics);
        const bool proxy_ok
            = candidate_proxy + 1e-3 < anchor_proxy * 1.030
              || (wl_gain >= min_wl_gain * 2L && via_gain >= min_via_gain * 2L);
        return proxy_ok;
      };

      const ScenarioResult* strong_wl_upgrade = nullptr;
      for (const ScenarioResult* candidate : std::array<const ScenarioResult*, 9>{
               budgeted_lift_ptr,
               layer_floor_wl_ptr,
               layer_lift_wl_ptr,
               absolute_wl_ptr,
               absolute_wl_deep_ptr,
               min_wl_wide_ptr,
               cugr_sp_balance_ptr,
               radical_wl_fusion_ptr,
               hyper_absolute_wl_ptr}) {
        if (!candidate_is_viable_wl_upgrade(candidate)) {
          continue;
        }
        if (strong_wl_upgrade == nullptr
            || wirelength_first_better(*candidate, *strong_wl_upgrade)) {
          strong_wl_upgrade = candidate;
        }
      }

      if (strong_wl_upgrade != nullptr) {
        logger_->info(
            GNR,
            6035,
            "NEWGR WL-first override selecting '{}' over anchor '{}' "
            "(wl gain {}, via gain {}, detour delta {}, high-layer delta {}).",
            strong_wl_upgrade->name,
            wl_anchor->name,
            wl_anchor->metrics.wirelength_dbu
                - strong_wl_upgrade->metrics.wirelength_dbu,
            wl_anchor->metrics.via_count - strong_wl_upgrade->metrics.via_count,
            strong_wl_upgrade->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            strong_wl_upgrade->metrics.high_layer_dbu
                - wl_anchor->metrics.high_layer_dbu);
        forced_wl_ptr = strong_wl_upgrade;
      }
    }
    if (forced_wl_ptr == wl_anchor && wl_anchor != nullptr
        && via_floor_ptr != nullptr) {
      // Via-elastic WL promotion: allow a bounded via increase if it buys
      // meaningful WL reduction while keeping structural DR proxies stable.
      const long wl_gain
          = wl_anchor->metrics.wirelength_dbu - via_floor_ptr->metrics.wirelength_dbu;
      const long min_wl_gain = std::max<long>(
          45L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00018)));
      const long via_rise
          = static_cast<long>(via_floor_ptr->metrics.via_count)
            - static_cast<long>(wl_anchor->metrics.via_count);
      const long max_via_rise = std::max<long>(
          240L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0032)));
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long detour_guard = std::max<long>(tile_size * 14L, 7000L);
      const long high_layer_guard = std::max<long>(
          tile_size * 24L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.055)));
      const int hotspot_guard = std::max<int>(5, wl_anchor->metrics.near_capacity_edges / 6);
      const bool structural_guard
          = via_floor_ptr->metrics.detour_dbu
                 <= wl_anchor->metrics.detour_dbu + detour_guard
            && via_floor_ptr->metrics.high_layer_dbu
                   <= wl_anchor->metrics.high_layer_dbu + high_layer_guard
            && via_floor_ptr->metrics.near_capacity_edges
                   <= wl_anchor->metrics.near_capacity_edges + hotspot_guard;
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
      const double via_floor_proxy
          = estimateDetailedRouteProxyCost(via_floor_ptr->metrics);
      const bool proxy_guard = via_floor_proxy + 1e-3 < anchor_proxy * 1.012;
      if (wl_gain >= min_wl_gain && via_rise >= 0 && via_rise <= max_via_rise
          && structural_guard && proxy_guard) {
        logger_->info(
            GNR,
            6037,
            "NEWGR via-elastic WL promotion selecting '{}' over '{}' "
            "(wl gain {}, via rise {}, proxy ratio {:.3f}).",
            via_floor_ptr->name,
            wl_anchor->name,
            wl_gain,
            via_rise,
            anchor_proxy > 1e-9 ? via_floor_proxy / anchor_proxy : 1.0);
        forced_wl_ptr = via_floor_ptr;
      }
    }

    // Keep FastRoute-like short guides as the baseline final choice; only
    // upgrade to patched guides when patching stays near anchor WL and wins
    // the DR-proxy objective.
    if (forced_wl_ptr != nullptr && patched_ptr != nullptr) {
      const long patch_wl_guard = std::max<long>(
          80,
          static_cast<long>(std::ceil(
              static_cast<double>(forced_wl_ptr->metrics.wirelength_dbu)
              * 0.00055)));
      const bool patched_within_guard
          = patched_ptr->metrics.wirelength_dbu
            <= forced_wl_ptr->metrics.wirelength_dbu + patch_wl_guard;
      if (patched_within_guard) {
        const double anchor_proxy
            = estimateDetailedRouteProxyCost(forced_wl_ptr->metrics);
        const double patched_proxy
            = estimateDetailedRouteProxyCost(patched_ptr->metrics);
        if (patched_proxy + 1e-3 < anchor_proxy) {
          forced_wl_ptr = patched_ptr;
        }
      }
    }
    if (wl_anchor != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long min_wl_gain = std::max<long>(
          40L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00016)));
      const long via_rise_cap = std::max<long>(
          300L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0032)));
      const long via_drop_cap = std::max<long>(
          600L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0060)));
      const long detour_guard = std::max<long>(tile_size * 22L, 11000L);
      const long high_layer_guard = std::max<long>(
          tile_size * 36L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.11)));
      const int hotspot_guard = std::max<int>(8, wl_anchor->metrics.near_capacity_edges / 4);
      const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);

      const ScenarioResult* elastic_wl_ptr = nullptr;
      auto consider_elastic_wl = [&](const ScenarioResult* candidate) {
        if (candidate == nullptr || candidate == wl_anchor) {
          return;
        }
        if (candidate->metrics.overflow_edges > wl_anchor->metrics.overflow_edges) {
          return;
        }

        const long wl_gain
            = wl_anchor->metrics.wirelength_dbu - candidate->metrics.wirelength_dbu;
        if (wl_gain < min_wl_gain) {
          return;
        }

        const long via_delta
            = static_cast<long>(candidate->metrics.via_count)
              - static_cast<long>(wl_anchor->metrics.via_count);
        if (via_delta > via_rise_cap || via_delta < -via_drop_cap) {
          return;
        }

        const bool structural_guard
            = candidate->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
              && candidate->metrics.high_layer_dbu
                     <= wl_anchor->metrics.high_layer_dbu + high_layer_guard
              && candidate->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard;
        if (!structural_guard) {
          return;
        }

        const double candidate_proxy
            = estimateDetailedRouteProxyCost(candidate->metrics);
        // Via-drop candidates are allowed a looser proxy cap only when they
        // deliver larger WL gains; via-rise candidates stay stricter.
        const double proxy_cap
            = via_delta < 0 ? (wl_gain >= min_wl_gain * 2L ? 1.024 : 1.018) : 1.012;
        if (candidate_proxy + 1e-3 >= anchor_proxy * proxy_cap) {
          return;
        }

        if (elastic_wl_ptr == nullptr
            || wirelength_with_dr_proxy_tie_better(*candidate, *elastic_wl_ptr)) {
          elastic_wl_ptr = candidate;
        }
      };

      consider_elastic_wl(via_floor_ptr);
      consider_elastic_wl(length_adaptive_ptr);
      consider_elastic_wl(layer_floor_wl_ptr);
      consider_elastic_wl(drt_elastic_ptr);
      consider_elastic_wl(budgeted_lift_ptr);
      consider_elastic_wl(layer_lift_wl_ptr);
      consider_elastic_wl(wl_corridor_ptr);
      consider_elastic_wl(cugr_sp_balance_ptr);
      consider_elastic_wl(wl_layerbudget_ptr);
      consider_elastic_wl(radical_wl_fusion_ptr);
      consider_elastic_wl(hyper_absolute_wl_ptr);
      consider_elastic_wl(absolute_wl_deep_ptr);
      consider_elastic_wl(absolute_via_recover_ptr);

      if (elastic_wl_ptr != nullptr
          && (forced_wl_ptr == nullptr
              || wirelength_with_dr_proxy_tie_better(*elastic_wl_ptr, *forced_wl_ptr))) {
        bool preserve_layerbudget_mix = false;
        if (forced_wl_ptr == wl_layerbudget_ptr && wl_layerbudget_ptr != nullptr
            && elastic_wl_ptr != wl_layerbudget_ptr && wl_anchor != nullptr) {
          const long extra_vias
              = static_cast<long>(elastic_wl_ptr->metrics.via_count)
                - static_cast<long>(wl_layerbudget_ptr->metrics.via_count);
          const long extra_wl_gain
              = wl_layerbudget_ptr->metrics.wirelength_dbu
                - elastic_wl_ptr->metrics.wirelength_dbu;
          const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
          const long via_to_wl_exchange = std::max<long>(tile_size / 10L, 400L);
          preserve_layerbudget_mix
              = extra_vias > 0
                && extra_wl_gain <= extra_vias * via_to_wl_exchange
                && wl_layerbudget_ptr->metrics.high_layer_dbu
                       >= wl_anchor->metrics.high_layer_dbu;
        }

        if (!preserve_layerbudget_mix) {
          logger_->info(
              GNR,
              6036,
              "NEWGR elastic WL promotion selecting '{}' over anchor '{}' "
              "(wl gain {}, via delta {}, detour delta {}, high-layer delta {}).",
              elastic_wl_ptr->name,
              wl_anchor->name,
              wl_anchor->metrics.wirelength_dbu - elastic_wl_ptr->metrics.wirelength_dbu,
              elastic_wl_ptr->metrics.via_count - wl_anchor->metrics.via_count,
              elastic_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
              elastic_wl_ptr->metrics.high_layer_dbu
                  - wl_anchor->metrics.high_layer_dbu);
          forced_wl_ptr = elastic_wl_ptr;
        } else {
          logger_->info(
              GNR,
              7347,
              "NEWGR preserving layer-budget '{}' over elastic '{}' "
              "(extra vias {}, extra wl gain {}, exchange threshold {}).",
              wl_layerbudget_ptr->name,
              elastic_wl_ptr->name,
              static_cast<long>(elastic_wl_ptr->metrics.via_count)
                  - static_cast<long>(wl_layerbudget_ptr->metrics.via_count),
              wl_layerbudget_ptr->metrics.wirelength_dbu
                  - elastic_wl_ptr->metrics.wirelength_dbu,
              std::max<long>(std::max(grouter_->grid_->getTileSize(), 1), 400L));
        }
      }
    }

    // Compact overflow-free runs are highly deterministic on this benchmark.
    // Keep the final pick anchored to the stable hybrid-netmix-wl route unless
    // a challenger is a strict WL+via structural improvement.
    if (compact_exploration_mode && wl_anchor != nullptr) {
      if (forced_wl_ptr != nullptr && forced_wl_ptr != wl_anchor) {
        const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
        const long wl_gain
            = wl_anchor->metrics.wirelength_dbu - forced_wl_ptr->metrics.wirelength_dbu;
        const long via_gain
            = wl_anchor->metrics.via_count - forced_wl_ptr->metrics.via_count;
        const long min_wl_gain = std::max<long>(
            28L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.wirelength_dbu) * 0.00010)));
        const long min_via_gain = std::max<long>(
            80L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.00070)));
        const long detour_guard = std::max<long>(tile_size * 28L, 14000L);
        const long high_layer_guard = std::max<long>(
            tile_size * 40L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.16)));
        const long high_layer_floor = std::max<long>(
            tile_size * 16L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.90)));
        const long relaxed_high_layer_floor = std::max<long>(
            tile_size * 12L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.80)));
        const int hotspot_guard
            = std::max<int>(6, wl_anchor->metrics.near_capacity_edges / 5);
        const bool strict_wl_via_improvement
            = forced_wl_ptr->metrics.wirelength_dbu <= wl_anchor->metrics.wirelength_dbu
              && forced_wl_ptr->metrics.via_count <= wl_anchor->metrics.via_count
              && wl_gain >= min_wl_gain && via_gain >= min_via_gain;
        const bool structural_guard
            = forced_wl_ptr->metrics.detour_dbu
                   <= wl_anchor->metrics.detour_dbu + detour_guard
              && forced_wl_ptr->metrics.high_layer_dbu
                     <= wl_anchor->metrics.high_layer_dbu + high_layer_guard
              && forced_wl_ptr->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard;
        const bool layer_balance_guard
            = forced_wl_ptr->metrics.high_layer_dbu >= high_layer_floor;
        const double anchor_proxy = estimateDetailedRouteProxyCost(wl_anchor->metrics);
        const double challenger_proxy
            = estimateDetailedRouteProxyCost(forced_wl_ptr->metrics);
        const bool proxy_guard = challenger_proxy + 1e-3 < anchor_proxy * 1.015;
        const bool layer_floor_unlock
            = forced_wl_ptr == layer_floor_wl_ptr && wl_gain >= min_wl_gain * 2L
              && via_gain >= 0;
        const long elastic_via_rise_guard = std::max<long>(
            220L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.0022)));
        const bool drt_elastic_unlock
            = forced_wl_ptr == drt_elastic_ptr && wl_gain >= min_wl_gain
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + elastic_via_rise_guard;
        const long budgeted_via_rise_guard = std::max<long>(
            120L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.0012)));
        const bool budgeted_lift_unlock
            = forced_wl_ptr == budgeted_lift_ptr && wl_gain >= min_wl_gain
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + budgeted_via_rise_guard;
        const long layerbudget_via_rise_guard = std::max<long>(
            60L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.00060)));
        const bool layerbudget_unlock
            = forced_wl_ptr == wl_layerbudget_ptr && wl_gain >= min_wl_gain
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + layerbudget_via_rise_guard;
        const bool radical_wl_mix_candidate
            = forced_wl_ptr == absolute_wl_ptr
              || forced_wl_ptr == absolute_wl_deep_ptr
              || forced_wl_ptr == min_wl_wide_ptr
              || forced_wl_ptr == cugr_sp_balance_ptr
              || forced_wl_ptr == wl_layerbudget_ptr
              || forced_wl_ptr == radical_wl_fusion_ptr
              || forced_wl_ptr == hyper_absolute_wl_ptr;
        const long radical_high_layer_floor = std::max<long>(
            tile_size * 10L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.72)));
        const long radical_detour_guard = std::max<long>(detour_guard, tile_size * 34L);
        const bool radical_wl_unlock
            = radical_wl_mix_candidate && wl_gain >= min_wl_gain * 3L
              && via_gain >= min_via_gain
              && forced_wl_ptr->metrics.detour_dbu
                     <= wl_anchor->metrics.detour_dbu + radical_detour_guard
              && forced_wl_ptr->metrics.high_layer_dbu >= radical_high_layer_floor
              && challenger_proxy + 1e-3 < anchor_proxy * 0.995;
        const bool layer_guard_ok
            = layer_balance_guard
              || (forced_wl_ptr == budgeted_lift_ptr
                  && forced_wl_ptr->metrics.high_layer_dbu
                         >= relaxed_high_layer_floor);
        const bool conservative_keep
            = (strict_wl_via_improvement || layer_floor_unlock
               || drt_elastic_unlock || budgeted_lift_unlock
               || layerbudget_unlock)
              && structural_guard && layer_guard_ok && proxy_guard;
        const bool layerbudget_dominance_unlock
            = forced_wl_ptr == wl_layerbudget_ptr && wl_gain >= min_wl_gain
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + layerbudget_via_rise_guard
              && forced_wl_ptr->metrics.detour_dbu
                     <= wl_anchor->metrics.detour_dbu + detour_guard
              && forced_wl_ptr->metrics.high_layer_dbu >= high_layer_floor
              && forced_wl_ptr->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard
              && forced_wl_ptr->metrics.overflow_edges
                     <= wl_anchor->metrics.overflow_edges;
        const long via_recover_via_rise_guard = std::max<long>(
            280L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.0028)));
        const bool via_recover_unlock
            = forced_wl_ptr == absolute_via_recover_ptr
              && wl_gain >= min_wl_gain * 3L
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + via_recover_via_rise_guard
              && forced_wl_ptr->metrics.detour_dbu
                     <= wl_anchor->metrics.detour_dbu + detour_guard
              && forced_wl_ptr->metrics.high_layer_dbu >= high_layer_floor
              && forced_wl_ptr->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard
              && challenger_proxy + 1e-3 < anchor_proxy * 1.020;
        const long hyper_via_rise_guard = std::max<long>(
            320L,
            static_cast<long>(std::ceil(
                static_cast<double>(wl_anchor->metrics.via_count) * 0.0032)));
        const bool hyper_absolute_unlock
            = forced_wl_ptr == hyper_absolute_wl_ptr
              && wl_gain >= min_wl_gain * 2L
              && forced_wl_ptr->metrics.via_count
                     <= wl_anchor->metrics.via_count + hyper_via_rise_guard
              && forced_wl_ptr->metrics.detour_dbu
                     <= wl_anchor->metrics.detour_dbu + detour_guard
              && forced_wl_ptr->metrics.high_layer_dbu >= high_layer_floor
              && forced_wl_ptr->metrics.near_capacity_edges
                     <= wl_anchor->metrics.near_capacity_edges + hotspot_guard
              && challenger_proxy + 1e-3 < anchor_proxy * 1.015;
        const bool keep_challenger
            = conservative_keep || radical_wl_unlock
              || layerbudget_dominance_unlock || via_recover_unlock
              || hyper_absolute_unlock;
        if (keep_challenger) {
          logger_->info(
              GNR,
              6040,
              "NEWGR compact dominance unlock keeping '{}' over anchor '{}' "
              "(wl gain {}, via gain {}, detour delta {}, high-layer delta {}, "
              "high-layer floor {}, radical floor {}, proxy ratio {:.3f}, "
              "radical unlock {}, layerbudget unlock {}, hyper unlock {}).",
              forced_wl_ptr->name,
              wl_anchor->name,
              wl_gain,
              via_gain,
              forced_wl_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
              forced_wl_ptr->metrics.high_layer_dbu
                  - wl_anchor->metrics.high_layer_dbu,
              high_layer_floor,
              radical_high_layer_floor,
              anchor_proxy > 1e-9 ? challenger_proxy / anchor_proxy : 1.0,
              radical_wl_unlock ? "yes" : "no",
              layerbudget_dominance_unlock ? "yes" : "no",
              hyper_absolute_unlock ? "yes" : "no");
        } else {
          logger_->info(
              GNR,
              6038,
              "NEWGR compact anchor-lock keeping '{}' over '{}' "
              "(wl delta {}, via delta {}, detour delta {}, high-layer delta {}).",
              wl_anchor->name,
              forced_wl_ptr->name,
              wl_anchor->metrics.wirelength_dbu
                  - forced_wl_ptr->metrics.wirelength_dbu,
              wl_anchor->metrics.via_count - forced_wl_ptr->metrics.via_count,
              wl_anchor->metrics.detour_dbu - forced_wl_ptr->metrics.detour_dbu,
              wl_anchor->metrics.high_layer_dbu
                  - forced_wl_ptr->metrics.high_layer_dbu);
          forced_wl_ptr = wl_anchor;
        }
      }
      if (forced_wl_ptr == nullptr) {
        forced_wl_ptr = wl_anchor;
      }
    }

    // Final overflow-free WL lock:
    // choose the shortest-wirelength scenario that stays within a bounded
    // structural/via envelope around the wl anchor.
    if (wl_anchor != nullptr) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long detour_guard = std::max<long>(tile_size * 34L, 16000L);
      const long high_layer_guard = std::max<long>(
          tile_size * 52L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.24)));
      const long high_layer_floor = std::max<long>(
          tile_size * 16L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.high_layer_dbu) * 0.90)));
      const int hotspot_guard
          = std::max<int>(12, wl_anchor->metrics.near_capacity_edges / 3);
      const long via_rise_cap = std::max<long>(
          140L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0014)));
      const long via_drop_cap = std::max<long>(
          900L,
          static_cast<long>(std::ceil(
              static_cast<double>(wl_anchor->metrics.via_count) * 0.0090)));
      const long via_trade_wl_factor = std::max<long>(40L, tile_size * 2L);
      const long via_drop_trade_wl_factor = std::max<long>(65L, tile_size * 3L);
      auto is_extreme_wl_candidate = [](const std::string& name) {
        return name == "hybrid-netmix-absolute-wl"
               || name == "hybrid-netmix-absolute-wl-deep"
               || name == "hybrid-netmix-min-wl-wide"
               || name == "hybrid-netmix-cugr-sp-balance-wl"
               || name == "hybrid-netmix-radical-wl-fusion";
      };

      const ScenarioResult* wl_locked_ptr = forced_wl_ptr;
      for (const ScenarioResult& candidate : scenario_results) {
        if (candidate.metrics.overflow_edges > wl_anchor->metrics.overflow_edges) {
          continue;
        }
        const long wl_gain
            = wl_anchor->metrics.wirelength_dbu - candidate.metrics.wirelength_dbu;
        const long via_rise = static_cast<long>(candidate.metrics.via_count)
                              - static_cast<long>(wl_anchor->metrics.via_count);
        const long via_drop = static_cast<long>(wl_anchor->metrics.via_count)
                              - static_cast<long>(candidate.metrics.via_count);
        const bool minwl_extreme_candidate
            = is_extreme_wl_candidate(candidate.name);
        if (via_rise > via_rise_cap && wl_gain < via_rise * via_trade_wl_factor) {
          continue;
        }
        if (!minwl_extreme_candidate && via_drop > via_drop_cap
            && wl_gain < via_drop * via_drop_trade_wl_factor) {
          continue;
        }
        if (candidate.metrics.detour_dbu > wl_anchor->metrics.detour_dbu + detour_guard
            || candidate.metrics.high_layer_dbu
                   > wl_anchor->metrics.high_layer_dbu + high_layer_guard
            || candidate.metrics.near_capacity_edges
                   > wl_anchor->metrics.near_capacity_edges + hotspot_guard) {
          continue;
        }
        // Preserve enough upper-layer guide structure (FastRoute/CUGR-style
        // routability guard) so extreme low-layer collapses do not degrade
        // detailed routing despite lower global-route WL.
        if (candidate.metrics.high_layer_dbu < high_layer_floor
            && candidate.name != "hybrid-netmix-wl-layerbudget"
            && candidate.name != "hybrid-netmix-budgeted-lift-wl") {
          continue;
        }
        if (wl_locked_ptr == nullptr) {
          wl_locked_ptr = &candidate;
          continue;
        }

        bool better_locked = false;
        if (candidate.metrics.wirelength_dbu < wl_locked_ptr->metrics.wirelength_dbu) {
          const long wl_gain_over_locked
              = wl_locked_ptr->metrics.wirelength_dbu
                - candidate.metrics.wirelength_dbu;
          const long via_rise_over_locked
              = static_cast<long>(candidate.metrics.via_count)
                - static_cast<long>(wl_locked_ptr->metrics.via_count);
          const long via_rise_guard = std::max<long>(
              72L,
              static_cast<long>(std::ceil(
                  static_cast<double>(wl_locked_ptr->metrics.via_count) * 0.0007)));
          const bool candidate_extreme = is_extreme_wl_candidate(candidate.name);
          const bool locked_is_via_recover
              = wl_locked_ptr->name == "hybrid-netmix-absolute-via-recover";
          const long min_extreme_gain = std::max<long>(tile_size * 2L, 900L);
          if (candidate_extreme && wl_gain_over_locked >= min_extreme_gain) {
            better_locked = true;
          } else {
            const long via_rise_trade_wl_factor
                = candidate_extreme ? std::max<long>(18L, tile_size / 2L)
                                    : std::max<long>(420L, tile_size / 3L);
            if (via_rise_over_locked > via_rise_guard
                && !locked_is_via_recover
                && wl_gain_over_locked
                       < via_rise_over_locked * via_rise_trade_wl_factor) {
              better_locked = false;
            } else {
              better_locked = true;
            }
          }
        } else if (candidate.metrics.wirelength_dbu
                       == wl_locked_ptr->metrics.wirelength_dbu
                   && candidate.metrics.via_count < wl_locked_ptr->metrics.via_count) {
          better_locked = true;
        }

        if (better_locked) {
          wl_locked_ptr = &candidate;
        }
      }

      if (wl_locked_ptr != nullptr
          && (forced_wl_ptr == nullptr
              || wirelength_first_better(*wl_locked_ptr, *forced_wl_ptr))) {
        logger_->info(
            GNR,
            7356,
            "NEWGR overflow-free WL lock selecting '{}' over '{}' "
            "(wl gain {}, via delta {}, detour delta {}, high-layer delta {}).",
            wl_locked_ptr->name,
            forced_wl_ptr != nullptr ? forced_wl_ptr->name : std::string("none"),
            wl_anchor->metrics.wirelength_dbu - wl_locked_ptr->metrics.wirelength_dbu,
            wl_locked_ptr->metrics.via_count - wl_anchor->metrics.via_count,
            wl_locked_ptr->metrics.detour_dbu - wl_anchor->metrics.detour_dbu,
            wl_locked_ptr->metrics.high_layer_dbu
                - wl_anchor->metrics.high_layer_dbu);
        forced_wl_ptr = wl_locked_ptr;
      }
    }

    // SPRoute-style deterministic WL-plateau collapse with CUGR-style
    // routability preservation:
    // keep the shortest overflow-free WL plateau, but avoid low-layer
    // collapses by preferring lower detour and stronger upper-layer
    // guide elasticity before via minimization.
    int min_overflow_edges = std::numeric_limits<int>::max();
    long min_plateau_wl = std::numeric_limits<long>::max();
    for (const ScenarioResult& candidate : scenario_results) {
      if (candidate.metrics.overflow_edges < min_overflow_edges) {
        min_overflow_edges = candidate.metrics.overflow_edges;
        min_plateau_wl = candidate.metrics.wirelength_dbu;
      } else if (candidate.metrics.overflow_edges == min_overflow_edges) {
        min_plateau_wl
            = std::min(min_plateau_wl, candidate.metrics.wirelength_dbu);
      }
    }
    if (min_overflow_edges != std::numeric_limits<int>::max()
        && min_plateau_wl != std::numeric_limits<long>::max()) {
      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const long wl_plateau_band
          = std::max<long>(4L, static_cast<long>(std::ceil(
                                   static_cast<double>(min_plateau_wl) * 0.00001)));
      const long detour_rise_guard = std::max<long>(tile_size * 14L, 5600L);
      long via_rise_guard = std::numeric_limits<long>::max();
      long high_layer_drop_guard = std::numeric_limits<long>::max();
      if (forced_wl_ptr != nullptr) {
        via_rise_guard = std::max<long>(
            240L,
            static_cast<long>(std::ceil(
                static_cast<double>(forced_wl_ptr->metrics.via_count) * 0.0025)));
        high_layer_drop_guard = std::max<long>(
            tile_size * 44L,
            static_cast<long>(std::ceil(
                static_cast<double>(forced_wl_ptr->metrics.high_layer_dbu) * 0.10)));
      }

      const ScenarioResult* plateau_pick = forced_wl_ptr;
      auto plateau_better = [&](const ScenarioResult& lhs,
                                const ScenarioResult& rhs) {
        if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
          return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
        }
        if (lhs.metrics.detour_dbu != rhs.metrics.detour_dbu) {
          return lhs.metrics.detour_dbu < rhs.metrics.detour_dbu;
        }
        if (lhs.metrics.high_layer_dbu != rhs.metrics.high_layer_dbu) {
          return lhs.metrics.high_layer_dbu > rhs.metrics.high_layer_dbu;
        }
        const double lhs_proxy = estimateDetailedRouteProxyCost(lhs.metrics);
        const double rhs_proxy = estimateDetailedRouteProxyCost(rhs.metrics);
        if (std::abs(lhs_proxy - rhs_proxy) > 1e-3) {
          return lhs_proxy < rhs_proxy;
        }
        if (lhs.metrics.via_count != rhs.metrics.via_count) {
          return lhs.metrics.via_count < rhs.metrics.via_count;
        }
        return lhs.name < rhs.name;
      };

      for (const ScenarioResult& candidate : scenario_results) {
        if (candidate.metrics.overflow_edges != min_overflow_edges) {
          continue;
        }
        if (candidate.metrics.wirelength_dbu > min_plateau_wl + wl_plateau_band) {
          continue;
        }
        if (forced_wl_ptr != nullptr) {
          const long via_rise
              = static_cast<long>(candidate.metrics.via_count)
                - static_cast<long>(forced_wl_ptr->metrics.via_count);
          const long detour_rise
              = candidate.metrics.detour_dbu - forced_wl_ptr->metrics.detour_dbu;
          const long high_layer_drop
              = forced_wl_ptr->metrics.high_layer_dbu
                - candidate.metrics.high_layer_dbu;
          const bool strong_wl_gain = candidate.metrics.wirelength_dbu + wl_plateau_band / 2
                                      < forced_wl_ptr->metrics.wirelength_dbu;
          if (!strong_wl_gain && via_rise > via_rise_guard) {
            continue;
          }
          if (!strong_wl_gain && detour_rise > detour_rise_guard) {
            continue;
          }
          if (!strong_wl_gain && high_layer_drop > high_layer_drop_guard) {
            continue;
          }
        }
        if (plateau_pick == nullptr || plateau_better(candidate, *plateau_pick)) {
          plateau_pick = &candidate;
        }
      }

      if (plateau_pick != nullptr
          && (forced_wl_ptr == nullptr
              || plateau_better(*plateau_pick, *forced_wl_ptr))) {
        logger_->info(
            GNR,
            7364,
            "NEWGR WL-plateau lock selecting '{}' (min WL {}, band {}, "
            "detour {}, high-layer {}, vias {}) over '{}'.",
            plateau_pick->name,
            min_plateau_wl,
            wl_plateau_band,
            plateau_pick->metrics.detour_dbu,
            plateau_pick->metrics.high_layer_dbu,
            plateau_pick->metrics.via_count,
            forced_wl_ptr != nullptr ? forced_wl_ptr->name : std::string("none"));
        forced_wl_ptr = plateau_pick;
      }
    }
  }

  const ScenarioResult* best_ptr = forced_wl_ptr;
  if (best_ptr == nullptr) {
    best_ptr = *std::min_element(
        shortlist.begin(),
        shortlist.end(),
        [&](const ScenarioResult* lhs, const ScenarioResult* rhs) {
          if (overflow_free_sweep) {
            return wirelength_with_dr_proxy_tie_better(*lhs, *rhs);
          }
          return robust_better(*lhs, *rhs);
        });
  }
  auto best_iter = scenario_results.begin()
                   + static_cast<std::ptrdiff_t>(best_ptr
                                                 - &scenario_results.front());
  ScenarioResult final_result = *best_iter;

  logger_->info(GNR,
                6007,
                "NEWGR best scenario '{}': wirelength {:.0f} um, vias {}, "
                "overflow edges {}, hot edges {}, max ratio {:.2f}, detour {}, "
                "high-layer {}",
                final_result.name,
                final_result.metrics.wirelength_um,
                final_result.metrics.via_count,
                final_result.metrics.overflow_edges,
                final_result.metrics.near_capacity_edges,
                final_result.metrics.max_usage_ratio,
                final_result.metrics.detour_dbu,
                final_result.metrics.high_layer_dbu);

  return std::move(final_result.routes);
}

}  // namespace grt
