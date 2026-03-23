#include "NEWGR/NewGR.h"

#include <algorithm>
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

bool isPreferredWirelengthScenario(const std::string& name)
{
  // Favor DR-stable wirelength hybrids over ultra-min guide-WL variants.
  return name == "hybrid-netmix-wl" || name == "hybrid-netmix-wl-safe"
         || name == "hybrid-netmix-dr-stable";
}

bool isAggressiveWirelengthScenario(const std::string& name)
{
  return name == "hybrid-netmix-hpwl-lock"
         || name == "hybrid-netmix-absolute-wl"
         || name == "hybrid-netmix-min-wl-extreme"
         || name == "hybrid-netmix-min-wl-wide"
         || name == "hybrid-netmix-min-wl";
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
      if (refinement_defs.size() >= 12) {
        return;
      }
      if (scenario_names.insert(def.name).second) {
        refinement_defs.push_back(std::move(def));
      }
    };

    const int elite_count = std::min<int>(4, ranked.size());
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
    const int layer_compact_source_count = std::min<int>(20, ranked.size());
    const int balanced_source_count = std::min<int>(8, ranked.size());
    const int softcap_source_count = std::min<int>(16, ranked.size());
    const int consensus_source_count = std::min<int>(22, ranked.size());
    const int detour_source_count = std::min<int>(24, ranked.size());
    const int pareto_softcap_source_count = std::min<int>(28, ranked.size());
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
    const long min_wl_via_weight = 0L;
    append_hybrid("hybrid-netmix-wl", wl_source_count, 0, 0.24, 0.50, 6010);
    append_hybrid(
        "hybrid-netmix-ultra-wl", ultra_wl_source_count, 0, 0.40, 0.95, 6013);
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
  const long tie_preferred_wl_band
      = std::max<long>(150, static_cast<long>(std::ceil(shortest_wl * 0.0034)));
  const long tie_proxy_wl_band
      = std::max<long>(220, static_cast<long>(std::ceil(shortest_wl * 0.0052)));
  auto wirelength_with_dr_proxy_tie_better = [&](const ScenarioResult& lhs,
                                                  const ScenarioResult& rhs) {
    const long wl_gap = std::llabs(lhs.metrics.wirelength_dbu
                                   - rhs.metrics.wirelength_dbu);
    const bool lhs_pref = isPreferredWirelengthScenario(lhs.name);
    const bool rhs_pref = isPreferredWirelengthScenario(rhs.name);
    if (lhs_pref != rhs_pref && wl_gap <= tie_preferred_wl_band) {
      return lhs_pref;
    }
    if (wl_gap <= tie_proxy_wl_band) {
      const bool lhs_aggressive = isAggressiveWirelengthScenario(lhs.name);
      const bool rhs_aggressive = isAggressiveWirelengthScenario(rhs.name);
      if (lhs_aggressive != rhs_aggressive) {
        return !lhs_aggressive;
      }
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
    for (const ScenarioResult& result : scenario_results) {
      if (result.name == "hybrid-netmix-wl") {
        forced_wl_ptr = &result;
        break;
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

  const ScenarioDefinition* replay_def = nullptr;
  if (best_iter->name != scenario_results.back().name) {
    if (best_iter->name == "baseline") {
      replay_def = &baseline_def;
    } else {
      for (const ScenarioDefinition& def : scenario_defs) {
        if (def.name == best_iter->name) {
          replay_def = &def;
          break;
        }
      }
    }
    if (replay_def != nullptr) {
      final_result = run_scenario(*replay_def, snapshot);
    }
  }

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
