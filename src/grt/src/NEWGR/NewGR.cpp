#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "grt/Rudy.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouteMetrics
{
  long wirelength_dbu = 0;
  long via_count = 0;
  long overflow = 0;
  double max_utilization = 0.0;
  double avg_utilization = 0.0;
  double stress_cost = 0.0;
  double reserve_score = 0.0;
  double wirelength_um = 0.0;
  double score = 0.0;
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
struct RudyStats
{
  float mean = 0.0f;
  float p50 = 0.0f;
  float p80 = 0.0f;
  float p90 = 0.0f;
  float max = 0.0f;
};

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
  std::vector<float> values;
  values.reserve(static_cast<size_t>(x_tiles * y_tiles));
  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      const float val = rudy->getTile(x, y).getRudy();
      normalized[x][y] = val;
      max_value = std::max(max_value, val);
      values.push_back(val);
    }
  }

  if (values.empty() || max_value <= std::numeric_limits<float>::epsilon()) {
    return normalized;
  }

  std::vector<float> sorted_values = values;
  std::sort(sorted_values.begin(), sorted_values.end());
  const size_t p95_idx = static_cast<size_t>(std::round(
      0.95 * static_cast<double>(std::max<size_t>(sorted_values.size() - 1, 0))));
  const float p95 = sorted_values[p95_idx];
  const float scale
      = std::max({p95, max_value * 0.6f, std::numeric_limits<float>::epsilon()});

  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      normalized[x][y] = std::clamp(normalized[x][y] / scale, 0.0f, 1.0f);
    }
  }

  RudyGrid smoothed = normalized;
  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      float accum = 0.0f;
      float weight_sum = 0.0f;
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          const int nx = x + dx;
          const int ny = y + dy;
          if (nx < 0 || ny < 0 || nx >= x_tiles || ny >= y_tiles) {
            continue;
          }
          const bool is_center = dx == 0 && dy == 0;
          const bool is_axis = (dx == 0) != (dy == 0);
          const float weight = is_center ? 4.0f : (is_axis ? 2.0f : 1.0f);
          accum += normalized[nx][ny] * weight;
          weight_sum += weight;
        }
      }
      if (weight_sum > std::numeric_limits<float>::epsilon()) {
        smoothed[x][y] = std::clamp(accum / weight_sum, 0.0f, 1.0f);
      }
    }
  }

  return smoothed;
}

RudyStats computeRudyStats(const RudyGrid& normalized_rudy)
{
  RudyStats stats;
  if (normalized_rudy.empty()) {
    return stats;
  }

  std::vector<float> values;
  values.reserve(normalized_rudy.size() * normalized_rudy.front().size());
  for (const auto& column : normalized_rudy) {
    for (float value : column) {
      stats.max = std::max(stats.max, value);
      stats.mean += value;
      values.push_back(value);
    }
  }

  if (values.empty()) {
    return stats;
  }

  const float inv_size = 1.0f / static_cast<float>(values.size());
  stats.mean *= inv_size;

  std::sort(values.begin(), values.end());
  auto percentile = [&](float pct) {
    if (values.empty()) {
      return 0.0f;
    }
    pct = std::clamp(pct, 0.0f, 1.0f);
    const size_t idx = std::min(
        static_cast<size_t>(pct * static_cast<float>(values.size() - 1)),
        values.size() - 1);
    return std::clamp(values[idx], 0.0f, 1.0f);
  };

  stats.p50 = percentile(0.50f);
  stats.p80 = percentile(0.80f);
  stats.p90 = percentile(0.90f);
  return stats;
}

float summarizeHotspotScore(const std::vector<Hotspot>& hotspots)
{
  if (hotspots.empty()) {
    return 0.0f;
  }

  float max_severity = 0.0f;
  float accum_severity = 0.0f;
  for (const Hotspot& hotspot : hotspots) {
    max_severity = std::max(max_severity, hotspot.severity);
    accum_severity += hotspot.severity;
  }

  const float avg_severity
      = accum_severity / std::max(static_cast<int>(hotspots.size()), 1);
  const float density_term = 0.12f
                             * static_cast<float>(
                                 std::log1p(static_cast<double>(hotspots.size())));
  const float score
      = 0.35f * max_severity + 0.25f * avg_severity + density_term;
  return std::clamp(score, 0.0f, 1.5f);
}

float computeCongestionSeverity(const RudyStats& stats,
                                const std::vector<Hotspot>& hotspots)
{
  const float hotspot_score = summarizeHotspotScore(hotspots);
  const float combined = 0.55f * stats.p80 + 0.20f * stats.mean
                         + 0.15f * stats.p90 + 0.10f * hotspot_score;
  return std::clamp(combined, 0.0f, 1.0f);
}

void adjustEdgeCapacity(GlobalRouter* grouter,
                        int x1,
                        int y1,
                        int x2,
                        int y2,
                        int layer,
                        float ratio)
{
  // Allow limited capacity boosts to open low-congestion corridors.
  const float max_ratio_cap = 1.12f;
  ratio = std::clamp(ratio, 0.05f, max_ratio_cap);
  FastRouteCore* core = grouter->fastroute();
  if (core == nullptr) {
    return;
  }
  const int current_cap = core->getEdgeCapacity(x1, y1, x2, y2, layer);
  if (current_cap <= 0) {
    return;
  }
  int new_cap
      = std::max(1, static_cast<int>(std::floor(current_cap * ratio)));
  if (ratio > 1.01f && new_cap == current_cap) {
    const int boosted_cap = static_cast<int>(std::ceil(current_cap * max_ratio_cap));
    new_cap = std::min(current_cap + 1, boosted_cap);
  }
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
                              float midpoint = 0.45f)
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
                                 float max_ratio) {
    normalized = std::clamp(normalized, 0.0f, 1.0f);
    const float exponent = slope * (normalized - midpoint);
    const float logistic = 1.0f / (1.0f + std::exp(exponent));
    const float blend = min_ratio + (max_ratio - min_ratio) * logistic;
    return std::clamp(blend, 0.05f, 0.99f);
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
            = logistic_ratio(normalized, slope, midpoint, min_ratio, max_ratio);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        const float ratio
            = logistic_ratio(normalized, slope, midpoint, min_ratio, max_ratio);
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

void applySelectiveRelief(GlobalRouter* grouter,
                          const RudyGrid& normalized_rudy,
                          const std::vector<Hotspot>& hotspots,
                          int min_layer,
                          int max_layer,
                          float threshold,
                          float min_ratio,
                          float hotspot_push,
                          int halo,
                          float cool_threshold = 0.0f,
                          float boost = 0.0f,
                          float boost_limit = 1.0f,
                          float layer_boost_falloff = 0.0f)
{
  Grid* grid = grouter->grid();
  if (grid == nullptr || (normalized_rudy.empty() && hotspots.empty())) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return;
  }

  threshold = std::clamp(threshold, 0.0f, 1.0f);
  min_ratio = std::clamp(min_ratio, 0.70f, 0.99f);
  hotspot_push = std::clamp(hotspot_push, 0.0f, 1.0f);
  halo = std::max(0, halo);

  const int x_tiles
      = normalized_rudy.empty() ? 0 : static_cast<int>(normalized_rudy.size());
  const int y_tiles = (normalized_rudy.empty() || normalized_rudy.front().empty())
                          ? 0
                          : static_cast<int>(normalized_rudy.front().size());

  const auto getNormalized = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= x_tiles || y >= y_tiles || x_tiles == 0
        || y_tiles == 0) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  std::vector<std::vector<float>> horiz_influence(
      std::max(x_grids, 1), std::vector<float>(std::max(y_grids, 1), 0.0f));
  std::vector<std::vector<float>> vert_influence(
      std::max(x_grids, 1), std::vector<float>(std::max(y_grids, 1), 0.0f));

  if (!hotspots.empty()) {
    for (const Hotspot& hotspot : hotspots) {
      for (int dx = -halo; dx <= halo; ++dx) {
        for (int dy = -halo; dy <= halo; ++dy) {
          const int gx = hotspot.gx + dx;
          const int gy = hotspot.gy + dy;
          if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
            continue;
          }
          const float distance = static_cast<float>(std::abs(dx) + std::abs(dy));
          const float decay = 1.0f / (1.0f + distance);
          const float influence = hotspot.severity * decay;
          if (hotspot.affect_horizontal && gx < x_grids - 1) {
            horiz_influence[gx][gy]
                = std::max(horiz_influence[gx][gy], influence);
          }
          if (hotspot.affect_vertical && gy < y_grids - 1) {
            vert_influence[gx][gy]
                = std::max(vert_influence[gx][gy], influence);
          }
        }
      }
    }
  }

  const int layer_span = std::max(max_layer - min_layer, 1);
  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float layer_min_ratio
        = std::clamp(min_ratio - 0.04f * layer_factor, 0.70f, 0.99f);

    for (int y = 0; y < y_grids; ++y) {
      for (int x = 0; x < x_grids - 1; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        float ratio = 1.0f;
        if (normalized > threshold) {
          const float t = (normalized - threshold)
                          / std::max(0.001f, 1.0f - threshold);
          ratio -= t * (1.0f - layer_min_ratio);
        }
        ratio = std::clamp(ratio, layer_min_ratio, 1.0f);

        const float influence = horiz_influence[x][y];
        if (influence > 0.0f && hotspot_push > 0.0f) {
          float dampen
              = 1.0f - hotspot_push * (influence - 1.0f) * 0.35f;
          dampen = std::clamp(dampen, 0.70f, 1.0f);
          ratio = std::clamp(ratio * dampen, layer_min_ratio * 0.92f, 1.0f);
        }

        if (boost > 0.0f && boost_limit > 1.0f && cool_threshold > 0.0f) {
          const float cool_score = std::clamp(
              (cool_threshold - normalized) / std::max(cool_threshold, 1e-3f),
              0.0f,
              1.0f);
          if (cool_score > 0.0f) {
            float layer_scale
                = 1.0f - layer_boost_falloff * std::clamp(layer_factor, 0.0f, 1.0f);
            layer_scale = std::clamp(layer_scale, 0.5f, 1.0f);
            float influence_guard
                = 1.0f - std::clamp(influence, 0.0f, 1.0f) * 0.65f;
            influence_guard = std::clamp(influence_guard, 0.35f, 1.0f);
            const float candidate
                = 1.0f + boost * cool_score * layer_scale * influence_guard;
            ratio = std::max(ratio, std::min(candidate, boost_limit));
          }
        }

        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < y_grids - 1; ++y) {
      for (int x = 0; x < x_grids; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        float ratio = 1.0f;
        if (normalized > threshold) {
          const float t = (normalized - threshold)
                          / std::max(0.001f, 1.0f - threshold);
          ratio -= t * (1.0f - layer_min_ratio);
        }
        ratio = std::clamp(ratio, layer_min_ratio, 1.0f);

        const float influence = vert_influence[x][y];
        if (influence > 0.0f && hotspot_push > 0.0f) {
          float dampen
              = 1.0f - hotspot_push * (influence - 1.0f) * 0.35f;
          dampen = std::clamp(dampen, 0.70f, 1.0f);
          ratio = std::clamp(ratio * dampen, layer_min_ratio * 0.92f, 1.0f);
        }

        if (boost > 0.0f && boost_limit > 1.0f && cool_threshold > 0.0f) {
          const float cool_score = std::clamp(
              (cool_threshold - normalized) / std::max(cool_threshold, 1e-3f),
              0.0f,
              1.0f);
          if (cool_score > 0.0f) {
            float layer_scale
                = 1.0f - layer_boost_falloff * std::clamp(layer_factor, 0.0f, 1.0f);
            layer_scale = std::clamp(layer_scale, 0.5f, 1.0f);
            float influence_guard
                = 1.0f - std::clamp(influence, 0.0f, 1.0f) * 0.65f;
            influence_guard = std::clamp(influence_guard, 0.35f, 1.0f);
            const float candidate
                = 1.0f + boost * cool_score * layer_scale * influence_guard;
            ratio = std::max(ratio, std::min(candidate, boost_limit));
          }
        }

        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
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
    const int tile_size
        = grouter_ != nullptr && grouter_->grid_ != nullptr
              ? std::max(grouter_->grid_->getTileSize(), 1)
              : 1;

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

    if (grouter_ != nullptr && grouter_->fastroute_ != nullptr) {
      FastRouteCore* core = grouter_->fastroute_;
      metrics.overflow = core->totalOverflow();

      const auto& usage = core->getTotalUsagePerLayer();
      const auto& capacity = core->getTotalCapacityPerLayer();
      const size_t layer_count = std::min(usage.size(), capacity.size());
      double stress = 0.0;
      double reserve = 0.0;
      for (size_t i = 0; i < layer_count; ++i) {
        const double cap = static_cast<double>(capacity[i]);
        if (cap <= std::numeric_limits<double>::epsilon()) {
          continue;
        }
        const double util = static_cast<double>(usage[i]) / cap;
        metrics.max_utilization = std::max(metrics.max_utilization, util);
        metrics.avg_utilization += util;
        const double overload = std::max(0.0, util - 0.86);
        stress += overload * overload;
        reserve += std::max(0.0, 1.0 - util);
      }

      if (layer_count > 0) {
        metrics.avg_utilization
            /= std::max<double>(static_cast<double>(layer_count), 1.0);
        metrics.stress_cost = stress * static_cast<double>(tile_size);
        metrics.reserve_score = reserve;
      }
    }

    if (metrics.wirelength_dbu > 0 && grouter_->db_ != nullptr
        && grouter_->db_->getTech() != nullptr) {
      metrics.wirelength_um
          = metrics.wirelength_dbu
            / static_cast<double>(
                grouter_->db_->getTech()->getDbUnitsPerMicron());
    }

    const double via_weight = static_cast<double>(tile_size) * 3.0;
    const double overflow_weight = static_cast<double>(tile_size) * 12.0;
    metrics.score = static_cast<double>(metrics.wirelength_dbu)
                    + via_weight * static_cast<double>(metrics.via_count)
                    + overflow_weight * static_cast<double>(metrics.overflow);
    return metrics;
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
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow {}, "
                  "max util {:.2f}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow,
                  metrics.max_utilization);
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
    logger_->info(GNR,
                  6006,
                  "NEWGR scenario {}: wirelength {:.0f} um, vias {}, overflow {}, "
                  "max util {:.2f}",
                  scenario.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow,
                  metrics.max_utilization);
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

  const RudyStats rudy_stats = computeRudyStats(normalized_rudy);
  const float congestion_severity
      = computeCongestionSeverity(rudy_stats, hotspots);
  const bool force_routability = baseline.metrics.overflow > 0;
  const bool has_congestion_data
      = (!normalized_rudy.empty() || !hotspots.empty());
  const float hotspot_bias = std::clamp(
      static_cast<float>(hotspots.size()) / 12.0f, 0.0f, 0.6f);
  const bool run_soft
      = force_routability
        || (has_congestion_data
            && (congestion_severity > 0.72f || hotspots.size() > 2
                || rudy_stats.p80 > 0.84f))
        || hotspots.size() > 5;
  const bool run_aggressive_soft
      = force_routability
        || (has_congestion_data
            && (congestion_severity > 0.88f || hotspots.size() > 4));
  const bool allow_seed_sweep = !force_routability && hotspots.size() <= 4;

  logger_->info(GNR,
                6008,
                "NEWGR congestion severity {:.2f} (RUDY mean {:.2f}, p80 "
                "{:.2f}, hotspots {}, baseline overflow {})",
                congestion_severity,
                rudy_stats.mean,
                rudy_stats.p80,
                hotspots.size(),
                baseline.metrics.overflow);

  auto make_soft_config
      = [&](const std::string& name,
            float min_base,
            float max_base,
            float slope,
            float midpoint,
            int halo,
            float hotspot_ratio,
            float severity_weight,
            float perturb_pct,
            int seed,
            float critical_pct) {
          ScenarioDefinition def;
          def.name = name;
          def.pre_init = [this, perturb_pct, seed, critical_pct]() {
            grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
            grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(seed);
            grouter_->setAllowCongestion(false);
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
                                           midpoint);
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

  const float wl_critical_pct
      = std::clamp(5.5f + 3.5f * (0.6f - congestion_severity), 4.5f, 11.0f);
  const float wl_perturb_pct
      = congestion_severity > 0.45f ? 0.22f : 0.0f;
  const int wl_seed = snapshot.seed + 5;
  ScenarioDefinition wl_variation;
  wl_variation.name = "wl-variation";
  wl_variation.pre_init
      = [this, wl_perturb_pct, wl_seed, wl_critical_pct]() {
          grouter_->setCapacitiesPerturbationPercentage(wl_perturb_pct);
          grouter_->setPerturbationAmount(wl_perturb_pct > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_critical_pct);
        };
  wl_variation.post_init
      = [this, &hotspots, min_routing_layer, max_routing_layer]() {
          if (!hotspots.empty()) {
            applyHotspotPenalties(grouter_,
                                  hotspots,
                                  min_routing_layer,
                                  max_routing_layer,
                                  1,
                                  0.985f,
                                  0.12f);
          }
        };
  scenario_defs.push_back(wl_variation);

  const float wl_greedy_perturb = congestion_severity > 0.30f ? 0.12f : 0.0f;
  const float wl_greedy_critical
      = std::clamp(wl_critical_pct - 0.8f, 3.5f, 10.0f);
  const int wl_greedy_seed = snapshot.seed + 37;
  ScenarioDefinition wl_greedy;
  wl_greedy.name = "wl-greedy";
  wl_greedy.pre_init
      = [this, wl_greedy_perturb, wl_greedy_seed, wl_greedy_critical]() {
          grouter_->setCapacitiesPerturbationPercentage(wl_greedy_perturb);
          grouter_->setPerturbationAmount(wl_greedy_perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_greedy_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_greedy_critical);
        };
  wl_greedy.post_init = []() {};
  scenario_defs.push_back(wl_greedy);

  const float wl_refine_perturb = wl_greedy_perturb * 0.6f;
  const float wl_refine_critical
      = std::clamp(wl_greedy_critical - 0.6f, 3.0f, 10.0f);
  const int wl_refine_seed = snapshot.seed + 71;
  ScenarioDefinition wl_refine;
  wl_refine.name = "wl-refine";
  wl_refine.pre_init
      = [this, wl_refine_perturb, wl_refine_seed, wl_refine_critical]() {
          const float perturb = std::clamp(wl_refine_perturb, 0.0f, 0.40f);
          grouter_->setCapacitiesPerturbationPercentage(perturb);
          grouter_->setPerturbationAmount(perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_refine_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_refine_critical);
        };
  wl_refine.post_init
      = [this, &hotspots, min_routing_layer, max_routing_layer]() {
          if (!hotspots.empty()) {
            applyHotspotPenalties(grouter_,
                                  hotspots,
                                  min_routing_layer,
                                  max_routing_layer,
                                  1,
                                  0.992f,
                                  0.10f);
          }
        };
  scenario_defs.push_back(wl_refine);

  auto make_wl_greedy_seed
      = [this, wl_greedy_perturb, wl_greedy_critical, &snapshot](
            const std::string& name, int seed_offset, float perturb_scale) {
          ScenarioDefinition def;
          def.name = name;
          const float perturb
              = std::clamp(wl_greedy_perturb * perturb_scale, 0.0f, 0.60f);
          const int seed = snapshot.seed + seed_offset;
          def.pre_init
              = [this, perturb, seed, wl_greedy_critical]() {
                  grouter_->setCapacitiesPerturbationPercentage(perturb);
                  grouter_->setPerturbationAmount(perturb > 0.0f ? 1 : 0);
                  grouter_->setSeed(seed);
                  grouter_->setAllowCongestion(false);
                  grouter_->fastroute_->setCriticalNetsPercentage(
                      wl_greedy_critical);
                };
          def.post_init = []() {};
          return def;
        };

  if (allow_seed_sweep) {
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s1", 73, 0.8f));
    if (hotspots.size() <= 2) {
      scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s2", 109, 1.1f));
    }
  }

  if (has_congestion_data && run_soft && !normalized_rudy.empty()) {
    const float contour_strength
        = std::clamp(congestion_severity * 0.70f + hotspot_bias * 0.40f,
                     0.0f,
                     1.0f);
    const float contour_min_base
        = std::clamp(0.88f - 0.06f * contour_strength, 0.80f, 0.90f);
    const float contour_max_base
        = std::clamp(0.97f - 0.03f * contour_strength,
                     contour_min_base + 0.02f,
                     0.985f);
    const float contour_slope = 2.6f + 1.0f * contour_strength;
    const float contour_midpoint
        = std::clamp(0.54f - 0.05f * contour_strength, 0.46f, 0.54f);
    const float contour_threshold
        = std::clamp(0.36f + 0.18f * (1.0f - congestion_severity),
                     0.32f,
                     0.60f);
    const float contour_min_ratio
        = std::clamp(0.84f + 0.08f * (1.0f - congestion_severity),
                     0.84f,
                     0.94f);
    const float contour_push
        = std::clamp(0.10f + 0.16f * contour_strength, 0.10f, 0.22f);
    const float contour_hotspot_ratio
        = std::clamp(0.98f - 0.08f * contour_strength, 0.86f, 0.99f);
    const float contour_severity_weight
        = std::clamp(0.12f + 0.22f * contour_strength, 0.12f, 0.38f);
    const int contour_halo = contour_strength > 0.55f ? 2 : 1;
    const float contour_perturb
        = std::clamp(0.06f + 0.50f * contour_strength, 0.04f, 0.55f);
    const float contour_critical
        = std::clamp(5.5f + 3.5f * (0.65f - congestion_severity),
                     4.5f,
                     10.0f);
    const int contour_seed = snapshot.seed + 13;

    ScenarioDefinition contour_soft;
    contour_soft.name = "contour-lite";
    contour_soft.pre_init
        = [this, contour_perturb, contour_seed, contour_critical]() {
            grouter_->setCapacitiesPerturbationPercentage(contour_perturb);
            grouter_->setPerturbationAmount(contour_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(contour_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(contour_critical);
          };
    contour_soft.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           contour_min_base,
           contour_max_base,
           contour_slope,
           contour_midpoint,
           contour_threshold,
           contour_min_ratio,
           contour_push,
           contour_hotspot_ratio,
           contour_severity_weight,
           contour_halo]() {
            applySoftCapacityScaling(grouter_,
                                     normalized_rudy,
                                     min_routing_layer,
                                     max_routing_layer,
                                     contour_min_base,
                                     contour_max_base,
                                     contour_slope,
                                     contour_midpoint);
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 contour_threshold,
                                 contour_min_ratio,
                                 contour_push,
                                 contour_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    contour_halo,
                                    contour_hotspot_ratio,
                                    contour_severity_weight);
            }
          };
    scenario_defs.push_back(contour_soft);
  }

  if (has_congestion_data) {
    const float corridor_threshold
        = std::clamp(0.34f + 0.14f * (1.0f - congestion_severity),
                     0.30f,
                     0.50f);
    const float corridor_min_ratio
        = std::clamp(0.88f + 0.05f * (1.0f - congestion_severity),
                     0.88f,
                     0.95f);
    const float corridor_hotspot_push
        = std::clamp(0.10f + 0.10f * congestion_severity + 0.04f * hotspot_bias,
                     0.10f,
                     0.22f);
    const int corridor_halo = congestion_severity > 0.62f ? 2 : 1;
    const float corridor_cool_threshold
        = std::clamp(corridor_threshold * 0.70f, 0.22f, 0.42f);
    const float corridor_boost
        = std::clamp(0.07f + 0.05f * (1.0f - congestion_severity),
                     0.06f,
                     0.13f);
    const float corridor_boost_limit
        = std::clamp(1.07f + 0.03f * (1.0f - congestion_severity),
                     1.07f,
                     1.11f);
    const float corridor_layer_falloff
        = std::clamp(0.16f + 0.08f * hotspot_bias, 0.12f, 0.30f);
    const float corridor_perturb
        = std::clamp(0.05f + 0.28f * congestion_severity, 0.05f, 0.32f);
    const float corridor_critical
        = std::clamp(5.0f + 2.5f * (0.55f - congestion_severity),
                     4.0f,
                     8.5f);
    const int corridor_seed = snapshot.seed + 131;

    ScenarioDefinition corridor;
    corridor.name = "cool-corridors";
    corridor.pre_init
        = [this, corridor_perturb, corridor_seed, corridor_critical]() {
            grouter_->setCapacitiesPerturbationPercentage(corridor_perturb);
            grouter_->setPerturbationAmount(corridor_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(corridor_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(corridor_critical);
          };
    corridor.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           corridor_threshold,
           corridor_min_ratio,
           corridor_hotspot_push,
           corridor_halo,
           corridor_cool_threshold,
           corridor_boost,
           corridor_boost_limit,
           corridor_layer_falloff]() {
           applySelectiveRelief(grouter_,
                                normalized_rudy,
                                hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                corridor_threshold,
                                 corridor_min_ratio,
                                 corridor_hotspot_push,
                                 corridor_halo,
                                 corridor_cool_threshold,
                                 corridor_boost,
                                 corridor_boost_limit,
                                 corridor_layer_falloff);
          };
    scenario_defs.push_back(corridor);
  }

  if (has_congestion_data) {
    const float boost_strength
        = std::clamp(0.08f + 0.10f * (0.65f - congestion_severity),
                     0.06f,
                     0.16f);
    const float boost_limit
        = std::clamp(1.08f + 0.08f * (0.50f - congestion_severity),
                     1.06f,
                     1.12f);
    const float threshold
        = std::clamp(0.44f + 0.18f * congestion_severity, 0.44f, 0.65f);
    const float min_ratio
        = std::clamp(0.90f - 0.06f * congestion_severity, 0.82f, 0.92f);
    const float hotspot_push
        = std::clamp(0.08f + 0.16f * congestion_severity + 0.05f * hotspot_bias,
                     0.08f,
                     0.26f);
    const int halo
        = (congestion_severity > 0.55f || hotspots.size() > 3) ? 2 : 1;
    const float cool_threshold
        = std::clamp(0.28f + 0.16f * (1.0f - congestion_severity),
                     0.24f,
                     0.44f);
    const float layer_falloff
        = std::clamp(0.14f + 0.10f * hotspot_bias, 0.12f, 0.28f);
    const float perturb_pct
        = std::clamp(0.03f + 0.30f * (0.60f - congestion_severity),
                     0.02f,
                     0.22f);
    const float critical_pct
        = std::clamp(4.2f + 3.5f * (0.70f - congestion_severity),
                     3.8f,
                     9.5f);
    const int seed = snapshot.seed + 191;

    ScenarioDefinition wl_coolboost;
    wl_coolboost.name = "wl-coolboost";
    wl_coolboost.pre_init
        = [this, perturb_pct, seed, critical_pct]() {
            grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
            grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
          };
    wl_coolboost.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           threshold,
           min_ratio,
           hotspot_push,
           halo,
           cool_threshold,
           boost_strength,
           boost_limit,
           layer_falloff]() {
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 threshold,
                                 min_ratio,
                                 hotspot_push,
                                 halo,
                                 cool_threshold,
                                 boost_strength,
                                 boost_limit,
                                 layer_falloff);
          };
    scenario_defs.push_back(wl_coolboost);
  }

  if (has_congestion_data && run_aggressive_soft) {
    const float selective_threshold = std::clamp(
        0.32f + 0.22f * (1.0f - congestion_severity), 0.26f, 0.58f);
    const float selective_min_ratio
        = std::clamp(0.78f + 0.10f * (1.0f - congestion_severity),
                     0.78f,
                     0.90f);
    const float selective_hotspot_push
        = std::clamp(0.10f + 0.22f * congestion_severity + 0.04f * hotspot_bias,
                     0.10f,
                     0.32f);
    const int selective_halo = congestion_severity > 0.55f ? 2 : 1;
    const float selective_perturb_pct = std::clamp(
        0.08f + 0.55f * congestion_severity, 0.05f, 0.55f);
    const float selective_critical_pct
        = std::clamp(6.0f + 4.0f * (congestion_severity + hotspot_bias),
                     5.5f,
                     11.5f);
    const int selective_seed = snapshot.seed + 23;

    ScenarioDefinition selective_relief;
    selective_relief.name = "focused-soft";
    selective_relief.pre_init
        = [this, selective_perturb_pct, selective_seed, selective_critical_pct]() {
            grouter_->setCapacitiesPerturbationPercentage(
                selective_perturb_pct);
            grouter_->setPerturbationAmount(selective_perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(selective_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(
                selective_critical_pct);
          };
    selective_relief.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           selective_threshold,
           selective_min_ratio,
           selective_hotspot_push,
           selective_halo]() {
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 selective_threshold,
                                 selective_min_ratio,
                                 selective_hotspot_push,
                                 selective_halo);
          };
    scenario_defs.push_back(selective_relief);
  }

  if (has_congestion_data && run_aggressive_soft) {
    const float tuned = std::clamp(congestion_severity * 0.65f
                                       + hotspot_bias * 0.35f
                                       + (force_routability ? 0.15f : 0.0f),
                                   0.0f,
                                   1.0f);
    const float min_base
        = std::clamp(0.80f - 0.09f * tuned, 0.72f, 0.88f);
    const float max_base
        = std::clamp(0.97f - 0.05f * tuned, min_base + 0.035f, 0.99f);
    const float slope = 2.2f + 1.0f * tuned;
    const float midpoint
        = std::clamp(0.50f - 0.04f * tuned, 0.42f, 0.50f);
    const int halo = tuned > 0.65f ? 2 : 1;
    const float hotspot_ratio
        = std::clamp(0.92f - 0.14f * tuned, 0.74f, 0.95f);
    const float severity_weight
        = std::clamp(0.30f + 0.40f * tuned, 0.30f, 0.85f);
    const float perturb_pct
        = std::clamp(0.35f + 1.6f * tuned, 0.20f, 2.40f);
    const float critical_pct
        = std::clamp(7.0f + 5.5f * tuned, 6.5f, 14.0f);
    const int soft_seed = snapshot.seed + 11;
    scenario_defs.push_back(make_soft_config("soft-relief",
                                             min_base,
                                             max_base,
                                             slope,
                                             midpoint,
                                             halo,
                                             hotspot_ratio,
                                             severity_weight,
                                             perturb_pct,
                                             soft_seed,
                                             critical_pct));
  }

  for (const ScenarioDefinition& def : scenario_defs) {
    ScenarioResult result = run_scenario(def, snapshot);
    scenario_results.push_back(std::move(result));
  }

  auto better_result = [&](const ScenarioResult& lhs,
                           const ScenarioResult& rhs) {
    if (lhs.metrics.overflow != rhs.metrics.overflow) {
      return lhs.metrics.overflow < rhs.metrics.overflow;
    }

    const double wl_a = static_cast<double>(lhs.metrics.wirelength_dbu);
    const double wl_b = static_cast<double>(rhs.metrics.wirelength_dbu);
    const double wl_den = std::max(std::max(wl_a, wl_b), 1.0);
    const double wl_rel = std::abs(wl_a - wl_b) / wl_den;
    const double wl_priority = 0.00035;   // ~0.035% difference
    const double wl_via_tie = 0.00090;    // ~0.090% difference

    if (wl_a != wl_b && wl_rel > wl_priority) {
      // Wirelength dominates until differences are small.
      return wl_a < wl_b;
    }

    if (lhs.metrics.via_count != rhs.metrics.via_count && wl_rel < wl_via_tie) {
      // With comparable wirelengths, prefer fewer vias.
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }

    const double util_diff
        = std::abs(lhs.metrics.max_utilization - rhs.metrics.max_utilization);
    if (wl_rel < 0.0015 && util_diff > 0.005) {
      // Allow small wirelength trade-offs when congestion relief is clear.
      return lhs.metrics.max_utilization < rhs.metrics.max_utilization;
    }

    auto relative_gap = [](double a, double b) {
      const double denom = std::max({std::abs(a), std::abs(b), 1e-9});
      return (a - b) / denom;
    };

    if (wl_rel < 0.0015) {
      const double stress_rel
          = relative_gap(lhs.metrics.stress_cost, rhs.metrics.stress_cost);
      if (std::abs(stress_rel) > 0.05) {
        return stress_rel < 0.0;
      }

      const double reserve_rel
          = relative_gap(lhs.metrics.reserve_score, rhs.metrics.reserve_score);
      if (std::abs(reserve_rel) > 0.05) {
        // Prefer solutions that leave a little more slack.
        return reserve_rel > 0.0;
      }
    }

    if (wl_a != wl_b) {
      return wl_a < wl_b;
    }

    if (std::abs(lhs.metrics.max_utilization - rhs.metrics.max_utilization)
        > 1e-4) {
      return lhs.metrics.max_utilization < rhs.metrics.max_utilization;
    }
    if (std::abs(lhs.metrics.stress_cost - rhs.metrics.stress_cost) > 1e-4) {
      return lhs.metrics.stress_cost < rhs.metrics.stress_cost;
    }
    return lhs.metrics.score < rhs.metrics.score;
  };

  auto best_iter = std::min_element(
      scenario_results.begin(), scenario_results.end(), better_result);
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
                "overflow {}, max util {:.2f}",
                final_result.name,
                final_result.metrics.wirelength_um,
                final_result.metrics.via_count,
                final_result.metrics.overflow,
                final_result.metrics.max_utilization);

  return std::move(final_result.routes);
}

}  // namespace grt
