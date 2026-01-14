#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_set>
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
  float via_cost_scale = 1.0f;
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

struct PatchSummary
{
  int nets_touched = 0;
  int segments_added = 0;
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

void applyCoolCapacityBoost(GlobalRouter* grouter,
                            const RudyGrid& normalized_rudy,
                            const std::vector<Hotspot>& hotspots,
                            int min_layer,
                            int max_layer,
                            float cool_threshold,
                            float base_boost,
                            float max_boost,
                            float layer_decay,
                            int halo,
                            float hotspot_guard = 0.6f)
{
  Grid* grid = grouter->grid();
  if (grid == nullptr || normalized_rudy.empty()) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return;
  }

  const int x_tiles = static_cast<int>(normalized_rudy.size());
  const int y_tiles
      = (normalized_rudy.empty() || normalized_rudy.front().empty())
            ? 0
            : static_cast<int>(normalized_rudy.front().size());

  const auto getNormalized = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= x_tiles || y >= y_tiles || x_tiles == 0
        || y_tiles == 0) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  std::vector<std::vector<float>> boost_mask(
      std::max(x_grids, 1), std::vector<float>(std::max(y_grids, 1), 1.0f));
  halo = std::max(0, halo);
  hotspot_guard = std::clamp(hotspot_guard, 0.3f, 1.0f);
  if (!hotspots.empty() && halo > 0) {
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
          const float guard = 1.0f - (1.0f - hotspot_guard) * decay;
          boost_mask[gx][gy] = std::min(boost_mask[gx][gy], guard);
        }
      }
    }
  }

  cool_threshold = std::clamp(cool_threshold, 0.05f, 1.0f);
  base_boost = std::clamp(base_boost, 0.0f, 1.0f);
  max_boost = std::clamp(max_boost, 1.0f, 1.30f);
  layer_decay = std::clamp(layer_decay, 0.0f, 0.8f);

  const int layer_span = std::max(max_layer - min_layer, 1);
  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    float layer_scale = 1.0f - layer_decay * layer_factor;
    layer_scale = std::clamp(layer_scale, 0.55f, 1.0f);

    for (int y = 0; y < y_grids; ++y) {
      for (int x = 0; x < x_grids - 1; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        if (normalized >= cool_threshold) {
          continue;
        }
        const float cool_score = (cool_threshold - normalized)
                                 / std::max(cool_threshold, 1e-3f);
        const float max_ratio = 1.0f + (max_boost - 1.0f) * layer_scale;
        float ratio = 1.0f + base_boost * cool_score * layer_scale;
        ratio = std::clamp(ratio, 1.0f, max_ratio);
        ratio = 1.0f + (ratio - 1.0f) * boost_mask[x][y];
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < y_grids - 1; ++y) {
      for (int x = 0; x < x_grids; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        if (normalized >= cool_threshold) {
          continue;
        }
        const float cool_score = (cool_threshold - normalized)
                                 / std::max(cool_threshold, 1e-3f);
        const float max_ratio = 1.0f + (max_boost - 1.0f) * layer_scale;
        float ratio = 1.0f + base_boost * cool_score * layer_scale;
        ratio = std::clamp(ratio, 1.0f, max_ratio);
        ratio = 1.0f + (ratio - 1.0f) * boost_mask[x][y];
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

void applyTopLayerBias(GlobalRouter* grouter,
                       const RudyGrid& normalized_rudy,
                       const std::vector<Hotspot>& hotspots,
                       int min_layer,
                       int max_layer,
                       int top_k,
                       float cool_threshold,
                       float base_boost,
                       float max_boost,
                       float hotspot_guard,
                       int halo)
{
  Grid* grid = grouter->grid();
  if (grid == nullptr || normalized_rudy.empty() || top_k <= 0) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return;
  }

  const int layer_start = std::max(min_layer, max_layer - top_k + 1);
  if (layer_start > max_layer) {
    return;
  }

  const int x_tiles = static_cast<int>(normalized_rudy.size());
  const int y_tiles
      = (normalized_rudy.empty() || normalized_rudy.front().empty())
            ? 0
            : static_cast<int>(normalized_rudy.front().size());

  const auto getNormalized = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= x_tiles || y >= y_tiles || x_tiles == 0
        || y_tiles == 0) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  std::vector<std::vector<float>> boost_mask(
      std::max(x_grids, 1), std::vector<float>(std::max(y_grids, 1), 1.0f));
  halo = std::max(halo, 0);
  hotspot_guard = std::clamp(hotspot_guard, 0.4f, 1.0f);
  if (!hotspots.empty() && halo > 0) {
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
          const float guard = 1.0f - (1.0f - hotspot_guard) * decay;
          boost_mask[gx][gy] = std::min(boost_mask[gx][gy], guard);
        }
      }
    }
  }

  cool_threshold = std::clamp(cool_threshold, 0.30f, 0.95f);
  base_boost = std::clamp(base_boost, 0.0f, 0.3f);
  max_boost = std::clamp(max_boost, 1.0f, 1.35f);

  const int layer_span = std::max(max_layer - layer_start, 1);
  for (int layer = layer_start; layer <= max_layer; ++layer) {
    const float layer_rel
        = static_cast<float>(layer - layer_start) / static_cast<float>(layer_span);
    float layer_scale = 0.90f + 0.18f * layer_rel;
    layer_scale = std::clamp(layer_scale, 0.85f, 1.12f);

    for (int y = 0; y < y_grids; ++y) {
      for (int x = 0; x < x_grids - 1; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        float cool_score = (cool_threshold - normalized)
                           / std::max(cool_threshold, 1e-3f);
        cool_score = std::clamp(cool_score, 0.0f, 1.0f);
        if (cool_score <= 0.0f) {
          continue;
        }

        float boost = base_boost * cool_score
                      + (max_boost - 1.0f - base_boost) * cool_score
                            * cool_score;
        float ratio = 1.0f + layer_scale * boost;
        const float max_layer_ratio
            = 1.0f + (max_boost - 1.0f) * layer_scale;
        ratio = std::clamp(ratio, 1.0f, max_layer_ratio);
        ratio = 1.0f + (ratio - 1.0f) * boost_mask[x][y];
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < y_grids - 1; ++y) {
      for (int x = 0; x < x_grids; ++x) {
        const float normalized
            = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        float cool_score = (cool_threshold - normalized)
                           / std::max(cool_threshold, 1e-3f);
        cool_score = std::clamp(cool_score, 0.0f, 1.0f);
        if (cool_score <= 0.0f) {
          continue;
        }

        float boost = base_boost * cool_score
                      + (max_boost - 1.0f - base_boost) * cool_score
                            * cool_score;
        float ratio = 1.0f + layer_scale * boost;
        const float max_layer_ratio
            = 1.0f + (max_boost - 1.0f) * layer_scale;
        ratio = std::clamp(ratio, 1.0f, max_layer_ratio);
        ratio = 1.0f + (ratio - 1.0f) * boost_mask[x][y];
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

bool segmentCoversHotspot(const GSegment& segment,
                          const Hotspot& hotspot,
                          int tile_size,
                          int x_min,
                          int y_min)
{
  if (tile_size <= 0) {
    return false;
  }

  const int min_x = std::min(segment.init_x, segment.final_x);
  const int max_x = std::max(segment.init_x, segment.final_x);
  const int min_y = std::min(segment.init_y, segment.final_y);
  const int max_y = std::max(segment.init_y, segment.final_y);

  const int gx0 = (min_x - x_min) / tile_size;
  const int gx1 = (max_x - x_min) / tile_size;
  const int gy0 = (min_y - y_min) / tile_size;
  const int gy1 = (max_y - y_min) / tile_size;

  return hotspot.gx >= gx0 && hotspot.gx <= gx1 && hotspot.gy >= gy0
         && hotspot.gy <= gy1;
}

PatchSummary applyHotspotPatches(NetRouteMap& routes,
                                 const std::vector<Hotspot>& hotspots,
                                 Grid* grid,
                                 int min_layer,
                                 int max_layer)
{
  PatchSummary summary;
  if (hotspots.empty() || grid == nullptr) {
    return summary;
  }

  const int tile_size = grid->getTileSize();
  if (tile_size <= 0) {
    return summary;
  }
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();

  // Focus on the worst congestion cells first.
  std::vector<Hotspot> prioritized = hotspots;
  std::sort(prioritized.begin(),
            prioritized.end(),
            [](const Hotspot& a, const Hotspot& b) {
              return a.severity > b.severity;
            });
  const float severity_floor = 1.05f;
  const size_t max_hotspots = 20;
  const size_t capped = std::min(max_hotspots, prioritized.size());

  size_t hotspot_limit = 0;
  while (hotspot_limit < capped
         && prioritized[hotspot_limit].severity >= severity_floor) {
    hotspot_limit++;
  }
  if (hotspot_limit > 0) {
    prioritized.resize(hotspot_limit);
  } else {
    prioritized.resize(std::min<size_t>(8, capped));
  }

  const int per_net_budget = 2;

  for (auto& [db_net, segments] : routes) {
    std::unordered_set<GSegment, GSegmentHash> seen(segments.begin(),
                                                    segments.end());
    int patches_used = 0;
    bool net_marked = false;

    for (const Hotspot& hotspot : prioritized) {
      if (patches_used >= per_net_budget) {
        break;
      }

      if (hotspot.severity < 1.02f) {
        continue;
      }

      int target_layer = -1;
      for (const GSegment& segment : segments) {
        if (segmentCoversHotspot(
                segment, hotspot, tile_size, x_min, y_min)) {
          target_layer = segment.init_layer;
          if (segment.init_layer != segment.final_layer) {
            target_layer = std::min(segment.init_layer, segment.final_layer);
          }
          break;
        }
      }

      if (target_layer < 0) {
        continue;
      }

      target_layer = std::clamp(target_layer, min_layer, max_layer);
      const int alt_layer
          = (target_layer < max_layer)
                ? target_layer + 1
                : (target_layer > min_layer ? target_layer - 1 : target_layer);

      const int x0 = x_min + hotspot.gx * tile_size;
      const int x1 = x_min + (hotspot.gx + 1) * tile_size;
      const int y0 = y_min + hotspot.gy * tile_size;
      const int y1 = y_min + (hotspot.gy + 1) * tile_size;

      bool added_patch_segments = false;
      auto add_segment = [&](int ix,
                             int iy,
                             int il,
                             int fx,
                             int fy,
                             int fl) {
        GSegment seg(ix, iy, il, fx, fy, fl);
        if (seen.insert(seg).second) {
          segments.push_back(seg);
          summary.segments_added++;
          added_patch_segments = true;
          if (!net_marked) {
            summary.nets_touched++;
            net_marked = true;
          }
        }
      };

      add_segment(x0, y0, target_layer, x1, y0, target_layer);
      add_segment(x0, y0, target_layer, x0, y1, target_layer);
      if (alt_layer != target_layer) {
        add_segment(x0, y0, target_layer, x0, y0, alt_layer);
        add_segment(x0, y0, alt_layer, x1, y0, alt_layer);
      }
      if (added_patch_segments) {
        patches_used++;
      }
    }
  }

  return summary;
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
      core->computeCongestionInformation();
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

    // Lean harder on wirelength for scenario ranking; still guard overflow.
    const double via_weight = static_cast<double>(tile_size) * 1.2;
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
    if (grouter_->fastroute_ != nullptr) {
      snapshot.critical_percentage
          = grouter_->fastroute_->getCriticalNetsPercentage();
      snapshot.via_cost_scale = grouter_->fastroute_->getViaCostScale();
    }
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    if (grouter_->fastroute_ != nullptr) {
      grouter_->fastroute_->setCriticalNetsPercentage(
          snapshot.critical_percentage);
      grouter_->fastroute_->setViaCostScale(snapshot.via_cost_scale);
    }
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
        if (usage_ratio < 0.68f) {
          continue;
        }
        const int gx
            = std::clamp((info.segment.init_x - x_min) / tile_size, 0, x_grids);
        const int gy
            = std::clamp((info.segment.init_y - y_min) / tile_size, 0, y_grids);
        Hotspot hotspot;
        hotspot.gx = std::clamp(gx, 0, std::max(x_grids - 1, 0));
        hotspot.gy = std::clamp(gy, 0, std::max(y_grids - 1, 0));
        hotspot.severity = std::clamp(usage_ratio, 0.7f, 2.4f);
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
  float congestion_severity
      = computeCongestionSeverity(rudy_stats, hotspots);
  if (baseline.metrics.overflow == 0) {
    const float relief_from_rudy
        = std::clamp(0.22f + 0.20f * std::max(0.0f, 0.90f - rudy_stats.p80)
                         + 0.10f * std::max(0.0f, 0.75f - rudy_stats.mean),
                     0.0f,
                     0.55f);
    const float hotspot_relief
        = 0.06f * std::max(0, 3 - static_cast<int>(hotspots.size()));
    congestion_severity
        = std::clamp(congestion_severity - relief_from_rudy - hotspot_relief,
                     0.0f,
                     1.0f);
  }
  const bool force_routability = baseline.metrics.overflow > 0;
  const bool has_congestion_data
      = (!normalized_rudy.empty() || !hotspots.empty());
  const float hotspot_bias = std::clamp(
      static_cast<float>(hotspots.size()) / 12.0f, 0.0f, 0.6f);
  const bool relaxed_utilization = baseline.metrics.max_utilization < 0.70f;
  const bool light_congestion = !force_routability
                                && congestion_severity < 0.70f
                                && hotspots.size() <= 3
                                && rudy_stats.p80 < 0.92f;
  const bool run_soft
      = force_routability
        || (has_congestion_data && !light_congestion
            && (congestion_severity > 0.70f || hotspots.size() > 3
                || rudy_stats.p80 > 0.90f))
        || hotspots.size() > 5;
  const bool run_aggressive_soft
      = force_routability
        || (has_congestion_data && !light_congestion
            && (congestion_severity > 0.86f || hotspots.size() > 4));
  const bool allow_seed_sweep = !force_routability && hotspots.size() <= 6
                                && congestion_severity > 0.38f;
  const bool allow_light_seed = !force_routability && light_congestion
                                && congestion_severity > 0.34f;
  const bool ultra_light = light_congestion
                           && baseline.metrics.overflow == 0
                           && hotspots.size() <= 2
                           && rudy_stats.p80 < 0.90f
                           && baseline.metrics.max_utilization < 0.70f;

  logger_->info(GNR,
                6008,
                "NEWGR congestion severity {:.2f} (RUDY mean {:.2f}, p80 "
                "{:.2f}, hotspots {}, baseline overflow {}, light {})",
                congestion_severity,
                rudy_stats.mean,
                rudy_stats.p80,
                hotspots.size(),
                baseline.metrics.overflow,
                light_congestion);

  float wl_via_scale = 1.0f;
  if (!force_routability) {
    float base_via_scale = 0.60f + 0.25f * congestion_severity
                           - 0.12f * hotspot_bias;
    if (relaxed_utilization) {
      base_via_scale -= 0.05f;
    }
    if (light_congestion) {
      base_via_scale -= 0.06f;
    }
    wl_via_scale = std::clamp(base_via_scale, 0.58f, 1.0f);
    if (baseline.metrics.overflow == 0 && light_congestion
        && hotspot_bias < 0.22f) {
      wl_via_scale = std::clamp(wl_via_scale - 0.10f, 0.50f, 0.96f);
    }
  } else {
    const float base_via_scale
        = 1.0f + 0.12f * congestion_severity + 0.08f * hotspot_bias;
    wl_via_scale = std::clamp(base_via_scale, 1.0f, 1.20f);
  }
  const auto apply_wl_via_scale = [this, wl_via_scale]() {
    if (grouter_->fastroute_ != nullptr) {
      grouter_->fastroute_->setViaCostScale(wl_via_scale);
    }
  };

  const float via_trim_scale_base
      = std::max(wl_via_scale * 1.05f,
                 0.98f + 0.30f * std::max(0.0f, 0.70f - congestion_severity)
                     - 0.10f * hotspot_bias);
  const float via_trim_scale = std::clamp(
      via_trim_scale_base,
      force_routability ? 0.95f : 1.02f,
      force_routability ? 1.32f : 1.28f);
  const auto apply_via_trim_scale = [this, via_trim_scale]() {
    if (grouter_->fastroute_ != nullptr) {
      grouter_->fastroute_->setViaCostScale(via_trim_scale);
    }
  };

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
          def.pre_init = [this, perturb_pct, seed, critical_pct, apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
            grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
            apply_wl_via_scale();
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
  float wl_perturb_pct
      = congestion_severity > 0.45f ? 0.22f : 0.0f;
  if (!force_routability && wl_perturb_pct > 0.0f
      && congestion_severity < 0.70f) {
    const float taper = std::clamp(
        0.75f + 0.25f * (congestion_severity / 0.70f), 0.70f, 0.98f);
    wl_perturb_pct *= taper;
  }
  const int wl_seed = snapshot.seed + 5;
  ScenarioDefinition wl_variation;
  wl_variation.name = "wl-variation";
  wl_variation.pre_init
      = [this,
         wl_perturb_pct,
         wl_seed,
         wl_critical_pct,
         apply_wl_via_scale]() {
          grouter_->setCapacitiesPerturbationPercentage(wl_perturb_pct);
          grouter_->setPerturbationAmount(wl_perturb_pct > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_critical_pct);
          apply_wl_via_scale();
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

  float wl_greedy_perturb = congestion_severity > 0.30f ? 0.12f : 0.0f;
  if (relaxed_utilization && wl_greedy_perturb > 0.0f) {
    wl_greedy_perturb *= 0.55f;
  }
  if (!force_routability && wl_greedy_perturb > 0.0f
      && congestion_severity < 0.70f) {
    const float taper = std::clamp(
        0.58f + 0.42f * (congestion_severity / 0.70f), 0.58f, 1.0f);
    wl_greedy_perturb *= taper;
  }
  if (!force_routability && light_congestion && wl_greedy_perturb > 0.0f) {
    wl_greedy_perturb *= 0.85f;
  }
  float wl_greedy_critical
      = std::clamp(wl_critical_pct - 0.8f, 3.5f, 10.0f);
  if (relaxed_utilization) {
    wl_greedy_critical = std::clamp(wl_greedy_critical + 0.6f, 3.8f, 11.0f);
  }
  const float wl_greedy_hotspot_ratio
      = std::clamp(0.995f - 0.04f * hotspot_bias, 0.97f, 0.995f);
  const float wl_greedy_hotspot_weight
      = std::clamp(0.06f + 0.12f * hotspot_bias, 0.06f, 0.14f);
  const int wl_greedy_hotspot_halo = hotspots.size() > 1 ? 2 : 1;
  const bool wl_greedy_skip_hotspot_penalties
      = hotspots.empty()
        || (!force_routability && baseline.metrics.overflow == 0
            && light_congestion && hotspot_bias < 0.20f
            && hotspots.size() <= 2 && congestion_severity < 0.62f);
  const bool wl_greedy_hotspot_guarded
      = force_routability || baseline.metrics.overflow > 0
        || congestion_severity > 0.70f || hotspot_bias > 0.22f
        || hotspots.size() > 3;
  const int wl_greedy_seed = snapshot.seed + 37;
  ScenarioDefinition wl_greedy;
  wl_greedy.name = "wl-greedy";
  wl_greedy.pre_init
      = [this,
         wl_greedy_perturb,
         wl_greedy_seed,
         wl_greedy_critical,
         apply_wl_via_scale]() {
          grouter_->setCapacitiesPerturbationPercentage(wl_greedy_perturb);
          grouter_->setPerturbationAmount(wl_greedy_perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_greedy_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_greedy_critical);
          apply_wl_via_scale();
        };
  wl_greedy.post_init
      = [this,
         &hotspots,
         min_routing_layer,
         max_routing_layer,
         wl_greedy_hotspot_ratio,
         wl_greedy_hotspot_weight,
         wl_greedy_hotspot_halo,
         wl_greedy_skip_hotspot_penalties,
         wl_greedy_hotspot_guarded,
         light_congestion]() {
          if (wl_greedy_skip_hotspot_penalties) {
            return;
          }

          const float hotspot_ratio
              = wl_greedy_hotspot_guarded
                    ? wl_greedy_hotspot_ratio
                    : std::clamp(wl_greedy_hotspot_ratio + 0.006f,
                                 0.0f,
                                 0.999f);
          const float hotspot_weight
              = wl_greedy_hotspot_weight
                * (wl_greedy_hotspot_guarded
                       ? 1.0f
                       : (light_congestion ? 0.55f : 0.70f));

          applyHotspotPenalties(grouter_,
                                hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                wl_greedy_hotspot_halo,
                                hotspot_ratio,
                                hotspot_weight);
        };
  scenario_defs.push_back(wl_greedy);

  if (ultra_light) {
    const float feather_perturb
        = std::clamp(wl_greedy_perturb * 0.65f, 0.0f, 0.08f);
    const float feather_critical
        = std::clamp(wl_greedy_critical - 0.5f, 3.2f, 8.5f);
    const int feather_seed = snapshot.seed + 977;
    const float feather_via_scale
        = std::clamp(wl_via_scale * 0.88f, 0.36f, 0.92f);
    const int feather_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float feather_top_threshold
        = std::clamp(0.54f + 0.04f * (congestion_severity - 0.40f),
                     0.50f,
                     0.66f);
    const float feather_top_base = std::clamp(
        0.05f + 0.03f * (0.60f - congestion_severity), 0.04f, 0.09f);
    const float feather_top_max = std::clamp(
        1.10f + 0.04f * (0.55f - congestion_severity), 1.06f, 1.16f);
    const float feather_top_guard
        = std::clamp(0.82f - 0.10f * hotspot_bias, 0.70f, 0.90f);
    const int feather_top_halo = 1;
    const float feather_cool_threshold
        = std::clamp(feather_top_threshold * 0.78f, 0.46f, 0.64f);
    const float feather_cool_base = std::clamp(
        0.05f + 0.03f * (0.60f - congestion_severity), 0.04f, 0.09f);
    const float feather_cool_max = std::clamp(
        1.08f + 0.04f * (0.55f - congestion_severity), 1.06f, 1.14f);
    const float feather_cool_decay
        = std::clamp(0.04f + 0.05f * hotspot_bias, 0.03f, 0.10f);
    const float feather_hotspot_guard
        = std::clamp(0.82f - 0.10f * hotspot_bias, 0.72f, 0.90f);
    const float feather_hotspot_ratio
        = std::clamp(0.996f - 0.02f * hotspot_bias, 0.98f, 0.998f);
    const float feather_hotspot_weight
        = std::clamp(0.04f + 0.06f * hotspot_bias, 0.04f, 0.08f);

    ScenarioDefinition wl_feather;
    wl_feather.name = "wl-feather";
    wl_feather.pre_init
        = [this,
           feather_perturb,
           feather_seed,
           feather_critical,
           feather_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(feather_perturb);
            grouter_->setPerturbationAmount(feather_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(feather_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(feather_critical);
            if (grouter_->fastroute_ != nullptr) {
              grouter_->fastroute_->setViaCostScale(feather_via_scale);
            }
          };
    wl_feather.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           feather_top_k,
           feather_top_threshold,
           feather_top_base,
           feather_top_max,
           feather_top_guard,
           feather_top_halo,
           feather_cool_threshold,
           feather_cool_base,
           feather_cool_max,
           feather_cool_decay,
           feather_hotspot_guard,
           feather_hotspot_ratio,
           feather_hotspot_weight]() {
            if (!normalized_rudy.empty()) {
              applyTopLayerBias(grouter_,
                                normalized_rudy,
                                hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                feather_top_k,
                                feather_top_threshold,
                                feather_top_base,
                                feather_top_max,
                                feather_top_guard,
                                feather_top_halo);
              applyCoolCapacityBoost(grouter_,
                                     normalized_rudy,
                                     hotspots,
                                     min_routing_layer,
                                     max_routing_layer,
                                     feather_cool_threshold,
                                     feather_cool_base,
                                     feather_cool_max,
                                     feather_cool_decay,
                                     feather_top_halo,
                                     feather_hotspot_guard);
            }
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    feather_top_halo,
                                    feather_hotspot_ratio,
                                    feather_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_feather);
  }

  const float wl_refine_perturb = wl_greedy_perturb * 0.6f;
  const float wl_refine_critical
      = std::clamp(wl_greedy_critical - 0.6f, 3.0f, 10.0f);
  const int wl_refine_seed = snapshot.seed + 71;
  ScenarioDefinition wl_refine;
  wl_refine.name = "wl-refine";
  wl_refine.pre_init
      = [this,
         wl_refine_perturb,
         wl_refine_seed,
         wl_refine_critical,
         apply_wl_via_scale]() {
          const float perturb = std::clamp(wl_refine_perturb, 0.0f, 0.40f);
          grouter_->setCapacitiesPerturbationPercentage(perturb);
          grouter_->setPerturbationAmount(perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(wl_refine_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(wl_refine_critical);
          apply_wl_via_scale();
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

  const float lean_strength
      = std::clamp(0.35f + 0.40f * congestion_severity
                       + 0.20f * hotspot_bias,
                   0.30f,
                   0.80f);
  const float lean_perturb
      = std::clamp(0.025f + 0.12f * congestion_severity, 0.02f, 0.16f);
  const float lean_critical = std::clamp(
      wl_greedy_critical - 0.3f + 0.4f * hotspot_bias, 3.2f, 10.5f);
  const int lean_seed = snapshot.seed + 311;
  const float lean_min_base = std::clamp(
      0.92f + 0.03f * (1.0f - congestion_severity), 0.91f, 0.96f);
  const float lean_max_base = std::clamp(
      0.985f + 0.01f * (0.65f - congestion_severity),
      lean_min_base + 0.01f,
      0.99f);
  const float lean_slope = 1.6f + 0.9f * lean_strength;
  const float lean_midpoint
      = std::clamp(0.50f - 0.02f * lean_strength, 0.45f, 0.52f);
  const int lean_top_k = std::max(
      1, std::min(2, max_routing_layer - min_routing_layer + 1));
  const float lean_top_threshold = std::clamp(
      0.60f + 0.06f * (congestion_severity - 0.50f), 0.56f, 0.72f);
  const float lean_top_base = std::clamp(
      0.05f + 0.04f * (0.60f - congestion_severity), 0.04f, 0.11f);
  const float lean_top_max = std::clamp(
      1.11f + 0.04f * (0.55f - congestion_severity), 1.08f, 1.18f);
  const float lean_top_guard
      = std::clamp(0.72f - 0.14f * hotspot_bias, 0.62f, 0.82f);
  const int lean_top_halo = hotspots.size() > 2 ? 2 : 1;
  const float lean_hotspot_ratio
      = std::clamp(0.988f - 0.04f * hotspot_bias, 0.95f, 0.992f);
  const float lean_hotspot_weight
      = std::clamp(0.08f + 0.10f * hotspot_bias, 0.06f, 0.16f);

  ScenarioDefinition wl_lean;
  wl_lean.name = "wl-lean";
  wl_lean.pre_init
      = [this, lean_perturb, lean_seed, lean_critical, apply_wl_via_scale]() {
          grouter_->setCapacitiesPerturbationPercentage(lean_perturb);
          grouter_->setPerturbationAmount(lean_perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(lean_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(lean_critical);
          apply_wl_via_scale();
        };
  wl_lean.post_init
      = [this,
         &normalized_rudy,
         &hotspots,
         min_routing_layer,
         max_routing_layer,
         lean_min_base,
         lean_max_base,
         lean_slope,
         lean_midpoint,
         lean_top_k,
         lean_top_threshold,
         lean_top_base,
         lean_top_max,
         lean_top_guard,
         lean_top_halo,
         lean_hotspot_ratio,
          lean_hotspot_weight]() {
          if (!normalized_rudy.empty()) {
            applySoftCapacityScaling(grouter_,
                                     normalized_rudy,
                                     min_routing_layer,
                                     max_routing_layer,
                                     lean_min_base,
                                     lean_max_base,
                                     lean_slope,
                                     lean_midpoint);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              lean_top_k,
                              lean_top_threshold,
                              lean_top_base,
                              lean_top_max,
                              lean_top_guard,
                              lean_top_halo);
          }
          if (!hotspots.empty()) {
            applyHotspotPenalties(grouter_,
                                  hotspots,
                                  min_routing_layer,
                                  max_routing_layer,
                                  lean_top_halo,
                                  lean_hotspot_ratio,
                                  lean_hotspot_weight);
          }
        };
  scenario_defs.push_back(wl_lean);

  const float balance_perturb
      = std::clamp(0.01f + 0.06f * congestion_severity
                       + 0.02f * hotspot_bias,
                   0.0f,
                   0.10f);
  const float balance_critical = std::clamp(
      6.2f + 2.2f * (0.60f - congestion_severity), 4.8f, 9.5f);
  const int balance_seed = snapshot.seed + 179;
  const float balance_min_base
      = std::clamp(0.945f + 0.02f * (0.55f - congestion_severity),
                   0.93f,
                   0.97f);
  const float balance_max_base
      = std::clamp(0.987f - 0.01f * congestion_severity,
                   balance_min_base + 0.01f,
                   0.995f);
  const float balance_slope = 1.35f + 0.55f * congestion_severity;
  const float balance_midpoint
      = std::clamp(0.49f - 0.01f * (1.0f - congestion_severity), 0.47f, 0.50f);
  const float balance_threshold = std::clamp(
      0.64f + 0.04f * congestion_severity, 0.60f, 0.72f);
  const float balance_min_ratio
      = std::clamp(0.965f - 0.01f * (1.0f - congestion_severity),
                   0.95f,
                   0.985f);
  const float balance_hotspot_push
      = std::clamp(0.05f + 0.06f * hotspot_bias, 0.04f, 0.12f);
  const int balance_halo = hotspots.size() > 3 ? 2 : 1;
  const int balance_top_k = std::max(
      1, std::min(2, max_routing_layer - min_routing_layer + 1));
  const float balance_top_threshold
      = std::clamp(0.60f + 0.06f * congestion_severity, 0.58f, 0.70f);
  const float balance_top_base = std::clamp(
      0.05f + 0.03f * (0.60f - congestion_severity), 0.04f, 0.10f);
  const float balance_top_max = std::clamp(
      1.10f + 0.03f * (0.55f - congestion_severity), 1.06f, 1.14f);
  const float balance_top_guard
      = std::clamp(0.78f - 0.12f * hotspot_bias, 0.66f, 0.82f);
  const int balance_top_halo = hotspots.size() > 2 ? 2 : 1;
  const float balance_hotspot_ratio
      = std::clamp(0.994f - 0.04f * hotspot_bias, 0.96f, 0.995f);
  const float balance_hotspot_weight
      = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);
  const bool balance_skip_hotspots
      = baseline.metrics.overflow == 0 && light_congestion
        && hotspot_bias < 0.24f;

  ScenarioDefinition wl_balance;
  wl_balance.name = "wl-balance";
  wl_balance.pre_init
      = [this,
         balance_perturb,
         balance_seed,
         balance_critical,
         apply_wl_via_scale]() {
          grouter_->setCapacitiesPerturbationPercentage(balance_perturb);
          grouter_->setPerturbationAmount(balance_perturb > 0.0f ? 1 : 0);
          grouter_->setSeed(balance_seed);
          grouter_->setAllowCongestion(false);
          grouter_->fastroute_->setCriticalNetsPercentage(balance_critical);
          apply_wl_via_scale();
        };
  wl_balance.post_init
      = [this,
         &normalized_rudy,
         &hotspots,
         min_routing_layer,
         max_routing_layer,
         balance_min_base,
         balance_max_base,
         balance_slope,
         balance_midpoint,
         balance_threshold,
         balance_min_ratio,
         balance_hotspot_push,
         balance_halo,
         balance_top_k,
         balance_top_threshold,
         balance_top_base,
         balance_top_max,
         balance_top_guard,
         balance_top_halo,
         balance_hotspot_ratio,
         balance_hotspot_weight,
         balance_skip_hotspots]() {
          if (!normalized_rudy.empty()) {
            applySoftCapacityScaling(grouter_,
                                     normalized_rudy,
                                     min_routing_layer,
                                     max_routing_layer,
                                     balance_min_base,
                                     balance_max_base,
                                     balance_slope,
                                     balance_midpoint);
          }
          applySelectiveRelief(grouter_,
                               normalized_rudy,
                               hotspots,
                               min_routing_layer,
                               max_routing_layer,
                               balance_threshold,
                               balance_min_ratio,
                               balance_hotspot_push,
                               balance_halo);
          applyTopLayerBias(grouter_,
                            normalized_rudy,
                            hotspots,
                            min_routing_layer,
                            max_routing_layer,
                            balance_top_k,
                            balance_top_threshold,
                            balance_top_base,
                            balance_top_max,
                            balance_top_guard,
                            balance_top_halo);
          if (!hotspots.empty() && !balance_skip_hotspots) {
            applyHotspotPenalties(grouter_,
                                  hotspots,
                                  min_routing_layer,
                                  max_routing_layer,
                                  balance_halo,
                                  balance_hotspot_ratio,
                                  balance_hotspot_weight);
          }
        };
  scenario_defs.push_back(wl_balance);

  if (baseline.metrics.overflow == 0 && has_congestion_data
      && baseline.metrics.max_utilization < 0.72f) {
    const float directlite_perturb = std::clamp(
        0.01f + 0.08f * congestion_severity, 0.01f, 0.10f);
    const float directlite_critical = std::clamp(
        6.5f + 3.0f * (0.60f - congestion_severity), 5.0f, 10.5f);
    const int directlite_seed = snapshot.seed + 287;
    const float directlite_cool_threshold
        = std::clamp(0.60f - 0.10f * hotspot_bias, 0.48f, 0.64f);
    const float directlite_base_boost = std::clamp(
        0.09f + 0.04f * (0.65f - congestion_severity), 0.06f, 0.12f);
    const float directlite_max_boost = std::clamp(
        1.10f + 0.04f * (0.55f - congestion_severity), 1.08f, 1.16f);
    const float directlite_layer_decay
        = std::clamp(0.06f + 0.08f * hotspot_bias, 0.05f, 0.18f);
    const int directlite_halo = hotspots.size() > 3 ? 2 : 1;
    const float directlite_hotspot_guard
        = std::clamp(0.70f - 0.10f * hotspot_bias, 0.60f, 0.74f);
    const int directlite_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float directlite_top_threshold
        = std::clamp(0.62f + 0.08f * congestion_severity, 0.60f, 0.78f);
    const float directlite_top_base = std::clamp(
        0.08f + 0.04f * (0.60f - congestion_severity), 0.06f, 0.12f);
    const float directlite_top_max = std::clamp(
        1.12f + 0.05f * (0.55f - congestion_severity), 1.10f, 1.20f);
    const float directlite_top_guard
        = std::clamp(0.70f - 0.16f * hotspot_bias, 0.60f, 0.80f);
    const int directlite_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float directlite_hotspot_ratio
        = std::clamp(0.992f - 0.05f * hotspot_bias, 0.95f, 0.995f);
    const float directlite_hotspot_weight
        = std::clamp(0.08f + 0.10f * hotspot_bias, 0.06f, 0.14f);

    ScenarioDefinition wl_directlite;
    wl_directlite.name = "wl-direct-lite";
    wl_directlite.pre_init
        = [this,
           directlite_perturb,
           directlite_seed,
           directlite_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(directlite_perturb);
            grouter_->setPerturbationAmount(directlite_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(directlite_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(
                directlite_critical);
            apply_wl_via_scale();
          };
    wl_directlite.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           directlite_cool_threshold,
           directlite_base_boost,
           directlite_max_boost,
           directlite_layer_decay,
           directlite_halo,
           directlite_hotspot_guard,
           directlite_top_k,
           directlite_top_threshold,
           directlite_top_base,
           directlite_top_max,
           directlite_top_guard,
           directlite_top_halo,
           directlite_hotspot_ratio,
           directlite_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   directlite_cool_threshold,
                                   directlite_base_boost,
                                   directlite_max_boost,
                                   directlite_layer_decay,
                                   directlite_halo,
                                   directlite_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              directlite_top_k,
                              directlite_top_threshold,
                              directlite_top_base,
                              directlite_top_max,
                              directlite_top_guard,
                              directlite_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    directlite_halo,
                                    directlite_hotspot_ratio,
                                    directlite_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_directlite);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data && light_congestion) {
    const float openlane_perturb
        = std::clamp(0.01f + 0.10f * congestion_severity, 0.01f, 0.14f);
    const float openlane_critical = std::clamp(
        5.5f + 2.0f * (0.58f - congestion_severity), 4.5f, 9.5f);
    const int openlane_seed = snapshot.seed + 641;
    const float openlane_cool_threshold
        = std::clamp(0.55f + 0.08f * congestion_severity, 0.50f, 0.68f);
    const float openlane_base_boost = std::clamp(
        0.05f + 0.05f * (0.60f - congestion_severity), 0.04f, 0.11f);
    const float openlane_max_boost = std::clamp(
        1.08f + 0.05f * (0.55f - congestion_severity), 1.06f, 1.16f);
    const float openlane_layer_decay
        = std::clamp(0.04f + 0.06f * hotspot_bias, 0.04f, 0.12f);
    const int openlane_halo = hotspots.size() > 3 ? 2 : 1;
    const float openlane_hotspot_guard = std::clamp(
        0.68f - 0.12f * hotspot_bias, 0.60f, 0.76f);
    const int openlane_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float openlane_top_threshold
        = std::clamp(0.58f + 0.06f * congestion_severity, 0.56f, 0.72f);
    const float openlane_top_base = std::clamp(
        0.05f + 0.05f * (0.58f - congestion_severity), 0.04f, 0.11f);
    const float openlane_top_max = std::clamp(
        1.10f + 0.05f * (0.55f - congestion_severity), 1.08f, 1.18f);
    const float openlane_top_guard
        = std::clamp(0.70f - 0.14f * hotspot_bias, 0.60f, 0.78f);
    const int openlane_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float openlane_hotspot_ratio
        = std::clamp(0.992f - 0.03f * hotspot_bias, 0.96f, 0.995f);
    const float openlane_hotspot_weight
        = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_openlane;
    wl_openlane.name = "wl-openlane";
    wl_openlane.pre_init
        = [this,
           openlane_perturb,
           openlane_seed,
           openlane_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(openlane_perturb);
            grouter_->setPerturbationAmount(openlane_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(openlane_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(openlane_critical);
            apply_wl_via_scale();
          };
    wl_openlane.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           openlane_cool_threshold,
           openlane_base_boost,
           openlane_max_boost,
           openlane_layer_decay,
           openlane_halo,
           openlane_hotspot_guard,
           openlane_top_k,
           openlane_top_threshold,
           openlane_top_base,
           openlane_top_max,
           openlane_top_guard,
           openlane_top_halo,
           openlane_hotspot_ratio,
           openlane_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   openlane_cool_threshold,
                                   openlane_base_boost,
                                   openlane_max_boost,
                                   openlane_layer_decay,
                                   openlane_halo,
                                   openlane_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              openlane_top_k,
                              openlane_top_threshold,
                              openlane_top_base,
                              openlane_top_max,
                              openlane_top_guard,
                              openlane_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    openlane_halo,
                                    openlane_hotspot_ratio,
                                    openlane_hotspot_weight);
           }
          };
    scenario_defs.push_back(wl_openlane);

    const float shortcut_via_scale = std::clamp(
        wl_via_scale * (0.60f + 0.12f * (0.70f - congestion_severity)),
        0.36f,
        0.82f);
    const float shortcut_perturb = std::clamp(
        0.008f + 0.10f * (0.65f - congestion_severity), 0.0f, 0.10f);
    const float shortcut_critical = std::clamp(
        4.2f + 2.0f * (0.60f - congestion_severity), 3.6f, 8.2f);
    const int shortcut_seed = snapshot.seed + 911;
    const int shortcut_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float shortcut_top_threshold
        = std::clamp(0.60f + 0.06f * congestion_severity, 0.58f, 0.72f);
    const float shortcut_top_base = std::clamp(
        0.09f + 0.05f * (0.60f - congestion_severity), 0.08f, 0.14f);
    const float shortcut_top_max = std::clamp(
        1.16f + 0.05f * (0.55f - congestion_severity), 1.12f, 1.26f);
    const float shortcut_top_guard
        = std::clamp(0.74f - 0.12f * hotspot_bias, 0.64f, 0.82f);
    const int shortcut_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float shortcut_cool_threshold
        = std::clamp(0.62f + 0.06f * (0.60f - congestion_severity),
                     0.58f,
                     0.70f);
    const float shortcut_base_boost = std::clamp(
        0.10f + 0.05f * (0.55f - congestion_severity), 0.08f, 0.16f);
    const float shortcut_max_boost = std::clamp(
        1.18f + 0.06f * (0.55f - congestion_severity), 1.12f, 1.26f);
    const float shortcut_layer_decay
        = std::clamp(0.05f + 0.05f * hotspot_bias, 0.04f, 0.12f);
    const float shortcut_hotspot_guard
        = std::clamp(0.72f - 0.12f * hotspot_bias, 0.62f, 0.82f);
    const float shortcut_hotspot_ratio
        = std::clamp(0.992f - 0.04f * hotspot_bias, 0.96f, 0.995f);
    const float shortcut_hotspot_weight
        = std::clamp(0.06f + 0.08f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_shortcut;
    wl_shortcut.name = "wl-shortcut";
    wl_shortcut.pre_init
        = [this,
           shortcut_perturb,
           shortcut_seed,
           shortcut_critical,
           shortcut_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(shortcut_perturb);
            grouter_->setPerturbationAmount(shortcut_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(shortcut_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(shortcut_critical);
            if (grouter_->fastroute_ != nullptr) {
              grouter_->fastroute_->setViaCostScale(shortcut_via_scale);
            }
          };
    wl_shortcut.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           shortcut_cool_threshold,
           shortcut_base_boost,
           shortcut_max_boost,
           shortcut_layer_decay,
           shortcut_hotspot_guard,
           shortcut_top_k,
           shortcut_top_threshold,
           shortcut_top_base,
           shortcut_top_max,
           shortcut_top_guard,
           shortcut_top_halo,
           shortcut_hotspot_ratio,
           shortcut_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   shortcut_cool_threshold,
                                   shortcut_base_boost,
                                   shortcut_max_boost,
                                   shortcut_layer_decay,
                                   shortcut_top_halo,
                                   shortcut_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              shortcut_top_k,
                              shortcut_top_threshold,
                              shortcut_top_base,
                              shortcut_top_max,
                              shortcut_top_guard,
                              shortcut_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    shortcut_top_halo,
                                    shortcut_hotspot_ratio,
                                    shortcut_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_shortcut);

    const float skim_strength
        = std::clamp(0.65f - congestion_severity, 0.0f, 0.35f);
    const float skim_perturb = std::clamp(
        0.015f + 0.08f * congestion_severity - 0.02f * hotspot_bias,
        0.0f,
        0.10f);
    const float skim_critical = std::clamp(
        4.5f + 2.0f * (0.60f - congestion_severity), 3.8f, 9.0f);
    const int skim_seed = snapshot.seed + 569;
    const float skim_via_scale = std::clamp(
        wl_via_scale * (0.78f + 0.14f * (0.65f - congestion_severity)),
        0.42f,
        0.86f);
    const float skim_cool_threshold = std::clamp(
        0.56f + 0.10f * (0.60f - congestion_severity), 0.48f, 0.70f);
    const float skim_base_boost
        = std::clamp(0.09f + 0.05f * skim_strength, 0.07f, 0.15f);
    const float skim_max_boost
        = std::clamp(1.14f + 0.08f * skim_strength, 1.12f, 1.30f);
    const float skim_layer_decay
        = std::clamp(0.04f + 0.08f * hotspot_bias, 0.04f, 0.14f);
    const int skim_halo = hotspots.size() > 2 ? 2 : 1;
    const float skim_hotspot_guard
        = std::clamp(0.70f - 0.12f * hotspot_bias, 0.60f, 0.82f);
    const int skim_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float skim_top_threshold
        = std::clamp(0.52f + 0.08f * congestion_severity, 0.50f, 0.68f);
    const float skim_top_base
        = std::clamp(0.08f + 0.05f * skim_strength, 0.06f, 0.15f);
    const float skim_top_max
        = std::clamp(1.15f + 0.09f * skim_strength, 1.12f, 1.30f);
    const float skim_top_guard
        = std::clamp(0.74f - 0.12f * hotspot_bias, 0.64f, 0.84f);
    const int skim_top_halo = hotspots.size() > 1 ? 2 : 1;
    const float skim_hotspot_ratio
        = std::clamp(0.994f - 0.04f * hotspot_bias, 0.96f, 0.996f);
    const float skim_hotspot_weight
        = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_skim;
    wl_skim.name = "wl-skim";
    wl_skim.pre_init
        = [this,
           skim_perturb,
           skim_seed,
           skim_critical,
           skim_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(skim_perturb);
            grouter_->setPerturbationAmount(skim_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(skim_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(skim_critical);
            if (grouter_->fastroute_ != nullptr) {
              grouter_->fastroute_->setViaCostScale(skim_via_scale);
            }
          };
    wl_skim.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           skim_cool_threshold,
           skim_base_boost,
           skim_max_boost,
           skim_layer_decay,
           skim_halo,
           skim_hotspot_guard,
           skim_top_k,
           skim_top_threshold,
           skim_top_base,
           skim_top_max,
           skim_top_guard,
           skim_top_halo,
           skim_hotspot_ratio,
           skim_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   skim_cool_threshold,
                                   skim_base_boost,
                                   skim_max_boost,
                                   skim_layer_decay,
                                   skim_halo,
                                   skim_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              skim_top_k,
                              skim_top_threshold,
                              skim_top_base,
                              skim_top_max,
                              skim_top_guard,
                              skim_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    skim_halo,
                                    skim_hotspot_ratio,
                                    skim_hotspot_weight);
           }
          };
    scenario_defs.push_back(wl_skim);

    const float direct_top_perturb
        = std::clamp(0.01f + 0.06f * congestion_severity, 0.0f, 0.10f);
    const float direct_top_critical = std::clamp(
        4.8f + 2.0f * (0.60f - congestion_severity), 4.0f, 9.0f);
    const int direct_top_seed = snapshot.seed + 887;
    const float direct_top_via_scale
        = std::clamp(wl_via_scale * 0.72f, 0.38f, 0.82f);
    const int direct_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float direct_top_threshold
        = std::clamp(0.55f + 0.08f * congestion_severity, 0.52f, 0.70f);
    const float direct_top_base = std::clamp(
        0.12f + 0.05f * (0.60f - congestion_severity), 0.10f, 0.18f);
    const float direct_top_max = std::clamp(
        1.18f + 0.06f * (0.55f - congestion_severity), 1.14f, 1.26f);
    const float direct_top_guard
        = std::clamp(0.74f - 0.12f * hotspot_bias, 0.64f, 0.82f);
    const int direct_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float direct_cool_threshold = std::clamp(
        0.62f - 0.06f * hotspot_bias, 0.50f, 0.68f);
    const float direct_base_boost = std::clamp(
        0.11f + 0.05f * (0.60f - congestion_severity), 0.09f, 0.16f);
    const float direct_max_boost = std::clamp(
        1.17f + 0.05f * (0.55f - congestion_severity), 1.12f, 1.26f);
    const float direct_layer_decay
        = std::clamp(0.05f + 0.08f * hotspot_bias, 0.04f, 0.14f);
    const int direct_halo = hotspots.size() > 3 ? 2 : 1;
    const float direct_hotspot_guard
        = std::clamp(0.70f - 0.10f * hotspot_bias, 0.60f, 0.76f);
    const float direct_hotspot_ratio
        = std::clamp(0.994f - 0.04f * hotspot_bias, 0.96f, 0.995f);
    const float direct_hotspot_weight
        = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_direct_top;
    wl_direct_top.name = "wl-direct-top";
    wl_direct_top.pre_init
        = [this,
           direct_top_perturb,
           direct_top_seed,
           direct_top_critical,
           direct_top_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(direct_top_perturb);
            grouter_->setPerturbationAmount(direct_top_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(direct_top_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(
                direct_top_critical);
            if (grouter_->fastroute_ != nullptr) {
              grouter_->fastroute_->setViaCostScale(direct_top_via_scale);
            }
          };
    wl_direct_top.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           direct_top_k,
           direct_top_threshold,
           direct_top_base,
           direct_top_max,
           direct_top_guard,
           direct_top_halo,
           direct_cool_threshold,
           direct_base_boost,
           direct_max_boost,
           direct_layer_decay,
           direct_halo,
           direct_hotspot_guard,
           direct_hotspot_ratio,
           direct_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   direct_cool_threshold,
                                   direct_base_boost,
                                   direct_max_boost,
                                   direct_layer_decay,
                                   direct_halo,
                                   direct_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              direct_top_k,
                              direct_top_threshold,
                              direct_top_base,
                              direct_top_max,
                              direct_top_guard,
                              direct_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    direct_halo,
                                    direct_hotspot_ratio,
                                    direct_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_direct_top);

    const float prime_perturb
        = std::clamp(wl_greedy_perturb * 0.35f, 0.0f, 0.08f);
    const float prime_critical
        = std::clamp(wl_greedy_critical - 0.3f, 3.2f, 8.5f);
    const int prime_seed = snapshot.seed + 457;
    const int prime_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float prime_top_threshold
        = std::clamp(0.60f + 0.06f * congestion_severity, 0.58f, 0.74f);
    const float prime_top_base = std::clamp(
        0.05f + 0.04f * (0.55f - congestion_severity), 0.04f, 0.10f);
    const float prime_top_max = std::clamp(
        1.08f + 0.05f * (0.55f - congestion_severity), 1.06f, 1.16f);
    const float prime_top_guard
        = std::clamp(0.78f - 0.10f * hotspot_bias, 0.70f, 0.84f);
    const int prime_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float prime_hotspot_ratio
        = std::clamp(0.994f - 0.03f * hotspot_bias, 0.97f, 0.995f);
    const float prime_hotspot_weight
        = std::clamp(0.05f + 0.08f * hotspot_bias, 0.04f, 0.12f);

    ScenarioDefinition wl_prime;
    wl_prime.name = "wl-direct-prime";
    wl_prime.pre_init
        = [this,
           prime_perturb,
           prime_seed,
           prime_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(prime_perturb);
            grouter_->setPerturbationAmount(prime_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(prime_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(prime_critical);
            apply_wl_via_scale();
          };
    wl_prime.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           prime_top_k,
           prime_top_threshold,
           prime_top_base,
           prime_top_max,
           prime_top_guard,
           prime_top_halo,
           prime_hotspot_ratio,
           prime_hotspot_weight]() {
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              prime_top_k,
                              prime_top_threshold,
                              prime_top_base,
                              prime_top_max,
                              prime_top_guard,
                              prime_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    prime_top_halo,
                                    prime_hotspot_ratio,
                                    prime_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_prime);

    const float compact_perturb
        = std::clamp(0.015f + 0.12f * congestion_severity, 0.0f, 0.16f);
    const float compact_critical = std::clamp(
        5.0f + 2.5f * (0.60f - congestion_severity), 4.2f, 9.0f);
    const int compact_seed = snapshot.seed + 503;
    const float compact_min_base = 0.94f;
    const float compact_max_base = 0.985f;
    const float compact_slope
        = 1.6f + 0.7f * std::clamp(congestion_severity, 0.0f, 1.0f);
    const float compact_midpoint = 0.48f;
    const float compact_threshold = std::clamp(
        0.64f + 0.08f * congestion_severity, 0.60f, 0.74f);
    const float compact_min_ratio
        = std::clamp(0.96f + 0.03f * (0.60f - congestion_severity),
                     0.95f,
                     0.985f);
    const float compact_hotspot_push
        = std::clamp(0.05f + 0.05f * hotspot_bias, 0.04f, 0.12f);
    const int compact_halo = hotspots.size() > 3 ? 2 : 1;
    const float compact_cool_threshold
        = std::clamp(compact_threshold * 0.72f, 0.42f, 0.56f);
    const float compact_boost = std::clamp(
        0.08f + 0.05f * (0.55f - congestion_severity), 0.06f, 0.14f);
    const float compact_boost_limit = std::clamp(
        1.07f + 0.05f * (0.55f - congestion_severity), 1.05f, 1.14f);
    const float compact_layer_falloff
        = std::clamp(0.08f + 0.08f * hotspot_bias, 0.07f, 0.16f);
    const int compact_top_k = std::max(
        1, std::min(2, max_routing_layer - min_routing_layer + 1));
    const float compact_top_threshold
        = std::clamp(0.58f + 0.10f * congestion_severity, 0.56f, 0.74f);
    const float compact_top_base = std::clamp(
        0.06f + 0.05f * (0.55f - congestion_severity), 0.04f, 0.12f);
    const float compact_top_max = std::clamp(
        1.09f + 0.05f * (0.55f - congestion_severity), 1.06f, 1.16f);
    const float compact_top_guard
        = std::clamp(0.70f - 0.20f * hotspot_bias, 0.60f, 0.74f);
    const int compact_top_halo = hotspots.size() > 2 ? 2 : 1;
    const float compact_hotspot_ratio
        = std::clamp(0.992f - 0.03f * hotspot_bias, 0.96f, 0.995f);
    const float compact_hotspot_weight
        = std::clamp(0.08f + 0.08f * hotspot_bias, 0.06f, 0.14f);

    ScenarioDefinition wl_compact;
    wl_compact.name = "wl-compact";
    wl_compact.pre_init
        = [this,
           compact_perturb,
           compact_seed,
           compact_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(compact_perturb);
            grouter_->setPerturbationAmount(compact_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(compact_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(compact_critical);
            apply_wl_via_scale();
          };
    wl_compact.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           compact_min_base,
           compact_max_base,
           compact_slope,
           compact_midpoint,
           compact_threshold,
           compact_min_ratio,
           compact_hotspot_push,
           compact_halo,
           compact_cool_threshold,
           compact_boost,
           compact_boost_limit,
           compact_layer_falloff,
           compact_top_k,
           compact_top_threshold,
           compact_top_base,
           compact_top_max,
           compact_top_guard,
           compact_top_halo,
           compact_hotspot_ratio,
           compact_hotspot_weight]() {
            if (!normalized_rudy.empty()) {
              applySoftCapacityScaling(grouter_,
                                       normalized_rudy,
                                       min_routing_layer,
                                       max_routing_layer,
                                       compact_min_base,
                                       compact_max_base,
                                       compact_slope,
                                       compact_midpoint);
            }

            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 compact_threshold,
                                 compact_min_ratio,
                                 compact_hotspot_push,
                                 compact_halo,
                                 compact_cool_threshold,
                                 compact_boost,
                                 compact_boost_limit,
                                 compact_layer_falloff);

            if (!normalized_rudy.empty()) {
              applyTopLayerBias(grouter_,
                                normalized_rudy,
                                hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                compact_top_k,
                                compact_top_threshold,
                                compact_top_base,
                                compact_top_max,
                                compact_top_guard,
                                compact_top_halo);
            }

            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    compact_halo,
                                    compact_hotspot_ratio,
                                    compact_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_compact);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data
      && baseline.metrics.max_utilization < 0.82f) {
    float via_trim_perturb = std::clamp(
        0.01f + 0.10f * congestion_severity + 0.05f * hotspot_bias,
        0.0f,
        0.20f);
    const float via_trim_critical = std::clamp(
        5.4f + 2.2f * (0.60f - congestion_severity), 4.0f, 10.5f);
    const int via_trim_seed = snapshot.seed + 823;
    const float via_trim_min_base = std::clamp(
        0.965f - 0.01f * hotspot_bias, 0.95f, 0.975f);
    const float via_trim_max_base = std::clamp(
        0.995f - 0.01f * congestion_severity,
        via_trim_min_base + 0.01f,
        0.997f);
    const float via_trim_slope = 1.3f + 0.5f * congestion_severity;
    const float via_trim_midpoint = 0.48f;
    const float via_trim_threshold = std::clamp(
        0.60f + 0.06f * congestion_severity, 0.58f, 0.72f);
    const float via_trim_min_ratio = std::clamp(
        0.97f - 0.015f * hotspot_bias, 0.95f, 0.985f);
    float via_trim_hotspot_push
        = std::clamp(0.06f + 0.08f * hotspot_bias, 0.05f, 0.16f);
    const int via_trim_halo = hotspots.size() > 3 ? 2 : 1;
    const float via_trim_cool_threshold
        = std::clamp(via_trim_threshold * 0.70f, 0.44f, 0.60f);
    float via_trim_boost = std::clamp(
        0.09f + 0.05f * (0.60f - congestion_severity), 0.06f, 0.16f);
    float via_trim_boost_limit = std::clamp(
        1.10f + 0.05f * (0.55f - congestion_severity), 1.08f, 1.18f);
    const float via_trim_layer_falloff
        = std::clamp(0.10f + 0.10f * hotspot_bias, 0.08f, 0.22f);
    const float via_trim_hotspot_ratio
        = std::clamp(0.993f - 0.04f * hotspot_bias, 0.96f, 0.995f);
    const float via_trim_hotspot_weight
        = std::clamp(0.08f + 0.10f * hotspot_bias, 0.06f, 0.16f);

    if (ultra_light) {
      via_trim_perturb = std::clamp(via_trim_perturb * 0.55f, 0.0f, 0.20f);
      via_trim_hotspot_push
          = std::clamp(via_trim_hotspot_push * 0.75f, 0.04f, 0.16f);
      via_trim_boost = std::clamp(via_trim_boost * 0.82f, 0.05f, 0.16f);
      via_trim_boost_limit
          = std::clamp(via_trim_boost_limit * 0.96f, 1.05f, 1.18f);
    }

    ScenarioDefinition via_trim;
    via_trim.name = "via-trim";
    via_trim.pre_init
        = [this,
           via_trim_perturb,
           via_trim_seed,
           via_trim_critical,
           apply_via_trim_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(via_trim_perturb);
            grouter_->setPerturbationAmount(via_trim_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(via_trim_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(via_trim_critical);
            apply_via_trim_scale();
          };
    via_trim.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           via_trim_min_base,
           via_trim_max_base,
           via_trim_slope,
           via_trim_midpoint,
           via_trim_threshold,
           via_trim_min_ratio,
           via_trim_hotspot_push,
           via_trim_halo,
           via_trim_cool_threshold,
           via_trim_boost,
           via_trim_boost_limit,
           via_trim_layer_falloff,
           via_trim_hotspot_ratio,
           via_trim_hotspot_weight]() {
            if (!normalized_rudy.empty()) {
              applySoftCapacityScaling(grouter_,
                                       normalized_rudy,
                                       min_routing_layer,
                                       max_routing_layer,
                                       via_trim_min_base,
                                       via_trim_max_base,
                                       via_trim_slope,
                                       via_trim_midpoint);
            }

            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 via_trim_threshold,
                                 via_trim_min_ratio,
                                 via_trim_hotspot_push,
                                 via_trim_halo,
                                 via_trim_cool_threshold,
                                 via_trim_boost,
                                 via_trim_boost_limit,
                                 via_trim_layer_falloff);

            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    via_trim_halo,
                                    via_trim_hotspot_ratio,
                                    via_trim_hotspot_weight);
            }
          };
    scenario_defs.push_back(via_trim);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data) {
    const float direct_threshold = std::clamp(
        0.60f + 0.12f * (congestion_severity - 0.40f), 0.54f, 0.80f);
    const float direct_min_ratio = std::clamp(
        0.94f + 0.06f * (0.55f - congestion_severity), 0.93f, 0.985f);
    const float direct_hotspot_push
        = std::clamp(0.05f + 0.08f * hotspot_bias, 0.04f, 0.14f);
    const int direct_halo = hotspots.size() > 4 ? 2 : 1;
    const float direct_cool_threshold
        = std::clamp(direct_threshold * 0.82f, 0.44f, 0.66f);
    const float direct_boost = std::clamp(
        0.12f + 0.07f * (0.55f - congestion_severity), 0.10f, 0.20f);
    const float direct_boost_limit = std::clamp(
        1.09f + 0.05f * (0.55f - congestion_severity), 1.08f, 1.16f);
    const float direct_layer_falloff
        = std::clamp(0.09f + 0.10f * hotspot_bias, 0.08f, 0.18f);
    const float direct_perturb = std::clamp(
        0.04f + 0.16f * (0.60f - congestion_severity), 0.02f, 0.16f);
    const float direct_critical = std::clamp(
        6.0f + 2.6f * (0.55f - congestion_severity), 5.2f, 8.8f);
    const int direct_seed = snapshot.seed + 389;
    const float direct_hotspot_ratio
        = std::clamp(0.992f - 0.03f * hotspot_bias, 0.96f, 0.995f);
    const float direct_hotspot_weight
        = std::clamp(0.10f + 0.08f * hotspot_bias, 0.08f, 0.16f);

    ScenarioDefinition wl_direct;
    wl_direct.name = "wl-direct";
    wl_direct.pre_init
        = [this,
           direct_perturb,
           direct_seed,
           direct_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(direct_perturb);
            grouter_->setPerturbationAmount(direct_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(direct_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(direct_critical);
            apply_wl_via_scale();
          };
    wl_direct.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           direct_threshold,
           direct_min_ratio,
           direct_hotspot_push,
           direct_halo,
           direct_cool_threshold,
           direct_boost,
           direct_boost_limit,
           direct_layer_falloff,
           direct_hotspot_ratio,
           direct_hotspot_weight]() {
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 direct_threshold,
                                 direct_min_ratio,
                                 direct_hotspot_push,
                                 direct_halo,
                                 direct_cool_threshold,
                                 direct_boost,
                                 direct_boost_limit,
                                 direct_layer_falloff);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    direct_halo,
                                    direct_hotspot_ratio,
                                    direct_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_direct);
  }

  auto make_wl_greedy_seed
      = [this,
         wl_greedy_perturb,
         wl_greedy_critical,
         &snapshot,
         apply_wl_via_scale](
            const std::string& name, int seed_offset, float perturb_scale) {
          ScenarioDefinition def;
          def.name = name;
          const float perturb
              = std::clamp(wl_greedy_perturb * perturb_scale, 0.0f, 0.60f);
          const int seed = snapshot.seed + seed_offset;
          def.pre_init
              = [this, perturb, seed, wl_greedy_critical, apply_wl_via_scale]() {
                  grouter_->setCapacitiesPerturbationPercentage(perturb);
                  grouter_->setPerturbationAmount(perturb > 0.0f ? 1 : 0);
                  grouter_->setSeed(seed);
                  grouter_->setAllowCongestion(false);
                  grouter_->fastroute_->setCriticalNetsPercentage(
                      wl_greedy_critical);
                  apply_wl_via_scale();
                };
          def.post_init = []() {};
          return def;
        };

  if (allow_light_seed) {
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s1", 73, 0.85f));
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-l2", 137, 1.05f));
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s2", 109, 0.70f));
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s3", 177, 1.20f));
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-zero", 251, 0.0f));
  } else if (allow_seed_sweep) {
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s1", 73, 0.8f));
    if (hotspots.size() <= 2 || congestion_severity > 0.78f) {
      scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s2", 109, 1.05f));
    }
    scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s3", 157, 0.6f));
    if (hotspots.size() <= 3 && congestion_severity > 0.78f) {
      scenario_defs.push_back(make_wl_greedy_seed("wl-greedy-s4", 211, 1.15f));
    }
  }

  if (has_congestion_data && light_congestion) {
    const float polish_threshold
        = std::clamp(0.58f + 0.12f * (1.0f - congestion_severity),
                     0.54f,
                     0.70f);
    const float polish_min_ratio
        = std::clamp(0.92f + 0.05f * (1.0f - congestion_severity),
                     0.90f,
                     0.97f);
    const float polish_hotspot_push
        = std::clamp(0.07f + 0.06f * hotspot_bias, 0.05f, 0.14f);
    const int polish_halo = hotspots.size() > 4 ? 2 : 1;
    const float polish_cool_threshold
        = std::clamp(polish_threshold * 0.75f, 0.42f, 0.58f);
    const float polish_boost = std::clamp(
        0.11f + 0.06f * (0.55f - congestion_severity), 0.08f, 0.16f);
    const float polish_boost_limit = std::clamp(
        1.10f + 0.04f * (0.55f - congestion_severity), 1.06f, 1.14f);
    const float polish_layer_falloff
        = std::clamp(0.10f + 0.10f * hotspot_bias, 0.08f, 0.20f);
    const float polish_perturb
        = std::clamp(0.02f + 0.10f * congestion_severity, 0.02f, 0.12f);
    const float polish_critical = std::clamp(
        8.0f + 3.5f * (0.60f - congestion_severity), 6.5f, 12.5f);
    const int polish_seed = snapshot.seed + 233;
    const float polish_hotspot_ratio = std::clamp(
        0.99f - 0.05f * congestion_severity, 0.94f, 0.99f);
    const float polish_hotspot_weight
        = std::clamp(0.08f + 0.20f * hotspot_bias, 0.08f, 0.18f);

    ScenarioDefinition wl_polish;
    wl_polish.name = "wl-polish-lite";
    wl_polish.pre_init
        = [this,
           polish_perturb,
           polish_seed,
           polish_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(polish_perturb);
            grouter_->setPerturbationAmount(polish_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(polish_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(polish_critical);
            apply_wl_via_scale();
          };
    wl_polish.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           polish_threshold,
           polish_min_ratio,
           polish_hotspot_push,
           polish_halo,
           polish_cool_threshold,
           polish_boost,
           polish_boost_limit,
           polish_layer_falloff,
           polish_hotspot_ratio,
           polish_hotspot_weight]() {
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 polish_threshold,
                                 polish_min_ratio,
                                 polish_hotspot_push,
                                 polish_halo,
                                 polish_cool_threshold,
                                 polish_boost,
                                 polish_boost_limit,
                                 polish_layer_falloff);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    polish_halo,
                                    polish_hotspot_ratio,
                                    polish_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_polish);
  }

  if (has_congestion_data && light_congestion && baseline.metrics.overflow == 0) {
    const float streamline_threshold = std::clamp(
        0.46f + 0.10f * congestion_severity, 0.44f, 0.62f);
    const float streamline_min_ratio = std::clamp(
        0.955f + 0.03f * (1.0f - congestion_severity), 0.94f, 0.985f);
    const float streamline_hotspot_push
        = std::clamp(0.06f + 0.05f * hotspot_bias, 0.05f, 0.13f);
    const int streamline_halo = hotspots.size() > 4 ? 2 : 1;
    const float streamline_cool_threshold = std::clamp(
        streamline_threshold * 0.62f, 0.34f, 0.50f);
    const float streamline_boost = std::clamp(
        0.12f + 0.05f * (0.55f - congestion_severity), 0.09f, 0.18f);
    const float streamline_boost_limit = std::clamp(
        1.10f + 0.06f * (0.55f - congestion_severity), 1.08f, 1.18f);
    const float streamline_layer_falloff
        = std::clamp(0.07f + 0.08f * hotspot_bias, 0.06f, 0.16f);
    const float streamline_perturb
        = std::clamp(0.01f + 0.06f * congestion_severity, 0.0f, 0.12f);
    const float streamline_critical = std::clamp(
        4.2f + 2.4f * (0.60f - congestion_severity), 3.8f, 8.5f);
    const int streamline_seed = snapshot.seed + 319;
    const float streamline_hotspot_ratio
        = std::clamp(0.992f - 0.03f * hotspot_bias, 0.96f, 0.995f);
    const float streamline_hotspot_weight
        = std::clamp(0.06f + 0.14f * hotspot_bias, 0.06f, 0.16f);

    ScenarioDefinition wl_streamline;
    wl_streamline.name = "wl-streamline";
    wl_streamline.pre_init
        = [this,
           streamline_perturb,
           streamline_seed,
           streamline_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(streamline_perturb);
            grouter_->setPerturbationAmount(streamline_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(streamline_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(
                streamline_critical);
            apply_wl_via_scale();
          };
    wl_streamline.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           streamline_threshold,
           streamline_min_ratio,
           streamline_hotspot_push,
           streamline_halo,
           streamline_cool_threshold,
           streamline_boost,
           streamline_boost_limit,
           streamline_layer_falloff,
           streamline_hotspot_ratio,
           streamline_hotspot_weight]() {
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 streamline_threshold,
                                 streamline_min_ratio,
                                 streamline_hotspot_push,
                                 streamline_halo,
                                 streamline_cool_threshold,
                                 streamline_boost,
                                 streamline_boost_limit,
                                 streamline_layer_falloff);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    streamline_halo,
                                    streamline_hotspot_ratio,
                                    streamline_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_streamline);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data
      && !normalized_rudy.empty()) {
    const int top_k
        = std::max(1, std::min(3, max_routing_layer - min_routing_layer + 1));
    const float top_threshold = std::clamp(
        0.60f + 0.18f * congestion_severity, 0.58f, 0.82f);
    const float top_base_boost = std::clamp(
        0.09f + 0.06f * (0.65f - congestion_severity), 0.06f, 0.16f);
    const float top_max_boost = std::clamp(
        1.14f + 0.05f * (0.55f - congestion_severity), 1.10f, 1.22f);
    const float top_guard
        = std::clamp(0.70f - 0.25f * hotspot_bias, 0.52f, 0.74f);
    const int top_halo = hotspots.size() > 2 ? 2 : 1;
    const float top_perturb = std::clamp(
        0.04f + 0.18f * (0.65f - congestion_severity), 0.02f, 0.22f);
    const float top_critical = std::clamp(
        7.5f + 3.0f * (0.55f - congestion_severity), 6.0f, 11.0f);
    const int top_seed = snapshot.seed + 427;
    const float top_hotspot_ratio
        = std::clamp(0.995f - 0.04f * hotspot_bias, 0.96f, 0.995f);
    const float top_hotspot_weight
        = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_toplane;
    wl_toplane.name = "wl-toplane";
    wl_toplane.pre_init
        = [this,
           top_perturb,
           top_seed,
           top_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(top_perturb);
            grouter_->setPerturbationAmount(top_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(top_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(top_critical);
            apply_wl_via_scale();
          };
    wl_toplane.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           top_k,
           top_threshold,
           top_base_boost,
           top_max_boost,
           top_guard,
           top_halo,
           top_hotspot_ratio,
           top_hotspot_weight]() {
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              top_k,
                              top_threshold,
                              top_base_boost,
                              top_max_boost,
                              top_guard,
                              top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    top_halo,
                                    top_hotspot_ratio,
                                    top_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_toplane);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data) {
    const float squeeze_perturb = std::clamp(
        0.02f + 0.12f * (0.55f - congestion_severity), 0.0f, 0.14f);
    const float squeeze_critical = std::clamp(
        wl_critical_pct + 0.8f * (0.60f - congestion_severity), 4.0f, 11.0f);
    const int squeeze_seed = snapshot.seed + 709;
    const float squeeze_via_scale = std::clamp(
        wl_via_scale * (0.82f + 0.10f * (0.65f - congestion_severity)),
        0.52f,
        0.98f);
    const float squeeze_cool_threshold = std::clamp(
        0.55f + 0.08f * (0.65f - congestion_severity), 0.50f, 0.70f);
    const float squeeze_base_boost
        = std::clamp(0.08f + 0.05f * (0.60f - congestion_severity),
                     0.06f,
                     0.14f);
    const float squeeze_max_boost
        = std::clamp(1.15f + 0.06f * (0.55f - congestion_severity),
                     1.10f,
                     1.24f);
    const float squeeze_layer_decay
        = std::clamp(0.05f + 0.06f * hotspot_bias, 0.04f, 0.14f);
    const int squeeze_halo = hotspots.size() > 2 ? 2 : 1;
    const float squeeze_hotspot_guard
        = std::clamp(0.70f - 0.14f * hotspot_bias, 0.60f, 0.78f);
    const int squeeze_top_k = std::max(
        2, std::min(3, max_routing_layer - min_routing_layer + 1));
    const float squeeze_top_threshold
        = std::clamp(0.60f + 0.10f * congestion_severity, 0.58f, 0.78f);
    const float squeeze_top_base = std::clamp(
        0.09f + 0.05f * (0.58f - congestion_severity), 0.08f, 0.16f);
    const float squeeze_top_max = std::clamp(
        1.16f + 0.06f * (0.55f - congestion_severity), 1.12f, 1.28f);
    const float squeeze_top_guard
        = std::clamp(0.74f - 0.18f * hotspot_bias, 0.62f, 0.84f);
    const int squeeze_top_halo = hotspots.size() > 1 ? 2 : 1;
    const float squeeze_hotspot_ratio
        = std::clamp(0.994f - 0.03f * hotspot_bias, 0.97f, 0.995f);
    const float squeeze_hotspot_weight
        = std::clamp(0.06f + 0.10f * hotspot_bias, 0.05f, 0.14f);

    ScenarioDefinition wl_squeeze;
    wl_squeeze.name = "wl-squeeze";
    wl_squeeze.pre_init
        = [this,
           squeeze_perturb,
           squeeze_seed,
           squeeze_critical,
           squeeze_via_scale]() {
            const float perturb = std::clamp(squeeze_perturb, 0.0f, 0.20f);
            grouter_->setCapacitiesPerturbationPercentage(perturb);
            grouter_->setPerturbationAmount(perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(squeeze_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(squeeze_critical);
            if (grouter_->fastroute_ != nullptr) {
              grouter_->fastroute_->setViaCostScale(squeeze_via_scale);
            }
          };
    wl_squeeze.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           squeeze_cool_threshold,
           squeeze_base_boost,
           squeeze_max_boost,
           squeeze_layer_decay,
           squeeze_halo,
           squeeze_hotspot_guard,
           squeeze_top_k,
           squeeze_top_threshold,
           squeeze_top_base,
           squeeze_top_max,
           squeeze_top_guard,
           squeeze_top_halo,
           squeeze_hotspot_ratio,
           squeeze_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   squeeze_cool_threshold,
                                   squeeze_base_boost,
                                   squeeze_max_boost,
                                   squeeze_layer_decay,
                                   squeeze_halo,
                                   squeeze_hotspot_guard);
            applyTopLayerBias(grouter_,
                              normalized_rudy,
                              hotspots,
                              min_routing_layer,
                              max_routing_layer,
                              squeeze_top_k,
                              squeeze_top_threshold,
                              squeeze_top_base,
                              squeeze_top_max,
                              squeeze_top_guard,
                              squeeze_top_halo);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    squeeze_halo,
                                    squeeze_hotspot_ratio,
                                    squeeze_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_squeeze);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data) {
    const float cool_threshold = std::clamp(
        0.42f + 0.10f * (0.65f - congestion_severity), 0.32f, 0.52f);
    const float cool_base_boost = std::clamp(
        0.05f + 0.12f * (0.70f - congestion_severity), 0.03f, 0.16f);
    const float cool_max_boost = std::clamp(
        1.07f + 0.10f * (0.60f - congestion_severity), 1.06f, 1.18f);
    const float cool_layer_decay
        = std::clamp(0.05f + 0.08f * hotspot_bias, 0.03f, 0.16f);
    const int cool_halo = hotspots.size() > 2 ? 2 : 1;
    const float cool_hotspot_guard = std::clamp(
        0.66f - 0.10f * hotspot_bias, 0.55f, 0.70f);
    const float cool_perturb = std::clamp(
        0.02f + 0.14f * congestion_severity, 0.015f, 0.24f);
    const float cool_critical = std::clamp(
        4.6f + 2.6f * (0.60f - congestion_severity), 3.8f, 9.0f);
    const int cool_seed = snapshot.seed + 59;
    const float cool_hotspot_ratio
        = std::clamp(0.995f - 0.05f * hotspot_bias, 0.96f, 0.995f);
    const float cool_hotspot_weight
        = std::clamp(0.06f + 0.12f * hotspot_bias, 0.05f, 0.16f);

    ScenarioDefinition wl_coolcorr;
    wl_coolcorr.name = "wl-coolcorr";
    wl_coolcorr.pre_init
        = [this,
           cool_perturb,
           cool_seed,
           cool_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(cool_perturb);
            grouter_->setPerturbationAmount(cool_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(cool_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(cool_critical);
            apply_wl_via_scale();
          };
    wl_coolcorr.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           cool_threshold,
           cool_base_boost,
           cool_max_boost,
           cool_layer_decay,
           cool_halo,
           cool_hotspot_guard,
           cool_hotspot_ratio,
           cool_hotspot_weight]() {
            applyCoolCapacityBoost(grouter_,
                                   normalized_rudy,
                                   hotspots,
                                   min_routing_layer,
                                   max_routing_layer,
                                   cool_threshold,
                                   cool_base_boost,
                                   cool_max_boost,
                                   cool_layer_decay,
                                   cool_halo,
                                   cool_hotspot_guard);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    cool_halo,
                                    cool_hotspot_ratio,
                                    cool_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_coolcorr);
  }

  if (baseline.metrics.overflow == 0 && has_congestion_data) {
    const float smooth_strength
        = std::clamp(0.32f + 0.48f * congestion_severity
                         + 0.20f * hotspot_bias,
                     0.30f,
                     0.85f);
    const float smooth_min_base
        = std::clamp(0.90f + 0.03f * (1.0f - congestion_severity),
                     0.90f,
                     0.95f);
    const float smooth_max_base
        = std::clamp(0.975f + 0.01f * (1.0f - congestion_severity),
                     smooth_min_base + 0.02f,
                     0.99f);
    const float smooth_slope = 2.0f + 1.4f * smooth_strength;
    const float smooth_midpoint
        = std::clamp(0.50f - 0.03f * smooth_strength, 0.44f, 0.52f);
    const float smooth_threshold
        = std::clamp(0.50f + 0.12f * congestion_severity, 0.48f, 0.64f);
    const float smooth_min_ratio
        = std::clamp(0.93f + 0.05f * (1.0f - congestion_severity),
                     0.92f,
                     0.985f);
    const float smooth_hotspot_push
        = std::clamp(0.07f + 0.06f * hotspot_bias, 0.06f, 0.16f);
    const int smooth_halo = hotspots.size() > 6 ? 2 : 1;
    const float smooth_cool_threshold
        = std::clamp(smooth_threshold * 0.70f, 0.36f, 0.52f);
    const float smooth_boost = std::clamp(
        0.10f + 0.08f * (0.55f - congestion_severity), 0.08f, 0.18f);
    const float smooth_boost_limit = std::clamp(
        1.09f + 0.06f * (0.55f - congestion_severity), 1.08f, 1.17f);
    const float smooth_layer_falloff
        = std::clamp(0.09f + 0.12f * hotspot_bias, 0.08f, 0.22f);
    const float smooth_perturb
        = std::clamp(0.02f + 0.12f * congestion_severity,
                     0.015f,
                     0.14f);
    const float smooth_critical = std::clamp(
        6.0f + 3.0f * (0.55f - congestion_severity), 5.5f, 10.5f);
    const int smooth_seed = snapshot.seed + 143;
    const float smooth_hotspot_ratio
        = std::clamp(0.99f - 0.05f * smooth_strength, 0.93f, 0.99f);
    const float smooth_hotspot_weight
        = std::clamp(0.10f + 0.18f * hotspot_bias, 0.08f, 0.22f);

    ScenarioDefinition wl_smooth;
    wl_smooth.name = "wl-smooth";
    wl_smooth.pre_init
        = [this,
           smooth_perturb,
           smooth_seed,
           smooth_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(smooth_perturb);
            grouter_->setPerturbationAmount(smooth_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(smooth_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(smooth_critical);
            apply_wl_via_scale();
          };
    wl_smooth.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           smooth_min_base,
           smooth_max_base,
           smooth_slope,
           smooth_midpoint,
           smooth_threshold,
           smooth_min_ratio,
           smooth_hotspot_push,
           smooth_halo,
           smooth_cool_threshold,
           smooth_boost,
           smooth_boost_limit,
           smooth_layer_falloff,
           smooth_hotspot_ratio,
           smooth_hotspot_weight]() {
            applySoftCapacityScaling(grouter_,
                                     normalized_rudy,
                                     min_routing_layer,
                                     max_routing_layer,
                                     smooth_min_base,
                                     smooth_max_base,
                                     smooth_slope,
                                     smooth_midpoint);
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 smooth_threshold,
                                 smooth_min_ratio,
                                 smooth_hotspot_push,
                                 smooth_halo,
                                 smooth_cool_threshold,
                                 smooth_boost,
                                 smooth_boost_limit,
                                 smooth_layer_falloff);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    smooth_halo,
                                    smooth_hotspot_ratio,
                                    smooth_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_smooth);
  }

  if (baseline.metrics.overflow == 0) {
    const float stability_threshold
        = std::clamp(0.52f + 0.10f * hotspot_bias, 0.50f, 0.62f);
    const float stability_min_ratio
        = std::clamp(0.93f + 0.04f * (1.0f - congestion_severity),
                     0.93f,
                     0.98f);
    const float stability_hotspot_push
        = std::clamp(0.12f + 0.08f * hotspot_bias, 0.10f, 0.20f);
    const int stability_halo = hotspots.size() > 5 ? 2 : 1;
    const float stability_cool_threshold
        = std::clamp(stability_threshold * 0.78f, 0.40f, 0.54f);
    const float stability_boost = std::clamp(
        0.12f + 0.10f * (0.55f - congestion_severity), 0.10f, 0.18f);
    const float stability_boost_limit = std::clamp(
        1.11f + 0.06f * (0.55f - congestion_severity), 1.09f, 1.18f);
    const float stability_layer_falloff
        = std::clamp(0.10f + 0.12f * hotspot_bias, 0.10f, 0.22f);
    const float stability_perturb
        = std::clamp(0.015f + 0.10f * congestion_severity, 0.01f, 0.14f);
    const float stability_critical = std::clamp(
        6.2f + 2.8f * (0.55f - congestion_severity), 5.2f, 9.5f);
    const int stability_seed = snapshot.seed + 277;
    const float stability_min_base = 0.90f;
    const float stability_max_base = 0.985f;
    const float stability_slope = 2.1f;
    const float stability_midpoint = 0.52f;
    const float stability_hotspot_ratio
        = std::clamp(0.985f - 0.03f * hotspot_bias, 0.96f, 0.99f);
    const float stability_hotspot_weight
        = std::clamp(0.10f + 0.12f * hotspot_bias, 0.10f, 0.20f);

    ScenarioDefinition wl_stability;
    wl_stability.name = "wl-stability";
    wl_stability.pre_init
        = [this,
           stability_perturb,
           stability_seed,
           stability_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(stability_perturb);
            grouter_->setPerturbationAmount(stability_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(stability_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(stability_critical);
            apply_wl_via_scale();
          };
    wl_stability.post_init
        = [this,
           &normalized_rudy,
           &hotspots,
           min_routing_layer,
           max_routing_layer,
           stability_min_base,
           stability_max_base,
           stability_slope,
           stability_midpoint,
           stability_threshold,
           stability_min_ratio,
           stability_hotspot_push,
           stability_halo,
           stability_cool_threshold,
           stability_boost,
           stability_boost_limit,
           stability_layer_falloff,
           stability_hotspot_ratio,
           stability_hotspot_weight]() {
            applySoftCapacityScaling(grouter_,
                                     normalized_rudy,
                                     min_routing_layer,
                                     max_routing_layer,
                                     stability_min_base,
                                     stability_max_base,
                                     stability_slope,
                                     stability_midpoint);
            applySelectiveRelief(grouter_,
                                 normalized_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 stability_threshold,
                                 stability_min_ratio,
                                 stability_hotspot_push,
                                 stability_halo,
                                 stability_cool_threshold,
                                 stability_boost,
                                 stability_boost_limit,
                                 stability_layer_falloff);
            if (!hotspots.empty()) {
              applyHotspotPenalties(grouter_,
                                    hotspots,
                                    min_routing_layer,
                                    max_routing_layer,
                                    stability_halo,
                                    stability_hotspot_ratio,
                                    stability_hotspot_weight);
            }
          };
    scenario_defs.push_back(wl_stability);
  }

  if (has_congestion_data && run_soft && !normalized_rudy.empty()
      && !light_congestion) {
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
        = [this,
           contour_perturb,
           contour_seed,
           contour_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(contour_perturb);
            grouter_->setPerturbationAmount(contour_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(contour_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(contour_critical);
            apply_wl_via_scale();
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

  if (has_congestion_data && !light_congestion) {
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
        = [this,
           corridor_perturb,
           corridor_seed,
           corridor_critical,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(corridor_perturb);
            grouter_->setPerturbationAmount(corridor_perturb > 0.0f ? 1 : 0);
            grouter_->setSeed(corridor_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(corridor_critical);
            apply_wl_via_scale();
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

  if (has_congestion_data && !light_congestion) {
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
        = [this, perturb_pct, seed, critical_pct, apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
            grouter_->setPerturbationAmount(perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(critical_pct);
            apply_wl_via_scale();
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

  if (has_congestion_data && run_aggressive_soft && !light_congestion) {
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
        = [this,
           selective_perturb_pct,
           selective_seed,
           selective_critical_pct,
           apply_wl_via_scale]() {
            grouter_->setCapacitiesPerturbationPercentage(
                selective_perturb_pct);
            grouter_->setPerturbationAmount(selective_perturb_pct > 0.0f ? 1 : 0);
            grouter_->setSeed(selective_seed);
            grouter_->setAllowCongestion(false);
            grouter_->fastroute_->setCriticalNetsPercentage(
                selective_critical_pct);
            apply_wl_via_scale();
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

  if (has_congestion_data && run_aggressive_soft && !light_congestion) {
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

  if (light_congestion && baseline.metrics.overflow == 0) {
    static const std::unordered_set<std::string> skip_names{
        "contour-lite",
        "cool-corridors",
        "wl-coolboost",
        "focused-soft",
        "wl-coolcorr"};
    scenario_defs.erase(
        std::remove_if(scenario_defs.begin(),
                       scenario_defs.end(),
                       [&](const ScenarioDefinition& def) {
                         return skip_names.find(def.name) != skip_names.end();
                       }),
        scenario_defs.end());
  }

  if (ultra_light) {
    static const std::unordered_set<std::string> ultra_skip{
        "wl-compact",
        "wl-squeeze",
        "wl-smooth",
        "wl-stability",
        "wl-coolcorr",
        "wl-coolboost",
        "contour-lite",
        "cool-corridors",
        "focused-soft",
        "soft-relief"};
    scenario_defs.erase(
        std::remove_if(scenario_defs.begin(),
                       scenario_defs.end(),
                       [&](const ScenarioDefinition& def) {
                         return ultra_skip.find(def.name) != ultra_skip.end();
                       }),
        scenario_defs.end());
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
    const double wl_primary = 0.00035;    // ~0.035% difference
    const double wl_tie = 0.00015;        // ~0.015% difference

    const double util_gap
        = lhs.metrics.max_utilization - rhs.metrics.max_utilization;
    const double reserve_gap
        = lhs.metrics.reserve_score - rhs.metrics.reserve_score;
    const double util_guard = 0.08;
    const double via_gap = static_cast<double>(lhs.metrics.via_count)
                           - static_cast<double>(rhs.metrics.via_count);
    const double via_rel
        = std::abs(via_gap)
          / std::max<double>(
              std::min(lhs.metrics.via_count, rhs.metrics.via_count), 1.0);

    // Prefer shorter wirelength while allowing a modest utilization cushion.
    if (wl_a != wl_b && wl_rel > wl_primary) {
      if (wl_a < wl_b
          && lhs.metrics.max_utilization
                 <= rhs.metrics.max_utilization + util_guard) {
        return true;
      }
      if (wl_b < wl_a
          && rhs.metrics.max_utilization
                 <= lhs.metrics.max_utilization + util_guard) {
        return false;
      }
      return wl_a < wl_b;
    }

    if (wl_rel < 0.0025 && std::abs(util_gap) > 0.025) {
      return util_gap < 0.0;
    }

    if (wl_rel < 0.0022 && std::abs(reserve_gap) > 0.05
        && std::abs(util_gap) < 0.05) {
      return reserve_gap > 0.0;
    }

    if (wl_rel < 0.0018) {
      const bool wl_close = wl_rel < 0.0010;
      const bool via_meaningful
          = via_rel > 0.0025 || std::abs(via_gap) > 80.0;
      const bool util_safe
          = lhs.metrics.max_utilization
            <= rhs.metrics.max_utilization + 0.02;
      if (via_meaningful && (wl_close || util_safe)) {
        return via_gap < 0.0;
      }
    }

    if (wl_rel > wl_tie) {
      return wl_a < wl_b;
    }

    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }

    if (std::abs(util_gap) > 1e-4) {
      return util_gap < 0.0;
    }

    auto relative_gap = [](double a, double b) {
      const double denom = std::max({std::abs(a), std::abs(b), 1e-9});
      return (a - b) / denom;
    };

    if (wl_rel < 0.0025) {
      const double stress_rel
          = relative_gap(lhs.metrics.stress_cost, rhs.metrics.stress_cost);
      if (std::abs(stress_rel) > 0.025) {
        return stress_rel < 0.0;
      }

      const double reserve_rel
          = relative_gap(lhs.metrics.reserve_score, rhs.metrics.reserve_score);
      if (std::abs(reserve_rel) > 0.025) {
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

  const size_t best_index
      = static_cast<size_t>(best_iter - scenario_results.begin());
  size_t wl_pref_index = best_index;
  const long target_overflow = scenario_results[best_index].metrics.overflow;
  const double wl_gain_threshold = 0.00012;
  const double util_soft_guard = 0.10;
  const long via_preference_limit
      = baseline.metrics.overflow == 0 && light_congestion ? 1400 : 800;

  for (size_t i = 0; i < scenario_results.size(); ++i) {
    const ScenarioResult& candidate = scenario_results[i];
    ScenarioResult& current = scenario_results[wl_pref_index];
    if (candidate.metrics.overflow != target_overflow) {
      continue;
    }
    if (candidate.metrics.max_utilization
        > current.metrics.max_utilization + util_soft_guard) {
      continue;
    }
    if (candidate.metrics.max_utilization > 1.08) {
      continue;
    }

    const double wl_diff = static_cast<double>(current.metrics.wirelength_dbu)
                           - static_cast<double>(candidate.metrics.wirelength_dbu);
    const double wl_rel = wl_diff
                          / std::max<double>(
                              static_cast<double>(current.metrics.wirelength_dbu),
                              1.0);
    if (wl_diff <= 0.0 || wl_rel <= wl_gain_threshold) {
      continue;
    }

    const long via_gap
        = candidate.metrics.via_count - current.metrics.via_count;
    if (via_gap > via_preference_limit) {
      continue;
    }

    wl_pref_index = i;
  }

  if (wl_pref_index != best_index) {
    ScenarioResult& chosen = scenario_results[wl_pref_index];
    const ScenarioResult& original_best = scenario_results[best_index];
    best_iter = scenario_results.begin()
                + static_cast<std::vector<ScenarioResult>::difference_type>(
                    wl_pref_index);
    logger_->info(GNR,
                  6009,
                  "NEWGR wirelength preference picked scenario {} over {} "
                  "(WL {:.0f} -> {:.0f} um, vias {} -> {}, max util {:.2f} -> {:.2f})",
                  chosen.name,
                  original_best.name,
                  original_best.metrics.wirelength_um,
                  chosen.metrics.wirelength_um,
                  original_best.metrics.via_count,
                  chosen.metrics.via_count,
                  original_best.metrics.max_utilization,
                  chosen.metrics.max_utilization);
  }

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

  const bool dense_hotspots = hotspots.size() > 4;
  const bool heavy_congestion = congestion_severity > 0.75f;
  const bool near_overflow = final_result.metrics.max_utilization > 0.92;
  const bool overflowing = final_result.metrics.overflow > 0;
  const bool severe_hotspots
      = congestion_severity > 0.68f && hotspots.size() > 2;
  const bool should_patch
      = overflowing || near_overflow
        || (heavy_congestion && dense_hotspots)
        || (severe_hotspots && final_result.metrics.max_utilization > 0.88);

  if (should_patch) {
    PatchSummary patch_summary
        = applyHotspotPatches(final_result.routes,
                              hotspots,
                              grouter_->grid(),
                              min_routing_layer,
                              max_routing_layer);
    if (patch_summary.segments_added > 0) {
      RouteMetrics patched_metrics = compute_metrics(final_result.routes);
      logger_->info(
          GNR,
          6010,
          "NEWGR applied {} hotspot patches on {} nets. "
          "Patched wirelength {:.0f} um, vias {}, overflow {}, max util {:.2f}",
          patch_summary.segments_added,
          patch_summary.nets_touched,
          patched_metrics.wirelength_um,
          patched_metrics.via_count,
          patched_metrics.overflow,
          patched_metrics.max_utilization);
    }
  }

  return std::move(final_result.routes);
}

}  // namespace grt
