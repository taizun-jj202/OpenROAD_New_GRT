#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
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
  double wirelength_um = 0.0;
  double score = 0.0;
};

struct ScenarioResult
{
  std::string name;
  RouteMetrics metrics;
  NetRouteMap routes;
  int overflow = 0;
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
using GridMatrix = std::vector<std::vector<float>>;

struct PlanarEdgeUsage
{
  GridMatrix horizontal;
  GridMatrix vertical;
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

PlanarEdgeUsage computeNormalizedBackboneUsage(GlobalRouter* grouter,
                                               const NetRouteMap& routes)
{
  PlanarEdgeUsage usage;
  if (grouter == nullptr || grouter->grid() == nullptr) {
    return usage;
  }

  Grid* grid = grouter->grid();
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return usage;
  }

  usage.horizontal.assign(
      std::max(x_grids - 1, 0), std::vector<float>(y_grids, 0.0f));
  usage.vertical.assign(
      x_grids, std::vector<float>(std::max(y_grids - 1, 0), 0.0f));

  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  float max_usage = 0.0f;

  auto grid_x = [&](int x_dbu) {
    return std::clamp((x_dbu - x_min) / tile_size, 0, std::max(x_grids - 1, 0));
  };
  auto grid_y = [&](int y_dbu) {
    return std::clamp((y_dbu - y_min) / tile_size, 0, std::max(y_grids - 1, 0));
  };

  for (const auto& [db_net, segments] : routes) {
    static_cast<void>(db_net);
    for (const GSegment& segment : segments) {
      if (segment.isVia()) {
        continue;
      }
      const int gx0 = grid_x(segment.init_x);
      const int gy0 = grid_y(segment.init_y);
      const int gx1 = grid_x(segment.final_x);
      const int gy1 = grid_y(segment.final_y);

      if (gy0 == gy1 && gx0 != gx1 && !usage.horizontal.empty()) {
        const int y = std::clamp(gy0, 0, std::max(y_grids - 1, 0));
        const int start = std::min(gx0, gx1);
        const int end = std::max(gx0, gx1);
        for (int x = start; x < end; ++x) {
          if (x < 0 || x >= static_cast<int>(usage.horizontal.size())) {
            continue;
          }
          usage.horizontal[x][y] += 1.0f;
          max_usage = std::max(max_usage, usage.horizontal[x][y]);
        }
      } else if (gx0 == gx1 && gy0 != gy1 && !usage.vertical.empty()) {
        const int x = std::clamp(gx0, 0, std::max(x_grids - 1, 0));
        const int start = std::min(gy0, gy1);
        const int end = std::max(gy0, gy1);
        for (int y = start; y < end; ++y) {
          if (y < 0 || y >= static_cast<int>(usage.vertical[x].size())) {
            continue;
          }
          usage.vertical[x][y] += 1.0f;
          max_usage = std::max(max_usage, usage.vertical[x][y]);
        }
      }
    }
  }

  if (max_usage <= std::numeric_limits<float>::epsilon()) {
    return usage;
  }

  for (auto& col : usage.horizontal) {
    for (float& val : col) {
      val = std::clamp(val / max_usage, 0.0f, 1.0f);
    }
  }
  for (auto& col : usage.vertical) {
    for (float& val : col) {
      val = std::clamp(val / max_usage, 0.0f, 1.0f);
    }
  }
  return usage;
}

void adjustEdgeCapacity(GlobalRouter* grouter,
                        int x1,
                        int y1,
                        int x2,
                        int y2,
                        int layer,
                        float ratio)
{
  ratio = std::clamp(ratio, 0.05f, 1.45f);
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

void applyAggressiveCapacityField(GlobalRouter* grouter,
                                  const RudyGrid& normalized_rudy,
                                  const std::vector<Hotspot>& hotspots,
                                  int min_layer,
                                  int max_layer,
                                  float high_rudy_threshold,
                                  float low_rudy_threshold,
                                  float high_rudy_ratio,
                                  float low_rudy_ratio,
                                  float orientation_bias,
                                  bool horizontal_preference)
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
  if (usable_x < 2 || usable_y < 2) {
    return;
  }

  high_rudy_threshold = std::clamp(high_rudy_threshold, 0.05f, 0.95f);
  low_rudy_threshold = std::clamp(low_rudy_threshold, 0.01f, 0.80f);
  if (low_rudy_threshold > high_rudy_threshold) {
    std::swap(low_rudy_threshold, high_rudy_threshold);
  }

  const int layer_span = std::max(max_layer - min_layer, 1);
  orientation_bias = std::clamp(orientation_bias, 0.0f, 0.8f);

  std::vector<float> hotspot_energy(usable_x * usable_y, 0.0f);
  auto idx = [usable_x](int x, int y) { return y * usable_x + x; };
  for (const Hotspot& hotspot : hotspots) {
    const int halo = 3;
    for (int dx = -halo; dx <= halo; ++dx) {
      for (int dy = -halo; dy <= halo; ++dy) {
        const int gx = hotspot.gx + dx;
        const int gy = hotspot.gy + dy;
        if (gx < 0 || gy < 0 || gx >= usable_x || gy >= usable_y) {
          continue;
        }
        const float distance = static_cast<float>(std::abs(dx) + std::abs(dy));
        const float decay = 1.0f / (1.0f + distance);
        hotspot_energy[idx(gx, gy)] += hotspot.severity * decay;
      }
    }
  }

  float max_hotspot_energy = 0.0f;
  for (const float energy : hotspot_energy) {
    max_hotspot_energy = std::max(max_hotspot_energy, energy);
  }

  const auto get_rudy = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };
  const auto get_hotspot = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y
        || max_hotspot_energy <= std::numeric_limits<float>::epsilon()) {
      return 0.0f;
    }
    return hotspot_energy[idx(x, y)] / max_hotspot_energy;
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float layer_relax = 1.0f + 0.18f * layer_factor;
    const float orient_boost = 1.0f + orientation_bias * (1.0f - 0.35f * layer_factor);
    const float orient_shrink = std::max(0.2f, 1.0f - orientation_bias);

    auto edge_ratio = [&](float local_rudy, float local_hotspot, bool is_horizontal) {
      float ratio = 1.0f;
      if (local_rudy >= high_rudy_threshold) {
        ratio = high_rudy_ratio;
      } else if (local_rudy <= low_rudy_threshold) {
        ratio = low_rudy_ratio;
      } else {
        const float blend
            = (local_rudy - low_rudy_threshold)
              / std::max(high_rudy_threshold - low_rudy_threshold, 0.01f);
        ratio = low_rudy_ratio + (high_rudy_ratio - low_rudy_ratio) * blend;
      }

      // Carve out hard "no-fly" channels around persistent hotspots.
      ratio *= std::clamp(1.0f - 0.65f * local_hotspot, 0.15f, 1.0f);
      ratio *= layer_relax;

      if (horizontal_preference) {
        ratio *= is_horizontal ? orient_boost : orient_shrink;
      } else {
        ratio *= is_horizontal ? orient_shrink : orient_boost;
      }
      return std::clamp(ratio, 0.08f, 1.45f);
    };

    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float local_rudy = 0.5f * (get_rudy(x, y) + get_rudy(x + 1, y));
        const float local_hotspot
            = 0.5f * (get_hotspot(x, y) + get_hotspot(x + 1, y));
        adjustEdgeCapacity(grouter,
                           x,
                           y,
                           x + 1,
                           y,
                           layer,
                           edge_ratio(local_rudy, local_hotspot, true));
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float local_rudy = 0.5f * (get_rudy(x, y) + get_rudy(x, y + 1));
        const float local_hotspot
            = 0.5f * (get_hotspot(x, y) + get_hotspot(x, y + 1));
        adjustEdgeCapacity(grouter,
                           x,
                           y,
                           x,
                           y + 1,
                           layer,
                           edge_ratio(local_rudy, local_hotspot, false));
      }
    }
  }
}

void applyBackboneCapacityReinforcement(GlobalRouter* grouter,
                                        const RudyGrid& normalized_rudy,
                                        const PlanarEdgeUsage& usage,
                                        int min_layer,
                                        int max_layer,
                                        float low_usage_ratio,
                                        float high_usage_ratio,
                                        float usage_gamma)
{
  Grid* grid = grouter->grid();
  if (grid == nullptr) {
    return;
  }
  if (usage.horizontal.empty() && usage.vertical.empty()) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return;
  }

  low_usage_ratio = std::clamp(low_usage_ratio, 0.10f, 1.20f);
  high_usage_ratio = std::clamp(high_usage_ratio, low_usage_ratio, 1.60f);
  usage_gamma = std::clamp(usage_gamma, 0.20f, 2.20f);

  const int layer_span = std::max(max_layer - min_layer, 1);
  const auto sample_rudy = [&](int x, int y) {
    if (normalized_rudy.empty()) {
      return 0.0f;
    }
    if (x < 0 || y < 0 || x >= static_cast<int>(normalized_rudy.size())
        || y >= static_cast<int>(normalized_rudy.front().size())) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  auto edge_ratio = [&](float normalized_usage, float local_rudy) {
    normalized_usage = std::clamp(normalized_usage, 0.0f, 1.0f);
    local_rudy = std::clamp(local_rudy, 0.0f, 1.0f);
    float ratio = low_usage_ratio
                  + (high_usage_ratio - low_usage_ratio)
                        * std::pow(normalized_usage, usage_gamma);
    // Reserve resources in high-RUDY regions and open low-RUDY channels.
    ratio *= std::clamp(1.18f - 0.72f * local_rudy, 0.30f, 1.35f);
    return std::clamp(ratio, 0.10f, 1.60f);
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float layer_relax = 0.90f + 0.22f * layer_factor;

    for (int x = 0; x < x_grids - 1; ++x) {
      for (int y = 0; y < y_grids; ++y) {
        if (x >= static_cast<int>(usage.horizontal.size())
            || y >= static_cast<int>(usage.horizontal[x].size())) {
          continue;
        }
        const float local_usage = usage.horizontal[x][y];
        const float local_rudy
            = 0.5f * (sample_rudy(x, y) + sample_rudy(x + 1, y));
        const float ratio = std::clamp(
            edge_ratio(local_usage, local_rudy) * layer_relax, 0.10f, 1.60f);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int x = 0; x < x_grids; ++x) {
      for (int y = 0; y < y_grids - 1; ++y) {
        if (x >= static_cast<int>(usage.vertical.size())
            || y >= static_cast<int>(usage.vertical[x].size())) {
          continue;
        }
        const float local_usage = usage.vertical[x][y];
        const float local_rudy
            = 0.5f * (sample_rudy(x, y) + sample_rudy(x, y + 1));
        const float ratio = std::clamp(
            edge_ratio(local_usage, local_rudy) * layer_relax, 0.10f, 1.60f);
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
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
    const float exponent = -slope * (normalized - midpoint);
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

void applyLayerPolarityField(GlobalRouter* grouter,
                             const RudyGrid& normalized_rudy,
                             int min_layer,
                             int max_layer,
                             float favored_ratio,
                             float suppressed_ratio,
                             float transition_rudy)
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
  if (usable_x < 2 || usable_y < 2) {
    return;
  }

  favored_ratio = std::clamp(favored_ratio, 1.0f, 1.60f);
  suppressed_ratio = std::clamp(suppressed_ratio, 0.08f, 0.98f);
  transition_rudy = std::clamp(transition_rudy, 0.10f, 0.90f);

  const auto sample = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const bool favor_horizontal = ((layer - min_layer) % 2) == 0;
    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float local_rudy = 0.5f * (sample(x, y) + sample(x + 1, y));
        const float polarity
            = std::clamp((local_rudy - transition_rudy) / 0.30f, -1.0f, 1.0f);
        float ratio = favor_horizontal
                          ? favored_ratio - 0.15f * polarity
                          : suppressed_ratio + 0.30f * local_rudy;
        ratio = std::clamp(ratio, 0.08f, 1.60f);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }
    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float local_rudy = 0.5f * (sample(x, y) + sample(x, y + 1));
        const float polarity
            = std::clamp((local_rudy - transition_rudy) / 0.30f, -1.0f, 1.0f);
        float ratio = favor_horizontal
                          ? suppressed_ratio + 0.30f * local_rudy
                          : favored_ratio - 0.15f * polarity;
        ratio = std::clamp(ratio, 0.08f, 1.60f);
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

float sampleRudyAt(const RudyGrid& normalized_rudy, int gx, int gy)
{
  if (normalized_rudy.empty()) {
    return 0.0f;
  }
  if (gx < 0 || gy < 0 || gx >= static_cast<int>(normalized_rudy.size())
      || gy >= static_cast<int>(normalized_rudy.front().size())) {
    return 0.0f;
  }
  return normalized_rudy[gx][gy];
}

void appendSegment(std::vector<GSegment>& route,
                   int x0,
                   int y0,
                   int l0,
                   int x1,
                   int y1,
                   int l1)
{
  if (x0 == x1 && y0 == y1 && l0 == l1) {
    return;
  }
  route.emplace_back(x0, y0, l0, x1, y1, l1);
}

[[maybe_unused]] void applyBraidedDetourWeave(GlobalRouter* grouter,
                                              NetRouteMap& routes,
                                              const RudyGrid& normalized_rudy)
{
  if (grouter == nullptr || grouter->grid() == nullptr || routes.empty()) {
    return;
  }

  Grid* grid = grouter->grid();
  const int tile = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int x_max = grid->getXMax();
  const int y_max = grid->getYMax();
  const int min_segment_len = tile * 6;
  const int detour_step = tile * 2;
  const int x_tiles = normalized_rudy.empty() ? 0 : normalized_rudy.size();
  const int y_tiles
      = normalized_rudy.empty() ? 0 : normalized_rudy.front().size();

  auto to_grid_x = [&](int x) {
    if (x_tiles <= 0) {
      return 0;
    }
    return std::clamp((x - x_min) / tile, 0, x_tiles - 1);
  };
  auto to_grid_y = [&](int y) {
    if (y_tiles <= 0) {
      return 0;
    }
    return std::clamp((y - y_min) / tile, 0, y_tiles - 1);
  };
  auto choose_detour_coord = [&](int base_coord,
                                 int fixed_coord,
                                 bool horizontal,
                                 std::uint64_t key) {
    const int lower = horizontal ? y_min + tile : x_min + tile;
    const int upper = horizontal ? y_max - tile : x_max - tile;
    if (lower >= upper) {
      return base_coord;
    }

    const int pos = std::clamp(base_coord + detour_step, lower, upper);
    const int neg = std::clamp(base_coord - detour_step, lower, upper);
    if (pos == base_coord && neg == base_coord) {
      return base_coord;
    }

    const int gx_pos = horizontal ? to_grid_x(fixed_coord) : to_grid_x(pos);
    const int gy_pos = horizontal ? to_grid_y(pos) : to_grid_y(fixed_coord);
    const int gx_neg = horizontal ? to_grid_x(fixed_coord) : to_grid_x(neg);
    const int gy_neg = horizontal ? to_grid_y(neg) : to_grid_y(fixed_coord);

    const float pos_rudy = sampleRudyAt(normalized_rudy, gx_pos, gy_pos);
    const float neg_rudy = sampleRudyAt(normalized_rudy, gx_neg, gy_neg);

    if (std::fabs(pos_rudy - neg_rudy) < 0.04f) {
      const bool prefer_pos = (key & 1ULL) == 0ULL;
      return prefer_pos ? pos : neg;
    }
    return (pos_rudy < neg_rudy) ? pos : neg;
  };

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    std::uint64_t net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if ((net_key % 100ULL) >= 60ULL) {
      continue;
    }

    std::vector<GSegment> detoured;
    detoured.reserve(route.size() * 3);
    for (size_t seg_idx = 0; seg_idx < route.size(); ++seg_idx) {
      const GSegment& segment = route[seg_idx];
      const bool is_via = segment.isVia();
      const bool horizontal
          = segment.init_layer == segment.final_layer
            && segment.init_y == segment.final_y
            && segment.init_x != segment.final_x;
      const bool vertical = segment.init_layer == segment.final_layer
                            && segment.init_x == segment.final_x
                            && segment.init_y != segment.final_y;
      const long length = segment.length();
      const std::uint64_t key = net_key + static_cast<std::uint64_t>(seg_idx) * 131ULL;

      const bool detour_candidate
          = !is_via && (horizontal || vertical) && length >= min_segment_len
            && ((key % 4ULL) == 0ULL);
      if (!detour_candidate) {
        detoured.push_back(segment);
        continue;
      }

      if (horizontal) {
        int mid_x = segment.init_x + (segment.final_x - segment.init_x) / 2;
        if (mid_x == segment.init_x || mid_x == segment.final_x) {
          mid_x = segment.init_x + (segment.final_x - segment.init_x) / 3;
        }
        const int detour_y
            = choose_detour_coord(segment.init_y, mid_x, true, key);
        if (detour_y == segment.init_y) {
          detoured.push_back(segment);
          continue;
        }

        appendSegment(detoured,
                      segment.init_x,
                      segment.init_y,
                      segment.init_layer,
                      mid_x,
                      segment.init_y,
                      segment.init_layer);
        appendSegment(detoured,
                      mid_x,
                      segment.init_y,
                      segment.init_layer,
                      mid_x,
                      detour_y,
                      segment.init_layer);
        appendSegment(detoured,
                      mid_x,
                      detour_y,
                      segment.init_layer,
                      segment.final_x,
                      detour_y,
                      segment.init_layer);
        appendSegment(detoured,
                      segment.final_x,
                      detour_y,
                      segment.init_layer,
                      segment.final_x,
                      segment.final_y,
                      segment.final_layer);
        continue;
      }

      int mid_y = segment.init_y + (segment.final_y - segment.init_y) / 2;
      if (mid_y == segment.init_y || mid_y == segment.final_y) {
        mid_y = segment.init_y + (segment.final_y - segment.init_y) / 3;
      }
      const int detour_x = choose_detour_coord(segment.init_x, mid_y, false, key);
      if (detour_x == segment.init_x) {
        detoured.push_back(segment);
        continue;
      }

      appendSegment(detoured,
                    segment.init_x,
                    segment.init_y,
                    segment.init_layer,
                    segment.init_x,
                    mid_y,
                    segment.init_layer);
      appendSegment(detoured,
                    segment.init_x,
                    mid_y,
                    segment.init_layer,
                    detour_x,
                    mid_y,
                    segment.init_layer);
      appendSegment(detoured,
                    detour_x,
                    mid_y,
                    segment.init_layer,
                    detour_x,
                    segment.final_y,
                    segment.init_layer);
      appendSegment(detoured,
                    detour_x,
                    segment.final_y,
                    segment.init_layer,
                    segment.final_x,
                    segment.final_y,
                    segment.final_layer);
    }
    if (!detoured.empty()) {
      route.swap(detoured);
    }
  }
}

[[maybe_unused]] void applyLayerHoppingDetours(GlobalRouter* grouter,
                                               NetRouteMap& routes,
                                               const RudyGrid& normalized_rudy,
                                               int min_layer,
                                               int max_layer)
{
  if (grouter == nullptr || grouter->grid() == nullptr || routes.empty()) {
    return;
  }

  Grid* grid = grouter->grid();
  const int tile = std::max(grid->getTileSize(), 1);
  const int min_segment_len = tile * 3;
  const int detour_step = tile * 2;
  const int x_min = grid->getXMin();
  const int x_max = grid->getXMax();
  const int y_min = grid->getYMin();
  const int y_max = grid->getYMax();
  const int x_tiles = normalized_rudy.empty() ? 0 : normalized_rudy.size();
  const int y_tiles
      = normalized_rudy.empty() ? 0 : normalized_rudy.front().size();

  auto to_grid_x = [&](int x) {
    if (x_tiles <= 0) {
      return 0;
    }
    return std::clamp((x - x_min) / tile, 0, x_tiles - 1);
  };
  auto to_grid_y = [&](int y) {
    if (y_tiles <= 0) {
      return 0;
    }
    return std::clamp((y - y_min) / tile, 0, y_tiles - 1);
  };

  auto choose_detour = [&](int straight_coord,
                           int fixed_coord,
                           bool horizontal,
                           std::uint64_t key) {
    const int lower = horizontal ? y_min + tile : x_min + tile;
    const int upper = horizontal ? y_max - tile : x_max - tile;
    if (lower >= upper) {
      return straight_coord;
    }

    const int pos = std::clamp(straight_coord + detour_step, lower, upper);
    const int neg = std::clamp(straight_coord - detour_step, lower, upper);
    if (pos == straight_coord && neg == straight_coord) {
      return straight_coord;
    }

    const int gx_pos = horizontal ? to_grid_x(fixed_coord) : to_grid_x(pos);
    const int gy_pos = horizontal ? to_grid_y(pos) : to_grid_y(fixed_coord);
    const int gx_neg = horizontal ? to_grid_x(fixed_coord) : to_grid_x(neg);
    const int gy_neg = horizontal ? to_grid_y(neg) : to_grid_y(fixed_coord);
    const float pos_rudy = sampleRudyAt(normalized_rudy, gx_pos, gy_pos);
    const float neg_rudy = sampleRudyAt(normalized_rudy, gx_neg, gy_neg);

    if (std::fabs(pos_rudy - neg_rudy) < 0.03f) {
      return ((key >> 1U) & 1ULL) == 0ULL ? pos : neg;
    }
    return pos_rudy < neg_rudy ? pos : neg;
  };

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const std::uint64_t net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    std::vector<GSegment> rewritten;
    rewritten.reserve(route.size() * 5);

    for (size_t seg_idx = 0; seg_idx < route.size(); ++seg_idx) {
      const GSegment& segment = route[seg_idx];
      const bool planar = !segment.isVia()
                          && segment.init_layer == segment.final_layer;
      const bool horizontal = planar && segment.init_y == segment.final_y
                              && segment.init_x != segment.final_x;
      const bool vertical = planar && segment.init_x == segment.final_x
                            && segment.init_y != segment.final_y;
      const long length = segment.length();
      const std::uint64_t key = net_key + static_cast<std::uint64_t>(seg_idx) * 97ULL;
      const bool rewrite_candidate
          = (horizontal || vertical) && length >= min_segment_len
            && (key % 3ULL != 0ULL);
      if (!rewrite_candidate) {
        rewritten.push_back(segment);
        continue;
      }

      const int base_layer = segment.init_layer;
      int hop_layer = ((key & 1ULL) == 0ULL) ? base_layer + 1 : base_layer - 1;
      if (hop_layer < min_layer || hop_layer > max_layer) {
        hop_layer = (hop_layer < min_layer) ? base_layer + 1 : base_layer - 1;
      }
      if (hop_layer < min_layer || hop_layer > max_layer || hop_layer == base_layer) {
        rewritten.push_back(segment);
        continue;
      }

      if (horizontal) {
        const int x0 = segment.init_x;
        const int x1 = segment.final_x;
        const int y = segment.init_y;
        const int mid_x = x0 + (x1 - x0) / 2;
        if (mid_x == x0 || mid_x == x1) {
          rewritten.push_back(segment);
          continue;
        }

        const int detour_y = choose_detour(y, mid_x, true, key);
        if (detour_y == y) {
          rewritten.push_back(segment);
          continue;
        }

        appendSegment(rewritten, x0, y, base_layer, mid_x, y, base_layer);
        appendSegment(rewritten, mid_x, y, base_layer, mid_x, y, hop_layer);
        appendSegment(rewritten, mid_x, y, hop_layer, mid_x, detour_y, hop_layer);
        appendSegment(rewritten, mid_x, detour_y, hop_layer, x1, detour_y, hop_layer);
        appendSegment(rewritten, x1, detour_y, hop_layer, x1, y, hop_layer);
        appendSegment(rewritten, x1, y, hop_layer, x1, y, base_layer);
        continue;
      }

      const int x = segment.init_x;
      const int y0 = segment.init_y;
      const int y1 = segment.final_y;
      const int mid_y = y0 + (y1 - y0) / 2;
      if (mid_y == y0 || mid_y == y1) {
        rewritten.push_back(segment);
        continue;
      }

      const int detour_x = choose_detour(x, mid_y, false, key);
      if (detour_x == x) {
        rewritten.push_back(segment);
        continue;
      }

      appendSegment(rewritten, x, y0, base_layer, x, mid_y, base_layer);
      appendSegment(rewritten, x, mid_y, base_layer, x, mid_y, hop_layer);
      appendSegment(rewritten, x, mid_y, hop_layer, detour_x, mid_y, hop_layer);
      appendSegment(rewritten, detour_x, mid_y, hop_layer, detour_x, y1, hop_layer);
      appendSegment(rewritten, detour_x, y1, hop_layer, x, y1, hop_layer);
      appendSegment(rewritten, x, y1, hop_layer, x, y1, base_layer);
    }

    if (!rewritten.empty()) {
      route.swap(rewritten);
    }
  }
}

std::vector<Hotspot> extractTopRudyHotspots(const RudyGrid& normalized_rudy,
                                            int max_hotspots)
{
  std::vector<Hotspot> hotspots;
  if (normalized_rudy.empty() || normalized_rudy.front().empty()
      || max_hotspots <= 0) {
    return hotspots;
  }

  struct RudyPoint
  {
    int gx;
    int gy;
    float rudy;
  };

  std::vector<RudyPoint> candidates;
  candidates.reserve(normalized_rudy.size() * normalized_rudy.front().size());
  for (int gx = 0; gx < static_cast<int>(normalized_rudy.size()); ++gx) {
    for (int gy = 0; gy < static_cast<int>(normalized_rudy[gx].size()); ++gy) {
      const float rudy = normalized_rudy[gx][gy];
      if (rudy < 0.65f) {
        continue;
      }
      candidates.push_back(RudyPoint{gx, gy, rudy});
    }
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const RudyPoint& lhs, const RudyPoint& rhs) {
              return lhs.rudy > rhs.rudy;
            });

  const int picked = std::min(max_hotspots, static_cast<int>(candidates.size()));
  hotspots.reserve(picked);
  for (int i = 0; i < picked; ++i) {
    const RudyPoint& point = candidates[i];
    Hotspot hotspot;
    hotspot.gx = point.gx;
    hotspot.gy = point.gy;
    hotspot.severity = 1.0f + 1.6f * point.rudy;
    hotspot.affect_horizontal = true;
    hotspot.affect_vertical = true;
    hotspots.push_back(hotspot);
  }

  return hotspots;
}

bool isPlanarHorizontal(const GSegment& segment)
{
  return !segment.isVia() && segment.init_layer == segment.final_layer
         && segment.init_y == segment.final_y && segment.init_x != segment.final_x;
}

bool isPlanarVertical(const GSegment& segment)
{
  return !segment.isVia() && segment.init_layer == segment.final_layer
         && segment.init_x == segment.final_x && segment.init_y != segment.final_y;
}

bool samePoint(const GSegment& lhs, bool lhs_final, const GSegment& rhs, bool rhs_init)
{
  const int lhs_x = lhs_final ? lhs.final_x : lhs.init_x;
  const int lhs_y = lhs_final ? lhs.final_y : lhs.init_y;
  const int rhs_x = rhs_init ? rhs.init_x : rhs.final_x;
  const int rhs_y = rhs_init ? rhs.init_y : rhs.final_y;
  return lhs_x == rhs_x && lhs_y == rhs_y;
}

void appendCompressedSegment(std::vector<GSegment>& compressed,
                             const GSegment& segment)
{
  const bool degenerate_planar
      = !segment.isVia()
        && segment.init_x == segment.final_x
        && segment.init_y == segment.final_y;
  if (degenerate_planar) {
    return;
  }

  if (compressed.empty()) {
    compressed.push_back(segment);
    return;
  }

  GSegment& prev = compressed.back();
  const bool reverse_cancel = prev.init_x == segment.final_x
                              && prev.init_y == segment.final_y
                              && prev.final_x == segment.init_x
                              && prev.final_y == segment.init_y
                              && prev.init_layer == segment.final_layer
                              && prev.final_layer == segment.init_layer;
  if (reverse_cancel) {
    compressed.pop_back();
    return;
  }

  const bool contiguous = prev.final_x == segment.init_x
                          && prev.final_y == segment.init_y
                          && prev.final_layer == segment.init_layer;
  if (contiguous) {
    const bool same_horizontal = isPlanarHorizontal(prev)
                                 && isPlanarHorizontal(segment)
                                 && prev.init_y == segment.init_y;
    const bool same_vertical = isPlanarVertical(prev)
                               && isPlanarVertical(segment)
                               && prev.init_x == segment.init_x;
    if (same_horizontal || same_vertical) {
      prev.final_x = segment.final_x;
      prev.final_y = segment.final_y;
      prev.final_layer = segment.final_layer;
      return;
    }
  }

  compressed.push_back(segment);
}

bool compressJogTail(std::vector<GSegment>& route, int jog_limit)
{
  if (route.size() < 3) {
    return false;
  }

  const int n = route.size();
  const GSegment& first = route[n - 3];
  const GSegment& middle = route[n - 2];
  const GSegment& last = route[n - 1];

  if (!samePoint(first, true, middle, true)
      || !samePoint(middle, false, last, true)) {
    return false;
  }
  if (first.init_layer != first.final_layer
      || middle.init_layer != middle.final_layer
      || last.init_layer != last.final_layer) {
    return false;
  }
  if (first.final_layer != middle.init_layer
      || middle.final_layer != last.init_layer) {
    return false;
  }

  if (isPlanarHorizontal(first) && isPlanarVertical(middle)
      && isPlanarHorizontal(last) && first.init_y == last.final_y
      && std::abs(middle.final_y - middle.init_y) <= jog_limit) {
    GSegment merged(first.init_x,
                    first.init_y,
                    first.init_layer,
                    last.final_x,
                    last.final_y,
                    last.final_layer);
    route.resize(n - 3);
    appendCompressedSegment(route, merged);
    return true;
  }

  if (isPlanarVertical(first) && isPlanarHorizontal(middle)
      && isPlanarVertical(last) && first.init_x == last.final_x
      && std::abs(middle.final_x - middle.init_x) <= jog_limit) {
    GSegment merged(first.init_x,
                    first.init_y,
                    first.init_layer,
                    last.final_x,
                    last.final_y,
                    last.final_layer);
    route.resize(n - 3);
    appendCompressedSegment(route, merged);
    return true;
  }
  return false;
}

void applyGuideCompression(NetRouteMap& routes, int jog_limit)
{
  jog_limit = std::max(jog_limit, 1);
  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty()) {
      continue;
    }

    std::vector<GSegment> compressed;
    compressed.reserve(route.size());
    for (const GSegment& segment : route) {
      appendCompressedSegment(compressed, segment);
      while (compressJogTail(compressed, jog_limit)) {
      }
    }

    std::vector<GSegment> second_pass;
    second_pass.reserve(compressed.size());
    for (const GSegment& segment : compressed) {
      appendCompressedSegment(second_pass, segment);
    }

    if (!second_pass.empty()) {
      route.swap(second_pass);
    }
  }
}

bool samePointAndLayer(const GSegment& lhs,
                       bool lhs_final,
                       const GSegment& rhs,
                       bool rhs_init)
{
  const int lhs_layer = lhs_final ? lhs.final_layer : lhs.init_layer;
  const int rhs_layer = rhs_init ? rhs.init_layer : rhs.final_layer;
  return samePoint(lhs, lhs_final, rhs, rhs_init) && lhs_layer == rhs_layer;
}

void applyViaExcursionCollapse(NetRouteMap& routes, int max_planar_excursion)
{
  max_planar_excursion = std::max(max_planar_excursion, 1);
  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.size() < 3) {
      continue;
    }

    std::vector<GSegment> simplified;
    simplified.reserve(route.size());

    for (size_t idx = 0; idx < route.size();) {
      if (idx + 2 < route.size()) {
        const GSegment& first = route[idx];
        const GSegment& middle = route[idx + 1];
        const GSegment& last = route[idx + 2];

        const bool first_via = first.isVia() && first.init_x == first.final_x
                               && first.init_y == first.final_y;
        const bool last_via = last.isVia() && last.init_x == last.final_x
                              && last.init_y == last.final_y;
        const bool middle_planar
            = !middle.isVia() && middle.init_layer == middle.final_layer;
        const long middle_len = std::abs(middle.final_x - middle.init_x)
                                + std::abs(middle.final_y - middle.init_y);

        const bool contiguous = samePointAndLayer(first, true, middle, true)
                                && samePointAndLayer(middle, false, last, true);
        const bool returns_to_origin_layer
            = first.init_layer == last.final_layer
              && first.final_layer == middle.init_layer
              && middle.final_layer == last.init_layer
              && first.final_layer == last.init_layer;

        if (first_via && middle_planar && last_via && contiguous
            && returns_to_origin_layer
            && middle_len <= max_planar_excursion) {
          GSegment collapsed(first.init_x,
                             first.init_y,
                             first.init_layer,
                             last.final_x,
                             last.final_y,
                             last.final_layer);
          appendCompressedSegment(simplified, collapsed);
          idx += 3;
          continue;
        }
      }

      appendCompressedSegment(simplified, route[idx]);
      ++idx;
    }

    std::vector<GSegment> second_pass;
    second_pass.reserve(simplified.size());
    for (const GSegment& segment : simplified) {
      appendCompressedSegment(second_pass, segment);
    }

    if (!second_pass.empty()) {
      route.swap(second_pass);
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

    const double via_weight
        = static_cast<double>(std::max(grouter_->grid_->getTileSize(), 1))
          * 3.0;
    metrics.score = static_cast<double>(metrics.wirelength_dbu)
                    + via_weight * static_cast<double>(metrics.via_count);
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
    const int overflow
        = grouter_->fastroute() != nullptr ? grouter_->fastroute()->totalOverflow()
                                           : 0;
    RouteMetrics metrics = compute_metrics(routes);
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow {}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  overflow);
    return ScenarioResult{name, metrics, std::move(routes), overflow};
  };

  auto compute_current_rudy = [&]() {
    Rudy* rudy = grouter_->getRudy();
    if (rudy == nullptr) {
      return RudyGrid{};
    }
    rudy->calculateRudy();
    return computeNormalizedRudyGrid(rudy);
  };

  RouterSnapshot snapshot = capture_snapshot();
  ScenarioResult baseline = run_existing_state("baseline", nets);
  if (baseline.routes.empty()) {
    restore_snapshot(snapshot);
    return {};
  }

  const int tile_size = grouter_->grid() != nullptr
                            ? std::max(grouter_->grid()->getTileSize(), 1)
                            : 1;
  const RudyGrid baseline_rudy = compute_current_rudy();

  // Candidate A: compress baseline routes instead of geometric detours.
  ScenarioResult compact = baseline;
  compact.name = "baseline_compact";
  applyViaExcursionCollapse(compact.routes, std::max(4 * tile_size, 1));
  applyGuideCompression(compact.routes, std::max(8 * tile_size, 1));
  compact.metrics = compute_metrics(compact.routes);

  // Candidate B: reroute with strong but wirelength-oriented capacity sculpting.
  ScenarioResult sculpted = compact;
  sculpted.name = "field_sculpted";
  bool sculpted_available = false;
  long nets_taken_from_sculpted = 0;
  try {
    for (Net* net : nets) {
      if (net != nullptr && net->getDbNet() != nullptr) {
        grouter_->fastroute()->clearNetRoute(net->getDbNet());
      }
    }

    const auto hotspots = extractTopRudyHotspots(baseline_rudy, 96);
    const PlanarEdgeUsage baseline_usage
        = computeNormalizedBackboneUsage(grouter_, compact.routes);
    const bool prefer_horizontal
        = grouter_->grid() != nullptr
          && grouter_->grid()->getXGrids() >= grouter_->grid()->getYGrids();
    applyAggressiveCapacityField(grouter_,
                                 baseline_rudy,
                                 hotspots,
                                 min_routing_layer,
                                 max_routing_layer,
                                 0.62f,
                                 0.18f,
                                 0.40f,
                                 1.18f,
                                 0.22f,
                                 prefer_horizontal);
    applyLayerPolarityField(grouter_,
                            baseline_rudy,
                            min_routing_layer,
                            max_routing_layer,
                            1.22f,
                            0.52f,
                            0.48f);
    applyBackboneCapacityReinforcement(grouter_,
                                       baseline_rudy,
                                       baseline_usage,
                                       min_routing_layer,
                                       max_routing_layer,
                                       0.70f,
                                       1.25f,
                                       1.05f);
    applyHotspotPenalties(grouter_,
                          hotspots,
                          min_routing_layer,
                          max_routing_layer,
                          2,
                          0.62f,
                          0.35f);
    applySoftCapacityScaling(grouter_,
                             baseline_rudy,
                             min_routing_layer,
                             max_routing_layer,
                             0.48f,
                             0.95f,
                             5.0f,
                             0.50f);

    sculpted = run_existing_state("field_sculpted", nets);
    applyViaExcursionCollapse(sculpted.routes, std::max(4 * tile_size, 1));
    applyGuideCompression(sculpted.routes, std::max(8 * tile_size, 1));
    sculpted.metrics = compute_metrics(sculpted.routes);
    sculpted_available = !sculpted.routes.empty();
  } catch (...) {
    logger_->warn(GNR,
                  6019,
                  "NEWGR field_sculpted candidate failed; reverting to "
                  "baseline_compact.");
    sculpted = compact;
    sculpted.name = "field_sculpted_failed";
  }

  ScenarioResult selected = compact;
  selected.name = "netblend_compact";
  auto has_planar_guide = [](const GRoute& route) {
    for (const GSegment& segment : route) {
      if (!segment.isVia()
          && (segment.init_x != segment.final_x
              || segment.init_y != segment.final_y)) {
        return true;
      }
    }
    return false;
  };
  auto route_score = [&](const GRoute& route, int scenario_overflow) {
    long route_wl = 0;
    long route_vias = 0;
    for (const GSegment& segment : route) {
      if (segment.isVia()) {
        route_vias++;
      } else {
        route_wl += std::abs(segment.final_x - segment.init_x)
                    + std::abs(segment.final_y - segment.init_y);
      }
    }
    const double via_weight = static_cast<double>(tile_size) * 2.2;
    const double overflow_penalty
        = static_cast<double>(std::max(scenario_overflow, 0))
          * static_cast<double>(tile_size) * 8.0;
    return static_cast<double>(route_wl)
           + via_weight * static_cast<double>(route_vias) + overflow_penalty;
  };

  if (sculpted_available) {
    selected.name = "netblend_field_sculpted";
    for (const auto& [db_net, compact_route] : compact.routes) {
      auto sculpted_it = sculpted.routes.find(db_net);
      if (sculpted_it == sculpted.routes.end()) {
        continue;
      }
      const GRoute& sculpted_route = sculpted_it->second;
      const bool compact_valid = has_planar_guide(compact_route);
      const bool sculpted_valid = has_planar_guide(sculpted_route);
      if (!sculpted_valid && compact_valid) {
        continue;
      }
      if (sculpted_valid && !compact_valid) {
        selected.routes[db_net] = sculpted_route;
        nets_taken_from_sculpted++;
        continue;
      }
      const double compact_score = route_score(compact_route, compact.overflow);
      const double sculpted_score = route_score(sculpted_route, sculpted.overflow);
      if (sculpted_score + 1e-3 < compact_score) {
        selected.routes[db_net] = sculpted_route;
        nets_taken_from_sculpted++;
      }
    }
  }

  for (const auto& [db_net, baseline_route] : compact.routes) {
    auto selected_it = selected.routes.find(db_net);
    if (selected_it == selected.routes.end()
        || !has_planar_guide(selected_it->second)) {
      selected.routes[db_net] = baseline_route;
    }
  }
  selected.metrics = compute_metrics(selected.routes);

  const double compact_delta_wl
      = compact.metrics.wirelength_um - baseline.metrics.wirelength_um;
  const long compact_delta_vias
      = compact.metrics.via_count - baseline.metrics.via_count;
  const double sculpted_delta_wl
      = sculpted.metrics.wirelength_um - baseline.metrics.wirelength_um;
  const long sculpted_delta_vias
      = sculpted.metrics.via_count - baseline.metrics.via_count;
  const double selected_delta_wl
      = selected.metrics.wirelength_um - baseline.metrics.wirelength_um;
  const long selected_delta_vias
      = selected.metrics.via_count - baseline.metrics.via_count;

  logger_->warn(GNR,
                6016,
                "NEWGR candidate {}: wl {:.0f} um vias {} (delta wl {:+.0f} "
                "um, delta vias {:+d}, overflow {}).",
                compact.name,
                compact.metrics.wirelength_um,
                compact.metrics.via_count,
                compact_delta_wl,
                compact_delta_vias,
                compact.overflow);
  logger_->warn(GNR,
                6017,
                "NEWGR candidate {}: wl {:.0f} um vias {} (delta wl {:+.0f} "
                "um, delta vias {:+d}, overflow {}).",
                sculpted.name,
                sculpted.metrics.wirelength_um,
                sculpted.metrics.via_count,
                sculpted_delta_wl,
                sculpted_delta_vias,
                sculpted.overflow);
  logger_->warn(GNR,
                6018,
                "NEWGR selected {} over baseline {:.0f} um vias {} -> {:.0f} "
                "um vias {} (delta wl {:+.0f} um, delta vias {:+d}).",
                selected.name,
                baseline.metrics.wirelength_um,
                baseline.metrics.via_count,
                selected.metrics.wirelength_um,
                selected.metrics.via_count,
                selected_delta_wl,
                selected_delta_vias);
  logger_->warn(GNR,
                6020,
                "NEWGR blended {} nets from {} into final selection.",
                nets_taken_from_sculpted,
                sculpted.name);

  restore_snapshot(snapshot);
  return selected.routes;
}

}  // namespace grt
