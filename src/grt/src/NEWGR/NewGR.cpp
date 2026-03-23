#include "NEWGR/NewGR.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
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

void applyAggressiveDoglegShortcuts(NetRouteMap& routes,
                                    int max_middle_len,
                                    int min_wirelength_gain)
{
  max_middle_len = std::max(max_middle_len, 1);
  min_wirelength_gain = std::max(min_wirelength_gain, 1);

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.size() < 3) {
      continue;
    }

    std::vector<GSegment> current = route;
    bool changed = false;

    for (int round = 0; round < 4; ++round) {
      if (current.size() < 3) {
        break;
      }

      bool round_changed = false;
      std::vector<GSegment> rewritten;
      rewritten.reserve(current.size());

      for (size_t idx = 0; idx < current.size();) {
        if (idx + 2 < current.size()) {
          const GSegment& first = current[idx];
          const GSegment& middle = current[idx + 1];
          const GSegment& last = current[idx + 2];

          const bool contiguous = samePointAndLayer(first, true, middle, true)
                                  && samePointAndLayer(middle, false, last, true);
          const bool same_layer = !first.isVia() && !middle.isVia() && !last.isVia()
                                  && first.init_layer == first.final_layer
                                  && middle.init_layer == middle.final_layer
                                  && last.init_layer == last.final_layer
                                  && first.init_layer == middle.init_layer
                                  && middle.init_layer == last.init_layer;

          if (contiguous && same_layer) {
            const bool hvh = isPlanarHorizontal(first) && isPlanarVertical(middle)
                             && isPlanarHorizontal(last)
                             && first.init_y == last.final_y;
            const bool vhv = isPlanarVertical(first) && isPlanarHorizontal(middle)
                             && isPlanarVertical(last)
                             && first.init_x == last.final_x;

            if (hvh || vhv) {
              const long middle_len = std::abs(middle.final_x - middle.init_x)
                                      + std::abs(middle.final_y - middle.init_y);
              const long old_len = first.length() + middle.length() + last.length();
              const GSegment direct(first.init_x,
                                    first.init_y,
                                    first.init_layer,
                                    last.final_x,
                                    last.final_y,
                                    last.final_layer);
              const long new_len = std::abs(direct.final_x - direct.init_x)
                                   + std::abs(direct.final_y - direct.init_y);
              const long gain = old_len - new_len;

              if (middle_len <= max_middle_len && gain >= min_wirelength_gain) {
                appendCompressedSegment(rewritten, direct);
                idx += 3;
                round_changed = true;
                continue;
              }
            }
          }
        }

        appendCompressedSegment(rewritten, current[idx]);
        ++idx;
      }

      if (!rewritten.empty()) {
        current.swap(rewritten);
      }
      if (!round_changed) {
        break;
      }
      changed = true;
    }

    if (changed && !current.empty()) {
      route.swap(current);
    }
  }
}

struct RouteNode
{
  int x = 0;
  int y = 0;
  int layer = 0;
};

std::vector<RouteNode> collectUniqueRouteNodes(const GRoute& route)
{
  std::vector<RouteNode> nodes;
  nodes.reserve(route.size() * 2);
  for (const GSegment& segment : route) {
    nodes.push_back(RouteNode{segment.init_x, segment.init_y, segment.init_layer});
    nodes.push_back(
        RouteNode{segment.final_x, segment.final_y, segment.final_layer});
  }
  std::sort(nodes.begin(),
            nodes.end(),
            [](const RouteNode& lhs, const RouteNode& rhs) {
              return std::tie(lhs.x, lhs.y, lhs.layer)
                     < std::tie(rhs.x, rhs.y, rhs.layer);
            });
  nodes.erase(
      std::unique(nodes.begin(),
                  nodes.end(),
                  [](const RouteNode& lhs, const RouteNode& rhs) {
                    return lhs.x == rhs.x && lhs.y == rhs.y
                           && lhs.layer == rhs.layer;
                  }),
      nodes.end());
  return nodes;
}

int chooseDominantLayer(const std::vector<RouteNode>& nodes,
                        int min_layer,
                        int max_layer)
{
  std::vector<int> layer_votes(std::max(max_layer + 1, 0), 0);
  for (const RouteNode& node : nodes) {
    if (node.layer < min_layer || node.layer > max_layer
        || node.layer >= static_cast<int>(layer_votes.size())) {
      continue;
    }
    layer_votes[node.layer]++;
  }

  int best_layer = std::clamp(min_layer, min_layer, max_layer);
  int best_votes = -1;
  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const int votes = layer < static_cast<int>(layer_votes.size())
                          ? layer_votes[layer]
                          : 0;
    if (votes > best_votes) {
      best_votes = votes;
      best_layer = layer;
    }
  }
  return best_layer;
}

struct GlobalPortal
{
  int x = 0;
  int y = 0;
  bool valid = false;
};

GlobalPortal computeLowRudyPortal(const RudyGrid& normalized_rudy, Grid* grid)
{
  GlobalPortal portal;
  if (grid == nullptr) {
    return portal;
  }

  const int tile = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  portal.x = x_min + (grid->getXGrids() / 2) * tile;
  portal.y = y_min + (grid->getYGrids() / 2) * tile;
  portal.valid = true;

  if (normalized_rudy.empty() || normalized_rudy.front().empty()) {
    return portal;
  }

  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  std::vector<float> col_cost(x_tiles, 0.0f);
  std::vector<float> row_cost(y_tiles, 0.0f);
  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      const float val = normalized_rudy[x][y];
      col_cost[x] += val;
      row_cost[y] += val;
    }
  }

  for (float& cost : col_cost) {
    cost /= static_cast<float>(std::max(y_tiles, 1));
  }
  for (float& cost : row_cost) {
    cost /= static_cast<float>(std::max(x_tiles, 1));
  }

  int best_x = 0;
  int best_y = 0;
  float best_x_score = std::numeric_limits<float>::max();
  float best_y_score = std::numeric_limits<float>::max();
  const float center_x = static_cast<float>(x_tiles - 1) * 0.5f;
  const float center_y = static_cast<float>(y_tiles - 1) * 0.5f;
  for (int x = 0; x < x_tiles; ++x) {
    const float center_penalty = 0.10f * std::abs(static_cast<float>(x) - center_x)
                                 / std::max(center_x, 1.0f);
    const float score = col_cost[x] + center_penalty;
    if (score < best_x_score) {
      best_x_score = score;
      best_x = x;
    }
  }
  for (int y = 0; y < y_tiles; ++y) {
    const float center_penalty = 0.10f * std::abs(static_cast<float>(y) - center_y)
                                 / std::max(center_y, 1.0f);
    const float score = row_cost[y] + center_penalty;
    if (score < best_y_score) {
      best_y_score = score;
      best_y = y;
    }
  }

  portal.x = x_min + best_x * tile;
  portal.y = y_min + best_y * tile;
  return portal;
}

float estimatePathRudy(const RudyGrid& normalized_rudy,
                       int gx0,
                       int gy0,
                       int gx1,
                       int gy1,
                       bool horizontal_first);

GlobalPortal computeQuadrantLowRudyPortal(const RudyGrid& normalized_rudy,
                                          Grid* grid,
                                          bool east_half,
                                          bool north_half)
{
  GlobalPortal portal;
  if (grid == nullptr) {
    return portal;
  }

  const int tile = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int x_grids = std::max(grid->getXGrids(), 1);
  const int y_grids = std::max(grid->getYGrids(), 1);
  const int mid_x = x_grids / 2;
  const int mid_y = y_grids / 2;
  const int fallback_gx = east_half ? std::max(mid_x + x_grids / 4, 0)
                                    : std::max(mid_x - x_grids / 4, 0);
  const int fallback_gy = north_half ? std::max(mid_y + y_grids / 4, 0)
                                     : std::max(mid_y - y_grids / 4, 0);

  if (normalized_rudy.empty() || normalized_rudy.front().empty()) {
    portal.x = x_min + std::clamp(fallback_gx, 0, x_grids - 1) * tile;
    portal.y = y_min + std::clamp(fallback_gy, 0, y_grids - 1) * tile;
    portal.valid = true;
    return portal;
  }

  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  const int x_mid = x_tiles / 2;
  const int y_mid = y_tiles / 2;
  const int x_start = east_half ? x_mid : 0;
  const int x_end = east_half ? x_tiles : std::max(x_mid, 1);
  const int y_start = north_half ? y_mid : 0;
  const int y_end = north_half ? y_tiles : std::max(y_mid, 1);

  if (x_start >= x_end || y_start >= y_end) {
    portal.x = x_min + std::clamp(fallback_gx, 0, x_grids - 1) * tile;
    portal.y = y_min + std::clamp(fallback_gy, 0, y_grids - 1) * tile;
    portal.valid = true;
    return portal;
  }

  const float quadrant_cx
      = 0.5f * static_cast<float>(x_start + std::max(x_end - 1, x_start));
  const float quadrant_cy
      = 0.5f * static_cast<float>(y_start + std::max(y_end - 1, y_start));
  int best_x = x_start;
  int best_y = y_start;
  float best_score = std::numeric_limits<float>::max();
  const float norm_x = std::max(static_cast<float>(x_end - x_start), 1.0f);
  const float norm_y = std::max(static_cast<float>(y_end - y_start), 1.0f);

  for (int x = x_start; x < x_end; ++x) {
    for (int y = y_start; y < y_end; ++y) {
      const float rudy = normalized_rudy[x][y];
      const float dx = std::abs(static_cast<float>(x) - quadrant_cx) / norm_x;
      const float dy = std::abs(static_cast<float>(y) - quadrant_cy) / norm_y;
      const float score = rudy + 0.09f * (dx + dy);
      if (score < best_score) {
        best_score = score;
        best_x = x;
        best_y = y;
      }
    }
  }

  portal.x = x_min + std::clamp(best_x, 0, x_grids - 1) * tile;
  portal.y = y_min + std::clamp(best_y, 0, y_grids - 1) * tile;
  portal.valid = true;
  return portal;
}

void applyQuadrantPortalHypergraphRebuild(GlobalRouter* grouter,
                                          NetRouteMap& routes,
                                          const RudyGrid& normalized_rudy,
                                          int min_unique_nodes,
                                          int max_unique_nodes,
                                          int coverage_percent,
                                          int min_layer,
                                          int max_layer)
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
  const int x_tiles
      = normalized_rudy.empty() ? 0 : static_cast<int>(normalized_rudy.size());
  const int y_tiles = normalized_rudy.empty() ? 0
                                              : static_cast<int>(
                                                    normalized_rudy.front().size());

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

  min_unique_nodes = std::max(min_unique_nodes, 3);
  max_unique_nodes = std::max(max_unique_nodes, min_unique_nodes);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);

  std::array<GlobalPortal, 4> portals{
      computeQuadrantLowRudyPortal(normalized_rudy, grid, false, false),
      computeQuadrantLowRudyPortal(normalized_rudy, grid, true, false),
      computeQuadrantLowRudyPortal(normalized_rudy, grid, false, true),
      computeQuadrantLowRudyPortal(normalized_rudy, grid, true, true)};

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const auto net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(net_key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes
        || static_cast<int>(nodes.size()) > max_unique_nodes) {
      continue;
    }

    long cx_acc = 0;
    long cy_acc = 0;
    int min_x = nodes.front().x;
    int max_x = nodes.front().x;
    int min_y = nodes.front().y;
    int max_y = nodes.front().y;
    for (const RouteNode& node : nodes) {
      cx_acc += node.x;
      cy_acc += node.y;
      min_x = std::min(min_x, node.x);
      max_x = std::max(max_x, node.x);
      min_y = std::min(min_y, node.y);
      max_y = std::max(max_y, node.y);
    }
    const int node_count = std::max(static_cast<int>(nodes.size()), 1);
    const int centroid_x = static_cast<int>(cx_acc / node_count);
    const int centroid_y = static_cast<int>(cy_acc / node_count);
    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(nodes.size() * 5 + 16);
    std::array<bool, 4> portal_used{false, false, false, false};

    for (size_t idx = 0; idx < nodes.size(); ++idx) {
      int cur_x = std::clamp(nodes[idx].x, x_min, x_max);
      int cur_y = std::clamp(nodes[idx].y, y_min, y_max);
      int cur_layer = nodes[idx].layer;

      if (cur_layer != trunk_layer) {
        const int step = (trunk_layer > cur_layer) ? 1 : -1;
        while (cur_layer != trunk_layer) {
          const int next_layer = cur_layer + step;
          appendSegment(
              rebuilt, cur_x, cur_y, cur_layer, cur_x, cur_y, next_layer);
          cur_layer = next_layer;
        }
      }

      const bool east = cur_x >= centroid_x;
      const bool north = cur_y >= centroid_y;
      const int preferred_portal = (east ? 1 : 0) + (north ? 2 : 0);
      int selected_portal = preferred_portal;
      double best_cost = std::numeric_limits<double>::max();
      bool best_horizontal_first = true;
      for (int p = 0; p < static_cast<int>(portals.size()); ++p) {
        const GlobalPortal& portal = portals[p];
        if (!portal.valid) {
          continue;
        }
        const long dist = std::abs(cur_x - portal.x) + std::abs(cur_y - portal.y);
        const float h_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(portal.x),
                                              to_grid_y(portal.y),
                                              true);
        const float v_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(portal.x),
                                              to_grid_y(portal.y),
                                              false);
        const float edge_rudy = std::min(h_rudy, v_rudy);
        double cost = static_cast<double>(dist)
                      + static_cast<double>(tile) * 3.2
                            * static_cast<double>(edge_rudy);
        if (p == preferred_portal) {
          cost -= static_cast<double>(tile) * 2.2;
        }
        if (cost + 1e-9 < best_cost) {
          best_cost = cost;
          selected_portal = p;
          best_horizontal_first = h_rudy <= v_rudy;
        }
      }

      const GlobalPortal& portal = portals[selected_portal];
      if (!portal.valid) {
        continue;
      }

      portal_used[selected_portal] = true;
      const bool horizontal_first
          = ((net_key + static_cast<std::uint64_t>(idx) * 41ULL) & 1ULL) == 0ULL
                ? best_horizontal_first
                : !best_horizontal_first;
      if (horizontal_first) {
        appendSegment(
            rebuilt, cur_x, cur_y, trunk_layer, portal.x, cur_y, trunk_layer);
        appendSegment(rebuilt,
                      portal.x,
                      cur_y,
                      trunk_layer,
                      portal.x,
                      portal.y,
                      trunk_layer);
      } else {
        appendSegment(
            rebuilt, cur_x, cur_y, trunk_layer, cur_x, portal.y, trunk_layer);
        appendSegment(rebuilt,
                      cur_x,
                      portal.y,
                      trunk_layer,
                      portal.x,
                      portal.y,
                      trunk_layer);
      }
    }

    std::vector<int> used_x;
    std::vector<int> used_y;
    used_x.reserve(portals.size());
    used_y.reserve(portals.size());
    for (int p = 0; p < static_cast<int>(portals.size()); ++p) {
      if (portal_used[p] && portals[p].valid) {
        used_x.push_back(portals[p].x);
        used_y.push_back(portals[p].y);
      }
    }
    if (used_x.empty()) {
      continue;
    }

    std::nth_element(used_x.begin(), used_x.begin() + used_x.size() / 2, used_x.end());
    std::nth_element(used_y.begin(), used_y.begin() + used_y.size() / 2, used_y.end());
    const int switch_x = std::clamp(used_x[used_x.size() / 2], x_min, x_max);
    const int switch_y = std::clamp(used_y[used_y.size() / 2], y_min, y_max);

    for (int p = 0; p < static_cast<int>(portals.size()); ++p) {
      if (!portal_used[p] || !portals[p].valid) {
        continue;
      }
      const GlobalPortal& portal = portals[p];
      const float h_rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(portal.x),
                                            to_grid_y(portal.y),
                                            to_grid_x(switch_x),
                                            to_grid_y(switch_y),
                                            true);
      const float v_rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(portal.x),
                                            to_grid_y(portal.y),
                                            to_grid_x(switch_x),
                                            to_grid_y(switch_y),
                                            false);
      if (h_rudy <= v_rudy) {
        appendSegment(rebuilt,
                      portal.x,
                      portal.y,
                      trunk_layer,
                      switch_x,
                      portal.y,
                      trunk_layer);
        appendSegment(rebuilt,
                      switch_x,
                      portal.y,
                      trunk_layer,
                      switch_x,
                      switch_y,
                      trunk_layer);
      } else {
        appendSegment(rebuilt,
                      portal.x,
                      portal.y,
                      trunk_layer,
                      portal.x,
                      switch_y,
                      trunk_layer);
        appendSegment(rebuilt,
                      portal.x,
                      switch_y,
                      trunk_layer,
                      switch_x,
                      switch_y,
                      trunk_layer);
      }
    }

    appendSegment(rebuilt, min_x, switch_y, trunk_layer, max_x, switch_y, trunk_layer);
    appendSegment(rebuilt, switch_x, min_y, trunk_layer, switch_x, max_y, trunk_layer);

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }
    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

void applyGlobalPortalRebuild(GlobalRouter* grouter,
                              NetRouteMap& routes,
                              const RudyGrid& normalized_rudy,
                              int min_unique_nodes,
                              int coverage_percent,
                              int min_layer,
                              int max_layer)
{
  if (grouter == nullptr || grouter->grid() == nullptr || routes.empty()) {
    return;
  }

  Grid* grid = grouter->grid();
  min_unique_nodes = std::max(min_unique_nodes, 3);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int x_max = grid->getXMax();
  const int y_max = grid->getYMax();
  const GlobalPortal global_portal = computeLowRudyPortal(normalized_rudy, grid);

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty()) {
      continue;
    }

    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes) {
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(nodes.size());
    ys.reserve(nodes.size());
    int min_x = nodes.front().x;
    int max_x = nodes.front().x;
    int min_y = nodes.front().y;
    int max_y = nodes.front().y;
    for (const RouteNode& node : nodes) {
      xs.push_back(node.x);
      ys.push_back(node.y);
      min_x = std::min(min_x, node.x);
      max_x = std::max(max_x, node.x);
      min_y = std::min(min_y, node.y);
      max_y = std::max(max_y, node.y);
    }
    std::nth_element(xs.begin(), xs.begin() + xs.size() / 2, xs.end());
    std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
    const int local_x = xs[xs.size() / 2];
    const int local_y = ys[ys.size() / 2];
    const int portal_x = std::clamp(
        (global_portal.x + 2 * local_x) / 3, x_min, x_max);
    const int portal_y = std::clamp(
        (2 * global_portal.y + local_y) / 3, y_min, y_max);
    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(nodes.size() * 4 + 4);
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
      const RouteNode& node = nodes[idx];
      int cur_x = node.x;
      int cur_y = node.y;
      int cur_layer = node.layer;

      if (cur_layer != trunk_layer) {
        const int step = (trunk_layer > cur_layer) ? 1 : -1;
        while (cur_layer != trunk_layer) {
          const int next_layer = cur_layer + step;
          appendSegment(
              rebuilt, cur_x, cur_y, cur_layer, cur_x, cur_y, next_layer);
          cur_layer = next_layer;
        }
      }

      const bool horizontal_first
          = ((key + static_cast<std::uint64_t>(idx) * 17ULL) & 1ULL) == 0ULL;
      if (horizontal_first) {
        appendSegment(rebuilt,
                      cur_x,
                      cur_y,
                      trunk_layer,
                      portal_x,
                      cur_y,
                      trunk_layer);
        appendSegment(rebuilt,
                      portal_x,
                      cur_y,
                      trunk_layer,
                      portal_x,
                      portal_y,
                      trunk_layer);
      } else {
        appendSegment(rebuilt,
                      cur_x,
                      cur_y,
                      trunk_layer,
                      cur_x,
                      portal_y,
                      trunk_layer);
        appendSegment(rebuilt,
                      cur_x,
                      portal_y,
                      trunk_layer,
                      portal_x,
                      portal_y,
                      trunk_layer);
      }
    }

    appendSegment(
        rebuilt, min_x, portal_y, trunk_layer, max_x, portal_y, trunk_layer);
    appendSegment(
        rebuilt, portal_x, min_y, trunk_layer, portal_x, max_y, trunk_layer);

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }
    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

float estimatePathRudy(const RudyGrid& normalized_rudy,
                       int gx0,
                       int gy0,
                       int gx1,
                       int gy1,
                       bool horizontal_first)
{
  if (normalized_rudy.empty() || normalized_rudy.front().empty()) {
    return 0.0f;
  }

  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  gx0 = std::clamp(gx0, 0, x_tiles - 1);
  gx1 = std::clamp(gx1, 0, x_tiles - 1);
  gy0 = std::clamp(gy0, 0, y_tiles - 1);
  gy1 = std::clamp(gy1, 0, y_tiles - 1);

  auto sample = [&](int gx, int gy) {
    gx = std::clamp(gx, 0, x_tiles - 1);
    gy = std::clamp(gy, 0, y_tiles - 1);
    return normalized_rudy[gx][gy];
  };

  int cur_x = gx0;
  int cur_y = gy0;
  float rudy_acc = sample(cur_x, cur_y);
  int samples = 1;

  auto walk_x = [&]() {
    while (cur_x != gx1) {
      cur_x += (gx1 > cur_x) ? 1 : -1;
      rudy_acc += sample(cur_x, cur_y);
      samples++;
    }
  };
  auto walk_y = [&]() {
    while (cur_y != gy1) {
      cur_y += (gy1 > cur_y) ? 1 : -1;
      rudy_acc += sample(cur_x, cur_y);
      samples++;
    }
  };

  if (horizontal_first) {
    walk_x();
    walk_y();
  } else {
    walk_y();
    walk_x();
  }

  return rudy_acc / static_cast<float>(std::max(samples, 1));
}

void applyRmstTrunkRebuild(GlobalRouter* grouter,
                           NetRouteMap& routes,
                           const RudyGrid& normalized_rudy,
                           int min_unique_nodes,
                           int max_unique_nodes,
                           int coverage_percent,
                           int min_layer,
                           int max_layer)
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
  const int x_tiles
      = normalized_rudy.empty() ? 0 : static_cast<int>(normalized_rudy.size());
  const int y_tiles = normalized_rudy.empty() ? 0
                                              : static_cast<int>(
                                                    normalized_rudy.front().size());

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

  min_unique_nodes = std::max(min_unique_nodes, 3);
  max_unique_nodes = std::max(max_unique_nodes, min_unique_nodes);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const auto net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(net_key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes
        || static_cast<int>(nodes.size()) > max_unique_nodes) {
      continue;
    }

    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);

    std::vector<RouteNode> planar_nodes;
    planar_nodes.reserve(nodes.size());
    for (const RouteNode& node : nodes) {
      const int clamped_x = std::clamp(node.x, x_min, x_max);
      const int clamped_y = std::clamp(node.y, y_min, y_max);
      planar_nodes.push_back(RouteNode{clamped_x, clamped_y, trunk_layer});
    }
    std::sort(planar_nodes.begin(),
              planar_nodes.end(),
              [](const RouteNode& lhs, const RouteNode& rhs) {
                return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y);
              });
    planar_nodes.erase(
        std::unique(planar_nodes.begin(),
                    planar_nodes.end(),
                    [](const RouteNode& lhs, const RouteNode& rhs) {
                      return lhs.x == rhs.x && lhs.y == rhs.y;
                    }),
        planar_nodes.end());
    if (static_cast<int>(planar_nodes.size()) < min_unique_nodes) {
      continue;
    }

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(planar_nodes.size() * 3);
    const int n = planar_nodes.size();
    std::vector<int> in_tree(n, 0);

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(n);
    ys.reserve(n);
    for (const RouteNode& node : planar_nodes) {
      xs.push_back(node.x);
      ys.push_back(node.y);
    }
    std::nth_element(xs.begin(), xs.begin() + xs.size() / 2, xs.end());
    std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
    const int seed_x = xs[xs.size() / 2];
    const int seed_y = ys[ys.size() / 2];

    int root = 0;
    long best_seed_dist = std::numeric_limits<long>::max();
    for (int i = 0; i < n; ++i) {
      const long dist = std::abs(planar_nodes[i].x - seed_x)
                        + std::abs(planar_nodes[i].y - seed_y);
      if (dist < best_seed_dist) {
        best_seed_dist = dist;
        root = i;
      }
    }
    in_tree[root] = 1;

    for (int connected = 1; connected < n; ++connected) {
      int best_u = -1;
      int best_v = -1;
      double best_cost = std::numeric_limits<double>::max();
      long best_dist = std::numeric_limits<long>::max();

      for (int u = 0; u < n; ++u) {
        if (!in_tree[u]) {
          continue;
        }
        for (int v = 0; v < n; ++v) {
          if (in_tree[v]) {
            continue;
          }
          const long dist = std::abs(planar_nodes[u].x - planar_nodes[v].x)
                            + std::abs(planar_nodes[u].y - planar_nodes[v].y);
          const int ux = to_grid_x(planar_nodes[u].x);
          const int uy = to_grid_y(planar_nodes[u].y);
          const int vx = to_grid_x(planar_nodes[v].x);
          const int vy = to_grid_y(planar_nodes[v].y);
          const float h_rudy
              = estimatePathRudy(normalized_rudy, ux, uy, vx, vy, true);
          const float v_rudy
              = estimatePathRudy(normalized_rudy, ux, uy, vx, vy, false);
          const float edge_rudy = std::min(h_rudy, v_rudy);
          const double cost
              = static_cast<double>(dist)
                + static_cast<double>(tile) * 2.0 * static_cast<double>(edge_rudy);
          if (cost + 1e-9 < best_cost
              || (std::abs(cost - best_cost) <= 1e-9 && dist < best_dist)) {
            best_cost = cost;
            best_dist = dist;
            best_u = u;
            best_v = v;
          }
        }
      }

      if (best_u < 0 || best_v < 0) {
        break;
      }

      const RouteNode& from = planar_nodes[best_u];
      const RouteNode& to = planar_nodes[best_v];
      const int from_gx = to_grid_x(from.x);
      const int from_gy = to_grid_y(from.y);
      const int to_gx = to_grid_x(to.x);
      const int to_gy = to_grid_y(to.y);
      const float h_rudy
          = estimatePathRudy(normalized_rudy, from_gx, from_gy, to_gx, to_gy, true);
      const float v_rudy
          = estimatePathRudy(normalized_rudy, from_gx, from_gy, to_gx, to_gy, false);
      bool horizontal_first = h_rudy + 1e-4f < v_rudy;
      if (std::fabs(h_rudy - v_rudy) <= 0.015f) {
        horizontal_first
            = ((net_key + static_cast<std::uint64_t>(best_u) * 31ULL
                + static_cast<std::uint64_t>(best_v) * 17ULL)
               & 1ULL)
              == 0ULL;
      }

      if (horizontal_first) {
        appendSegment(rebuilt,
                      from.x,
                      from.y,
                      trunk_layer,
                      to.x,
                      from.y,
                      trunk_layer);
        appendSegment(
            rebuilt, to.x, from.y, trunk_layer, to.x, to.y, trunk_layer);
      } else {
        appendSegment(rebuilt,
                      from.x,
                      from.y,
                      trunk_layer,
                      from.x,
                      to.y,
                      trunk_layer);
        appendSegment(
            rebuilt, from.x, to.y, trunk_layer, to.x, to.y, trunk_layer);
      }
      in_tree[best_v] = 1;
    }

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }
    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

void applyBipolarPortalBackboneRebuild(GlobalRouter* grouter,
                                       NetRouteMap& routes,
                                       const RudyGrid& normalized_rudy,
                                       int min_unique_nodes,
                                       int max_unique_nodes,
                                       int coverage_percent,
                                       int min_layer,
                                       int max_layer)
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
  const int x_tiles
      = normalized_rudy.empty() ? 0 : static_cast<int>(normalized_rudy.size());
  const int y_tiles = normalized_rudy.empty() ? 0
                                              : static_cast<int>(
                                                    normalized_rudy.front().size());

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

  min_unique_nodes = std::max(min_unique_nodes, 3);
  max_unique_nodes = std::max(max_unique_nodes, min_unique_nodes);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const auto net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(net_key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes
        || static_cast<int>(nodes.size()) > max_unique_nodes) {
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(nodes.size());
    ys.reserve(nodes.size());
    int min_x = nodes.front().x;
    int max_x = nodes.front().x;
    int min_y = nodes.front().y;
    int max_y = nodes.front().y;
    for (const RouteNode& node : nodes) {
      const int clamped_x = std::clamp(node.x, x_min, x_max);
      const int clamped_y = std::clamp(node.y, y_min, y_max);
      xs.push_back(clamped_x);
      ys.push_back(clamped_y);
      min_x = std::min(min_x, clamped_x);
      max_x = std::max(max_x, clamped_x);
      min_y = std::min(min_y, clamped_y);
      max_y = std::max(max_y, clamped_y);
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());
    if (xs.empty() || ys.empty()) {
      continue;
    }

    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);
    const int width = max_x - min_x;
    const int height = max_y - min_y;
    const bool horizontal_backbone = width >= height;

    int anchor_a_x = xs[(xs.size() - 1) / 3];
    int anchor_b_x = xs[(2 * (xs.size() - 1)) / 3];
    int anchor_a_y = ys[(ys.size() - 1) / 3];
    int anchor_b_y = ys[(2 * (ys.size() - 1)) / 3];
    int backbone_x = xs[xs.size() / 2];
    int backbone_y = ys[ys.size() / 2];

    if (horizontal_backbone) {
      float best_score = std::numeric_limits<float>::max();
      int best_y = backbone_y;
      for (int offset = -4; offset <= 4; ++offset) {
        const int cand_y = std::clamp(backbone_y + offset * tile, y_min, y_max);
        const float rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(anchor_a_x),
                                            to_grid_y(cand_y),
                                            to_grid_x(anchor_b_x),
                                            to_grid_y(cand_y),
                                            true);
        const float score = rudy + 0.03f * static_cast<float>(std::abs(offset));
        if (score < best_score) {
          best_score = score;
          best_y = cand_y;
        }
      }
      backbone_y = best_y;
      if (anchor_a_x == anchor_b_x) {
        anchor_a_x = std::clamp(anchor_a_x - tile, x_min, x_max);
        anchor_b_x = std::clamp(anchor_b_x + tile, x_min, x_max);
      }
    } else {
      float best_score = std::numeric_limits<float>::max();
      int best_x = backbone_x;
      for (int offset = -4; offset <= 4; ++offset) {
        const int cand_x = std::clamp(backbone_x + offset * tile, x_min, x_max);
        const float rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(cand_x),
                                            to_grid_y(anchor_a_y),
                                            to_grid_x(cand_x),
                                            to_grid_y(anchor_b_y),
                                            false);
        const float score = rudy + 0.03f * static_cast<float>(std::abs(offset));
        if (score < best_score) {
          best_score = score;
          best_x = cand_x;
        }
      }
      backbone_x = best_x;
      if (anchor_a_y == anchor_b_y) {
        anchor_a_y = std::clamp(anchor_a_y - tile, y_min, y_max);
        anchor_b_y = std::clamp(anchor_b_y + tile, y_min, y_max);
      }
    }

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(nodes.size() * 4 + 4);
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
      int cur_x = std::clamp(nodes[idx].x, x_min, x_max);
      int cur_y = std::clamp(nodes[idx].y, y_min, y_max);
      int cur_layer = nodes[idx].layer;

      if (cur_layer != trunk_layer) {
        const int step = (trunk_layer > cur_layer) ? 1 : -1;
        while (cur_layer != trunk_layer) {
          const int next_layer = cur_layer + step;
          appendSegment(
              rebuilt, cur_x, cur_y, cur_layer, cur_x, cur_y, next_layer);
          cur_layer = next_layer;
        }
      }

      if (horizontal_backbone) {
        const int dist_a = std::abs(cur_x - anchor_a_x);
        const int dist_b = std::abs(cur_x - anchor_b_x);
        const int target_x = dist_a <= dist_b ? anchor_a_x : anchor_b_x;
        const float h_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(target_x),
                                              to_grid_y(backbone_y),
                                              true);
        const float v_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(target_x),
                                              to_grid_y(backbone_y),
                                              false);
        const bool horizontal_first = h_rudy <= v_rudy;
        if (horizontal_first) {
          appendSegment(
              rebuilt, cur_x, cur_y, trunk_layer, target_x, cur_y, trunk_layer);
          appendSegment(rebuilt,
                        target_x,
                        cur_y,
                        trunk_layer,
                        target_x,
                        backbone_y,
                        trunk_layer);
        } else {
          appendSegment(
              rebuilt, cur_x, cur_y, trunk_layer, cur_x, backbone_y, trunk_layer);
          appendSegment(rebuilt,
                        cur_x,
                        backbone_y,
                        trunk_layer,
                        target_x,
                        backbone_y,
                        trunk_layer);
        }
      } else {
        const int dist_a = std::abs(cur_y - anchor_a_y);
        const int dist_b = std::abs(cur_y - anchor_b_y);
        const int target_y = dist_a <= dist_b ? anchor_a_y : anchor_b_y;
        const float h_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(backbone_x),
                                              to_grid_y(target_y),
                                              true);
        const float v_rudy = estimatePathRudy(normalized_rudy,
                                              to_grid_x(cur_x),
                                              to_grid_y(cur_y),
                                              to_grid_x(backbone_x),
                                              to_grid_y(target_y),
                                              false);
        const bool horizontal_first = h_rudy <= v_rudy;
        if (horizontal_first) {
          appendSegment(
              rebuilt, cur_x, cur_y, trunk_layer, backbone_x, cur_y, trunk_layer);
          appendSegment(rebuilt,
                        backbone_x,
                        cur_y,
                        trunk_layer,
                        backbone_x,
                        target_y,
                        trunk_layer);
        } else {
          appendSegment(
              rebuilt, cur_x, cur_y, trunk_layer, cur_x, target_y, trunk_layer);
          appendSegment(rebuilt,
                        cur_x,
                        target_y,
                        trunk_layer,
                        backbone_x,
                        target_y,
                        trunk_layer);
        }
      }
    }

    if (horizontal_backbone) {
      appendSegment(rebuilt,
                    anchor_a_x,
                    backbone_y,
                    trunk_layer,
                    anchor_b_x,
                    backbone_y,
                    trunk_layer);
      appendSegment(
          rebuilt, min_x, backbone_y, trunk_layer, max_x, backbone_y, trunk_layer);
    } else {
      appendSegment(rebuilt,
                    backbone_x,
                    anchor_a_y,
                    trunk_layer,
                    backbone_x,
                    anchor_b_y,
                    trunk_layer);
      appendSegment(
          rebuilt, backbone_x, min_y, trunk_layer, backbone_x, max_y, trunk_layer);
    }

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }
    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

void connectNodeToMedianSpine(std::vector<GSegment>& rebuilt,
                              const RouteNode& node,
                              int spine_x,
                              int spine_y,
                              int trunk_layer)
{
  int cur_x = node.x;
  int cur_y = node.y;
  int cur_layer = node.layer;
  if (cur_layer != trunk_layer) {
    const int step = (trunk_layer > cur_layer) ? 1 : -1;
    while (cur_layer != trunk_layer) {
      const int next_layer = cur_layer + step;
      appendSegment(
          rebuilt, cur_x, cur_y, cur_layer, cur_x, cur_y, next_layer);
      cur_layer = next_layer;
    }
  }
  appendSegment(rebuilt, cur_x, cur_y, cur_layer, spine_x, cur_y, trunk_layer);
  appendSegment(rebuilt, spine_x, cur_y, trunk_layer, spine_x, spine_y, trunk_layer);
}

void applyMedianSpineRebuild(NetRouteMap& routes,
                             int min_unique_nodes,
                             int coverage_percent,
                             int min_layer,
                             int max_layer)
{
  min_unique_nodes = std::max(min_unique_nodes, 3);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty()) {
      continue;
    }

    const auto net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(net_key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes) {
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(nodes.size());
    ys.reserve(nodes.size());
    int min_x = nodes.front().x;
    int max_x = nodes.front().x;
    int min_y = nodes.front().y;
    int max_y = nodes.front().y;
    for (const RouteNode& node : nodes) {
      xs.push_back(node.x);
      ys.push_back(node.y);
      min_x = std::min(min_x, node.x);
      max_x = std::max(max_x, node.x);
      min_y = std::min(min_y, node.y);
      max_y = std::max(max_y, node.y);
    }

    std::nth_element(xs.begin(), xs.begin() + xs.size() / 2, xs.end());
    std::nth_element(ys.begin(), ys.begin() + ys.size() / 2, ys.end());
    const int spine_x = xs[xs.size() / 2];
    const int spine_y = ys[ys.size() / 2];
    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(nodes.size() * 4);
    for (const RouteNode& node : nodes) {
      connectNodeToMedianSpine(rebuilt, node, spine_x, spine_y, trunk_layer);
    }

    // Reinforce central trunks so detailed routing has a strong connected backbone.
    appendSegment(rebuilt, min_x, spine_y, trunk_layer, max_x, spine_y, trunk_layer);
    appendSegment(rebuilt, spine_x, min_y, trunk_layer, spine_x, max_y, trunk_layer);

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }

    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

void applyDualBackboneWarp(GlobalRouter* grouter,
                           NetRouteMap& routes,
                           const RudyGrid& normalized_rudy,
                           int min_unique_nodes,
                           int coverage_percent,
                           int min_layer,
                           int max_layer)
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
  const int x_grids = std::max(grid->getXGrids(), 1);
  const int y_grids = std::max(grid->getYGrids(), 1);
  const int x_tiles = normalized_rudy.empty()
                          ? x_grids
                          : static_cast<int>(normalized_rudy.size());
  const int y_tiles = normalized_rudy.empty()
                          ? y_grids
                          : static_cast<int>(normalized_rudy.front().size());

  min_unique_nodes = std::max(min_unique_nodes, 3);
  coverage_percent = std::clamp(coverage_percent, 1, 100);
  min_layer = std::max(min_layer, 0);
  max_layer = std::max(max_layer, min_layer);

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

  std::vector<float> col_cost(std::max(x_tiles, 1), 0.0f);
  std::vector<float> row_cost(std::max(y_tiles, 1), 0.0f);
  if (!normalized_rudy.empty() && !normalized_rudy.front().empty()) {
    for (int x = 0; x < x_tiles; ++x) {
      for (int y = 0; y < y_tiles; ++y) {
        const float rudy = normalized_rudy[x][y];
        col_cost[x] += rudy;
        row_cost[y] += rudy;
      }
    }
    for (float& val : col_cost) {
      val /= static_cast<float>(std::max(y_tiles, 1));
    }
    for (float& val : row_cost) {
      val /= static_cast<float>(std::max(x_tiles, 1));
    }
  }

  auto pick_anchor = [](const std::vector<float>& cost,
                        int target,
                        int window) -> int {
    if (cost.empty()) {
      return 0;
    }
    const int size = static_cast<int>(cost.size());
    const int clamped_target = std::clamp(target, 0, size - 1);
    const int search = std::max(window, 1);
    const int lo = std::max(0, clamped_target - search);
    const int hi = std::min(size - 1, clamped_target + search);

    int best_idx = clamped_target;
    float best_score = std::numeric_limits<float>::max();
    for (int i = lo; i <= hi; ++i) {
      const float proximity_penalty
          = 0.10f
            * static_cast<float>(std::abs(i - clamped_target))
            / static_cast<float>(search + 1);
      const float score = cost[i] + proximity_penalty;
      if (score < best_score) {
        best_score = score;
        best_idx = i;
      }
    }
    return best_idx;
  };

  const int west_anchor_tile
      = pick_anchor(col_cost, x_tiles / 4, std::max(2, x_tiles / 5));
  const int east_anchor_tile = pick_anchor(
      col_cost, std::max((3 * x_tiles) / 4, 0), std::max(2, x_tiles / 5));
  const int south_anchor_tile
      = pick_anchor(row_cost, y_tiles / 4, std::max(2, y_tiles / 5));
  const int north_anchor_tile = pick_anchor(
      row_cost, std::max((3 * y_tiles) / 4, 0), std::max(2, y_tiles / 5));

  const int west_anchor_x
      = std::clamp(x_min + west_anchor_tile * tile, x_min, x_max);
  const int east_anchor_x
      = std::clamp(x_min + east_anchor_tile * tile, x_min, x_max);
  const int south_anchor_y
      = std::clamp(y_min + south_anchor_tile * tile, y_min, y_max);
  const int north_anchor_y
      = std::clamp(y_min + north_anchor_tile * tile, y_min, y_max);

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (static_cast<int>(key % 100ULL) >= coverage_percent) {
      continue;
    }

    const std::vector<RouteNode> nodes = collectUniqueRouteNodes(route);
    if (static_cast<int>(nodes.size()) < min_unique_nodes) {
      continue;
    }

    int min_x = nodes.front().x;
    int max_x = nodes.front().x;
    int min_y = nodes.front().y;
    int max_y = nodes.front().y;
    long centroid_x_acc = 0;
    long centroid_y_acc = 0;
    for (const RouteNode& node : nodes) {
      min_x = std::min(min_x, node.x);
      max_x = std::max(max_x, node.x);
      min_y = std::min(min_y, node.y);
      max_y = std::max(max_y, node.y);
      centroid_x_acc += node.x;
      centroid_y_acc += node.y;
    }
    const int node_count = std::max(static_cast<int>(nodes.size()), 1);
    const int centroid_x = static_cast<int>(centroid_x_acc / node_count);
    const int centroid_y = static_cast<int>(centroid_y_acc / node_count);

    const int backbone_x = (std::abs(centroid_x - west_anchor_x)
                            <= std::abs(centroid_x - east_anchor_x))
                               ? west_anchor_x
                               : east_anchor_x;
    const int backbone_y = (std::abs(centroid_y - south_anchor_y)
                            <= std::abs(centroid_y - north_anchor_y))
                               ? south_anchor_y
                               : north_anchor_y;
    const int trunk_layer = chooseDominantLayer(nodes, min_layer, max_layer);

    std::vector<GSegment> rebuilt;
    rebuilt.reserve(nodes.size() * 5 + 8);
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
      int cur_x = std::clamp(nodes[idx].x, x_min, x_max);
      int cur_y = std::clamp(nodes[idx].y, y_min, y_max);
      int cur_layer = nodes[idx].layer;

      if (cur_layer != trunk_layer) {
        const int step = (trunk_layer > cur_layer) ? 1 : -1;
        while (cur_layer != trunk_layer) {
          const int next_layer = cur_layer + step;
          appendSegment(
              rebuilt, cur_x, cur_y, cur_layer, cur_x, cur_y, next_layer);
          cur_layer = next_layer;
        }
      }

      const float h_rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(cur_x),
                                            to_grid_y(cur_y),
                                            to_grid_x(backbone_x),
                                            to_grid_y(backbone_y),
                                            true);
      const float v_rudy = estimatePathRudy(normalized_rudy,
                                            to_grid_x(cur_x),
                                            to_grid_y(cur_y),
                                            to_grid_x(backbone_x),
                                            to_grid_y(backbone_y),
                                            false);
      bool horizontal_first = h_rudy <= v_rudy;
      if (std::fabs(h_rudy - v_rudy) < 0.02f) {
        horizontal_first
            = ((key + static_cast<std::uint64_t>(idx) * 13ULL) & 1ULL) == 0ULL;
      }
      if (horizontal_first) {
        appendSegment(rebuilt,
                      cur_x,
                      cur_y,
                      trunk_layer,
                      backbone_x,
                      cur_y,
                      trunk_layer);
        appendSegment(rebuilt,
                      backbone_x,
                      cur_y,
                      trunk_layer,
                      backbone_x,
                      backbone_y,
                      trunk_layer);
      } else {
        appendSegment(rebuilt,
                      cur_x,
                      cur_y,
                      trunk_layer,
                      cur_x,
                      backbone_y,
                      trunk_layer);
        appendSegment(rebuilt,
                      cur_x,
                      backbone_y,
                      trunk_layer,
                      backbone_x,
                      backbone_y,
                      trunk_layer);
      }
    }

    appendSegment(rebuilt,
                  std::clamp(min_x, x_min, x_max),
                  backbone_y,
                  trunk_layer,
                  std::clamp(max_x, x_min, x_max),
                  backbone_y,
                  trunk_layer);
    appendSegment(rebuilt,
                  backbone_x,
                  std::clamp(min_y, y_min, y_max),
                  trunk_layer,
                  backbone_x,
                  std::clamp(max_y, y_min, y_max),
                  trunk_layer);
    appendSegment(rebuilt,
                  west_anchor_x,
                  backbone_y,
                  trunk_layer,
                  east_anchor_x,
                  backbone_y,
                  trunk_layer);
    appendSegment(rebuilt,
                  backbone_x,
                  south_anchor_y,
                  trunk_layer,
                  backbone_x,
                  north_anchor_y,
                  trunk_layer);

    std::vector<GSegment> compressed;
    compressed.reserve(rebuilt.size());
    for (const GSegment& segment : rebuilt) {
      appendCompressedSegment(compressed, segment);
    }
    if (!compressed.empty()) {
      route.swap(compressed);
    }
  }
}

void applyWavefrontDetours(GlobalRouter* grouter,
                           NetRouteMap& routes,
                           const RudyGrid& normalized_rudy,
                           int min_segment_len,
                           int detour_amplitude,
                           int coverage_percent)
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
  const int x_tiles = normalized_rudy.empty() ? 0 : normalized_rudy.size();
  const int y_tiles
      = normalized_rudy.empty() ? 0 : normalized_rudy.front().size();

  min_segment_len = std::max(min_segment_len, tile * 2);
  detour_amplitude = std::max(detour_amplitude, tile);
  coverage_percent = std::clamp(coverage_percent, 1, 100);

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

  auto choose_detour_coord = [&](int straight_coord,
                                 int sample_coord_a,
                                 int sample_coord_b,
                                 bool horizontal,
                                 std::uint64_t key) {
    const int low = horizontal ? y_min + tile : x_min + tile;
    const int high = horizontal ? y_max - tile : x_max - tile;
    if (low >= high) {
      return straight_coord;
    }

    const int pos = std::clamp(straight_coord + detour_amplitude, low, high);
    const int neg = std::clamp(straight_coord - detour_amplitude, low, high);
    if (pos == straight_coord && neg == straight_coord) {
      return straight_coord;
    }

    const int mid_coord = (sample_coord_a + sample_coord_b) / 2;
    const int gx_pos = horizontal ? to_grid_x(mid_coord) : to_grid_x(pos);
    const int gy_pos = horizontal ? to_grid_y(pos) : to_grid_y(mid_coord);
    const int gx_neg = horizontal ? to_grid_x(mid_coord) : to_grid_x(neg);
    const int gy_neg = horizontal ? to_grid_y(neg) : to_grid_y(mid_coord);
    const float pos_rudy = sampleRudyAt(normalized_rudy, gx_pos, gy_pos);
    const float neg_rudy = sampleRudyAt(normalized_rudy, gx_neg, gy_neg);

    if (std::fabs(pos_rudy - neg_rudy) < 0.04f) {
      return ((key >> 2U) & 1ULL) == 0ULL ? pos : neg;
    }
    return pos_rudy < neg_rudy ? pos : neg;
  };

  for (auto& [db_net, route] : routes) {
    if (route.empty()) {
      continue;
    }

    const std::uint64_t net_key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    std::vector<GSegment> transformed;
    transformed.reserve(route.size() * 4);

    for (size_t seg_idx = 0; seg_idx < route.size(); ++seg_idx) {
      const GSegment& segment = route[seg_idx];
      const bool planar
          = !segment.isVia() && segment.init_layer == segment.final_layer;
      const bool horizontal = planar && segment.init_y == segment.final_y
                              && segment.init_x != segment.final_x;
      const bool vertical = planar && segment.init_x == segment.final_x
                            && segment.init_y != segment.final_y;
      const long length = segment.length();
      const std::uint64_t key
          = net_key + static_cast<std::uint64_t>(seg_idx) * 131ULL;
      const bool selected = static_cast<int>(key % 100ULL) < coverage_percent;
      const bool candidate
          = selected && (horizontal || vertical) && length >= min_segment_len;
      if (!candidate) {
        transformed.push_back(segment);
        continue;
      }

      if (horizontal) {
        const int x0 = segment.init_x;
        const int x1 = segment.final_x;
        const int y = segment.init_y;
        const int p1 = x0 + (x1 - x0) / 3;
        const int p2 = x0 + (2 * (x1 - x0)) / 3;
        if (p1 == x0 || p2 == x1 || p1 == p2) {
          transformed.push_back(segment);
          continue;
        }

        const int detour_y = choose_detour_coord(y, p1, p2, true, key);
        if (detour_y == y) {
          transformed.push_back(segment);
          continue;
        }

        appendSegment(transformed, x0, y, segment.init_layer, p1, y, segment.init_layer);
        appendSegment(
            transformed, p1, y, segment.init_layer, p1, detour_y, segment.init_layer);
        appendSegment(transformed,
                      p1,
                      detour_y,
                      segment.init_layer,
                      p2,
                      detour_y,
                      segment.init_layer);
        appendSegment(
            transformed, p2, detour_y, segment.init_layer, p2, y, segment.init_layer);
        appendSegment(transformed, p2, y, segment.init_layer, x1, y, segment.init_layer);
        continue;
      }

      const int x = segment.init_x;
      const int y0 = segment.init_y;
      const int y1 = segment.final_y;
      const int p1 = y0 + (y1 - y0) / 3;
      const int p2 = y0 + (2 * (y1 - y0)) / 3;
      if (p1 == y0 || p2 == y1 || p1 == p2) {
        transformed.push_back(segment);
        continue;
      }

      const int detour_x = choose_detour_coord(x, p1, p2, false, key);
      if (detour_x == x) {
        transformed.push_back(segment);
        continue;
      }

      appendSegment(transformed, x, y0, segment.init_layer, x, p1, segment.init_layer);
      appendSegment(
          transformed, x, p1, segment.init_layer, detour_x, p1, segment.init_layer);
      appendSegment(transformed,
                    detour_x,
                    p1,
                    segment.init_layer,
                    detour_x,
                    p2,
                    segment.init_layer);
      appendSegment(
          transformed, detour_x, p2, segment.init_layer, x, p2, segment.init_layer);
      appendSegment(transformed, x, p2, segment.init_layer, x, y1, segment.init_layer);
    }

    if (!transformed.empty()) {
      route.swap(transformed);
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
  long nets_taken_from_shortcuts = 0;
  long nets_taken_from_anisotropic = 0;
  long nets_taken_from_wirelength_hunter = 0;
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
  auto route_stats = [](const GRoute& route) {
    std::pair<long, long> stats{0, 0};
    for (const GSegment& segment : route) {
      if (segment.isVia()) {
        stats.second++;
      } else {
        stats.first += std::abs(segment.final_x - segment.init_x)
                       + std::abs(segment.final_y - segment.init_y);
      }
    }
    return stats;
  };
  auto route_score = [&](const GRoute& route, int scenario_overflow) {
    const auto [route_wl, route_vias] = route_stats(route);
    // Keep this score wirelength-forward so net blending can escape
    // FastRoute-equivalent local minima.
    const double via_weight = static_cast<double>(tile_size) * 0.85;
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

  ScenarioResult shortcut = selected;
  shortcut.name = "ortholine_shortcut";
  applyAggressiveDoglegShortcuts(shortcut.routes,
                                 std::max(30 * tile_size, 1),
                                 std::max(2 * tile_size, 1));
  applyViaExcursionCollapse(shortcut.routes, std::max(14 * tile_size, 1));
  applyGuideCompression(shortcut.routes, std::max(24 * tile_size, 1));
  shortcut.metrics = compute_metrics(shortcut.routes);

  for (const auto& [db_net, route] : selected.routes) {
    auto shortcut_it = shortcut.routes.find(db_net);
    if (shortcut_it == shortcut.routes.end()) {
      continue;
    }
    const GRoute& shortcut_route = shortcut_it->second;

    const bool base_valid = has_planar_guide(route);
    const bool shortcut_valid = has_planar_guide(shortcut_route);
    if (!shortcut_valid && base_valid) {
      continue;
    }
    if (shortcut_valid && !base_valid) {
      selected.routes[db_net] = shortcut_route;
      nets_taken_from_shortcuts++;
      continue;
    }

    const double base_score = route_score(route, selected.overflow);
    const double shortcut_score = route_score(shortcut_route, selected.overflow);
    const auto [base_wl, base_vias] = route_stats(route);
    const auto [shortcut_wl, shortcut_vias] = route_stats(shortcut_route);
    bool use_shortcut = shortcut_score + 1e-3 < base_score
                        || (shortcut_wl < base_wl && shortcut_vias <= base_vias * 2);
    if (!use_shortcut && shortcut_wl <= base_wl && shortcut_vias <= base_vias * 3) {
      use_shortcut = true;
    }
    if (!use_shortcut && shortcut_score <= base_score * 1.06) {
      const auto key
          = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
      use_shortcut = (key % 3ULL) == 0ULL;
    }
    if (use_shortcut) {
      selected.routes[db_net] = shortcut_route;
      nets_taken_from_shortcuts++;
    }
  }
  selected.name += "+ortholine";
  selected.metrics = compute_metrics(selected.routes);

  ScenarioResult anisotropic = selected;
  anisotropic.name = "anisotropic_escape";
  applyLayerHoppingDetours(grouter_,
                           anisotropic.routes,
                           baseline_rudy,
                           min_routing_layer,
                           max_routing_layer);
  applyBraidedDetourWeave(grouter_, anisotropic.routes, baseline_rudy);
  applyWavefrontDetours(grouter_,
                        anisotropic.routes,
                        baseline_rudy,
                        std::max(12 * tile_size, 1),
                        std::max(4 * tile_size, 1),
                        68);
  applyAggressiveDoglegShortcuts(anisotropic.routes,
                                 std::max(28 * tile_size, 1),
                                 std::max(2 * tile_size, 1));
  applyViaExcursionCollapse(anisotropic.routes, std::max(16 * tile_size, 1));
  applyGuideCompression(anisotropic.routes, std::max(20 * tile_size, 1));
  anisotropic.metrics = compute_metrics(anisotropic.routes);

  for (const auto& [db_net, route] : selected.routes) {
    auto anisotropic_it = anisotropic.routes.find(db_net);
    if (anisotropic_it == anisotropic.routes.end()) {
      continue;
    }
    const GRoute& anisotropic_route = anisotropic_it->second;

    const bool base_valid = has_planar_guide(route);
    const bool anisotropic_valid = has_planar_guide(anisotropic_route);
    if (!anisotropic_valid && base_valid) {
      continue;
    }
    if (anisotropic_valid && !base_valid) {
      selected.routes[db_net] = anisotropic_route;
      nets_taken_from_anisotropic++;
      continue;
    }

    const double base_score = route_score(route, selected.overflow);
    const double anisotropic_score
        = route_score(anisotropic_route, selected.overflow);
    const auto [base_wl, base_vias] = route_stats(route);
    const auto [anisotropic_wl, anisotropic_vias] = route_stats(anisotropic_route);
    bool use_anisotropic
        = anisotropic_score + 1e-3 < base_score
          || (anisotropic_wl < base_wl && anisotropic_vias <= base_vias * 2);
    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (!use_anisotropic && anisotropic_wl <= base_wl
        && anisotropic_vias <= base_vias * 3) {
      use_anisotropic = true;
    }
    if (!use_anisotropic && anisotropic_score <= base_score * 1.12) {
      use_anisotropic = (key % 2ULL) == 0ULL;
    }
    if (use_anisotropic) {
      selected.routes[db_net] = anisotropic_route;
      nets_taken_from_anisotropic++;
    }
  }
  selected.name += "+anisotropic";
  selected.metrics = compute_metrics(selected.routes);

  // Radical candidate: rebuild most eligible nets through
  // low-congestion portal/spine trunks, then re-select with
  // wirelength-priority acceptance.
  ScenarioResult wirelength_hunter = selected;
  wirelength_hunter.name = "wirelength_hunter";
  applyMedianSpineRebuild(wirelength_hunter.routes,
                          5,
                          94,
                          min_routing_layer,
                          max_routing_layer);
  applyGlobalPortalRebuild(grouter_,
                           wirelength_hunter.routes,
                           baseline_rudy,
                           4,
                           86,
                           min_routing_layer,
                           max_routing_layer);
  applyRmstTrunkRebuild(grouter_,
                        wirelength_hunter.routes,
                        baseline_rudy,
                        3,
                        74,
                        90,
                        min_routing_layer,
                        max_routing_layer);
  applyBipolarPortalBackboneRebuild(grouter_,
                                    wirelength_hunter.routes,
                                    baseline_rudy,
                                    5,
                                    180,
                                    92,
                                    min_routing_layer,
                                    max_routing_layer);
  applyAggressiveDoglegShortcuts(wirelength_hunter.routes,
                                 std::max(36 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyWavefrontDetours(grouter_,
                        wirelength_hunter.routes,
                        baseline_rudy,
                        std::max(22 * tile_size, 1),
                        std::max(2 * tile_size, 1),
                        61);
  applyGuideCompression(wirelength_hunter.routes, std::max(10 * tile_size, 1));
  applyViaExcursionCollapse(wirelength_hunter.routes, std::max(5 * tile_size, 1));
  wirelength_hunter.metrics = compute_metrics(wirelength_hunter.routes);

  for (const auto& [db_net, route] : selected.routes) {
    auto hunter_it = wirelength_hunter.routes.find(db_net);
    if (hunter_it == wirelength_hunter.routes.end()) {
      continue;
    }
    const GRoute& hunter_route = hunter_it->second;

    const bool base_valid = has_planar_guide(route);
    const bool hunter_valid = has_planar_guide(hunter_route);
    if (!hunter_valid && base_valid) {
      continue;
    }
    if (hunter_valid && !base_valid) {
      selected.routes[db_net] = hunter_route;
      nets_taken_from_wirelength_hunter++;
      continue;
    }

    const auto [base_wl, base_vias] = route_stats(route);
    const auto [hunter_wl, hunter_vias] = route_stats(hunter_route);
    const double base_score = route_score(route, selected.overflow);
    const double hunter_score = route_score(hunter_route, selected.overflow);
    bool use_hunter = false;
    if (hunter_wl < static_cast<long>(base_wl * 0.97)) {
      use_hunter = hunter_vias <= static_cast<long>(base_vias * 1.9 + 3);
    }
    if (!use_hunter && hunter_wl <= base_wl) {
      use_hunter = hunter_vias <= static_cast<long>(base_vias * 1.35 + 2);
    }
    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    if (!use_hunter && hunter_score <= base_score * 1.12) {
      use_hunter = (key % 4ULL) == 0ULL;
    }
    if (use_hunter) {
      selected.routes[db_net] = hunter_route;
      nets_taken_from_wirelength_hunter++;
    }
  }
  selected.name += "+wirehunter";
  selected.metrics = compute_metrics(selected.routes);

  // Radical move: replace chained unconditional post-rewrites with a
  // deterministic candidate tournament and strict acceptance gates.
  ScenarioResult radical_hyper = compact;
  radical_hyper.name = "radical_hyper";
  applyQuadrantPortalHypergraphRebuild(grouter_,
                                       radical_hyper.routes,
                                       baseline_rudy,
                                       4,
                                       220,
                                       96,
                                       min_routing_layer,
                                       max_routing_layer);
  applyBipolarPortalBackboneRebuild(grouter_,
                                    radical_hyper.routes,
                                    baseline_rudy,
                                    4,
                                    220,
                                    96,
                                    min_routing_layer,
                                    max_routing_layer);
  applyRmstTrunkRebuild(grouter_,
                        radical_hyper.routes,
                        baseline_rudy,
                        3,
                        90,
                        96,
                        min_routing_layer,
                        max_routing_layer);
  applyAggressiveDoglegShortcuts(radical_hyper.routes,
                                 std::max(28 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(radical_hyper.routes, std::max(10 * tile_size, 1));
  applyViaExcursionCollapse(radical_hyper.routes, std::max(5 * tile_size, 1));
  radical_hyper.metrics = compute_metrics(radical_hyper.routes);

  // Axial-force scenario: rebuild almost every eligible net through
  // median/portal/trunk structures to escape compact-route fixed points.
  ScenarioResult axial_force = compact;
  axial_force.name = "axial_force";
  applyMedianSpineRebuild(
      axial_force.routes, 3, 100, min_routing_layer, max_routing_layer);
  applyQuadrantPortalHypergraphRebuild(grouter_,
                                       axial_force.routes,
                                       baseline_rudy,
                                       3,
                                       4096,
                                       100,
                                       min_routing_layer,
                                       max_routing_layer);
  applyGlobalPortalRebuild(grouter_,
                           axial_force.routes,
                           baseline_rudy,
                           3,
                           100,
                           min_routing_layer,
                           max_routing_layer);
  applyRmstTrunkRebuild(grouter_,
                        axial_force.routes,
                        baseline_rudy,
                        3,
                        4096,
                        100,
                        min_routing_layer,
                        max_routing_layer);
  applyBipolarPortalBackboneRebuild(grouter_,
                                    axial_force.routes,
                                    baseline_rudy,
                                    3,
                                    4096,
                                    100,
                                    min_routing_layer,
                                    max_routing_layer);
  applyAggressiveDoglegShortcuts(axial_force.routes,
                                 std::max(42 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(axial_force.routes, std::max(8 * tile_size, 1));
  applyViaExcursionCollapse(axial_force.routes, std::max(4 * tile_size, 1));
  axial_force.metrics = compute_metrics(axial_force.routes);

  // Shockwave scenario: aggressively collapse most nets onto global
  // portal/spine topologies so we force a topology break from compact routes.
  ScenarioResult shockwave = compact;
  shockwave.name = "shockwave";
  applyMedianSpineRebuild(
      shockwave.routes, 2, 100, min_routing_layer, max_routing_layer);
  applyQuadrantPortalHypergraphRebuild(grouter_,
                                       shockwave.routes,
                                       baseline_rudy,
                                       2,
                                       4096,
                                       100,
                                       min_routing_layer,
                                       max_routing_layer);
  applyGlobalPortalRebuild(grouter_,
                           shockwave.routes,
                           baseline_rudy,
                           2,
                           100,
                           min_routing_layer,
                           max_routing_layer);
  applyBipolarPortalBackboneRebuild(grouter_,
                                    shockwave.routes,
                                    baseline_rudy,
                                    2,
                                    4096,
                                    100,
                                    min_routing_layer,
                                    max_routing_layer);
  applyRmstTrunkRebuild(grouter_,
                        shockwave.routes,
                        baseline_rudy,
                        2,
                        4096,
                        100,
                        min_routing_layer,
                        max_routing_layer);
  applyAggressiveDoglegShortcuts(shockwave.routes,
                                 std::max(48 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(shockwave.routes, std::max(7 * tile_size, 1));
  applyViaExcursionCollapse(shockwave.routes, std::max(3 * tile_size, 1));
  shockwave.metrics = compute_metrics(shockwave.routes);

  // Monorail scenario: force most large nets onto a shared median spine and
  // dual portal hubs so we intentionally break residual compact-route geometry.
  ScenarioResult monorail = compact;
  monorail.name = "monorail_spine";
  applyMedianSpineRebuild(
      monorail.routes, 2, 100, min_routing_layer, max_routing_layer);
  applyGlobalPortalRebuild(grouter_,
                           monorail.routes,
                           baseline_rudy,
                           2,
                           100,
                           min_routing_layer,
                           max_routing_layer);
  applyQuadrantPortalHypergraphRebuild(grouter_,
                                       monorail.routes,
                                       baseline_rudy,
                                       2,
                                       4096,
                                       100,
                                       min_routing_layer,
                                       max_routing_layer);
  applyAggressiveDoglegShortcuts(monorail.routes,
                                 std::max(56 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(monorail.routes, std::max(6 * tile_size, 1));
  applyViaExcursionCollapse(monorail.routes, std::max(3 * tile_size, 1));
  monorail.metrics = compute_metrics(monorail.routes);

  // Mesh-warp scenario: force eligible nets through low-RUDY dual backbones.
  ScenarioResult meshwarp = compact;
  meshwarp.name = "dual_backbone_warp";
  applyDualBackboneWarp(grouter_,
                        meshwarp.routes,
                        baseline_rudy,
                        3,
                        100,
                        min_routing_layer,
                        max_routing_layer);
  applyAggressiveDoglegShortcuts(meshwarp.routes,
                                 std::max(20 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(meshwarp.routes, std::max(9 * tile_size, 1));
  applyViaExcursionCollapse(meshwarp.routes, std::max(5 * tile_size, 1));
  meshwarp.metrics = compute_metrics(meshwarp.routes);

  struct CandidateEntry
  {
    const char* name;
    const ScenarioResult* scenario;
  };

  std::vector<CandidateEntry> candidates{
      {"compact", &compact},
      {"shortcut", &shortcut},
      {"anisotropic", &anisotropic},
      {"wirehunter", &wirelength_hunter},
      {"radical_hyper", &radical_hyper},
      {"axial_force", &axial_force},
      {"shockwave", &shockwave},
      {"monorail", &monorail},
      {"meshwarp", &meshwarp}};
  if (sculpted_available) {
    candidates.push_back({"sculpted", &sculpted});
  }

  auto objective = [&](const RouteMetrics& metrics, int overflow) {
    const double via_weight = static_cast<double>(tile_size) * 0.72;
    const double overflow_penalty
        = static_cast<double>(std::max(overflow, 0))
          * static_cast<double>(tile_size) * 18.0;
    return static_cast<double>(metrics.wirelength_dbu)
           + via_weight * static_cast<double>(metrics.via_count)
           + overflow_penalty;
  };

  auto scenario_admissible = [&](const ScenarioResult& scenario) {
    if (scenario.routes.empty()) {
      return false;
    }
    const bool wl_guard
        = static_cast<double>(scenario.metrics.wirelength_dbu)
          <= static_cast<double>(baseline.metrics.wirelength_dbu) * 3.00;
    const bool via_guard
        = static_cast<double>(scenario.metrics.via_count)
          <= static_cast<double>(baseline.metrics.via_count) * 4.00;
    return wl_guard && via_guard;
  };

  auto route_admissible = [&](const GRoute& candidate, const GRoute& reference) {
    const auto [cand_wl, cand_vias] = route_stats(candidate);
    const auto [ref_wl, ref_vias] = route_stats(reference);
    if (ref_wl <= 0) {
      return true;
    }
    const long wl_cap
        = std::max(ref_wl + static_cast<long>(8 * tile_size),
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_wl) * 2.80)));
    const long via_cap
        = std::max(ref_vias + 8L,
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_vias) * 5.00 + 8.0)));
    return cand_wl <= wl_cap && cand_vias <= via_cap;
  };

  auto route_admissible_relaxed = [&](const GRoute& candidate,
                                      const GRoute& reference) {
    const auto [cand_wl, cand_vias] = route_stats(candidate);
    const auto [ref_wl, ref_vias] = route_stats(reference);
    if (ref_wl <= 0) {
      return true;
    }
    const long wl_cap
        = std::max(ref_wl + static_cast<long>(12 * tile_size),
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_wl) * 1.85)));
    const long via_cap
        = std::max(ref_vias + 16L,
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_vias) * 3.40 + 8.0)));
    return cand_wl <= wl_cap && cand_vias <= via_cap;
  };

  auto route_admissible_loose = [&](const GRoute& candidate,
                                    const GRoute& reference) {
    const auto [cand_wl, cand_vias] = route_stats(candidate);
    const auto [ref_wl, ref_vias] = route_stats(reference);
    if (ref_wl <= 0) {
      return true;
    }
    const long wl_cap
        = std::max(ref_wl + static_cast<long>(20 * tile_size),
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_wl) * 2.40)));
    const long via_cap
        = std::max(ref_vias + 20L,
                   static_cast<long>(
                       std::ceil(static_cast<double>(ref_vias) * 4.20 + 16.0)));
    return cand_wl <= wl_cap && cand_vias <= via_cap;
  };

  auto route_objective = [&](const GRoute& route, int overflow) {
    const auto [route_wl, route_vias] = route_stats(route);
    const double via_weight = static_cast<double>(tile_size) * 0.72;
    const double overflow_penalty
        = static_cast<double>(std::max(overflow, 0))
          * static_cast<double>(tile_size) * 18.0;
    return static_cast<double>(route_wl)
           + via_weight * static_cast<double>(route_vias) + overflow_penalty;
  };

  selected = compact;
  selected.name = "deterministic_tournament";
  std::vector<long> tournament_picks(candidates.size(), 0);
  long axial_forced_nets = 0;
  long shock_forced_nets = 0;
  long monorail_forced_nets = 0;
  long meshwarp_forced_nets = 0;

  for (const auto& [db_net, base_route] : compact.routes) {
    const GRoute* best_route = &base_route;
    int best_index = 0;
    double best_score = route_objective(base_route, compact.overflow);
    bool best_planar = has_planar_guide(base_route);

    for (size_t idx = 1; idx < candidates.size(); ++idx) {
      const ScenarioResult& scenario = *candidates[idx].scenario;
      if (!scenario_admissible(scenario)) {
        continue;
      }

      const auto it = scenario.routes.find(db_net);
      if (it == scenario.routes.end()) {
        continue;
      }

      const GRoute& candidate_route = it->second;
      const bool candidate_planar = has_planar_guide(candidate_route);
      if (!candidate_planar && best_planar) {
        continue;
      }
      if (candidate_planar && !best_planar) {
        best_route = &candidate_route;
        best_index = static_cast<int>(idx);
        best_score = route_objective(candidate_route, scenario.overflow);
        best_planar = true;
        continue;
      }
      if (!route_admissible(candidate_route, base_route)) {
        continue;
      }

      const double candidate_score
          = route_objective(candidate_route, scenario.overflow);
      if (candidate_score + 1e-3 < best_score) {
        best_route = &candidate_route;
        best_index = static_cast<int>(idx);
        best_score = candidate_score;
        best_planar = candidate_planar;
      }
    }

    selected.routes[db_net] = *best_route;
    if (best_index >= 0
        && best_index < static_cast<int>(tournament_picks.size())) {
      tournament_picks[best_index]++;
    }
  }
  selected.metrics = compute_metrics(selected.routes);

  // Force a second-pass axial replacement on complex nets so the run
  // can break out of compact-route local minima.
  for (const auto& [db_net, current_route] : selected.routes) {
    const auto axial_it = axial_force.routes.find(db_net);
    if (axial_it == axial_force.routes.end()) {
      continue;
    }
    const GRoute& axial_route = axial_it->second;
    if (!has_planar_guide(axial_route)) {
      continue;
    }
    const int node_count
        = static_cast<int>(collectUniqueRouteNodes(current_route).size());
    if (node_count < 2) {
      continue;
    }
    if (!route_admissible_relaxed(axial_route, current_route)) {
      continue;
    }

    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    const auto [base_wl, base_vias] = route_stats(current_route);
    const auto [axial_wl, axial_vias] = route_stats(axial_route);
    const double base_score = route_objective(current_route, selected.overflow);
    const double axial_score
        = route_objective(axial_route, axial_force.overflow);

    bool use_axial = axial_wl < static_cast<long>(base_wl * 0.985);
    if (!use_axial
        && axial_score <= base_score * (node_count >= 14 ? 1.10 : 1.06)) {
      use_axial = (key % 3ULL) == 0ULL;
    }
    if (!use_axial && node_count >= 16
        && axial_wl <= static_cast<long>(base_wl * 1.22)
        && axial_vias <= static_cast<long>(base_vias * 2.90 + 6)) {
      use_axial = (key % 2ULL) == 0ULL;
    }
    if (!use_axial && node_count >= 2
        && axial_wl <= static_cast<long>(base_wl * 2.50)
        && axial_vias <= static_cast<long>(base_vias * 6.00 + 20)) {
      use_axial = (key % 5ULL) == 0ULL;
    }

    if (use_axial) {
      selected.routes[db_net] = axial_route;
      axial_forced_nets++;
    }
  }
  if (axial_forced_nets > 0) {
    selected.name += "+axial_force";
    selected.metrics = compute_metrics(selected.routes);
  }

  // Final shock injection: force topology perturbation on a subset of
  // complex nets so NEWGR keeps escaping compact-route fixed points.
  for (const auto& [db_net, current_route] : selected.routes) {
    const auto shock_it = shockwave.routes.find(db_net);
    if (shock_it == shockwave.routes.end()) {
      continue;
    }
    const GRoute& shock_route = shock_it->second;
    if (!has_planar_guide(shock_route)) {
      continue;
    }

    const int node_count
        = static_cast<int>(collectUniqueRouteNodes(current_route).size());
    if (node_count < 7) {
      continue;
    }
    if (!route_admissible_loose(shock_route, current_route)) {
      continue;
    }

    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));
    const auto [base_wl, base_vias] = route_stats(current_route);
    const auto [shock_wl, shock_vias] = route_stats(shock_route);
    const double base_score = route_objective(current_route, selected.overflow);
    const double shock_score = route_objective(shock_route, shockwave.overflow);

    bool use_shock = false;
    if (node_count >= 14 && (key % 2ULL) == 1ULL) {
      use_shock = true;
    }
    if (!use_shock
        && shock_wl <= static_cast<long>(base_wl * 1.05)
        && shock_vias <= static_cast<long>(base_vias * 1.30 + 6)) {
      use_shock = true;
    }
    if (!use_shock && shock_score <= base_score * 1.18) {
      use_shock = (key % 3ULL) == 1ULL;
    }
    if (!use_shock && node_count >= 10
        && shock_wl <= static_cast<long>(base_wl * 1.35)
        && shock_vias <= static_cast<long>(base_vias * 3.80 + 12)) {
      use_shock = (key % 5ULL) == 2ULL;
    }
    if (use_shock) {
      selected.routes[db_net] = shock_route;
      shock_forced_nets++;
    }
  }
  if (shock_forced_nets > 0) {
    selected.name += "+shockwave";
    selected.metrics = compute_metrics(selected.routes);
  }

  // Monorail injection: apply a whole-net spine rewrite on complex nets.
  for (const auto& [db_net, current_route] : selected.routes) {
    const auto monorail_it = monorail.routes.find(db_net);
    if (monorail_it == monorail.routes.end()) {
      continue;
    }
    const GRoute& monorail_route = monorail_it->second;
    if (!has_planar_guide(monorail_route)) {
      continue;
    }

    const int node_count
        = static_cast<int>(collectUniqueRouteNodes(current_route).size());
    if (node_count < 8) {
      continue;
    }
    if (!route_admissible_loose(monorail_route, current_route)) {
      continue;
    }

    const auto [base_wl, base_vias] = route_stats(current_route);
    const auto [mono_wl, mono_vias] = route_stats(monorail_route);
    const double base_score = route_objective(current_route, selected.overflow);
    const double mono_score = route_objective(monorail_route, monorail.overflow);
    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));

    bool use_monorail = false;
    if (mono_wl < static_cast<long>(base_wl * 0.98)) {
      use_monorail = mono_vias <= static_cast<long>(base_vias * 2.40 + 8);
    }
    if (!use_monorail && node_count >= 14
        && mono_wl <= static_cast<long>(base_wl * 1.25)
        && mono_vias <= static_cast<long>(base_vias * 4.20 + 16)) {
      use_monorail = (key % 2ULL) == 0ULL;
    }
    if (!use_monorail && mono_score <= base_score * 1.16) {
      use_monorail = (key % 3ULL) == 1ULL;
    }
    if (use_monorail) {
      selected.routes[db_net] = monorail_route;
      monorail_forced_nets++;
    }
  }
  if (monorail_forced_nets > 0) {
    selected.name += "+monorail";
    selected.metrics = compute_metrics(selected.routes);
  }

  // Mesh-warp injection: force a dual-backbone topology shift on a subset
  // of medium/large nets to break repeated compact-route fixed points.
  for (const auto& [db_net, current_route] : selected.routes) {
    const auto mesh_it = meshwarp.routes.find(db_net);
    if (mesh_it == meshwarp.routes.end()) {
      continue;
    }
    const GRoute& mesh_route = mesh_it->second;
    if (!has_planar_guide(mesh_route)) {
      continue;
    }

    const int node_count
        = static_cast<int>(collectUniqueRouteNodes(current_route).size());
    if (node_count < 6) {
      continue;
    }
    if (!route_admissible_loose(mesh_route, current_route)) {
      continue;
    }

    const auto [base_wl, base_vias] = route_stats(current_route);
    const auto [mesh_wl, mesh_vias] = route_stats(mesh_route);
    const double base_score = route_objective(current_route, selected.overflow);
    const double mesh_score = route_objective(mesh_route, meshwarp.overflow);
    const auto key
        = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(db_net));

    bool use_mesh = false;
    if (mesh_wl < static_cast<long>(base_wl * 0.99)) {
      use_mesh = mesh_vias <= static_cast<long>(base_vias * 2.80 + 10);
    }
    if (!use_mesh && node_count >= 12
        && mesh_wl <= static_cast<long>(base_wl * 1.30)
        && mesh_vias <= static_cast<long>(base_vias * 4.40 + 16)) {
      use_mesh = (key % 2ULL) == 0ULL;
    }
    if (!use_mesh && mesh_score <= base_score * 1.20) {
      use_mesh = (key % 4ULL) == 1ULL;
    }
    if (!use_mesh && node_count >= 8
        && mesh_wl <= static_cast<long>(base_wl * 1.60)
        && mesh_vias <= static_cast<long>(base_vias * 5.20 + 20)) {
      use_mesh = (key % 7ULL) == 3ULL;
    }
    if (use_mesh) {
      selected.routes[db_net] = mesh_route;
      meshwarp_forced_nets++;
    }
  }
  if (meshwarp_forced_nets > 0) {
    selected.name += "+meshwarp";
    selected.metrics = compute_metrics(selected.routes);
  }

  ScenarioResult stabilized = selected;
  stabilized.name = "stabilized_radical_hyper";
  applyQuadrantPortalHypergraphRebuild(grouter_,
                                       stabilized.routes,
                                       baseline_rudy,
                                       5,
                                       180,
                                       62,
                                       min_routing_layer,
                                       max_routing_layer);
  applyAggressiveDoglegShortcuts(stabilized.routes,
                                 std::max(16 * tile_size, 1),
                                 std::max(tile_size, 1));
  applyGuideCompression(stabilized.routes, std::max(10 * tile_size, 1));
  applyViaExcursionCollapse(stabilized.routes, std::max(6 * tile_size, 1));
  stabilized.metrics = compute_metrics(stabilized.routes);

  if (scenario_admissible(stabilized)
      && objective(stabilized.metrics, selected.overflow) + 1e-3
             < objective(selected.metrics, selected.overflow)) {
    selected = stabilized;
  }

  // Final radical move for this iteration: force a global wave detour pass on
  // the selected routes so wirelength escapes compact-route fixed points.
  applyWavefrontDetours(grouter_,
                        selected.routes,
                        baseline_rudy,
                        std::max(2 * tile_size, 1),
                        std::max(6 * tile_size, 1),
                        100);
  applyViaExcursionCollapse(selected.routes, std::max(2 * tile_size, 1));
  selected.name += "+final_waveforce";
  selected.metrics = compute_metrics(selected.routes);

  // Final safeguard to prevent catastrophic regressions.
  const bool catastrophic
      = static_cast<double>(selected.metrics.wirelength_dbu)
            > static_cast<double>(baseline.metrics.wirelength_dbu) * 3.50
        || static_cast<double>(selected.metrics.via_count)
               > static_cast<double>(baseline.metrics.via_count) * 6.00;
  if (catastrophic) {
    selected = compact;
    selected.name = "catastrophic_fallback_compact";
  }
  selected.metrics = compute_metrics(selected.routes);

  const double selected_delta_wl
      = selected.metrics.wirelength_um - baseline.metrics.wirelength_um;
  const long selected_delta_vias
      = selected.metrics.via_count - baseline.metrics.via_count;
  logger_->warn(GNR,
                6032,
                "NEWGR deterministic tournament selected {}: baseline {:.0f} um/"
                "{} vias -> {:.0f} um/{} vias (delta wl {:+.0f} um, delta "
                "vias {:+d}).",
                selected.name,
                baseline.metrics.wirelength_um,
                baseline.metrics.via_count,
                selected.metrics.wirelength_um,
                selected.metrics.via_count,
                selected_delta_wl,
                selected_delta_vias);
  logger_->warn(GNR,
                6033,
                "NEWGR tournament picks: compact {} shortcut {} anisotropic {} "
                "wirehunter {} radical_hyper {} axial_force {} shockwave {} "
                "monorail {} meshwarp {} sculpted {}.",
                tournament_picks.size() > 0 ? tournament_picks[0] : 0,
                tournament_picks.size() > 1 ? tournament_picks[1] : 0,
                tournament_picks.size() > 2 ? tournament_picks[2] : 0,
                tournament_picks.size() > 3 ? tournament_picks[3] : 0,
                tournament_picks.size() > 4 ? tournament_picks[4] : 0,
                tournament_picks.size() > 5 ? tournament_picks[5] : 0,
                tournament_picks.size() > 6 ? tournament_picks[6] : 0,
                tournament_picks.size() > 7 ? tournament_picks[7] : 0,
                tournament_picks.size() > 8 ? tournament_picks[8] : 0,
                tournament_picks.size() > 9 ? tournament_picks[9] : 0);
  logger_->warn(GNR,
                6034,
                "NEWGR blend counters: sculpted {} shortcuts {} anisotropic {} "
                "wirehunter {} axial_force {} shockwave {} monorail {} "
                "meshwarp {}.",
                nets_taken_from_sculpted,
                nets_taken_from_shortcuts,
                nets_taken_from_anisotropic,
                nets_taken_from_wirelength_hunter,
                axial_forced_nets,
                shock_forced_nets,
                monorail_forced_nets,
                meshwarp_forced_nets);

  restore_snapshot(snapshot);
  return selected.routes;
}

}  // namespace grt
