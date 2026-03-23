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
    RouteMetrics metrics = compute_metrics(routes);
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count);
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
                  "NEWGR scenario {}: wirelength {:.0f} um, vias {}",
                  scenario.name,
                  metrics.wirelength_um,
                  metrics.via_count);
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
  PlanarEdgeUsage baseline_usage
      = computeNormalizedBackboneUsage(grouter_, baseline.routes);

  ScenarioDefinition backbone_def;
  backbone_def.name = "backbone-channelized";
  backbone_def.pre_init = [this]() {
    grouter_->setCapacitiesPerturbationPercentage(10.0f);
    grouter_->setPerturbationAmount(1);
    grouter_->setSeed(29);
    grouter_->setAllowCongestion(false);
    grouter_->fastroute_->setCriticalNetsPercentage(28.0f);
  };
  backbone_def.post_init
      = [this,
         &normalized_rudy,
         &hotspots,
         &baseline_usage,
         min_routing_layer,
         max_routing_layer]() {
          applySoftCapacityScaling(grouter_,
                                   normalized_rudy,
                                   min_routing_layer,
                                   max_routing_layer,
                                   0.56f,
                                   0.97f,
                                   4.4f,
                                   0.43f);
          applyBackboneCapacityReinforcement(grouter_,
                                             normalized_rudy,
                                             baseline_usage,
                                             min_routing_layer,
                                             max_routing_layer,
                                             0.34f,
                                             1.55f,
                                             0.62f);
          applyAggressiveCapacityField(grouter_,
                                       normalized_rudy,
                                       hotspots,
                                       min_routing_layer,
                                       max_routing_layer,
                                       0.76f,
                                       0.12f,
                                       0.23f,
                                       1.38f,
                                       0.35f,
                                       true);
          applyHotspotPenalties(grouter_,
                                hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                4,
                                0.50f,
                                0.90f);
        };

  ScenarioResult channelized = run_scenario(backbone_def, snapshot);
  const int channelized_overflow
      = grouter_->fastroute() != nullptr ? grouter_->fastroute()->totalOverflow()
                                         : 0;
  const bool channelized_congested = channelized_overflow > 0;

  auto better_result = [](const ScenarioResult& lhs,
                          const ScenarioResult& rhs) {
    if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
      return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
    }
    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    return lhs.metrics.score < rhs.metrics.score;
  };

  const ScenarioResult& better
      = better_result(channelized, baseline) ? channelized : baseline;
  logger_->info(GNR,
                6007,
                "NEWGR diagnostic best between baseline/backbone is '{}': "
                "wirelength {:.0f} um, vias {}",
                better.name,
                better.metrics.wirelength_um,
                better.metrics.via_count);

  if (!channelized_congested && !channelized.routes.empty()) {
    restore_snapshot(snapshot);
    return std::move(channelized.routes);
  }

  logger_->warn(
      GNR,
      6008,
      "NEWGR backbone-channelized result has overflow {} ; replaying baseline "
      "to guarantee routable guides.",
      channelized_overflow);
  ScenarioDefinition baseline_replay{"baseline-replay", nullptr, nullptr};
  ScenarioResult fallback = run_scenario(baseline_replay, snapshot);
  restore_snapshot(snapshot);
  if (!fallback.routes.empty()) {
    return std::move(fallback.routes);
  }
  return std::move(baseline.routes);
}

}  // namespace grt
