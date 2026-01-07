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
  bool affect_horizontal = false;
  bool affect_vertical = false;
};

using RudyGrid = std::vector<std::vector<float>>;

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
  ratio = std::clamp(ratio, 0.05f, 1.0f);
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
                              int max_layer)
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
    const float min_ratio = 0.50f + 0.15f * layer_factor;
    const float max_ratio = 0.92f + 0.05f * layer_factor;

    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        const float ratio = logistic_ratio(normalized, 6.0f, 0.45f, min_ratio, max_ratio);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        const float ratio = logistic_ratio(normalized, 6.0f, 0.45f, min_ratio, max_ratio);
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
                           float base_ratio)
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

    for (const Hotspot& hotspot : hotspots) {
      for (int dx = -halo; dx <= halo; ++dx) {
        for (int dy = -halo; dy <= halo; ++dy) {
          const int gx = hotspot.gx + dx;
          const int gy = hotspot.gy + dy;
          if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
            continue;
          }
          if (hotspot.affect_horizontal && gx < x_grids - 1) {
            adjustEdgeCapacity(grouter, gx, gy, gx + 1, gy, layer, layer_ratio);
          }
          if (hotspot.affect_vertical && gy < y_grids - 1) {
            adjustEdgeCapacity(grouter, gx, gy, gx, gy + 1, layer, layer_ratio);
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
          * 5.0;
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
        const int gx
            = std::clamp((info.segment.init_x - x_min) / tile_size, 0, x_grids);
        const int gy
            = std::clamp((info.segment.init_y - y_min) / tile_size, 0, y_grids);
        Hotspot hotspot;
        hotspot.gx = std::clamp(gx, 0, std::max(x_grids - 1, 0));
        hotspot.gy = std::clamp(gy, 0, std::max(y_grids - 1, 0));
        hotspot.affect_vertical = is_vertical;
        hotspot.affect_horizontal = !is_vertical;
        hotspots.push_back(hotspot);
      }
    };

    append_hotspots(horizontal, false);
    append_hotspots(vertical, true);
    return hotspots;
  };

  auto adjust_edge_capacity
      = [&](int x1, int y1, int x2, int y2, int layer, float ratio) {
          ratio = std::clamp(ratio, 0.05f, 1.0f);
          const int current_cap
              = grouter_->fastroute_->getEdgeCapacity(x1, y1, x2, y2, layer);
          if (current_cap <= 0) {
            return;
          }
          const int new_cap = std::max(
              1, static_cast<int>(std::floor(current_cap * ratio)));
          if (new_cap == current_cap) {
            return;
          }
          const bool is_reduce = new_cap < current_cap;
          grouter_->fastroute_->addAdjustment(
              x1, y1, x2, y2, layer, new_cap, is_reduce);
        };

  auto apply_soft_capacity = [&](const RudyGrid& normalized_rudy) {
    if (normalized_rudy.empty() || grouter_->grid_ == nullptr) {
      return;
    }
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();
    const int usable_x = std::min<int>(x_grids, normalized_rudy.size());
    const int usable_y
        = std::min<int>(y_grids, normalized_rudy.front().size());
    if (usable_x == 0 || usable_y == 0) {
      return;
    }
    const int layer_span = std::max(max_routing_layer - min_routing_layer, 1);
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
    const auto normalized_at = [&](int x, int y) {
      if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
        return 0.0f;
      }
      return normalized_rudy[x][y];
    };

    for (int layer = min_routing_layer; layer <= max_routing_layer; ++layer) {
      const float layer_factor = static_cast<float>(layer - min_routing_layer)
                                 / static_cast<float>(layer_span);
      const float min_ratio = 0.50f + 0.15f * layer_factor;
      const float max_ratio = 0.92f + 0.05f * layer_factor;

      for (int y = 0; y < usable_y; ++y) {
        for (int x = 0; x < usable_x - 1; ++x) {
          const float normalized
              = 0.5f * (normalized_at(x, y) + normalized_at(x + 1, y));
          const float ratio
              = logistic_ratio(normalized, 6.0f, 0.45f, min_ratio, max_ratio);
          adjust_edge_capacity(x, y, x + 1, y, layer, ratio);
        }
      }
      for (int y = 0; y < usable_y - 1; ++y) {
        for (int x = 0; x < usable_x; ++x) {
          const float normalized
              = 0.5f * (normalized_at(x, y) + normalized_at(x, y + 1));
          const float ratio
              = logistic_ratio(normalized, 6.0f, 0.45f, min_ratio, max_ratio);
          adjust_edge_capacity(x, y, x, y + 1, layer, ratio);
        }
      }
    }
  };

  auto apply_hotspot_penalties = [&](const std::vector<Hotspot>& hotspots,
                                     int halo,
                                     float base_ratio) {
    if (hotspots.empty() || grouter_->grid_ == nullptr) {
      return;
    }
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();
    const int layer_span = std::max(max_routing_layer - min_routing_layer, 1);
    halo = std::max(0, halo);

    for (int layer = min_routing_layer; layer <= max_routing_layer; ++layer) {
      const float layer_factor
          = static_cast<float>(layer - min_routing_layer)
            / static_cast<float>(layer_span);
      const float layer_ratio
          = std::clamp(base_ratio + 0.1f * layer_factor, 0.4f, 0.95f);

      for (const Hotspot& hotspot : hotspots) {
        for (int dx = -halo; dx <= halo; ++dx) {
          for (int dy = -halo; dy <= halo; ++dy) {
            const int gx = hotspot.gx + dx;
            const int gy = hotspot.gy + dy;
            if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
              continue;
            }
            if (hotspot.affect_horizontal && gx < x_grids - 1) {
              adjust_edge_capacity(gx, gy, gx + 1, gy, layer, layer_ratio);
            }
            if (hotspot.affect_vertical && gy < y_grids - 1) {
              adjust_edge_capacity(gx, gy, gx, gy + 1, layer, layer_ratio);
            }
          }
        }
      }
    }
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

  Rudy* rudy = grouter_->getRudy();
  rudy->calculateRudy();
  RudyGrid normalized_rudy = computeNormalizedRudyGrid(rudy);

  std::vector<ScenarioResult> scenario_results;
  scenario_results.push_back(baseline);

  ScenarioDefinition baseline_def{"baseline", nullptr, nullptr};
  std::vector<ScenarioDefinition> scenario_defs;

  if (!normalized_rudy.empty()) {
    ScenarioDefinition soft_cap;
    soft_cap.name = "soft-cap";
    soft_cap.post_init = [this,
                          &apply_soft_capacity,
                          &apply_hotspot_penalties,
                          &normalized_rudy,
                          &hotspots,
                          min_routing_layer,
                          max_routing_layer]() {
      apply_soft_capacity(normalized_rudy);
      apply_hotspot_penalties(
          hotspots, 1, 0.65f);
    };
    scenario_defs.push_back(soft_cap);
  }

  auto make_random_def = [&](int seed, float perturb_pct) {
    ScenarioDefinition def;
    def.name = "perturb-seed" + std::to_string(seed);
    def.pre_init = [this, seed, perturb_pct]() {
      grouter_->setCapacitiesPerturbationPercentage(perturb_pct);
      grouter_->setPerturbationAmount(1);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(5.0f);
    };
    return def;
  };

  scenario_defs.push_back(make_random_def(11, 6.0f));
  scenario_defs.push_back(make_random_def(29, 4.0f));

  for (const ScenarioDefinition& def : scenario_defs) {
    ScenarioResult result = run_scenario(def, snapshot);
    scenario_results.push_back(std::move(result));
  }

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
                "NEWGR best scenario '{}': wirelength {:.0f} um, vias {}",
                final_result.name,
                final_result.metrics.wirelength_um,
                final_result.metrics.via_count);

  return std::move(final_result.routes);
}

}  // namespace grt
