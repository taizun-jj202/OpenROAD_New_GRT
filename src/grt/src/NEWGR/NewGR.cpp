#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstddef>
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
                  "hot edges {}, max ratio {:.2f}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow_edges,
                  metrics.near_capacity_edges,
                  metrics.max_usage_ratio);
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
                  "edges {}, hot edges {}, max ratio {:.2f}",
                  scenario.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.overflow_edges,
                  metrics.near_capacity_edges,
                  metrics.max_usage_ratio);
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

  if (!normalized_rudy.empty()) {
    scenario_defs.push_back(make_soft_config("cugr-prob-strong",
                                             0.38f,
                                             0.86f,
                                             8.0f,
                                             0.52f,
                                             2,
                                             0.52f,
                                             0.85f,
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
                                             2.5f,
                                             5,
                                             8.0f));
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

  for (const ScenarioDefinition& def : scenario_defs) {
    ScenarioResult result = run_scenario(def, snapshot);
    scenario_results.push_back(std::move(result));
  }

  auto robust_better = [](const ScenarioResult& lhs,
                          const ScenarioResult& rhs) {
    if (lhs.metrics.overflow_edges != rhs.metrics.overflow_edges) {
      return lhs.metrics.overflow_edges < rhs.metrics.overflow_edges;
    }
    if (lhs.metrics.near_capacity_edges != rhs.metrics.near_capacity_edges) {
      return lhs.metrics.near_capacity_edges < rhs.metrics.near_capacity_edges;
    }
    if (lhs.metrics.max_usage_ratio != rhs.metrics.max_usage_ratio) {
      return lhs.metrics.max_usage_ratio < rhs.metrics.max_usage_ratio;
    }
    if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
      return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
    }
    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    return lhs.metrics.score < rhs.metrics.score;
  };

  const auto shortest_wl_iter = std::min_element(
      scenario_results.begin(),
      scenario_results.end(),
      [](const ScenarioResult& lhs, const ScenarioResult& rhs) {
        return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
      });
  const long shortest_wl = shortest_wl_iter->metrics.wirelength_dbu;
  const long wl_guard_band
      = std::max<long>(200, static_cast<long>(std::ceil(0.0125 * shortest_wl)));

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

  auto best_ptr = *std::min_element(
      shortlist.begin(),
      shortlist.end(),
      [&](const ScenarioResult* lhs, const ScenarioResult* rhs) {
        return robust_better(*lhs, *rhs);
      });
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
                "overflow edges {}, hot edges {}, max ratio {:.2f}",
                final_result.name,
                final_result.metrics.wirelength_um,
                final_result.metrics.via_count,
                final_result.metrics.overflow_edges,
                final_result.metrics.near_capacity_edges,
                final_result.metrics.max_usage_ratio);

  return std::move(final_result.routes);
}

}  // namespace grt
