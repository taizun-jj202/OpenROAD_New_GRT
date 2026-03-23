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
    // Even in low-overflow designs, run light RUDY-guided variants to expose
    // alternative low-detour topologies for detailed route.
    scenario_defs.push_back(make_soft_config("rudy-light-direct",
                                             0.70f,
                                             1.00f,
                                             3.2f,
                                             0.30f,
                                             1,
                                             0.88f,
                                             0.15f,
                                             0.10f,
                                             0.0f,
                                             snapshot.seed,
                                             2.0f));
    scenario_defs.push_back(make_soft_config("rudy-light-balanced",
                                             0.62f,
                                             0.98f,
                                             3.8f,
                                             0.34f,
                                             1,
                                             0.82f,
                                             0.20f,
                                             0.09f,
                                             0.0f,
                                             snapshot.seed,
                                             3.0f));
    logger_->info(GNR,
                  6008,
                  "NEWGR enabling light soft-capacity scenarios in low "
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

    const auto get_route_metrics = [](const GRoute& route) {
      long wl = 0;
      int vias = 0;
      for (const GSegment& segment : route) {
        if (segment.isVia()) {
          vias++;
        } else {
          wl += std::abs(segment.final_x - segment.init_x)
                + std::abs(segment.final_y - segment.init_y);
        }
      }
      return std::pair<long, int>{wl, vias};
    };

    auto append_hybrid = [&](const std::string& hybrid_name,
                             int source_count,
                             long via_weight,
                             int logger_code) {
      source_count = std::max(1, std::min(source_count, static_cast<int>(ranked.size())));
      ScenarioResult hybrid_result;
      hybrid_result.name = hybrid_name;

      for (odb::dbNet* db_net : hybrid_nets) {
        if (db_net == nullptr) {
          continue;
        }
        const GRoute* best_route = nullptr;
        long best_cost = std::numeric_limits<long>::max();
        long best_wl = std::numeric_limits<long>::max();
        int best_via = std::numeric_limits<int>::max();

        for (int idx = 0; idx < source_count; ++idx) {
          const ScenarioResult* source = ranked[idx];
          if (source == nullptr) {
            continue;
          }
          const auto route_it = source->routes.find(db_net);
          if (route_it == source->routes.end()) {
            continue;
          }
          const GRoute& source_route = route_it->second;
          const auto [wl, vias] = get_route_metrics(source_route);
          const long cost = wl + via_weight * static_cast<long>(vias);
          if (cost < best_cost || (cost == best_cost && wl < best_wl)
              || (cost == best_cost && wl == best_wl && vias < best_via)) {
            best_cost = cost;
            best_wl = wl;
            best_via = vias;
            best_route = &source_route;
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
                      "vias {}, routed nets {}/{}",
                      hybrid_name,
                      source_count,
                      hybrid_result.metrics.wirelength_um,
                      hybrid_result.metrics.via_count,
                      hybrid_result.routes.size(),
                      expected_net_count);
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

    // Drastic recombination:
    // 1) a wirelength-first hybrid (FastRoute shortest-path intent),
    // 2) a balanced hybrid that softly penalizes vias (SPRoute/CUGR flavor).
    const int wl_source_count = std::min<int>(14, ranked.size());
    const int balanced_source_count = std::min<int>(8, ranked.size());
    const long balanced_via_weight
        = std::max<long>(1, static_cast<long>(std::max(grouter_->grid_->getTileSize(), 1)) / 3L);
    append_hybrid("hybrid-netmix-wl", wl_source_count, 0, 6010);
    append_hybrid(
        "hybrid-netmix-balanced", balanced_source_count, balanced_via_weight, 6012);
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
  const double wl_guard_ratio = overflow_free_sweep ? 0.0035 : 0.0125;
  const long wl_guard_band
      = std::max<long>(200, static_cast<long>(std::ceil(wl_guard_ratio * shortest_wl)));

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
      = std::max<long>(24, static_cast<long>(std::ceil(shortest_wl * 0.00008)));
  auto wirelength_with_via_tie_better = [&](const ScenarioResult& lhs,
                                            const ScenarioResult& rhs) {
    const long wl_gap = std::llabs(lhs.metrics.wirelength_dbu
                                   - rhs.metrics.wirelength_dbu);
    if (wl_gap <= tie_via_wl_band && lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    return wirelength_first_better(lhs, rhs);
  };

  auto best_ptr = *std::min_element(
      shortlist.begin(),
      shortlist.end(),
      [&](const ScenarioResult* lhs, const ScenarioResult* rhs) {
        if (overflow_free_sweep) {
          return wirelength_with_via_tie_better(*lhs, *rhs);
        }
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
