#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "RoutingTracks.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
  int congestion_iterations = 0;
  float global_adjustment = 0.0f;
  std::vector<RegionAdjustment> region_adjustments;
};

struct CandidateConfig
{
  std::string name;
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  int seed = 0;
  float critical_percentage = 0.0f;
  int congestion_iterations = 0;
};

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

  struct RouteMetrics
  {
    long wirelength_dbu = 0;
    long via_count = 0;
    double wirelength_um = 0.0;
  };

  auto compute_metrics = [&](const NetRouteMap& routes) -> RouteMetrics {
    RouteMetrics metrics;
    for (const auto& [db_net, segments] : routes) {
      static_cast<void>(db_net);
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          metrics.via_count++;
        } else {
          metrics.wirelength_dbu += segment.length();
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
    snapshot.congestion_iterations = grouter_->congestion_iterations_;
    snapshot.global_adjustment = grouter_->adjustment_;
    snapshot.region_adjustments = grouter_->region_adjustments_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
    grouter_->setCongestionIterations(snapshot.congestion_iterations);
    grouter_->adjustment_ = snapshot.global_adjustment;
    grouter_->region_adjustments_ = snapshot.region_adjustments;
  };

  const auto prepare_fastroute = [&](const CandidateConfig& config) {
    // Re-initialize FastRoute state and grid/capacities without rebuilding
    // Net objects/pins (which is expensive). We reuse the `nets` list built
    // by the caller's initFastRoute() and only rebuild the FastRoute core
    // data structures, capacities, and adjustments.
    //
    // This keeps multi-candidate exploration affordable while still allowing
    // capacity perturbations/seed changes to take effect.

    grouter_->pad_pins_connections_.clear();

    grouter_->fastroute_->clear();
    grouter_->sproute_grid_ready_ = false;
    grouter_->sproute_nets_ready_ = false;
    grouter_->sproute_grid_data_ = SprouteGridData();
    grouter_->sproute_nets_.clear();
    grouter_->sproute_total_overflow_ = 0;
    grouter_->h_nets_in_pos_.clear();
    grouter_->v_nets_in_pos_.clear();

    grouter_->setCongestionIterations(config.congestion_iterations);
    grouter_->setCapacitiesPerturbationPercentage(config.caps_percentage);
    grouter_->setPerturbationAmount(config.perturbation_amount);
    grouter_->setSeed(config.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(config.critical_percentage);

    grouter_->routing_tracks_.clear();
    grouter_->routing_layers_.clear();
    grouter_->grid_->clear();
    grouter_->vertical_capacities_.clear();
    grouter_->horizontal_capacities_.clear();

    grouter_->ensureLayerForGuideDimension(max_routing_layer);
    grouter_->configFastRoute();
    grouter_->initRoutingLayers(min_routing_layer, max_routing_layer);
    grouter_->initRoutingTracks(max_routing_layer);
    grouter_->initCoreGrid(max_routing_layer);
    grouter_->setCapacities(min_routing_layer, max_routing_layer);
    grouter_->captureSprouteGridData();
    grouter_->applyAdjustments(min_routing_layer, max_routing_layer);
    grouter_->perturbCapacities();

    // Init the data structures to monitor 3D capacity during 2D phases
    grouter_->fastroute_->initEdgesCapacityPerLayer();

    // Rebuild the FastRoute netlist from the already-created Net objects.
    grouter_->initFastRouteIncr(nets);
    grouter_->initialized_ = true;
  };

  RouterSnapshot snapshot = capture_snapshot();
  // Strategy:
  // - Explore a small set of FastRoute configurations (seed/perturbation/
  //   critical nets %) and pick the best by (global) wirelength first, then
  //   by via count.
  // - Re-initialize FastRoute per candidate without rebuilding Net objects,
  //   keeping runtime reasonable.

  std::vector<CandidateConfig> candidates;
  candidates.push_back({"baseline",
                        snapshot.caps_percentage,
                        snapshot.perturbation_amount,
                        snapshot.seed,
                        snapshot.critical_percentage,
                        snapshot.congestion_iterations});
  candidates.push_back({"perturb-seed11-crit0",
                        6.0f,
                        1,
                        11,
                        0.0f,
                        snapshot.congestion_iterations});
  candidates.push_back({"perturb-seed11-crit5",
                        6.0f,
                        1,
                        11,
                        5.0f,
                        snapshot.congestion_iterations});
  candidates.push_back({"perturb-seed17-crit0",
                        6.0f,
                        1,
                        17,
                        0.0f,
                        snapshot.congestion_iterations});

  const auto better = [](const RouteMetrics& lhs, const RouteMetrics& rhs) {
    if (lhs.wirelength_dbu != rhs.wirelength_dbu) {
      return lhs.wirelength_dbu < rhs.wirelength_dbu;
    }
    return lhs.via_count < rhs.via_count;
  };

  NetRouteMap chosen_routes;
  RouteMetrics chosen_metrics;
  std::string chosen_name;
  bool have_choice = false;

  for (const auto& candidate : candidates) {
    restore_snapshot(snapshot);
    prepare_fastroute(candidate);

    NetRouteMap routes
        = grouter_->findRouting(nets, min_routing_layer, max_routing_layer);
    const RouteMetrics metrics = compute_metrics(routes);
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}",
                  candidate.name,
                  metrics.wirelength_um,
                  metrics.via_count);

    if (!have_choice || better(metrics, chosen_metrics)) {
      chosen_routes = std::move(routes);
      chosen_metrics = metrics;
      chosen_name = candidate.name;
      have_choice = true;
    }
  }

  logger_->info(GNR, 6007, "NEWGR picked {}", chosen_name);

  restore_snapshot(snapshot);
  return chosen_routes;
}

}  // namespace grt
