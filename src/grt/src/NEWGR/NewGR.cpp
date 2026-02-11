#include "NEWGR/NewGR.h"

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "FastRoute.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouteMetrics
{
  long wirelength_dbu = 0;
  long via_count = 0;
  double wirelength_um = 0.0;
};

struct CandidateSettings
{
  std::string name;
  int seed = 0;
  float caps_perturbation_percentage = 0.0f;
  int perturbation_amount = 1;
  float critical_nets_percentage = 10.0f;
};

struct CandidateResult
{
  CandidateSettings settings;
  RouteMetrics metrics;
  int overflow = std::numeric_limits<int>::max();
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
  static_cast<void>(nets);

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

  RouterSnapshot snapshot = capture_snapshot();

  auto run_candidate = [&](const CandidateSettings& settings,
                           bool keep_routes) -> CandidateResult {
    CandidateResult result;
    result.settings = settings;

    restore_snapshot(snapshot);
    grouter_->setAllowCongestion(false);
    grouter_->setCapacitiesPerturbationPercentage(
        settings.caps_perturbation_percentage);
    grouter_->setPerturbationAmount(settings.perturbation_amount);
    grouter_->setSeed(settings.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        settings.critical_nets_percentage);

    std::vector<Net*> candidate_nets
        = grouter_->initFastRoute(min_routing_layer, max_routing_layer);
    NetRouteMap routes;
    if (!candidate_nets.empty()) {
      routes = grouter_->findRouting(
          candidate_nets, min_routing_layer, max_routing_layer);
    }

    result.metrics = compute_metrics(routes);
    result.overflow = (grouter_->fastroute_ != nullptr)
                          ? grouter_->fastroute_->totalOverflow()
                          : std::numeric_limits<int>::max();
    if (keep_routes) {
      result.routes = std::move(routes);
    }

    logger_->info(GNR,
                  6006,
                  "NEWGR candidate {}: wl {:.0f} um, vias {}, overflow {}",
                  result.settings.name,
                  result.metrics.wirelength_um,
                  result.metrics.via_count,
                  result.overflow);

    return result;
  };

  CandidateSettings baseline;
  baseline.name = "baseline";
  baseline.seed = snapshot.seed;
  baseline.caps_perturbation_percentage = snapshot.caps_percentage;
  baseline.perturbation_amount = snapshot.perturbation_amount;
  baseline.critical_nets_percentage = snapshot.critical_percentage;

  CandidateSettings perturb_11;
  perturb_11.name = "perturb-s11";
  perturb_11.seed = 11;
  perturb_11.caps_perturbation_percentage = 1.5f;
  perturb_11.perturbation_amount = 1;
  perturb_11.critical_nets_percentage = snapshot.critical_percentage;

  CandidateSettings perturb_29;
  perturb_29.name = "perturb-s29";
  perturb_29.seed = 29;
  perturb_29.caps_perturbation_percentage = 1.5f;
  perturb_29.perturbation_amount = 1;
  perturb_29.critical_nets_percentage = snapshot.critical_percentage;

  std::vector<CandidateSettings> candidates;
  candidates.push_back(baseline);
  candidates.push_back(perturb_11);
  candidates.push_back(perturb_29);

  CandidateResult best;
  for (const CandidateSettings& candidate : candidates) {
    CandidateResult current = run_candidate(candidate, /*keep_routes=*/false);
    const bool current_ok = current.overflow == 0;
    const bool best_ok = best.overflow == 0;

    const bool better = [&]() {
      if (current_ok != best_ok) {
        return current_ok;  // prefer overflow-free
      }
      if (!current_ok && current.overflow != best.overflow) {
        return current.overflow < best.overflow;
      }
      if (current.metrics.wirelength_dbu != best.metrics.wirelength_dbu) {
        return current.metrics.wirelength_dbu < best.metrics.wirelength_dbu;
      }
      return current.metrics.via_count < best.metrics.via_count;
    }();

    if (best.settings.name.empty() || better) {
      best = std::move(current);
    }
  }

  // Replay the winner so GlobalRouter state (congestion DB, etc.) matches the
  // returned routes.
  CandidateResult replay = run_candidate(best.settings, /*keep_routes=*/true);
  logger_->info(GNR,
                6005,
                "NEWGR selected {}: wl {:.0f} um, vias {}, overflow {}",
                replay.settings.name,
                replay.metrics.wirelength_um,
                replay.metrics.via_count,
                replay.overflow);

  return std::move(replay.routes);
}

}  // namespace grt
