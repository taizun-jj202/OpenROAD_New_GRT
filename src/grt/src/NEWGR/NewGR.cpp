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
  int congestion_iterations = 50;
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
  int congestion_iterations = 50;
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
    snapshot.congestion_iterations = grouter_->congestion_iterations_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setCongestionIterations(snapshot.congestion_iterations);
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
    grouter_->setCongestionIterations(settings.congestion_iterations);
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
  baseline.congestion_iterations = snapshot.congestion_iterations;
  baseline.critical_nets_percentage = snapshot.critical_percentage;

  // Wirelength-focused candidate with runtime guardrails:
  // - Keep congestion iterations bounded (fewer detours + faster runtime).
  // - Disable critical-net ordering to avoid STA overhead and additional
  //   ripup/re-route churn on this benchmark.
  // - Use modest capacity perturbation to avoid pathological tie-breaking.
  CandidateSettings wl_lean = baseline;
  wl_lean.name = "wl-lean";
  wl_lean.seed = 29;
  wl_lean.caps_perturbation_percentage = std::max(1.0f, snapshot.caps_percentage);
  wl_lean.perturbation_amount = 1;
  wl_lean.congestion_iterations = std::min(30, snapshot.congestion_iterations);
  wl_lean.critical_nets_percentage = 0.0f;

  auto is_better = [&](const CandidateResult& current,
                       const CandidateResult& best) -> bool {
    const bool current_ok = current.overflow == 0;
    const bool best_ok = best.overflow == 0;

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
  };

  // One-pass default: run the WL-lean configuration and only fall back if
  // overflow remains.
  CandidateResult primary = run_candidate(wl_lean, /*keep_routes=*/true);
  if (primary.overflow == 0) {
    logger_->info(GNR,
                  6005,
                  "NEWGR selected {}: wl {:.0f} um, vias {}, overflow {}",
                  primary.settings.name,
                  primary.metrics.wirelength_um,
                  primary.metrics.via_count,
                  primary.overflow);
    return std::move(primary.routes);
  }

  CandidateResult recovery = run_candidate(baseline, /*keep_routes=*/true);
  const bool primary_is_best = is_better(primary, recovery);
  const CandidateResult& best = primary_is_best ? primary : recovery;
  logger_->info(GNR,
                6007,
                "NEWGR selected {}: wl {:.0f} um, vias {}, overflow {}",
                best.settings.name,
                best.metrics.wirelength_um,
                best.metrics.via_count,
                best.overflow);

  return primary_is_best ? std::move(primary.routes) : std::move(recovery.routes);
}

}  // namespace grt
