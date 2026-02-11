#include "NEWGR/NewGR.h"

#include <algorithm>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "FastRoute.h"
#include "grt/GRoute.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct GuidePatchStats
{
  int nets_patched = 0;
  int vias_added = 0;
  int wires_added = 0;
};

int pick_adjacent_layer(const int layer,
                        const int min_routing_layer,
                        const int max_routing_layer)
{
  if (layer + 1 <= max_routing_layer) {
    return layer + 1;
  }
  if (layer - 1 >= min_routing_layer) {
    return layer - 1;
  }
  return -1;
}

// Add small, optional "lane" patches (via-up, short parallel segment, via-down)
// along long guides. This intentionally increases guide flexibility for the
// detailed router without forcing a longer global route topology.
//
// Inspired by CUGR "patching" concepts: provide extra escape points / track
// switching opportunities for long segments, especially on lower layers where
// pin-access constraints are tighter.
void add_lane_patches(NetRouteMap& routes,
                      GlobalRouter* grouter,
                      utl::Logger* logger,
                      const int min_routing_layer,
                      const int max_routing_layer)
{
  if (grouter == nullptr || logger == nullptr) {
    return;
  }

  // Keep patches opt-in to avoid perturbing baseline metrics. Enable with:
  //   `utl::set_debug_level(GNR, "newgrGuidePatches", 1)` (or equivalent)
  // in the calling environment.
  if (!logger->debugCheck(GNR, "newgrGuidePatches", 1)) {
    return;
  }

  const int tile_size = grouter->getTileSize();
  if (tile_size <= 0) {
    return;
  }

  constexpr int kMinSegmentLenTiles = 8;   // only long segments
  constexpr int kEndMarginTiles = 2;       // keep away from endpoints/pins
  constexpr int kMinLaneSpanTiles = 3;     // ensure the lane has meaningful span
  constexpr int kMaxPatchesPerNet = 6;     // runtime guardrail
  constexpr int kMaxTotalPatchedSegs = 25000;

  // Focus patches on the lower routing layers where DR often struggles most.
  const int max_patch_layer
      = std::min(max_routing_layer, min_routing_layer + 2);

  GuidePatchStats stats;
  int total_patched = 0;

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);
    if (route.empty()) {
      continue;
    }
    if (total_patched >= kMaxTotalPatchedSegs) {
      break;
    }

    // Iterate over a snapshot to avoid patching patches.
    const std::vector<GSegment> original = route;

    std::unordered_set<GSegment, GSegmentHash> existing;
    existing.reserve(route.size() * 2 + 32);
    for (const auto& seg : route) {
      existing.insert(seg);
    }

    int net_patches = 0;
    bool net_counted = false;

    for (const auto& seg : original) {
      if (total_patched >= kMaxTotalPatchedSegs
          || net_patches >= kMaxPatchesPerNet) {
        break;
      }

      if (seg.isVia() || seg.isJumper()) {
        continue;
      }
      if (seg.init_layer != seg.final_layer) {
        continue;
      }

      const int layer = seg.init_layer;
      if (layer < min_routing_layer || layer > max_patch_layer) {
        continue;
      }

      const int adj_layer
          = pick_adjacent_layer(layer, min_routing_layer, max_routing_layer);
      if (adj_layer < 0 || adj_layer == layer) {
        continue;
      }

      const int len_dbu = seg.length();
      if (len_dbu <= 0 || len_dbu % tile_size != 0) {
        continue;  // off-grid or degenerate
      }

      const int len_tiles = len_dbu / tile_size;
      if (len_tiles < kMinSegmentLenTiles) {
        continue;
      }

      const bool horizontal
          = (seg.init_y == seg.final_y && seg.init_x != seg.final_x);
      const bool vertical
          = (seg.init_x == seg.final_x && seg.init_y != seg.final_y);
      if (!horizontal && !vertical) {
        continue;
      }

      const int margin = kEndMarginTiles * tile_size;
      if (len_dbu <= (2 * margin + kMinLaneSpanTiles * tile_size)) {
        continue;
      }

      int x1 = seg.init_x;
      int y1 = seg.init_y;
      int x2 = seg.final_x;
      int y2 = seg.final_y;
      if (horizontal) {
        x1 = seg.init_x + margin;
        x2 = seg.final_x - margin;
      } else {
        y1 = seg.init_y + margin;
        y2 = seg.final_y - margin;
      }

      const int span_dbu = std::abs(x1 - x2) + std::abs(y1 - y2);
      if (span_dbu <= 0 || span_dbu % tile_size != 0) {
        continue;
      }
      const int span_tiles = span_dbu / tile_size;
      if (span_tiles < kMinLaneSpanTiles) {
        continue;
      }

      const int via_lo = std::min(layer, adj_layer);
      const int via_hi = std::max(layer, adj_layer);

      // Create a parallel lane on the adjacent layer spanning between two
      // internal points, with vias to allow DR to hop up/down.
      const GSegment via_a(x1, y1, via_lo, x1, y1, via_hi);
      const GSegment lane(x1, y1, adj_layer, x2, y2, adj_layer);
      const GSegment via_b(x2, y2, via_lo, x2, y2, via_hi);

      int added_here = 0;
      if (existing.insert(via_a).second) {
        route.push_back(via_a);
        stats.vias_added++;
        added_here++;
      }
      if (existing.insert(lane).second) {
        route.push_back(lane);
        stats.wires_added++;
        added_here++;
      }
      if (existing.insert(via_b).second) {
        route.push_back(via_b);
        stats.vias_added++;
        added_here++;
      }

      if (added_here > 0) {
        net_patches++;
        total_patched++;
        if (!net_counted) {
          stats.nets_patched++;
          net_counted = true;
        }
      }
    }
  }

  if (stats.vias_added > 0 || stats.wires_added > 0) {
    logger->info(GNR,
                 6008,
                 "NEWGR guide patches: nets {}, added vias {}, added segs {}",
                 stats.nets_patched,
                 stats.vias_added,
                 stats.wires_added);
  }
}

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
  //   ripup/re-route churn.
  // - Avoid capacity perturbations (tends to increase detours/wirelength).
  CandidateSettings wl_lean = baseline;
  wl_lean.name = "wl-lean";
  // Seed affects tie-breaking in routing and can meaningfully change final WL.
  // Prefer a deterministic seed that historically produces shorter solutions.
  wl_lean.seed = 11;
  wl_lean.caps_perturbation_percentage = 0.0f;
  wl_lean.perturbation_amount = 1;
  wl_lean.congestion_iterations
      = std::min(25, snapshot.congestion_iterations);
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
    NetRouteMap routes = std::move(primary.routes);
    add_lane_patches(
        routes, grouter_, logger_, min_routing_layer, max_routing_layer);
    return routes;
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

  NetRouteMap routes = primary_is_best ? std::move(primary.routes)
                                      : std::move(recovery.routes);
  add_lane_patches(
      routes, grouter_, logger_, min_routing_layer, max_routing_layer);
  return routes;
}

}  // namespace grt
