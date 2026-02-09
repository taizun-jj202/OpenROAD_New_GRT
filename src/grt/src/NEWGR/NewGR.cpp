#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "RoutingTracks.h"
#include "grt/Rudy.h"
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
  int rudy_hotspots = 0;            // number of hotspots to patch (0 disables)
  int rudy_expand_tiles = 0;        // expand each hotspot by N tiles
  float rudy_adjustment = 1.0f;     // capacity multiplier (e.g. 0.8 = -20%)
  int rudy_layers = 0;              // apply to [min_layer, min_layer + L)
};

struct CongestionScore
{
  // 0.0 is ideal; higher means "tighter" routing.
  double score = 0.0;
  double max_util = 0.0;
  double frac_util_90 = 0.0;
  double frac_util_80 = 0.0;
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

	  const std::vector<Net*> canonical_nets = nets;

	  struct RouteMetrics
	  {
	    long wirelength_dbu = 0;
	    long via_count = 0;
	    double wirelength_um = 0.0;
	    int total_overflow = 0;
      CongestionScore congestion;
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

  const auto compute_congestion_score = [&]() -> CongestionScore {
    CongestionScore out;
    const auto& h_edges = grouter_->fastroute_->getHorizontalEdges3D();
    const auto& v_edges = grouter_->fastroute_->getVerticalEdges3D();

    auto scan_edges = [&](const auto& edges) {
      std::size_t edge_count = 0;
      std::size_t util_80 = 0;
      std::size_t util_90 = 0;

      for (const auto& layer : edges) {
        for (const auto& row : layer) {
          for (const auto& edge : row) {
            const int eff_cap = static_cast<int>(edge.cap)
                                - static_cast<int>(edge.red);
            if (eff_cap <= 0) {
              continue;
            }
            const double util = static_cast<double>(edge.usage)
                                / static_cast<double>(eff_cap);
            out.max_util = std::max(out.max_util, util);
            if (util >= 0.8) {
              util_80++;
            }
            if (util >= 0.9) {
              util_90++;
            }
            edge_count++;
          }
        }
      }

      if (edge_count > 0) {
        out.frac_util_80
            = std::max(out.frac_util_80,
                       static_cast<double>(util_80) / edge_count);
        out.frac_util_90
            = std::max(out.frac_util_90,
                       static_cast<double>(util_90) / edge_count);
      }
    };

    scan_edges(h_edges);
    scan_edges(v_edges);

    // A light "DR-friendliness" proxy: peak utilization dominates, with small
    // penalties when lots of edges are near-full.
    out.score = out.max_util + 0.10 * out.frac_util_90 + 0.03 * out.frac_util_80;
    return out;
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

  const RouterSnapshot snapshot = capture_snapshot();

  // Build a deterministic list of "hotspot" regions from Rudy (wire density)
  // to optionally reserve extra soft-capacity around congested areas, similar
  // in spirit to SPRoute soft capacity / CUGR patching.
  std::vector<odb::Rect> rudy_hotspot_regions;
  {
    Rudy* rudy = grouter_->getRudy();
    if (rudy != nullptr && grouter_->grid_ != nullptr) {
      const int x_grids = grouter_->grid_->getXGrids();
      const int y_grids = grouter_->grid_->getYGrids();
      if (x_grids > 0 && y_grids > 0) {
        rudy->setGridConfig(grouter_->grid_->getGridArea(), x_grids, y_grids);
        rudy->calculateRudy();

        struct TileRudy
        {
          int x = 0;
          int y = 0;
          float rudy = 0.0f;
        };
        std::vector<TileRudy> tiles;
        tiles.reserve(static_cast<std::size_t>(x_grids) * y_grids);

        for (int x = 0; x < x_grids; x++) {
          for (int y = 0; y < y_grids; y++) {
            const float value = rudy->getTile(x, y).getRudy();
            if (value <= 0.0f) {
              continue;
            }
            tiles.push_back({x, y, value});
          }
        }

        std::sort(tiles.begin(),
                  tiles.end(),
                  [](const TileRudy& a, const TileRudy& b) {
                    if (a.rudy != b.rudy) {
                      return a.rudy > b.rudy;
                    }
                    if (a.x != b.x) {
                      return a.x < b.x;
                    }
                    return a.y < b.y;
                  });

        // Keep the list small and stable; we'll take prefixes for candidates.
        const int max_regions
            = std::min<int>(80, static_cast<int>(tiles.size()));
        rudy_hotspot_regions.reserve(max_regions);
        for (int i = 0; i < max_regions; i++) {
          const auto& t = tiles[i];
          rudy_hotspot_regions.push_back(rudy->getTile(t.x, t.y).getRect());
        }
      }
    }
  }

  const auto add_rudy_hotspot_adjustments = [&](const CandidateConfig& config) {
    if (config.rudy_hotspots <= 0 || config.rudy_adjustment >= 1.0f
        || config.rudy_layers <= 0 || rudy_hotspot_regions.empty()) {
      return;
    }

    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();
    const int tile = grouter_->grid_->getTileSize();
    if (x_grids <= 0 || y_grids <= 0 || tile <= 0) {
      return;
    }

    const int max_hotspots
        = std::min<int>(config.rudy_hotspots,
                        static_cast<int>(rudy_hotspot_regions.size()));
    const int min_layer = min_routing_layer;
    const int max_layer = std::min(max_routing_layer,
                                   min_routing_layer + config.rudy_layers - 1);

    for (int i = 0; i < max_hotspots; i++) {
      odb::Rect rect = rudy_hotspot_regions[i];
      if (config.rudy_expand_tiles > 0) {
        const int expand = config.rudy_expand_tiles * tile;
        rect = odb::Rect(rect.xMin() - expand,
                         rect.yMin() - expand,
                         rect.xMax() + expand,
                         rect.yMax() + expand);
      }
      for (int layer = min_layer; layer <= max_layer; layer++) {
        grouter_->region_adjustments_.emplace_back(rect.xMin(),
                                                  rect.yMin(),
                                                  rect.xMax(),
                                                  rect.yMax(),
                                                  layer,
                                                  config.rudy_adjustment);
      }
    }
  };

	  const auto prepare_fastroute = [&](const CandidateConfig& config,
	                                     std::vector<Net*>& candidate_nets) {
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

    // Start from the snapshot's adjustments and optionally add Rudy-based
    // hotspot "soft-capacity" reservations for DR-friendliness.
    grouter_->region_adjustments_ = snapshot.region_adjustments;

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

    add_rudy_hotspot_adjustments(config);

    grouter_->applyAdjustments(min_routing_layer, max_routing_layer);
    grouter_->perturbCapacities();

	    // Init the data structures to monitor 3D capacity during 2D phases
	    grouter_->fastroute_->initEdgesCapacityPerLayer();

	    // Rebuild the FastRoute netlist from the already-created Net objects.
	    grouter_->initFastRouteIncr(candidate_nets);
	    grouter_->initialized_ = true;
	  };

  logger_->info(GNR,
                6004,
                "NEWGR snapshot: caps% {:.1f} pert_amt {} seed {} crit% {:.1f} cong_iters {} global_adj {:.2f}",
                snapshot.caps_percentage,
                snapshot.perturbation_amount,
                snapshot.seed,
                snapshot.critical_percentage,
                snapshot.congestion_iterations,
                snapshot.global_adjustment);

  const auto better = [](const RouteMetrics& lhs, const RouteMetrics& rhs) {
    const bool lhs_routable = lhs.total_overflow == 0;
    const bool rhs_routable = rhs.total_overflow == 0;
    if (lhs_routable != rhs_routable) {
      return lhs_routable;
    }
    if (!lhs_routable && lhs.total_overflow != rhs.total_overflow) {
      return lhs.total_overflow < rhs.total_overflow;
    }

    // Both solutions are routable. Global wirelength deltas across small
    // perturbations are often tiny. For downstream detailed routing, we use
    // a small tolerance window where we prefer lower "tightness" (peak
    // utilization) and fewer vias as a proxy for DR-friendliness.
    // Keep the window tight: we only allow tie-breaks (tightness/vias) when
    // the global wirelength difference is extremely small.
    constexpr double kWirelengthEps = 0.00008;  // 0.008%
    if (lhs.wirelength_dbu < rhs.wirelength_dbu * (1.0 - kWirelengthEps)) {
      return true;
    }
    if (rhs.wirelength_dbu < lhs.wirelength_dbu * (1.0 - kWirelengthEps)) {
      return false;
    }

    if (lhs.congestion.score != rhs.congestion.score) {
      return lhs.congestion.score < rhs.congestion.score;
    }

    if (lhs.via_count != rhs.via_count) {
      return lhs.via_count < rhs.via_count;
    }
    return lhs.wirelength_dbu < rhs.wirelength_dbu;
  };

  // NEWGR strategy (wirelength-first, via-second):
  // - Explore a small, deterministic set of perturbation/seed/critical-net
  //   configurations and pick the best by (1) routability, (2) global
  //   wirelength, and (3) via count.
  //
  // Rationale: for our regression, small deterministic perturbations can
  // improve downstream detailed-routing metrics by nudging congestion away
  // from hard-to-route pin-access regions, even if global cost differences
  // are small.
  const CandidateConfig baseline{"baseline",
                                 snapshot.caps_percentage,
                                 snapshot.perturbation_amount,
                                 snapshot.seed,
                                 snapshot.critical_percentage,
                                 snapshot.congestion_iterations,
                                 0,
                                 0,
                                 1.0f,
                                 0};

	  // Small candidate set tuned to be:
	  // - deterministic across runs (fixed seeds),
	  // - reasonably cheap (single-digit candidates),
	  // - wirelength-first, with DR-friendliness tie-breaks.
  const std::vector<CandidateConfig> candidates = {
      baseline,
      {"perturb3-seed11-crit0", 3.0f, 1, 11, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed11-crit0", 6.0f, 1, 11, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed17-crit0", 6.0f, 1, 17, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed23-crit0", 6.0f, 1, 23, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed29-crit0", 6.0f, 1, 29, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed31-crit0", 6.0f, 1, 31, 0.0f, snapshot.congestion_iterations},
      {"perturb6-seed37-crit0", 6.0f, 1, 37, 0.0f, snapshot.congestion_iterations},
      // Keep critical nets enabled while still perturbing capacities. The goal
      // is to preserve "straight" routes for a small subset of high-slack nets
      // while nudging congestion away from pin-access hotspots.
      {"perturb6-seed29-crit10", 6.0f, 1, 29, 10.0f, snapshot.congestion_iterations},
      // More perturbation magnitude can sometimes help DR escape from local
      // pin-access deadlocks at the cost of minor global detours.
      {"perturb6-amt2-seed29", 6.0f, 2, 29, 0.0f, snapshot.congestion_iterations},
      // Slightly fewer overflow iterations can reduce late-stage detours.
      {"perturb6-seed29-it40", 6.0f, 1, 29, 0.0f, 40},
      {"perturb9-seed11-crit0", 9.0f, 1, 11, 0.0f, snapshot.congestion_iterations},
      // A small "soft-capacity" reservation around the hottest Rudy tiles to
      // improve detailed-routability; chosen only if it stays competitive in
      // wirelength within tolerance and reduces tightness/vias.
      {"perturb6-seed29-rudy", 6.0f, 1, 29, 0.0f, snapshot.congestion_iterations, 25, 1, 0.85f, 2},
  };

	  const auto run_candidate = [&](const CandidateConfig& candidate) {
	    restore_snapshot(snapshot);
	    grouter_->setAllowCongestion(snapshot.allow_congestion);
	    std::vector<Net*> candidate_nets = canonical_nets;
	    prepare_fastroute(candidate, candidate_nets);

	    NetRouteMap routes = grouter_->findRouting(
	        candidate_nets, min_routing_layer, max_routing_layer);
	    RouteMetrics metrics = compute_metrics(routes);
	    metrics.total_overflow = grouter_->fastroute_->totalOverflow();
      metrics.congestion = compute_congestion_score();
	    logger_->info(GNR,
	                  6005,
	                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow {}, cong {:.3f} (max util {:.3f})",
                  candidate.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.total_overflow,
                  metrics.congestion.score,
                  metrics.congestion.max_util);
    return std::make_pair(std::move(routes), metrics);
  };

  auto finalize = [&](int overflow) {
    restore_snapshot(snapshot);
    if (overflow > 0 && !snapshot.allow_congestion) {
      logger_->warn(GNR,
                    6008,
                    "NEWGR finished with overflow ({}); enabling allow_congestion to avoid abort.",
                    overflow);
      grouter_->setAllowCongestion(true);
    }
  };

  CandidateConfig best_candidate = baseline;
  RouteMetrics best_metrics;
  bool have_best = false;

  for (const auto& candidate : candidates) {
    auto candidate_result = run_candidate(candidate);
    const RouteMetrics& metrics = candidate_result.second;
    if (!have_best || better(metrics, best_metrics)) {
      best_candidate = candidate;
      best_metrics = metrics;
      have_best = true;
    }
  }

  // Re-run the winner so `grouter_->fastroute_` internal state matches the
  // returned routes (important for subsequent incremental calls in the flow).
  auto [winner_routes, winner_metrics] = run_candidate(best_candidate);

  logger_->info(GNR,
                6007,
                "NEWGR picked {} (overflow {}, wl {:.0f} um, vias {})",
                best_candidate.name,
                winner_metrics.total_overflow,
                winner_metrics.wirelength_um,
                winner_metrics.via_count);
  finalize(winner_metrics.total_overflow);
  return winner_routes;
}

}  // namespace grt
