#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_set>
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
  bool override_global_adjustment = false;
  float global_adjustment = 0.0f;
};

struct GuideInflationConfig
{
  // Expand each (non-via) guide by adding parallel segments one tile away.
  // This is intended to give DRT more local flexibility and reduce detours
  // that increase final wirelength.
  int radius_tiles = 1;
  int min_length_tiles = 3;
  // Edges with small remaining resources are a proxy for regions where DRT
  // tends to detour. Treat edges with availability <= threshold as "risky"
  // and only inflate those to avoid guide bloat.
  int risky_avail_threshold = 2;
  int max_layer_inflate = std::numeric_limits<int>::max();
};

static void inflate_guides(NetRouteMap& routes,
                           const Grid* grid,
                           FastRouteCore* fastroute,
                           int min_routing_layer,
                           int max_routing_layer,
                           const GuideInflationConfig& cfg)
{
  if (grid == nullptr || fastroute == nullptr) {
    return;
  }
  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int x_grids = std::max(grid->getXGrids(), 1);
  const int y_grids = std::max(grid->getYGrids(), 1);

  auto to_grid_index = [&](int x_dbu, int y_dbu) -> std::pair<int, int> {
    const int x_idx = (x_dbu - x_min) / tile_size;
    const int y_idx = (y_dbu - y_min) / tile_size;
    return {std::clamp(x_idx, 0, x_grids - 1),
            std::clamp(y_idx, 0, y_grids - 1)};
  };

  auto to_dbu = [&](int x_idx, int y_idx) -> std::pair<int, int> {
    // GSegment coordinates are expected to be at GCELL centers (see
    // GlobalRouter::globalRoutingToBox()).
    const int x_dbu = x_min + x_idx * tile_size + (tile_size / 2);
    const int y_dbu = y_min + y_idx * tile_size + (tile_size / 2);
    return {x_dbu, y_dbu};
  };

  for (auto& [db_net, route] : routes) {
    static_cast<void>(db_net);

    std::unordered_set<GSegment, GSegmentHash> uniq;
    uniq.reserve(route.size() * 3 + 8);
    std::vector<GSegment> deduped;
    deduped.reserve(route.size());
    for (const auto& seg : route) {
      if (uniq.insert(seg).second) {
        deduped.push_back(seg);
      }
    }

    std::vector<GSegment> extra;
    extra.reserve(deduped.size() * 2);

    for (const auto& seg : deduped) {
      if (seg.isVia()) {
        continue;
      }
      if (seg.init_layer != seg.final_layer) {
        continue;
      }
      const int layer = seg.init_layer;
      if (layer < min_routing_layer || layer > max_routing_layer) {
        continue;
      }
      if (layer > cfg.max_layer_inflate) {
        continue;
      }

      const auto [gx0, gy0] = to_grid_index(seg.init_x, seg.init_y);
      const auto [gx1, gy1] = to_grid_index(seg.final_x, seg.final_y);
      const int len_tiles = std::abs(gx1 - gx0) + std::abs(gy1 - gy0);
      if (len_tiles < cfg.min_length_tiles) {
        continue;
      }

      const bool is_horizontal = (gy0 == gy1) && (gx0 != gx1);
      const bool is_vertical = (gx0 == gx1) && (gy0 != gy1);
      if (!is_horizontal && !is_vertical) {
        continue;
      }

      const int xs = std::min(gx0, gx1);
      const int xe = std::max(gx0, gx1);
      const int ys = std::min(gy0, gy1);
      const int ye = std::max(gy0, gy1);

      // Only widen guides when the current segment runs along an edge with low
      // remaining resources (proxy for DRT detours).
      bool is_risky = false;
      if (is_horizontal) {
        for (int x = xs; x < xe; x++) {
          const int avail = fastroute->getAvailableResources(x, gy0, x + 1, gy0, layer);
          if (avail <= cfg.risky_avail_threshold) {
            is_risky = true;
            break;
          }
        }
      } else {  // vertical
        for (int y = ys; y < ye; y++) {
          const int avail = fastroute->getAvailableResources(gx0, y, gx0, y + 1, layer);
          if (avail <= cfg.risky_avail_threshold) {
            is_risky = true;
            break;
          }
        }
      }
      if (!is_risky) {
        continue;
      }

      // Pick the "better" side (more remaining resources) and add a single
      // connected parallel guide there to avoid guide bloat.
      struct SideScore
      {
        bool valid = false;
        int min_avail = std::numeric_limits<int>::min();
        long sum_avail = std::numeric_limits<long>::min();
        int offset_idx = 0;
      };

      auto score_horizontal = [&](int gy_off) -> SideScore {
        SideScore s;
        if (gy_off < 0 || gy_off >= y_grids) {
          return s;
        }
        s.valid = true;
        s.offset_idx = gy_off;
        s.min_avail = std::numeric_limits<int>::max();
        s.sum_avail = 0;
        for (int x = xs; x < xe; x++) {
          const int avail = fastroute->getAvailableResources(x, gy_off, x + 1, gy_off, layer);
          s.min_avail = std::min(s.min_avail, avail);
          s.sum_avail += avail;
        }
        return s;
      };

      auto score_vertical = [&](int gx_off) -> SideScore {
        SideScore s;
        if (gx_off < 0 || gx_off >= x_grids) {
          return s;
        }
        s.valid = true;
        s.offset_idx = gx_off;
        s.min_avail = std::numeric_limits<int>::max();
        s.sum_avail = 0;
        for (int y = ys; y < ye; y++) {
          const int avail = fastroute->getAvailableResources(gx_off, y, gx_off, y + 1, layer);
          s.min_avail = std::min(s.min_avail, avail);
          s.sum_avail += avail;
        }
        return s;
      };

      SideScore best;
      if (is_horizontal) {
        const SideScore down = score_horizontal(gy0 - 1);
        const SideScore up = score_horizontal(gy0 + 1);
        best = down;
        if (!best.valid || (up.valid && (up.min_avail > best.min_avail
                                         || (up.min_avail == best.min_avail
                                             && up.sum_avail > best.sum_avail)))) {
          best = up;
        }
        if (!best.valid) {
          continue;
        }
        const int gy_off = best.offset_idx;
        const auto [x0_off_dbu, y_off_dbu] = to_dbu(gx0, gy_off);
        const auto [x1_off_dbu, y_off_dbu2] = to_dbu(gx1, gy_off);
        static_cast<void>(y_off_dbu2);
        const auto [x0_dbu, y0_dbu] = to_dbu(gx0, gy0);
        const auto [x1_dbu, y1_dbu] = to_dbu(gx1, gy0);
        static_cast<void>(y1_dbu);

        GSegment inflated(
            x0_off_dbu, y_off_dbu, layer, x1_off_dbu, y_off_dbu, layer);
        if (uniq.insert(inflated).second) {
          extra.push_back(inflated);
        }
        GSegment conn0(x0_dbu, y0_dbu, layer, x0_dbu, y_off_dbu, layer);
        if (uniq.insert(conn0).second) {
          extra.push_back(conn0);
        }
        GSegment conn1(x1_dbu, y0_dbu, layer, x1_dbu, y_off_dbu, layer);
        if (uniq.insert(conn1).second) {
          extra.push_back(conn1);
        }
      } else {  // vertical
        const SideScore left = score_vertical(gx0 - 1);
        const SideScore right = score_vertical(gx0 + 1);
        best = left;
        if (!best.valid || (right.valid && (right.min_avail > best.min_avail
                                            || (right.min_avail == best.min_avail
                                                && right.sum_avail
                                                       > best.sum_avail)))) {
          best = right;
        }
        if (!best.valid) {
          continue;
        }
        const int gx_off = best.offset_idx;
        const auto [x_off_dbu, y0_off_dbu] = to_dbu(gx_off, gy0);
        const auto [x_off_dbu2, y1_off_dbu] = to_dbu(gx_off, gy1);
        static_cast<void>(x_off_dbu2);
        const auto [x0_dbu, y0_dbu] = to_dbu(gx0, gy0);
        const auto [x1_dbu, y1_dbu] = to_dbu(gx0, gy1);
        static_cast<void>(x1_dbu);

        GSegment inflated(
            x_off_dbu, y0_off_dbu, layer, x_off_dbu, y1_off_dbu, layer);
        if (uniq.insert(inflated).second) {
          extra.push_back(inflated);
        }
        GSegment conn0(x0_dbu, y0_dbu, layer, x_off_dbu, y0_dbu, layer);
        if (uniq.insert(conn0).second) {
          extra.push_back(conn0);
        }
        GSegment conn1(x0_dbu, y1_dbu, layer, x_off_dbu, y1_dbu, layer);
        if (uniq.insert(conn1).second) {
          extra.push_back(conn1);
        }
      }
    }

    route = std::move(deduped);
    route.insert(route.end(), extra.begin(), extra.end());
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

  struct RouteMetrics
  {
    long wirelength_dbu = 0;
    long via_count = 0;
    double wirelength_um = 0.0;
    int total_overflow = 0;
    long saturated_edges = 0;
    long near_saturated_edges = 0;
    long risk_cost_dbu = 0;
    long via_cost_dbu = 0;
    long overflow_cost_dbu = 0;
    long score_dbu = 0;
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

  auto compute_risk_cost = [&](RouteMetrics& metrics, const NetRouteMap& routes) {
    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = std::max(grouter_->grid_->getXGrids(), 1);
    const int y_grids = std::max(grouter_->grid_->getYGrids(), 1);

    auto to_grid_index = [&](int x_dbu, int y_dbu) -> std::pair<int, int> {
      const int x_idx = (x_dbu - x_min) / tile_size;
      const int y_idx = (y_dbu - y_min) / tile_size;
      return {std::clamp(x_idx, 0, x_grids - 1),
              std::clamp(y_idx, 0, y_grids - 1)};
    };

    long saturated = 0;
    long near_saturated = 0;

    for (const auto& [db_net, segments] : routes) {
      static_cast<void>(db_net);
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          continue;
        }
        const int layer = segment.init_layer;
        const auto [x0, y0] = to_grid_index(segment.init_x, segment.init_y);
        const auto [x1, y1] = to_grid_index(segment.final_x, segment.final_y);

        if (y0 == y1 && x0 != x1) {
          const int y = y0;
          const int xs = std::min(x0, x1);
          const int xe = std::max(x0, x1);
          for (int x = xs; x < xe; x++) {
            const int avail
                = grouter_->fastroute_->getAvailableResources(x, y, x + 1, y, layer);
            if (avail <= 0) {
              saturated++;
            } else if (avail <= 2) {
              near_saturated++;
            }
          }
        } else if (x0 == x1 && y0 != y1) {
          const int x = x0;
          const int ys = std::min(y0, y1);
          const int ye = std::max(y0, y1);
          for (int y = ys; y < ye; y++) {
            const int avail
                = grouter_->fastroute_->getAvailableResources(x, y, x, y + 1, layer);
            if (avail <= 0) {
              saturated++;
            } else if (avail <= 2) {
              near_saturated++;
            }
          }
        }
      }
    }

    metrics.saturated_edges = saturated;
    metrics.near_saturated_edges = near_saturated;

    // Heuristic: near-saturated edges are a proxy for detailed-router detours.
    // Keep wirelength primary, but steer away from routes that hug capacity.
    // NEWGR adds a post-pass that widens guides; that extra flexibility tends
    // to make DRT less sensitive to "near-saturated" edges. So keep risk as a
    // tie-breaker rather than a dominant term.
    // Risk weights are intentionally low to keep wirelength dominant while
    // still biasing away from heavily saturated edges (which often correlate
    // with DRT detours).
    // Empirically, DRT wirelength tends to correlate more with "near-sat"
    // edges than with small differences in 2D global wirelength. Make risk
    // large enough to influence candidate choice, but still keep wirelength
    // as the dominant term.
    metrics.risk_cost_dbu = (saturated * tile_size) * 4
                            + (near_saturated * tile_size) / 8;
    metrics.via_cost_dbu = metrics.via_count * tile_size / 50;
    metrics.overflow_cost_dbu = static_cast<long>(metrics.total_overflow) * tile_size * 200;
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
  // - Use a known-good perturbation configuration (from prior iterations)
  //   to target wirelength improvements.
  // - Post-process the chosen global route to *widen* guides (add parallel
  //   guide segments) so DRT has more flexibility and is less likely to detour.
  logger_->info(GNR,
                6004,
                "NEWGR snapshot: caps% {:.1f} pert_amt {} seed {} crit% {:.1f} cong_iters {} global_adj {:.2f}",
                snapshot.caps_percentage,
                snapshot.perturbation_amount,
                snapshot.seed,
                snapshot.critical_percentage,
                snapshot.congestion_iterations,
                snapshot.global_adjustment);

  std::vector<CandidateConfig> candidates;
  candidates.push_back({"baseline",
                        snapshot.caps_percentage,
                        snapshot.perturbation_amount,
                        snapshot.seed,
                        snapshot.critical_percentage,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});

  // Candidate pool:
  // - Keep the "known-good" perturbation profile from earlier iterations.
  // - Explore a handful of seeds (bounded) and a couple critical-net settings.
  candidates.push_back({"perturb6-seed11-crit0",
                        6.0f,
                        1,
                        11,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb6-seed11-crit5",
                        6.0f,
                        1,
                        11,
                        5.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb6-seed17-crit0",
                        6.0f,
                        1,
                        17,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb6-seed23-crit0",
                        6.0f,
                        1,
                        23,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb6-seed29-crit0",
                        6.0f,
                        1,
                        29,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb4-seed11-crit0",
                        4.0f,
                        1,
                        11,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});
  candidates.push_back({"perturb8-seed11-crit0",
                        8.0f,
                        1,
                        11,
                        0.0f,
                        snapshot.congestion_iterations,
                        false,
                        snapshot.global_adjustment});

  // Occasionally, removing the flow-level global adjustment (capacity derate)
  // yields shorter routes that DRT can still realize with the extra guide
  // flexibility. Try it when the flow has a non-zero adjustment.
  if (snapshot.global_adjustment > 0.0f) {
    candidates.push_back({"perturb6-seed17-crit0-adj0",
                          6.0f,
                          1,
                          17,
                          0.0f,
                          snapshot.congestion_iterations,
                          true,
                          0.0f});
  }

  const auto better = [](const RouteMetrics& lhs, const RouteMetrics& rhs) {
    if (lhs.total_overflow != rhs.total_overflow) {
      return lhs.total_overflow < rhs.total_overflow;
    }
    if (lhs.score_dbu != rhs.score_dbu) {
      return lhs.score_dbu < rhs.score_dbu;
    }
    if (lhs.via_count != rhs.via_count) {
      return lhs.via_count < rhs.via_count;
    }
    return lhs.wirelength_dbu < rhs.wirelength_dbu;
  };

  NetRouteMap chosen_routes;
  RouteMetrics chosen_metrics;
  std::string chosen_name;
  CandidateConfig chosen_candidate;
  bool have_choice = false;

  for (const auto& candidate : candidates) {
    restore_snapshot(snapshot);
    if (candidate.override_global_adjustment) {
      grouter_->adjustment_ = candidate.global_adjustment;
    }
    prepare_fastroute(candidate);

    NetRouteMap routes
        = grouter_->findRouting(nets, min_routing_layer, max_routing_layer);
    RouteMetrics metrics = compute_metrics(routes);
    metrics.total_overflow = grouter_->fastroute_->totalOverflow();
    compute_risk_cost(metrics, routes);
    metrics.score_dbu = metrics.wirelength_dbu + metrics.via_cost_dbu
                        + metrics.risk_cost_dbu + metrics.overflow_cost_dbu;
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}, overflow {}, sat_edges {}, near_sat {}, score {}",
                  candidate.name,
                  metrics.wirelength_um,
                  metrics.via_count,
                  metrics.total_overflow,
                  metrics.saturated_edges,
                  metrics.near_saturated_edges,
                  metrics.score_dbu);

    if (!have_choice || better(metrics, chosen_metrics)) {
      chosen_routes = std::move(routes);
      chosen_metrics = metrics;
      chosen_name = candidate.name;
      chosen_candidate = candidate;
      have_choice = true;
    }
  }

  logger_->info(GNR, 6007, "NEWGR picked {}", chosen_name);

  // Re-run the chosen candidate so `fastroute_` resources match the returned
  // routes (the multi-candidate loop ends with FastRoute state from the last
  // explored candidate).
  restore_snapshot(snapshot);
  if (chosen_candidate.override_global_adjustment) {
    grouter_->adjustment_ = chosen_candidate.global_adjustment;
  }
  prepare_fastroute(chosen_candidate);
  NetRouteMap final_routes
      = grouter_->findRouting(nets, min_routing_layer, max_routing_layer);

  // Post-process: selectively widen guides around risky segments to give DRT
  // flexibility and reduce detours (primary objective: final wirelength).
  GuideInflationConfig inflate_cfg;
  inflate_cfg.radius_tiles = 1;
  inflate_cfg.min_length_tiles = 3;
  inflate_cfg.risky_avail_threshold = 2;
  inflate_cfg.max_layer_inflate = max_routing_layer;
  inflate_guides(final_routes,
                 grouter_->grid_,
                 grouter_->fastroute_,
                 min_routing_layer,
                 max_routing_layer,
                 inflate_cfg);

  // Leave GlobalRouter knobs as they were when entering NEWGR, but do not
  // rebuild the routing structures again; the returned guides are derived from
  // `final_routes` and are independent of these settings.
  restore_snapshot(snapshot);
  return final_routes;
}

}  // namespace grt
