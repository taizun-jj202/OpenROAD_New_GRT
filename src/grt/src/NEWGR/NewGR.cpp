#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "NEWGR/src/NewgrEngine.h"
#include "Net.h"
#include "Pin.h"
#include "fastroute/include/FastRoute.h"
#include "utl/Logger.h"

namespace grt {

namespace {

struct RouteScore
{
  uint64_t wirelength{0};
  uint64_t vias{0};
  uint64_t routed_nets{0};
  uint64_t segments{0};
};

RouteScore computeRouteScore(const GRoute& route)
{
  RouteScore score;
  for (const GSegment& segment : route) {
    score.wirelength += static_cast<uint64_t>(segment.length());
    ++score.segments;
    if (segment.isVia()) {
      score.vias += static_cast<uint64_t>(
          std::abs(segment.final_layer - segment.init_layer));
    }
  }
  return score;
}

RouteScore computeRouteScore(const NetRouteMap& routes)
{
  RouteScore score;
  for (const auto& [db_net, route] : routes) {
    (void) db_net;
    const RouteScore net_score = computeRouteScore(route);
    if (net_score.segments == 0) {
      continue;
    }
    score.wirelength += net_score.wirelength;
    score.vias += net_score.vias;
    score.segments += net_score.segments;
    ++score.routed_nets;
  }
  return score;
}

uint64_t routeCostForPins(const RouteScore& score, int pin_count)
{
  // Prioritize wirelength while keeping via count under control.
  const uint64_t via_weight
      = (pin_count <= 4) ? 10 : ((pin_count <= 12) ? 14 : 18);
  return score.wirelength + score.vias * via_weight;
}

uint64_t minWirelengthGainForPins(int pin_count)
{
  if (pin_count <= 4) {
    return 32;
  }
  if (pin_count <= 12) {
    return 120;
  }
  if (pin_count <= 24) {
    return 320;
  }
  return 800;
}

int64_t viaBudgetForPins(int pin_count)
{
  if (pin_count <= 4) {
    return 2;
  }
  if (pin_count <= 12) {
    return 4;
  }
  if (pin_count <= 24) {
    return 6;
  }
  return 8;
}

bool shouldUseNewgrOverFastRoute(const RouteScore& newgr_score,
                                 const RouteScore& fastroute_score,
                                 int pin_count)
{
  if (newgr_score.segments == 0) {
    return false;
  }
  if (fastroute_score.segments == 0) {
    return true;
  }

  // Keep high-fanout nets on FastRoute for topology stability.
  if (pin_count > 24) {
    return false;
  }

  const int64_t wl_gain = static_cast<int64_t>(fastroute_score.wirelength)
                          - static_cast<int64_t>(newgr_score.wirelength);
  const int64_t via_delta = static_cast<int64_t>(newgr_score.vias)
                            - static_cast<int64_t>(fastroute_score.vias);
  const int64_t via_budget = viaBudgetForPins(pin_count);
  const int64_t min_wl_gain
      = static_cast<int64_t>(minWirelengthGainForPins(pin_count));

  // Primary gate: meaningful wirelength reduction with controlled via impact.
  if (wl_gain >= min_wl_gain && via_delta <= via_budget) {
    return true;
  }

  // Allow near-neutral wirelength moves when NEWGR removes many vias.
  if (via_delta <= -3 && wl_gain >= -(min_wl_gain / 4)) {
    return true;
  }

  const uint64_t newgr_cost = routeCostForPins(newgr_score, pin_count);
  const uint64_t fastroute_cost = routeCostForPins(fastroute_score, pin_count);
  if (newgr_cost + minWirelengthGainForPins(pin_count) < fastroute_cost
      && via_delta <= via_budget * 2) {
    return true;
  }

  // Strong wirelength win override.
  if (wl_gain >= std::max<int64_t>(
                     min_wl_gain * 4,
                     static_cast<int64_t>(fastroute_score.wirelength / 10))
      && via_delta <= via_budget * 2) {
    return true;
  }

  return false;
}

struct HybridStats
{
  uint64_t replaced_with_newgr{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
};

NetRouteMap buildHybridRoutes(const NetRouteMap& newgr_routes,
                              const NetRouteMap& fastroute_routes,
                              const std::map<odb::dbNet*, Net*>& db_net_map,
                              HybridStats& stats)
{
  // FastRoute-first backbone with selective NEWGR replacements.
  NetRouteMap hybrid_routes = fastroute_routes;
  for (const auto& [db_net, newgr_route] : newgr_routes) {
    if (newgr_route.empty()) {
      continue;
    }
    const auto net_it = db_net_map.find(db_net);
    if (net_it == db_net_map.end() || net_it->second == nullptr) {
      continue;
    }

    const int pin_count = std::max(1, net_it->second->getNumPins());
    auto hybrid_it = hybrid_routes.find(db_net);
    if (hybrid_it == hybrid_routes.end() || hybrid_it->second.empty()) {
      hybrid_routes[db_net] = newgr_route;
      ++stats.added_missing_nets;
      continue;
    }

    const RouteScore fastroute_score = computeRouteScore(hybrid_it->second);
    const RouteScore newgr_score = computeRouteScore(newgr_route);
    ++stats.considered_nets;
    if (shouldUseNewgrOverFastRoute(newgr_score, fastroute_score, pin_count)) {
      hybrid_it->second = newgr_route;
      ++stats.replaced_with_newgr;
    }
  }

  return hybrid_routes;
}

bool isBetterRoute(int overflow_a,
                   const RouteScore& score_a,
                   int overflow_b,
                   const RouteScore& score_b)
{
  if (overflow_a != overflow_b) {
    return overflow_a < overflow_b;
  }
  if (score_a.routed_nets != score_b.routed_nets) {
    return score_a.routed_nets > score_b.routed_nets;
  }
  if (score_a.wirelength != score_b.wirelength) {
    return score_a.wirelength < score_b.wirelength;
  }
  if (score_a.vias != score_b.vias) {
    return score_a.vias < score_b.vias;
  }
  return false;
}

}  // namespace

NewGR::NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NewGR::~NewGR() = default;

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  if (nets.empty()) {
    return {};
  }

  if (!grouter_->hasSprouteGridData() || !grouter_->hasSprouteNetData()) {
    logger_->error(utl::GRT,
                   6003,
                   "NEWGR router selected, but grid/net data is not initialized.");
  }

  if (!engine_) {
    engine_ = std::make_unique<NewgrEngine>(logger_);
  }

  engine_->init(grouter_->sproute_grid_data_, grouter_->sproute_nets_);
  NetRouteMap routes = engine_->run();
  last_total_overflow_ = engine_->getTotalOverflow();
  used_fastroute_last_run_ = false;

  const bool disable_fastroute_graft
      = std::getenv("NEWGR_DISABLE_FASTROUTE_GRAFT") != nullptr;
  if (!disable_fastroute_graft && grouter_->fastroute() != nullptr) {
    NetRouteMap fastroute_routes = grouter_->fastroute()->run();
    const int fastroute_overflow = grouter_->fastroute()->totalOverflow();
    const RouteScore newgr_score = computeRouteScore(routes);
    const RouteScore fastroute_score = computeRouteScore(fastroute_routes);
    HybridStats hybrid_stats;
    NetRouteMap hybrid_routes
        = buildHybridRoutes(routes, fastroute_routes, grouter_->db_net_map_, hybrid_stats);
    const RouteScore hybrid_score = computeRouteScore(hybrid_routes);

    logger_->info(utl::GRT,
                  6006,
                  "NEWGR candidate summary: NEWGR(ofl={}, wl={}, vias={}, nets={}), "
                  "HYBRID_FR_BASE(wl={}, vias={}, nets={}, adopt_newgr={}, add={}), "
                  "FastRoute(ofl={}, wl={}, vias={}, nets={})",
                  last_total_overflow_,
                  newgr_score.wirelength,
                  newgr_score.vias,
                  newgr_score.routed_nets,
                  hybrid_score.wirelength,
                  hybrid_score.vias,
                  hybrid_score.routed_nets,
                  hybrid_stats.replaced_with_newgr,
                  hybrid_stats.added_missing_nets,
                  fastroute_overflow,
                  fastroute_score.wirelength,
                  fastroute_score.vias,
                  fastroute_score.routed_nets);

    RouteScore best_score = newgr_score;
    int best_overflow = last_total_overflow_;
    const char* selected_label = "NEWGR-core";
    bool selected_hybrid = false;

    if (isBetterRoute(last_total_overflow_,
                      hybrid_score,
                      best_overflow,
                      best_score)) {
      routes = std::move(hybrid_routes);
      best_score = hybrid_score;
      selected_label = "NEWGR+FastRoute net-graft";
      selected_hybrid = true;
    }

    if (isBetterRoute(fastroute_overflow,
                      fastroute_score,
                      best_overflow,
                      best_score)) {
      routes = std::move(fastroute_routes);
      best_overflow = fastroute_overflow;
      best_score = fastroute_score;
      selected_label = "FastRoute";
      used_fastroute_last_run_ = true;
    }

    last_total_overflow_ = best_overflow;
    if (used_fastroute_last_run_) {
      logger_->info(utl::GRT,
                    6004,
                    "NEWGR selected {} path: overflow={}, route_wl={}, "
                    "route_vias={}, routed_nets={}",
                    selected_label,
                    last_total_overflow_,
                    best_score.wirelength,
                    best_score.vias,
                    best_score.routed_nets);
    } else if (selected_hybrid) {
      logger_->info(utl::GRT,
                    6007,
                    "NEWGR selected {} path: overflow={}, route_wl={}, "
                    "route_vias={}, routed_nets={}, swapped_nets={}, added_nets={}",
                    selected_label,
                    last_total_overflow_,
                    best_score.wirelength,
                    best_score.vias,
                    best_score.routed_nets,
                    hybrid_stats.replaced_with_newgr,
                    hybrid_stats.added_missing_nets);
    } else {
      logger_->info(utl::GRT,
                    6005,
                    "NEWGR kept NEWGR-core path: overflow={}, route_wl={}, "
                    "route_vias={}, routed_nets={}",
                    last_total_overflow_,
                    best_score.wirelength,
                    best_score.vias,
                    best_score.routed_nets);
    }
  } else if (disable_fastroute_graft && grouter_->fastroute() != nullptr) {
    logger_->info(utl::GRT,
                  6008,
                  "NEWGR FastRoute net-graft is disabled via "
                  "NEWGR_DISABLE_FASTROUTE_GRAFT.");
  }

  grouter_->addRemainingGuides(routes, nets, min_routing_layer, max_routing_layer);
  grouter_->connectPadPins(routes);
  for (auto& net_route : routes) {
    std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
    GRoute& route = net_route.second;
    grouter_->mergeSegments(pins, route);
  }

  return routes;
}

int NewGR::getTotalOverflow() const
{
  return last_total_overflow_;
}

void NewGR::updateDbCongestion(odb::dbBlock* block)
{
  if (used_fastroute_last_run_ && grouter_->fastroute() != nullptr) {
    grouter_->fastroute()->updateDbCongestion(grouter_->getMinRoutingLayer(),
                                              grouter_->getMaxRoutingLayer());
    return;
  }
  if (engine_ == nullptr) {
    return;
  }
  engine_->updateDbCongestion(block);
}

}  // namespace grt
