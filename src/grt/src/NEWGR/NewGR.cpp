#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <memory>
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
  int64_t wirelength{0};
  int64_t vias{0};
  int64_t weighted_cost{0};
};

enum class RouteSource
{
  kFastRoute,
  kNewgrBalanced,
  kNewgrWirelength
};

RouteScore scoreRoute(const GRoute& route, int via_penalty)
{
  RouteScore score;
  for (const GSegment& segment : route) {
    score.wirelength += segment.length();
    score.vias += std::abs(segment.init_layer - segment.final_layer);
  }
  score.weighted_cost = score.wirelength + score.vias * via_penalty;
  return score;
}

int64_t effectiveWirelengthCost(const RouteScore& baseline_score,
                                const RouteScore& candidate_score,
                                int via_tradeoff)
{
  const int64_t extra_vias
      = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
  return candidate_score.wirelength + extra_vias * via_tradeoff;
}

bool betterCandidate(const RouteScore& baseline_score,
                     const RouteScore& current_best_score,
                     const RouteScore& candidate_score,
                     int via_tradeoff)
{
  const int64_t current_effective_cost
      = effectiveWirelengthCost(baseline_score, current_best_score, via_tradeoff);
  const int64_t candidate_effective_cost
      = effectiveWirelengthCost(baseline_score, candidate_score, via_tradeoff);
  if (candidate_effective_cost != current_effective_cost) {
    return candidate_effective_cost < current_effective_cost;
  }
  if (candidate_score.wirelength != current_best_score.wirelength) {
    return candidate_score.wirelength < current_best_score.wirelength;
  }
  if (candidate_score.vias != current_best_score.vias) {
    return candidate_score.vias < current_best_score.vias;
  }
  return candidate_score.weighted_cost < current_best_score.weighted_cost;
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

  NetRouteMap routes = grouter_->fastroute()->run();

  engine_->init(grouter_->sproute_grid_data_, grouter_->sproute_nets_);
  NetRouteMap balanced_routes = engine_->run();
  NetRouteMap wirelength_routes = engine_->runWirelengthFirst();

  const int tile_size = std::max(1, grouter_->getTileSize());
  const int via_penalty = std::max(1, tile_size / 16);
  const int via_tradeoff = std::max(1, tile_size / 20);

  int selected_from_balanced = 0;
  int selected_from_wl = 0;
  int kept_fastroute = 0;
  int inserted_from_balanced = 0;
  int inserted_from_wl = 0;

  for (auto& [db_net, route] : routes) {
    const RouteScore baseline_score = scoreRoute(route, via_penalty);
    RouteScore best_score = baseline_score;
    const GRoute* selected_route = &route;
    RouteSource selected_source = RouteSource::kFastRoute;

    auto consider = [&](const NetRouteMap& candidate_routes, RouteSource source) {
      auto candidate_it = candidate_routes.find(db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      const RouteScore candidate_score
          = scoreRoute(candidate_it->second, via_penalty);
      if (betterCandidate(
              baseline_score, best_score, candidate_score, via_tradeoff)) {
        best_score = candidate_score;
        selected_route = &candidate_it->second;
        selected_source = source;
      }
    };

    consider(balanced_routes, RouteSource::kNewgrBalanced);
    consider(wirelength_routes, RouteSource::kNewgrWirelength);

    if (selected_route != &route) {
      route = *selected_route;
    }

    switch (selected_source) {
      case RouteSource::kFastRoute:
        kept_fastroute++;
        break;
      case RouteSource::kNewgrBalanced:
        selected_from_balanced++;
        break;
      case RouteSource::kNewgrWirelength:
        selected_from_wl++;
        break;
    }
  }

  for (const auto& [db_net, route] : balanced_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_balanced++;
    }
  }
  for (const auto& [db_net, route] : wirelength_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_wl++;
    }
  }

  logger_->info(utl::GRT,
                6004,
                "NEWGR hybrid selected {} balanced routes, {} WL-first routes "
                "(kept {} FastRoute, +{} balanced-only, +{} WL-only) out of {} "
                "total.",
                selected_from_balanced,
                selected_from_wl,
                kept_fastroute,
                inserted_from_balanced,
                inserted_from_wl,
                routes.size());

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
  if (!engine_) {
    return 0;
  }
  return engine_->getTotalOverflow();
}

void NewGR::updateDbCongestion(odb::dbBlock* block)
{
  if (!engine_) {
    return;
  }
  engine_->updateDbCongestion(block);
}

}  // namespace grt
