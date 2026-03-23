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
  int64_t bends{0};
};

enum class RouteSource
{
  kFastRoute,
  kNewgrBalanced,
  kNewgrWirelength
};

enum class SegmentDirection
{
  kNone,
  kHorizontal,
  kVertical,
  kVia
};

SegmentDirection getSegmentDirection(const GSegment& segment)
{
  if (segment.init_layer != segment.final_layer) {
    return SegmentDirection::kVia;
  }
  if (segment.init_x != segment.final_x) {
    return SegmentDirection::kHorizontal;
  }
  if (segment.init_y != segment.final_y) {
    return SegmentDirection::kVertical;
  }
  return SegmentDirection::kNone;
}

RouteScore scoreRoute(const GRoute& route)
{
  RouteScore score;
  SegmentDirection previous_planar_direction = SegmentDirection::kNone;
  for (const GSegment& segment : route) {
    score.wirelength += segment.length();
    score.vias += std::abs(segment.init_layer - segment.final_layer);
    const SegmentDirection direction = getSegmentDirection(segment);
    if (direction == SegmentDirection::kHorizontal
        || direction == SegmentDirection::kVertical) {
      if (previous_planar_direction != SegmentDirection::kNone
          && previous_planar_direction != direction) {
        score.bends++;
      }
      previous_planar_direction = direction;
    }
  }
  return score;
}

int64_t viaGuardForNet(const RouteScore& baseline_score, int tile_size)
{
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 24;
  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 12;
  if (baseline_score.wirelength >= long_net_threshold) {
    return std::max<int64_t>(10, baseline_score.vias / 2 + 4);
  }
  if (baseline_score.wirelength >= medium_net_threshold) {
    return std::max<int64_t>(4, baseline_score.vias / 3 + 2);
  }
  return std::max<int64_t>(1, baseline_score.vias / 6);
}

int64_t effectiveWirelengthCost(const RouteScore& baseline_score,
                                const RouteScore& candidate_score,
                                int via_tradeoff,
                                int bend_tradeoff,
                                int64_t via_guard)
{
  const int64_t extra_vias
      = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
  const int64_t penalized_vias = std::max<int64_t>(0, extra_vias - via_guard);
  const int64_t extra_bends
      = std::max<int64_t>(0, candidate_score.bends - baseline_score.bends);
  return candidate_score.wirelength + penalized_vias * via_tradeoff
         + extra_bends * bend_tradeoff;
}

bool betterCandidate(const RouteScore& baseline_score,
                     const RouteScore& current_best_score,
                     const RouteScore& candidate_score,
                     int via_tradeoff,
                     int bend_tradeoff,
                     int64_t via_guard)
{
  const int64_t current_effective_cost
      = effectiveWirelengthCost(
          baseline_score, current_best_score, via_tradeoff, bend_tradeoff, via_guard);
  const int64_t candidate_effective_cost
      = effectiveWirelengthCost(
          baseline_score, candidate_score, via_tradeoff, bend_tradeoff, via_guard);
  if (candidate_effective_cost != current_effective_cost) {
    return candidate_effective_cost < current_effective_cost;
  }
  if (candidate_score.wirelength != current_best_score.wirelength) {
    return candidate_score.wirelength < current_best_score.wirelength;
  }
  if (candidate_score.vias != current_best_score.vias) {
    return candidate_score.vias < current_best_score.vias;
  }
  return candidate_score.bends < current_best_score.bends;
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
  const int via_tradeoff = std::max(1, tile_size / 30);
  const int bend_tradeoff = std::max(1, tile_size / 48);

  int selected_from_balanced = 0;
  int selected_from_wl = 0;
  int kept_fastroute = 0;
  int inserted_from_balanced = 0;
  int inserted_from_wl = 0;

  for (auto& [db_net, route] : routes) {
    const RouteScore baseline_score = scoreRoute(route);
    const int64_t via_guard = viaGuardForNet(baseline_score, tile_size);
    RouteScore best_score = baseline_score;
    const GRoute* selected_route = &route;
    RouteSource selected_source = RouteSource::kFastRoute;

    auto consider = [&](const NetRouteMap& candidate_routes, RouteSource source) {
      auto candidate_it = candidate_routes.find(db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      if (betterCandidate(
              baseline_score,
              best_score,
              candidate_score,
              via_tradeoff,
              bend_tradeoff,
              via_guard)) {
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
                "NEWGR hybrid selected {} balanced and {} WL-first routes (kept "
                "{} FastRoute, +{} balanced-only, +{} WL-only) out of {} total.",
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
