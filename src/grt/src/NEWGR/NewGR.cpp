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

struct SelectionPolicy
{
  int via_tradeoff{1};
  int bend_tradeoff{1};
  int64_t via_guard{0};
  int64_t hard_via_guard{0};
  int64_t min_wl_improve{1};
  bool medium_net{false};
  bool long_net{false};
};

enum class RouteSource
{
  kFastRoute,
  kNewgrBalanced,
  kNewgrWirelength,
  kNewgrRegion,
  kNewgrRegularRegion,
  kNewgrFineGrain
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

SelectionPolicy buildSelectionPolicy(const RouteScore& baseline_score, int tile_size)
{
  SelectionPolicy policy;
  policy.via_tradeoff = 1;
  policy.bend_tradeoff = std::max(1, tile_size / 72);
  policy.via_guard = viaGuardForNet(baseline_score, tile_size) + 2;
  policy.hard_via_guard = policy.via_guard * 3 + 2;
  policy.min_wl_improve = 1;

  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 8;
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 16;
  if (baseline_score.wirelength >= long_net_threshold) {
    policy.long_net = true;
    policy.via_tradeoff = 1;
    policy.bend_tradeoff = std::max(1, tile_size / 180);
    policy.via_guard = policy.via_guard * 4 + 12;
    policy.hard_via_guard = policy.via_guard * 2 + 16;
    policy.min_wl_improve = 1;
  } else if (baseline_score.wirelength >= medium_net_threshold) {
    policy.medium_net = true;
    policy.via_tradeoff = 1;
    policy.bend_tradeoff = std::max(1, tile_size / 120);
    policy.via_guard += 8;
    policy.hard_via_guard = policy.via_guard * 3 + 12;
    policy.min_wl_improve = 1;
  }

  return policy;
}

int64_t effectiveWirelengthCost(const RouteScore& baseline_score,
                                const RouteScore& candidate_score,
                                const SelectionPolicy& policy)
{
  const int64_t extra_vias
      = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
  const int64_t penalized_vias
      = std::max<int64_t>(0, extra_vias - policy.via_guard);
  const int64_t extra_bends
      = std::max<int64_t>(0, candidate_score.bends - baseline_score.bends);
  return candidate_score.wirelength + penalized_vias * policy.via_tradeoff
         + extra_bends * policy.bend_tradeoff;
}

bool betterCandidate(const RouteScore& baseline_score,
                     const RouteScore& current_best_score,
                     const RouteScore& candidate_score,
                     const SelectionPolicy& policy)
{
  const int64_t current_extra_vias
      = std::max<int64_t>(0, current_best_score.vias - baseline_score.vias);
  const int64_t candidate_extra_vias
      = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);

  if (candidate_extra_vias > policy.hard_via_guard
      && candidate_extra_vias > current_extra_vias) {
    return false;
  }

  if (candidate_score.wirelength < current_best_score.wirelength) {
    const int64_t allowed_extra_vs_current
        = policy.long_net ? (policy.via_guard + 10)
                          : (policy.medium_net ? (policy.via_guard / 2 + 5) : 4);
    return candidate_extra_vias <= current_extra_vias + allowed_extra_vs_current;
  }
  if (candidate_score.wirelength > current_best_score.wirelength) {
    return false;
  }

  const int64_t current_effective_cost
      = effectiveWirelengthCost(baseline_score, current_best_score, policy);
  const int64_t candidate_effective_cost
      = effectiveWirelengthCost(baseline_score, candidate_score, policy);
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
  NetRouteMap region_routes = engine_->runRegionAware();
  NetRouteMap regular_region_routes = engine_->runRegularRegionAware();
  NetRouteMap finegrain_routes = engine_->runFineGrainRefine();

  int selected_from_balanced = 0;
  int selected_from_wl = 0;
  int selected_from_region = 0;
  int selected_from_regular_region = 0;
  int selected_from_finegrain = 0;
  int kept_fastroute = 0;
  int inserted_from_balanced = 0;
  int inserted_from_wl = 0;
  int inserted_from_region = 0;
  int inserted_from_regular_region = 0;
  int inserted_from_finegrain = 0;

  for (auto& [db_net, route] : routes) {
    const int tile_size = std::max(1, grouter_->getTileSize());
    const RouteScore baseline_score = scoreRoute(route);
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
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
              policy)) {
        best_score = candidate_score;
        selected_route = &candidate_it->second;
        selected_source = source;
      }
    };

    consider(balanced_routes, RouteSource::kNewgrBalanced);
    consider(wirelength_routes, RouteSource::kNewgrWirelength);
    consider(region_routes, RouteSource::kNewgrRegion);
    consider(regular_region_routes, RouteSource::kNewgrRegularRegion);
    consider(finegrain_routes, RouteSource::kNewgrFineGrain);

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
      case RouteSource::kNewgrRegion:
        selected_from_region++;
        break;
      case RouteSource::kNewgrRegularRegion:
        selected_from_regular_region++;
        break;
      case RouteSource::kNewgrFineGrain:
        selected_from_finegrain++;
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
  for (const auto& [db_net, route] : region_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_region++;
    }
  }
  for (const auto& [db_net, route] : regular_region_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_regular_region++;
    }
  }
  for (const auto& [db_net, route] : finegrain_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_finegrain++;
    }
  }

  logger_->info(utl::GRT,
                6004,
                "NEWGR hybrid selected {} balanced, {} WL-first, {} region, {} "
                "regular-region, and {} fine-grain routes (kept {} FastRoute, "
                "+{} balanced-only, +{} WL-only, +{} region-only, +{} "
                "regular-region-only, +{} fine-grain-only) out of {} total.",
                selected_from_balanced,
                selected_from_wl,
                selected_from_region,
                selected_from_regular_region,
                selected_from_finegrain,
                kept_fastroute,
                inserted_from_balanced,
                inserted_from_wl,
                inserted_from_region,
                inserted_from_regular_region,
                inserted_from_finegrain,
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
