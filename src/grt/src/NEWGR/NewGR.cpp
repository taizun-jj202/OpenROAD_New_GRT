#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include <utility>
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
  int64_t wl_per_extra_via{8};
  int64_t aggressive_wl_gain{1};
  bool medium_net{false};
  bool long_net{false};
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

enum class EdgeAxis : uint8_t
{
  kHorizontal = 0,
  kVertical = 1
};

struct RouteEdgeKey
{
  int x{0};
  int y{0};
  int layer{0};  // 0-based routing layer.
  EdgeAxis axis{EdgeAxis::kHorizontal};

  bool operator==(const RouteEdgeKey& rhs) const
  {
    return x == rhs.x && y == rhs.y && layer == rhs.layer && axis == rhs.axis;
  }
};

struct RouteEdgeKeyHash
{
  std::size_t operator()(const RouteEdgeKey& key) const
  {
    std::size_t hash = static_cast<std::size_t>(key.x * 1315423911u);
    hash ^= static_cast<std::size_t>(key.y + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    hash ^= static_cast<std::size_t>(key.layer + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    hash ^= static_cast<std::size_t>(static_cast<uint8_t>(key.axis)
                                     + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    return hash;
  }
};

using EdgeUsageMap = std::unordered_map<RouteEdgeKey, int, RouteEdgeKeyHash>;

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

int gridIndexFromCoord(int coord, int origin, int tile_size)
{
  if (tile_size <= 0) {
    return 0;
  }
  return (coord - origin) / tile_size;
}

template <typename Fn>
void forEachUnitPlanarEdge(const GRoute& route,
                           int origin_x,
                           int origin_y,
                           int tile_size,
                           Fn&& fn)
{
  if (tile_size <= 0) {
    return;
  }

  for (const GSegment& segment : route) {
    if (segment.init_layer != segment.final_layer) {
      continue;
    }

    const int layer = std::max(0, segment.init_layer - 1);
    const int x0 = gridIndexFromCoord(segment.init_x, origin_x, tile_size);
    const int y0 = gridIndexFromCoord(segment.init_y, origin_y, tile_size);
    const int x1 = gridIndexFromCoord(segment.final_x, origin_x, tile_size);
    const int y1 = gridIndexFromCoord(segment.final_y, origin_y, tile_size);

    if (x0 == x1 && y0 == y1) {
      continue;
    }

    if (y0 == y1) {
      const int step = (x1 > x0) ? 1 : -1;
      for (int x = x0; x != x1; x += step) {
        fn(RouteEdgeKey{
            std::min(x, x + step), y0, layer, EdgeAxis::kHorizontal});
      }
      continue;
    }

    if (x0 == x1) {
      const int step = (y1 > y0) ? 1 : -1;
      for (int y = y0; y != y1; y += step) {
        fn(RouteEdgeKey{
            x0, std::min(y, y + step), layer, EdgeAxis::kVertical});
      }
      continue;
    }

    // Defensive Manhattan decomposition for non-orthogonal merged segments.
    const int x_step = (x1 > x0) ? 1 : -1;
    for (int x = x0; x != x1; x += x_step) {
      fn(RouteEdgeKey{
          std::min(x, x + x_step), y0, layer, EdgeAxis::kHorizontal});
    }
    const int y_step = (y1 > y0) ? 1 : -1;
    for (int y = y0; y != y1; y += y_step) {
      fn(RouteEdgeKey{
          x1, std::min(y, y + y_step), layer, EdgeAxis::kVertical});
    }
  }
}

int edgeHardCapacity(const RouteEdgeKey& edge, const SprouteGridData& grid)
{
  if (edge.layer < 0) {
    return 1;
  }

  if (edge.axis == EdgeAxis::kHorizontal) {
    if (edge.layer < static_cast<int>(grid.h_capacities.size())) {
      return std::max<int>(1, grid.h_capacities[edge.layer]);
    }
    return 1;
  }

  if (edge.layer < static_cast<int>(grid.v_capacities.size())) {
    return std::max<int>(1, grid.v_capacities[edge.layer]);
  }
  return 1;
}

double softRatioFromCongestion(const RouteEdgeKey& edge, double congestion)
{
  // SPRoute-style logistic soft capacity reservation, tighter on lower layers.
  double min_ratio = 0.65;
  double max_ratio = 0.92;
  double cong_mid = 0.85;
  double slope = 3.0;
  if (edge.layer <= 1) {
    min_ratio = 0.55;
    max_ratio = 0.86;
    cong_mid = 0.72;
    slope = 3.8;
  } else if (edge.layer <= 3) {
    min_ratio = 0.60;
    max_ratio = 0.90;
    cong_mid = 0.80;
    slope = 3.4;
  }
  return min_ratio
         + (max_ratio - min_ratio) / (1.0 + std::exp((congestion - cong_mid) * slope));
}

int softCapacityFromDemand(const RouteEdgeKey& edge,
                           int hard_capacity,
                           int baseline_demand)
{
  if (hard_capacity <= 1) {
    return 1;
  }
  const double congestion
      = static_cast<double>(baseline_demand) / static_cast<double>(hard_capacity);
  const double ratio = softRatioFromCongestion(edge, congestion);
  const int soft_capacity = static_cast<int>(std::floor(hard_capacity * ratio));
  return std::max(1, std::min(hard_capacity, soft_capacity));
}

void addRouteToUsage(const GRoute& route,
                     int origin_x,
                     int origin_y,
                     int tile_size,
                     EdgeUsageMap& usage)
{
  forEachUnitPlanarEdge(
      route, origin_x, origin_y, tile_size, [&](const RouteEdgeKey& edge) {
        usage[edge]++;
      });
}

int64_t routeCongestionPenalty(const GRoute& route,
                               const SprouteGridData& grid,
                               int origin_x,
                               int origin_y,
                               int tile_size,
                               const EdgeUsageMap& selected_usage,
                               const EdgeUsageMap& soft_capacities)
{
  double penalty = 0.0;
  forEachUnitPlanarEdge(
      route, origin_x, origin_y, tile_size, [&](const RouteEdgeKey& edge) {
        const auto usage_it = selected_usage.find(edge);
        const int current_usage = usage_it == selected_usage.end() ? 0 : usage_it->second;

        int soft_capacity = edgeHardCapacity(edge, grid);
        const auto soft_it = soft_capacities.find(edge);
        if (soft_it != soft_capacities.end()) {
          soft_capacity = soft_it->second;
        } else {
          soft_capacity = softCapacityFromDemand(edge, soft_capacity, 0);
        }

        const int projected_usage = current_usage + 1;
        const int overflow = std::max(0, projected_usage - soft_capacity);
        const double logistic_cost
            = 1.0 / (1.0 + std::exp(1.6 * static_cast<double>(soft_capacity - projected_usage)));

        // CUGR-like expected overflow + logistic routability pressure.
        penalty += static_cast<double>(overflow) * 28.0 + logistic_cost * 4.0;
        if (projected_usage >= soft_capacity) {
          penalty += 1.0;
        }
      });
  return static_cast<int64_t>(std::llround(penalty));
}

int64_t viaGuardForNet(const RouteScore& baseline_score, int tile_size)
{
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 20;
  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 10;
  if (baseline_score.wirelength >= long_net_threshold) {
    return std::max<int64_t>(12, baseline_score.vias / 2 + 6);
  }
  if (baseline_score.wirelength >= medium_net_threshold) {
    return std::max<int64_t>(5, baseline_score.vias / 3 + 3);
  }
  return std::max<int64_t>(1, baseline_score.vias / 8);
}

SelectionPolicy buildSelectionPolicy(const RouteScore& baseline_score, int tile_size)
{
  SelectionPolicy policy;
  policy.via_tradeoff = 1;
  policy.bend_tradeoff = std::max(1, tile_size / 160);
  policy.via_guard = viaGuardForNet(baseline_score, tile_size) + 4;
  policy.hard_via_guard = policy.via_guard * 4 + 8;
  policy.min_wl_improve = 1;
  policy.wl_per_extra_via = 12;
  policy.aggressive_wl_gain = std::max<int64_t>(4, tile_size / 3);

  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 8;
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 16;
  if (baseline_score.wirelength >= long_net_threshold) {
    policy.long_net = true;
    policy.via_tradeoff = 1;
    policy.bend_tradeoff = std::max(1, tile_size / 220);
    policy.via_guard = policy.via_guard * 5 + 20;
    policy.hard_via_guard = policy.via_guard * 2 + 28;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 4;
    policy.aggressive_wl_gain = std::max<int64_t>(2, tile_size / 4);
  } else if (baseline_score.wirelength >= medium_net_threshold) {
    policy.medium_net = true;
    policy.via_tradeoff = 1;
    policy.bend_tradeoff = std::max(1, tile_size / 180);
    policy.via_guard = policy.via_guard * 2 + 12;
    policy.hard_via_guard = policy.via_guard * 3 + 18;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 8;
    policy.aggressive_wl_gain = std::max<int64_t>(3, tile_size / 3);
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
      = std::max<int64_t>(0, extra_vias - policy.via_guard / 2);
  const int64_t extra_bends
      = std::max<int64_t>(0, candidate_score.bends - baseline_score.bends);
  return candidate_score.wirelength + penalized_vias * policy.via_tradeoff
         + extra_bends * policy.bend_tradeoff;
}

int64_t extraViaBudget(const SelectionPolicy& policy,
                       int64_t current_extra_vias,
                       int64_t wl_gain_vs_current)
{
  const int64_t wl_credit
      = wl_gain_vs_current / std::max<int64_t>(1, policy.wl_per_extra_via);
  return std::min<int64_t>(
      policy.hard_via_guard,
      current_extra_vias + policy.via_guard + wl_credit);
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

  const int64_t wl_gain_vs_current = std::max<int64_t>(
      0, current_best_score.wirelength - candidate_score.wirelength);
  if (wl_gain_vs_current >= policy.min_wl_improve) {
    if (candidate_extra_vias
        <= extraViaBudget(policy, current_extra_vias, wl_gain_vs_current)) {
      return true;
    }
    const int64_t aggressive_via_budget
        = policy.hard_via_guard
          + wl_gain_vs_current
                / std::max<int64_t>(1, policy.wl_per_extra_via * 2);
    return wl_gain_vs_current >= policy.aggressive_wl_gain
           && candidate_extra_vias <= aggressive_via_budget;
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

  const SprouteGridData& grid = grouter_->sproute_grid_data_;
  const int origin_x = grid.origin.x();
  const int origin_y = grid.origin.y();
  const int tile_size = std::max(1, grouter_->getTileSize());

  EdgeUsageMap baseline_demand;
  for (const auto& [db_net, route] : routes) {
    (void) db_net;
    addRouteToUsage(route, origin_x, origin_y, tile_size, baseline_demand);
  }

  EdgeUsageMap soft_capacities;
  soft_capacities.reserve(baseline_demand.size());
  for (const auto& [edge, demand] : baseline_demand) {
    const int hard_capacity = edgeHardCapacity(edge, grid);
    soft_capacities.emplace(
        edge, softCapacityFromDemand(edge, hard_capacity, demand));
  }

  struct OrderedNet
  {
    odb::dbNet* db_net{nullptr};
    RouteScore baseline_score;
  };

  std::vector<OrderedNet> ordered_nets;
  ordered_nets.reserve(routes.size());
  for (const auto& [db_net, route] : routes) {
    ordered_nets.push_back({db_net, scoreRoute(route)});
  }
  std::sort(ordered_nets.begin(),
            ordered_nets.end(),
            [](const OrderedNet& lhs, const OrderedNet& rhs) {
              if (lhs.baseline_score.wirelength != rhs.baseline_score.wirelength) {
                return lhs.baseline_score.wirelength > rhs.baseline_score.wirelength;
              }
              if (lhs.baseline_score.vias != rhs.baseline_score.vias) {
                return lhs.baseline_score.vias > rhs.baseline_score.vias;
              }
              return lhs.db_net < rhs.db_net;
            });

  int selected_from_balanced = 0;
  int selected_from_wl = 0;
  int kept_fastroute = 0;
  int inserted_from_balanced = 0;
  int inserted_from_wl = 0;

  EdgeUsageMap selected_usage;
  selected_usage.reserve(baseline_demand.size());

  for (const OrderedNet& ordered_net : ordered_nets) {
    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      continue;
    }

    GRoute& route = route_it->second;
    const RouteScore baseline_score = ordered_net.baseline_score;
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
    const int64_t congestion_tradeoff = policy.long_net ? 1 : (policy.medium_net ? 2 : 3);

    RouteScore best_score = baseline_score;
    int64_t best_congestion_cost
        = routeCongestionPenalty(route,
                                 grid,
                                 origin_x,
                                 origin_y,
                                 tile_size,
                                 selected_usage,
                                 soft_capacities);
    int64_t best_total_cost = effectiveWirelengthCost(
                                  baseline_score, best_score, policy)
                              + congestion_tradeoff * best_congestion_cost;
    const GRoute* selected_route = &route;
    RouteSource selected_source = RouteSource::kFastRoute;

    auto consider = [&](const NetRouteMap& candidate_routes, RouteSource source) {
      auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      if (!betterCandidate(
              baseline_score, best_score, candidate_score, policy)) {
        return;
      }

      const int64_t candidate_congestion_cost
          = routeCongestionPenalty(candidate_it->second,
                                   grid,
                                   origin_x,
                                   origin_y,
                                   tile_size,
                                   selected_usage,
                                   soft_capacities);
      const int64_t candidate_total_cost
          = effectiveWirelengthCost(baseline_score, candidate_score, policy)
            + congestion_tradeoff * candidate_congestion_cost;

      if (candidate_total_cost < best_total_cost
          || (candidate_total_cost == best_total_cost
              && candidate_score.wirelength < best_score.wirelength)
          || (candidate_total_cost == best_total_cost
              && candidate_score.wirelength == best_score.wirelength
              && candidate_score.vias < best_score.vias)) {
        best_score = candidate_score;
        best_total_cost = candidate_total_cost;
        selected_route = &candidate_it->second;
        selected_source = source;
      }
    };

    consider(balanced_routes, RouteSource::kNewgrBalanced);
    consider(wirelength_routes, RouteSource::kNewgrWirelength);

    if (selected_route != &route) {
      route = *selected_route;
    }
    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);

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
                "NEWGR congestion-aware hybrid selected {} balanced and {} "
                "WL-first routes (kept {} FastRoute, +{} balanced-only, +{} "
                "WL-only) out of {} total.",
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
