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
  kNewgrWirelength,
  kNewgrDataWirelength,
  kNewgrRegionAware,
  kNewgrRegularRegion,
  kNewgrFineGrain,
  kNewgrSmallNet,
  kNewgrAstar,
  kNewgrRudy,
  kNewgrAstarEarly,
  kNewgrDetPartClassic,
  kNewgrNonDetHybrid,
  kNewgrRudyPartition
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

void removeRouteFromUsage(const GRoute& route,
                          int origin_x,
                          int origin_y,
                          int tile_size,
                          EdgeUsageMap& usage)
{
  forEachUnitPlanarEdge(
      route, origin_x, origin_y, tile_size, [&](const RouteEdgeKey& edge) {
        const auto usage_it = usage.find(edge);
        if (usage_it == usage.end()) {
          return;
        }
        if (usage_it->second <= 1) {
          usage.erase(usage_it);
          return;
        }
        usage_it->second--;
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

double congestionRiskDensity(const RouteScore& score,
                             int64_t congestion_penalty,
                             int tile_size)
{
  const double norm = static_cast<double>(std::max(1, tile_size));
  const double wl_units
      = static_cast<double>(std::max<int64_t>(1, score.wirelength)) / norm;
  return static_cast<double>(congestion_penalty) / std::max(1.0, wl_units);
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
    policy.via_guard = policy.via_guard * 7 + 28;
    policy.hard_via_guard = policy.via_guard * 2 + 36;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 3;
    policy.aggressive_wl_gain = std::max<int64_t>(1, tile_size / 5);
  } else if (baseline_score.wirelength >= medium_net_threshold) {
    policy.medium_net = true;
    policy.via_tradeoff = 1;
    policy.bend_tradeoff = std::max(1, tile_size / 180);
    policy.via_guard = policy.via_guard * 3 + 16;
    policy.hard_via_guard = policy.via_guard * 3 + 24;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 6;
    policy.aggressive_wl_gain = std::max<int64_t>(2, tile_size / 4);
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
  NetRouteMap data_wirelength_routes = engine_->runDataDrivenWirelength();
  NetRouteMap region_aware_routes = engine_->runRegionAware();
  NetRouteMap regular_region_routes = engine_->runRegularRegionAware();
  NetRouteMap finegrain_routes = engine_->runFineGrainRefine();
  NetRouteMap smallnet_routes = engine_->runSmallNetAware();
  NetRouteMap astar_routes = engine_->runAstarClassic();
  NetRouteMap rudy_routes = engine_->runRudyDriven();
  NetRouteMap astar_early_routes = engine_->runAstarEarly();
  NetRouteMap detpart_routes = engine_->runDetPartClassic();
  NetRouteMap nondet_routes = engine_->runNonDetHybrid();
  NetRouteMap rudy_partition_routes = engine_->runRudyPartition();

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

  int64_t baseline_total_vias = 0;
  for (const OrderedNet& ordered_net : ordered_nets) {
    baseline_total_vias += ordered_net.baseline_score.vias;
  }

  const int64_t global_base_via_budget
      = std::max<int64_t>(16, baseline_total_vias / 4000);
  const int64_t global_wl_to_via_credit = std::max<int64_t>(8, tile_size / 2);

  int selected_from_balanced = 0;
  int selected_from_wl = 0;
  int selected_from_data_wl = 0;
  int selected_from_region = 0;
  int selected_from_regular_region = 0;
  int selected_from_finegrain = 0;
  int selected_from_smallnet = 0;
  int selected_from_astar = 0;
  int selected_from_rudy = 0;
  int selected_from_astar_early = 0;
  int selected_from_detpart = 0;
  int selected_from_nondet = 0;
  int selected_from_rudy_partition = 0;
  int kept_fastroute = 0;
  int inserted_from_balanced = 0;
  int inserted_from_wl = 0;
  int inserted_from_data_wl = 0;
  int inserted_from_region = 0;
  int inserted_from_regular_region = 0;
  int inserted_from_finegrain = 0;
  int inserted_from_smallnet = 0;
  int inserted_from_astar = 0;
  int inserted_from_rudy = 0;
  int inserted_from_astar_early = 0;
  int inserted_from_detpart = 0;
  int inserted_from_nondet = 0;
  int inserted_from_rudy_partition = 0;

  EdgeUsageMap selected_usage;
  selected_usage.reserve(baseline_demand.size());
  int64_t cumulative_wl_gain = 0;
  int64_t cumulative_extra_vias = 0;
  const int top_trunk_nets
      = std::max<int>(1, static_cast<int>(ordered_nets.size() / 6));
  int net_rank = 0;

  for (const OrderedNet& ordered_net : ordered_nets) {
    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      continue;
    }

    GRoute& route = route_it->second;
    const RouteScore baseline_score = ordered_net.baseline_score;
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
    const bool ultra_wl_mode = net_rank < top_trunk_nets;
    net_rank++;
    const int64_t congestion_tradeoff = policy.long_net ? 0 : (policy.medium_net ? 0 : 1);

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

    auto consider = [&](const NetRouteMap& candidate_routes,
                        RouteSource source,
                        int64_t min_wl_drop_required,
                        bool exploratory_mode) {
      auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      if (min_wl_drop_required > 0
          && candidate_score.wirelength + min_wl_drop_required
                 > best_score.wirelength) {
        return;
      }
      if (!betterCandidate(
              baseline_score, best_score, candidate_score, policy)) {
        return;
      }

      const int64_t current_net_extra_vias
          = std::max<int64_t>(0, best_score.vias - baseline_score.vias);
      const int64_t candidate_net_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t current_net_wl_gain
          = std::max<int64_t>(0, baseline_score.wirelength - best_score.wirelength);
      const int64_t candidate_net_wl_gain
          = std::max<int64_t>(0, baseline_score.wirelength - candidate_score.wirelength);

      const int64_t prospective_total_extra_vias
          = cumulative_extra_vias - current_net_extra_vias + candidate_net_extra_vias;
      const int64_t prospective_total_wl_gain
          = cumulative_wl_gain - current_net_wl_gain + candidate_net_wl_gain;
      int64_t global_via_budget = global_base_via_budget;
      if (prospective_total_wl_gain > 0) {
        global_via_budget += prospective_total_wl_gain / global_wl_to_via_credit;
      }

      if (prospective_total_extra_vias > global_via_budget
          && candidate_net_extra_vias > current_net_extra_vias) {
        const int64_t wl_drop_vs_best = std::max<int64_t>(
            0, best_score.wirelength - candidate_score.wirelength);
        const int64_t emergency_budget = global_via_budget + policy.via_guard / 2;
        if (wl_drop_vs_best < std::max<int64_t>(tile_size / 2, 1)
            || prospective_total_extra_vias > emergency_budget) {
          return;
        }
      }

      const int64_t candidate_congestion_cost
          = routeCongestionPenalty(candidate_it->second,
                                   grid,
                                   origin_x,
                                   origin_y,
                                   tile_size,
                                   selected_usage,
                                   soft_capacities);
      if (exploratory_mode && candidate_congestion_cost > best_congestion_cost) {
        const int64_t wl_drop_vs_best
            = std::max<int64_t>(0, best_score.wirelength - candidate_score.wirelength);
        const int64_t congestion_delta = candidate_congestion_cost - best_congestion_cost;
        const int64_t allowed_delta = policy.long_net
                                          ? std::max<int64_t>(10, best_congestion_cost / 7)
                                          : std::max<int64_t>(6, best_congestion_cost / 12);
        const int64_t required_wl_drop = policy.long_net
                                             ? std::max<int64_t>(1, tile_size / 8)
                                             : std::max<int64_t>(1, tile_size / 5);
        if (congestion_delta > allowed_delta && wl_drop_vs_best < required_wl_drop) {
          return;
        }
      }
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
        best_congestion_cost = candidate_congestion_cost;
        best_total_cost = candidate_total_cost;
        selected_route = &candidate_it->second;
        selected_source = source;
      }
    };

    int64_t exploratory_min_wl_drop
        = policy.long_net
              ? std::max<int64_t>(1, tile_size / 10)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 8)
                                   : std::max<int64_t>(1, tile_size / 6));
    int64_t aggressive_min_wl_drop
        = policy.long_net
              ? std::max<int64_t>(1, tile_size / 12)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 10)
                                   : std::max<int64_t>(1, tile_size / 7));
    if (ultra_wl_mode) {
      exploratory_min_wl_drop = std::max<int64_t>(1, tile_size / 16);
      aggressive_min_wl_drop = std::max<int64_t>(1, tile_size / 18);
    }
    consider(balanced_routes, RouteSource::kNewgrBalanced, 0, false);
    consider(wirelength_routes, RouteSource::kNewgrWirelength, 0, false);
    consider(data_wirelength_routes,
             RouteSource::kNewgrDataWirelength,
             exploratory_min_wl_drop,
             true);
    consider(region_aware_routes,
             RouteSource::kNewgrRegionAware,
             exploratory_min_wl_drop,
             true);
    consider(regular_region_routes,
             RouteSource::kNewgrRegularRegion,
             exploratory_min_wl_drop,
             true);
    consider(finegrain_routes,
             RouteSource::kNewgrFineGrain,
             aggressive_min_wl_drop,
             true);
    consider(smallnet_routes,
             RouteSource::kNewgrSmallNet,
             aggressive_min_wl_drop,
             true);
    consider(astar_routes,
             RouteSource::kNewgrAstar,
             aggressive_min_wl_drop,
             true);
    consider(rudy_routes, RouteSource::kNewgrRudy, aggressive_min_wl_drop, true);
    consider(astar_early_routes,
             RouteSource::kNewgrAstarEarly,
             aggressive_min_wl_drop,
             true);
    consider(detpart_routes,
             RouteSource::kNewgrDetPartClassic,
             aggressive_min_wl_drop,
             true);
    consider(nondet_routes,
             RouteSource::kNewgrNonDetHybrid,
             aggressive_min_wl_drop,
             true);
    consider(rudy_partition_routes,
             RouteSource::kNewgrRudyPartition,
             aggressive_min_wl_drop,
             true);

    // Wirelength champion pass: if one candidate has a material WL gain and
    // stays within manageable congestion/via deltas, force-select it.
    const GRoute* wl_champion_route = selected_route;
    RouteScore wl_champion_score = best_score;
    int64_t wl_champion_congestion_cost = best_congestion_cost;
    RouteSource wl_champion_source = selected_source;
    auto maybeUpdateChampion = [&](const NetRouteMap& candidate_routes,
                                   RouteSource source) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      if (candidate_extra_vias > policy.hard_via_guard * 2 + 24) {
        return;
      }
      if (candidate_score.wirelength >= wl_champion_score.wirelength) {
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
      wl_champion_route = &candidate_it->second;
      wl_champion_score = candidate_score;
      wl_champion_congestion_cost = candidate_congestion_cost;
      wl_champion_source = source;
    };

    maybeUpdateChampion(balanced_routes, RouteSource::kNewgrBalanced);
    maybeUpdateChampion(wirelength_routes, RouteSource::kNewgrWirelength);
    maybeUpdateChampion(data_wirelength_routes, RouteSource::kNewgrDataWirelength);
    maybeUpdateChampion(region_aware_routes, RouteSource::kNewgrRegionAware);
    maybeUpdateChampion(regular_region_routes, RouteSource::kNewgrRegularRegion);
    maybeUpdateChampion(finegrain_routes, RouteSource::kNewgrFineGrain);
    maybeUpdateChampion(smallnet_routes, RouteSource::kNewgrSmallNet);
    maybeUpdateChampion(astar_routes, RouteSource::kNewgrAstar);
    maybeUpdateChampion(rudy_routes, RouteSource::kNewgrRudy);
    maybeUpdateChampion(astar_early_routes, RouteSource::kNewgrAstarEarly);
    maybeUpdateChampion(detpart_routes, RouteSource::kNewgrDetPartClassic);
    maybeUpdateChampion(nondet_routes, RouteSource::kNewgrNonDetHybrid);
    maybeUpdateChampion(rudy_partition_routes, RouteSource::kNewgrRudyPartition);

    if (wl_champion_route != selected_route) {
      const int64_t wl_drop_vs_best
          = std::max<int64_t>(0, best_score.wirelength - wl_champion_score.wirelength);
      const int64_t congestion_delta
          = wl_champion_congestion_cost - best_congestion_cost;
      const int64_t allowed_congestion_delta
          = policy.long_net ? std::max<int64_t>(14, best_congestion_cost / 4)
                            : (policy.medium_net
                                   ? std::max<int64_t>(10, best_congestion_cost / 5)
                                   : std::max<int64_t>(6, best_congestion_cost / 8));
      const int64_t champion_min_wl_gain
          = policy.long_net ? std::max<int64_t>(1, tile_size / 8)
                            : (policy.medium_net
                                   ? std::max<int64_t>(1, tile_size / 6)
                                   : std::max<int64_t>(1, tile_size / 5));
      const int64_t champion_force_gain
          = policy.long_net ? std::max<int64_t>(2, tile_size / 3)
                            : std::max<int64_t>(2, tile_size / 2);
      const double best_congestion_density
          = congestionRiskDensity(best_score, best_congestion_cost, tile_size);
      const double champion_congestion_density
          = congestionRiskDensity(
              wl_champion_score, wl_champion_congestion_cost, tile_size);
      const double allowed_density_ratio = policy.long_net ? 1.18 : 1.24;
      const bool hotspot_density_safe
          = champion_congestion_density
                <= best_congestion_density * allowed_density_ratio
            || wl_drop_vs_best >= champion_force_gain;

      const int64_t current_net_extra_vias
          = std::max<int64_t>(0, best_score.vias - baseline_score.vias);
      const int64_t champion_net_extra_vias
          = std::max<int64_t>(0, wl_champion_score.vias - baseline_score.vias);
      const int64_t current_net_wl_gain
          = std::max<int64_t>(0, baseline_score.wirelength - best_score.wirelength);
      const int64_t champion_net_wl_gain
          = std::max<int64_t>(0, baseline_score.wirelength - wl_champion_score.wirelength);
      const int64_t prospective_total_extra_vias
          = cumulative_extra_vias - current_net_extra_vias + champion_net_extra_vias;
      const int64_t prospective_total_wl_gain
          = cumulative_wl_gain - current_net_wl_gain + champion_net_wl_gain;
      int64_t global_via_budget = global_base_via_budget;
      if (prospective_total_wl_gain > 0) {
        global_via_budget += prospective_total_wl_gain / global_wl_to_via_credit;
      }

      // Long trunk override: for high-WL nets, prioritize shorter trunks even
      // when congestion/via projections worsen moderately.
      const int64_t long_override_gain = std::max<int64_t>(1, tile_size / 8);
      const int64_t trunk_override_gain = std::max<int64_t>(1, tile_size / 4);
      const int64_t champion_via_cap
          = policy.hard_via_guard * 2 + (policy.long_net ? 48 : 24);
      const bool long_wl_override
          = policy.long_net && wl_drop_vs_best >= long_override_gain
            && champion_net_extra_vias <= champion_via_cap;
      const bool medium_trunk_override
          = policy.medium_net && wl_drop_vs_best >= trunk_override_gain
            && champion_net_extra_vias <= champion_via_cap
            && congestion_delta <= allowed_congestion_delta;
      const int64_t ultra_min_gain = std::max<int64_t>(1, tile_size / 12);
      const int64_t ultra_force_gain = std::max<int64_t>(2, tile_size / 5);
      const int64_t ultra_via_cap = policy.hard_via_guard * 3 + 80;
      const bool ultra_trunk_override
          = ultra_wl_mode && wl_drop_vs_best >= ultra_min_gain
            && champion_net_extra_vias <= ultra_via_cap
            && (congestion_delta <= allowed_congestion_delta * 2
                || wl_drop_vs_best >= ultra_force_gain);

      if ((wl_drop_vs_best >= champion_min_wl_gain
          && (congestion_delta <= allowed_congestion_delta
              || wl_drop_vs_best >= champion_force_gain)
          && hotspot_density_safe
          && (prospective_total_extra_vias
                  <= global_via_budget + policy.via_guard
              || wl_drop_vs_best >= champion_force_gain))
          || long_wl_override || medium_trunk_override || ultra_trunk_override) {
        best_score = wl_champion_score;
        best_congestion_cost = wl_champion_congestion_cost;
        best_total_cost = effectiveWirelengthCost(
                              baseline_score, wl_champion_score, policy)
                          + congestion_tradeoff * wl_champion_congestion_cost;
        selected_route = wl_champion_route;
        selected_source = wl_champion_source;
      }
    }

    if (selected_route != &route) {
      route = *selected_route;
    }
    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);
    cumulative_wl_gain += std::max<int64_t>(
        0, baseline_score.wirelength - best_score.wirelength);
    cumulative_extra_vias += std::max<int64_t>(
        0, best_score.vias - baseline_score.vias);

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
      case RouteSource::kNewgrDataWirelength:
        selected_from_data_wl++;
        break;
      case RouteSource::kNewgrRegionAware:
        selected_from_region++;
        break;
      case RouteSource::kNewgrRegularRegion:
        selected_from_regular_region++;
        break;
      case RouteSource::kNewgrFineGrain:
        selected_from_finegrain++;
        break;
      case RouteSource::kNewgrSmallNet:
        selected_from_smallnet++;
        break;
      case RouteSource::kNewgrAstar:
        selected_from_astar++;
        break;
      case RouteSource::kNewgrRudy:
        selected_from_rudy++;
        break;
      case RouteSource::kNewgrAstarEarly:
        selected_from_astar_early++;
        break;
      case RouteSource::kNewgrDetPartClassic:
        selected_from_detpart++;
        break;
      case RouteSource::kNewgrNonDetHybrid:
        selected_from_nondet++;
        break;
      case RouteSource::kNewgrRudyPartition:
        selected_from_rudy_partition++;
        break;
    }
  }

  int second_pass_swaps = 0;
  int64_t second_pass_wl_gain = 0;
  int64_t second_pass_via_delta = 0;
  for (const OrderedNet& ordered_net : ordered_nets) {
    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      continue;
    }

    GRoute& route = route_it->second;
    const RouteScore baseline_score = ordered_net.baseline_score;
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
    const RouteScore current_score = scoreRoute(route);

    removeRouteFromUsage(route, origin_x, origin_y, tile_size, selected_usage);
    const int64_t current_congestion_cost
        = routeCongestionPenalty(route,
                                 grid,
                                 origin_x,
                                 origin_y,
                                 tile_size,
                                 selected_usage,
                                 soft_capacities);

    const GRoute* best_refine_route = &route;
    RouteScore best_refine_score = current_score;
    int64_t best_refine_congestion_cost = current_congestion_cost;

    auto consider_refine = [&](const NetRouteMap& candidate_routes) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }

      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      const int64_t wl_drop_vs_current = std::max<int64_t>(
          0, current_score.wirelength - candidate_score.wirelength);
      if (wl_drop_vs_current == 0) {
        return;
      }
      if (candidate_score.wirelength >= best_refine_score.wirelength) {
        return;
      }

      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t allowed_total_extra_vias
          = policy.hard_via_guard * 2 + (policy.long_net ? 48 : 24);
      if (candidate_extra_vias > allowed_total_extra_vias) {
        return;
      }

      const int64_t via_delta = candidate_score.vias - current_score.vias;
      const int64_t min_wl_gain
          = policy.long_net ? std::max<int64_t>(1, tile_size / 9)
                            : (policy.medium_net
                                   ? std::max<int64_t>(1, tile_size / 7)
                                   : std::max<int64_t>(1, tile_size / 5));
      const int64_t force_wl_gain
          = policy.long_net ? std::max<int64_t>(2, tile_size / 2)
                            : std::max<int64_t>(2, tile_size / 3);
      const bool via_safe
          = via_delta <= std::max<int64_t>(2, policy.via_guard / 2)
            || wl_drop_vs_current
                   >= std::max<int64_t>(
                       1, via_delta * std::max<int64_t>(1, policy.wl_per_extra_via / 2));
      if (wl_drop_vs_current < min_wl_gain && !via_safe) {
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
      const int64_t congestion_delta
          = candidate_congestion_cost - current_congestion_cost;
      const int64_t allowed_congestion_delta
          = policy.long_net ? std::max<int64_t>(16, current_congestion_cost / 3)
                            : (policy.medium_net
                                   ? std::max<int64_t>(12, current_congestion_cost / 4)
                                   : std::max<int64_t>(8, current_congestion_cost / 5));
      if (congestion_delta > allowed_congestion_delta
          && wl_drop_vs_current < force_wl_gain) {
        return;
      }

      if (candidate_score.wirelength < best_refine_score.wirelength
          || (candidate_score.wirelength == best_refine_score.wirelength
              && candidate_congestion_cost < best_refine_congestion_cost)
          || (candidate_score.wirelength == best_refine_score.wirelength
              && candidate_congestion_cost == best_refine_congestion_cost
              && candidate_score.vias < best_refine_score.vias)) {
        best_refine_route = &candidate_it->second;
        best_refine_score = candidate_score;
        best_refine_congestion_cost = candidate_congestion_cost;
      }
    };

    consider_refine(balanced_routes);
    consider_refine(wirelength_routes);
    consider_refine(data_wirelength_routes);
    consider_refine(region_aware_routes);
    consider_refine(regular_region_routes);
    consider_refine(finegrain_routes);
    consider_refine(smallnet_routes);
    consider_refine(astar_routes);
    consider_refine(rudy_routes);
    consider_refine(astar_early_routes);
    consider_refine(detpart_routes);
    consider_refine(nondet_routes);
    consider_refine(rudy_partition_routes);

    if (best_refine_route != &route) {
      second_pass_swaps++;
      second_pass_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_refine_score.wirelength);

      const int64_t current_extra_vias
          = std::max<int64_t>(0, current_score.vias - baseline_score.vias);
      const int64_t refined_extra_vias
          = std::max<int64_t>(0, best_refine_score.vias - baseline_score.vias);
      const int64_t extra_via_delta = refined_extra_vias - current_extra_vias;
      second_pass_via_delta += extra_via_delta;
      cumulative_extra_vias += extra_via_delta;
      cumulative_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_refine_score.wirelength);

      route = *best_refine_route;
    }

    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);
  }

  int third_pass_swaps = 0;
  int64_t third_pass_wl_gain = 0;
  int64_t third_pass_via_delta = 0;
  int third_pass_rank = 0;
  for (const OrderedNet& ordered_net : ordered_nets) {
    if (third_pass_rank >= top_trunk_nets) {
      break;
    }
    third_pass_rank++;

    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      continue;
    }

    GRoute& route = route_it->second;
    const RouteScore baseline_score = ordered_net.baseline_score;
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
    const RouteScore current_score = scoreRoute(route);

    removeRouteFromUsage(route, origin_x, origin_y, tile_size, selected_usage);
    const int64_t current_congestion_cost
        = routeCongestionPenalty(route,
                                 grid,
                                 origin_x,
                                 origin_y,
                                 tile_size,
                                 selected_usage,
                                 soft_capacities);

    const GRoute* best_rescue_route = &route;
    RouteScore best_rescue_score = current_score;
    auto consider_rescue = [&](const NetRouteMap& candidate_routes) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }

      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      const int64_t wl_drop_vs_current = std::max<int64_t>(
          0, current_score.wirelength - candidate_score.wirelength);
      if (wl_drop_vs_current < std::max<int64_t>(1, tile_size / 14)) {
        return;
      }
      if (candidate_score.wirelength >= best_rescue_score.wirelength) {
        return;
      }

      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t rescue_via_cap = policy.hard_via_guard * 3 + 80;
      if (candidate_extra_vias > rescue_via_cap) {
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
      const int64_t congestion_delta
          = candidate_congestion_cost - current_congestion_cost;
      const int64_t allowed_congestion_delta
          = std::max<int64_t>(20, current_congestion_cost / 2);
      const int64_t force_wl_gain = std::max<int64_t>(2, tile_size / 5);
      if (congestion_delta > allowed_congestion_delta
          && wl_drop_vs_current < force_wl_gain) {
        return;
      }

      best_rescue_route = &candidate_it->second;
      best_rescue_score = candidate_score;
    };

    consider_rescue(balanced_routes);
    consider_rescue(wirelength_routes);
    consider_rescue(data_wirelength_routes);
    consider_rescue(region_aware_routes);
    consider_rescue(regular_region_routes);
    consider_rescue(finegrain_routes);
    consider_rescue(smallnet_routes);
    consider_rescue(astar_routes);
    consider_rescue(rudy_routes);
    consider_rescue(astar_early_routes);
    consider_rescue(detpart_routes);
    consider_rescue(nondet_routes);
    consider_rescue(rudy_partition_routes);

    if (best_rescue_route != &route) {
      third_pass_swaps++;
      third_pass_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_rescue_score.wirelength);

      const int64_t current_extra_vias
          = std::max<int64_t>(0, current_score.vias - baseline_score.vias);
      const int64_t rescue_extra_vias
          = std::max<int64_t>(0, best_rescue_score.vias - baseline_score.vias);
      const int64_t extra_via_delta = rescue_extra_vias - current_extra_vias;
      third_pass_via_delta += extra_via_delta;
      cumulative_extra_vias += extra_via_delta;
      cumulative_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_rescue_score.wirelength);

      route = *best_rescue_route;
    }

    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);
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
  for (const auto& [db_net, route] : data_wirelength_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_data_wl++;
    }
  }
  for (const auto& [db_net, route] : region_aware_routes) {
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
  for (const auto& [db_net, route] : smallnet_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_smallnet++;
    }
  }
  for (const auto& [db_net, route] : astar_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_astar++;
    }
  }
  for (const auto& [db_net, route] : rudy_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_rudy++;
    }
  }
  for (const auto& [db_net, route] : astar_early_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_astar_early++;
    }
  }
  for (const auto& [db_net, route] : detpart_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_detpart++;
    }
  }
  for (const auto& [db_net, route] : nondet_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_nondet++;
    }
  }
  for (const auto& [db_net, route] : rudy_partition_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_rudy_partition++;
    }
  }
  logger_->info(utl::GRT,
                6004,
                "NEWGR WL-priority hybrid selected balanced={} wl={} data={} "
                "region={} regular={} fine={} small={} astar={} rudy={} "
                "astarEarly={} detpart={} nondet={} rudyPart={} "
                "(kept FR={}; +balanced={} "
                "+wl={} +data={} +region={} +regular={} +fine={} +small={} "
                "+astar={} +rudy={} +astarEarly={} +detpart={} +nondet={} "
                "+rudyPart={}). "
                "Refine swaps={} wl-gain={} via-delta={}. "
                "Rescue swaps={} wl-gain={} via-delta={}. "
                "Global WL gain={} "
                "extra-vias={} (base via budget={} + gain/{}) out of {} total.",
                selected_from_balanced,
                selected_from_wl,
                selected_from_data_wl,
                selected_from_region,
                selected_from_regular_region,
                selected_from_finegrain,
                selected_from_smallnet,
                selected_from_astar,
                selected_from_rudy,
                selected_from_astar_early,
                selected_from_detpart,
                selected_from_nondet,
                selected_from_rudy_partition,
                kept_fastroute,
                inserted_from_balanced,
                inserted_from_wl,
                inserted_from_data_wl,
                inserted_from_region,
                inserted_from_regular_region,
                inserted_from_finegrain,
                inserted_from_smallnet,
                inserted_from_astar,
                inserted_from_rudy,
                inserted_from_astar_early,
                inserted_from_detpart,
                inserted_from_nondet,
                inserted_from_rudy_partition,
                second_pass_swaps,
                second_pass_wl_gain,
                second_pass_via_delta,
                third_pass_swaps,
                third_pass_wl_gain,
                third_pass_via_delta,
                cumulative_wl_gain,
                cumulative_extra_vias,
                global_base_via_budget,
                global_wl_to_via_credit,
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
