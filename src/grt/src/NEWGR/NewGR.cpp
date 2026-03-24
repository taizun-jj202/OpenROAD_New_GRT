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
  kNewgrRudyPartition,
  kNewgrLegacySqueeze,
  kNewgrDirectWirelength,
  kNewgrUltraDirectWirelength
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
using EdgeHistoryMap = std::unordered_map<RouteEdgeKey, double, RouteEdgeKeyHash>;

uint64_t routeFingerprint(const GRoute& route)
{
  uint64_t hash = 1469598103934665603ull;
  for (const GSegment& segment : route) {
    const uint64_t values[] = {
        static_cast<uint64_t>(static_cast<uint32_t>(segment.init_x)),
        static_cast<uint64_t>(static_cast<uint32_t>(segment.init_y)),
        static_cast<uint64_t>(static_cast<uint32_t>(segment.init_layer)),
        static_cast<uint64_t>(static_cast<uint32_t>(segment.final_x)),
        static_cast<uint64_t>(static_cast<uint32_t>(segment.final_y)),
        static_cast<uint64_t>(static_cast<uint32_t>(segment.final_layer))};
    for (uint64_t value : values) {
      hash ^= value + 0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2);
      hash *= 1099511628211ull;
    }
  }
  return hash;
}

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

double routeHistoryPenalty(const GRoute& route,
                           int origin_x,
                           int origin_y,
                           int tile_size,
                           const EdgeHistoryMap& history)
{
  double penalty = 0.0;
  forEachUnitPlanarEdge(
      route, origin_x, origin_y, tile_size, [&](const RouteEdgeKey& edge) {
        const auto history_it = history.find(edge);
        if (history_it != history.end()) {
          penalty += history_it->second;
        }
      });
  return penalty;
}

int64_t totalSoftOverflow(const EdgeUsageMap& usage,
                          const SprouteGridData& grid,
                          const EdgeUsageMap& soft_capacities)
{
  int64_t total_overflow = 0;
  for (const auto& [edge, edge_usage] : usage) {
    int soft_capacity = edgeHardCapacity(edge, grid);
    const auto soft_it = soft_capacities.find(edge);
    if (soft_it != soft_capacities.end()) {
      soft_capacity = soft_it->second;
    } else {
      soft_capacity = softCapacityFromDemand(edge, soft_capacity, 0);
    }
    total_overflow += std::max(0, edge_usage - soft_capacity);
  }
  return total_overflow;
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
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 12;
  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 5;
  if (baseline_score.wirelength >= long_net_threshold) {
    return std::max<int64_t>(20, baseline_score.vias / 2 + 14);
  }
  if (baseline_score.wirelength >= medium_net_threshold) {
    return std::max<int64_t>(9, baseline_score.vias / 3 + 7);
  }
  return std::max<int64_t>(2, baseline_score.vias / 6 + 2);
}

SelectionPolicy buildSelectionPolicy(const RouteScore& baseline_score, int tile_size)
{
  SelectionPolicy policy;
  policy.via_tradeoff = 0;
  policy.bend_tradeoff = std::max(1, tile_size / 420);
  policy.via_guard = viaGuardForNet(baseline_score, tile_size) + 10;
  policy.hard_via_guard = policy.via_guard * 8 + 48;
  policy.min_wl_improve = 1;
  policy.wl_per_extra_via = 4;
  policy.aggressive_wl_gain = std::max<int64_t>(2, tile_size / 8);

  const int64_t medium_net_threshold = static_cast<int64_t>(tile_size) * 5;
  const int64_t long_net_threshold = static_cast<int64_t>(tile_size) * 10;
  if (baseline_score.wirelength >= long_net_threshold) {
    policy.long_net = true;
    policy.via_tradeoff = 0;
    policy.bend_tradeoff = std::max(1, tile_size / 620);
    policy.via_guard = policy.via_guard * 10 + 90;
    policy.hard_via_guard = policy.via_guard * 3 + 160;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 1;
    policy.aggressive_wl_gain = std::max<int64_t>(1, tile_size / 18);
  } else if (baseline_score.wirelength >= medium_net_threshold) {
    policy.medium_net = true;
    policy.via_tradeoff = 0;
    policy.bend_tradeoff = std::max(1, tile_size / 520);
    policy.via_guard = policy.via_guard * 6 + 48;
    policy.hard_via_guard = policy.via_guard * 4 + 96;
    policy.min_wl_improve = 1;
    policy.wl_per_extra_via = 2;
    policy.aggressive_wl_gain = std::max<int64_t>(1, tile_size / 14);
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
  NetRouteMap rudy_classic_routes = engine_->runRudyClassic();
  NetRouteMap pin_density_routes = engine_->runPinDensityClassic();
  NetRouteMap rudy_pin_hybrid_routes = engine_->runRudyPinHybrid();
  NetRouteMap astar_early_routes = engine_->runAstarEarly();
  NetRouteMap detpart_routes = engine_->runDetPartClassic();
  NetRouteMap nondet_routes = engine_->runNonDetHybrid();
  NetRouteMap rudy_partition_routes = engine_->runRudyPartition();
  NetRouteMap legacy_squeeze_routes = engine_->runLegacySqueeze();
  NetRouteMap direct_wl_routes = engine_->runDirectWirelength();
  NetRouteMap ultra_direct_wl_routes = engine_->runUltraDirectWirelength();

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
      = std::max<int64_t>(36, baseline_total_vias / 1800);
  const int64_t global_wl_to_via_credit = std::max<int64_t>(2, tile_size / 8);

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
  int selected_from_legacy_squeeze = 0;
  int selected_from_direct_wl = 0;
  int selected_from_ultra_direct_wl = 0;
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
  int inserted_from_legacy_squeeze = 0;
  int inserted_from_direct_wl = 0;
  int inserted_from_ultra_direct_wl = 0;

  EdgeUsageMap selected_usage;
  selected_usage.reserve(baseline_demand.size());
  int64_t cumulative_wl_gain = 0;
  int64_t cumulative_extra_vias = 0;
  const int top_trunk_nets
      = std::max<int>(1, static_cast<int>(ordered_nets.size() / 2));
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
    const int64_t congestion_tradeoff = 0;

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
      if (ultra_wl_mode) {
        global_via_budget += std::max<int64_t>(policy.via_guard * 2, tile_size / 2);
      }

      if (prospective_total_extra_vias > global_via_budget
          && candidate_net_extra_vias > current_net_extra_vias) {
        const int64_t wl_drop_vs_best = std::max<int64_t>(
            0, best_score.wirelength - candidate_score.wirelength);
        const int64_t emergency_budget = global_via_budget + policy.via_guard / 2;
        const bool ultra_override
            = ultra_wl_mode
              && wl_drop_vs_best >= std::max<int64_t>(1, tile_size / 14)
              && prospective_total_extra_vias
                     <= emergency_budget + std::max<int64_t>(policy.via_guard, 16);
        if (ultra_override) {
          // Keep trunk routes short even when projected via usage rises.
        } else if (wl_drop_vs_best < std::max<int64_t>(tile_size / 20, 1)
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
        const int64_t allowed_delta = ultra_wl_mode
                                          ? std::max<int64_t>(20, best_congestion_cost / 3)
                                      : policy.long_net
                                          ? std::max<int64_t>(10, best_congestion_cost / 7)
                                          : std::max<int64_t>(6, best_congestion_cost / 12);
        const int64_t required_wl_drop = policy.long_net
                                             ? std::max<int64_t>(1, tile_size / 24)
                                         : ultra_wl_mode
                                             ? std::max<int64_t>(1, tile_size / 42)
                                             : std::max<int64_t>(1, tile_size / 20);
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
              ? std::max<int64_t>(1, tile_size / 56)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 52)
                                   : std::max<int64_t>(1, tile_size / 48));
    int64_t aggressive_min_wl_drop
        = policy.long_net
              ? std::max<int64_t>(1, tile_size / 68)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 64)
                                   : std::max<int64_t>(1, tile_size / 60));
    int64_t direct_min_wl_drop
        = policy.long_net
              ? std::max<int64_t>(1, tile_size / 96)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 90)
                                   : std::max<int64_t>(1, tile_size / 84));
    int64_t ultra_direct_min_wl_drop
        = policy.long_net
              ? std::max<int64_t>(1, tile_size / 120)
              : (policy.medium_net ? std::max<int64_t>(1, tile_size / 112)
                                   : std::max<int64_t>(1, tile_size / 104));
    if (ultra_wl_mode) {
      exploratory_min_wl_drop = std::max<int64_t>(1, tile_size / 120);
      aggressive_min_wl_drop = std::max<int64_t>(1, tile_size / 140);
      direct_min_wl_drop = std::max<int64_t>(1, tile_size / 180);
      ultra_direct_min_wl_drop = std::max<int64_t>(1, tile_size / 220);
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
    consider(pin_density_routes,
             RouteSource::kNewgrDataWirelength,
             aggressive_min_wl_drop,
             true);
    consider(astar_routes,
             RouteSource::kNewgrAstar,
             aggressive_min_wl_drop,
             true);
    consider(rudy_routes, RouteSource::kNewgrRudy, aggressive_min_wl_drop, true);
    consider(rudy_classic_routes,
             RouteSource::kNewgrRudy,
             aggressive_min_wl_drop,
             true);
    consider(rudy_pin_hybrid_routes,
             RouteSource::kNewgrRudy,
             direct_min_wl_drop,
             true);
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
    consider(legacy_squeeze_routes,
             RouteSource::kNewgrLegacySqueeze,
             exploratory_min_wl_drop,
             true);
    consider(direct_wl_routes,
             RouteSource::kNewgrDirectWirelength,
             direct_min_wl_drop,
             true);
    consider(ultra_direct_wl_routes,
             RouteSource::kNewgrUltraDirectWirelength,
             ultra_direct_min_wl_drop,
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
    maybeUpdateChampion(pin_density_routes, RouteSource::kNewgrDataWirelength);
    maybeUpdateChampion(astar_routes, RouteSource::kNewgrAstar);
    maybeUpdateChampion(rudy_routes, RouteSource::kNewgrRudy);
    maybeUpdateChampion(rudy_classic_routes, RouteSource::kNewgrRudy);
    maybeUpdateChampion(rudy_pin_hybrid_routes, RouteSource::kNewgrRudy);
    maybeUpdateChampion(astar_early_routes, RouteSource::kNewgrAstarEarly);
    maybeUpdateChampion(detpart_routes, RouteSource::kNewgrDetPartClassic);
    maybeUpdateChampion(nondet_routes, RouteSource::kNewgrNonDetHybrid);
    maybeUpdateChampion(rudy_partition_routes, RouteSource::kNewgrRudyPartition);
    maybeUpdateChampion(legacy_squeeze_routes, RouteSource::kNewgrLegacySqueeze);
    maybeUpdateChampion(direct_wl_routes, RouteSource::kNewgrDirectWirelength);
    maybeUpdateChampion(ultra_direct_wl_routes,
                        RouteSource::kNewgrUltraDirectWirelength);

    if (wl_champion_route != selected_route) {
      const int64_t wl_drop_vs_best
          = std::max<int64_t>(0, best_score.wirelength - wl_champion_score.wirelength);
      const int64_t congestion_delta
          = wl_champion_congestion_cost - best_congestion_cost;
      const int64_t allowed_congestion_delta
          = ultra_wl_mode ? std::max<int64_t>(26, best_congestion_cost / 2)
                          : policy.long_net
                                ? std::max<int64_t>(14, best_congestion_cost / 4)
                                : (policy.medium_net
                                       ? std::max<int64_t>(10, best_congestion_cost / 5)
                                       : std::max<int64_t>(6, best_congestion_cost / 8));
      const int64_t champion_min_wl_gain
          = policy.long_net ? std::max<int64_t>(1, tile_size / 36)
                            : (policy.medium_net
                                   ? std::max<int64_t>(1, tile_size / 44)
                                   : std::max<int64_t>(1, tile_size / 52));
      const int64_t champion_force_gain
          = ultra_wl_mode ? std::max<int64_t>(1, tile_size / 16)
                          : policy.long_net ? std::max<int64_t>(1, tile_size / 14)
                                            : std::max<int64_t>(1, tile_size / 12);
      const double best_congestion_density
          = congestionRiskDensity(best_score, best_congestion_cost, tile_size);
      const double champion_congestion_density
          = congestionRiskDensity(
              wl_champion_score, wl_champion_congestion_cost, tile_size);
      const double allowed_density_ratio
          = ultra_wl_mode ? 1.80 : (policy.long_net ? 1.35 : 1.50);
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
      if (ultra_wl_mode) {
        global_via_budget += std::max<int64_t>(policy.via_guard * 2, tile_size / 2);
      }

      // Long trunk override: for high-WL nets, prioritize shorter trunks even
      // when congestion/via projections worsen moderately.
      const int64_t long_override_gain = std::max<int64_t>(1, tile_size / 8);
      const int64_t trunk_override_gain = std::max<int64_t>(1, tile_size / 4);
      const int64_t champion_via_cap
          = policy.hard_via_guard * 3 + (policy.long_net ? 80 : 36);
      const bool long_wl_override
          = policy.long_net && wl_drop_vs_best >= long_override_gain
            && champion_net_extra_vias <= champion_via_cap;
      const bool medium_trunk_override
          = policy.medium_net && wl_drop_vs_best >= trunk_override_gain
            && champion_net_extra_vias <= champion_via_cap
            && congestion_delta <= allowed_congestion_delta;
      const int64_t ultra_min_gain = std::max<int64_t>(1, tile_size / 12);
      const int64_t ultra_force_gain = std::max<int64_t>(2, tile_size / 5);
      const int64_t ultra_via_cap = policy.hard_via_guard * 4 + 180;
      const bool ultra_trunk_override
          = ultra_wl_mode && wl_drop_vs_best >= ultra_min_gain
            && champion_net_extra_vias <= ultra_via_cap
            && (congestion_delta <= allowed_congestion_delta * 2
                || wl_drop_vs_best >= ultra_force_gain);
      const bool ultra_direct_override
          = ultra_wl_mode
            && wl_drop_vs_best >= std::max<int64_t>(1, tile_size / 16)
            && champion_net_extra_vias <= policy.hard_via_guard * 6 + 240
            && (congestion_delta <= std::max<int64_t>(32, best_congestion_cost)
                || wl_drop_vs_best >= std::max<int64_t>(2, tile_size / 6));

      if ((wl_drop_vs_best >= champion_min_wl_gain
          && (congestion_delta <= allowed_congestion_delta
              || wl_drop_vs_best >= champion_force_gain)
          && hotspot_density_safe
          && (prospective_total_extra_vias
                  <= global_via_budget + policy.via_guard
              || wl_drop_vs_best >= champion_force_gain))
          || long_wl_override || medium_trunk_override || ultra_trunk_override
          || ultra_direct_override) {
        best_score = wl_champion_score;
        best_congestion_cost = wl_champion_congestion_cost;
        best_total_cost = effectiveWirelengthCost(
                              baseline_score, wl_champion_score, policy)
                          + congestion_tradeoff * wl_champion_congestion_cost;
        selected_route = wl_champion_route;
        selected_source = wl_champion_source;
      }
    }

    // Hard WL clamp for high-priority trunks: for long and top-ranked nets,
    // force-pick the shortest direct-style candidate if it is meaningfully
    // shorter and within a very loose via cap.
    if (ultra_wl_mode || policy.long_net || policy.medium_net) {
      const GRoute* forced_route = selected_route;
      RouteScore forced_score = best_score;
      RouteSource forced_source = selected_source;

      auto consider_forced = [&](const NetRouteMap& candidate_routes,
                                 RouteSource source) {
        const auto candidate_it = candidate_routes.find(ordered_net.db_net);
        if (candidate_it == candidate_routes.end()) {
          return;
        }
        const RouteScore candidate_score = scoreRoute(candidate_it->second);
        if (candidate_score.wirelength >= forced_score.wirelength) {
          return;
        }
        const int64_t candidate_extra_vias
            = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
        const int64_t forced_via_cap
            = policy.hard_via_guard * 10 + (ultra_wl_mode ? 560 : 360);
        if (candidate_extra_vias > forced_via_cap) {
          return;
        }
        forced_route = &candidate_it->second;
        forced_score = candidate_score;
        forced_source = source;
      };

      consider_forced(direct_wl_routes, RouteSource::kNewgrDirectWirelength);
      consider_forced(ultra_direct_wl_routes,
                      RouteSource::kNewgrUltraDirectWirelength);
      consider_forced(rudy_classic_routes, RouteSource::kNewgrRudy);
      consider_forced(rudy_pin_hybrid_routes, RouteSource::kNewgrRudy);

      const int64_t forced_wl_gain
          = std::max<int64_t>(0, best_score.wirelength - forced_score.wirelength);
      const int64_t forced_min_gain
          = ultra_wl_mode ? std::max<int64_t>(1, tile_size / 280)
                          : std::max<int64_t>(1, tile_size / 180);
      if (forced_route != selected_route && forced_wl_gain >= forced_min_gain) {
        selected_route = forced_route;
        selected_source = forced_source;
        best_score = forced_score;
        best_congestion_cost = routeCongestionPenalty(*forced_route,
                                                      grid,
                                                      origin_x,
                                                      origin_y,
                                                      tile_size,
                                                      selected_usage,
                                                      soft_capacities);
        best_total_cost = effectiveWirelengthCost(
                              baseline_score, forced_score, policy)
                          + congestion_tradeoff * best_congestion_cost;
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
      case RouteSource::kNewgrLegacySqueeze:
        selected_from_legacy_squeeze++;
        break;
      case RouteSource::kNewgrDirectWirelength:
        selected_from_direct_wl++;
        break;
      case RouteSource::kNewgrUltraDirectWirelength:
        selected_from_ultra_direct_wl++;
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
          = policy.long_net ? std::max<int64_t>(1, tile_size / 32)
                            : (policy.medium_net
                                   ? std::max<int64_t>(1, tile_size / 38)
                                   : std::max<int64_t>(1, tile_size / 44));
      const int64_t force_wl_gain
          = policy.long_net ? std::max<int64_t>(1, tile_size / 11)
                            : std::max<int64_t>(1, tile_size / 10);
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
    consider_refine(pin_density_routes);
    consider_refine(region_aware_routes);
    consider_refine(regular_region_routes);
    consider_refine(finegrain_routes);
    consider_refine(smallnet_routes);
    consider_refine(astar_routes);
    consider_refine(rudy_routes);
    consider_refine(rudy_classic_routes);
    consider_refine(rudy_pin_hybrid_routes);
    consider_refine(astar_early_routes);
    consider_refine(detpart_routes);
    consider_refine(nondet_routes);
    consider_refine(rudy_partition_routes);
    consider_refine(legacy_squeeze_routes);
    consider_refine(direct_wl_routes);
    consider_refine(ultra_direct_wl_routes);

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
      if (wl_drop_vs_current < std::max<int64_t>(1, tile_size / 80)) {
        return;
      }
      if (candidate_score.wirelength >= best_rescue_score.wirelength) {
        return;
      }

      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t rescue_via_cap = policy.hard_via_guard * 5 + 220;
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
      const int64_t force_wl_gain = std::max<int64_t>(1, tile_size / 16);
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
    consider_rescue(pin_density_routes);
    consider_rescue(region_aware_routes);
    consider_rescue(regular_region_routes);
    consider_rescue(finegrain_routes);
    consider_rescue(smallnet_routes);
    consider_rescue(astar_routes);
    consider_rescue(rudy_routes);
    consider_rescue(rudy_classic_routes);
    consider_rescue(rudy_pin_hybrid_routes);
    consider_rescue(astar_early_routes);
    consider_rescue(detpart_routes);
    consider_rescue(nondet_routes);
    consider_rescue(rudy_partition_routes);
    consider_rescue(legacy_squeeze_routes);
    consider_rescue(direct_wl_routes);
    consider_rescue(ultra_direct_wl_routes);

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

  int trunk_polish_swaps = 0;
  int64_t trunk_polish_wl_gain = 0;
  int64_t trunk_polish_via_delta = 0;
  int trunk_polish_rank = 0;
  for (const OrderedNet& ordered_net : ordered_nets) {
    if (trunk_polish_rank >= top_trunk_nets) {
      break;
    }
    trunk_polish_rank++;

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

    const GRoute* best_polish_route = &route;
    RouteScore best_polish_score = current_score;
    int64_t best_polish_congestion_cost = current_congestion_cost;

    auto consider_polish = [&](const NetRouteMap& candidate_routes) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }

      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      if (candidate_score.wirelength >= best_polish_score.wirelength) {
        return;
      }
      const int64_t wl_drop_vs_current
          = std::max<int64_t>(0, current_score.wirelength - candidate_score.wirelength);
      if (wl_drop_vs_current == 0) {
        return;
      }

      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t polish_via_cap = policy.hard_via_guard * 4 + 140;
      if (candidate_extra_vias > polish_via_cap) {
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
      const int64_t congestion_delta = candidate_congestion_cost - current_congestion_cost;
      const int64_t allowed_congestion_delta
          = std::max<int64_t>(30, current_congestion_cost * 2 + 12);
      if (congestion_delta > allowed_congestion_delta
          && wl_drop_vs_current < std::max<int64_t>(1, tile_size / 22)) {
        return;
      }

      if (candidate_score.wirelength < best_polish_score.wirelength
          || (candidate_score.wirelength == best_polish_score.wirelength
              && candidate_congestion_cost < best_polish_congestion_cost)
          || (candidate_score.wirelength == best_polish_score.wirelength
              && candidate_congestion_cost == best_polish_congestion_cost
              && candidate_score.vias < best_polish_score.vias)) {
        best_polish_route = &candidate_it->second;
        best_polish_score = candidate_score;
        best_polish_congestion_cost = candidate_congestion_cost;
      }
    };

    consider_polish(balanced_routes);
    consider_polish(wirelength_routes);
    consider_polish(data_wirelength_routes);
    consider_polish(pin_density_routes);
    consider_polish(region_aware_routes);
    consider_polish(regular_region_routes);
    consider_polish(finegrain_routes);
    consider_polish(smallnet_routes);
    consider_polish(astar_routes);
    consider_polish(rudy_routes);
    consider_polish(rudy_classic_routes);
    consider_polish(rudy_pin_hybrid_routes);
    consider_polish(astar_early_routes);
    consider_polish(detpart_routes);
    consider_polish(nondet_routes);
    consider_polish(rudy_partition_routes);
    consider_polish(legacy_squeeze_routes);
    consider_polish(direct_wl_routes);
    consider_polish(ultra_direct_wl_routes);

    if (best_polish_route != &route) {
      trunk_polish_swaps++;
      trunk_polish_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_polish_score.wirelength);

      const int64_t current_extra_vias
          = std::max<int64_t>(0, current_score.vias - baseline_score.vias);
      const int64_t polish_extra_vias
          = std::max<int64_t>(0, best_polish_score.vias - baseline_score.vias);
      const int64_t extra_via_delta = polish_extra_vias - current_extra_vias;
      trunk_polish_via_delta += extra_via_delta;
      cumulative_extra_vias += extra_via_delta;
      cumulative_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_polish_score.wirelength);

      route = *best_polish_route;
    }

    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);
  }

  int wl_crush_swaps = 0;
  int64_t wl_crush_gain = 0;
  int64_t wl_crush_via_delta = 0;
  int wl_crush_rank = 0;
  const int wl_crush_nets = std::max<int>(
      1, static_cast<int>(ordered_nets.size()));
  for (const OrderedNet& ordered_net : ordered_nets) {
    if (wl_crush_rank >= wl_crush_nets) {
      break;
    }
    wl_crush_rank++;

    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      continue;
    }

    GRoute& route = route_it->second;
    const RouteScore baseline_score = ordered_net.baseline_score;
    const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
    const RouteScore current_score = scoreRoute(route);
    const int64_t current_extra_vias
        = std::max<int64_t>(0, current_score.vias - baseline_score.vias);

    removeRouteFromUsage(route, origin_x, origin_y, tile_size, selected_usage);
    const int64_t current_congestion_cost
        = routeCongestionPenalty(route,
                                 grid,
                                 origin_x,
                                 origin_y,
                                 tile_size,
                                 selected_usage,
                                 soft_capacities);

    const GRoute* best_crush_route = &route;
    RouteScore best_crush_score = current_score;

    auto consider_crush = [&](const NetRouteMap& candidate_routes) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }

      const RouteScore candidate_score = scoreRoute(candidate_it->second);
      if (candidate_score.wirelength >= best_crush_score.wirelength) {
        return;
      }

      const int64_t wl_drop_vs_best = std::max<int64_t>(
          0, best_crush_score.wirelength - candidate_score.wirelength);
      if (wl_drop_vs_best == 0) {
        return;
      }

      const int64_t candidate_extra_vias
          = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
      const int64_t crush_via_cap = policy.hard_via_guard * 6
                                    + (policy.long_net ? 280 : 180);
      if (candidate_extra_vias > crush_via_cap) {
        return;
      }

      const int64_t extra_vias_vs_best
          = std::max<int64_t>(0, candidate_score.vias - best_crush_score.vias);
      if (extra_vias_vs_best > 0
          && wl_drop_vs_best
                 < std::max<int64_t>(1, extra_vias_vs_best / 4)) {
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
          = policy.long_net ? std::max<int64_t>(42, current_congestion_cost * 2 + 16)
                            : (policy.medium_net
                                   ? std::max<int64_t>(28,
                                                       current_congestion_cost * 3 / 2 + 12)
                                   : std::max<int64_t>(16, current_congestion_cost + 8));
      const int64_t force_wl_gain
          = policy.long_net ? std::max<int64_t>(1, tile_size / 18)
                            : std::max<int64_t>(1, tile_size / 20);
      if (congestion_delta > allowed_congestion_delta
          && wl_drop_vs_best < force_wl_gain) {
        return;
      }

      best_crush_route = &candidate_it->second;
      best_crush_score = candidate_score;
    };

    consider_crush(balanced_routes);
    consider_crush(wirelength_routes);
    consider_crush(data_wirelength_routes);
    consider_crush(pin_density_routes);
    consider_crush(region_aware_routes);
    consider_crush(regular_region_routes);
    consider_crush(finegrain_routes);
    consider_crush(smallnet_routes);
    consider_crush(astar_routes);
    consider_crush(rudy_routes);
    consider_crush(rudy_classic_routes);
    consider_crush(rudy_pin_hybrid_routes);
    consider_crush(astar_early_routes);
    consider_crush(detpart_routes);
    consider_crush(nondet_routes);
    consider_crush(rudy_partition_routes);
    consider_crush(legacy_squeeze_routes);
    consider_crush(direct_wl_routes);
    consider_crush(ultra_direct_wl_routes);

    if (best_crush_route != &route) {
      wl_crush_swaps++;
      wl_crush_gain += std::max<int64_t>(
          0, current_score.wirelength - best_crush_score.wirelength);

      const int64_t crush_extra_vias
          = std::max<int64_t>(0, best_crush_score.vias - baseline_score.vias);
      const int64_t extra_via_delta = crush_extra_vias - current_extra_vias;
      wl_crush_via_delta += extra_via_delta;
      cumulative_extra_vias += extra_via_delta;
      cumulative_wl_gain += std::max<int64_t>(
          0, current_score.wirelength - best_crush_score.wirelength);

      route = *best_crush_route;
    }

    addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);
  }

  // Portfolio re-optimization pass:
  // Mix-and-match candidates from all routing modes and run a negotiated
  // selection loop that favors shorter trunks while progressively penalizing
  // persistent hotspot edges (FastRoute-style history + SPRoute soft-cap view).
  struct PortfolioCandidate
  {
    const GRoute* route{nullptr};
    RouteScore score;
    RouteSource source{RouteSource::kFastRoute};
    uint64_t fingerprint{0};
    int support{1};
  };
  using CandidateList = std::vector<PortfolioCandidate>;

  std::unordered_map<odb::dbNet*, GRoute> frozen_current_routes;
  frozen_current_routes.reserve(ordered_nets.size());
  std::unordered_map<odb::dbNet*, CandidateList> candidate_portfolios;
  candidate_portfolios.reserve(ordered_nets.size());

  auto build_portfolio = [&](const OrderedNet& ordered_net) {
    auto route_it = routes.find(ordered_net.db_net);
    if (route_it == routes.end()) {
      return;
    }

    const auto inserted
        = frozen_current_routes.emplace(ordered_net.db_net, route_it->second);
    const GRoute* frozen_route = &inserted.first->second;

    CandidateList candidates;
    candidates.reserve(24);
    std::unordered_map<uint64_t, size_t> fingerprint_to_index;
    fingerprint_to_index.reserve(32);

    auto append_candidate = [&](const GRoute* candidate_route, RouteSource source) {
      if (candidate_route == nullptr) {
        return;
      }
      const uint64_t fingerprint = routeFingerprint(*candidate_route);
      const auto seen_it = fingerprint_to_index.find(fingerprint);
      if (seen_it != fingerprint_to_index.end()) {
        candidates[seen_it->second].support++;
        return;
      }
      candidates.push_back(
          PortfolioCandidate{
              candidate_route, scoreRoute(*candidate_route), source, fingerprint, 1});
      fingerprint_to_index.emplace(fingerprint, candidates.size() - 1);
    };

    auto append_from_map = [&](const NetRouteMap& candidate_routes, RouteSource source) {
      const auto candidate_it = candidate_routes.find(ordered_net.db_net);
      if (candidate_it == candidate_routes.end()) {
        return;
      }
      append_candidate(&candidate_it->second, source);
    };

    append_candidate(frozen_route, RouteSource::kFastRoute);
    append_from_map(balanced_routes, RouteSource::kNewgrBalanced);
    append_from_map(wirelength_routes, RouteSource::kNewgrWirelength);
    append_from_map(data_wirelength_routes, RouteSource::kNewgrDataWirelength);
    append_from_map(pin_density_routes, RouteSource::kNewgrDataWirelength);
    append_from_map(region_aware_routes, RouteSource::kNewgrRegionAware);
    append_from_map(regular_region_routes, RouteSource::kNewgrRegularRegion);
    append_from_map(finegrain_routes, RouteSource::kNewgrFineGrain);
    append_from_map(smallnet_routes, RouteSource::kNewgrSmallNet);
    append_from_map(astar_routes, RouteSource::kNewgrAstar);
    append_from_map(rudy_routes, RouteSource::kNewgrRudy);
    append_from_map(rudy_classic_routes, RouteSource::kNewgrRudy);
    append_from_map(rudy_pin_hybrid_routes, RouteSource::kNewgrRudy);
    append_from_map(astar_early_routes, RouteSource::kNewgrAstarEarly);
    append_from_map(detpart_routes, RouteSource::kNewgrDetPartClassic);
    append_from_map(nondet_routes, RouteSource::kNewgrNonDetHybrid);
    append_from_map(rudy_partition_routes, RouteSource::kNewgrRudyPartition);
    append_from_map(legacy_squeeze_routes, RouteSource::kNewgrLegacySqueeze);
    append_from_map(direct_wl_routes, RouteSource::kNewgrDirectWirelength);
    append_from_map(ultra_direct_wl_routes, RouteSource::kNewgrUltraDirectWirelength);

    std::sort(candidates.begin(),
              candidates.end(),
              [](const PortfolioCandidate& lhs, const PortfolioCandidate& rhs) {
                if (lhs.score.wirelength != rhs.score.wirelength) {
                  return lhs.score.wirelength < rhs.score.wirelength;
                }
                if (lhs.score.vias != rhs.score.vias) {
                  return lhs.score.vias < rhs.score.vias;
                }
                if (lhs.score.bends != rhs.score.bends) {
                  return lhs.score.bends < rhs.score.bends;
                }
                if (lhs.support != rhs.support) {
                  return lhs.support > rhs.support;
                }
                return static_cast<int>(lhs.source) < static_cast<int>(rhs.source);
              });

    // Keep a compact portfolio for speed while preserving the shortest options.
    constexpr size_t kMaxPortfolioPerNet = 12;
    if (candidates.size() > kMaxPortfolioPerNet) {
      candidates.resize(kMaxPortfolioPerNet);
    }

    candidate_portfolios.emplace(ordered_net.db_net, std::move(candidates));
  };

  for (const OrderedNet& ordered_net : ordered_nets) {
    build_portfolio(ordered_net);
  }

  int negotiated_swaps = 0;
  int64_t negotiated_wl_gain = 0;
  int64_t negotiated_via_delta = 0;
  int64_t negotiated_final_overflow = totalSoftOverflow(
      selected_usage, grid, soft_capacities);
  EdgeHistoryMap edge_history;
  edge_history.reserve(std::max<size_t>(selected_usage.size(), 4096));
  double congestion_weight = 0.12;
  double history_weight = 0.08;
  constexpr int kNegotiationRounds = 3;

  for (int round = 0; round < kNegotiationRounds; round++) {
    int round_swaps = 0;
    int64_t round_wl_gain = 0;
    int64_t round_via_delta = 0;
    int round_rank = 0;

    for (const OrderedNet& ordered_net : ordered_nets) {
      auto route_it = routes.find(ordered_net.db_net);
      if (route_it == routes.end()) {
        continue;
      }
      const auto portfolio_it = candidate_portfolios.find(ordered_net.db_net);
      if (portfolio_it == candidate_portfolios.end()
          || portfolio_it->second.empty()) {
        continue;
      }

      const bool trunk_priority = round_rank < top_trunk_nets;
      round_rank++;

      GRoute& route = route_it->second;
      const RouteScore baseline_score = ordered_net.baseline_score;
      const SelectionPolicy policy = buildSelectionPolicy(baseline_score, tile_size);
      const RouteScore current_score = scoreRoute(route);

      removeRouteFromUsage(route, origin_x, origin_y, tile_size, selected_usage);

      const GRoute* best_route = &route;
      RouteScore best_score = current_score;
      RouteSource best_source = RouteSource::kFastRoute;
      int best_support = 1;
      const int64_t current_extra_vias
          = std::max<int64_t>(0, current_score.vias - baseline_score.vias);
      const int64_t best_congestion_cost_initial
          = routeCongestionPenalty(route,
                                   grid,
                                   origin_x,
                                   origin_y,
                                   tile_size,
                                   selected_usage,
                                   soft_capacities);
      int64_t best_congestion_cost = best_congestion_cost_initial;
      const double best_history_cost_initial = routeHistoryPenalty(route,
                                                                   origin_x,
                                                                   origin_y,
                                                                   tile_size,
                                                                   edge_history);
      double best_history_cost = best_history_cost_initial;
      const double via_weight = trunk_priority ? 0.08 : (policy.long_net ? 0.12 : 0.16);
      const double bend_weight = trunk_priority ? 0.02 : 0.04;

      auto objective = [&](const RouteScore& candidate_score,
                           int64_t candidate_congestion_cost,
                           double candidate_history_cost,
                           int candidate_support) {
        const int64_t candidate_extra_vias
            = std::max<int64_t>(0, candidate_score.vias - baseline_score.vias);
        const double support_bonus = trunk_priority ? 2.0 : 1.2;
        return static_cast<double>(candidate_extra_vias) * via_weight
               + static_cast<double>(candidate_score.bends) * bend_weight
               + static_cast<double>(candidate_congestion_cost) * congestion_weight
               + candidate_history_cost * history_weight
               - static_cast<double>(candidate_support) * support_bonus;
      };

      double best_objective = objective(
          current_score, best_congestion_cost, best_history_cost, best_support);
      for (const PortfolioCandidate& candidate : portfolio_it->second) {
        if (candidate.route == nullptr) {
          continue;
        }

        const int64_t candidate_extra_vias
            = std::max<int64_t>(0, candidate.score.vias - baseline_score.vias);
        const int64_t via_cap = policy.hard_via_guard * 7
                                + (trunk_priority ? 320 : (policy.long_net ? 240 : 160));
        if (candidate_extra_vias > via_cap) {
          continue;
        }
        if (candidate.score.wirelength > best_score.wirelength) {
          continue;
        }

        const int64_t wl_gain_vs_current
            = std::max<int64_t>(0, current_score.wirelength - candidate.score.wirelength);
        const int64_t extra_vias_vs_current
            = std::max<int64_t>(0, candidate.score.vias - current_score.vias);
        if (extra_vias_vs_current > 0
            && wl_gain_vs_current
                   < std::max<int64_t>(1, extra_vias_vs_current / 3)) {
          continue;
        }

        const int64_t candidate_congestion_cost
            = routeCongestionPenalty(*candidate.route,
                                     grid,
                                     origin_x,
                                     origin_y,
                                     tile_size,
                                     selected_usage,
                                     soft_capacities);
        const double candidate_history_cost = routeHistoryPenalty(*candidate.route,
                                                                  origin_x,
                                                                  origin_y,
                                                                  tile_size,
                                                                  edge_history);

        const int64_t wl_gain_vs_best = std::max<int64_t>(
            0, best_score.wirelength - candidate.score.wirelength);
        const int64_t congestion_delta_vs_best
            = candidate_congestion_cost - best_congestion_cost;
        if (wl_gain_vs_best > 0) {
          const int64_t allowed_congestion_delta
              = trunk_priority
                    ? std::max<int64_t>(48,
                                        best_congestion_cost * 2 + wl_gain_vs_best / 2)
                : policy.long_net
                    ? std::max<int64_t>(28,
                                        best_congestion_cost + wl_gain_vs_best / 3)
                    : std::max<int64_t>(16,
                                        best_congestion_cost / 2 + wl_gain_vs_best / 4);
          const int64_t force_wl_gain = trunk_priority
                                            ? std::max<int64_t>(1, tile_size / 22)
                                        : policy.long_net
                                            ? std::max<int64_t>(1, tile_size / 18)
                                            : std::max<int64_t>(1, tile_size / 14);
          if (congestion_delta_vs_best > allowed_congestion_delta
              && wl_gain_vs_best < force_wl_gain) {
            continue;
          }
        } else {
          const int64_t allowed_tie_congestion_delta
              = std::max<int64_t>(10, best_congestion_cost / 3);
          if (candidate_congestion_cost
                  > best_congestion_cost + allowed_tie_congestion_delta
              && candidate.support <= best_support) {
            continue;
          }
        }

        const double candidate_objective = objective(candidate.score,
                                                     candidate_congestion_cost,
                                                     candidate_history_cost,
                                                     candidate.support);

        bool take_candidate = false;
        if (candidate.score.wirelength < best_score.wirelength) {
          take_candidate = true;
        } else if (candidate.support != best_support) {
          take_candidate = candidate.support > best_support;
        } else if (std::abs(candidate_objective - best_objective) > 1e-9) {
          take_candidate = candidate_objective < best_objective;
        } else if (candidate.score.vias != best_score.vias) {
          take_candidate = candidate.score.vias < best_score.vias;
        } else if (candidate_congestion_cost != best_congestion_cost) {
          take_candidate = candidate_congestion_cost < best_congestion_cost;
        } else if (candidate.score.bends != best_score.bends) {
          take_candidate = candidate.score.bends < best_score.bends;
        }

        if (take_candidate) {
          best_objective = candidate_objective;
          best_route = candidate.route;
          best_score = candidate.score;
          best_source = candidate.source;
          best_support = candidate.support;
          best_congestion_cost = candidate_congestion_cost;
          best_history_cost = candidate_history_cost;
        }
      }

      if (best_route != &route) {
        route = *best_route;
        round_swaps++;
        round_wl_gain += std::max<int64_t>(
            0, current_score.wirelength - best_score.wirelength);
        const int64_t refined_extra_vias
            = std::max<int64_t>(0, best_score.vias - baseline_score.vias);
        round_via_delta += refined_extra_vias - current_extra_vias;
      }

      addRouteToUsage(route, origin_x, origin_y, tile_size, selected_usage);

      switch (best_source) {
        case RouteSource::kFastRoute:
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
        case RouteSource::kNewgrLegacySqueeze:
          selected_from_legacy_squeeze++;
          break;
        case RouteSource::kNewgrDirectWirelength:
          selected_from_direct_wl++;
          break;
        case RouteSource::kNewgrUltraDirectWirelength:
          selected_from_ultra_direct_wl++;
          break;
      }
    }

    int64_t round_overflow = 0;
    for (const auto& [edge, edge_usage] : selected_usage) {
      int soft_capacity = edgeHardCapacity(edge, grid);
      const auto soft_it = soft_capacities.find(edge);
      if (soft_it != soft_capacities.end()) {
        soft_capacity = soft_it->second;
      } else {
        soft_capacity = softCapacityFromDemand(edge, soft_capacity, 0);
      }

      const int overflow = std::max(0, edge_usage - soft_capacity);
      if (overflow > 0) {
        round_overflow += overflow;
        edge_history[edge] += static_cast<double>(overflow) * 0.70;
      } else {
        const auto history_it = edge_history.find(edge);
        if (history_it != edge_history.end()) {
          history_it->second *= 0.94;
        }
      }
    }

    negotiated_swaps += round_swaps;
    negotiated_wl_gain += round_wl_gain;
    negotiated_via_delta += round_via_delta;
    cumulative_wl_gain += round_wl_gain;
    cumulative_extra_vias += round_via_delta;
    negotiated_final_overflow = round_overflow;

    if (round_overflow > 0) {
      congestion_weight = std::min(0.72, congestion_weight * 1.12 + 0.02);
      history_weight = std::min(0.56, history_weight * 1.10 + 0.015);
    } else {
      congestion_weight = std::max(0.08, congestion_weight * 0.95);
      history_weight = std::max(0.05, history_weight * 0.93);
    }

    if (round_swaps == 0) {
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
  for (const auto& [db_net, route] : data_wirelength_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_data_wl++;
    }
  }
  for (const auto& [db_net, route] : pin_density_routes) {
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
  for (const auto& [db_net, route] : rudy_classic_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_rudy++;
    }
  }
  for (const auto& [db_net, route] : rudy_pin_hybrid_routes) {
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
  for (const auto& [db_net, route] : legacy_squeeze_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_legacy_squeeze++;
    }
  }
  for (const auto& [db_net, route] : direct_wl_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_direct_wl++;
    }
  }
  for (const auto& [db_net, route] : ultra_direct_wl_routes) {
    if (routes.find(db_net) == routes.end()) {
      routes.emplace(db_net, route);
      inserted_from_ultra_direct_wl++;
    }
  }
  logger_->info(utl::GRT,
                6004,
                "NEWGR WL-priority hybrid selected balanced={} wl={} data={} "
                "region={} regular={} fine={} small={} astar={} rudy={} "
                "astarEarly={} detpart={} nondet={} rudyPart={} "
                "legacySqz={} directWl={} ultraDirectWl={} "
                "(kept FR={}; +balanced={} "
                "+wl={} +data={} +region={} +regular={} +fine={} +small={} "
                "+astar={} +rudy={} +astarEarly={} +detpart={} +nondet={} "
                "+rudyPart={} +legacySqz={} +directWl={} +ultraDirectWl={}). "
                "Refine swaps={} wl-gain={} via-delta={}. "
                "Rescue swaps={} wl-gain={} via-delta={}. "
                "Trunk polish swaps={} wl-gain={} via-delta={}. "
                "WL crush swaps={} wl-gain={} via-delta={}. "
                "Negotiated swaps={} wl-gain={} via-delta={} "
                "final-soft-overflow={}. "
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
                selected_from_legacy_squeeze,
                selected_from_direct_wl,
                selected_from_ultra_direct_wl,
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
                inserted_from_legacy_squeeze,
                inserted_from_direct_wl,
                inserted_from_ultra_direct_wl,
                second_pass_swaps,
                second_pass_wl_gain,
                second_pass_via_delta,
                third_pass_swaps,
                third_pass_wl_gain,
                third_pass_via_delta,
                trunk_polish_swaps,
                trunk_polish_wl_gain,
                trunk_polish_via_delta,
                wl_crush_swaps,
                wl_crush_gain,
                wl_crush_via_delta,
                negotiated_swaps,
                negotiated_wl_gain,
                negotiated_via_delta,
                negotiated_final_overflow,
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
