#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <limits>
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

struct RouteLayerUsage
{
  uint64_t low_layer_wl{0};
  int min_x{std::numeric_limits<int>::max()};
  int max_x{std::numeric_limits<int>::min()};
  bool has_segment{false};
};

RouteLayerUsage analyzeRouteLayerUsage(const GRoute& route)
{
  RouteLayerUsage usage;
  for (const GSegment& segment : route) {
    usage.has_segment = true;
    usage.min_x = std::min(usage.min_x, std::min(segment.init_x, segment.final_x));
    usage.max_x = std::max(usage.max_x, std::max(segment.init_x, segment.final_x));
    if (segment.isVia()) {
      continue;
    }
    const int low_layer = std::min(segment.init_layer, segment.final_layer);
    if (low_layer <= 4) {
      usage.low_layer_wl += static_cast<uint64_t>(segment.length());
    }
  }
  if (!usage.has_segment) {
    usage.min_x = 0;
    usage.max_x = 0;
  }
  return usage;
}

uint64_t computeLowLayerWirelength(const NetRouteMap& routes)
{
  uint64_t low_layer_wl = 0;
  for (const auto& [db_net, route] : routes) {
    (void) db_net;
    low_layer_wl += analyzeRouteLayerUsage(route).low_layer_wl;
  }
  return low_layer_wl;
}

uint64_t routeCostForPins(const RouteScore& score, int pin_count)
{
  // Prioritize wirelength while keeping via count under control.
  const uint64_t via_weight
      = (pin_count <= 4) ? 10 : ((pin_count <= 12) ? 12 : 16);
  return score.wirelength + score.vias * via_weight;
}

uint64_t minSwapWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 16;
  }
  if (pin_count <= 12) {
    return 64;
  }
  if (pin_count <= 24) {
    return 180;
  }
  return 360;
}

int64_t maxViaIncreaseForSwap(int pin_count)
{
  if (pin_count <= 4) {
    return 8;
  }
  if (pin_count <= 12) {
    return 12;
  }
  if (pin_count <= 24) {
    return 18;
  }
  return 24;
}

int64_t maxViaDropForNewgrSwap(int pin_count)
{
  if (pin_count <= 4) {
    return 6;
  }
  if (pin_count <= 12) {
    return 12;
  }
  if (pin_count <= 24) {
    return 20;
  }
  return 28;
}

uint64_t minNewgrGraftWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 24;
  }
  if (pin_count <= 12) {
    return 96;
  }
  if (pin_count <= 24) {
    return 260;
  }
  return 700;
}

uint64_t minHighFanoutNewgrWirelengthGain(int pin_count)
{
  if (pin_count <= 32) {
    return 2600;
  }
  if (pin_count <= 48) {
    return 5200;
  }
  return 8600;
}

bool shouldSwapToFastRoute(const RouteScore& newgr_score,
                           const RouteScore& fastroute_score,
                           int pin_count)
{
  if (fastroute_score.segments == 0) {
    return false;
  }
  if (newgr_score.segments == 0) {
    return true;
  }

  // Keep the largest high-fanout nets on NEWGR-core to avoid unstable
  // topology flips.
  if (pin_count > 36) {
    return false;
  }

  // Positive means FastRoute is shorter.
  const int64_t wl_gain = static_cast<int64_t>(newgr_score.wirelength)
                          - static_cast<int64_t>(fastroute_score.wirelength);
  const int64_t via_increase = static_cast<int64_t>(fastroute_score.vias)
                               - static_cast<int64_t>(newgr_score.vias);
  if (wl_gain <= 0) {
    return false;
  }

  const int64_t via_budget = maxViaIncreaseForSwap(pin_count);
  const int64_t min_wl_gain
      = static_cast<int64_t>(minSwapWirelengthGain(pin_count));

  // Primary gate: meaningful wirelength reduction with controlled via impact.
  if (wl_gain >= min_wl_gain && via_increase <= via_budget) {
    return true;
  }

  const uint64_t newgr_cost = routeCostForPins(newgr_score, pin_count);
  const uint64_t fastroute_cost = routeCostForPins(fastroute_score, pin_count);
  if (fastroute_cost + minSwapWirelengthGain(pin_count) < newgr_cost
      && via_increase <= via_budget * 2) {
    return true;
  }

  // Strong wirelength win override.
  if (wl_gain >= std::max<int64_t>(
                     min_wl_gain * 3,
                     static_cast<int64_t>(newgr_score.wirelength / 12))
      && via_increase <= via_budget * 3) {
    return true;
  }

  return false;
}

bool shouldSwapToNewgr(const RouteScore& fastroute_score,
                       const RouteScore& newgr_score,
                       int pin_count)
{
  if (newgr_score.segments == 0) {
    return false;
  }
  if (fastroute_score.segments == 0) {
    return true;
  }

  const int64_t wl_gain = static_cast<int64_t>(fastroute_score.wirelength)
                          - static_cast<int64_t>(newgr_score.wirelength);
  if (wl_gain <= 0) {
    return false;
  }
  const bool high_fanout = pin_count > 28;
  if (high_fanout && pin_count > 80) {
    return false;
  }

  const int64_t via_increase = static_cast<int64_t>(newgr_score.vias)
                               - static_cast<int64_t>(fastroute_score.vias);
  const int64_t via_drop = static_cast<int64_t>(fastroute_score.vias)
                           - static_cast<int64_t>(newgr_score.vias);
  if (via_increase > maxViaIncreaseForSwap(pin_count) / 2) {
    return false;
  }
  if (via_drop > maxViaDropForNewgrSwap(pin_count)) {
    return false;
  }

  const uint64_t min_gain = minNewgrGraftWirelengthGain(pin_count);
  if (static_cast<uint64_t>(wl_gain) < min_gain) {
    return false;
  }
  if (high_fanout) {
    const uint64_t high_fanout_min_gain
        = std::max<uint64_t>(min_gain * 2,
                             minHighFanoutNewgrWirelengthGain(pin_count));
    if (static_cast<uint64_t>(wl_gain) < high_fanout_min_gain) {
      return false;
    }
    if (via_increase > 3 || via_drop > 10) {
      return false;
    }
  }

  const uint64_t newgr_cost = routeCostForPins(newgr_score, pin_count);
  const uint64_t fastroute_cost = routeCostForPins(fastroute_score, pin_count);
  if (high_fanout) {
    return newgr_cost + min_gain < fastroute_cost;
  }
  return newgr_cost + min_gain / 2 < fastroute_cost;
}

struct NewgrBackboneStats
{
  uint64_t replaced_with_fastroute{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
};

NetRouteMap buildNewgrBackboneHybrid(const NetRouteMap& newgr_routes,
                                     const NetRouteMap& fastroute_routes,
                                     const std::map<odb::dbNet*, Net*>& db_net_map,
                                     NewgrBackboneStats& stats)
{
  // NEWGR-core backbone with FastRoute swap-in for clear WL wins.
  NetRouteMap hybrid_routes = newgr_routes;
  for (const auto& [db_net, fastroute_route] : fastroute_routes) {
    if (fastroute_route.empty()) {
      continue;
    }
    const auto net_it = db_net_map.find(db_net);
    if (net_it == db_net_map.end() || net_it->second == nullptr) {
      continue;
    }

    const int pin_count = std::max(1, net_it->second->getNumPins());
    auto hybrid_it = hybrid_routes.find(db_net);
    if (hybrid_it == hybrid_routes.end() || hybrid_it->second.empty()) {
      hybrid_routes[db_net] = fastroute_route;
      ++stats.added_missing_nets;
      continue;
    }

    const RouteScore newgr_score = computeRouteScore(hybrid_it->second);
    const RouteScore fastroute_score = computeRouteScore(fastroute_route);
    ++stats.considered_nets;
    if (shouldSwapToFastRoute(newgr_score, fastroute_score, pin_count)) {
      hybrid_it->second = fastroute_route;
      ++stats.replaced_with_fastroute;
    }
  }

  return hybrid_routes;
}

struct FastRouteBackboneStats
{
  uint64_t replaced_with_newgr{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_swap_limit{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildFastRouteBackboneHybrid(const NetRouteMap& fastroute_routes,
                                         const NetRouteMap& newgr_routes,
                                         const std::map<odb::dbNet*, Net*>& db_net_map,
                                         uint64_t total_fastroute_vias,
                                         uint64_t total_fastroute_wirelength,
                                         FastRouteBackboneStats& stats)
{
  NetRouteMap hybrid_routes = fastroute_routes;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* newgr_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_drop{0};
    int64_t via_increase{0};
    int64_t low_layer_delta{0};
    int center_x{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(newgr_routes.size());
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
    if (!shouldSwapToNewgr(fastroute_score, newgr_score, pin_count)) {
      continue;
    }

    const RouteLayerUsage fastroute_usage = analyzeRouteLayerUsage(hybrid_it->second);
    const RouteLayerUsage newgr_usage = analyzeRouteLayerUsage(newgr_route);
    const int64_t wl_gain = static_cast<int64_t>(fastroute_score.wirelength)
                            - static_cast<int64_t>(newgr_score.wirelength);
    const int64_t via_drop = static_cast<int64_t>(fastroute_score.vias)
                             - static_cast<int64_t>(newgr_score.vias);
    const int64_t via_increase = static_cast<int64_t>(newgr_score.vias)
                                 - static_cast<int64_t>(fastroute_score.vias);
    const int64_t low_layer_delta
        = static_cast<int64_t>(newgr_usage.low_layer_wl)
          - static_cast<int64_t>(fastroute_usage.low_layer_wl);
    const int center_x = (newgr_usage.min_x + newgr_usage.max_x) / 2;
    const int64_t via_bonus = std::max<int64_t>(0, via_drop) * 20;
    const int64_t via_penalty = std::max<int64_t>(0, via_increase) * 280;
    const int64_t hotspot_penalty = std::max<int64_t>(0, low_layer_delta) / 4;
    const int64_t hotspot_bonus = std::max<int64_t>(0, -low_layer_delta) / 8;
    const int64_t fanout_penalty = (pin_count > 24) ? pin_count * 3 : 0;
    const int64_t priority
        = wl_gain + via_bonus - via_penalty - hotspot_penalty + hotspot_bonus
          - fanout_penalty;
    candidates.push_back({db_net,
                          &newgr_route,
                          pin_count,
                          wl_gain,
                          via_drop,
                          via_increase,
                          low_layer_delta,
                          center_x,
                          priority});
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              const int64_t lhs_via_inc = std::max<int64_t>(0, lhs.via_increase);
              const int64_t rhs_via_inc = std::max<int64_t>(0, rhs.via_increase);
              const int64_t lhs_via_drop = std::max<int64_t>(0, lhs.via_drop);
              const int64_t rhs_via_drop = std::max<int64_t>(0, rhs.via_drop);
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     -std::max<int64_t>(0, lhs.low_layer_delta),
                                     lhs_via_drop,
                                     -lhs_via_inc)
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       -std::max<int64_t>(0, rhs.low_layer_delta),
                                       rhs_via_drop,
                                       -rhs_via_inc);
            });

  stats.candidate_pool_size = candidates.size();
  const size_t swap_limit = std::min<size_t>(
      1200, std::max<size_t>(96, fastroute_routes.size() / 75));
  const int64_t total_via_drop_budget
      = std::max<int64_t>(6400, static_cast<int64_t>(total_fastroute_vias / 16));
  const int64_t total_via_increase_budget
      = std::max<int64_t>(640, static_cast<int64_t>(total_fastroute_vias / 160));
  const uint64_t wl_gain_target = std::max<uint64_t>(
      7200000, static_cast<uint64_t>(total_fastroute_wirelength / 70));
  const size_t min_swaps_before_stop = std::max<size_t>(64, (swap_limit * 3) / 5);
  int64_t consumed_via_drop = 0;
  int64_t consumed_via_increase = 0;
  uint64_t consumed_wl_gain = 0;

  int min_center_x = 0;
  int max_center_x = 0;
  if (!candidates.empty()) {
    min_center_x = candidates.front().center_x;
    max_center_x = candidates.front().center_x;
    for (const Candidate& candidate : candidates) {
      min_center_x = std::min(min_center_x, candidate.center_x);
      max_center_x = std::max(max_center_x, candidate.center_x);
    }
  }
  const int bucket_count = std::min<int>(
      16, std::max<int>(6, static_cast<int>(candidates.size() / 80) + 6));
  std::vector<std::vector<size_t>> buckets(bucket_count);
  for (size_t idx = 0; idx < candidates.size(); ++idx) {
    int bucket = 0;
    if (max_center_x > min_center_x) {
      const int64_t numer = static_cast<int64_t>(candidates[idx].center_x - min_center_x)
                            * bucket_count;
      bucket = static_cast<int>(numer / (max_center_x - min_center_x + 1));
      if (bucket >= bucket_count) {
        bucket = bucket_count - 1;
      }
    }
    buckets[bucket].push_back(idx);
  }
  std::vector<size_t> cursor(bucket_count, 0);
  int start_bucket = 0;
  while (stats.replaced_with_newgr < swap_limit) {
    bool consumed_any = false;
    for (int shift = 0; shift < bucket_count; ++shift) {
      const int bucket = (start_bucket + shift) % bucket_count;
      if (cursor[bucket] >= buckets[bucket].size()) {
        continue;
      }
      consumed_any = true;
      const Candidate& candidate = candidates[buckets[bucket][cursor[bucket]++]];
      const int64_t via_drop = std::max<int64_t>(0, candidate.via_drop);
      const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
      const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
      if (low_layer_delta
          > std::max<int64_t>(500, candidate.wl_gain / 2 + candidate.pin_count * 18)) {
        ++stats.skipped_by_layer_guard;
        continue;
      }
      if (consumed_via_drop + via_drop > total_via_drop_budget
          || consumed_via_increase + via_increase > total_via_increase_budget) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      // Accept positive-via swaps only when WL improvement is very strong and
      // low-layer demand does not increase significantly.
      if (via_increase > 0
          && candidate.wl_gain
                 < via_increase * 300 + low_layer_delta / 3
                       + static_cast<int64_t>(120)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      hybrid_routes[candidate.db_net] = *candidate.newgr_route;
      consumed_via_drop += via_drop;
      consumed_via_increase += via_increase;
      consumed_wl_gain
          += static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
      ++stats.replaced_with_newgr;
      if (consumed_wl_gain >= wl_gain_target
          && stats.replaced_with_newgr >= min_swaps_before_stop) {
        break;
      }
    }
    if (consumed_wl_gain >= wl_gain_target
        && stats.replaced_with_newgr >= min_swaps_before_stop) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    start_bucket = (start_bucket + 1) % bucket_count;
  }
  if (stats.replaced_with_newgr >= swap_limit
      && stats.candidate_pool_size > stats.replaced_with_newgr
                                       + stats.skipped_by_via_guard
                                       + stats.skipped_by_layer_guard) {
    stats.skipped_by_swap_limit
        = stats.candidate_pool_size - stats.replaced_with_newgr
          - stats.skipped_by_via_guard - stats.skipped_by_layer_guard;
  }
  stats.consumed_wl_gain = consumed_wl_gain;

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

bool shouldPreferFastRouteBackboneHybrid(int fr_overflow,
                                         const RouteScore& fastroute_score,
                                         uint64_t fastroute_low_layer_wl,
                                         int hybrid_overflow,
                                         const RouteScore& hybrid_score,
                                         uint64_t hybrid_low_layer_wl)
{
  if (hybrid_overflow > fr_overflow) {
    return false;
  }
  if (hybrid_overflow < fr_overflow) {
    return true;
  }
  if (hybrid_score.routed_nets < fastroute_score.routed_nets) {
    return false;
  }
  if (hybrid_score.wirelength >= fastroute_score.wirelength) {
    return false;
  }

  const int64_t via_increase = static_cast<int64_t>(hybrid_score.vias)
                               - static_cast<int64_t>(fastroute_score.vias);
  const int64_t via_drop = static_cast<int64_t>(fastroute_score.vias)
                           - static_cast<int64_t>(hybrid_score.vias);
  const int64_t via_increase_budget = std::max<int64_t>(
      360, static_cast<int64_t>(fastroute_score.vias / 220));
  const int64_t via_drop_budget = std::max<int64_t>(
      4200, static_cast<int64_t>(fastroute_score.vias / 20));
  if (via_increase > via_increase_budget || via_drop > via_drop_budget) {
    return false;
  }

  const uint64_t wl_gain = fastroute_score.wirelength - hybrid_score.wirelength;
  const int64_t low_layer_delta = static_cast<int64_t>(hybrid_low_layer_wl)
                                  - static_cast<int64_t>(fastroute_low_layer_wl);
  const int64_t low_layer_budget = std::max<int64_t>(
      2400000, static_cast<int64_t>(fastroute_low_layer_wl / 16));
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta * 2 + static_cast<int64_t>(260000))) {
    return false;
  }
  const uint64_t min_gain = std::max<uint64_t>(
      140000, fastroute_score.wirelength / 4600);
  return wl_gain >= min_gain;
}

bool shouldPreferNewgrBackboneHybrid(int incumbent_overflow,
                                     const RouteScore& incumbent_score,
                                     uint64_t incumbent_low_layer_wl,
                                     int candidate_overflow,
                                     const RouteScore& candidate_score,
                                     uint64_t candidate_low_layer_wl)
{
  if (candidate_overflow > incumbent_overflow) {
    return false;
  }
  if (candidate_score.routed_nets < incumbent_score.routed_nets) {
    return false;
  }
  if (candidate_overflow < incumbent_overflow) {
    return true;
  }
  if (candidate_score.wirelength >= incumbent_score.wirelength) {
    return false;
  }

  const uint64_t wl_gain = incumbent_score.wirelength - candidate_score.wirelength;
  const int64_t via_increase = static_cast<int64_t>(candidate_score.vias)
                               - static_cast<int64_t>(incumbent_score.vias);
  const int64_t low_layer_delta = static_cast<int64_t>(candidate_low_layer_wl)
                                  - static_cast<int64_t>(incumbent_low_layer_wl);
  const int64_t via_increase_budget
      = std::max<int64_t>(520, static_cast<int64_t>(incumbent_score.vias / 140));
  const int64_t low_layer_budget = std::max<int64_t>(
      3000000, static_cast<int64_t>(incumbent_low_layer_wl / 10));
  const uint64_t min_wl_gain = std::max<uint64_t>(
      420000, incumbent_score.wirelength / 1200);
  if (wl_gain < min_wl_gain) {
    return false;
  }
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 1200 + static_cast<int64_t>(300000))) {
    return false;
  }
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta * 2 + static_cast<int64_t>(360000))) {
    return false;
  }
  return true;
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
    const uint64_t newgr_low_layer_wl = computeLowLayerWirelength(routes);
    const uint64_t fastroute_low_layer_wl = computeLowLayerWirelength(fastroute_routes);
    NewgrBackboneStats newgr_backbone_stats;
    NetRouteMap newgr_backbone_hybrid = buildNewgrBackboneHybrid(
        routes, fastroute_routes, grouter_->db_net_map_, newgr_backbone_stats);
    const RouteScore newgr_backbone_score
        = computeRouteScore(newgr_backbone_hybrid);
    const uint64_t newgr_backbone_low_layer_wl
        = computeLowLayerWirelength(newgr_backbone_hybrid);
    FastRouteBackboneStats fastroute_backbone_stats;
    NetRouteMap fastroute_backbone_hybrid
        = buildFastRouteBackboneHybrid(fastroute_routes,
                                       routes,
                                       grouter_->db_net_map_,
                                       fastroute_score.vias,
                                       fastroute_score.wirelength,
                                       fastroute_backbone_stats);
    const RouteScore fastroute_backbone_score
        = computeRouteScore(fastroute_backbone_hybrid);
    const uint64_t fastroute_backbone_low_layer_wl
        = computeLowLayerWirelength(fastroute_backbone_hybrid);

    logger_->info(utl::GRT,
                  6006,
                  "NEWGR candidate summary: NEWGR(ofl={}, wl={}, vias={}, low_wl={}, nets={}), "
                  "NEWGR_BB_HYBRID(wl={}, vias={}, low_wl={}, nets={}, fr_swap={}, add={}), "
                  "FR_BB_HYBRID(wl={}, vias={}, low_wl={}, nets={}, ng_swap={}, add={}, "
                  "cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "swap_limit_skip={}, wl_gain={}), "
                  "FastRoute(ofl={}, wl={}, vias={}, low_wl={}, nets={})",
                  last_total_overflow_,
                  newgr_score.wirelength,
                  newgr_score.vias,
                  newgr_low_layer_wl,
                  newgr_score.routed_nets,
                  newgr_backbone_score.wirelength,
                  newgr_backbone_score.vias,
                  newgr_backbone_low_layer_wl,
                  newgr_backbone_score.routed_nets,
                  newgr_backbone_stats.replaced_with_fastroute,
                  newgr_backbone_stats.added_missing_nets,
                  fastroute_backbone_score.wirelength,
                  fastroute_backbone_score.vias,
                  fastroute_backbone_low_layer_wl,
                  fastroute_backbone_score.routed_nets,
                  fastroute_backbone_stats.replaced_with_newgr,
                  fastroute_backbone_stats.added_missing_nets,
                  fastroute_backbone_stats.candidate_pool_size,
                  fastroute_backbone_stats.skipped_by_via_guard,
                  fastroute_backbone_stats.skipped_by_layer_guard,
                  fastroute_backbone_stats.skipped_by_swap_limit,
                  fastroute_backbone_stats.consumed_wl_gain,
                  fastroute_overflow,
                  fastroute_score.wirelength,
                  fastroute_score.vias,
                  fastroute_low_layer_wl,
                  fastroute_score.routed_nets);

    // Detailed-route QoR has been more stable when FastRoute is used as the
    // default backbone, and NEWGR/hybrid are only used as overflow fallback.
    RouteScore best_score = fastroute_score;
    int best_overflow = fastroute_overflow;
    uint64_t best_low_layer_wl = fastroute_low_layer_wl;
    const char* selected_label = "FastRoute";
    bool selected_hybrid = false;
    uint64_t selected_swapped_nets = 0;
    uint64_t selected_added_nets = 0;
    used_fastroute_last_run_ = true;
    routes = std::move(fastroute_routes);

    if (shouldPreferFastRouteBackboneHybrid(fastroute_overflow,
                                            fastroute_score,
                                            fastroute_low_layer_wl,
                                            last_total_overflow_,
                                            fastroute_backbone_score,
                                            fastroute_backbone_low_layer_wl)) {
      routes = std::move(fastroute_backbone_hybrid);
      best_score = fastroute_backbone_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = fastroute_backbone_low_layer_wl;
      selected_label = "FastRoute+NEWGR targeted-graft";
      selected_hybrid = true;
      selected_swapped_nets = fastroute_backbone_stats.replaced_with_newgr;
      selected_added_nets = fastroute_backbone_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferNewgrBackboneHybrid(best_overflow,
                                        best_score,
                                        best_low_layer_wl,
                                        last_total_overflow_,
                                        newgr_backbone_score,
                                        newgr_backbone_low_layer_wl)) {
      routes = std::move(newgr_backbone_hybrid);
      best_score = newgr_backbone_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = newgr_backbone_low_layer_wl;
      selected_label = "NEWGR+FastRoute net-graft";
      selected_hybrid = true;
      selected_swapped_nets = newgr_backbone_stats.replaced_with_fastroute;
      selected_added_nets = newgr_backbone_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    } else if (isBetterRoute(last_total_overflow_,
                             newgr_score,
                             best_overflow,
                             best_score)
               && last_total_overflow_ < best_overflow) {
      best_score = newgr_score;
      best_overflow = last_total_overflow_;
      selected_label = "NEWGR-core";
      used_fastroute_last_run_ = false;
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
                    selected_swapped_nets,
                    selected_added_nets);
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
