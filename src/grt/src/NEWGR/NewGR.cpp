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
  int min_y{std::numeric_limits<int>::max()};
  int max_y{std::numeric_limits<int>::min()};
  bool has_segment{false};
};

RouteLayerUsage analyzeRouteLayerUsage(const GRoute& route)
{
  RouteLayerUsage usage;
  for (const GSegment& segment : route) {
    usage.has_segment = true;
    usage.min_x = std::min(usage.min_x, std::min(segment.init_x, segment.final_x));
    usage.max_x = std::max(usage.max_x, std::max(segment.init_x, segment.final_x));
    usage.min_y = std::min(usage.min_y, std::min(segment.init_y, segment.final_y));
    usage.max_y = std::max(usage.max_y, std::max(segment.init_y, segment.final_y));
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
    usage.min_y = 0;
    usage.max_y = 0;
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

uint64_t envelopeRouteValue(const RouteScore& score,
                            const RouteLayerUsage& usage,
                            int pin_count)
{
  // Envelope pass is wirelength-dominant; vias and low-layer pressure are
  // secondary guards to avoid pathological detailed-routing outcomes.
  const uint64_t via_weight
      = (pin_count <= 8) ? 4 : ((pin_count <= 24) ? 5 : 7);
  const uint64_t low_layer_divisor
      = (pin_count <= 8) ? 44 : ((pin_count <= 24) ? 56 : 72);
  return score.wirelength + score.vias * via_weight
         + usage.low_layer_wl / low_layer_divisor;
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

uint64_t minMicroGraftWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 20;
  }
  if (pin_count <= 12) {
    return 36;
  }
  if (pin_count <= 24) {
    return 72;
  }
  return 140;
}

uint64_t minInterleavedWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 12;
  }
  if (pin_count <= 12) {
    return 48;
  }
  if (pin_count <= 24) {
    return 140;
  }
  return 280;
}

uint64_t minWirelengthSweepGain(int pin_count)
{
  if (pin_count <= 4) {
    return 8;
  }
  if (pin_count <= 12) {
    return 24;
  }
  if (pin_count <= 24) {
    return 64;
  }
  return 120;
}

uint64_t minClosureWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 6;
  }
  if (pin_count <= 12) {
    return 18;
  }
  if (pin_count <= 24) {
    return 44;
  }
  return 92;
}

int64_t maxInterleavedViaIncrease(int pin_count)
{
  if (pin_count <= 4) {
    return 2;
  }
  if (pin_count <= 12) {
    return 3;
  }
  if (pin_count <= 24) {
    return 5;
  }
  return 7;
}

uint64_t minRadicalRefineWirelengthGain(int pin_count)
{
  if (pin_count <= 4) {
    return 36;
  }
  if (pin_count <= 12) {
    return 120;
  }
  if (pin_count <= 24) {
    return 300;
  }
  return 720;
}

int64_t maxRadicalRefineViaIncrease(int pin_count)
{
  if (pin_count <= 4) {
    return 2;
  }
  if (pin_count <= 12) {
    return 3;
  }
  if (pin_count <= 24) {
    return 4;
  }
  return 6;
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
  uint64_t polish_swaps{0};
  uint64_t polish_wl_gain{0};
  uint64_t micro_swaps{0};
  uint64_t micro_wl_gain{0};
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
    int center_y{0};
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
    const int center_y = (newgr_usage.min_y + newgr_usage.max_y) / 2;
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
                          center_y,
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
      1800, std::max<size_t>(160, fastroute_routes.size() / 45));
  const size_t polish_swap_budget
      = std::min<size_t>(700, std::max<size_t>(64, swap_limit / 2));
  const int64_t total_via_drop_budget
      = std::max<int64_t>(7600, static_cast<int64_t>(total_fastroute_vias / 14));
  const int64_t total_via_increase_budget
      = std::max<int64_t>(560, static_cast<int64_t>(total_fastroute_vias / 180));
  const uint64_t wl_gain_target = std::max<uint64_t>(
      11200000, static_cast<uint64_t>(total_fastroute_wirelength / 48));
  const uint64_t polish_wl_gain_target = std::max<uint64_t>(
      2600000, static_cast<uint64_t>(total_fastroute_wirelength / 180));
  const size_t min_swaps_before_stop = std::max<size_t>(120, (swap_limit * 3) / 4);
  int64_t consumed_via_drop = 0;
  int64_t consumed_via_increase = 0;
  uint64_t consumed_wl_gain = 0;
  std::vector<bool> selected_candidate(candidates.size(), false);

  int min_center_x = 0;
  int max_center_x = 0;
  int min_center_y = 0;
  int max_center_y = 0;
  if (!candidates.empty()) {
    min_center_x = candidates.front().center_x;
    max_center_x = candidates.front().center_x;
    min_center_y = candidates.front().center_y;
    max_center_y = candidates.front().center_y;
    for (const Candidate& candidate : candidates) {
      min_center_x = std::min(min_center_x, candidate.center_x);
      max_center_x = std::max(max_center_x, candidate.center_x);
      min_center_y = std::min(min_center_y, candidate.center_y);
      max_center_y = std::max(max_center_y, candidate.center_y);
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
      const size_t candidate_idx = buckets[bucket][cursor[bucket]++];
      const Candidate& candidate = candidates[candidate_idx];
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
                 < via_increase * 420 + low_layer_delta / 2
                       + static_cast<int64_t>(220)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      hybrid_routes[candidate.db_net] = *candidate.newgr_route;
      selected_candidate[candidate_idx] = true;
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

  // SPRoute-style deterministic second pass: consume remaining candidates in
  // priority order, but only when NEWGR strictly improves via and low-layer
  // demand while still providing WL gain.
  uint64_t polish_wl_gain = 0;
  for (size_t idx = 0;
       idx < candidates.size()
       && stats.polish_swaps < polish_swap_budget
       && stats.replaced_with_newgr < swap_limit + polish_swap_budget;
       ++idx) {
    if (selected_candidate[idx]) {
      continue;
    }
    const Candidate& candidate = candidates[idx];
    const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
    const int64_t via_drop = std::max<int64_t>(0, candidate.via_drop);
    if (via_increase > 0 || candidate.low_layer_delta > 0) {
      continue;
    }
    const int64_t min_polish_gain
        = std::max<int64_t>(96, candidate.pin_count <= 12 ? 72 : 180);
    if (candidate.wl_gain < min_polish_gain) {
      continue;
    }
    if (consumed_via_drop + via_drop > total_via_drop_budget) {
      continue;
    }
    hybrid_routes[candidate.db_net] = *candidate.newgr_route;
    selected_candidate[idx] = true;
    consumed_via_drop += via_drop;
    const uint64_t wl_gain = static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
    consumed_wl_gain += wl_gain;
    polish_wl_gain += wl_gain;
    ++stats.replaced_with_newgr;
    ++stats.polish_swaps;
    if (polish_wl_gain >= polish_wl_gain_target && stats.polish_swaps >= 24) {
      break;
    }
  }

  // Deterministic micro-graft pass:
  // consume smaller WL wins while keeping via non-increasing and bounding
  // low-layer growth. A 2D checkerboard schedule spreads updates spatially.
  const size_t micro_swap_budget = std::min<size_t>(
      520, std::max<size_t>(56, candidates.size() / 3));
  const uint64_t micro_wl_gain_target = std::max<uint64_t>(
      1300000, static_cast<uint64_t>(total_fastroute_wirelength / 420));
  const int micro_x_bucket_count = std::min<int>(
      12, std::max<int>(4, static_cast<int>(bucket_count / 2) + 1));
  const int micro_y_bucket_count = std::min<int>(
      12, std::max<int>(4, static_cast<int>(candidates.size() / 140) + 4));
  std::vector<std::vector<size_t>> micro_buckets(
      micro_x_bucket_count * micro_y_bucket_count);
  for (size_t idx = 0; idx < candidates.size(); ++idx) {
    if (selected_candidate[idx]) {
      continue;
    }
    int bucket_x = 0;
    if (max_center_x > min_center_x) {
      const int64_t numer
          = static_cast<int64_t>(candidates[idx].center_x - min_center_x)
            * micro_x_bucket_count;
      bucket_x = static_cast<int>(numer / (max_center_x - min_center_x + 1));
      if (bucket_x >= micro_x_bucket_count) {
        bucket_x = micro_x_bucket_count - 1;
      }
    }
    int bucket_y = 0;
    if (max_center_y > min_center_y) {
      const int64_t numer
          = static_cast<int64_t>(candidates[idx].center_y - min_center_y)
            * micro_y_bucket_count;
      bucket_y = static_cast<int>(numer / (max_center_y - min_center_y + 1));
      if (bucket_y >= micro_y_bucket_count) {
        bucket_y = micro_y_bucket_count - 1;
      }
    }
    const int bucket_id = bucket_y * micro_x_bucket_count + bucket_x;
    micro_buckets[bucket_id].push_back(idx);
  }

  std::vector<size_t> micro_cursor(micro_buckets.size(), 0);
  uint64_t micro_wl_gain = 0;
  int parity_seed = 0;
  while (stats.micro_swaps < micro_swap_budget
         && stats.replaced_with_newgr < swap_limit + polish_swap_budget + micro_swap_budget) {
    bool consumed_any = false;
    for (int phase = 0; phase < 2; ++phase) {
      const int parity = (parity_seed + phase) & 1;
      bool phase_consumed = false;
      for (int y = 0; y < micro_y_bucket_count; ++y) {
        for (int x = 0; x < micro_x_bucket_count; ++x) {
          if (((x + y) & 1) != parity) {
            continue;
          }
          const int bucket_id = y * micro_x_bucket_count + x;
          auto& bucket = micro_buckets[bucket_id];
          while (micro_cursor[bucket_id] < bucket.size()
                 && selected_candidate[bucket[micro_cursor[bucket_id]]]) {
            ++micro_cursor[bucket_id];
          }
          if (micro_cursor[bucket_id] >= bucket.size()) {
            continue;
          }
          const size_t candidate_idx = bucket[micro_cursor[bucket_id]++];
          const Candidate& candidate = candidates[candidate_idx];
          const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
          if (via_increase > 0) {
            continue;
          }
          const int64_t via_drop = std::max<int64_t>(0, candidate.via_drop);
          const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
          if (low_layer_delta > 160) {
            continue;
          }
          if (candidate.wl_gain
              < static_cast<int64_t>(minMicroGraftWirelengthGain(candidate.pin_count))) {
            continue;
          }
          if (low_layer_delta > 0
              && candidate.wl_gain
                     < low_layer_delta * 2 + static_cast<int64_t>(100)) {
            continue;
          }
          if (consumed_via_drop + via_drop > total_via_drop_budget) {
            continue;
          }

          hybrid_routes[candidate.db_net] = *candidate.newgr_route;
          selected_candidate[candidate_idx] = true;
          consumed_via_drop += via_drop;
          const uint64_t wl_gain
              = static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
          consumed_wl_gain += wl_gain;
          micro_wl_gain += wl_gain;
          ++stats.replaced_with_newgr;
          ++stats.micro_swaps;
          consumed_any = true;
          phase_consumed = true;
          if (stats.micro_swaps >= micro_swap_budget
              || (micro_wl_gain >= micro_wl_gain_target && stats.micro_swaps >= 28)) {
            break;
          }
        }
        if (stats.micro_swaps >= micro_swap_budget
            || (micro_wl_gain >= micro_wl_gain_target && stats.micro_swaps >= 28)) {
          break;
        }
      }
      if (!phase_consumed) {
        continue;
      }
      if (stats.micro_swaps >= micro_swap_budget
          || (micro_wl_gain >= micro_wl_gain_target && stats.micro_swaps >= 28)) {
        break;
      }
    }
    if (stats.micro_swaps >= micro_swap_budget
        || (micro_wl_gain >= micro_wl_gain_target && stats.micro_swaps >= 28)) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    parity_seed = (parity_seed + 1) & 1;
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
  stats.polish_wl_gain = polish_wl_gain;
  stats.micro_wl_gain = micro_wl_gain;

  return hybrid_routes;
}

struct InterleavedHybridStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildInterleavedBackboneHybrid(
    const NetRouteMap& base_routes,
    const NetRouteMap& donor_routes,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t base_total_vias,
    uint64_t base_low_layer_wl,
    InterleavedHybridStats& stats)
{
  NetRouteMap hybrid_routes = base_routes;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_drop{0};
    int64_t via_increase{0};
    int64_t low_layer_delta{0};
    int center_x{0};
    int center_y{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(donor_routes.size());
  for (const auto& [db_net, donor_route] : donor_routes) {
    if (donor_route.empty()) {
      continue;
    }
    const auto net_it = db_net_map.find(db_net);
    if (net_it == db_net_map.end() || net_it->second == nullptr) {
      continue;
    }

    const int pin_count = std::max(1, net_it->second->getNumPins());
    auto hybrid_it = hybrid_routes.find(db_net);
    if (hybrid_it == hybrid_routes.end() || hybrid_it->second.empty()) {
      hybrid_routes[db_net] = donor_route;
      ++stats.added_missing_nets;
      continue;
    }

    const RouteScore base_score = computeRouteScore(hybrid_it->second);
    const RouteScore donor_score = computeRouteScore(donor_route);
    ++stats.considered_nets;
    if (donor_score.segments == 0 || base_score.segments == 0) {
      continue;
    }

    const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                            - static_cast<int64_t>(donor_score.wirelength);
    if (wl_gain <= 0) {
      continue;
    }
    const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                 - static_cast<int64_t>(base_score.vias);
    const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                             - static_cast<int64_t>(donor_score.vias);
    const int64_t via_increase_limit
        = maxInterleavedViaIncrease(pin_count)
          + ((pin_count <= 20) ? 1 : ((pin_count <= 40) ? 2 : 3));
    if (via_increase > via_increase_limit) {
      ++stats.skipped_by_via_guard;
      continue;
    }

    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(hybrid_it->second);
    const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(donor_route);
    const int64_t low_layer_delta
        = static_cast<int64_t>(donor_usage.low_layer_wl)
          - static_cast<int64_t>(base_usage.low_layer_wl);
    if (low_layer_delta
        > std::max<int64_t>(320, wl_gain + static_cast<int64_t>(pin_count * 10))) {
      ++stats.skipped_by_layer_guard;
      continue;
    }
    const uint64_t min_gain = minInterleavedWirelengthGain(pin_count);
    if (static_cast<uint64_t>(wl_gain) < min_gain
        && !(via_drop > 0 && low_layer_delta <= 0)) {
      continue;
    }

    const int64_t via_bonus = std::max<int64_t>(0, via_drop) * 18;
    const int64_t via_penalty = std::max<int64_t>(0, via_increase) * 320;
    const int64_t low_layer_penalty = std::max<int64_t>(0, low_layer_delta) / 3;
    const int64_t low_layer_bonus = std::max<int64_t>(0, -low_layer_delta) / 6;
    const int64_t priority
        = wl_gain + via_bonus - via_penalty - low_layer_penalty + low_layer_bonus;
    candidates.push_back({db_net,
                          &donor_route,
                          pin_count,
                          wl_gain,
                          via_drop,
                          via_increase,
                          low_layer_delta,
                          (donor_usage.min_x + donor_usage.max_x) / 2,
                          (donor_usage.min_y + donor_usage.max_y) / 2,
                          priority});
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     -std::max<int64_t>(0, lhs.low_layer_delta),
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       -std::max<int64_t>(0, rhs.low_layer_delta),
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase));
            });

  stats.candidate_pool_size = candidates.size();
  const size_t swap_limit = std::min<size_t>(
      840, std::max<size_t>(120, base_routes.size() / 24));
  const int64_t total_via_increase_budget = std::max<int64_t>(
      420, static_cast<int64_t>(base_total_vias / 220));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      2600000, static_cast<int64_t>(base_low_layer_wl / 14));
  const uint64_t base_wirelength = computeRouteScore(base_routes).wirelength;
  const uint64_t wl_gain_target = std::max<uint64_t>(
      2200000, base_wirelength / 250);
  const size_t min_swaps_before_stop = std::max<size_t>(42, swap_limit / 4);

  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;

  int min_center_x = 0;
  int max_center_x = 0;
  int min_center_y = 0;
  int max_center_y = 0;
  if (!candidates.empty()) {
    min_center_x = candidates.front().center_x;
    max_center_x = candidates.front().center_x;
    min_center_y = candidates.front().center_y;
    max_center_y = candidates.front().center_y;
    for (const Candidate& candidate : candidates) {
      min_center_x = std::min(min_center_x, candidate.center_x);
      max_center_x = std::max(max_center_x, candidate.center_x);
      min_center_y = std::min(min_center_y, candidate.center_y);
      max_center_y = std::max(max_center_y, candidate.center_y);
    }
  }
  const int x_bucket_count = std::min<int>(
      18, std::max<int>(6, static_cast<int>(candidates.size() / 90) + 6));
  const int y_bucket_count = std::min<int>(
      16, std::max<int>(5, static_cast<int>(candidates.size() / 120) + 5));
  std::vector<std::vector<size_t>> buckets(x_bucket_count * y_bucket_count);
  for (size_t idx = 0; idx < candidates.size(); ++idx) {
    int bucket_x = 0;
    if (max_center_x > min_center_x) {
      const int64_t numer = static_cast<int64_t>(candidates[idx].center_x - min_center_x)
                            * x_bucket_count;
      bucket_x = static_cast<int>(numer / (max_center_x - min_center_x + 1));
      if (bucket_x >= x_bucket_count) {
        bucket_x = x_bucket_count - 1;
      }
    }
    int bucket_y = 0;
    if (max_center_y > min_center_y) {
      const int64_t numer = static_cast<int64_t>(candidates[idx].center_y - min_center_y)
                            * y_bucket_count;
      bucket_y = static_cast<int>(numer / (max_center_y - min_center_y + 1));
      if (bucket_y >= y_bucket_count) {
        bucket_y = y_bucket_count - 1;
      }
    }
    const int bucket_id = bucket_y * x_bucket_count + bucket_x;
    buckets[bucket_id].push_back(idx);
  }

  std::vector<size_t> cursor(buckets.size(), 0);
  std::vector<bool> selected_candidate(candidates.size(), false);
  int parity_seed = 0;
  while (stats.replaced_with_donor < swap_limit) {
    bool consumed_any = false;
    for (int phase = 0; phase < 2; ++phase) {
      const int parity = (parity_seed + phase) & 1;
      bool phase_consumed = false;
      for (int y = 0; y < y_bucket_count; ++y) {
        for (int x = 0; x < x_bucket_count; ++x) {
          if (((x + y) & 1) != parity) {
            continue;
          }
          const int bucket_id = y * x_bucket_count + x;
          auto& bucket = buckets[bucket_id];
          while (cursor[bucket_id] < bucket.size()) {
            const size_t candidate_idx = bucket[cursor[bucket_id]++];
            if (selected_candidate[candidate_idx]) {
              continue;
            }
            const Candidate& candidate = candidates[candidate_idx];
            const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
            const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
            if (consumed_via_increase + via_increase > total_via_increase_budget
                || consumed_low_layer_growth + low_layer_delta > low_layer_growth_budget) {
              ++stats.skipped_by_budget_guard;
              continue;
            }
            if (via_increase > 0
                && candidate.wl_gain
                       < via_increase * 180 + low_layer_delta / 5
                             + static_cast<int64_t>(80)) {
              ++stats.skipped_by_via_guard;
              continue;
            }
            if (low_layer_delta > 0
                && candidate.wl_gain
                       < low_layer_delta + static_cast<int64_t>(64)) {
              ++stats.skipped_by_layer_guard;
              continue;
            }
            if (candidate.pin_count > 48 && low_layer_delta > 0
                && candidate.wl_gain
                       < low_layer_delta * 6 / 5 + static_cast<int64_t>(180)) {
              ++stats.skipped_by_layer_guard;
              continue;
            }

            hybrid_routes[candidate.db_net] = *candidate.donor_route;
            selected_candidate[candidate_idx] = true;
            consumed_via_increase += via_increase;
            consumed_low_layer_growth += low_layer_delta;
            consumed_wl_gain += static_cast<uint64_t>(
                std::max<int64_t>(0, candidate.wl_gain));
            ++stats.replaced_with_donor;
            consumed_any = true;
            phase_consumed = true;
            break;
          }
          if (stats.replaced_with_donor >= swap_limit
              || (consumed_wl_gain >= wl_gain_target
                  && stats.replaced_with_donor >= min_swaps_before_stop)) {
            break;
          }
        }
        if (stats.replaced_with_donor >= swap_limit
            || (consumed_wl_gain >= wl_gain_target
                && stats.replaced_with_donor >= min_swaps_before_stop)) {
          break;
        }
      }
      if (!phase_consumed) {
        continue;
      }
      if (stats.replaced_with_donor >= swap_limit
          || (consumed_wl_gain >= wl_gain_target
              && stats.replaced_with_donor >= min_swaps_before_stop)) {
        break;
      }
    }
    if (stats.replaced_with_donor >= swap_limit
        || (consumed_wl_gain >= wl_gain_target
            && stats.replaced_with_donor >= min_swaps_before_stop)) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    parity_seed = (parity_seed + 1) & 1;
  }

  // FastRoute-style wirelength closure after spatial interleaving:
  // consume residual donor routes with small, bounded via/layer risk.
  const size_t closure_swap_budget = std::min<size_t>(
      560, std::max<size_t>(64, candidates.size() / 3));
  const int64_t closure_via_budget = std::max<int64_t>(
      320, static_cast<int64_t>(base_total_vias / 420));
  const int64_t closure_low_layer_budget = 950000;
  const uint64_t closure_wl_gain_target = std::max<uint64_t>(
      1400000, base_wirelength / 320);
  uint64_t closure_wl_gain = 0;
  size_t closure_swaps = 0;

  for (size_t idx = 0;
       idx < candidates.size()
       && closure_swaps < closure_swap_budget
       && stats.replaced_with_donor < swap_limit + closure_swap_budget;
       ++idx) {
    if (selected_candidate[idx]) {
      continue;
    }
    const Candidate& candidate = candidates[idx];
    const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
    const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
    const int64_t min_closure_gain = std::max<int64_t>(
        48, static_cast<int64_t>(minInterleavedWirelengthGain(candidate.pin_count) / 2));

    if (candidate.wl_gain < min_closure_gain) {
      continue;
    }
    if (via_increase > 3) {
      continue;
    }
    if (low_layer_delta > 0
        && candidate.wl_gain
               < low_layer_delta * 8 / 5 + static_cast<int64_t>(120)) {
      continue;
    }
    if (consumed_via_increase + via_increase
        > total_via_increase_budget + closure_via_budget) {
      continue;
    }
    if (consumed_low_layer_growth + low_layer_delta
        > low_layer_growth_budget + closure_low_layer_budget) {
      continue;
    }

    hybrid_routes[candidate.db_net] = *candidate.donor_route;
    selected_candidate[idx] = true;
    consumed_via_increase += via_increase;
    consumed_low_layer_growth += low_layer_delta;
    const uint64_t wl_gain
        = static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
    consumed_wl_gain += wl_gain;
    closure_wl_gain += wl_gain;
    ++stats.replaced_with_donor;
    ++closure_swaps;
    if (closure_wl_gain >= closure_wl_gain_target && closure_swaps >= 24) {
      break;
    }
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct WirelengthSweepStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildWirelengthSweepHybrid(
    const NetRouteMap& base_routes,
    const NetRouteMap& fastroute_routes,
    const NetRouteMap& newgr_routes,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t base_total_vias,
    uint64_t base_low_layer_wl,
    WirelengthSweepStats& stats)
{
  NetRouteMap hybrid_routes = base_routes;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_increase{0};
    int64_t via_drop{0};
    int64_t low_layer_delta{0};
    int center_x{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(db_net_map.size());

  for (const auto& [db_net, net] : db_net_map) {
    if (db_net == nullptr || net == nullptr) {
      continue;
    }
    const int pin_count = std::max(1, net->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    auto fr_it = fastroute_routes.find(db_net);
    auto ng_it = newgr_routes.find(db_net);
    const GRoute* fr_route = (fr_it == fastroute_routes.end()) ? nullptr : &fr_it->second;
    const GRoute* ng_route = (ng_it == newgr_routes.end()) ? nullptr : &ng_it->second;

    if ((base_it == hybrid_routes.end() || base_it->second.empty())
        && ((fr_route == nullptr || fr_route->empty())
            && (ng_route == nullptr || ng_route->empty()))) {
      continue;
    }

    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      const GRoute* donor = nullptr;
      RouteScore donor_score;
      if (fr_route != nullptr && !fr_route->empty()) {
        donor = fr_route;
        donor_score = computeRouteScore(*fr_route);
      }
      if (ng_route != nullptr && !ng_route->empty()) {
        const RouteScore ng_score = computeRouteScore(*ng_route);
        if (donor == nullptr
            || std::make_tuple(ng_score.wirelength, ng_score.vias)
                   < std::make_tuple(donor_score.wirelength, donor_score.vias)) {
          donor = ng_route;
          donor_score = ng_score;
        }
      }
      if (donor != nullptr) {
        hybrid_routes[db_net] = *donor;
        ++stats.added_missing_nets;
      }
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    if (base_score.segments == 0) {
      continue;
    }

    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);
    const GRoute* donor = nullptr;
    RouteScore donor_score;
    RouteLayerUsage donor_usage;

    auto consider_donor = [&](const GRoute* route) {
      if (route == nullptr || route->empty()) {
        return;
      }
      const RouteScore score = computeRouteScore(*route);
      if (score.segments == 0 || score.wirelength >= base_score.wirelength) {
        return;
      }
      const RouteLayerUsage usage = analyzeRouteLayerUsage(*route);
      if (donor == nullptr
          || std::make_tuple(score.wirelength, score.vias)
                 < std::make_tuple(donor_score.wirelength, donor_score.vias)) {
        donor = route;
        donor_score = score;
        donor_usage = usage;
      }
    };

    consider_donor(fr_route);
    consider_donor(ng_route);
    ++stats.considered_nets;
    if (donor == nullptr) {
      continue;
    }

    const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                            - static_cast<int64_t>(donor_score.wirelength);
    const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                 - static_cast<int64_t>(base_score.vias);
    const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                             - static_cast<int64_t>(donor_score.vias);
    const int64_t low_layer_delta
        = static_cast<int64_t>(donor_usage.low_layer_wl)
          - static_cast<int64_t>(base_usage.low_layer_wl);
    const int64_t via_increase_limit
        = maxInterleavedViaIncrease(pin_count)
          + ((pin_count <= 16) ? 3 : ((pin_count <= 40) ? 5 : 6));
    if (via_increase > via_increase_limit) {
      ++stats.skipped_by_via_guard;
      continue;
    }
    if (low_layer_delta
        > std::max<int64_t>(360, wl_gain + static_cast<int64_t>(pin_count * 11))) {
      ++stats.skipped_by_layer_guard;
      continue;
    }
    if (wl_gain < static_cast<int64_t>(minWirelengthSweepGain(pin_count))
        && !(via_drop > 0 && low_layer_delta <= 0)) {
      continue;
    }
    if (via_increase > 0
        && wl_gain
               < via_increase * 75 + low_layer_delta / 10
                     + static_cast<int64_t>(30)) {
      ++stats.skipped_by_via_guard;
      continue;
    }

    const int64_t via_bonus = std::max<int64_t>(0, via_drop) * 12;
    const int64_t via_penalty = std::max<int64_t>(0, via_increase) * 140;
    const int64_t low_layer_penalty = std::max<int64_t>(0, low_layer_delta) / 6;
    const int64_t low_layer_bonus = std::max<int64_t>(0, -low_layer_delta) / 10;
    const int64_t priority
        = wl_gain * 2 + via_bonus - via_penalty - low_layer_penalty
          + low_layer_bonus;
    candidates.push_back({db_net,
                          donor,
                          pin_count,
                          wl_gain,
                          via_increase,
                          via_drop,
                          low_layer_delta,
                          (donor_usage.min_x + donor_usage.max_x) / 2,
                          priority});
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase),
                                     -std::max<int64_t>(0, lhs.low_layer_delta))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase),
                                       -std::max<int64_t>(0, rhs.low_layer_delta));
            });

  stats.candidate_pool_size = candidates.size();
  const size_t swap_limit = std::min<size_t>(
      1320, std::max<size_t>(170, base_routes.size() / 15));
  const int64_t via_increase_budget = std::max<int64_t>(
      220, static_cast<int64_t>(base_total_vias / 560));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      2400000, static_cast<int64_t>(base_low_layer_wl / 16));
  const uint64_t base_wirelength = computeRouteScore(base_routes).wirelength;
  const uint64_t wl_gain_target = std::max<uint64_t>(
      3000000, base_wirelength / 200);
  const size_t min_swaps_before_stop = std::max<size_t>(72, swap_limit / 4);

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
      14, std::max<int>(5, static_cast<int>(candidates.size() / 120) + 5));
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
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;
  int start_bucket = 0;
  while (stats.replaced_with_donor < swap_limit) {
    bool consumed_any = false;
    for (int shift = 0; shift < bucket_count; ++shift) {
      const int bucket = (start_bucket + shift) % bucket_count;
      if (cursor[bucket] >= buckets[bucket].size()) {
        continue;
      }
      consumed_any = true;
      const Candidate& candidate = candidates[buckets[bucket][cursor[bucket]++]];
      const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
      const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
      if (via_increase > 0
          && candidate.wl_gain
                 < via_increase * 58 + low_layer_delta / 10
                       + static_cast<int64_t>(30)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      if (consumed_via_increase + via_increase > via_increase_budget
          || consumed_low_layer_growth + low_layer_delta > low_layer_growth_budget) {
        ++stats.skipped_by_budget_guard;
        continue;
      }
      hybrid_routes[candidate.db_net] = *candidate.donor_route;
      consumed_via_increase += via_increase;
      consumed_low_layer_growth += low_layer_delta;
      consumed_wl_gain += static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
      ++stats.replaced_with_donor;
      if (stats.replaced_with_donor >= swap_limit
          || (consumed_wl_gain >= wl_gain_target
              && stats.replaced_with_donor >= min_swaps_before_stop)) {
        break;
      }
    }
    if (stats.replaced_with_donor >= swap_limit
        || (consumed_wl_gain >= wl_gain_target
            && stats.replaced_with_donor >= min_swaps_before_stop)) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    start_bucket = (start_bucket + 1) % bucket_count;
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct RadicalRefineStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildRadicalRefineHybrid(
    const NetRouteMap& base_routes,
    const NetRouteMap& donor_routes,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t base_total_vias,
    uint64_t base_low_layer_wl,
    RadicalRefineStats& stats)
{
  NetRouteMap hybrid_routes = base_routes;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_increase{0};
    int64_t via_drop{0};
    int64_t low_layer_delta{0};
    int center_x{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(donor_routes.size());

  for (const auto& [db_net, donor_route] : donor_routes) {
    if (donor_route.empty()) {
      continue;
    }
    const auto net_it = db_net_map.find(db_net);
    if (net_it == db_net_map.end() || net_it->second == nullptr) {
      continue;
    }

    const int pin_count = std::max(1, net_it->second->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      hybrid_routes[db_net] = donor_route;
      ++stats.added_missing_nets;
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    const RouteScore donor_score = computeRouteScore(donor_route);
    ++stats.considered_nets;
    if (base_score.segments == 0 || donor_score.segments == 0
        || donor_score.wirelength >= base_score.wirelength) {
      continue;
    }

    const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                            - static_cast<int64_t>(donor_score.wirelength);
    const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                 - static_cast<int64_t>(base_score.vias);
    const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                             - static_cast<int64_t>(donor_score.vias);
    const int64_t via_increase_limit
        = maxRadicalRefineViaIncrease(pin_count)
          + ((pin_count <= 20) ? 1 : ((pin_count <= 40) ? 2 : 3));
    if (via_increase > via_increase_limit) {
      ++stats.skipped_by_via_guard;
      continue;
    }

    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);
    const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(donor_route);
    const int64_t low_layer_delta
        = static_cast<int64_t>(donor_usage.low_layer_wl)
          - static_cast<int64_t>(base_usage.low_layer_wl);
    if (low_layer_delta
        > std::max<int64_t>(
               820, wl_gain * 4 + static_cast<int64_t>(pin_count * 36))) {
      ++stats.skipped_by_layer_guard;
      continue;
    }
    const uint64_t min_gain = minRadicalRefineWirelengthGain(pin_count);
    if (wl_gain < static_cast<int64_t>(min_gain)
        && !(via_drop > 0 && via_increase <= 0 && low_layer_delta <= 0)) {
      continue;
    }
    if (via_increase > 0
        && wl_gain
               < via_increase * 150 + std::max<int64_t>(0, low_layer_delta) / 8
                     + static_cast<int64_t>(90)) {
      ++stats.skipped_by_via_guard;
      continue;
    }

    const int64_t via_bonus = std::max<int64_t>(0, via_drop) * 50;
    const int64_t via_penalty = std::max<int64_t>(0, via_increase) * 320;
    const int64_t low_layer_penalty = std::max<int64_t>(0, low_layer_delta) / 5;
    const int64_t low_layer_bonus = std::max<int64_t>(0, -low_layer_delta) / 12;
    const int64_t fanout_penalty = (pin_count > 32) ? pin_count * 4 : 0;
    const int64_t priority
        = wl_gain * 3 + via_bonus - via_penalty - low_layer_penalty
          + low_layer_bonus - fanout_penalty;
    candidates.push_back({db_net,
                          &donor_route,
                          pin_count,
                          wl_gain,
                          via_increase,
                          via_drop,
                          low_layer_delta,
                          (donor_usage.min_x + donor_usage.max_x) / 2,
                          priority});
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase),
                                     -std::max<int64_t>(0, lhs.low_layer_delta))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase),
                                       -std::max<int64_t>(0, rhs.low_layer_delta));
            });

  stats.candidate_pool_size = candidates.size();
  const size_t swap_limit = std::min<size_t>(
      1900, std::max<size_t>(220, base_routes.size() / 9));
  const int64_t via_increase_budget = std::max<int64_t>(
      360, static_cast<int64_t>(base_total_vias / 300));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      8600000, static_cast<int64_t>(base_low_layer_wl / 6));
  const uint64_t base_wirelength = computeRouteScore(base_routes).wirelength;
  const uint64_t wl_gain_target = std::max<uint64_t>(
      6400000, base_wirelength / 120);
  const size_t min_swaps_before_stop = std::max<size_t>(96, swap_limit / 5);

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
      18, std::max<int>(6, static_cast<int>(candidates.size() / 140) + 6));
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
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;
  int start_bucket = 0;
  while (stats.replaced_with_donor < swap_limit) {
    bool consumed_any = false;
    for (int shift = 0; shift < bucket_count; ++shift) {
      const int bucket = (start_bucket + shift) % bucket_count;
      if (cursor[bucket] >= buckets[bucket].size()) {
        continue;
      }
      consumed_any = true;
      const Candidate& candidate = candidates[buckets[bucket][cursor[bucket]++]];
      const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
      const int64_t low_layer_delta = std::max<int64_t>(0, candidate.low_layer_delta);
      if (via_increase > 0
          && candidate.wl_gain
                 < via_increase * 120 + low_layer_delta / 7
                       + static_cast<int64_t>(70)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      if (low_layer_delta > 0
          && candidate.wl_gain
                 < low_layer_delta / 7 + static_cast<int64_t>(120)) {
        ++stats.skipped_by_layer_guard;
        continue;
      }
      if (consumed_via_increase + via_increase > via_increase_budget
          || consumed_low_layer_growth + low_layer_delta > low_layer_growth_budget) {
        ++stats.skipped_by_budget_guard;
        continue;
      }
      hybrid_routes[candidate.db_net] = *candidate.donor_route;
      consumed_via_increase += via_increase;
      consumed_low_layer_growth += low_layer_delta;
      consumed_wl_gain += static_cast<uint64_t>(std::max<int64_t>(0, candidate.wl_gain));
      ++stats.replaced_with_donor;
      if (stats.replaced_with_donor >= swap_limit
          || (consumed_wl_gain >= wl_gain_target
              && stats.replaced_with_donor >= min_swaps_before_stop)) {
        break;
      }
    }
    if (stats.replaced_with_donor >= swap_limit
        || (consumed_wl_gain >= wl_gain_target
            && stats.replaced_with_donor >= min_swaps_before_stop)) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    start_bucket = (start_bucket + 1) % bucket_count;
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct EnvelopeHybridStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

const GRoute* findNetRoute(const NetRouteMap& routes, odb::dbNet* db_net)
{
  const auto it = routes.find(db_net);
  if (it == routes.end() || it->second.empty()) {
    return nullptr;
  }
  return &it->second;
}

NetRouteMap buildRadicalEnvelopeHybrid(
    const NetRouteMap& seed_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t seed_total_vias,
    uint64_t seed_low_layer_wl,
    EnvelopeHybridStats& stats)
{
  NetRouteMap hybrid_routes = seed_routes;
  const int64_t via_increase_budget = std::max<int64_t>(
      720, static_cast<int64_t>(seed_total_vias / 110));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      8800000, static_cast<int64_t>(seed_low_layer_wl / 5));
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;

  for (const auto& [db_net, net] : db_net_map) {
    if (net == nullptr) {
      continue;
    }
    const int pin_count = std::max(1, net->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      const GRoute* best_missing_route = nullptr;
      RouteScore best_missing_score;
      RouteLayerUsage best_missing_usage;
      uint64_t best_missing_cost = std::numeric_limits<uint64_t>::max();
      for (const NetRouteMap* donor_set : donor_route_sets) {
        if (donor_set == nullptr) {
          continue;
        }
        const GRoute* donor_route = findNetRoute(*donor_set, db_net);
        if (donor_route == nullptr) {
          continue;
        }
        ++stats.candidate_pool_size;
        const RouteScore donor_score = computeRouteScore(*donor_route);
        if (donor_score.segments == 0) {
          continue;
        }
        const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*donor_route);
        const uint64_t donor_cost
            = envelopeRouteValue(donor_score, donor_usage, pin_count);
        if (donor_cost < best_missing_cost
            || (donor_cost == best_missing_cost
                && donor_score.wirelength < best_missing_score.wirelength)) {
          best_missing_route = donor_route;
          best_missing_score = donor_score;
          best_missing_usage = donor_usage;
          best_missing_cost = donor_cost;
        }
      }
      if (best_missing_route != nullptr) {
        (void) best_missing_usage;
        hybrid_routes[db_net] = *best_missing_route;
        ++stats.added_missing_nets;
      }
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    if (base_score.segments == 0) {
      continue;
    }
    ++stats.considered_nets;
    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);

    const GRoute* best_route = nullptr;
    RouteScore best_score = base_score;
    RouteLayerUsage best_usage = base_usage;
    uint64_t best_value = envelopeRouteValue(base_score, base_usage, pin_count);

    for (const NetRouteMap* donor_set : donor_route_sets) {
      if (donor_set == nullptr) {
        continue;
      }
      const GRoute* donor_route = findNetRoute(*donor_set, db_net);
      if (donor_route == nullptr) {
        continue;
      }
      ++stats.candidate_pool_size;
      const RouteScore donor_score = computeRouteScore(*donor_route);
      if (donor_score.segments == 0) {
        continue;
      }

      const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                              - static_cast<int64_t>(donor_score.wirelength);
      if (wl_gain <= 0) {
        continue;
      }
      const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                   - static_cast<int64_t>(base_score.vias);
      const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*donor_route);
      const int64_t low_layer_delta
          = static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(base_usage.low_layer_wl);

      const int64_t via_limit
          = maxRadicalRefineViaIncrease(pin_count)
            + ((pin_count <= 20) ? 2 : ((pin_count <= 40) ? 3 : 4));
      if (via_increase > via_limit
          && wl_gain
                 < via_increase * 150 + std::max<int64_t>(0, low_layer_delta) / 8
                       + static_cast<int64_t>(90)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      if (low_layer_delta > 0
          && wl_gain
                 < low_layer_delta / 10 + static_cast<int64_t>(80)) {
        ++stats.skipped_by_layer_guard;
        continue;
      }

      const uint64_t donor_value
          = envelopeRouteValue(donor_score, donor_usage, pin_count);
      const int64_t wl_vs_best = static_cast<int64_t>(best_score.wirelength)
                                 - static_cast<int64_t>(donor_score.wirelength);
      const int64_t via_vs_best = std::max<int64_t>(
          0,
          static_cast<int64_t>(donor_score.vias)
              - static_cast<int64_t>(best_score.vias));
      const int64_t low_layer_vs_best = std::max<int64_t>(
          0,
          static_cast<int64_t>(donor_usage.low_layer_wl)
              - static_cast<int64_t>(best_usage.low_layer_wl));
      const int64_t via_tradeoff_penalty
          = via_vs_best * ((pin_count <= 16) ? 45 : ((pin_count <= 32) ? 60 : 80));
      const int64_t low_layer_tradeoff_penalty
          = low_layer_vs_best / ((pin_count <= 16) ? 22 : 30);
      const bool better_wirelength
          = donor_score.wirelength + minWirelengthSweepGain(pin_count)
            < best_score.wirelength;
      const bool better_cost
          = donor_value + minInterleavedWirelengthGain(pin_count) / 2 < best_value;
      const bool stronger_tradeoff
          = wl_vs_best
            > via_tradeoff_penalty + low_layer_tradeoff_penalty
                  + static_cast<int64_t>(20);
      if (better_wirelength || better_cost || stronger_tradeoff) {
        best_route = donor_route;
        best_score = donor_score;
        best_usage = donor_usage;
        best_value = donor_value;
      }
    }

    if (best_route == nullptr) {
      continue;
    }

    const int64_t final_via_increase
        = static_cast<int64_t>(best_score.vias)
          - static_cast<int64_t>(base_score.vias);
    const int64_t final_low_layer_delta
        = static_cast<int64_t>(best_usage.low_layer_wl)
          - static_cast<int64_t>(base_usage.low_layer_wl);
    const int64_t via_growth = std::max<int64_t>(0, final_via_increase);
    const int64_t low_layer_growth = std::max<int64_t>(0, final_low_layer_delta);
    if (consumed_via_increase + via_growth > via_increase_budget
        || consumed_low_layer_growth + low_layer_growth > low_layer_growth_budget) {
      ++stats.skipped_by_budget_guard;
      continue;
    }

    base_it->second = *best_route;
    consumed_via_increase += via_growth;
    consumed_low_layer_growth += low_layer_growth;
    consumed_wl_gain
        += static_cast<uint64_t>(
            std::max<int64_t>(0,
                              static_cast<int64_t>(base_score.wirelength)
                                  - static_cast<int64_t>(best_score.wirelength)));
    ++stats.replaced_with_donor;
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct WirelengthOracleStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildWirelengthOracleHybrid(
    const NetRouteMap& seed_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t seed_total_vias,
    uint64_t seed_low_layer_wl,
    WirelengthOracleStats& stats)
{
  NetRouteMap hybrid_routes = seed_routes;
  const uint64_t seed_total_wirelength = computeRouteScore(seed_routes).wirelength;
  // CUGR-style detailed-routability guard:
  // keep oracle swaps from consuming too much low-layer resource for small WL wins.
  const int64_t via_increase_budget = std::max<int64_t>(
      1200, static_cast<int64_t>(seed_total_vias / 95));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      1300000, static_cast<int64_t>(seed_low_layer_wl / 320));
  const uint64_t wl_gain_target = std::max<uint64_t>(
      3000000, seed_total_wirelength / 130);
  const size_t min_swaps_before_stop
      = std::max<size_t>(120, db_net_map.size() / 140);
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_increase{0};
    int64_t via_drop{0};
    int64_t low_layer_delta{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(db_net_map.size());

  for (const auto& [db_net, net] : db_net_map) {
    if (db_net == nullptr || net == nullptr) {
      continue;
    }
    const int pin_count = std::max(1, net->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      const GRoute* donor = nullptr;
      RouteScore donor_score;
      for (const NetRouteMap* donor_set : donor_route_sets) {
        if (donor_set == nullptr) {
          continue;
        }
        const GRoute* donor_route = findNetRoute(*donor_set, db_net);
        if (donor_route == nullptr) {
          continue;
        }
        ++stats.candidate_pool_size;
        const RouteScore candidate_score = computeRouteScore(*donor_route);
        if (candidate_score.segments == 0) {
          continue;
        }
        if (donor == nullptr
            || std::make_tuple(candidate_score.wirelength, candidate_score.vias)
                   < std::make_tuple(donor_score.wirelength, donor_score.vias)) {
          donor = donor_route;
          donor_score = candidate_score;
        }
      }
      if (donor != nullptr) {
        hybrid_routes[db_net] = *donor;
        ++stats.added_missing_nets;
      }
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    if (base_score.segments == 0) {
      continue;
    }
    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);
    ++stats.considered_nets;

    for (const NetRouteMap* donor_set : donor_route_sets) {
      if (donor_set == nullptr) {
        continue;
      }
      const GRoute* donor_route = findNetRoute(*donor_set, db_net);
      if (donor_route == nullptr) {
        continue;
      }
      ++stats.candidate_pool_size;
      const RouteScore donor_score = computeRouteScore(*donor_route);
      if (donor_score.segments == 0) {
        continue;
      }

      const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                              - static_cast<int64_t>(donor_score.wirelength);
      if (wl_gain <= 0) {
        continue;
      }

      const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                   - static_cast<int64_t>(base_score.vias);
      const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                               - static_cast<int64_t>(donor_score.vias);
      const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*donor_route);
      const int64_t low_layer_delta
          = static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(base_usage.low_layer_wl);

      const int64_t via_soft_limit
          = maxInterleavedViaIncrease(pin_count)
            + ((pin_count <= 12) ? 2 : ((pin_count <= 32) ? 4 : 5));
      if (via_increase > via_soft_limit
          && wl_gain
                 < via_increase * 150 + std::max<int64_t>(0, low_layer_delta) / 5
                       + static_cast<int64_t>(90)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      if (low_layer_delta > 0
          && wl_gain
                 < low_layer_delta / 3 + via_increase * 90
                       + static_cast<int64_t>(120)) {
        ++stats.skipped_by_layer_guard;
        continue;
      }
      const int64_t via_bonus = std::max<int64_t>(0, via_drop) * 14;
      const int64_t via_penalty = std::max<int64_t>(0, via_increase) * 180;
      const int64_t low_layer_penalty = std::max<int64_t>(0, low_layer_delta) / 3;
      const int64_t low_layer_bonus = std::max<int64_t>(0, -low_layer_delta) / 16;
      const int64_t priority
          = wl_gain * 4 + via_bonus - via_penalty - low_layer_penalty
            + low_layer_bonus;

      candidates.push_back({db_net,
                            donor_route,
                            pin_count,
                            wl_gain,
                            via_increase,
                            via_drop,
                            low_layer_delta,
                            priority});
    }
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase),
                                     -std::max<int64_t>(0, lhs.low_layer_delta))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase),
                                       -std::max<int64_t>(0, rhs.low_layer_delta));
            });

  const size_t swap_limit = std::min<size_t>(
      3600, std::max<size_t>(420, db_net_map.size() / 5));
  for (const Candidate& candidate : candidates) {
    if (stats.replaced_with_donor >= swap_limit) {
      break;
    }
    auto base_it = hybrid_routes.find(candidate.db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      continue;
    }
    const RouteScore live_base_score = computeRouteScore(base_it->second);
    if (live_base_score.segments == 0) {
      continue;
    }
    const RouteScore donor_score = computeRouteScore(*candidate.donor_route);
    if (donor_score.segments == 0 || donor_score.wirelength >= live_base_score.wirelength) {
      continue;
    }
    const RouteLayerUsage live_base_usage = analyzeRouteLayerUsage(base_it->second);
    const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*candidate.donor_route);
    const int64_t wl_gain = static_cast<int64_t>(live_base_score.wirelength)
                            - static_cast<int64_t>(donor_score.wirelength);
    const int64_t via_increase = std::max<int64_t>(
        0, static_cast<int64_t>(donor_score.vias) - static_cast<int64_t>(live_base_score.vias));
    const int64_t low_layer_growth = std::max<int64_t>(
        0,
        static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(live_base_usage.low_layer_wl));
    if (via_increase > 0
        && wl_gain
               < via_increase * 120 + low_layer_growth / 4
                     + static_cast<int64_t>(90)) {
      ++stats.skipped_by_via_guard;
      continue;
    }
    if (low_layer_growth > 0
        && wl_gain
               < low_layer_growth / 4 + static_cast<int64_t>(120)) {
      ++stats.skipped_by_layer_guard;
      continue;
    }
    if (consumed_via_increase + via_increase > via_increase_budget
        || consumed_low_layer_growth + low_layer_growth > low_layer_growth_budget) {
      ++stats.skipped_by_budget_guard;
      continue;
    }

    base_it->second = *candidate.donor_route;
    consumed_via_increase += via_increase;
    consumed_low_layer_growth += low_layer_growth;
    consumed_wl_gain += static_cast<uint64_t>(std::max<int64_t>(0, wl_gain));
    ++stats.replaced_with_donor;
    if (consumed_wl_gain >= wl_gain_target
        && stats.replaced_with_donor >= min_swaps_before_stop) {
      break;
    }
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct WirelengthClosureStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildWirelengthClosureHybrid(
    const NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t base_total_vias,
    uint64_t base_low_layer_wl,
    WirelengthClosureStats& stats)
{
  NetRouteMap hybrid_routes = base_routes;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_increase{0};
    int64_t via_drop{0};
    int64_t low_layer_delta{0};
    int center_x{0};
    int center_y{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(db_net_map.size());

  for (const auto& [db_net, net] : db_net_map) {
    if (db_net == nullptr || net == nullptr) {
      continue;
    }

    const int pin_count = std::max(1, net->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      const GRoute* donor = nullptr;
      RouteScore donor_score;
      for (const NetRouteMap* donor_set : donor_route_sets) {
        if (donor_set == nullptr) {
          continue;
        }
        const GRoute* donor_route = findNetRoute(*donor_set, db_net);
        if (donor_route == nullptr) {
          continue;
        }
        ++stats.candidate_pool_size;
        const RouteScore candidate_score = computeRouteScore(*donor_route);
        if (candidate_score.segments == 0) {
          continue;
        }
        if (donor == nullptr
            || std::make_tuple(candidate_score.wirelength, candidate_score.vias)
                   < std::make_tuple(donor_score.wirelength, donor_score.vias)) {
          donor = donor_route;
          donor_score = candidate_score;
        }
      }
      if (donor != nullptr) {
        hybrid_routes[db_net] = *donor;
        ++stats.added_missing_nets;
      }
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    if (base_score.segments == 0) {
      continue;
    }
    ++stats.considered_nets;
    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);

    const GRoute* best_donor = nullptr;
    RouteLayerUsage best_usage;
    int64_t best_priority = std::numeric_limits<int64_t>::min();
    int64_t best_wl_gain = 0;
    int64_t best_via_increase = 0;
    int64_t best_via_drop = 0;
    int64_t best_low_layer_delta = 0;

    for (const NetRouteMap* donor_set : donor_route_sets) {
      if (donor_set == nullptr) {
        continue;
      }
      const GRoute* donor_route = findNetRoute(*donor_set, db_net);
      if (donor_route == nullptr) {
        continue;
      }
      ++stats.candidate_pool_size;
      const RouteScore donor_score = computeRouteScore(*donor_route);
      if (donor_score.segments == 0) {
        continue;
      }

      const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                              - static_cast<int64_t>(donor_score.wirelength);
      if (wl_gain <= 0) {
        continue;
      }

      const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                   - static_cast<int64_t>(base_score.vias);
      const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                               - static_cast<int64_t>(donor_score.vias);
      const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*donor_route);
      const int64_t low_layer_delta
          = static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(base_usage.low_layer_wl);

      const int64_t via_increase_limit
          = (pin_count <= 8) ? 1 : ((pin_count <= 20) ? 2 : 4);
      if (via_increase > via_increase_limit) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      const bool via_safe_swap = via_increase <= 0;
      if (low_layer_delta
          > std::max<int64_t>(180, static_cast<int64_t>(pin_count * 8 + wl_gain / 3))
          && !(via_safe_swap
               && wl_gain
                      > std::max<int64_t>(
                          static_cast<int64_t>(120),
                          low_layer_delta / 2 + static_cast<int64_t>(32)))) {
        ++stats.skipped_by_layer_guard;
        continue;
      }
      if (wl_gain < static_cast<int64_t>(minClosureWirelengthGain(pin_count))
          && !(via_drop > 0 && low_layer_delta <= 0)) {
        continue;
      }
      if (via_increase > 0
          && wl_gain
                 < via_increase * 80 + std::max<int64_t>(0, low_layer_delta) / 4
                       + static_cast<int64_t>(40)) {
        ++stats.skipped_by_via_guard;
        continue;
      }

      const int64_t priority
          = wl_gain * 5 + std::max<int64_t>(0, via_drop) * 42
            - std::max<int64_t>(0, via_increase) * 180
            - std::max<int64_t>(0, low_layer_delta) / 3
            + std::max<int64_t>(0, -low_layer_delta) / 3;
      if (best_donor == nullptr
          || std::make_tuple(priority,
                             wl_gain,
                             std::max<int64_t>(0, via_drop),
                             -std::max<int64_t>(0, via_increase),
                             -std::max<int64_t>(0, low_layer_delta))
                 > std::make_tuple(best_priority,
                                   best_wl_gain,
                                   std::max<int64_t>(0, best_via_drop),
                                   -std::max<int64_t>(0, best_via_increase),
                                   -std::max<int64_t>(0, best_low_layer_delta))) {
        best_donor = donor_route;
        best_usage = donor_usage;
        best_priority = priority;
        best_wl_gain = wl_gain;
        best_via_increase = via_increase;
        best_via_drop = via_drop;
        best_low_layer_delta = low_layer_delta;
      }
    }

    if (best_donor == nullptr) {
      continue;
    }

    candidates.push_back({db_net,
                          best_donor,
                          pin_count,
                          best_wl_gain,
                          best_via_increase,
                          best_via_drop,
                          best_low_layer_delta,
                          (best_usage.min_x + best_usage.max_x) / 2,
                          (best_usage.min_y + best_usage.max_y) / 2,
                          best_priority});
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase),
                                     -std::max<int64_t>(0, lhs.low_layer_delta))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase),
                                       -std::max<int64_t>(0, rhs.low_layer_delta));
            });

  if (candidates.empty()) {
    return hybrid_routes;
  }

  const size_t swap_limit = std::min<size_t>(
      1600, std::max<size_t>(220, db_net_map.size() / 10));
  const int64_t via_increase_budget = std::max<int64_t>(
      420, static_cast<int64_t>(base_total_vias / 280));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      900000, static_cast<int64_t>(base_low_layer_wl / 520));
  const uint64_t base_wirelength = computeRouteScore(base_routes).wirelength;
  const uint64_t wl_gain_target = std::max<uint64_t>(
      2500000, base_wirelength / 220);
  const size_t min_swaps_before_stop = std::max<size_t>(56, swap_limit / 5);

  int min_center_x = candidates.front().center_x;
  int max_center_x = candidates.front().center_x;
  int min_center_y = candidates.front().center_y;
  int max_center_y = candidates.front().center_y;
  for (const Candidate& candidate : candidates) {
    min_center_x = std::min(min_center_x, candidate.center_x);
    max_center_x = std::max(max_center_x, candidate.center_x);
    min_center_y = std::min(min_center_y, candidate.center_y);
    max_center_y = std::max(max_center_y, candidate.center_y);
  }

  const int x_bucket_count = std::min<int>(
      16, std::max<int>(6, static_cast<int>(candidates.size() / 120) + 6));
  const int y_bucket_count = std::min<int>(
      14, std::max<int>(5, static_cast<int>(candidates.size() / 150) + 5));
  std::vector<std::vector<size_t>> buckets(x_bucket_count * y_bucket_count);
  for (size_t idx = 0; idx < candidates.size(); ++idx) {
    int bucket_x = 0;
    if (max_center_x > min_center_x) {
      const int64_t numer = static_cast<int64_t>(candidates[idx].center_x - min_center_x)
                            * x_bucket_count;
      bucket_x = static_cast<int>(numer / (max_center_x - min_center_x + 1));
      if (bucket_x >= x_bucket_count) {
        bucket_x = x_bucket_count - 1;
      }
    }
    int bucket_y = 0;
    if (max_center_y > min_center_y) {
      const int64_t numer = static_cast<int64_t>(candidates[idx].center_y - min_center_y)
                            * y_bucket_count;
      bucket_y = static_cast<int>(numer / (max_center_y - min_center_y + 1));
      if (bucket_y >= y_bucket_count) {
        bucket_y = y_bucket_count - 1;
      }
    }
    const int bucket_id = bucket_y * x_bucket_count + bucket_x;
    buckets[bucket_id].push_back(idx);
  }

  std::vector<size_t> cursor(buckets.size(), 0);
  std::vector<bool> selected(candidates.size(), false);
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;
  int parity_seed = 0;
  while (stats.replaced_with_donor < swap_limit) {
    bool consumed_any = false;
    for (int phase = 0; phase < 2; ++phase) {
      const int parity = (parity_seed + phase) & 1;
      bool phase_consumed = false;
      for (int y = 0; y < y_bucket_count; ++y) {
        for (int x = 0; x < x_bucket_count; ++x) {
          if (((x + y) & 1) != parity) {
            continue;
          }
          const int bucket_id = y * x_bucket_count + x;
          auto& bucket = buckets[bucket_id];
          while (cursor[bucket_id] < bucket.size()
                 && selected[bucket[cursor[bucket_id]]]) {
            ++cursor[bucket_id];
          }
          if (cursor[bucket_id] >= bucket.size()) {
            continue;
          }
          const size_t candidate_idx = bucket[cursor[bucket_id]++];
          const Candidate& candidate = candidates[candidate_idx];
          const int64_t via_increase = std::max<int64_t>(0, candidate.via_increase);
          const int64_t low_layer_growth = std::max<int64_t>(0, candidate.low_layer_delta);
          if (via_increase > 0
              && candidate.wl_gain
                     < via_increase * 90 + low_layer_growth * 2 / 5
                           + static_cast<int64_t>(50)) {
            ++stats.skipped_by_via_guard;
            continue;
          }
          if (low_layer_growth > 0) {
            const int64_t low_layer_threshold
                = (via_increase > 0)
                      ? (low_layer_growth + static_cast<int64_t>(40))
                      : (low_layer_growth * 3 / 4 + static_cast<int64_t>(24));
            if (candidate.wl_gain < low_layer_threshold) {
              ++stats.skipped_by_layer_guard;
              continue;
            }
          }
          if (consumed_via_increase + via_increase > via_increase_budget
              || consumed_low_layer_growth + low_layer_growth
                     > low_layer_growth_budget) {
            ++stats.skipped_by_budget_guard;
            continue;
          }

          hybrid_routes[candidate.db_net] = *candidate.donor_route;
          selected[candidate_idx] = true;
          consumed_via_increase += via_increase;
          consumed_low_layer_growth += low_layer_growth;
          consumed_wl_gain += static_cast<uint64_t>(
              std::max<int64_t>(0, candidate.wl_gain));
          ++stats.replaced_with_donor;
          consumed_any = true;
          phase_consumed = true;
          break;
        }
        if (stats.replaced_with_donor >= swap_limit
            || (consumed_wl_gain >= wl_gain_target
                && stats.replaced_with_donor >= min_swaps_before_stop)) {
          break;
        }
      }
      if (!phase_consumed) {
        continue;
      }
      if (stats.replaced_with_donor >= swap_limit
          || (consumed_wl_gain >= wl_gain_target
              && stats.replaced_with_donor >= min_swaps_before_stop)) {
        break;
      }
    }
    if (stats.replaced_with_donor >= swap_limit
        || (consumed_wl_gain >= wl_gain_target
            && stats.replaced_with_donor >= min_swaps_before_stop)) {
      break;
    }
    if (!consumed_any) {
      break;
    }
    parity_seed = (parity_seed + 1) & 1;
  }

  stats.consumed_wl_gain = consumed_wl_gain;
  return hybrid_routes;
}

struct WirelengthFusionStats
{
  uint64_t replaced_with_donor{0};
  uint64_t added_missing_nets{0};
  uint64_t considered_nets{0};
  uint64_t candidate_pool_size{0};
  uint64_t skipped_by_via_guard{0};
  uint64_t skipped_by_layer_guard{0};
  uint64_t skipped_by_budget_guard{0};
  uint64_t consumed_wl_gain{0};
};

NetRouteMap buildWirelengthFusionHybrid(
    const NetRouteMap& seed_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    const std::map<odb::dbNet*, Net*>& db_net_map,
    uint64_t seed_total_vias,
    uint64_t seed_low_layer_wl,
    WirelengthFusionStats& stats)
{
  NetRouteMap hybrid_routes = seed_routes;
  const uint64_t seed_total_wirelength = computeRouteScore(seed_routes).wirelength;
  // Aggressive wirelength-majority fusion:
  // combine all donor families while using soft-cap-like budgets to cap risk.
  const int64_t via_increase_budget = std::max<int64_t>(
      2200, static_cast<int64_t>(seed_total_vias / 55));
  const int64_t low_layer_growth_budget = std::max<int64_t>(
      2200000, static_cast<int64_t>(seed_low_layer_wl / 180));
  const uint64_t wl_gain_target = std::max<uint64_t>(
      4200000, seed_total_wirelength / 90);
  const size_t min_swaps_before_stop
      = std::max<size_t>(180, db_net_map.size() / 120);
  int64_t consumed_via_increase = 0;
  int64_t consumed_low_layer_growth = 0;
  uint64_t consumed_wl_gain = 0;

  struct Candidate
  {
    odb::dbNet* db_net{nullptr};
    const GRoute* donor_route{nullptr};
    int pin_count{0};
    int64_t wl_gain{0};
    int64_t via_increase{0};
    int64_t via_drop{0};
    int64_t low_layer_delta{0};
    int64_t priority{0};
  };

  std::vector<Candidate> candidates;
  candidates.reserve(db_net_map.size());

  for (const auto& [db_net, net] : db_net_map) {
    if (db_net == nullptr || net == nullptr) {
      continue;
    }
    const int pin_count = std::max(1, net->getNumPins());
    auto base_it = hybrid_routes.find(db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      const GRoute* donor = nullptr;
      RouteScore donor_score;
      for (const NetRouteMap* donor_set : donor_route_sets) {
        if (donor_set == nullptr) {
          continue;
        }
        const GRoute* donor_route = findNetRoute(*donor_set, db_net);
        if (donor_route == nullptr) {
          continue;
        }
        ++stats.candidate_pool_size;
        const RouteScore candidate_score = computeRouteScore(*donor_route);
        if (candidate_score.segments == 0) {
          continue;
        }
        if (donor == nullptr
            || std::make_tuple(candidate_score.wirelength, candidate_score.vias)
                   < std::make_tuple(donor_score.wirelength, donor_score.vias)) {
          donor = donor_route;
          donor_score = candidate_score;
        }
      }
      if (donor != nullptr) {
        hybrid_routes[db_net] = *donor;
        ++stats.added_missing_nets;
      }
      continue;
    }

    const RouteScore base_score = computeRouteScore(base_it->second);
    if (base_score.segments == 0) {
      continue;
    }
    const RouteLayerUsage base_usage = analyzeRouteLayerUsage(base_it->second);
    ++stats.considered_nets;

    for (const NetRouteMap* donor_set : donor_route_sets) {
      if (donor_set == nullptr) {
        continue;
      }
      const GRoute* donor_route = findNetRoute(*donor_set, db_net);
      if (donor_route == nullptr) {
        continue;
      }
      ++stats.candidate_pool_size;
      const RouteScore donor_score = computeRouteScore(*donor_route);
      if (donor_score.segments == 0) {
        continue;
      }

      const int64_t wl_gain = static_cast<int64_t>(base_score.wirelength)
                              - static_cast<int64_t>(donor_score.wirelength);
      if (wl_gain <= 0) {
        continue;
      }

      const int64_t via_increase = static_cast<int64_t>(donor_score.vias)
                                   - static_cast<int64_t>(base_score.vias);
      const int64_t via_drop = static_cast<int64_t>(base_score.vias)
                               - static_cast<int64_t>(donor_score.vias);
      const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*donor_route);
      const int64_t low_layer_delta
          = static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(base_usage.low_layer_wl);
      const int64_t low_layer_release = std::max<int64_t>(0, -low_layer_delta);
      const int64_t low_layer_growth = std::max<int64_t>(0, low_layer_delta);

      const int64_t via_soft_limit
          = maxInterleavedViaIncrease(pin_count)
            + ((pin_count <= 12) ? 3 : ((pin_count <= 32) ? 6 : 8));
      const int64_t effective_gain
          = wl_gain + low_layer_release / 3 - std::max<int64_t>(0, via_increase) * 90
            - low_layer_growth / 8;
      if (via_increase > via_soft_limit
          && effective_gain
                 < via_increase * 140 + static_cast<int64_t>(100)) {
        ++stats.skipped_by_via_guard;
        continue;
      }
      if (low_layer_growth > 0
          && wl_gain
                 < low_layer_growth / 3 + std::max<int64_t>(0, via_increase) * 70
                       + static_cast<int64_t>(50)) {
        ++stats.skipped_by_layer_guard;
        continue;
      }
      if (effective_gain <= 0) {
        continue;
      }

      const int64_t priority
          = effective_gain * 6 + std::max<int64_t>(0, via_drop) * 20
            + low_layer_release / 4 - low_layer_growth / 10;
      candidates.push_back({db_net,
                            donor_route,
                            pin_count,
                            wl_gain,
                            via_increase,
                            via_drop,
                            low_layer_delta,
                            priority});
    }
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              return std::make_tuple(lhs.priority,
                                     lhs.wl_gain,
                                     std::max<int64_t>(0, lhs.via_drop),
                                     -std::max<int64_t>(0, lhs.via_increase),
                                     -std::max<int64_t>(0, lhs.low_layer_delta))
                     > std::make_tuple(rhs.priority,
                                       rhs.wl_gain,
                                       std::max<int64_t>(0, rhs.via_drop),
                                       -std::max<int64_t>(0, rhs.via_increase),
                                       -std::max<int64_t>(0, rhs.low_layer_delta));
            });

  const size_t swap_limit = std::min<size_t>(
      4200, std::max<size_t>(520, db_net_map.size() / 4));
  for (const Candidate& candidate : candidates) {
    if (stats.replaced_with_donor >= swap_limit) {
      break;
    }
    auto base_it = hybrid_routes.find(candidate.db_net);
    if (base_it == hybrid_routes.end() || base_it->second.empty()) {
      continue;
    }
    const RouteScore live_base_score = computeRouteScore(base_it->second);
    if (live_base_score.segments == 0) {
      continue;
    }
    const RouteScore donor_score = computeRouteScore(*candidate.donor_route);
    if (donor_score.segments == 0 || donor_score.wirelength >= live_base_score.wirelength) {
      continue;
    }

    const RouteLayerUsage live_base_usage = analyzeRouteLayerUsage(base_it->second);
    const RouteLayerUsage donor_usage = analyzeRouteLayerUsage(*candidate.donor_route);
    const int64_t wl_gain = static_cast<int64_t>(live_base_score.wirelength)
                            - static_cast<int64_t>(donor_score.wirelength);
    const int64_t via_increase = std::max<int64_t>(
        0, static_cast<int64_t>(donor_score.vias) - static_cast<int64_t>(live_base_score.vias));
    const int64_t low_layer_growth = std::max<int64_t>(
        0,
        static_cast<int64_t>(donor_usage.low_layer_wl)
            - static_cast<int64_t>(live_base_usage.low_layer_wl));
    const int64_t low_layer_release = std::max<int64_t>(
        0,
        static_cast<int64_t>(live_base_usage.low_layer_wl)
            - static_cast<int64_t>(donor_usage.low_layer_wl));
    const int64_t effective_gain = wl_gain + low_layer_release / 3;

    if (via_increase > 0
        && effective_gain
               < via_increase * 120 + low_layer_growth / 4
                     + static_cast<int64_t>(80)) {
      ++stats.skipped_by_via_guard;
      continue;
    }
    if (low_layer_growth > 0
        && wl_gain
               < low_layer_growth / 3 + via_increase * 60
                     + static_cast<int64_t>(60)) {
      ++stats.skipped_by_layer_guard;
      continue;
    }
    if (consumed_via_increase + via_increase > via_increase_budget
        || consumed_low_layer_growth + low_layer_growth > low_layer_growth_budget) {
      ++stats.skipped_by_budget_guard;
      continue;
    }

    base_it->second = *candidate.donor_route;
    consumed_via_increase += via_increase;
    consumed_low_layer_growth += low_layer_growth;
    consumed_wl_gain += static_cast<uint64_t>(std::max<int64_t>(0, wl_gain));
    ++stats.replaced_with_donor;
    if (consumed_wl_gain >= wl_gain_target
        && stats.replaced_with_donor >= min_swaps_before_stop) {
      break;
    }
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
      = std::max<int64_t>(420, static_cast<int64_t>(incumbent_score.vias / 170));
  const int64_t low_layer_budget = std::max<int64_t>(
      1800000, static_cast<int64_t>(incumbent_low_layer_wl / 24));
  const uint64_t min_wl_gain = std::max<uint64_t>(
      520000, incumbent_score.wirelength / 900);
  if (wl_gain < min_wl_gain) {
    return false;
  }
  if (low_layer_delta > low_layer_budget) {
    return false;
  }
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 1200 + static_cast<int64_t>(300000))) {
    return false;
  }
  if (low_layer_delta > 0
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta / 2 + static_cast<int64_t>(520000))) {
    return false;
  }
  return true;
}

bool shouldPreferInterleavedHybrid(int incumbent_overflow,
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
  const uint64_t min_wl_gain
      = std::max<uint64_t>(90000, incumbent_score.wirelength / 14000);
  if (wl_gain < min_wl_gain) {
    return false;
  }
  const int64_t via_increase_budget
      = std::max<int64_t>(150, static_cast<int64_t>(incumbent_score.vias / 700));
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 700 + static_cast<int64_t>(160000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      900000, static_cast<int64_t>(incumbent_low_layer_wl / 35));
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta / 2 + static_cast<int64_t>(180000))) {
    return false;
  }
  return true;
}

bool shouldPreferWirelengthSweepHybrid(int incumbent_overflow,
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
  const uint64_t min_wl_gain
      = std::max<uint64_t>(70000, incumbent_score.wirelength / 18000);
  if (wl_gain < min_wl_gain) {
    return false;
  }
  const int64_t via_increase_budget
      = std::max<int64_t>(280, static_cast<int64_t>(incumbent_score.vias / 390));
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 520 + static_cast<int64_t>(120000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      1200000, static_cast<int64_t>(incumbent_low_layer_wl / 24));
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta / 3 + static_cast<int64_t>(180000))) {
    return false;
  }
  return true;
}

bool shouldPreferRadicalRefineHybrid(int incumbent_overflow,
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
  const int64_t via_drop = static_cast<int64_t>(incumbent_score.vias)
                           - static_cast<int64_t>(candidate_score.vias);
  const int64_t low_layer_delta = static_cast<int64_t>(candidate_low_layer_wl)
                                  - static_cast<int64_t>(incumbent_low_layer_wl);
  const uint64_t min_wl_gain
      = std::max<uint64_t>(120000, incumbent_score.wirelength / 9600);
  if (wl_gain < min_wl_gain) {
    return false;
  }

  const int64_t via_increase_budget
      = std::max<int64_t>(360, static_cast<int64_t>(incumbent_score.vias / 300));
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 460 + static_cast<int64_t>(180000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      2200000, static_cast<int64_t>(incumbent_low_layer_wl / 12));
  if (via_increase <= 0 && via_drop > 0) {
    if (low_layer_delta > low_layer_budget * 2
        && wl_gain
               < static_cast<uint64_t>(low_layer_delta / 4 + static_cast<int64_t>(240000))) {
      return false;
    }
    return true;
  }
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta / 4 + static_cast<int64_t>(220000))) {
    return false;
  }
  return true;
}

bool shouldPreferEnvelopeHybrid(int incumbent_overflow,
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
  const int64_t via_increase_budget = std::max<int64_t>(
      980, static_cast<int64_t>(incumbent_score.vias / 95));
  if (via_increase > via_increase_budget
      && wl_gain
             < static_cast<uint64_t>(via_increase * 240 + static_cast<int64_t>(140000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      3600000, static_cast<int64_t>(incumbent_low_layer_wl / 8));
  if (low_layer_delta > low_layer_budget
      && wl_gain
             < static_cast<uint64_t>(low_layer_delta / 5 + static_cast<int64_t>(180000))) {
    return false;
  }
  return true;
}

bool shouldPreferWirelengthOracleHybrid(int incumbent_overflow,
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
  // SPRoute/CUGR-style soft-cap guard:
  // allow modest via growth when candidate releases significant low-layer
  // demand (proxy for reducing pin-access/resource risk).
  const uint64_t low_layer_release_credit = static_cast<uint64_t>(
      std::max<int64_t>(0, -low_layer_delta) / 4);
  const uint64_t effective_wl_gain = wl_gain + low_layer_release_credit;
  const bool strong_risk_drop
      = low_layer_delta <= -80000
        && (via_increase <= 0 || via_increase <= 220);
  const uint64_t min_wl_gain
      = std::max<uint64_t>(22000, incumbent_score.wirelength / 32000);
  if (effective_wl_gain < min_wl_gain
      && !(strong_risk_drop && wl_gain >= min_wl_gain / 2)) {
    return false;
  }
  const int64_t via_increase_budget
      = std::max<int64_t>(520, static_cast<int64_t>(incumbent_score.vias / 180));
  if (via_increase > via_increase_budget
      && effective_wl_gain
             < static_cast<uint64_t>(via_increase * 300 + static_cast<int64_t>(50000))) {
    return false;
  }
  if (via_increase > 0
      && effective_wl_gain
             < std::max<uint64_t>(
                 static_cast<uint64_t>(via_increase * 260
                                       + static_cast<int64_t>(50000)),
                 static_cast<uint64_t>(80000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      1500000, static_cast<int64_t>(incumbent_low_layer_wl / 180));
  if (low_layer_delta > low_layer_budget
      || (low_layer_delta > 0
          && effective_wl_gain
                 < static_cast<uint64_t>(low_layer_delta / 3
                                         + static_cast<int64_t>(50000)))) {
    return false;
  }
  return true;
}

bool shouldPreferWirelengthClosureHybrid(int incumbent_overflow,
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
  // Promote routes that release low-layer resources while keeping via growth
  // bounded; this mirrors soft-cap balancing from SPRoute and CUGR.
  const uint64_t low_layer_release_credit = static_cast<uint64_t>(
      std::max<int64_t>(0, -low_layer_delta) / 4);
  const uint64_t effective_wl_gain = wl_gain + low_layer_release_credit;
  const bool strong_risk_drop
      = low_layer_delta <= -70000
        && (via_increase <= 0 || via_increase <= 220);
  const uint64_t min_wl_gain
      = std::max<uint64_t>(18000, incumbent_score.wirelength / 36000);
  if (effective_wl_gain < min_wl_gain
      && !(strong_risk_drop && wl_gain >= min_wl_gain / 2)) {
    return false;
  }
  const int64_t via_increase_budget
      = std::max<int64_t>(380, static_cast<int64_t>(incumbent_score.vias / 250));
  if (via_increase > via_increase_budget
      && effective_wl_gain
             < static_cast<uint64_t>(via_increase * 320 + static_cast<int64_t>(60000))) {
    return false;
  }
  if (via_increase > 0
      && effective_wl_gain
             < std::max<uint64_t>(
                 static_cast<uint64_t>(via_increase * 260
                                       + static_cast<int64_t>(60000)),
                 static_cast<uint64_t>(90000))) {
    return false;
  }
  const int64_t low_layer_budget = std::max<int64_t>(
      1200000, static_cast<int64_t>(incumbent_low_layer_wl / 150));
  if (low_layer_delta > low_layer_budget
      && effective_wl_gain
             < static_cast<uint64_t>(low_layer_delta / 2
                                     + static_cast<int64_t>(60000))) {
    return false;
  }
  if (low_layer_delta > 0 && via_increase > 0
      && effective_wl_gain
             < static_cast<uint64_t>(low_layer_delta / 2
                                     + via_increase * 320
                                     + static_cast<int64_t>(60000))) {
    return false;
  }
  return true;
}

bool shouldPreferWirelengthFusionHybrid(int incumbent_overflow,
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
  const int64_t via_growth = std::max<int64_t>(0, via_increase);
  const int64_t low_layer_delta = static_cast<int64_t>(candidate_low_layer_wl)
                                  - static_cast<int64_t>(incumbent_low_layer_wl);
  const uint64_t low_layer_release_credit = static_cast<uint64_t>(
      std::max<int64_t>(0, -low_layer_delta) / 3);
  const uint64_t effective_wl_gain = wl_gain + low_layer_release_credit;
  const uint64_t min_wl_gain
      = std::max<uint64_t>(16000, incumbent_score.wirelength / 50000);
  if (effective_wl_gain < min_wl_gain) {
    return false;
  }

  const int64_t via_increase_budget
      = std::max<int64_t>(520, static_cast<int64_t>(incumbent_score.vias / 220));
  if (via_growth > via_increase_budget
      && effective_wl_gain
             < static_cast<uint64_t>(via_growth * 240 + static_cast<int64_t>(60000))) {
    return false;
  }
  if (via_growth > 0
      && effective_wl_gain
             < std::max<uint64_t>(
                 static_cast<uint64_t>(
                     via_growth * 170
                     + std::max<int64_t>(0, low_layer_delta) / 4
                     + static_cast<int64_t>(40000)),
                 static_cast<uint64_t>(70000))) {
    return false;
  }

  const int64_t low_layer_budget = std::max<int64_t>(
      1300000, static_cast<int64_t>(incumbent_low_layer_wl / 130));
  if (low_layer_delta > low_layer_budget
      && effective_wl_gain
             < static_cast<uint64_t>(low_layer_delta / 3
                                     + via_growth * 180
                                     + static_cast<int64_t>(80000))) {
    return false;
  }
  if (via_growth > 0 && low_layer_delta > 0
      && effective_wl_gain
             < static_cast<uint64_t>(via_growth * 210
                                     + low_layer_delta / 3
                                     + static_cast<int64_t>(90000))) {
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
    InterleavedHybridStats interleaved_stats;
    NetRouteMap interleaved_hybrid = buildInterleavedBackboneHybrid(
        fastroute_backbone_hybrid,
        newgr_backbone_hybrid,
        grouter_->db_net_map_,
        fastroute_backbone_score.vias,
        fastroute_backbone_low_layer_wl,
        interleaved_stats);
    const RouteScore interleaved_score = computeRouteScore(interleaved_hybrid);
    const uint64_t interleaved_low_layer_wl
        = computeLowLayerWirelength(interleaved_hybrid);
    WirelengthSweepStats sweep_stats;
    NetRouteMap wirelength_sweep_hybrid = buildWirelengthSweepHybrid(
        interleaved_hybrid,
        fastroute_routes,
        routes,
        grouter_->db_net_map_,
        interleaved_score.vias,
        interleaved_low_layer_wl,
        sweep_stats);
    const RouteScore sweep_score = computeRouteScore(wirelength_sweep_hybrid);
    const uint64_t sweep_low_layer_wl
        = computeLowLayerWirelength(wirelength_sweep_hybrid);
    WirelengthSweepStats extreme_sweep_stats;
    NetRouteMap extreme_sweep_hybrid = buildWirelengthSweepHybrid(
        wirelength_sweep_hybrid,
        newgr_backbone_hybrid,
        routes,
        grouter_->db_net_map_,
        sweep_score.vias,
        sweep_low_layer_wl,
        extreme_sweep_stats);
    const RouteScore extreme_sweep_score
        = computeRouteScore(extreme_sweep_hybrid);
    const uint64_t extreme_sweep_low_layer_wl
        = computeLowLayerWirelength(extreme_sweep_hybrid);
    RadicalRefineStats radical_refine_stats;
    NetRouteMap radical_refine_hybrid = buildRadicalRefineHybrid(
        extreme_sweep_hybrid,
        newgr_backbone_hybrid,
        grouter_->db_net_map_,
        extreme_sweep_score.vias,
        extreme_sweep_low_layer_wl,
        radical_refine_stats);
    const RouteScore radical_refine_score
        = computeRouteScore(radical_refine_hybrid);
    const uint64_t radical_refine_low_layer_wl
        = computeLowLayerWirelength(radical_refine_hybrid);
    EnvelopeHybridStats envelope_stats;
    std::vector<const NetRouteMap*> envelope_donors{
        &fastroute_routes,
        &routes,
        &newgr_backbone_hybrid,
        &fastroute_backbone_hybrid,
        &interleaved_hybrid,
        &wirelength_sweep_hybrid,
        &extreme_sweep_hybrid,
        &radical_refine_hybrid};
    NetRouteMap envelope_hybrid = buildRadicalEnvelopeHybrid(
        radical_refine_hybrid,
        envelope_donors,
        grouter_->db_net_map_,
        radical_refine_score.vias,
        radical_refine_low_layer_wl,
        envelope_stats);
    const RouteScore envelope_score = computeRouteScore(envelope_hybrid);
    const uint64_t envelope_low_layer_wl
        = computeLowLayerWirelength(envelope_hybrid);
    InterleavedHybridStats envelope_newgr_stats;
    NetRouteMap envelope_newgr_hybrid = buildInterleavedBackboneHybrid(
        envelope_hybrid,
        newgr_backbone_hybrid,
        grouter_->db_net_map_,
        envelope_score.vias,
        envelope_low_layer_wl,
        envelope_newgr_stats);
    const RouteScore envelope_newgr_score
        = computeRouteScore(envelope_newgr_hybrid);
    const uint64_t envelope_newgr_low_layer_wl
        = computeLowLayerWirelength(envelope_newgr_hybrid);
    WirelengthOracleStats wirelength_oracle_stats;
    std::vector<const NetRouteMap*> wirelength_oracle_donors{
        &fastroute_routes,
        &routes,
        &newgr_backbone_hybrid,
        &fastroute_backbone_hybrid,
        &interleaved_hybrid,
        &wirelength_sweep_hybrid,
        &extreme_sweep_hybrid,
        &radical_refine_hybrid,
        &envelope_hybrid};
    NetRouteMap wirelength_oracle_hybrid = buildWirelengthOracleHybrid(
        envelope_hybrid,
        wirelength_oracle_donors,
        grouter_->db_net_map_,
        envelope_score.vias,
        envelope_low_layer_wl,
        wirelength_oracle_stats);
    const RouteScore wirelength_oracle_score
        = computeRouteScore(wirelength_oracle_hybrid);
    const uint64_t wirelength_oracle_low_layer_wl
        = computeLowLayerWirelength(wirelength_oracle_hybrid);
    WirelengthClosureStats closure_stats;
    std::vector<const NetRouteMap*> closure_donors{
        &fastroute_routes,
        &routes,
        &newgr_backbone_hybrid,
        &fastroute_backbone_hybrid,
        &interleaved_hybrid,
        &wirelength_sweep_hybrid,
        &extreme_sweep_hybrid,
        &radical_refine_hybrid,
        &envelope_hybrid,
        &wirelength_oracle_hybrid};
    NetRouteMap wirelength_closure_hybrid = buildWirelengthClosureHybrid(
        wirelength_oracle_hybrid,
        closure_donors,
        grouter_->db_net_map_,
        wirelength_oracle_score.vias,
        wirelength_oracle_low_layer_wl,
        closure_stats);
    const RouteScore wirelength_closure_score
        = computeRouteScore(wirelength_closure_hybrid);
    const uint64_t wirelength_closure_low_layer_wl
        = computeLowLayerWirelength(wirelength_closure_hybrid);
    WirelengthFusionStats fusion_stats;
    std::vector<const NetRouteMap*> fusion_donors{
        &fastroute_routes,
        &routes,
        &newgr_backbone_hybrid,
        &fastroute_backbone_hybrid,
        &interleaved_hybrid,
        &wirelength_sweep_hybrid,
        &extreme_sweep_hybrid,
        &radical_refine_hybrid,
        &envelope_hybrid,
        &wirelength_oracle_hybrid,
        &wirelength_closure_hybrid};
    NetRouteMap wirelength_fusion_hybrid = buildWirelengthFusionHybrid(
        wirelength_closure_hybrid,
        fusion_donors,
        grouter_->db_net_map_,
        wirelength_closure_score.vias,
        wirelength_closure_low_layer_wl,
        fusion_stats);
    const RouteScore wirelength_fusion_score
        = computeRouteScore(wirelength_fusion_hybrid);
    const uint64_t wirelength_fusion_low_layer_wl
        = computeLowLayerWirelength(wirelength_fusion_hybrid);

    logger_->info(utl::GRT,
                  6006,
                  "NEWGR candidate summary: NEWGR(ofl={}, wl={}, vias={}, low_wl={}, nets={}), "
                  "NEWGR_BB_HYBRID(wl={}, vias={}, low_wl={}, nets={}, fr_swap={}, add={}), "
                  "FR_BB_HYBRID(wl={}, vias={}, low_wl={}, nets={}, ng_swap={}, add={}, "
                  "cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "swap_limit_skip={}, wl_gain={}, polish_swap={}, polish_wl_gain={}, "
                  "micro_swap={}, micro_wl_gain={}), "
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
                  fastroute_backbone_stats.polish_swaps,
                  fastroute_backbone_stats.polish_wl_gain,
                  fastroute_backbone_stats.micro_swaps,
                  fastroute_backbone_stats.micro_wl_gain,
                  fastroute_overflow,
                  fastroute_score.wirelength,
                  fastroute_score.vias,
                  fastroute_low_layer_wl,
                  fastroute_score.routed_nets);
    logger_->info(utl::GRT,
                  6010,
                  "NEWGR interleaved hybrid summary: "
                  "INTERLEAVED(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  interleaved_score.wirelength,
                  interleaved_score.vias,
                  interleaved_low_layer_wl,
                  interleaved_score.routed_nets,
                  interleaved_stats.replaced_with_donor,
                  interleaved_stats.added_missing_nets,
                  interleaved_stats.candidate_pool_size,
                  interleaved_stats.skipped_by_via_guard,
                  interleaved_stats.skipped_by_layer_guard,
                  interleaved_stats.skipped_by_budget_guard,
                  interleaved_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6011,
                  "NEWGR wirelength sweep summary: "
                  "SWEEP(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  sweep_score.wirelength,
                  sweep_score.vias,
                  sweep_low_layer_wl,
                  sweep_score.routed_nets,
                  sweep_stats.replaced_with_donor,
                  sweep_stats.added_missing_nets,
                  sweep_stats.candidate_pool_size,
                  sweep_stats.skipped_by_via_guard,
                  sweep_stats.skipped_by_layer_guard,
                  sweep_stats.skipped_by_budget_guard,
                  sweep_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6012,
                  "NEWGR extreme sweep summary: "
                  "EXTREME_SWEEP(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  extreme_sweep_score.wirelength,
                  extreme_sweep_score.vias,
                  extreme_sweep_low_layer_wl,
                  extreme_sweep_score.routed_nets,
                  extreme_sweep_stats.replaced_with_donor,
                  extreme_sweep_stats.added_missing_nets,
                  extreme_sweep_stats.candidate_pool_size,
                  extreme_sweep_stats.skipped_by_via_guard,
                  extreme_sweep_stats.skipped_by_layer_guard,
                  extreme_sweep_stats.skipped_by_budget_guard,
                  extreme_sweep_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6013,
                  "NEWGR radical refine summary: "
                  "RADICAL_REFINE(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  radical_refine_score.wirelength,
                  radical_refine_score.vias,
                  radical_refine_low_layer_wl,
                  radical_refine_score.routed_nets,
                  radical_refine_stats.replaced_with_donor,
                  radical_refine_stats.added_missing_nets,
                  radical_refine_stats.candidate_pool_size,
                  radical_refine_stats.skipped_by_via_guard,
                  radical_refine_stats.skipped_by_layer_guard,
                  radical_refine_stats.skipped_by_budget_guard,
                  radical_refine_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6014,
                  "NEWGR envelope summary: "
                  "ENVELOPE(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  envelope_score.wirelength,
                  envelope_score.vias,
                  envelope_low_layer_wl,
                  envelope_score.routed_nets,
                  envelope_stats.replaced_with_donor,
                  envelope_stats.added_missing_nets,
                  envelope_stats.candidate_pool_size,
                  envelope_stats.skipped_by_via_guard,
                  envelope_stats.skipped_by_layer_guard,
                  envelope_stats.skipped_by_budget_guard,
                  envelope_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6018,
                  "NEWGR envelope+newgr graft summary: "
                  "ENVELOPE_NEWGR(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  envelope_newgr_score.wirelength,
                  envelope_newgr_score.vias,
                  envelope_newgr_low_layer_wl,
                  envelope_newgr_score.routed_nets,
                  envelope_newgr_stats.replaced_with_donor,
                  envelope_newgr_stats.added_missing_nets,
                  envelope_newgr_stats.candidate_pool_size,
                  envelope_newgr_stats.skipped_by_via_guard,
                  envelope_newgr_stats.skipped_by_layer_guard,
                  envelope_newgr_stats.skipped_by_budget_guard,
                  envelope_newgr_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6015,
                  "NEWGR wirelength oracle summary: "
                  "WL_ORACLE(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  wirelength_oracle_score.wirelength,
                  wirelength_oracle_score.vias,
                  wirelength_oracle_low_layer_wl,
                  wirelength_oracle_score.routed_nets,
                  wirelength_oracle_stats.replaced_with_donor,
                  wirelength_oracle_stats.added_missing_nets,
                  wirelength_oracle_stats.candidate_pool_size,
                  wirelength_oracle_stats.skipped_by_via_guard,
                  wirelength_oracle_stats.skipped_by_layer_guard,
                  wirelength_oracle_stats.skipped_by_budget_guard,
                  wirelength_oracle_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6016,
                  "NEWGR wirelength closure summary: "
                  "WL_CLOSURE(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  wirelength_closure_score.wirelength,
                  wirelength_closure_score.vias,
                  wirelength_closure_low_layer_wl,
                  wirelength_closure_score.routed_nets,
                  closure_stats.replaced_with_donor,
                  closure_stats.added_missing_nets,
                  closure_stats.candidate_pool_size,
                  closure_stats.skipped_by_via_guard,
                  closure_stats.skipped_by_layer_guard,
                  closure_stats.skipped_by_budget_guard,
                  closure_stats.consumed_wl_gain);
    logger_->info(utl::GRT,
                  6017,
                  "NEWGR wirelength fusion summary: "
                  "WL_FUSION(wl={}, vias={}, low_wl={}, nets={}, donor_swap={}, "
                  "add={}, cand={}, via_guard_skip={}, layer_guard_skip={}, "
                  "budget_skip={}, wl_gain={})",
                  wirelength_fusion_score.wirelength,
                  wirelength_fusion_score.vias,
                  wirelength_fusion_low_layer_wl,
                  wirelength_fusion_score.routed_nets,
                  fusion_stats.replaced_with_donor,
                  fusion_stats.added_missing_nets,
                  fusion_stats.candidate_pool_size,
                  fusion_stats.skipped_by_via_guard,
                  fusion_stats.skipped_by_layer_guard,
                  fusion_stats.skipped_by_budget_guard,
                  fusion_stats.consumed_wl_gain);

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

    if (shouldPreferInterleavedHybrid(best_overflow,
                                      best_score,
                                      best_low_layer_wl,
                                      last_total_overflow_,
                                      interleaved_score,
                                      interleaved_low_layer_wl)) {
      routes = std::move(interleaved_hybrid);
      best_score = interleaved_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = interleaved_low_layer_wl;
      selected_label = "FastRoute+NEWGR interleaved-graft";
      selected_hybrid = true;
      selected_swapped_nets = interleaved_stats.replaced_with_donor;
      selected_added_nets = interleaved_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferWirelengthSweepHybrid(best_overflow,
                                          best_score,
                                          best_low_layer_wl,
                                          last_total_overflow_,
                                          sweep_score,
                                          sweep_low_layer_wl)) {
      routes = std::move(wirelength_sweep_hybrid);
      best_score = sweep_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = sweep_low_layer_wl;
      selected_label = "FastRoute+NEWGR wirelength-sweep";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferWirelengthSweepHybrid(best_overflow,
                                          best_score,
                                          best_low_layer_wl,
                                          last_total_overflow_,
                                          extreme_sweep_score,
                                          extreme_sweep_low_layer_wl)) {
      routes = std::move(extreme_sweep_hybrid);
      best_score = extreme_sweep_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = extreme_sweep_low_layer_wl;
      selected_label = "FastRoute+NEWGR extreme-sweep";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferRadicalRefineHybrid(best_overflow,
                                        best_score,
                                        best_low_layer_wl,
                                        last_total_overflow_,
                                        radical_refine_score,
                                        radical_refine_low_layer_wl)) {
      routes = std::move(radical_refine_hybrid);
      best_score = radical_refine_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = radical_refine_low_layer_wl;
      selected_label = "FastRoute+NEWGR radical-refine";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferEnvelopeHybrid(best_overflow,
                                   best_score,
                                   best_low_layer_wl,
                                   last_total_overflow_,
                                   envelope_score,
                                   envelope_low_layer_wl)) {
      routes = std::move(envelope_hybrid);
      best_score = envelope_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = envelope_low_layer_wl;
      selected_label = "FastRoute+NEWGR multi-router-envelope";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor
                              + envelope_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets
                            + envelope_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (shouldPreferInterleavedHybrid(best_overflow,
                                      best_score,
                                      best_low_layer_wl,
                                      last_total_overflow_,
                                      envelope_newgr_score,
                                      envelope_newgr_low_layer_wl)) {
      routes = std::move(envelope_newgr_hybrid);
      best_score = envelope_newgr_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = envelope_newgr_low_layer_wl;
      selected_label = "FastRoute+NEWGR envelope-newgr-graft";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor
                              + envelope_stats.replaced_with_donor
                              + envelope_newgr_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets
                            + envelope_stats.added_missing_nets
                            + envelope_newgr_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (false
        && shouldPreferWirelengthOracleHybrid(best_overflow,
                                              best_score,
                                              best_low_layer_wl,
                                              last_total_overflow_,
                                              wirelength_oracle_score,
                                              wirelength_oracle_low_layer_wl)) {
      routes = std::move(wirelength_oracle_hybrid);
      best_score = wirelength_oracle_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = wirelength_oracle_low_layer_wl;
      selected_label = "FastRoute+NEWGR wirelength-oracle";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor
                              + envelope_stats.replaced_with_donor
                              + wirelength_oracle_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets
                            + envelope_stats.added_missing_nets
                            + wirelength_oracle_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (false
        && shouldPreferWirelengthClosureHybrid(best_overflow,
                                               best_score,
                                               best_low_layer_wl,
                                               last_total_overflow_,
                                               wirelength_closure_score,
                                               wirelength_closure_low_layer_wl)) {
      routes = std::move(wirelength_closure_hybrid);
      best_score = wirelength_closure_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = wirelength_closure_low_layer_wl;
      selected_label = "FastRoute+NEWGR wirelength-closure";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor
                              + envelope_stats.replaced_with_donor
                              + wirelength_oracle_stats.replaced_with_donor
                              + closure_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets
                            + envelope_stats.added_missing_nets
                            + wirelength_oracle_stats.added_missing_nets
                            + closure_stats.added_missing_nets;
      used_fastroute_last_run_ = false;
    }

    if (false
        && shouldPreferWirelengthFusionHybrid(best_overflow,
                                              best_score,
                                              best_low_layer_wl,
                                              last_total_overflow_,
                                              wirelength_fusion_score,
                                              wirelength_fusion_low_layer_wl)) {
      routes = std::move(wirelength_fusion_hybrid);
      best_score = wirelength_fusion_score;
      best_overflow = last_total_overflow_;
      best_low_layer_wl = wirelength_fusion_low_layer_wl;
      selected_label = "FastRoute+NEWGR wirelength-fusion";
      selected_hybrid = true;
      selected_swapped_nets = sweep_stats.replaced_with_donor
                              + extreme_sweep_stats.replaced_with_donor
                              + radical_refine_stats.replaced_with_donor
                              + envelope_stats.replaced_with_donor
                              + wirelength_oracle_stats.replaced_with_donor
                              + closure_stats.replaced_with_donor
                              + fusion_stats.replaced_with_donor;
      selected_added_nets = sweep_stats.added_missing_nets
                            + extreme_sweep_stats.added_missing_nets
                            + radical_refine_stats.added_missing_nets
                            + envelope_stats.added_missing_nets
                            + wirelength_oracle_stats.added_missing_nets
                            + closure_stats.added_missing_nets
                            + fusion_stats.added_missing_nets;
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
