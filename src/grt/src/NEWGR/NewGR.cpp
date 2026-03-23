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

bool preferNewgr(const RouteScore& fastroute_score,
                 const RouteScore& newgr_score,
                 int tile_size)
{
  // Prefer the NEWGR segment set if it is clearly better in WL and does not
  // regress vias too much. Otherwise use weighted score as a tie-breaker.
  if (newgr_score.wirelength + tile_size <= fastroute_score.wirelength
      && newgr_score.vias <= fastroute_score.vias + 2) {
    return true;
  }
  if (newgr_score.wirelength <= fastroute_score.wirelength
      && newgr_score.vias < fastroute_score.vias) {
    return true;
  }
  return newgr_score.weighted_cost + tile_size < fastroute_score.weighted_cost;
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
  NetRouteMap newgr_routes = engine_->run();

  const int tile_size = std::max(1, grouter_->getTileSize());
  const int via_penalty = std::max(1, tile_size / 4);
  int swapped_to_newgr = 0;
  int inserted_from_newgr = 0;

  for (auto& [db_net, newgr_route] : newgr_routes) {
    auto route_it = routes.find(db_net);
    if (route_it == routes.end()) {
      routes.emplace(db_net, std::move(newgr_route));
      inserted_from_newgr++;
      continue;
    }

    const RouteScore fastroute_score = scoreRoute(route_it->second, via_penalty);
    const RouteScore newgr_score = scoreRoute(newgr_route, via_penalty);
    if (preferNewgr(fastroute_score, newgr_score, tile_size)) {
      route_it->second = std::move(newgr_route);
      swapped_to_newgr++;
    }
  }

  logger_->info(utl::GRT,
                6004,
                "NEWGR hybrid selected {} NEWGR net routes (+{} NEWGR-only) "
                "out of {} total.",
                swapped_to_newgr,
                inserted_from_newgr,
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
