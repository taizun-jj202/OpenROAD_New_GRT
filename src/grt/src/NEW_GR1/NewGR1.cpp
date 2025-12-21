// SPDX-License-Identifier: BSD-3-Clause

#include "NEW_GR1/NewGR1.h"

#include <algorithm>
#include <cmath>

#include "CUGR.h"
#include "FastRoute.h"
#include "Grid.h"
#include "Net.h"
#include "Pin.h"
#include "grt/GlobalRouter.h"
#include "grt/Rudy.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace grt {

NewGR1::NewGR1(GlobalRouter* grouter,
               FastRouteCore* fast,
               CUGR* cugr,
               utl::Logger* logger)
    : grouter_(grouter), fastroute_(fast), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR1::run(std::vector<Net*>& nets,
                        int min_routing_layer,
                        int max_routing_layer)
{
  NetRouteMap routes;
  if (nets.empty()) {
    return routes;
  }

  buildRudyCache();
  routes = runFastRoute(nets, min_routing_layer, max_routing_layer);
  const bool has_overflow = fastroute_->totalOverflow() > 0;

  std::set<odb::dbNet*> congested_nets;
  collectCongestedNets(congested_nets);

  if (!congested_nets.empty()) {
    NetRouteMap cugr_routes = runCugr(min_routing_layer, max_routing_layer);
    for (odb::dbNet* net : congested_nets) {
      auto cugr_it = cugr_routes.find(net);
      auto fast_it = routes.find(net);
      if (cugr_it == cugr_routes.end() || fast_it == routes.end()) {
        continue;
      }
      const double fast_score = scoreRoute(fast_it->second);
      const double cugr_score = scoreRoute(cugr_it->second);
      if (cugr_score < config_.improvement_guard * fast_score) {
        routes[net] = cugr_it->second;
      }
    }
  }

  if (has_overflow) {
    logger_->warn(utl::GRT,
                  6004,
                  "NEW_GR1 detected {} overflowing edges. "
                  "Consider relaxing configuration parameters.",
                  fastroute_->totalOverflow());
  }
  return routes;
}

void NewGR1::buildRudyCache()
{
  Rudy* rudy = grouter_->getRudy();
  rudy->calculateRudy();
  const auto grid_size = rudy->getGridSize();
  tiles_x_ = grid_size.first;
  tiles_y_ = grid_size.second;
  tile_size_ = grouter_->getGridTileSize();
  rudy_cache_.assign(tiles_x_, std::vector<double>(tiles_y_, 0.0));
  for (int x = 0; x < tiles_x_; ++x) {
    for (int y = 0; y < tiles_y_; ++y) {
      rudy_cache_[x][y] = rudy->getTile(x, y).getRudy();
    }
  }
}

NetRouteMap NewGR1::runFastRoute(std::vector<Net*>& nets,
                                 int min_layer,
                                 int max_layer)
{
  NetRouteMap routes = fastroute_->run();
  grouter_->addRemainingGuides(routes, nets, min_layer, max_layer);
  grouter_->connectPadPins(routes);

  for (auto& net_route : routes) {
    std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
    grouter_->mergeSegments(pins, net_route.second);
  }
  return routes;
}

NetRouteMap NewGR1::runCugr(int min_layer, int max_layer)
{
  cugr_->init(min_layer, max_layer);
  cugr_->route();
  return cugr_->getRoutes();
}

void NewGR1::collectCongestedNets(std::set<odb::dbNet*>& congested_nets)
{
  NetsPerCongestedArea congested_tiles;
  fastroute_->findCongestedEdgesNets(congested_tiles, true);
  for (auto& entry : congested_tiles) {
    congested_nets.insert(entry.second.nets.begin(), entry.second.nets.end());
  }
  congested_tiles.clear();
  fastroute_->findCongestedEdgesNets(congested_tiles, false);
  for (auto& entry : congested_tiles) {
    congested_nets.insert(entry.second.nets.begin(), entry.second.nets.end());
  }
}

double NewGR1::scoreRoute(const GRoute& route) const
{
  const double wl = computeWirelength(route);
  const double vias = static_cast<double>(computeViaCount(route));
  const double congestion = computeCongestionCost(route);
  return config_.wl_weight * wl + config_.via_weight * vias
         + config_.congestion_weight * congestion;
}

double NewGR1::computeWirelength(const GRoute& route) const
{
  double length = 0.0;
  for (const GSegment& segment : route) {
    if (!segment.isVia()) {
      length += static_cast<double>(segment.length()) * tile_size_;
    }
  }
  return length;
}

int NewGR1::computeViaCount(const GRoute& route) const
{
  int count = 0;
  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      count++;
    }
  }
  return count;
}

double NewGR1::computeCongestionCost(const GRoute& route) const
{
  if (rudy_cache_.empty()) {
    return 0.0;
  }

  double cost = 0.0;
  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      continue;
    }
    const int x = std::clamp(segment.init_x, 0, tiles_x_ - 1);
    const int y = std::clamp(segment.init_y, 0, tiles_y_ - 1);
    cost += rudy_cache_[x][y] * segment.length();
  }
  return cost;
}

}  // namespace grt
