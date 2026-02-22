#include "NEWGR/NewGR.h"

#include <algorithm>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "FastRoute.h"
#include "src/NewgrEngine.h"
#include "Net.h"
#include "Pin.h"
#include "utl/Logger.h"

namespace grt {

namespace {

using SegmentKey = std::tuple<int, int, int, int, int, int>;

SegmentKey canonicalKey(const GSegment& segment)
{
  std::tuple<int, int, int> first{
      segment.init_x, segment.init_y, segment.init_layer};
  std::tuple<int, int, int> second{
      segment.final_x, segment.final_y, segment.final_layer};
  if (second < first) {
    std::swap(first, second);
  }
  return SegmentKey{std::get<0>(first),
                    std::get<1>(first),
                    std::get<2>(first),
                    std::get<0>(second),
                    std::get<1>(second),
                    std::get<2>(second)};
}

void removeDuplicateSegments(NetRouteMap& routes)
{
  for (auto& [ignored_db_net, route] : routes) {
    static_cast<void>(ignored_db_net);
    GRoute deduped;
    deduped.reserve(route.size());
    std::set<SegmentKey> seen;
    for (const GSegment& segment : route) {
      const SegmentKey key = canonicalKey(segment);
      if (seen.insert(key).second) {
        deduped.push_back(segment);
      }
    }
    route.swap(deduped);
  }
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
  active_backend_ = Backend::None;
  last_total_overflow_ = 0;
  if (nets.empty()) {
    return {};
  }

  // Quality-oriented pass: leverage FastRoute's lower wirelength/via tendency.
  NetRouteMap routes
      = grouter_->findRouting(nets, min_routing_layer, max_routing_layer);
  if (!routes.empty()) {
    removeDuplicateSegments(routes);
    active_backend_ = Backend::FastRoute;
    last_total_overflow_ = grouter_->fastroute_->totalOverflow();
    return routes;
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
  routes = engine_->run();

  grouter_->addRemainingGuides(routes, nets, min_routing_layer, max_routing_layer);
  grouter_->connectPadPins(routes);
  for (auto& net_route : routes) {
    std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
    GRoute& route = net_route.second;
    grouter_->mergeSegments(pins, route);
  }
  removeDuplicateSegments(routes);
  active_backend_ = Backend::NewgrEngine;
  last_total_overflow_ = engine_->getTotalOverflow();

  return routes;
}

int NewGR::getTotalOverflow() const
{
  if (active_backend_ == Backend::NewgrEngine && engine_) {
    return engine_->getTotalOverflow();
  }
  return last_total_overflow_;
}

void NewGR::updateDbCongestion(odb::dbBlock* block)
{
  if (active_backend_ == Backend::FastRoute) {
    int min_routing_layer = 0;
    int max_routing_layer = 0;
    grouter_->getMinMaxLayer(min_routing_layer, max_routing_layer);
    grouter_->fastroute_->updateDbCongestion(min_routing_layer,
                                             max_routing_layer);
  } else if (active_backend_ == Backend::NewgrEngine && engine_) {
    engine_->updateDbCongestion(block);
  }
}

}  // namespace grt
