#include "NEWGR/NewGR.h"

#include <algorithm>
#include <memory>
#include <unordered_set>
#include <vector>

#include "FastRoute.h"
#include "Net.h"
#include "Pin.h"
#include "src/NewgrEngine.h"
#include "utl/Logger.h"

namespace grt {

namespace {

GSegment normalizeSegment(const GSegment& segment)
{
  GSegment normalized = segment;
  if (normalized.isVia() && normalized.init_layer > normalized.final_layer) {
    std::swap(normalized.init_layer, normalized.final_layer);
  }
  return normalized;
}

void cleanupRouteSegments(NetRouteMap& routes)
{
  for (auto& [ignored_db_net, route] : routes) {
    static_cast<void>(ignored_db_net);
    GRoute filtered;
    filtered.reserve(route.size());
    std::unordered_set<GSegment, GSegmentHash> seen;
    seen.reserve(route.size());
    for (const GSegment& segment : route) {
      const bool is_zero_length_stub
          = segment.init_x == segment.final_x
            && segment.init_y == segment.final_y
            && segment.init_layer == segment.final_layer;
      if (is_zero_length_stub) {
        continue;
      }
      const GSegment normalized_segment = normalizeSegment(segment);
      if (seen.insert(normalized_segment).second) {
        filtered.push_back(normalized_segment);
      }
    }
    route.swap(filtered);
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

  // Start from FastRoute guides and run NEWGR-specific cleanup.
  NetRouteMap routes = grouter_->fastroute_->run();
  if (!routes.empty()) {
    grouter_->addRemainingGuides(
        routes, nets, min_routing_layer, max_routing_layer);
    grouter_->connectPadPins(routes);
    for (auto& net_route : routes) {
      std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
      GRoute& route = net_route.second;
      grouter_->mergeSegments(pins, route);
    }
    cleanupRouteSegments(routes);
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
  cleanupRouteSegments(routes);
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
