#include "NEWGR/NewGR.h"

#include <algorithm>
#include <memory>
#include <tuple>
#include <vector>

#include "NEWGR/src/NewgrEngine.h"
#include "Net.h"
#include "Pin.h"
#include "utl/Logger.h"

namespace grt {

namespace {

bool isViaSegment(const GSegment& segment)
{
  return segment.init_x == segment.final_x && segment.init_y == segment.final_y
         && segment.init_layer != segment.final_layer;
}

void sortRouteForMerge(GRoute& route)
{
  std::sort(route.begin(), route.end(), [](const GSegment& lhs, const GSegment& rhs) {
    const bool lhs_is_via = isViaSegment(lhs);
    const bool rhs_is_via = isViaSegment(rhs);
    if (lhs_is_via != rhs_is_via) {
      return !lhs_is_via;
    }

    if (lhs_is_via) {
      const auto lhs_key = std::make_tuple(std::min(lhs.init_layer, lhs.final_layer),
                                           std::max(lhs.init_layer, lhs.final_layer),
                                           lhs.init_x,
                                           lhs.init_y);
      const auto rhs_key = std::make_tuple(std::min(rhs.init_layer, rhs.final_layer),
                                           std::max(rhs.init_layer, rhs.final_layer),
                                           rhs.init_x,
                                           rhs.init_y);
      return lhs_key < rhs_key;
    }

    const bool lhs_vertical = lhs.init_x == lhs.final_x;
    const bool rhs_vertical = rhs.init_x == rhs.final_x;
    if (lhs.init_layer != rhs.init_layer) {
      return lhs.init_layer < rhs.init_layer;
    }
    if (lhs_vertical != rhs_vertical) {
      return lhs_vertical < rhs_vertical;
    }

    if (lhs_vertical) {
      const auto lhs_key
          = std::make_tuple(lhs.init_x,
                            std::min(lhs.init_y, lhs.final_y),
                            std::max(lhs.init_y, lhs.final_y));
      const auto rhs_key
          = std::make_tuple(rhs.init_x,
                            std::min(rhs.init_y, rhs.final_y),
                            std::max(rhs.init_y, rhs.final_y));
      return lhs_key < rhs_key;
    }

    const auto lhs_key = std::make_tuple(lhs.init_y,
                                         std::min(lhs.init_x, lhs.final_x),
                                         std::max(lhs.init_x, lhs.final_x));
    const auto rhs_key = std::make_tuple(rhs.init_y,
                                         std::min(rhs.init_x, rhs.final_x),
                                         std::max(rhs.init_x, rhs.final_x));
    return lhs_key < rhs_key;
  });
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

  grouter_->addRemainingGuides(routes, nets, min_routing_layer, max_routing_layer);
  grouter_->connectPadPins(routes);
  for (auto& net_route : routes) {
    std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
    GRoute& route = net_route.second;
    sortRouteForMerge(route);
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
