// SPDX-License-Identifier: BSD-3-Clause

#include "NEW_GR1/NewGR1.h"

#include "FastRoute.h"
#include "Net.h"
#include "Pin.h"
#include "grt/GlobalRouter.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace grt {

using utl::GRT;

NewGR1::NewGR1(GlobalRouter* grouter,
               FastRouteCore* fast,
               CUGR* cugr,
               utl::Logger* logger)
    : grouter_(grouter),
      fastroute_(fast),
      cugr_(cugr),
      logger_(logger),
      config_(newgr::getConfig())
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

  routes = runFastRoute(nets, min_routing_layer, max_routing_layer);

  if (config_.log_summary) {
    reportBaselineMetrics(routes);
  }

  return routes;
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

void NewGR1::reportBaselineMetrics(const NetRouteMap& routes) const
{
  if (routes.empty()) {
    return;
  }

  double wl_dbu = 0.0;
  double via_count = 0.0;
  const int tile_size = grouter_->getGridTileSize();

  for (const auto& [net, route] : routes) {
    for (const GSegment& segment : route) {
      if (segment.isVia()) {
        via_count += 1.0;
      } else {
        wl_dbu += static_cast<double>(segment.length()) * tile_size;
      }
    }
  }

  const double wl_um
      = wl_dbu / grouter_->db()->getTech()->getDbUnitsPerMicron();
  logger_->info(GRT,
                6005,
                "NEW_GR1 baseline: wirelength {:.0f} um, vias {:.0f}",
                wl_um,
                via_count);
}

}  // namespace grt
