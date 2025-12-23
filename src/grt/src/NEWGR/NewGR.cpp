#include "NEWGR/NewGR.h"

#include <cmath>

#include "utl/Logger.h"

namespace grt {

using utl::GNR;

NewGR::NewGR(GlobalRouter* grouter,
             FastRouteCore* fastroute,
             CUGR* cugr,
             utl::Logger* logger)
    : grouter_(grouter), fastroute_(fastroute), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  NetRouteMap routes;
  if (!nets.empty()) {
    routes = grouter_->findRouting(nets, min_routing_layer, max_routing_layer);
  }

  long wirelength_dbu = 0;
  long via_count = 0;
  for (auto& net_route : routes) {
    for (const GSegment& segment : net_route.second) {
      if (segment.isVia()) {
        via_count++;
      } else {
        wirelength_dbu += std::abs(segment.final_x - segment.init_x)
                          + std::abs(segment.final_y - segment.init_y);
      }
    }
  }

  double wirelength_um = 0.0;
  if (wirelength_dbu > 0 && grouter_->db_ != nullptr
      && grouter_->db_->getTech() != nullptr) {
    wirelength_um = wirelength_dbu
                    / static_cast<double>(
                        grouter_->db_->getTech()->getDbUnitsPerMicron());
  }

  logger_->info(GNR,
                6005,
                "NEWGR baseline: wirelength {:.0f} um, vias {:.0f}",
                wirelength_um,
                static_cast<double>(via_count));

  return routes;
}

}  // namespace grt
