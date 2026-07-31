#include "NEWGR/NewGR.h"

#include "utl/Logger.h"

namespace grt {

using utl::GNR;

NewGR::NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  if (nets.empty()) {
    return {};
  }

  logger_->info(GNR, 6005, "NEWGR wirelength-first single-pass routing");
  return grouter_->findRouting(
      nets, min_routing_layer, max_routing_layer);
}

}  // namespace grt
