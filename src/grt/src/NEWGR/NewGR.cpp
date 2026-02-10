#include "NEWGR/NewGR.h"

#include "newgr/CUGR.h"
#include "utl/Logger.h"

namespace grt {

using utl::GRT;

NewGR::NewGR(GlobalRouter* grouter, newgr::CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  static_cast<void>(grouter_);
  static_cast<void>(nets);

  if (cugr_ == nullptr) {
    logger_->error(GRT, 6003, "NEWGR selected but CUGR engine is null.");
  }

  cugr_->init(min_routing_layer, max_routing_layer);
  cugr_->route();
  return cugr_->getRoutes();
}

}  // namespace grt
