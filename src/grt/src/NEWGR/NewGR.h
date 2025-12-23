#pragma once

#include "grt/GlobalRouter.h"

namespace grt {

class NewGR
{
 public:
  NewGR(GlobalRouter* grouter,
        FastRouteCore* fastroute,
        CUGR* cugr,
        utl::Logger* logger);

  NetRouteMap run(std::vector<Net*>& nets,
                  int min_routing_layer,
                  int max_routing_layer);

 private:
  GlobalRouter* grouter_;
  FastRouteCore* fastroute_;
  CUGR* cugr_;
  utl::Logger* logger_;
};

}  // namespace grt
