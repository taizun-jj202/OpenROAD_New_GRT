// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "NEW_GR1/NEW_GR_config.h"
#include "grt/GRoute.h"

namespace utl {
class Logger;
}  // namespace utl

namespace grt {

class CUGR;
class FastRouteCore;
class GlobalRouter;
class Net;

// NEW_GR1 baseline: run FastRoute and leave hooks for future enhancements.
class NewGR1
{
 public:
  NewGR1(GlobalRouter* grouter,
         FastRouteCore* fast,
         CUGR* cugr,
         utl::Logger* logger);

  NetRouteMap run(std::vector<Net*>& nets,
                  int min_routing_layer,
                  int max_routing_layer);

 private:
  NetRouteMap runFastRoute(std::vector<Net*>& nets,
                           int min_layer,
                           int max_layer);
  void reportBaselineMetrics(const NetRouteMap& routes) const;

  GlobalRouter* grouter_;
  FastRouteCore* fastroute_;
  CUGR* cugr_;
  utl::Logger* logger_;
  const newgr::Config& config_;
};

}  // namespace grt
