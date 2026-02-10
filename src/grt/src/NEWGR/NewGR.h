#pragma once

#include <vector>

// Keep this header lightweight to avoid pulling the full GlobalRouter
// dependency graph into every translation unit.
#include "grt/GRoute.h"

namespace utl {
class Logger;
}  // namespace utl

namespace grt {

class GlobalRouter;
class CUGR;
class Net;
namespace newgr {
class CUGR;
}  // namespace newgr

class NewGR
{
 public:
  NewGR(GlobalRouter* grouter, newgr::CUGR* cugr, utl::Logger* logger);

  NetRouteMap run(std::vector<Net*>& nets,
                  int min_routing_layer,
                  int max_routing_layer);

 private:
  GlobalRouter* grouter_;
  newgr::CUGR* cugr_;
  utl::Logger* logger_;
};

}  // namespace grt
