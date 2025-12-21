// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025

#pragma once

#include <map>
#include <set>
#include <vector>

#include "grt/GRoute.h"

namespace utl {
class Logger;
}  // namespace utl

namespace grt {

class CUGR;
class FastRouteCore;
class GlobalRouter;
class Net;

// NEW_GR1 blends FastRoute's multi-stage routing, SPRoute's soft capacity
// heuristics, and CUGR's probabilistic awareness to target better congestion,
// via count, and wirelength simultaneously.
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
  struct Config
  {
    double wl_weight = 1.0;
    double via_weight = 1.08;
    double congestion_weight = 0.0012;
    double improvement_guard = 0.99;  // require >=1% score reduction
  };

  void buildRudyCache();
  NetRouteMap runFastRoute(std::vector<Net*>& nets,
                           int min_layer,
                           int max_layer);
  NetRouteMap runCugr(int min_layer, int max_layer);
  void collectCongestedNets(std::set<odb::dbNet*>& congested_nets);
  double scoreRoute(const GRoute& route) const;
  double computeWirelength(const GRoute& route) const;
  int computeViaCount(const GRoute& route) const;
  double computeCongestionCost(const GRoute& route) const;

  GlobalRouter* grouter_;
  FastRouteCore* fastroute_;
  CUGR* cugr_;
  utl::Logger* logger_;
  Config config_;

  std::vector<std::vector<double>> rudy_cache_;
  int tiles_x_{0};
  int tiles_y_{0};
  int tile_size_{0};
};

}  // namespace grt
