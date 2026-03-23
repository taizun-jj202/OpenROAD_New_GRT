#pragma once

#include <memory>

#include "grt/GlobalRouter.h"

namespace grt {

class NewgrEngine;

class NewGR
{
 public:
  NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger);
  ~NewGR();

  NetRouteMap run(std::vector<Net*>& nets,
                  int min_routing_layer,
                  int max_routing_layer);

  int getTotalOverflow() const;
  void updateDbCongestion(odb::dbBlock* block);

 private:
  GlobalRouter* grouter_;
  CUGR* cugr_;
  utl::Logger* logger_;
  std::unique_ptr<NewgrEngine> engine_;
  int last_total_overflow_{0};
  bool used_fastroute_last_run_{false};
};

}  // namespace grt
