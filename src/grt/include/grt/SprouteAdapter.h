#pragma once

#include <memory>
#include <vector>

#include "grt/GlobalRouter.h"

namespace utl {
class Logger;
}

namespace grt {

class SprouteEngine;

class SprouteAdapter
{
 public:
  explicit SprouteAdapter(utl::Logger* logger);
  ~SprouteAdapter();

  void initialize(const SprouteGridData& grid,
                  const std::vector<SprouteNetData>& nets);
  NetRouteMap run();
  bool ready() const { return initialized_; }
  int getTotalOverflow() const;
  void updateDbCongestion(odb::dbBlock* block);

 private:
  std::unique_ptr<SprouteEngine> engine_;
  bool initialized_{false};
};

}  // namespace grt
