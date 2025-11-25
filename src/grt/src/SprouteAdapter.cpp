#include "grt/SprouteAdapter.h"

#include "SprouteEngine.h"
#include "utl/Logger.h"

namespace grt {

SprouteAdapter::SprouteAdapter(utl::Logger* logger)
    : engine_(std::make_unique<SprouteEngine>(logger))
{
}

SprouteAdapter::~SprouteAdapter() = default;

void SprouteAdapter::initialize(const SprouteGridData& grid,
                                const std::vector<SprouteNetData>& nets)
{
  engine_->init(grid, nets);
  initialized_ = true;
}

NetRouteMap SprouteAdapter::run()
{
  if (!initialized_) {
    return {};
  }
  return engine_->run();
}

int SprouteAdapter::getTotalOverflow() const
{
  if (!engine_) {
    return 0;
  }
  return engine_->getTotalOverflow();
}

void SprouteAdapter::updateDbCongestion(odb::dbBlock* block)
{
  if (!engine_) {
    return;
  }
  engine_->updateDbCongestion(block);
}

}  // namespace grt
