#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "grt/GRoute.h"
#include "grt/GlobalRouter.h"

namespace parser {
class grGenerator;
}

namespace utl {
class Logger;
}

namespace grt {

struct NewgrInputNet
{
  std::string name;
  odb::dbNet* db_net{nullptr};
  bool is_clock{false};
  int min_layer{0};
  int max_layer{0};
  int root_pin_index{-1};
  std::vector<RoutePt> pins;
};

struct NewgrInput
{
  SprouteGridData grid;
  std::vector<NewgrInputNet> nets;
};

// Thin wrapper that will eventually host the full SPRoute implementation.
class NewgrEngine
{
 public:
  explicit NewgrEngine(utl::Logger* logger);

  void init(const SprouteGridData& grid, const std::vector<SprouteNetData>& nets);
  NetRouteMap run();
  NetRouteMap runWirelengthFirst();
  NetRouteMap runDataDrivenWirelength();
  NetRouteMap runRegionAware();
  NetRouteMap runRegularRegionAware();
  NetRouteMap runFineGrainRefine();

  const NewgrInput& getInput() const { return input_; }
  int getTotalOverflow() const { return last_total_overflow_; }
  void updateDbCongestion(odb::dbBlock* block);

 private:
  void buildInput();
  void prepareLefDefMetadata();
  parser::grGenerator buildGenerator() const;
  int coordFromIndex(int index, bool is_y = false) const;
  std::string routingLayerName(int layer_index) const;
  std::string cutLayerName(int layer_index) const;
  NetRouteMap extractRoutes() const;
  void appendRouteSegments(int net_id, grt::GRoute& route) const;
  void addSegment(grt::GRoute& route,
                  int grid_x0,
                  int grid_y0,
                  int grid_l0,
                  int grid_x1,
                  int grid_y1,
                  int grid_l1) const;
  NetRouteMap runWithConfig(int max_maze_round, int algo_id, int warn_id);
  int toDbLayer(int sproute_layer) const;
  std::string sanitizeNetName(const std::string& name) const;

  utl::Logger* logger_;
  SprouteGridData grid_;
  std::vector<SprouteNetData> nets_;
  NewgrInput input_;
  bool input_ready_{false};
  int last_total_overflow_{0};
};

}  // namespace grt
