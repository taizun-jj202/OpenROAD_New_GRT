#pragma once

#include <csignal>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "grt/GRoute.h"

namespace odb {
class dbDatabase;
}  // namespace odb

namespace stt {
class SteinerTreeBuilder;
}  // namespace stt

namespace utl {
class Logger;
}  // namespace utl

namespace grt::newgr {

class Design;
class GridGraph;
class GRNet;
class BoxT;

struct Constants
{
  double weight_wire_length = 0.5;
  double weight_via_number = 6.3;
  double weight_short_area = 500.0;

  int min_routing_layer = 1;

  double cost_logistic_slope = 1.0;

  // allowed stem length increase to trunk length ratio
  double max_detour_ratio = 0.06;
  int target_detour_count = 2;
  int detour_overflow_threshold = 2;

  double via_multiplier = 2.25;
  // Require a meaningful gain before switching a child branch to a higher
  // layer during layer assignment.
  double layer_assignment_hysteresis_ratio = 0.30;
  // Discourage unnecessary use of higher layers even when congestion is low.
  double layer_usage_penalty_ratio = 0.048;

  double maze_logistic_slope = 0.5;
  int maze_overflow_threshold = 6;

  // Apply guide patching conservatively for low-overflow nets to avoid extra
  // detailed-route detours and vias.
  int guide_patch_overflow_threshold = 2;
  double pin_patch_threshold = 18.0;
  int pin_patch_padding = 1;
  double wire_patch_threshold = 1.2;
  double wire_patch_inflation_rate = 1.03;

  bool write_heatmap = false;
};

class CUGR
{
 public:
  CUGR(odb::dbDatabase* db,
       utl::Logger* log,
       stt::SteinerTreeBuilder* stt_builder);
  ~CUGR();
  void init(int min_routing_layer, int max_routing_layer);
  void route();
  void write(const std::string& guide_file);
  NetRouteMap getRoutes();
  void updateDbCongestion();

 private:
  void updateOverflowNets(std::vector<int>& netIndices);
  void patternRoute(std::vector<int>& netIndices);
  void patternRouteWithDetours(std::vector<int>& netIndices);
  void mazeRoute(std::vector<int>& netIndices);
  void sortNetIndices(std::vector<int>& netIndices, bool large_first) const;
  void getGuides(const GRNet* net,
                 std::vector<std::pair<int, BoxT>>& guides);
  void printStatistics() const;

  std::unique_ptr<Design> design_;
  std::unique_ptr<GridGraph> grid_graph_;
  std::vector<std::unique_ptr<GRNet>> gr_nets_;

  odb::dbDatabase* db_;
  utl::Logger* logger_;
  stt::SteinerTreeBuilder* stt_builder_;

  Constants constants_;

  int area_of_pin_patches_ = 0;
  int area_of_wire_patches_ = 0;
};

}  // namespace grt::newgr
