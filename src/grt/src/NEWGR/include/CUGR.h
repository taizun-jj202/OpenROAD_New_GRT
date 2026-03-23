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
  double weight_wire_length = 0.9;
  double weight_via_number = 2.0;
  double weight_short_area = 320.0;

  int min_routing_layer = 1;

  double cost_logistic_slope = 0.75;

  // allowed stem length increase to trunk length ratio
  double max_detour_ratio = 0.18;
  int target_detour_count = 12;

  double via_multiplier = 1.4;

  double maze_logistic_slope = 0.38;

  // Hybrid FastRoute/CUGR cost shaping:
  // keep under-utilized edges close to pure wirelength, ramp near saturation,
  // then switch to strong linear overflow cost.
  double wl_relaxation_util_threshold = 0.82;
  double wl_relaxation_penalty_floor = 0.08;
  double overflow_linear_penalty = 6.0;

  // SPRoute-style soft capacity model (utilization driven)
  bool use_soft_capacity = false;
  double soft_cap_min_ratio = 0.78;
  double soft_cap_max_ratio = 0.96;
  double soft_cap_mid_util = 0.72;
  double soft_cap_slope = 7.0;

  // FastRoute-style critical-net refinement schedule
  bool wirelength_first_refinement = true;
  int refinement_overflow_threshold = 1;
  int refinement_hpwl_threshold = 80;
  double detour_refine_ratio = 0.90;
  double maze_refine_ratio = 0.72;

  // FastRoute-style post-processing pass: reroute a subset of long,
  // overflow-free nets and keep only wirelength-improving solutions.
  bool enable_wirelength_recovery = true;
  bool recovery_use_maze = true;
  int recovery_hpwl_threshold = 72;
  double recovery_refine_ratio = 0.70;
  int recovery_maze_sparse_x = 5;
  int recovery_maze_sparse_y = 5;

  double pin_patch_threshold = 20.0;
  int pin_patch_padding = 1;
  double wire_patch_threshold = 2.0;
  double wire_patch_inflation_rate = 1.2;

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
  void wirelengthRecovery(const std::vector<int>& netIndices);
  void sortNetIndices(std::vector<int>& netIndices) const;
  void getGuides(const GRNet* net,
                 std::vector<std::pair<int, BoxT>>& guides);
  std::vector<int> selectCriticalNets(const std::vector<int>& candidates,
                                      double reroute_ratio) const;
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
