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
  double weight_via_number = 1.6;
  double weight_short_area = 320.0;

  int min_routing_layer = 1;

  double cost_logistic_slope = 0.75;

  // allowed stem length increase to trunk length ratio
  double max_detour_ratio = 0.18;
  int target_detour_count = 12;

  double via_multiplier = 1.0;

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
  int refinement_overflow_threshold = 2;
  int refinement_hpwl_threshold = 120;
  double refinement_stretch_threshold = 1.08;
  int refinement_max_selected_nets = 384;
  double detour_refine_ratio = 0.78;
  double maze_refine_ratio = 0.52;

  // Stage-3 (maze) exploration controls.
  int stage3_dense_sparse_x = 2;
  int stage3_dense_sparse_y = 2;
  bool stage3_try_offset = true;
  bool stage3_full_offset_sweep = false;
  int stage3_full_offset_hpwl_threshold = 180;
  int stage3_max_maze_configs = 8;
  bool stage3_use_wirelength_maze = true;
  int stage3_wl_only_hpwl_threshold = 140;
  int stage3_wl_config_limit = 3;
  double stage3_wl_via_cost_scale = 0.35;
  bool stage3_use_full_grid_wl_maze = true;
  int stage3_full_grid_hpwl_threshold = 260;
  int stage3_full_grid_pin_limit = 16;
  int stage3_full_grid_overflow_threshold = 1;
  double stage3_full_grid_via_cost_scale = 0.02;
  double stage3_full_grid_min_stretch = 1.30;

  // FastRoute-style post-processing pass: reroute a subset of long,
  // overflow-free nets and keep only wirelength-improving solutions.
  bool enable_wirelength_recovery = true;
  bool recovery_use_maze = true;
  bool recovery_use_wirelength_maze = true;
  double recovery_wl_maze_via_cost_scale = 0.30;
  int recovery_wl_only_hpwl_threshold = 140;
  int recovery_hpwl_threshold = 96;
  double recovery_min_stretch = 1.10;
  int recovery_candidate_cap = 64;
  int recovery_deep_search_cap = 24;
  int recovery_shallow_maze_config_limit = 2;
  int recovery_deep_maze_config_limit = 4;
  int recovery_shallow_wl_config_limit = 1;
  int recovery_deep_wl_config_limit = 2;
  double recovery_refine_ratio = 0.45;
  int recovery_max_passes = 1;
  double recovery_pass_decay = 0.80;
  int recovery_maze_sparse_x = 3;
  int recovery_maze_sparse_y = 3;
  double recovery_deep_ratio = 0.65;
  int recovery_deep_hpwl_threshold = 140;
  int recovery_aniso_sparse_long = 7;
  int recovery_aniso_sparse_short = 3;
  bool recovery_try_offset = false;
  bool recovery_full_offset_sweep = false;
  int recovery_full_offset_hpwl_threshold = 220;
  int recovery_dense_sparse_x = 2;
  int recovery_dense_sparse_y = 2;
  int recovery_max_maze_configs = 6;
  int recovery_max_via_increase = 6;
  bool recovery_use_full_grid_maze = false;
  int recovery_full_grid_top_n = 8;
  int recovery_full_grid_hpwl_threshold = 160;
  double recovery_full_grid_via_cost_scale = 0.03;

  double pin_patch_threshold = 1000000.0;
  int pin_patch_padding = 3;
  double wire_patch_threshold = 1000000.0;
  double wire_patch_inflation_rate = 1.0;

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
