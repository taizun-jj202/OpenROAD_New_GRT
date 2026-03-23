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
  double weight_wire_length = 1.0;
  double weight_via_number = 0.8;
  double weight_short_area = 240.0;

  int min_routing_layer = 1;

  double cost_logistic_slope = 1.0;

  // allowed stem length increase to trunk length ratio
  double max_detour_ratio = 0.25;
  int target_detour_count = 20;

  double via_multiplier = 0.85;

  int mst_topology_pin_threshold = 7;
  double mst_bend_penalty = 3.0;
  int hub_topology_pin_threshold = 20;
  int pin_access_layer_extension = 1;
  bool enable_early_detour_stage = false;

  double maze_logistic_slope = 0.5;
  int maze_base_interval = 5;
  int maze_min_interval = 2;
  int global_rebalance_rounds = 1;
  bool enable_critical_compaction_stage = true;
  int critical_compaction_rounds = 2;
  int critical_compaction_interval = 1;
  double critical_compaction_net_ratio = 0.12;
  double critical_compaction_overflow_weight = 180.0;

  // SPRoute-inspired soft-capacity shaping.
  bool enable_soft_capacity = true;
  double soft_capacity_ratio_min = 0.82;
  double soft_capacity_ratio_max = 1.00;
  double soft_capacity_slope = 4.0;
  double soft_capacity_congestion_mid = 0.75;
  double soft_capacity_rudy_weight = 1.0;
  double soft_capacity_lower_layer_boost = 0.10;
  double soft_capacity_min_absolute = 0.75;

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
  void globalRebalanceRoute(const std::vector<int>& allNetIndices);
  void criticalCompactionRoute(const std::vector<int>& allNetIndices);
  void sortNetIndices(std::vector<int>& netIndices) const;
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
