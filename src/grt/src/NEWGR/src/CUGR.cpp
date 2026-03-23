#include "CUGR.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

#include "Design.h"
#include "GRNet.h"
#include "GRTree.h"
#include "GridGraph.h"
#include "Layers.h"
#include "MazeRoute.h"
#include "Netlist.h"
#include "PatternRoute.h"
#include "geo.h"
#include "grt/GRoute.h"
#include "odb/db.h"
#include "stt/SteinerTreeBuilder.h"
#include "utl/Logger.h"

namespace grt::newgr {

namespace {

struct RouteStats
{
  uint64_t wirelength = 0;
  int vias = 0;
  int overflow = 0;
  CapacityT total_overflow = std::numeric_limits<CapacityT>::max();
};

RouteStats measureRouteStats(const GridGraph* grid_graph,
                             const std::shared_ptr<GRTreeNode>& tree)
{
  RouteStats stats;
  if (!tree) {
    return stats;
  }

  GRTreeNode::preorder(
      tree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            const int direction
                = grid_graph->getLayerDirection(node->getLayerIdx());
            const int l = std::min((*node)[direction], (*child)[direction]);
            const int h = std::max((*node)[direction], (*child)[direction]);
            for (int c = l; c < h; c++) {
              stats.wirelength += grid_graph->getEdgeLength(direction, c);
            }
          } else {
            stats.vias += abs(node->getLayerIdx() - child->getLayerIdx());
          }
        }
      });
  stats.overflow = grid_graph->checkOverflow(tree);
  return stats;
}

bool isBetterStage3Candidate(const RouteStats& candidate,
                             const RouteStats& current_best,
                             const int baseline_overflow)
{
  constexpr double kOverflowEpsilon = 1e-6;
  constexpr double kAllowedOverflowIncreaseForWlGain = 16.0;
  constexpr double kStrongOverflowDropThreshold = 20.0;

  // Wirelength-first objective:
  // keep shorter candidates as long as they don't cause a large overflow jump.
  if (candidate.wirelength < current_best.wirelength) {
    return candidate.total_overflow
           <= current_best.total_overflow + kAllowedOverflowIncreaseForWlGain;
  }

  const double overflow_drop
      = current_best.total_overflow - candidate.total_overflow;
  if (overflow_drop > kStrongOverflowDropThreshold) {
    // Accept a large overflow reduction even without immediate WL gain.
    constexpr int64_t kMaxWirelengthTradeoff = 40;
    if (candidate.wirelength
        > current_best.wirelength + kMaxWirelengthTradeoff) {
      return false;
    }
    if (candidate.vias != current_best.vias) {
      return candidate.vias < current_best.vias;
    }
    return true;
  }

  if (std::abs(overflow_drop) > kOverflowEpsilon) {
    if (baseline_overflow > 0) {
      return overflow_drop > 0.0;
    }
    return overflow_drop > 0.0 && candidate.wirelength <= current_best.wirelength;
  }

  if (candidate.vias != current_best.vias) {
    return candidate.vias < current_best.vias;
  }
  return candidate.overflow < current_best.overflow;
}

}  // namespace

CUGR::CUGR(odb::dbDatabase* db,
           utl::Logger* log,
           stt::SteinerTreeBuilder* stt_builder)
    : db_(db), logger_(log), stt_builder_(stt_builder)
{
}

CUGR::~CUGR() = default;

void CUGR::init(const int min_routing_layer, const int max_routing_layer)
{
  design_ = std::make_unique<Design>(
      db_, logger_, constants_, min_routing_layer, max_routing_layer);
  grid_graph_ = std::make_unique<GridGraph>(design_.get(), constants_, logger_);
  // Instantiate the global routing netlist
  const std::vector<CUGRNet>& baseNets = design_->getAllNets();
  gr_nets_.reserve(baseNets.size());
  for (const CUGRNet& baseNet : baseNets) {
    gr_nets_.push_back(std::make_unique<GRNet>(baseNet, grid_graph_.get()));
  }
}

void CUGR::updateOverflowNets(std::vector<int>& netIndices)
{
  netIndices.clear();
  for (const auto& net : gr_nets_) {
    if (grid_graph_->checkOverflow(net->getRoutingTree()) > 0) {
      netIndices.push_back(net->getIndex());
    }
  }
  const int num_nets = gr_nets_.size();
  logger_->report("{} / {} nets have overflow.", netIndices.size(), num_nets);
}

void CUGR::patternRoute(std::vector<int>& netIndices)
{
  logger_->report("stage 1: pattern routing");
  sortNetIndices(netIndices);
  for (const int netIndex : netIndices) {
    PatternRoute patternRoute(gr_nets_[netIndex].get(),
                              grid_graph_.get(),
                              stt_builder_,
                              constants_,
                              logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();
    grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree());
  }

  updateOverflowNets(netIndices);
}

void CUGR::patternRouteWithDetours(std::vector<int>& netIndices)
{
  if (netIndices.empty()) {
    return;
  }
  logger_->report("stage 2: pattern routing with possible detours");
  // (2d) direction -> x -> y -> has overflow?
  GridGraphView<bool> congestionView;
  grid_graph_->extractCongestionView(congestionView);
  sortNetIndices(netIndices);
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    grid_graph_->commitTree(net->getRoutingTree(), /*ripup*/ true);
    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    // KEY DIFFERENCE compared to stage 1 (patternRoute)
    patternRoute.constructDetours(congestionView);
    patternRoute.run();
    grid_graph_->commitTree(net->getRoutingTree());
  }

  updateOverflowNets(netIndices);
}

void CUGR::mazeRoute(std::vector<int>& netIndices)
{
  if (netIndices.empty()) {
    return;
  }
  logger_->report("stage 3: maze routing on sparsified routing graph");

  std::vector<std::shared_ptr<GRTreeNode>> baseline_trees(gr_nets_.size(),
                                                           nullptr);
  std::vector<RouteStats> baseline_stats(gr_nets_.size());
  for (const int netIndex : netIndices) {
    baseline_trees[netIndex] = gr_nets_[netIndex]->getRoutingTree();
    baseline_stats[netIndex]
        = measureRouteStats(grid_graph_.get(), baseline_trees[netIndex]);
    if (baseline_trees[netIndex]) {
      grid_graph_->commitTree(baseline_trees[netIndex], /*ripup*/ true);
    }
  }

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  sortNetIndices(netIndices);

  struct MazeConfig
  {
    int sparse_x;
    int sparse_y;
    int offset_x;
    int offset_y;
  };

  int improved_nets = 0;
  int fallback_nets = 0;
  int evaluated_candidates = 0;

  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const auto baseline_tree = baseline_trees[netIndex];
    const RouteStats original_stats = baseline_stats[netIndex];

    RouteStats best_stats;
    std::shared_ptr<GRTreeNode> best_tree = nullptr;
    bool best_is_baseline = false;
    int baseline_overflow = original_stats.overflow;

    auto evaluateTree = [&](const std::shared_ptr<GRTreeNode>& tree,
                            RouteStats& stats) {
      if (!tree) {
        return false;
      }
      grid_graph_->commitTree(tree);
      stats = measureRouteStats(grid_graph_.get(), tree);
      stats.total_overflow = grid_graph_->getTotalOverflow();
      grid_graph_->commitTree(tree, /*ripup*/ true);
      return true;
    };

    if (baseline_tree && evaluateTree(baseline_tree, best_stats)) {
      best_tree = baseline_tree;
      best_is_baseline = true;
      baseline_overflow = best_stats.overflow;
    }

    std::vector<MazeConfig> maze_configs;
    const int max_maze_configs = std::max(6, constants_.stage3_max_maze_configs);
    maze_configs.reserve(max_maze_configs);
    auto addMazeConfig = [&](int sparse_x,
                             int sparse_y,
                             int offset_x,
                             int offset_y) {
      if (static_cast<int>(maze_configs.size()) >= max_maze_configs) {
        return;
      }
      sparse_x = std::max(2, sparse_x);
      sparse_y = std::max(2, sparse_y);
      offset_x = std::clamp(offset_x, 0, sparse_x - 1);
      offset_y = std::clamp(offset_y, 0, sparse_y - 1);
      for (const auto& cfg : maze_configs) {
        if (cfg.sparse_x == sparse_x && cfg.sparse_y == sparse_y
            && cfg.offset_x == offset_x && cfg.offset_y == offset_y) {
          return;
        }
      }
      maze_configs.push_back({sparse_x, sparse_y, offset_x, offset_y});
    };

    const BoxT& bbox = net->getBoundingBox();
    const int hpwl = bbox.hp();
    const bool wide_bbox = bbox.width() >= bbox.height();
    const int base_sparse = std::clamp(hpwl >= 240 ? 8 : (hpwl >= 120 ? 7 : 6),
                                       4,
                                       9);
    const int dense_sparse_x = std::max(2, constants_.stage3_dense_sparse_x);
    const int dense_sparse_y = std::max(2, constants_.stage3_dense_sparse_y);
    const int anis_long = std::min(10, base_sparse + 2);
    const int anis_short = std::max(3, base_sparse - 2);

    // Mix FastRoute-style shifted sparse grids with SPRoute-style denser local
    // search and anisotropic grids for elongated nets.
    addMazeConfig(base_sparse, base_sparse, 0, 0);
    addMazeConfig(dense_sparse_x, dense_sparse_y, 0, 0);
    if (wide_bbox) {
      addMazeConfig(anis_long, anis_short, 0, 0);
      addMazeConfig(anis_long, anis_short, anis_long / 2, anis_short / 2);
    } else {
      addMazeConfig(anis_short, anis_long, 0, 0);
      addMazeConfig(anis_short, anis_long, anis_short / 2, anis_long / 2);
    }
    if (constants_.stage3_try_offset) {
      addMazeConfig(base_sparse, base_sparse, base_sparse / 2, base_sparse / 2);
      addMazeConfig(
          dense_sparse_x, dense_sparse_y, dense_sparse_x / 2, dense_sparse_y / 2);
    }
    if (hpwl >= 220) {
      addMazeConfig(base_sparse - 1,
                    base_sparse + 1,
                    (base_sparse - 1) / 2,
                    (base_sparse + 1) / 2);
      addMazeConfig(base_sparse + 1,
                    base_sparse - 1,
                    (base_sparse + 1) / 2,
                    (base_sparse - 1) / 2);
    }
    if (constants_.stage3_full_offset_sweep
        && hpwl >= constants_.stage3_full_offset_hpwl_threshold) {
      for (int offset_x = 0;
           offset_x < base_sparse
           && static_cast<int>(maze_configs.size()) < max_maze_configs;
           offset_x++) {
        for (int offset_y = 0;
             offset_y < base_sparse
             && static_cast<int>(maze_configs.size()) < max_maze_configs;
             offset_y++) {
          addMazeConfig(base_sparse, base_sparse, offset_x, offset_y);
        }
      }
      for (int offset_x = 0;
           offset_x < dense_sparse_x
           && static_cast<int>(maze_configs.size()) < max_maze_configs;
           offset_x++) {
        for (int offset_y = 0;
             offset_y < dense_sparse_y
             && static_cast<int>(maze_configs.size()) < max_maze_configs;
             offset_y++) {
          addMazeConfig(dense_sparse_x, dense_sparse_y, offset_x, offset_y);
        }
      }
    }

    auto considerCandidate = [&](const std::shared_ptr<GRTreeNode>& tree,
                                 const bool is_baseline_candidate) {
      RouteStats candidate_stats;
      if (!evaluateTree(tree, candidate_stats)) {
        return;
      }
      if (!best_tree
          || isBetterStage3Candidate(
              candidate_stats, best_stats, baseline_overflow)) {
        best_tree = tree;
        best_stats = candidate_stats;
        best_is_baseline = is_baseline_candidate;
      }
    };

    for (const auto& cfg : maze_configs) {
      MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
      SparseGrid sparse_grid(cfg.sparse_x, cfg.sparse_y, cfg.offset_x, cfg.offset_y);
      mazeRoute.constructSparsifiedGraph(wireCostView, sparse_grid);
      mazeRoute.run();
      const std::shared_ptr<SteinerTreeNode> tree = mazeRoute.getSteinerTree();
      if (!tree) {
        continue;
      }

      PatternRoute patternRoute(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      patternRoute.setSteinerTree(tree);
      patternRoute.constructRoutingDAG();
      patternRoute.run();
      considerCandidate(net->getRoutingTree(), /*is_baseline_candidate*/ false);
      evaluated_candidates++;
    }

    if (!best_tree && baseline_tree) {
      best_tree = baseline_tree;
      best_stats = original_stats;
      best_stats.total_overflow = grid_graph_->getTotalOverflow();
      best_is_baseline = true;
    }

    if (best_tree) {
      net->setRoutingTree(best_tree);
      grid_graph_->commitTree(best_tree);
      grid_graph_->updateWireCostView(wireCostView, best_tree);
      if (!best_is_baseline && best_stats.wirelength < original_stats.wirelength) {
        improved_nets++;
      }
    } else {
      fallback_nets++;
      net->clearRoutingTree();
    }
  }

  logger_->report("stage 3 evaluated {} maze candidates; improved {} nets; {} "
                  "nets had no candidate tree.",
                  evaluated_candidates,
                  improved_nets,
                  fallback_nets);
  updateOverflowNets(netIndices);
}

void CUGR::wirelengthRecovery(const std::vector<int>& netIndices)
{
  if (!constants_.enable_wirelength_recovery || netIndices.empty()) {
    return;
  }

  struct Candidate
  {
    int index;
    int hpwl;
  };

  auto measureRoute = [&](const std::shared_ptr<GRTreeNode>& tree) {
    std::pair<uint64_t, int> stats{0, 0};
    if (!tree) {
      return stats;
    }
    GRTreeNode::preorder(
        tree, [&](const std::shared_ptr<GRTreeNode>& node) {
          for (const auto& child : node->getChildren()) {
            if (node->getLayerIdx() == child->getLayerIdx()) {
              const int direction
                  = grid_graph_->getLayerDirection(node->getLayerIdx());
              const int l = std::min((*node)[direction], (*child)[direction]);
              const int h = std::max((*node)[direction], (*child)[direction]);
              for (int c = l; c < h; c++) {
                stats.first += grid_graph_->getEdgeLength(direction, c);
              }
            } else {
              stats.second += abs(node->getLayerIdx() - child->getLayerIdx());
            }
          }
        });
    return stats;
  };

  auto isImprovement = [](const std::pair<uint64_t, int>& candidate_stats,
                          const std::pair<uint64_t, int>& baseline_stats) {
    return candidate_stats.first < baseline_stats.first
           || (candidate_stats.first == baseline_stats.first
               && candidate_stats.second < baseline_stats.second);
  };

  int total_accepted = 0;
  const int max_passes = std::max(1, constants_.recovery_max_passes);
  const double pass_decay = std::clamp(constants_.recovery_pass_decay, 0.25, 1.0);
  for (int pass = 0; pass < max_passes; pass++) {
    std::vector<Candidate> candidates;
    candidates.reserve(netIndices.size());
    const double pass_scale = std::pow(pass_decay, pass);
    const int hpwl_threshold = std::max(
        24, static_cast<int>(std::ceil(constants_.recovery_hpwl_threshold * pass_scale)));
    for (const int netIndex : netIndices) {
      const auto& net = gr_nets_[netIndex];
      const auto& tree = net->getRoutingTree();
      if (!tree) {
        continue;
      }
      const int hpwl = net->getBoundingBox().hp();
      if (hpwl < hpwl_threshold) {
        continue;
      }
      if (grid_graph_->checkOverflow(tree) > 0) {
        continue;
      }
      candidates.push_back({netIndex, hpwl});
    }

    if (candidates.empty()) {
      break;
    }

    std::sort(candidates.begin(),
              candidates.end(),
              [](const Candidate& lhs, const Candidate& rhs) {
                return lhs.hpwl > rhs.hpwl;
              });

    const double pass_ratio = std::clamp(
        constants_.recovery_refine_ratio
            * (pass == 0 ? 1.0 : std::max(0.45, pass_scale)),
        0.2,
        1.0);
    int keep
        = static_cast<int>(std::ceil(candidates.size() * pass_ratio));
    keep = std::max(1, std::min(keep, static_cast<int>(candidates.size())));
    int deep_keep = static_cast<int>(std::ceil(
        keep * std::clamp(constants_.recovery_deep_ratio, 0.0, 1.0)));
    deep_keep = std::max(1, std::min(deep_keep, keep));

    logger_->report("stage 4.{}: wirelength recovery on {} / {} nets (deep "
                    "search on {} nets)",
                    pass + 1,
                    keep,
                    candidates.size(),
                    deep_keep);

    int accepted_in_pass = 0;
    for (int candidateIndex = 0; candidateIndex < keep; candidateIndex++) {
      const int netIndex = candidates[candidateIndex].index;
      GRNet* net = gr_nets_[netIndex].get();
      const std::shared_ptr<GRTreeNode> original_tree = net->getRoutingTree();
      if (!original_tree) {
        continue;
      }
      const bool deep_search = candidateIndex < deep_keep
                               || candidates[candidateIndex].hpwl
                                      >= constants_.recovery_deep_hpwl_threshold
                               || pass > 0;
      const auto original_stats = measureRoute(original_tree);
      const int original_tree_overflow = grid_graph_->checkOverflow(original_tree);
      const int max_via_allowed
          = original_stats.second + std::max(0, constants_.recovery_max_via_increase);

      grid_graph_->commitTree(original_tree, /*ripup*/ true);

      std::shared_ptr<GRTreeNode> best_tree = nullptr;
      std::pair<uint64_t, int> best_stats = original_stats;

      auto considerCandidate = [&](const std::shared_ptr<GRTreeNode>& tree) {
        if (!tree) {
          return;
        }
        constexpr int kRecoveryOverflowSlack = 2;
        // Evaluate overflow after adding the candidate tree back to the live
        // graph. This avoids selecting routes that appear legal in rip-up mode
        // but create new overflows once committed.
        grid_graph_->commitTree(tree);
        const int overflow_after_commit = grid_graph_->checkOverflow(tree);
        grid_graph_->commitTree(tree, /*ripup*/ true);
        if (overflow_after_commit
            > original_tree_overflow + kRecoveryOverflowSlack) {
          return;
        }
        const auto stats = measureRoute(tree);
        if (stats.second > max_via_allowed) {
          return;
        }
        if (!best_tree || isImprovement(stats, best_stats)) {
          best_tree = tree;
          best_stats = stats;
        }
      };

      // Candidate A: CUGR DAG reroute.
      PatternRoute patternRoute(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      patternRoute.constructSteinerTree();
      patternRoute.constructRoutingDAG();
      patternRoute.run();
      considerCandidate(net->getRoutingTree());

      // Candidate B: Sparsified maze reroute with FastRoute-style shifted
      // coarse grids and SPRoute-style dense local search.
      if (constants_.recovery_use_maze) {
        GridGraphView<CostT> recoveryWireCostView;
        grid_graph_->extractWireCostView(recoveryWireCostView);

        struct MazeConfig
        {
          int sparse_x;
          int sparse_y;
          int offset_x;
          int offset_y;
        };

        std::vector<MazeConfig> maze_configs;
        maze_configs.reserve(64);
        const int max_maze_configs
            = std::max(2, constants_.recovery_max_maze_configs);
        auto addMazeConfig = [&](int sparse_x,
                                 int sparse_y,
                                 int offset_x,
                                 int offset_y) {
          if (static_cast<int>(maze_configs.size()) >= max_maze_configs) {
            return;
          }
          sparse_x = std::max(2, sparse_x);
          sparse_y = std::max(2, sparse_y);
          offset_x = std::clamp(offset_x, 0, sparse_x - 1);
          offset_y = std::clamp(offset_y, 0, sparse_y - 1);
          for (const auto& cfg : maze_configs) {
            if (cfg.sparse_x == sparse_x && cfg.sparse_y == sparse_y
                && cfg.offset_x == offset_x && cfg.offset_y == offset_y) {
              return;
            }
          }
          maze_configs.push_back({sparse_x, sparse_y, offset_x, offset_y});
        };

        const int base_sparse_x = std::max(2, constants_.recovery_maze_sparse_x);
        const int base_sparse_y = std::max(2, constants_.recovery_maze_sparse_y);
        addMazeConfig(base_sparse_x, base_sparse_y, 0, 0);
        addMazeConfig(
            constants_.recovery_dense_sparse_x, constants_.recovery_dense_sparse_y, 0, 0);

        if (deep_search) {
          const auto& bbox = net->getBoundingBox();
          const bool is_wide = bbox.width() >= bbox.height();
          const int long_sparse
              = std::max(2, constants_.recovery_aniso_sparse_long);
          const int short_sparse = std::max(
              2, std::min(long_sparse, constants_.recovery_aniso_sparse_short));

          if (is_wide) {
            addMazeConfig(long_sparse, short_sparse, 0, 0);
          } else {
            addMazeConfig(short_sparse, long_sparse, 0, 0);
          }

          addMazeConfig(std::max(2, base_sparse_x - 1),
                        std::max(2, base_sparse_y - 1),
                        0,
                        0);
          addMazeConfig(base_sparse_x + 1, base_sparse_y + 1, 0, 0);

          if (constants_.recovery_try_offset) {
            const int n_base_configs
                = std::min(4, static_cast<int>(maze_configs.size()));
            for (int cfg_idx = 0; cfg_idx < n_base_configs; cfg_idx++) {
              const auto cfg = maze_configs[cfg_idx];
              addMazeConfig(
                  cfg.sparse_x, cfg.sparse_y, cfg.sparse_x / 2, cfg.sparse_y / 2);
            }
          }

          const bool enable_full_offset_sweep
              = constants_.recovery_full_offset_sweep
                && candidates[candidateIndex].hpwl
                       >= constants_.recovery_full_offset_hpwl_threshold;
          if (enable_full_offset_sweep) {
            for (int offset_x = 0;
                 offset_x < base_sparse_x
                 && static_cast<int>(maze_configs.size()) < max_maze_configs;
                 offset_x++) {
              for (int offset_y = 0;
                   offset_y < base_sparse_y
                   && static_cast<int>(maze_configs.size()) < max_maze_configs;
                   offset_y++) {
                addMazeConfig(base_sparse_x, base_sparse_y, offset_x, offset_y);
              }
            }
            const int dense_sparse_x
                = std::max(2, constants_.recovery_dense_sparse_x);
            const int dense_sparse_y
                = std::max(2, constants_.recovery_dense_sparse_y);
            for (int offset_x = 0;
                 offset_x < dense_sparse_x
                 && static_cast<int>(maze_configs.size()) < max_maze_configs;
                 offset_x++) {
              for (int offset_y = 0;
                   offset_y < dense_sparse_y
                   && static_cast<int>(maze_configs.size()) < max_maze_configs;
                   offset_y++) {
                addMazeConfig(dense_sparse_x, dense_sparse_y, offset_x, offset_y);
              }
            }
          }
        }

        for (const auto& cfg : maze_configs) {
          MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
          SparseGrid recoveryGrid(
              cfg.sparse_x, cfg.sparse_y, cfg.offset_x, cfg.offset_y);
          mazeRoute.constructSparsifiedGraph(recoveryWireCostView, recoveryGrid);
          mazeRoute.run();
          if (const std::shared_ptr<SteinerTreeNode> maze_tree
              = mazeRoute.getSteinerTree()) {
            PatternRoute mazePatternRoute(
                net, grid_graph_.get(), stt_builder_, constants_, logger_);
            mazePatternRoute.setSteinerTree(maze_tree);
            mazePatternRoute.constructRoutingDAG();
            mazePatternRoute.run();
            considerCandidate(net->getRoutingTree());
          }
        }

        // FastRoute-style wirelength-first recovery: rerun sparse maze using
        // pure wirelength edge costs and keep only legal improvements.
        const bool enable_wl_only_maze
            = constants_.recovery_use_wirelength_maze
              && (deep_search || candidateIndex < keep / 2)
              && candidates[candidateIndex].hpwl
                     >= constants_.recovery_wl_only_hpwl_threshold;
        if (enable_wl_only_maze) {
          GridGraphView<CostT> recoveryWlOnlyView;
          grid_graph_->extractWireLengthCostView(recoveryWlOnlyView);
          const double wl_via_cost_scale
              = std::clamp(constants_.recovery_wl_maze_via_cost_scale, 0.0, 1.0);
          const int wl_config_count = std::min(static_cast<int>(maze_configs.size()),
                                               deep_search ? 10 : 4);
          for (int cfg_index = 0; cfg_index < wl_config_count; cfg_index++) {
            const auto& cfg = maze_configs[cfg_index];
            MazeRoute wlMazeRoute(net, grid_graph_.get(), logger_);
            SparseGrid wlGrid(
                cfg.sparse_x, cfg.sparse_y, cfg.offset_x, cfg.offset_y);
            wlMazeRoute.constructSparsifiedGraph(
                recoveryWlOnlyView, wlGrid, wl_via_cost_scale);
            wlMazeRoute.run();
            if (const std::shared_ptr<SteinerTreeNode> wl_tree
                = wlMazeRoute.getSteinerTree()) {
              PatternRoute wlPatternRoute(
                  net, grid_graph_.get(), stt_builder_, constants_, logger_);
              wlPatternRoute.setSteinerTree(wl_tree);
              wlPatternRoute.constructRoutingDAG();
              wlPatternRoute.run();
              considerCandidate(net->getRoutingTree());
            }
          }
        }
      }

      const bool accept = best_tree && isImprovement(best_stats, original_stats);
      if (accept) {
        net->setRoutingTree(best_tree);
        grid_graph_->commitTree(best_tree);
        accepted_in_pass++;
      } else {
        net->setRoutingTree(original_tree);
        grid_graph_->commitTree(original_tree);
      }
    }

    total_accepted += accepted_in_pass;
    logger_->report("stage 4.{} accepted {} reroutes.",
                    pass + 1,
                    accepted_in_pass);
    if (accepted_in_pass == 0) {
      break;
    }
  }

  logger_->report("wirelength recovery accepted {} reroutes.",
                  total_accepted);
}

void CUGR::route()
{
  std::vector<int> allNetIndices;
  allNetIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    allNetIndices.push_back(net->getIndex());
  }

  std::vector<int> netIndices = allNetIndices;
  patternRoute(netIndices);

  std::vector<int> detourIndices = netIndices;
  if (constants_.wirelength_first_refinement) {
    detourIndices
        = selectCriticalNets(netIndices, constants_.detour_refine_ratio);
  }
  patternRouteWithDetours(detourIndices);

  std::vector<int> mazeIndices = detourIndices;
  if (constants_.wirelength_first_refinement) {
    mazeIndices
        = selectCriticalNets(detourIndices, constants_.maze_refine_ratio);
  }
  mazeRoute(mazeIndices);
  wirelengthRecovery(allNetIndices);

  printStatistics();
  if (constants_.write_heatmap) {
    grid_graph_->write();
  }
}

std::vector<int> CUGR::selectCriticalNets(
    const std::vector<int>& candidates,
    const double reroute_ratio) const
{
  if (candidates.empty()) {
    return {};
  }

  struct ScoredNet
  {
    int index;
    int overflow;
    int hpwl;
    bool critical;
  };

  std::vector<ScoredNet> scored;
  scored.reserve(candidates.size());
  for (const int netIndex : candidates) {
    const auto& net = gr_nets_[netIndex];
    const int overflow = grid_graph_->checkOverflow(net->getRoutingTree());
    const int hpwl = net->getBoundingBox().hp();
    const bool critical
        = overflow >= constants_.refinement_overflow_threshold
          || hpwl >= constants_.refinement_hpwl_threshold;
    scored.push_back({netIndex, overflow, hpwl, critical});
  }

  std::sort(scored.begin(),
            scored.end(),
            [](const ScoredNet& lhs, const ScoredNet& rhs) {
              if (lhs.critical != rhs.critical) {
                return lhs.critical > rhs.critical;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              return lhs.hpwl > rhs.hpwl;
            });

  int keep = static_cast<int>(std::ceil(scored.size() * reroute_ratio));
  keep = std::max(1, std::min(keep, static_cast<int>(scored.size())));

  // Always keep all critical nets to avoid starvation.
  int critical_count = 0;
  while (critical_count < static_cast<int>(scored.size())
         && scored[critical_count].critical) {
    critical_count++;
  }
  keep = std::max(keep, critical_count);

  std::vector<int> selected;
  selected.reserve(keep);
  for (int i = 0; i < keep; i++) {
    selected.push_back(scored[i].index);
  }
  logger_->report("wirelength-first refinement: {} -> {} nets",
                  candidates.size(),
                  selected.size());
  return selected;
}

void CUGR::write(const std::string& guide_file)
{
  area_of_pin_patches_ = 0;
  area_of_wire_patches_ = 0;
  std::stringstream ss;
  for (const auto& net : gr_nets_) {
    std::vector<std::pair<int, BoxT>> guides;
    getGuides(net.get(), guides);

    ss << net->getName() << '\n';
    ss << "(\n";
    for (const auto& guide : guides) {
      ss << grid_graph_->getGridline(0, guide.second.lx()) << " "
         << grid_graph_->getGridline(1, guide.second.ly()) << " "
         << grid_graph_->getGridline(0, guide.second.hx() + 1) << " "
         << grid_graph_->getGridline(1, guide.second.hy() + 1) << " "
         << grid_graph_->getLayerName(guide.first) << "\n";
    }
    ss << ")\n";
  }
  logger_->report("total area of pin access patches: {}", area_of_pin_patches_);
  logger_->report("total area of wire segment patches: {}",
                  area_of_wire_patches_);
  std::ofstream fout(guide_file);
  fout << ss.str();
  fout.close();
}

NetRouteMap CUGR::getRoutes()
{
  NetRouteMap routes;
  for (const auto& net : gr_nets_) {
    if (net->getNumPins() < 2) {
      continue;
    }
    odb::dbNet* db_net = net->getDbNet();
    GRoute& route = routes[db_net];

    const int half_gcell = design_->getGridlineSize() / 2;

    auto& routing_tree = net->getRoutingTree();
    if (!routing_tree) {
      continue;
    }
    GRTreeNode::preorder(
        routing_tree, [&](const std::shared_ptr<GRTreeNode>& node) {
          for (const auto& child : node->getChildren()) {
            if (node->getLayerIdx() == child->getLayerIdx()) {
              auto [min_x, max_x] = std::minmax({node->x(), child->x()});
              auto [min_y, max_y] = std::minmax({node->y(), child->y()});

              // convert to dbu
              min_x = grid_graph_->getGridline(0, min_x) + half_gcell;
              min_y = grid_graph_->getGridline(1, min_y) + half_gcell;
              max_x = grid_graph_->getGridline(0, max_x) + half_gcell;
              max_y = grid_graph_->getGridline(1, max_y) + half_gcell;

              route.emplace_back(min_x,
                                 min_y,
                                 node->getLayerIdx() + 1,
                                 max_x,
                                 max_y,
                                 child->getLayerIdx() + 1,
                                 false);
            } else {
              const auto [bottom_layer, top_layer]
                  = std::minmax({node->getLayerIdx(), child->getLayerIdx()});
              for (int layer_idx = bottom_layer; layer_idx < top_layer;
                   layer_idx++) {
                const int x
                    = grid_graph_->getGridline(0, node->x()) + half_gcell;
                const int y
                    = grid_graph_->getGridline(1, node->y()) + half_gcell;

                route.emplace_back(
                    x, y, layer_idx + 1, x, y, layer_idx + 2, true);
              }
            }
          }
        });
  }

  return routes;
}

void CUGR::sortNetIndices(std::vector<int>& netIndices) const
{
  std::vector<int> halfParameters(gr_nets_.size());
  for (int netIndex : netIndices) {
    auto& net = gr_nets_[netIndex];
    halfParameters[netIndex] = net->getBoundingBox().hp();
  }
  sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    return halfParameters[lhs] > halfParameters[rhs];
  });
}

void CUGR::getGuides(const GRNet* net,
                     std::vector<std::pair<int, BoxT>>& guides)
{
  auto& routingTree = net->getRoutingTree();
  if (!routingTree) {
    return;
  }
  // 0. Basic guides
  GRTreeNode::preorder(
      routingTree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            guides.emplace_back(node->getLayerIdx(),
                                BoxT(std::min(node->x(), child->x()),
                                     std::min(node->y(), child->y()),
                                     std::max(node->x(), child->x()),
                                     std::max(node->y(), child->y())));
          } else {
            const int maxLayerIndex
                = std::max(node->getLayerIdx(), child->getLayerIdx());
            for (int layerIdx
                 = std::min(node->getLayerIdx(), child->getLayerIdx());
                 layerIdx <= maxLayerIndex;
                 layerIdx++) {
              guides.emplace_back(layerIdx, BoxT(node->x(), node->y()));
            }
          }
        }
      });

  auto getSpareResource = [&](const GRPoint& point) {
    double resource = std::numeric_limits<double>::max();
    const int direction = grid_graph_->getLayerDirection(point.getLayerIdx());
    if (point[direction] + 1 < grid_graph_->getSize(direction)) {
      resource = std::min(
          resource,
          grid_graph_->getEdge(point.getLayerIdx(), point.x(), point.y())
              .getResource());
    }
    if (point[direction] > 0) {
      GRPoint lower = point;
      lower[direction] -= 1;
      resource = std::min(
          resource,
          grid_graph_->getEdge(lower.getLayerIdx(), point.x(), point.y())
              .getResource());
    }
    return resource;
  };

  // 1. Pin access patches
  assert(constants_.min_routing_layer + 1 < grid_graph_->getNumLayers());
  for (auto& gpts : net->getPinAccessPoints()) {
    for (auto& gpt : gpts) {
      if (gpt.getLayerIdx() < constants_.min_routing_layer) {
        int padding = 0;
        if (getSpareResource({constants_.min_routing_layer, gpt.x(), gpt.y()})
            < constants_.pin_patch_threshold) {
          padding = constants_.pin_patch_padding;
        }
        for (int layerIdx = gpt.getLayerIdx();
             layerIdx <= constants_.min_routing_layer + 1;
             layerIdx++) {
          guides.emplace_back(
              layerIdx,
              BoxT(std::max(gpt.x() - padding, 0),
                   std::max(gpt.y() - padding, 0),
                   std::min(gpt.x() + padding,
                            (int) grid_graph_->getSize(0) - 1),
                   std::min(gpt.y() + padding,
                            (int) grid_graph_->getSize(1) - 1)));
          area_of_pin_patches_ += (guides.back().second.x().range() + 1)
                                  * (guides.back().second.y().range() + 1);
        }
      }
    }
  }

  // 2. Wire segment patches
  GRTreeNode::preorder(
      routingTree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            double wire_patch_threshold = constants_.wire_patch_threshold;
            const int direction
                = grid_graph_->getLayerDirection(node->getLayerIdx());
            const int l = std::min((*node)[direction], (*child)[direction]);
            const int h = std::max((*node)[direction], (*child)[direction]);
            const int r = (*node)[1 - direction];
            for (int c = l; c <= h; c++) {
              bool patched = false;
              const GRPoint point = (direction == MetalLayer::H
                                         ? GRPoint(node->getLayerIdx(), c, r)
                                         : GRPoint(node->getLayerIdx(), r, c));
              if (getSpareResource(point) < wire_patch_threshold) {
                for (int layerIndex = node->getLayerIdx() - 1;
                     layerIndex <= node->getLayerIdx() + 1;
                     layerIndex += 2) {
                  if (layerIndex < constants_.min_routing_layer
                      || layerIndex >= grid_graph_->getNumLayers()) {
                    continue;
                  }
                  if (getSpareResource({layerIndex, point.x(), point.y()})
                      >= 1.0) {
                    guides.emplace_back(layerIndex, BoxT(point.x(), point.y()));
                    area_of_wire_patches_ += 1;
                    patched = true;
                  }
                }
              }
              if (patched) {
                wire_patch_threshold = constants_.wire_patch_threshold;
              } else {
                wire_patch_threshold *= constants_.wire_patch_inflation_rate;
              }
            }
          }
        }
      });
}

void CUGR::printStatistics() const
{
  logger_->report("routing statistics");

  // wire length and via count
  uint64_t wireLength = 0;
  int viaCount = 0;
  std::vector<std::vector<std::vector<int>>> wireUsage;
  wireUsage.assign(grid_graph_->getNumLayers(),
                   std::vector<std::vector<int>>(
                       grid_graph_->getSize(0),
                       std::vector<int>(grid_graph_->getSize(1), 0)));
  for (const auto& net : gr_nets_) {
    GRTreeNode::preorder(
        net->getRoutingTree(), [&](const std::shared_ptr<GRTreeNode>& node) {
          for (const auto& child : node->getChildren()) {
            if (node->getLayerIdx() == child->getLayerIdx()) {
              const int direction
                  = grid_graph_->getLayerDirection(node->getLayerIdx());
              const int l = std::min((*node)[direction], (*child)[direction]);
              const int h = std::max((*node)[direction], (*child)[direction]);
              const int r = (*node)[1 - direction];
              for (int c = l; c < h; c++) {
                wireLength += grid_graph_->getEdgeLength(direction, c);
                const int x = direction == MetalLayer::H ? c : r;
                const int y = direction == MetalLayer::H ? r : c;
                wireUsage[node->getLayerIdx()][x][y] += 1;
              }
            } else {
              viaCount += abs(node->getLayerIdx() - child->getLayerIdx());
            }
          }
        });
  }

  // resource
  CapacityT overflow = 0;

  CapacityT minResource = std::numeric_limits<CapacityT>::max();
  GRPoint bottleneck(-1, -1, -1);
  for (int layerIndex = constants_.min_routing_layer;
       layerIndex < grid_graph_->getNumLayers();
       layerIndex++) {
    const int direction = grid_graph_->getLayerDirection(layerIndex);
    for (int x = 0; x < grid_graph_->getSize(0) - 1 + direction; x++) {
      for (int y = 0; y < grid_graph_->getSize(1) - direction; y++) {
        const CapacityT resource
            = grid_graph_->getEdge(layerIndex, x, y).getResource();
        if (resource < minResource) {
          minResource = resource;
          bottleneck = {layerIndex, x, y};
        }
        const CapacityT usage = wireUsage[layerIndex][x][y];
        const CapacityT capacity
            = std::max(grid_graph_->getEdge(layerIndex, x, y).capacity, 0.0);
        if (usage > 0.0 && usage > capacity) {
          overflow += usage - capacity;
        }
      }
    }
  }

  logger_->report("wire length (metric):  {}",
                  wireLength / grid_graph_->getM2Pitch());
  logger_->report("total via count:       {}", viaCount);
  logger_->report("total wire overflow:   {}", (int) overflow);

  logger_->report("min resource: {}", minResource);
  logger_->report("bottleneck:   {}", bottleneck);
}

void CUGR::updateDbCongestion()
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::dbGCellGrid* db_gcell = block->getGCellGrid();
  if (db_gcell == nullptr) {
    db_gcell = odb::dbGCellGrid::create(block);
  } else {
    db_gcell->resetGrid();
  }

  const int x_corner_ = design_->getDieRegion().lx();
  const int y_corner_ = design_->getDieRegion().ly();
  const int x_size_ = grid_graph_->getXSize();
  const int y_size_ = grid_graph_->getYSize();
  const int gridline_size = design_->getGridlineSize();
  db_gcell->addGridPatternX(x_corner_, x_size_, gridline_size);
  db_gcell->addGridPatternY(y_corner_, y_size_, gridline_size);

  odb::dbTech* db_tech = db_->getTech();
  for (int layer = 0; layer < grid_graph_->getNumLayers(); layer++) {
    odb::dbTechLayer* db_layer = db_tech->findRoutingLayer(layer + 1);
    if (db_layer == nullptr) {
      continue;
    }

    for (int y = 0; y < y_size_; y++) {
      for (int x = 0; x < x_size_; x++) {
        const GraphEdge& edge = grid_graph_->getEdge(layer, x, y);
        db_gcell->setCapacity(db_layer, x, y, edge.capacity);
        db_gcell->setUsage(db_layer, x, y, edge.demand);
      }
    }
  }
}

}  // namespace grt::newgr
