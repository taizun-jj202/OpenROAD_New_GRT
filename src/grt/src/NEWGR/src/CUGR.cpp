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

struct TreeStats
{
  uint64_t wire_length{0};
  int via_count{0};
};

TreeStats getTreeStats(const std::shared_ptr<GRTreeNode>& tree,
                       const GridGraph* grid_graph)
{
  TreeStats stats;
  if (!tree) {
    return stats;
  }

  GRTreeNode::preorder(tree, [&](const std::shared_ptr<GRTreeNode>& node) {
    for (const auto& child : node->getChildren()) {
      if (node->getLayerIdx() == child->getLayerIdx()) {
        const int direction = grid_graph->getLayerDirection(node->getLayerIdx());
        const int l = std::min((*node)[direction], (*child)[direction]);
        const int h = std::max((*node)[direction], (*child)[direction]);
        for (int edge = l; edge < h; edge++) {
          stats.wire_length += grid_graph->getEdgeLength(direction, edge);
        }
      } else {
        stats.via_count += abs(node->getLayerIdx() - child->getLayerIdx());
      }
    }
  });

  return stats;
}

}  // namespace

CUGR::CUGR(odb::dbDatabase* db,
           utl::Logger* log,
           stt::SteinerTreeBuilder* stt_builder)
    : db_(db), logger_(log), stt_builder_(stt_builder)
{
  // Radical wirelength-first policy:
  // permit extra vias and reduce short-area pressure so long trunks stay direct.
  constants_.weight_wire_length = 2.4;
  constants_.weight_via_number = 1.4;
  constants_.weight_short_area = 150.0;
  constants_.cost_logistic_slope = 0.22;
  constants_.maze_logistic_slope = 0.20;
  constants_.max_detour_ratio = 0.04;
  constants_.target_detour_count = 4;
  constants_.via_multiplier = 0.7;
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

void CUGR::patternRouteSubset(const std::vector<int>& netIndices,
                              const char* stage_name)
{
  if (netIndices.empty()) {
    return;
  }

  logger_->report("{} ({}, wire_scale {:.3f}, maze_scale {:.3f})",
                  stage_name,
                  netIndices.size(),
                  grid_graph_->getWireCongestionScale(),
                  grid_graph_->getMazeCongestionScale());
  std::vector<int> ordered = netIndices;
  sortNetIndices(ordered);
  for (const int netIndex : ordered) {
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
  for (const int netIndex : netIndices) {
    grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree(),
                            /*ripup*/ true);
  }
  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  sortNetIndices(netIndices);
  // Coarser sparse grid to speed up maze rerouting.
  SparseGrid grid(16, 16, 0, 0);
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
    mazeRoute.constructSparsifiedGraph(wireCostView, grid);
    mazeRoute.run();
    std::shared_ptr<SteinerTreeNode> tree = mazeRoute.getSteinerTree();
    assert(tree != nullptr);

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.setSteinerTree(tree);
    patternRoute.constructRoutingDAG();
    patternRoute.run();

    grid_graph_->commitTree(net->getRoutingTree());
    grid_graph_->updateWireCostView(wireCostView, net->getRoutingTree());
    grid.step();
  }

  updateOverflowNets(netIndices);
}

void CUGR::wirelengthRefine()
{
  constexpr int kMaxRefineNets = 1400;
  constexpr int kViaMargin = 10;
  const uint64_t kMinWireImprovement = 0;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();
  // Collapse detours first, then keep only overflow-safe improvements.
  grid_graph_->setCongestionPenaltyScales(0.01, original_maze_scale);

  std::vector<int> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    if (net->getNumPins() >= 2 && net->getRoutingTree()) {
      candidates.push_back(net->getIndex());
    }
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(), candidates.end(), [&](const int lhs, const int rhs) {
    const GRNet* lhsNet = gr_nets_[lhs].get();
    const GRNet* rhsNet = gr_nets_[rhs].get();
    const int lhsHpwl = lhsNet->getBoundingBox().hp();
    const int rhsHpwl = rhsNet->getBoundingBox().hp();
    if (lhsHpwl != rhsHpwl) {
      return lhsHpwl > rhsHpwl;
    }
    if (lhsNet->getNumPins() != rhsNet->getNumPins()) {
      return lhsNet->getNumPins() > rhsNet->getNumPins();
    }
    return lhs < rhs;
  });

  if ((int) candidates.size() > kMaxRefineNets) {
    candidates.resize(kMaxRefineNets);
  }

  int accepted = 0;
  int attempted = 0;
  logger_->report("stage 4: wirelength rescue on {} long nets",
                  candidates.size());

  for (const int netIndex : candidates) {
    GRNet* net = gr_nets_[netIndex].get();
    const std::shared_ptr<GRTreeNode> oldTree = net->getRoutingTree();
    if (!oldTree) {
      continue;
    }

    attempted++;
    const TreeStats oldStats = getTreeStats(oldTree, grid_graph_.get());

    // Evaluate old overflow with existing committed demands.
    const int oldOverflow = grid_graph_->checkOverflow(oldTree);

    grid_graph_->commitTree(oldTree, /*rip_up*/ true);

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();

    const std::shared_ptr<GRTreeNode> candidateTree = net->getRoutingTree();
    if (!candidateTree) {
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      continue;
    }

    const TreeStats candidateStats = getTreeStats(candidateTree, grid_graph_.get());

    // Temporarily commit candidate to evaluate overflow on its used edges.
    grid_graph_->commitTree(candidateTree);
    const int candidateOverflow = grid_graph_->checkOverflow(candidateTree);

    const bool improveWire = candidateStats.wire_length + kMinWireImprovement
                             < oldStats.wire_length;
    const bool improveComposite
        = (candidateStats.wire_length <= oldStats.wire_length
           && candidateStats.via_count + kViaMargin < oldStats.via_count);
    const bool keepOverflow = candidateOverflow <= oldOverflow;

    if (keepOverflow && (improveWire || improveComposite)) {
      accepted++;
      continue;
    }

    // Reject candidate and restore previous tree/demand.
    grid_graph_->commitTree(candidateTree, /*rip_up*/ true);
    net->setRoutingTree(oldTree);
    grid_graph_->commitTree(oldTree);
  }

  logger_->report("wirelength rescue accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::mazeWirelengthCollapse()
{
  constexpr int kMaxCriticalNets = 480;
  constexpr int kViaGrowthLimit = 120;
  constexpr uint64_t kMinWireImprovement = 1;
  constexpr int kOverflowPriorityWeight = 900;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical phase: rebuild only critical long/overflowing nets on a very fine
  // sparse maze to aggressively collapse trunks while preserving overflow.
  grid_graph_->setCongestionPenaltyScales(0.008, 0.020);

  struct CriticalCandidate
  {
    int net_index;
    int hpwl;
    int overflow;
    int score;
  };

  std::vector<CriticalCandidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    if (net->getNumPins() >= 3 && net->getRoutingTree()) {
      const int hpwl = net->getBoundingBox().hp();
      const int overflow = grid_graph_->checkOverflow(net->getRoutingTree());
      const int score = hpwl + overflow * kOverflowPriorityWeight;
      candidates.push_back({net->getIndex(), hpwl, overflow, score});
    }
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(
      candidates.begin(),
      candidates.end(),
      [](const CriticalCandidate& lhs, const CriticalCandidate& rhs) {
        if (lhs.score != rhs.score) {
          return lhs.score > rhs.score;
        }
        if (lhs.overflow != rhs.overflow) {
          return lhs.overflow > rhs.overflow;
        }
        if (lhs.hpwl != rhs.hpwl) {
          return lhs.hpwl > rhs.hpwl;
        }
        return lhs.net_index < rhs.net_index;
      });

  if ((int) candidates.size() > kMaxCriticalNets) {
    candidates.resize(kMaxCriticalNets);
  }

  std::vector<int> orderedIndices;
  orderedIndices.reserve(candidates.size());
  for (const auto& candidate : candidates) {
    orderedIndices.push_back(candidate.net_index);
  }

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  // Finer than the default maze stage to allow radical topology movement.
  SparseGrid sparseGrid(4, 4, 0, 0);
  int accepted = 0;
  int attempted = 0;
  logger_->report("stage 5: critical-net spine rebuild on {} nets",
                  orderedIndices.size());

  for (const int netIndex : orderedIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const std::shared_ptr<GRTreeNode> oldTree = net->getRoutingTree();
    if (!oldTree) {
      sparseGrid.step();
      continue;
    }

    attempted++;
    const TreeStats oldStats = getTreeStats(oldTree, grid_graph_.get());
    const int oldOverflow = grid_graph_->checkOverflow(oldTree);

    grid_graph_->commitTree(oldTree, /*rip_up*/ true);
    grid_graph_->updateWireCostView(wireCostView, oldTree);

    MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
    mazeRoute.constructSparsifiedGraph(wireCostView, sparseGrid);
    mazeRoute.run();
    std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();

    if (!steinerTree) {
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      grid_graph_->updateWireCostView(wireCostView, oldTree);
      sparseGrid.step();
      continue;
    }

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.setSteinerTree(steinerTree);
    patternRoute.constructRoutingDAG();
    patternRoute.run();

    const std::shared_ptr<GRTreeNode> candidateTree = net->getRoutingTree();
    if (!candidateTree) {
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      grid_graph_->updateWireCostView(wireCostView, oldTree);
      sparseGrid.step();
      continue;
    }

    const TreeStats candidateStats = getTreeStats(candidateTree, grid_graph_.get());

    grid_graph_->commitTree(candidateTree);
    grid_graph_->updateWireCostView(wireCostView, candidateTree);
    const int candidateOverflow = grid_graph_->checkOverflow(candidateTree);

    const bool improveWire = candidateStats.wire_length + kMinWireImprovement
                             < oldStats.wire_length;
    const bool viaGrowthBound
        = candidateStats.via_count <= oldStats.via_count + kViaGrowthLimit;
    const bool keepOverflow = candidateOverflow <= oldOverflow;
    const bool relieveOverflowWithoutWireRegression
        = (candidateOverflow + 2 < oldOverflow
           && candidateStats.wire_length <= oldStats.wire_length);

    if (viaGrowthBound
        && ((improveWire && keepOverflow) || relieveOverflowWithoutWireRegression)) {
      accepted++;
    } else {
      grid_graph_->commitTree(candidateTree, /*rip_up*/ true);
      grid_graph_->updateWireCostView(wireCostView, candidateTree);
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      grid_graph_->updateWireCostView(wireCostView, oldTree);
    }

    sparseGrid.step();
  }

  logger_->report("critical-net spine rebuild accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::hybridTopologySurgery()
{
  constexpr int kMaxSurgeryNets = 180;
  constexpr int kViaGrowthLimit = 96;
  constexpr uint64_t kMinWireImprovement = 8;
  constexpr uint64_t kOverflowWeight = 2200000;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical policy: dual-route each critical net with two distinct engines
  // (pattern-first vs maze-first) and keep only overflow-safe wirelength wins.
  grid_graph_->setCongestionPenaltyScales(0.010, 0.018);

  struct SurgeryCandidate
  {
    int net_index;
    uint64_t score;
    int overflow;
    int hpwl;
  };

  std::vector<SurgeryCandidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const auto routingTree = net->getRoutingTree();
    if (net->getNumPins() < 3 || !routingTree) {
      continue;
    }
    const TreeStats stats = getTreeStats(routingTree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routingTree);
    const int hpwl = net->getBoundingBox().hp();
    const uint64_t score
        = stats.wire_length + (uint64_t) overflow * kOverflowWeight;
    candidates.push_back({net->getIndex(), score, overflow, hpwl});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(
      candidates.begin(),
      candidates.end(),
      [](const SurgeryCandidate& lhs, const SurgeryCandidate& rhs) {
        if (lhs.score != rhs.score) {
          return lhs.score > rhs.score;
        }
        if (lhs.overflow != rhs.overflow) {
          return lhs.overflow > rhs.overflow;
        }
        if (lhs.hpwl != rhs.hpwl) {
          return lhs.hpwl > rhs.hpwl;
        }
        return lhs.net_index < rhs.net_index;
      });

  if ((int) candidates.size() > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  SparseGrid sparseGrid(6, 6, 0, 0);

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 6: hybrid topology surgery on {} critical nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> oldTree = net->getRoutingTree();
    if (!oldTree) {
      sparseGrid.step();
      continue;
    }

    attempted++;
    const TreeStats oldStats = getTreeStats(oldTree, grid_graph_.get());
    const int oldOverflow = grid_graph_->checkOverflow(oldTree);

    auto evaluateTree = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      grid_graph_->updateWireCostView(wireCostView, tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      grid_graph_->updateWireCostView(wireCostView, tree);
      return result;
    };

    grid_graph_->commitTree(oldTree, /*rip_up*/ true);
    grid_graph_->updateWireCostView(wireCostView, oldTree);

    // Candidate A: direct pattern rebuild.
    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();
    RerouteResult patternResult = evaluateTree(net->getRoutingTree());

    // Candidate B: maze-driven topology rebuild, then pattern legalization.
    RerouteResult mazeResult;
    MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
    mazeRoute.constructSparsifiedGraph(wireCostView, sparseGrid);
    mazeRoute.run();
    std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();
    if (steinerTree) {
      PatternRoute patternFromMaze(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      patternFromMaze.setSteinerTree(steinerTree);
      patternFromMaze.constructRoutingDAG();
      patternFromMaze.run();
      mazeResult = evaluateTree(net->getRoutingTree());
    }

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improveWire
          = result.stats.wire_length + kMinWireImprovement <= oldStats.wire_length;
      const bool keepOverflow = result.overflow <= oldOverflow;
      const bool boundedVia
          = result.stats.via_count <= oldStats.via_count + kViaGrowthLimit;
      return improveWire && keepOverflow && boundedVia;
    };

    const bool patternAcceptable = isAcceptable(patternResult);
    const bool mazeAcceptable = isAcceptable(mazeResult);

    auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
      if (lhs.stats.wire_length != rhs.stats.wire_length) {
        return lhs.stats.wire_length < rhs.stats.wire_length;
      }
      if (lhs.overflow != rhs.overflow) {
        return lhs.overflow < rhs.overflow;
      }
      return lhs.stats.via_count < rhs.stats.via_count;
    };

    const RerouteResult* best = nullptr;
    if (patternAcceptable) {
      best = &patternResult;
    }
    if (mazeAcceptable && (!best || isBetter(mazeResult, *best))) {
      best = &mazeResult;
    }

    if (best) {
      net->setRoutingTree(best->tree);
      grid_graph_->commitTree(best->tree);
      grid_graph_->updateWireCostView(wireCostView, best->tree);
      accepted++;
    } else {
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      grid_graph_->updateWireCostView(wireCostView, oldTree);
    }

    sparseGrid.step();
  }

  logger_->report("hybrid topology surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::route()
{
  constexpr int kMaxDetourNets = 2200;
  constexpr int kMaxMazeNets = 760;
  constexpr double kLongNetRatio = 0.72;

  auto trimOverflowSet = [&](std::vector<int>& netIndices,
                             const int limit,
                             const char* stage) {
    if ((int) netIndices.size() <= limit) {
      return;
    }

    struct OverflowStat
    {
      int net_index;
      int overflow;
      int hpwl;
    };

    std::vector<OverflowStat> overflowStats;
    overflowStats.reserve(netIndices.size());
    for (const int netIndex : netIndices) {
      const int overflow
          = grid_graph_->checkOverflow(gr_nets_[netIndex]->getRoutingTree());
      if (overflow <= 0) {
        continue;
      }
      overflowStats.push_back(
          {netIndex, overflow, gr_nets_[netIndex]->getBoundingBox().hp()});
    }

    if ((int) overflowStats.size() <= limit) {
      netIndices.clear();
      netIndices.reserve(overflowStats.size());
      for (const auto& stat : overflowStats) {
        netIndices.push_back(stat.net_index);
      }
      return;
    }

    std::sort(overflowStats.begin(),
              overflowStats.end(),
              [](const OverflowStat& lhs, const OverflowStat& rhs) {
                if (lhs.overflow != rhs.overflow) {
                  return lhs.overflow > rhs.overflow;
                }
                if (lhs.hpwl != rhs.hpwl) {
                  return lhs.hpwl > rhs.hpwl;
                }
                return lhs.net_index < rhs.net_index;
              });

    netIndices.clear();
    netIndices.reserve(limit);
    for (int i = 0; i < limit; i++) {
      netIndices.push_back(overflowStats[i].net_index);
    }

    logger_->report("runtime guard [{}]: rerouting {} / {} overflow nets",
                    stage,
                    netIndices.size(),
                    overflowStats.size());
  };

  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }
  if (netIndices.empty()) {
    return;
  }

  sortNetIndices(netIndices);
  int longNetCount = std::max(1, (int) (netIndices.size() * kLongNetRatio));
  longNetCount = std::min(longNetCount, (int) netIndices.size());
  std::vector<int> longNetSet(netIndices.begin(), netIndices.begin() + longNetCount);
  std::vector<int> remainingNetSet(netIndices.begin() + longNetCount,
                                   netIndices.end());

  grid_graph_->setCongestionPenaltyScales(0.04, 0.12);
  patternRouteSubset(longNetSet, "stage 1a: long-net trunk route");
  grid_graph_->setCongestionPenaltyScales(0.12, 0.16);
  patternRouteSubset(remainingNetSet, "stage 1b: congestion-aware fill route");
  updateOverflowNets(netIndices);

  grid_graph_->setCongestionPenaltyScales(0.22, 0.20);
  trimOverflowSet(netIndices, kMaxDetourNets, "detour");
  patternRouteWithDetours(netIndices);

  grid_graph_->setCongestionPenaltyScales(0.14, 0.26);
  trimOverflowSet(netIndices, kMaxMazeNets, "maze");
  mazeRoute(netIndices);

  hybridTopologySurgery();
  mazeWirelengthCollapse();
  grid_graph_->setCongestionPenaltyScales(0.05, 0.10);
  wirelengthRefine();

  printStatistics();
  if (constants_.write_heatmap) {
    grid_graph_->write();
  }
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
  std::sort(netIndices.begin(),
            netIndices.end(),
            [&](const int lhs, const int rhs) {
    const GRNet* lhsNet = gr_nets_[lhs].get();
    const GRNet* rhsNet = gr_nets_[rhs].get();
    const int lhsHpwl = lhsNet->getBoundingBox().hp();
    const int rhsHpwl = rhsNet->getBoundingBox().hp();
    if (lhsHpwl != rhsHpwl) {
      // Route long nets first so they lock in near-Manhattan trunks.
      return lhsHpwl > rhsHpwl;
    }
    if (lhsNet->getNumPins() != rhsNet->getNumPins()) {
      return lhsNet->getNumPins() > rhsNet->getNumPins();
    }
    return lhs < rhs;
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
