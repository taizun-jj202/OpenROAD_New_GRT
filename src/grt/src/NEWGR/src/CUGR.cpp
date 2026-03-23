#include "CUGR.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
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

struct RouteScore
{
  double wire_length = 0.0;
  int via_count = 0;
  int overflow_edges = 0;
  double objective = std::numeric_limits<double>::max();
};

RouteScore scoreRouteTree(const std::shared_ptr<GRTreeNode>& tree,
                          const GridGraph& graph,
                          const double wireWeight,
                          const double viaWeight,
                          const double overflowWeight)
{
  RouteScore score;
  if (!tree) {
    return score;
  }

  uint64_t wire_length_dbu = 0;
  GRTreeNode::preorder(tree, [&](const std::shared_ptr<GRTreeNode>& node) {
    for (const auto& child : node->getChildren()) {
      if (node->getLayerIdx() == child->getLayerIdx()) {
        const int direction = graph.getLayerDirection(node->getLayerIdx());
        if (direction == MetalLayer::H) {
          const auto [l, h] = std::minmax({node->x(), child->x()});
          for (int x = l; x < h; x++) {
            wire_length_dbu += graph.getEdgeLength(direction, x);
          }
        } else {
          const auto [l, h] = std::minmax({node->y(), child->y()});
          for (int y = l; y < h; y++) {
            wire_length_dbu += graph.getEdgeLength(direction, y);
          }
        }
      } else {
        score.via_count += std::abs(node->getLayerIdx() - child->getLayerIdx());
      }
    }
  });

  score.wire_length
      = static_cast<double>(wire_length_dbu) / std::max(graph.getM2Pitch(), 1);
  score.overflow_edges = graph.checkOverflow(tree);
  score.objective = wireWeight * score.wire_length + viaWeight * score.via_count
                    + overflowWeight * score.overflow_edges;
  return score;
}

RouteScore scoreRouteTree(const std::shared_ptr<GRTreeNode>& tree,
                          const GridGraph& graph,
                          const Constants& constants)
{
  return scoreRouteTree(tree,
                        graph,
                        constants.weight_wire_length,
                        constants.weight_via_number,
                        constants.weight_short_area);
}

struct CandidateGrid
{
  int interval;
  int x_offset;
  int y_offset;
};

void appendWireSegment(GRoute& route,
                       int x0,
                       int y0,
                       int layer,
                       int x1,
                       int y1,
                       uint64_t& totalRouteSegments)
{
  if (x0 == x1 && y0 == y1) {
    return;
  }
  route.emplace_back(x0, y0, layer, x1, y1, layer, false);
  totalRouteSegments++;
}

void appendViaStack(GRoute& route,
                    int x,
                    int y,
                    int fromLayer,
                    int toLayer,
                    uint64_t& totalRouteSegments)
{
  if (fromLayer == toLayer) {
    return;
  }
  const int lowLayer = std::min(fromLayer, toLayer);
  const int highLayer = std::max(fromLayer, toLayer);
  for (int layer = lowLayer; layer < highLayer; layer++) {
    route.emplace_back(x, y, layer, x, y, layer + 1, true);
    totalRouteSegments++;
  }
}

void pushGridCandidate(std::vector<CandidateGrid>& candidates,
                       int interval,
                       int xOffset,
                       int yOffset)
{
  if (interval <= 0) {
    return;
  }
  xOffset = ((xOffset % interval) + interval) % interval;
  yOffset = ((yOffset % interval) + interval) % interval;
  for (const auto& candidate : candidates) {
    if (candidate.interval == interval && candidate.x_offset == xOffset
        && candidate.y_offset == yOffset) {
      return;
    }
  }
  candidates.push_back({interval, xOffset, yOffset});
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
  grid_graph_->applySoftCapacities(gr_nets_);
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
  for (const int netIndex : netIndices) {
    grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree(),
                            /*ripup*/ true);
  }
  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  sortNetIndices(netIndices);
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const int hp = std::max(net->getBoundingBox().hp(), 1);
    const int maxShrink = std::max(constants_.maze_base_interval
                                       - constants_.maze_min_interval,
                                   0);
    const int adaptiveInterval = std::max(
        constants_.maze_min_interval,
        constants_.maze_base_interval - std::min(maxShrink, hp / 30));
    const int denseInterval
        = std::max(constants_.maze_min_interval, adaptiveInterval - 2);
    const int coarseInterval = adaptiveInterval + 2;

    std::vector<CandidateGrid> candidateGrids;
    pushGridCandidate(
        candidateGrids, adaptiveInterval, netIndex + hp, netIndex / 2 + hp);
    pushGridCandidate(
        candidateGrids, denseInterval, netIndex * 3 + hp, netIndex * 7 + hp * 2);
    pushGridCandidate(candidateGrids,
                      constants_.maze_min_interval,
                      netIndex * 11 + hp * 5,
                      netIndex * 13 + hp);
    pushGridCandidate(
        candidateGrids, coarseInterval, netIndex * 5 + hp * 3, netIndex + hp * 4);

    RouteScore bestScore;
    std::shared_ptr<GRTreeNode> bestTree = nullptr;
    for (const auto& candidate : candidateGrids) {
      SparseGrid grid(candidate.interval,
                      candidate.interval,
                      candidate.x_offset,
                      candidate.y_offset);
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

      const auto& candidateTree = net->getRoutingTree();
      const int pinCount = net->getNumPins();
      const double hp = std::max(net->getBoundingBox().hp(), 1);
      double wireWeight = constants_.weight_wire_length;
      double viaWeight = constants_.weight_via_number;
      double overflowWeight = constants_.weight_short_area;
      // Large nets dominate total wirelength; bias candidate selection to
      // compact geometry even at moderate via cost.
      if (pinCount >= 8) {
        wireWeight *= 1.45;
        viaWeight *= 1.35;
      }
      if (hp >= 120.0) {
        wireWeight *= 1.35;
        overflowWeight *= 1.20;
      }
      const RouteScore candidateScore = scoreRouteTree(candidateTree,
                                                       *grid_graph_,
                                                       wireWeight,
                                                       viaWeight,
                                                       overflowWeight);
      if (candidateScore.objective < bestScore.objective) {
        bestScore = candidateScore;
        bestTree = candidateTree;
      }
    }

    assert(bestTree != nullptr);
    net->setRoutingTree(bestTree);
    grid_graph_->commitTree(bestTree);
    grid_graph_->updateWireCostView(wireCostView, bestTree);
  }

  updateOverflowNets(netIndices);
}

void CUGR::wirelengthPulseRoute(const std::vector<int>& allNetIndices)
{
  if (!constants_.enable_wirelength_pulse_stage || allNetIndices.empty()
      || constants_.wirelength_pulse_rounds <= 0) {
    return;
  }

  const int totalNets = static_cast<int>(allNetIndices.size());
  const int selectedCount = std::clamp(
      static_cast<int>(std::round(constants_.wirelength_pulse_net_ratio
                                  * static_cast<double>(totalNets))),
      1,
      totalNets);

  logger_->report(
      "stage 4: wirelength pulse reroute ({} rounds, {} nets/round)",
      constants_.wirelength_pulse_rounds,
      selectedCount);

  for (int round = 0; round < constants_.wirelength_pulse_rounds; round++) {
    std::vector<std::pair<double, int>> rankedNets;
    rankedNets.reserve(allNetIndices.size());
    for (const int netIndex : allNetIndices) {
      const GRNet* net = gr_nets_[netIndex].get();
      const int hp = std::max(net->getBoundingBox().hp(), 1);
      const auto& tree = net->getRoutingTree();
      const RouteScore score = scoreRouteTree(tree, *grid_graph_, constants_);
      const double impact = score.wire_length + 1.8 * score.via_count
                            + 9.0 * hp + 180.0 * score.overflow_edges;
      rankedNets.emplace_back(impact, netIndex);
    }
    std::sort(rankedNets.begin(),
              rankedNets.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

    std::vector<int> rerouteIndices;
    rerouteIndices.reserve(selectedCount);
    const int shift = (round * std::max(selectedCount / 3, 1)) % totalNets;
    for (int i = 0; i < selectedCount; i++) {
      rerouteIndices.push_back(rankedNets[(shift + i) % totalNets].second);
    }

    std::sort(rerouteIndices.begin(), rerouteIndices.end(), [&](int lhs, int rhs) {
      const int lhsHp = gr_nets_[lhs]->getBoundingBox().hp();
      const int rhsHp = gr_nets_[rhs]->getBoundingBox().hp();
      if (lhsHp != rhsHp) {
        return lhsHp > rhsHp;
      }
      return gr_nets_[lhs]->getNumPins() > gr_nets_[rhs]->getNumPins();
    });

    for (const int netIndex : rerouteIndices) {
      grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree(), true);
    }

    GridGraphView<CostT> wireCostView;
    grid_graph_->extractWireCostView(wireCostView);
    int order = 0;
    for (const int netIndex : rerouteIndices) {
      GRNet* net = gr_nets_[netIndex].get();
      const int hp = std::max(net->getBoundingBox().hp(), 1);
      const int denseInterval
          = std::max(constants_.maze_min_interval,
                     constants_.wirelength_pulse_dense_interval);
      const int relaxedInterval = std::max(
          denseInterval + 1,
          std::max(constants_.maze_min_interval,
                   constants_.wirelength_pulse_relaxed_interval));
      std::vector<CandidateGrid> candidateGrids;
      pushGridCandidate(candidateGrids,
                        denseInterval,
                        netIndex * 5 + hp + round,
                        netIndex * 9 + order + hp);
      pushGridCandidate(candidateGrids,
                        denseInterval + 1,
                        netIndex * 11 + hp * 3,
                        netIndex * 13 + round);
      pushGridCandidate(candidateGrids,
                        relaxedInterval,
                        netIndex * 17 + order,
                        netIndex * 19 + hp);
      pushGridCandidate(candidateGrids,
                        relaxedInterval + 1,
                        netIndex * 23 + hp + round,
                        netIndex * 29 + order);

      RouteScore bestScore;
      std::shared_ptr<GRTreeNode> bestTree = nullptr;
      auto scoreCandidate = [&](const std::shared_ptr<GRTreeNode>& candidateTree) {
        double wireWeight = constants_.weight_wire_length * 2.3;
        double viaWeight = constants_.weight_via_number * 1.7;
        const double overflowWeight = constants_.weight_short_area * 2.3;
        if (hp >= 120 || net->getNumPins() >= 10) {
          wireWeight *= 1.10;
          viaWeight *= 1.35;
        }
        return scoreRouteTree(candidateTree,
                              *grid_graph_,
                              wireWeight,
                              viaWeight,
                              overflowWeight);
      };

      {
        PatternRoute patternRoute(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        patternRoute.constructSteinerTree();
        patternRoute.constructRoutingDAG();
        patternRoute.run();
        const auto& candidateTree = net->getRoutingTree();
        const RouteScore candidateScore = scoreCandidate(candidateTree);
        if (candidateScore.objective < bestScore.objective) {
          bestScore = candidateScore;
          bestTree = candidateTree;
        }
      }

      for (const auto& candidate : candidateGrids) {
        SparseGrid grid(candidate.interval,
                        candidate.interval,
                        candidate.x_offset,
                        candidate.y_offset);
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

        const auto& candidateTree = net->getRoutingTree();
        const RouteScore candidateScore = scoreCandidate(candidateTree);
        if (candidateScore.objective < bestScore.objective) {
          bestScore = candidateScore;
          bestTree = candidateTree;
        }
      }

      assert(bestTree != nullptr);
      net->setRoutingTree(bestTree);
      grid_graph_->commitTree(bestTree);
      grid_graph_->updateWireCostView(wireCostView, bestTree);
      order++;
    }

    std::vector<int> overflowIndices;
    updateOverflowNets(overflowIndices);
    logger_->report("wirelength pulse round {} complete, {} overflow nets remain",
                    round + 1,
                    overflowIndices.size());
  }
}

void CUGR::globalRebalanceRoute(const std::vector<int>& allNetIndices)
{
  if (allNetIndices.empty() || constants_.global_rebalance_rounds <= 0) {
    return;
  }

  logger_->report(
      "stage 5: global full-net annealed rebalance ({} rounds)",
      constants_.global_rebalance_rounds);
  std::vector<int> rerouteIndices = allNetIndices;
  for (int round = 0; round < constants_.global_rebalance_rounds; round++) {
    const bool longNetsFirst = (round % 2 == 0);
    sort(rerouteIndices.begin(), rerouteIndices.end(), [&](int lhs, int rhs) {
      const int lhsHp = gr_nets_[lhs]->getBoundingBox().hp();
      const int rhsHp = gr_nets_[rhs]->getBoundingBox().hp();
      if (lhsHp != rhsHp) {
        return longNetsFirst ? (lhsHp > rhsHp) : (lhsHp < rhsHp);
      }
      const int lhsPins = gr_nets_[lhs]->getNumPins();
      const int rhsPins = gr_nets_[rhs]->getNumPins();
      return longNetsFirst ? (lhsPins > rhsPins) : (lhsPins < rhsPins);
    });

    for (const int netIndex : rerouteIndices) {
      grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree(), true);
    }

    GridGraphView<CostT> wireCostView;
    grid_graph_->extractWireCostView(wireCostView);
    int order = 0;
    for (const int netIndex : rerouteIndices) {
      GRNet* net = gr_nets_[netIndex].get();
      const int hp = std::max(net->getBoundingBox().hp(), 1);
      const int maxShrink = std::max(constants_.maze_base_interval
                                         - constants_.maze_min_interval,
                                     0);
      const int adaptiveInterval = std::max(
          constants_.maze_min_interval,
          constants_.maze_base_interval - std::min(maxShrink, hp / 20));
      const int denseInterval = std::max(constants_.maze_min_interval,
                                         adaptiveInterval - 1);
      const int coarseInterval = adaptiveInterval + 1;

      std::vector<CandidateGrid> candidateGrids;
      pushGridCandidate(candidateGrids,
                        denseInterval,
                        round + order + hp,
                        netIndex * 3 + hp);
      pushGridCandidate(candidateGrids,
                        adaptiveInterval,
                        netIndex + hp * 5,
                        round * 11 + order);
      pushGridCandidate(candidateGrids,
                        coarseInterval,
                        netIndex * 7 + order,
                        hp * 13 + round);
      pushGridCandidate(candidateGrids,
                        constants_.maze_min_interval,
                        netIndex * 17 + hp,
                        round * 19 + netIndex);

      RouteScore bestScore;
      std::shared_ptr<GRTreeNode> bestTree = nullptr;
      for (const auto& candidate : candidateGrids) {
        SparseGrid grid(candidate.interval,
                        candidate.interval,
                        candidate.x_offset,
                        candidate.y_offset);
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

        double wireWeight = constants_.weight_wire_length
                            * (longNetsFirst ? 2.4 : 2.1);
        double viaWeight = constants_.weight_via_number * 1.25;
        if (net->getNumPins() >= 12 || hp >= 140) {
          wireWeight *= 1.15;
          viaWeight *= 1.35;
        }
        const double overflowWeight
            = constants_.weight_short_area * (1.3 + 0.25 * round);
        const auto& candidateTree = net->getRoutingTree();
        const RouteScore candidateScore = scoreRouteTree(candidateTree,
                                                         *grid_graph_,
                                                         wireWeight,
                                                         viaWeight,
                                                         overflowWeight);
        if (candidateScore.objective < bestScore.objective) {
          bestScore = candidateScore;
          bestTree = candidateTree;
        }
      }

      assert(bestTree != nullptr);
      net->setRoutingTree(bestTree);
      grid_graph_->commitTree(bestTree);
      grid_graph_->updateWireCostView(wireCostView, bestTree);
      order++;
    }

    std::vector<int> overflowIndices;
    updateOverflowNets(overflowIndices);
    logger_->report(
        "annealed rebalance round {} ({}) complete, {} overflow nets remain",
        round + 1,
        longNetsFirst ? "long-first" : "short-first",
        overflowIndices.size());
  }
}

void CUGR::criticalCompactionRoute(const std::vector<int>& allNetIndices)
{
  if (!constants_.enable_critical_compaction_stage || allNetIndices.empty()
      || constants_.critical_compaction_rounds <= 0) {
    return;
  }

  const int totalNets = static_cast<int>(allNetIndices.size());
  const int selectedCount = std::clamp(
      static_cast<int>(std::round(constants_.critical_compaction_net_ratio
                                  * static_cast<double>(totalNets))),
      1,
      totalNets);

  logger_->report(
      "stage 6: selective critical-net compaction ({} rounds, {} nets/round)",
      constants_.critical_compaction_rounds,
      selectedCount);
  for (int round = 0; round < constants_.critical_compaction_rounds; round++) {
    std::vector<std::pair<double, int>> rankedNets;
    rankedNets.reserve(allNetIndices.size());
    for (const int netIndex : allNetIndices) {
      const GRNet* net = gr_nets_[netIndex].get();
      const int hp = std::max(net->getBoundingBox().hp(), 1);
      const int pinDegree = std::max(net->getNumPins(), 2);
      // Emphasize long multi-pin nets, which dominate total wire/via impact.
      const double impact
          = static_cast<double>(hp) * pinDegree
            + 0.35 * static_cast<double>(hp) * static_cast<double>(hp);
      rankedNets.emplace_back(impact, netIndex);
    }
    std::sort(rankedNets.begin(),
              rankedNets.end(),
              [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

    std::vector<int> rerouteIndices;
    rerouteIndices.reserve(selectedCount);
    const int shift = (round * std::max(selectedCount / 2, 1)) % totalNets;
    for (int i = 0; i < selectedCount; i++) {
      rerouteIndices.push_back(rankedNets[(shift + i) % totalNets].second);
    }

    std::sort(rerouteIndices.begin(),
              rerouteIndices.end(),
              [&](int lhs, int rhs) {
                const int lhsHp = gr_nets_[lhs]->getBoundingBox().hp();
                const int rhsHp = gr_nets_[rhs]->getBoundingBox().hp();
                if (lhsHp != rhsHp) {
                  return (round % 2 == 0) ? (lhsHp > rhsHp) : (lhsHp < rhsHp);
                }
                return gr_nets_[lhs]->getNumPins() > gr_nets_[rhs]->getNumPins();
              });

    for (const int netIndex : rerouteIndices) {
      grid_graph_->commitTree(gr_nets_[netIndex]->getRoutingTree(), true);
    }

    GridGraphView<CostT> wireCostView;
    grid_graph_->extractWireCostView(wireCostView);
    int order = 0;
    for (const int netIndex : rerouteIndices) {
      GRNet* net = gr_nets_[netIndex].get();
      const int hp = std::max(net->getBoundingBox().hp(), 1);
      const int denseInterval = std::max(constants_.critical_compaction_interval
                                             + std::min(hp / 80, 2),
                                         1);
      std::vector<CandidateGrid> candidateGrids;
      pushGridCandidate(candidateGrids,
                        denseInterval,
                        round + order + hp,
                        netIndex + hp * 3);
      pushGridCandidate(
          candidateGrids, denseInterval + 1, netIndex + round, order + hp);
      pushGridCandidate(candidateGrids,
                        std::max(constants_.maze_min_interval, 1),
                        netIndex * 5 + hp,
                        netIndex * 11 + round);

      RouteScore bestScore;
      std::shared_ptr<GRTreeNode> bestTree = nullptr;
      for (const auto& candidate : candidateGrids) {
        SparseGrid grid(candidate.interval,
                        candidate.interval,
                        candidate.x_offset,
                        candidate.y_offset);
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

        const auto& candidateTree = net->getRoutingTree();
        const RouteScore candidateScore
            = scoreRouteTree(candidateTree,
                             *grid_graph_,
                             4.6,
                             4.4,
                             constants_.critical_compaction_overflow_weight);
        if (candidateScore.objective < bestScore.objective) {
          bestScore = candidateScore;
          bestTree = candidateTree;
        }
      }

      assert(bestTree != nullptr);
      net->setRoutingTree(bestTree);
      grid_graph_->commitTree(bestTree);
      grid_graph_->updateWireCostView(wireCostView, bestTree);
      order++;
    }

    std::vector<int> overflowIndices;
    updateOverflowNets(overflowIndices);
    logger_->report("critical compaction round {} complete, {} overflow nets remain",
                    round + 1,
                    overflowIndices.size());
  }
}

void CUGR::route()
{
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }
  const std::vector<int> allNetIndices = netIndices;

  patternRoute(netIndices);

  if (constants_.enable_early_detour_stage) {
    patternRouteWithDetours(netIndices);
  }

  mazeRoute(netIndices);
  wirelengthPulseRoute(allNetIndices);
  globalRebalanceRoute(allNetIndices);
  criticalCompactionRoute(allNetIndices);

  updateOverflowNets(netIndices);
  if (constants_.enable_early_detour_stage && !netIndices.empty()) {
    patternRouteWithDetours(netIndices);
    updateOverflowNets(netIndices);
  }
  logger_->report("completed pulse + rebalance + compaction stages");

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
  uint64_t totalRouteSegments = 0;
  for (const auto& net : gr_nets_) {
    if (net->getDbNet()->getTermCount() < 2) {
      continue;
    }
    odb::dbNet* db_net = net->getDbNet();
    GRoute& route = routes[db_net];

    const int half_gcell = design_->getGridlineSize() / 2;
    auto toDbuX = [&](int grid_x) {
      return grid_graph_->getGridline(0, grid_x) + half_gcell;
    };
    auto toDbuY = [&](int grid_y) {
      return grid_graph_->getGridline(1, grid_y) + half_gcell;
    };

    std::vector<GRPoint> pinAnchors;
    pinAnchors.reserve(net->getPinAccessPoints().size());
    for (const auto& pinAps : net->getPinAccessPoints()) {
      if (pinAps.empty()) {
        continue;
      }
      const GRPoint* best = &pinAps.front();
      for (const auto& ap : pinAps) {
        if (ap.getLayerIdx() < best->getLayerIdx()) {
          best = &ap;
        }
      }
      pinAnchors.push_back(*best);
    }

    auto& routing_tree = net->getRoutingTree();
    if (routing_tree) {
      GRTreeNode::preorder(
          routing_tree, [&](const std::shared_ptr<GRTreeNode>& node) {
            for (const auto& child : node->getChildren()) {
              if (node->getLayerIdx() == child->getLayerIdx()) {
                auto [min_x, max_x] = std::minmax({node->x(), child->x()});
                auto [min_y, max_y] = std::minmax({node->y(), child->y()});

                route.emplace_back(toDbuX(min_x),
                                   toDbuY(min_y),
                                   node->getLayerIdx() + 1,
                                   toDbuX(max_x),
                                   toDbuY(max_y),
                                   child->getLayerIdx() + 1,
                                   false);
                totalRouteSegments++;
              } else {
                const auto [bottom_layer, top_layer]
                    = std::minmax({node->getLayerIdx(), child->getLayerIdx()});
                const int x = toDbuX(node->x());
                const int y = toDbuY(node->y());
                for (int layer_idx = bottom_layer; layer_idx < top_layer;
                     layer_idx++) {
                  route.emplace_back(
                      x, y, layer_idx + 1, x, y, layer_idx + 2, true);
                  totalRouteSegments++;
                }
              }
            }
          });
    }

    // Ensure every multi-terminal net has at least one connected route pattern,
    // even when tree construction/routing fails for that net.
    if (route.empty() && pinAnchors.size() >= 2) {
      const GRPoint& hub = pinAnchors.front();
      const int hubX = toDbuX(hub.x());
      const int hubY = toDbuY(hub.y());
      const int hubLayer = hub.getLayerIdx() + 1;

      for (size_t i = 1; i < pinAnchors.size(); i++) {
        const GRPoint& pin = pinAnchors[i];
        const int pinX = toDbuX(pin.x());
        const int pinY = toDbuY(pin.y());
        const int pinLayer = pin.getLayerIdx() + 1;

        appendViaStack(
            route, pinX, pinY, pinLayer, hubLayer, totalRouteSegments);
        appendWireSegment(
            route, pinX, pinY, hubLayer, hubX, pinY, totalRouteSegments);
        appendWireSegment(
            route, hubX, pinY, hubLayer, hubX, hubY, totalRouteSegments);
      }
    }

    // Add per-pin anchor points so pins are explicitly represented in exported
    // segments for downstream parasitic/antenna stages.
    std::set<std::tuple<int, int, int>> emittedAnchors;
    for (const auto& pin : pinAnchors) {
      const int x = toDbuX(pin.x());
      const int y = toDbuY(pin.y());
      const int layer = pin.getLayerIdx() + 1;
      if (emittedAnchors.insert({x, y, layer}).second) {
        route.emplace_back(x, y, layer, x, y, layer, false);
        totalRouteSegments++;
      }
    }

    if (route.empty() && !pinAnchors.empty()) {
      const auto& pin = pinAnchors.front();
      route.emplace_back(toDbuX(pin.x()),
                         toDbuY(pin.y()),
                         pin.getLayerIdx() + 1,
                         toDbuX(pin.x()),
                         toDbuY(pin.y()),
                         pin.getLayerIdx() + 1,
                         false);
      totalRouteSegments++;
    }
  }

  logger_->report("exported {} global route segments", totalRouteSegments);

  return routes;
}

void CUGR::sortNetIndices(std::vector<int>& netIndices) const
{
  std::vector<int> halfParameters(gr_nets_.size(), 0);
  for (int netIndex : netIndices) {
    auto& net = gr_nets_[netIndex];
    halfParameters[netIndex] = net->getBoundingBox().hp();
  }
  sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (halfParameters[lhs] != halfParameters[rhs]) {
      return halfParameters[lhs] > halfParameters[rhs];
    }
    return gr_nets_[lhs]->getNumPins() > gr_nets_[rhs]->getNumPins();
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
