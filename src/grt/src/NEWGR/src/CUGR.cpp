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

struct RouteScore
{
  int overflow_edges = std::numeric_limits<int>::max();
  uint64_t wire_length = std::numeric_limits<uint64_t>::max();
  int via_count = std::numeric_limits<int>::max();
};

RouteScore evaluateRouteScore(const std::shared_ptr<GRTreeNode>& tree,
                              const GridGraph* grid_graph)
{
  if (!tree) {
    return {};
  }

  RouteScore score;
  score.overflow_edges = 0;
  score.wire_length = 0;
  score.via_count = 0;

  GRTreeNode::preorder(
      tree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            const int layer_idx = node->getLayerIdx();
            const int direction = grid_graph->getLayerDirection(layer_idx);
            const int l = std::min((*node)[direction], (*child)[direction]);
            const int h = std::max((*node)[direction], (*child)[direction]);
            const int r = (*node)[1 - direction];
            for (int c = l; c < h; c++) {
              score.wire_length += grid_graph->getEdgeLength(direction, c);
              const int x = direction == MetalLayer::H ? c : r;
              const int y = direction == MetalLayer::H ? r : c;
              if (grid_graph->checkOverflow(layer_idx, x, y)) {
                score.overflow_edges += 1;
              }
            }
          } else {
            score.via_count
                += std::abs(node->getLayerIdx() - child->getLayerIdx());
          }
        }
      });

  return score;
}

bool isBetterScore(const RouteScore& candidate, const RouteScore& baseline)
{
  if (baseline.overflow_edges == std::numeric_limits<int>::max()) {
    return true;
  }
  if (candidate.overflow_edges < baseline.overflow_edges) {
    const int overflowGain = baseline.overflow_edges - candidate.overflow_edges;
    // Keep overflow reduction, but tightly cap wire inflation. This borrows
    // FastRoute-style RRR acceptance pressure to prevent long detours from
    // being admitted early and then requiring expensive cleanup later.
    const uint64_t allowedIncrease
        = static_cast<uint64_t>(overflowGain) * 5ULL;
    const uint64_t allowedWireLength
        = baseline.wire_length + allowedIncrease;
    if (candidate.wire_length > allowedWireLength) {
      return false;
    }
    const int viaSlack = std::max(1, overflowGain / 7);
    return candidate.via_count <= baseline.via_count + viaSlack;
  }
  if (candidate.overflow_edges > baseline.overflow_edges) {
    return false;
  }
  if (candidate.wire_length != baseline.wire_length) {
    return candidate.wire_length < baseline.wire_length;
  }
  return candidate.via_count < baseline.via_count;
}

bool isRecoveryScoreBetter(const RouteScore& candidate,
                           const RouteScore& baseline)
{
  if (baseline.overflow_edges == std::numeric_limits<int>::max()) {
    return true;
  }
  if (candidate.overflow_edges > baseline.overflow_edges) {
    return false;
  }
  if (candidate.overflow_edges < baseline.overflow_edges) {
    const int overflowGain = baseline.overflow_edges - candidate.overflow_edges;
    const uint64_t allowedIncrease
        = static_cast<uint64_t>(overflowGain) * 3ULL;
    const uint64_t allowedWireLength
        = baseline.wire_length + allowedIncrease;
    if (candidate.wire_length > allowedWireLength) {
      return false;
    }
    const int viaSlack = std::max(1, overflowGain / 8);
    return candidate.via_count <= baseline.via_count + viaSlack;
  }
  if (candidate.wire_length < baseline.wire_length) {
    return true;
  }
  if (candidate.wire_length > baseline.wire_length) {
    return false;
  }
  return candidate.via_count < baseline.via_count;
}

bool isTightenScoreBetter(const RouteScore& candidate, const RouteScore& baseline)
{
  if (baseline.overflow_edges == std::numeric_limits<int>::max()) {
    return true;
  }
  if (candidate.overflow_edges > baseline.overflow_edges) {
    return false;
  }
  if (candidate.overflow_edges < baseline.overflow_edges) {
    const int overflowGain = baseline.overflow_edges - candidate.overflow_edges;
    // Tightening is wirelength-focused. Overflow reduction can grow wire only
    // slightly to avoid reintroducing long detours late in the flow.
    const uint64_t allowedIncrease
        = static_cast<uint64_t>(overflowGain) * 2ULL;
    const uint64_t allowedWireLength
        = baseline.wire_length + allowedIncrease;
    if (candidate.wire_length > allowedWireLength) {
      return false;
    }
    return candidate.via_count <= baseline.via_count + 1;
  }
  if (candidate.wire_length < baseline.wire_length) {
    return true;
  }
  if (candidate.wire_length > baseline.wire_length) {
    return false;
  }
  return candidate.via_count < baseline.via_count;
}

bool isCompactionScoreBetter(const RouteScore& candidate,
                             const RouteScore& baseline)
{
  if (baseline.overflow_edges == std::numeric_limits<int>::max()) {
    return true;
  }
  if (candidate.overflow_edges > baseline.overflow_edges) {
    return false;
  }
  if (candidate.overflow_edges < baseline.overflow_edges) {
    const int overflowGain = baseline.overflow_edges - candidate.overflow_edges;
    // Keep overflow relief but enforce near wirelength neutrality.
    const uint64_t allowedIncrease
        = std::max<uint64_t>(1ULL, static_cast<uint64_t>(overflowGain) / 3ULL);
    if (candidate.wire_length > baseline.wire_length + allowedIncrease) {
      return false;
    }
    return candidate.via_count <= baseline.via_count + 1;
  }
  if (candidate.wire_length < baseline.wire_length) {
    if (candidate.via_count <= baseline.via_count + 2) {
      return true;
    }
    const uint64_t wireGain = baseline.wire_length - candidate.wire_length;
    const int viaIncrease = candidate.via_count - baseline.via_count;
    return viaIncrease > 0
           && wireGain >= static_cast<uint64_t>(viaIncrease) * 10ULL;
  }
  if (candidate.wire_length > baseline.wire_length) {
    return false;
  }
  return candidate.via_count < baseline.via_count;
}

bool isStrictWirelengthScoreBetter(const RouteScore& candidate,
                                   const RouteScore& baseline)
{
  if (baseline.overflow_edges == std::numeric_limits<int>::max()) {
    return true;
  }
  if (candidate.overflow_edges > baseline.overflow_edges) {
    return false;
  }
  if (candidate.overflow_edges < baseline.overflow_edges) {
    const int overflowGain = baseline.overflow_edges - candidate.overflow_edges;
    // Preserve overflow relief while keeping the route close to wirelength
    // neutral.
    const uint64_t allowedIncrease
        = std::max<uint64_t>(1ULL, static_cast<uint64_t>(overflowGain) / 4ULL);
    if (candidate.wire_length > baseline.wire_length + allowedIncrease) {
      return false;
    }
    return candidate.via_count <= baseline.via_count;
  }
  // For overflow-neutral updates, require strict wirelength reduction.
  if (candidate.wire_length >= baseline.wire_length) {
    return false;
  }
  if (candidate.via_count <= baseline.via_count + 2) {
    return true;
  }
  const uint64_t wireGain = baseline.wire_length - candidate.wire_length;
  const int viaIncrease = candidate.via_count - baseline.via_count;
  return viaIncrease > 0
         && wireGain >= static_cast<uint64_t>(viaIncrease) * 24ULL;
}

std::vector<SparseGrid> buildMazeCandidateGrids(int base_interval,
                                                 int rank,
                                                 int hp,
                                                 int pins,
                                                 int max_candidates)
{
  const bool criticalLongNet
      = (pins >= 8 || hp >= 120) && rank < 2048;
  const bool extremeLongNet = (pins >= 12 || hp >= 180) && rank < 1536;
  std::vector<int> intervals{
      base_interval,
      std::max(3, base_interval - 1),
      std::max(3, base_interval - 2),
      std::min(12, base_interval + 1)};
  if (criticalLongNet) {
    // Add a dense sparse-grid option only for critical long nets to expose
    // shorter reconnection opportunities without applying the runtime hit to
    // all nets.
    intervals.emplace_back(2);
  }
  if (pins >= 12 || hp >= 160) {
    intervals.emplace_back(3);
  }
  if (extremeLongNet) {
    intervals.emplace_back(2);
    intervals.emplace_back(3);
  }
  if (pins <= 3 && hp <= 60) {
    intervals.emplace_back(std::min(12, base_interval + 2));
  }

  std::vector<SparseGrid> grids;
  grids.reserve(max_candidates);
  auto addGrid = [&](int interval, int x_offset, int y_offset) {
    interval = std::clamp(interval, 2, 12);
    x_offset = ((x_offset % interval) + interval) % interval;
    y_offset = ((y_offset % interval) + interval) % interval;
    for (const auto& grid : grids) {
      if (grid.interval.x() == interval && grid.interval.y() == interval
          && grid.offset.x() == x_offset && grid.offset.y() == y_offset) {
        return;
      }
    }
    if (grids.size() < static_cast<size_t>(max_candidates)) {
      grids.emplace_back(interval, interval, x_offset, y_offset);
    }
  };

  for (size_t idx = 0;
       idx < intervals.size()
       && grids.size() < static_cast<size_t>(max_candidates);
       idx++) {
    const int interval = std::clamp(intervals[idx], 2, 12);
    const int idx_int = static_cast<int>(idx);
    const int x_offset
        = (rank * (3 + idx_int * 2) + hp + idx_int * 5) % interval;
    const int y_offset
        = (rank * (5 + idx_int * 2) + pins + idx_int * 7) % interval;
    addGrid(interval, x_offset, y_offset);
    addGrid(interval, interval - 1 - x_offset, interval - 1 - y_offset);
  }
  if (grids.empty()) {
    grids.emplace_back(4, 4, 0, 0);
  }
  return grids;
}

std::vector<int> buildSpatialCompactionOrder(
    const std::vector<int>& prioritized_nets,
    const std::vector<std::unique_ptr<GRNet>>& gr_nets,
    int max_count,
    bool use_x_axis)
{
  if (prioritized_nets.empty() || max_count <= 0) {
    return {};
  }

  constexpr int bucket_count = 32;
  std::vector<std::vector<int>> buckets(bucket_count);
  for (const int netIndex : prioritized_nets) {
    const auto& box = gr_nets[netIndex]->getBoundingBox();
    int coordinate = use_x_axis ? box.cx() : box.cy();
    int bucket = coordinate % bucket_count;
    if (bucket < 0) {
      bucket += bucket_count;
    }
    buckets[bucket].push_back(netIndex);
  }

  std::vector<int> order;
  order.reserve(std::min(max_count, static_cast<int>(prioritized_nets.size())));
  for (size_t offset = 0; order.size() < static_cast<size_t>(max_count);
       offset++) {
    bool inserted = false;
    for (int bucket = 0; bucket < bucket_count; bucket++) {
      if (offset < buckets[bucket].size()) {
        order.push_back(buckets[bucket][offset]);
        inserted = true;
        if (order.size() == static_cast<size_t>(max_count)) {
          break;
        }
      }
    }
    if (!inserted) {
      break;
    }
  }
  return order;
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
  std::vector<int> overflowEdges(gr_nets_.size(), 0);
  for (const int netIndex : netIndices) {
    overflowEdges[netIndex]
        = grid_graph_->checkOverflow(gr_nets_[netIndex]->getRoutingTree());
  }
  std::stable_sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (overflowEdges[lhs] != overflowEdges[rhs]) {
      return overflowEdges[lhs] > overflowEdges[rhs];
    }
    return lhs < rhs;
  });
  const int detourBudget = std::max(128, static_cast<int>(netIndices.size() / 4));
  const int minOverflowForDetour = 2;
  if (overflowEdges[netIndices.front()] < minOverflowForDetour) {
    logger_->report("stage 2 skipped: overflow is below detour trigger.");
    updateOverflowNets(netIndices);
    return;
  }
  int considered = 0;
  int accepted = 0;
  for (int rank = 0; rank < netIndices.size(); rank++) {
    const int netIndex = netIndices[rank];
    if (rank >= detourBudget || overflowEdges[netIndex] < minOverflowForDetour) {
      continue;
    }
    considered++;
    GRNet* net = gr_nets_[netIndex].get();
    const auto oldTree = net->getRoutingTree();
    grid_graph_->commitTree(oldTree, /*ripup*/ true);

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    // KEY DIFFERENCE compared to stage 1 (patternRoute)
    patternRoute.constructDetours(congestionView);
    patternRoute.run();
    const auto candidateTree = net->getRoutingTree();

    if (!oldTree || !candidateTree) {
      if (candidateTree) {
        grid_graph_->commitTree(candidateTree);
      } else if (oldTree) {
        net->setRoutingTree(oldTree);
        grid_graph_->commitTree(oldTree);
      }
      continue;
    }

    grid_graph_->commitTree(oldTree);
    const RouteScore oldScore = evaluateRouteScore(oldTree, grid_graph_.get());
    grid_graph_->commitTree(oldTree, /*ripup*/ true);

    grid_graph_->commitTree(candidateTree);
    const RouteScore candidateScore
        = evaluateRouteScore(candidateTree, grid_graph_.get());

    if (isBetterScore(candidateScore, oldScore)) {
      accepted++;
      continue;
    }

    grid_graph_->commitTree(candidateTree, /*ripup*/ true);
    net->setRoutingTree(oldTree);
    grid_graph_->commitTree(oldTree);
  }
  logger_->report("stage 2 detour candidates accepted {} / {} nets.",
                  accepted,
                  considered);

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
  const int base_sparse_interval = netIndices.size() < 2000 ? 6 : 10;
  int rank = 0;
  int accepted = 0;
  int totalCandidates = 0;
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const auto oldTree = net->getRoutingTree();

    // Mix FastRoute/SPRoute ideas: use denser sparse graphs on long/high-fanout
    // nets and coarser sampling on simpler nets to control runtime.
    const int hp = net->getBoundingBox().hp();
    const int pins = net->getNumPins();
    int interval = base_sparse_interval;
    if (pins >= 16 || hp >= 180) {
      interval = std::max(3, base_sparse_interval - 3);
    } else if (pins >= 8 || hp >= 110) {
      interval = std::max(4, base_sparse_interval - 2);
    } else if (pins >= 4 || hp >= 70) {
      interval = std::max(5, base_sparse_interval - 1);
    }
    const int max_candidates = (pins >= 10 || hp >= 130) ? 5 : 4;
    const auto candidateGrids
        = buildMazeCandidateGrids(interval, rank, hp, pins, max_candidates);
    totalCandidates += static_cast<int>(candidateGrids.size());

    std::shared_ptr<GRTreeNode> bestTree = oldTree;
    RouteScore bestScore;
    if (oldTree) {
      grid_graph_->commitTree(oldTree);
      bestScore = evaluateRouteScore(oldTree, grid_graph_.get());
      grid_graph_->commitTree(oldTree, /*ripup*/ true);
    } else {
      bestScore = RouteScore{};
    }

    for (const auto& grid : candidateGrids) {
      MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
      mazeRoute.constructSparsifiedGraph(wireCostView, grid);
      mazeRoute.run();
      std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();
      if (!steinerTree) {
        continue;
      }

      PatternRoute patternRoute(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      patternRoute.setSteinerTree(steinerTree);
      patternRoute.constructRoutingDAG();
      patternRoute.run();
      const auto candidateTree = net->getRoutingTree();
      if (!candidateTree) {
        continue;
      }

      grid_graph_->commitTree(candidateTree);
      const RouteScore candidateScore
          = evaluateRouteScore(candidateTree, grid_graph_.get());
      grid_graph_->commitTree(candidateTree, /*ripup*/ true);
      if (isBetterScore(candidateScore, bestScore)) {
        bestTree = candidateTree;
        bestScore = candidateScore;
      }
    }

    if (!bestTree && oldTree) {
      bestTree = oldTree;
    }

    if (bestTree) {
      if (bestTree != oldTree) {
        accepted++;
      }
      net->setRoutingTree(bestTree);
      grid_graph_->commitTree(bestTree);
      grid_graph_->updateWireCostView(wireCostView, bestTree);
    } else if (oldTree) {
      net->setRoutingTree(oldTree);
      grid_graph_->commitTree(oldTree);
      grid_graph_->updateWireCostView(wireCostView, oldTree);
    }
    rank++;
  }

  logger_->report("stage 3 accepted {} / {} nets ({} candidates tested).",
                  accepted,
                  netIndices.size(),
                  totalCandidates);
  updateOverflowNets(netIndices);
}

void CUGR::wirelengthRecovery()
{
  logger_->report(
      "stage 4: wirelength recovery with multi-candidate reroute");
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }
  // Prioritize nets by current committed wire length contribution so the
  // largest detours are tightened first.
  std::vector<uint64_t> routedWireLength(gr_nets_.size(), 0);
  for (const int netIndex : netIndices) {
    const auto tree = gr_nets_[netIndex]->getRoutingTree();
    routedWireLength[netIndex]
        = evaluateRouteScore(tree, grid_graph_.get()).wire_length;
  }
  std::sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (routedWireLength[lhs] != routedWireLength[rhs]) {
      return routedWireLength[lhs] > routedWireLength[rhs];
    }
    return lhs < rhs;
  });

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);

  const int mazeCandidateBudget
      = std::max(128, static_cast<int>(netIndices.size() / 4));
  const int denseMazeBudget
      = std::max(64, static_cast<int>(netIndices.size() / 20));
  const int longNetRank
      = std::min(static_cast<int>(netIndices.size()) - 1,
                 std::max(0, static_cast<int>(netIndices.size() / 8)));
  const uint64_t longWireThreshold
      = routedWireLength[netIndices[longNetRank]];
  int accepted = 0;
  int acceptedFromMaze = 0;
  int rank = 0;
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const auto oldTree = net->getRoutingTree();
    if (!oldTree) {
      rank++;
      continue;
    }

    const RouteScore oldScore = evaluateRouteScore(oldTree, grid_graph_.get());

    // Remove current net and keep sparse maze costs synchronized.
    grid_graph_->commitTree(oldTree, /*ripup*/ true);
    grid_graph_->updateWireCostView(wireCostView, oldTree);

    std::shared_ptr<GRTreeNode> bestTree = oldTree;
    RouteScore bestScore = oldScore;
    bool bestFromMaze = false;

    auto tryCandidate = [&](const std::shared_ptr<GRTreeNode>& candidateTree,
                            const bool fromMaze) {
      if (!candidateTree) {
        return;
      }
      grid_graph_->commitTree(candidateTree);
      const RouteScore candidateScore
          = evaluateRouteScore(candidateTree, grid_graph_.get());
      grid_graph_->commitTree(candidateTree, /*ripup*/ true);
      if (isRecoveryScoreBetter(candidateScore, bestScore)) {
        bestTree = candidateTree;
        bestScore = candidateScore;
        bestFromMaze = fromMaze;
      }
    };

    // Candidate A: CUGR-style pattern reroute.
    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();
    tryCandidate(net->getRoutingTree(), /*fromMaze*/ false);

    // Candidate B: SPRoute-like selective sparse maze refinement for
    // critical/overflow nets.
    const bool runMazeCandidate
        = (rank < mazeCandidateBudget || oldScore.overflow_edges > 0
           || oldScore.wire_length >= longWireThreshold);
    if (runMazeCandidate) {
      auto runSparseMazeCandidate = [&](int interval, int xOffset, int yOffset) {
        SparseGrid recoveryGrid(interval, interval, xOffset, yOffset);

        MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
        mazeRoute.constructSparsifiedGraph(wireCostView, recoveryGrid);
        mazeRoute.run();
        std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();
        if (!steinerTree) {
          return;
        }
        PatternRoute mazeRefineRoute(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        mazeRefineRoute.setSteinerTree(steinerTree);
        mazeRefineRoute.constructRoutingDAG();
        mazeRefineRoute.run();
        tryCandidate(net->getRoutingTree(), /*fromMaze*/ true);
      };

      int interval = oldScore.overflow_edges > 0 ? 4 : 5;
      if (rank < denseMazeBudget || oldScore.wire_length >= longWireThreshold) {
        interval = 3;
      }
      runSparseMazeCandidate(
          interval,
          (rank * 3 + oldScore.overflow_edges) % interval,
          (rank * 5 + net->getNumPins()) % interval);
    }

    if (bestTree != oldTree) {
      accepted++;
      if (bestFromMaze) {
        acceptedFromMaze++;
      }
    }

    net->setRoutingTree(bestTree);
    grid_graph_->commitTree(bestTree);
    grid_graph_->updateWireCostView(wireCostView, bestTree);
    rank++;
  }

  logger_->report("wirelength recovery accepted {} net updates ({} from maze).",
                  accepted,
                  acceptedFromMaze);
}

void CUGR::finalPatternTighten()
{
  logger_->report("stage 5: final pattern tightening on long nets");
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }

  std::vector<uint64_t> routedWireLength(gr_nets_.size(), 0);
  for (const int netIndex : netIndices) {
    routedWireLength[netIndex]
        = evaluateRouteScore(gr_nets_[netIndex]->getRoutingTree(),
                             grid_graph_.get())
              .wire_length;
  }
  std::sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (routedWireLength[lhs] != routedWireLength[rhs]) {
      return routedWireLength[lhs] > routedWireLength[rhs];
    }
    return lhs < rhs;
  });

  const int tightenBudget
      = std::min(static_cast<int>(netIndices.size()),
                 std::max(4096, static_cast<int>(netIndices.size() / 2)));
  const int mazeTightenBudget
      = std::min(tightenBudget, std::max(512, tightenBudget / 12));
  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);
  int accepted = 0;
  int acceptedFromMaze = 0;
  for (int rank = 0; rank < tightenBudget; rank++) {
    const int netIndex = netIndices[rank];
    GRNet* net = gr_nets_[netIndex].get();
    const auto oldTree = net->getRoutingTree();
    if (!oldTree) {
      continue;
    }

    const RouteScore oldScore = evaluateRouteScore(oldTree, grid_graph_.get());
    grid_graph_->commitTree(oldTree, /*ripup*/ true);
    grid_graph_->updateWireCostView(wireCostView, oldTree);

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();
    const auto patternTree = net->getRoutingTree();

    std::shared_ptr<GRTreeNode> bestTree = oldTree;
    RouteScore bestScore = oldScore;
    if (patternTree) {
      grid_graph_->commitTree(patternTree);
      const RouteScore patternScore
          = evaluateRouteScore(patternTree, grid_graph_.get());
      grid_graph_->commitTree(patternTree, /*ripup*/ true);
      if (isTightenScoreBetter(patternScore, bestScore)) {
        bestTree = patternTree;
        bestScore = patternScore;
      }
    }

    if (rank < mazeTightenBudget) {
      const int interval = rank < mazeTightenBudget / 4 ? 3 : 4;
      const int xOffset = (rank * 5) % interval;
      const int yOffset = (rank * 7) % interval;
      SparseGrid tightenGrid(interval, interval, xOffset, yOffset);
      MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
      mazeRoute.constructSparsifiedGraph(wireCostView, tightenGrid);
      mazeRoute.run();
      std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();
      if (steinerTree) {
        PatternRoute mazeRefineRoute(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        mazeRefineRoute.setSteinerTree(steinerTree);
        mazeRefineRoute.constructRoutingDAG();
        mazeRefineRoute.run();
        const auto mazeTree = net->getRoutingTree();
        if (mazeTree) {
          grid_graph_->commitTree(mazeTree);
          const RouteScore mazeScore
              = evaluateRouteScore(mazeTree, grid_graph_.get());
          grid_graph_->commitTree(mazeTree, /*ripup*/ true);
          if (isTightenScoreBetter(mazeScore, bestScore)) {
            bestTree = mazeTree;
            bestScore = mazeScore;
            acceptedFromMaze++;
          }
        }
      }
    }

    if (bestTree != oldTree) {
      accepted++;
    }
    net->setRoutingTree(bestTree);
    grid_graph_->commitTree(bestTree);
    grid_graph_->updateWireCostView(wireCostView, bestTree);
  }

  logger_->report(
      "final pattern tightening accepted {} net updates ({} from maze).",
      accepted,
      acceptedFromMaze);
}

void CUGR::globalCompaction()
{
  logger_->report("stage 6: multi-round global compaction with strict score "
                  "gating");
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }
  if (netIndices.empty()) {
    return;
  }

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);

  enum class CandidateSource
  {
    Original,
    Pattern,
    Maze
  };

  struct CompactionRound
  {
    int budget_divisor;
    int min_budget;
    bool use_x_axis;
    int maze_candidate_base;
  };

  const std::vector<CompactionRound> rounds{
      {3, 4096, true, 2},
      {6, 2048, false, 2},
      {12, 1024, true, 1}};

  int totalAccepted = 0;
  int totalAcceptedPattern = 0;
  int totalAcceptedMaze = 0;

  std::vector<RouteScore> routedScores(gr_nets_.size());
  for (size_t roundIdx = 0; roundIdx < rounds.size(); roundIdx++) {
    const auto& round = rounds[roundIdx];
    for (const int netIndex : netIndices) {
      routedScores[netIndex]
          = evaluateRouteScore(gr_nets_[netIndex]->getRoutingTree(),
                               grid_graph_.get());
    }
    std::sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
      if (routedScores[lhs].overflow_edges != routedScores[rhs].overflow_edges) {
        return routedScores[lhs].overflow_edges > routedScores[rhs].overflow_edges;
      }
      if (routedScores[lhs].wire_length != routedScores[rhs].wire_length) {
        return routedScores[lhs].wire_length > routedScores[rhs].wire_length;
      }
      if (routedScores[lhs].via_count != routedScores[rhs].via_count) {
        return routedScores[lhs].via_count > routedScores[rhs].via_count;
      }
      return lhs < rhs;
    });

    const int totalNets = static_cast<int>(netIndices.size());
    const int compactionBudget = std::min(totalNets,
                                          std::max(round.min_budget,
                                                   totalNets
                                                       / round.budget_divisor));
    if (compactionBudget <= 0) {
      continue;
    }

    const int longNetRank
        = std::min(compactionBudget - 1, std::max(0, compactionBudget / 6));
    const uint64_t longWireThreshold
        = routedScores[netIndices[longNetRank]].wire_length;
    const int denseMazeBudget
        = std::min(compactionBudget,
                   std::max(256,
                            compactionBudget
                                / std::max(2, static_cast<int>(roundIdx) + 2)));
    std::vector<int> scheduledNetIndices = buildSpatialCompactionOrder(
        netIndices, gr_nets_, compactionBudget, round.use_x_axis);
    if (scheduledNetIndices.empty()) {
      scheduledNetIndices.assign(netIndices.begin(),
                                 netIndices.begin() + compactionBudget);
    }

    int accepted = 0;
    int acceptedPattern = 0;
    int acceptedMaze = 0;
    const int scheduledCount = static_cast<int>(scheduledNetIndices.size());
    for (int rank = 0; rank < scheduledCount; rank++) {
      const int netIndex = scheduledNetIndices[rank];
      GRNet* net = gr_nets_[netIndex].get();
      const auto oldTree = net->getRoutingTree();
      if (!oldTree) {
        continue;
      }
      const RouteScore oldScore
          = evaluateRouteScore(oldTree, grid_graph_.get());

      grid_graph_->commitTree(oldTree, /*ripup*/ true);
      grid_graph_->updateWireCostView(wireCostView, oldTree);

      std::shared_ptr<GRTreeNode> bestTree = oldTree;
      RouteScore bestScore = oldScore;
      CandidateSource bestSource = CandidateSource::Original;

      auto tryCandidate = [&](const std::shared_ptr<GRTreeNode>& candidateTree,
                              const CandidateSource source) {
        if (!candidateTree) {
          return;
        }
        grid_graph_->commitTree(candidateTree);
        const RouteScore candidateScore
            = evaluateRouteScore(candidateTree, grid_graph_.get());
        grid_graph_->commitTree(candidateTree, /*ripup*/ true);
        if (isCompactionScoreBetter(candidateScore, bestScore)) {
          bestTree = candidateTree;
          bestScore = candidateScore;
          bestSource = source;
        }
      };

      PatternRoute patternRoute(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      patternRoute.constructSteinerTree();
      patternRoute.constructRoutingDAG();
      patternRoute.run();
      tryCandidate(net->getRoutingTree(), CandidateSource::Pattern);

      const bool runMazeCandidate
          = rank < denseMazeBudget || oldScore.overflow_edges > 0
            || oldScore.wire_length >= longWireThreshold
            || oldScore.via_count >= (roundIdx == 0 ? 8 : 6);
      if (runMazeCandidate) {
        const int hp = net->getBoundingBox().hp();
        int interval = 4;
        if (rank < denseMazeBudget / 3 || oldScore.overflow_edges > 0) {
          interval = 3;
        } else if (roundIdx > 0 && net->getNumPins() <= 3 && hp <= 70) {
          interval = 5;
        }
        const int maxMazeCandidates = std::clamp(
            round.maze_candidate_base
                + ((rank < denseMazeBudget / 3
                    || oldScore.wire_length >= longWireThreshold)
                       ? 1
                       : 0),
            1,
            4);
        const int seed = rank + oldScore.via_count
                         + static_cast<int>(roundIdx) * 97;
        const auto candidateGrids
            = buildMazeCandidateGrids(
                interval, seed, hp, net->getNumPins(), maxMazeCandidates);
        for (const auto& compactionGrid : candidateGrids) {
          MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
          mazeRoute.constructSparsifiedGraph(wireCostView, compactionGrid);
          mazeRoute.run();
          std::shared_ptr<SteinerTreeNode> steinerTree
              = mazeRoute.getSteinerTree();
          if (!steinerTree) {
            continue;
          }
          PatternRoute mazeRefineRoute(
              net, grid_graph_.get(), stt_builder_, constants_, logger_);
          mazeRefineRoute.setSteinerTree(steinerTree);
          mazeRefineRoute.constructRoutingDAG();
          mazeRefineRoute.run();
          tryCandidate(net->getRoutingTree(), CandidateSource::Maze);
        }
      }

      if (bestTree != oldTree) {
        accepted++;
        if (bestSource == CandidateSource::Pattern) {
          acceptedPattern++;
        } else if (bestSource == CandidateSource::Maze) {
          acceptedMaze++;
        }
      }

      net->setRoutingTree(bestTree);
      grid_graph_->commitTree(bestTree);
      grid_graph_->updateWireCostView(wireCostView, bestTree);
    }

    totalAccepted += accepted;
    totalAcceptedPattern += acceptedPattern;
    totalAcceptedMaze += acceptedMaze;
    logger_->report("stage 6 round {} scheduled {} nets (axis={}): accepted "
                    "{} ({} pattern / {} maze).",
                    roundIdx + 1,
                    scheduledCount,
                    round.use_x_axis ? "x" : "y",
                    accepted,
                    acceptedPattern,
                    acceptedMaze);
  }

  logger_->report("global compaction accepted {} net updates total ({} pattern "
                  "/ {} maze).",
                  totalAccepted,
                  totalAcceptedPattern,
                  totalAcceptedMaze);
}

void CUGR::strictWirelengthCompaction()
{
  logger_->report("stage 7: strict wirelength compaction on detour-heavy nets");
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }
  if (netIndices.empty()) {
    return;
  }

  std::vector<RouteScore> routedScores(gr_nets_.size());
  for (const int netIndex : netIndices) {
    routedScores[netIndex]
        = evaluateRouteScore(gr_nets_[netIndex]->getRoutingTree(),
                             grid_graph_.get());
  }
  std::sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (routedScores[lhs].wire_length != routedScores[rhs].wire_length) {
      return routedScores[lhs].wire_length > routedScores[rhs].wire_length;
    }
    if (routedScores[lhs].overflow_edges != routedScores[rhs].overflow_edges) {
      return routedScores[lhs].overflow_edges > routedScores[rhs].overflow_edges;
    }
    if (routedScores[lhs].via_count != routedScores[rhs].via_count) {
      return routedScores[lhs].via_count > routedScores[rhs].via_count;
    }
    return lhs < rhs;
  });

  const int totalNets = static_cast<int>(netIndices.size());
  static int strictCallCount = 0;
  const bool useXAxisWavefront = (strictCallCount % 2 == 0);
  strictCallCount++;
  const int compactionBudget
      = std::min(totalNets, std::max(6144, (totalNets * 3) / 5));
  if (compactionBudget <= 0) {
    return;
  }

  const int longNetRank
      = std::min(compactionBudget - 1, std::max(0, compactionBudget / 6));
  const uint64_t longWireThreshold
      = routedScores[netIndices[longNetRank]].wire_length;
  const int denseMazeBudget
      = std::min(compactionBudget, std::max(1536, compactionBudget / 2));
  std::vector<int> scheduledNetIndices = buildSpatialCompactionOrder(
      netIndices, gr_nets_, compactionBudget, useXAxisWavefront);
  if (scheduledNetIndices.empty()) {
    scheduledNetIndices.assign(netIndices.begin(),
                               netIndices.begin() + compactionBudget);
  }

  GridGraphView<CostT> wireCostView;
  grid_graph_->extractWireCostView(wireCostView);

  int accepted = 0;
  int acceptedPattern = 0;
  int acceptedMaze = 0;
  const int scheduledCount = static_cast<int>(scheduledNetIndices.size());
  for (int rank = 0; rank < scheduledCount; rank++) {
    const int netIndex = scheduledNetIndices[rank];
    GRNet* net = gr_nets_[netIndex].get();
    const auto oldTree = net->getRoutingTree();
    if (!oldTree) {
      continue;
    }

    const RouteScore oldScore = evaluateRouteScore(oldTree, grid_graph_.get());
    grid_graph_->commitTree(oldTree, /*ripup*/ true);
    grid_graph_->updateWireCostView(wireCostView, oldTree);

    std::shared_ptr<GRTreeNode> bestTree = oldTree;
    RouteScore bestScore = oldScore;
    bool bestFromMaze = false;

    auto tryCandidate = [&](const std::shared_ptr<GRTreeNode>& candidateTree,
                            const bool fromMaze) {
      if (!candidateTree) {
        return;
      }
      grid_graph_->commitTree(candidateTree);
      const RouteScore candidateScore
          = evaluateRouteScore(candidateTree, grid_graph_.get());
      grid_graph_->commitTree(candidateTree, /*ripup*/ true);
      if (isStrictWirelengthScoreBetter(candidateScore, bestScore)) {
        bestTree = candidateTree;
        bestScore = candidateScore;
        bestFromMaze = fromMaze;
      }
    };

    PatternRoute patternRoute(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    patternRoute.constructSteinerTree();
    patternRoute.constructRoutingDAG();
    patternRoute.run();
    tryCandidate(net->getRoutingTree(), /*fromMaze*/ false);

    const bool runMazeCandidate
        = rank < denseMazeBudget || oldScore.overflow_edges > 0
          || oldScore.wire_length >= longWireThreshold
          || oldScore.via_count >= 8;
    if (runMazeCandidate) {
      const int hp = net->getBoundingBox().hp();
      const int pins = net->getNumPins();
      int interval = 4;
      if (rank < denseMazeBudget * 2 / 3 || oldScore.overflow_edges > 0) {
        interval = 3;
      }
      if (rank < denseMazeBudget / 6
          || (oldScore.wire_length >= longWireThreshold && pins >= 6)) {
        interval = 2;
      } else if (pins <= 3 && hp <= 80) {
        interval = 5;
      }
      const int maxMazeCandidates = rank < denseMazeBudget / 6
                                        ? 6
                                        : (rank < denseMazeBudget / 2 ? 5 : 4);
      const auto candidateGrids
          = buildMazeCandidateGrids(interval,
                                    rank + oldScore.via_count * 3
                                        + strictCallCount * 31,
                                    hp,
                                    pins,
                                    maxMazeCandidates);
      for (const auto& compactionGrid : candidateGrids) {
        MazeRoute mazeRoute(net, grid_graph_.get(), logger_);
        mazeRoute.constructSparsifiedGraph(wireCostView, compactionGrid);
        mazeRoute.run();
        std::shared_ptr<SteinerTreeNode> steinerTree = mazeRoute.getSteinerTree();
        if (!steinerTree) {
          continue;
        }
        PatternRoute mazeRefineRoute(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        mazeRefineRoute.setSteinerTree(steinerTree);
        mazeRefineRoute.constructRoutingDAG();
        mazeRefineRoute.run();
        tryCandidate(net->getRoutingTree(), /*fromMaze*/ true);
      }
    }

    if (bestTree != oldTree) {
      accepted++;
      if (bestFromMaze) {
        acceptedMaze++;
      } else {
        acceptedPattern++;
      }
    }

    net->setRoutingTree(bestTree);
    grid_graph_->commitTree(bestTree);
    grid_graph_->updateWireCostView(wireCostView, bestTree);
  }

  logger_->report("stage 7 scheduled {} nets (axis={}): accepted {} updates "
                  "({} pattern / {} maze).",
                  scheduledCount,
                  useXAxisWavefront ? "x" : "y",
                  accepted,
                  acceptedPattern,
                  acceptedMaze);
}

void CUGR::route()
{
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }

  grid_graph_->setSoftCapacityEnabled(false);

  // FastRoute-style adaptive emphasis: start wirelength-first and raise
  // congestion pressure in the middle RRR passes.
  grid_graph_->setStageCostScales(0.75, 0.75, 1.20);
  patternRoute(netIndices);

  // Run maze reroute before detours so most overflow repairs come from a
  // shortest-path engine instead of detour inflation.
  grid_graph_->setStageCostScales(1.10, 1.15, 1.10);
  mazeRoute(netIndices);

  // Keep detours as a final cleanup pass for residual difficult hotspots.
  grid_graph_->setStageCostScales(1.22, 1.24, 1.00);
  patternRouteWithDetours(netIndices);

  // FastRoute-style final RRR cleanup: re-run maze search to pull inflated
  // detours back to shorter legal paths after hotspot repair.
  grid_graph_->setStageCostScales(1.32, 1.35, 1.05);
  mazeRoute(netIndices);

  // FastRoute-inspired post-congestion tightening: re-run pure pattern
  // routing and accept only net-level improvements in
  // overflow/wirelength/via score.
  grid_graph_->setSoftCapacityEnabled(false);
  grid_graph_->setStageCostScales(0.46, 0.50, 1.20);
  wirelengthRecovery();
  grid_graph_->setStageCostScales(0.36, 0.38, 1.15);
  finalPatternTighten();
  grid_graph_->setStageCostScales(0.18, 0.20, 1.05);
  globalCompaction();
  grid_graph_->setStageCostScales(0.09, 0.10, 0.95);
  strictWirelengthCompaction();
  grid_graph_->setStageCostScales(0.04, 0.05, 0.90);
  strictWirelengthCompaction();
  grid_graph_->setStageCostScales(0.01, 0.02, 0.88);
  strictWirelengthCompaction();
  grid_graph_->setStageCostScales(0.0, 0.01, 0.92);
  strictWirelengthCompaction();

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
  std::vector<double> priorities(gr_nets_.size(), 0.0);
  std::vector<int> halfParameters(gr_nets_.size(), 0);
  for (int netIndex : netIndices) {
    auto& net = gr_nets_[netIndex];
    const int hp = net->getBoundingBox().hp();
    const int pins = std::max(net->getNumPins(), 2);
    // FastRoute-inspired priority: protect large/high-fanout nets first.
    priorities[netIndex] = static_cast<double>(hp) * (1.0 + std::log2(pins));
    halfParameters[netIndex] = hp;
  }
  sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (priorities[lhs] != priorities[rhs]) {
      return priorities[lhs] > priorities[rhs];
    }
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
