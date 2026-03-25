#include "MazeRoute.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "Design.h"
#include "GridGraph.h"
#include "Layers.h"
#include "PatternRoute.h"
#include "geo.h"
#include "robin_hood.h"
#include "utl/Logger.h"

namespace grt::newgr {

void SparseGraph::init(const GridGraphView<CostT>& wire_cost_view,
                       const SparseGrid& grid)
{
  pseudo_pins_.clear();
  xs_.clear();
  ys_.clear();
  vertices_.clear();
  edges_.clear();
  costs_.clear();
  vertex_pins_.clear();
  pin_vertex_.clear();

  const int xSize = grid_graph_->getSize(0);
  const int ySize = grid_graph_->getSize(1);

  // 0. Create pseudo pins
  const auto selectedAccessPoints = grid_graph_->selectAccessPoints(net_);
  pseudo_pins_.reserve(selectedAccessPoints.size());
  for (const auto& selectedPoint : selectedAccessPoints) {
    pseudo_pins_.push_back(selectedPoint);
  }

  // 1. Collect additional routing grid lines
  std::vector<int> pxs;
  std::vector<int> pys;
  pxs.reserve(net_->getNumPins() * 2);
  pys.reserve(net_->getNumPins() * 2);
  for (const auto& pin : pseudo_pins_) {
    pxs.emplace_back(pin.point.x());
    pys.emplace_back(pin.point.y());
  }

  const auto& box = net_->getBoundingBox();
  const int pins = std::max(2, net_->getNumPins());
  const int hp = std::max(1, box.hp());
  const bool shortestTopologyMode
      = hp >= 120 || pins >= 8
        || std::max(grid.interval.x(), grid.interval.y()) <= 3;
  const auto& oldTree = net_->getRoutingTree();
  const int oldTreeOverflow = oldTree ? grid_graph_->checkOverflow(oldTree) : 0;
  const bool preserveOldTopologyAnchors = oldTreeOverflow > 0;
  // Keep only nearby old-tree anchors so maze reroute can aggressively compact
  // long detours instead of preserving the previous expanded topology.
  int anchorMargin = shortestTopologyMode
                         ? std::clamp(hp / (pins >= 12 ? 20 : 18), 1, 8)
                         : std::clamp(hp / (pins >= 12 ? 10 : 9), 3, 22);
  if (hp >= 220 || pins >= 18) {
    anchorMargin = shortestTopologyMode ? std::min(anchorMargin + 1, 12)
                                        : std::min(anchorMargin + 2, 24);
  }
  const int anchorXLow = std::max(0, box.lx() - anchorMargin);
  const int anchorXHigh = std::min(xSize - 1, box.hx() + anchorMargin);
  const int anchorYLow = std::max(0, box.ly() - anchorMargin);
  const int anchorYHigh = std::min(ySize - 1, box.hy() + anchorMargin);
  if (oldTree && anchorMargin > 0 && preserveOldTopologyAnchors) {
    GRTreeNode::preorder(
        oldTree, [&](const std::shared_ptr<GRTreeNode>& node) {
          if (node->x() >= anchorXLow && node->x() <= anchorXHigh
              && node->y() >= anchorYLow && node->y() <= anchorYHigh) {
            pxs.emplace_back(node->x());
            pys.emplace_back(node->y());
          }
        });
  }

  pxs.emplace_back(std::clamp(box.cx(), 0, xSize - 1));
  pys.emplace_back(std::clamp(box.cy(), 0, ySize - 1));
  pxs.emplace_back(std::clamp(box.lx(), 0, xSize - 1));
  pxs.emplace_back(std::clamp(box.hx(), 0, xSize - 1));
  pys.emplace_back(std::clamp(box.ly(), 0, ySize - 1));
  pys.emplace_back(std::clamp(box.hy(), 0, ySize - 1));

  // CUGR-style coarse-to-fine corridor search: bound sparse-graph expansion to
  // each net neighborhood to suppress long global detours.
  int margin = shortestTopologyMode
                   ? std::clamp(hp / (pins >= 12 ? 15 : 14), 3, 28)
                   : std::clamp(hp / (pins >= 12 ? 8 : 7), 5, 52);
  if (pins <= 3) {
    margin = std::max(margin, 8);
  }
  if (hp >= 220 || pins >= 18) {
    margin += shortestTopologyMode ? std::max(1, hp / 90)
                                   : std::max(2, hp / 65);
  }
  if (shortestTopologyMode
      && std::max(grid.interval.x(), grid.interval.y()) <= 3) {
    margin = std::max(4, margin - 2);
  }
  // If the old tree is overflow-clean, avoid preserving its topology too
  // aggressively so sparse maze can collapse detours toward shorter trees.
  if (!preserveOldTopologyAnchors) {
    margin = shortestTopologyMode ? std::max(2, margin - 5)
                                  : std::max(4, margin - 2);
  }
  if (grid.interval.x() == 1 || grid.interval.y() == 1) {
    // Full-resolution sparse grid is already expressive; shrink the corridor to
    // aggressively collapse legacy detours toward bbox-local shortest paths.
    const int compactMargin
        = std::clamp(hp / (pins >= 10 ? 20 : 16), 2, 12);
    margin = std::min(margin, compactMargin);
    if (!preserveOldTopologyAnchors) {
      margin = std::max(2, margin - 1);
    }
  }
  margin += std::max(1, std::max(grid.interval.x(), grid.interval.y()) / 2);
  int xLow = std::max(0, box.lx() - margin);
  int xHigh = std::min(xSize - 1, box.hx() + margin);
  int yLow = std::max(0, box.ly() - margin);
  int yHigh = std::min(ySize - 1, box.hy() + margin);

  if (xLow == xHigh) {
    if (xHigh + 1 < xSize) {
      xHigh += 1;
    } else if (xLow > 0) {
      xLow -= 1;
    }
  }
  if (yLow == yHigh) {
    if (yHigh + 1 < ySize) {
      yHigh += 1;
    } else if (yLow > 0) {
      yLow -= 1;
    }
  }

  auto buildSparseAxis = [](std::vector<int>& axis,
                            const std::vector<int>& pinCoords,
                            int low,
                            int high,
                            int interval,
                            int offset) {
    std::vector<int> coords;
    coords.reserve(pinCoords.size() + 16);
    for (const int c : pinCoords) {
      if (c >= low && c <= high) {
        coords.push_back(c);
      }
    }
    coords.push_back(low);
    coords.push_back(high);
    coords.push_back((low + high) / 2);
    if (high - low > 2 * interval) {
      coords.push_back(std::clamp(low + interval, low, high));
      coords.push_back(std::clamp(high - interval, low, high));
    }

    const int first = (low - offset + interval - 1) / interval;
    const int last = (high - offset) / interval;
    for (int i = first; i <= last; i++) {
      const int c = i * interval + offset;
      if (c >= low && c <= high) {
        coords.push_back(c);
      }
    }

    std::sort(coords.begin(), coords.end());
    coords.erase(std::unique(coords.begin(), coords.end()), coords.end());
    axis.swap(coords);
  };

  buildSparseAxis(
      xs_, pxs, xLow, xHigh, grid.interval.x(), grid.offset.x());
  buildSparseAxis(
      ys_, pys, yLow, yHigh, grid.interval.y(), grid.offset.y());

  // 2. Add vertices
  vertices_.reserve(2 * xs_.size() * ys_.size());
  for (int direction = 0; direction < 2; direction++) {
    for (auto& y : ys_) {
      for (auto& x : xs_) {
        vertices_.emplace_back(direction, x, y);
      }
    }
  }

  // 3. Add same-layer connections
  edges_.resize(vertices_.size(), {-1, -1, -1});
  costs_.resize(vertices_.size(), {-1, -1, -1});
  auto addSameLayerEdge = [&](const int direction, const int xi, const int yi) {
    const int u = getVertexIndex(direction, xi, yi);
    const int v = direction == MetalLayer::H ? u + 1 : u + xs_.size();
    const PointT U(xs_[xi], ys_[yi]);
    const PointT V(xs_[xi + 1 - direction], ys_[yi + direction]);

    edges_[u][0] = v;
    edges_[v][1] = u;
    costs_[u][0] = costs_[v][1] = wire_cost_view.sum(U, V);
  };

  for (int direction = 0; direction < 2; direction++) {
    if (direction == MetalLayer::H) {
      for (int yi = 0; yi < ys_.size(); yi++) {
        for (int xi = 0; xi + 1 < xs_.size(); xi++) {
          addSameLayerEdge(direction, xi, yi);
        }
      }
    } else {
      for (int xi = 0; xi < xs_.size(); xi++) {
        for (int yi = 0; yi + 1 < ys_.size(); yi++) {
          addSameLayerEdge(direction, xi, yi);
        }
      }
    }
  }

  // 4. Add diff-layer connections
  auto addDiffLayerEdge = [&](const int xi, const int yi) {
    const int u = getVertexIndex(0, xi, yi);
    const int v = u + xs_.size() * ys_.size();

    edges_[u][2] = v;
    edges_[v][2] = u;
    costs_[u][2] = costs_[v][2] = grid_graph_->getUnitViaCost();
  };

  for (int xi = 0; xi < xs_.size(); xi++) {
    for (int yi = 0; yi < ys_.size(); yi++) {
      addDiffLayerEdge(xi, yi);
    }
  }

  // 5. Add pseudo pin locations
  robin_hood::unordered_map<int, int> xtoxi;
  robin_hood::unordered_map<int, int> ytoyi;
  for (int xi = 0; xi < xs_.size(); xi++) {
    xtoxi.emplace(xs_[xi], xi);
  }
  for (int yi = 0; yi < ys_.size(); yi++) {
    ytoyi.emplace(ys_[yi], yi);
  }

  pin_vertex_.resize(pseudo_pins_.size(), -1);
  for (int pinIndex = 0; pinIndex < pseudo_pins_.size(); pinIndex++) {
    const auto& pin = pseudo_pins_[pinIndex];
    const int xi = xtoxi[pin.point.x()];
    const int yi = ytoyi[pin.point.y()];
    const int u = getVertexIndex(0, xi, yi);
    vertex_pins_[u].push_back(pinIndex);
    pin_vertex_[pinIndex] = u;
    // Set the cost of the diff-layer connection at u to be 0
    costs_[u][2] = 0;
    costs_[u + xs_.size() * ys_.size()][2] = 0;
  }
}

void MazeRoute::run()
{
  solutions_.clear();
  const int numPseudoPins = graph_.getNumPseudoPins();
  if (numPseudoPins <= 0) {
    return;
  }
  const auto& box = net_->getBoundingBox();
  const PointT center(box.cx(), box.cy());
  std::vector<int64_t> centerDistances(numPseudoPins, 0);

  int startPinIndex = -1;
  int64_t bestTotalDistance = std::numeric_limits<int64_t>::max();
  int64_t bestCenterDistance = std::numeric_limits<int64_t>::max();
  for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
    const auto& pseudoPin = graph_.getPseudoPin(pinIndex);
    int64_t totalDistance = 0;
    for (int otherPin = 0; otherPin < numPseudoPins; otherPin++) {
      if (otherPin == pinIndex) {
        continue;
      }
      const auto& other = graph_.getPseudoPin(otherPin);
      totalDistance += std::llabs(static_cast<int64_t>(pseudoPin.point.x())
                                  - other.point.x());
      totalDistance += std::llabs(static_cast<int64_t>(pseudoPin.point.y())
                                  - other.point.y());
    }
    const int64_t centerDistance
        = std::llabs(static_cast<int64_t>(pseudoPin.point.x()) - center.x())
          + std::llabs(static_cast<int64_t>(pseudoPin.point.y()) - center.y());
    centerDistances[pinIndex] = centerDistance;
    if (totalDistance < bestTotalDistance
        || (totalDistance == bestTotalDistance
            && centerDistance < bestCenterDistance)) {
      bestTotalDistance = totalDistance;
      bestCenterDistance = centerDistance;
      startPinIndex = pinIndex;
    }
  }

  auto addStartCandidate = [](std::vector<int>& starts, int candidate) {
    if (candidate < 0) {
      return;
    }
    if (std::find(starts.begin(), starts.end(), candidate) == starts.end()) {
      starts.push_back(candidate);
    }
  };

  std::vector<int> startCandidates;
  addStartCandidate(startCandidates, startPinIndex);

  if (numPseudoPins >= 4) {
    int farthestFromPrimary = -1;
    int64_t farthestDist = -1;
    const auto& primaryPin = graph_.getPseudoPin(startPinIndex);
    for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
      const auto& pseudoPin = graph_.getPseudoPin(pinIndex);
      const int64_t dist
          = std::llabs(static_cast<int64_t>(pseudoPin.point.x())
                       - primaryPin.point.x())
            + std::llabs(static_cast<int64_t>(pseudoPin.point.y())
                         - primaryPin.point.y());
      if (dist > farthestDist) {
        farthestDist = dist;
        farthestFromPrimary = pinIndex;
      }
    }
    addStartCandidate(startCandidates, farthestFromPrimary);
  }
  if (numPseudoPins >= 6) {
    int farthestFromCenter = -1;
    int64_t farthestCenterDistance = -1;
    for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
      if (centerDistances[pinIndex] > farthestCenterDistance) {
        farthestCenterDistance = centerDistances[pinIndex];
        farthestFromCenter = pinIndex;
      }
    }
    addStartCandidate(startCandidates, farthestFromCenter);
  }
  if (numPseudoPins >= 8) {
    int minXPin = -1;
    int maxXPin = -1;
    int minYPin = -1;
    int maxYPin = -1;
    for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
      const auto& pseudoPin = graph_.getPseudoPin(pinIndex);
      if (minXPin == -1
          || pseudoPin.point.x() < graph_.getPseudoPin(minXPin).point.x()) {
        minXPin = pinIndex;
      }
      if (maxXPin == -1
          || pseudoPin.point.x() > graph_.getPseudoPin(maxXPin).point.x()) {
        maxXPin = pinIndex;
      }
      if (minYPin == -1
          || pseudoPin.point.y() < graph_.getPseudoPin(minYPin).point.y()) {
        minYPin = pinIndex;
      }
      if (maxYPin == -1
          || pseudoPin.point.y() > graph_.getPseudoPin(maxYPin).point.y()) {
        maxYPin = pinIndex;
      }
    }
    addStartCandidate(startCandidates, minXPin);
    addStartCandidate(startCandidates, maxXPin);
    addStartCandidate(startCandidates, minYPin);
    addStartCandidate(startCandidates, maxYPin);
  }

  const int maxStartCandidates
      = numPseudoPins >= 14 ? 8 : (numPseudoPins >= 8 ? 5 : 2);
  if (startCandidates.size() > static_cast<size_t>(maxStartCandidates)) {
    startCandidates.resize(maxStartCandidates);
  }

  auto runFromStartPin = [&](int startPin,
                             std::vector<std::shared_ptr<Solution>>& output) {
    output.clear();
    output.reserve(numPseudoPins);
    std::vector<CostT> minCosts(graph_.getNumVertices(),
                                std::numeric_limits<CostT>::max());
    auto compareSolution = [&](const std::shared_ptr<Solution>& lhs,
                               const std::shared_ptr<Solution>& rhs) {
      return lhs->cost > rhs->cost;
    };
    std::priority_queue<std::shared_ptr<Solution>,
                        std::vector<std::shared_ptr<Solution>>,
                        decltype(compareSolution)>
        queue(compareSolution);

    auto updateSolution = [&](const std::shared_ptr<Solution>& solution) {
      queue.push(solution);
      if (solution->cost < minCosts[solution->vertex]) {
        minCosts[solution->vertex] = solution->cost;
      }
    };

    std::vector<bool> visited(numPseudoPins, false);
    int numDetached = numPseudoPins;
    auto markVisitedPinsAtVertex = [&](int vertex) {
      int newlyVisited = 0;
      const auto& pinsAtVertex = graph_.getVertexPins(vertex);
      for (const int pinIndex : pinsAtVertex) {
        if (!visited[pinIndex]) {
          visited[pinIndex] = true;
          newlyVisited++;
        }
      }
      return newlyVisited;
    };
    numDetached -= markVisitedPinsAtVertex(graph_.getPinVertex(startPin));
    updateSolution(
        std::make_shared<Solution>(0, graph_.getPinVertex(startPin), nullptr));

    while (numDetached > 0) {
      std::shared_ptr<Solution> foundSolution;
      int foundPinIndex = -1;
      while (!queue.empty()) {
        auto solution = queue.top();
        queue.pop();
        const auto& pinsAtVertex = graph_.getVertexPins(solution->vertex);
        for (const int pinIndex : pinsAtVertex) {
          if (!visited[pinIndex]) {
            foundPinIndex = pinIndex;
            foundSolution = std::move(solution);
            break;
          }
        }
        if (foundSolution) {
          break;
        }
        if (solution->cost > minCosts[solution->vertex]) {
          continue;
        }
        for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
          const int nextVertex
              = graph_.getNextVertex(solution->vertex, edgeIndex);
          if (nextVertex == -1
              || (solution->prev && nextVertex == solution->prev->vertex)) {
            continue;
          }
          const CostT nextCost = solution->cost
                                 + graph_.getEdgeCost(solution->vertex,
                                                      edgeIndex);
          if (nextCost < minCosts[nextVertex]) {
            updateSolution(
                std::make_shared<Solution>(nextCost, nextVertex, solution));
          }
        }
      }

      if (!foundSolution || foundPinIndex < 0) {
        break;
      }

      output.emplace_back(foundSolution);
      numDetached -= markVisitedPinsAtVertex(foundSolution->vertex);

      // Multi-source expansion: every accepted path is added as zero-cost
      // frontier so the next connection can attach to any routed branch.
      std::shared_ptr<Solution> temp = std::move(foundSolution);
      while (temp && temp->cost != 0) {
        updateSolution(std::make_shared<Solution>(0, temp->vertex, temp->prev));
        temp = temp->prev;
      }
    }

    return numDetached == 0;
  };

  struct CandidateScore
  {
    uint64_t unique_wire_length = std::numeric_limits<uint64_t>::max();
    CostT total_path_cost = std::numeric_limits<CostT>::max();
    int via_steps = std::numeric_limits<int>::max();
  };

  auto scoreSolutions = [&](const std::vector<std::shared_ptr<Solution>>& sols) {
    CandidateScore score;
    score.unique_wire_length = 0;
    score.total_path_cost = 0;
    score.via_steps = 0;
    robin_hood::unordered_set<uint64_t> visitedEdges;
    visitedEdges.reserve(sols.size() * 16);
    for (const auto& solution : sols) {
      if (!solution) {
        continue;
      }
      score.total_path_cost += solution->cost;
      std::shared_ptr<Solution> temp = solution;
      while (temp && temp->prev) {
        const int u = temp->vertex;
        const int v = temp->prev->vertex;
        const uint32_t lo = static_cast<uint32_t>(std::min(u, v));
        const uint32_t hi = static_cast<uint32_t>(std::max(u, v));
        const uint64_t key = (static_cast<uint64_t>(lo) << 32)
                             | static_cast<uint64_t>(hi);
        if (visitedEdges.insert(key).second) {
          const auto p = graph_.getPoint(u);
          const auto q = graph_.getPoint(v);
          const int dx = std::abs(p.x() - q.x());
          const int dy = std::abs(p.y() - q.y());
          if (dx == 0 && dy == 0) {
            score.via_steps += 1;
          } else {
            const int direction = dx > 0 ? MetalLayer::H : MetalLayer::V;
            const int l = std::min(p[direction], q[direction]);
            const int h = std::max(p[direction], q[direction]);
            for (int edge = l; edge < h; edge++) {
              score.unique_wire_length += static_cast<uint64_t>(
                  grid_graph_->getEdgeLength(direction, edge));
            }
          }
        }
        temp = temp->prev;
      }
    }
    return score;
  };

  auto isBetterCandidate = [&](const CandidateScore& lhs,
                               const CandidateScore& rhs) {
    if (rhs.unique_wire_length == std::numeric_limits<uint64_t>::max()) {
      return true;
    }
    if (lhs.unique_wire_length != rhs.unique_wire_length) {
      return lhs.unique_wire_length < rhs.unique_wire_length;
    }
    if (lhs.total_path_cost != rhs.total_path_cost) {
      return lhs.total_path_cost < rhs.total_path_cost;
    }
    return lhs.via_steps < rhs.via_steps;
  };

  const CostT kInfCost = std::numeric_limits<CostT>::max();
  const CostT kCostEps = 1e-9;
  auto edgeCostBetweenVertices = [&](int u, int v) {
    for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
      if (graph_.getNextVertex(u, edgeIndex) == v) {
        return graph_.getEdgeCost(u, edgeIndex);
      }
    }
    for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
      if (graph_.getNextVertex(v, edgeIndex) == u) {
        return graph_.getEdgeCost(v, edgeIndex);
      }
    }
    return kInfCost;
  };

  auto runMstCandidate = [&](std::vector<std::shared_ptr<Solution>>& output) {
    output.clear();

    const int numVertices = graph_.getNumVertices();
    const auto& box = net_->getBoundingBox();
    const bool criticalNet = numPseudoPins >= 6 || box.hp() >= 90;
    if (numPseudoPins < 3 || !criticalNet || numPseudoPins > 48
        || numVertices <= 0 || numVertices > 28000) {
      return false;
    }

    std::vector<int> pinVertices(numPseudoPins, -1);
    for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
      pinVertices[pinIndex] = graph_.getPinVertex(pinIndex);
      if (pinVertices[pinIndex] < 0) {
        return false;
      }
    }

    auto runDijkstra = [&](int sourceVertex,
                           std::vector<CostT>& dist,
                           std::vector<int>& prev) {
      dist.assign(numVertices, kInfCost);
      prev.assign(numVertices, -1);
      using QueueNode = std::pair<CostT, int>;
      std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<>>
          queue;
      dist[sourceVertex] = 0;
      queue.emplace(0, sourceVertex);
      while (!queue.empty()) {
        const auto [cost, vertex] = queue.top();
        queue.pop();
        if (cost > dist[vertex] + kCostEps) {
          continue;
        }
        for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++) {
          const int nextVertex = graph_.getNextVertex(vertex, edgeIndex);
          if (nextVertex == -1) {
            continue;
          }
          const CostT edgeCost = graph_.getEdgeCost(vertex, edgeIndex);
          if (edgeCost >= kInfCost / 4) {
            continue;
          }
          const CostT nextCost = cost + edgeCost;
          if (nextCost + kCostEps < dist[nextVertex]) {
            dist[nextVertex] = nextCost;
            prev[nextVertex] = vertex;
            queue.emplace(nextCost, nextVertex);
          }
        }
      }
    };

    std::vector<std::vector<CostT>> pinDistances(
        numPseudoPins,
        std::vector<CostT>(numPseudoPins, kInfCost));
    for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
      pinDistances[pinIndex][pinIndex] = 0;
    }
    for (int sourcePin = 0; sourcePin < numPseudoPins; sourcePin++) {
      std::vector<CostT> dist;
      std::vector<int> prev;
      runDijkstra(pinVertices[sourcePin], dist, prev);
      for (int targetPin = sourcePin + 1; targetPin < numPseudoPins;
           targetPin++) {
        const CostT distance = dist[pinVertices[targetPin]];
        pinDistances[sourcePin][targetPin] = distance;
        pinDistances[targetPin][sourcePin] = distance;
      }
    }

    // FastRoute-style multi-source reconnection converted to an MST over
    // pseudo pins using sparse-graph shortest-path costs.
    std::vector<CostT> bestEdge(numPseudoPins, kInfCost);
    std::vector<int> parent(numPseudoPins, -1);
    std::vector<bool> inMst(numPseudoPins, false);
    bestEdge[0] = 0;
    for (int iter = 0; iter < numPseudoPins; iter++) {
      int bestPin = -1;
      CostT minCost = kInfCost;
      for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
        if (!inMst[pinIndex] && bestEdge[pinIndex] < minCost) {
          minCost = bestEdge[pinIndex];
          bestPin = pinIndex;
        }
      }
      if (bestPin < 0 || minCost >= kInfCost / 4) {
        return false;
      }
      inMst[bestPin] = true;
      for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
        const CostT edgeCost = pinDistances[bestPin][pinIndex];
        if (!inMst[pinIndex] && edgeCost + kCostEps < bestEdge[pinIndex]) {
          bestEdge[pinIndex] = edgeCost;
          parent[pinIndex] = bestPin;
        }
      }
    }

    std::vector<bool> cachedSource(numPseudoPins, false);
    std::vector<std::vector<CostT>> cachedDist(numPseudoPins);
    std::vector<std::vector<int>> cachedPrev(numPseudoPins);
    auto ensureSourceData = [&](int sourcePin) {
      if (!cachedSource[sourcePin]) {
        runDijkstra(pinVertices[sourcePin],
                    cachedDist[sourcePin],
                    cachedPrev[sourcePin]);
        cachedSource[sourcePin] = true;
      }
      return true;
    };

    auto buildSolutionFromPinPair = [&](int sourcePin, int targetPin) {
      if (!ensureSourceData(sourcePin)) {
        return std::shared_ptr<Solution>();
      }
      const int sourceVertex = pinVertices[sourcePin];
      const int targetVertex = pinVertices[targetPin];
      if (cachedDist[sourcePin][targetVertex] >= kInfCost / 4) {
        return std::shared_ptr<Solution>();
      }

      std::vector<int> reversedPath;
      reversedPath.reserve(64);
      int current = targetVertex;
      while (current != -1
             && reversedPath.size() <= static_cast<size_t>(numVertices + 4)) {
        reversedPath.push_back(current);
        if (current == sourceVertex) {
          break;
        }
        current = cachedPrev[sourcePin][current];
      }
      if (reversedPath.empty() || reversedPath.back() != sourceVertex) {
        return std::shared_ptr<Solution>();
      }
      std::reverse(reversedPath.begin(), reversedPath.end());

      CostT cumulativeCost = 0;
      std::shared_ptr<Solution> pathSolution
          = std::make_shared<Solution>(0, reversedPath.front(), nullptr);
      for (size_t idx = 1; idx < reversedPath.size(); idx++) {
        const CostT edgeCost
            = edgeCostBetweenVertices(reversedPath[idx - 1], reversedPath[idx]);
        if (edgeCost >= kInfCost / 4) {
          return std::shared_ptr<Solution>();
        }
        cumulativeCost += edgeCost;
        pathSolution = std::make_shared<Solution>(
            cumulativeCost, reversedPath[idx], pathSolution);
      }
      return pathSolution;
    };

    output.reserve(numPseudoPins - 1);
    for (int pinIndex = 1; pinIndex < numPseudoPins; pinIndex++) {
      if (parent[pinIndex] < 0) {
        return false;
      }
      auto pairSolution = buildSolutionFromPinPair(parent[pinIndex], pinIndex);
      if (!pairSolution) {
        return false;
      }
      output.emplace_back(pairSolution);
    }
    return !output.empty();
  };

  std::vector<std::shared_ptr<Solution>> bestSolutions;
  CandidateScore bestScore;
  for (const int candidateStart : startCandidates) {
    std::vector<std::shared_ptr<Solution>> candidateSolutions;
    if (!runFromStartPin(candidateStart, candidateSolutions)) {
      continue;
    }
    const CandidateScore candidateScore = scoreSolutions(candidateSolutions);
    if (isBetterCandidate(candidateScore, bestScore)) {
      bestScore = candidateScore;
      bestSolutions = std::move(candidateSolutions);
    }
  }
  std::vector<std::shared_ptr<Solution>> mstSolutions;
  if (runMstCandidate(mstSolutions)) {
    const CandidateScore mstScore = scoreSolutions(mstSolutions);
    if (isBetterCandidate(mstScore, bestScore)) {
      bestScore = mstScore;
      bestSolutions = std::move(mstSolutions);
    }
  }

  if (bestSolutions.empty()) {
    logger_->warn(utl::GRT,
                  7002,
                  "sparse maze candidate skipped: failed to connect all pins.");
    return;
  }

  solutions_ = std::move(bestSolutions);
}

std::shared_ptr<SteinerTreeNode> MazeRoute::getSteinerTree() const
{
  std::shared_ptr<SteinerTreeNode> tree = nullptr;
  if (graph_.getNumPseudoPins() == 1) {
    const auto& pseudoPin = graph_.getPseudoPin(0);
    tree = std::make_shared<SteinerTreeNode>(pseudoPin.point, pseudoPin.layers);
    return tree;
  }
  if (solutions_.empty()) {
    return nullptr;
  }

  robin_hood::unordered_map<int, std::shared_ptr<SteinerTreeNode>> created;
  auto mergeVertexPinLayers = [&](int vertex,
                                  const std::shared_ptr<SteinerTreeNode>& node) {
    const auto& pinsAtVertex = graph_.getVertexPins(vertex);
    if (pinsAtVertex.empty()) {
      return;
    }
    bool hasLayers = false;
    IntervalT mergedLayers;
    for (const int pinIndex : pinsAtVertex) {
      const auto& layers = graph_.getPseudoPin(pinIndex).layers;
      if (!hasLayers) {
        mergedLayers = layers;
        hasLayers = true;
      } else {
        mergedLayers.UnionWith(layers);
      }
    }
    if (!hasLayers) {
      return;
    }
    if (node->getFixedLayers().IsValid()) {
      node->getFixedLayers().UnionWith(mergedLayers);
    } else {
      node->setFixedLayers(mergedLayers);
    }
  };
  for (auto& solution : solutions_) {
    std::shared_ptr<Solution> temp = solution;
    std::shared_ptr<SteinerTreeNode> lastNode = nullptr;
    while (temp) {
      auto it = created.find(temp->vertex);
      if (it == created.end()) {
        const PointT point = graph_.getPoint(temp->vertex);
        auto node = std::make_shared<SteinerTreeNode>(point);
        created.emplace(temp->vertex, node);
        if (lastNode) {
          node->addChild(lastNode);
        }
        if (!temp->prev) {
          tree = node;
        }
        if (!lastNode || !temp->prev) {
          // Both the start and the end of a path represent pseudo pins.
          // Merge all colocated pin layer constraints at this vertex.
          mergeVertexPinLayers(temp->vertex, node);
        }
        lastNode = std::move(node);
        temp = temp->prev;
      } else {
        if (!lastNode || !temp->prev) {
          mergeVertexPinLayers(temp->vertex, it->second);
        }
        if (lastNode) {
          it->second->addChild(lastNode);
        }
        break;
      }
    }
  }
  if (!tree) {
    return nullptr;
  }

  // Remove redundant tree nodes
  SteinerTreeNode::preorder(
      tree, [&](const std::shared_ptr<SteinerTreeNode>& node) {
        for (int childIndex = 0; childIndex < node->getNumChildren();
             childIndex++) {
          const std::shared_ptr<SteinerTreeNode> child
              = node->getChildren()[childIndex];
          if (node->x() == child->x() && node->y() == child->y()) {
            for (const auto& gradchild : child->getChildren()) {
              node->addChild(gradchild);
            }
            if (child->getFixedLayers().IsValid()) {
              if (node->getFixedLayers().IsValid()) {
                node->getFixedLayers().UnionWith(child->getFixedLayers());
              } else {
                node->setFixedLayers(child->getFixedLayers());
              }
            }
            node->removeChild(childIndex);
            childIndex -= 1;
          }
        }
      });

  // Remove intermediate tree nodes
  SteinerTreeNode::preorder(
      tree, [&](const std::shared_ptr<SteinerTreeNode>& node) {
        for (std::shared_ptr<SteinerTreeNode>& child : node->getChildren()) {
          const int direction
              = (node->y() == child->y() ? MetalLayer::H : MetalLayer::V);
          std::shared_ptr<SteinerTreeNode> temp = child;
          while (!temp->getFixedLayers().IsValid()
                 && temp->getNumChildren() == 1
                 && (*temp)[1 - direction]
                        == (*(temp->getChildren()[0]))[1 - direction]) {
            temp = temp->getChildren()[0];
          }
          child = std::move(temp);
        }
      });

  // Check duplicate tree nodes
  SteinerTreeNode::preorder(
      tree, [&](const std::shared_ptr<SteinerTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->x() == child->x() && node->y() == child->y()) {
            logger_->error(utl::GRT,
                           7003,
                           "duplicate tree nodes encountered.");
          }
        }
      });
  return tree;
}

}  // namespace grt::newgr
