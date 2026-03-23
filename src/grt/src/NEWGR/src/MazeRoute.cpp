#include "MazeRoute.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
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
  vertex_pin_.clear();
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

  // Reuse current routed trunks as sparse-graph anchors (FastRoute-like
  // topology preservation). This keeps maze search near short existing trees
  // while still enabling congestion repairs.
  const auto& oldTree = net_->getRoutingTree();
  if (oldTree) {
    GRTreeNode::preorder(
        oldTree, [&](const std::shared_ptr<GRTreeNode>& node) {
          pxs.emplace_back(node->x());
          pys.emplace_back(node->y());
        });
  }

  const auto& box = net_->getBoundingBox();
  pxs.emplace_back(std::clamp(box.cx(), 0, xSize - 1));
  pys.emplace_back(std::clamp(box.cy(), 0, ySize - 1));

  // CUGR-style coarse-to-fine corridor search: bound sparse-graph expansion to
  // each net neighborhood to suppress long global detours.
  const int pins = std::max(2, net_->getNumPins());
  const int hp = std::max(1, box.hp());
  int margin = std::clamp(hp / (pins >= 12 ? 5 : 4), 12, 96);
  if (pins <= 3) {
    margin = std::max(margin, 18);
  }
  if (hp >= 180 || pins >= 16) {
    margin += std::max(6, hp / 30);
  }
  margin += std::max(grid.interval.x(), grid.interval.y());
  int xLow = std::max(0, box.lx() - margin);
  int xHigh = std::min(xSize - 1, box.hx() + margin);
  int yLow = std::max(0, box.ly() - margin);
  int yHigh = std::min(ySize - 1, box.hy() + margin);

  if (oldTree) {
    int treeXL = xSize - 1;
    int treeXH = 0;
    int treeYL = ySize - 1;
    int treeYH = 0;
    GRTreeNode::preorder(
        oldTree, [&](const std::shared_ptr<GRTreeNode>& node) {
          treeXL = std::min(treeXL, node->x());
          treeXH = std::max(treeXH, node->x());
          treeYL = std::min(treeYL, node->y());
          treeYH = std::max(treeYH, node->y());
        });
    const int treePadding = std::max(6, margin / 2);
    xLow = std::min(xLow, std::max(0, treeXL - treePadding));
    xHigh = std::max(xHigh, std::min(xSize - 1, treeXH + treePadding));
    yLow = std::min(yLow, std::max(0, treeYL - treePadding));
    yHigh = std::max(yHigh, std::min(ySize - 1, treeYH + treePadding));
  }

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
    vertex_pin_.emplace(u, pinIndex);
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

  std::vector<CostT> minCosts(graph_.getNumVertices(),
                              std::numeric_limits<CostT>::max());

  // lambda to compare solutions
  auto compareSolution = [&](const std::shared_ptr<Solution>& lhs,
                             const std::shared_ptr<Solution>& rhs) {
    return lhs->cost > rhs->cost;
  };

  std::priority_queue<std::shared_ptr<Solution>,
                      std::vector<std::shared_ptr<Solution>>,
                      decltype(compareSolution)>
      queue(compareSolution);

  // lambda to update solution
  auto updateSolution = [&](const std::shared_ptr<Solution>& solution) {
    queue.push(solution);
    if (solution->cost < minCosts[solution->vertex]) {
      minCosts[solution->vertex] = solution->cost;
    }
  };

  solutions_.reserve(numPseudoPins);

  const auto& box = net_->getBoundingBox();
  const PointT center(box.cx(), box.cy());
  int startPinIndex = 0;
  int64_t bestCenterDistance = std::numeric_limits<int64_t>::max();
  for (int pinIndex = 0; pinIndex < numPseudoPins; pinIndex++) {
    const auto& pseudoPin = graph_.getPseudoPin(pinIndex);
    const int64_t distance = std::llabs(static_cast<int64_t>(pseudoPin.point.x())
                                        - center.x())
                             + std::llabs(static_cast<int64_t>(pseudoPin.point.y())
                                          - center.y());
    if (distance < bestCenterDistance) {
      bestCenterDistance = distance;
      startPinIndex = pinIndex;
    }
  }

  std::vector<bool> visited(numPseudoPins, false);
  visited[startPinIndex] = true;
  int numDetached = numPseudoPins - 1;
  updateSolution(std::make_shared<Solution>(
      0, graph_.getPinVertex(startPinIndex), nullptr));

  while (numDetached > 0) {
    std::shared_ptr<Solution> foundSolution;
    int foundPinIndex = 0;
    while (!queue.empty()) {
      auto solution = queue.top();
      queue.pop();
      foundPinIndex = graph_.getVertexPin(solution->vertex);
      if (foundPinIndex != -1 && !visited[foundPinIndex]) {
        foundSolution = std::move(solution);
        break;
      }
      // Pruning
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
        const CostT nextCost
            = solution->cost + graph_.getEdgeCost(solution->vertex, edgeIndex);
        if (nextCost < minCosts[nextVertex]) {
          updateSolution(
              std::make_shared<Solution>(nextCost, nextVertex, solution));
        }
      }
    }

    solutions_.emplace_back(foundSolution);
    assert(foundPinIndex >= 0);
    visited[foundPinIndex] = true;
    numDetached -= 1;

    // Update the cost of the vertices_ on the path
    std::shared_ptr<Solution> temp = std::move(foundSolution);
    while (temp && temp->cost != 0) {
      updateSolution(std::make_shared<Solution>(0, temp->vertex, temp->prev));
      temp = temp->prev;
    }
  }

  if (numDetached != 0) {
    logger_->error(utl::GRT, 7002, "failed to connect all pins.");
  }
}

std::shared_ptr<SteinerTreeNode> MazeRoute::getSteinerTree() const
{
  std::shared_ptr<SteinerTreeNode> tree = nullptr;
  if (graph_.getNumPseudoPins() == 1) {
    const auto& pseudoPin = graph_.getPseudoPin(0);
    tree = std::make_shared<SteinerTreeNode>(pseudoPin.point, pseudoPin.layers);
    return tree;
  }

  std::vector<bool> visited(net_->getNumPins(), false);
  robin_hood::unordered_map<int, std::shared_ptr<SteinerTreeNode>> created;
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
          // Both the start and the end of the path should contain pins
          const int pinIndex = graph_.getVertexPin(temp->vertex);
          assert(pinIndex != -1);
          node->setFixedLayers(graph_.getPseudoPin(pinIndex).layers);
        }
        lastNode = std::move(node);
        temp = temp->prev;
      } else {
        if (lastNode) {
          it->second->addChild(lastNode);
        }
        break;
      }
    }
  }
  assert(tree);

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
