#include "MazeRoute.h"

#include <algorithm>
#include <cassert>
#include <cmath>
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
                       const SparseGrid& grid,
                       const double via_cost_scale)
{
  via_cost_scale_ = std::max(0.0, via_cost_scale);
  // 0. Create pseudo pins
  const auto selectedAccessPoints = grid_graph_->selectAccessPoints(net_);
  pseudo_pins_.reserve(selectedAccessPoints.size());
  for (const auto& selectedPoint : selectedAccessPoints) {
    pseudo_pins_.push_back(selectedPoint);
  }

  // 1. Collect additional routing grid lines
  std::vector<int> pxs;
  std::vector<int> pys;
  pxs.reserve(net_->getNumPins());
  pys.reserve(net_->getNumPins());
  for (const auto& pin : pseudo_pins_) {
    pxs.emplace_back(pin.point.x());
    pys.emplace_back(pin.point.y());
  }
  std::sort(pxs.begin(), pxs.end());
  std::sort(pys.begin(), pys.end());

  // Add Steiner-like anchors so sparse routing can form shorter trunk lines
  // than a pin-only Hanan subset.
  auto addAnchor = [](std::vector<int>& anchors, int value, int upper_bound) {
    if (upper_bound <= 0) {
      return;
    }
    anchors.push_back(std::clamp(value, 0, upper_bound - 1));
  };
  std::vector<int> x_guides = pxs;
  std::vector<int> y_guides = pys;

  const int xSize = grid_graph_->getSize(0);
  const int ySize = grid_graph_->getSize(1);
  if (!pxs.empty()) {
    const int x_min = pxs.front();
    const int x_max = pxs.back();
    const int y_min = pys.front();
    const int y_max = pys.back();
    addAnchor(x_guides, (x_min + x_max) / 2, xSize);
    addAnchor(y_guides, (y_min + y_max) / 2, ySize);

    int sum_x = 0;
    int sum_y = 0;
    for (const int x : pxs) {
      sum_x += x;
    }
    for (const int y : pys) {
      sum_y += y;
    }
    addAnchor(x_guides, sum_x / static_cast<int>(pxs.size()), xSize);
    addAnchor(y_guides, sum_y / static_cast<int>(pys.size()), ySize);

    if (pxs.size() >= 4) {
      const int n = static_cast<int>(pxs.size()) - 1;
      addAnchor(x_guides, pxs[n / 4], xSize);
      addAnchor(y_guides, pys[n / 4], ySize);
      addAnchor(x_guides, pxs[(3 * n) / 4], xSize);
      addAnchor(y_guides, pys[(3 * n) / 4], ySize);
    }
  }

  std::sort(x_guides.begin(), x_guides.end());
  x_guides.erase(std::unique(x_guides.begin(), x_guides.end()), x_guides.end());
  std::sort(y_guides.begin(), y_guides.end());
  y_guides.erase(std::unique(y_guides.begin(), y_guides.end()), y_guides.end());

  xs_.reserve(xSize / grid.interval.x() + x_guides.size());
  ys_.reserve(ySize / grid.interval.y() + y_guides.size());
  for (int i = 0, j = 0; true; i++) {
    const int x = i * grid.interval.x() + grid.offset.x();
    for (; j < x_guides.size() && x_guides[j] <= x; j++) {
      if ((!xs_.empty() && x_guides[j] == xs_.back()) || x_guides[j] == x) {
        continue;
      }
      xs_.emplace_back(x_guides[j]);
    }
    if (x < xSize) {
      xs_.emplace_back(x);
    } else {
      break;
    }
  }
  for (int i = 0, j = 0; true; i++) {
    const int y = i * grid.interval.y() + grid.offset.y();
    for (; j < y_guides.size() && y_guides[j] <= y; j++) {
      if ((!ys_.empty() && y_guides[j] == ys_.back()) || y_guides[j] == y) {
        continue;
      }
      ys_.emplace_back(y_guides[j]);
    }
    if (y < ySize) {
      ys_.emplace_back(y);
    } else {
      break;
    }
  }

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
    costs_[u][2] = costs_[v][2]
                   = via_cost_scale_ * grid_graph_->getUnitViaCost();
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
  struct RunResult
  {
    std::vector<std::shared_ptr<Solution>> solutions;
    CostT total_cost = std::numeric_limits<CostT>::max();
    int via_steps = std::numeric_limits<int>::max();
    uint64_t unique_wirelength = std::numeric_limits<uint64_t>::max();
    int unique_vias = std::numeric_limits<int>::max();
    bool valid = false;
  };

  constexpr CostT kCostEpsilon = static_cast<CostT>(1e-6);
  auto computeRouteGeometry
      = [&](const std::vector<std::shared_ptr<Solution>>& solutions) {
          uint64_t wirelength = 0;
          int vias = 0;
          robin_hood::unordered_set<uint64_t> seen_edges;
          for (const auto& solution : solutions) {
            std::shared_ptr<Solution> temp = solution;
            while (temp && temp->prev) {
              const int a = std::min(temp->vertex, temp->prev->vertex);
              const int b = std::max(temp->vertex, temp->prev->vertex);
              const uint64_t key
                  = (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
              if (seen_edges.emplace(key).second) {
                const auto lhs = graph_.getPoint(temp->vertex);
                const auto rhs = graph_.getPoint(temp->prev->vertex);
                if (lhs.getLayerIdx() == rhs.getLayerIdx()) {
                  const int direction = lhs.y() == rhs.y() ? MetalLayer::H
                                                           : MetalLayer::V;
                  const int lo = std::min(lhs[direction], rhs[direction]);
                  const int hi = std::max(lhs[direction], rhs[direction]);
                  for (int edge_index = lo; edge_index < hi; edge_index++) {
                    wirelength += static_cast<uint64_t>(
                        grid_graph_->getEdgeLength(direction, edge_index));
                  }
                } else {
                  vias += std::abs(lhs.getLayerIdx() - rhs.getLayerIdx());
                }
              }
              temp = temp->prev;
            }
          }
          return std::pair<uint64_t, int>{wirelength, vias};
        };

  auto finalizeResult = [&](RunResult& result) {
    if (!result.valid) {
      return;
    }
    const auto [wirelength, vias] = computeRouteGeometry(result.solutions);
    result.unique_wirelength = wirelength;
    result.unique_vias = vias;
  };

  auto isBetterResult = [&](const RunResult& candidate,
                            const RunResult& current_best) {
    if (!candidate.valid) {
      return false;
    }
    if (!current_best.valid) {
      return true;
    }
    // Strong wirelength-first tie breaking:
    // this is intentionally more aggressive than congestion-cost wins.
    constexpr uint64_t kStrongWireGain = 5;
    constexpr uint64_t kModerateWireGain = 1;
    constexpr double kCostSlackForWireGain = 1.30;
    constexpr double kCostSlackForModerateWireGain = 1.45;
    if (candidate.unique_wirelength + kStrongWireGain
            < current_best.unique_wirelength
        && candidate.total_cost <= current_best.total_cost * kCostSlackForWireGain) {
      return true;
    }
    if (candidate.unique_wirelength + kModerateWireGain
            < current_best.unique_wirelength
        && candidate.total_cost
               <= current_best.total_cost * kCostSlackForModerateWireGain
        && candidate.unique_vias <= current_best.unique_vias + 6) {
      return true;
    }
    if (candidate.unique_wirelength < current_best.unique_wirelength
        && candidate.unique_vias <= current_best.unique_vias + 8
        && candidate.total_cost <= current_best.total_cost * 1.55) {
      return true;
    }
    // Allow pure congestion-cost wins only when geometric quality is preserved.
    if (candidate.total_cost + kCostEpsilon < current_best.total_cost
        && candidate.unique_wirelength <= current_best.unique_wirelength
        && candidate.unique_vias <= current_best.unique_vias + 4) {
      return true;
    }
    if (std::abs(candidate.total_cost - current_best.total_cost) <= kCostEpsilon
        && candidate.unique_wirelength < current_best.unique_wirelength) {
      return true;
    }
    if (std::abs(candidate.total_cost - current_best.total_cost) <= kCostEpsilon
        && candidate.unique_wirelength == current_best.unique_wirelength
        && candidate.unique_vias < current_best.unique_vias) {
      return true;
    }
    if (std::abs(candidate.total_cost - current_best.total_cost) <= kCostEpsilon
        && candidate.unique_wirelength == current_best.unique_wirelength
        && candidate.unique_vias == current_best.unique_vias
        && candidate.via_steps < current_best.via_steps) {
      return true;
    }
    return false;
  };

  auto runFromSeed = [&](const int start_pin_index) -> RunResult {
    RunResult result;
    const int num_vertices = graph_.getNumVertices();
    const int num_pins = graph_.getNumPseudoPins();
    if (num_vertices == 0 || num_pins == 0 || start_pin_index < 0
        || start_pin_index >= num_pins) {
      return result;
    }

    std::vector<CostT> min_costs(num_vertices,
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
      if (solution->cost < min_costs[solution->vertex]) {
        min_costs[solution->vertex] = solution->cost;
      }
    };

    std::vector<bool> visited(num_pins, false);
    const int start_vertex = graph_.getPinVertex(start_pin_index);
    if (start_vertex < 0) {
      return result;
    }
    visited[start_pin_index] = true;
    int num_detached = num_pins - 1;
    updateSolution(std::make_shared<Solution>(0, start_vertex, nullptr));

    result.solutions.reserve(num_pins);
    result.total_cost = 0;
    result.via_steps = 0;
    while (num_detached > 0) {
      std::shared_ptr<Solution> found_solution;
      int found_pin_index = -1;
      while (!queue.empty()) {
        auto solution = queue.top();
        queue.pop();
        found_pin_index = graph_.getVertexPin(solution->vertex);
        if (found_pin_index != -1 && !visited[found_pin_index]) {
          found_solution = std::move(solution);
          break;
        }
        if (solution->cost > min_costs[solution->vertex]) {
          continue;
        }
        for (int edge_index = 0; edge_index < 3; edge_index++) {
          const int next_vertex
              = graph_.getNextVertex(solution->vertex, edge_index);
          if (next_vertex == -1
              || (solution->prev && next_vertex == solution->prev->vertex)) {
            continue;
          }
          const CostT next_cost
              = solution->cost + graph_.getEdgeCost(solution->vertex, edge_index);
          if (next_cost + kCostEpsilon < min_costs[next_vertex]) {
            updateSolution(
                std::make_shared<Solution>(next_cost, next_vertex, solution));
          }
        }
      }

      if (!found_solution || found_pin_index < 0) {
        result.valid = false;
        return result;
      }

      result.total_cost += found_solution->cost;
      std::shared_ptr<Solution> temp = found_solution;
      while (temp && temp->prev) {
        const auto curr_point = graph_.getPoint(temp->vertex);
        const auto prev_point = graph_.getPoint(temp->prev->vertex);
        if (curr_point.getLayerIdx() != prev_point.getLayerIdx()) {
          result.via_steps++;
        }
        temp = temp->prev;
      }

      result.solutions.emplace_back(found_solution);
      visited[found_pin_index] = true;
      num_detached -= 1;

      temp = std::move(found_solution);
      while (temp && temp->cost != 0) {
        updateSolution(std::make_shared<Solution>(0, temp->vertex, temp->prev));
        temp = temp->prev;
      }
    }

    result.valid = true;
    finalizeResult(result);
    return result;
  };

  // Trunk-first variant:
  // connect a distant pin pair first, then grow the tree from the trunk.
  auto runFromPinPair = [&](const int start_pin_index,
                            const int second_pin_index) -> RunResult {
    if (start_pin_index == second_pin_index) {
      return runFromSeed(start_pin_index);
    }

    RunResult result;
    const int num_vertices = graph_.getNumVertices();
    const int num_pins = graph_.getNumPseudoPins();
    if (num_vertices == 0 || num_pins == 0 || start_pin_index < 0
        || start_pin_index >= num_pins || second_pin_index < 0
        || second_pin_index >= num_pins) {
      return result;
    }

    const int source_vertex = graph_.getPinVertex(start_pin_index);
    const int target_vertex = graph_.getPinVertex(second_pin_index);
    if (source_vertex < 0 || source_vertex >= num_vertices || target_vertex < 0
        || target_vertex >= num_vertices) {
      return result;
    }

    std::vector<CostT> min_costs(num_vertices,
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
      if (solution->cost < min_costs[solution->vertex]) {
        min_costs[solution->vertex] = solution->cost;
      }
    };

    // 1) Build a shortest trunk between the chosen pin pair.
    std::vector<CostT> dist(num_vertices, std::numeric_limits<CostT>::max());
    std::vector<int> parent(num_vertices, -1);
    using QueueItem = std::pair<CostT, int>;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>>
        dijkstra_queue;
    dist[source_vertex] = 0;
    dijkstra_queue.emplace(0, source_vertex);

    while (!dijkstra_queue.empty()) {
      const auto [cost, vertex] = dijkstra_queue.top();
      dijkstra_queue.pop();
      if (cost > dist[vertex] + kCostEpsilon) {
        continue;
      }
      if (vertex == target_vertex) {
        break;
      }
      for (int edge_index = 0; edge_index < 3; edge_index++) {
        const int next_vertex = graph_.getNextVertex(vertex, edge_index);
        if (next_vertex < 0) {
          continue;
        }
        const CostT next_cost = cost + graph_.getEdgeCost(vertex, edge_index);
        if (next_cost + kCostEpsilon < dist[next_vertex]
            || (std::abs(next_cost - dist[next_vertex]) <= kCostEpsilon
                && (parent[next_vertex] < 0 || vertex < parent[next_vertex]))) {
          dist[next_vertex] = next_cost;
          parent[next_vertex] = vertex;
          dijkstra_queue.emplace(next_cost, next_vertex);
        }
      }
    }

    if (!std::isfinite(dist[target_vertex])) {
      return result;
    }

    std::vector<int> reverse_path;
    reverse_path.reserve(64);
    int vertex = target_vertex;
    reverse_path.push_back(vertex);
    const int max_steps = num_vertices + 1;
    int steps = 0;
    while (vertex != source_vertex && steps < max_steps) {
      vertex = parent[vertex];
      if (vertex < 0) {
        return RunResult{};
      }
      reverse_path.push_back(vertex);
      steps++;
    }
    if (vertex != source_vertex) {
      return RunResult{};
    }

    std::shared_ptr<Solution> trunk_solution = nullptr;
    for (auto it = reverse_path.rbegin(); it != reverse_path.rend(); ++it) {
      trunk_solution = std::make_shared<Solution>(0, *it, trunk_solution);
    }
    if (!trunk_solution) {
      return RunResult{};
    }

    result.solutions.reserve(num_pins);
    result.solutions.emplace_back(trunk_solution);
    result.total_cost = dist[target_vertex];
    result.via_steps = 0;
    for (int i = 1; i < static_cast<int>(reverse_path.size()); i++) {
      const auto lhs = graph_.getPoint(reverse_path[i - 1]);
      const auto rhs = graph_.getPoint(reverse_path[i]);
      if (lhs.getLayerIdx() != rhs.getLayerIdx()) {
        result.via_steps++;
      }
    }

    std::vector<bool> visited(num_pins, false);
    visited[start_pin_index] = true;
    visited[second_pin_index] = true;
    int num_detached = num_pins - 2;

    // Add all trunk vertices as 0-cost frontier.
    std::shared_ptr<Solution> temp = trunk_solution;
    while (temp) {
      updateSolution(std::make_shared<Solution>(0, temp->vertex, temp->prev));
      temp = temp->prev;
    }

    // 2) Grow the remaining tree from the trunk frontier.
    while (num_detached > 0) {
      std::shared_ptr<Solution> found_solution;
      int found_pin_index = -1;
      while (!queue.empty()) {
        auto solution = queue.top();
        queue.pop();
        found_pin_index = graph_.getVertexPin(solution->vertex);
        if (found_pin_index != -1 && !visited[found_pin_index]) {
          found_solution = std::move(solution);
          break;
        }
        if (solution->cost > min_costs[solution->vertex]) {
          continue;
        }
        for (int edge_index = 0; edge_index < 3; edge_index++) {
          const int next_vertex
              = graph_.getNextVertex(solution->vertex, edge_index);
          if (next_vertex == -1
              || (solution->prev && next_vertex == solution->prev->vertex)) {
            continue;
          }
          const CostT next_cost
              = solution->cost + graph_.getEdgeCost(solution->vertex, edge_index);
          if (next_cost + kCostEpsilon < min_costs[next_vertex]) {
            updateSolution(
                std::make_shared<Solution>(next_cost, next_vertex, solution));
          }
        }
      }

      if (!found_solution || found_pin_index < 0) {
        result.valid = false;
        return result;
      }

      result.total_cost += found_solution->cost;
      temp = found_solution;
      while (temp && temp->prev) {
        const auto curr_point = graph_.getPoint(temp->vertex);
        const auto prev_point = graph_.getPoint(temp->prev->vertex);
        if (curr_point.getLayerIdx() != prev_point.getLayerIdx()) {
          result.via_steps++;
        }
        temp = temp->prev;
      }

      result.solutions.emplace_back(found_solution);
      visited[found_pin_index] = true;
      num_detached -= 1;

      temp = std::move(found_solution);
      while (temp && temp->cost != 0) {
        updateSolution(std::make_shared<Solution>(0, temp->vertex, temp->prev));
        temp = temp->prev;
      }
    }

    result.valid = true;
    finalizeResult(result);
    return result;
  };

  auto runMetricClosureMst = [&]() -> RunResult {
    RunResult result;
    const int num_vertices = graph_.getNumVertices();
    const int num_pins = graph_.getNumPseudoPins();
    if (num_vertices == 0 || num_pins <= 1) {
      return result;
    }

    // Runtime guardrail: use the metric closure only for small/medium nets.
    constexpr int kMaxMetricClosurePins = 56;
    if (num_pins > kMaxMetricClosurePins) {
      return result;
    }

    const CostT inf = std::numeric_limits<CostT>::max();
    std::vector<std::vector<CostT>> pin_pair_cost(
        num_pins, std::vector<CostT>(num_pins, inf));
    std::vector<std::vector<int>> source_parents(
        num_pins, std::vector<int>(num_vertices, -1));

    using QueueItem = std::pair<CostT, int>;
    auto runSingleSource = [&](const int source_pin) -> bool {
      const int source_vertex = graph_.getPinVertex(source_pin);
      if (source_vertex < 0 || source_vertex >= num_vertices) {
        return false;
      }

      std::vector<CostT> dist(num_vertices, inf);
      std::vector<int> parent(num_vertices, -1);
      std::priority_queue<QueueItem,
                          std::vector<QueueItem>,
                          std::greater<QueueItem>>
          queue;

      dist[source_vertex] = 0;
      queue.emplace(0, source_vertex);
      while (!queue.empty()) {
        const auto [cost, vertex] = queue.top();
        queue.pop();
        if (cost > dist[vertex] + kCostEpsilon) {
          continue;
        }
        for (int edge_index = 0; edge_index < 3; edge_index++) {
          const int next_vertex = graph_.getNextVertex(vertex, edge_index);
          if (next_vertex < 0) {
            continue;
          }
          const CostT next_cost = cost + graph_.getEdgeCost(vertex, edge_index);
          if (next_cost + kCostEpsilon < dist[next_vertex]
              || (std::abs(next_cost - dist[next_vertex]) <= kCostEpsilon
                  && (parent[next_vertex] < 0 || vertex < parent[next_vertex]))) {
            dist[next_vertex] = next_cost;
            parent[next_vertex] = vertex;
            queue.emplace(next_cost, next_vertex);
          }
        }
      }

      source_parents[source_pin] = std::move(parent);
      for (int target_pin = 0; target_pin < num_pins; target_pin++) {
        const int target_vertex = graph_.getPinVertex(target_pin);
        if (target_vertex < 0 || target_vertex >= num_vertices) {
          return false;
        }
        pin_pair_cost[source_pin][target_pin] = dist[target_vertex];
      }
      return true;
    };

    for (int source_pin = 0; source_pin < num_pins; source_pin++) {
      if (!runSingleSource(source_pin)) {
        return result;
      }
    }

    std::vector<bool> used(num_pins, false);
    std::vector<CostT> min_edge_cost(num_pins, inf);
    std::vector<int> parent_pin(num_pins, -1);
    min_edge_cost[0] = 0;

    std::vector<std::pair<int, int>> mst_edges;
    mst_edges.reserve(num_pins - 1);
    for (int iter = 0; iter < num_pins; iter++) {
      int next_pin = -1;
      CostT best_cost = inf;
      for (int pin_index = 0; pin_index < num_pins; pin_index++) {
        if (!used[pin_index] && min_edge_cost[pin_index] < best_cost) {
          best_cost = min_edge_cost[pin_index];
          next_pin = pin_index;
        }
      }
      if (next_pin < 0 || !std::isfinite(best_cost)) {
        return result;
      }
      used[next_pin] = true;
      if (parent_pin[next_pin] >= 0) {
        mst_edges.emplace_back(parent_pin[next_pin], next_pin);
      }

      for (int pin_index = 0; pin_index < num_pins; pin_index++) {
        if (used[pin_index]) {
          continue;
        }
        const CostT candidate_cost
            = std::min(pin_pair_cost[next_pin][pin_index],
                       pin_pair_cost[pin_index][next_pin]);
        if (candidate_cost + kCostEpsilon < min_edge_cost[pin_index]
            || (std::abs(candidate_cost - min_edge_cost[pin_index])
                    <= kCostEpsilon
                && (parent_pin[pin_index] < 0 || next_pin < parent_pin[pin_index]))) {
          min_edge_cost[pin_index] = candidate_cost;
          parent_pin[pin_index] = next_pin;
        }
      }
    }

    result.total_cost = 0;
    result.via_steps = 0;
    result.solutions.reserve(mst_edges.size());
    for (const auto& [source_pin, target_pin] : mst_edges) {
      const int source_vertex = graph_.getPinVertex(source_pin);
      const int target_vertex = graph_.getPinVertex(target_pin);
      if (source_vertex < 0 || target_vertex < 0) {
        return RunResult{};
      }
      const auto& parent = source_parents[source_pin];

      std::vector<int> reverse_path;
      reverse_path.reserve(32);
      int vertex = target_vertex;
      reverse_path.push_back(vertex);
      const int max_steps = num_vertices + 1;
      int steps = 0;
      while (vertex != source_vertex && steps < max_steps) {
        vertex = parent[vertex];
        if (vertex < 0) {
          return RunResult{};
        }
        reverse_path.push_back(vertex);
        steps++;
      }
      if (vertex != source_vertex) {
        return RunResult{};
      }

      std::shared_ptr<Solution> chain = nullptr;
      for (auto it = reverse_path.rbegin(); it != reverse_path.rend(); ++it) {
        chain = std::make_shared<Solution>(0, *it, chain);
      }
      if (!chain) {
        return RunResult{};
      }
      result.solutions.emplace_back(chain);

      for (int i = 1; i < reverse_path.size(); i++) {
        const auto lhs = graph_.getPoint(reverse_path[i - 1]);
        const auto rhs = graph_.getPoint(reverse_path[i]);
        if (lhs.getLayerIdx() != rhs.getLayerIdx()) {
          result.via_steps++;
        }
      }

      const CostT edge_cost = pin_pair_cost[source_pin][target_pin];
      if (!std::isfinite(edge_cost)) {
        return RunResult{};
      }
      result.total_cost += edge_cost;
    }
    result.valid = true;
    finalizeResult(result);
    return result;
  };

  // SPRoute/FastRoute hybrid candidate:
  // build a Manhattan MST over pins, then embed each MST edge on the sparse
  // graph using shortest paths. This preserves a low-wirelength topology while
  // still honoring congestion-aware edge costs in path embedding.
  auto runGeometricMstEmbedding = [&](const int root_pin) -> RunResult {
    RunResult result;
    const int num_vertices = graph_.getNumVertices();
    const int num_pins = graph_.getNumPseudoPins();
    if (num_vertices == 0 || num_pins <= 1 || root_pin < 0
        || root_pin >= num_pins) {
      return result;
    }

    constexpr int kMaxGeometricMstPins = 96;
    if (num_pins > kMaxGeometricMstPins) {
      return result;
    }

    std::vector<bool> used(num_pins, false);
    std::vector<int> parent_pin(num_pins, -1);
    std::vector<int> min_manhattan(num_pins, std::numeric_limits<int>::max());
    min_manhattan[root_pin] = 0;

    for (int iter = 0; iter < num_pins; iter++) {
      int next_pin = -1;
      int best_dist = std::numeric_limits<int>::max();
      for (int pin = 0; pin < num_pins; pin++) {
        if (!used[pin] && min_manhattan[pin] < best_dist) {
          best_dist = min_manhattan[pin];
          next_pin = pin;
        }
      }
      if (next_pin < 0) {
        return result;
      }
      used[next_pin] = true;

      const PointT next_point = graph_.getPseudoPin(next_pin).point;
      for (int pin = 0; pin < num_pins; pin++) {
        if (used[pin]) {
          continue;
        }
        const PointT point = graph_.getPseudoPin(pin).point;
        const int manhattan_dist
            = std::abs(next_point.x() - point.x())
              + std::abs(next_point.y() - point.y());
        if (manhattan_dist < min_manhattan[pin]
            || (manhattan_dist == min_manhattan[pin]
                && (parent_pin[pin] < 0 || next_pin < parent_pin[pin]))) {
          min_manhattan[pin] = manhattan_dist;
          parent_pin[pin] = next_pin;
        }
      }
    }

    struct SingleSourceData
    {
      bool ready = false;
      std::vector<CostT> dist;
      std::vector<int> parent;
    };
    std::vector<SingleSourceData> source_cache(num_pins);
    using QueueItem = std::pair<CostT, int>;
    const CostT inf = std::numeric_limits<CostT>::max();
    auto buildSingleSource = [&](const int source_pin) -> bool {
      if (source_pin < 0 || source_pin >= num_pins) {
        return false;
      }
      auto& cache = source_cache[source_pin];
      if (cache.ready) {
        return true;
      }
      const int source_vertex = graph_.getPinVertex(source_pin);
      if (source_vertex < 0 || source_vertex >= num_vertices) {
        return false;
      }

      cache.dist.assign(num_vertices, inf);
      cache.parent.assign(num_vertices, -1);
      std::priority_queue<QueueItem,
                          std::vector<QueueItem>,
                          std::greater<QueueItem>>
          queue;
      cache.dist[source_vertex] = 0;
      queue.emplace(0, source_vertex);
      while (!queue.empty()) {
        const auto [cost, vertex] = queue.top();
        queue.pop();
        if (cost > cache.dist[vertex] + kCostEpsilon) {
          continue;
        }
        for (int edge_index = 0; edge_index < 3; edge_index++) {
          const int next_vertex = graph_.getNextVertex(vertex, edge_index);
          if (next_vertex < 0) {
            continue;
          }
          const CostT next_cost = cost + graph_.getEdgeCost(vertex, edge_index);
          if (next_cost + kCostEpsilon < cache.dist[next_vertex]
              || (std::abs(next_cost - cache.dist[next_vertex]) <= kCostEpsilon
                  && (cache.parent[next_vertex] < 0
                      || vertex < cache.parent[next_vertex]))) {
            cache.dist[next_vertex] = next_cost;
            cache.parent[next_vertex] = vertex;
            queue.emplace(next_cost, next_vertex);
          }
        }
      }
      cache.ready = true;
      return true;
    };

    result.total_cost = 0;
    result.via_steps = 0;
    result.solutions.reserve(num_pins - 1);
    for (int target_pin = 0; target_pin < num_pins; target_pin++) {
      const int source_pin = parent_pin[target_pin];
      if (source_pin < 0) {
        continue;
      }
      if (!buildSingleSource(source_pin)) {
        return RunResult{};
      }
      const auto& cache = source_cache[source_pin];
      const int source_vertex = graph_.getPinVertex(source_pin);
      const int target_vertex = graph_.getPinVertex(target_pin);
      if (source_vertex < 0 || target_vertex < 0 || target_vertex >= num_vertices) {
        return RunResult{};
      }
      if (!std::isfinite(cache.dist[target_vertex])) {
        return RunResult{};
      }

      std::vector<int> reverse_path;
      reverse_path.reserve(32);
      int vertex = target_vertex;
      reverse_path.push_back(vertex);
      const int max_steps = num_vertices + 1;
      int steps = 0;
      while (vertex != source_vertex && steps < max_steps) {
        vertex = cache.parent[vertex];
        if (vertex < 0) {
          return RunResult{};
        }
        reverse_path.push_back(vertex);
        steps++;
      }
      if (vertex != source_vertex) {
        return RunResult{};
      }

      std::shared_ptr<Solution> chain = nullptr;
      for (auto it = reverse_path.rbegin(); it != reverse_path.rend(); ++it) {
        chain = std::make_shared<Solution>(0, *it, chain);
      }
      if (!chain) {
        return RunResult{};
      }
      result.solutions.emplace_back(chain);
      result.total_cost += cache.dist[target_vertex];

      for (int i = 1; i < static_cast<int>(reverse_path.size()); i++) {
        const auto lhs = graph_.getPoint(reverse_path[i - 1]);
        const auto rhs = graph_.getPoint(reverse_path[i]);
        if (lhs.getLayerIdx() != rhs.getLayerIdx()) {
          result.via_steps++;
        }
      }
    }
    result.valid = true;
    finalizeResult(result);
    return result;
  };

  solutions_.clear();
  const int num_pins = graph_.getNumPseudoPins();
  if (num_pins <= 1) {
    return;
  }

  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int min_y = std::numeric_limits<int>::max();
  int max_y = std::numeric_limits<int>::min();
  for (int pin_index = 0; pin_index < num_pins; pin_index++) {
    const PointT point = graph_.getPseudoPin(pin_index).point;
    min_x = std::min(min_x, point.x());
    max_x = std::max(max_x, point.x());
    min_y = std::min(min_y, point.y());
    max_y = std::max(max_y, point.y());
  }
  const PointT center((min_x + max_x) / 2, (min_y + max_y) / 2);

  int center_seed = 0;
  int far_seed = 0;
  int min_x_seed = 0;
  int max_x_seed = 0;
  int min_y_seed = 0;
  int max_y_seed = 0;
  int min_sum_seed = 0;
  int max_sum_seed = 0;
  int min_diff_seed = 0;
  int max_diff_seed = 0;
  int far_pair_lhs = -1;
  int far_pair_rhs = -1;
  int best_center_dist = std::numeric_limits<int>::max();
  int best_far_dist = std::numeric_limits<int>::min();
  int min_sum = std::numeric_limits<int>::max();
  int max_sum = std::numeric_limits<int>::min();
  int min_diff = std::numeric_limits<int>::max();
  int max_diff = std::numeric_limits<int>::min();
  for (int pin_index = 0; pin_index < num_pins; pin_index++) {
    const PointT point = graph_.getPseudoPin(pin_index).point;
    const int center_dist
        = std::abs(point.x() - center.x()) + std::abs(point.y() - center.y());
    if (center_dist < best_center_dist) {
      best_center_dist = center_dist;
      center_seed = pin_index;
    }
    if (center_dist > best_far_dist) {
      best_far_dist = center_dist;
      far_seed = pin_index;
    }
    if (point.x() < graph_.getPseudoPin(min_x_seed).point.x()) {
      min_x_seed = pin_index;
    }
    if (point.x() > graph_.getPseudoPin(max_x_seed).point.x()) {
      max_x_seed = pin_index;
    }
    if (point.y() < graph_.getPseudoPin(min_y_seed).point.y()) {
      min_y_seed = pin_index;
    }
    if (point.y() > graph_.getPseudoPin(max_y_seed).point.y()) {
      max_y_seed = pin_index;
    }
    const int sum = point.x() + point.y();
    if (sum < min_sum) {
      min_sum = sum;
      min_sum_seed = pin_index;
    }
    if (sum > max_sum) {
      max_sum = sum;
      max_sum_seed = pin_index;
    }
    const int diff = point.x() - point.y();
    if (diff < min_diff) {
      min_diff = diff;
      min_diff_seed = pin_index;
    }
    if (diff > max_diff) {
      max_diff = diff;
      max_diff_seed = pin_index;
    }
  }
  std::vector<std::pair<int, int>> far_pairs;
  if (num_pins <= 96) {
    struct PairCandidate
    {
      int lhs;
      int rhs;
      int dist;
    };
    std::vector<PairCandidate> ranked_pairs;
    ranked_pairs.reserve((num_pins * (num_pins - 1)) / 2);
    for (int lhs = 0; lhs < num_pins; lhs++) {
      const PointT lhs_point = graph_.getPseudoPin(lhs).point;
      for (int rhs = lhs + 1; rhs < num_pins; rhs++) {
        const PointT rhs_point = graph_.getPseudoPin(rhs).point;
        const int dist = std::abs(lhs_point.x() - rhs_point.x())
                         + std::abs(lhs_point.y() - rhs_point.y());
        ranked_pairs.push_back({lhs, rhs, dist});
      }
    }
    std::sort(ranked_pairs.begin(),
              ranked_pairs.end(),
              [](const PairCandidate& lhs, const PairCandidate& rhs) {
                if (lhs.dist != rhs.dist) {
                  return lhs.dist > rhs.dist;
                }
                if (lhs.lhs != rhs.lhs) {
                  return lhs.lhs < rhs.lhs;
                }
                return lhs.rhs < rhs.rhs;
              });

    int far_pair_budget = 3;
    if (num_pins <= 16) {
      far_pair_budget = 10;
    } else if (num_pins <= 32) {
      far_pair_budget = 8;
    } else if (num_pins <= 64) {
      far_pair_budget = 6;
    }

    std::vector<bool> endpoint_used(num_pins, false);
    for (const auto& pair : ranked_pairs) {
      if (static_cast<int>(far_pairs.size()) >= far_pair_budget) {
        break;
      }
      if (!endpoint_used[pair.lhs] && !endpoint_used[pair.rhs]) {
        far_pairs.emplace_back(pair.lhs, pair.rhs);
        endpoint_used[pair.lhs] = true;
        endpoint_used[pair.rhs] = true;
      }
    }
    for (const auto& pair : ranked_pairs) {
      if (static_cast<int>(far_pairs.size()) >= far_pair_budget) {
        break;
      }
      far_pairs.emplace_back(pair.lhs, pair.rhs);
    }

    if (!far_pairs.empty()) {
      far_pair_lhs = far_pairs.front().first;
      far_pair_rhs = far_pairs.front().second;
    }
  }

  std::vector<int> seeds;
  seeds.reserve(6);
  auto addSeed = [&](const int seed) {
    if (seed < 0 || seed >= num_pins) {
      return;
    }
    if (std::find(seeds.begin(), seeds.end(), seed) == seeds.end()) {
      seeds.push_back(seed);
    }
  };
  addSeed(center_seed);
  addSeed(far_seed);
  addSeed(min_x_seed);
  addSeed(max_x_seed);
  addSeed(min_y_seed);
  addSeed(max_y_seed);
  addSeed(min_sum_seed);
  addSeed(max_sum_seed);
  addSeed(min_diff_seed);
  addSeed(max_diff_seed);

  std::vector<int> sorted_by_center;
  sorted_by_center.reserve(num_pins);
  for (int pin_index = 0; pin_index < num_pins; pin_index++) {
    sorted_by_center.push_back(pin_index);
  }
  std::sort(sorted_by_center.begin(),
            sorted_by_center.end(),
            [&](const int lhs, const int rhs) {
              const PointT lhs_point = graph_.getPseudoPin(lhs).point;
              const PointT rhs_point = graph_.getPseudoPin(rhs).point;
              const int lhs_dist = std::abs(lhs_point.x() - center.x())
                                   + std::abs(lhs_point.y() - center.y());
              const int rhs_dist = std::abs(rhs_point.x() - center.x())
                                   + std::abs(rhs_point.y() - center.y());
              if (lhs_dist != rhs_dist) {
                return lhs_dist < rhs_dist;
              }
              return lhs < rhs;
            });
  for (const int seed : sorted_by_center) {
    addSeed(seed);
  }

  std::vector<int> sorted_by_far = sorted_by_center;
  std::reverse(sorted_by_far.begin(), sorted_by_far.end());
  for (const int seed : sorted_by_far) {
    addSeed(seed);
  }

  int max_seeds = 1;
  if (num_pins <= 16) {
    max_seeds = num_pins;
  } else if (num_pins <= 32) {
    max_seeds = 16;
  } else if (num_pins <= 64) {
    max_seeds = 14;
  } else if (num_pins <= 96) {
    max_seeds = 10;
  } else {
    max_seeds = 6;
  }
  if (max_seeds < static_cast<int>(seeds.size())) {
    seeds.resize(max_seeds);
  }

  RunResult best_result = runMetricClosureMst();
  auto runGeometricCandidate = [&](const int root_pin) {
    RunResult geometric_mst_result = runGeometricMstEmbedding(root_pin);
    if (isBetterResult(geometric_mst_result, best_result)) {
      best_result = std::move(geometric_mst_result);
    }
  };
  runGeometricCandidate(center_seed);
  if (num_pins <= 128) {
    runGeometricCandidate(far_seed);
  }
  if (num_pins <= 64) {
    runGeometricCandidate(min_x_seed);
    runGeometricCandidate(max_x_seed);
    runGeometricCandidate(min_y_seed);
    runGeometricCandidate(max_y_seed);
    runGeometricCandidate(min_sum_seed);
    runGeometricCandidate(max_sum_seed);
    runGeometricCandidate(min_diff_seed);
    runGeometricCandidate(max_diff_seed);
  }

  auto runPairCandidate = [&](const int lhs, const int rhs) {
    if (lhs < 0 || rhs < 0 || lhs >= num_pins || rhs >= num_pins || lhs == rhs) {
      return;
    }
    RunResult pair_result = runFromPinPair(lhs, rhs);
    if (isBetterResult(pair_result, best_result)) {
      best_result = std::move(pair_result);
    }
  };
  runPairCandidate(min_x_seed, max_x_seed);
  runPairCandidate(min_y_seed, max_y_seed);
  runPairCandidate(center_seed, far_seed);
  if (num_pins <= 128 && far_pair_lhs >= 0 && far_pair_rhs >= 0) {
    runPairCandidate(far_pair_lhs, far_pair_rhs);
  }
  if (num_pins <= 96) {
    for (const auto& pair : far_pairs) {
      runPairCandidate(pair.first, pair.second);
    }
  }
  if (num_pins <= 64) {
    runPairCandidate(min_x_seed, max_y_seed);
    runPairCandidate(max_x_seed, min_y_seed);
    runPairCandidate(min_sum_seed, max_sum_seed);
    runPairCandidate(min_diff_seed, max_diff_seed);
    runPairCandidate(center_seed, min_x_seed);
    runPairCandidate(center_seed, max_x_seed);
  }

  for (const int seed : seeds) {
    RunResult candidate = runFromSeed(seed);
    if (isBetterResult(candidate, best_result)) {
      best_result = std::move(candidate);
    }
  }

  if (!best_result.valid) {
    logger_->error(utl::GRT, 7002, "failed to connect all pins.");
    return;
  }

  solutions_ = std::move(best_result.solutions);
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
