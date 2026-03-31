#include "CUGR.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
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

int getMedian(std::vector<int> values)
{
  if (values.empty()) {
    return 0;
  }
  const size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  return values[middle];
}

void mergeFixedLayers(const IntervalT& layers,
                      std::shared_ptr<SteinerTreeNode>& node)
{
  if (!layers.IsValid()) {
    return;
  }
  if (node->getFixedLayers().IsValid()) {
    node->setFixedLayers(node->getFixedLayers().UnionWith(layers));
  } else {
    node->setFixedLayers(layers);
  }
}

std::vector<AccessPoint> materializeAccessPoints(
    const AccessPointSet& access_points)
{
  std::vector<AccessPoint> points;
  points.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    points.push_back(access_point);
  }
  std::sort(points.begin(),
            points.end(),
            [](const AccessPoint& lhs, const AccessPoint& rhs) {
              if (lhs.point.x() != rhs.point.x()) {
                return lhs.point.x() < rhs.point.x();
              }
              if (lhs.point.y() != rhs.point.y()) {
                return lhs.point.y() < rhs.point.y();
              }
              if (lhs.layers.low() != rhs.layers.low()) {
                return lhs.layers.low() < rhs.layers.low();
              }
              return lhs.layers.high() < rhs.layers.high();
            });
  return points;
}

std::shared_ptr<SteinerTreeNode> buildMedianTrunkSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const bool horizontal_trunk)
{
  if (access_points.empty()) {
    return nullptr;
  }
  if (access_points.size() == 1) {
    return std::make_shared<SteinerTreeNode>(access_points.front().point,
                                             access_points.front().layers);
  }

  std::vector<int> primary_coords;
  std::vector<int> secondary_coords;
  primary_coords.reserve(access_points.size());
  secondary_coords.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    if (horizontal_trunk) {
      primary_coords.push_back(access_point.point.x());
      secondary_coords.push_back(access_point.point.y());
    } else {
      primary_coords.push_back(access_point.point.y());
      secondary_coords.push_back(access_point.point.x());
    }
  }

  const int trunk_axis = getMedian(secondary_coords);
  const int root_coord = getMedian(primary_coords);
  primary_coords.push_back(root_coord);
  std::sort(primary_coords.begin(), primary_coords.end());
  primary_coords.erase(
      std::unique(primary_coords.begin(), primary_coords.end()),
      primary_coords.end());

  if (primary_coords.empty()) {
    return nullptr;
  }

  const int root_index
      = std::distance(primary_coords.begin(),
                      std::find(primary_coords.begin(),
                                primary_coords.end(),
                                root_coord));
  std::vector<std::shared_ptr<SteinerTreeNode>> trunk_nodes;
  trunk_nodes.reserve(primary_coords.size());
  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> anchor_nodes;
  anchor_nodes.reserve(primary_coords.size());
  for (const int coordinate : primary_coords) {
    const PointT loc
        = horizontal_trunk ? PointT(coordinate, trunk_axis)
                           : PointT(trunk_axis, coordinate);
    std::shared_ptr<SteinerTreeNode> node
        = std::make_shared<SteinerTreeNode>(loc);
    anchor_nodes.emplace(coordinate, node);
    trunk_nodes.push_back(node);
  }

  const int trunk_count = static_cast<int>(trunk_nodes.size());
  for (int index = root_index - 1; index >= 0; index--) {
    trunk_nodes[index + 1]->addChild(trunk_nodes[index]);
  }
  for (int index = root_index + 1; index < trunk_count; index++) {
    trunk_nodes[index - 1]->addChild(trunk_nodes[index]);
  }

  for (const AccessPoint& access_point : access_points) {
    const int coordinate
        = horizontal_trunk ? access_point.point.x() : access_point.point.y();
    auto it = anchor_nodes.find(coordinate);
    if (it == anchor_nodes.end()) {
      continue;
    }
    std::shared_ptr<SteinerTreeNode>& anchor_node = it->second;
    const bool on_trunk
        = horizontal_trunk ? access_point.point.y() == trunk_axis
                           : access_point.point.x() == trunk_axis;
    if (on_trunk) {
      mergeFixedLayers(access_point.layers, anchor_node);
    } else {
      anchor_node->addChild(
          std::make_shared<SteinerTreeNode>(access_point.point,
                                            access_point.layers));
    }
  }

  return trunk_nodes[root_index];
}

uint64_t pointKey(const PointT& point)
{
  return (static_cast<uint64_t>(static_cast<uint32_t>(point.x())) << 32)
         | static_cast<uint32_t>(point.y());
}

void attachChildUnique(const std::shared_ptr<SteinerTreeNode>& parent,
                       const std::shared_ptr<SteinerTreeNode>& child)
{
  if (!parent || !child || parent == child) {
    return;
  }
  for (const auto& existing : parent->getChildren()) {
    if (existing == child) {
      return;
    }
  }
  parent->addChild(child);
}

std::shared_ptr<SteinerTreeNode> buildHubSpineSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const PointT& hub_point,
    const bool horizontal_trunk)
{
  if (access_points.empty()) {
    return nullptr;
  }
  if (access_points.size() == 1) {
    return std::make_shared<SteinerTreeNode>(access_points.front().point,
                                             access_points.front().layers);
  }

  std::vector<int> primary_coords;
  primary_coords.reserve(access_points.size() + 1);
  for (const AccessPoint& access_point : access_points) {
    primary_coords.push_back(horizontal_trunk ? access_point.point.x()
                                              : access_point.point.y());
  }
  const int root_coord = horizontal_trunk ? hub_point.x() : hub_point.y();
  const int trunk_axis = horizontal_trunk ? hub_point.y() : hub_point.x();
  primary_coords.push_back(root_coord);
  std::sort(primary_coords.begin(), primary_coords.end());
  primary_coords.erase(
      std::unique(primary_coords.begin(), primary_coords.end()),
      primary_coords.end());

  if (primary_coords.empty()) {
    return nullptr;
  }

  const int root_index
      = std::distance(primary_coords.begin(),
                      std::find(primary_coords.begin(),
                                primary_coords.end(),
                                root_coord));
  std::vector<std::shared_ptr<SteinerTreeNode>> trunk_nodes;
  trunk_nodes.reserve(primary_coords.size());
  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> anchor_nodes;
  anchor_nodes.reserve(primary_coords.size());
  for (const int coordinate : primary_coords) {
    const PointT loc
        = horizontal_trunk ? PointT(coordinate, trunk_axis)
                           : PointT(trunk_axis, coordinate);
    std::shared_ptr<SteinerTreeNode> node
        = std::make_shared<SteinerTreeNode>(loc);
    anchor_nodes.emplace(coordinate, node);
    trunk_nodes.push_back(node);
  }

  const int trunk_count = static_cast<int>(trunk_nodes.size());
  for (int index = root_index - 1; index >= 0; index--) {
    trunk_nodes[index + 1]->addChild(trunk_nodes[index]);
  }
  for (int index = root_index + 1; index < trunk_count; index++) {
    trunk_nodes[index - 1]->addChild(trunk_nodes[index]);
  }

  std::unordered_set<uint64_t> emitted_stems;
  emitted_stems.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    const int coordinate
        = horizontal_trunk ? access_point.point.x() : access_point.point.y();
    auto it = anchor_nodes.find(coordinate);
    if (it == anchor_nodes.end()) {
      continue;
    }
    std::shared_ptr<SteinerTreeNode>& anchor_node = it->second;
    const bool on_trunk
        = horizontal_trunk ? access_point.point.y() == trunk_axis
                           : access_point.point.x() == trunk_axis;
    if (on_trunk) {
      mergeFixedLayers(access_point.layers, anchor_node);
      continue;
    }
    const uint64_t key = pointKey(access_point.point);
    if (!emitted_stems.insert(key).second) {
      continue;
    }
    anchor_node->addChild(
        std::make_shared<SteinerTreeNode>(access_point.point,
                                          access_point.layers));
  }

  return trunk_nodes[root_index];
}

std::vector<PointT> buildHubCandidates(const std::vector<AccessPoint>& access_points,
                                       const BoxT& bbox,
                                       const int max_candidates)
{
  std::vector<PointT> candidates;
  if (access_points.empty() || max_candidates <= 0) {
    return candidates;
  }

  std::vector<int> xs;
  std::vector<int> ys;
  xs.reserve(access_points.size());
  ys.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    xs.push_back(access_point.point.x());
    ys.push_back(access_point.point.y());
  }

  std::vector<int> xs_sorted = xs;
  std::vector<int> ys_sorted = ys;
  std::sort(xs_sorted.begin(), xs_sorted.end());
  std::sort(ys_sorted.begin(), ys_sorted.end());
  const int median_x = getMedian(xs);
  const int median_y = getMedian(ys);
  const int q1_x = xs_sorted[xs_sorted.size() / 4];
  const int q3_x = xs_sorted[(3 * xs_sorted.size()) / 4];
  const int q1_y = ys_sorted[ys_sorted.size() / 4];
  const int q3_y = ys_sorted[(3 * ys_sorted.size()) / 4];

  std::unordered_set<uint64_t> seen;
  seen.reserve(max_candidates * 2);
  auto pushCandidate = [&](PointT point) {
    point[0] = std::clamp(point.x(), bbox.lx(), bbox.hx());
    point[1] = std::clamp(point.y(), bbox.ly(), bbox.hy());
    const uint64_t key = pointKey(point);
    if (seen.insert(key).second) {
      candidates.push_back(point);
    }
  };

  pushCandidate({median_x, median_y});
  pushCandidate({bbox.cx(), bbox.cy()});
  pushCandidate({q1_x, median_y});
  pushCandidate({q3_x, median_y});
  pushCandidate({median_x, q1_y});
  pushCandidate({median_x, q3_y});

  struct FarPin
  {
    int distance;
    PointT point;
  };
  std::vector<FarPin> far_pins;
  far_pins.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    const int distance = std::abs(access_point.point.x() - median_x)
                         + std::abs(access_point.point.y() - median_y);
    far_pins.push_back({distance, access_point.point});
  }
  std::sort(far_pins.begin(),
            far_pins.end(),
            [](const FarPin& lhs, const FarPin& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance > rhs.distance;
              }
              if (lhs.point.x() != rhs.point.x()) {
                return lhs.point.x() < rhs.point.x();
              }
              return lhs.point.y() < rhs.point.y();
            });

  const int far_limit = std::min(6, static_cast<int>(far_pins.size()));
  for (int idx = 0; idx < far_limit; idx++) {
    pushCandidate(far_pins[idx].point);
    pushCandidate({far_pins[idx].point.x(), median_y});
    pushCandidate({median_x, far_pins[idx].point.y()});
  }

  if (static_cast<int>(candidates.size()) > max_candidates) {
    candidates.resize(max_candidates);
  }
  return candidates;
}

int manhattanDistance(const PointT& lhs, const PointT& rhs)
{
  return std::abs(lhs.x() - rhs.x()) + std::abs(lhs.y() - rhs.y());
}

struct UnionFind
{
  explicit UnionFind(const int size) : parent(size), rank(size, 0)
  {
    for (int i = 0; i < size; i++) {
      parent[i] = i;
    }
  }

  int find(int value)
  {
    if (parent[value] == value) {
      return value;
    }
    parent[value] = find(parent[value]);
    return parent[value];
  }

  bool unite(const int lhs, const int rhs)
  {
    int lhs_root = find(lhs);
    int rhs_root = find(rhs);
    if (lhs_root == rhs_root) {
      return false;
    }
    if (rank[lhs_root] < rank[rhs_root]) {
      std::swap(lhs_root, rhs_root);
    }
    parent[rhs_root] = lhs_root;
    if (rank[lhs_root] == rank[rhs_root]) {
      rank[lhs_root]++;
    }
    return true;
  }

  std::vector<int> parent;
  std::vector<uint8_t> rank;
};

struct MstVertexSpec
{
  PointT point;
  IntervalT layers;
  bool has_layers{false};
};

void mergeVertexLayers(MstVertexSpec& vertex, const IntervalT& layers)
{
  if (!layers.IsValid()) {
    return;
  }
  if (vertex.has_layers) {
    vertex.layers = vertex.layers.UnionWith(layers);
  } else {
    vertex.layers = layers;
    vertex.has_layers = true;
  }
}

std::shared_ptr<SteinerTreeNode> buildRectilinearMstSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const PointT& root_hint,
    const bool include_hub,
    const PointT& hub_point,
    const bool prefer_horizontal_first)
{
  if (access_points.empty()) {
    return nullptr;
  }

  std::vector<MstVertexSpec> vertices;
  vertices.reserve(access_points.size() + 1);
  std::unordered_map<uint64_t, int> point_to_index;
  point_to_index.reserve(access_points.size() + 1);

  auto addVertex = [&](const PointT& point, const IntervalT& layers) {
    const uint64_t key = pointKey(point);
    auto existing = point_to_index.find(key);
    if (existing != point_to_index.end()) {
      mergeVertexLayers(vertices[existing->second], layers);
      return existing->second;
    }
    MstVertexSpec vertex;
    vertex.point = point;
    mergeVertexLayers(vertex, layers);
    const int new_index = vertices.size();
    vertices.push_back(vertex);
    point_to_index.emplace(key, new_index);
    return new_index;
  };

  for (const AccessPoint& access_point : access_points) {
    addVertex(access_point.point, access_point.layers);
  }

  int root_index = 0;
  if (include_hub) {
    root_index = addVertex(hub_point, {});
  } else {
    int best_dist = std::numeric_limits<int>::max();
    for (int i = 0; i < static_cast<int>(vertices.size()); i++) {
      const int dist = manhattanDistance(vertices[i].point, root_hint);
      if (dist < best_dist) {
        best_dist = dist;
        root_index = i;
      }
    }
  }

  if (vertices.size() == 1) {
    if (vertices[0].has_layers) {
      return std::make_shared<SteinerTreeNode>(vertices[0].point,
                                               vertices[0].layers);
    }
    return std::make_shared<SteinerTreeNode>(vertices[0].point);
  }

  struct CandidateEdge
  {
    int u;
    int v;
    int distance;
    int tie_breaker;
  };

  std::vector<CandidateEdge> edges;
  edges.reserve((vertices.size() * (vertices.size() - 1)) / 2);
  for (int lhs = 0; lhs < static_cast<int>(vertices.size()); lhs++) {
    for (int rhs = lhs + 1; rhs < static_cast<int>(vertices.size()); rhs++) {
      const int distance
          = manhattanDistance(vertices[lhs].point, vertices[rhs].point);
      const int tie_breaker
          = std::min(manhattanDistance(vertices[lhs].point, root_hint),
                     manhattanDistance(vertices[rhs].point, root_hint));
      edges.push_back({lhs, rhs, distance, tie_breaker});
    }
  }
  std::sort(edges.begin(),
            edges.end(),
            [](const CandidateEdge& lhs, const CandidateEdge& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
              }
              if (lhs.tie_breaker != rhs.tie_breaker) {
                return lhs.tie_breaker < rhs.tie_breaker;
              }
              if (lhs.u != rhs.u) {
                return lhs.u < rhs.u;
              }
              return lhs.v < rhs.v;
            });

  UnionFind union_find(vertices.size());
  std::vector<std::vector<int>> adjacency(vertices.size());
  int selected_edges = 0;
  for (const CandidateEdge& edge : edges) {
    if (!union_find.unite(edge.u, edge.v)) {
      continue;
    }
    adjacency[edge.u].push_back(edge.v);
    adjacency[edge.v].push_back(edge.u);
    selected_edges++;
    if (selected_edges + 1 == static_cast<int>(vertices.size())) {
      break;
    }
  }
  if (selected_edges + 1 != static_cast<int>(vertices.size())) {
    return nullptr;
  }

  std::vector<std::shared_ptr<SteinerTreeNode>> mst_nodes(vertices.size(),
                                                           nullptr);
  auto getOrCreateNode = [&](const int index) {
    if (!mst_nodes[index]) {
      if (vertices[index].has_layers) {
        mst_nodes[index]
            = std::make_shared<SteinerTreeNode>(vertices[index].point,
                                                vertices[index].layers);
      } else {
        mst_nodes[index] = std::make_shared<SteinerTreeNode>(vertices[index].point);
      }
    }
    return mst_nodes[index];
  };

  std::function<void(int, int)> buildTree = [&](const int current,
                                                 const int parent) {
    std::shared_ptr<SteinerTreeNode> current_node = getOrCreateNode(current);
    for (const int next : adjacency[current]) {
      if (next == parent) {
        continue;
      }
      std::shared_ptr<SteinerTreeNode> next_node = getOrCreateNode(next);
      const PointT current_point = vertices[current].point;
      const PointT next_point = vertices[next].point;

      if (current_point.x() == next_point.x() || current_point.y() == next_point.y()) {
        current_node->addChild(next_node);
      } else {
        const PointT bend_horizontal_first(next_point.x(), current_point.y());
        const PointT bend_vertical_first(current_point.x(), next_point.y());
        const int horizontal_score
            = manhattanDistance(bend_horizontal_first, root_hint);
        const int vertical_score
            = manhattanDistance(bend_vertical_first, root_hint);
        const bool choose_horizontal
            = horizontal_score < vertical_score
              || (horizontal_score == vertical_score && prefer_horizontal_first);
        const PointT bend = choose_horizontal ? bend_horizontal_first
                                              : bend_vertical_first;
        std::shared_ptr<SteinerTreeNode> bend_node
            = std::make_shared<SteinerTreeNode>(bend);
        current_node->addChild(bend_node);
        bend_node->addChild(next_node);
      }

      buildTree(next, current);
    }
  };

  std::shared_ptr<SteinerTreeNode> root = getOrCreateNode(root_index);
  buildTree(root_index, -1);
  return root;
}

std::shared_ptr<SteinerTreeNode> buildDualHubBackboneSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const PointT& hub_a,
    const PointT& hub_b,
    const bool horizontal_first)
{
  if (access_points.empty()) {
    return nullptr;
  }

  std::shared_ptr<SteinerTreeNode> hub_a_node
      = std::make_shared<SteinerTreeNode>(hub_a);
  std::shared_ptr<SteinerTreeNode> hub_b_node = hub_a_node;

  const bool same_hub = hub_a.x() == hub_b.x() && hub_a.y() == hub_b.y();
  if (!same_hub) {
    if (hub_a.x() == hub_b.x() || hub_a.y() == hub_b.y()) {
      hub_b_node = std::make_shared<SteinerTreeNode>(hub_b);
      hub_a_node->addChild(hub_b_node);
    } else {
      const PointT pivot
          = horizontal_first ? PointT(hub_b.x(), hub_a.y())
                             : PointT(hub_a.x(), hub_b.y());
      std::shared_ptr<SteinerTreeNode> pivot_node
          = std::make_shared<SteinerTreeNode>(pivot);
      hub_b_node = std::make_shared<SteinerTreeNode>(hub_b);
      hub_a_node->addChild(pivot_node);
      pivot_node->addChild(hub_b_node);
    }
  }

  for (const AccessPoint& access_point : access_points) {
    if (access_point.point == hub_a) {
      mergeFixedLayers(access_point.layers, hub_a_node);
      continue;
    }
    if (access_point.point == hub_b) {
      mergeFixedLayers(access_point.layers, hub_b_node);
      continue;
    }

    const int dist_to_a = manhattanDistance(access_point.point, hub_a);
    const int dist_to_b = manhattanDistance(access_point.point, hub_b);
    const std::shared_ptr<SteinerTreeNode> attach_node
        = dist_to_a <= dist_to_b ? hub_a_node : hub_b_node;
    attach_node->addChild(std::make_shared<SteinerTreeNode>(access_point.point,
                                                             access_point.layers));
  }

  return hub_a_node;
}

std::vector<int> buildAxisCandidates(const std::vector<AccessPoint>& access_points,
                                     const BoxT& bbox,
                                     const bool use_x_axis,
                                     const int max_candidates)
{
  std::vector<int> candidates;
  if (access_points.empty() || max_candidates <= 0) {
    return candidates;
  }

  const int lower = use_x_axis ? bbox.lx() : bbox.ly();
  const int upper = use_x_axis ? bbox.hx() : bbox.hy();

  std::vector<int> coords;
  coords.reserve(access_points.size());
  std::unordered_map<int, int> frequency;
  frequency.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    const int coordinate = use_x_axis ? access_point.point.x()
                                      : access_point.point.y();
    coords.push_back(coordinate);
    frequency[coordinate]++;
  }

  if (coords.empty()) {
    return candidates;
  }

  std::vector<int> sorted = coords;
  std::sort(sorted.begin(), sorted.end());
  const int median = getMedian(coords);
  const int q1 = sorted[sorted.size() / 4];
  const int q3 = sorted[(3 * sorted.size()) / 4];
  const int center = use_x_axis ? bbox.cx() : bbox.cy();
  const int span = std::max(1, upper - lower);

  auto pushCandidate = [&](const int coordinate) {
    candidates.push_back(std::clamp(coordinate, lower, upper));
  };

  pushCandidate(median);
  pushCandidate(center);
  pushCandidate(q1);
  pushCandidate(q3);
  pushCandidate(lower);
  pushCandidate(upper);
  pushCandidate(lower + span / 3);
  pushCandidate(upper - span / 3);

  std::vector<std::pair<int, int>> ranked_frequency;
  ranked_frequency.reserve(frequency.size());
  for (const auto& [coordinate, count] : frequency) {
    ranked_frequency.push_back({coordinate, count});
  }
  std::sort(
      ranked_frequency.begin(),
      ranked_frequency.end(),
      [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) {
        if (lhs.second != rhs.second) {
          return lhs.second > rhs.second;
        }
        return lhs.first < rhs.first;
      });
  for (int i = 0; i < static_cast<int>(ranked_frequency.size()) && i < 4; i++) {
    pushCandidate(ranked_frequency[i].first);
  }

  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  if (!sorted.empty()) {
    pushCandidate(sorted.front());
    pushCandidate(sorted.back());
    pushCandidate(sorted[sorted.size() / 2]);
  }

  std::vector<int> unique_candidates;
  unique_candidates.reserve(candidates.size());
  std::unordered_set<int> seen;
  seen.reserve(candidates.size());
  for (const int coordinate : candidates) {
    if (seen.insert(coordinate).second) {
      unique_candidates.push_back(coordinate);
    }
  }
  if (static_cast<int>(unique_candidates.size()) > max_candidates) {
    unique_candidates.resize(max_candidates);
  }
  return unique_candidates;
}

std::shared_ptr<SteinerTreeNode> buildCrossbarSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const PointT& spine_origin,
    const bool prefer_vertical)
{
  if (access_points.empty()) {
    return nullptr;
  }
  if (access_points.size() == 1) {
    return std::make_shared<SteinerTreeNode>(access_points.front().point,
                                             access_points.front().layers);
  }

  const int spine_x = spine_origin.x();
  const int spine_y = spine_origin.y();
  std::shared_ptr<SteinerTreeNode> root
      = std::make_shared<SteinerTreeNode>(spine_origin);

  std::vector<int> x_coords;
  std::vector<int> y_coords;
  x_coords.reserve(access_points.size() + 1);
  y_coords.reserve(access_points.size() + 1);
  for (const AccessPoint& access_point : access_points) {
    x_coords.push_back(access_point.point.x());
    y_coords.push_back(access_point.point.y());
  }
  x_coords.push_back(spine_x);
  y_coords.push_back(spine_y);
  std::sort(x_coords.begin(), x_coords.end());
  std::sort(y_coords.begin(), y_coords.end());
  x_coords.erase(std::unique(x_coords.begin(), x_coords.end()), x_coords.end());
  y_coords.erase(std::unique(y_coords.begin(), y_coords.end()), y_coords.end());

  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> horizontal_nodes;
  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> vertical_nodes;
  horizontal_nodes.reserve(x_coords.size());
  vertical_nodes.reserve(y_coords.size());
  for (const int x : x_coords) {
    std::shared_ptr<SteinerTreeNode> node
        = (x == spine_x ? root
                        : std::make_shared<SteinerTreeNode>(PointT(x, spine_y)));
    horizontal_nodes.emplace(x, node);
  }
  for (const int y : y_coords) {
    std::shared_ptr<SteinerTreeNode> node
        = (y == spine_y ? root
                        : std::make_shared<SteinerTreeNode>(PointT(spine_x, y)));
    vertical_nodes.emplace(y, node);
  }

  auto connectLine = [](const std::vector<int>& coordinates,
                        const int root_coordinate,
                        std::unordered_map<int, std::shared_ptr<SteinerTreeNode>>&
                            nodes) {
    const auto root_it
        = std::find(coordinates.begin(), coordinates.end(), root_coordinate);
    if (root_it == coordinates.end()) {
      return;
    }
    const int root_index = std::distance(coordinates.begin(), root_it);
    for (int index = root_index - 1; index >= 0; index--) {
      nodes[coordinates[index + 1]]->addChild(nodes[coordinates[index]]);
    }
    for (int index = root_index + 1; index < static_cast<int>(coordinates.size());
         index++) {
      nodes[coordinates[index - 1]]->addChild(nodes[coordinates[index]]);
    }
  };

  connectLine(x_coords, spine_x, horizontal_nodes);
  connectLine(y_coords, spine_y, vertical_nodes);

  std::unordered_map<uint64_t, std::shared_ptr<SteinerTreeNode>> node_by_point;
  node_by_point.reserve(horizontal_nodes.size() + vertical_nodes.size()
                        + access_points.size());
  node_by_point.emplace(pointKey(spine_origin), root);
  for (const auto& [x, node] : horizontal_nodes) {
    node_by_point.emplace(pointKey(PointT(x, spine_y)), node);
  }
  for (const auto& [y, node] : vertical_nodes) {
    node_by_point.emplace(pointKey(PointT(spine_x, y)), node);
  }

  for (const AccessPoint& access_point : access_points) {
    const uint64_t point_key = pointKey(access_point.point);
    auto existing = node_by_point.find(point_key);
    if (existing != node_by_point.end()) {
      mergeFixedLayers(access_point.layers, existing->second);
      continue;
    }

    const int dx = std::abs(access_point.point.x() - spine_x);
    const int dy = std::abs(access_point.point.y() - spine_y);
    const bool attach_vertical = dx < dy || (dx == dy && prefer_vertical);
    const PointT anchor_point
        = attach_vertical ? PointT(spine_x, access_point.point.y())
                          : PointT(access_point.point.x(), spine_y);
    auto anchor_it = node_by_point.find(pointKey(anchor_point));
    if (anchor_it == node_by_point.end()) {
      continue;
    }

    std::shared_ptr<SteinerTreeNode> pin_node
        = std::make_shared<SteinerTreeNode>(access_point.point,
                                            access_point.layers);
    anchor_it->second->addChild(pin_node);
    node_by_point.emplace(point_key, pin_node);
  }

  return root;
}

std::shared_ptr<SteinerTreeNode> buildQuadrantHubSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const PointT& root_point)
{
  if (access_points.empty()) {
    return nullptr;
  }
  if (access_points.size() == 1) {
    return std::make_shared<SteinerTreeNode>(access_points.front().point,
                                             access_points.front().layers);
  }

  std::shared_ptr<SteinerTreeNode> root
      = std::make_shared<SteinerTreeNode>(root_point);
  std::unordered_map<uint64_t, std::shared_ptr<SteinerTreeNode>> node_by_point;
  node_by_point.reserve(access_points.size() + 8);
  node_by_point.emplace(pointKey(root_point), root);

  std::array<std::vector<const AccessPoint*>, 4> clusters;
  for (const AccessPoint& access_point : access_points) {
    const int x_side = access_point.point.x() >= root_point.x() ? 1 : 0;
    const int y_side = access_point.point.y() >= root_point.y() ? 1 : 0;
    const int cluster_index = (x_side << 1) | y_side;
    clusters[cluster_index].push_back(&access_point);
  }

  for (const auto& cluster : clusters) {
    if (cluster.empty()) {
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(cluster.size());
    ys.reserve(cluster.size());
    for (const AccessPoint* access_point : cluster) {
      xs.push_back(access_point->point.x());
      ys.push_back(access_point->point.y());
    }
    const PointT hub_point(getMedian(xs), getMedian(ys));

    const uint64_t hub_key = pointKey(hub_point);
    auto [hub_it, inserted]
        = node_by_point.emplace(hub_key,
                                std::make_shared<SteinerTreeNode>(hub_point));
    std::shared_ptr<SteinerTreeNode> hub = hub_it->second;
    if (inserted) {
      attachChildUnique(root, hub);
    } else if (hub != root) {
      attachChildUnique(root, hub);
    }

    for (const AccessPoint* access_point : cluster) {
      const uint64_t pin_key = pointKey(access_point->point);
      auto existing = node_by_point.find(pin_key);
      if (existing != node_by_point.end()) {
        mergeFixedLayers(access_point->layers, existing->second);
        continue;
      }

      std::shared_ptr<SteinerTreeNode> pin_node
          = std::make_shared<SteinerTreeNode>(access_point->point,
                                              access_point->layers);
      attachChildUnique(hub, pin_node);
      node_by_point.emplace(pin_key, pin_node);
    }
  }

  return root;
}

std::shared_ptr<SteinerTreeNode> buildLadderBackboneSteinerTree(
    const std::vector<AccessPoint>& access_points,
    const int spine_coord,
    const int lane_a,
    const int lane_b,
    const bool vertical_spine)
{
  if (access_points.empty()) {
    return nullptr;
  }
  if (access_points.size() == 1) {
    return std::make_shared<SteinerTreeNode>(access_points.front().point,
                                             access_points.front().layers);
  }

  if (vertical_spine) {
    std::shared_ptr<SteinerTreeNode> root_a
        = std::make_shared<SteinerTreeNode>(PointT(spine_coord, lane_a));
    std::shared_ptr<SteinerTreeNode> root_b = root_a;
    if (lane_b != lane_a) {
      root_b = std::make_shared<SteinerTreeNode>(PointT(spine_coord, lane_b));
      attachChildUnique(root_a, root_b);
    }

    std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> lane_a_nodes;
    std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> lane_b_nodes;
    lane_a_nodes.reserve(access_points.size() + 1);
    lane_b_nodes.reserve(access_points.size() + 1);
    lane_a_nodes.emplace(spine_coord, root_a);
    lane_b_nodes.emplace(spine_coord, root_b);

    std::vector<bool> use_lane_a;
    use_lane_a.reserve(access_points.size());
    for (const AccessPoint& access_point : access_points) {
      const int dist_a = std::abs(access_point.point.y() - lane_a);
      const int dist_b = std::abs(access_point.point.y() - lane_b);
      const bool pick_a = dist_a <= dist_b;
      use_lane_a.push_back(pick_a);
      auto& lane_nodes = pick_a ? lane_a_nodes : lane_b_nodes;
      if (lane_nodes.find(access_point.point.x()) == lane_nodes.end()) {
        lane_nodes.emplace(
            access_point.point.x(),
            std::make_shared<SteinerTreeNode>(
                PointT(access_point.point.x(), pick_a ? lane_a : lane_b)));
      }
    }

    auto connectLane = [&](std::unordered_map<int,
                                               std::shared_ptr<SteinerTreeNode>>&
                               lane_nodes) {
      std::vector<int> coords;
      coords.reserve(lane_nodes.size());
      for (const auto& [coord, _] : lane_nodes) {
        coords.push_back(coord);
      }
      std::sort(coords.begin(), coords.end());
      const auto root_it = std::find(coords.begin(), coords.end(), spine_coord);
      if (root_it == coords.end()) {
        return;
      }
      const int root_idx = std::distance(coords.begin(), root_it);
      for (int i = root_idx - 1; i >= 0; i--) {
        attachChildUnique(lane_nodes[coords[i + 1]], lane_nodes[coords[i]]);
      }
      for (int i = root_idx + 1; i < static_cast<int>(coords.size()); i++) {
        attachChildUnique(lane_nodes[coords[i - 1]], lane_nodes[coords[i]]);
      }
    };
    connectLane(lane_a_nodes);
    connectLane(lane_b_nodes);

    for (int pin_idx = 0; pin_idx < static_cast<int>(access_points.size());
         pin_idx++) {
      const AccessPoint& access_point = access_points[pin_idx];
      const bool pick_a = use_lane_a[pin_idx];
      auto& lane_nodes = pick_a ? lane_a_nodes : lane_b_nodes;
      auto anchor_it = lane_nodes.find(access_point.point.x());
      if (anchor_it == lane_nodes.end()) {
        continue;
      }

      std::shared_ptr<SteinerTreeNode> anchor = anchor_it->second;
      if (anchor->x() == access_point.point.x()
          && anchor->y() == access_point.point.y()) {
        mergeFixedLayers(access_point.layers, anchor);
        continue;
      }

      anchor->addChild(std::make_shared<SteinerTreeNode>(access_point.point,
                                                          access_point.layers));
    }

    return root_a;
  }

  std::shared_ptr<SteinerTreeNode> root_a
      = std::make_shared<SteinerTreeNode>(PointT(lane_a, spine_coord));
  std::shared_ptr<SteinerTreeNode> root_b = root_a;
  if (lane_b != lane_a) {
    root_b = std::make_shared<SteinerTreeNode>(PointT(lane_b, spine_coord));
    attachChildUnique(root_a, root_b);
  }

  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> lane_a_nodes;
  std::unordered_map<int, std::shared_ptr<SteinerTreeNode>> lane_b_nodes;
  lane_a_nodes.reserve(access_points.size() + 1);
  lane_b_nodes.reserve(access_points.size() + 1);
  lane_a_nodes.emplace(spine_coord, root_a);
  lane_b_nodes.emplace(spine_coord, root_b);

  std::vector<bool> use_lane_a;
  use_lane_a.reserve(access_points.size());
  for (const AccessPoint& access_point : access_points) {
    const int dist_a = std::abs(access_point.point.x() - lane_a);
    const int dist_b = std::abs(access_point.point.x() - lane_b);
    const bool pick_a = dist_a <= dist_b;
    use_lane_a.push_back(pick_a);
    auto& lane_nodes = pick_a ? lane_a_nodes : lane_b_nodes;
    if (lane_nodes.find(access_point.point.y()) == lane_nodes.end()) {
      lane_nodes.emplace(
          access_point.point.y(),
          std::make_shared<SteinerTreeNode>(
              PointT(pick_a ? lane_a : lane_b, access_point.point.y())));
    }
  }

  auto connectLane = [&](std::unordered_map<int, std::shared_ptr<SteinerTreeNode>>&
                             lane_nodes) {
    std::vector<int> coords;
    coords.reserve(lane_nodes.size());
    for (const auto& [coord, _] : lane_nodes) {
      coords.push_back(coord);
    }
    std::sort(coords.begin(), coords.end());
    const auto root_it = std::find(coords.begin(), coords.end(), spine_coord);
    if (root_it == coords.end()) {
      return;
    }
    const int root_idx = std::distance(coords.begin(), root_it);
    for (int i = root_idx - 1; i >= 0; i--) {
      attachChildUnique(lane_nodes[coords[i + 1]], lane_nodes[coords[i]]);
    }
    for (int i = root_idx + 1; i < static_cast<int>(coords.size()); i++) {
      attachChildUnique(lane_nodes[coords[i - 1]], lane_nodes[coords[i]]);
    }
  };
  connectLane(lane_a_nodes);
  connectLane(lane_b_nodes);

  for (int pin_idx = 0; pin_idx < static_cast<int>(access_points.size());
       pin_idx++) {
    const AccessPoint& access_point = access_points[pin_idx];
    const bool pick_a = use_lane_a[pin_idx];
    auto& lane_nodes = pick_a ? lane_a_nodes : lane_b_nodes;
    auto anchor_it = lane_nodes.find(access_point.point.y());
    if (anchor_it == lane_nodes.end()) {
      continue;
    }

    std::shared_ptr<SteinerTreeNode> anchor = anchor_it->second;
    if (anchor->x() == access_point.point.x()
        && anchor->y() == access_point.point.y()) {
      mergeFixedLayers(access_point.layers, anchor);
      continue;
    }

    anchor->addChild(std::make_shared<SteinerTreeNode>(access_point.point,
                                                        access_point.layers));
  }

  return root_a;
}

}  // namespace

CUGR::CUGR(odb::dbDatabase* db,
           utl::Logger* log,
           stt::SteinerTreeBuilder* stt_builder)
    : db_(db), logger_(log), stt_builder_(stt_builder)
{
  // Radical wirelength-first policy:
  // keep soft-cap shaping but bias it much closer to hard capacity so the
  // search can keep direct trunks instead of spilling into long detours.
  constants_.weight_wire_length = 2.4;
  constants_.weight_via_number = 1.4;
  constants_.weight_short_area = 150.0;
  constants_.cost_logistic_slope = 0.22;
  constants_.maze_logistic_slope = 0.20;
  constants_.soft_cap_min_ratio = 0.82;
  constants_.soft_cap_max_ratio = 0.99;
  constants_.soft_cap_mid_util = 0.91;
  constants_.soft_cap_slope = 4.8;
  constants_.soft_cap_neighbor_weight = 0.22;
  constants_.maze_bbox_penalty = 1.85;
  constants_.maze_bbox_padding = 10;
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
    MazeRoute mazeRoute(net, grid_graph_.get(), constants_, logger_);
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
  constexpr int kMaxRefineNets = 1200;
  constexpr int kViaGrowthLimit = 360;
  constexpr int kOverflowSlack = 3;
  constexpr uint64_t kMinWireImprovement = 3;
  constexpr uint64_t kStrongWireImprovement = 12;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();
  // Radical pass: test multiple topology families and keep only local wins.
  grid_graph_->setCongestionPenaltyScales(0.006, original_maze_scale);

  struct Candidate
  {
    int net_index;
    int hpwl;
    int overflow;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (net->getNumPins() >= 3 && routing_tree) {
      const int hpwl = std::max(1, net->getBoundingBox().hp());
      const int overflow = grid_graph_->checkOverflow(routing_tree);
      const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
      const double stretch = static_cast<double>(stats.wire_length) / hpwl;
      candidates.push_back(
          {net->getIndex(), hpwl, overflow, net->getNumPins(), stretch});
    }
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxRefineNets) {
    candidates.resize(kMaxRefineNets);
  }

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  int accepted = 0;
  int attempted = 0;
  logger_->report("stage 4: multi-topology wirelength collapse on {} nets",
                  candidates.size());

  auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
    if (lhs.stats.wire_length != rhs.stats.wire_length) {
      return lhs.stats.wire_length < rhs.stats.wire_length;
    }
    if (lhs.overflow != rhs.overflow) {
      return lhs.overflow < rhs.overflow;
    }
    return lhs.stats.via_count < rhs.stats.via_count;
  };

  for (const Candidate& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> oldTree = net->getRoutingTree();
    if (!oldTree) {
      continue;
    }

    attempted++;
    const TreeStats oldStats = getTreeStats(oldTree, grid_graph_.get());

    const int oldOverflow = grid_graph_->checkOverflow(oldTree);
    grid_graph_->commitTree(oldTree, /*rip_up*/ true);

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improveWire
          = result.stats.wire_length + kMinWireImprovement
            <= oldStats.wire_length;
      const bool strongWireImprovement
          = result.stats.wire_length + kStrongWireImprovement
            <= oldStats.wire_length;
      const bool relieveOverflow
          = result.overflow + 2 < oldOverflow
            && result.stats.wire_length <= oldStats.wire_length + 6;
      const bool keepOverflow = result.overflow <= oldOverflow + kOverflowSlack;
      const bool boundedVia
          = result.stats.via_count <= oldStats.via_count + kViaGrowthLimit;
      return boundedVia
             && ((improveWire && keepOverflow) || strongWireImprovement
                 || relieveOverflow);
    };

    RerouteResult bestResult;
    auto considerCurrent = [&]() {
      const RerouteResult result = evaluate(net->getRoutingTree());
      if (isAcceptable(result)
          && (!bestResult.valid || isBetter(result, bestResult))) {
        bestResult = result;
      }
    };

    PatternRoute flutePattern(
        net, grid_graph_.get(), stt_builder_, constants_, logger_);
    flutePattern.constructSteinerTree();
    flutePattern.constructRoutingDAG();
    flutePattern.run();
    considerCurrent();

    const AccessPointSet selectedAccessPoints = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> accessPoints
        = materializeAccessPoints(selectedAccessPoints);
    if (accessPoints.size() >= 2) {
      auto trySteinerTemplate
          = [&](const std::shared_ptr<SteinerTreeNode>& steinerTree) {
              if (!steinerTree) {
                return;
              }
              PatternRoute pattern(
                  net, grid_graph_.get(), stt_builder_, constants_, logger_);
              pattern.setSteinerTree(steinerTree);
              pattern.constructRoutingDAG();
              pattern.run();
              considerCurrent();
            };

      trySteinerTemplate(buildMedianTrunkSteinerTree(accessPoints, true));
      trySteinerTemplate(buildMedianTrunkSteinerTree(accessPoints, false));

      const PointT center(net->getBoundingBox().cx(), net->getBoundingBox().cy());
      trySteinerTemplate(
          buildCrossbarSteinerTree(accessPoints, center, true));
      trySteinerTemplate(
          buildCrossbarSteinerTree(accessPoints, center, false));
    }

    if (bestResult.valid) {
      net->setRoutingTree(bestResult.tree);
      grid_graph_->commitTree(bestResult.tree);
      accepted++;
      continue;
    }

    net->setRoutingTree(oldTree);
    grid_graph_->commitTree(oldTree);
  }

  logger_->report("multi-topology collapse accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::mazeWirelengthCollapse()
{
  constexpr int kMaxCriticalNets = 640;
  constexpr int kViaGrowthLimit = 220;
  constexpr uint64_t kMinWireImprovement = 1;
  constexpr uint64_t kStrongWireImprovement = 10;
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

    MazeRoute mazeRoute(net, grid_graph_.get(), constants_, logger_);
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
    const bool strongWireImprove
        = candidateStats.wire_length + kStrongWireImprovement
          < oldStats.wire_length;
    const bool viaGrowthBound
        = candidateStats.via_count <= oldStats.via_count + kViaGrowthLimit;
    const bool keepOverflow = candidateOverflow <= oldOverflow + 2;
    const bool relieveOverflowWithoutWireRegression
        = (candidateOverflow + 2 < oldOverflow
           && candidateStats.wire_length <= oldStats.wire_length);

    if (viaGrowthBound
        && (((improveWire && keepOverflow) || strongWireImprove)
            || relieveOverflowWithoutWireRegression)) {
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
  constexpr int kMaxSurgeryNets = 120;
  constexpr int kViaGrowthLimit = 120;
  constexpr uint64_t kMinWireImprovement = 12;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical policy: replace FLUTE-first topology with median-spine templates
  // and keep only wirelength wins that preserve or improve overflow.
  grid_graph_->setCongestionPenaltyScales(0.008, 0.016);

  struct SurgeryCandidate
  {
    int net_index;
    int overflow;
    int hpwl;
    double stretch;
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
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    candidates.push_back({net->getIndex(), overflow, hpwl, stretch});
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
        if (lhs.overflow != rhs.overflow) {
          return lhs.overflow > rhs.overflow;
        }
        if (lhs.stretch != rhs.stretch) {
          return lhs.stretch > rhs.stretch;
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
  SparseGrid sparseGrid(5, 5, 0, 0);

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 6: median-spine topology surgery on {} critical nets",
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

    const AccessPointSet selectedAccessPoints = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> accessPoints
        = materializeAccessPoints(selectedAccessPoints);

    RerouteResult spineHorizontalResult;
    RerouteResult spineVerticalResult;
    if (accessPoints.size() >= 2) {
      auto runSpineCandidate = [&](const bool horizontal, RerouteResult& result) {
        std::shared_ptr<SteinerTreeNode> spineTree
            = buildMedianTrunkSteinerTree(accessPoints, horizontal);
        if (!spineTree) {
          return;
        }
        PatternRoute spinePattern(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        spinePattern.setSteinerTree(spineTree);
        spinePattern.constructRoutingDAG();
        spinePattern.run();
        result = evaluateTree(net->getRoutingTree());
      };
      runSpineCandidate(/*horizontal=*/true, spineHorizontalResult);
      runSpineCandidate(/*horizontal=*/false, spineVerticalResult);
    }

    // Candidate C: maze-driven topology rebuild, then pattern legalization.
    RerouteResult mazeResult;
    MazeRoute mazeRoute(net, grid_graph_.get(), constants_, logger_);
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
      const bool overflowRescue
          = result.overflow + 3 < oldOverflow
            && result.stats.wire_length <= oldStats.wire_length;
      return boundedVia && ((improveWire && keepOverflow) || overflowRescue);
    };

    const bool spineHorizontalAcceptable = isAcceptable(spineHorizontalResult);
    const bool spineVerticalAcceptable = isAcceptable(spineVerticalResult);
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
    if (spineHorizontalAcceptable) {
      best = &spineHorizontalResult;
    }
    if (spineVerticalAcceptable
        && (!best || isBetter(spineVerticalResult, *best))) {
      best = &spineVerticalResult;
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

  logger_->report("median-spine topology surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::hubTopologySurgery()
{
  constexpr int kMaxSurgeryNets = 28;
  constexpr int kMaxHubCandidates = 10;
  constexpr int kViaGrowthLimit = 160;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 4;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: sweep multiple hub locations and rebuild large nets using
  // hub-centered shared spines; keep only wirelength wins under tight
  // overflow/via constraints.
  grid_graph_->setCongestionPenaltyScales(0.008, 0.016);

  struct Candidate
  {
    int net_index;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const auto routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 4) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 40) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    candidates.push_back(
        {net->getIndex(), overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 7: hub-sweep topology surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);
    const std::vector<PointT> hub_candidates
        = buildHubCandidates(access_points,
                             net->getBoundingBox(),
                             kMaxHubCandidates);

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool preserve_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 2 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length;
      return bounded_via
             && ((improve_wire && preserve_overflow) || overflow_rescue);
    };

    auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
      if (lhs.stats.wire_length != rhs.stats.wire_length) {
        return lhs.stats.wire_length < rhs.stats.wire_length;
      }
      if (lhs.overflow != rhs.overflow) {
        return lhs.overflow < rhs.overflow;
      }
      return lhs.stats.via_count < rhs.stats.via_count;
    };

    RerouteResult best_result;
    for (const PointT& hub : hub_candidates) {
      for (const bool horizontal : {true, false}) {
        std::shared_ptr<SteinerTreeNode> spine_tree
            = buildHubSpineSteinerTree(access_points, hub, horizontal);
        if (!spine_tree) {
          continue;
        }
        PatternRoute pattern_route(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        pattern_route.setSteinerTree(spine_tree);
        pattern_route.constructRoutingDAG();
        pattern_route.run();
        RerouteResult result = evaluate(net->getRoutingTree());
        if (isAcceptable(result)
            && (!best_result.valid || isBetter(result, best_result))) {
          best_result = result;
        }
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("hub-sweep topology surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::dualHubBackboneSurgery()
{
  constexpr int kMaxSurgeryNets = 36;
  constexpr int kMaxHubCandidates = 8;
  constexpr int kMaxHubPairs = 14;
  constexpr int kViaGrowthLimit = 240;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 4;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: build two-hub shared backbones so high-degree nets can
  // collapse multiple long branches into a small number of long trunks.
  grid_graph_->setCongestionPenaltyScales(0.007, 0.015);

  struct Candidate
  {
    int net_index;
    int score;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 6) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 80) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    const int score = overflow * 2200 + hpwl + net->getNumPins() * 10
                      + static_cast<int>(stretch * 100.0);
    candidates.push_back(
        {net->getIndex(), score, overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct HubPair
  {
    PointT hub_a;
    PointT hub_b;
    int span;
  };

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 8: dual-hub backbone surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);
    const std::vector<PointT> hub_candidates
        = buildHubCandidates(access_points,
                             net->getBoundingBox(),
                             kMaxHubCandidates);
    if (hub_candidates.size() < 2) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::vector<HubPair> hub_pairs;
    hub_pairs.reserve(hub_candidates.size() * hub_candidates.size());
    for (int lhs = 0; lhs < static_cast<int>(hub_candidates.size()); lhs++) {
      for (int rhs = lhs + 1; rhs < static_cast<int>(hub_candidates.size());
           rhs++) {
        const int span
            = manhattanDistance(hub_candidates[lhs], hub_candidates[rhs]);
        if (span == 0) {
          continue;
        }
        hub_pairs.push_back(
            {hub_candidates[lhs], hub_candidates[rhs], span});
      }
    }
    if (hub_pairs.empty()) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::sort(hub_pairs.begin(),
              hub_pairs.end(),
              [](const HubPair& lhs, const HubPair& rhs) {
                if (lhs.span != rhs.span) {
                  return lhs.span > rhs.span;
                }
                if (lhs.hub_a.x() != rhs.hub_a.x()) {
                  return lhs.hub_a.x() < rhs.hub_a.x();
                }
                if (lhs.hub_a.y() != rhs.hub_a.y()) {
                  return lhs.hub_a.y() < rhs.hub_a.y();
                }
                if (lhs.hub_b.x() != rhs.hub_b.x()) {
                  return lhs.hub_b.x() < rhs.hub_b.x();
                }
                return lhs.hub_b.y() < rhs.hub_b.y();
              });
    if (static_cast<int>(hub_pairs.size()) > kMaxHubPairs) {
      hub_pairs.resize(kMaxHubPairs);
    }

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool preserve_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 3 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length;
      return bounded_via
             && ((improve_wire && preserve_overflow) || overflow_rescue);
    };

    auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
      if (lhs.stats.wire_length != rhs.stats.wire_length) {
        return lhs.stats.wire_length < rhs.stats.wire_length;
      }
      if (lhs.overflow != rhs.overflow) {
        return lhs.overflow < rhs.overflow;
      }
      return lhs.stats.via_count < rhs.stats.via_count;
    };

    RerouteResult best_result;
    for (const HubPair& hub_pair : hub_pairs) {
      for (const bool horizontal_first : {true, false}) {
        std::shared_ptr<SteinerTreeNode> dual_hub_tree
            = buildDualHubBackboneSteinerTree(access_points,
                                              hub_pair.hub_a,
                                              hub_pair.hub_b,
                                              horizontal_first);
        if (!dual_hub_tree) {
          continue;
        }

        PatternRoute pattern_route(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        pattern_route.setSteinerTree(dual_hub_tree);
        pattern_route.constructRoutingDAG();
        pattern_route.run();

        RerouteResult result = evaluate(net->getRoutingTree());
        if (isAcceptable(result)
            && (!best_result.valid || isBetter(result, best_result))) {
          best_result = result;
        }
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("dual-hub backbone surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::crossbarBackboneSurgery()
{
  constexpr int kMaxSurgeryNets = 52;
  constexpr int kMaxAxisCandidates = 10;
  constexpr int kMaxCrossbarPairs = 18;
  constexpr int kViaGrowthLimit = 260;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 6;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: rebuild high-degree nets with a full X/Y shared backbone.
  // This topology class differs from point-hub trees by forcing trunk reuse
  // along two global axes, often reducing repeated branch overlap.
  grid_graph_->setCongestionPenaltyScales(0.008, 0.017);

  struct Candidate
  {
    int net_index;
    int score;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 5) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 60) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    const int score = overflow * 3400 + hpwl + net->getNumPins() * 10
                      + static_cast<int>(stretch * 120.0);
    candidates.push_back(
        {net->getIndex(), score, overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct CrossbarPair
  {
    int x;
    int y;
    int score;
  };

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
    if (lhs.stats.wire_length != rhs.stats.wire_length) {
      return lhs.stats.wire_length < rhs.stats.wire_length;
    }
    if (lhs.overflow != rhs.overflow) {
      return lhs.overflow < rhs.overflow;
    }
    return lhs.stats.via_count < rhs.stats.via_count;
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 9: crossbar backbone surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points
        = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);

    const std::vector<int> x_candidates
        = buildAxisCandidates(access_points,
                              net->getBoundingBox(),
                              /*use_x_axis=*/true,
                              kMaxAxisCandidates);
    const std::vector<int> y_candidates
        = buildAxisCandidates(access_points,
                              net->getBoundingBox(),
                              /*use_x_axis=*/false,
                              kMaxAxisCandidates);
    if (x_candidates.empty() || y_candidates.empty()) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::vector<CrossbarPair> pairs;
    pairs.reserve(x_candidates.size() * y_candidates.size());
    const BoxT& bbox = net->getBoundingBox();
    for (const int x : x_candidates) {
      for (const int y : y_candidates) {
        int estimated_stem = 0;
        for (const AccessPoint& access_point : access_points) {
          const int dx = std::abs(access_point.point.x() - x);
          const int dy = std::abs(access_point.point.y() - y);
          estimated_stem += std::min(dx, dy);
        }
        const int center_bias
            = std::abs(x - bbox.cx()) + std::abs(y - bbox.cy());
        pairs.push_back({x, y, estimated_stem * 2 + center_bias});
      }
    }
    if (pairs.empty()) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::sort(pairs.begin(),
              pairs.end(),
              [](const CrossbarPair& lhs, const CrossbarPair& rhs) {
                if (lhs.score != rhs.score) {
                  return lhs.score < rhs.score;
                }
                if (lhs.x != rhs.x) {
                  return lhs.x < rhs.x;
                }
                return lhs.y < rhs.y;
              });
    if (static_cast<int>(pairs.size()) > kMaxCrossbarPairs) {
      pairs.resize(kMaxCrossbarPairs);
    }

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool keep_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 3 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length;
      return bounded_via
             && ((improve_wire && keep_overflow) || overflow_rescue);
    };

    RerouteResult best_result;
    for (const CrossbarPair& pair : pairs) {
      for (const bool prefer_vertical : {true, false}) {
        std::shared_ptr<SteinerTreeNode> crossbar_tree
            = buildCrossbarSteinerTree(
                access_points, PointT(pair.x, pair.y), prefer_vertical);
        if (!crossbar_tree) {
          continue;
        }

        PatternRoute pattern_route(
            net, grid_graph_.get(), stt_builder_, constants_, logger_);
        pattern_route.setSteinerTree(crossbar_tree);
        pattern_route.constructRoutingDAG();
        pattern_route.run();

        RerouteResult result = evaluate(net->getRoutingTree());
        if (isAcceptable(result)
            && (!best_result.valid || isBetter(result, best_result))) {
          best_result = result;
        }
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("crossbar backbone surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::quadrantHubHierarchySurgery()
{
  constexpr int kMaxSurgeryNets = 64;
  constexpr int kMaxRootCandidates = 6;
  constexpr int kViaGrowthLimit = 280;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 4;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: rebuild high-pin nets using a two-level hierarchy
  // (quadrant hubs + root trunk) to force branch sharing and collapse stems.
  grid_graph_->setCongestionPenaltyScales(0.006, 0.014);

  struct Candidate
  {
    int net_index;
    int score;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 8) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 90) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    const int score = overflow * 3600 + hpwl + net->getNumPins() * 12
                      + static_cast<int>(stretch * 140.0);
    candidates.push_back(
        {net->getIndex(), score, overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
    if (lhs.stats.wire_length != rhs.stats.wire_length) {
      return lhs.stats.wire_length < rhs.stats.wire_length;
    }
    if (lhs.overflow != rhs.overflow) {
      return lhs.overflow < rhs.overflow;
    }
    return lhs.stats.via_count < rhs.stats.via_count;
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 10: quadrant-hub hierarchy surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points
        = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);
    if (access_points.size() < 4) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::vector<PointT> root_candidates
        = buildHubCandidates(access_points,
                             net->getBoundingBox(),
                             kMaxRootCandidates);
    if (root_candidates.empty()) {
      root_candidates.push_back(
          PointT(net->getBoundingBox().cx(), net->getBoundingBox().cy()));
    }

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool keep_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 3 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length;
      return bounded_via
             && ((improve_wire && keep_overflow) || overflow_rescue);
    };

    RerouteResult best_result;
    for (const PointT& root_point : root_candidates) {
      std::shared_ptr<SteinerTreeNode> hierarchy_tree
          = buildQuadrantHubSteinerTree(access_points, root_point);
      if (!hierarchy_tree) {
        continue;
      }

      PatternRoute pattern_route(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      pattern_route.setSteinerTree(hierarchy_tree);
      pattern_route.constructRoutingDAG();
      pattern_route.run();

      RerouteResult result = evaluate(net->getRoutingTree());
      if (isAcceptable(result)
          && (!best_result.valid || isBetter(result, best_result))) {
        best_result = result;
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("quadrant-hub hierarchy surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::ladderBackboneSurgery()
{
  constexpr int kMaxSurgeryNets = 48;
  constexpr int kMaxConfigs = 14;
  constexpr int kViaGrowthLimit = 300;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 2;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: rebuild selected high-pin nets with a two-lane ladder
  // backbone (parallel trunks + connecting spine). This is intentionally
  // different from hub/crossbar/MST topologies and targets stem collapse.
  grid_graph_->setCongestionPenaltyScales(0.006, 0.013);

  struct Candidate
  {
    int net_index;
    int score;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 7) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 70) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    const int score = overflow * 4200 + hpwl + net->getNumPins() * 18
                      + static_cast<int>(stretch * 200.0);
    candidates.push_back(
        {net->getIndex(), score, overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct LadderConfig
  {
    bool vertical_spine;
    int spine_coord;
    int lane_a;
    int lane_b;
    int estimate;
  };

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
    if (lhs.stats.wire_length != rhs.stats.wire_length) {
      return lhs.stats.wire_length < rhs.stats.wire_length;
    }
    if (lhs.overflow != rhs.overflow) {
      return lhs.overflow < rhs.overflow;
    }
    return lhs.stats.via_count < rhs.stats.via_count;
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 11: ladder-backbone surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points
        = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);
    if (access_points.size() < 4) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(access_points.size());
    ys.reserve(access_points.size());
    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    for (const AccessPoint& access_point : access_points) {
      xs.push_back(access_point.point.x());
      ys.push_back(access_point.point.y());
      min_x = std::min(min_x, access_point.point.x());
      max_x = std::max(max_x, access_point.point.x());
      min_y = std::min(min_y, access_point.point.y());
      max_y = std::max(max_y, access_point.point.y());
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());

    const int q1_x = xs[xs.size() / 4];
    const int median_x = xs[xs.size() / 2];
    const int q3_x = xs[(3 * xs.size()) / 4];
    const int q1_y = ys[ys.size() / 4];
    const int median_y = ys[ys.size() / 2];
    const int q3_y = ys[(3 * ys.size()) / 4];
    const BoxT& bbox = net->getBoundingBox();

    std::vector<int> x_spines{median_x, bbox.cx(), q1_x, q3_x};
    std::vector<int> y_spines{median_y, bbox.cy(), q1_y, q3_y};
    std::vector<std::pair<int, int>> horizontal_lanes{
        {q1_y, q3_y}, {bbox.ly(), bbox.hy()}, {q1_y, median_y}, {median_y, q3_y}};
    std::vector<std::pair<int, int>> vertical_lanes{
        {q1_x, q3_x}, {bbox.lx(), bbox.hx()}, {q1_x, median_x}, {median_x, q3_x}};

    auto normalizeLanes = [](std::vector<std::pair<int, int>>& lanes) {
      for (auto& lane_pair : lanes) {
        if (lane_pair.first > lane_pair.second) {
          std::swap(lane_pair.first, lane_pair.second);
        }
      }
      std::sort(lanes.begin(), lanes.end());
      lanes.erase(std::unique(lanes.begin(), lanes.end()), lanes.end());
    };
    normalizeLanes(horizontal_lanes);
    normalizeLanes(vertical_lanes);

    std::vector<LadderConfig> configs;
    configs.reserve(48);
    auto appendConfigs = [&](const bool vertical_spine) {
      const auto& spines = vertical_spine ? x_spines : y_spines;
      const auto& lane_pairs = vertical_spine ? horizontal_lanes : vertical_lanes;
      for (const int spine_coord : spines) {
        for (const auto& [lane_a, lane_b] : lane_pairs) {
          int stem_cost = 0;
          for (const AccessPoint& access_point : access_points) {
            if (vertical_spine) {
              stem_cost += std::min(std::abs(access_point.point.y() - lane_a),
                                    std::abs(access_point.point.y() - lane_b));
            } else {
              stem_cost += std::min(std::abs(access_point.point.x() - lane_a),
                                    std::abs(access_point.point.x() - lane_b));
            }
          }
          const int trunk_span = vertical_spine
                                     ? ((lane_a == lane_b ? 1 : 2) * (max_x - min_x)
                                        + std::abs(lane_b - lane_a))
                                     : ((lane_a == lane_b ? 1 : 2) * (max_y - min_y)
                                        + std::abs(lane_b - lane_a));
          const int spine_bias = vertical_spine ? std::abs(spine_coord - bbox.cx())
                                                : std::abs(spine_coord - bbox.cy());
          configs.push_back(
              {vertical_spine, spine_coord, lane_a, lane_b, stem_cost * 2 + trunk_span
                                                             + spine_bias});
        }
      }
    };
    appendConfigs(/*vertical_spine=*/true);
    appendConfigs(/*vertical_spine=*/false);
    if (configs.empty()) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::sort(configs.begin(),
              configs.end(),
              [](const LadderConfig& lhs, const LadderConfig& rhs) {
                if (lhs.estimate != rhs.estimate) {
                  return lhs.estimate < rhs.estimate;
                }
                if (lhs.vertical_spine != rhs.vertical_spine) {
                  return lhs.vertical_spine < rhs.vertical_spine;
                }
                if (lhs.spine_coord != rhs.spine_coord) {
                  return lhs.spine_coord < rhs.spine_coord;
                }
                if (lhs.lane_a != rhs.lane_a) {
                  return lhs.lane_a < rhs.lane_a;
                }
                return lhs.lane_b < rhs.lane_b;
              });

    std::vector<LadderConfig> filtered_configs;
    filtered_configs.reserve(kMaxConfigs);
    std::unordered_set<uint64_t> seen;
    seen.reserve(configs.size());
    for (const LadderConfig& config : configs) {
      const uint64_t key = (static_cast<uint64_t>(config.vertical_spine) << 48)
                           | (static_cast<uint64_t>(static_cast<uint16_t>(
                                  config.spine_coord))
                              << 32)
                           | (static_cast<uint64_t>(static_cast<uint16_t>(
                                  config.lane_a))
                              << 16)
                           | static_cast<uint16_t>(config.lane_b);
      if (!seen.insert(key).second) {
        continue;
      }
      filtered_configs.push_back(config);
      if (static_cast<int>(filtered_configs.size()) >= kMaxConfigs) {
        break;
      }
    }

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool keep_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 3 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length;
      return bounded_via
             && ((improve_wire && keep_overflow) || overflow_rescue);
    };

    RerouteResult best_result;
    for (const LadderConfig& config : filtered_configs) {
      std::shared_ptr<SteinerTreeNode> ladder_tree
          = buildLadderBackboneSteinerTree(access_points,
                                           config.spine_coord,
                                           config.lane_a,
                                           config.lane_b,
                                           config.vertical_spine);
      if (!ladder_tree) {
        continue;
      }

      PatternRoute pattern_route(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      pattern_route.setSteinerTree(ladder_tree);
      pattern_route.constructRoutingDAG();
      pattern_route.run();

      RerouteResult result = evaluate(net->getRoutingTree());
      if (isAcceptable(result)
          && (!best_result.valid || isBetter(result, best_result))) {
        best_result = result;
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("ladder-backbone surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::mstBackboneSurgery()
{
  constexpr int kMaxSurgeryNets = 40;
  constexpr int kViaGrowthLimit = 320;
  constexpr int kOverflowSlack = 1;
  constexpr uint64_t kMinWireImprovement = 4;
  const double original_wire_scale = grid_graph_->getWireCongestionScale();
  const double original_maze_scale = grid_graph_->getMazeCongestionScale();

  // Radical move: replace tree templates with Manhattan MST backbones, with an
  // optional synthetic hub to force long trunk sharing.
  grid_graph_->setCongestionPenaltyScales(0.007, 0.014);

  struct Candidate
  {
    int net_index;
    int score;
    int overflow;
    int hpwl;
    int pins;
    double stretch;
  };

  std::vector<Candidate> candidates;
  candidates.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    const std::shared_ptr<GRTreeNode> routing_tree = net->getRoutingTree();
    if (!routing_tree || net->getNumPins() < 6) {
      continue;
    }
    const int hpwl = std::max(1, net->getBoundingBox().hp());
    if (hpwl < 70) {
      continue;
    }
    const TreeStats stats = getTreeStats(routing_tree, grid_graph_.get());
    const int overflow = grid_graph_->checkOverflow(routing_tree);
    const double stretch = static_cast<double>(stats.wire_length) / hpwl;
    const int score = overflow * 3800 + hpwl + net->getNumPins() * 16
                      + static_cast<int>(stretch * 160.0);
    candidates.push_back(
        {net->getIndex(), score, overflow, hpwl, net->getNumPins(), stretch});
  }

  if (candidates.empty()) {
    grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                            original_maze_scale);
    return;
  }

  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate& lhs, const Candidate& rhs) {
              if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
              }
              if (lhs.overflow != rhs.overflow) {
                return lhs.overflow > rhs.overflow;
              }
              if (lhs.stretch != rhs.stretch) {
                return lhs.stretch > rhs.stretch;
              }
              if (lhs.pins != rhs.pins) {
                return lhs.pins > rhs.pins;
              }
              if (lhs.hpwl != rhs.hpwl) {
                return lhs.hpwl > rhs.hpwl;
              }
              return lhs.net_index < rhs.net_index;
            });

  if (static_cast<int>(candidates.size()) > kMaxSurgeryNets) {
    candidates.resize(kMaxSurgeryNets);
  }

  struct RerouteResult
  {
    std::shared_ptr<GRTreeNode> tree;
    TreeStats stats;
    int overflow{std::numeric_limits<int>::max()};
    bool valid{false};
  };

  auto isBetter = [](const RerouteResult& lhs, const RerouteResult& rhs) {
    if (lhs.stats.wire_length != rhs.stats.wire_length) {
      return lhs.stats.wire_length < rhs.stats.wire_length;
    }
    if (lhs.overflow != rhs.overflow) {
      return lhs.overflow < rhs.overflow;
    }
    return lhs.stats.via_count < rhs.stats.via_count;
  };

  int attempted = 0;
  int accepted = 0;
  logger_->report("stage 12: Manhattan-MST backbone surgery on {} nets",
                  candidates.size());

  for (const auto& candidate : candidates) {
    GRNet* net = gr_nets_[candidate.net_index].get();
    const std::shared_ptr<GRTreeNode> old_tree = net->getRoutingTree();
    if (!old_tree) {
      continue;
    }

    attempted++;
    const TreeStats old_stats = getTreeStats(old_tree, grid_graph_.get());
    const int old_overflow = grid_graph_->checkOverflow(old_tree);
    grid_graph_->commitTree(old_tree, /*rip_up*/ true);

    const AccessPointSet selected_access_points
        = grid_graph_->selectAccessPoints(net);
    const std::vector<AccessPoint> access_points
        = materializeAccessPoints(selected_access_points);
    if (access_points.size() < 3) {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
      continue;
    }

    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(access_points.size());
    ys.reserve(access_points.size());
    for (const AccessPoint& access_point : access_points) {
      xs.push_back(access_point.point.x());
      ys.push_back(access_point.point.y());
    }

    const PointT center(net->getBoundingBox().cx(), net->getBoundingBox().cy());
    const PointT median_hub(getMedian(xs), getMedian(ys));

    auto evaluate = [&](const std::shared_ptr<GRTreeNode>& tree) {
      RerouteResult result;
      if (!tree) {
        return result;
      }
      result.tree = tree;
      result.stats = getTreeStats(tree, grid_graph_.get());
      grid_graph_->commitTree(tree);
      result.overflow = grid_graph_->checkOverflow(tree);
      result.valid = true;
      grid_graph_->commitTree(tree, /*rip_up*/ true);
      return result;
    };

    auto isAcceptable = [&](const RerouteResult& result) {
      if (!result.valid) {
        return false;
      }
      const bool improve_wire
          = result.stats.wire_length + kMinWireImprovement <= old_stats.wire_length;
      const bool keep_overflow = result.overflow <= old_overflow + kOverflowSlack;
      const bool bounded_via
          = result.stats.via_count <= old_stats.via_count + kViaGrowthLimit;
      const bool overflow_rescue
          = result.overflow + 3 < old_overflow
            && result.stats.wire_length <= old_stats.wire_length + 2;
      return bounded_via
             && ((improve_wire && keep_overflow) || overflow_rescue);
    };

    RerouteResult best_result;
    auto tryCandidate = [&](const bool include_hub,
                            const PointT& hub,
                            const bool prefer_horizontal_first) {
      std::shared_ptr<SteinerTreeNode> mst_tree
          = buildRectilinearMstSteinerTree(access_points,
                                           center,
                                           include_hub,
                                           hub,
                                           prefer_horizontal_first);
      if (!mst_tree) {
        return;
      }
      PatternRoute pattern_route(
          net, grid_graph_.get(), stt_builder_, constants_, logger_);
      pattern_route.setSteinerTree(mst_tree);
      pattern_route.constructRoutingDAG();
      pattern_route.run();

      RerouteResult result = evaluate(net->getRoutingTree());
      if (isAcceptable(result)
          && (!best_result.valid || isBetter(result, best_result))) {
        best_result = result;
      }
    };

    for (const bool prefer_horizontal_first : {true, false}) {
      tryCandidate(/*include_hub=*/false, center, prefer_horizontal_first);
      tryCandidate(/*include_hub=*/true, median_hub, prefer_horizontal_first);
      if (median_hub != center) {
        tryCandidate(/*include_hub=*/true, center, prefer_horizontal_first);
      }
    }

    if (best_result.valid) {
      net->setRoutingTree(best_result.tree);
      grid_graph_->commitTree(best_result.tree);
      accepted++;
    } else {
      net->setRoutingTree(old_tree);
      grid_graph_->commitTree(old_tree);
    }
  }

  logger_->report("Manhattan-MST backbone surgery accepted {} / {} nets",
                  accepted,
                  attempted);
  grid_graph_->setCongestionPenaltyScales(original_wire_scale,
                                          original_maze_scale);
}

void CUGR::route()
{
  constexpr int kMaxDetourNets = 2000;
  constexpr int kMaxMazeNets = 680;
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

  // Backbone cascade strategy:
  //   1) median-spine collapse
  //   2) hub-sweep trunk sharing
  //   3) maze collapse + MST rebuild
  //   4) two wirelength-only cleanup pulses
  hybridTopologySurgery();
  hubTopologySurgery();
  mazeWirelengthCollapse();
  mstBackboneSurgery();
  grid_graph_->setCongestionPenaltyScales(0.05, 0.10);
  wirelengthRefine();
  grid_graph_->setCongestionPenaltyScales(0.04, 0.08);
  mazeWirelengthCollapse();
  grid_graph_->setCongestionPenaltyScales(0.04, 0.08);
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
