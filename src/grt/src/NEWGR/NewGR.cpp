#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Net.h"
#include "Pin.h"
#include "src/NewgrEngine.h"
#include "utl/Logger.h"

namespace grt {

namespace {

constexpr int kViaPenalty = 1000;

struct RouteStats
{
  int64_t wirelength{0};
  int via_count{0};
};

GSegment normalizeSegment(const GSegment& segment)
{
  GSegment normalized = segment;
  if (normalized.isVia() && normalized.init_layer > normalized.final_layer) {
    std::swap(normalized.init_layer, normalized.final_layer);
  }
  return normalized;
}

struct NodeKey
{
  int x;
  int y;
  int layer;

  bool operator==(const NodeKey& other) const
  {
    return x == other.x && y == other.y && layer == other.layer;
  }
};

struct NodeKeyHash
{
  size_t operator()(const NodeKey& key) const
  {
    size_t seed = 0;
    seed ^= std::hash<int>{}(key.x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>{}(key.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>{}(key.layer) + 0x9e3779b9 + (seed << 6)
            + (seed >> 2);
    return seed;
  }
};

struct EdgeKey
{
  int u;
  int v;

  bool operator==(const EdgeKey& other) const
  {
    return u == other.u && v == other.v;
  }
};

struct EdgeKeyHash
{
  size_t operator()(const EdgeKey& key) const
  {
    size_t seed = 0;
    seed ^= std::hash<int>{}(key.u) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>{}(key.v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct EdgeInfo
{
  int u;
  int v;
  int weight;
  GSegment segment;
};

class DisjointSet
{
 public:
  explicit DisjointSet(int size) : parent_(size), rank_(size, 0)
  {
    for (int i = 0; i < size; ++i) {
      parent_[i] = i;
    }
  }

  int find(int node)
  {
    if (parent_[node] != node) {
      parent_[node] = find(parent_[node]);
    }
    return parent_[node];
  }

  bool unite(int lhs, int rhs)
  {
    lhs = find(lhs);
    rhs = find(rhs);
    if (lhs == rhs) {
      return false;
    }
    if (rank_[lhs] < rank_[rhs]) {
      std::swap(lhs, rhs);
    }
    parent_[rhs] = lhs;
    if (rank_[lhs] == rank_[rhs]) {
      rank_[lhs]++;
    }
    return true;
  }

 private:
  std::vector<int> parent_;
  std::vector<int> rank_;
};

uint64_t xyKey(int x, int y)
{
  return (static_cast<uint64_t>(static_cast<uint32_t>(x)) << 32)
         | static_cast<uint32_t>(y);
}

int segmentWeight(const GSegment& segment)
{
  const int via_cost = std::abs(segment.final_layer - segment.init_layer)
                       * kViaPenalty;
  return segment.length() + via_cost;
}

int segmentViaSpan(const GSegment& segment)
{
  return std::abs(segment.final_layer - segment.init_layer);
}

bool isBetterSegmentForEqualWeight(const GSegment& candidate,
                                   const GSegment& incumbent)
{
  const int candidate_len = candidate.length();
  const int incumbent_len = incumbent.length();
  if (candidate_len != incumbent_len) {
    return candidate_len < incumbent_len;
  }

  const int candidate_vias = segmentViaSpan(candidate);
  const int incumbent_vias = segmentViaSpan(incumbent);
  if (candidate_vias != incumbent_vias) {
    return candidate_vias < incumbent_vias;
  }

  if (candidate.init_layer != incumbent.init_layer) {
    return candidate.init_layer < incumbent.init_layer;
  }
  if (candidate.final_layer != incumbent.final_layer) {
    return candidate.final_layer < incumbent.final_layer;
  }
  if (candidate.init_x != incumbent.init_x) {
    return candidate.init_x < incumbent.init_x;
  }
  if (candidate.init_y != incumbent.init_y) {
    return candidate.init_y < incumbent.init_y;
  }
  if (candidate.final_x != incumbent.final_x) {
    return candidate.final_x < incumbent.final_x;
  }
  return candidate.final_y < incumbent.final_y;
}

bool isBetterEdgeForMst(const EdgeInfo& lhs, const EdgeInfo& rhs)
{
  if (lhs.weight != rhs.weight) {
    return lhs.weight < rhs.weight;
  }

  const GSegment& left_seg = lhs.segment;
  const GSegment& right_seg = rhs.segment;
  if (isBetterSegmentForEqualWeight(left_seg, right_seg)) {
    return true;
  }
  if (isBetterSegmentForEqualWeight(right_seg, left_seg)) {
    return false;
  }

  if (lhs.u != rhs.u) {
    return lhs.u < rhs.u;
  }
  return lhs.v < rhs.v;
}

RouteStats computeRouteStats(const GRoute& route)
{
  RouteStats stats;
  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      stats.via_count += std::abs(segment.final_layer - segment.init_layer);
    } else {
      stats.wirelength += segment.length();
    }
  }
  return stats;
}

void dedupeAndDropStubs(GRoute& route)
{
  GRoute filtered;
  filtered.reserve(route.size());
  std::unordered_set<GSegment, GSegmentHash> seen;
  seen.reserve(route.size());
  for (const GSegment& segment : route) {
    const bool is_zero_length_stub
        = segment.init_x == segment.final_x && segment.init_y == segment.final_y
          && segment.init_layer == segment.final_layer;
    if (is_zero_length_stub) {
      continue;
    }
    const GSegment normalized_segment = normalizeSegment(segment);
    if (seen.insert(normalized_segment).second) {
      filtered.push_back(normalized_segment);
    }
  }
  route.swap(filtered);
}

void optimizeRouteTopology(const std::vector<Pin>& pins, GRoute& route)
{
  dedupeAndDropStubs(route);
  if (pins.size() < 3 || route.size() < 3) {
    return;
  }
  const RouteStats original_stats = computeRouteStats(route);

  std::unordered_map<NodeKey, int, NodeKeyHash> node_to_id;
  std::vector<NodeKey> nodes;
  nodes.reserve(route.size() * 2);
  auto getNodeId = [&](const NodeKey& node) -> int {
    const auto [it, inserted] = node_to_id.emplace(node, nodes.size());
    if (inserted) {
      nodes.push_back(node);
    }
    return it->second;
  };

  std::unordered_map<EdgeKey, EdgeInfo, EdgeKeyHash> edge_table;
  edge_table.reserve(route.size());

  for (const GSegment& raw_segment : route) {
    const GSegment segment = normalizeSegment(raw_segment);
    NodeKey src{segment.init_x, segment.init_y, segment.init_layer};
    NodeKey dst{segment.final_x, segment.final_y, segment.final_layer};
    int src_id = getNodeId(src);
    int dst_id = getNodeId(dst);
    if (src_id == dst_id) {
      continue;
    }
    if (src_id > dst_id) {
      std::swap(src_id, dst_id);
    }
    EdgeKey key{src_id, dst_id};
    const int weight = segmentWeight(segment);
    auto found = edge_table.find(key);
    if (found == edge_table.end()
        || weight < found->second.weight
        || (weight == found->second.weight
            && isBetterSegmentForEqualWeight(segment, found->second.segment))) {
      edge_table[key] = EdgeInfo{src_id, dst_id, weight, segment};
    }
  }

  if (edge_table.empty()) {
    return;
  }

  std::vector<EdgeInfo> edges;
  edges.reserve(edge_table.size());
  for (const auto& [ignored_key, edge] : edge_table) {
    static_cast<void>(ignored_key);
    edges.push_back(edge);
  }

  const int node_count = nodes.size();
  DisjointSet components(node_count);
  for (const EdgeInfo& edge : edges) {
    components.unite(edge.u, edge.v);
  }

  std::unordered_map<uint64_t, std::vector<int>> xy_to_nodes;
  xy_to_nodes.reserve(nodes.size());
  for (int node_id = 0; node_id < node_count; ++node_id) {
    xy_to_nodes[xyKey(nodes[node_id].x, nodes[node_id].y)].push_back(node_id);
  }

  std::vector<bool> required(node_count, false);
  for (const Pin& pin : pins) {
    const odb::Point& pin_pos = pin.getOnGridPosition();
    const uint64_t key = xyKey(pin_pos.x(), pin_pos.y());
    auto matching_nodes = xy_to_nodes.find(key);
    if (matching_nodes != xy_to_nodes.end() && !matching_nodes->second.empty()) {
      // Anchor each pin to one best-fit node to avoid preserving
      // unnecessary vertical stacks at the same (x, y).
      int best_node = -1;
      int best_score = std::numeric_limits<int>::max();
      const int pin_layer = pin.getConnectionLayer();
      for (int node_id : matching_nodes->second) {
        const NodeKey& node = nodes[node_id];
        const int layer_distance = std::min(std::abs(node.layer - pin_layer),
                                            std::abs(node.layer - (pin_layer + 1)));
        if (layer_distance < best_score) {
          best_score = layer_distance;
          best_node = node_id;
        }
      }
      if (best_node >= 0) {
        required[best_node] = true;
      }
      continue;
    }

    int nearest_node = -1;
    int best_distance = std::numeric_limits<int>::max();
    const int pin_layer = pin.getConnectionLayer();
    for (int node_id = 0; node_id < node_count; ++node_id) {
      const NodeKey& node = nodes[node_id];
      const int distance = std::abs(node.x - pin_pos.x())
                           + std::abs(node.y - pin_pos.y())
                           + 100 * std::min(std::abs(node.layer - pin_layer),
                                            std::abs(node.layer - (pin_layer + 1)));
      if (distance < best_distance) {
        best_distance = distance;
        nearest_node = node_id;
      }
    }
    if (nearest_node >= 0) {
      required[nearest_node] = true;
    }
  }

  std::unordered_map<int, std::vector<int>> component_edges;
  component_edges.reserve(edges.size());
  std::unordered_map<int, int> component_required_count;
  component_required_count.reserve(nodes.size());

  for (int edge_id = 0; edge_id < static_cast<int>(edges.size()); ++edge_id) {
    const int root = components.find(edges[edge_id].u);
    component_edges[root].push_back(edge_id);
  }
  for (int node_id = 0; node_id < node_count; ++node_id) {
    if (!required[node_id]) {
      continue;
    }
    const int root = components.find(node_id);
    component_required_count[root]++;
  }

  std::vector<bool> keep_edge(edges.size(), false);
  for (const auto& [root, comp_edge_ids] : component_edges) {
    const auto it = component_required_count.find(root);
    const int required_count = (it == component_required_count.end()) ? 0
                                                                       : it->second;
    if (required_count <= 1) {
      for (int edge_id : comp_edge_ids) {
        keep_edge[edge_id] = true;
      }
      continue;
    }

    std::vector<int> comp_nodes;
    comp_nodes.reserve(comp_edge_ids.size() * 2);
    std::unordered_map<int, int> local_id;
    local_id.reserve(comp_edge_ids.size() * 2);

    for (int edge_id : comp_edge_ids) {
      const EdgeInfo& edge = edges[edge_id];
      if (local_id.find(edge.u) == local_id.end()) {
        local_id[edge.u] = comp_nodes.size();
        comp_nodes.push_back(edge.u);
      }
      if (local_id.find(edge.v) == local_id.end()) {
        local_id[edge.v] = comp_nodes.size();
        comp_nodes.push_back(edge.v);
      }
    }

    std::vector<int> sorted_edges = comp_edge_ids;
    std::sort(sorted_edges.begin(),
              sorted_edges.end(),
              [&](int lhs, int rhs) {
                return isBetterEdgeForMst(edges[lhs], edges[rhs]);
              });

    DisjointSet mst_sets(comp_nodes.size());
    std::vector<int> mst_edges;
    mst_edges.reserve(comp_nodes.size() > 0 ? comp_nodes.size() - 1 : 0);

    for (int edge_id : sorted_edges) {
      const EdgeInfo& edge = edges[edge_id];
      const int u_local = local_id[edge.u];
      const int v_local = local_id[edge.v];
      if (mst_sets.unite(u_local, v_local)) {
        mst_edges.push_back(edge_id);
      }
    }

    if (mst_edges.empty()) {
      for (int edge_id : comp_edge_ids) {
        keep_edge[edge_id] = true;
      }
      continue;
    }

    std::vector<std::vector<int>> incident(comp_nodes.size());
    std::vector<int> degree(comp_nodes.size(), 0);
    for (int mst_id = 0; mst_id < static_cast<int>(mst_edges.size()); ++mst_id) {
      const EdgeInfo& edge = edges[mst_edges[mst_id]];
      const int u_local = local_id[edge.u];
      const int v_local = local_id[edge.v];
      incident[u_local].push_back(mst_id);
      incident[v_local].push_back(mst_id);
      degree[u_local]++;
      degree[v_local]++;
    }

    std::vector<bool> mst_active(mst_edges.size(), true);
    std::deque<int> leaves;
    for (int node_local = 0; node_local < static_cast<int>(comp_nodes.size());
         ++node_local) {
      const int node_global = comp_nodes[node_local];
      if (degree[node_local] == 1 && !required[node_global]) {
        leaves.push_back(node_local);
      }
    }

    while (!leaves.empty()) {
      const int node_local = leaves.front();
      leaves.pop_front();
      if (degree[node_local] != 1) {
        continue;
      }
      const int node_global = comp_nodes[node_local];
      if (required[node_global]) {
        continue;
      }

      for (int mst_id : incident[node_local]) {
        if (!mst_active[mst_id]) {
          continue;
        }
        mst_active[mst_id] = false;
        const EdgeInfo& edge = edges[mst_edges[mst_id]];
        const int u_local = local_id[edge.u];
        const int v_local = local_id[edge.v];
        const int other_local = (u_local == node_local) ? v_local : u_local;
        degree[node_local]--;
        degree[other_local]--;
        const int other_global = comp_nodes[other_local];
        if (degree[other_local] == 1 && !required[other_global]) {
          leaves.push_back(other_local);
        }
        break;
      }
    }

    int active_edges = 0;
    for (int mst_id = 0; mst_id < static_cast<int>(mst_edges.size()); ++mst_id) {
      if (mst_active[mst_id]) {
        keep_edge[mst_edges[mst_id]] = true;
        active_edges++;
      }
    }
    if (active_edges == 0) {
      for (int edge_id : comp_edge_ids) {
        keep_edge[edge_id] = true;
      }
    }
  }

  GRoute optimized;
  optimized.reserve(route.size());
  std::unordered_set<GSegment, GSegmentHash> seen;
  seen.reserve(route.size());
  for (int edge_id = 0; edge_id < static_cast<int>(edges.size()); ++edge_id) {
    if (!keep_edge[edge_id]) {
      continue;
    }
    const GSegment normalized_segment = normalizeSegment(edges[edge_id].segment);
    if (seen.insert(normalized_segment).second) {
      optimized.push_back(normalized_segment);
    }
  }

  const RouteStats optimized_stats = computeRouteStats(optimized);
  const bool improves_wirelength
      = optimized_stats.wirelength < original_stats.wirelength;
  const bool same_wire_and_no_more_vias
      = optimized_stats.wirelength == original_stats.wirelength
        && optimized_stats.via_count <= original_stats.via_count;

  if (!optimized.empty() && (improves_wirelength || same_wire_and_no_more_vias)) {
    route.swap(optimized);
  }
}

void cleanupRouteSegments(
    NetRouteMap& routes,
    const std::unordered_map<odb::dbNet*, const std::vector<Pin>*>* net_pins)
{
  for (auto& [db_net, route] : routes) {
    if (net_pins != nullptr) {
      auto pins_it = net_pins->find(db_net);
      if (pins_it != net_pins->end() && pins_it->second != nullptr) {
        optimizeRouteTopology(*pins_it->second, route);
        continue;
      }
    }
    dedupeAndDropStubs(route);
  }
}

}  // namespace

NewGR::NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NewGR::~NewGR() = default;

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  active_backend_ = Backend::None;
  last_total_overflow_ = 0;
  if (nets.empty()) {
    return {};
  }

  // Start from FastRoute guides and run NEWGR-specific cleanup.
  NetRouteMap routes = grouter_->fastroute_->run();
  if (!routes.empty()) {
    grouter_->addRemainingGuides(
        routes, nets, min_routing_layer, max_routing_layer);
    grouter_->connectPadPins(routes);
    std::unordered_map<odb::dbNet*, const std::vector<Pin>*> net_pins;
    net_pins.reserve(routes.size());
    for (auto& net_route : routes) {
      std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
      GRoute& route = net_route.second;
      grouter_->mergeSegments(pins, route);
      net_pins.emplace(net_route.first, &pins);
    }
    cleanupRouteSegments(routes, &net_pins);
    active_backend_ = Backend::FastRoute;
    last_total_overflow_ = grouter_->fastroute_->totalOverflow();
    return routes;
  }

  if (!grouter_->hasSprouteGridData() || !grouter_->hasSprouteNetData()) {
    logger_->error(utl::GRT,
                   6003,
                   "NEWGR router selected, but grid/net data is not initialized.");
  }
  if (!engine_) {
    engine_ = std::make_unique<NewgrEngine>(logger_);
  }

  engine_->init(grouter_->sproute_grid_data_, grouter_->sproute_nets_);
  routes = engine_->run();

  grouter_->addRemainingGuides(routes, nets, min_routing_layer, max_routing_layer);
  grouter_->connectPadPins(routes);
  std::unordered_map<odb::dbNet*, const std::vector<Pin>*> net_pins;
  net_pins.reserve(routes.size());
  for (auto& net_route : routes) {
    std::vector<Pin>& pins = grouter_->db_net_map_[net_route.first]->getPins();
    GRoute& route = net_route.second;
    grouter_->mergeSegments(pins, route);
    net_pins.emplace(net_route.first, &pins);
  }
  cleanupRouteSegments(routes, &net_pins);
  active_backend_ = Backend::NewgrEngine;
  last_total_overflow_ = engine_->getTotalOverflow();

  return routes;
}

int NewGR::getTotalOverflow() const
{
  if (active_backend_ == Backend::NewgrEngine && engine_) {
    return engine_->getTotalOverflow();
  }
  return last_total_overflow_;
}

void NewGR::updateDbCongestion(odb::dbBlock* block)
{
  if (active_backend_ == Backend::FastRoute) {
    int min_routing_layer = 0;
    int max_routing_layer = 0;
    grouter_->getMinMaxLayer(min_routing_layer, max_routing_layer);
    grouter_->fastroute_->updateDbCongestion(min_routing_layer,
                                             max_routing_layer);
  } else if (active_backend_ == Backend::NewgrEngine && engine_) {
    engine_->updateDbCongestion(block);
  }
}

}  // namespace grt
