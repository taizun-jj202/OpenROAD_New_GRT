#include "newgr/CUGR.h"

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

namespace grt {
namespace newgr {

CUGR::CUGR(odb::dbDatabase* db,
           utl::Logger* log,
           stt::SteinerTreeBuilder* stt_builder)
    : db_(db), logger_(log), stt_builder_(stt_builder)
{
}

CUGR::~CUGR() = default;

void CUGR::init(const int min_routing_layer, const int max_routing_layer)
{
  // The OpenROAD API passes tech routing "levels" (1-based). Internally this
  // router uses 0-based indices, so we translate here and keep the rest of the
  // implementation consistent with the design/grid representation.
  constants_.min_routing_layer = std::max(0, min_routing_layer - 1);

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
    const int overflow_edges = grid_graph_->checkOverflow(net->getRoutingTree());
    if (overflow_edges >= constants_.reroute_overflow_edge_threshold) {
      netIndices.push_back(net->getIndex());
    }
  }
  const int num_nets = gr_nets_.size();
  logger_->report("{} / {} nets have overflow (threshold >= {}).",
                  netIndices.size(),
                  num_nets,
                  constants_.reroute_overflow_edge_threshold);
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
  // Slightly denser sparse grid to improve solution quality for overflow nets
  // while keeping runtime close to the fast baseline.
  SparseGrid grid(8, 8, 0, 0);
  for (const int netIndex : netIndices) {
    GRNet* net = gr_nets_[netIndex].get();
    const int hp = net->getBoundingBox().hp();
    if (hp >= 80) {
      grid.reset(5, 5);
    } else if (hp >= 40) {
      grid.reset(6, 6);
    } else {
      grid.reset(7, 7);
    }
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
  }

  updateOverflowNets(netIndices);
}

void CUGR::route()
{
  std::vector<int> netIndices;
  netIndices.reserve(gr_nets_.size());
  for (const auto& net : gr_nets_) {
    netIndices.push_back(net->getIndex());
  }

  patternRoute(netIndices);

  if (constants_.enable_detours) {
    patternRouteWithDetours(netIndices);
  }

  if (constants_.enable_maze) {
    mazeRoute(netIndices);
  }

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

    // `check_antennas` requires that every multi-term net has at least one
    // guide. Some nets can collapse to a single routing point at this coarse
    // gcell granularity, resulting in no edges (and thus no guides) unless we
    // emit a degenerate segment.
    if (route.empty()) {
      const auto& pin_aps = net->getPinAccessPoints();
      GRPoint anchor(constants_.min_routing_layer, 0, 0);
      bool found_anchor = false;
      for (const auto& aps : pin_aps) {
        if (!aps.empty()) {
          anchor = aps.front();
          found_anchor = true;
          break;
        }
      }
      if (found_anchor) {
        int layer_idx = std::max(constants_.min_routing_layer,
                                 anchor.getLayerIdx());
        layer_idx = std::min(layer_idx, grid_graph_->getNumLayers() - 1);
        const int x = grid_graph_->getGridline(0, anchor.x()) + half_gcell;
        const int y = grid_graph_->getGridline(1, anchor.y()) + half_gcell;
        route.emplace_back(x, y, layer_idx + 1, x, y, layer_idx + 1, false);
      }
    }
  }

  return routes;
}

void CUGR::sortNetIndices(std::vector<int>& netIndices) const
{
  // Route "harder" nets first to reduce later-stage overflows/detours.
  // In this implementation, bounding-box HPWL is a good proxy for difficulty.
  std::vector<int> halfParameters(gr_nets_.size());
  std::vector<int> pinCounts(gr_nets_.size());
  for (int netIndex : netIndices) {
    auto& net = gr_nets_[netIndex];
    halfParameters[netIndex] = net->getBoundingBox().hp();
    pinCounts[netIndex] = net->getNumPins();
  }
  sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
    if (halfParameters[lhs] != halfParameters[rhs]) {
      return halfParameters[lhs] > halfParameters[rhs];
    }
    return pinCounts[lhs] > pinCounts[rhs];
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

}  // namespace newgr
}  // namespace grt
