#include "GridGraph.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "CUGR.h"
#include "Design.h"
#include "GRNet.h"
#include "GRTree.h"
#include "geo.h"
#include "robin_hood.h"
#include "utl/Logger.h"

namespace grt::newgr {

GridGraph::GridGraph(const Design* design,
                     const Constants& constants,
                     utl::Logger* logger)
    : logger_(logger),
      gridlines_(design->getGridlines()),
      lib_dbu_(design->getLibDBU()),
      m2_pitch_(design->getLayer(1).getPitch()),
      num_layers_(design->getNumLayers()),
      x_size_(gridlines_[0].size() - 1),
      y_size_(gridlines_[1].size() - 1),
      constants_(constants)
{
  grid_centers_.resize(2);
  for (int dimension = 0; dimension <= 1; dimension++) {
    grid_centers_[dimension].resize(gridlines_[dimension].size() - 1);
    for (int gridIndex = 0; gridIndex < gridlines_[dimension].size() - 1;
         gridIndex++) {
      grid_centers_[dimension][gridIndex]
          = (gridlines_[dimension][gridIndex]
             + gridlines_[dimension][gridIndex + 1])
            / 2;
    }
  }

  layer_names_.resize(num_layers_);
  layer_directions_.resize(num_layers_);
  layer_min_lengths_.resize(num_layers_);
  for (int layer_index = 0; layer_index < num_layers_; layer_index++) {
    const auto& layer = design->getLayer(layer_index);
    layer_names_[layer_index] = layer.getName();
    layer_directions_[layer_index] = layer.getDirection();
    layer_min_lengths_[layer_index] = layer.getMinLength();
  }

  unit_length_wire_cost_ = design->getUnitLengthWireCost();
  unit_via_cost_ = design->getUnitViaCost();
  unit_length_short_costs_.resize(num_layers_);
  for (int layer_index = 0; layer_index < num_layers_; layer_index++) {
    unit_length_short_costs_[layer_index]
        = design->getUnitLengthShortCost(layer_index);
  }

  // Init grid graph edges
  std::vector<std::vector<int>> gridTracks(num_layers_);
  graph_edges_.assign(num_layers_,
                      std::vector<std::vector<GraphEdge>>(
                          x_size_, std::vector<GraphEdge>(y_size_)));
  for (int layer_index = 0; layer_index < num_layers_; layer_index++) {
    const MetalLayer& layer = design->getLayer(layer_index);
    const int direction = layer.getDirection();

    const int nGrids = gridlines_[1 - direction].size() - 1;
    gridTracks[layer_index].resize(nGrids);
    for (size_t gridIndex = 0; gridIndex < nGrids; gridIndex++) {
      IntervalT locRange(gridlines_[1 - direction][gridIndex],
                         gridlines_[1 - direction][gridIndex + 1]);
      auto trackRange = layer.rangeSearchTracks(locRange);
      if (trackRange.IsValid()) {
        gridTracks[layer_index][gridIndex] = trackRange.range() + 1;
        // exclude the track on the higher gridline
        if (gridIndex != nGrids - 1
            && layer.getTrackLocation(trackRange.high()) == locRange.high()) {
          gridTracks[layer_index][gridIndex]--;
        }
      } else {
        gridTracks[layer_index][gridIndex] = 0;
      }
    }

    // Initialize edges' capacity to the number of tracks
    if (direction == MetalLayer::V) {
      for (size_t x = 0; x < x_size_; x++) {
        const CapacityT nTracks = gridTracks[layer_index][x];
        for (size_t y = 0; y + 1 < y_size_; y++) {
          graph_edges_[layer_index][x][y].capacity = nTracks;
        }
      }
    } else {
      for (size_t y = 0; y < y_size_; y++) {
        const CapacityT nTracks = gridTracks[layer_index][y];
        for (size_t x = 0; x + 1 < x_size_; x++) {
          graph_edges_[layer_index][x][y].capacity = nTracks;
        }
      }
    }
  }

  // Deduct obstacles usage for layers EXCEPT Metal 1
  std::vector<std::vector<BoxT>> obstacles(num_layers_);
  design->getAllObstacles(obstacles, true);
  for (int layer_index = 1; layer_index < num_layers_; layer_index++) {
    const MetalLayer& layer = design->getLayer(layer_index);
    const int direction = layer.getDirection();
    const int nGrids = gridlines_[1 - direction].size() - 1;
    const int nEdges = gridlines_[direction].size() - 2;
    int minEdgeLength = std::numeric_limits<int>::max();
    for (int edge_index = 0; edge_index < nEdges; edge_index++) {
      minEdgeLength = std::min(minEdgeLength,
                               grid_centers_[direction][edge_index + 1]
                                   - grid_centers_[direction][edge_index]);
    }
    std::vector<std::vector<std::shared_ptr<std::pair<BoxT, IntervalT>>>>
        obstaclesInGrid(nGrids);  // obstacle indices sorted in track grids
    // Sort obstacles in track grids
    for (auto& obs : obstacles[layer_index]) {
      const int width = std::min(obs.x().range(), obs.y().range());
      const int spacing
          = layer.getParallelSpacing(
                width, std::min(minEdgeLength, obs[direction].range()))
            + layer.getWidth() / 2 - 1;
      PointT margin(0, 0);
      margin[1 - direction] = spacing;
      const BoxT obsBox(obs.lx() - margin.x(),
                        obs.ly() - margin.y(),
                        obs.hx() + margin.x(),
                        obs.hy() + margin.y());  // enlarged obstacle box
      const IntervalT trackRange
          = layer.rangeSearchTracks(obsBox[1 - direction]);
      const std::shared_ptr<std::pair<BoxT, IntervalT>> obstacle
          = std::make_shared<std::pair<BoxT, IntervalT>>(obsBox, trackRange);
      // Get grid range
      const IntervalT gridRange
          = rangeSearchRows(1 - direction, obsBox[1 - direction]);
      for (int gridIndex = gridRange.low(); gridIndex <= gridRange.high();
           gridIndex++) {
        obstaclesInGrid[gridIndex].push_back(obstacle);
      }
    }
    // Handle each track grid
    IntervalT gridTrackRange;
    for (int gridIndex = 0; gridIndex < nGrids; gridIndex++) {
      if (gridIndex == 0) {
        gridTrackRange.Set(0, gridTracks[layer_index][gridIndex] - 1);
      } else {
        gridTrackRange.SetLow(gridTrackRange.high() + 1);
        gridTrackRange.addToHigh(gridTracks[layer_index][gridIndex]);
      }
      if (!gridTrackRange.IsValid()) {
        continue;
      }
      if (obstaclesInGrid[gridIndex].empty()) {
        continue;
      }
      std::vector<std::vector<std::shared_ptr<std::pair<BoxT, IntervalT>>>>
          obstaclesAtEdge(nEdges);
      for (auto& obstacle : obstaclesInGrid[gridIndex]) {
        const IntervalT gridlineRange
            = rangeSearchGridlines(direction, obstacle->first[direction]);
        const IntervalT edgeRange(std::max(gridlineRange.low() - 2, 0),
                                  std::min(gridlineRange.high(), nEdges - 1));
        for (int edge_index = edgeRange.low(); edge_index <= edgeRange.high();
             edge_index++) {
          obstaclesAtEdge[edge_index].emplace_back(obstacle);
        }
      }
      for (int edge_index = 0; edge_index < nEdges; edge_index++) {
        if (obstaclesAtEdge[edge_index].empty()) {
          continue;
        }
        const int gridline = gridlines_[direction][edge_index + 1];
        const IntervalT edgeInterval(grid_centers_[direction][edge_index],
                                     grid_centers_[direction][edge_index + 1]);
        // Update cpacity
        std::vector<IntervalT> usable_intervals(gridTrackRange.range() + 1,
                                                edgeInterval);
        for (auto& obstacle : obstaclesAtEdge[edge_index]) {
          const IntervalT affectedTrackRange
              = gridTrackRange.IntersectWith(obstacle->second);
          if (!affectedTrackRange.IsValid()) {
            continue;
          }
          for (int trackIndex = affectedTrackRange.low();
               trackIndex <= affectedTrackRange.high();
               trackIndex++) {
            const int tIdx = trackIndex - gridTrackRange.low();
            if (obstacle->first[direction].low() <= gridline
                && obstacle->first[direction].high() >= gridline) {
              // Completely blocked
              usable_intervals[tIdx] = {gridline, gridline};
            } else if (obstacle->first[direction].high() < gridline) {
              usable_intervals[tIdx].SetLow(
                  std::max(usable_intervals[tIdx].low(),
                           obstacle->first[direction].high()));
            } else if (obstacle->first[direction].low() > gridline) {
              usable_intervals[tIdx].SetHigh(
                  std::min(usable_intervals[tIdx].high(),
                           obstacle->first[direction].low()));
            }
          }
        }
        CapacityT capacity = 0;
        for (IntervalT& usable_interval : usable_intervals) {
          capacity
              += (CapacityT) usable_interval.range() / edgeInterval.range();
        }
        if (direction == MetalLayer::V) {
          graph_edges_[layer_index][gridIndex][edge_index].capacity = capacity;
        } else {
          graph_edges_[layer_index][edge_index][gridIndex].capacity = capacity;
        }
      }
    }
  }

  // Apply user-defined capacity adjustment
  for (int layer_index = 0; layer_index < num_layers_; layer_index++) {
    const float adjustment = design->getLayer(layer_index).getAdjustment();
    if (adjustment != 0.0) {
      for (size_t x = 0; x < x_size_; x++) {
        for (size_t y = 0; y < y_size_; y++) {
          graph_edges_[layer_index][x][y].capacity *= (1.0 - adjustment);
        }
      }
    }
  }
}

IntervalT GridGraph::rangeSearchGridlines(const int dimension,
                                          const IntervalT& loc_interval) const
{
  IntervalT range;
  range.Set(lower_bound(gridlines_[dimension].begin(),
                        gridlines_[dimension].end(),
                        loc_interval.low())
                - gridlines_[dimension].begin(),
            lower_bound(gridlines_[dimension].begin(),
                        gridlines_[dimension].end(),
                        loc_interval.high())
                - gridlines_[dimension].begin());
  if (range.high() >= gridlines_[dimension].size()) {
    range.SetHigh(gridlines_[dimension].size() - 1);
  } else if (gridlines_[dimension][range.high()] > loc_interval.high()) {
    range.addToHigh(-1);
  }
  return range;
}

IntervalT GridGraph::rangeSearchRows(const int dimension,
                                     const IntervalT& loc_interval) const
{
  const auto& lineRange = rangeSearchGridlines(dimension, loc_interval);
  return {gridlines_[dimension][lineRange.low()] == loc_interval.low()
              ? lineRange.low()
              : std::max(lineRange.low() - 1, 0),
          gridlines_[dimension][lineRange.high()] == loc_interval.high()
              ? lineRange.high() - 1
              : std::min(lineRange.high(),
                         static_cast<int>(getSize(dimension)) - 1)};
}

BoxT GridGraph::getCellBox(PointT point) const
{
  return {getGridline(0, point.x()),
          getGridline(1, point.y()),
          getGridline(0, point.x() + 1),
          getGridline(1, point.y() + 1)};
}

BoxT GridGraph::rangeSearchCells(const BoxT& box) const
{
  return {rangeSearchRows(0, box[0]), rangeSearchRows(1, box[1])};
}

int GridGraph::getEdgeLength(int direction, int edge_index) const
{
  return grid_centers_[direction][edge_index + 1]
         - grid_centers_[direction][edge_index];
}

double GridGraph::logistic(const CapacityT& input, const double slope) const
{
  return 1.0 / (1.0 + exp(input * slope));
}

CapacityT GridGraph::getSoftCapacity(const GraphEdge& edge) const
{
  if (!constants_.use_soft_capacity) {
    return edge.capacity;
  }
  if (edge.capacity < 1.0) {
    return edge.capacity;
  }
  const double util = edge.demand / std::max(edge.capacity, 1.0);
  const double ratio = constants_.soft_cap_min_ratio
                       + (constants_.soft_cap_max_ratio
                          - constants_.soft_cap_min_ratio)
                             / (1.0 + exp((util - constants_.soft_cap_mid_util)
                                          * constants_.soft_cap_slope));
  return edge.capacity * ratio;
}

CostT GridGraph::getCongestionPenalty(const CapacityT reference_capacity,
                                      const CapacityT demand,
                                      const double slope) const
{
  if (reference_capacity < 1.0) {
    return 1.0;
  }

  const double util = demand / std::max(reference_capacity, 1.0);
  const double threshold = std::clamp(
      constants_.wl_relaxation_util_threshold, 0.0, 0.99);
  const double floor
      = std::clamp(constants_.wl_relaxation_penalty_floor, 0.0, 1.0);

  if (demand <= reference_capacity) {
    if (util <= threshold) {
      return floor * logistic(reference_capacity - demand, slope);
    }
    const double span = std::max(1e-6, 1.0 - threshold);
    const double ramp = std::clamp((util - threshold) / span, 0.0, 1.0);
    return floor + (1.0 - floor) * ramp * ramp;
  }

  const double overflow = demand - reference_capacity;
  return 1.0 + constants_.overflow_linear_penalty * overflow;
}

CostT GridGraph::getCongestionPenalty(const GraphEdge& edge,
                                      const double slope) const
{
  const CapacityT reference_capacity = getSoftCapacity(edge);
  return getCongestionPenalty(reference_capacity, edge.demand, slope);
}

CostT GridGraph::getWireCost(const int layer_index,
                             const PointT lower,
                             const CapacityT demand) const
{
  const int direction = layer_directions_[layer_index];
  const int edgeLength = getEdgeLength(direction, lower[direction]);
  const int demandLength = demand * edgeLength;
  const auto& edge = graph_edges_[layer_index][lower.x()][lower.y()];
  CostT cost = demandLength * unit_length_wire_cost_;
  cost += demandLength * unit_length_short_costs_[layer_index]
          * getCongestionPenalty(edge, constants_.cost_logistic_slope);
  return cost;
}

CostT GridGraph::getWireCost(const int layer_index,
                             const PointT u,
                             const PointT v) const
{
  const int direction = layer_directions_[layer_index];
  assert(u[1 - direction] == v[1 - direction]);
  CostT cost = 0;
  if (direction == MetalLayer::H) {
    const auto [l, h] = std::minmax({u.x(), v.x()});
    for (int x = l; x < h; x++) {
      cost += getWireCost(layer_index, {x, u.y()});
    }
  } else {
    const auto [l, h] = std::minmax({u.y(), v.y()});
    for (int y = l; y < h; y++) {
      cost += getWireCost(layer_index, {u.x(), y});
    }
  }
  return cost;
}

CostT GridGraph::getViaCost(const int layer_index, const PointT loc) const
{
  assert(layer_index + 1 < num_layers_);
  CostT cost = unit_via_cost_;
  // Estimated wire cost to satisfy min-area
  for (int l = layer_index; l <= layer_index + 1; l++) {
    const int direction = layer_directions_[l];
    PointT lowerLoc = loc;
    lowerLoc[direction] -= 1;
    const int lowerEdgeLength
        = loc[direction] > 0 ? getEdgeLength(direction, lowerLoc[direction])
                             : 0;
    const int higherEdgeLength = loc[direction] < getSize(direction) - 1
                                     ? getEdgeLength(direction, loc[direction])
                                     : 0;
    const CapacityT demand = (CapacityT) layer_min_lengths_[l]
                             / (lowerEdgeLength + higherEdgeLength)
                             * constants_.via_multiplier;
    if (lowerEdgeLength > 0) {
      cost += getWireCost(l, lowerLoc, demand);
    }
    if (higherEdgeLength > 0) {
      cost += getWireCost(l, loc, demand);
    }
  }
  return cost;
}

AccessPointSet GridGraph::selectAccessPoints(const GRNet* net) const
{
  AccessPointHash hasher(y_size_);
  AccessPointSet selected_access_points(0, hasher);
  // cell hash (2d) -> access point, fixed layer interval
  selected_access_points.reserve(net->getNumPins());
  const auto& pin_access_points = net->getPinAccessPoints();
  if (pin_access_points.empty()) {
    return selected_access_points;
  }

  const auto& bounding_box = net->getBoundingBox();
  const PointT net_center(bounding_box.cx(), bounding_box.cy());

  auto getAccessibility = [&](const GRPoint& point) {
    int accessibility = 0;
    if (point.getLayerIdx() >= constants_.min_routing_layer) {
      const int direction = getLayerDirection(point.getLayerIdx());
      accessibility
          += getEdge(point.getLayerIdx(), point.x(), point.y()).capacity >= 1;
      if (point[direction] > 0) {
        auto lower = point;
        lower[direction] -= 1;
        accessibility
            += getEdge(lower.getLayerIdx(), lower.x(), lower.y()).capacity >= 1;
      }
    } else {
      accessibility = 1;
    }
    return accessibility;
  };

  std::vector<int> selected_indices(pin_access_points.size(), -1);
  std::vector<PointT> selected_points(pin_access_points.size(), PointT(0, 0));
  std::vector<int> min_allowed_accessibility(pin_access_points.size(), 0);
  for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
    int max_pin_accessibility = -1;
    for (const auto& point : pin_access_points[pin_index]) {
      max_pin_accessibility = std::max(max_pin_accessibility, getAccessibility(point));
    }
    min_allowed_accessibility[pin_index]
        = std::max(0, max_pin_accessibility - 1);
  }

  auto evaluateHpwl = [&](const int pin_to_replace, const GRPoint& replacement) {
    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    bool has_point = false;
    for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
      if (pin_index == pin_to_replace) {
        min_x = std::min(min_x, replacement.x());
        max_x = std::max(max_x, replacement.x());
        min_y = std::min(min_y, replacement.y());
        max_y = std::max(max_y, replacement.y());
        has_point = true;
        continue;
      }
      const int selected_index = selected_indices[pin_index];
      if (selected_index < 0) {
        continue;
      }
      const auto& point = pin_access_points[pin_index][selected_index];
      min_x = std::min(min_x, point.x());
      max_x = std::max(max_x, point.x());
      min_y = std::min(min_y, point.y());
      max_y = std::max(max_y, point.y());
      has_point = true;
    }
    if (!has_point) {
      return 0;
    }
    return (max_x - min_x) + (max_y - min_y);
  };

  auto evaluateMstWire = [&](const int pin_to_replace,
                             const GRPoint& replacement) {
    std::vector<PointT> points;
    points.reserve(pin_access_points.size());
    for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
      if (pin_index == pin_to_replace) {
        points.emplace_back(replacement.x(), replacement.y());
        continue;
      }
      const int selected_index = selected_indices[pin_index];
      if (selected_index < 0) {
        continue;
      }
      const auto& point = pin_access_points[pin_index][selected_index];
      points.emplace_back(point.x(), point.y());
    }
    if (points.size() <= 1) {
      return 0;
    }

    const int n = points.size();
    std::vector<int> best(n, std::numeric_limits<int>::max());
    std::vector<bool> used(n, false);
    best[0] = 0;
    int total = 0;
    for (int iter = 0; iter < n; iter++) {
      int v = -1;
      int best_cost = std::numeric_limits<int>::max();
      for (int i = 0; i < n; i++) {
        if (!used[i] && best[i] < best_cost) {
          best_cost = best[i];
          v = i;
        }
      }
      if (v < 0) {
        break;
      }
      used[v] = true;
      total += best[v];
      for (int u = 0; u < n; u++) {
        if (used[u]) {
          continue;
        }
        const int dist = std::abs(points[v].x() - points[u].x())
                         + std::abs(points[v].y() - points[u].y());
        if (dist < best[u]) {
          best[u] = dist;
        }
      }
    }
    return total;
  };

  // First pass: keep accessibility but bias toward compact net shape.
  for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
    const auto& access_points = pin_access_points[pin_index];
    int best_index = -1;
    int best_accessibility = -1;
    int64_t best_score = std::numeric_limits<int64_t>::max();
    for (int index = 0; index < access_points.size(); index++) {
      const GRPoint& point = access_points[index];
      const int accessibility = getAccessibility(point);
      const int center_dist
          = std::abs(net_center.x() - point.x()) + std::abs(net_center.y() - point.y());
      const int hpwl = evaluateHpwl(pin_index, point);
      const int mst_wire = evaluateMstWire(pin_index, point);
      const int64_t score = static_cast<int64_t>(hpwl) * 4096
                            + static_cast<int64_t>(mst_wire) * 1024
                            + static_cast<int64_t>(center_dist) * 16;
      if (accessibility > best_accessibility
          || (accessibility == best_accessibility && score < best_score)) {
        best_index = index;
        best_accessibility = accessibility;
        best_score = score;
      }
    }

    if (best_accessibility <= 0) {
      logger_->warn(utl::GRT, 7001, "pin is hard to access.");
    }
    if (best_index >= 0) {
      selected_indices[pin_index] = best_index;
      selected_points[pin_index] = access_points[best_index];
    }
  }

  // FastRoute-style topology shaping with SPRoute-like accessibility guardband:
  // allow at most one accessibility level drop to obtain shorter trees.
  const int max_refine_passes = 7;
  for (int pass = 0; pass < max_refine_passes; pass++) {
    std::vector<int> xs;
    std::vector<int> ys;
    xs.reserve(selected_points.size());
    ys.reserve(selected_points.size());
    for (int pin_index = 0; pin_index < selected_indices.size(); pin_index++) {
      if (selected_indices[pin_index] < 0) {
        continue;
      }
      xs.push_back(selected_points[pin_index].x());
      ys.push_back(selected_points[pin_index].y());
    }
    if (xs.empty()) {
      break;
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());
    const int median_x = xs[xs.size() / 2];
    const int median_y = ys[ys.size() / 2];

    bool changed = false;
    for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
      const auto& access_points = pin_access_points[pin_index];
      if (access_points.empty()) {
        continue;
      }

      int best_index = selected_indices[pin_index];
      int best_accessibility = std::numeric_limits<int>::min();
      int best_hpwl = std::numeric_limits<int>::max();
      int best_mst_wire = std::numeric_limits<int>::max();
      int64_t best_score = std::numeric_limits<int64_t>::max();
      for (int index = 0; index < access_points.size(); index++) {
        const GRPoint& point = access_points[index];
        const int accessibility = getAccessibility(point);
        if (accessibility < min_allowed_accessibility[pin_index]) {
          continue;
        }
        const int hpwl = evaluateHpwl(pin_index, point);
        const int mst_wire = evaluateMstWire(pin_index, point);
        const int median_dist
            = std::abs(median_x - point.x()) + std::abs(median_y - point.y());
        const int center_dist
            = std::abs(net_center.x() - point.x()) + std::abs(net_center.y() - point.y());
        const int64_t score = static_cast<int64_t>(hpwl) * 4096
                              + static_cast<int64_t>(mst_wire) * 1024
                              + static_cast<int64_t>(median_dist) * 32
                              + static_cast<int64_t>(center_dist) * 8;
        if (score < best_score
            || (score == best_score && accessibility > best_accessibility)
            || (score == best_score && accessibility == best_accessibility
                && hpwl < best_hpwl)
            || (score == best_score && accessibility == best_accessibility
                && hpwl == best_hpwl && mst_wire < best_mst_wire)) {
          best_index = index;
          best_accessibility = accessibility;
          best_hpwl = hpwl;
          best_mst_wire = mst_wire;
          best_score = score;
        }
      }

      if (best_index < 0) {
        for (int index = 0; index < access_points.size(); index++) {
          const GRPoint& point = access_points[index];
          const int accessibility = getAccessibility(point);
          const int hpwl = evaluateHpwl(pin_index, point);
          if (accessibility > best_accessibility
              || (accessibility == best_accessibility && hpwl < best_hpwl)) {
            best_index = index;
            best_accessibility = accessibility;
            best_hpwl = hpwl;
          }
        }
      }

      if (best_index >= 0 && best_index != selected_indices[pin_index]) {
        selected_indices[pin_index] = best_index;
        selected_points[pin_index] = access_points[best_index];
        changed = true;
      }
    }
    if (!changed) {
      break;
    }
  }

  auto computeGlobalShapeScore = [&](const std::vector<PointT>& points) {
    std::vector<PointT> active_points;
    active_points.reserve(points.size());
    int min_x = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int min_y = std::numeric_limits<int>::max();
    int max_y = std::numeric_limits<int>::min();
    int64_t center_dist_sum = 0;
    for (int pin_index = 0; pin_index < selected_indices.size(); pin_index++) {
      if (selected_indices[pin_index] < 0) {
        continue;
      }
      const auto& point = points[pin_index];
      active_points.push_back(point);
      min_x = std::min(min_x, point.x());
      max_x = std::max(max_x, point.x());
      min_y = std::min(min_y, point.y());
      max_y = std::max(max_y, point.y());
      center_dist_sum += std::abs(net_center.x() - point.x())
                         + std::abs(net_center.y() - point.y());
    }
    if (active_points.size() <= 1) {
      return static_cast<int64_t>(0);
    }

    const int hpwl = (max_x - min_x) + (max_y - min_y);
    const int n = active_points.size();
    std::vector<int> best(n, std::numeric_limits<int>::max());
    std::vector<bool> used(n, false);
    best[0] = 0;
    int mst_wire = 0;
    for (int iter = 0; iter < n; iter++) {
      int v = -1;
      int best_cost = std::numeric_limits<int>::max();
      for (int i = 0; i < n; i++) {
        if (!used[i] && best[i] < best_cost) {
          best_cost = best[i];
          v = i;
        }
      }
      if (v < 0) {
        break;
      }
      used[v] = true;
      mst_wire += best[v];
      for (int u = 0; u < n; u++) {
        if (used[u]) {
          continue;
        }
        const int dist = std::abs(active_points[v].x() - active_points[u].x())
                         + std::abs(active_points[v].y() - active_points[u].y());
        if (dist < best[u]) {
          best[u] = dist;
        }
      }
    }

    return static_cast<int64_t>(hpwl) * 4096
           + static_cast<int64_t>(mst_wire) * 1024 + center_dist_sum * 8;
  };

  // Pairwise local search on the most topologically influential pins.
  // This escapes one-pin-at-a-time local minima without exploding runtime.
  if (pin_access_points.size() >= 4) {
    struct PinPriority
    {
      int pin_index;
      int center_dist;
      int num_options;
    };

    std::vector<PinPriority> variable_pins;
    variable_pins.reserve(pin_access_points.size());
    for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
      if (selected_indices[pin_index] < 0 || pin_access_points[pin_index].size() <= 1) {
        continue;
      }
      const auto& point = selected_points[pin_index];
      const int center_dist
          = std::abs(net_center.x() - point.x()) + std::abs(net_center.y() - point.y());
      variable_pins.push_back(
          {pin_index, center_dist, static_cast<int>(pin_access_points[pin_index].size())});
    }

    std::sort(variable_pins.begin(),
              variable_pins.end(),
              [](const PinPriority& lhs, const PinPriority& rhs) {
                if (lhs.center_dist != rhs.center_dist) {
                  return lhs.center_dist > rhs.center_dist;
                }
                if (lhs.num_options != rhs.num_options) {
                  return lhs.num_options > rhs.num_options;
                }
                return lhs.pin_index < rhs.pin_index;
              });

    const int max_pair_pins
        = std::min(static_cast<int>(variable_pins.size()), 8);
    if (max_pair_pins >= 2) {
      std::vector<int> focus_pins;
      focus_pins.reserve(max_pair_pins);
      for (int i = 0; i < max_pair_pins; i++) {
        focus_pins.push_back(variable_pins[i].pin_index);
      }

      auto getCandidateIndices = [&](const int pin_index) {
        struct RankedOption
        {
          int index;
          int accessibility;
          int center_dist;
        };
        std::vector<RankedOption> ranked;
        ranked.reserve(pin_access_points[pin_index].size());
        for (int index = 0; index < pin_access_points[pin_index].size(); index++) {
          const auto& point = pin_access_points[pin_index][index];
          const int accessibility = getAccessibility(point);
          if (accessibility < min_allowed_accessibility[pin_index]) {
            continue;
          }
          const int center_dist = std::abs(net_center.x() - point.x())
                                  + std::abs(net_center.y() - point.y());
          ranked.push_back({index, accessibility, center_dist});
        }
        std::sort(ranked.begin(),
                  ranked.end(),
                  [](const RankedOption& lhs, const RankedOption& rhs) {
                    if (lhs.accessibility != rhs.accessibility) {
                      return lhs.accessibility > rhs.accessibility;
                    }
                    if (lhs.center_dist != rhs.center_dist) {
                      return lhs.center_dist < rhs.center_dist;
                    }
                    return lhs.index < rhs.index;
                  });

        std::vector<int> indices;
        constexpr int kMaxCandidatesPerPin = 4;
        indices.reserve(kMaxCandidatesPerPin + 1);
        for (int i = 0;
             i < static_cast<int>(ranked.size()) && i < kMaxCandidatesPerPin;
             i++) {
          indices.push_back(ranked[i].index);
        }
        const int selected_index = selected_indices[pin_index];
        if (selected_index >= 0
            && std::find(indices.begin(), indices.end(), selected_index)
                   == indices.end()) {
          indices.push_back(selected_index);
        }
        return indices;
      };

      int64_t global_best_score = computeGlobalShapeScore(selected_points);
      bool improved = false;
      constexpr int kPairBudget = 12;
      for (int pass = 0; pass < 2; pass++) {
        int pair_trials = 0;
        bool improved_in_pass = false;
        for (int lhs_idx = 0; lhs_idx < static_cast<int>(focus_pins.size());
             lhs_idx++) {
          for (int rhs_idx = lhs_idx + 1;
               rhs_idx < static_cast<int>(focus_pins.size());
               rhs_idx++) {
            if (pair_trials >= kPairBudget) {
              break;
            }
            pair_trials++;

            const int lhs_pin = focus_pins[lhs_idx];
            const int rhs_pin = focus_pins[rhs_idx];
            const auto lhs_options = getCandidateIndices(lhs_pin);
            const auto rhs_options = getCandidateIndices(rhs_pin);
            if (lhs_options.empty() || rhs_options.empty()) {
              continue;
            }

            int best_lhs_index = selected_indices[lhs_pin];
            int best_rhs_index = selected_indices[rhs_pin];
            int best_access_sum = -1;
            int64_t best_pair_score = global_best_score;
            const PointT orig_lhs = selected_points[lhs_pin];
            const PointT orig_rhs = selected_points[rhs_pin];

            for (const int lhs_option : lhs_options) {
              const auto& lhs_point = pin_access_points[lhs_pin][lhs_option];
              const int lhs_access = getAccessibility(lhs_point);
              for (const int rhs_option : rhs_options) {
                if (lhs_option == selected_indices[lhs_pin]
                    && rhs_option == selected_indices[rhs_pin]) {
                  continue;
                }
                const auto& rhs_point = pin_access_points[rhs_pin][rhs_option];
                const int rhs_access = getAccessibility(rhs_point);
                selected_points[lhs_pin] = {lhs_point.x(), lhs_point.y()};
                selected_points[rhs_pin] = {rhs_point.x(), rhs_point.y()};
                const int64_t score = computeGlobalShapeScore(selected_points);
                const int access_sum = lhs_access + rhs_access;
                if (score < best_pair_score
                    || (score == best_pair_score && access_sum > best_access_sum)) {
                  best_pair_score = score;
                  best_access_sum = access_sum;
                  best_lhs_index = lhs_option;
                  best_rhs_index = rhs_option;
                }
              }
            }
            selected_points[lhs_pin] = orig_lhs;
            selected_points[rhs_pin] = orig_rhs;

            if ((best_lhs_index != selected_indices[lhs_pin]
                 || best_rhs_index != selected_indices[rhs_pin])
                && best_pair_score <= global_best_score) {
              selected_indices[lhs_pin] = best_lhs_index;
              selected_indices[rhs_pin] = best_rhs_index;
              const auto& best_lhs_point = pin_access_points[lhs_pin][best_lhs_index];
              const auto& best_rhs_point = pin_access_points[rhs_pin][best_rhs_index];
              selected_points[lhs_pin]
                  = {best_lhs_point.x(), best_lhs_point.y()};
              selected_points[rhs_pin]
                  = {best_rhs_point.x(), best_rhs_point.y()};
              global_best_score = best_pair_score;
              improved = true;
              improved_in_pass = true;
            }
          }
          if (pair_trials >= kPairBudget) {
            break;
          }
        }
        if (!improved_in_pass) {
          break;
        }
      }

      if (improved) {
        for (const int pin_index : focus_pins) {
          const int selected_index = selected_indices[pin_index];
          if (selected_index < 0) {
            continue;
          }
          const auto& point = pin_access_points[pin_index][selected_index];
          selected_points[pin_index] = {point.x(), point.y()};
        }
      }
    }
  }

  for (int pin_index = 0; pin_index < pin_access_points.size(); pin_index++) {
    const auto& access_points = pin_access_points[pin_index];
    const int selected_index = selected_indices[pin_index];
    if (selected_index < 0 || selected_index >= access_points.size()) {
      continue;
    }
    const PointT selected_point = access_points[selected_index];
    const AccessPoint ap{selected_point, {}};
    auto it = selected_access_points.emplace(ap).first;
    IntervalT& fixed_layer_interval = it->layers;
    for (const auto& point : access_points) {
      if (point.x() == selected_point.x() && point.y() == selected_point.y()) {
        fixed_layer_interval.Update(point.getLayerIdx());
      }
    }
  }

  // Extend the fixed layers to 2 layers higher to facilitate track switching
  for (auto& accessPoint : selected_access_points) {
    IntervalT& fixedLayers = accessPoint.layers;
    fixedLayers.SetHigh(
        std::min(fixedLayers.high() + 2, (int) getNumLayers() - 1));
  }
  return selected_access_points;
}

void GridGraph::commit(const int layer_index,
                       const PointT lower,
                       const CapacityT demand)
{
  graph_edges_[layer_index][lower.x()][lower.y()].demand += demand;
}

void GridGraph::commitWire(const int layer_index,
                           const PointT lower,
                           const bool rip_up)
{
  const int direction = layer_directions_[layer_index];
  const int edgeLength = getEdgeLength(direction, lower[direction]);
  if (rip_up) {
    commit(layer_index, lower, -1);
    total_length_ -= edgeLength;
  } else {
    commit(layer_index, lower, 1);
    total_length_ += edgeLength;
  }
}

void GridGraph::commitVia(const int layer_index,
                          const PointT loc,
                          const bool rip_up)
{
  assert(layer_index + 1 < num_layers_);
  for (int l = layer_index; l <= layer_index + 1; l++) {
    const int direction = layer_directions_[l];
    PointT lowerLoc = loc;
    lowerLoc[direction] -= 1;
    const int lowerEdgeLength
        = loc[direction] > 0 ? getEdgeLength(direction, lowerLoc[direction])
                             : 0;
    const int higherEdgeLength = loc[direction] < getSize(direction) - 1
                                     ? getEdgeLength(direction, loc[direction])
                                     : 0;
    const CapacityT demand = (CapacityT) layer_min_lengths_[l]
                             / (lowerEdgeLength + higherEdgeLength)
                             * constants_.via_multiplier;
    if (lowerEdgeLength > 0) {
      commit(l, lowerLoc, (rip_up ? -demand : demand));
    }
    if (higherEdgeLength > 0) {
      commit(l, loc, (rip_up ? -demand : demand));
    }
  }
  if (rip_up) {
    total_num_vias_ -= 1;
  } else {
    total_num_vias_ += 1;
  }
}

void GridGraph::commitTree(const std::shared_ptr<GRTreeNode>& tree,
                           const bool rip_up)
{
  GRTreeNode::preorder(tree, [&](const std::shared_ptr<GRTreeNode>& node) {
    for (const auto& child : node->getChildren()) {
      if (node->getLayerIdx() == child->getLayerIdx()) {
        const int direction = layer_directions_[node->getLayerIdx()];
        if (direction == MetalLayer::H) {
          assert(node->y() == child->y());
          const auto [l, h] = std::minmax({node->x(), child->x()});
          for (int x = l; x < h; x++) {
            commitWire(node->getLayerIdx(), {x, node->y()}, rip_up);
          }
        } else {
          assert(node->x() == child->x());
          const auto [l, h] = std::minmax({node->y(), child->y()});
          for (int y = l; y < h; y++) {
            commitWire(node->getLayerIdx(), {node->x(), y}, rip_up);
          }
        }
      } else {
        const int maxLayerIndex
            = std::max(node->getLayerIdx(), child->getLayerIdx());
        for (int layerIdx = std::min(node->getLayerIdx(), child->getLayerIdx());
             layerIdx < maxLayerIndex;
             layerIdx++) {
          commitVia(layerIdx, {node->x(), node->y()}, rip_up);
        }
      }
    }
  });
}

int GridGraph::checkOverflow(const int layer_index,
                             const PointT u,
                             const PointT v) const
{
  int num = 0;
  const int direction = layer_directions_[layer_index];
  if (direction == MetalLayer::H) {
    assert(u.y() == v.y());
    const auto [l, h] = std::minmax({u.x(), v.x()});
    for (int x = l; x < h; x++) {
      if (checkOverflow(layer_index, x, u.y())) {
        num++;
      }
    }
  } else {
    assert(u.x() == v.x());
    const auto [l, h] = std::minmax({u.y(), v.y()});
    for (int y = l; y < h; y++) {
      if (checkOverflow(layer_index, u.x(), y)) {
        num++;
      }
    }
  }
  return num;
}

int GridGraph::checkOverflow(const std::shared_ptr<GRTreeNode>& tree) const
{
  if (!tree) {
    return 0;
  }
  int num = 0;
  GRTreeNode::preorder(tree, [&](const std::shared_ptr<GRTreeNode>& node) {
    for (auto& child : node->getChildren()) {
      // Only check wires
      if (node->getLayerIdx() == child->getLayerIdx()) {
        num += checkOverflow(
            node->getLayerIdx(), (PointT) *node, (PointT) *child);
      }
    }
  });
  return num;
}

CapacityT GridGraph::getTotalOverflow() const
{
  CapacityT overflow = 0.0;
  for (int layer_index = constants_.min_routing_layer;
       layer_index < num_layers_;
       layer_index++) {
    const int direction = layer_directions_[layer_index];
    for (int x = 0; x < x_size_ - 1 + direction; x++) {
      for (int y = 0; y < y_size_ - direction; y++) {
        const GraphEdge& edge = graph_edges_[layer_index][x][y];
        if (edge.demand > edge.capacity) {
          overflow += edge.demand - edge.capacity;
        }
      }
    }
  }
  return overflow;
}

std::string GridGraph::getPythonString(
    const std::shared_ptr<GRTreeNode>& routing_tree) const
{
  std::vector<std::tuple<PointT, PointT, bool>> edges;
  GRTreeNode::preorder(
      routing_tree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            const int direction = getLayerDirection(node->getLayerIdx());
            const int r = (*node)[1 - direction];
            const int l = std::min((*node)[direction], (*child)[direction]);
            const int h = std::max((*node)[direction], (*child)[direction]);
            if (l == h) {
              continue;
            }
            PointT lpoint
                = (direction == MetalLayer::H ? PointT(l, r) : PointT(r, l));
            const PointT hpoint
                = (direction == MetalLayer::H ? PointT(h, r) : PointT(r, h));
            bool congested = false;
            for (int c = l; c < h; c++) {
              const PointT cpoint
                  = (direction == MetalLayer::H ? PointT(c, r) : PointT(r, c));
              if (checkOverflow(node->getLayerIdx(), cpoint.x(), cpoint.y())
                  != congested) {
                if (lpoint != cpoint) {
                  edges.emplace_back(lpoint, cpoint, congested);
                  lpoint = cpoint;
                }
                congested = !congested;
              }
            }
            if (lpoint != hpoint) {
              edges.emplace_back(lpoint, hpoint, congested);
            }
          }
        }
      });
  std::stringstream ss;
  ss << "[";
  for (int i = 0; i < edges.size(); i++) {
    const auto& edge = edges[i];
    ss << "[" << std::get<0>(edge) << ", " << std::get<1>(edge) << ", "
       << (std::get<2>(edge) ? 1 : 0) << "]";
    ss << (i < edges.size() - 1 ? ", " : "]");
  }
  return ss.str();
}

void GridGraph::extractBlockageView(GridGraphView<bool>& view) const
{
  view.assign(2,
              std::vector<std::vector<bool>>(x_size_,
                                             std::vector<bool>(y_size_, true)));
  for (int layer_index = constants_.min_routing_layer;
       layer_index < num_layers_;
       layer_index++) {
    const int direction = getLayerDirection(layer_index);
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        if (getEdge(layer_index, x, y).capacity >= 1.0) {
          view[direction][x][y] = false;
        }
      }
    }
  }
}

void GridGraph::extractCongestionView(GridGraphView<bool>& view) const
{
  view.assign(2,
              std::vector<std::vector<bool>>(
                  x_size_, std::vector<bool>(y_size_, false)));
  for (int layer_index = constants_.min_routing_layer;
       layer_index < num_layers_;
       layer_index++) {
    const int direction = getLayerDirection(layer_index);
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        if (checkOverflow(layer_index, x, y)) {
          view[direction][x][y] = true;
        }
      }
    }
  }
}

void GridGraph::extractWireCostView(GridGraphView<CostT>& view) const
{
  view.assign(
      2,
      std::vector<std::vector<CostT>>(
          x_size_,
          std::vector<CostT>(y_size_, std::numeric_limits<CostT>::max())));
  for (int direction = 0; direction < 2; direction++) {
    std::vector<int> layerIndices;
    CostT unitLengthShortCost = std::numeric_limits<CostT>::max();
    for (int layer_index = constants_.min_routing_layer;
         layer_index < getNumLayers();
         layer_index++) {
      if (getLayerDirection(layer_index) == direction) {
        layerIndices.emplace_back(layer_index);
        unitLengthShortCost = std::min(unitLengthShortCost,
                                       getUnitLengthShortCost(layer_index));
      }
    }
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        const int edge_index = direction == MetalLayer::H ? x : y;
        if (edge_index >= getSize(direction) - 1) {
          continue;
        }
        CapacityT capacity = 0;
        CapacityT demand = 0;
        for (int layer_index : layerIndices) {
          const auto& edge = getEdge(layer_index, x, y);
          capacity += getSoftCapacity(edge);
          demand += edge.demand;
        }
        const int length = getEdgeLength(direction, edge_index);
        view[direction][x][y]
            = length
              * (unit_length_wire_cost_
                 + unitLengthShortCost
                       * getCongestionPenalty(capacity,
                                              demand,
                                              constants_.maze_logistic_slope));
      }
    }
  }
}

void GridGraph::extractWireLengthCostView(GridGraphView<CostT>& view) const
{
  view.assign(
      2,
      std::vector<std::vector<CostT>>(
          x_size_,
          std::vector<CostT>(y_size_, std::numeric_limits<CostT>::max())));
  for (int direction = 0; direction < 2; direction++) {
    for (int x = 0; x < x_size_; x++) {
      for (int y = 0; y < y_size_; y++) {
        const int edge_index = direction == MetalLayer::H ? x : y;
        if (edge_index >= getSize(direction) - 1) {
          continue;
        }
        const int length = getEdgeLength(direction, edge_index);
        view[direction][x][y] = length * unit_length_wire_cost_;
      }
    }
  }
}

void GridGraph::updateWireCostView(
    GridGraphView<CostT>& view,
    const std::shared_ptr<GRTreeNode>& routing_tree) const
{
  std::vector<std::vector<int>> sameDirectionLayers(2);
  std::vector<CostT> unitLengthShortCost(2, std::numeric_limits<CostT>::max());
  for (int layer_index = constants_.min_routing_layer;
       layer_index < getNumLayers();
       layer_index++) {
    const int direction = getLayerDirection(layer_index);
    sameDirectionLayers[direction].emplace_back(layer_index);
    unitLengthShortCost[direction] = std::min(
        unitLengthShortCost[direction], getUnitLengthShortCost(layer_index));
  }
  auto update = [&](int direction, int x, int y) {
    const int edge_index = direction == MetalLayer::H ? x : y;
    if (edge_index >= getSize(direction) - 1) {
      return;
    }
    CapacityT capacity = 0;
    CapacityT demand = 0;
    for (int layer_index : sameDirectionLayers[direction]) {
      if (getLayerDirection(layer_index) != direction) {
        continue;
      }
      const auto& edge = getEdge(layer_index, x, y);
      capacity += getSoftCapacity(edge);
      demand += edge.demand;
    }
    const int length = getEdgeLength(direction, edge_index);
    view[direction][x][y]
        = length
          * (unit_length_wire_cost_
             + unitLengthShortCost[direction]
                   * getCongestionPenalty(capacity,
                                          demand,
                                          constants_.maze_logistic_slope));
  };
  GRTreeNode::preorder(
      routing_tree, [&](const std::shared_ptr<GRTreeNode>& node) {
        for (const auto& child : node->getChildren()) {
          if (node->getLayerIdx() == child->getLayerIdx()) {
            const int direction = getLayerDirection(node->getLayerIdx());
            if (direction == MetalLayer::H) {
              assert(node->y() == child->y());
              const int l = std::min(node->x(), child->x()),
                        h = std::max(node->x(), child->x());
              for (int x = l; x < h; x++) {
                update(direction, x, node->y());
              }
            } else {
              assert(node->x() == child->x());
              const int l = std::min(node->y(), child->y()),
                        h = std::max(node->y(), child->y());
              for (int y = l; y < h; y++) {
                update(direction, node->x(), y);
              }
            }
          } else {
            const int maxLayerIndex
                = std::max(node->getLayerIdx(), child->getLayerIdx());
            for (int layerIdx
                 = std::min(node->getLayerIdx(), child->getLayerIdx());
                 layerIdx < maxLayerIndex;
                 layerIdx++) {
              const int direction = getLayerDirection(layerIdx);
              update(direction, node->x(), node->y());
              if ((*node)[direction] > 0) {
                update(direction,
                       node->x() - 1 + direction,
                       node->y() - direction);
              }
            }
          }
        }
      });
}

void GridGraph::write(const std::string& heatmap_file) const
{
  logger_->report("writing heatmap to file...");
  std::stringstream ss;

  ss << num_layers_ << " " << x_size_ << " " << y_size_ << '\n';
  for (int layer_index = 0; layer_index < num_layers_; layer_index++) {
    ss << layer_names_[layer_index] << '\n';
    for (int y = 0; y < y_size_; y++) {
      for (int x = 0; x < x_size_; x++) {
        ss << (graph_edges_[layer_index][x][y].capacity
               - graph_edges_[layer_index][x][y].demand)
           << (x == x_size_ - 1 ? "" : " ");
      }
      ss << '\n';
    }
  }
  std::ofstream fout(heatmap_file);
  fout << ss.str();
  fout.close();
}

}  // namespace grt::newgr
