#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FastRoute.h"
#include "Grid.h"
#include "grt/Rudy.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
  int congestion_iterations = 0;
};

using RudyGrid = std::vector<std::vector<float>>;

RudyGrid computeNormalizedRudyGrid(Rudy* rudy)
{
  RudyGrid normalized;
  if (rudy == nullptr) {
    return normalized;
  }
  const auto [x_tiles, y_tiles] = rudy->getGridSize();
  if (x_tiles == 0 || y_tiles == 0) {
    return normalized;
  }

  normalized.resize(x_tiles, std::vector<float>(y_tiles, 0.0f));
  float max_value = 0.0f;
  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      const float val = rudy->getTile(x, y).getRudy();
      normalized[x][y] = val;
      max_value = std::max(max_value, val);
    }
  }

  if (max_value <= std::numeric_limits<float>::epsilon()) {
    return normalized;
  }

  for (int x = 0; x < x_tiles; ++x) {
    for (int y = 0; y < y_tiles; ++y) {
      normalized[x][y] = std::clamp(normalized[x][y] / max_value, 0.0f, 1.0f);
    }
  }

  return normalized;
}

struct HotspotTile
{
  int gx = 0;
  int gy = 0;
  float score = 0.0f;
};

std::vector<HotspotTile> collectRudyHotspots(const RudyGrid& normalized_rudy,
                                             Grid* grid,
                                             float threshold,
                                             int max_hotspots)
{
  std::vector<HotspotTile> hotspots;
  if (normalized_rudy.empty() || grid == nullptr) {
    return hotspots;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return hotspots;
  }

  const int usable_x = std::min<int>(x_grids, normalized_rudy.size());
  const int usable_y = std::min<int>(y_grids, normalized_rudy.front().size());
  threshold = std::clamp(threshold, 0.0f, 1.0f);
  max_hotspots = std::max(0, max_hotspots);

  for (int x = 0; x < usable_x; ++x) {
    for (int y = 0; y < usable_y; ++y) {
      const float val = normalized_rudy[x][y];
      if (val >= threshold) {
        hotspots.push_back({x, y, val});
      }
    }
  }

  std::sort(hotspots.begin(),
            hotspots.end(),
            [](const HotspotTile& a, const HotspotTile& b) {
              return a.score > b.score;
            });

  if (max_hotspots > 0 && static_cast<int>(hotspots.size()) > max_hotspots) {
    hotspots.resize(max_hotspots);
  }
  return hotspots;
}

std::vector<std::vector<uint8_t>> buildHotspotMask(
    const std::vector<HotspotTile>& hotspots,
    Grid* grid)
{
  std::vector<std::vector<uint8_t>> mask;
  if (grid == nullptr) {
    return mask;
  }
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 0 || y_grids <= 0) {
    return mask;
  }
  mask.assign(x_grids, std::vector<uint8_t>(y_grids, 0));
  for (const HotspotTile& h : hotspots) {
    if (h.gx >= 0 && h.gy >= 0 && h.gx < x_grids && h.gy < y_grids) {
      mask[h.gx][h.gy] = 1;
    }
  }
  return mask;
}

odb::Point tileCenterDbu(Grid* grid, int gx, int gy)
{
  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int cx = tile_size * (gx + 0.5) + x_min;
  const int cy = tile_size * (gy + 0.5) + y_min;
  return odb::Point(cx, cy);
}

std::pair<int, int> dbuToTileIndex(Grid* grid, int x_dbu, int y_dbu)
{
  const int tile_size = std::max(grid->getTileSize(), 1);
  const int x_min = grid->getXMin();
  const int y_min = grid->getYMin();
  const int gx = (x_dbu - x_min) / tile_size;
  const int gy = (y_dbu - y_min) / tile_size;
  return {gx, gy};
}

bool addUniqueSegment(GRoute& route,
                      std::unordered_set<GSegment, GSegmentHash>& seen,
                      const GSegment& segment)
{
  if (seen.insert(segment).second) {
    route.push_back(segment);
    return true;
  }
  return false;
}

void addCrossPatchGuides(Grid* grid,
                         GRoute& route,
                         std::unordered_set<GSegment, GSegmentHash>& seen,
                         int gx,
                         int gy,
                         int layer,
                         int min_layer,
                         int max_layer)
{
  if (grid == nullptr) {
    return;
  }
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
    return;
  }

  layer = std::clamp(layer, min_layer, max_layer);

  const odb::Point center = tileCenterDbu(grid, gx, gy);
  const int cx = center.getX();
  const int cy = center.getY();
  const int tile_size = std::max(grid->getTileSize(), 1);

  auto add_wire_to_neighbor = [&](int ngx, int ngy) {
    if (ngx < 0 || ngy < 0 || ngx >= x_grids || ngy >= y_grids) {
      return;
    }
    const odb::Point neigh = tileCenterDbu(grid, ngx, ngy);
    addUniqueSegment(route,
                     seen,
                     GSegment(cx, cy, layer, neigh.getX(), neigh.getY(), layer));
  };

  add_wire_to_neighbor(gx - 1, gy);
  add_wire_to_neighbor(gx + 1, gy);
  add_wire_to_neighbor(gx, gy - 1);
  add_wire_to_neighbor(gx, gy + 1);

  const int up_layer = layer + 1;
  const int down_layer = layer - 1;
  if (up_layer <= max_layer) {
    addUniqueSegment(route, seen, GSegment(cx, cy, layer, cx, cy, up_layer));
    // Also add a short segment on the adjacent layer to make the via usable.
    if (gx + 1 < x_grids) {
      addUniqueSegment(route,
                       seen,
                       GSegment(cx,
                                cy,
                                up_layer,
                                cx + tile_size,
                                cy,
                                up_layer));
    }
  } else if (down_layer >= min_layer) {
    addUniqueSegment(route, seen, GSegment(cx, cy, layer, cx, cy, down_layer));
    if (gx + 1 < x_grids) {
      addUniqueSegment(route,
                       seen,
                       GSegment(cx,
                                cy,
                                down_layer,
                                cx + tile_size,
                                cy,
                                down_layer));
    }
  }
}

}  // namespace

NewGR::NewGR(GlobalRouter* grouter, CUGR* cugr, utl::Logger* logger)
    : grouter_(grouter), cugr_(cugr), logger_(logger)
{
}

NetRouteMap NewGR::run(std::vector<Net*>& nets,
                       int min_routing_layer,
                       int max_routing_layer)
{
  if (nets.empty()) {
    return {};
  }

  auto capture_snapshot = [&]() -> RouterSnapshot {
    RouterSnapshot snapshot;
    snapshot.caps_percentage = grouter_->caps_perturbation_percentage_;
    snapshot.perturbation_amount = grouter_->perturbation_amount_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    snapshot.congestion_iterations = grouter_->congestion_iterations_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
    grouter_->setCongestionIterations(snapshot.congestion_iterations);
  };

  RouterSnapshot snapshot = capture_snapshot();
  // NEWGR objective for this branch: provide the detailed router with more
  // freedom near predicted hotspots without paying for multiple full GR runs.
  //
  // We keep the routing engine (FastRoute) unchanged, but post-process the
  // guides by adding small "patch" guides in high-RUDY regions that the net
  // already traverses. This aims to reduce detailed-routing detours, improving
  // final wirelength while keeping runtime close to FastRoute.
  grouter_->setCapacitiesPerturbationPercentage(0.0f);
  grouter_->setPerturbationAmount(1);
  grouter_->fastroute_->setCriticalNetsPercentage(
      snapshot.critical_percentage);
  grouter_->setCongestionIterations(snapshot.congestion_iterations);

  Rudy* rudy = grouter_->getRudy();
  if (rudy != nullptr) {
    rudy->calculateRudy();
  }
  const RudyGrid normalized_rudy = computeNormalizedRudyGrid(rudy);
  const std::vector<HotspotTile> rudy_hotspots
      = collectRudyHotspots(normalized_rudy, grouter_->grid_, 0.85f, 350);
  const auto hotspot_mask = buildHotspotMask(rudy_hotspots, grouter_->grid_);

  NetRouteMap routes = grouter_->findRouting(
      nets, min_routing_layer, max_routing_layer);

  if (grouter_->grid_ != nullptr && !hotspot_mask.empty()) {
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();
    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);

    int patched_nets = 0;
    int patched_points = 0;

    for (auto& [db_net, route] : routes) {
      if (route.empty() || db_net == nullptr) {
        continue;
      }

      std::unordered_set<GSegment, GSegmentHash> seen(route.begin(),
                                                      route.end());
      int patch_budget = 2;
      bool patched_this_net = false;

      for (const GSegment& seg : route) {
        if (patch_budget == 0) {
          break;
        }
        if (seg.isVia()) {
          continue;
        }
        if (seg.init_layer != seg.final_layer) {
          continue;
        }

        auto [gx1, gy1]
            = dbuToTileIndex(grouter_->grid_, seg.init_x, seg.init_y);
        auto [gx2, gy2]
            = dbuToTileIndex(grouter_->grid_, seg.final_x, seg.final_y);
        gx1 = std::clamp(gx1, 0, std::max(x_grids - 1, 0));
        gx2 = std::clamp(gx2, 0, std::max(x_grids - 1, 0));
        gy1 = std::clamp(gy1, 0, std::max(y_grids - 1, 0));
        gy2 = std::clamp(gy2, 0, std::max(y_grids - 1, 0));

        const int layer = seg.init_layer;
        if (gx1 == gx2 && gy1 == gy2) {
          continue;
        }

        if (gy1 == gy2) {
          const int y = gy1;
          const int start = std::min(gx1, gx2);
          const int end = std::max(gx1, gx2);
          for (int x = start; x <= end; ++x) {
            if (hotspot_mask[x][y]) {
              addCrossPatchGuides(grouter_->grid_,
                                  route,
                                  seen,
                                  x,
                                  y,
                                  layer,
                                  min_routing_layer,
                                  max_routing_layer);
              patch_budget--;
              patched_points++;
              patched_this_net = true;
              break;
            }
          }
        } else if (gx1 == gx2) {
          const int x = gx1;
          const int start = std::min(gy1, gy2);
          const int end = std::max(gy1, gy2);
          for (int y = start; y <= end; ++y) {
            if (hotspot_mask[x][y]) {
              addCrossPatchGuides(grouter_->grid_,
                                  route,
                                  seen,
                                  x,
                                  y,
                                  layer,
                                  min_routing_layer,
                                  max_routing_layer);
              patch_budget--;
              patched_points++;
              patched_this_net = true;
              break;
            }
          }
        } else {
          // Should not happen for rectilinear segments, but be robust.
          const int steps = std::max(std::abs(gx1 - gx2), std::abs(gy1 - gy2));
          for (int i = 0; i <= steps; ++i) {
            const int gx
                = gx1 + static_cast<int>(std::round((gx2 - gx1) * (i / static_cast<double>(steps))));
            const int gy
                = gy1 + static_cast<int>(std::round((gy2 - gy1) * (i / static_cast<double>(steps))));
            if (gx >= 0 && gy >= 0 && gx < x_grids && gy < y_grids
                && hotspot_mask[gx][gy]) {
              addCrossPatchGuides(grouter_->grid_,
                                  route,
                                  seen,
                                  gx,
                                  gy,
                                  layer,
                                  min_routing_layer,
                                  max_routing_layer);
              patch_budget--;
              patched_points++;
              patched_this_net = true;
              break;
            }
          }
        }
      }

      if (patched_this_net) {
        patched_nets++;
      }
    }

    if (patched_points > 0) {
      logger_->info(GNR,
                    6010,
                    "NEWGR: added {} patch points across {} nets using {} RUDY hotspots",
                    patched_points,
                    patched_nets,
                    rudy_hotspots.size());
    }
    static_cast<void>(tile_size);
  }

  restore_snapshot(snapshot);
  return routes;
}

}  // namespace grt
