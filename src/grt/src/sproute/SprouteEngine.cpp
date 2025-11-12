#include "SprouteEngine.h"

#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <utility>

#include "global.h"
#include "fastroute.h"
#include "grGen.h"
#include "utl/Logger.h"

#ifdef HORIZONTAL
#undef HORIZONTAL
#endif
#ifdef VERTICAL
#undef VERTICAL
#endif

namespace grt {

namespace {

const Edge3D dummy_edge{};

const Edge3D& horizontalEdge(const SprouteGridData& grid,
                             int layer,
                             int y,
                             int x)
{
  if (grid.x_grids <= 1) {
    return dummy_edge;
  }
  const int per_layer = (grid.x_grids - 1) * grid.y_grids;
  const int idx = layer * per_layer + y * (grid.x_grids - 1) + x;
  return h_edges3D[idx];
}

const Edge3D& verticalEdge(const SprouteGridData& grid,
                           int layer,
                           int y,
                           int x)
{
  if (grid.y_grids <= 1) {
    return dummy_edge;
  }
  const int per_layer = grid.x_grids * (grid.y_grids - 1);
  const int idx = layer * per_layer + y * grid.x_grids + x;
  return v_edges3D[idx];
}

}  // namespace

SprouteEngine::SprouteEngine(utl::Logger* logger) : logger_(logger)
{
}

void SprouteEngine::init(const SprouteGridData& grid,
                         const std::vector<SprouteNetData>& nets)
{
  grid_ = grid;
  nets_ = nets;
  buildInput();
}

NetRouteMap SprouteEngine::run()
{
  if (!input_ready_) {
    buildInput();
  }

  prepareLefDefMetadata();
  parser::grGenerator generator = buildGenerator();

  galois::SharedMemSys galois_runtime;
  if (numThreads <= 0) {
    numThreads = 1;
  }
  galois::preAlloc(numThreads * 2);
  numThreads = galois::setActiveThreads(numThreads);

  parser::CongestionMap congestion_map(
      generator.grid.z, generator.grid.x, generator.grid.y);
  galois::StatTimer timer("sproute_timer");
  timer.start();
  runFastRoute(generator,
               /*benchFile=*/"",
               /*OutFileName=*/"",
               congestion_map,
               timer,
               /*maxMazeRound=*/500,
               Algo::FineGrain);
  timer.stop();
  last_total_overflow_ = totalOverflow;

  return extractRoutes();
}

void SprouteEngine::buildInput()
{
  input_ = SprouteInput();
  input_.grid = grid_;
  input_.nets.reserve(nets_.size());

  for (const auto& net : nets_) {
    SprouteInputNet input_net;
    input_net.db_net = net.db_net;
    if (net.db_net != nullptr) {
      input_net.name = net.db_net->getConstName();
    }
    input_net.is_clock = net.is_clock;
    input_net.min_layer = net.min_layer;
    input_net.max_layer = net.max_layer;
    input_net.root_pin_index = net.root_pin_index;
    input_net.pins = net.pins;
    input_.nets.push_back(std::move(input_net));
  }

  input_ready_ = true;
}

void SprouteEngine::prepareLefDefMetadata()
{
  defDB.xGcellBoundaries.clear();
  defDB.yGcellBoundaries.clear();
  defDB.layerName2trackidx.clear();
  defDB.tracks.clear();
  defDB.netName2netidx.clear();
  defDB.nets.clear();

  const int origin_x = grid_.origin.x();
  const int origin_y = grid_.origin.y();
  const int tile = grid_.tile_size;

  for (int x = 0; x <= grid_.x_grids; ++x) {
    defDB.xGcellBoundaries.push_back(origin_x + x * tile);
  }

  for (int y = 0; y <= grid_.y_grids; ++y) {
    defDB.yGcellBoundaries.push_back(origin_y + y * tile);
  }

  defDB.gcellGridDim.x = grid_.x_grids;
  defDB.gcellGridDim.y = grid_.y_grids;
  defDB.gcellGridDim.z = grid_.num_layers;
  defDB.dbuPerMicro = grid_.dbu_per_micron;

  lefDB.layers.clear();
  lefDB.layer2idx.clear();

  for (int layer = 0; layer < grid_.num_layers; ++layer) {
    parser::Layer routing_layer;
    routing_layer.reset();
    routing_layer.name = routingLayerName(layer);
    routing_layer.type = "ROUTING";
    routing_layer.direction
        = (grid_.layer_directions[layer] == odb::dbTechLayerDir::Value::VERTICAL)
              ? "VERTICAL"
              : "HORIZONTAL";
    routing_layer.idx = static_cast<int>(lefDB.layers.size());
    lefDB.layer2idx[routing_layer.name] = routing_layer.idx;
    lefDB.layers.push_back(routing_layer);

    parser::Layer cut_layer;
    cut_layer.reset();
    cut_layer.name = cutLayerName(layer);
    cut_layer.type = "CUT";
    cut_layer.idx = static_cast<int>(lefDB.layers.size());
    lefDB.layer2idx[cut_layer.name] = cut_layer.idx;
    lefDB.layers.push_back(cut_layer);

    parser::Track track;
    if (routing_layer.direction == "HORIZONTAL") {
      track.direction = "Y";
      track.start = origin_y;
      track.numTracks = grid_.y_grids;
      track.step = tile;
    } else {
      track.direction = "X";
      track.start = origin_x;
      track.numTracks = grid_.x_grids;
      track.step = tile;
    }
    track.layerNames.push_back(routing_layer.name);
    int track_index = static_cast<int>(defDB.tracks.size());
    defDB.tracks.push_back(track);
    defDB.layerName2trackidx[routing_layer.name] = track_index;
  }
}

parser::grGenerator SprouteEngine::buildGenerator() const
{
  parser::grGenerator generator;
  generator.grid.x = grid_.x_grids;
  generator.grid.y = grid_.y_grids;
  generator.grid.z = grid_.num_layers;

  generator.gcellBegin.set(grid_.origin.x(), grid_.origin.y());
  generator.gcellSize.set(grid_.tile_size, grid_.tile_size);

  generator.vCap = grid_.v_capacities;
  generator.hCap = grid_.h_capacities;

  generator.numNets = static_cast<int>(input_.nets.size());
  generator.grnets.reserve(generator.numNets);

  int net_idx = 0;
  for (const auto& net : input_.nets) {
    parser::grNet generator_net;
    generator_net.name = net.name;
    generator_net.idx = net_idx;
    generator_net.numPins = static_cast<int>(net.pins.size());
    generator_net.pins.reserve(net.pins.size());

    for (const RoutePt& pin : net.pins) {
      parser::Point3D<int> point;
      point.x = coordFromIndex(pin.x());
      point.y = coordFromIndex(pin.y(), /*is_y=*/true);
      point.z = pin.layer();
      generator_net.pins.push_back(point);
    }

    generator.grnets.push_back(generator_net);
    generator.netName2grnetidx[generator_net.name] = net_idx;
    if (static_cast<int>(defDB.nets.size()) <= net_idx) {
      defDB.nets.resize(net_idx + 1);
    }
    defDB.nets[net_idx].reset();
    defDB.nets[net_idx].name = generator_net.name;
    defDB.netName2netidx[generator_net.name] = net_idx;
    ++net_idx;
  }

  return generator;
}

int SprouteEngine::coordFromIndex(int index, bool is_y) const
{
  const int origin = is_y ? grid_.origin.y() : grid_.origin.x();
  return origin + index * grid_.tile_size + grid_.tile_size / 2;
}

std::string SprouteEngine::routingLayerName(int layer_index) const
{
  return "metal" + std::to_string(layer_index + 1);
}

std::string SprouteEngine::cutLayerName(int layer_index) const
{
  return "via" + std::to_string(layer_index + 1);
}

NetRouteMap SprouteEngine::extractRoutes() const
{
  NetRouteMap routes;
  if (sttrees == nullptr || grid_.tile_size == 0) {
    return routes;
  }

  std::unordered_map<std::string, odb::dbNet*> net_lookup;
  net_lookup.reserve(input_.nets.size());
  for (const auto& net : input_.nets) {
    if (net.db_net != nullptr) {
      net_lookup.emplace(net.name, net.db_net);
    }
  }

  for (int net_id = 0; net_id < numValidNets; ++net_id) {
    const std::string net_name = nets[net_id]->name;
    auto it = net_lookup.find(net_name);
    if (it == net_lookup.end()) {
      continue;
    }
    GRoute route_segments;
    appendRouteSegments(net_id, route_segments);
    if (!route_segments.empty()) {
      routes[it->second] = std::move(route_segments);
    }
  }

  return routes;
}

void SprouteEngine::appendRouteSegments(int net_id, GRoute& route) const
{
  const StTree& tree = sttrees[net_id];
  if (tree.deg == 0 || tree.edges == nullptr) {
    return;
  }

  const int edge_count = 2 * tree.deg - 3;
  for (int edge_id = 0; edge_id < edge_count; ++edge_id) {
    const TreeEdge& tree_edge = tree.edges[edge_id];
    const Route& edge_route = tree_edge.route;
    if (edge_route.routelen <= 0 || edge_route.gridsX == nullptr
        || edge_route.gridsY == nullptr || edge_route.gridsL == nullptr) {
      continue;
    }
    for (int i = 0; i < edge_route.routelen; ++i) {
      const int x0 = edge_route.gridsX[i];
      const int y0 = edge_route.gridsY[i];
      const int l0 = edge_route.gridsL[i];
      const int x1 = edge_route.gridsX[i + 1];
      const int y1 = edge_route.gridsY[i + 1];
      const int l1 = edge_route.gridsL[i + 1];
      addSegment(route, x0, y0, l0, x1, y1, l1);
    }
  }
}

void SprouteEngine::addSegment(GRoute& route,
                               int grid_x0,
                               int grid_y0,
                               int grid_l0,
                               int grid_x1,
                               int grid_y1,
                               int grid_l1) const
{
  const int x0 = coordFromIndex(grid_x0, /*is_y=*/false);
  const int y0 = coordFromIndex(grid_y0, /*is_y=*/true);
  const int x1 = coordFromIndex(grid_x1, /*is_y=*/false);
  const int y1 = coordFromIndex(grid_y1, /*is_y=*/true);
  const int l0 = toDbLayer(grid_l0);
  const int l1 = toDbLayer(grid_l1);

  route.emplace_back(x0, y0, l0, x1, y1, l1);
}

int SprouteEngine::toDbLayer(int sproute_layer) const
{
  int layer = sproute_layer + 1;
  if (layer < 1) {
    layer = 1;
  }
  if (layer > grid_.num_layers) {
    layer = grid_.num_layers;
  }
  return layer;
}

void SprouteEngine::updateDbCongestion(odb::dbBlock* block)
{
  if (block == nullptr || grid_.x_grids == 0 || grid_.y_grids == 0) {
    return;
  }
  odb::dbGCellGrid* db_gcell = block->getGCellGrid();
  if (db_gcell) {
    db_gcell->resetGrid();
  } else {
    db_gcell = odb::dbGCellGrid::create(block);
  }

  db_gcell->addGridPatternX(
      grid_.origin.x(), grid_.x_grids, grid_.tile_size);
  db_gcell->addGridPatternY(
      grid_.origin.y(), grid_.y_grids, grid_.tile_size);

  odb::dbTech* tech = block->getDb()->getTech();
  const int x_grids = grid_.x_grids;
  const int y_grids = grid_.y_grids;

  for (int layer = 0; layer < grid_.num_layers; ++layer) {
    odb::dbTechLayer* tech_layer = tech->findRoutingLayer(layer + 1);
    if (tech_layer == nullptr) {
      continue;
    }

    const bool is_horizontal
        = layer < static_cast<int>(grid_.layer_directions.size())
          && grid_.layer_directions[layer]
                 == odb::dbTechLayerDir::Value::HORIZONTAL;
    const uint16_t capH
        = layer < static_cast<int>(grid_.h_capacities.size())
              ? grid_.h_capacities[layer]
              : 0;
    const uint16_t capV
        = layer < static_cast<int>(grid_.v_capacities.size())
              ? grid_.v_capacities[layer]
              : 0;

    for (int y = 0; y < y_grids; ++y) {
      for (int x = 0; x < x_grids; ++x) {
        uint16_t usageH = 0;
        uint16_t usageV = 0;

        if (x_grids > 1 && capH > 0) {
          const int hx = std::max(
              0, std::min(x, x_grids - 2));
          const Edge3D& hedge = horizontalEdge(grid_, layer, y, hx);
          const uint16_t blockageH
              = capH > hedge.cap ? capH - hedge.cap : 0;
          usageH = hedge.usage + blockageH;
        }

        if (y_grids > 1 && capV > 0) {
          const int vy = std::max(
              0, std::min(y, y_grids - 2));
          const Edge3D& vedge = verticalEdge(grid_, layer, vy, x);
          const uint16_t blockageV
              = capV > vedge.cap ? capV - vedge.cap : 0;
          usageV = vedge.usage + blockageV;
        }

        const uint16_t capacity = is_horizontal ? capH : capV;
        db_gcell->setCapacity(tech_layer, x, y, capacity);
        db_gcell->setUsage(tech_layer, x, y, usageH + usageV);
      }
    }
  }
}

}  // namespace grt
