#include "NewgrEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
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

void ensureGaloisRuntime()
{
  static std::unique_ptr<galois::SharedMemSys> runtime
      = std::make_unique<galois::SharedMemSys>();
  (void) runtime;
}

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

galois::InsertBag<parser::CapReduction>& newgrCapReductions()
{
  static galois::InsertBag<parser::CapReduction> cap_reductions;
  return cap_reductions;
}

int cellIndex(int x, int y, int x_grids)
{
  return y * x_grids + x;
}

struct RouteCost
{
  int64_t wirelength{0};
  int64_t vias{0};
};

RouteCost estimateRouteCost(const GRoute& route)
{
  RouteCost cost;
  for (const GSegment& segment : route) {
    cost.wirelength += std::llabs(static_cast<long long>(segment.init_x)
                                  - static_cast<long long>(segment.final_x));
    cost.wirelength += std::llabs(static_cast<long long>(segment.init_y)
                                  - static_cast<long long>(segment.final_y));
    cost.vias += std::llabs(static_cast<long long>(segment.init_layer)
                            - static_cast<long long>(segment.final_layer));
  }
  return cost;
}

double combinedRouteScore(const RouteCost& cost, int tile_size, int pin_count)
{
  const double tile_scale
      = std::sqrt(static_cast<double>(std::max(1, tile_size)));
  const double via_weight_scale = (pin_count >= 24) ? 0.95
                                  : (pin_count >= 12) ? 1.35
                                                      : 1.75;
  const double via_weight = tile_scale * via_weight_scale;
  return static_cast<double>(cost.wirelength)
         + static_cast<double>(cost.vias) * via_weight;
}

float reductionRatio(int base_cap,
                     float local_pressure,
                     float mean_pressure,
                     float stdev_pressure,
                     float distance_to_center,
                     float normalized_layer,
                     bool horizontal)
{
  constexpr float kEpsilon = 1e-3f;
  const float pressure_sigma_raw
      = (local_pressure - mean_pressure) / std::max(stdev_pressure, kEpsilon);
  const float pressure_sigma = std::clamp(pressure_sigma_raw, -3.0f, 4.0f);
  const float pressure_term
      = 1.0f / (1.0f + std::exp(-(pressure_sigma + 0.15f) * 2.2f));
  const float center_term
      = std::clamp(1.0f - distance_to_center * 1.35f, 0.0f, 1.0f);

  // Radical experiment: create a strong "highway core" by preserving center
  // capacity and sharply reducing off-core resources.
  float ratio = 0.22f;
  ratio += 0.58f * center_term;
  ratio += 0.20f * pressure_term;
  ratio += 0.10f * normalized_layer;
  if (horizontal) {
    ratio += 0.04f;
  }

  if (local_pressure < mean_pressure - 0.40f * std::max(stdev_pressure, kEpsilon)) {
    ratio -= 0.10f;
  }
  if (base_cap <= 2) {
    ratio += 0.08f;
  }
  return std::clamp(ratio, 0.10f, 1.00f);
}

int buildLocalizedCapacityReductions(const NewgrInput& input,
                                     const SprouteGridData& grid,
                                     const parser::grGenerator& generator,
                                     galois::InsertBag<parser::CapReduction>&
                                         cap_reductions)
{
  const int x_grids = grid.x_grids;
  const int y_grids = grid.y_grids;
  const int num_layers = grid.num_layers;
  if (x_grids <= 1 || y_grids <= 1 || num_layers <= 0) {
    return 0;
  }
  const int total_edges = num_layers
                          * ((x_grids - 1) * y_grids + x_grids * (y_grids - 1));
  // Radical experiment: aggressively reshape most of the resource graph.
  const int max_reductions = std::max(4096, (total_edges * 3) / 4);

  std::vector<float> pin_pressure(x_grids * y_grids, 0.0f);
  for (const auto& net : input.nets) {
    const float net_weight
        = 1.0f + std::min(6.0f, std::sqrt(static_cast<float>(net.pins.size())));
    for (const auto& pin : net.pins) {
      const int x = std::clamp(pin.x(), 0, x_grids - 1);
      const int y = std::clamp(pin.y(), 0, y_grids - 1);
      pin_pressure[cellIndex(x, y, x_grids)] += net_weight;
    }
  }

  std::vector<float> smooth_pressure(x_grids * y_grids, 0.0f);
  float mean_pressure = 0.0f;
  for (int y = 0; y < y_grids; ++y) {
    for (int x = 0; x < x_grids; ++x) {
      float sum = 0.0f;
      int count = 0;
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          const int nx = x + dx;
          const int ny = y + dy;
          if (nx < 0 || nx >= x_grids || ny < 0 || ny >= y_grids) {
            continue;
          }
          sum += pin_pressure[cellIndex(nx, ny, x_grids)];
          ++count;
        }
      }
      const float pressure = (count > 0) ? (sum / count) : 0.0f;
      smooth_pressure[cellIndex(x, y, x_grids)] = pressure;
      mean_pressure += pressure;
    }
  }
  mean_pressure /= static_cast<float>(x_grids * y_grids);

  float variance = 0.0f;
  for (float pressure : smooth_pressure) {
    const float delta = pressure - mean_pressure;
    variance += delta * delta;
  }
  variance /= static_cast<float>(x_grids * y_grids);
  const float stdev_pressure = std::sqrt(std::max(variance, 0.0f));

  const float center_x = (x_grids - 1) * 0.5f;
  const float center_y = (y_grids - 1) * 0.5f;
  const float max_center_dist
      = std::max(1.0f, std::sqrt(center_x * center_x + center_y * center_y));

  int reduction_count = 0;
  for (int layer = 0; layer < num_layers; ++layer) {
    const bool horizontal
        = layer < static_cast<int>(grid.layer_directions.size())
          && grid.layer_directions[layer]
                 == odb::dbTechLayerDir::Value::HORIZONTAL;
    const int base_hcap
        = layer < static_cast<int>(generator.hCap.size()) ? generator.hCap[layer]
                                                           : 0;
    const int base_vcap
        = layer < static_cast<int>(generator.vCap.size()) ? generator.vCap[layer]
                                                           : 0;

    if (horizontal) {
      if (base_hcap <= 1) {
        continue;
      }
      for (int y = 0; y < y_grids; ++y) {
        for (int x = 0; x < x_grids - 1; ++x) {
          const float p0 = smooth_pressure[cellIndex(x, y, x_grids)];
          const float p1 = smooth_pressure[cellIndex(x + 1, y, x_grids)];
          const float local_pressure = 0.5f * (p0 + p1);
          const float ratio = reductionRatio(base_hcap,
                                             local_pressure,
                                             mean_pressure,
                                             stdev_pressure,
                                             std::sqrt((x + 0.5f - center_x)
                                                           * (x + 0.5f - center_x)
                                                       + (y - center_y)
                                                             * (y - center_y))
                                                 / max_center_dist,
                                             static_cast<float>(layer)
                                                 / std::max(1, num_layers - 1),
                                             /*horizontal=*/true);
          if (reduction_count >= max_reductions) {
            return reduction_count;
          }
          const int new_cap
              = std::max(1, static_cast<int>(std::lround(base_hcap * ratio)));
          if (new_cap >= base_hcap) {
            continue;
          }
          parser::CapReduction reduction;
          reduction.x = x;
          reduction.y = y;
          reduction.z = layer;
          reduction.newCap = new_cap;
          cap_reductions.push(reduction);
          ++reduction_count;
        }
      }
    } else {
      if (base_vcap <= 1) {
        continue;
      }
      for (int y = 0; y < y_grids - 1; ++y) {
        for (int x = 0; x < x_grids; ++x) {
          const float p0 = smooth_pressure[cellIndex(x, y, x_grids)];
          const float p1 = smooth_pressure[cellIndex(x, y + 1, x_grids)];
          const float local_pressure = 0.5f * (p0 + p1);
          const float ratio = reductionRatio(base_vcap,
                                             local_pressure,
                                             mean_pressure,
                                             stdev_pressure,
                                             std::sqrt((x - center_x)
                                                           * (x - center_x)
                                                       + (y + 0.5f - center_y)
                                                             * (y + 0.5f - center_y))
                                                 / max_center_dist,
                                             static_cast<float>(layer)
                                                 / std::max(1, num_layers - 1),
                                             /*horizontal=*/false);
          if (reduction_count >= max_reductions) {
            return reduction_count;
          }
          const int new_cap
              = std::max(1, static_cast<int>(std::lround(base_vcap * ratio)));
          if (new_cap >= base_vcap) {
            continue;
          }
          parser::CapReduction reduction;
          reduction.x = x;
          reduction.y = y;
          reduction.z = layer;
          reduction.newCap = new_cap;
          cap_reductions.push(reduction);
          ++reduction_count;
        }
      }
    }
  }

  return reduction_count;
}

}  // namespace

NewgrEngine::NewgrEngine(utl::Logger* logger) : logger_(logger)
{
}

void NewgrEngine::init(const SprouteGridData& grid,
                         const std::vector<SprouteNetData>& nets)
{
  grid_ = grid;
  nets_ = nets;
  buildInput();
}

NetRouteMap NewgrEngine::run()
{
  if (!input_ready_) {
    buildInput();
  }

  prepareLefDefMetadata();
  ensureGaloisRuntime();
  if (numThreads <= 0) {
    numThreads = 1;
  }
  galois::preAlloc(numThreads * 2);
  numThreads = galois::setActiveThreads(numThreads);
  parser::grGenerator generator = buildGenerator();

  parser::CongestionMap congestion_map(
      generator.grid.z, generator.grid.x, generator.grid.y);
  galois::StatTimer timer("newgr_timer");
  timer.start();
  if (generator.capReductions_p == nullptr) {
    logger_->warn(utl::GRT,
                  401,
                  "NEWGR generator has no localized capacity reductions; "
                  "continuing without adjustments.");
  }
  runFastRoute(generator,
               /*benchFile=*/"",
               /*OutFileName=*/"",
               congestion_map,
               timer,
               /*maxMazeRound=*/20,
               Algo::DetPart_Astar_RUDY);
  timer.stop();
  last_total_overflow_ = totalOverflow;

  NetRouteMap extracted_routes = extractRoutes();
  return extractPinFallbackRoutes(extracted_routes);
}

void NewgrEngine::buildInput()
{
  input_ = NewgrInput();
  input_.grid = grid_;
  input_.nets.reserve(nets_.size());

  for (const auto& net : nets_) {
    NewgrInputNet input_net;
    input_net.db_net = net.db_net;
    std::string base_name;
    if (net.db_net != nullptr) {
      base_name = net.db_net->getConstName();
    } else {
      base_name = "newgr_net_" + std::to_string(input_.nets.size());
    }
    input_net.name = sanitizeNetName(base_name);
    input_net.is_clock = net.is_clock;
    input_net.min_layer = net.min_layer;
    input_net.max_layer = net.max_layer;
    input_net.root_pin_index = net.root_pin_index;
    input_net.pins = net.pins;
    input_.nets.push_back(std::move(input_net));
  }

  input_ready_ = true;
}

void NewgrEngine::prepareLefDefMetadata()
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

parser::grGenerator NewgrEngine::buildGenerator() const
{
  parser::grGenerator generator;
  generator.grid.x = grid_.x_grids;
  generator.grid.y = grid_.y_grids;
  generator.grid.z = grid_.num_layers;

  generator.gcellBegin.set(grid_.origin.x(), grid_.origin.y());
  generator.gcellSize.set(grid_.tile_size, grid_.tile_size);

  generator.vCap = grid_.v_capacities;
  generator.hCap = grid_.h_capacities;
  generator.capReductions_p = nullptr;

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

  auto& cap_reductions = newgrCapReductions();
  cap_reductions.clear();
  const int reduction_count
      = buildLocalizedCapacityReductions(input_, grid_, generator, cap_reductions);
  if (reduction_count > 0) {
    generator.capReductions_p = &cap_reductions;
    logger_->report("NEWGR applied {} localized capacity reductions.",
                    reduction_count);
  } else {
    generator.capReductions_p = nullptr;
    logger_->warn(utl::GRT,
                  402,
                  "NEWGR generated zero localized capacity reductions.");
  }

  return generator;
}

int NewgrEngine::coordFromIndex(int index, bool is_y) const
{
  const int origin = is_y ? grid_.origin.y() : grid_.origin.x();
  return origin + index * grid_.tile_size + grid_.tile_size / 2;
}

std::string NewgrEngine::routingLayerName(int layer_index) const
{
  return "metal" + std::to_string(layer_index + 1);
}

std::string NewgrEngine::cutLayerName(int layer_index) const
{
  return "via" + std::to_string(layer_index + 1);
}

NetRouteMap NewgrEngine::extractRoutes() const
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
    if (appendRouteSegments(net_id, route_segments) && !route_segments.empty()) {
      routes[it->second] = std::move(route_segments);
    }
  }

  return routes;
}

NetRouteMap NewgrEngine::extractPinFallbackRoutes(
    const NetRouteMap& seeded_routes) const
{
  NetRouteMap routes = seeded_routes;
  int missing_net_fallback_count = 0;
  int replaced_net_count = 0;
  for (const auto& net : input_.nets) {
    if (net.db_net == nullptr || net.pins.size() < 2) {
      continue;
    }

    GRoute mst_route;
    const bool has_mst_route
        = appendFallbackMstRoute(net, mst_route) && !mst_route.empty();
    GRoute backbone_route;
    const bool has_backbone_route
        = appendFallbackBackboneRoute(net, backbone_route) && !backbone_route.empty();
    GRoute cross_route;
    const bool has_cross_route
        = appendFallbackCrossRoute(net, cross_route) && !cross_route.empty();
    GRoute dual_hub_route;
    const bool has_dual_hub_route
        = appendFallbackDualHubRoute(net, dual_hub_route)
          && !dual_hub_route.empty();
    if (!has_mst_route && !has_backbone_route && !has_cross_route
        && !has_dual_hub_route) {
      continue;
    }

    const int pin_count = static_cast<int>(net.pins.size());
    auto selectBestFallback = [&]() -> std::pair<const GRoute*, RouteCost> {
      const GRoute* selected = nullptr;
      RouteCost selected_cost;
      double selected_score = std::numeric_limits<double>::max();

      if (has_mst_route) {
        const RouteCost mst_cost = estimateRouteCost(mst_route);
        const double mst_score
            = combinedRouteScore(mst_cost, grid_.tile_size, pin_count);
        selected = &mst_route;
        selected_cost = mst_cost;
        selected_score = mst_score;
      }
      if (has_backbone_route) {
        const RouteCost backbone_cost = estimateRouteCost(backbone_route);
        const double backbone_score
            = combinedRouteScore(backbone_cost, grid_.tile_size, pin_count);
        if (selected == nullptr || backbone_score < selected_score) {
          selected = &backbone_route;
          selected_cost = backbone_cost;
          selected_score = backbone_score;
        }
      }
      if (has_cross_route) {
        const RouteCost cross_cost = estimateRouteCost(cross_route);
        const double cross_score
            = combinedRouteScore(cross_cost, grid_.tile_size, pin_count);
        if (selected == nullptr || cross_score < selected_score) {
          selected = &cross_route;
          selected_cost = cross_cost;
          selected_score = cross_score;
        }
      }
      if (has_dual_hub_route) {
        const RouteCost dual_hub_cost = estimateRouteCost(dual_hub_route);
        const double dual_hub_score
            = combinedRouteScore(dual_hub_cost, grid_.tile_size, pin_count);
        if (selected == nullptr || dual_hub_score < selected_score) {
          selected = &dual_hub_route;
          selected_cost = dual_hub_cost;
          selected_score = dual_hub_score;
        }
      }
      return {selected, selected_cost};
    };
    const auto [fallback_route, fallback_cost] = selectBestFallback();
    if (fallback_route == nullptr) {
      continue;
    }

    auto route_it = routes.find(net.db_net);
    if (route_it == routes.end()) {
      routes[net.db_net] = *fallback_route;
      ++missing_net_fallback_count;
      continue;
    }

    const RouteCost seeded_cost = estimateRouteCost(route_it->second);
    const double seeded_score
        = combinedRouteScore(seeded_cost, grid_.tile_size, pin_count);
    const double fallback_score
        = combinedRouteScore(fallback_cost, grid_.tile_size, pin_count);
    const bool strong_wire_improved
        = fallback_cost.wirelength
          < static_cast<double>(seeded_cost.wirelength) * 0.94;
    const bool balanced_improved
        = fallback_cost.wirelength
              < static_cast<double>(seeded_cost.wirelength) * 0.98
          && fallback_cost.vias
                 < static_cast<double>(seeded_cost.vias) * 0.92;
    const bool via_strongly_improved
        = fallback_cost.vias < static_cast<double>(seeded_cost.vias) * 0.78
          && fallback_cost.wirelength
                 < static_cast<double>(seeded_cost.wirelength) * 1.10;
    const bool giant_net_wire_bias
        = pin_count >= 24
          && fallback_cost.wirelength
                 < static_cast<double>(seeded_cost.wirelength) * 0.99
          && fallback_cost.vias
                 <= static_cast<double>(seeded_cost.vias) * 1.02;
    const bool force_geometric_route = pin_count >= 2;
    const bool should_replace
        = force_geometric_route || (fallback_score < seeded_score * 0.90)
          || strong_wire_improved || balanced_improved || via_strongly_improved
          || giant_net_wire_bias;

    if (should_replace) {
      route_it->second = *fallback_route;
      ++replaced_net_count;
    }
  }
  if (missing_net_fallback_count > 0 || replaced_net_count > 0) {
    logger_->report(
        "NEWGR fallback routing: {} unrouted nets filled, {} FastRoute nets "
        "replaced with wire-focused geometric guides.",
        missing_net_fallback_count,
        replaced_net_count);
  }
  return routes;
}

bool NewgrEngine::appendRouteSegments(int net_id, GRoute& route) const
{
  if (net_id < 0 || net_id >= numValidNets) {
    return false;
  }
  const StTree& tree = sttrees[net_id];
  if (tree.deg == 0 || tree.edges == nullptr) {
    return true;
  }

  const int max_reasonable_routelen
      = std::max(1024, grid_.x_grids * grid_.y_grids * std::max(1, grid_.num_layers));
  const int edge_count = 2 * tree.deg - 3;
  for (int edge_id = 0; edge_id < edge_count; ++edge_id) {
    const TreeEdge& tree_edge = tree.edges[edge_id];
    const Route& edge_route = tree_edge.route;
    if (edge_route.routelen <= 0 || edge_route.gridsX == nullptr
        || edge_route.gridsY == nullptr || edge_route.gridsL == nullptr) {
      continue;
    }
    if (edge_route.routelen > max_reasonable_routelen) {
      return false;
    }
    const int route_steps = std::max(0, edge_route.routelen - 1);
    for (int i = 0; i < route_steps; ++i) {
      const int x0 = edge_route.gridsX[i];
      const int y0 = edge_route.gridsY[i];
      const int l0 = edge_route.gridsL[i];
      const int x1 = edge_route.gridsX[i + 1];
      const int y1 = edge_route.gridsY[i + 1];
      const int l1 = edge_route.gridsL[i + 1];
      if (!isGridPointValid(x0, y0, l0) || !isGridPointValid(x1, y1, l1)) {
        return false;
      }

      const int manhattan_delta
          = std::abs(x1 - x0) + std::abs(y1 - y0) + std::abs(l1 - l0);
      if (manhattan_delta == 0) {
        continue;
      }
      if (manhattan_delta == 1) {
        addSegment(route, x0, y0, l0, x1, y1, l1);
        continue;
      }
      if (!appendManhattanBridge(x0, y0, l0, x1, y1, l1, route)) {
        return false;
      }
    }
  }

  return true;
}

bool NewgrEngine::appendFallbackMstRoute(const NewgrInputNet& net,
                                         GRoute& route) const
{
  if (net.pins.size() < 2) {
    return false;
  }

  const int pin_count = static_cast<int>(net.pins.size());
  int seed_pin = net.root_pin_index;
  if (seed_pin < 0 || seed_pin >= pin_count) {
    seed_pin = 0;
  }
  std::vector<int> pin_layers;
  pin_layers.reserve(pin_count);
  for (const auto& pin : net.pins) {
    pin_layers.push_back(pin.layer());
  }
  std::sort(pin_layers.begin(), pin_layers.end());
  int preferred_layer = pin_layers[pin_layers.size() / 2];
  preferred_layer = std::clamp(preferred_layer, 0, std::max(0, grid_.num_layers - 1));
  const int layer_mismatch_weight = (pin_count >= 8) ? 4 : 3;
  const bool use_layer_spine = pin_count >= 5;

  std::vector<bool> in_tree(pin_count, false);
  std::vector<bool> pin_access_added(pin_count, false);
  in_tree[seed_pin] = true;
  int connected_count = 1;

  auto ensurePinAccess = [&](int pin_index, int spine_layer) -> bool {
    if (pin_index < 0 || pin_index >= pin_count) {
      return false;
    }
    if (pin_access_added[pin_index]) {
      return true;
    }
    const RoutePt& pin = net.pins[pin_index];
    if (!appendManhattanBridge(pin.x(),
                               pin.y(),
                               pin.layer(),
                               pin.x(),
                               pin.y(),
                               spine_layer,
                               route)) {
      return false;
    }
    pin_access_added[pin_index] = true;
    return true;
  };

  auto weightedDistance = [&](int lhs_idx, int rhs_idx) {
    const RoutePt& lhs = net.pins[lhs_idx];
    const RoutePt& rhs = net.pins[rhs_idx];
    const int lhs_bias = std::abs(lhs.layer() - preferred_layer);
    const int rhs_bias = std::abs(rhs.layer() - preferred_layer);
    return std::abs(lhs.x() - rhs.x()) + std::abs(lhs.y() - rhs.y())
           + layer_mismatch_weight * std::abs(lhs.layer() - rhs.layer())
           + lhs_bias + rhs_bias;
  };

  while (connected_count < pin_count) {
    int best_u = -1;
    int best_v = -1;
    int best_cost = std::numeric_limits<int>::max();
    for (int u = 0; u < pin_count; ++u) {
      if (!in_tree[u]) {
        continue;
      }
      for (int v = 0; v < pin_count; ++v) {
        if (in_tree[v]) {
          continue;
        }
        const int cost = weightedDistance(u, v);
        if (cost < best_cost) {
          best_cost = cost;
          best_u = u;
          best_v = v;
        }
      }
    }

    if (best_u < 0 || best_v < 0) {
      return false;
    }

    const RoutePt& src = net.pins[best_u];
    const RoutePt& dst = net.pins[best_v];
    const int src_spine_layer = use_layer_spine ? preferred_layer : src.layer();
    const int dst_spine_layer = use_layer_spine ? preferred_layer : dst.layer();

    if (!ensurePinAccess(best_u, src_spine_layer)
        || !ensurePinAccess(best_v, dst_spine_layer)
        || !appendManhattanBridge(src.x(),
                                  src.y(),
                                  src_spine_layer,
                                  dst.x(),
                                  dst.y(),
                                  dst_spine_layer,
                                  route)) {
      return false;
    }
    in_tree[best_v] = true;
    ++connected_count;
  }

  return true;
}

bool NewgrEngine::appendFallbackBackboneRoute(const NewgrInputNet& net,
                                              GRoute& route) const
{
  if (net.pins.size() < 2) {
    return false;
  }

  std::vector<int> xs;
  std::vector<int> ys;
  std::vector<int> ls;
  xs.reserve(net.pins.size());
  ys.reserve(net.pins.size());
  ls.reserve(net.pins.size());

  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int min_y = std::numeric_limits<int>::max();
  int max_y = std::numeric_limits<int>::min();
  for (const auto& pin : net.pins) {
    xs.push_back(pin.x());
    ys.push_back(pin.y());
    ls.push_back(pin.layer());
    min_x = std::min(min_x, pin.x());
    max_x = std::max(max_x, pin.x());
    min_y = std::min(min_y, pin.y());
    max_y = std::max(max_y, pin.y());
  }

  std::sort(xs.begin(), xs.end());
  std::sort(ys.begin(), ys.end());
  std::sort(ls.begin(), ls.end());
  const int median_x = xs[xs.size() / 2];
  const int median_y = ys[ys.size() / 2];
  const int preferred_layer
      = std::clamp(ls[ls.size() / 2], 0, std::max(0, grid_.num_layers - 1));

  int64_t branch_cost_vertical = static_cast<int64_t>(max_y) - min_y;
  int64_t branch_cost_horizontal = static_cast<int64_t>(max_x) - min_x;
  for (const auto& pin : net.pins) {
    branch_cost_vertical += std::llabs(static_cast<long long>(pin.x()) - median_x);
    branch_cost_horizontal += std::llabs(static_cast<long long>(pin.y()) - median_y);
  }
  const bool vertical_backbone = branch_cost_vertical <= branch_cost_horizontal;

  if (vertical_backbone) {
    if (!appendManhattanBridge(
            median_x, min_y, preferred_layer, median_x, max_y, preferred_layer, route)) {
      return false;
    }
  } else {
    if (!appendManhattanBridge(
            min_x, median_y, preferred_layer, max_x, median_y, preferred_layer, route)) {
      return false;
    }
  }

  for (const auto& pin : net.pins) {
    if (!appendManhattanBridge(pin.x(),
                               pin.y(),
                               pin.layer(),
                               pin.x(),
                               pin.y(),
                               preferred_layer,
                               route)) {
      return false;
    }

    const int target_x = vertical_backbone ? median_x : pin.x();
    const int target_y = vertical_backbone ? pin.y() : median_y;
    if (!appendManhattanBridge(
            pin.x(), pin.y(), preferred_layer, target_x, target_y, preferred_layer, route)) {
      return false;
    }
  }

  return true;
}

bool NewgrEngine::appendFallbackCrossRoute(const NewgrInputNet& net,
                                           GRoute& route) const
{
  if (net.pins.size() < 2) {
    return false;
  }
  if (net.pins.size() < 4) {
    return appendFallbackBackboneRoute(net, route);
  }

  std::vector<int> xs;
  std::vector<int> ys;
  std::vector<int> ls;
  xs.reserve(net.pins.size());
  ys.reserve(net.pins.size());
  ls.reserve(net.pins.size());

  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int min_y = std::numeric_limits<int>::max();
  int max_y = std::numeric_limits<int>::min();
  for (const auto& pin : net.pins) {
    xs.push_back(pin.x());
    ys.push_back(pin.y());
    ls.push_back(pin.layer());
    min_x = std::min(min_x, pin.x());
    max_x = std::max(max_x, pin.x());
    min_y = std::min(min_y, pin.y());
    max_y = std::max(max_y, pin.y());
  }

  std::sort(xs.begin(), xs.end());
  std::sort(ys.begin(), ys.end());
  std::sort(ls.begin(), ls.end());
  const int median_x = xs[xs.size() / 2];
  const int median_y = ys[ys.size() / 2];
  const int preferred_layer
      = std::clamp(ls[ls.size() / 2], 0, std::max(0, grid_.num_layers - 1));

  if (!appendManhattanBridge(
          median_x, min_y, preferred_layer, median_x, max_y, preferred_layer, route)
      || !appendManhattanBridge(
          min_x, median_y, preferred_layer, max_x, median_y, preferred_layer, route)) {
    return false;
  }

  for (const auto& pin : net.pins) {
    if (!appendManhattanBridge(pin.x(),
                               pin.y(),
                               pin.layer(),
                               pin.x(),
                               pin.y(),
                               preferred_layer,
                               route)) {
      return false;
    }

    const int target_vx = median_x;
    const int target_vy = std::clamp(pin.y(), min_y, max_y);
    const int target_hx = std::clamp(pin.x(), min_x, max_x);
    const int target_hy = median_y;
    const int vertical_distance = std::abs(pin.x() - target_vx);
    const int horizontal_distance = std::abs(pin.y() - target_hy);

    const int target_x
        = (vertical_distance <= horizontal_distance) ? target_vx : target_hx;
    const int target_y
        = (vertical_distance <= horizontal_distance) ? target_vy : target_hy;

    if (!appendManhattanBridge(pin.x(),
                               pin.y(),
                               preferred_layer,
                               target_x,
                               target_y,
                               preferred_layer,
                               route)) {
      return false;
    }
  }

  return true;
}

bool NewgrEngine::appendFallbackDualHubRoute(const NewgrInputNet& net,
                                             GRoute& route) const
{
  const int pin_count = static_cast<int>(net.pins.size());
  if (pin_count < 4) {
    return false;
  }

  std::vector<int> xs;
  std::vector<int> ys;
  std::vector<int> ls;
  xs.reserve(pin_count);
  ys.reserve(pin_count);
  ls.reserve(pin_count);

  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  int min_y = std::numeric_limits<int>::max();
  int max_y = std::numeric_limits<int>::min();
  for (const auto& pin : net.pins) {
    xs.push_back(pin.x());
    ys.push_back(pin.y());
    ls.push_back(pin.layer());
    min_x = std::min(min_x, pin.x());
    max_x = std::max(max_x, pin.x());
    min_y = std::min(min_y, pin.y());
    max_y = std::max(max_y, pin.y());
  }

  std::sort(xs.begin(), xs.end());
  std::sort(ys.begin(), ys.end());
  std::sort(ls.begin(), ls.end());

  const int q1 = (pin_count - 1) / 4;
  const int q3 = ((pin_count - 1) * 3) / 4;
  const int median = pin_count / 2;
  const int preferred_layer
      = std::clamp(ls[median], 0, std::max(0, grid_.num_layers - 1));
  const bool dominant_x = (max_x - min_x) >= (max_y - min_y);

  if (dominant_x) {
    int hub1_x = xs[q1];
    int hub2_x = xs[q3];
    if (hub1_x > hub2_x) {
      std::swap(hub1_x, hub2_x);
    }
    const int hub_y = ys[median];

    if (!appendManhattanBridge(
            hub1_x, hub_y, preferred_layer, hub2_x, hub_y, preferred_layer, route)) {
      return false;
    }

    for (const auto& pin : net.pins) {
      if (!appendManhattanBridge(pin.x(),
                                 pin.y(),
                                 pin.layer(),
                                 pin.x(),
                                 pin.y(),
                                 preferred_layer,
                                 route)) {
        return false;
      }

      const int d1 = std::abs(pin.x() - hub1_x) + std::abs(pin.y() - hub_y);
      const int d2 = std::abs(pin.x() - hub2_x) + std::abs(pin.y() - hub_y);
      const int hub_x = (d1 <= d2) ? hub1_x : hub2_x;

      if (!appendManhattanBridge(pin.x(),
                                 pin.y(),
                                 preferred_layer,
                                 hub_x,
                                 pin.y(),
                                 preferred_layer,
                                 route)
          || !appendManhattanBridge(hub_x,
                                    pin.y(),
                                    preferred_layer,
                                    hub_x,
                                    hub_y,
                                    preferred_layer,
                                    route)) {
        return false;
      }
    }
  } else {
    int hub1_y = ys[q1];
    int hub2_y = ys[q3];
    if (hub1_y > hub2_y) {
      std::swap(hub1_y, hub2_y);
    }
    const int hub_x = xs[median];

    if (!appendManhattanBridge(
            hub_x, hub1_y, preferred_layer, hub_x, hub2_y, preferred_layer, route)) {
      return false;
    }

    for (const auto& pin : net.pins) {
      if (!appendManhattanBridge(pin.x(),
                                 pin.y(),
                                 pin.layer(),
                                 pin.x(),
                                 pin.y(),
                                 preferred_layer,
                                 route)) {
        return false;
      }

      const int d1 = std::abs(pin.x() - hub_x) + std::abs(pin.y() - hub1_y);
      const int d2 = std::abs(pin.x() - hub_x) + std::abs(pin.y() - hub2_y);
      const int hub_y = (d1 <= d2) ? hub1_y : hub2_y;

      if (!appendManhattanBridge(pin.x(),
                                 pin.y(),
                                 preferred_layer,
                                 pin.x(),
                                 hub_y,
                                 preferred_layer,
                                 route)
          || !appendManhattanBridge(pin.x(),
                                    hub_y,
                                    preferred_layer,
                                    hub_x,
                                    hub_y,
                                    preferred_layer,
                                    route)) {
        return false;
      }
    }
  }

  return true;
}

bool NewgrEngine::appendManhattanBridge(int grid_x0,
                                        int grid_y0,
                                        int grid_l0,
                                        int grid_x1,
                                        int grid_y1,
                                        int grid_l1,
                                        GRoute& route) const
{
  if (!isGridPointValid(grid_x0, grid_y0, grid_l0)
      || !isGridPointValid(grid_x1, grid_y1, grid_l1)) {
    return false;
  }

  int x = grid_x0;
  int y = grid_y0;
  int l = grid_l0;
  const int max_steps = std::max(32,
                                 2 * (std::abs(grid_x1 - grid_x0)
                                      + std::abs(grid_y1 - grid_y0)
                                      + std::abs(grid_l1 - grid_l0))
                                     + 8);
  int step_count = 0;

  while (l != grid_l1) {
    if (++step_count > max_steps) {
      return false;
    }
    const int next_l = l + ((grid_l1 > l) ? 1 : -1);
    if (!isGridPointValid(x, y, next_l)) {
      return false;
    }
    addSegment(route, x, y, l, x, y, next_l);
    l = next_l;
  }
  while (x != grid_x1) {
    if (++step_count > max_steps) {
      return false;
    }
    const int next_x = x + ((grid_x1 > x) ? 1 : -1);
    if (!isGridPointValid(next_x, y, l)) {
      return false;
    }
    addSegment(route, x, y, l, next_x, y, l);
    x = next_x;
  }
  while (y != grid_y1) {
    if (++step_count > max_steps) {
      return false;
    }
    const int next_y = y + ((grid_y1 > y) ? 1 : -1);
    if (!isGridPointValid(x, next_y, l)) {
      return false;
    }
    addSegment(route, x, y, l, x, next_y, l);
    y = next_y;
  }

  return true;
}

bool NewgrEngine::isGridPointValid(int grid_x, int grid_y, int grid_l) const
{
  return grid_x >= 0 && grid_x < grid_.x_grids && grid_y >= 0
         && grid_y < grid_.y_grids && grid_l >= 0 && grid_l < grid_.num_layers;
}

void NewgrEngine::addSegment(GRoute& route,
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

int NewgrEngine::toDbLayer(int sproute_layer) const
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

std::string NewgrEngine::sanitizeNetName(const std::string& name) const
{
  constexpr size_t kLimit = 96;  // SPRoute internal buffers default to 100.
  if (name.size() < kLimit) {
    return name;
  }
  const size_t hash_value = std::hash<std::string>{}(name);
  std::ostringstream oss;
  oss << std::hex << hash_value;
  const std::string hash = oss.str();
  const size_t suffix_len = hash.size() + 1;  // include '_' separator.
  const size_t prefix_len = (suffix_len < kLimit) ? (kLimit - suffix_len) : 0;
  std::string result = name.substr(0, prefix_len);
  result.push_back('_');
  result += hash;
  return result;
}

void NewgrEngine::updateDbCongestion(odb::dbBlock* block)
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
