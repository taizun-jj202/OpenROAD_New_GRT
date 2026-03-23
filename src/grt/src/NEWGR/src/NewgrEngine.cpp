#include "NewgrEngine.h"

#include <algorithm>
#include <cstdlib>
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

struct GridPoint
{
  int x{0};
  int y{0};
  int l{0};
};

struct RouteMetrics
{
  uint64_t wirelength{0};
  uint64_t vias{0};
  uint64_t congestion_risk{0};
  uint64_t proxy_cost{0};
};

struct CandidateResult
{
  Algo algo{Algo::Astar};
  std::string mode_name;
  int capacity_profile{NEWGR_CAP_PROFILE_BALANCED};
  int maze_rounds{0};
  int overflow{std::numeric_limits<int>::max()};
  RouteMetrics metrics{};
  NetRouteMap routes;
};

uint64_t computeCongestionRisk(const SprouteGridData& grid)
{
  uint64_t risk = 0;
  for (int layer = 0; layer < grid.num_layers; ++layer) {
    const uint64_t layer_weight = (layer <= 2) ? 9 : ((layer <= 4) ? 6 : 4);
    if (grid.x_grids > 1) {
      for (int y = 0; y < grid.y_grids; ++y) {
        for (int x = 0; x + 1 < grid.x_grids; ++x) {
          const Edge3D& edge = horizontalEdge(grid, layer, y, x);
          const int cap = std::max(1, static_cast<int>(edge.cap));
          const int usage = static_cast<int>(edge.usage);
          const int util_permil = (usage * 1000) / cap;
          if (util_permil > 700) {
            const uint64_t pressure = static_cast<uint64_t>(util_permil - 700);
            risk += pressure * pressure * layer_weight;
          }
          if (usage > cap) {
            risk += static_cast<uint64_t>(usage - cap) * 3000 * layer_weight;
          }
        }
      }
    }
    if (grid.y_grids > 1) {
      for (int y = 0; y + 1 < grid.y_grids; ++y) {
        for (int x = 0; x < grid.x_grids; ++x) {
          const Edge3D& edge = verticalEdge(grid, layer, y, x);
          const int cap = std::max(1, static_cast<int>(edge.cap));
          const int usage = static_cast<int>(edge.usage);
          const int util_permil = (usage * 1000) / cap;
          if (util_permil > 700) {
            const uint64_t pressure = static_cast<uint64_t>(util_permil - 700);
            risk += pressure * pressure * layer_weight;
          }
          if (usage > cap) {
            risk += static_cast<uint64_t>(usage - cap) * 3000 * layer_weight;
          }
        }
      }
    }
  }
  return risk;
}

RouteMetrics computeRouteMetrics(const SprouteGridData& grid,
                                 const NetRouteMap& routes)
{
  RouteMetrics metrics;
  for (const auto& [net, route] : routes) {
    (void) net;
    for (const GSegment& segment : route) {
      metrics.wirelength += static_cast<uint64_t>(segment.length());
      if (segment.isVia()) {
        metrics.vias += static_cast<uint64_t>(
            std::abs(segment.final_layer - segment.init_layer));
      }
    }
  }
  metrics.congestion_risk = computeCongestionRisk(grid);
  // Proxy is intentionally risk-heavy because detailed-route wirelength
  // regresses when global-route hotspots are left unresolved.
  metrics.proxy_cost = metrics.wirelength + metrics.vias * 12
                       + metrics.congestion_risk / 24;
  return metrics;
}

uint64_t selectionCost(const RouteMetrics& metrics)
{
  // Wirelength-first ranking with only a very light risk guard.
  return metrics.wirelength + metrics.vias * 8 + metrics.congestion_risk / 3000;
}

bool isBetterCandidate(const CandidateResult& lhs, const CandidateResult& rhs)
{
  if (lhs.overflow != rhs.overflow) {
    return lhs.overflow < rhs.overflow;
  }

  const uint64_t lhs_cost = selectionCost(lhs.metrics);
  const uint64_t rhs_cost = selectionCost(rhs.metrics);
  if (lhs_cost != rhs_cost) {
    return lhs_cost < rhs_cost;
  }
  if (lhs.metrics.wirelength != rhs.metrics.wirelength) {
    return lhs.metrics.wirelength < rhs.metrics.wirelength;
  }
  if (lhs.metrics.vias != rhs.metrics.vias) {
    return lhs.metrics.vias < rhs.metrics.vias;
  }
  if (lhs.metrics.congestion_risk != rhs.metrics.congestion_risk) {
    return lhs.metrics.congestion_risk < rhs.metrics.congestion_risk;
  }
  return false;
}

const char* algoName(Algo algo)
{
  switch (algo) {
    case Algo::DetPart_Astar_Local:
      return "DetPart_Astar_Local";
    case Algo::Astar:
      return "Astar";
    case Algo::FineGrain:
      return "FineGrain";
    default:
      return "Other";
  }
}

const char* capacityProfileName(int profile)
{
  switch (profile) {
    case NEWGR_CAP_PROFILE_3D_SHORT:
      return "3D_SHORT";
    case NEWGR_CAP_PROFILE_ULTRA_WL:
      return "ULTRA_WL";
    case NEWGR_CAP_PROFILE_RADICAL_WL:
      return "RADICAL_WL";
    case NEWGR_CAP_PROFILE_RADICAL_MIX:
      return "RADICAL_MIX";
    case NEWGR_CAP_PROFILE_WL_FOCUSED:
      return "WL_FOCUSED";
    case NEWGR_CAP_PROFILE_DR_FOCUSED:
      return "DR_FOCUSED";
    case NEWGR_CAP_PROFILE_BALANCED:
    default:
      return "BALANCED";
  }
}

std::string candidateName(const CandidateResult& candidate)
{
  if (!candidate.mode_name.empty()) {
    return candidate.mode_name;
  }
  return std::string(algoName(candidate.algo));
}

uint64_t packPoint(const GridPoint& point)
{
  constexpr uint64_t mask = 0x1FFFFF;
  return (static_cast<uint64_t>(point.x) & mask)
         | ((static_cast<uint64_t>(point.y) & mask) << 21)
         | ((static_cast<uint64_t>(point.l) & mask) << 42);
}

bool canMergeMiddlePoint(const GridPoint& a,
                         const GridPoint& b,
                         const GridPoint& c)
{
  if (a.l == b.l && b.l == c.l) {
    if (a.x == b.x && b.x == c.x) {
      return true;
    }
    if (a.y == b.y && b.y == c.y) {
      return true;
    }
  }
  if (a.x == b.x && b.x == c.x && a.y == b.y && b.y == c.y) {
    return std::abs(a.l - c.l) <= 1;
  }
  return false;
}

std::vector<GridPoint> simplifyGridPath(const Route& edge_route)
{
  std::vector<GridPoint> path;
  path.reserve(edge_route.routelen + 1);
  for (int i = 0; i <= edge_route.routelen; ++i) {
    path.push_back(
        {edge_route.gridsX[i], edge_route.gridsY[i], edge_route.gridsL[i]});
  }

  std::vector<GridPoint> simple_path;
  simple_path.reserve(path.size());
  std::unordered_map<uint64_t, size_t> last_seen;
  last_seen.reserve(path.size() * 2);

  for (const GridPoint& point : path) {
    const uint64_t key = packPoint(point);
    const auto found = last_seen.find(key);
    if (found != last_seen.end()) {
      const size_t keep = found->second + 1;
      for (size_t i = keep; i < simple_path.size(); ++i) {
        last_seen.erase(packPoint(simple_path[i]));
      }
      simple_path.resize(keep);
      continue;
    }
    last_seen[key] = simple_path.size();
    simple_path.push_back(point);
  }

  std::vector<GridPoint> compact_path;
  compact_path.reserve(simple_path.size());
  for (const GridPoint& point : simple_path) {
    compact_path.push_back(point);
    while (compact_path.size() >= 3) {
      const size_t n = compact_path.size();
      const GridPoint& a = compact_path[n - 3];
      const GridPoint& b = compact_path[n - 2];
      const GridPoint& c = compact_path[n - 1];
      if (!canMergeMiddlePoint(a, b, c)) {
        break;
      }
      compact_path.erase(compact_path.end() - 2);
    }
  }

  return compact_path;
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

  auto run_candidate = [&](Algo algo,
                           int maze_rounds,
                           int capacity_profile,
                           const char* mode_name) {
    parser::grGenerator generator = buildGenerator();
    if (generator.capReductions_p == nullptr) {
      logger_->warn(utl::GRT,
                    401,
                    "NEWGR generator has no localized capacity reductions; "
                    "continuing without adjustments.");
    }

    CandidateResult candidate;
    candidate.algo = algo;
    candidate.mode_name = mode_name;
    candidate.capacity_profile = capacity_profile;
    candidate.maze_rounds = maze_rounds;
    newgr_capacity_profile = capacity_profile;

    parser::CongestionMap congestion_map(
        generator.grid.z, generator.grid.x, generator.grid.y);
    galois::StatTimer timer("newgr_timer");
    timer.start();
    runFastRoute(generator,
                 /*benchFile=*/"",
                 /*OutFileName=*/"",
                 congestion_map,
                 timer,
                 maze_rounds,
                 algo);
    timer.stop();

    candidate.overflow = totalOverflow;
    candidate.routes = extractRoutes();
    candidate.metrics = computeRouteMetrics(grid_, candidate.routes);

    logger_->info(utl::GRT,
                  402,
                  "NEWGR candidate {} [{}]: overflow={}, route_wl={}, "
                  "route_vias={}, cong_risk={}, proxy_cost={}",
                  candidateName(candidate),
                  capacityProfileName(candidate.capacity_profile),
                  candidate.overflow,
                  candidate.metrics.wirelength,
                  candidate.metrics.vias,
                  candidate.metrics.congestion_risk,
                  candidate.metrics.proxy_cost);
    return candidate;
  };

  // NEWGR radical ensemble:
  // - RADICAL_WL keeps shortest-path pressure high.
  // - RADICAL_MIX adds deterministic hotspot diffusion from SPRoute ideas.
  // - ULTRA_WL provides a tighter-box shortest-path alternative.
  CandidateResult best = run_candidate(Algo::Astar,
                                       520,
                                       NEWGR_CAP_PROFILE_RADICAL_WL,
                                       "Astar_RadicalWL");
  CandidateResult radical_mix = run_candidate(Algo::Astar,
                                              500,
                                              NEWGR_CAP_PROFILE_RADICAL_MIX,
                                              "Astar_RadicalMix");
  if (isBetterCandidate(radical_mix, best)) {
    best = std::move(radical_mix);
  }
  CandidateResult ultra_wl = run_candidate(Algo::Astar,
                                           460,
                                           NEWGR_CAP_PROFILE_ULTRA_WL,
                                           "Astar_UltraWL");
  if (isBetterCandidate(ultra_wl, best)) {
    best = std::move(ultra_wl);
  }

  if (best.overflow > 0) {
    CandidateResult short3d = run_candidate(Algo::Astar,
                                            620,
                                            NEWGR_CAP_PROFILE_3D_SHORT,
                                            "Astar_3DShort_Overflow");
    if (isBetterCandidate(short3d, best)) {
      best = std::move(short3d);
    }
    CandidateResult fallback = run_candidate(Algo::DetPart_Astar_Local,
                                             560,
                                             NEWGR_CAP_PROFILE_DR_FOCUSED,
                                             "DetPart_DRFallback");
    if (isBetterCandidate(fallback, best)) {
      best = std::move(fallback);
    }
  }

  newgr_capacity_profile = best.capacity_profile;
  last_total_overflow_ = best.overflow;
  logger_->info(utl::GRT,
                404,
                "NEWGR selected candidate {} [{}]: overflow={}, route_wl={}, "
                "route_vias={}, cong_risk={}, proxy_cost={}",
                candidateName(best),
                capacityProfileName(best.capacity_profile),
                best.overflow,
                best.metrics.wirelength,
                best.metrics.vias,
                best.metrics.congestion_risk,
                best.metrics.proxy_cost);
  return best.routes;
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
    appendRouteSegments(net_id, route_segments);
    if (!route_segments.empty()) {
      routes[it->second] = std::move(route_segments);
    }
  }

  return routes;
}

void NewgrEngine::appendRouteSegments(int net_id, GRoute& route) const
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
    const std::vector<GridPoint> simplified_path = simplifyGridPath(edge_route);
    if (simplified_path.size() < 2) {
      continue;
    }
    for (size_t i = 0; i + 1 < simplified_path.size(); ++i) {
      const GridPoint& p0 = simplified_path[i];
      const GridPoint& p1 = simplified_path[i + 1];
      if (p0.x == p1.x && p0.y == p1.y && p0.l == p1.l) {
        continue;
      }
      addSegment(route, p0.x, p0.y, p0.l, p1.x, p1.y, p1.l);
    }
  }
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
