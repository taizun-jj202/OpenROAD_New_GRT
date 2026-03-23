#include "NEWGR/NewGR.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "CUGR.h"
#include "FastRoute.h"
#include "Grid.h"
#include "Net.h"
#include "grt/Rudy.h"
#include "grt/SprouteAdapter.h"
#include "utl/Logger.h"

namespace grt {

using utl::GNR;

namespace {

struct RouteMetrics
{
  long wirelength_dbu = 0;
  long via_count = 0;
  double wirelength_um = 0.0;
  double score = 0.0;
};

struct ScenarioResult
{
  std::string name;
  RouteMetrics metrics;
  NetRouteMap routes;
};

struct RouterSnapshot
{
  float caps_percentage = 0.0f;
  int perturbation_amount = 0;
  float critical_percentage = 0.0f;
  bool allow_congestion = false;
  int seed = 0;
};

struct ScenarioDefinition
{
  std::string name;
  std::function<void()> pre_init;
  std::function<void()> post_init;
  std::function<void(std::vector<Net*>&)> order_nets;
  bool aggressive = false;
};

struct Hotspot
{
  int gx = 0;
  int gy = 0;
  float severity = 1.0f;
  bool affect_horizontal = false;
  bool affect_vertical = false;
};

struct NetOrderFeatures
{
  long hpwl = 0;
  int x_span = 0;
  int y_span = 0;
  int center_x = 0;
  int center_y = 0;
  int pin_count = 0;
  int layer_span = 1;
  int port_count = 0;
  float slack = 0.0f;
  double score = 0.0;
};

struct NetRouteCost
{
  long wirelength_dbu = 0;
  long via_count = 0;
  double hotspot_exposure = 0.0;
  double score = 0.0;
};

NetOrderFeatures getNetOrderFeatures(Net* net)
{
  NetOrderFeatures features;
  if (net == nullptr) {
    return features;
  }

  const auto& pins = net->getPins();
  features.pin_count = std::max(1, static_cast<int>(pins.size()));
  features.slack = net->getSlack();

  int xmin = std::numeric_limits<int>::max();
  int xmax = std::numeric_limits<int>::min();
  int ymin = std::numeric_limits<int>::max();
  int ymax = std::numeric_limits<int>::min();
  int min_layer = std::numeric_limits<int>::max();
  int max_layer = std::numeric_limits<int>::min();

  for (const Pin& pin : pins) {
    const odb::Point& pos = pin.getPosition();
    xmin = std::min(xmin, pos.x());
    xmax = std::max(xmax, pos.x());
    ymin = std::min(ymin, pos.y());
    ymax = std::max(ymax, pos.y());

    if (pin.isPort()) {
      features.port_count++;
    }

    const auto& pin_layers = pin.getLayers();
    for (const int layer : pin_layers) {
      min_layer = std::min(min_layer, layer);
      max_layer = std::max(max_layer, layer);
    }
  }

  if (xmin <= xmax && ymin <= ymax) {
    features.hpwl = static_cast<long>(xmax - xmin + ymax - ymin);
    features.x_span = std::max(0, xmax - xmin);
    features.y_span = std::max(0, ymax - ymin);
    features.center_x = xmin + (features.x_span / 2);
    features.center_y = ymin + (features.y_span / 2);
  }
  if (min_layer <= max_layer) {
    features.layer_span = std::max(1, max_layer - min_layer + 1);
  }

  const double hpwl_term = static_cast<double>(std::max<long>(features.hpwl, 1));
  const double pin_term
      = 1.0 + 0.28 * std::log1p(static_cast<double>(features.pin_count));
  const double layer_term = 1.0 + 0.16 * std::max(features.layer_span - 1, 0);
  const double port_term = 1.0 + 0.08 * features.port_count;
  const double slack_term = features.slack < 0.0f
                                ? 1.0 + std::min(1.2, std::abs(static_cast<double>(features.slack)))
                                : 1.0;
  features.score = hpwl_term * pin_term * layer_term * port_term * slack_term;

  return features;
}

void reorderNetsByWirelengthPriority(std::vector<Net*>& nets)
{
  if (nets.size() < 2) {
    return;
  }

  std::stable_sort(nets.begin(), nets.end(), [](Net* lhs, Net* rhs) {
    const NetOrderFeatures lhs_features = getNetOrderFeatures(lhs);
    const NetOrderFeatures rhs_features = getNetOrderFeatures(rhs);
    if (lhs_features.score != rhs_features.score) {
      return lhs_features.score > rhs_features.score;
    }
    if (lhs_features.hpwl != rhs_features.hpwl) {
      return lhs_features.hpwl > rhs_features.hpwl;
    }
    if (lhs_features.pin_count != rhs_features.pin_count) {
      return lhs_features.pin_count > rhs_features.pin_count;
    }
    return lhs < rhs;
  });
}

void reorderNetsBySpatialRoundRobin(std::vector<Net*>& nets)
{
  if (nets.size() < 2) {
    return;
  }

  struct NetOrderEntry
  {
    Net* net = nullptr;
    NetOrderFeatures features;
  };

  std::vector<NetOrderEntry> entries;
  entries.reserve(nets.size());
  int min_center_x = std::numeric_limits<int>::max();
  int max_center_x = std::numeric_limits<int>::min();

  for (Net* net : nets) {
    NetOrderEntry entry;
    entry.net = net;
    entry.features = getNetOrderFeatures(net);
    min_center_x = std::min(min_center_x, entry.features.center_x);
    max_center_x = std::max(max_center_x, entry.features.center_x);
    entries.push_back(entry);
  }

  std::stable_sort(entries.begin(),
                   entries.end(),
                   [](const NetOrderEntry& lhs, const NetOrderEntry& rhs) {
                     const int lhs_span = std::max(lhs.features.x_span, lhs.features.y_span);
                     const int rhs_span = std::max(rhs.features.x_span, rhs.features.y_span);
                     if (lhs_span != rhs_span) {
                       return lhs_span > rhs_span;
                     }
                     if (lhs.features.score != rhs.features.score) {
                       return lhs.features.score > rhs.features.score;
                     }
                     if (lhs.features.hpwl != rhs.features.hpwl) {
                       return lhs.features.hpwl > rhs.features.hpwl;
                     }
                     return lhs.net < rhs.net;
                   });

  const int bucket_count
      = std::clamp(static_cast<int>(nets.size() / 600), 4, 12);
  std::vector<std::vector<NetOrderEntry>> buckets(bucket_count);
  const int x_span = std::max(1, max_center_x - min_center_x + 1);

  for (const NetOrderEntry& entry : entries) {
    int bucket_id
        = ((entry.features.center_x - min_center_x) * bucket_count) / x_span;
    bucket_id = std::clamp(bucket_id, 0, bucket_count - 1);
    buckets[bucket_id].push_back(entry);
  }

  for (std::vector<NetOrderEntry>& bucket : buckets) {
    std::stable_sort(
        bucket.begin(), bucket.end(), [](const NetOrderEntry& lhs, const NetOrderEntry& rhs) {
          if (lhs.features.score != rhs.features.score) {
            return lhs.features.score > rhs.features.score;
          }
          if (lhs.features.center_y != rhs.features.center_y) {
            return lhs.features.center_y < rhs.features.center_y;
          }
          return lhs.net < rhs.net;
        });
  }

  std::vector<size_t> idx(bucket_count, 0);
  size_t emitted = 0;
  nets.clear();
  nets.reserve(entries.size());

  while (emitted < entries.size()) {
    bool progress = false;
    for (int b = 0; b < bucket_count; ++b) {
      if (idx[b] >= buckets[b].size()) {
        continue;
      }
      nets.push_back(buckets[b][idx[b]].net);
      idx[b]++;
      emitted++;
      progress = true;
    }
    if (!progress) {
      break;
    }
  }
}

void reorderNetsByBspScheduler(std::vector<Net*>& nets)
{
  if (nets.size() < 2) {
    return;
  }

  struct NetOrderEntry
  {
    Net* net = nullptr;
    NetOrderFeatures features;
  };

  std::vector<NetOrderEntry> entries;
  entries.reserve(nets.size());

  int min_center_x = std::numeric_limits<int>::max();
  int max_center_x = std::numeric_limits<int>::min();
  int min_center_y = std::numeric_limits<int>::max();
  int max_center_y = std::numeric_limits<int>::min();

  for (Net* net : nets) {
    NetOrderEntry entry;
    entry.net = net;
    entry.features = getNetOrderFeatures(net);
    min_center_x = std::min(min_center_x, entry.features.center_x);
    max_center_x = std::max(max_center_x, entry.features.center_x);
    min_center_y = std::min(min_center_y, entry.features.center_y);
    max_center_y = std::max(max_center_y, entry.features.center_y);
    entries.push_back(entry);
  }

  const bool sort_by_x = (max_center_x - min_center_x) >= (max_center_y - min_center_y);
  std::stable_sort(entries.begin(),
                   entries.end(),
                   [sort_by_x](const NetOrderEntry& lhs, const NetOrderEntry& rhs) {
                     if (sort_by_x) {
                       if (lhs.features.center_x != rhs.features.center_x) {
                         return lhs.features.center_x < rhs.features.center_x;
                       }
                     } else {
                       if (lhs.features.center_y != rhs.features.center_y) {
                         return lhs.features.center_y < rhs.features.center_y;
                       }
                     }
                     if (lhs.features.score != rhs.features.score) {
                       return lhs.features.score > rhs.features.score;
                     }
                     if (lhs.features.hpwl != rhs.features.hpwl) {
                       return lhs.features.hpwl > rhs.features.hpwl;
                     }
                     return lhs.net < rhs.net;
                   });

  // SPRoute-style deterministic round-robin binning: nearby nets are spread
  // across bins to reduce direct contention during reroute.
  const int bin_count = std::clamp(static_cast<int>(entries.size() / 420), 6, 22);
  std::vector<std::vector<NetOrderEntry>> bins(bin_count);
  for (size_t i = 0; i < entries.size(); ++i) {
    bins[i % bin_count].push_back(entries[i]);
  }

  for (int b = 0; b < bin_count; ++b) {
    std::stable_sort(
        bins[b].begin(),
        bins[b].end(),
        [b](const NetOrderEntry& lhs, const NetOrderEntry& rhs) {
          if (lhs.features.score != rhs.features.score) {
            return lhs.features.score > rhs.features.score;
          }
          if ((b % 2) == 0) {
            if (lhs.features.center_y != rhs.features.center_y) {
              return lhs.features.center_y < rhs.features.center_y;
            }
          } else {
            if (lhs.features.center_y != rhs.features.center_y) {
              return lhs.features.center_y > rhs.features.center_y;
            }
          }
          return lhs.net < rhs.net;
        });
  }

  std::vector<size_t> idx(bin_count, 0);
  nets.clear();
  nets.reserve(entries.size());
  size_t emitted = 0;
  while (emitted < entries.size()) {
    bool progressed = false;
    for (int b = 0; b < bin_count; ++b) {
      if (idx[b] >= bins[b].size()) {
        continue;
      }
      nets.push_back(bins[b][idx[b]].net);
      idx[b]++;
      emitted++;
      progressed = true;
    }
    if (!progressed) {
      break;
    }
  }
}

void reorderNetsByBspThenHpwlBurst(std::vector<Net*>& nets)
{
  if (nets.size() < 2) {
    return;
  }

  reorderNetsByBspScheduler(nets);

  const size_t burst_count
      = std::min(nets.size(),
                 static_cast<size_t>(
                     std::clamp<int>(nets.size() / 5, 72, 1800)));
  std::stable_sort(
      nets.begin(), nets.begin() + burst_count, [](Net* lhs, Net* rhs) {
        const NetOrderFeatures lhs_features = getNetOrderFeatures(lhs);
        const NetOrderFeatures rhs_features = getNetOrderFeatures(rhs);
        if (lhs_features.score != rhs_features.score) {
          return lhs_features.score > rhs_features.score;
        }
        if (lhs_features.hpwl != rhs_features.hpwl) {
          return lhs_features.hpwl > rhs_features.hpwl;
        }
        return lhs < rhs;
      });
}

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

void adjustEdgeCapacity(GlobalRouter* grouter,
                        int x1,
                        int y1,
                        int x2,
                        int y2,
                        int layer,
                        float ratio)
{
  ratio = std::clamp(ratio, 0.05f, 1.35f);
  FastRouteCore* core = grouter->fastroute();
  if (core == nullptr) {
    return;
  }
  const int current_cap = core->getEdgeCapacity(x1, y1, x2, y2, layer);
  if (current_cap <= 0) {
    return;
  }
  const int new_cap
      = std::max(1, static_cast<int>(std::floor(current_cap * ratio)));
  if (new_cap == current_cap) {
    return;
  }
  const bool is_reduce = new_cap < current_cap;
  core->addAdjustment(x1, y1, x2, y2, layer, new_cap, is_reduce);
}

void applyUniformCapacityBoost(GlobalRouter* grouter,
                               int min_layer,
                               int max_layer,
                               float ratio)
{
  Grid* grid = grouter->grid();
  if (grid == nullptr) {
    return;
  }
  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  if (x_grids <= 1 || y_grids <= 1) {
    return;
  }

  ratio = std::clamp(ratio, 1.0f, 1.35f);
  for (int layer = min_layer; layer <= max_layer; ++layer) {
    for (int y = 0; y < y_grids; ++y) {
      for (int x = 0; x < x_grids - 1; ++x) {
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }
    for (int y = 0; y < y_grids - 1; ++y) {
      for (int x = 0; x < x_grids; ++x) {
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

void applySoftCapacityScaling(GlobalRouter* grouter,
                              const RudyGrid& normalized_rudy,
                              int min_layer,
                              int max_layer,
                              float min_ratio_base = 0.50f,
                              float max_ratio_base = 0.92f,
                              float slope = 6.0f,
                              float midpoint = 0.45f)
{
  Grid* grid = grouter->grid();
  if (normalized_rudy.empty() || grid == nullptr) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  const int usable_x = std::min(x_grids, x_tiles);
  const int usable_y = std::min(y_grids, y_tiles);
  if (usable_x == 0 || usable_y == 0) {
    return;
  }

  const int layer_span = std::max(max_layer - min_layer, 1);
  const auto logistic_ratio = [](float normalized,
                                 float slope,
                                 float midpoint,
                                 float min_ratio,
                                 float max_ratio) {
    normalized = std::clamp(normalized, 0.0f, 1.0f);
    const float exponent = -slope * (normalized - midpoint);
    const float logistic = 1.0f / (1.0f + std::exp(exponent));
    const float blend = min_ratio + (max_ratio - min_ratio) * logistic;
    return std::clamp(blend, 0.05f, 0.99f);
  };

  const auto getNormalized = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float min_ratio = std::clamp(
        min_ratio_base + 0.12f * layer_factor, 0.05f, 0.99f);
    const float max_ratio = std::clamp(
        max_ratio_base + 0.04f * layer_factor, min_ratio, 0.995f);

    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x + 1, y));
        const float ratio
            = logistic_ratio(normalized, slope, midpoint, min_ratio, max_ratio);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float normalized = 0.5f * (getNormalized(x, y) + getNormalized(x, y + 1));
        const float ratio
            = logistic_ratio(normalized, slope, midpoint, min_ratio, max_ratio);
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

void applyHybridCapacityRemap(GlobalRouter* grouter,
                              const RudyGrid& normalized_rudy,
                              const std::vector<Hotspot>& hotspots,
                              int min_layer,
                              int max_layer,
                              float min_ratio_base = 0.48f,
                              float max_ratio_base = 0.96f,
                              float slope = 6.0f,
                              float midpoint = 0.42f,
                              float cool_boost = 0.45f,
                              int hotspot_halo = 2)
{
  Grid* grid = grouter->grid();
  if (normalized_rudy.empty() || grid == nullptr) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  const int x_tiles = normalized_rudy.size();
  const int y_tiles = normalized_rudy.front().size();
  const int usable_x = std::min(x_grids, x_tiles);
  const int usable_y = std::min(y_grids, y_tiles);
  if (usable_x == 0 || usable_y == 0) {
    return;
  }

  std::vector<std::vector<float>> hotspot_map(
      usable_x, std::vector<float>(usable_y, 0.0f));
  hotspot_halo = std::max(hotspot_halo, 0);
  for (const Hotspot& hotspot : hotspots) {
    for (int dx = -hotspot_halo; dx <= hotspot_halo; ++dx) {
      for (int dy = -hotspot_halo; dy <= hotspot_halo; ++dy) {
        const int x = hotspot.gx + dx;
        const int y = hotspot.gy + dy;
        if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
          continue;
        }
        const int manhattan = std::abs(dx) + std::abs(dy);
        const float decay = 1.0f / static_cast<float>(1 + manhattan);
        const float severity
            = std::max(0.0f, (hotspot.severity - 1.0f) * decay);
        hotspot_map[x][y] = std::max(hotspot_map[x][y], severity);
      }
    }
  }

  const int layer_span = std::max(max_layer - min_layer, 1);
  const auto logistic_ratio = [](float normalized,
                                 float slope,
                                 float midpoint,
                                 float min_ratio,
                                 float max_ratio) {
    normalized = std::clamp(normalized, 0.0f, 1.0f);
    const float exponent = -slope * (normalized - midpoint);
    const float logistic = 1.0f / (1.0f + std::exp(exponent));
    const float blend = min_ratio + (max_ratio - min_ratio) * logistic;
    return std::clamp(blend, 0.10f, 1.25f);
  };

  const auto getRudy = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return normalized_rudy[x][y];
  };
  const auto getHotspot = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= usable_x || y >= usable_y) {
      return 0.0f;
    }
    return hotspot_map[x][y];
  };

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float min_ratio = std::clamp(
        min_ratio_base + 0.08f * layer_factor, 0.10f, 1.10f);
    const float max_ratio = std::clamp(
        max_ratio_base + 0.06f * layer_factor, min_ratio, 1.25f);
    const float layer_cool_boost = cool_boost * (1.0f - 0.30f * layer_factor);

    for (int y = 0; y < usable_y; ++y) {
      for (int x = 0; x < usable_x - 1; ++x) {
        const float norm = 0.5f * (getRudy(x, y) + getRudy(x + 1, y));
        const float hotspot = 0.5f * (getHotspot(x, y) + getHotspot(x + 1, y));
        const float soft_ratio
            = logistic_ratio(norm, slope, midpoint, min_ratio, max_ratio);
        const float coolness = std::max(0.0f, 0.55f - norm);
        const float boost = 1.0f + layer_cool_boost * coolness;
        const float hotspot_shrink
            = std::clamp(1.0f - 0.35f * hotspot, 0.55f, 1.0f);
        const float ratio = std::clamp(soft_ratio * boost * hotspot_shrink,
                                       0.20f,
                                       1.30f);
        adjustEdgeCapacity(grouter, x, y, x + 1, y, layer, ratio);
      }
    }

    for (int y = 0; y < usable_y - 1; ++y) {
      for (int x = 0; x < usable_x; ++x) {
        const float norm = 0.5f * (getRudy(x, y) + getRudy(x, y + 1));
        const float hotspot = 0.5f * (getHotspot(x, y) + getHotspot(x, y + 1));
        const float soft_ratio
            = logistic_ratio(norm, slope, midpoint, min_ratio, max_ratio);
        const float coolness = std::max(0.0f, 0.55f - norm);
        const float boost = 1.0f + layer_cool_boost * coolness;
        const float hotspot_shrink
            = std::clamp(1.0f - 0.35f * hotspot, 0.55f, 1.0f);
        const float ratio = std::clamp(soft_ratio * boost * hotspot_shrink,
                                       0.20f,
                                       1.30f);
        adjustEdgeCapacity(grouter, x, y, x, y + 1, layer, ratio);
      }
    }
  }
}

void applyHotspotPenalties(GlobalRouter* grouter,
                           const std::vector<Hotspot>& hotspots,
                           int min_layer,
                           int max_layer,
                           int halo,
                           float base_ratio,
                           float severity_weight = 0.5f)
{
  Grid* grid = grouter->grid();
  if (hotspots.empty() || grid == nullptr) {
    return;
  }

  const int x_grids = grid->getXGrids();
  const int y_grids = grid->getYGrids();
  const int layer_span = std::max(max_layer - min_layer, 1);
  halo = std::max(0, halo);

  for (int layer = min_layer; layer <= max_layer; ++layer) {
    const float layer_factor
        = static_cast<float>(layer - min_layer) / static_cast<float>(layer_span);
    const float layer_ratio
        = std::clamp(base_ratio + 0.1f * layer_factor, 0.4f, 0.95f);
    const float scaled_severity_weight = std::clamp(severity_weight, 0.0f, 1.0f);

    for (const Hotspot& hotspot : hotspots) {
      for (int dx = -halo; dx <= halo; ++dx) {
        for (int dy = -halo; dy <= halo; ++dy) {
          const int gx = hotspot.gx + dx;
          const int gy = hotspot.gy + dy;
          if (gx < 0 || gy < 0 || gx >= x_grids || gy >= y_grids) {
            continue;
          }
          const float ratio_scale = std::clamp(
              1.0f - scaled_severity_weight * (hotspot.severity - 1.0f), 0.5f, 1.5f);
          const float adjusted_ratio
              = std::clamp(layer_ratio * ratio_scale, 0.25f, 0.98f);
          if (hotspot.affect_horizontal && gx < x_grids - 1) {
            adjustEdgeCapacity(
                grouter, gx, gy, gx + 1, gy, layer, adjusted_ratio);
          }
          if (hotspot.affect_vertical && gy < y_grids - 1) {
            adjustEdgeCapacity(
                grouter, gx, gy, gx, gy + 1, layer, adjusted_ratio);
          }
        }
      }
    }
  }
}

void appendUniqueRouteSegment(GRoute& route, const GSegment& segment)
{
  if (segment.init_x == segment.final_x && segment.init_y == segment.final_y
      && segment.init_layer == segment.final_layer) {
    return;
  }
  if (std::find(route.begin(), route.end(), segment) != route.end()) {
    return;
  }
  route.push_back(segment);
}

std::int64_t makeGridKey(int gx, int gy)
{
  return (static_cast<std::int64_t>(gx) << 32)
         | static_cast<std::uint32_t>(gy);
}

std::map<std::int64_t, float> buildHotspotSeverityMap(const std::vector<Hotspot>& hotspots)
{
  std::map<std::int64_t, float> severity_map;
  for (const Hotspot& hotspot : hotspots) {
    const std::int64_t key = makeGridKey(hotspot.gx, hotspot.gy);
    auto [it, inserted] = severity_map.emplace(key, hotspot.severity);
    if (!inserted) {
      it->second = std::max(it->second, hotspot.severity);
    }
  }
  return severity_map;
}

NetRouteCost computeNetRouteCost(const GRoute& route,
                                 int tile_size,
                                 int x_min,
                                 int y_min,
                                 int x_grids,
                                 int y_grids,
                                 const std::map<std::int64_t, float>& hotspot_map)
{
  NetRouteCost cost;
  tile_size = std::max(tile_size, 1);
  x_grids = std::max(x_grids, 1);
  y_grids = std::max(y_grids, 1);

  const auto to_grid = [&](int x, int y) {
    const int gx = std::clamp((x - x_min) / tile_size, 0, x_grids - 1);
    const int gy = std::clamp((y - y_min) / tile_size, 0, y_grids - 1);
    return std::pair<int, int>{gx, gy};
  };
  const auto sample_hotspot = [&](int x, int y) {
    const auto [gx, gy] = to_grid(x, y);
    auto it = hotspot_map.find(makeGridKey(gx, gy));
    if (it != hotspot_map.end()) {
      cost.hotspot_exposure += static_cast<double>(it->second);
    }
  };

  for (const GSegment& segment : route) {
    if (segment.isVia()) {
      cost.via_count++;
      sample_hotspot(segment.init_x, segment.init_y);
      continue;
    }

    const long seg_length = static_cast<long>(std::abs(segment.final_x - segment.init_x)
                                              + std::abs(segment.final_y - segment.init_y));
    cost.wirelength_dbu += seg_length;

    sample_hotspot(segment.init_x, segment.init_y);
    sample_hotspot(segment.final_x, segment.final_y);
    sample_hotspot((segment.init_x + segment.final_x) / 2,
                   (segment.init_y + segment.final_y) / 2);
  }

  const double via_weight = static_cast<double>(tile_size) * 2.8;
  const double hotspot_weight = static_cast<double>(tile_size) * 4.0;
  cost.score = static_cast<double>(cost.wirelength_dbu)
               + via_weight * static_cast<double>(cost.via_count)
               + hotspot_weight * cost.hotspot_exposure;
  return cost;
}

int applySelectiveNetRouteGrafting(NetRouteMap& base_routes,
                                   const NetRouteMap& alternate_routes,
                                   int tile_size,
                                   int x_min,
                                   int y_min,
                                   int x_grids,
                                   int y_grids,
                                   const std::map<std::int64_t, float>& hotspot_map)
{
  int replaced_nets = 0;
  const long min_wl_gain = std::max(2 * tile_size, 1);

  for (auto& [db_net, base_route] : base_routes) {
    auto alt_it = alternate_routes.find(db_net);
    if (alt_it == alternate_routes.end()) {
      continue;
    }

    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const NetRouteCost alt_cost = computeNetRouteCost(alt_it->second,
                                                      tile_size,
                                                      x_min,
                                                      y_min,
                                                      x_grids,
                                                      y_grids,
                                                      hotspot_map);

    // Accept alternate routes only when they strictly improve local length
    // without worsening via count and while preserving hotspot pressure.
    const bool wl_better
        = (base_cost.wirelength_dbu - alt_cost.wirelength_dbu) >= min_wl_gain;
    const bool via_not_worse = alt_cost.via_count <= base_cost.via_count;
    const bool hotspot_not_worse
        = alt_cost.hotspot_exposure <= (base_cost.hotspot_exposure + 0.20);
    const bool score_better = alt_cost.score + static_cast<double>(tile_size)
                              < base_cost.score;

    if (wl_better && via_not_worse && hotspot_not_worse && score_better) {
      base_route = alt_it->second;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

double computeWirelengthDominantObjective(const NetRouteCost& cost, int tile_size)
{
  const double via_weight = static_cast<double>(tile_size) * 1.65;
  const double hotspot_weight = static_cast<double>(tile_size) * 2.10;
  return static_cast<double>(cost.wirelength_dbu)
         + via_weight * static_cast<double>(cost.via_count)
         + hotspot_weight * cost.hotspot_exposure;
}

int applyMultiScenarioWirelengthFusion(
    NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    int tile_size,
    int x_min,
    int y_min,
    int x_grids,
    int y_grids,
    const std::map<std::int64_t, float>& hotspot_map,
    long min_wl_gain,
    int via_slack,
    double hotspot_slack)
{
  int replaced_nets = 0;
  min_wl_gain = std::max(min_wl_gain, 1L);
  via_slack = std::max(via_slack, 0);
  hotspot_slack = std::max(hotspot_slack, 0.0);

  for (auto& [db_net, base_route] : base_routes) {
    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const GRoute* chosen_route = &base_route;
    NetRouteCost chosen_cost = base_cost;
    double chosen_objective
        = computeWirelengthDominantObjective(base_cost, tile_size);

    for (const NetRouteMap* donor_routes : donor_route_sets) {
      if (donor_routes == nullptr) {
        continue;
      }

      const auto donor_it = donor_routes->find(db_net);
      if (donor_it == donor_routes->end()) {
        continue;
      }

      const NetRouteCost donor_cost = computeNetRouteCost(donor_it->second,
                                                          tile_size,
                                                          x_min,
                                                          y_min,
                                                          x_grids,
                                                          y_grids,
                                                          hotspot_map);
      if (donor_cost.via_count > (base_cost.via_count + via_slack)) {
        continue;
      }
      if (donor_cost.hotspot_exposure > (base_cost.hotspot_exposure + hotspot_slack)) {
        continue;
      }

      const double donor_objective
          = computeWirelengthDominantObjective(donor_cost, tile_size);
      const bool wl_strictly_better
          = donor_cost.wirelength_dbu
            < (chosen_cost.wirelength_dbu - min_wl_gain);
      const bool wl_tie_or_better
          = donor_cost.wirelength_dbu <= chosen_cost.wirelength_dbu;
      const bool objective_better
          = donor_objective + (0.30 * static_cast<double>(tile_size))
            < chosen_objective;

      if (!wl_strictly_better && !(wl_tie_or_better && objective_better)) {
        continue;
      }

      chosen_route = &donor_it->second;
      chosen_cost = donor_cost;
      chosen_objective = donor_objective;
    }

    if (chosen_route != &base_route) {
      base_route = *chosen_route;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

int applyAdaptiveWirelengthFusion(
    NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    int tile_size,
    int x_min,
    int y_min,
    int x_grids,
    int y_grids,
    const std::map<std::int64_t, float>& hotspot_map,
    long min_wl_gain,
    int base_via_slack,
    double base_hotspot_slack)
{
  int replaced_nets = 0;
  min_wl_gain = std::max(min_wl_gain, 1L);
  base_via_slack = std::max(base_via_slack, 0);
  base_hotspot_slack = std::max(base_hotspot_slack, 0.0);
  tile_size = std::max(tile_size, 1);

  const auto compute_fusion_objective = [tile_size](const NetRouteCost& cost) {
    // Strongly prioritize WL. Allow moderate via increase when donor routes
    // give meaningful Manhattan path shortening.
    const double via_weight = static_cast<double>(tile_size) * 0.85;
    const double hotspot_weight = static_cast<double>(tile_size) * 1.40;
    return static_cast<double>(cost.wirelength_dbu)
           + via_weight * static_cast<double>(cost.via_count)
           + hotspot_weight * cost.hotspot_exposure;
  };

  for (auto& [db_net, base_route] : base_routes) {
    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const GRoute* chosen_route = &base_route;
    NetRouteCost chosen_cost = base_cost;
    double chosen_objective = compute_fusion_objective(base_cost);

    for (const NetRouteMap* donor_routes : donor_route_sets) {
      if (donor_routes == nullptr) {
        continue;
      }
      const auto donor_it = donor_routes->find(db_net);
      if (donor_it == donor_routes->end()) {
        continue;
      }

      const NetRouteCost donor_cost = computeNetRouteCost(donor_it->second,
                                                          tile_size,
                                                          x_min,
                                                          y_min,
                                                          x_grids,
                                                          y_grids,
                                                          hotspot_map);
      const long wl_gain = base_cost.wirelength_dbu - donor_cost.wirelength_dbu;
      if (wl_gain < min_wl_gain) {
        continue;
      }

      const int wl_bonus_via_slack = static_cast<int>(
          wl_gain / std::max(3 * tile_size, 1));
      const int allowed_vias
          = base_cost.via_count + base_via_slack + wl_bonus_via_slack;
      if (donor_cost.via_count > allowed_vias) {
        continue;
      }

      const double wl_bonus_hotspot
          = static_cast<double>(wl_gain)
            / static_cast<double>(std::max(8 * tile_size, 1));
      const double allowed_hotspot
          = base_cost.hotspot_exposure + base_hotspot_slack + wl_bonus_hotspot;
      if (donor_cost.hotspot_exposure > allowed_hotspot) {
        continue;
      }

      const double donor_objective = compute_fusion_objective(donor_cost);
      const bool wl_better_than_chosen
          = donor_cost.wirelength_dbu
            < (chosen_cost.wirelength_dbu - std::max<long>(1, min_wl_gain / 2));
      const bool objective_better
          = donor_objective + (0.15 * static_cast<double>(tile_size))
            < chosen_objective;
      if (!wl_better_than_chosen && !objective_better) {
        continue;
      }

      chosen_route = &donor_it->second;
      chosen_cost = donor_cost;
      chosen_objective = donor_objective;
    }

    if (chosen_route != &base_route) {
      base_route = *chosen_route;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

int applyExtremeWirelengthFusion(
    NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    int tile_size,
    int x_min,
    int y_min,
    int x_grids,
    int y_grids,
    const std::map<std::int64_t, float>& hotspot_map,
    long min_wl_gain,
    double min_wl_gain_ratio,
    int base_via_slack,
    double via_gain_scale,
    double base_hotspot_slack)
{
  int replaced_nets = 0;
  tile_size = std::max(tile_size, 1);
  min_wl_gain = std::max(min_wl_gain, 1L);
  min_wl_gain_ratio = std::clamp(min_wl_gain_ratio, 0.0, 1.0);
  base_via_slack = std::max(base_via_slack, 0);
  via_gain_scale = std::max(via_gain_scale, 0.4);
  base_hotspot_slack = std::max(base_hotspot_slack, 0.0);

  const auto compute_extreme_objective = [tile_size](const NetRouteCost& cost) {
    // Extreme mode: keep WL as the dominant signal but still discourage
    // pathological donor swaps with excessive vias/hotspot exposure.
    const double via_weight = static_cast<double>(tile_size) * 0.30;
    const double hotspot_weight = static_cast<double>(tile_size) * 0.50;
    return static_cast<double>(cost.wirelength_dbu)
           + via_weight * static_cast<double>(cost.via_count)
           + hotspot_weight * cost.hotspot_exposure;
  };

  for (auto& [db_net, base_route] : base_routes) {
    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const GRoute* chosen_route = &base_route;
    NetRouteCost chosen_cost = base_cost;
    double chosen_objective = compute_extreme_objective(base_cost);

    for (const NetRouteMap* donor_routes : donor_route_sets) {
      if (donor_routes == nullptr) {
        continue;
      }
      const auto donor_it = donor_routes->find(db_net);
      if (donor_it == donor_routes->end()) {
        continue;
      }

      const NetRouteCost donor_cost = computeNetRouteCost(donor_it->second,
                                                          tile_size,
                                                          x_min,
                                                          y_min,
                                                          x_grids,
                                                          y_grids,
                                                          hotspot_map);
      const long wl_gain_from_base
          = base_cost.wirelength_dbu - donor_cost.wirelength_dbu;
      if (wl_gain_from_base < min_wl_gain) {
        continue;
      }

      const double wl_gain_ratio
          = static_cast<double>(wl_gain_from_base)
            / static_cast<double>(std::max(base_cost.wirelength_dbu, 1L));
      const bool has_meaningful_gain
          = wl_gain_ratio >= min_wl_gain_ratio
            || wl_gain_from_base >= (2 * min_wl_gain);
      if (!has_meaningful_gain) {
        continue;
      }

      const int wl_bonus_via_slack = static_cast<int>(
          static_cast<double>(wl_gain_from_base)
          / static_cast<double>(std::max(
              static_cast<int>(std::round(via_gain_scale * tile_size)), 1)));
      const int large_gain_bonus = (wl_gain_from_base >= (12 * tile_size)) ? 2 : 0;
      const int allowed_vias
          = base_cost.via_count + base_via_slack + wl_bonus_via_slack + large_gain_bonus;
      if (donor_cost.via_count > allowed_vias) {
        continue;
      }

      const double wl_bonus_hotspot
          = static_cast<double>(wl_gain_from_base)
            / static_cast<double>(std::max(10 * tile_size, 1));
      const double allowed_hotspot
          = base_cost.hotspot_exposure + base_hotspot_slack + wl_bonus_hotspot;
      if (donor_cost.hotspot_exposure > allowed_hotspot) {
        continue;
      }

      const double donor_objective = compute_extreme_objective(donor_cost);
      const bool wl_strictly_better
          = donor_cost.wirelength_dbu
            < (chosen_cost.wirelength_dbu - std::max<long>(1, min_wl_gain / 2));
      const bool wl_not_worse
          = donor_cost.wirelength_dbu <= chosen_cost.wirelength_dbu;
      const bool objective_better
          = donor_objective + (0.05 * static_cast<double>(tile_size))
            < chosen_objective;
      if (!wl_strictly_better && !(wl_not_worse && objective_better)) {
        continue;
      }

      chosen_route = &donor_it->second;
      chosen_cost = donor_cost;
      chosen_objective = donor_objective;
    }

    if (chosen_route != &base_route) {
      base_route = *chosen_route;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

int applyLongNetPriorityFusion(
    NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    int tile_size,
    int x_min,
    int y_min,
    int x_grids,
    int y_grids,
    const std::map<std::int64_t, float>& hotspot_map,
    long min_wl_gain,
    int base_via_slack,
    double base_hotspot_slack)
{
  int replaced_nets = 0;
  if (base_routes.empty()) {
    return replaced_nets;
  }

  tile_size = std::max(tile_size, 1);
  min_wl_gain = std::max(min_wl_gain, 1L);
  base_via_slack = std::max(base_via_slack, 0);
  base_hotspot_slack = std::max(base_hotspot_slack, 0.0);

  std::vector<long> route_lengths;
  route_lengths.reserve(base_routes.size());
  for (const auto& [db_net, route] : base_routes) {
    static_cast<void>(db_net);
    const NetRouteCost cost = computeNetRouteCost(route,
                                                  tile_size,
                                                  x_min,
                                                  y_min,
                                                  x_grids,
                                                  y_grids,
                                                  hotspot_map);
    route_lengths.push_back(std::max(cost.wirelength_dbu, 0L));
  }
  if (route_lengths.empty()) {
    return replaced_nets;
  }
  std::sort(route_lengths.begin(), route_lengths.end());

  const auto percentile_value = [&route_lengths](double q) {
    q = std::clamp(q, 0.0, 1.0);
    const size_t idx = static_cast<size_t>(
        std::floor(q * static_cast<double>(route_lengths.size() - 1)));
    return route_lengths[idx];
  };

  const long long_net_threshold = percentile_value(0.70);
  const long xlong_net_threshold = percentile_value(0.90);

  const auto compute_tiered_objective = [tile_size](const NetRouteCost& cost,
                                                    bool is_long,
                                                    bool is_xlong) {
    if (is_xlong) {
      const double via_weight = static_cast<double>(tile_size) * 0.14;
      const double hotspot_weight = static_cast<double>(tile_size) * 0.24;
      return static_cast<double>(cost.wirelength_dbu)
             + via_weight * static_cast<double>(cost.via_count)
             + hotspot_weight * cost.hotspot_exposure;
    }
    if (is_long) {
      const double via_weight = static_cast<double>(tile_size) * 0.22;
      const double hotspot_weight = static_cast<double>(tile_size) * 0.38;
      return static_cast<double>(cost.wirelength_dbu)
             + via_weight * static_cast<double>(cost.via_count)
             + hotspot_weight * cost.hotspot_exposure;
    }
    const double via_weight = static_cast<double>(tile_size) * 0.40;
    const double hotspot_weight = static_cast<double>(tile_size) * 0.65;
    return static_cast<double>(cost.wirelength_dbu)
           + via_weight * static_cast<double>(cost.via_count)
           + hotspot_weight * cost.hotspot_exposure;
  };

  for (auto& [db_net, base_route] : base_routes) {
    static_cast<void>(db_net);
    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const bool is_xlong = base_cost.wirelength_dbu >= xlong_net_threshold;
    const bool is_long = is_xlong || base_cost.wirelength_dbu >= long_net_threshold;
    const long tier_min_gain = is_xlong ? std::max(1L, min_wl_gain / 3)
                             : is_long ? std::max(1L, min_wl_gain / 2)
                                       : min_wl_gain;

    const GRoute* chosen_route = &base_route;
    NetRouteCost chosen_cost = base_cost;
    double chosen_objective
        = compute_tiered_objective(base_cost, is_long, is_xlong);

    for (const NetRouteMap* donor_routes : donor_route_sets) {
      if (donor_routes == nullptr) {
        continue;
      }
      const auto donor_it = donor_routes->find(db_net);
      if (donor_it == donor_routes->end()) {
        continue;
      }

      const NetRouteCost donor_cost = computeNetRouteCost(donor_it->second,
                                                          tile_size,
                                                          x_min,
                                                          y_min,
                                                          x_grids,
                                                          y_grids,
                                                          hotspot_map);

      const long wl_gain_from_base
          = base_cost.wirelength_dbu - donor_cost.wirelength_dbu;
      if (wl_gain_from_base < tier_min_gain) {
        continue;
      }

      const int via_scale = is_xlong ? 2 : (is_long ? 3 : 4);
      const int wl_bonus_via = static_cast<int>(
          wl_gain_from_base / std::max(via_scale * tile_size, 1));
      const int tier_via_bonus = is_xlong ? 10 : (is_long ? 5 : 1);
      const int allowed_vias = base_cost.via_count + base_via_slack
                               + tier_via_bonus + wl_bonus_via;
      if (donor_cost.via_count > allowed_vias) {
        continue;
      }

      const int hotspot_scale = is_xlong ? 8 : (is_long ? 10 : 12);
      const double wl_bonus_hotspot
          = static_cast<double>(wl_gain_from_base)
            / static_cast<double>(std::max(hotspot_scale * tile_size, 1));
      const double tier_hotspot_bonus = is_xlong ? 1.5 : (is_long ? 0.8 : 0.2);
      const double allowed_hotspot = base_cost.hotspot_exposure
                                     + base_hotspot_slack
                                     + tier_hotspot_bonus
                                     + wl_bonus_hotspot;
      if (donor_cost.hotspot_exposure > allowed_hotspot) {
        continue;
      }

      const double donor_objective
          = compute_tiered_objective(donor_cost, is_long, is_xlong);
      const bool wl_strictly_better
          = donor_cost.wirelength_dbu
            < (chosen_cost.wirelength_dbu - std::max<long>(1, tier_min_gain / 2));
      const bool wl_not_worse
          = donor_cost.wirelength_dbu <= chosen_cost.wirelength_dbu;
      const bool objective_better
          = donor_objective + (0.05 * static_cast<double>(tile_size))
            < chosen_objective;

      if (!wl_strictly_better && !(wl_not_worse && objective_better)) {
        continue;
      }

      chosen_route = &donor_it->second;
      chosen_cost = donor_cost;
      chosen_objective = donor_objective;
    }

    if (chosen_route != &base_route) {
      base_route = *chosen_route;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

int applyViaAwareStabilizationFusion(
    NetRouteMap& base_routes,
    const std::vector<const NetRouteMap*>& donor_route_sets,
    int tile_size,
    int x_min,
    int y_min,
    int x_grids,
    int y_grids,
    const std::map<std::int64_t, float>& hotspot_map)
{
  int replaced_nets = 0;
  tile_size = std::max(tile_size, 1);

  for (auto& [db_net, base_route] : base_routes) {
    const NetRouteCost base_cost = computeNetRouteCost(base_route,
                                                       tile_size,
                                                       x_min,
                                                       y_min,
                                                       x_grids,
                                                       y_grids,
                                                       hotspot_map);
    const GRoute* chosen_route = &base_route;
    double best_utility = 0.0;

    for (const NetRouteMap* donor_routes : donor_route_sets) {
      if (donor_routes == nullptr) {
        continue;
      }
      const auto donor_it = donor_routes->find(db_net);
      if (donor_it == donor_routes->end()) {
        continue;
      }

      const NetRouteCost donor_cost = computeNetRouteCost(donor_it->second,
                                                          tile_size,
                                                          x_min,
                                                          y_min,
                                                          x_grids,
                                                          y_grids,
                                                          hotspot_map);
      const long wl_delta
          = donor_cost.wirelength_dbu - base_cost.wirelength_dbu;
      const long wl_loss = std::max(0L, wl_delta);

      const long allowed_wl_loss = std::max<long>(
          tile_size,
          std::min<long>(base_cost.wirelength_dbu / 55, 3L * tile_size));
      if (wl_loss > allowed_wl_loss) {
        continue;
      }

      const long via_gain = base_cost.via_count - donor_cost.via_count;
      const double hotspot_gain
          = base_cost.hotspot_exposure - donor_cost.hotspot_exposure;
      const long extra_vias
          = std::max(0L, donor_cost.via_count - base_cost.via_count);

      // Prefer route simplification for detailed routing while keeping WL
      // almost unchanged. This pass trims donor swaps that are too "spiky".
      const double utility
          = static_cast<double>(via_gain) * (0.90 * static_cast<double>(tile_size))
            + hotspot_gain * (1.80 * static_cast<double>(tile_size))
            - static_cast<double>(wl_loss) * 1.90
            - static_cast<double>(extra_vias)
                  * (0.65 * static_cast<double>(tile_size));

      const bool improves_geometry = via_gain >= 3 || hotspot_gain >= 0.35;
      const bool has_wl_benefit = wl_delta < 0;
      if (!has_wl_benefit && !improves_geometry) {
        continue;
      }

      if (utility > best_utility) {
        best_utility = utility;
        chosen_route = &donor_it->second;
      }
    }

    if (chosen_route != &base_route) {
      base_route = *chosen_route;
      replaced_nets++;
    }
  }

  return replaced_nets;
}

void applyCugrStyleGuidePatching(GlobalRouter* grouter,
                                 NetRouteMap& routes,
                                 int min_routing_layer,
                                 int max_routing_layer,
                                 utl::Logger* logger)
{
  if (grouter == nullptr || routes.empty()) {
    return;
  }

  FastRouteCore* core = grouter->fastroute();
  if (core == nullptr || grouter->grid() == nullptr) {
    return;
  }
  const int tile_size = std::max(grouter->grid()->getTileSize(), 1);

  core->computeCongestionInformation();

  std::vector<CongestionInformation> vertical;
  std::vector<CongestionInformation> horizontal;
  core->getCongestionGrid(vertical, horizontal);

  struct CongestionPatchCandidate
  {
    const CongestionInformation* info = nullptr;
    float severity = 0.0f;
  };

  std::vector<CongestionPatchCandidate> candidates;
  candidates.reserve(vertical.size() + horizontal.size());

  const auto append_candidates = [&](const std::vector<CongestionInformation>& edges) {
    for (const CongestionInformation& info : edges) {
      if (info.sources.empty()) {
        continue;
      }

      const int capacity = std::max(info.congestion.capacity, 1);
      const int usage = std::max(info.congestion.usage, 0);
      const int overflow = std::max(usage - capacity, 0);
      const float usage_ratio
          = static_cast<float>(usage) / static_cast<float>(capacity);

      // Keep patching focused on severe regions. Patching low-pressure edges
      // can broaden guides and invite unnecessary detailed-route detours.
      if (overflow == 0 && usage_ratio < 0.90f) {
        continue;
      }

      const float severity
          = usage_ratio + 0.22f * static_cast<float>(overflow);
      candidates.push_back({&info, severity});
    }
  };

  append_candidates(horizontal);
  append_candidates(vertical);

  if (candidates.empty()) {
    return;
  }

  std::stable_sort(
      candidates.begin(),
      candidates.end(),
      [](const CongestionPatchCandidate& lhs, const CongestionPatchCandidate& rhs) {
        return lhs.severity > rhs.severity;
      });

  const int edge_budget = 64;
  const int sources_per_edge = 1;
  const int patches_per_net = 4;
  const int total_patch_budget = 780;

  int patched_edges = 0;
  int added_segments = 0;
  int congestion_added_segments = 0;
  std::map<odb::dbNet*, int> per_net_patch_count;

  for (const CongestionPatchCandidate& candidate : candidates) {
    if (patched_edges >= edge_budget || added_segments >= total_patch_budget) {
      break;
    }
    if (candidate.info == nullptr) {
      continue;
    }

    const GSegment& base = candidate.info->segment;
    const int base_layer = base.init_layer;
    if (base_layer < min_routing_layer || base_layer > max_routing_layer) {
      continue;
    }

    const int alt_layer = (base_layer < max_routing_layer) ? base_layer + 1
                                                            : base_layer - 1;
    if (alt_layer < min_routing_layer || alt_layer > max_routing_layer
        || alt_layer == base_layer) {
      continue;
    }

    bool edge_patched = false;
    int source_count = 0;
    const int low_layer = std::min(base_layer, alt_layer);
    const int high_layer = std::max(base_layer, alt_layer);

    for (odb::dbNet* db_net : candidate.info->sources) {
      if (source_count >= sources_per_edge
          || added_segments >= total_patch_budget) {
        break;
      }
      auto route_it = routes.find(db_net);
      if (route_it == routes.end()) {
        continue;
      }

      int& net_budget = per_net_patch_count[db_net];
      if (net_budget >= patches_per_net) {
        continue;
      }

      GRoute& route = route_it->second;
      const size_t route_size_before = route.size();

      appendUniqueRouteSegment(route,
                               GSegment(base.init_x,
                                        base.init_y,
                                        alt_layer,
                                        base.final_x,
                                        base.final_y,
                                        alt_layer));
      appendUniqueRouteSegment(route,
                               GSegment(base.init_x,
                                        base.init_y,
                                        low_layer,
                                        base.init_x,
                                        base.init_y,
                                        high_layer));
      appendUniqueRouteSegment(route,
                               GSegment(base.final_x,
                                        base.final_y,
                                        low_layer,
                                        base.final_x,
                                        base.final_y,
                                        high_layer));

      const int new_segments
          = static_cast<int>(route.size() - route_size_before);
      if (new_segments > 0) {
        added_segments += new_segments;
        congestion_added_segments += new_segments;
        net_budget++;
        source_count++;
        edge_patched = true;
      }
    }

    if (edge_patched) {
      patched_edges++;
    }
  }

  // Add long-segment alternatives only when there are enough severe hotspots
  // to justify extra guide flexibility.
  const bool enable_long_seg_patching = patched_edges >= (edge_budget / 3);
  const int long_seg_threshold = 4 * tile_size;
  const int long_seg_patches_per_net = 3;
  const int long_seg_segment_budget = 360;
  int long_seg_added_segments = 0;
  std::map<odb::dbNet*, int> long_seg_patch_count;

  if (enable_long_seg_patching) {
    for (auto& [db_net, route] : routes) {
      if (added_segments >= total_patch_budget
          || long_seg_added_segments >= long_seg_segment_budget
          ) {
        break;
      }
      int& net_budget = long_seg_patch_count[db_net];
      if (net_budget >= long_seg_patches_per_net) {
        continue;
      }

      const GRoute original_route = route;
      for (const GSegment& base : original_route) {
        if (added_segments >= total_patch_budget
            || long_seg_added_segments >= long_seg_segment_budget
            || net_budget >= long_seg_patches_per_net) {
          break;
        }
        if (base.isVia() || base.length() < long_seg_threshold) {
          continue;
        }
        if (base.init_layer != base.final_layer) {
          continue;
        }

        const int base_layer = base.init_layer;
        if (base_layer < min_routing_layer || base_layer > max_routing_layer) {
          continue;
        }
        const int alt_layer = (base_layer < max_routing_layer) ? base_layer + 1
                                                                : base_layer - 1;
        if (alt_layer < min_routing_layer || alt_layer > max_routing_layer
            || alt_layer == base_layer) {
          continue;
        }

        const int low_layer = std::min(base_layer, alt_layer);
        const int high_layer = std::max(base_layer, alt_layer);
        const size_t route_size_before = route.size();

        appendUniqueRouteSegment(route,
                                 GSegment(base.init_x,
                                          base.init_y,
                                          alt_layer,
                                          base.final_x,
                                          base.final_y,
                                          alt_layer));
        appendUniqueRouteSegment(route,
                                 GSegment(base.init_x,
                                          base.init_y,
                                          low_layer,
                                          base.init_x,
                                          base.init_y,
                                          high_layer));
        appendUniqueRouteSegment(route,
                                 GSegment(base.final_x,
                                          base.final_y,
                                          low_layer,
                                          base.final_x,
                                          base.final_y,
                                          high_layer));

        const int new_segments
            = static_cast<int>(route.size() - route_size_before);
        if (new_segments > 0) {
          added_segments += new_segments;
          long_seg_added_segments += new_segments;
          net_budget++;
        }
      }
    }
  }

  if (logger != nullptr && added_segments > 0) {
    logger->info(GNR,
                 6009,
                 "NEWGR patching: added {} guide segments (congestion {}, long-segment {}) across {} hotspots",
                 added_segments,
                 congestion_added_segments,
                 long_seg_added_segments,
                 patched_edges);
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

  auto compute_metrics = [&](const NetRouteMap& routes) -> RouteMetrics {
    RouteMetrics metrics;
    for (const auto& [db_net, segments] : routes) {
      static_cast<void>(db_net);
      for (const GSegment& segment : segments) {
        if (segment.isVia()) {
          metrics.via_count++;
        } else {
          metrics.wirelength_dbu
              += std::abs(segment.final_x - segment.init_x)
                 + std::abs(segment.final_y - segment.init_y);
        }
      }
    }

    if (metrics.wirelength_dbu > 0 && grouter_->db_ != nullptr
        && grouter_->db_->getTech() != nullptr) {
      metrics.wirelength_um
          = metrics.wirelength_dbu
            / static_cast<double>(
                grouter_->db_->getTech()->getDbUnitsPerMicron());
    }

    const double via_weight
        = static_cast<double>(std::max(grouter_->grid_->getTileSize(), 1))
          * 3.0;
    metrics.score = static_cast<double>(metrics.wirelength_dbu)
                    + via_weight * static_cast<double>(metrics.via_count);
    return metrics;
  };

  auto capture_snapshot = [&]() -> RouterSnapshot {
    RouterSnapshot snapshot;
    snapshot.caps_percentage = grouter_->caps_perturbation_percentage_;
    snapshot.perturbation_amount = grouter_->perturbation_amount_;
    snapshot.critical_percentage
        = grouter_->fastroute_->getCriticalNetsPercentage();
    snapshot.allow_congestion = grouter_->allow_congestion_;
    snapshot.seed = grouter_->seed_;
    return snapshot;
  };

  auto restore_snapshot = [&](const RouterSnapshot& snapshot) {
    grouter_->setCapacitiesPerturbationPercentage(snapshot.caps_percentage);
    grouter_->setPerturbationAmount(snapshot.perturbation_amount);
    grouter_->setAllowCongestion(snapshot.allow_congestion);
    grouter_->setSeed(snapshot.seed);
    grouter_->fastroute_->setCriticalNetsPercentage(
        snapshot.critical_percentage);
  };

  auto run_existing_state = [&](const std::string& name,
                                std::vector<Net*>& state_nets) {
    NetRouteMap routes;
    if (!state_nets.empty()) {
      routes = grouter_->findRouting(
          state_nets, min_routing_layer, max_routing_layer);
    }
    RouteMetrics metrics = compute_metrics(routes);
    logger_->info(GNR,
                  6005,
                  "NEWGR {}: wirelength {:.0f} um, vias {}",
                  name,
                  metrics.wirelength_um,
                  metrics.via_count);
    return ScenarioResult{name, metrics, std::move(routes)};
  };

  auto collect_hotspots = [&]() -> std::vector<Hotspot> {
    std::vector<Hotspot> hotspots;
    if (grouter_->fastroute_ == nullptr || grouter_->grid_ == nullptr) {
      return hotspots;
    }

    grouter_->fastroute_->computeCongestionInformation();

    std::vector<CongestionInformation> vertical;
    std::vector<CongestionInformation> horizontal;
    grouter_->fastroute_->getCongestionGrid(vertical, horizontal);

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    auto append_hotspots = [&](const std::vector<CongestionInformation>& edges,
                               bool is_vertical) {
      for (const auto& info : edges) {
        const int capacity = std::max(info.congestion.capacity, 1);
        const float usage_ratio
            = static_cast<float>(info.congestion.usage)
              / static_cast<float>(capacity);
        if (usage_ratio < 0.6f) {
          continue;
        }
        const int gx
            = std::clamp((info.segment.init_x - x_min) / tile_size, 0, x_grids);
        const int gy
            = std::clamp((info.segment.init_y - y_min) / tile_size, 0, y_grids);
        Hotspot hotspot;
        hotspot.gx = std::clamp(gx, 0, std::max(x_grids - 1, 0));
        hotspot.gy = std::clamp(gy, 0, std::max(y_grids - 1, 0));
        hotspot.severity = std::clamp(usage_ratio, 0.6f, 3.0f);
        hotspot.affect_vertical = is_vertical;
        hotspot.affect_horizontal = !is_vertical;
        hotspots.push_back(hotspot);
      }
    };

    append_hotspots(horizontal, false);
    append_hotspots(vertical, true);
    return hotspots;
  };

  auto focus_hotspots = [](const std::vector<Hotspot>& all_hotspots,
                           size_t max_count,
                           float min_severity) {
    std::vector<Hotspot> filtered;
    filtered.reserve(all_hotspots.size());
    for (const Hotspot& hotspot : all_hotspots) {
      if (hotspot.severity >= min_severity) {
        filtered.push_back(hotspot);
      }
    }

    std::stable_sort(filtered.begin(),
                     filtered.end(),
                     [](const Hotspot& lhs, const Hotspot& rhs) {
                       return lhs.severity > rhs.severity;
                     });
    if (filtered.size() > max_count) {
      filtered.resize(max_count);
    }
    return filtered;
  };

  auto run_scenario = [&](const ScenarioDefinition& scenario,
                          const RouterSnapshot& snapshot) {
    restore_snapshot(snapshot);
    if (scenario.pre_init) {
      scenario.pre_init();
    }

    std::vector<Net*> scenario_nets
        = grouter_->initFastRoute(min_routing_layer, max_routing_layer);
    if (scenario.order_nets) {
      scenario.order_nets(scenario_nets);
    }
    if (scenario.post_init) {
      scenario.post_init();
    }

    NetRouteMap routes;
    if (!scenario_nets.empty()) {
      routes = grouter_->findRouting(
          scenario_nets, min_routing_layer, max_routing_layer);
    }
    RouteMetrics metrics = compute_metrics(routes);
    logger_->info(GNR,
                  6006,
                  "NEWGR scenario {} [{}]: wirelength {:.0f} um, vias {}",
                  scenario.name,
                  scenario.aggressive ? "aggressive" : "conservative",
                  metrics.wirelength_um,
                  metrics.via_count);
    return ScenarioResult{scenario.name, metrics, std::move(routes)};
  };

  RouterSnapshot snapshot = capture_snapshot();

  ScenarioResult baseline
      = run_existing_state("baseline", nets);
  std::vector<Hotspot> hotspots = collect_hotspots();
  const std::vector<Hotspot> focused_hotspots
      = focus_hotspots(hotspots, hotspots.size(), 0.6f);
  const std::map<std::int64_t, float> hotspot_map
      = buildHotspotSeverityMap(focused_hotspots);

  RudyGrid normalized_rudy;
  if (Rudy* rudy = grouter_->getRudy()) {
    rudy->calculateRudy();
    normalized_rudy = computeNormalizedRudyGrid(rudy);
  }

  std::vector<ScenarioResult> scenario_results;
  scenario_results.push_back(baseline);

  ScenarioDefinition baseline_def{"baseline", nullptr, nullptr, nullptr, false};
  std::vector<ScenarioDefinition> scenario_defs;

  ScenarioDefinition sporder_shortest_def;
  sporder_shortest_def.name = "sporder-shortest";
  sporder_shortest_def.pre_init = [this, seed = 31]() {
    grouter_->setCapacitiesPerturbationPercentage(0.0f);
    grouter_->setPerturbationAmount(0);
    grouter_->setAllowCongestion(true);
    grouter_->setSeed(seed);
    grouter_->fastroute_->setCriticalNetsPercentage(24.0f);
  };
  sporder_shortest_def.order_nets = [](std::vector<Net*>& scenario_nets) {
    reorderNetsByWirelengthPriority(scenario_nets);
  };
  sporder_shortest_def.post_init
      = [this, &focused_hotspots, min_routing_layer, max_routing_layer]() {
          applyUniformCapacityBoost(
              grouter_, min_routing_layer, max_routing_layer, 1.02f);
          applyHotspotPenalties(grouter_,
                                focused_hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                0,
                                0.99f,
                                0.03f);
        };
  sporder_shortest_def.aggressive = false;
  scenario_defs.push_back(std::move(sporder_shortest_def));

  ScenarioDefinition bsp_scheduler_def;
  bsp_scheduler_def.name = "bsp-scheduler";
  bsp_scheduler_def.pre_init = [this, seed = 53]() {
    grouter_->setCapacitiesPerturbationPercentage(0.0f);
    grouter_->setPerturbationAmount(0);
    grouter_->setAllowCongestion(true);
    grouter_->setSeed(seed);
    grouter_->fastroute_->setCriticalNetsPercentage(20.0f);
  };
  bsp_scheduler_def.order_nets = [](std::vector<Net*>& scenario_nets) {
    reorderNetsByBspScheduler(scenario_nets);
  };
  bsp_scheduler_def.aggressive = false;
  scenario_defs.push_back(std::move(bsp_scheduler_def));

  ScenarioDefinition spatial_round_robin_def;
  spatial_round_robin_def.name = "spatial-roundrobin-turbo";
  spatial_round_robin_def.pre_init = [this, seed = 47]() {
    grouter_->setCapacitiesPerturbationPercentage(0.0f);
    grouter_->setPerturbationAmount(0);
    grouter_->setAllowCongestion(true);
    grouter_->setSeed(seed);
    grouter_->fastroute_->setCriticalNetsPercentage(26.0f);
  };
  spatial_round_robin_def.order_nets = [](std::vector<Net*>& scenario_nets) {
    reorderNetsBySpatialRoundRobin(scenario_nets);
  };
  spatial_round_robin_def.post_init
      = [this, &focused_hotspots, min_routing_layer, max_routing_layer]() {
          applyUniformCapacityBoost(
              grouter_, min_routing_layer, max_routing_layer, 1.06f);
          applyHotspotPenalties(grouter_,
                                focused_hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                0,
                                0.985f,
                                0.06f);
        };
  spatial_round_robin_def.aggressive = true;
  scenario_defs.push_back(std::move(spatial_round_robin_def));

  ScenarioDefinition ultra_compact_bsp_def;
  ultra_compact_bsp_def.name = "ultra-compact-bsp";
  ultra_compact_bsp_def.pre_init = [this, seed = 71]() {
    grouter_->setCapacitiesPerturbationPercentage(0.0f);
    grouter_->setPerturbationAmount(0);
    grouter_->setAllowCongestion(true);
    grouter_->setSeed(seed);
    // Treat most nets as critical so maze routing shrinks detour windows and
    // strongly prefers straight reconnection paths.
    grouter_->fastroute_->setCriticalNetsPercentage(88.0f);
  };
  ultra_compact_bsp_def.order_nets = [](std::vector<Net*>& scenario_nets) {
    reorderNetsByBspThenHpwlBurst(scenario_nets);
  };
  ultra_compact_bsp_def.post_init
      = [this, &focused_hotspots, min_routing_layer, max_routing_layer]() {
          // Keep global resources slightly relaxed to preserve shortest
          // Manhattan paths unless heavy overflow appears.
          applyUniformCapacityBoost(
              grouter_, min_routing_layer, max_routing_layer, 1.12f);
          applyHotspotPenalties(grouter_,
                                focused_hotspots,
                                min_routing_layer,
                                max_routing_layer,
                                0,
                                0.995f,
                                0.02f);
        };
  ultra_compact_bsp_def.aggressive = true;
  scenario_defs.push_back(std::move(ultra_compact_bsp_def));

  if (!normalized_rudy.empty()) {
    ScenarioDefinition wl_direct_focus_def;
    wl_direct_focus_def.name = "wl-direct-focused";
    wl_direct_focus_def.pre_init = [this, seed = 83]() {
      grouter_->setCapacitiesPerturbationPercentage(0.0f);
      grouter_->setPerturbationAmount(0);
      grouter_->setAllowCongestion(true);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(72.0f);
    };
    wl_direct_focus_def.order_nets = [](std::vector<Net*>& scenario_nets) {
      reorderNetsByBspThenHpwlBurst(scenario_nets);
    };
    wl_direct_focus_def.post_init = [this,
                                     &normalized_rudy,
                                     &focused_hotspots,
                                     min_routing_layer,
                                     max_routing_layer]() {
      // FastRoute-style shortest-path bias with SPRoute-style deterministic
      // ordering and a light CUGR-inspired soft-cap remap.
      applyUniformCapacityBoost(
          grouter_, min_routing_layer, max_routing_layer, 1.10f);
      applyHybridCapacityRemap(grouter_,
                               normalized_rudy,
                               focused_hotspots,
                               min_routing_layer,
                               max_routing_layer,
                               0.96f,
                               1.16f,
                               2.2f,
                               0.78f,
                               0.42f,
                               0);
      applyHotspotPenalties(grouter_,
                            focused_hotspots,
                            min_routing_layer,
                            max_routing_layer,
                            0,
                            0.97f,
                            0.18f);
    };
    wl_direct_focus_def.aggressive = true;
    scenario_defs.push_back(std::move(wl_direct_focus_def));

    ScenarioDefinition cugr_softcap_wl_def;
    cugr_softcap_wl_def.name = "cugr-softcap-wirelength";
    cugr_softcap_wl_def.pre_init = [this, seed = 67]() {
      grouter_->setCapacitiesPerturbationPercentage(0.0f);
      grouter_->setPerturbationAmount(0);
      grouter_->setAllowCongestion(true);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(24.0f);
    };
    cugr_softcap_wl_def.order_nets = [](std::vector<Net*>& scenario_nets) {
      reorderNetsByWirelengthPriority(scenario_nets);
    };
    cugr_softcap_wl_def.post_init = [this,
                                     &normalized_rudy,
                                     &focused_hotspots,
                                     min_routing_layer,
                                     max_routing_layer]() {
      applyHybridCapacityRemap(grouter_,
                               normalized_rudy,
                               focused_hotspots,
                               min_routing_layer,
                               max_routing_layer,
                               0.94f,
                               1.08f,
                               3.1f,
                               0.75f,
                               0.28f,
                               0);
      applySoftCapacityScaling(grouter_,
                               normalized_rudy,
                               min_routing_layer,
                               max_routing_layer,
                               0.90f,
                               0.98f,
                               2.2f,
                               0.84f);
      applyUniformCapacityBoost(
          grouter_, min_routing_layer, max_routing_layer, 1.025f);
      applyHotspotPenalties(grouter_,
                            focused_hotspots,
                            min_routing_layer,
                            max_routing_layer,
                            0,
                            0.99f,
                            0.05f);
    };
    cugr_softcap_wl_def.aggressive = true;
    scenario_defs.push_back(std::move(cugr_softcap_wl_def));

    ScenarioDefinition wl_squeeze_hybrid_def;
    wl_squeeze_hybrid_def.name = "wl-squeeze-hybrid";
    wl_squeeze_hybrid_def.pre_init = [this, seed = 79]() {
      grouter_->setCapacitiesPerturbationPercentage(0.0f);
      grouter_->setPerturbationAmount(0);
      grouter_->setAllowCongestion(true);
      grouter_->setSeed(seed);
      grouter_->fastroute_->setCriticalNetsPercentage(32.0f);
    };
    wl_squeeze_hybrid_def.order_nets = [](std::vector<Net*>& scenario_nets) {
      reorderNetsByBspThenHpwlBurst(scenario_nets);
    };
    wl_squeeze_hybrid_def.post_init = [this,
                                       &normalized_rudy,
                                       &focused_hotspots,
                                       min_routing_layer,
                                       max_routing_layer]() {
      applyHybridCapacityRemap(grouter_,
                               normalized_rudy,
                               focused_hotspots,
                               min_routing_layer,
                               max_routing_layer,
                               0.82f,
                               1.18f,
                               2.8f,
                               0.68f,
                               0.55f,
                               1);
      applySoftCapacityScaling(grouter_,
                               normalized_rudy,
                               min_routing_layer,
                               max_routing_layer,
                               0.78f,
                               1.04f,
                               2.3f,
                               0.72f);
      applyUniformCapacityBoost(
          grouter_, min_routing_layer, max_routing_layer, 1.04f);
      applyHotspotPenalties(grouter_,
                            focused_hotspots,
                            min_routing_layer,
                            max_routing_layer,
                            1,
                            0.95f,
                            0.12f);
    };
    wl_squeeze_hybrid_def.aggressive = true;
    scenario_defs.push_back(std::move(wl_squeeze_hybrid_def));
  }

  for (const ScenarioDefinition& def : scenario_defs) {
    ScenarioResult result = run_scenario(def, snapshot);
    scenario_results.push_back(std::move(result));
  }

  auto find_scenario_result = [&](const std::string& name) -> ScenarioResult* {
    for (ScenarioResult& result : scenario_results) {
      if (result.name == name) {
        return &result;
      }
    }
    return nullptr;
  };

  // Cross-router donor collection: pull full-route candidates from CUGR and
  // SPRoute so per-net fusion can use different router strengths.
  if (cugr_ != nullptr) {
    try {
      ScenarioResult cugr_donor;
      cugr_donor.name = "cugr-router-donor";
      cugr_->init(min_routing_layer, max_routing_layer);
      cugr_->route();
      cugr_donor.routes = cugr_->getRoutes();
      cugr_donor.metrics = compute_metrics(cugr_donor.routes);
      logger_->info(GNR,
                    6012,
                    "NEWGR scenario {} [router-donor]: wirelength {:.0f} um, vias {}",
                    cugr_donor.name,
                    cugr_donor.metrics.wirelength_um,
                    cugr_donor.metrics.via_count);
      scenario_results.push_back(std::move(cugr_donor));
    } catch (...) {
      logger_->warn(GNR, 6013, "NEWGR: CUGR donor routing failed; continuing without donor.");
    }
  }

  if (grouter_->sproute_adapter_ != nullptr && grouter_->hasSprouteGridData()
      && grouter_->hasSprouteNetData()) {
    try {
      ScenarioResult sproute_donor;
      sproute_donor.name = "sproute-router-donor";
      grouter_->sproute_adapter_->initialize(grouter_->sproute_grid_data_,
                                             grouter_->sproute_nets_);
      sproute_donor.routes = grouter_->sproute_adapter_->run();
      sproute_donor.metrics = compute_metrics(sproute_donor.routes);
      logger_->info(GNR,
                    6014,
                    "NEWGR scenario {} [router-donor]: wirelength {:.0f} um, vias {}",
                    sproute_donor.name,
                    sproute_donor.metrics.wirelength_um,
                    sproute_donor.metrics.via_count);
      scenario_results.push_back(std::move(sproute_donor));
    } catch (...) {
      logger_->warn(
          GNR, 6015, "NEWGR: SPRoute donor routing failed; continuing without donor.");
    }
  }

  if (ScenarioResult* sporder = find_scenario_result("sporder-shortest")) {
    if (ScenarioResult* spatial
        = find_scenario_result("spatial-roundrobin-turbo")) {
      ScenarioResult fusion;
      fusion.name = "spatial-wirelength-grafting";
      fusion.routes = sporder->routes;

      const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
      const int x_min = grouter_->grid_->getXMin();
      const int y_min = grouter_->grid_->getYMin();
      const int x_grids = grouter_->grid_->getXGrids();
      const int y_grids = grouter_->grid_->getYGrids();

      const int spatial_swaps = applySelectiveNetRouteGrafting(fusion.routes,
                                                               spatial->routes,
                                                               tile_size,
                                                               x_min,
                                                               y_min,
                                                               x_grids,
                                                               y_grids,
                                                               hotspot_map);
      if (ScenarioResult* direct_focus = find_scenario_result("wl-direct-focused")) {
        applySelectiveNetRouteGrafting(fusion.routes,
                                       direct_focus->routes,
                                       tile_size,
                                       x_min,
                                       y_min,
                                       x_grids,
                                       y_grids,
                                       hotspot_map);
      }

      fusion.metrics = compute_metrics(fusion.routes);
      logger_->info(GNR,
                    6008,
                    "NEWGR scenario {} [hybrid]: wirelength {:.0f} um, vias {}, grafted nets {}",
                    fusion.name,
                    fusion.metrics.wirelength_um,
                    fusion.metrics.via_count,
                    spatial_swaps);
      scenario_results.push_back(std::move(fusion));
    }
  }

  if (ScenarioResult* sporder = find_scenario_result("sporder-shortest")) {
    ScenarioResult fusion;
    fusion.name = "multi-router-wirelength-fusion";
    fusion.routes = sporder->routes;

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> donors;
    donors.reserve(10);
    const std::vector<std::string> donor_names{
        "spatial-wirelength-grafting",
        "wl-direct-focused",
        "cugr-softcap-wirelength",
        "wl-squeeze-hybrid",
        "ultra-compact-bsp",
        "spatial-roundrobin-turbo",
        "bsp-scheduler",
        "cugr-router-donor",
        "sproute-router-donor",
        "baseline"};
    for (const std::string& donor_name : donor_names) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        donors.push_back(&donor->routes);
      }
    }

    const int fused_swaps = applyMultiScenarioWirelengthFusion(
        fusion.routes,
        donors,
        tile_size,
        x_min,
        y_min,
        x_grids,
        y_grids,
        hotspot_map,
        std::max(tile_size / 2, 1),
        2,
        0.45);
    fusion.metrics = compute_metrics(fusion.routes);
    logger_->info(GNR,
                  6010,
                  "NEWGR scenario {} [hybrid]: wirelength {:.0f} um, vias {}, fused nets {}",
                  fusion.name,
                  fusion.metrics.wirelength_um,
                  fusion.metrics.via_count,
                  fused_swaps);
    scenario_results.push_back(std::move(fusion));
  }

  auto best_wirelength_iter
      = std::min_element(scenario_results.begin(),
                         scenario_results.end(),
                         [](const ScenarioResult& lhs, const ScenarioResult& rhs) {
                           if (lhs.metrics.wirelength_dbu != rhs.metrics.wirelength_dbu) {
                             return lhs.metrics.wirelength_dbu < rhs.metrics.wirelength_dbu;
                           }
                           return lhs.metrics.via_count < rhs.metrics.via_count;
                         });

  const bool has_best_wirelength = best_wirelength_iter != scenario_results.end();
  std::string best_wirelength_name;
  NetRouteMap best_wirelength_routes;
  if (has_best_wirelength) {
    best_wirelength_name = best_wirelength_iter->name;
    best_wirelength_routes = best_wirelength_iter->routes;
  }

  if (has_best_wirelength) {
    ScenarioResult fusion;
    fusion.name = "cross-router-wirelength-fusion";
    fusion.routes = best_wirelength_routes;

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> donors;
    donors.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      if (result.name == fusion.name || result.name == best_wirelength_name) {
        continue;
      }
      donors.push_back(&result.routes);
    }

    int fused_swaps = applyMultiScenarioWirelengthFusion(fusion.routes,
                                                         donors,
                                                         tile_size,
                                                         x_min,
                                                         y_min,
                                                         x_grids,
                                                         y_grids,
                                                         hotspot_map,
                                                         std::max(tile_size / 4, 1),
                                                         3,
                                                         0.70);

    std::vector<const NetRouteMap*> radical_donors;
    radical_donors.reserve(4);
    for (const char* donor_name : {"cugr-router-donor",
                                   "sproute-router-donor",
                                   "multi-router-wirelength-fusion",
                                   "spatial-wirelength-grafting"}) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        radical_donors.push_back(&donor->routes);
      }
    }
    fused_swaps += applyMultiScenarioWirelengthFusion(fusion.routes,
                                                      radical_donors,
                                                      tile_size,
                                                      x_min,
                                                      y_min,
                                                      x_grids,
                                                      y_grids,
                                                      hotspot_map,
                                                      1,
                                                      4,
                                                      0.90);

    fusion.metrics = compute_metrics(fusion.routes);
    logger_->info(
        GNR,
        6016,
        "NEWGR scenario {} [cross-router]: wirelength {:.0f} um, vias {}, fused nets {}",
        fusion.name,
        fusion.metrics.wirelength_um,
        fusion.metrics.via_count,
        fused_swaps);
    scenario_results.push_back(std::move(fusion));
  }

  if (has_best_wirelength) {
    ScenarioResult radical;
    radical.name = "radical-shortpath-fusion";
    radical.routes = best_wirelength_routes;

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> donors;
    donors.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      if (result.name == radical.name || result.name == best_wirelength_name) {
        continue;
      }
      donors.push_back(&result.routes);
    }

    int fused_swaps = applyAdaptiveWirelengthFusion(radical.routes,
                                                    donors,
                                                    tile_size,
                                                    x_min,
                                                    y_min,
                                                    x_grids,
                                                    y_grids,
                                                    hotspot_map,
                                                    std::max(tile_size / 6, 1),
                                                    10,
                                                    1.10);

    std::vector<const NetRouteMap*> radical_donors;
    radical_donors.reserve(5);
    for (const char* donor_name : {"cugr-router-donor",
                                   "sproute-router-donor",
                                   "cross-router-wirelength-fusion",
                                   "multi-router-wirelength-fusion",
                                   "spatial-wirelength-grafting"}) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        radical_donors.push_back(&donor->routes);
      }
    }

    fused_swaps += applyAdaptiveWirelengthFusion(radical.routes,
                                                 radical_donors,
                                                 tile_size,
                                                 x_min,
                                                 y_min,
                                                 x_grids,
                                                 y_grids,
                                                 hotspot_map,
                                                 1,
                                                 16,
                                                 2.20);

    radical.metrics = compute_metrics(radical.routes);
    logger_->info(
        GNR,
        6017,
        "NEWGR scenario {} [radical]: wirelength {:.0f} um, vias {}, fused nets {}",
        radical.name,
        radical.metrics.wirelength_um,
        radical.metrics.via_count,
        fused_swaps);
    scenario_results.push_back(std::move(radical));
  }

  if (has_best_wirelength) {
    ScenarioResult extreme;
    extreme.name = "extreme-wirelength-stitch";
    if (ScenarioResult* radical = find_scenario_result("radical-shortpath-fusion")) {
      extreme.routes = radical->routes;
    } else if (ScenarioResult* cross = find_scenario_result("cross-router-wirelength-fusion")) {
      extreme.routes = cross->routes;
    } else {
      extreme.routes = best_wirelength_routes;
    }

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> donors;
    donors.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      if (result.name == extreme.name) {
        continue;
      }
      donors.push_back(&result.routes);
    }

    int fused_swaps = applyExtremeWirelengthFusion(extreme.routes,
                                                   donors,
                                                   tile_size,
                                                   x_min,
                                                   y_min,
                                                   x_grids,
                                                   y_grids,
                                                   hotspot_map,
                                                   std::max(tile_size / 2, 1),
                                                   0.0150,
                                                   3,
                                                   3.40,
                                                   0.45);

    std::vector<const NetRouteMap*> specialist_donors;
    specialist_donors.reserve(6);
    for (const char* donor_name : {"radical-shortpath-fusion",
                                   "cross-router-wirelength-fusion",
                                   "multi-router-wirelength-fusion",
                                   "spatial-wirelength-grafting"}) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        specialist_donors.push_back(&donor->routes);
      }
    }

    fused_swaps += applyExtremeWirelengthFusion(extreme.routes,
                                                specialist_donors,
                                                tile_size,
                                                x_min,
                                                y_min,
                                                x_grids,
                                                y_grids,
                                                hotspot_map,
                                                std::max(tile_size / 3, 1),
                                                0.0100,
                                                4,
                                                2.80,
                                                0.80);

    extreme.metrics = compute_metrics(extreme.routes);
    logger_->info(GNR,
                  6018,
                  "NEWGR scenario {} [radical]: wirelength {:.0f} um, vias {}, fused nets {}",
                  extreme.name,
                  extreme.metrics.wirelength_um,
                  extreme.metrics.via_count,
                  fused_swaps);
    scenario_results.push_back(std::move(extreme));
  }

  if (has_best_wirelength) {
    ScenarioResult stabilized;
    stabilized.name = "stabilized-dr-fusion";
    if (ScenarioResult* extreme = find_scenario_result("extreme-wirelength-stitch")) {
      stabilized.routes = extreme->routes;
    } else if (ScenarioResult* radical = find_scenario_result("radical-shortpath-fusion")) {
      stabilized.routes = radical->routes;
    } else {
      stabilized.routes = best_wirelength_routes;
    }

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> stabilization_donors;
    stabilization_donors.reserve(8);
    for (const char* donor_name : {"radical-shortpath-fusion",
                                   "cross-router-wirelength-fusion",
                                   "multi-router-wirelength-fusion",
                                   "spatial-wirelength-grafting",
                                   "ultra-compact-bsp",
                                   "sporder-shortest",
                                   "bsp-scheduler",
                                   "baseline"}) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        stabilization_donors.push_back(&donor->routes);
      }
    }

    const int stabilized_swaps = applyViaAwareStabilizationFusion(
        stabilized.routes,
        stabilization_donors,
        tile_size,
        x_min,
        y_min,
        x_grids,
        y_grids,
        hotspot_map);

    stabilized.metrics = compute_metrics(stabilized.routes);
    logger_->info(
        GNR,
        6020,
        "NEWGR scenario {} [stabilized]: wirelength {:.0f} um, vias {}, stabilized nets {}",
        stabilized.name,
        stabilized.metrics.wirelength_um,
        stabilized.metrics.via_count,
        stabilized_swaps);
    scenario_results.push_back(std::move(stabilized));
  }

  if (has_best_wirelength) {
    ScenarioResult longnet_fusion;
    longnet_fusion.name = "longnet-priority-fusion";
    if (ScenarioResult* extreme = find_scenario_result("extreme-wirelength-stitch")) {
      longnet_fusion.routes = extreme->routes;
    } else if (ScenarioResult* radical = find_scenario_result("radical-shortpath-fusion")) {
      longnet_fusion.routes = radical->routes;
    } else {
      longnet_fusion.routes = best_wirelength_routes;
    }

    const int tile_size = std::max(grouter_->grid_->getTileSize(), 1);
    const int x_min = grouter_->grid_->getXMin();
    const int y_min = grouter_->grid_->getYMin();
    const int x_grids = grouter_->grid_->getXGrids();
    const int y_grids = grouter_->grid_->getYGrids();

    std::vector<const NetRouteMap*> donors;
    donors.reserve(scenario_results.size());
    for (const ScenarioResult& result : scenario_results) {
      if (result.name == longnet_fusion.name) {
        continue;
      }
      donors.push_back(&result.routes);
    }

    int fused_swaps = applyLongNetPriorityFusion(
        longnet_fusion.routes,
        donors,
        tile_size,
        x_min,
        y_min,
        x_grids,
        y_grids,
        hotspot_map,
        std::max(tile_size / 3, 1),
        6,
        1.20);

    std::vector<const NetRouteMap*> specialist_donors;
    specialist_donors.reserve(7);
    for (const char* donor_name : {"radical-shortpath-fusion",
                                   "cross-router-wirelength-fusion",
                                   "multi-router-wirelength-fusion",
                                   "spatial-wirelength-grafting",
                                   "wl-direct-focused",
                                   "sporder-shortest"}) {
      if (ScenarioResult* donor = find_scenario_result(donor_name)) {
        specialist_donors.push_back(&donor->routes);
      }
    }

    fused_swaps += applyLongNetPriorityFusion(longnet_fusion.routes,
                                              specialist_donors,
                                              tile_size,
                                              x_min,
                                              y_min,
                                              x_grids,
                                              y_grids,
                                              hotspot_map,
                                              1,
                                              12,
                                              2.00);

    longnet_fusion.metrics = compute_metrics(longnet_fusion.routes);
    logger_->info(
        GNR,
        6019,
        "NEWGR scenario {} [radical]: wirelength {:.0f} um, vias {}, fused nets {}",
        longnet_fusion.name,
        longnet_fusion.metrics.wirelength_um,
        longnet_fusion.metrics.via_count,
        fused_swaps);
    scenario_results.push_back(std::move(longnet_fusion));
  }

  const long baseline_vias = baseline.metrics.via_count;
  auto better_result = [baseline_vias](const ScenarioResult& lhs,
                                       const ScenarioResult& rhs) {
    const long lhs_wl = lhs.metrics.wirelength_dbu;
    const long rhs_wl = rhs.metrics.wirelength_dbu;
    const long wl_tie_window
        = std::max<long>(4500, std::max(lhs_wl, rhs_wl) / 70000);
    const long wl_delta = lhs_wl > rhs_wl ? lhs_wl - rhs_wl : rhs_wl - lhs_wl;
    if (wl_delta > wl_tie_window) {
      return lhs_wl < rhs_wl;
    }

    if (baseline_vias > 0) {
      const long max_reasonable_via
          = static_cast<long>(std::ceil(1.32 * baseline_vias));
      const bool lhs_via_high = lhs.metrics.via_count > max_reasonable_via;
      const bool rhs_via_high = rhs.metrics.via_count > max_reasonable_via;
      if (lhs_via_high != rhs_via_high) {
        if (lhs_via_high) {
          return lhs_wl + wl_tie_window < rhs_wl;
        }
        return !(rhs_wl + wl_tie_window < lhs_wl);
      }
    }

    if (lhs.metrics.via_count != rhs.metrics.via_count) {
      return lhs.metrics.via_count < rhs.metrics.via_count;
    }
    if (lhs_wl != rhs_wl) {
      return lhs_wl < rhs_wl;
    }
    return lhs.metrics.score < rhs.metrics.score;
  };

  auto best_iter
      = std::min_element(scenario_results.begin(), scenario_results.end(), better_result);

  if (best_iter == scenario_results.end()) {
    return {};
  }

  ScenarioResult final_result = *best_iter;

  const ScenarioDefinition* replay_def = nullptr;
  if (best_iter->name != scenario_results.back().name) {
    if (best_iter->name == "baseline") {
      replay_def = &baseline_def;
    } else {
      for (const ScenarioDefinition& def : scenario_defs) {
        if (def.name == best_iter->name) {
          replay_def = &def;
          break;
        }
      }
    }
    if (replay_def != nullptr) {
      final_result = run_scenario(*replay_def, snapshot);
    }
  }

  logger_->info(GNR,
                6007,
                "NEWGR best scenario '{}': wirelength {:.0f} um, vias {}",
                final_result.name,
                final_result.metrics.wirelength_um,
                final_result.metrics.via_count);

  const bool apply_patching
      = final_result.name == "cugr-softcap-wirelength"
        || final_result.name.find("fusion") != std::string::npos
        || final_result.name.find("stitch") != std::string::npos;
  if (apply_patching) {
    applyCugrStyleGuidePatching(grouter_,
                                final_result.routes,
                                min_routing_layer,
                                max_routing_layer,
                                logger_);
  }

  return std::move(final_result.routes);
}

}  // namespace grt
