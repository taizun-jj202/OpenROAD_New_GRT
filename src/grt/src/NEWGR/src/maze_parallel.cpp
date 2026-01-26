// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "FastRoute.h"

namespace grt::newgr {

namespace {

static int parent_index(int i)
{
  return (i - 1) / 2;
}

static int left_index(int i)
{
  return 2 * i + 1;
}

static int right_index(int i)
{
  return 2 * i + 2;
}

static inline int heap_elem_index(const double* base, const double* elem)
{
  return static_cast<int>(elem - base);
}

static void heapify(std::vector<double*>& array,
                    std::vector<int>& heap_pos,
                    double* base)
{
  bool stop = false;
  const int heapSize = static_cast<int>(array.size());
  int i = 0;

  double* tmp = array[i];
  do {
    const int l = left_index(i);
    const int r = right_index(i);

    int smallest;
    if (l < heapSize && *(array[l]) < *tmp) {
      smallest = l;
      if (r < heapSize && *(array[r]) < *(array[l])) {
        smallest = r;
      }
    } else {
      smallest = i;
      if (r < heapSize && *(array[r]) < *tmp) {
        smallest = r;
      }
    }
    if (smallest != i) {
      array[i] = array[smallest];
      heap_pos[heap_elem_index(base, array[i])] = i;
      i = smallest;
    } else {
      array[i] = tmp;
      heap_pos[heap_elem_index(base, tmp)] = i;
      stop = true;
    }
  } while (!stop);
}

static void updateHeap(std::vector<double*>& array,
                       std::vector<int>& heap_pos,
                       double* base,
                       int i)
{
  double* tmpi = array[i];
  const int tmp_index = heap_elem_index(base, tmpi);
  while (i > 0 && *(array[parent_index(i)]) > *tmpi) {
    const int parent = parent_index(i);
    array[i] = array[parent];
    heap_pos[heap_elem_index(base, array[i])] = i;
    i = parent;
  }
  array[i] = tmpi;
  heap_pos[tmp_index] = i;
}

static void removeMin(std::vector<double*>& array,
                      std::vector<int>& heap_pos,
                      double* base)
{
  if (array.empty()) {
    return;
  }

  heap_pos[heap_elem_index(base, array[0])] = -1;

  if (array.size() == 1) {
    array.pop_back();
    return;
  }

  array[0] = array.back();
  array.pop_back();
  heap_pos[heap_elem_index(base, array[0])] = 0;
  heapify(array, heap_pos, base);
}

struct UsageUpdate
{
  bool vertical = false;
  int16_t x = 0;
  int16_t y = 0;
  int delta = 0;
  FrNet* net = nullptr;
};

static void recordUsageH(void* ctx, int x, int y, FrNet* net, int delta)
{
  auto* updates = static_cast<std::vector<UsageUpdate>*>(ctx);
  updates->push_back(
      UsageUpdate{false, static_cast<int16_t>(x), static_cast<int16_t>(y), delta, net});
}

static void recordUsageV(void* ctx, int x, int y, FrNet* net, int delta)
{
  auto* updates = static_cast<std::vector<UsageUpdate>*>(ctx);
  updates->push_back(
      UsageUpdate{true, static_cast<int16_t>(x), static_cast<int16_t>(y), delta, net});
}

struct Maze2DScratch
{
  int x_range = 0;
  int y_range = 0;
  int x_grid = 0;
  int y_grid = 0;

  multi_array<short, 2> parent_x1;
  multi_array<short, 2> parent_y1;
  multi_array<short, 2> parent_x3;
  multi_array<short, 2> parent_y3;
  multi_array<bool, 2> hv;
  multi_array<bool, 2> hyper_v;
  multi_array<bool, 2> hyper_h;
  multi_array<bool, 2> in_region;
  multi_array<int, 2> corr_edge;
  multi_array<double, 2> d1;
  multi_array<double, 2> d2;

  std::vector<bool> pop_heap2;
  std::vector<int> src_heap_pos;
  std::vector<int> src_heap_touched;
  std::vector<double*> src_heap;
  std::vector<double*> dest_heap;

  void ensure(int x_range_in, int y_range_in, int x_grid_in, int y_grid_in)
  {
    if (x_range == x_range_in && y_range == y_range_in && x_grid == x_grid_in
        && y_grid == y_grid_in) {
      return;
    }
    x_range = x_range_in;
    y_range = y_range_in;
    x_grid = x_grid_in;
    y_grid = y_grid_in;

    parent_x1.resize(boost::extents[y_range][x_range]);
    parent_y1.resize(boost::extents[y_range][x_range]);
    parent_x3.resize(boost::extents[y_range][x_range]);
    parent_y3.resize(boost::extents[y_range][x_range]);
    hv.resize(boost::extents[y_range][x_range]);
    hyper_v.resize(boost::extents[y_range][x_range]);
    hyper_h.resize(boost::extents[y_range][x_range]);
    in_region.resize(boost::extents[y_range][x_range]);
    corr_edge.resize(boost::extents[y_range][x_range]);
    d1.resize(boost::extents[y_range][x_range]);
    d2.resize(boost::extents[y_range][x_range]);

    const int heap_map_size = y_grid * x_range;
    pop_heap2.assign(heap_map_size, false);
    src_heap_pos.assign(heap_map_size, -1);
    src_heap_touched.clear();
    src_heap_touched.reserve(1024);
    src_heap.clear();
    dest_heap.clear();
    src_heap.reserve(static_cast<size_t>(y_grid) * static_cast<size_t>(x_range));
    dest_heap.reserve(static_cast<size_t>(y_grid) * static_cast<size_t>(x_range));
  }
};

}  // namespace

void FastRouteCore::mazeRouteMSMDParallel(const int iter,
                                         const int expand,
                                         const int ripup_threshold,
                                         const int maze_edge_threshold,
                                         const bool ordering,
                                         const int via,
                                         const int L,
                                         const CostParams& cost_params,
                                         float& slack_th)
{
#ifndef _OPENMP
  mazeRouteMSMD(iter,
                expand,
                ripup_threshold,
                maze_edge_threshold,
                ordering,
                via,
                L,
                cost_params,
                slack_th);
  return;
#else
  const int max_threads = omp_get_max_threads();
  if (max_threads <= 1 || static_cast<int>(net_ids_.size()) < 512) {
    mazeRouteMSMD(iter,
                  expand,
                  ripup_threshold,
                  maze_edge_threshold,
                  ordering,
                  via,
                  L,
                  cost_params,
                  slack_th);
    return;
  }

  int threads = std::min(max_threads, 16);
  if (threads < 2) {
    threads = 2;
  }

  int batch_size
      = static_cast<int>(net_ids_.size()) / std::max(1, threads * 2);
  batch_size = std::clamp(batch_size, 256, 2048);
  if (iter > 10) {
    batch_size = std::max(128, batch_size / 2);
  }
  if (iter > 20) {
    batch_size = std::max(96, batch_size / 2);
  }

  auto route_one = [&](int nidRPC, const UsageUpdateCallbacks* updates_cb) {
    static thread_local Maze2DScratch scratch;
    scratch.ensure(x_range_, y_range_, x_grid_, y_grid_);

    auto& src_heap = scratch.src_heap;
    auto& dest_heap = scratch.dest_heap;
    auto& d1 = scratch.d1;
    auto& d2 = scratch.d2;
    auto& pop_heap2 = scratch.pop_heap2;
    auto& src_heap_pos = scratch.src_heap_pos;
    auto& src_heap_touched = scratch.src_heap_touched;
    auto& parent_x1_ = scratch.parent_x1;
    auto& parent_y1_ = scratch.parent_y1;
    auto& parent_x3_ = scratch.parent_x3;
    auto& parent_y3_ = scratch.parent_y3;
    auto& hv_ = scratch.hv;
    auto& hyper_h_ = scratch.hyper_h;
    auto& hyper_v_ = scratch.hyper_v;
    auto& in_region_ = scratch.in_region;
    auto& corr_edge_ = scratch.corr_edge;

    auto update_usage_h = [&](int x, int y, FrNet* net, int delta) {
      if (updates_cb != nullptr && updates_cb->enabled()) {
        updates_cb->updateH(updates_cb->ctx, x, y, net, delta);
      } else {
        graph2d_.updateUsageH(x, y, net, delta);
      }
    };
    auto update_usage_v = [&](int x, int y, FrNet* net, int delta) {
      if (updates_cb != nullptr && updates_cb->enabled()) {
        updates_cb->updateV(updates_cb->ctx, x, y, net, delta);
      } else {
        graph2d_.updateUsageV(x, y, net, delta);
      }
    };

    const int netID
        = ordering ? tree_order_cong_[nidRPC].treeIndex : net_ids_[nidRPC];
    const int num_terminals = sttrees_[netID].num_terminals;
    const int origENG = expand;

    std::vector<OrderNetEdge> net_eo;
    net_eo.reserve(2ul * max_degree_);
    netedgeOrderDec(netID, net_eo);

    auto& treeedges = sttrees_[netID].edges;
    auto& treenodes = sttrees_[netID].nodes;
    const int num_edges = sttrees_[netID].num_edges();

    for (int edgeREC = 0; edgeREC < num_edges; edgeREC++) {
      const int edgeID = net_eo[edgeREC].edgeID;
      TreeEdge* treeedge = &(treeedges[edgeID]);

      int n1 = treeedge->n1;
      int n2 = treeedge->n2;
      const int n1x = treenodes[n1].x;
      const int n1y = treenodes[n1].y;
      const int n2x = treenodes[n2].x;
      const int n2y = treenodes[n2].y;
      treeedge->len = abs(n2x - n1x) + abs(n2y - n1y);

      if (treeedge->len <= maze_edge_threshold) {
        continue;
      }

      const bool enter = newRipupCheck(treeedge,
                                       n1x,
                                       n1y,
                                       n2x,
                                       n2y,
                                       ripup_threshold,
                                       slack_th,
                                       netID,
                                       edgeID,
                                       updates_cb);

      if (!enter) {
        continue;
      }

      const auto [ymin, ymax] = std::minmax(n1y, n2y);
      const auto [xmin, xmax] = std::minmax(n1x, n2x);

      const int manhattan_len = treeedge->len;
      const int min_local_expand = 3;
      const double expand_ratio = 0.35;
      const int dynamic_cap
          = min_local_expand
            + static_cast<int>(std::round(manhattan_len * expand_ratio));
      const int local_enlarge
          = std::min(origENG, (iter / 6 + 3) * treeedge->route.routelen);
      const int effective_enlarge
          = std::max(min_local_expand, std::min(local_enlarge, dynamic_cap));

      int decrease = 0;
      if (nets_[netID]->isCritical()) {
        decrease = std::min((iter / 7) * 5, effective_enlarge / 2);
      }
      const int regionX1 = std::max(xmin - effective_enlarge + decrease, 0);
      const int regionX2
          = std::min(xmax + effective_enlarge - decrease, x_grid_ - 1);
      const int regionY1 = std::max(ymin - effective_enlarge + decrease, 0);
      const int regionY2
          = std::min(ymax + effective_enlarge - decrease, y_grid_ - 1);

      const int region_area
          = (regionY2 - regionY1 + 1) * (regionX2 - regionX1 + 1);
      for (int i = regionY1; i <= regionY2; i++) {
        for (int j = regionX1; j <= regionX2; j++) {
          d1[i][j] = BIG_INT;
          d2[i][j] = BIG_INT;
          hyper_h_[i][j] = false;
          hyper_v_[i][j] = false;
        }
      }

      for (const int idx : src_heap_touched) {
        src_heap_pos[idx] = -1;
      }
      src_heap_touched.clear();

      setupHeap(netID,
                edgeID,
                src_heap,
                dest_heap,
                src_heap_pos,
                src_heap_touched,
                d1,
                d2,
                in_region_,
                corr_edge_,
                regionX1,
                regionX2,
                regionY1,
                regionY2);

      double* d1_base = d1.data();

      auto updateAdjacent = [&](const int cur_x,
                                const int cur_y,
                                const int adj_x,
                                const int adj_y,
                                double cost) {
        double adj_cost = d1[adj_y][adj_x];
        if (adj_cost <= cost) {
          return;
        }

        d1[adj_y][adj_x] = cost;

        if (cur_x != adj_x) {
          parent_x3_[adj_y][adj_x] = cur_x;
          parent_y3_[adj_y][adj_x] = cur_y;
          hv_[adj_y][adj_x] = false;
        } else {
          parent_x1_[adj_y][adj_x] = cur_x;
          parent_y1_[adj_y][adj_x] = cur_y;
          hv_[adj_y][adj_x] = true;
        }

        if (adj_cost >= BIG_INT) {
          src_heap.push_back(&d1[adj_y][adj_x]);
          const int idx = adj_y * x_range_ + adj_x;
          if (src_heap_pos[idx] == -1) {
            src_heap_touched.push_back(idx);
          }
          src_heap_pos[idx] = static_cast<int>(src_heap.size()) - 1;
          updateHeap(src_heap, src_heap_pos, d1_base, src_heap.size() - 1);
        } else if (adj_cost > cost) {
          const int idx = adj_y * x_range_ + adj_x;
          const int pos = src_heap_pos[idx];
          if (pos >= 0) {
            updateHeap(src_heap, src_heap_pos, d1_base, pos);
          }
        }
      };

      auto relaxAdjacent = [&](const int cur_x,
                               const int cur_y,
                               const int d_x,
                               const int d_y,
                               const bool add_via,
                               const bool maybe_hyper) {
        const bool is_horizontal = d_x != 0;
        auto& hyper = is_horizontal ? hyper_h_ : hyper_v_;

        const int p1_x = cur_x - (d_x == -1);
        const int p1_y = cur_y - (d_y == -1);
        const int p2_x = cur_x - (d_x == 1);
        const int p2_y = cur_y - (d_y == 1);

        auto usage_red
            = is_horizontal ? &Graph2D::getUsageRedH : &Graph2D::getUsageRedV;
        auto last_usage
            = is_horizontal ? &Graph2D::getLastUsageH : &Graph2D::getLastUsageV;
        const int pos1 = (graph2d_.*usage_red)(p1_x, p1_y)
                         + L * (graph2d_.*last_usage)(p1_x, p1_y);

        const double cost1 = getCost(pos1, is_horizontal, cost_params);
        double tmp = d1[cur_y][cur_x] + cost1;

        if (add_via && d1[cur_y][cur_x] != 0) {
          tmp += via;

          if (maybe_hyper) {
            const int pos2 = (graph2d_.*usage_red)(p2_x, p2_y)
                             + L * (graph2d_.*last_usage)(p2_x, p2_y);
            const double cost2 = getCost(pos2, is_horizontal, cost_params);
            const int tmp_cost = d1[cur_y - d_y][cur_x - d_x] + cost2;
            if (tmp_cost < d1[cur_y][cur_x] + via) {
              hyper[cur_y][cur_x] = true;
            }
          }
        }

        updateAdjacent(cur_x, cur_y, cur_x + d_x, cur_y + d_y, tmp);
      };

      int ind1 = static_cast<int>(src_heap[0] - d1_base);
      for (int i = 0; i < static_cast<int>(dest_heap.size()); i++) {
        pop_heap2[static_cast<size_t>(dest_heap[i] - &d2[0][0])] = true;
      }

      while (pop_heap2[ind1] == false) {
        const int curX = ind1 % x_range_;
        const int curY = ind1 / x_range_;

        int preX = curX;
        int preY = curY;
        if (d1[curY][curX] != 0) {
          preX = hv_[curY][curX] ? parent_x1_[curY][curX]
                                 : parent_x3_[curY][curX];
          preY = hv_[curY][curX] ? parent_y1_[curY][curX]
                                 : parent_y3_[curY][curX];
        }

        removeMin(src_heap, src_heap_pos, d1_base);

        if (curX > regionX1) {
          relaxAdjacent(curX, curY, -1, 0, preY != curY, curX < regionX2 - 1);
        }
        if (curX < regionX2) {
          relaxAdjacent(curX, curY, 1, 0, preY != curY, curX > regionX1 + 1);
        }
        if (curY > regionY1) {
          relaxAdjacent(curX, curY, 0, -1, preX != curX, curY < regionY2 - 1);
        }
        if (curY < regionY2) {
          relaxAdjacent(curX, curY, 0, 1, preX != curX, curY > regionY1 + 1);
        }

        ind1 = static_cast<int>(src_heap[0] - d1_base);
      }

      for (int i = 0; i < static_cast<int>(dest_heap.size()); i++) {
        pop_heap2[static_cast<size_t>(dest_heap[i] - &d2[0][0])] = false;
      }

      const int16_t crossX = ind1 % x_range_;
      const int16_t crossY = ind1 / x_range_;

      int tmpX = 0;
      int tmpY = 0;
      int cnt = 0;
      int16_t curX = crossX;
      int16_t curY = crossY;
      std::vector<GPoint3D> tmp_grids;
      while (d1[curY][curX] != 0) {
        bool hypered = false;
        if (cnt != 0) {
          if (curX != tmpX && hyper_h_[curY][curX]) {
            curX = 2 * curX - tmpX;
            hypered = true;
          }

          if (curY != tmpY && hyper_v_[curY][curX]) {
            curY = 2 * curY - tmpY;
            hypered = true;
          }
        }
        tmpX = curX;
        tmpY = curY;
        if (!hypered) {
          if (hv_[tmpY][tmpX]) {
            curY = parent_y1_[tmpY][tmpX];
          } else {
            curX = parent_x3_[tmpY][tmpX];
          }
        }
        tmp_grids.push_back({curX, curY, -1});
        cnt++;
      }

      std::vector<GPoint3D> grids(tmp_grids.rbegin(), tmp_grids.rend());
      grids.push_back({crossX, crossY, -1});
      cnt++;

      const int cnt_n1n2 = cnt;
      const int E1x = grids[0].x;
      const int E1y = grids[0].y;
      const int E2x = grids.back().x;
      const int E2y = grids.back().y;
      const int edge_n1n2 = edgeID;

      if (n1 < num_terminals && (E1x != n1x || E1y != n1y)) {
        n1 = splitEdge(treeedges, treenodes, n2, n1, edgeID);
      }

      if (n1 >= num_terminals && (E1x != n1x || E1y != n1y)) {
        const int endpt1 = treeedges[corr_edge_[E1y][E1x]].n1;
        const int endpt2 = treeedges[corr_edge_[E1y][E1x]].n2;

        int A1, A2;
        int edge_n1A1, edge_n1A2;
        if (treenodes[n1].nbr[0] == n2) {
          A1 = treenodes[n1].nbr[1];
          A2 = treenodes[n1].nbr[2];
          edge_n1A1 = treenodes[n1].edge[1];
          edge_n1A2 = treenodes[n1].edge[2];
        } else if (treenodes[n1].nbr[1] == n2) {
          A1 = treenodes[n1].nbr[0];
          A2 = treenodes[n1].nbr[2];
          edge_n1A1 = treenodes[n1].edge[0];
          edge_n1A2 = treenodes[n1].edge[2];
        } else {
          A1 = treenodes[n1].nbr[0];
          A2 = treenodes[n1].nbr[1];
          edge_n1A1 = treenodes[n1].edge[0];
          edge_n1A2 = treenodes[n1].edge[1];
        }

        if (endpt1 == n1 || endpt2 == n1) {
          if (endpt1 == A2 || endpt2 == A2) {
            std::swap(A1, A2);
            std::swap(edge_n1A1, edge_n1A2);
          }

          bool route_ok = updateRouteType1(netID,
                                           treenodes,
                                           n1,
                                           A1,
                                           A2,
                                           E1x,
                                           E1y,
                                           treeedges,
                                           edge_n1A1,
                                           edge_n1A2);
          if (!route_ok) {
            return;
          }

          treenodes[n1].x = E1x;
          treenodes[n1].y = E1y;
        } else {
          const int C1 = endpt1;
          const int C2 = endpt2;
          const int edge_C1C2 = corr_edge_[E1y][E1x];

          bool route_ok = updateRouteType2(netID,
                                           treenodes,
                                           n1,
                                           A1,
                                           A2,
                                           C1,
                                           C2,
                                           E1x,
                                           E1y,
                                           treeedges,
                                           edge_n1A1,
                                           edge_n1A2,
                                           edge_C1C2);
          if (!route_ok) {
            return;
          }

          treenodes[n1].x = E1x;
          treenodes[n1].y = E1y;

          const int edge_n1C1 = edge_n1A1;
          treeedges[edge_n1C1].n1 = C1;
          treeedges[edge_n1C1].n2 = n1;
          const int edge_n1C2 = edge_n1A2;
          treeedges[edge_n1C2].n1 = n1;
          treeedges[edge_n1C2].n2 = C2;
          const int edge_A1A2 = edge_C1C2;
          treeedges[edge_A1A2].n1 = A1;
          treeedges[edge_A1A2].n2 = A2;

          treenodes[n1].nbr[0] = n2;
          treenodes[n1].edge[0] = edge_n1n2;
          treenodes[n1].nbr[1] = C1;
          treenodes[n1].edge[1] = edge_n1C1;
          treenodes[n1].nbr[2] = C2;
          treenodes[n1].edge[2] = edge_n1C2;
        }
      }

      if (n2 < num_terminals && (E2x != n2x || E2y != n2y)) {
        n2 = splitEdge(treeedges, treenodes, n1, n2, edgeID);
      }

      if (n2 >= num_terminals && (E2x != n2x || E2y != n2y)) {
        const int endpt1 = treeedges[corr_edge_[E2y][E2x]].n1;
        const int endpt2 = treeedges[corr_edge_[E2y][E2x]].n2;

        int B1, B2;
        int edge_n2B1, edge_n2B2;
        if (treenodes[n2].nbr[0] == n1) {
          B1 = treenodes[n2].nbr[1];
          B2 = treenodes[n2].nbr[2];
          edge_n2B1 = treenodes[n2].edge[1];
          edge_n2B2 = treenodes[n2].edge[2];
        } else if (treenodes[n2].nbr[1] == n1) {
          B1 = treenodes[n2].nbr[0];
          B2 = treenodes[n2].nbr[2];
          edge_n2B1 = treenodes[n2].edge[0];
          edge_n2B2 = treenodes[n2].edge[2];
        } else {
          B1 = treenodes[n2].nbr[0];
          B2 = treenodes[n2].nbr[1];
          edge_n2B1 = treenodes[n2].edge[0];
          edge_n2B2 = treenodes[n2].edge[1];
        }

        if (endpt1 == n2 || endpt2 == n2) {
          if (endpt1 == B2 || endpt2 == B2) {
            std::swap(B1, B2);
            std::swap(edge_n2B1, edge_n2B2);
          }

          bool route_ok = updateRouteType1(netID,
                                           treenodes,
                                           n2,
                                           B1,
                                           B2,
                                           E2x,
                                           E2y,
                                           treeedges,
                                           edge_n2B1,
                                           edge_n2B2);
          if (!route_ok) {
            return;
          }

          treenodes[n2].x = E2x;
          treenodes[n2].y = E2y;
        } else {
          const int D1 = endpt1;
          const int D2 = endpt2;
          const int edge_D1D2 = corr_edge_[E2y][E2x];

          bool route_ok = updateRouteType2(netID,
                                           treenodes,
                                           n2,
                                           B1,
                                           B2,
                                           D1,
                                           D2,
                                           E2x,
                                           E2y,
                                           treeedges,
                                           edge_n2B1,
                                           edge_n2B2,
                                           edge_D1D2);
          if (!route_ok) {
            return;
          }

          treenodes[n2].x = E2x;
          treenodes[n2].y = E2y;
        }
      }

      if (treeedges[edge_n1n2].route.type == RouteType::MazeRoute) {
        treeedges[edge_n1n2].route.grids.clear();
      }
      treeedges[edge_n1n2].route.grids.resize(cnt_n1n2);
      treeedges[edge_n1n2].route.type = RouteType::MazeRoute;
      treeedges[edge_n1n2].route.routelen = cnt_n1n2 - 1;
      treeedges[edge_n1n2].len = abs(E1x - E2x) + abs(E1y - E2y);

      for (int i = 0; i < cnt_n1n2; i++) {
        treeedges[edge_n1n2].route.grids[i].x = grids[i].x;
        treeedges[edge_n1n2].route.grids[i].y = grids[i].y;
      }

      FrNet* net = nets_[netID];
      const int8_t edgeCost = net->getEdgeCost();
      for (int i = 0; i < cnt_n1n2 - 1; i++) {
        if (grids[i].x == grids[i + 1].x) {
          const int min_y = std::min(grids[i].y, grids[i + 1].y);
          update_usage_v(grids[i].x, min_y, net, edgeCost);
        } else {
          const int min_x = std::min(grids[i].x, grids[i + 1].x);
          update_usage_h(min_x, grids[i].y, net, edgeCost);
        }
      }
    }
  };

  for (int batch_begin = 0; batch_begin < static_cast<int>(net_ids_.size());
       batch_begin += batch_size) {
    const int batch_end
        = std::min(static_cast<int>(net_ids_.size()), batch_begin + batch_size);
    std::vector<std::vector<UsageUpdate>> updates(
        static_cast<size_t>(batch_end - batch_begin));

#pragma omp parallel for schedule(static) num_threads(threads)
    for (int idx = batch_begin; idx < batch_end; ++idx) {
      auto& net_updates = updates[static_cast<size_t>(idx - batch_begin)];
      net_updates.reserve(256);
      UsageUpdateCallbacks cb;
      cb.ctx = &net_updates;
      cb.updateH = recordUsageH;
      cb.updateV = recordUsageV;
      route_one(idx, &cb);
    }

    for (int idx = batch_begin; idx < batch_end; ++idx) {
      const auto& net_updates = updates[static_cast<size_t>(idx - batch_begin)];
      for (const UsageUpdate& upd : net_updates) {
        if (upd.vertical) {
          graph2d_.updateUsageV(upd.x, upd.y, upd.net, upd.delta);
        } else {
          graph2d_.updateUsageH(upd.x, upd.y, upd.net, upd.delta);
        }
      }
    }
  }
#endif
}

}  // namespace grt::newgr

