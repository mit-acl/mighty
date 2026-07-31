#include "mighty/frontier_detector.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <queue>
#include <unordered_set>

#include "mighty/visited_map.hpp"

namespace {

// 8-connected neighbor offsets.
constexpr int kDx8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
constexpr int kDy8[8] = {-1, -1, -1,  0, 0,  1, 1, 1};

inline size_t flatIdx(int ix, int iy, int width) {
  return static_cast<size_t>(iy) * width + ix;
}

// Spiral search outward from (cx, cy) for the nearest FREE cell within
// `max_radius_cells`. Returns true on success and writes the result into
// (*ox, *oy). Naive O(R²) scan — fine for the small radii we use here.
bool snapToFree(const OccGrid2D& grid, int cx, int cy, int max_radius_cells,
                int* ox, int* oy) {
  if (grid.isFree(cx, cy)) {
    *ox = cx;
    *oy = cy;
    return true;
  }
  for (int r = 1; r <= max_radius_cells; ++r) {
    for (int dy = -r; dy <= r; ++dy) {
      for (int dx = -r; dx <= r; ++dx) {
        // Only consider cells on the ring of radius r.
        if (std::max(std::abs(dx), std::abs(dy)) != r) continue;
        int nx = cx + dx;
        int ny = cy + dy;
        if (grid.isFree(nx, ny)) {
          *ox = nx;
          *oy = ny;
          return true;
        }
      }
    }
  }
  return false;
}

}  // namespace

std::vector<FrontierCluster> FrontierDetector::detect(
    const OccGrid2D& grid, const Eigen::Vector2d& robot_xy,
    const VisitedMap* visited_map) const {
  std::vector<FrontierCluster> out;

  const int W = grid.width();
  const int H = grid.height();
  if (W <= 0 || H <= 0) return out;

  // ---- Snap robot pose to nearest FREE cell ----
  int rx, ry;
  grid.worldToGrid(robot_xy.x(), robot_xy.y(), rx, ry);
  const int snap_cells = std::max(
      0, static_cast<int>(std::ceil(params_.robot_snap_radius_m
                                    / grid.resolution())));
  int seed_x = -1, seed_y = -1;
  if (!snapToFree(grid, rx, ry, snap_cells, &seed_x, &seed_y)) {
    // Robot not in or near any free cell — bail out silently. Caller may log.
    return out;
  }

  // ---- Pass A: BFS over known-free cells, tag frontier seeds ----
  const size_t N = static_cast<size_t>(W) * H;
  std::vector<uint8_t> visited(N, 0);
  std::vector<uint8_t> is_frontier(N, 0);
  // Remaining UNKNOWN-bridging budget each cell was reached with. Stepping onto
  // a known-free cell refills it; stepping onto an UNKNOWN cell spends one.
  // With unknown_bridge_cells == 0 this degenerates to the free-only BFS.
  const int bridge = std::max(0, params_.unknown_bridge_cells);
  std::vector<int8_t> budget(N, 0);

  std::queue<std::pair<int, int>> q;
  q.push({seed_x, seed_y});
  visited[flatIdx(seed_x, seed_y, W)] = 1;
  budget[flatIdx(seed_x, seed_y, W)] = static_cast<int8_t>(std::min(bridge, 127));

  while (!q.empty()) {
    auto [cx, cy] = q.front();
    q.pop();

    // Tag this cell as a frontier seed if it matches the shared per-cell
    // frontier predicate. Keeping the definition in one place (isFrontierCell)
    // guarantees record re-validation (isStillFrontier) can never drift from
    // what detection considers a frontier. Bridged UNKNOWN cells fail the
    // predicate's isFree() gate, so they can never become seeds.
    const size_t cidx = flatIdx(cx, cy, W);
    const bool cur_free = grid.isFree(cx, cy);
    if (cur_free && isFrontierCell(grid, cx, cy, visited_map)) {
      is_frontier[cidx] = 1;
    }

    // Expand BFS through known-free cells, plus short UNKNOWN seams when
    // bridging is enabled. OCCUPIED always blocks, so bridging can only
    // reconnect free space the sensor simply hasn't filled in yet.
    const int cur_budget = cur_free ? bridge : budget[cidx];
    for (int i = 0; i < 8; ++i) {
      const int nx = cx + kDx8[i];
      const int ny = cy + kDy8[i];
      if (!grid.inBounds(nx, ny)) continue;
      if (grid.isOccupied(nx, ny)) continue;
      const size_t nidx = flatIdx(nx, ny, W);
      const bool nfree = grid.isFree(nx, ny);
      const int nbudget = nfree ? bridge : cur_budget - 1;
      if (nbudget < 0) continue;  // out of bridging budget — seam too wide
      // Revisit only when we arrive with strictly more budget to spend, so a
      // cell first reached mid-seam can still expand once free space is found.
      if (visited[nidx] && budget[nidx] >= nbudget) continue;
      visited[nidx] = 1;
      budget[nidx] = static_cast<int8_t>(std::min(nbudget, 127));
      q.push({nx, ny});
    }
  }

  // ---- Pass B: BFS over frontier seeds to form clusters ----
  std::vector<uint8_t> assigned(N, 0);
  for (int sy = 0; sy < H; ++sy) {
    for (int sx = 0; sx < W; ++sx) {
      const size_t sidx = flatIdx(sx, sy, W);
      if (!is_frontier[sidx] || assigned[sidx]) continue;

      FrontierCluster cluster;
      double sum_wx = 0.0, sum_wy = 0.0;
      double min_wx = std::numeric_limits<double>::infinity();
      double min_wy = std::numeric_limits<double>::infinity();
      double max_wx = -std::numeric_limits<double>::infinity();
      double max_wy = -std::numeric_limits<double>::infinity();

      std::queue<std::pair<int, int>> cq;
      cq.push({sx, sy});
      assigned[sidx] = 1;

      while (!cq.empty()) {
        auto [cx, cy] = cq.front();
        cq.pop();

        double wx, wy;
        grid.gridToWorld(cx, cy, wx, wy);
        cluster.cells.emplace_back(wx, wy);
        sum_wx += wx;
        sum_wy += wy;
        min_wx = std::min(min_wx, wx);
        min_wy = std::min(min_wy, wy);
        max_wx = std::max(max_wx, wx);
        max_wy = std::max(max_wy, wy);

        for (int i = 0; i < 8; ++i) {
          const int nx = cx + kDx8[i];
          const int ny = cy + kDy8[i];
          if (!grid.inBounds(nx, ny)) continue;
          const size_t nidx = flatIdx(nx, ny, W);
          if (assigned[nidx]) continue;
          if (!is_frontier[nidx]) continue;
          assigned[nidx] = 1;
          cq.push({nx, ny});
        }
      }

      cluster.size_cells = static_cast<int>(cluster.cells.size());
      if (cluster.size_cells < params_.cluster_min_cells) continue;

      cluster.size_m2 = cluster.size_cells
                        * grid.resolution() * grid.resolution();
      cluster.centroid = Eigen::Vector2d(sum_wx / cluster.size_cells,
                                         sum_wy / cluster.size_cells);
      cluster.aabb_min = Eigen::Vector2d(min_wx, min_wy);
      cluster.aabb_max = Eigen::Vector2d(max_wx, max_wy);
      out.push_back(std::move(cluster));
    }
  }

  // Sort by descending size for stable downstream behavior.
  std::sort(out.begin(), out.end(),
            [](const FrontierCluster& a, const FrontierCluster& b) {
              return a.size_cells > b.size_cells;
            });

  return out;
}

bool FrontierDetector::isFrontierCell(const OccGrid2D& grid, int cx, int cy,
                                      const VisitedMap* visited_map) const {
  // Must be a known-FREE cell (the detector's outer BFS only ever walks free
  // cells; enforce it here so the predicate is safe to call on any cell).
  if (!grid.isFree(cx, cy)) return false;

  // Not on the border ring. Cells near the window edge may have just slid in
  // and not been observed yet; their out-of-bounds neighbors would otherwise
  // look UNKNOWN and spuriously flag the whole edge as a frontier.
  const int m = params_.border_margin_cells;
  if (cx < m || cx >= grid.width() - m ||
      cy < m || cy >= grid.height() - m) {
    return false;
  }

  // At least one *in-bounds* UNKNOWN 8-neighbor. Out-of-bounds neighbors are
  // "outside the mapper window", handled by border_margin_cells, not counted
  // as unknown here. VisitedMap filter: an UNKNOWN neighbor already in the
  // persistent visited bitmap is stale slid-out area, not genuine unexplored
  // space, so it does not generate a frontier.
  bool has_unknown_nbr = false;
  for (int i = 0; i < 8; ++i) {
    const int nx = cx + kDx8[i];
    const int ny = cy + kDy8[i];
    if (!grid.inBounds(nx, ny)) continue;
    if (!grid.isUnknown(nx, ny)) continue;
    if (visited_map != nullptr) {
      double wx_n, wy_n;
      grid.gridToWorld(nx, ny, wx_n, wy_n);
      if (visited_map->isVisitedWorld(wx_n, wy_n)) continue;
    }
    has_unknown_nbr = true;
    break;
  }
  if (!has_unknown_nbr) return false;

  // Reject cells hugging an obstacle (usually wall-thickness slop / raycast
  // endpoints, not worth visiting).
  if (params_.obstacle_clearance_cells > 0) {
    const int rcells = params_.obstacle_clearance_cells;
    for (int dy = -rcells; dy <= rcells; ++dy) {
      for (int dx = -rcells; dx <= rcells; ++dx) {
        if (dx == 0 && dy == 0) continue;
        const int nx = cx + dx;
        const int ny = cy + dy;
        if (!grid.inBounds(nx, ny)) continue;
        if (grid.isOccupied(nx, ny)) return false;
      }
    }
  }

  // Drop cells outside the optional exploration bounds box (world frame).
  if (params_.bounds_enabled) {
    double wx, wy;
    grid.gridToWorld(cx, cy, wx, wy);
    if (wx < params_.bounds_min_x || wx > params_.bounds_max_x ||
        wy < params_.bounds_min_y || wy > params_.bounds_max_y) {
      return false;
    }
  }

  return true;
}

bool FrontierDetector::isStillFrontier(const OccGrid2D& grid,
                                       const Eigen::Vector2d& world_xy,
                                       int min_cells, int search_radius_cells,
                                       const VisitedMap* visited_map) const {
  const int W = grid.width();
  const int H = grid.height();
  if (W <= 0 || H <= 0) return false;
  if (min_cells < 1) min_cells = 1;
  if (search_radius_cells < 0) search_radius_cells = 0;

  int cx, cy;
  grid.worldToGrid(world_xy.x(), world_xy.y(), cx, cy);
  if (!grid.inBounds(cx, cy)) return false;

  // Find a frontier-cell seed within search_radius_cells, scanning outward ring
  // by ring so a centroid nudged off the frontier by EMA blending still binds
  // to its frontier. First match wins (nearest seed).
  int sx = -1, sy = -1;
  for (int r = 0; r <= search_radius_cells && sx < 0; ++r) {
    for (int dy = -r; dy <= r && sx < 0; ++dy) {
      for (int dx = -r; dx <= r; ++dx) {
        if (std::max(std::abs(dx), std::abs(dy)) != r) continue;  // ring only
        const int nx = cx + dx;
        const int ny = cy + dy;
        if (isFrontierCell(grid, nx, ny, visited_map)) {
          sx = nx;
          sy = ny;
          break;
        }
      }
    }
  }
  if (sx < 0) return false;  // no frontier cell nearby -> no longer a frontier

  // BFS the connected frontier-cell component from the seed, early-exiting the
  // moment it reaches min_cells. A local hash set bounds memory: the walk stops
  // at min_cells (still a real frontier) or exhausts a sub-threshold remnant.
  std::unordered_set<size_t> seen;
  std::queue<std::pair<int, int>> q;
  q.push({sx, sy});
  seen.insert(flatIdx(sx, sy, W));
  int count = 0;
  while (!q.empty()) {
    auto [x, y] = q.front();
    q.pop();
    if (++count >= min_cells) return true;  // enough cells -> still a frontier
    for (int i = 0; i < 8; ++i) {
      const int nx = x + kDx8[i];
      const int ny = y + kDy8[i];
      if (!grid.inBounds(nx, ny)) continue;
      const size_t nidx = flatIdx(nx, ny, W);
      if (seen.count(nidx)) continue;
      if (!isFrontierCell(grid, nx, ny, visited_map)) continue;
      seen.insert(nidx);
      q.push({nx, ny});
    }
  }
  return count >= min_cells;  // component smaller than min_cells -> retire it
}
