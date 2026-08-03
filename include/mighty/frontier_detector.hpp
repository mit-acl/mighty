// Wavefront Frontier Detection (Keidar & Kaminka 2014).
//
// A "frontier cell" is a known-FREE cell adjacent to at least one UNKNOWN cell.
// This detector runs two BFS passes seeded at the robot pose:
//
//   Pass A — outer BFS over reachable known-free cells from the robot.
//            Frontier cells are tagged en passant.
//   Pass B — inner BFS clusters frontier cells via 8-connectivity.
//
// The two-pass structure visits only the reachable known-free region of the
// map, so cost is O(reachable cells) per cycle — much cheaper than naive
// scans, and naturally ignores known-free regions disconnected from the robot.
//
// Stateless: all tunables are passed in via FrontierDetectorParams.

#pragma once

#include <Eigen/Core>
#include <memory>
#include <vector>

#include "mighty/occ_grid_2d.hpp"

class VisitedMap;

struct FrontierDetectorParams {
  int    cluster_min_cells       = 6;   // ~0.24 m² at 0.2 m res
  int    border_margin_cells     = 2;   // exclude N-cell border ring from seeds
  int    obstacle_clearance_cells = 1;  // disqualify free cells within N cells
                                        // of an occupied cell (0 = disabled).
                                        // Filters noisy frontiers hugging walls.
  double robot_snap_radius_m     = 1.0; // spiral search if robot is not on FREE

  // Max run of consecutive UNKNOWN cells the reachability BFS may step through
  // to treat two known-free regions as connected (0 = strictly free-only, the
  // textbook WFD behavior). Sensors have a near-field blind zone: the mapper
  // clears a small cylinder around the robot, but the first ray-cleared cells
  // sit farther out, so the robot's own free pocket is ringed by UNKNOWN and
  // is not 8-connected to anything it can see. A free-only BFS then reaches
  // only that pocket, the sole cluster it finds is centred on the robot, the
  // dwell check retires it as VISITED within a second, and exploration
  // deadlocks before the robot ever moves. Bridging a couple of cells closes
  // that seam. OCCUPIED cells always block, so this never leaks through walls,
  // and bridged UNKNOWN cells are never tagged as frontier seeds themselves.
  int    unknown_bridge_cells    = 0;

  // Optional axis-aligned bounding box that confines exploration. When
  // bounds_enabled is true, frontier seeds outside [min_x,max_x]×[min_y,max_y]
  // are dropped, so the robot never receives exploration goals outside the
  // box. The BFS still expands through the entire reachable known-free area
  // (so the box does not artificially split clusters), only the per-cell
  // seed-tag is gated.
  bool   bounds_enabled = false;
  double bounds_min_x   = 0.0;
  double bounds_max_x   = 0.0;
  double bounds_min_y   = 0.0;
  double bounds_max_y   = 0.0;
};

struct FrontierCluster {
  Eigen::Vector2d centroid;             // world frame
  std::vector<Eigen::Vector2d> cells;   // world frame, one entry per cell
  int             size_cells = 0;
  double          size_m2    = 0.0;
  Eigen::Vector2d aabb_min = Eigen::Vector2d::Zero();
  Eigen::Vector2d aabb_max = Eigen::Vector2d::Zero();
};

class FrontierDetector {
 public:
  explicit FrontierDetector(FrontierDetectorParams p) : params_(p) {}

  /** @brief Run WFD on `grid` seeded at `robot_xy`.
   *  @param grid        Local 2D occupancy window from the mapper.
   *  @param robot_xy    Robot position in world frame.
   *  @param visited_map Optional persistent "previously observed" bitmap.
   *                     UNKNOWN neighbors that are already visited do *not*
   *                     count as frontier-generating — they are stale slid-out
   *                     cells, not genuine unexplored area.
   *  @return Vector of frontier clusters in world frame, sorted by descending
   *          size_cells. Empty if the robot pose can't be snapped to a free
   *          cell or no frontiers exist.
   */
  std::vector<FrontierCluster> detect(
      const OccGrid2D& grid,
      const Eigen::Vector2d& robot_xy,
      const VisitedMap* visited_map = nullptr) const;

  /** @brief Re-validate a persistent frontier record against the current grid.
   *  Returns true iff `world_xy` lies on — or within `search_radius_cells` of —
   *  a connected frontier-cell component of at least `min_cells` cells, using
   *  the SAME per-cell definition as detect(). The search radius absorbs the
   *  EMA drift of a record's centroid so a frontier that merely shifted a little
   *  is not falsely retired.
   *
   *  Lets the frontier manager retire records whose unknown area has since been
   *  explored (or shrank below the detection threshold, or is a permanently
   *  occluded speck) so they stop lingering as phantom frontiers on known cells.
   *  Walk early-exits once `min_cells` is reached, so cost is O(min_cells).
   */
  bool isStillFrontier(const OccGrid2D& grid, const Eigen::Vector2d& world_xy,
                       int min_cells, int search_radius_cells,
                       const VisitedMap* visited_map = nullptr) const;

  const FrontierDetectorParams& params() const { return params_; }

 private:
  // Shared per-cell frontier predicate: a known-FREE cell (off the border ring)
  // with at least one in-bounds, not-yet-visited UNKNOWN neighbor, that is not
  // hugging an obstacle and lies inside the optional exploration bounds. Single
  // source of truth for both detect() seed-tagging and isStillFrontier().
  bool isFrontierCell(const OccGrid2D& grid, int cx, int cy,
                      const VisitedMap* visited_map) const;

  FrontierDetectorParams params_;
};
