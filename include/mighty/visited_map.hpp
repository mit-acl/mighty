// Persistent tristate occupancy map for frontier exploration + visualization.
//
// The mapper publishes a *sliding* 2D occupancy window centered on the robot.
// When the robot moves, cells fall off the window and are re-marked UNKNOWN
// in the next message. This causes two visible problems:
//
//   1) The WFD detector would happily re-create frontiers along the seam
//      where local UNKNOWN re-appears, even though the area was already
//      fully observed in the past.
//   2) RViz visualization flickers UNKNOWN over revisited cells until the
//      lidar re-observes them.
//
// VisitedMap solves both with a fixed-size persistent grid, world-frame
// indexed, that records the last-known tristate value (-1 / 0 / 100, matching
// nav_msgs/OccupancyGrid codes) for every cell ever observed in the mission.
//
//   - The frontier detector consults isVisitedWorld(): if an UNKNOWN neighbor
//     is already in the persistent map, it's not really unknown — it just
//     slid out of the local window — so it does not count as a frontier.
//   - publishVisitedMap() in MIGHTY_NODE publishes the raw tristate buffer as
//     a nav_msgs/OccupancyGrid on /exploration/visited_map. RViz can layer
//     this behind the sliding occ_2d so revisited areas keep their old
//     FREE/OCCUPIED colors instead of flickering UNKNOWN.
//
// Memory: at 0.15 m resolution, a 100 x 100 m area is ~444 KB. The absorb()
// loop iterates the *local* grid (not this map), so per-frame CPU stays
// constant even if this map is sized in km². See Computational Cost section
// of the design plan for the full table.
//
// Threading: not internally synchronized. All access must come from the same
// callback group (in MIGHTY, occ2DCallback owns it).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "mighty/occ_grid_2d.hpp"

class VisitedMap {
 public:
  // Cell encoding (matches nav_msgs/OccupancyGrid).
  static constexpr int8_t kUnknown  = -1;
  static constexpr int8_t kFree     =  0;
  static constexpr int8_t kOccupied = 100;

  VisitedMap() = default;

  /** @brief Construct a fixed-size grid centered on (cx, cy) in world frame.
   *  @param cx,cy   Center of the persistent area, in world coordinates.
   *  @param width_m, height_m  Total extent of the grid.
   *  @param resolution         Cell size, normally matched to the local grid.
   */
  VisitedMap(double cx, double cy, double width_m, double height_m,
             double resolution)
      : origin_x_(cx - 0.5 * width_m),
        origin_y_(cy - 0.5 * height_m),
        resolution_(resolution),
        inv_resolution_(1.0 / resolution),
        width_(static_cast<int>(std::round(width_m / resolution))),
        height_(static_cast<int>(std::round(height_m / resolution))),
        data_(static_cast<size_t>(width_) * height_, kUnknown) {}

  /** @brief Copy every observed (non-unknown) cell of `grid` into the
   *  persistent buffer, recording its actual FREE/OCCUPIED state.
   *
   *  Observed cells overwrite whatever was previously stored at the same
   *  world location, so the persistent map always reflects the most-recent
   *  observation. Cells outside the persistent map's extent are silently
   *  dropped (worldToGrid returns false).
   */
  void absorb(const OccGrid2D& grid) {
    if (data_.empty()) return;
    const int W = grid.width();
    const int H = grid.height();
    for (int iy = 0; iy < H; ++iy) {
      for (int ix = 0; ix < W; ++ix) {
        if (grid.isUnknown(ix, iy)) continue;  // skip unobserved local cells
        const int8_t v = grid.isOccupied(ix, iy) ? kOccupied : kFree;
        double wx, wy;
        grid.gridToWorld(ix, iy, wx, wy);
        setStateWorld(wx, wy, v);
      }
    }
  }

  /** @brief Return the persistent state at (wx, wy), or kUnknown if outside
   *  the persistent map or never observed. */
  int8_t getStateWorld(double wx, double wy) const {
    if (data_.empty()) return kUnknown;
    int ix, iy;
    if (!worldToGrid(wx, wy, ix, iy)) return kUnknown;
    return data_[static_cast<size_t>(iy) * width_ + ix];
  }

  /** @brief Was this cell ever observed (FREE or OCCUPIED) before? Used by
   *  the frontier detector to filter stale-unknown neighbors. */
  bool isVisitedWorld(double wx, double wy) const {
    return getStateWorld(wx, wy) != kUnknown;
  }

  /** @brief Write `v` (kUnknown / kFree / kOccupied) at the given world coords.
   *  No-op if the coords fall outside the persistent map. */
  void setStateWorld(double wx, double wy, int8_t v) {
    int ix, iy;
    if (!worldToGrid(wx, wy, ix, iy)) return;
    data_[static_cast<size_t>(iy) * width_ + ix] = v;
  }

  /** @brief Merge a peer robot's visited map into this one.  Only overwrites
   *  cells that are still UNKNOWN locally — i.e. never overrides this robot's
   *  own observations.  Handles different grid origins/dimensions by converting
   *  through world coordinates.
   *
   *  @param data       Row-major tristate buffer (kUnknown / kFree / kOccupied).
   *  @param w,h        Grid width/height in cells.
   *  @param origin_x,origin_y  World-frame lower-left corner.
   *  @param resolution Cell size in meters.
   */
  void mergeFrom(const int8_t* data, int w, int h,
                 double origin_x, double origin_y, double resolution) {
    if (data_.empty()) return;
    for (int iy = 0; iy < h; ++iy) {
      for (int ix = 0; ix < w; ++ix) {
        const int8_t v = data[iy * w + ix];
        if (v == kUnknown) continue;  // nothing to learn
        const double wx = origin_x + (ix + 0.5) * resolution;
        const double wy = origin_y + (iy + 0.5) * resolution;
        if (getStateWorld(wx, wy) == kUnknown) {
          setStateWorld(wx, wy, v);
        }
      }
    }
  }

  void clear() { std::fill(data_.begin(), data_.end(), kUnknown); }

  // Geometry accessors (used by visualization & frontier filter).
  double originX() const { return origin_x_; }
  double originY() const { return origin_y_; }
  double resolution() const { return resolution_; }
  int width() const { return width_; }
  int height() const { return height_; }
  bool empty() const { return data_.empty(); }
  const std::vector<int8_t>& data() const { return data_; }

 private:
  bool worldToGrid(double wx, double wy, int& ix, int& iy) const {
    ix = static_cast<int>(std::floor((wx - origin_x_) * inv_resolution_));
    iy = static_cast<int>(std::floor((wy - origin_y_) * inv_resolution_));
    return ix >= 0 && ix < width_ && iy >= 0 && iy < height_;
  }

  double origin_x_ = 0.0;
  double origin_y_ = 0.0;
  double resolution_ = 0.0;
  double inv_resolution_ = 0.0;
  int width_ = 0;
  int height_ = 0;
  std::vector<int8_t> data_;  // -1 unknown / 0 free / 100 occupied
};
