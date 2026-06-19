/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

/**
 * @file hgp_planner.h
 * @brief HGPPlanner
 */
#pragma once

#include <cmath>

#include <hgp/data_type.hpp>
#include <hgp/graph_search.hpp>
#include <hgp/map_util.hpp>

/** @brief Global path planner using A* and JPS on a 3D voxel map.
 *
 *  HGPPlanner wraps GraphSearch and provides path planning with line-of-sight
 *  shortcutting, collinear/corner point removal, and heat-aware Laplacian
 *  smoothing. Supports both 3D UAV and 2D ground robot modes.
 */
class HGPPlanner {
 public:
  /** @brief Construct the planner with algorithm choice and dynamic constraints.
   *  @param global_planner Name of the search algorithm (e.g., "astar_heat", "sjps").
   *  @param verbose Enable verbose debug output.
   *  @param v_max Maximum velocity constraint.
   *  @param a_max Maximum acceleration constraint.
   *  @param j_max Maximum jerk constraint.
   *  @param hgp_timeout_duration_ms Planning timeout in milliseconds.
   *  @param w_unknown Cost weight for traversing unknown cells.
   *  @param w_align Alignment penalty weight for start-direction bias.
   *  @param decay_len_cells E-folding distance for alignment decay in cells.
   *  @param w_side Side tie-break weight for handedness preference.
   *  @param los_cells Inflation radius in cells for line-of-sight shortcutting.
   *  @param min_len Minimum edge length for path simplification.
   *  @param min_turn Minimum turn angle (degrees) for path simplification.
   */
  HGPPlanner(std::string global_planner, bool verbose, double v_max, double a_max, double j_max,
             int hgp_timeout_duration_ms, double w_unknown, double w_align, double decay_len_cells,
             double w_side, int los_cells = 3, double min_len = 0.5, double min_turn = 10.0);

  /** @brief Set the voxel map utility used for collision checking.
   *  @param map_util Shared pointer to a 3D voxel map utility.
   */
  void setMapUtil(const std::shared_ptr<mighty::MapUtil<3>>& map_util);
  /**
   * @brief Status of the planner
   *
   * 0 --- exit normally;
   * -1 --- no path found;
   * 1, 2 --- start or goal is not free.
   */
  int status();
  /** @brief Get the post-processed path.
   *  @return Cleaned-up waypoint path.
   */
  vec_Vecf<3> getPath();

  /** @brief Get the raw path directly from the graph search.
   *  @return Unmodified waypoint path.
   */
  vec_Vecf<3> getRawPath();

  /** @brief Set the planning timeout duration.
   *  @param timeout_duration_ms Timeout in milliseconds.
   */
  void setTimeOutDurationMs(int timeout_duration_ms);

  /** @brief Remove collinear intermediate waypoints from a path.
   *  @param path Input waypoint path.
   *  @return Path with redundant collinear points removed.
   */
  vec_Vecf<3> removeLinePts(const vec_Vecf<3>& path);

  /** @brief Remove unnecessary corner waypoints from a path.
   *  @param path Input waypoint path.
   *  @return Simplified path.
   */
  vec_Vecf<3> removeCornerPts(const vec_Vecf<3>& path);

  /** @brief Apply all path cleanup steps (collinear and corner removal) in place.
   *  @param path Waypoint path to clean up.
   */
  void cleanUpPath(vec_Vecf<3>& path);

  /** @brief Run the global path planner from start to goal.
   *  @param start Start position.
   *  @param start_vel Start velocity for alignment heuristic.
   *  @param goal Goal position.
   *  @param final_g Output cost of the found path.
   *  @param current_time Current simulation/wall time.
   *  @param eps Heuristic weight for A* search.
   *  @return True if a valid path was found.
   */
  bool plan(const Vecf<3>& start, const Vecf<3>& start_vel, const Vecf<3>& goal, double& final_g,
            double current_time, decimal_t eps = 1);

  /** @brief Get the nodes in the open set after planning.
   *  @return Positions of open-set nodes.
   */
  vec_Vecf<3> getOpenSet() const;

  /** @brief Get the nodes in the closed set after planning.
   *  @return Positions of closed-set nodes.
   */
  vec_Vecf<3> getCloseSet() const;

  /** @brief Get all explored nodes after planning.
   *  @return Positions of all explored nodes.
   */
  vec_Vecf<3> getAllSet() const;

  /** @brief Get the total initial-guess planning time.
   *  @return Time in seconds.
   */
  double getInitialGuessPlanningTime();

  /** @brief Get the static JPS planning time.
   *  @return Time in seconds.
   */
  double getStaticJPSPlanningTime();

  /** @brief Get the path-checking time.
   *  @return Time in seconds.
   */
  double getCheckPathTime();

  /** @brief Get the dynamic A* planning time.
   *  @return Time in seconds.
   */
  double getDynamicAstarTime();

  /** @brief Get the path recovery time.
   *  @return Time in seconds.
   */
  double getRecoverPathTime();

  /** @brief Update the maximum velocity used by the planner.
   *  @param v_max New maximum velocity.
   */
  void updateVmax(double v_max);

  /** @brief Set the maximum number of node expansions for graph search.
   *  @param max_expand Maximum node expansions.
   */
  void setMaxExpand(int max_expand) { max_expand_ = max_expand; }

  /** @brief Shorten a path using line-of-sight checks with inflated capsule collision tests.
   *  @param in Input waypoint path.
   *  @param inflate_radius_cells Inflation radius in voxel cells.
   *  @return Shortened path.
   */
  vec_Vecf<3> shortCutByLoS(const vec_Vecf<3>& in, int inflate_radius_cells) const;

  /** @brief Check line-of-sight between two points using an inflated capsule collision test.
   *  @param a Start point.
   *  @param b End point.
   *  @param inflate_radius_cells Inflation radius in voxel cells.
   *  @return True if the line segment is collision-free.
   */
  bool lineOfSightCapsule(const Vecf<3>& a, const Vecf<3>& b, int inflate_radius_cells) const;

  /** @brief Merge consecutive edges shorter than a minimum length.
   *  @param in Input waypoint path.
   *  @param min_len Minimum edge length.
   *  @return Path with short edges collapsed.
   */
  vec_Vecf<3> collapseShortEdges(const vec_Vecf<3>& in, double min_len) const;

  /** @brief Remove waypoints that create turns below a minimum angle threshold.
   *  @param in Input waypoint path.
   *  @param min_turn_deg Minimum turn angle in degrees.
   *  @param min_seg_len Minimum segment length to consider.
   *  @return Filtered path.
   */
  vec_Vecf<3> angleSpacingFilter(const vec_Vecf<3>& in, double min_turn_deg,
                                 double min_seg_len) const;

 protected:
  // Assume using 3D voxel map for all 2d and 3d planning
  std::shared_ptr<mighty::VoxelMapUtil> map_util_;
  // The planner
  std::shared_ptr<mighty::GraphSearch> graph_search_;
  // Raw path from planner
  vec_Vecf<3> raw_path_;
  // Modified path for future usage
  vec_Vecf<3> path_;
  // Flag indicating the success of planning
  int status_ = 0;
  // Enabled for printing info
  bool planner_verbose_;

  // Parameters
  std::string global_planner_;

  // For benchmarking time recording
  double global_planning_time_ = 0.0;
  double hgp_static_jps_time_ = 0.0;
  double hgp_check_path_time_ = 0.0;
  double hgp_dynamic_astar_time_ = 0.0;
  double hgp_recover_path_time_ = 0.0;

  // time out duration
  int hgp_timeout_duration_ms_ = 1000;

  // max node expansion
  int max_expand_ = 10000;

  // max values
  double v_max_;
  double a_max_;
  double j_max_;

  // alignment
  double w_align_ = 60.0;          // weight for alignment with start velocity
  double decay_len_cells_ = 20.0;  // decay length in cells for alignment
  double w_side_ = 0.2;            // weight for side tie-break

  // unknown cell cost weight
  double w_unknown_ = 1.0;

  // ESDF for ground robot A* cost
  std::shared_ptr<const EsdfGrid2D> esdf_grid_;
  double esdf_weight_astar_ = 0.0;
  double esdf_d_safe_astar_ = 0.0;

  // Corridor-center corner snap (ground robot only). Params are set via
  // setCornerSnapParams; snap_enabled_ defaults to false so UAV flow is
  // untouched.
  bool   snap_enabled_ = false;
  double snap_corner_angle_rad_ = 0.0;
  double snap_clearance_threshold_m_ = 0.0;
  double snap_max_ascent_m_ = 0.0;
  double snap_ascent_step_m_ = 0.0;
  int    snap_ascent_max_iters_ = 0;

  // LOS post processing
  int los_cells_ = 3;       // number of cells for inflation in LoS
  double min_len_ = 0.5;    // minimum length of edges
  double min_turn_ = 10.0;  // minimum turn angle in degrees

 public:
  /// Whether 2D ground robot planning mode is active.
  bool is_2d_mode_{false};

  /** @brief Enable or disable 2D ground robot planning mode.
   *  @param enabled True to restrict planning to 2D.
   */
  void set2DMode(bool enabled) { is_2d_mode_ = enabled; }

  /** @brief Set ESDF grid for distance-based A* cost (ground robot only). */
  void setEsdfGrid(std::shared_ptr<const EsdfGrid2D> grid, double weight, double d_safe) {
    esdf_grid_ = grid;
    esdf_weight_astar_ = weight;
    esdf_d_safe_astar_ = d_safe;
  }

  /** @brief Configure corridor-center corner snap post-processing (ground robot only).
   *
   *  At each sharp corner waypoint in the processed path, performs ESDF gradient
   *  ascent until either (a) the local clearance exceeds @p clearance_threshold
   *  (open-field / corridor-center reached) or (b) the accumulated displacement
   *  exceeds @p max_ascent. Corners already ≥ @p clearance_threshold from any
   *  obstacle are left alone. No-op when @p enabled is false or the ESDF grid
   *  is unset.
   */
  void setCornerSnapParams(bool enabled, double corner_angle_deg,
                           double clearance_threshold, double max_ascent,
                           double ascent_step, int ascent_max_iters) {
    snap_enabled_ = enabled;
    snap_corner_angle_rad_ = corner_angle_deg * M_PI / 180.0;
    snap_clearance_threshold_m_ = clearance_threshold;
    snap_max_ascent_m_ = max_ascent;
    snap_ascent_step_m_ = ascent_step;
    snap_ascent_max_iters_ = ascent_max_iters;
  }

  /** @brief Snap sharp corners in @p path toward local ESDF max-clearance points.
   *
   *  See setCornerSnapParams for behavior. Mutates @p path in place. Preserves
   *  z coordinates (operates in the xy plane only).
   */
  void snapCornersToClearance(vec_Vecf<3>& path) const;

  /// Whether all path smoothing is disabled (raw A* path is used directly).
  bool disable_all_smoothing_{false};

  /** @brief Enable or disable all path smoothing.
   *  @param enabled True to skip all smoothing and return the raw search path.
   */
  void setDisableAllSmoothing(bool enabled) { disable_all_smoothing_ = enabled; }

  /// Whether line-of-sight shortcutting is skipped in favor of heat-aware smoothing.
  bool skip_path_smoothing_{false};

  /** @brief Enable or disable skipping of LoS-based path shortcutting.
   *  @param enabled True to skip LoS shortcutting and use Laplacian smoothing instead.
   */
  void setSkipPathSmoothing(bool enabled) { skip_path_smoothing_ = enabled; }

  /// Number of iterations for heat-aware Laplacian smoothing.
  int smooth_iterations_{50};
  /// Blending weight for heat-aware Laplacian smoothing.
  double smooth_alpha_{0.3};

  /** @brief Set parameters for heat-aware Laplacian path smoothing.
   *  @param iterations Number of smoothing iterations.
   *  @param alpha Blending weight in [0,1]; higher values produce smoother paths.
   */
  void setSmoothParams(int iterations, double alpha) {
    smooth_iterations_ = iterations;
    smooth_alpha_ = alpha;
  }

  /// Maximum distance between consecutive path waypoints for 2D resampling [m].
  double max_dist_vertexes_2d_{1.0};

  /** @brief Set max distance between waypoints for 2D path resampling. */
  void setMaxDistVertexes2D(double d) { max_dist_vertexes_2d_ = d; }

  /** @brief Smooth a path using heat-aware Laplacian smoothing.
   *
   *  Iteratively moves interior waypoints toward the midpoint of their
   *  neighbors while biasing away from high-heat (obstacle proximity) regions.
   *
   *  @param path Input waypoint path.
   *  @param iterations Number of smoothing iterations.
   *  @param alpha Blending weight per iteration.
   *  @return Smoothed path.
   */
  vec_Vecf<3> smoothPathHeatAware(const vec_Vecf<3>& path, int iterations, double alpha) const;
};
