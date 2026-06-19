/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/StdVector>

#include <pcl/kdtree/kdtree.h>

#include <mighty/lbfgs_solver.hpp>
#include <mighty/utils.hpp>

#include "hgp/hgp_manager.hpp"
#include "hgp/termcolor.hpp"
#include "mighty/mighty_type.hpp"

#include "timer.hpp"

enum { MAP = 0, UNKNOWN_MAP = 1 };
enum { RETURN_LAST_VERTEX = 0, RETURN_INTERSECTION = 1 };

using Vec3 = Eigen::Vector3d;
using Vec3f = Eigen::Matrix<double, 3, 1>;
using MatXd = Eigen::MatrixXd;
using VecXd = Eigen::VectorXd;

// ------------------------------------------
// 2) The aligned‐allocator for local blocks:
// ------------------------------------------
// your PlaneBlock is still:
using PlaneBlock = std::pair<Eigen::Matrix<double, Eigen::Dynamic, 3>, Eigen::VectorXd>;

// replace your old ConstraintBlocks with:

// 1) an aligned‐allocator for a vector of PlaneBlock
using AlignedPlaneBlockVec = std::vector<PlaneBlock, Eigen::aligned_allocator<PlaneBlock>>;

// 2) then your top‐level blocks is a vector of those
using ConstraintBlocks =
    std::vector<AlignedPlaneBlockVec, Eigen::aligned_allocator<AlignedPlaneBlockVec>>;

enum DroneStatus { YAWING = 0, TRAVELING = 1, GOAL_SEEN = 2, GOAL_REACHED = 3 };

/** @brief Core trajectory planner using Hermite spline-based optimization.
 *
 *  Implements the MIGHTY planning algorithm including global path search,
 *  safe corridor computation, local trajectory optimization via L-BFGS,
 *  and a state machine (YAWING -> TRAVELING -> GOAL_SEEN -> GOAL_REACHED).
 */
class MIGHTY {
 public:
  /** @brief Construct the planner with the given parameter set.
   *  @param par Configuration parameters for dynamics, map, and planner.
   */
  MIGHTY(parameters par);

  /** @brief Check whether a replan is needed based on current and goal states.
   *  @param local_state Current robot state.
   *  @param local_G_term Terminal goal state.
   *  @param last_plan_state Last planned state.
   *  @return True if replanning is necessary.
   */
  bool needReplan(const state& local_state, const state& local_G_term,
                  const state& last_plan_state);

  /** @brief Find the planning start state A and its timestamp.
   *  @param A Output start state for the next plan.
   *  @param A_time Output timestamp of the start state.
   *  @param current_time Current wall/sim time.
   *  @param last_replaning_computation_time Computation time of the last replan.
   *  @return True if a valid start state was found.
   */
  bool findAandAtime(state& A, double& A_time, double current_time,
                     double last_replaning_computation_time);

  /** @brief Check if a point is in an occupied voxel.
   *  @param point Query position.
   *  @return True if occupied.
   */
  bool checkIfPointOccupied(const Vec3f& point);

  /** @brief Check if a point is in a free voxel.
   *  @param point Query position.
   *  @return True if free.
   */
  bool checkIfPointFree(const Vec3f& point);

  /** @brief Compute safe flight corridors (convex polytopes) along a global path.
   *  @param global_path Waypoint path for corridor computation.
   *  @param A Planning start state.
   *  @return True if safe corridors were successfully computed.
   */
  bool getSafeCorridor(const vec_Vecf<3>& global_path, const state& A);

  /** @brief Execute one full replan cycle (global path + local trajectory).
   *  @param last_replaning_computation_time Computation time of the previous replan.
   *  @param current_time Current wall/sim time.
   *  @return Tuple of (replan_success, trajectory_changed).
   */
  std::tuple<bool, bool> replan(double last_replaning_computation_time, double current_time);

  /** @brief Begin adaptive k-value estimation for replanning intervals. */
  void startAdaptKValue();

  /** @brief Get the terminal goal state (thread-safe).
   *  @param G_term Output terminal goal.
   */
  void getGterm(state& G_term);

  /** @brief Set the terminal goal state (thread-safe).
   *  @param G_term New terminal goal.
   */
  void setGterm(const state& G_term);

  /** @brief Get the projected in-map goal G (thread-safe).
   *  @param G Output projected goal.
   */
  void getG(state& G);

  /** @brief Get the local trajectory endpoint E (thread-safe).
   *  @param E Output endpoint state.
   */
  void getE(state& E);

  /** @brief Set the projected in-map goal G (thread-safe).
   *  @param G New projected goal.
   */
  void setG(const state& G);

  /** @brief Get the planning start state A (thread-safe).
   *  @param A Output start state.
   */
  void getA(state& A);

  /** @brief Set the planning start state A (thread-safe).
   *  @param A New start state.
   */
  void setA(const state& A);

  /** @brief Get the planning start time (thread-safe).
   *  @param A_time Output start time.
   */
  void getA_time(double& A_time);

  /** @brief Set the planning start time (thread-safe).
   *  @param A_time New start time.
   */
  void setA_time(double A_time);

  /** @brief Get the current robot state (thread-safe).
   *  @param state Output current state.
   */
  void getState(state& state);

  /** @brief Get all tracked dynamic obstacle trajectories.
   *  @return Vector of shared pointers to dynamic trajectories.
   */
  std::vector<std::shared_ptr<dynTraj>> getTrajs();

  /** @brief Get the last state from the current plan.
   *  @param state Output last planned state.
   */
  void getLastPlanState(state& state);

  /** @brief Remove expired dynamic obstacle trajectories.
   *  @param current_time Current time; trajectories older than this are removed.
   */
  void cleanUpOldTrajs(double current_time);

  /** @brief Add a new dynamic obstacle trajectory (thread-safe).
   *  @param new_traj Shared pointer to the new trajectory.
   *  @param current_time Current time for age tracking.
   */
  void addTraj(std::shared_ptr<dynTraj> new_traj, double current_time);

  /** @brief Update the current robot state (thread-safe).
   *  @param data New state data.
   */
  void updateState(state data);

  /** @brief Sample the next goal setpoint from the current plan.
   *  @param next_goal Output goal state.
   *  @return True if a valid goal was obtained.
   */
  bool getNextGoal(state& next_goal);

  /** @brief Check if the planner is ready to replan.
   *  @return True if state and terminal goal have been initialized.
   */
  bool checkReadyToReplan();

  /** @brief Set the terminal goal and reset the drone status to YAWING.
   *  @param term_goal New terminal goal state.
   */
  void setTerminalGoal(const state& term_goal);

  /** @brief If the goal lies in an occupied voxel, relocate it to the nearest
   *         free/unknown cell pushed outward by ||drone_bbox|| along the
   *         direction from the original (occupied) goal. The z component is
   *         left to the caller to re-clamp if needed. Goal is mutated in place.
   *  @param goal Terminal goal state to sanitize (modified in place).
   *  @return True if goal is safe to use (was already non-occupied or was
   *          successfully relocated). False if the map is initialized but no
   *          non-occupied cell could be found within the BFS budget — caller
   *          should drop the goal in that case.
   */
  bool sanitizeTerminalGoal(state& goal);

  /** @brief Change the drone state machine status.
   *  @param new_status One of DroneStatus enum values.
   */
  void changeDroneStatus(int new_status);

  /** @brief Compute the desired yaw for the next goal based on velocity direction.
   *  @param next_goal Goal state whose yaw field is updated.
   */
  void getDesiredYaw(state& next_goal);

  /** @brief Command a yaw rotation by the given angular difference.
   *  @param diff Yaw angle difference in radians.
   *  @param next_goal Output goal state with updated yaw.
   */
  void yaw(double diff, state& next_goal);

  /** @brief Compute the in-map goal G by projecting G_term into the known map.
   *  @param A Current start state.
   *  @param G_term Terminal goal.
   *  @param horizon Planning horizon distance.
   */
  void computeG(const state& A, const state& G_term, double horizon);

  /** @brief Corridor-center subgoal hopping / bend pre-alignment (ground robot only).
   *
   *  Detects sharp bends on the RAW A* path @p raw_global_path using a
   *  windowed direction average (robust to resample phase and A* diagonal
   *  grid noise). At the first qualifying bend, computes a subgoal backed
   *  off by corridor_backoff_m along the windowed incoming direction, sets
   *  G.pos to that point and G.yaw to the windowed outgoing direction, and
   *  replaces @p global_path with the raw path truncated + ending at the
   *  backed-off subgoal so the L-BFGS solver plans only up to the subgoal.
   *  Fall-through (no bend found) uses the path end as G.
   */
  void computeG_corridorHop(const state& A, const state& G_term,
                            vec_Vecf<3>& global_path,
                            const vec_Vecf<3>& raw_global_path);

  /** @brief Check if the robot has reached the terminal goal.
   *  @return True if the goal is reached.
   */
  bool goalReachedCheck();

  /** @brief Set the map bounding box dimensions.
   *  @param min_pos Minimum corner of the map.
   *  @param max_pos Maximum corner of the map.
   */
  void computeMapSize(const Eigen::Vector3d& min_pos, const Eigen::Vector3d& max_pos);

  /** @brief Check if a point lies within the current map bounds.
   *  @param point Query position.
   *  @return True if the point is inside the map.
   */
  bool checkPointWithinMap(const Eigen::Vector3d& point) const;

  /** @brief Get the static push points used for path deformation (thread-safe).
   *  @param static_push_points Output push point positions.
   */
  void getStaticPushPoints(vec_Vecf<3>& static_push_points);

  /** @brief Get the local segment of the global path before and after push deformation.
   *  @param local_global_path Output local path before push.
   *  @param local_global_path_after_push Output local path after push.
   */
  void getLocalGlobalPath(vec_Vecf<3>& local_global_path,
                          vec_Vecf<3>& local_global_path_after_push);

  /** @brief Get the current global path (thread-safe).
   *  @param global_path Output global waypoint path.
   */
  void getGlobalPath(vec_Vecf<3>& global_path);

  /** @brief Get the original global path before any modifications (thread-safe).
   *  @param original_global_path Output original path.
   */
  void getOriginalGlobalPath(vec_Vecf<3>& original_global_path);

  /** @brief Get the global path truncated to free space.
   *  @param free_global_path Output free-space path.
   */
  void getFreeGlobalPath(vec_Vecf<3>& free_global_path);

  /** @brief Generate a local trajectory along the global path using L-BFGS optimization.
   *  @param local_A Start state for the local trajectory.
   *  @param A_time Start timestamp.
   *  @param global_path Global waypoints guiding the trajectory.
   *  @param initial_guess_computation_time Output time for initial guess computation.
   *  @param local_traj_computation_time Output time for local trajectory optimization.
   *  @param whole_traj_solver_ptr L-BFGS solver instance.
   *  @return True if a feasible local trajectory was generated.
   */
  bool generateLocalTrajectory(const state& local_A, double A_time, vec_Vec3f& global_path,
                               double& initial_guess_computation_time,
                               double& local_traj_computation_time,
                               std::shared_ptr<lbfgs::SolverLBFGS>& whole_traj_solver_ptr);

  /** @brief Reset all internal timing and cost data to zero. */
  void resetData();

  /** @brief Retrieve computation time breakdowns from the last replan.
   *  @param final_g Final path cost.
   *  @param global_planning_time Total global planning time.
   *  @param hgp_static_jps_time Static JPS time.
   *  @param hgp_check_path_time Path checking time.
   *  @param hgp_dynamic_astar_time Dynamic A* time.
   *  @param hgp_recover_path_time Path recovery time.
   *  @param cvx_decomp_time Convex decomposition time.
   *  @param initial_guess_computation_time Initial guess computation time.
   *  @param local_traj_computatoin_time Local trajectory optimization time.
   *  @param safety_check_time Safety check time.
   *  @param safe_paths_time Safe paths computation time.
   *  @param yaw_sequence_time Yaw sequence computation time.
   *  @param yaw_fitting_time Yaw fitting time.
   */
  void retrieveData(double& final_g, double& global_planning_time, double& hgp_static_jps_time,
                    double& hgp_check_path_time, double& hgp_dynamic_astar_time,
                    double& hgp_recover_path_time, double& cvx_decomp_time,
                    double& initial_guess_computation_time, double& local_traj_computatoin_time,
                    double& safety_check_time, double& safe_paths_time, double& yaw_sequence_time,
                    double& yaw_fitting_time);

  /** @brief Retrieve polytopes from the last safe corridor computation.
   *  @param poly_out_whole Output whole-trajectory polytopes.
   *  @param poly_out_safe Output safe-path polytopes.
   */
  void retrievePolytopes(vec_E<Polyhedron<3>>& poly_out_whole, vec_E<Polyhedron<3>>& poly_out_safe);

  /** @brief Retrieve goal setpoints from the last trajectory optimization.
   *  @param goal_setpoints Output vector of goal states.
   */
  void retrieveGoalSetpoints(std::vector<state>& goal_setpoints);

  /** @brief Retrieve sub-optimal goal setpoint lists from multi-initial-guess optimization.
   *  @param list_subopt_goal_setpoints Output list of sub-optimal goal state vectors.
   */
  void retrieveListSubOptGoalSetpoints(std::vector<std::vector<state>>& list_subopt_goal_setpoints);

  /** @brief Retrieve control points from the last trajectory optimization.
   *  @param cps Output vector of 3x6 control point matrices per segment.
   */
  void retrieveCPs(std::vector<Eigen::Matrix<double, 3, 6>>& cps);

  /** @brief Generate a global path from the current state to the goal.
   *  @param global_path Output global waypoint path.
   *  @param current_time Current wall/sim time.
   *  @param last_replaning_computation_time Previous replan computation time.
   *  @return True if a valid global path was found.
   */
  bool generateGlobalPath(vec_Vecf<3>& global_path, double current_time,
                          double last_replaning_computation_time);

  /** @brief Push the global path waypoints into free space.
   *  @param global_path Input/output global path.
   *  @param free_global_path Output path truncated to free space.
   *  @param current_time Current wall/sim time.
   *  @return True if the pushed path is valid.
   */
  bool pushPath(vec_Vecf<3>& global_path, vec_Vecf<3>& free_global_path, double current_time);

  /** @brief Plan a local trajectory along the given global path.
   *  @param global_path Global waypoints guiding trajectory optimization.
   *  @return True if local trajectory planning succeeded.
   */
  bool planLocalTrajectory(vec_Vecf<3>& global_path);

  /** @brief Append the latest local trajectory to the current plan.
   *  @return True if the trajectory was successfully appended.
   */
  bool appendToPlan();

  /** @brief Set the initial pose transform for hardware frame alignment.
   *  @param init_pose Transform from planner frame to hardware frame.
   */
  void setInitialPose(const geometry_msgs::msg::TransformStamped& init_pose);

  /** @brief Set the 2D ESDF grid for ground robot obstacle cost. Thread-safe (immutable snapshot). */
  void setEsdfGrid(std::shared_ptr<const class EsdfGrid2D> grid) { esdf_grid_ = grid; }

  /** @brief Set the binary 2D occupancy grid for ground robot A* planning. */
  void setOccGrid2D(std::shared_ptr<const class OccGrid2D> grid) { occ_grid_2d_ = grid; }

  /** @brief Apply the initial pose transform to a piecewise polynomial trajectory.
   *  @param pwp Trajectory to transform in place.
   */
  void applyInitiPoseTransform(PieceWisePol& pwp);

  /** @brief Apply the inverse initial pose transform to a piecewise polynomial trajectory.
   *  @param pwp Trajectory to inverse-transform in place.
   */
  void applyInitiPoseInverseTransform(PieceWisePol& pwp);

  /** @brief Update both the occupancy and unknown maps from point clouds.
   *  @param pclptr_map Occupied point cloud.
   *  @param pclptr_unk Unknown-space point cloud.
   */
  void updateMap(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_map,
                 const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_unk);

  /** @brief Extract dynamic obstacle positions and predicted samples for heat map.
   *  @param trajs Dynamic trajectories.
   *  @param current_time Current timestamp.
   *  @param agent_pos Current agent position.
   *  @param obst_pos Output: current obstacle positions.
   *  @param obst_bbox Output: obstacle bounding box half-extents.
   *  @param traj_max_time Output: prediction horizon.
   */
  void extractDynamicHeatData(const std::vector<std::shared_ptr<dynTraj>>& trajs,
                              double current_time, const Eigen::Vector3d& agent_pos,
                              vec_Vecf<3>& obst_pos, vec_Vecf<3>& obst_bbox,
                              double& traj_max_time);

  /** @brief Update only the occupancy map from a point cloud.
   *  @param pclptr_map Occupied point cloud.
   */
  void updateOccupancyMap(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_map);

  /** @brief Update only the unknown-space cloud.
   *  @param pclptr_unk Unknown-space point cloud.
   */
  void updateUnknownCloud(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_unk);

  /** @brief Get the piecewise quintic polynomial representation of the current trajectory.
   *  @param pwp Output piecewise quintic polynomial.
   */
  void getPiecewiseQuinticPol(PieceWiseQuinticPol& pwp);

  /** @brief Get a shared pointer to the underlying voxel map utility.
   *  @return Shared pointer to VoxelMapUtil.
   */
  std::shared_ptr<mighty::VoxelMapUtil> getMapUtil() const { return hgp_manager_.map_util_; }

 private:
  // Parameters
  parameters par_;          // Parameters of the planner
  HGPManager hgp_manager_;  // HGP Manager
  std::vector<LinearConstraint3D>
      safe_corridor_polytopes_whole_;  // Polytope (Linear) constraints for whole trajectory
  std::shared_ptr<lbfgs::SolverLBFGS>
      whole_traj_solver_ptr_;                    // L-BFGS solver pointer for the whole trajectory
  std::shared_ptr<const class EsdfGrid2D> esdf_grid_;  // 2D ESDF grid (ground robot only)
  std::shared_ptr<const class OccGrid2D> occ_grid_2d_;  // Binary 2D occupancy (ground robot only)
  std::vector<std::shared_ptr<dynTraj>> trajs_;  // Dynamic trajectory
  Eigen::Vector3d v_max_3d_;                     // Maximum velocity
  Eigen::Vector3d a_max_3d_;                     // Maximum acceleration
  Eigen::Vector3d j_max_3d_;                     // Maximum jerk
  double v_max_;                                 // Maximum speed
  double max_dist_vertexes_;                     // Maximum velocity
  lbfgs::planner_params_t planner_params_;       // Planner parameters
  lbfgs::lbfgs_parameter_t lbfgs_params_;        // L-BFGS parameters

  // Flags
  bool state_initialized_ = false;          // State initialized
  bool terminal_goal_initialized_ = false;  // Terminal goal initialized
  bool use_adapt_k_value_ = false;          // Use adapt k value
  bool kdtree_map_initialized_ = false;     // Kd-tree for the map initialized
  bool kdtree_unk_initialized_ = false;     // Kd-tree for the map initialized

  // Data
  double final_g_ = 0.0;
  double global_planning_time_ = 0.0;
  double hgp_static_jps_time_ = 0.0;
  double hgp_check_path_time_ = 0.0;
  double hgp_dynamic_astar_time_ = 0.0;
  double hgp_recover_path_time_ = 0.0;
  double cvx_decomp_time_ = 0.0;
  double initial_guess_computation_time_ = 0.0;
  double local_traj_computation_time_ = 0.0;
  double safe_paths_time_ = 0.0;
  double safety_check_time_ = 0.0;
  double yaw_sequence_time_ = 0.0;
  double yaw_fitting_time_ = 0.0;
  vec_E<Polyhedron<3>> poly_out_whole_;
  vec_E<Polyhedron<3>> poly_out_safe_;
  std::vector<state> goal_setpoints_;
  std::vector<double> optimal_yaw_sequence_;
  std::vector<double> yaw_control_points_;
  std::vector<double> yaw_knots_;
  std::vector<Eigen::Matrix<double, 3, 6>> cps_;
  std::vector<Eigen::VectorXd> list_z_subopt_;
  std::vector<std::vector<Eigen::Vector3d>> list_initial_guess_wps_subopt_;
  std::vector<std::vector<state>> list_subopt_goal_setpoints_;

  // Basis
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_;
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_inverse_;

  // Replanning-related variables
  state state_;                                     // State for the drone
  state G_;                                         // This goal is always inside of the map
  state A_;                                         // Starting point of the drone
  double A_time_;                                   // Time of the starting point
  state E_;                                         // The goal point of actual trajectory
  state G_term_;                                    // Terminal goal
  std::deque<state> plan_;                          // Plan for the drone
  std::deque<std::vector<state>> plan_safe_paths_;  // Indicate if the state has a safe path
  double previous_yaw_ = 0.0;                       // Previous yaw
  double prev_dyaw_ = 0.0;                          // Previous dyaw
  double dyaw_filtered_ = 0.0;                      // Filtered dyaw
  PieceWiseQuinticPol pwp_to_share_;                // Piecewise polynomial to share

  // Drone status
  int drone_status_ =
      DroneStatus::GOAL_REACHED;  // status_ can be TRAVELING, GOAL_SEEN, GOAL_REACHED

  // Corridor-hop turn-in-place override (ground robot only). When true, the
  // YAWING branch in the yaw controller uses G_.yaw (the heading to the next
  // subgoal) as the target instead of the default "direction from current
  // pos toward G_term". Cleared on YAWING -> TRAVELING transition.
  bool corridor_hop_yawing_ = false;

  // Hysteresis for bend pre-alignment subgoal selection. Once a subgoal is
  // chosen for a leg, we hold it across replans (even if the HGP path wiggles
  // with sensor noise) until the robot arrives (goal_radius) or the held
  // subgoal is no longer near any waypoint on the new path. Cleared on
  // YAWING -> TRAVELING transition so the next leg re-picks.
  bool held_subgoal_valid_ = false;
  Eigen::Vector3d held_subgoal_pos_ = Eigen::Vector3d::Zero();
  double held_subgoal_yaw_ = 0.0;

  // Mutex
  std::mutex mtx_plan_;                  // Mutex for the plan_
  std::mutex mtx_state_;                 // Mutex for the state_
  std::mutex mtx_G_;                     // Mutex for the G_
  std::mutex mtx_A_;                     // Mutex for the A_
  std::mutex mtx_A_time_;                // Mutex for the A_time_
  std::mutex mtx_G_term_;                // Mutex for the G_term_
  std::mutex mtx_E_;                     // Mutex for the E_
  std::mutex mtx_trajs_;                 // Mutex for the trajs_
  std::mutex mtx_solve_hgp_;             // Mutex for the solveHGP
  std::mutex mtx_global_path_;           // Mutex for the global_path_
  std::mutex mtx_original_global_path_;  // Mutex for the original_global_path_
  std::mutex mtx_kdtree_map_;            // Mutex for the map_
  std::mutex mtx_kdtree_unk_;            // Mutex for the unknown map_
  pcl::PointCloud<pcl::PointXYZ>::ConstPtr pclptr_map_;
  pcl::PointCloud<pcl::PointXYZ>::ConstPtr pclptr_unk_;

  // Map resolution
  double map_res_;

  // Counter for replanning failure
  int replanning_failure_count_ = 0;

  // Map size
  bool map_size_initialized_ = false;  // Map size initialized
  int hgp_failure_count_ = 0;          // HGP failure count used for map size adaptation

  // Communication delay
  std::unordered_map<int, double> comm_delay_map_;

  // Dynamic k value
  int num_replanning_ = 0;                       // Number of replanning
  bool got_enough_replanning_ = false;           // Check if tot enough replanning
  int k_value_ = 0;                              // k value
  std::vector<double> store_computation_times_;  // Store computation times
  double est_comp_time_ = 0.0;                   // Computation time estimation

  // Map size
  double wdx_, wdy_, wdz_;         // Width of the map
  Vec3f map_center_;               // Center of the map
  double current_dynamic_buffer_;  // Current dynamic buffer

  // Global path
  vec_Vecf<3> global_path_;
  vec_Vecf<3> original_global_path_;  // For visualization
  vec_Vecf<3> free_global_path_;

  // Static push
  vec_Vecf<3> static_push_points_;
  vec_Vecf<3> p_points_;
  vec_Vecf<3> local_global_path_;
  vec_Vecf<3> local_global_path_after_push_;

  // Initial pose
  geometry_msgs::msg::TransformStamped init_pose_;
  Eigen::Matrix4d init_pose_transform_;
  Eigen::Matrix3d init_pose_transform_rotation_;
  Eigen::Matrix4d init_pose_transform_inv_;
  Eigen::Matrix3d init_pose_transform_rotation_inv_;
  double yaw_init_offset_ = 0.0;

  // Safe corridor
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3>,
              Eigen::aligned_allocator<Eigen::Matrix<double, Eigen::Dynamic, 3>>>
      A_stat_;
  std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> b_stat_;

  // kd-tree for the map
  pcl::KdTreeFLANN<pcl::PointXYZ> kdtree_map_;  // kdtree of the point cloud of the occupancy grid
  pcl::KdTreeFLANN<pcl::PointXYZ> kdtree_unk_;  // kdtree of the point cloud of the unknown grid

  // Store data
  Eigen::VectorXd zopt_;
  double fopt_;
};
