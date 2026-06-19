/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <fstream>
#include <optional>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <message_filters/time_synchronizer.h>
#include <pcl_conversions/pcl_conversions.h>
#include <tf2_eigen/tf2_eigen.h>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>

#include <dynus_interfaces/msg/dyn_traj.hpp>
#include <dynus_interfaces/msg/dyn_traj_array.hpp>
#include <dynus_interfaces/msg/frontier_list.hpp>
#include <geometry_msgs/msg/point_stamped.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/path.hpp>
#include <dynus_interfaces/msg/speedy_path.hpp>
#include <std_msgs/msg/empty.hpp>
#include <std_msgs/msg/float64.hpp>
#include <std_msgs/msg/string.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include <mighty/utils.hpp>

#include "hgp/utils.hpp"
#include "mighty/mighty.hpp"
#include "mighty/mighty_type.hpp"

// Frontier exploration (ground robot only). These are forward-declared in the
// header to keep build dependencies minimal — full headers are included in
// mighty_node.cpp.
class FrontierDetector;
class FrontierManager;
class VisitedMap;
struct FrontierRecord;

#include "mighty/peer_tracker.hpp"

#include "ament_index_cpp/get_package_share_directory.hpp"
#include "rclcpp/rclcpp.hpp"

#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/pn_adaptation.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include "dynus_interfaces/msg/trajectory.hpp"
#include "dynus_interfaces/msg/yaw_output.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored \
    "-Wdeprecated-declarations"  // pcl::SAC_SAMPLE_SIZE is protected since PCL 1.8.0 // NOLINT
#include <pcl/sample_consensus/model_types.h>
#pragma GCC diagnostic pop
#include <algorithm>
#include <execution>
#include <memory>
#include <string>

#include <pcl/filters/extract_indices.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/segmentation/sac_segmentation.h>

#include <sensor_msgs/point_cloud2_iterator.hpp>

#include "nav_msgs/msg/occupancy_grid.hpp"
#include "pcl_ros/transforms.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "std_msgs/msg/color_rgba.hpp"
#include "std_srvs/srv/empty.hpp"
#include "tf2_ros/create_timer_ros.h"
#include "tf2_ros/message_filter.h"

typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::msg::PointCloud2,
                                                        sensor_msgs::msg::PointCloud2>
    MySyncPolicy;
typedef message_filters::Synchronizer<MySyncPolicy> Sync;

namespace mighty {

using PCLPoint = pcl::PointXYZ;
using PCLPointCloud = pcl::PointCloud<pcl::PointXYZ>;

/** @brief Main ROS 2 node for the MIGHTY trajectory planner.
 *
 *  Manages subscriptions (point clouds, odometry, Vicon, inter-agent trajectories),
 *  replanning timer callbacks, trajectory publishing, and visualization. Orchestrates
 *  the MIGHTY planner, HGP global planner, and L-BFGS local optimizer.
 */
class MIGHTY_NODE : public rclcpp::Node {
 public:
  /** @brief Construct the node, declare parameters, set up publishers/subscribers and timers. */
  MIGHTY_NODE();

  /** @brief Destructor. Logs benchmarking data if enabled. */
  ~MIGHTY_NODE();

 private:
  // Callbacks
  void replanCallback();
  void trajCallback(const dynus_interfaces::msg::DynTraj::SharedPtr msg);
  void stateCallback(const dynus_interfaces::msg::State::SharedPtr msg);
  void terminalGoalCallback(const geometry_msgs::msg::PoseStamped& msg);
  // Internal entry point used by both the public terminalGoalCallback (sets
  // manual_goal_active_) and the exploration loop (does not).
  void terminalGoalCallbackImpl(const geometry_msgs::msg::PoseStamped& msg, bool from_user);
  // Shared swarm goal subscriber: applies formation_self_offset and forwards
  // through terminalGoalCallbackImpl. Only created when use_formation: true.
  void swarmGoalCallback(const geometry_msgs::msg::PoseStamped& msg);
  // Frontier exploration goal-selection loop, runs at expl_select_rate_hz.
  void exploreSelectCallback();
  void mapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& pcl2ptr_map_ros,
                   const sensor_msgs::msg::PointCloud2::ConstPtr& pcl2ptr_unk_ros);
  void occupancyMapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& map_msg);
  void unknownMapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& unk_msg);
  void esdfCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg);
  void occ2DCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg);
  void goalReachedCheckCallback();
  void convertDynTrajMsg2DynTraj(const dynus_interfaces::msg::DynTraj& msg,
                                 std::shared_ptr<dynTraj>& traj, double current_time);
  void cleanUpOldTrajsCallback();
  void getInitialPoseHwCallback();
  void frameAlignCallback(const geometry_msgs::msg::TransformStamped::SharedPtr msg, int agent_id);
  void applyFrameAlignTransform(std::shared_ptr<dynTraj>& traj, int sender_id);
  void applyTransformToTraj(std::shared_ptr<dynTraj>& traj, const Eigen::Matrix4d& T);

  // Others
  void declareParameters();
  void setParameters();
  void printParameters();
  void createMarkerArrayFromVec_Vec3f(const vec_Vec3f& occupied_cells,
                                      const std_msgs::msg::ColorRGBA& color, int namespace_id,
                                      double scale,
                                      visualization_msgs::msg::MarkerArray* marker_array);
  void clearMarkerArray(
      visualization_msgs::msg::MarkerArray& path_marker,
      rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr publisher);
  void runSim();
  void printComputationTime(bool result);
  void recordData(bool result);
  void logData();
  void setComputationTimesToZero();
  void constructFOVMarker();
  void retrieveData();
  void retrieveDataForVisualizationAndTopics();

  // Functions to publish
  void publishGlobalPath();
  void publishFreeGlobalPath();
  void publishPoly();
  void publishTraj();
  void publishOwnTraj();
  void publishActualTraj();
  void publishGoal();          // Publish the goal (trajectory points)
  void publishTrajectory();    // Publish the full trajectory with replan support
  void publishMpcPath();       // Publish smoothed global path as nav_msgs/Path for MPC
  void publishPointG() const;  // Publish the point G (projected terminal goal)
  void publishPointE() const;  // Publish the point E (Local trajectory's goal)
  void publishPointA() const;  // Publish the point A (starting point)
  void publishCurrentState(const state& state) const;
  void publishState(
      const state& data,
      const rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr& publisher) const;
  void publishFOV();
  void publishCps();
  void publishHeatCloud();
  void publishGround2DOccupied();
  void publishGround2DHeat();
  void publishStaticPushPoints();
  void publishLocalGlobalPath();
  void publishVelocityInText(const Eigen::Vector3d& position, double velocity);
  // Frontier exploration visualization
  void publishFrontierMarkers();
  void publishFrontierData();
  void publishExplorationCurrentGoal(const FrontierRecord& r);
  void publishVisitedMap();

  // Timers for callback
  rclcpp::TimerBase::SharedPtr timer_replanning_;
  rclcpp::TimerBase::SharedPtr timer_goal_;
  rclcpp::TimerBase::SharedPtr timer_update_tmap_;
  rclcpp::TimerBase::SharedPtr timer_goal_reached_check_;
  rclcpp::TimerBase::SharedPtr timer_cleanup_old_trajs_;
  rclcpp::TimerBase::SharedPtr timer_trajs_update_for_tmap_;
  rclcpp::TimerBase::SharedPtr timer_check_plan_safety_;
  rclcpp::TimerBase::SharedPtr timer_initial_pose_;

  // Callback groups
  std::vector<rclcpp::CallbackGroup::SharedPtr> cb_groups_mu_;  // mutually-exclusive groups
  std::vector<rclcpp::CallbackGroup::SharedPtr> cb_groups_re_;  // reentrant groups
  rclcpp::CallbackGroup::SharedPtr cb_group_map_;
  rclcpp::CallbackGroup::SharedPtr cb_group_replan_;
  rclcpp::CallbackGroup::SharedPtr cb_group_goal_;

  // Publishers
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_hgp_path_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_original_hgp_path_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_free_hgp_path_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_local_global_path_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr
      pub_local_global_path_after_push_marker_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_dynamic_map_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_free_map_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_unknown_map_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_static_map_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_dynamic_map_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_free_map_marker_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_unknown_map_marker_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_heat_cloud_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_ground_2d_occ_;
  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_ground_2d_heat_;
  rclcpp::Publisher<decomp_ros_msgs::msg::PolyhedronArray>::SharedPtr pub_poly_whole_;
  rclcpp::Publisher<decomp_ros_msgs::msg::PolyhedronArray>::SharedPtr pub_poly_safe_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_traj_committed_colored_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_traj_subopt_colored_;
  rclcpp::Publisher<dynus_interfaces::msg::DynTraj>::SharedPtr pub_own_traj_;
  rclcpp::Publisher<dynus_interfaces::msg::Goal>::SharedPtr pub_goal_;
  rclcpp::Publisher<dynus_interfaces::msg::Trajectory>::SharedPtr pub_trajectory_;
  rclcpp::Publisher<dynus_interfaces::msg::SpeedyPath>::SharedPtr pub_mpc_path_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_G_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_E_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_G_term_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_A_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_current_state_;
  rclcpp::Publisher<std_msgs::msg::Empty>::SharedPtr pub_goal_reached_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_setpoint_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_actual_traj_;
  rclcpp::Publisher<dynus_interfaces::msg::YawOutput>::SharedPtr pub_yaw_output_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_fov_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_cp_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_static_push_points_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_p_points_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_vel_text_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_traj_received_;
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_traj_transformed_;
  // Corridor-hop yaw-target arrow at G (ground robot only).
  rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_corridor_yaw_target_;

  // Subscribers
  rclcpp::Subscription<dynus_interfaces::msg::DynTraj>::SharedPtr sub_traj_;
  rclcpp::Subscription<dynus_interfaces::msg::DynTraj>::SharedPtr sub_predicted_traj_;
  rclcpp::Subscription<dynus_interfaces::msg::State>::SharedPtr sub_state_;
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr sub_terminal_goal_;
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr sub_swarm_goal_;
  rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr sub_fake_sim_occupancy_map_;

  // Time synchronizer
  message_filters::Subscriber<sensor_msgs::msg::PointCloud2> occup_grid_sub_;
  message_filters::Subscriber<sensor_msgs::msg::PointCloud2> unknown_grid_sub_;
  std::shared_ptr<Sync> sync_;

  // Independent map subscribers (fallback when sync fails, e.g. hardware)
  rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr sub_occupancy_grid_;
  rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr sub_unknown_grid_;

  // ESDF subscription (ground robot only)
  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr sub_esdf_2d_;
  std::shared_ptr<const class EsdfGrid2D> esdf_grid_;

  // Binary 2D occupancy subscription (ground robot only)
  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr sub_occ_2d_;
  std::shared_ptr<const class OccGrid2D> occ_grid_2d_;

  // Frontier exploration (ground robot only). Detector + persistent global
  // frontier database. See plan: /home/kkondo/.claude/plans/snazzy-moseying-donut.md
  std::unique_ptr<FrontierDetector> frontier_detector_;
  std::unique_ptr<FrontierManager>  frontier_manager_;
  // Persistent "ever-observed" bitmap. Suppresses re-detection of frontiers
  // along the seam where the sliding mapper window re-blanks revisited
  // areas to UNKNOWN. Owned and updated here, queried by frontier_detector_.
  std::unique_ptr<VisitedMap>       visited_map_;
  // Most recent grid the detector ran on. Aliases occ_grid_2d_ when
  // expl_detect_on_visited_map is false; otherwise it's a snapshot of the
  // fused (own + peers) visited_map. The manager and goal-selection paths
  // must read state from the same grid the detector saw, otherwise frontiers
  // cleared by peers stay stuck ACTIVE because the local sliding window
  // still shows UNKNOWN there.
  std::shared_ptr<const OccGrid2D>  current_detect_grid_;
  bool     exploration_active_     = false;  // we issued the current goal
  bool     manual_goal_active_     = false;  // user issued the current goal
  uint64_t current_explore_id_     = 0;
  int      unreachable_consec_count_ = 0;
  Eigen::Vector3d exploration_start_pos_{0.0, 0.0, 0.0};
  bool exploration_start_captured_ = false;  // sticky for the whole exploration session
  rclcpp::TimerBase::SharedPtr timer_explore_select_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_frontiers_;
  rclcpp::Publisher<dynus_interfaces::msg::FrontierList>::SharedPtr pub_frontier_data_;
  rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pub_explore_current_goal_;
  rclcpp::Publisher<nav_msgs::msg::OccupancyGrid>::SharedPtr pub_visited_map_;
  // MinPos peer tracking (multi-robot frontier allocation)
  PeerTracker peer_tracker_;
  rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pub_peer_pose_;
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr sub_peer_pose_;
  double last_peer_pose_publish_t_ = 0.0;
  // Visited map sharing (global topic, all agents pub+sub)
  rclcpp::Publisher<nav_msgs::msg::OccupancyGrid>::SharedPtr pub_peer_visited_map_;
  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr sub_peer_visited_map_;
  double last_peer_visited_publish_t_ = 0.0;
  // Global return-home trigger. A single Empty publish on /exploration/return_home
  // makes every agent issue a goal back to its captured exploration_start_pos_
  // and stop accepting new frontier goals (so it stays parked once arrived).
  rclcpp::Subscription<std_msgs::msg::Empty>::SharedPtr sub_return_home_;
  bool home_return_requested_ = false;
  // Wall-clock seconds of the last successful publishVisitedMap() call.
  // Used to throttle the (potentially large) tristate-grid publish to ~1 Hz
  // — RViz only needs occasional updates because the persistent map only
  // grows incrementally and `transient_local` QoS replays the latest snapshot
  // to late subscribers.
  double last_visited_publish_t_ = 0.0;

  // Latched origin.z of the most recent occ_2d_topic message, so the
  // visited_map we republish renders at the same ground plane as the live
  // occupancy grid (global_mapper sets it to z_ground). std::nullopt until
  // we've seen a real occ_2d — falls back to expl_default_goal_z then.
  std::optional<double> occ2d_origin_z_;

  // Wall-clock seconds of the last visualization publish in replanCallback.
  // The replan loop runs at 100 Hz which is fine for control but floods RViz
  // (especially for traj_committed_colored, hgp_path_marker, and the other
  // multi-KB MarkerArrays). RViz drops messages when the publish rate exceeds
  // its processing capacity, showing "some messages were lost" warnings.
  // We throttle the entire viz block to ~20 Hz (every 50 ms) — plenty smooth
  // visually, ~5x cheaper to render, no message loss.
  double last_replan_viz_publish_t_ = 0.0;

  // Same throttle for actual_traj. publishActualTraj() is called from
  // stateCallback which fires at the state publisher's rate (100 Hz from
  // fake_sim), and the persistent LINE_STRIP marker grows large enough that
  // RViz drops messages at that rate too.
  double last_actual_traj_publish_t_ = 0.0;

  // Visualization
  visualization_msgs::msg::MarkerArray hgp_path_marker_;
  visualization_msgs::msg::MarkerArray original_hgp_path_marker_;
  visualization_msgs::msg::MarkerArray hgp_free_path_marker_;
  visualization_msgs::msg::MarkerArray hgp_local_global_path_marker_;
  visualization_msgs::msg::MarkerArray hgp_local_global_path_after_push_marker_;
  visualization_msgs::msg::MarkerArray traj_committed_colored_;
  visualization_msgs::msg::MarkerArray traj_subopt_colored_;
  visualization_msgs::msg::Marker marker_fov_;

  // Mutex
  std::mutex cloud_callback_mutex_;  // Mutex for cloud callback

  // Parameters
  int id_;
  std::string ns_;
  std::string id_str_;
  parameters par_;
  uint32_t trajectory_id_ = 0;             // Trajectory ID for replan detection
  double final_g_ = 0.0;                   // only for debugging
  bool verbose_computation_time_ = false;  // only for debugging
  int marker_fov_id_ = 0;

  // DYNUS pointer
  std::shared_ptr<MIGHTY> mighty_ptr_;

  // Global Path Benchmarking
  std::string file_path_;         // only for benchmarking
  std::vector<std::tuple<bool,    // Result
                         double,  // Cost
                         double,  // Total replanning time
                         double,  // Global planning time
                         double,  // CVX decomposition time
                         double,  // Initial guess computation time
                         double,  // Local trajectory time
                         double,  // Safe paths time
                         double,  // Safety Check time
                         double,  // Yaw sequence time
                         double,  // Yaw fitting time
                         double,  // Static JPS time in HGP
                         double,  // Check path time in HGP
                         double,  // Dynamic A* time in HGP
                         double   // Recover path time in HGP
                         >>
      global_path_benchmark_;  // only for benchmarking

  // debug
  double replan_last_time_called_ = 0.0;

  // Record computation time
  double global_planning_time_ = 0.0;
  double hgp_static_jps_time_ = 0.0;
  double hgp_check_path_time_ = 0.0;
  double hgp_dynamic_astar_time_ = 0.0;
  double hgp_recover_path_time_ = 0.0;
  double cvx_decomp_time_ = 0.0;
  double initial_guess_computation_time_ = 0.0;  // Time for computing initial guess
  double local_traj_computation_time_ = 0.0;
  double safe_paths_time_ = 0.0;
  double safety_check_time_ = 0.0;
  double yaw_sequence_time_ = 0.0;
  double yaw_fitting_time_ = 0.0;
  double replanning_computation_time_ = 0.0;
  double current_time_for_debug_ = 0.0;

  // Visualization
  vec_E<Polyhedron<3>> poly_whole_;
  vec_E<Polyhedron<3>> poly_safe_;
  std::vector<state> goal_setpoints_;
  std::vector<std::vector<state>> list_subopt_goal_setpoints_;

  // Yaw Optimization Debugging
  std::vector<double> optimal_yaw_sequence_;
  std::vector<double> yaw_control_points_;
  std::vector<double> yaw_knots_;

  // Local trajectory debugging
  std::vector<Eigen::Matrix<double, 3, 6>> cps_;

  // Static push points and p points
  vec_Vecf<3> static_push_points_;
  vec_Vecf<3> p_points_;
  int static_push_points_id_ = 0;
  int p_points_id_ = 0;

  // Trajectory sharing
  PieceWiseQuinticPol pwp_to_share_;  // Piecewise polynomial

  // Flags
  bool state_initialized_ = false;           // State initialized
  bool replan_timer_started_ = false;        // Replan timer started
  bool use_benchmark_ = false;               // Use benchmark

  // Actual-trajectory history for visualization (sando-style: a single
  // persistent LINE_STRIP marker built from a bounded history of past states,
  // colored by speed). Replaces the legacy "publish one ARROW marker per
  // sample" approach which leaked thousands of markers into RViz and lost
  // them frequently because the topic was a single Marker (not MarkerArray).
  std::vector<state> actual_traj_hist_;
  bool actual_traj_initialized_ = false;
  Eigen::Vector3d actual_traj_prev_pos_{0.0, 0.0, 0.0};
  double actual_traj_prev_time_ = 0.0;
  size_t actual_traj_max_hist_       = 4000;  // max stored states
  size_t actual_traj_max_points_vis_ = 300;   // max points shown in RViz after downsampling
  double actual_traj_line_width_     = 0.15;  // meters

  int last_subopt_count_{0};  // tracks how many subopt strips we drew last time

  // D435 parameters
  std::string d435_depth_frame_id_;
  std::string lidar_frame_id_;
  std::string d435_camera_info_topic_;

  // TF2 buffer and listener
  std::shared_ptr<tf2_ros::Buffer> tf2_buffer_;
  std::shared_ptr<tf2_ros::TransformListener> tf2_listener_;

  // initial pose (for hardware)
  bool initial_pose_received_ = false;
  std::string initial_pose_topic_;
  geometry_msgs::msg::TransformStamped init_pose_transform_stamped_;

  // Frame alignment transforms (inter-agent)
  std::vector<rclcpp::Subscription<geometry_msgs::msg::TransformStamped>::SharedPtr>
      frame_align_subs_;
  std::unordered_map<int, Eigen::Matrix4d> frame_align_transforms_;
  std::unordered_map<int, bool> frame_align_received_;
  std::mutex frame_align_mutex_;

  // Simulated frame offset matrix (for testing frame alignment in fake_sim)
  Eigen::Matrix4d sim_frame_offset_{Eigen::Matrix4d::Identity()};
  Eigen::Matrix4d sim_frame_offset_inv_{Eigen::Matrix4d::Identity()};

  // Timer to make sure we don't sample point cloud too often
  rclcpp::Time last_lidar_callback_time_;
  rclcpp::Time last_depth_camera_callback_time_;

  // Command-to-execution timing (time from goal received to first trajectory)
  rclcpp::Time goal_received_time_;
  bool waiting_for_first_traj_ = false;
  rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr pub_command_to_exec_time_;
};

}  // namespace mighty