/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <mighty/mighty_node.hpp>
#include <mighty/esdf_grid_2d.hpp>
#include <mighty/occ_grid_2d.hpp>
#include <mighty/frontier_detector.hpp>
#include <mighty/frontier_manager.hpp>
#include <mighty/visited_map.hpp>

#include <fstream>
#include <iomanip>

using namespace std::chrono_literals;

namespace mighty {

// Debug log: writes every published /goal and every received /term_goal to a
// file so the dip-toward-origin bug can be analyzed offline. Set debug_log_file
// to "" in the YAML to disable.
namespace {
double& debug_log_t0() {
  static double t0 = 0.0;
  return t0;
}
std::ofstream& debug_log_stream() {
  static std::ofstream s;  // never closed; flushed on every write
  return s;
}
void debug_log_open(const std::string& path) {
  if (path.empty()) return;
  auto& s = debug_log_stream();
  if (!s.is_open()) {
    s.open(path, std::ios::out | std::ios::trunc);
    if (s.is_open()) {
      s << std::fixed << std::setprecision(4);
      s << "# t_rel_sec | event | px | py | pz | extra\n";
      s.flush();
    }
  }
}
double debug_log_t(double now_sec) {
  auto& t0 = debug_log_t0();
  if (t0 == 0.0) t0 = now_sec;
  return now_sec - t0;
}
}  // namespace

// ----------------------------------------------------------------------------

/**
 * @brief Constructor
 */
MIGHTY_NODE::MIGHTY_NODE() : Node("mighty_node") {
  // Get id from ns
  ns_ = this->get_namespace();
  ns_ = ns_.substr(ns_.find_last_of("/") + 1);
  id_str_ = ns_.substr(ns_.size() -
                       2);  // ns is like NX01, so we get the last two characters and convert to int
  id_ = std::stoi(id_str_);

  // Declare, set, and print parameters
  this->declareParameters();
  this->setParameters();
  this->printParameters();

  // Qos policy settings
  rclcpp::QoS critical_qos(rclcpp::KeepLast(10));
  critical_qos.reliable().durability_volatile();
  // Visualization: best-effort + depth 1. Markers/clouds for RViz only — losing
  // a frame is fine, and reliable+depth-10 was filling Zenoh's TX queue and
  // closing transports under multi-agent load.
  rclcpp::QoS viz_qos(rclcpp::KeepLast(1));
  viz_qos.best_effort().durability_volatile();

  // Initialize callback groups
  cb_groups_mu_.resize(9);
  cb_groups_re_.resize(9);
  for (int i = 0; i < 9; i++) {
    cb_groups_mu_[i] = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
    cb_groups_re_[i] = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  }
  // MutuallyExclusive (NOT Reentrant): every callback on this group mutates
  // shared map state — occ_grid_2d_, esdf_grid_, visited_map_, frontier_manager_
  // — without any internal locking. With Reentrant, two concurrent occ2DCallback
  // invocations would race the `occ_grid_2d_ = OccGrid2D::fromOccupancyGrid(...)`
  // assignment in occ2DCallback: the in-flight FrontierDetector::detect() in one
  // thread holds a raw reference to the previous OccGrid2D, and the second
  // thread's shared_ptr replacement drops its refcount to zero, freeing the
  // unknown_/occupied_ vector data the BFS is still reading. SIGSEGV in
  // isUnknown() inside the BFS expansion. Serializing the group fixes this and
  // also closes the analogous race on frontier_manager_::update()/evict.
  this->cb_group_map_ =
      this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_replan_ =
      this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_goal_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);

  // Options for callback group
  rclcpp::SubscriptionOptions options_re_1;
  options_re_1.callback_group = this->cb_groups_re_[0];
  rclcpp::SubscriptionOptions options_re_2;
  options_re_2.callback_group = this->cb_groups_re_[1];
  rclcpp::SubscriptionOptions options_map;
  options_map.callback_group = this->cb_group_map_;

  // Visulaization publishers
  pub_dynamic_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "dynamic_occupied_grid", viz_qos);
  pub_static_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "static_map_marker", viz_qos);
  pub_dynamic_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "dynamic_map_marker", viz_qos);
  pub_free_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "free_map_marker", viz_qos);
  pub_unknown_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "unknown_map_marker", viz_qos);
  pub_heat_cloud_ =
      this->create_publisher<sensor_msgs::msg::PointCloud2>("heat_cloud", viz_qos);
  pub_ground_2d_occ_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "ground_2d_occupied", viz_qos);
  pub_ground_2d_heat_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "ground_2d_heat", viz_qos);
  pub_free_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "free_grid", viz_qos);
  pub_unknown_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "unknown_grid", viz_qos);
  pub_hgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "hgp_path_marker", viz_qos);
  pub_original_hgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "original_hgp_path_marker", viz_qos);
  pub_free_hgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "free_hgp_path_marker", viz_qos);
  pub_local_global_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "local_global_path_marker", viz_qos);
  pub_local_global_path_after_push_marker_ =
      this->create_publisher<visualization_msgs::msg::MarkerArray>(
          "local_global_path_after_push_marker", viz_qos);
  pub_poly_whole_ = this->create_publisher<decomp_ros_msgs::msg::PolyhedronArray>(
      "poly_whole", viz_qos);
  pub_poly_safe_ = this->create_publisher<decomp_ros_msgs::msg::PolyhedronArray>(
      "poly_safe", viz_qos);
  pub_traj_committed_colored_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "traj_committed_colored", viz_qos);
  pub_traj_subopt_colored_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "traj_subopt_colored", viz_qos);
  pub_setpoint_ = this->create_publisher<geometry_msgs::msg::PointStamped>("setpoint_vis", viz_qos);
  pub_actual_traj_ =
      this->create_publisher<visualization_msgs::msg::MarkerArray>("actual_traj", viz_qos);
  pub_fov_ = this->create_publisher<visualization_msgs::msg::Marker>("fov", viz_qos);
  pub_cp_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("cp", viz_qos);
  pub_static_push_points_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "static_push_points", viz_qos);
  pub_p_points_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
      "p_points", viz_qos);
  pub_point_A_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_A", viz_qos);
  pub_point_G_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_G", viz_qos);
  pub_point_E_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_E", viz_qos);
  pub_point_G_term_ = this->create_publisher<geometry_msgs::msg::PointStamped>(
      "point_G_term", viz_qos);
  pub_current_state_ = this->create_publisher<geometry_msgs::msg::PointStamped>(
      "point_current_state", viz_qos);
  pub_vel_text_ = this->create_publisher<visualization_msgs::msg::Marker>("vel_text", viz_qos);
  pub_traj_received_ = this->create_publisher<visualization_msgs::msg::Marker>(
      "traj_received", viz_qos);
  pub_traj_transformed_ = this->create_publisher<visualization_msgs::msg::Marker>(
      "traj_transformed", viz_qos);
  pub_corridor_yaw_target_ = this->create_publisher<visualization_msgs::msg::Marker>(
      "corridor_yaw_target", viz_qos);

  // Debug publishers
  pub_yaw_output_ = this->create_publisher<dynus_interfaces::msg::YawOutput>("yaw_output", 10);

  // Essential publishers
  pub_own_traj_ = this->create_publisher<dynus_interfaces::msg::DynTraj>("/trajs", critical_qos);
  pub_goal_ = this->create_publisher<dynus_interfaces::msg::Goal>("goal", critical_qos);
  pub_trajectory_ =
      this->create_publisher<dynus_interfaces::msg::Trajectory>("trajectory", critical_qos);
  pub_mpc_path_ = this->create_publisher<dynus_interfaces::msg::SpeedyPath>("mpc_waypoints", 10);
  pub_goal_reached_ = this->create_publisher<std_msgs::msg::Empty>("goal_reached", critical_qos);
  pub_command_to_exec_time_ =
      this->create_publisher<std_msgs::msg::Float64>("command_to_exec_time", 10);

  // Subscribers
  sub_traj_ = this->create_subscription<dynus_interfaces::msg::DynTraj>(
      "/trajs", critical_qos, std::bind(&MIGHTY_NODE::trajCallback, this, std::placeholders::_1),
      options_re_1);
  sub_predicted_traj_ = this->create_subscription<dynus_interfaces::msg::DynTraj>(
      "predicted_trajs", critical_qos,
      std::bind(&MIGHTY_NODE::trajCallback, this, std::placeholders::_1), options_re_1);
  sub_state_ = this->create_subscription<dynus_interfaces::msg::State>(
      "state", critical_qos, std::bind(&MIGHTY_NODE::stateCallback, this, std::placeholders::_1),
      options_re_1);
  sub_terminal_goal_ = this->create_subscription<geometry_msgs::msg::PoseStamped>(
      "term_goal", critical_qos,
      std::bind(&MIGHTY_NODE::terminalGoalCallback, this, std::placeholders::_1));
  // Shared swarm goal: a single /swarm_goal publication is fanned out to each
  // agent's local terminal-goal pin via formation_self_offset, so the whole
  // swarm can be commanded with one PoseStamped message.
  if (par_.use_formation) {
    sub_swarm_goal_ = this->create_subscription<geometry_msgs::msg::PoseStamped>(
        "/swarm_goal", critical_qos,
        std::bind(&MIGHTY_NODE::swarmGoalCallback, this, std::placeholders::_1));
    RCLCPP_INFO(this->get_logger(),
                "Subscribing to /swarm_goal with self offset [%.2f %.2f %.2f]",
                par_.formation_self_offset[0], par_.formation_self_offset[1],
                par_.formation_self_offset[2]);
  }
  // Frame alignment subscriptions (inter-agent transforms)
  if (par_.use_frame_alignment) {
    for (int i = 1; i <= par_.num_agents; i++) {
      if (i == id_) continue;
      std::string prefix = ns_.substr(0, ns_.size() - 2);
      char other_name[16];
      std::snprintf(other_name, sizeof(other_name), "%s%02d", prefix.c_str(), i);
      std::string topic = "/frame_align/" + ns_ + "/" + std::string(other_name);
      frame_align_transforms_[i] = Eigen::Matrix4d::Identity();
      frame_align_received_[i] = false;
      auto sub = this->create_subscription<geometry_msgs::msg::TransformStamped>(
          topic, 10,
          [this, i](const geometry_msgs::msg::TransformStamped::SharedPtr msg) {
            this->frameAlignCallback(msg, i);
          },
          options_re_1);
      frame_align_subs_.push_back(sub);
      RCLCPP_INFO(this->get_logger(), "Frame align: subscribed to %s", topic.c_str());
    }
  }

  // Timer for callback
  timer_replanning_ = this->create_wall_timer(10ms, std::bind(&MIGHTY_NODE::replanCallback, this),
                                              this->cb_group_replan_);
  timer_goal_ =
      this->create_wall_timer(std::chrono::duration<double>(par_.dc),
                              std::bind(&MIGHTY_NODE::publishGoal, this), this->cb_group_goal_);
  // Goal-reached check is needed for both benchmark logging and exploration
  // (so the manager knows when the robot has reached a frontier).
  if (use_benchmark_ || par_.expl_enabled)
    timer_goal_reached_check_ = this->create_wall_timer(
        100ms, std::bind(&MIGHTY_NODE::goalReachedCheckCallback, this), this->cb_groups_re_[2]);
  timer_cleanup_old_trajs_ = this->create_wall_timer(
      500ms, std::bind(&MIGHTY_NODE::cleanUpOldTrajsCallback, this), this->cb_groups_mu_[4]);
  if (par_.use_hardware)
    timer_initial_pose_ = this->create_wall_timer(
        100ms, std::bind(&MIGHTY_NODE::getInitialPoseHwCallback, this), this->cb_groups_mu_[8]);

  // Stop the timer for callback
  if (timer_replanning_) timer_replanning_->cancel();
  if (timer_goal_) timer_goal_->cancel();
  if (!use_benchmark_ && !par_.expl_enabled && timer_goal_reached_check_)
    timer_goal_reached_check_->cancel();

  // Initialize the DYNUS object
  mighty_ptr_ = std::make_shared<MIGHTY>(par_);

  // Initialize the tf2 buffer and listener
  tf2_buffer_ = std::make_shared<tf2_ros::Buffer>(this->get_clock());
  tf2_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf2_buffer_);
  auto timer_interface = std::make_shared<tf2_ros::CreateTimerROS>(
      this->get_node_base_interface(), this->get_node_timers_interface());
  tf2_buffer_->setCreateTimerInterface(timer_interface);

  // Initialize the d435 depth frame ID and camera
  d435_depth_frame_id_ = ns_ + "/" + ns_ + "_d435_depth_optical_frame";
  lidar_frame_id_ = ns_ + "/NX01_livox";

  // Construct the FOV marker
  constructFOVMarker();

  // Initialize the initial pose topic name
  initial_pose_topic_ = ns_ + "/init_pose";

  if (par_.use_hardware) {
    // Hardware: subscribe independently (time-sync can fail due to timestamp mismatch)
    sub_occupancy_grid_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(
        "occupancy_grid",
        rclcpp::QoS(rclcpp::QoSInitialization::from_rmw(rmw_qos_profile_sensor_data)),
        std::bind(&MIGHTY_NODE::occupancyMapCallback, this, std::placeholders::_1), options_map);
    sub_unknown_grid_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(
        "unknown_grid",
        rclcpp::QoS(rclcpp::QoSInitialization::from_rmw(rmw_qos_profile_sensor_data)),
        std::bind(&MIGHTY_NODE::unknownMapCallback, this, std::placeholders::_1), options_map);
    RCLCPP_INFO(this->get_logger(),
                "Hardware mode: subscribing to occupancy_grid and unknown_grid independently");
  } else if (par_.sim_env != "fake_sim") {
    // Gazebo sim: synchronize the occupancy grid and unknown grid
    occup_grid_sub_.subscribe(this, "occupancy_grid", rmw_qos_profile_sensor_data, options_map);
    unknown_grid_sub_.subscribe(this, "unknown_grid", rmw_qos_profile_sensor_data, options_map);
    sync_.reset(new Sync(MySyncPolicy(10), occup_grid_sub_, unknown_grid_sub_));
    sync_->registerCallback(
        std::bind(&MIGHTY_NODE::mapCallback, this, std::placeholders::_1, std::placeholders::_2));
  } else {
    sub_fake_sim_occupancy_map_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(
        "sensor_point_cloud",
        rclcpp::QoS(rclcpp::QoSInitialization::from_rmw(rmw_qos_profile_sensor_data)),
        std::bind(&MIGHTY_NODE::occupancyMapCallback, this, std::placeholders::_1), options_map);
  }

  // ESDF subscription (ground robot only). global_mapper_ros publishes these
  // with BEST_EFFORT reliability — must match here or DDS silently drops
  // delivery and the exploration pipeline never sees an occupancy grid.
  if (par_.use_esdf_cost && par_.vehicle_type == "ground_robot") {
    // SensorDataQoS = BEST_EFFORT + VOLATILE + KEEP_LAST/5 — copies the whole
    // profile, unlike QoSInitialization::from_rmw which only sets history+depth.
    auto map_qos = rclcpp::SensorDataQoS();
    sub_esdf_2d_ = this->create_subscription<nav_msgs::msg::OccupancyGrid>(
        "esdf_2d_topic", map_qos,
        std::bind(&MIGHTY_NODE::esdfCallback, this, std::placeholders::_1), options_map);
    RCLCPP_INFO(this->get_logger(), "ESDF: Subscribed to esdf_2d_topic (d_safe=%.1f m, weight=%.0f)",
                par_.esdf_d_safe, par_.esdf_weight);

    // Also subscribe to binary 2D occupancy for A* planning
    sub_occ_2d_ = this->create_subscription<nav_msgs::msg::OccupancyGrid>(
        "occ_2d_topic", map_qos,
        std::bind(&MIGHTY_NODE::occ2DCallback, this, std::placeholders::_1), options_map);
    RCLCPP_INFO(this->get_logger(), "Occ2D: Subscribed to occ_2d_topic for ground robot A* planning");

    // Frontier-based exploration. The detector + persistent manager run inside
    // occ2DCallback; the explore-select timer issues exploration goals through
    // the same pathway as a manual term_goal.
    if (par_.expl_enabled) {
      FrontierDetectorParams dp;
      dp.cluster_min_cells       = par_.expl_cluster_min_cells;
      dp.border_margin_cells     = par_.expl_border_margin_cells;
      dp.obstacle_clearance_cells = par_.expl_obstacle_clearance_cells;
      dp.robot_snap_radius_m     = par_.expl_robot_snap_radius_m;
      dp.unknown_bridge_cells    = par_.expl_unknown_bridge_cells;
      dp.bounds_enabled          = par_.expl_bounds_enabled;
      dp.bounds_min_x            = par_.expl_bounds_min_x;
      dp.bounds_max_x            = par_.expl_bounds_max_x;
      dp.bounds_min_y            = par_.expl_bounds_min_y;
      dp.bounds_max_y            = par_.expl_bounds_max_y;
      frontier_detector_ = std::make_unique<FrontierDetector>(dp);

      FrontierManagerParams mp;
      mp.merge_radius_m            = par_.expl_merge_radius_m;
      mp.centroid_ema_alpha        = par_.expl_centroid_ema_alpha;
      mp.visit_radius_m            = par_.expl_visit_radius_m;
      mp.visit_dwell_sec           = par_.expl_visit_dwell_sec;
      mp.verify_radius_cells       = par_.expl_verify_radius_cells;
      mp.max_frontiers             = par_.expl_max_frontiers;
      mp.w_size     = par_.expl_w_size;
      mp.w_dist     = par_.expl_w_dist;
      mp.w_info     = par_.expl_w_info;
      mp.w_revisit  = par_.expl_w_revisit;
      mp.w_heading  = par_.expl_w_heading;
      mp.size_ref_m2     = par_.expl_size_ref_m2;
      mp.dist_ref_m      = par_.expl_dist_ref_m;
      mp.sensor_radius_m = par_.expl_sensor_radius_m;
      mp.goal_select_threshold = par_.expl_goal_select_threshold;
      mp.pursuit_timeout_factor  = par_.expl_pursuit_timeout_factor;
      mp.pursuit_timeout_v_ref   = par_.expl_pursuit_timeout_v_ref;
      mp.pursuit_timeout_min_sec = par_.expl_pursuit_timeout_min_sec;
      mp.pursuit_timeout_max_strikes = par_.expl_pursuit_timeout_max_strikes;
      mp.invalidation_keep_out_radius_m = par_.expl_invalidation_keep_out_radius_m;
      mp.invalidation_cooldown_sec      = par_.expl_invalidation_cooldown_sec;
      mp.peer_visit_radius_m            = par_.expl_peer_visit_radius_m;
      frontier_manager_ = std::make_unique<FrontierManager>(mp);

      // Persistent visited bitmap. Records every cell ever observed across
      // the mission so the detector can suppress re-detection along the
      // sliding-window seam when the robot revisits an area.
      visited_map_ = std::make_unique<VisitedMap>(
          par_.expl_visited_map_center_x,
          par_.expl_visited_map_center_y,
          par_.expl_visited_map_width_m,
          par_.expl_visited_map_height_m,
          par_.expl_visited_map_resolution_m);

      pub_frontiers_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
          "exploration/frontiers", 10);
      pub_frontier_data_ = this->create_publisher<dynus_interfaces::msg::FrontierList>(
          "exploration/frontiers_data", 10);
      pub_explore_current_goal_ = this->create_publisher<geometry_msgs::msg::PoseStamped>(
          "exploration/current_goal", 10);
      pub_visited_map_ = this->create_publisher<nav_msgs::msg::OccupancyGrid>(
          "exploration/visited_map", rclcpp::QoS(1).transient_local());

      // MinPos peer pose sharing (global topic, all agents pub+sub)
      if (par_.expl_use_minpos) {
        auto peer_qos = rclcpp::QoS(10).reliable();
        pub_peer_pose_ = this->create_publisher<geometry_msgs::msg::PoseStamped>(
            "/exploration/peer_poses", peer_qos);
        sub_peer_pose_ = this->create_subscription<geometry_msgs::msg::PoseStamped>(
            "/exploration/peer_poses", peer_qos,
            [this](const geometry_msgs::msg::PoseStamped::SharedPtr msg) {
              // Filter out own messages
              if (msg->header.frame_id == ns_) return;
              peer_tracker_.updatePeer(
                  msg->header.frame_id,
                  msg->pose.position.x,
                  msg->pose.position.y,
                  rclcpp::Time(msg->header.stamp).seconds());
            },
            options_re_1);
        RCLCPP_INFO(this->get_logger(),
                    "Exploration MinPos: enabled, peer_timeout=%.1fs, publish_rate=%.1fHz",
                    par_.expl_peer_timeout_sec, par_.expl_peer_publish_rate_hz);
      }

      // Visited-map sharing: broadcast our persistent map so peers can skip
      // frontiers in areas we've already explored. Subscriber on cb_group_map_
      // so mergeFrom() is serialized with occ2DCallback / frontier detection.
      if (par_.expl_use_minpos) {
        pub_peer_visited_map_ = this->create_publisher<nav_msgs::msg::OccupancyGrid>(
            "/exploration/visited_maps", rclcpp::QoS(1).reliable());
        sub_peer_visited_map_ = this->create_subscription<nav_msgs::msg::OccupancyGrid>(
            "/exploration/visited_maps", rclcpp::QoS(1).reliable(),
            [this](const nav_msgs::msg::OccupancyGrid::SharedPtr msg) {
              if (msg->header.frame_id == ns_ || !visited_map_) return;
              visited_map_->mergeFrom(
                  msg->data.data(),
                  static_cast<int>(msg->info.width),
                  static_cast<int>(msg->info.height),
                  msg->info.origin.position.x,
                  msg->info.origin.position.y,
                  msg->info.resolution);
            },
            options_map);
        RCLCPP_INFO(this->get_logger(), "Exploration MinPos: visited-map sharing enabled");
      }

      // Global return-home trigger. One publish on /exploration/return_home
      // makes every subscribing agent drop new frontier work and head back to
      // its captured start position. Reliable QoS so a single shot reaches
      // every agent. cb_group_map_ keeps it serialized with explore-select.
      sub_return_home_ = this->create_subscription<std_msgs::msg::Empty>(
          "/exploration/return_home", rclcpp::QoS(1).reliable(),
          [this](const std_msgs::msg::Empty::SharedPtr) {
            if (home_return_requested_) {
              RCLCPP_INFO(this->get_logger(),
                          "Return-home trigger received, but already returning home");
              return;
            }
            if (!exploration_start_captured_) {
              RCLCPP_WARN(this->get_logger(),
                          "Return-home trigger received before exploration started — "
                          "no start pose captured, ignoring");
              return;
            }
            home_return_requested_ = true;
            // External command: come home and stay parked (do not auto-resume).
            home_return_resume_on_arrival_ = false;
            geometry_msgs::msg::PoseStamped home;
            home.header.frame_id    = par_.map_frame_id;
            home.header.stamp       = this->now();
            home.pose.position.x    = exploration_start_pos_.x();
            home.pose.position.y    = exploration_start_pos_.y();
            home.pose.position.z    = par_.expl_default_goal_z;
            home.pose.orientation.w = 1.0;
            terminalGoalCallbackImpl(home, /*from_user=*/false);
            exploration_active_  = false;
            RCLCPP_INFO(this->get_logger(),
                        "Return-home: heading to (%.2f, %.2f, %.2f)",
                        exploration_start_pos_.x(), exploration_start_pos_.y(),
                        par_.expl_default_goal_z);
          },
          options_map);

      const double rate_hz = std::max(0.1, par_.expl_select_rate_hz);
      const auto period = std::chrono::duration<double>(1.0 / rate_hz);
      timer_explore_select_ = this->create_wall_timer(
          std::chrono::duration_cast<std::chrono::nanoseconds>(period),
          std::bind(&MIGHTY_NODE::exploreSelectCallback, this), this->cb_group_map_);

      RCLCPP_INFO(this->get_logger(),
                  "Exploration: enabled, select rate=%.1f Hz, min_cells=%d, "
                  "merge_radius=%.2fm, visit_radius=%.2fm",
                  rate_hz, par_.expl_cluster_min_cells,
                  par_.expl_merge_radius_m, par_.expl_visit_radius_m);
    }
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Destructor
 */
MIGHTY_NODE::~MIGHTY_NODE() {
  // release the memory
  mighty_ptr_.reset();
}

// ----------------------------------------------------------------------------

/**
 * @brief Declare the parameters
 */
void MIGHTY_NODE::declareParameters() {
  // Sim enviroment
  this->declare_parameter("sim_env", "fake_sim");

  // UAV or Ground robot
  this->declare_parameter("vehicle_type", "uav");
  this->declare_parameter("provide_goal_in_global_frame", false);
  this->declare_parameter("use_hardware", false);
  this->declare_parameter("map_frame_id", "map");
  this->declare_parameter("share_traj", true);
  this->declare_parameter("use_frame_alignment", false);
  this->declare_parameter("num_agents", 10);
  this->declare_parameter("sim_frame_offset_qx", 0.0);
  this->declare_parameter("sim_frame_offset_qy", 0.0);
  this->declare_parameter("sim_frame_offset_qz", 0.0);
  this->declare_parameter("sim_frame_offset_qw", 1.0);

  // Formation flight
  this->declare_parameter("use_formation", false);
  this->declare_parameter("formation_weight", 0.0);
  this->declare_parameter("formation_self_offset", std::vector<double>{0.0, 0.0, 0.0});
  this->declare_parameter("formation_neighbor_ids", std::vector<int64_t>{});
  this->declare_parameter("formation_neighbor_offsets", std::vector<double>{});

  // Flight mode
  this->declare_parameter("flight_mode", "terminal_goal");

  // Visual
  this->declare_parameter("visual_level", 1);

  // Global planner parameters
  this->declare_parameter("file_path", "/home/kkondo/code/dynus_ws/src/dynus/data/data.txt");
  this->declare_parameter("use_benchmark", false);
  this->declare_parameter("start_yaw", -90.0);
  this->declare_parameter("global_planner", "sjps");
  this->declare_parameter("global_planner_verbose", false);
  this->declare_parameter("global_planner_heuristic_weight", 1.0);
  this->declare_parameter("factor_hgp", 1.0);
  this->declare_parameter("inflation_hgp", 0.5);
  this->declare_parameter("x_min", -100.0);
  this->declare_parameter("x_max", 100.0);
  this->declare_parameter("y_min", -100.0);
  this->declare_parameter("y_max", 100.0);
  this->declare_parameter("z_min", 0.0);
  this->declare_parameter("z_max", 5.0);
  this->declare_parameter("hgp_timeout_duration_ms", 1000);
  this->declare_parameter("max_expand", 10000);
  this->declare_parameter("use_free_start", false);
  this->declare_parameter("free_start_factor", 1.0);
  this->declare_parameter("use_free_goal", false);
  this->declare_parameter("free_goal_factor", 1.0);
  this->declare_parameter("relocate_occupied_goal", true);
  this->declare_parameter("max_dist_vertexes", 5.0);
  this->declare_parameter("w_unknown", 1.0);
  this->declare_parameter("w_align", 60.0);
  this->declare_parameter("decay_len_cells", 20.0);
  this->declare_parameter("w_side", 0.2);

  // LOS post processing parameters
  this->declare_parameter("los_cells", 3);
  this->declare_parameter("min_len",
                          0.5);  // [m] minimum length between two waypoints after post processing
  this->declare_parameter("min_turn", 10.0);  // [deg] minimum turn angle after post processing
  this->declare_parameter("heat_cutoff_ratio", 0.5);
  this->declare_parameter("disable_all_smoothing", false);
  this->declare_parameter("skip_path_smoothing", false);
  this->declare_parameter("smooth_iterations", 50);
  this->declare_parameter("smooth_alpha", 0.3);

  // Path push visualization parameters
  this->declare_parameter("use_state_update", true);
  this->declare_parameter("use_random_color_for_global_path", false);
  this->declare_parameter("use_path_push_for_visualization", false);

  // Decomposition parameters
  this->declare_parameter("local_box_size", std::vector<float>{2.0, 2.0, 2.0});
  this->declare_parameter("min_dist_from_agent_to_traj", 6.0);
  this->declare_parameter("use_shrinked_box", false);
  this->declare_parameter("shrinked_box_size", 0.2);
  this->declare_parameter("sfc_use_unknown_as_obstacle", false);

  // Map parameters
  this->declare_parameter("map_buffer", 6.0);
  this->declare_parameter("center_shift_factor", 0.5);
  this->declare_parameter("initial_wdx", 1.0);
  this->declare_parameter("initial_wdy", 1.0);
  this->declare_parameter("initial_wdz", 1.0);
  this->declare_parameter("min_wdx", 10.0);
  this->declare_parameter("min_wdy", 10.0);
  this->declare_parameter("min_wdz", 2.0);
  this->declare_parameter("mighty_map_res", 0.1);

  // Heat map parameters
  this->declare_parameter("use_heat_map", false);
  this->declare_parameter("heat_weight", 1.0);
  this->declare_parameter("dynamic_heat_enabled", false);
  this->declare_parameter("dynamic_as_occupied_current", true);
  this->declare_parameter("dynamic_as_occupied_future", false);
  this->declare_parameter("heat_alpha0", 1.0);
  this->declare_parameter("heat_alpha1", 2.0);
  this->declare_parameter("heat_p", 2);
  this->declare_parameter("heat_q", 2);
  this->declare_parameter("heat_tau_ratio", 0.5);
  this->declare_parameter("heat_gamma", 0.0);
  this->declare_parameter("heat_Hmax", 10.0);
  this->declare_parameter("obst_max_vel", 1.0);
  this->declare_parameter("dyn_base_inflation_m", 0.5);
  this->declare_parameter("dyn_heat_tube_radius_m", 2.0);
  this->declare_parameter("heat_num_samples", 15);
  this->declare_parameter("prediction_horizon", 3.0);
  this->declare_parameter("prediction_mask_distance", 2.0);
  this->declare_parameter("static_heat_enabled", false);
  this->declare_parameter("static_heat_alpha", 2.0);
  this->declare_parameter("static_heat_p", 2);
  this->declare_parameter("static_heat_Hmax", 50.0);
  this->declare_parameter("static_heat_rmax_m", 1.0);
  this->declare_parameter("static_heat_default_radius_m", 0.5);
  this->declare_parameter("static_heat_boundary_only", true);
  this->declare_parameter("static_heat_apply_on_unknown", false);
  this->declare_parameter("static_heat_exclude_dynamic", true);

  // Soft-cost obstacle parameters
  this->declare_parameter("use_soft_cost_obstacles", false);
  this->declare_parameter("obstacle_soft_cost", 100.0);

  // Communication delay parameters
  this->declare_parameter("use_comm_delay_inflation", true);
  this->declare_parameter("comm_delay_inflation_alpha", 0.2);
  this->declare_parameter("comm_delay_inflation_max", 0.5);
  this->declare_parameter("comm_delay_filter_alpha", 0.8);

  // Simulation parameters
  this->declare_parameter("depth_camera_depth_max", 10.0);
  this->declare_parameter("fov_visual_depth", 10.0);
  this->declare_parameter("fov_visual_x_deg", 10.0);
  this->declare_parameter("fov_visual_y_deg", 10.0);
  this->declare_parameter("use_sphere_sensing", false);
  this->declare_parameter("sphere_sensing_radius", 5.0);

  // Initial guess parameters
  this->declare_parameter("use_multiple_initial_guesses", true);
  this->declare_parameter("num_perturbation_for_ig", 8);
  this->declare_parameter("r_max_for_ig", 1.0);

  // Optimization parameters
  this->declare_parameter("num_N", 5);
  this->declare_parameter("horizon", 20.0);
  this->declare_parameter("dc", 0.01);
  this->declare_parameter("v_max", 1.0);
  this->declare_parameter("a_max", 1.0);
  this->declare_parameter("j_max", 1.0);
  this->declare_parameter("closed_form_traj_verbose", false);
  this->declare_parameter("jerk_weight", 1.0);
  this->declare_parameter("dynamic_weight", 1.0);
  this->declare_parameter("time_weight", 1.0);
  this->declare_parameter("pos_anchor_weight", 1.0);
  this->declare_parameter("stat_weight", 1.0);
  this->declare_parameter("dyn_constr_bodyrate_weight", 1.0);
  this->declare_parameter("dyn_constr_tilt_weight", 1.0);
  this->declare_parameter("dyn_constr_thrust_weight", 1.0);
  this->declare_parameter("dyn_constr_vel_weight", 1.0);
  this->declare_parameter("dyn_constr_acc_weight", 1.0);
  this->declare_parameter("dyn_constr_jerk_weight", 1.0);
  this->declare_parameter("num_dyn_obst_samples", 10);
  this->declare_parameter("planner_Co", 0.5);
  this->declare_parameter("planner_Cw", 1.0);
  this->declare_parameter("verbose_computation_time", false);
  this->declare_parameter("drone_bbox", std::vector<double>{0.5, 0.5, 0.5});
  this->declare_parameter("goal_radius", 0.5);
  this->declare_parameter("goal_seen_radius", 2.0);
  this->declare_parameter("init_turn_bf", 15.0);
  this->declare_parameter("integral_resolution", 30);
  this->declare_parameter("hinge_mu", 1e-2);
  this->declare_parameter("omega_max", 1e-2);
  this->declare_parameter("tilt_max_rad", 0.6);
  this->declare_parameter("f_min", 0.0);
  this->declare_parameter("f_max", 20.0);
  this->declare_parameter("mass", 1.0);
  this->declare_parameter("g", 9.81);
  this->declare_parameter("fopt_threshold", 0.1);
  this->declare_parameter("spline_degree", 5);

  // L-BFGS parameters
  this->declare_parameter("f_dec_coeff", 1e-2);
  this->declare_parameter("cautious_factor", 0.0);
  this->declare_parameter("past", 5);
  this->declare_parameter("max_linesearch", 64);
  this->declare_parameter("max_iterations", 30);
  this->declare_parameter("g_epsilon", 0.0);
  this->declare_parameter("delta", 1e-6);

  // Dynamic obstacles parameters
  this->declare_parameter("traj_lifetime", 10.0);

  // Dynamic k_value parameters
  this->declare_parameter("num_replanning_before_adapt", 10);
  this->declare_parameter("default_k_value", 150);
  this->declare_parameter("alpha_k_value_filtering", 0.8);
  this->declare_parameter("k_value_factor", 1.2);

  // Yaw-related parameters
  this->declare_parameter("alpha_filter_dyaw", 0.8);
  this->declare_parameter("w_max", 0.5);
  this->declare_parameter("yaw_spinning_threshold", 10);
  this->declare_parameter("yaw_spinning_dyaw", 0.1);

  // Simulation env parameters
  this->declare_parameter("force_goal_z", true);
  this->declare_parameter("default_goal_z", 2.5);

  // Debug flags
  this->declare_parameter("debug_verbose", false);
  this->declare_parameter("debug_log_file", std::string(""));  // path to per-publish/per-term-goal log; "" disables

  // Ground robot control parameters
  this->declare_parameter("ground_robot_kx", 0.1);
  this->declare_parameter("ground_robot_ky", 0.1);
  this->declare_parameter("ground_robot_kyaw", 1.0);
  this->declare_parameter("ground_robot_eps", 0.1);
  this->declare_parameter("ground_robot_v_max", 1.0);
  this->declare_parameter("ground_robot_w_max", 1.5);
  this->declare_parameter("ground_robot_L_min", 1.0);

  // 2D ground robot planning parameters
  this->declare_parameter("use_2d_planning", false);
  this->declare_parameter("robot_height", 0.5);
  this->declare_parameter("obstacle_min_height", 0.3);
  this->declare_parameter("use_column_any_occupied", true);
  this->declare_parameter("column_min_z", 0.15);
  this->declare_parameter("terrain_cost_weight", 1.0);
  this->declare_parameter("terrain_cost_mode", std::string("max"));
  this->declare_parameter("ground_slab_margin", 0.3);

  // ESDF-based obstacle avoidance (ground robot only)
  this->declare_parameter("use_esdf_cost", false);
  this->declare_parameter("use_persistent_global_map", false);
  this->declare_parameter("reuse_global_path", false);
  this->declare_parameter("global_path_max_age_sec", 2.0);
  this->declare_parameter("global_path_max_deviation_m", 1.0);
  this->declare_parameter("esdf_weight", 1e+3);
  this->declare_parameter("esdf_d_safe", 1.0);
  this->declare_parameter("esdf_truncation_distance", 10);

  // Bend pre-alignment (ground robot only)
  this->declare_parameter("corridor_hop_enabled", false);
  this->declare_parameter("corridor_corner_angle_deg", 75.0);
  this->declare_parameter("corridor_detection_window_m", 0.8);
  this->declare_parameter("corridor_backoff_m", 0.4);
  this->declare_parameter("corridor_min_leg_m", 0.3);
  this->declare_parameter("corridor_clearance_threshold_m", 0.7);
  this->declare_parameter("corridor_max_ascent_m", 0.0);
  this->declare_parameter("corridor_ascent_step_m", 0.05);
  this->declare_parameter("corridor_ascent_max_iters", 0);

  this->declare_parameter("trajectory_downsample_points", 500);
  this->declare_parameter("mpc_path_spacing", 0.05);

  // Frontier-based exploration (ground robot only).
  this->declare_parameter("exploration.enabled", false);
  this->declare_parameter("exploration.select_rate_hz", 1.0);
  this->declare_parameter("exploration.default_goal_z", 0.0);
  this->declare_parameter("exploration.select_nearest", false);
  this->declare_parameter("exploration.select_commit_margin_m", 0.5);
  this->declare_parameter("exploration.home_grace_sec", 3.0);
  this->declare_parameter("exploration.stall_timeout_sec", 30.0);
  this->declare_parameter("exploration.stall_eps_m", 0.15);
  this->declare_parameter("exploration.detector.cluster_min_cells", 6);
  this->declare_parameter("exploration.detector.border_margin_cells", 2);
  this->declare_parameter("exploration.detector.obstacle_clearance_cells", 1);
  this->declare_parameter("exploration.detector.robot_snap_radius_m", 1.0);
  this->declare_parameter("exploration.detector.unknown_bridge_cells", 0);
  this->declare_parameter("exploration.detector.reconcile_stale", true);
  this->declare_parameter("exploration.bounds.enabled", false);
  this->declare_parameter("exploration.bounds.min_x", -50.0);
  this->declare_parameter("exploration.bounds.max_x",  50.0);
  this->declare_parameter("exploration.bounds.min_y", -50.0);
  this->declare_parameter("exploration.bounds.max_y",  50.0);
  this->declare_parameter("exploration.detector.min_obstacle_distance_m", 0.0);
  this->declare_parameter("exploration.utility.w_size", 1.0);
  this->declare_parameter("exploration.utility.w_dist", 2.0);
  this->declare_parameter("exploration.utility.w_info", 1.0);
  this->declare_parameter("exploration.utility.w_revisit", 0.5);
  this->declare_parameter("exploration.utility.w_heading", 0.3);
  this->declare_parameter("exploration.utility.size_ref_m2", 5.0);
  this->declare_parameter("exploration.utility.dist_ref_m", 25.0);
  this->declare_parameter("exploration.utility.sensor_radius_m", 5.0);
  this->declare_parameter("exploration.utility.goal_select_threshold", -1.0e9);
  this->declare_parameter("exploration.manager.merge_radius_m", 1.0);
  this->declare_parameter("exploration.manager.centroid_ema_alpha", 0.5);
  this->declare_parameter("exploration.manager.visit_radius_m", 2.0);
  this->declare_parameter("exploration.manager.visit_dwell_sec", 1.0);
  this->declare_parameter("exploration.manager.verify_radius_cells", 2);
  this->declare_parameter("exploration.manager.max_frontiers", 1000);
  this->declare_parameter("exploration.manager.unreachable_consec_thresh", 5);
  this->declare_parameter("exploration.manager.pursuit_timeout_factor", 10.0);
  this->declare_parameter("exploration.manager.pursuit_timeout_v_ref", 0.5);
  this->declare_parameter("exploration.manager.pursuit_timeout_min_sec", 10.0);
  this->declare_parameter("exploration.manager.pursuit_timeout_max_strikes", 3);
  this->declare_parameter("exploration.manager.invalidation_keep_out_radius_m", 1.5);
  this->declare_parameter("exploration.manager.invalidation_cooldown_sec", 30.0);
  this->declare_parameter("exploration.visited_map.center_x", 0.0);
  this->declare_parameter("exploration.visited_map.center_y", 0.0);
  this->declare_parameter("exploration.visited_map.width_m", 100.0);
  this->declare_parameter("exploration.visited_map.height_m", 100.0);
  this->declare_parameter("exploration.visited_map.resolution_m", 0.15);
  this->declare_parameter("exploration.visited_map.publish", true);
  this->declare_parameter("exploration.visited_map.fuse_into_local", true);
  this->declare_parameter("exploration.visited_map.detect_on_visited_map", true);
  this->declare_parameter("exploration.visualization.publish_markers", true);
  // MinPos multi-robot frontier allocation
  this->declare_parameter("exploration.minpos.enabled", false);
  this->declare_parameter("exploration.minpos.peer_timeout_sec", 5.0);
  this->declare_parameter("exploration.minpos.peer_publish_rate_hz", 5.0);
  this->declare_parameter("exploration.minpos.min_frontier_dist_to_peers_m", 0.0);
  this->declare_parameter("exploration.minpos.peer_visit_radius_m", 2.0);
  this->declare_parameter("exploration.external_selector", false);
}

// ----------------------------------------------------------------------------

/**
 * @brief Set the parameters
 */
void MIGHTY_NODE::setParameters() {
  // Set the parameters

  // Sim enviroment
  par_.sim_env = this->get_parameter("sim_env").as_string();

  // Vehicle type (UAV, Wheeled Robit, or Quadruped)
  par_.vehicle_type = this->get_parameter("vehicle_type").as_string();
  par_.provide_goal_in_global_frame = this->get_parameter("provide_goal_in_global_frame").as_bool();
  par_.use_hardware = this->get_parameter("use_hardware").as_bool();
  par_.map_frame_id = this->get_parameter("map_frame_id").as_string();
  par_.share_traj = this->get_parameter("share_traj").as_bool();
  par_.use_frame_alignment = this->get_parameter("use_frame_alignment").as_bool();
  par_.num_agents = this->get_parameter("num_agents").as_int();
  par_.sim_frame_offset_qx = this->get_parameter("sim_frame_offset_qx").as_double();
  par_.sim_frame_offset_qy = this->get_parameter("sim_frame_offset_qy").as_double();
  par_.sim_frame_offset_qz = this->get_parameter("sim_frame_offset_qz").as_double();
  par_.sim_frame_offset_qw = this->get_parameter("sim_frame_offset_qw").as_double();

  // Formation flight
  par_.use_formation = this->get_parameter("use_formation").as_bool();
  par_.formation_weight = this->get_parameter("formation_weight").as_double();
  par_.formation_self_offset = this->get_parameter("formation_self_offset").as_double_array();
  par_.formation_neighbor_ids = this->get_parameter("formation_neighbor_ids").as_integer_array();
  par_.formation_neighbor_offsets =
      this->get_parameter("formation_neighbor_offsets").as_double_array();
  if (par_.use_formation) {
    if (par_.formation_self_offset.size() != 3) {
      RCLCPP_FATAL(this->get_logger(),
                   "formation_self_offset must have length 3 (got %zu); shutting down",
                   par_.formation_self_offset.size());
      rclcpp::shutdown();
    }
    if (par_.formation_neighbor_offsets.size() != 3 * par_.formation_neighbor_ids.size()) {
      RCLCPP_FATAL(this->get_logger(),
                   "formation_neighbor_offsets length (%zu) must be 3 * "
                   "formation_neighbor_ids length (%zu); shutting down",
                   par_.formation_neighbor_offsets.size(),
                   par_.formation_neighbor_ids.size());
      rclcpp::shutdown();
    }
  }

  // Build the sim frame offset matrix
  {
    Eigen::Quaterniond q(par_.sim_frame_offset_qw, par_.sim_frame_offset_qx,
                         par_.sim_frame_offset_qy, par_.sim_frame_offset_qz);
    sim_frame_offset_ = Eigen::Matrix4d::Identity();
    sim_frame_offset_.block<3, 3>(0, 0) = q.normalized().toRotationMatrix();
    // Inverse: for pure rotation, R^-1 = R^T
    sim_frame_offset_inv_ = Eigen::Matrix4d::Identity();
    sim_frame_offset_inv_.block<3, 3>(0, 0) = q.normalized().toRotationMatrix().transpose();
  }

  // Flight mode
  par_.flight_mode = this->get_parameter("flight_mode").as_string();

  // Visual level
  par_.visual_level = this->get_parameter("visual_level").as_int();

  // Global Planner parameters
  file_path_ = this->get_parameter("file_path").as_string();
  use_benchmark_ = this->get_parameter("use_benchmark").as_bool();
  par_.global_planner = this->get_parameter("global_planner").as_string();
  par_.global_planner_verbose = this->get_parameter("global_planner_verbose").as_bool();
  par_.global_planner_heuristic_weight =
      this->get_parameter("global_planner_heuristic_weight").as_double();
  par_.factor_hgp = this->get_parameter("factor_hgp").as_double();
  par_.inflation_hgp = this->get_parameter("inflation_hgp").as_double();
  par_.x_min = this->get_parameter("x_min").as_double();
  par_.x_max = this->get_parameter("x_max").as_double();
  par_.y_min = this->get_parameter("y_min").as_double();
  par_.y_max = this->get_parameter("y_max").as_double();
  par_.z_min = this->get_parameter("z_min").as_double();
  par_.z_max = this->get_parameter("z_max").as_double();
  par_.hgp_timeout_duration_ms = this->get_parameter("hgp_timeout_duration_ms").as_int();
  par_.max_expand = this->get_parameter("max_expand").as_int();
  par_.max_num_expansion = par_.max_expand;

  par_.use_free_start = this->get_parameter("use_free_start").as_bool();
  par_.free_start_factor = this->get_parameter("free_start_factor").as_double();
  par_.use_free_goal = this->get_parameter("use_free_goal").as_bool();
  par_.free_goal_factor = this->get_parameter("free_goal_factor").as_double();
  par_.relocate_occupied_goal = this->get_parameter("relocate_occupied_goal").as_bool();
  par_.max_dist_vertexes = this->get_parameter("max_dist_vertexes").as_double();
  par_.w_unknown = this->get_parameter("w_unknown").as_double();
  par_.w_align = this->get_parameter("w_align").as_double();
  par_.decay_len_cells = this->get_parameter("decay_len_cells").as_double();
  par_.w_side = this->get_parameter("w_side").as_double();

  // LOS post processing parameters
  par_.los_cells = this->get_parameter("los_cells").as_int();
  par_.min_len = this->get_parameter("min_len").as_double();
  par_.min_turn = this->get_parameter("min_turn").as_double();
  par_.heat_cutoff_ratio = this->get_parameter("heat_cutoff_ratio").as_double();
  par_.disable_all_smoothing = this->get_parameter("disable_all_smoothing").as_bool();
  par_.skip_path_smoothing = this->get_parameter("skip_path_smoothing").as_bool();
  par_.smooth_iterations = this->get_parameter("smooth_iterations").as_int();
  par_.smooth_alpha = this->get_parameter("smooth_alpha").as_double();

  // Path push visualization parameters
  par_.use_state_update = this->get_parameter("use_state_update").as_bool();
  par_.use_random_color_for_global_path =
      this->get_parameter("use_random_color_for_global_path").as_bool();
  par_.use_path_push_for_visualization =
      this->get_parameter("use_path_push_for_visualization").as_bool();

  // Static obstacle push parameters

  // Decomposition parameters
  par_.local_box_size = this->get_parameter("local_box_size").as_double_array();
  par_.min_dist_from_agent_to_traj = this->get_parameter("min_dist_from_agent_to_traj").as_double();
  par_.use_shrinked_box = this->get_parameter("use_shrinked_box").as_bool();
  par_.shrinked_box_size = this->get_parameter("shrinked_box_size").as_double();
  par_.sfc_use_unknown_as_obstacle = this->get_parameter("sfc_use_unknown_as_obstacle").as_bool();

  // Map parameters
  par_.map_buffer = this->get_parameter("map_buffer").as_double();
  par_.center_shift_factor = this->get_parameter("center_shift_factor").as_double();
  par_.initial_wdx = this->get_parameter("initial_wdx").as_double();
  par_.initial_wdy = this->get_parameter("initial_wdy").as_double();
  par_.initial_wdz = this->get_parameter("initial_wdz").as_double();
  par_.min_wdx = this->get_parameter("min_wdx").as_double();
  par_.min_wdy = this->get_parameter("min_wdy").as_double();
  par_.min_wdz = this->get_parameter("min_wdz").as_double();
  par_.res = this->get_parameter("mighty_map_res").as_double();

  // Heat map parameters
  par_.use_heat_map = this->get_parameter("use_heat_map").as_bool();
  par_.heat_weight = this->get_parameter("heat_weight").as_double();
  par_.dynamic_heat_enabled = this->get_parameter("dynamic_heat_enabled").as_bool();
  par_.dynamic_as_occupied_current = this->get_parameter("dynamic_as_occupied_current").as_bool();
  par_.dynamic_as_occupied_future = this->get_parameter("dynamic_as_occupied_future").as_bool();
  par_.heat_alpha0 = this->get_parameter("heat_alpha0").as_double();
  par_.heat_alpha1 = this->get_parameter("heat_alpha1").as_double();
  par_.heat_p = this->get_parameter("heat_p").as_int();
  par_.heat_q = this->get_parameter("heat_q").as_int();
  par_.heat_tau_ratio = this->get_parameter("heat_tau_ratio").as_double();
  par_.heat_gamma = this->get_parameter("heat_gamma").as_double();
  par_.heat_Hmax = this->get_parameter("heat_Hmax").as_double();
  par_.obst_max_vel = this->get_parameter("obst_max_vel").as_double();
  par_.dyn_base_inflation_m = this->get_parameter("dyn_base_inflation_m").as_double();
  par_.dyn_heat_tube_radius_m = this->get_parameter("dyn_heat_tube_radius_m").as_double();
  par_.heat_num_samples = this->get_parameter("heat_num_samples").as_int();
  par_.prediction_horizon = this->get_parameter("prediction_horizon").as_double();
  par_.prediction_mask_distance = this->get_parameter("prediction_mask_distance").as_double();
  par_.static_heat_enabled = this->get_parameter("static_heat_enabled").as_bool();
  par_.static_heat_alpha = this->get_parameter("static_heat_alpha").as_double();
  par_.static_heat_p = this->get_parameter("static_heat_p").as_int();
  par_.static_heat_Hmax = this->get_parameter("static_heat_Hmax").as_double();
  par_.static_heat_rmax_m = this->get_parameter("static_heat_rmax_m").as_double();
  par_.static_heat_default_radius_m =
      this->get_parameter("static_heat_default_radius_m").as_double();
  par_.static_heat_boundary_only = this->get_parameter("static_heat_boundary_only").as_bool();
  par_.static_heat_apply_on_unknown = this->get_parameter("static_heat_apply_on_unknown").as_bool();
  par_.static_heat_exclude_dynamic = this->get_parameter("static_heat_exclude_dynamic").as_bool();

  // Soft-cost obstacle parameters
  par_.use_soft_cost_obstacles = this->get_parameter("use_soft_cost_obstacles").as_bool();
  par_.obstacle_soft_cost = this->get_parameter("obstacle_soft_cost").as_double();

  // Communication delay parameters
  par_.use_comm_delay_inflation = this->get_parameter("use_comm_delay_inflation").as_bool();
  par_.comm_delay_inflation_alpha = this->get_parameter("comm_delay_inflation_alpha").as_double();
  par_.comm_delay_inflation_max = this->get_parameter("comm_delay_inflation_max").as_double();
  par_.comm_delay_filter_alpha = this->get_parameter("comm_delay_filter_alpha").as_double();

  // Simulation parameters
  par_.depth_camera_depth_max = this->get_parameter("depth_camera_depth_max").as_double();
  par_.fov_visual_depth = this->get_parameter("fov_visual_depth").as_double();
  par_.fov_visual_x_deg = this->get_parameter("fov_visual_x_deg").as_double();
  par_.fov_visual_y_deg = this->get_parameter("fov_visual_y_deg").as_double();
  par_.use_sphere_sensing = this->get_parameter("use_sphere_sensing").as_bool();
  par_.sphere_sensing_radius = this->get_parameter("sphere_sensing_radius").as_double();

  // Initial guess parameters
  par_.use_multiple_initial_guesses = this->get_parameter("use_multiple_initial_guesses").as_bool();
  par_.num_perturbation_for_ig = this->get_parameter("num_perturbation_for_ig").as_int();
  par_.r_max_for_ig = this->get_parameter("r_max_for_ig").as_double();

  // Optimization parameters
  par_.num_N = this->get_parameter("num_N").as_int();
  par_.horizon = this->get_parameter("horizon").as_double();
  par_.dc = this->get_parameter("dc").as_double();
  par_.v_max = this->get_parameter("v_max").as_double();
  par_.a_max = this->get_parameter("a_max").as_double();
  par_.j_max = this->get_parameter("j_max").as_double();
  par_.closed_form_traj_verbose = this->get_parameter("closed_form_traj_verbose").as_bool();
  par_.jerk_weight = this->get_parameter("jerk_weight").as_double();
  par_.dynamic_weight = this->get_parameter("dynamic_weight").as_double();
  par_.time_weight = this->get_parameter("time_weight").as_double();
  par_.pos_anchor_weight = this->get_parameter("pos_anchor_weight").as_double();
  par_.stat_weight = this->get_parameter("stat_weight").as_double();
  par_.dyn_constr_bodyrate_weight = this->get_parameter("dyn_constr_bodyrate_weight").as_double();
  par_.dyn_constr_tilt_weight = this->get_parameter("dyn_constr_tilt_weight").as_double();
  par_.dyn_constr_thrust_weight = this->get_parameter("dyn_constr_thrust_weight").as_double();
  par_.dyn_constr_vel_weight = this->get_parameter("dyn_constr_vel_weight").as_double();
  par_.dyn_constr_acc_weight = this->get_parameter("dyn_constr_acc_weight").as_double();
  par_.dyn_constr_jerk_weight = this->get_parameter("dyn_constr_jerk_weight").as_double();
  par_.num_dyn_obst_samples = this->get_parameter("num_dyn_obst_samples").as_int();
  par_.planner_Co = this->get_parameter("planner_Co").as_double();
  par_.planner_Cw = this->get_parameter("planner_Cw").as_double();
  verbose_computation_time_ = this->get_parameter("verbose_computation_time").as_bool();
  par_.drone_bbox = this->get_parameter("drone_bbox").as_double_array();
  par_.drone_radius = par_.drone_bbox[0] / 2.0;
  par_.goal_radius = this->get_parameter("goal_radius").as_double();
  par_.goal_seen_radius = this->get_parameter("goal_seen_radius").as_double();
  par_.init_turn_bf = this->get_parameter("init_turn_bf").as_double();
  par_.integral_resolution = this->get_parameter("integral_resolution").as_int();
  par_.hinge_mu = this->get_parameter("hinge_mu").as_double();
  par_.omega_max = this->get_parameter("omega_max").as_double();
  par_.tilt_max_rad = this->get_parameter("tilt_max_rad").as_double();
  par_.f_min = this->get_parameter("f_min").as_double();
  par_.f_max = this->get_parameter("f_max").as_double();
  par_.mass = this->get_parameter("mass").as_double();
  par_.g = this->get_parameter("g").as_double();
  par_.fopt_threshold = this->get_parameter("fopt_threshold").as_double();
  par_.spline_degree = this->get_parameter("spline_degree").as_int();

  // L-BFGS parameters
  par_.f_dec_coeff = this->get_parameter("f_dec_coeff").as_double();
  par_.cautious_factor = this->get_parameter("cautious_factor").as_double();
  par_.past = this->get_parameter("past").as_int();
  par_.max_linesearch = this->get_parameter("max_linesearch").as_int();
  par_.max_iterations = this->get_parameter("max_iterations").as_int();
  par_.g_epsilon = this->get_parameter("g_epsilon").as_double();
  par_.delta = this->get_parameter("delta").as_double();

  // Dynamic obstacles parameters
  par_.traj_lifetime = this->get_parameter("traj_lifetime").as_double();

  // Dynamic k_value parameters
  par_.num_replanning_before_adapt = this->get_parameter("num_replanning_before_adapt").as_int();
  par_.default_k_value = this->get_parameter("default_k_value").as_int();
  par_.alpha_k_value_filtering = this->get_parameter("alpha_k_value_filtering").as_double();
  par_.k_value_factor = this->get_parameter("k_value_factor").as_double();

  // Yaw-related parameters
  par_.alpha_filter_dyaw = this->get_parameter("alpha_filter_dyaw").as_double();
  par_.w_max = this->get_parameter("w_max").as_double();
  par_.yaw_spinning_threshold = this->get_parameter("yaw_spinning_threshold").as_int();
  par_.yaw_spinning_dyaw = this->get_parameter("yaw_spinning_dyaw").as_double();

  // Simulation env parameters
  par_.force_goal_z = this->get_parameter("force_goal_z").as_bool();
  par_.default_goal_z = this->get_parameter("default_goal_z").as_double();

  if (par_.default_goal_z <= par_.z_min) {
    RCLCPP_ERROR(this->get_logger(), "Default goal z is lower than the ground level");
  }

  if (par_.default_goal_z >= par_.z_max) {
    RCLCPP_ERROR(this->get_logger(), "Default goal z is higher than the max level");
  }

  // Debug flags
  par_.debug_verbose = this->get_parameter("debug_verbose").as_bool();
  {
    auto path = this->get_parameter("debug_log_file").as_string();
    debug_log_open(path);
    if (!path.empty()) {
      RCLCPP_INFO(this->get_logger(), "Debug log file: %s", path.c_str());
    }
  }

  // Ground robot control parameters
  par_.ground_robot_kx = this->get_parameter("ground_robot_kx").as_double();
  par_.ground_robot_ky = this->get_parameter("ground_robot_ky").as_double();
  par_.ground_robot_kyaw = this->get_parameter("ground_robot_kyaw").as_double();
  par_.ground_robot_eps = this->get_parameter("ground_robot_eps").as_double();
  par_.ground_robot_v_max = this->get_parameter("ground_robot_v_max").as_double();
  par_.ground_robot_w_max = this->get_parameter("ground_robot_w_max").as_double();
  par_.ground_robot_L_min = this->get_parameter("ground_robot_L_min").as_double();

  // 2D ground robot planning parameters
  par_.use_2d_planning = this->get_parameter("use_2d_planning").as_bool();
  par_.robot_height = this->get_parameter("robot_height").as_double();
  par_.obstacle_min_height = this->get_parameter("obstacle_min_height").as_double();
  par_.use_column_any_occupied = this->get_parameter("use_column_any_occupied").as_bool();
  par_.column_min_z = this->get_parameter("column_min_z").as_double();
  par_.terrain_cost_weight = this->get_parameter("terrain_cost_weight").as_double();
  par_.terrain_cost_mode = this->get_parameter("terrain_cost_mode").as_string();
  par_.ground_slab_margin = this->get_parameter("ground_slab_margin").as_double();

  par_.use_esdf_cost = this->get_parameter("use_esdf_cost").as_bool();
  par_.use_persistent_global_map = this->get_parameter("use_persistent_global_map").as_bool();
  par_.reuse_global_path = this->get_parameter("reuse_global_path").as_bool();
  par_.global_path_max_age_sec = this->get_parameter("global_path_max_age_sec").as_double();
  par_.global_path_max_deviation_m = this->get_parameter("global_path_max_deviation_m").as_double();
  par_.esdf_weight = this->get_parameter("esdf_weight").as_double();
  par_.esdf_d_safe = this->get_parameter("esdf_d_safe").as_double();
  par_.esdf_truncation_distance = this->get_parameter("esdf_truncation_distance").as_int();

  par_.corridor_hop_enabled = this->get_parameter("corridor_hop_enabled").as_bool();
  par_.corridor_corner_angle_deg = this->get_parameter("corridor_corner_angle_deg").as_double();
  par_.corridor_detection_window_m = this->get_parameter("corridor_detection_window_m").as_double();
  par_.corridor_backoff_m = this->get_parameter("corridor_backoff_m").as_double();
  par_.corridor_min_leg_m = this->get_parameter("corridor_min_leg_m").as_double();
  par_.corridor_clearance_threshold_m = this->get_parameter("corridor_clearance_threshold_m").as_double();
  par_.corridor_max_ascent_m = this->get_parameter("corridor_max_ascent_m").as_double();
  par_.corridor_ascent_step_m = this->get_parameter("corridor_ascent_step_m").as_double();
  par_.corridor_ascent_max_iters = this->get_parameter("corridor_ascent_max_iters").as_int();

  par_.trajectory_downsample_points = this->get_parameter("trajectory_downsample_points").as_int();
  par_.mpc_path_spacing = this->get_parameter("mpc_path_spacing").as_double();

  // Frontier-based exploration
  par_.expl_enabled              = this->get_parameter("exploration.enabled").as_bool();
  par_.expl_select_rate_hz       = this->get_parameter("exploration.select_rate_hz").as_double();
  par_.expl_default_goal_z       = this->get_parameter("exploration.default_goal_z").as_double();
  par_.expl_select_nearest       = this->get_parameter("exploration.select_nearest").as_bool();
  par_.expl_select_commit_margin_m = this->get_parameter("exploration.select_commit_margin_m").as_double();
  par_.expl_home_grace_sec = this->get_parameter("exploration.home_grace_sec").as_double();
  par_.expl_stall_timeout_sec = this->get_parameter("exploration.stall_timeout_sec").as_double();
  par_.expl_stall_eps_m = this->get_parameter("exploration.stall_eps_m").as_double();
  par_.expl_cluster_min_cells    = this->get_parameter("exploration.detector.cluster_min_cells").as_int();
  par_.expl_border_margin_cells  = this->get_parameter("exploration.detector.border_margin_cells").as_int();
  par_.expl_obstacle_clearance_cells =
      this->get_parameter("exploration.detector.obstacle_clearance_cells").as_int();
  par_.expl_robot_snap_radius_m  = this->get_parameter("exploration.detector.robot_snap_radius_m").as_double();
  par_.expl_unknown_bridge_cells =
      this->get_parameter("exploration.detector.unknown_bridge_cells").as_int();
  par_.expl_reconcile_stale      = this->get_parameter("exploration.detector.reconcile_stale").as_bool();
  par_.expl_bounds_enabled       = this->get_parameter("exploration.bounds.enabled").as_bool();
  par_.expl_bounds_min_x         = this->get_parameter("exploration.bounds.min_x").as_double();
  par_.expl_bounds_max_x         = this->get_parameter("exploration.bounds.max_x").as_double();
  par_.expl_bounds_min_y         = this->get_parameter("exploration.bounds.min_y").as_double();
  par_.expl_bounds_max_y         = this->get_parameter("exploration.bounds.max_y").as_double();
  par_.expl_min_obstacle_distance_m =
      this->get_parameter("exploration.detector.min_obstacle_distance_m").as_double();
  par_.expl_w_size               = this->get_parameter("exploration.utility.w_size").as_double();
  par_.expl_w_dist               = this->get_parameter("exploration.utility.w_dist").as_double();
  par_.expl_w_info               = this->get_parameter("exploration.utility.w_info").as_double();
  par_.expl_w_revisit            = this->get_parameter("exploration.utility.w_revisit").as_double();
  par_.expl_w_heading            = this->get_parameter("exploration.utility.w_heading").as_double();
  par_.expl_size_ref_m2          = this->get_parameter("exploration.utility.size_ref_m2").as_double();
  par_.expl_dist_ref_m           = this->get_parameter("exploration.utility.dist_ref_m").as_double();
  par_.expl_sensor_radius_m      = this->get_parameter("exploration.utility.sensor_radius_m").as_double();
  par_.expl_goal_select_threshold = this->get_parameter("exploration.utility.goal_select_threshold").as_double();
  par_.expl_merge_radius_m       = this->get_parameter("exploration.manager.merge_radius_m").as_double();
  par_.expl_centroid_ema_alpha   = this->get_parameter("exploration.manager.centroid_ema_alpha").as_double();
  par_.expl_visit_radius_m       = this->get_parameter("exploration.manager.visit_radius_m").as_double();
  par_.expl_visit_dwell_sec      = this->get_parameter("exploration.manager.visit_dwell_sec").as_double();
  par_.expl_verify_radius_cells  = this->get_parameter("exploration.manager.verify_radius_cells").as_int();
  par_.expl_max_frontiers        = this->get_parameter("exploration.manager.max_frontiers").as_int();
  par_.expl_unreachable_consec_thresh =
      this->get_parameter("exploration.manager.unreachable_consec_thresh").as_int();
  par_.expl_pursuit_timeout_factor =
      this->get_parameter("exploration.manager.pursuit_timeout_factor").as_double();
  par_.expl_pursuit_timeout_v_ref =
      this->get_parameter("exploration.manager.pursuit_timeout_v_ref").as_double();
  par_.expl_pursuit_timeout_min_sec =
      this->get_parameter("exploration.manager.pursuit_timeout_min_sec").as_double();
  par_.expl_pursuit_timeout_max_strikes =
      this->get_parameter("exploration.manager.pursuit_timeout_max_strikes").as_int();
  par_.expl_invalidation_keep_out_radius_m =
      this->get_parameter("exploration.manager.invalidation_keep_out_radius_m").as_double();
  par_.expl_invalidation_cooldown_sec =
      this->get_parameter("exploration.manager.invalidation_cooldown_sec").as_double();
  par_.expl_visited_map_center_x   = this->get_parameter("exploration.visited_map.center_x").as_double();
  par_.expl_visited_map_center_y   = this->get_parameter("exploration.visited_map.center_y").as_double();
  par_.expl_visited_map_width_m    = this->get_parameter("exploration.visited_map.width_m").as_double();
  par_.expl_visited_map_height_m   = this->get_parameter("exploration.visited_map.height_m").as_double();
  par_.expl_visited_map_resolution_m = this->get_parameter("exploration.visited_map.resolution_m").as_double();
  par_.expl_publish_visited_map    = this->get_parameter("exploration.visited_map.publish").as_bool();
  par_.expl_fuse_persistent_into_local =
      this->get_parameter("exploration.visited_map.fuse_into_local").as_bool();
  par_.expl_detect_on_visited_map =
      this->get_parameter("exploration.visited_map.detect_on_visited_map").as_bool();
  par_.expl_publish_markers      = this->get_parameter("exploration.visualization.publish_markers").as_bool();
  par_.expl_use_minpos           = this->get_parameter("exploration.minpos.enabled").as_bool();
  par_.expl_peer_timeout_sec     = this->get_parameter("exploration.minpos.peer_timeout_sec").as_double();
  par_.expl_peer_publish_rate_hz = this->get_parameter("exploration.minpos.peer_publish_rate_hz").as_double();
  par_.expl_min_frontier_dist_to_peers_m =
      this->get_parameter("exploration.minpos.min_frontier_dist_to_peers_m").as_double();
  par_.expl_peer_visit_radius_m =
      this->get_parameter("exploration.minpos.peer_visit_radius_m").as_double();
  par_.expl_external_selector =
      this->get_parameter("exploration.external_selector").as_bool();
}

// ----------------------------------------------------------------------------

/**
 * @brief Print the parameters
 */
void MIGHTY_NODE::printParameters() {
  // Print the parameters

  // Sim enviroment
  RCLCPP_INFO(this->get_logger(), "Sim Enviroment: %s", par_.sim_env.c_str());

  // Vehicle type (UAV, Wheeled Robit, or Quadruped)
  RCLCPP_INFO(this->get_logger(), "Vehicle Type: %d", par_.vehicle_type);
  RCLCPP_INFO(this->get_logger(), "Provide Goal in Global Frame: %d",
              par_.provide_goal_in_global_frame);
  RCLCPP_INFO(this->get_logger(), "Use Hardware: %d", par_.use_hardware);
  RCLCPP_INFO(this->get_logger(), "Share Traj: %d", par_.share_traj);
  RCLCPP_INFO(this->get_logger(), "Use Frame Alignment: %d", par_.use_frame_alignment);
  RCLCPP_INFO(this->get_logger(), "Num Agents: %d", par_.num_agents);

  // Flight mode
  RCLCPP_INFO(this->get_logger(), "Flight Mode: %s", par_.flight_mode.c_str());

  // Visual
  RCLCPP_INFO(this->get_logger(), "Visual Level: %d", par_.visual_level);

  // HGP parameters
  RCLCPP_INFO(this->get_logger(), "File Path: %s", file_path_.c_str());
  RCLCPP_INFO(this->get_logger(), "Perform Benchmark?: %d", use_benchmark_);
  RCLCPP_INFO(this->get_logger(), "Initial Guess Planner: %s", par_.global_planner.c_str());
  RCLCPP_INFO(this->get_logger(), "HGP Planner Verbose: %d", par_.global_planner_verbose);
  RCLCPP_INFO(this->get_logger(), "Global Planner Heuristic Weight: %f",
              par_.global_planner_heuristic_weight);
  RCLCPP_INFO(this->get_logger(), "Factor HGP: %f", par_.factor_hgp);
  RCLCPP_INFO(this->get_logger(), "Inflation HGP: %f", par_.inflation_hgp);
  RCLCPP_INFO(this->get_logger(), "X Min: %f", par_.x_min);
  RCLCPP_INFO(this->get_logger(), "X Max: %f", par_.x_max);
  RCLCPP_INFO(this->get_logger(), "Y Min: %f", par_.y_min);
  RCLCPP_INFO(this->get_logger(), "Y Max: %f", par_.y_max);
  RCLCPP_INFO(this->get_logger(), "Z Ground: %f", par_.z_min);
  RCLCPP_INFO(this->get_logger(), "Z Max: %f", par_.z_max);
  RCLCPP_INFO(this->get_logger(), "HGP Timeout Duration: %d", par_.hgp_timeout_duration_ms);
  RCLCPP_INFO(this->get_logger(), "Use Free Start?: %d", par_.use_free_start);
  RCLCPP_INFO(this->get_logger(), "Free Start Factor: %f", par_.free_start_factor);
  RCLCPP_INFO(this->get_logger(), "Use Free Goal?: %d", par_.use_free_goal);
  RCLCPP_INFO(this->get_logger(), "Free Goal Factor: %f", par_.free_goal_factor);
  RCLCPP_INFO(this->get_logger(), "Relocate Occupied Goal?: %d", par_.relocate_occupied_goal);
  RCLCPP_INFO(this->get_logger(), "max_dist_vertexes: %f", par_.max_dist_vertexes);
  RCLCPP_INFO(this->get_logger(), "w_unknown: %f", par_.w_unknown);
  RCLCPP_INFO(this->get_logger(), "w_align: %f", par_.w_align);
  RCLCPP_INFO(this->get_logger(), "decay_len_cells: %f", par_.decay_len_cells);
  RCLCPP_INFO(this->get_logger(), "w_side: %f", par_.w_side);

  // LOS post processing parameters
  RCLCPP_INFO(this->get_logger(), "LOS Cells: %d", par_.los_cells);
  RCLCPP_INFO(this->get_logger(), "Min Len: %f", par_.min_len);
  RCLCPP_INFO(this->get_logger(), "Min Turn: %f", par_.min_turn);
  RCLCPP_INFO(this->get_logger(), "Skip Path Smoothing: %d", par_.skip_path_smoothing);

  // Path push visualization parameters
  RCLCPP_INFO(this->get_logger(), "Use State Update?: %d", par_.use_state_update);
  RCLCPP_INFO(this->get_logger(), "Use Random Color for Global Path?: %d",
              par_.use_random_color_for_global_path);
  RCLCPP_INFO(this->get_logger(), "Use Path Push for Paper?: %d",
              par_.use_path_push_for_visualization);

  // Static obstacle push parameters

  RCLCPP_INFO(this->get_logger(), "Local Box Size: (%f, %f, %f)", par_.local_box_size[0],
              par_.local_box_size[1], par_.local_box_size[2]);
  RCLCPP_INFO(this->get_logger(), "Min Dist from Agent to Traj: %f",
              par_.min_dist_from_agent_to_traj);
  RCLCPP_INFO(this->get_logger(), "Use Shrinked Box: %d", par_.use_shrinked_box);
  RCLCPP_INFO(this->get_logger(), "Shrinked Box Size: %f", par_.shrinked_box_size);
  RCLCPP_INFO(this->get_logger(), "SFC Use Unknown As Obstacle: %d",
              par_.sfc_use_unknown_as_obstacle);

  // Map parameters
  RCLCPP_INFO(this->get_logger(), "Local Map Buffer: %f", par_.map_buffer);
  RCLCPP_INFO(this->get_logger(), "Center Shift Factor: %f", par_.center_shift_factor);
  RCLCPP_INFO(this->get_logger(), "initial_wdx: %f", par_.initial_wdx);
  RCLCPP_INFO(this->get_logger(), "initial_wdy: %f", par_.initial_wdy);
  RCLCPP_INFO(this->get_logger(), "initial_wdz: %f", par_.initial_wdz);
  RCLCPP_INFO(this->get_logger(), "min_wdx: %f", par_.min_wdx);
  RCLCPP_INFO(this->get_logger(), "min_wdy: %f", par_.min_wdy);
  RCLCPP_INFO(this->get_logger(), "min_wdz: %f", par_.min_wdz);
  RCLCPP_INFO(this->get_logger(), "Res: %f", par_.res);

  // Communication delay parameters
  RCLCPP_INFO(this->get_logger(), "Use Comm Delay Inflation: %d", par_.use_comm_delay_inflation);
  RCLCPP_INFO(this->get_logger(), "Comm Delay Inflation Alpha: %f",
              par_.comm_delay_inflation_alpha);
  RCLCPP_INFO(this->get_logger(), "Comm Delay Inflation Max: %f", par_.comm_delay_inflation_max);
  RCLCPP_INFO(this->get_logger(), "Comm Delay Filter Alpha: %f", par_.comm_delay_filter_alpha);

  // Simulation parameters
  RCLCPP_INFO(this->get_logger(), "D435 Depth Max: %f", par_.depth_camera_depth_max);
  RCLCPP_INFO(this->get_logger(), "FOV Visual Depth: %f", par_.fov_visual_depth);
  RCLCPP_INFO(this->get_logger(), "FOV Visual X Deg: %f", par_.fov_visual_x_deg);
  RCLCPP_INFO(this->get_logger(), "FOV Visual Y Deg: %f", par_.fov_visual_y_deg);
  RCLCPP_INFO(this->get_logger(), "Use Sphere Sensing: %d", par_.use_sphere_sensing);
  RCLCPP_INFO(this->get_logger(), "Sphere Sensing Radius: %f", par_.sphere_sensing_radius);

  // Initial guess parameters
  RCLCPP_INFO(this->get_logger(), "Use Multiple Initial Guesses: %d",
              par_.use_multiple_initial_guesses);
  RCLCPP_INFO(this->get_logger(), "Num of Perturbations: %d", par_.num_perturbation_for_ig);
  RCLCPP_INFO(this->get_logger(), "r_max: %f", par_.r_max_for_ig);

  // Optimization parameters
  RCLCPP_INFO(this->get_logger(), "Horizon: %f", par_.horizon);
  RCLCPP_INFO(this->get_logger(), "DC: %f", par_.dc);
  RCLCPP_INFO(this->get_logger(), "V Max: %f", par_.v_max);
  RCLCPP_INFO(this->get_logger(), "A Max: %f", par_.a_max);
  RCLCPP_INFO(this->get_logger(), "J Max: %f", par_.j_max);
  RCLCPP_INFO(this->get_logger(), "Closed Form Verbose: %d", par_.closed_form_traj_verbose);
  RCLCPP_INFO(this->get_logger(), "Control Cost Weight: %f", par_.jerk_weight);
  RCLCPP_INFO(this->get_logger(), "Obstacles and Agents Distance Weight: %f", par_.dynamic_weight);
  RCLCPP_INFO(this->get_logger(), "Time Weight: %f", par_.time_weight);
  RCLCPP_INFO(this->get_logger(), "Position Anchor Weight: %f", par_.pos_anchor_weight);
  RCLCPP_INFO(this->get_logger(), "Static Obstacle Weight: %f", par_.stat_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Bodyrate Constr. Weight: %f",
              par_.dyn_constr_bodyrate_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Tilt Constr. Weight: %f",
              par_.dyn_constr_tilt_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Thrust Constr. Weight: %f",
              par_.dyn_constr_thrust_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Vel Constr. Weight: %f",
              par_.dyn_constr_vel_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Aeccel Constr. Weight: %f",
              par_.dyn_constr_acc_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Jerk Constr. Weight: %f",
              par_.dyn_constr_jerk_weight);
  RCLCPP_INFO(this->get_logger(), "Num Dynamic Obstacles Samples: %d", par_.num_dyn_obst_samples);
  RCLCPP_INFO(this->get_logger(), "Local Traj Co: %f", par_.planner_Co);
  RCLCPP_INFO(this->get_logger(), "Local Traj Cw: %f", par_.planner_Cw);
  RCLCPP_INFO(this->get_logger(), "Verbose Computation Time: %d", verbose_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Drone Bbox: (%f, %f, %f)", par_.drone_bbox[0],
              par_.drone_bbox[1], par_.drone_bbox[2]);
  RCLCPP_INFO(this->get_logger(), "Goal Radius: %f", par_.goal_radius);
  RCLCPP_INFO(this->get_logger(), "Goal Seen Radius: %f", par_.goal_seen_radius);
  RCLCPP_INFO(this->get_logger(), "Init Turn BF: %f", par_.init_turn_bf);
  RCLCPP_INFO(this->get_logger(), "Integral Resolution: %d", par_.integral_resolution);
  RCLCPP_INFO(this->get_logger(), "Hinge Mu: %f", par_.hinge_mu);
  RCLCPP_INFO(this->get_logger(), "Omega Max: %f", par_.omega_max);
  RCLCPP_INFO(this->get_logger(), "Tilt Max Rad: %f", par_.tilt_max_rad);
  RCLCPP_INFO(this->get_logger(), "F Min: %f", par_.f_min);
  RCLCPP_INFO(this->get_logger(), "F Max: %f", par_.f_max);
  RCLCPP_INFO(this->get_logger(), "Mass: %f", par_.mass);
  RCLCPP_INFO(this->get_logger(), "Gravity: %f", par_.g);
  RCLCPP_INFO(this->get_logger(), "Fopt Threshold: %f", par_.fopt_threshold);

  // L-BFGS parameters
  RCLCPP_INFO(this->get_logger(), "f_dec_coeff: %f", par_.f_dec_coeff);
  RCLCPP_INFO(this->get_logger(), "Cautious Factor: %f", par_.cautious_factor);
  RCLCPP_INFO(this->get_logger(), "Past: %d", par_.past);
  RCLCPP_INFO(this->get_logger(), "Max Linesearch: %d", par_.max_linesearch);
  RCLCPP_INFO(this->get_logger(), "Max Iterations: %d", par_.max_iterations);
  RCLCPP_INFO(this->get_logger(), "g_epsilon: %f", par_.g_epsilon);
  RCLCPP_INFO(this->get_logger(), "Delta: %f", par_.delta);

  // Dynamic obstacles parameters
  RCLCPP_INFO(this->get_logger(), "Traj Lifetime: %f", par_.traj_lifetime);

  // Dynamic k_value parameters
  RCLCPP_INFO(this->get_logger(), "Num Replanning Before Adapt: %d",
              par_.num_replanning_before_adapt);
  RCLCPP_INFO(this->get_logger(), "Default K Value End: %d", par_.default_k_value);
  RCLCPP_INFO(this->get_logger(), "Alpha K Value: %f", par_.alpha_k_value_filtering);
  RCLCPP_INFO(this->get_logger(), "K Value Inflation: %f", par_.k_value_factor);

  // Yaw-related parameters
  RCLCPP_INFO(this->get_logger(), "Alpha Filter Dyaw: %f", par_.alpha_filter_dyaw);
  RCLCPP_INFO(this->get_logger(), "W Max: %f", par_.w_max);
  RCLCPP_INFO(this->get_logger(), "Yaw Spinning Threshold: %d", par_.yaw_spinning_threshold);
  RCLCPP_INFO(this->get_logger(), "Yaw Spinning Dyaw: %f", par_.yaw_spinning_dyaw);

  // Simulation env parameters
  RCLCPP_INFO(this->get_logger(), "Force Goal Z: %d", par_.force_goal_z);
  RCLCPP_INFO(this->get_logger(), "Default Goal Z: %f", par_.default_goal_z);

  // Debug flag
  RCLCPP_INFO(this->get_logger(), "Debug Verbose: %d", par_.debug_verbose);
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function to clean up old trajs in DYNUS
 */
void MIGHTY_NODE::cleanUpOldTrajsCallback() {
  // Get current time
  double current_time = this->now().seconds();

  // Clean up old trajs
  mighty_ptr_->cleanUpOldTrajs(current_time);
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function to update the traj
 * @param msg Trajectory message
 */
void MIGHTY_NODE::trajCallback(const dynus_interfaces::msg::DynTraj::SharedPtr msg) {
  // Filter out its own traj. Only drop on id match when the sender is an agent;
  // the obstacle tracker reuses small ids (0,1,2,...) that can collide with our
  // own namespace-derived id_ (e.g. RR04 -> id_=4), which would silently drop
  // dynamic-obstacle predictions.
  if (msg->is_agent && msg->id == id_) return;

  // Get current time
  double current_time = this->now().seconds();

  // Get dynTraj from the message
  auto traj = std::make_shared<dynTraj>();
  convertDynTrajMsg2DynTraj(*msg, traj, current_time);

  // Helper to publish a trajectory as a LINE_STRIP marker
  auto publishTrajMarker = [&](const std::shared_ptr<dynTraj>& t,
                               rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr& pub,
                               float r, float g, float b, double z_offset) {
    visualization_msgs::msg::Marker m;
    m.header.frame_id = par_.map_frame_id;
    m.header.stamp = this->now();
    m.ns = "other_agent_traj_" + std::to_string(msg->id);
    m.id = msg->id;
    m.type = visualization_msgs::msg::Marker::LINE_STRIP;
    m.action = visualization_msgs::msg::Marker::ADD;
    m.scale.x = 0.08;
    m.color.r = r;
    m.color.g = g;
    m.color.b = b;
    m.color.a = 1.0;
    m.lifetime = rclcpp::Duration(1, 0);  // 1 second
    const int N = 30;
    // Sample over the trajectory's actual time range instead of a fixed window
    double t0 = current_time;
    double tf;
    if (t->mode == dynTraj::Mode::Piecewise && !t->pwp.times.empty()) {
      tf = t->pwp.getEndTime();
    } else if (t->mode == dynTraj::Mode::Quintic) {
      tf = t->poly_end_time;
    } else {
      tf = current_time + 2.0;  // fallback for analytic
    }
    if (tf <= t0) tf = t0 + 2.0;  // safety: ensure non-degenerate range
    // Skip if piecewise trajectory has no data (empty times → eval would crash)
    if (t->mode == dynTraj::Mode::Piecewise && t->pwp.times.empty()) return;
    for (int i = 0; i <= N; i++) {
      double ti = t0 + (tf - t0) * i / N;
      Eigen::Vector3d p = t->eval(ti);
      geometry_msgs::msg::Point pt;
      pt.x = p.x();
      pt.y = p.y();
      pt.z = p.z() + z_offset;
      m.points.push_back(pt);
    }
    pub->publish(m);
  };

  // Apply frame alignment transform if enabled
  if (par_.use_frame_alignment) {
    // Drop trajectory if we haven't received frame alignment for this sender yet
    {
      std::lock_guard<std::mutex> lock(frame_align_mutex_);
      if (!frame_align_received_.count(msg->id) || !frame_align_received_[msg->id]) {
        RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 5000,
                             "Dropping traj from agent %d — no frame alignment received yet", msg->id);
        return;
      }
    }

    // Publish BEFORE transform (red, slight z offset for debug)
    publishTrajMarker(traj, pub_traj_received_, 1.0, 0.2, 0.2, 0.1);

    // Frame alignment: other's frame → ego's frame
    applyFrameAlignTransform(traj, msg->id);

    // Undo ego's sim_frame_offset: ego's frame → world frame
    // (In fake_sim the agent operates in world frame, so we need to bring it back)
    if (!sim_frame_offset_inv_.isIdentity(1e-9)) applyTransformToTraj(traj, sim_frame_offset_inv_);
  }

  // Always publish other agents' trajectories in solid blue (after transform if applicable)
  publishTrajMarker(traj, pub_traj_transformed_, 0.1, 0.4, 1.0, 0.0);

  // Pass the dynTraj to mighty.cpp
  mighty_ptr_->addTraj(traj, current_time);
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function for the state of the agent
 * @param msg State message
 */
void MIGHTY_NODE::stateCallback(const dynus_interfaces::msg::State::SharedPtr msg) {
  if (par_.use_state_update) {
    state current_state;
    current_state.setPos(msg->pos.x, msg->pos.y, msg->pos.z);
    current_state.setVel(msg->vel.x, msg->vel.y, msg->vel.z);
    current_state.setAccel(0.0, 0.0, 0.0);
    double roll, pitch, yaw;
    quaternion2Euler(msg->quat, roll, pitch, yaw);
    current_state.setYaw(yaw);
    current_state.t = this->now().seconds();
    mighty_ptr_->updateState(current_state);

    // publish the state
    publishCurrentState(current_state);

    // Publish the velocity in text
    if (par_.visual_level >= 1) publishVelocityInText(current_state.pos, current_state.vel.norm());
  }

  if (!state_initialized_) {
    // If we don't use state update, we need to initialize the state
    if (!par_.use_state_update) {
      state current_state;
      current_state.setPos(msg->pos.x, msg->pos.y, msg->pos.z);
      current_state.setVel(msg->vel.x, msg->vel.y, msg->vel.z);
      current_state.setAccel(0.0, 0.0, 0.0);
      double roll, pitch, yaw;
      quaternion2Euler(msg->quat, roll, pitch, yaw);
      current_state.setYaw(yaw);
      current_state.t = this->now().seconds();
      mighty_ptr_->updateState(current_state);
    }
    RCLCPP_INFO(this->get_logger(), "State initialized");
    state_initialized_ = true;
    timer_goal_->reset();
  }

  // Throttle actual_traj viz to ~20 Hz to match the replan viz throttle.
  // stateCallback runs at 100 Hz (fake_sim publishes state every 10 ms) and
  // the persistent LINE_STRIP marker is large enough that RViz drops
  // messages with "some messages were lost" warnings without throttling.
  if (par_.visual_level >= 1) {
    const double t_now_actual = this->now().seconds();
    if (t_now_actual - last_actual_traj_publish_t_ >= 0.05) {
      publishActualTraj();
      last_actual_traj_publish_t_ = t_now_actual;
    }
  }

  // MinPos: broadcast own position to peers (throttled)
  if (pub_peer_pose_ && par_.expl_peer_publish_rate_hz > 0.0) {
    const double t_now_peer = this->now().seconds();
    const double period = 1.0 / par_.expl_peer_publish_rate_hz;
    if (t_now_peer - last_peer_pose_publish_t_ >= period) {
      geometry_msgs::msg::PoseStamped pose_msg;
      pose_msg.header.stamp = this->now();
      pose_msg.header.frame_id = ns_;
      pose_msg.pose.position.x = msg->pos.x;
      pose_msg.pose.position.y = msg->pos.y;
      pose_msg.pose.position.z = msg->pos.z;
      pub_peer_pose_->publish(pose_msg);
      last_peer_pose_publish_t_ = t_now_peer;
    }
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function for replanning
 */
void MIGHTY_NODE::replanCallback() {
  // Get the current time as double
  double current_time = this->now().seconds();

  // Set computation times to zero
  setComputationTimesToZero();

  // Pass current ESDF snapshot to planner (ground robot only)
  if (par_.use_esdf_cost && esdf_grid_) {
    mighty_ptr_->setEsdfGrid(esdf_grid_);
  }

  // Replan
  auto [replanning_result, hgp_result] =
      mighty_ptr_->replan(replanning_computation_time_, current_time);

  // Republish the terminal goal marker so RViz tracks any in-replan
  // relocation done by sanitizeTerminalGoal (e.g. when an obstacle gets
  // sensed after the original click and the goal jumps to a clear cell).
  if (par_.relocate_occupied_goal) {
    state gterm_now;
    mighty_ptr_->getGterm(gterm_now);
    publishState(gterm_now, pub_point_G_term_);
  }

  // Get computation time (used to find point A) - note this value is not updated in the replan
  // function
  if (replanning_result) {
    // Get the replanning computation time
    replanning_computation_time_ = this->now().seconds() - current_time;
    if (par_.debug_verbose)
      printf("Total Replanning: %f ms\n", replanning_computation_time_ * 1000.0);
  }

  // Debug log: dump global path + local trajectory endpoints on every replan
  // so backward-planning bugs can be diagnosed offline.
  if (auto& s = debug_log_stream(); s.is_open()) {
    state cur, gterm, gsub, lA;
    mighty_ptr_->getState(cur);
    mighty_ptr_->getGterm(gterm);
    mighty_ptr_->getG(gsub);
    mighty_ptr_->getA(lA);
    s << debug_log_t(this->now().seconds()) << " REPLAN"
      << " ok=" << (replanning_result ? 1 : 0)
      << " hgp=" << (hgp_result ? 1 : 0)
      << " state=" << cur.pos.x() << "," << cur.pos.y() << "," << cur.pos.z()
      << " A=" << lA.pos.x() << "," << lA.pos.y() << "," << lA.pos.z()
      << " G=" << gsub.pos.x() << "," << gsub.pos.y() << "," << gsub.pos.z()
      << " Gterm=" << gterm.pos.x() << "," << gterm.pos.y() << "," << gterm.pos.z()
      << "\n";
    if (hgp_result) {
      vec_Vecf<3> gp;
      mighty_ptr_->getGlobalPath(gp);
      s << "  HGP_PATH n=" << gp.size();
      auto dump = [&](size_t i) {
        s << " [" << i << "]=" << gp[i].x() << "," << gp[i].y() << "," << gp[i].z();
      };
      if (gp.size() >= 1) dump(0);
      if (gp.size() >= 2) dump(1);
      if (gp.size() >= 3) dump(2);
      if (gp.size() >= 4) dump(gp.size() - 1);
      s << "\n";
    }
    s.flush();
  }

  // Frontier-unreachable detection: if the global planner consistently fails
  // on an exploration goal, mark it INVALIDATED so the manager can pick the
  // next one. We use hgp_result rather than replanning_result because the
  // local L-BFGS may legitimately fail intermittently while HGP is fine.
  if (par_.expl_enabled && exploration_active_ && frontier_manager_) {
    if (!hgp_result) {
      if (++unreachable_consec_count_ >= par_.expl_unreachable_consec_thresh) {
        RCLCPP_WARN(this->get_logger(),
                    "Exploration: frontier %lu unreachable after %d HGP failures, "
                    "invalidating",
                    static_cast<unsigned long>(current_explore_id_),
                    unreachable_consec_count_);
        frontier_manager_->markInvalidated(current_explore_id_, this->now().seconds());
        exploration_active_ = false;
        unreachable_consec_count_ = 0;
      }
    } else {
      unreachable_consec_count_ = 0;
    }
  }

  // To share trajectory with other agents
  if (replanning_result && par_.share_traj) publishOwnTraj();

  // Publish trajectory for tracking (increments trajectory_id on replan)
  if (replanning_result) {
    publishTrajectory();
  }

  // Publish waypoint path for the MPC controller (ground robot only)
  if (replanning_result && par_.vehicle_type == "ground_robot") {
    publishMpcPath();
  }

  // Publish command-to-execution time (time from goal received to first trajectory)
  if (replanning_result && waiting_for_first_traj_) {
    double command_to_exec_time_ms = (this->now() - goal_received_time_).seconds() * 1000.0;
    std_msgs::msg::Float64 time_msg;
    time_msg.data = command_to_exec_time_ms;
    pub_command_to_exec_time_->publish(time_msg);
    waiting_for_first_traj_ = false;
    RCLCPP_INFO(this->get_logger(), "Command to execution time: %.2f ms", command_to_exec_time_ms);
  }

  // Throttle the visualization block to ~20 Hz. The replan loop is 100 Hz
  // which is fine for control, but publishing all of these MarkerArrays at
  // that rate (especially traj_committed_colored, hgp_path_marker,
  // original_hgp_path_marker — measured at 100-200 Hz) overwhelms RViz and it
  // drops messages with "some messages were lost" warnings. 20 Hz is visually
  // smooth and cuts marker traffic ~5x.
  bool do_viz = false;
  if (par_.visual_level >= 1) {
    const double t_now_viz = this->now().seconds();
    if (t_now_viz - last_replan_viz_publish_t_ >= 0.05) {
      do_viz = true;
      last_replan_viz_publish_t_ = t_now_viz;
    }
  }

  // For visualization of global path
  if (hgp_result && do_viz) publishGlobalPath();

  // For visualization of free global path
  if (hgp_result && do_viz) publishFreeGlobalPath();

  // For visualization of local_global_path and local_global_path_after_push_
  if (hgp_result && do_viz) publishLocalGlobalPath();

  // For visualization of the local trajectory
  if (replanning_result && do_viz) publishTraj();

  // For visualization of the safe corridor
  if (hgp_result && do_viz) publishPoly();

  // For visualization of point G and point A
  if (replanning_result && do_viz) {
    publishPointG();
    publishPointE();
    publishPointA();
  }

  // For visualization of control points
  if (replanning_result && do_viz) publishCps();

  // For visualization of static push points and P points
  if (replanning_result && do_viz) {
    mighty_ptr_->getStaticPushPoints(static_push_points_);
    publishStaticPushPoints();
  }

  // If verbose_computation_time_ or use_benchmark_ is true, we need to retrieve data from
  // mighty_ptr_
  if (verbose_computation_time_ || use_benchmark_) retrieveData();

  // Verbose computation time to the terminal
  if (verbose_computation_time_) printComputationTime(replanning_result);

  // Record the data
  // if (replanning_result && use_benchmark_)
  if (use_benchmark_) recordData(replanning_result);

  // Usually this is done is goal callback but becuase we don't call that in push path test, we need
  // to call it here
  if (par_.use_path_push_for_visualization) publishFOV();
}

// ----------------------------------------------------------------------------

/**
 * @brief Public callback for the terminal goal topic. A user-issued goal
 *        preempts any active frontier exploration.
 */
void MIGHTY_NODE::terminalGoalCallback(const geometry_msgs::msg::PoseStamped& msg) {
  terminalGoalCallbackImpl(msg, /*from_user=*/true);
}

/**
 * @brief Apply the per-agent formation offset to a shared /swarm_goal
 *        publication and forward to the standard goal-pin path.
 *
 *        Routing through terminalGoalCallbackImpl with from_user=true means
 *        the swarm goal preempts exploration just like a manual goal would,
 *        and the existing reconstruct() pin (P[M_]=xf_) automatically becomes
 *        G_swarm + δ_i for each agent.
 */
void MIGHTY_NODE::swarmGoalCallback(const geometry_msgs::msg::PoseStamped& msg) {
  geometry_msgs::msg::PoseStamped offset_msg = msg;
  offset_msg.pose.position.x += par_.formation_self_offset[0];
  offset_msg.pose.position.y += par_.formation_self_offset[1];
  offset_msg.pose.position.z += par_.formation_self_offset[2];
  terminalGoalCallbackImpl(offset_msg, /*from_user=*/true);
}

/**
 * @brief Internal goal-issuing entry point shared by the public callback and
 *        the frontier exploration loop. `from_user=true` records that a manual
 *        goal is now active (preempting exploration); `from_user=false` is
 *        used by the explore-select callback so the same routine doesn't
 *        clobber its own state.
 */
void MIGHTY_NODE::terminalGoalCallbackImpl(const geometry_msgs::msg::PoseStamped& msg,
                                           bool from_user) {
  // Record the time when goal is received (for command-to-execution timing)
  goal_received_time_ = this->now();
  waiting_for_first_traj_ = true;

  // Debug log: every received term_goal
  if (auto& s = debug_log_stream(); s.is_open()) {
    s << debug_log_t(this->now().seconds()) << " TERM_GOAL "
      << msg.pose.position.x << " " << msg.pose.position.y << " " << msg.pose.position.z
      << " from_user=" << (from_user ? 1 : 0)
      << " frame=" << msg.header.frame_id
      << "\n";
    s.flush();
  }

  if (from_user) {
    manual_goal_active_ = true;
    // A manual goal preempts any in-progress exploration goal.
    exploration_active_ = false;
    exploration_start_captured_ = false;
    unreachable_consec_count_ = 0;
    // Operator override of a return-home — start a fresh session.
    home_return_requested_ = false;
  }

  // Set the terminal goal
  state G_term;
  double goal_z;

  // If force_goal_z is true, set the goal_z to default_goal_z
  if (par_.force_goal_z)
    goal_z = par_.default_goal_z;
  else
    goal_z = msg.pose.position.z;

  // Check if the goal_z is within the limits
  if (goal_z < par_.z_min || goal_z > par_.z_max) {
    RCLCPP_ERROR(this->get_logger(), "Goal z is out of bounds: %f", goal_z);
    return;
  }

  // Set the terminal goal
  G_term.setPos(msg.pose.position.x, msg.pose.position.y, goal_z);

  // If the goal lies in an occupied voxel, relocate it to the nearest
  // free/unknown cell with ||drone_bbox|| clearance from the original point.
  // The mutated G_term then flows into both the planner and the RViz marker.
  if (par_.relocate_occupied_goal) {
    if (!mighty_ptr_->sanitizeTerminalGoal(G_term)) {
      // ERROR (not WARN) so it survives --log-level error in the launch file.
      RCLCPP_ERROR(this->get_logger(),
                   "Goal at (%.2f,%.2f,%.2f) is in occupied space and could not be "
                   "relocated; ignoring.",
                   msg.pose.position.x, msg.pose.position.y, goal_z);
      return;
    }
    // Re-clamp z in case the BFS pushed the relocated goal out of bounds.
    if (G_term.pos.z() < par_.z_min) G_term.pos.z() = par_.z_min;
    if (G_term.pos.z() > par_.z_max) G_term.pos.z() = par_.z_max;
  }

  // Update the terminal goal
  mighty_ptr_->setTerminalGoal(G_term);

  // Publish the term goal for visualization
  publishState(G_term, pub_point_G_term_);

  // Start replanning
  timer_replanning_->reset();
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish velocity in text
 * @param position Position to publish the text
 * @param velocity Velocity to publish
 */
void MIGHTY_NODE::publishVelocityInText(const Eigen::Vector3d& position, double velocity) {
  // Set velocity's precision to 2 decimal points
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << velocity;

  // Make a string
  std::string text = oss.str() + "m/s";

  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = par_.map_frame_id;
  marker.header.stamp = this->get_clock()->now();
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.ns = "velocity";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
  marker.scale.z = 1.0;
  marker.color.a = 1.0;
  marker.color.r = 0.0;
  marker.color.g = 0.0;
  marker.color.b = 0.0;
  marker.text = text;
  marker.pose.position.x = position.x();
  marker.pose.position.y = position.y();
  marker.pose.position.z = position.z() + 5.0;
  marker.pose.orientation.w = 1.0;
  pub_vel_text_->publish(marker);
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function to check if the goal is reached
 */
void MIGHTY_NODE::goalReachedCheckCallback() {
  if (!mighty_ptr_->goalReachedCheck()) return;

  if (use_benchmark_) {
    logData();
  }
  pub_goal_reached_->publish(std_msgs::msg::Empty());

  // If the robot reached an exploration goal, mark it visited and clear our
  // active flag so the next explore-select tick can pick the next frontier.
  if (par_.expl_enabled && exploration_active_ && frontier_manager_) {
    frontier_manager_->markVisited(current_explore_id_);
    exploration_active_ = false;
    unreachable_consec_count_ = 0;
  }
  // A manual goal that just completed releases the preemption.
  if (manual_goal_active_) {
    manual_goal_active_ = false;
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function to get the initial pose from tf (for hardware use case)
 */
void MIGHTY_NODE::getInitialPoseHwCallback() {
  // First find the transformation matrix from map to camera
  try {
    init_pose_transform_stamped_ =
        tf2_buffer_->lookupTransform("map", initial_pose_topic_, tf2::TimePointZero);

    // Print out the initial pose
    RCLCPP_INFO(this->get_logger(), "Initial pose received: (%f, %f, %f)",
                init_pose_transform_stamped_.transform.translation.x,
                init_pose_transform_stamped_.transform.translation.y,
                init_pose_transform_stamped_.transform.translation.z);

    // Push the initial pose to dynus
    mighty_ptr_->setInitialPose(init_pose_transform_stamped_);
  } catch (tf2::TransformException& ex) {
    RCLCPP_WARN(this->get_logger(), "Transform error: %s", ex.what());
    return;
  }

  // flag
  if (!initial_pose_received_) {
    initial_pose_received_ = true;
    timer_initial_pose_->cancel();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback for frame alignment transforms from lidar registration
 * @param msg TransformStamped message (T^{ego/map}_{other/map})
 * @param agent_id The ID of the other agent
 */
void MIGHTY_NODE::frameAlignCallback(const geometry_msgs::msg::TransformStamped::SharedPtr msg,
                                     int agent_id) {
  Eigen::Matrix4d T = mighty_utils::transformStampedToMatrix(*msg);
  std::lock_guard<std::mutex> lock(frame_align_mutex_);
  frame_align_transforms_[agent_id] = T;
  if (!frame_align_received_[agent_id]) {
    frame_align_received_[agent_id] = true;
    RCLCPP_INFO(this->get_logger(), "Received first frame alignment for agent %d", agent_id);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Apply frame alignment transform to a trajectory from another agent
 * @param traj The trajectory to transform (modified in place)
 * @param sender_id The ID of the sending agent
 */
void MIGHTY_NODE::applyFrameAlignTransform(std::shared_ptr<dynTraj>& traj, int sender_id) {
  Eigen::Matrix4d T;
  {
    std::lock_guard<std::mutex> lock(frame_align_mutex_);
    auto it = frame_align_transforms_.find(sender_id);
    if (it == frame_align_transforms_.end() || !frame_align_received_[sender_id]) {
      RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 5000,
                           "No frame alignment for agent %d yet", sender_id);
      return;
    }
    T = it->second;
  }
  applyTransformToTraj(traj, T);
}

// ----------------------------------------------------------------------------

/**
 * @brief Apply a 4x4 transform to a trajectory (modifies in place)
 * @param traj The trajectory to transform
 * @param T The 4x4 transformation matrix
 */
void MIGHTY_NODE::applyTransformToTraj(std::shared_ptr<dynTraj>& traj, const Eigen::Matrix4d& T) {
  Eigen::Matrix3d R = T.block<3, 3>(0, 0);

  // Transform PieceWiseQuinticPol
  if (traj->mode == dynTraj::Mode::Piecewise) {
    for (size_t i = 0; i < traj->pwp.coeff_x.size(); i++) {
      for (int j = 0; j < 5; j++)  // non-constant: rotate only
      {
        Eigen::Vector3d c(traj->pwp.coeff_x[i](j), traj->pwp.coeff_y[i](j),
                          traj->pwp.coeff_z[i](j));
        c = R * c;
        traj->pwp.coeff_x[i](j) = c.x();
        traj->pwp.coeff_y[i](j) = c.y();
        traj->pwp.coeff_z[i](j) = c.z();
      }
      // constant term (j=5): full transform
      Eigen::Vector4d h;
      h << traj->pwp.coeff_x[i](5), traj->pwp.coeff_y[i](5), traj->pwp.coeff_z[i](5), 1.0;
      h = T * h;
      traj->pwp.coeff_x[i](5) = h(0);
      traj->pwp.coeff_y[i](5) = h(1);
      traj->pwp.coeff_z[i](5) = h(2);
    }
  } else if (traj->mode == dynTraj::Mode::Quintic) {
    for (int j = 0; j < 5; j++) {
      Eigen::Vector3d c(traj->cx(j), traj->cy(j), traj->cz(j));
      c = R * c;
      traj->cx(j) = c.x();
      traj->cy(j) = c.y();
      traj->cz(j) = c.z();
    }
    Eigen::Vector4d h;
    h << traj->cx(5), traj->cy(5), traj->cz(5), 1.0;
    h = T * h;
    traj->cx(5) = h(0);
    traj->cy(5) = h(1);
    traj->cz(5) = h(2);
  }
  // Analytic mode: string expressions can't be transformed — agents use pwp/quintic

  // Transform goal
  if (traj->is_agent) {
    Eigen::Vector4d g;
    g << traj->goal.x(), traj->goal.y(), traj->goal.z(), 1.0;
    g = T * g;
    traj->goal = g.head<3>();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Convert dynTraj message to dynTraj
 * @param msg dynTraj message
 * @param traj dynTraj
 * @param current_time current time
 */
void MIGHTY_NODE::convertDynTrajMsg2DynTraj(const dynus_interfaces::msg::DynTraj& msg,
                                            std::shared_ptr<dynTraj>& traj, double current_time) {
  // Inflate bbox using drone_bbox
  // We need to use the obstacle's bbox as well as ego drone's bbox
  traj->bbox << msg.bbox[0] / 2.0 + par_.drone_bbox[0] / 2.0,
      msg.bbox[1] / 2.0 + par_.drone_bbox[1] / 2.0, msg.bbox[2] / 2.0 + par_.drone_bbox[2] / 2.0;

  // Get id
  traj->id = msg.id;

  // Get pwp — quintic (6-coeff) takes priority, fall back to cubic (4-coeff) promoted to quintic
  if (msg.mode == "pwp") {
    if (!msg.quintic_pwp.times.empty()) {
      traj->pwp = mighty_utils::convertPwpMsg2Pwp(msg.quintic_pwp);
    } else if (!msg.pwp.times.empty()) {
      // Cubic PWP from obstacle tracker — promote to quintic by zero-padding high-order coeffs
      auto cubic = mighty_utils::convertPwpMsg2Pwp(msg.pwp);
      traj->pwp.times = cubic.times;
      traj->pwp.coeff_x.resize(cubic.coeff_x.size());
      traj->pwp.coeff_y.resize(cubic.coeff_y.size());
      traj->pwp.coeff_z.resize(cubic.coeff_z.size());
      for (size_t i = 0; i < cubic.coeff_x.size(); ++i) {
        // Quintic: [a*u^5, b*u^4, c*u^3, d*u^2, e*u, f] — set a=0, b=0, then c,d,e,f from cubic
        traj->pwp.coeff_x[i] << 0.0, 0.0, cubic.coeff_x[i](0), cubic.coeff_x[i](1),
            cubic.coeff_x[i](2), cubic.coeff_x[i](3);
        traj->pwp.coeff_y[i] << 0.0, 0.0, cubic.coeff_y[i](0), cubic.coeff_y[i](1),
            cubic.coeff_y[i](2), cubic.coeff_y[i](3);
        traj->pwp.coeff_z[i] << 0.0, 0.0, cubic.coeff_z[i](0), cubic.coeff_z[i](1),
            cubic.coeff_z[i](2), cubic.coeff_z[i](3);
      }
    }
    traj->mode = dynTraj::Mode::Piecewise;
  }

  // Find quihtic coefficients from the given pwp
  if (msg.mode == "quintic") {
    traj->cx = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_x);
    traj->cy = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_y);
    traj->cz = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_z);
    traj->poly_start_time = msg.poly_start_time;
    traj->poly_end_time = msg.poly_end_time;
    traj->mode = dynTraj::Mode::Quintic;
  }

  if (msg.mode == "analytic") {
    traj->mode = dynTraj::Mode::Analytic;
  }

  // Get covariances
  if (!msg.is_agent) {
    if (msg.ekf_cov_p.size() != 0) {
      traj->ekf_cov_p = mighty_utils::convertCovMsg2Cov(msg.ekf_cov_p);  // ekf cov
    }

    if (msg.ekf_cov_q.size() != 0) {
      traj->ekf_cov_q = mighty_utils::convertCovMsg2Cov(msg.ekf_cov_q);  // ekf cov
    }

    if (msg.poly_cov.size() != 0) {
      traj->poly_cov = mighty_utils::convertCovMsg2Cov(msg.poly_cov);  // future traj cov
    }

    if (msg.function.size() == 3) {
      traj->traj_x = msg.function[0];
      traj->traj_y = msg.function[1];
      traj->traj_z = msg.function[2];
    }

    if (msg.velocity.size() == 3) {
      traj->traj_vx = msg.velocity[0];
      traj->traj_vy = msg.velocity[1];
      traj->traj_vz = msg.velocity[2];
    }

    if (msg.function.size() == 3 && msg.velocity.size() == 3) {
      if (traj->compileAnalytic()) {
        // Change the mode only when we successfully compiled the analytic trajectory
        traj->mode = dynTraj::Mode::Analytic;
        // printf("Successfully compiled analytic traj id=%d\n", traj->id);
      } else {
        RCLCPP_ERROR(this->get_logger(),
                     "Failed to compile analytic traj id=%d, falling back to zeros.", traj->id);
        // leave mode as whatever it was (Piecewise/Quintic),
        // or explicitly set a safe default here
      }
    }
  }

  // Record received time
  traj->time_received = current_time;

  // Get is_agent
  traj->is_agent = msg.is_agent;

  // Get terminal goal
  if (traj->is_agent) traj->goal << msg.goal[0], msg.goal[1], msg.goal[2];

  // Get communication delay.
  //
  // NOTE: this value is currently WRITTEN BUT NEVER READ — nothing in the
  // planner consumes traj->communication_delay, so comm-delay inflation is not
  // actually applied despite use_comm_delay_inflation /
  // comm_delay_inflation_alpha / comm_delay_inflation_max existing as
  // parameters. Fixing the arithmetic anyway so the field is meaningful if it is
  // ever wired up.
  //
  // The previous form ended the statement after now().seconds(), leaving
  // "-msg.header.stamp.sec + ..." as a discarded expression, so the field got
  // the full epoch (~1.79e9 s) instead of a delay of a few milliseconds.
  if (traj->is_agent && par_.use_comm_delay_inflation) {
    const double msg_stamp = static_cast<double>(msg.header.stamp.sec)
                           + static_cast<double>(msg.header.stamp.nanosec) * 1e-9;
    traj->communication_delay = this->now().seconds() - msg_stamp;

    // Sanity check - if the delay is negative, set it to 0 - send warning message: it's probably
    // due to the clock synchronization issue
    if (traj->communication_delay < 0) {
      traj->communication_delay = 0;
      RCLCPP_WARN(this->get_logger(), "Communication delay is negative. Setting it to 0.");
    }
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish control points
 */

void MIGHTY_NODE::publishCps() {
  // Retrieve control points
  mighty_ptr_->retrieveCPs(cps_);

  // Create a marker array
  visualization_msgs::msg::MarkerArray marker_array;
  marker_array.markers.resize(cps_.size());

  // Loop through the control points (std::vector<Eigen::Matrix<double, 3, 4>>)
  for (int seg = 0; seg < cps_.size(); seg++) {
    // Create a marker
    visualization_msgs::msg::Marker marker;
    marker.header.frame_id = par_.map_frame_id;
    marker.header.stamp = this->now();
    marker.ns = "cp";
    marker.id = seg;
    marker.type = visualization_msgs::msg::Marker::SPHERE_LIST;
    marker.action = visualization_msgs::msg::Marker::ADD;
    marker.pose.orientation.w = 1.0;

    marker.scale.x = 0.5;
    marker.scale.y = 0.5;
    marker.scale.z = 0.5;
    marker.color.a = 1.0;

    // Set different colors for different segments
    if (seg % 3 == 0) {
      marker.color.r = 1.0;
      marker.color.g = 0.0;
      marker.color.b = 0.0;
    } else if (seg % 3 == 1) {
      marker.color.r = 0.0;
      marker.color.g = 1.0;
      marker.color.b = 0.0;
    } else {
      marker.color.r = 0.0;
      marker.color.g = 0.0;
      marker.color.b = 1.0;
    }

    auto cp = cps_[seg];

    // Loop through the control points for each segment
    for (int i = 0; i < cp.cols(); i++) {
      geometry_msgs::msg::Point point;
      point.x = cp(0, i);
      point.y = cp(1, i);
      point.z = cp(2, i);
      marker.points.push_back(point);
    }

    // Add the marker to the marker array
    marker_array.markers[seg] = marker;
  }

  // Publish
  pub_cp_->publish(marker_array);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish static push points
 */
void MIGHTY_NODE::publishStaticPushPoints() {
  // Create a marker array
  visualization_msgs::msg::MarkerArray marker_array;
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = par_.map_frame_id;
  marker.header.stamp = this->now();
  marker.ns = "static_push_points";
  marker.id = static_push_points_id_;
  marker.type = visualization_msgs::msg::Marker::SPHERE_LIST;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 3.0;
  marker.scale.y = 3.0;
  marker.scale.z = 3.0;
  marker.color.a = 1.0;
  marker.color.r = 0.0;
  marker.color.g = 0.925;
  marker.color.b = 1.0;

  // Loop through the static push points
  for (int idx = 0; idx < static_push_points_.size(); idx++) {
    geometry_msgs::msg::Point point;
    point.x = static_push_points_[idx](0);
    point.y = static_push_points_[idx](1);
    point.z = static_push_points_[idx](2);
    marker.points.push_back(point);
  }

  // Add the single marker to the marker array
  marker_array.markers.push_back(marker);

  // Publish
  pub_static_push_points_->publish(marker_array);
}

// ----------------------------------------------------------------------------

/**
 * @brief Set computation times to zero
 */
void MIGHTY_NODE::setComputationTimesToZero() {
  final_g_ = 0.0;
  global_planning_time_ = 0.0;
  hgp_static_jps_time_ = 0.0;
  hgp_check_path_time_ = 0.0;
  hgp_dynamic_astar_time_ = 0.0;
  hgp_recover_path_time_ = 0.0;
  cvx_decomp_time_ = 0.0;
  initial_guess_computation_time_ = 0.0;
  local_traj_computation_time_ = 0.0;
  safe_paths_time_ = 0.0;
  safety_check_time_ = 0.0;
  yaw_sequence_time_ = 0.0;
  yaw_fitting_time_ = 0.0;
}

// ----------------------------------------------------------------------------

/**
 * @brief Retrive computation times from mighty_ptr_
 */
void MIGHTY_NODE::retrieveData() {
  mighty_ptr_->retrieveData(final_g_, global_planning_time_, hgp_static_jps_time_,
                            hgp_check_path_time_, hgp_dynamic_astar_time_, hgp_recover_path_time_,
                            cvx_decomp_time_, initial_guess_computation_time_,
                            local_traj_computation_time_, safety_check_time_, safe_paths_time_,
                            yaw_sequence_time_, yaw_fitting_time_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Print the computation times
 */
void MIGHTY_NODE::printComputationTime(bool result) {
  // Print the computation times
  RCLCPP_INFO(this->get_logger(), "Planner: %s", par_.global_planner.c_str());
  RCLCPP_INFO(this->get_logger(), "Result: %d", result);
  RCLCPP_INFO(this->get_logger(), "Cost (final node's g): %f", final_g_);
  RCLCPP_INFO(this->get_logger(), "Total replanning time [ms]: %f",
              replanning_computation_time_ * 1000.0);
  RCLCPP_INFO(this->get_logger(), "Global Planning Time [ms]: %f", global_planning_time_);
  RCLCPP_INFO(this->get_logger(), "CVX Decomposition Time [ms]: %f", cvx_decomp_time_);
  RCLCPP_INFO(this->get_logger(), "Initial Guess Time [ms]: %f", initial_guess_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Local Traj Time [ms]: %f", local_traj_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Safe Paths Time [ms]: %f", safe_paths_time_);
  RCLCPP_INFO(this->get_logger(), "Safety Check Time [ms]: %f", safety_check_time_);
  RCLCPP_INFO(this->get_logger(), "Yaw Sequence Time [ms]: %f", yaw_sequence_time_);
  RCLCPP_INFO(this->get_logger(), "Yaw Fitting Time [ms]: %f", yaw_fitting_time_);
  if (par_.global_planner == "hgp") {
    RCLCPP_INFO(this->get_logger(), "Static JPS Time [ms]: %f", hgp_static_jps_time_);
    RCLCPP_INFO(this->get_logger(), "Check Path Time [ms]: %f", hgp_check_path_time_);
    RCLCPP_INFO(this->get_logger(), "Dynamic A* Time [ms]: %f", hgp_dynamic_astar_time_);
    RCLCPP_INFO(this->get_logger(), "Recover Path Time [ms]: %f", hgp_recover_path_time_);
  }
  RCLCPP_INFO(this->get_logger(), "------------------------");
}

// ----------------------------------------------------------------------------

/**
 * @brief Record the data
 * @param result result of the replanning
 */
void MIGHTY_NODE::recordData(bool result) {
  // Record all the data into global_path_benchmark_
  std::tuple<bool, double, double, double, double, double, double, double, double, double, double,
             double, double, double, double>
      data;
  if (par_.global_planner == "hgp") {
    data = std::make_tuple(result, final_g_, replanning_computation_time_, global_planning_time_,
                           cvx_decomp_time_, initial_guess_computation_time_,
                           local_traj_computation_time_, safe_paths_time_, safety_check_time_,
                           yaw_sequence_time_, yaw_fitting_time_, hgp_static_jps_time_,
                           hgp_check_path_time_, hgp_dynamic_astar_time_, hgp_recover_path_time_);
  } else {
    data = std::make_tuple(result, final_g_, replanning_computation_time_, global_planning_time_,
                           cvx_decomp_time_, initial_guess_computation_time_,
                           local_traj_computation_time_, safe_paths_time_, safety_check_time_,
                           yaw_sequence_time_, yaw_fitting_time_, 0.0, 0.0, 0.0, 0.0);
  }

  global_path_benchmark_.push_back(data);
}

// ----------------------------------------------------------------------------

/**
 * @brief Log the data to a csv file
 */
void MIGHTY_NODE::logData() {
  // Loc the computation times to csv file
  // std::ofstream log_file(file_path_, std::ios_base::app); // Open the file in append mode
  std::ofstream log_file(file_path_);  // Open the file in overwrite mode
  if (log_file.is_open()) {
    if (par_.global_planner == "hgp") {
      // Header
      log_file << "Planner,Result,Cost (final node's g),Total replanning time [ms],Global Planning "
                  "Time [ms],CVX Decomposition Time [ms],Initial Guess Time [ms],Local Traj Time "
                  "[ms],Safe Paths Time [ms],Safety Check Time [ms],Yaw Sequence Time [ms],Yaw "
                  "Fitting Time [ms],Static JPS Time [ms],Check Path Time [ms],Dynamic A* Time "
                  "[ms],Recover Path Time [ms]\n";
    } else {
      // Header
      log_file << "Planner,Result,Cost (final node's g),Total replanning time [ms],Global Planning "
                  "Time [ms],CVX Decomposition Time [ms],Initial Guess Time [ms],Local Traj Time "
                  "[ms],Safe Paths Time [ms],Safety Check Time [ms],Yaw Sequence Time [ms],Yaw "
                  "Fitting Time [ms]\n";
    }

    // Data
    for (const auto& row : global_path_benchmark_) {
      if (par_.global_planner == "hgp") {
        log_file << par_.global_planner << "," << std::get<0>(row) << "," << std::get<1>(row) << ","
                 << std::get<2>(row) * 1000.0 << "," << std::get<3>(row) << "," << std::get<4>(row)
                 << "," << std::get<5>(row) << "," << std::get<6>(row) << "," << std::get<7>(row)
                 << "," << std::get<8>(row) << "," << std::get<9>(row) << "," << std::get<10>(row)
                 << "," << std::get<11>(row) << "," << std::get<12>(row) << "," << std::get<13>(row)
                 << std::get<14>(row) << "\n";
      } else {
        log_file << par_.global_planner << "," << std::get<0>(row) << "," << std::get<1>(row) << ","
                 << std::get<2>(row) * 1000.0 << "," << std::get<3>(row) << "," << std::get<4>(row)
                 << "," << std::get<5>(row) << "," << std::get<6>(row) << "," << std::get<7>(row)
                 << "," << std::get<8>(row) << "," << std::get<9>(row) << "," << std::get<10>(row)
                 << "\n";
      }
    }

    log_file.close();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the Point G (sub goal)
 */
void MIGHTY_NODE::publishPointG() const {
  // get projected goal (G)
  state G;
  mighty_ptr_->getG(G);

  // Publish the goal for visualization
  publishState(G, pub_point_G_);

  // Corridor-hop TIP target arrow (ground robot only).
  // Draws a yellow arrow at G pointing in direction G.yaw, which is the
  // heading the robot will turn-in-place to once it arrives at this
  // subgoal. Skips visualization entirely for UAVs or when the feature is
  // disabled so no markers clutter the scene.
  if (par_.vehicle_type == "ground_robot" && par_.corridor_hop_enabled) {
    visualization_msgs::msg::Marker arrow;
    arrow.header.frame_id = par_.map_frame_id;
    arrow.header.stamp = this->now();
    arrow.ns = "corridor_yaw_target";
    arrow.id = 0;
    arrow.type = visualization_msgs::msg::Marker::ARROW;
    arrow.action = visualization_msgs::msg::Marker::ADD;
    arrow.scale.x = 0.05;  // shaft diameter
    arrow.scale.y = 0.10;  // head diameter
    arrow.scale.z = 0.10;  // head length
    arrow.color.r = 1.0;
    arrow.color.g = 1.0;
    arrow.color.b = 0.0;
    arrow.color.a = 1.0;
    arrow.lifetime = rclcpp::Duration::from_seconds(1.0);
    geometry_msgs::msg::Point p_start, p_end;
    p_start.x = G.pos.x();
    p_start.y = G.pos.y();
    p_start.z = G.pos.z();
    const double arrow_len = 0.6;
    p_end.x = G.pos.x() + arrow_len * std::cos(G.yaw);
    p_end.y = G.pos.y() + arrow_len * std::sin(G.yaw);
    p_end.z = G.pos.z();
    arrow.points.push_back(p_start);
    arrow.points.push_back(p_end);
    pub_corridor_yaw_target_->publish(arrow);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the Point E (sub goal)
 */
void MIGHTY_NODE::publishPointE() const {
  // get projected goal (E)
  state E;
  mighty_ptr_->getE(E);

  // Publish the goal for visualization
  publishState(E, pub_point_E_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the Point A (trajectory start point)
 */
void MIGHTY_NODE::publishPointA() const {
  // get projected goal (A)
  state A;
  mighty_ptr_->getA(A);

  // Publish the goal for visualization
  publishState(A, pub_point_A_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the current state
 */
void MIGHTY_NODE::publishCurrentState(const state& state) const {
  // Publish the goal for visualization
  publishState(state, pub_current_state_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish state
 */
void MIGHTY_NODE::publishState(
    const state& data,
    const rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr& publisher) const {
  geometry_msgs::msg::PointStamped p;
  p.header.frame_id = par_.map_frame_id;
  p.header.stamp = this->now();
  p.point = eigen2point(data.pos);
  publisher->publish(p);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish its own trajectory for deconfliction
 */
void MIGHTY_NODE::publishOwnTraj() {
  // Get the piecewise quintic polynomial trajectory to share
  mighty_ptr_->getPiecewiseQuinticPol(pwp_to_share_);

  // Apply simulated frame offset (for testing frame alignment in fake_sim).
  // This pre-rotates the trajectory into the agent's own (simulated) frame,
  // so that the receiver's frame alignment can correct it back.
  if (par_.use_frame_alignment && !sim_frame_offset_.isIdentity(1e-9)) {
    Eigen::Matrix3d R = sim_frame_offset_.block<3, 3>(0, 0);
    for (size_t i = 0; i < pwp_to_share_.coeff_x.size(); i++) {
      for (int j = 0; j < 5; j++)  // non-constant: rotate only
      {
        Eigen::Vector3d c(pwp_to_share_.coeff_x[i](j), pwp_to_share_.coeff_y[i](j),
                          pwp_to_share_.coeff_z[i](j));
        c = R * c;
        pwp_to_share_.coeff_x[i](j) = c.x();
        pwp_to_share_.coeff_y[i](j) = c.y();
        pwp_to_share_.coeff_z[i](j) = c.z();
      }
      // constant term (j=5): full transform (rotation + translation)
      Eigen::Vector4d h;
      h << pwp_to_share_.coeff_x[i](5), pwp_to_share_.coeff_y[i](5), pwp_to_share_.coeff_z[i](5),
          1.0;
      h = sim_frame_offset_ * h;
      pwp_to_share_.coeff_x[i](5) = h(0);
      pwp_to_share_.coeff_y[i](5) = h(1);
      pwp_to_share_.coeff_z[i](5) = h(2);
    }
  }

  // Create the message
  dynus_interfaces::msg::DynTraj msg;
  msg.header.stamp = this->now();
  msg.header.frame_id = par_.map_frame_id;
  msg.bbox.push_back(par_.drone_bbox[0]);
  msg.bbox.push_back(par_.drone_bbox[1]);
  msg.bbox.push_back(par_.drone_bbox[2]);
  msg.id = id_;
  msg.mode = "pwp";
  msg.quintic_pwp = mighty_utils::convertPwp2PwpMsg(pwp_to_share_);
  msg.is_agent = true;

  // Get the terminal goal
  state G;
  mighty_ptr_->getG(G);

  // Apply sim frame offset to the goal as well
  if (par_.use_frame_alignment && !sim_frame_offset_.isIdentity(1e-9)) {
    Eigen::Vector4d g;
    g << G.pos(0), G.pos(1), G.pos(2), 1.0;
    g = sim_frame_offset_ * g;
    msg.goal.push_back(g(0));
    msg.goal.push_back(g(1));
    msg.goal.push_back(g(2));
  } else {
    msg.goal.push_back(G.pos(0));
    msg.goal.push_back(G.pos(1));
    msg.goal.push_back(G.pos(2));
  }

  // Publish the trajectory
  pub_own_traj_->publish(msg);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the trajectory the agent actually followed.
 *
 * Sando-style: keep a bounded history of past states and re-publish a single
 * persistent LINE_STRIP marker (id=1) colored by speed each call. RViz only
 * has to render one marker (downsampled to ~max_points_vis_), and there's no
 * marker accumulation between cycles.
 *
 * The previous implementation published a fresh ARROW Marker per sample with
 * an ever-incrementing id, which leaked thousands of markers into RViz, was
 * impossible to clean up reliably, and had a per-message Marker (not
 * MarkerArray) topic type — causing the rosbag "topic has more than one type"
 * error whenever a sando-style consumer was on the same graph.
 */
void MIGHTY_NODE::publishActualTraj() {
  if (!pub_actual_traj_) return;

  state current_state;
  mighty_ptr_->getState(current_state);
  const Eigen::Vector3d current_pos = current_state.pos;

  // Skip until the state is actually populated.
  if (current_pos.norm() < 1e-2) return;

  const auto now = this->now();
  const double tnow = now.seconds();

  // Initialize on first valid sample.
  if (!actual_traj_initialized_) {
    actual_traj_prev_pos_  = current_pos;
    actual_traj_prev_time_ = tnow;

    // Ground-robot velocity is unknown on the very first sample.
    if (par_.vehicle_type != "uav") current_state.vel.setZero();

    actual_traj_hist_.clear();
    actual_traj_hist_.push_back(current_state);
    actual_traj_initialized_ = true;
    return;  // wait for second sample to draw a line
  }

  // Velocity handling:
  //   UAV       -> trust current_state.vel from the estimator/sim
  //   Ground    -> approximate from position diff (TF-based state publisher)
  if (par_.vehicle_type != "uav") {
    const double dt = tnow - actual_traj_prev_time_;
    if (dt > 1e-3)
      current_state.vel = (current_pos - actual_traj_prev_pos_) / dt;
    else
      current_state.vel.setZero();
  }

  actual_traj_prev_pos_  = current_pos;
  actual_traj_prev_time_ = tnow;

  // Append to history only if it moved enough — prevents dense duplicates
  // when the robot is stationary. Still update the latest sample's velocity
  // so colors stay accurate.
  const double eps = 1e-3;
  if (!actual_traj_hist_.empty()) {
    const Eigen::Vector3d last_pos = actual_traj_hist_.back().pos;
    if ((current_pos - last_pos).norm() < eps) {
      actual_traj_hist_.back().vel = current_state.vel;
    } else {
      actual_traj_hist_.push_back(current_state);
    }
  } else {
    actual_traj_hist_.push_back(current_state);
  }

  // Bound history to keep RViz responsive.
  if (actual_traj_hist_.size() > actual_traj_max_hist_) {
    const size_t overflow = actual_traj_hist_.size() - actual_traj_max_hist_;
    actual_traj_hist_.erase(actual_traj_hist_.begin(),
                            actual_traj_hist_.begin() + overflow);
  }

  // Build a single persistent colored LINE_STRIP. par_.v_max is per-axis;
  // matching sando we use it directly as the jet-colormap max.
  const double vmax_for_color = par_.v_max;
  visualization_msgs::msg::MarkerArray ma = stateVector2ColoredLineStripMarkerArray(
      actual_traj_hist_,
      /*id=*/1,
      /*ns=*/"actual_traj_" + id_str_,
      /*max_value=*/vmax_for_color,
      /*stamp=*/now,
      /*line_width=*/actual_traj_line_width_,
      /*max_points_vis=*/actual_traj_max_points_vis_);

  pub_actual_traj_->publish(ma);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish goal (setpoint)
 */
void MIGHTY_NODE::publishGoal() {
  // Initialize the goal
  state next_goal;

  // Get the next goal
  if (mighty_ptr_->getNextGoal(next_goal) && par_.use_state_update) {
    // Publish the goal (actual setpoint)
    dynus_interfaces::msg::Goal quadGoal;
    quadGoal.header.stamp = this->now();
    quadGoal.header.frame_id = par_.map_frame_id;
    quadGoal.p = eigen2rosvector(next_goal.pos);
    quadGoal.v = eigen2rosvector(next_goal.vel);
    quadGoal.a = eigen2rosvector(next_goal.accel);
    quadGoal.j = eigen2rosvector(next_goal.jerk);
    quadGoal.yaw = next_goal.yaw;
    quadGoal.dyaw = next_goal.dyaw;
    pub_goal_->publish(quadGoal);

    // Debug log: every published setpoint, with current state and plan front
    // for context. Helps trace dip-toward-origin bugs.
    if (auto& s = debug_log_stream(); s.is_open()) {
      state cur, gterm, gsub;
      mighty_ptr_->getState(cur);
      mighty_ptr_->getGterm(gterm);
      mighty_ptr_->getG(gsub);
      s << debug_log_t(this->now().seconds()) << " GOAL "
        << next_goal.pos.x() << " " << next_goal.pos.y() << " " << next_goal.pos.z()
        << " state=" << cur.pos.x() << "," << cur.pos.y() << "," << cur.pos.z()
        << " G=" << gsub.pos.x() << "," << gsub.pos.y() << "," << gsub.pos.z()
        << " Gterm=" << gterm.pos.x() << "," << gterm.pos.y() << "," << gterm.pos.z()
        << "\n";
      s.flush();
    }

    // Publish the goal (setpoint) for visualization
    if (par_.visual_level >= 1) publishState(next_goal, pub_setpoint_);
  }

  // Publish FOV
  if (par_.visual_level >= 1) publishFOV();
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the full trajectory for robust tracking with replan support
 */
void MIGHTY_NODE::publishTrajectory() {
  // Retrieve the goal setpoints (full trajectory)
  mighty_ptr_->retrieveGoalSetpoints(goal_setpoints_);

  if (goal_setpoints_.empty()) {
    RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 1000,
                         "Cannot publish trajectory: goal_setpoints_ is empty");
    return;
  }

  // Increment trajectory_id for new trajectory (replan detected by change in ID)
  trajectory_id_++;

  // Downsample trajectory: publish ~50 points instead of all points
  const size_t target_points = par_.trajectory_downsample_points;
  size_t total_points = goal_setpoints_.size();
  size_t step = std::max(size_t(1), total_points / target_points);

  // Create trajectory message
  dynus_interfaces::msg::Trajectory traj_msg;
  traj_msg.header.stamp = this->now();
  traj_msg.header.frame_id = par_.map_frame_id;
  traj_msg.trajectory_id = trajectory_id_;
  traj_msg.dt = par_.dc * step;  // Adjusted time step for downsampled trajectory

  // Downsample: take every 'step' points
  traj_msg.goals.reserve(target_points);
  for (size_t i = 0; i < goal_setpoints_.size(); i += step) {
    const auto& state_point = goal_setpoints_[i];
    dynus_interfaces::msg::Goal goal;
    goal.p = eigen2rosvector(state_point.pos);
    goal.v = eigen2rosvector(state_point.vel);
    goal.a = eigen2rosvector(state_point.accel);
    goal.j = eigen2rosvector(state_point.jerk);
    goal.yaw = state_point.yaw;
    goal.dyaw = state_point.dyaw;
    traj_msg.goals.push_back(goal);
  }

  // Always include the last point
  if (!goal_setpoints_.empty() && (goal_setpoints_.size() - 1) % step != 0) {
    const auto& last_point = goal_setpoints_.back();
    dynus_interfaces::msg::Goal goal;
    goal.p = eigen2rosvector(last_point.pos);
    goal.v = eigen2rosvector(last_point.vel);
    goal.a = eigen2rosvector(last_point.accel);
    goal.j = eigen2rosvector(last_point.jerk);
    goal.yaw = last_point.yaw;
    goal.dyaw = last_point.dyaw;
    traj_msg.goals.push_back(goal);
  }

  // Publish trajectory
  pub_trajectory_->publish(traj_msg);
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::publishMpcPath() {
  // Get optimized trajectory setpoints (smooth quintic Hermite spline)
  mighty_ptr_->retrieveGoalSetpoints(goal_setpoints_);

  if (goal_setpoints_.size() < 2) return;

  dynus_interfaces::msg::SpeedyPath path_msg;
  path_msg.header.stamp = this->now();
  path_msg.header.frame_id = par_.map_frame_id;

  const double spacing = par_.mpc_path_spacing;

  // Collect speeds alongside poses
  std::vector<double> speeds;

  // Always include first point
  auto addPose = [&](const state& sp) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = path_msg.header;
    pose.pose.position.x = sp.pos.x();
    pose.pose.position.y = sp.pos.y();
    pose.pose.position.z = sp.pos.z();
    pose.pose.orientation.w = 1.0;
    path_msg.poses.push_back(pose);
    speeds.push_back(sp.vel.head<2>().norm());
  };

  addPose(goal_setpoints_.front());
  Eigen::Vector3d last_emitted = goal_setpoints_.front().pos;

  // Walk through setpoints, emit when distance threshold exceeded
  for (size_t i = 1; i < goal_setpoints_.size(); ++i) {
    double dist = (goal_setpoints_[i].pos - last_emitted).head<2>().norm();
    if (dist >= spacing) {
      addPose(goal_setpoints_[i]);
      last_emitted = goal_setpoints_[i].pos;
    }
  }

  // Always include last point
  const auto& last = goal_setpoints_.back();
  if ((last.pos - last_emitted).head<2>().norm() > 1e-6) addPose(last);

  // Need at least 2 points for MPC
  if (path_msg.poses.size() < 2) return;

  // Set yaw from direction to next waypoint
  for (size_t i = 0; i < path_msg.poses.size(); ++i) {
    double yaw = 0.0;
    if (i + 1 < path_msg.poses.size()) {
      double dx = path_msg.poses[i + 1].pose.position.x - path_msg.poses[i].pose.position.x;
      double dy = path_msg.poses[i + 1].pose.position.y - path_msg.poses[i].pose.position.y;
      yaw = std::atan2(dy, dx);
    } else if (i > 0) {
      double dx = path_msg.poses[i].pose.position.x - path_msg.poses[i - 1].pose.position.x;
      double dy = path_msg.poses[i].pose.position.y - path_msg.poses[i - 1].pose.position.y;
      yaw = std::atan2(dy, dx);
    }
    path_msg.poses[i].pose.orientation.x = 0.0;
    path_msg.poses[i].pose.orientation.y = 0.0;
    path_msg.poses[i].pose.orientation.z = std::sin(yaw / 2.0);
    path_msg.poses[i].pose.orientation.w = std::cos(yaw / 2.0);
  }

  // Assign per-waypoint speeds
  path_msg.speeds = speeds;

  pub_mpc_path_->publish(path_msg);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish Sefe Corridor Polyhedra
 */
void MIGHTY_NODE::publishPoly() {
  // retrieve the polyhedra
  mighty_ptr_->retrievePolytopes(poly_whole_, poly_safe_);

  // For whole trajectory
  if (!poly_whole_.empty()) {
    decomp_ros_msgs::msg::PolyhedronArray poly_whole_msg =
        DecompROS::polyhedron_array_to_ros(poly_whole_);
    poly_whole_msg.header.stamp = this->now();
    poly_whole_msg.header.frame_id = par_.map_frame_id;
    poly_whole_msg.lifetime = rclcpp::Duration::from_seconds(1.0);
    pub_poly_whole_->publish(poly_whole_msg);
  }

  // For safe trajectory
  if (!poly_safe_.empty()) {
    decomp_ros_msgs::msg::PolyhedronArray poly_safe_msg =
        DecompROS::polyhedron_array_to_ros(poly_safe_);
    poly_safe_msg.header.stamp = this->now();
    poly_safe_msg.header.frame_id = par_.map_frame_id;
    poly_safe_msg.lifetime = rclcpp::Duration::from_seconds(1.0);
    pub_poly_safe_->publish(poly_safe_msg);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the trajectory
 */
void MIGHTY_NODE::publishTraj() {
  auto now = this->now();

  // Build ONE MarkerArray that contains DELETEALL followed by the new markers,
  // for each topic. Publishing DELETEALL and the new markers as two separate
  // messages was racy: if the second message was dropped or arrived out of
  // order, RViz would show an empty topic. A single MarkerArray is processed
  // atomically — RViz applies DELETEALL first, then the ADDs that follow.
  auto deleteAll = [&]() {
    visualization_msgs::msg::Marker m;
    m.header.frame_id = par_.map_frame_id;
    m.header.stamp    = rclcpp::Time(0, 0, RCL_ROS_TIME);  // latest TF
    m.action          = visualization_msgs::msg::Marker::DELETEALL;
    return m;
  };

  // 1) Committed (best) trajectory
  mighty_ptr_->retrieveGoalSetpoints(goal_setpoints_);
  {
    visualization_msgs::msg::MarkerArray committed_ma;
    committed_ma.markers.push_back(deleteAll());
    auto added = stateVector2ColoredMarkerArray(
        goal_setpoints_, /*type=*/1, par_.v_max, now, par_.map_frame_id);
    committed_ma.markers.insert(committed_ma.markers.end(),
                                added.markers.begin(), added.markers.end());
    pub_traj_committed_colored_->publish(committed_ma);
  }

  // 2) Sub-optimal trajectories
  if (par_.use_multiple_initial_guesses) {
    mighty_ptr_->retrieveListSubOptGoalSetpoints(list_subopt_goal_setpoints_);

    visualization_msgs::msg::MarkerArray subopt_ma;
    subopt_ma.markers.push_back(deleteAll());
    for (int i = 0; i < (int)list_subopt_goal_setpoints_.size(); ++i) {
      auto single =
          stateVector2ColoredMarkerArray(list_subopt_goal_setpoints_[i],
                                         /*type=*/i + 2, par_.v_max, now, par_.map_frame_id);
      subopt_ma.markers.insert(subopt_ma.markers.end(), single.markers.begin(),
                               single.markers.end());
    }
    pub_traj_subopt_colored_->publish(subopt_ma);
  } else {
    // Even when subopt is disabled, still publish a DELETEALL so the topic
    // doesn't accumulate stale markers from a previous toggle.
    visualization_msgs::msg::MarkerArray clear_only;
    clear_only.markers.push_back(deleteAll());
    pub_traj_subopt_colored_->publish(clear_only);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the global path (that can go through unknown space)
 */
void MIGHTY_NODE::publishGlobalPath() {
  int global_path_color = RED;
  int original_global_path_color = ORANGE;

  // Generate random integer from 1 to 10 to generate random color
  if (par_.use_random_color_for_global_path) global_path_color = rand() % 10 + 1;

  // Get global_path
  vec_Vecf<3> global_path;
  mighty_ptr_->getGlobalPath(global_path);

  if (!global_path.empty()) {
    // Publish global_path (thin line + dots)
    clearMarkerArray(hgp_path_marker_, pub_hgp_path_marker_);

    pathLineDotsToMarkerArray(global_path, &hgp_path_marker_, color(global_path_color),
                              /*line_width=*/0.03,    // meters
                              /*dot_diameter=*/0.06,  // meters
                              /*base_id=*/50000,
                              /*frame_id=*/par_.map_frame_id,
                              // Long lifetime so the path stays visible across
                              // intermittent HGP failures (which are common
                              // when the robot is near the edge of the local
                              // mapper window). The next successful publish
                              // calls clearMarkerArray to wipe stale markers,
                              // so we don't accumulate.
                              /*lifetime_sec=*/10.0);

    pub_hgp_path_marker_->publish(hgp_path_marker_);
  }

  // Get the original global path
  vec_Vecf<3> original_global_path;
  mighty_ptr_->getOriginalGlobalPath(original_global_path);

  if (!original_global_path.empty()) {
    // Publish original_global_path (thin line + dots)
    clearMarkerArray(original_hgp_path_marker_, pub_original_hgp_path_marker_);

    pathLineDotsToMarkerArray(original_global_path, &original_hgp_path_marker_,
                              color(original_global_path_color),
                              /*line_width=*/0.03,    // meters
                              /*dot_diameter=*/0.06,  // meters
                              /*base_id=*/60000,
                              /*frame_id=*/par_.map_frame_id,
                              /*lifetime_sec=*/10.0);

    pub_original_hgp_path_marker_->publish(original_hgp_path_marker_);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the free global path (that only goes through free space)
 */
void MIGHTY_NODE::publishFreeGlobalPath() {
  // Get free_global_path
  vec_Vecf<3> free_global_path;
  mighty_ptr_->getFreeGlobalPath(free_global_path);

  if (free_global_path.empty()) return;

  // Publish free_global_path
  clearMarkerArray(hgp_free_path_marker_, pub_free_hgp_path_marker_);
  vectorOfVectors2MarkerArray(free_global_path, &hgp_free_path_marker_, color(GREEN),
                              visualization_msgs::msg::Marker::ARROW, {}, par_.map_frame_id);
  pub_free_hgp_path_marker_->publish(hgp_free_path_marker_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the local_global_path and local_global_path_after_push_
 */
void MIGHTY_NODE::publishLocalGlobalPath() {
  // Get the local global path and local global path after push
  vec_Vecf<3> local_global_path;
  vec_Vecf<3> local_global_path_after_push;
  mighty_ptr_->getLocalGlobalPath(local_global_path, local_global_path_after_push);

  if (!local_global_path.empty()) {
    // Publish local_global_path
    clearMarkerArray(hgp_local_global_path_marker_, pub_local_global_path_marker_);
    vectorOfVectors2MarkerArray(local_global_path, &hgp_local_global_path_marker_, color(BLUE),
                                visualization_msgs::msg::Marker::ARROW, {}, par_.map_frame_id);
    pub_local_global_path_marker_->publish(hgp_local_global_path_marker_);
  }

  if (!local_global_path_after_push.empty()) {
    // Publish local_global_path_after_push
    clearMarkerArray(hgp_local_global_path_after_push_marker_,
                     pub_local_global_path_after_push_marker_);
    vectorOfVectors2MarkerArray(local_global_path_after_push,
                                &hgp_local_global_path_after_push_marker_, color(ORANGE),
                                visualization_msgs::msg::Marker::ARROW, {}, par_.map_frame_id);
    pub_local_global_path_after_push_marker_->publish(hgp_local_global_path_after_push_marker_);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Create MarkerArray from vec_Vec3f
 */
void MIGHTY_NODE::createMarkerArrayFromVec_Vec3f(
    const vec_Vec3f& occupied_cells, const std_msgs::msg::ColorRGBA& color, int namespace_id,
    double scale, visualization_msgs::msg::MarkerArray* marker_array) {
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = par_.map_frame_id;
  marker.header.stamp = rclcpp::Clock().now();  // Use ROS2 clock
  marker.ns = "namespace_" + std::to_string(namespace_id);
  marker.id = 0;
  marker.type =
      visualization_msgs::msg::Marker::CUBE_LIST;  // Each point will be visualized as a cube
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.scale.x = par_.res;
  marker.scale.y = par_.res;
  marker.scale.z = par_.res;
  marker.color = color;

  for (const auto& cell : occupied_cells) {
    geometry_msgs::msg::Point point;
    point.x = cell(0);
    point.y = cell(1);
    point.z = cell(2);
    marker.points.push_back(point);
  }

  marker_array->markers.push_back(marker);
}

// ----------------------------------------------------------------------------

/**
 * @brief Clear any marker array
 */
void MIGHTY_NODE::clearMarkerArray(
    visualization_msgs::msg::MarkerArray& path_marker,
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr publisher) {
  // If the marker array is empty, return
  if (path_marker.markers.size() == 0) return;

  // Clear the marker array
  int id_begin = path_marker.markers[0].id;

  for (int i = 0; i < path_marker.markers.size(); i++) {
    visualization_msgs::msg::Marker m;
    m.type = visualization_msgs::msg::Marker::ARROW;
    m.action = visualization_msgs::msg::Marker::DELETE;
    m.id = i + id_begin;
    path_marker.markers[i] = m;
  }

  publisher->publish(path_marker);
  path_marker.markers.clear();
}

// ----------------------------------------------------------------------------

/**
 * @brief Construct the FOV marker for visualization
 */
void MIGHTY_NODE::constructFOVMarker() {
  marker_fov_.header.stamp = this->now();
  marker_fov_.header.frame_id = d435_depth_frame_id_;
  marker_fov_.ns = "marker_fov";
  marker_fov_.id = marker_fov_id_++;
  marker_fov_.frame_locked = true;
  marker_fov_.action = marker_fov_.ADD;
  marker_fov_.pose = mighty_utils::identityGeometryMsgsPose();

  // Sphere sensing mode: render a translucent sphere instead of the camera frustum.
  if (par_.use_sphere_sensing) {
    marker_fov_.type = marker_fov_.SPHERE;
    marker_fov_.points.clear();
    const double diameter = 2.0 * par_.sphere_sensing_radius;
    marker_fov_.scale.x = diameter;
    marker_fov_.scale.y = diameter;
    marker_fov_.scale.z = diameter;
    marker_fov_.color.a = 0.15;
    marker_fov_.color.r = 0.0;
    marker_fov_.color.g = 1.0;
    marker_fov_.color.b = 0.0;
    return;
  }

  marker_fov_.type = marker_fov_.LINE_LIST;

  double delta_y = par_.fov_visual_depth * fabs(tan((par_.fov_visual_x_deg * M_PI / 180) / 2.0));
  double delta_z = par_.fov_visual_depth * fabs(tan((par_.fov_visual_y_deg * M_PI / 180) / 2.0));

  geometry_msgs::msg::Point v0 = eigen2point(Eigen::Vector3d(0.0, 0.0, 0.0));
  geometry_msgs::msg::Point v1 =
      eigen2point(Eigen::Vector3d(-delta_y, delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v2 =
      eigen2point(Eigen::Vector3d(delta_y, delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v3 =
      eigen2point(Eigen::Vector3d(delta_y, -delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v4 =
      eigen2point(Eigen::Vector3d(-delta_y, -delta_z, par_.fov_visual_depth));

  marker_fov_.points.clear();

  // Line
  marker_fov_.points.push_back(v0);
  marker_fov_.points.push_back(v1);

  // Line
  marker_fov_.points.push_back(v0);
  marker_fov_.points.push_back(v2);

  // Line
  marker_fov_.points.push_back(v0);
  marker_fov_.points.push_back(v3);

  // Line
  marker_fov_.points.push_back(v0);
  marker_fov_.points.push_back(v4);

  // Line
  marker_fov_.points.push_back(v1);
  marker_fov_.points.push_back(v2);

  // Line
  marker_fov_.points.push_back(v2);
  marker_fov_.points.push_back(v3);

  // Line
  marker_fov_.points.push_back(v3);
  marker_fov_.points.push_back(v4);

  // Line
  marker_fov_.points.push_back(v4);
  marker_fov_.points.push_back(v1);

  marker_fov_.scale.x = 0.03;
  marker_fov_.scale.y = 0.00001;
  marker_fov_.scale.z = 0.00001;
  marker_fov_.color.a = 1.0;
  marker_fov_.color.r = 0.0;
  marker_fov_.color.g = 1.0;
  marker_fov_.color.b = 0.0;
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the FOV marker for visualization
 */
void MIGHTY_NODE::publishFOV() {
  marker_fov_.header.stamp = this->now();
  pub_fov_->publish(marker_fov_);
  return;
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::mapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& map_msg,
                              const sensor_msgs::msg::PointCloud2::ConstPtr& unk_msg) {
  RCLCPP_INFO_ONCE(this->get_logger(),
                   "mapCallback triggered — synced occupancy_grid + unknown_grid received");

  // use PCL's own Ptr (boost::shared_ptr)
  pcl::PointCloud<pcl::PointXYZ>::Ptr map_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*map_msg, *map_pc);

  pcl::PointCloud<pcl::PointXYZ>::Ptr unk_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*unk_msg, *unk_pc);

  mighty_ptr_->updateMap(map_pc, unk_pc);

  // Publish heat cloud visualization if enabled
  if (par_.use_heat_map) {
    publishHeatCloud();
  }
  publishGround2DOccupied();
  publishGround2DHeat();
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::occupancyMapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& map_msg) {
  // use PCL's own Ptr (boost::shared_ptr)
  pcl::PointCloud<pcl::PointXYZ>::Ptr map_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*map_msg, *map_pc);

  mighty_ptr_->updateOccupancyMap(map_pc);

  // Publish heat cloud visualization if enabled
  if (par_.use_heat_map) {
    publishHeatCloud();
  }
  publishGround2DOccupied();
  publishGround2DHeat();
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::unknownMapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr& unk_msg) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr unk_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*unk_msg, *unk_pc);
  mighty_ptr_->updateUnknownCloud(unk_pc);
}

void MIGHTY_NODE::esdfCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg) {
  esdf_grid_ = EsdfGrid2D::fromOccupancyGrid(*msg, par_.esdf_truncation_distance);
}

void MIGHTY_NODE::occ2DCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg) {
  // Remember the global mapper's ground-plane z so publishVisitedMap() can
  // render at the same height as the live occ_2d layer in RViz.
  occ2d_origin_z_ = msg->info.origin.position.z;

  // Persistent-map fusion: any cell that arrived UNKNOWN from the mapper but
  // was previously observed (FREE or OCCUPIED) gets restored from
  // visited_map_ before we build OccGrid2D. So when the robot revisits a
  // previously-explored region, the sliding window comes back instantly with
  // its last-known occupancy instead of one frame of UNKNOWN flicker. This
  // re-introduces stale OCCUPIED for moving obstacles that have since left,
  // so it's only safe in static environments — gated by a parameter.
  // Mutating msg->data in place is safe here because (a) we're the only
  // subscriber to occ_2d_topic in this process, (b) cb_group_map_ is now
  // MutuallyExclusive, and (c) inter-process delivery gives us a unique copy.
  if (visited_map_ && par_.expl_fuse_persistent_into_local && !msg->data.empty()) {
    const double res = msg->info.resolution;
    const double ox  = msg->info.origin.position.x;
    const double oy  = msg->info.origin.position.y;
    const unsigned W = msg->info.width;
    const unsigned H = msg->info.height;
    for (unsigned iy = 0; iy < H; ++iy) {
      for (unsigned ix = 0; ix < W; ++ix) {
        const size_t i = static_cast<size_t>(iy) * W + ix;
        if (i >= msg->data.size()) break;
        if (msg->data[i] >= 0) continue;  // already known (FREE or OCCUPIED)
        const double wx = ox + (ix + 0.5) * res;
        const double wy = oy + (iy + 0.5) * res;
        const int8_t v = visited_map_->getStateWorld(wx, wy);
        if (v != VisitedMap::kUnknown) {
          msg->data[i] = v;  // restore old persistent value
        }
      }
    }
  }

  occ_grid_2d_ = OccGrid2D::fromOccupancyGrid(*msg);

  // Choose the grid handed to the GLOBAL planner. By default it's the sliding
  // occ_2d window (legacy). When use_persistent_global_map is on (ground robot),
  // feed the accumulated visited_map_ instead, so the A* plans over the whole
  // mission history — a stable route that doesn't churn as the window slides.
  // The node keeps `occ_grid_2d_` (the sliding grid) for exploration/frontier
  // use; only the core planner's copy is swapped. Falls back to the sliding grid
  // until the persistent map exists and is non-empty (e.g. exploration disabled
  // -> no visited_map_). absorb() here folds in the newest observation first;
  // the exploration block below may absorb the same (small) sliding grid again,
  // which is idempotent and negligible.
  std::shared_ptr<const OccGrid2D> planner_grid = occ_grid_2d_;
  if (par_.use_persistent_global_map && par_.vehicle_type == "ground_robot"
      && visited_map_) {
    visited_map_->absorb(*occ_grid_2d_);
    if (!visited_map_->empty()) {
      planner_grid = OccGrid2D::fromTristate(
          visited_map_->width(), visited_map_->height(),
          visited_map_->resolution(),
          visited_map_->originX(), visited_map_->originY(),
          visited_map_->data());
    }
  }
  mighty_ptr_->setOccGrid2D(planner_grid);

  // Frontier-based exploration: detect frontiers in the new grid, update the
  // persistent global database, then immediately try to issue an exploration
  // goal so the robot starts moving the moment the first frontier appears
  // (instead of waiting up to 1 s for the explore-select timer tick).
  if (par_.expl_enabled && occ_grid_2d_ && frontier_detector_ && frontier_manager_
      && state_initialized_) {
    state cur;
    mighty_ptr_->getState(cur);
    Eigen::Vector3d robot_pose(cur.pos.x(), cur.pos.y(), cur.yaw);
    Eigen::Vector2d robot_xy(cur.pos.x(), cur.pos.y());

    // Absorb the freshly observed cells into the persistent visited bitmap
    // *before* running the detector, so cells we are observing right now are
    // already marked visited and never become "stale unknown" on the next
    // sliding step.
    if (visited_map_) visited_map_->absorb(*occ_grid_2d_);

    // Pick the grid the detector runs on. The persistent visited_map_ is the
    // entire mission history, so frontiers at the boundary of explored area
    // are always reachable via WFD even when the robot is far from them. The
    // local sliding window is the legacy path; both paths feed the visited
    // map as the suppression filter when on the local grid.
    if (par_.expl_detect_on_visited_map && visited_map_ && !visited_map_->empty()) {
      current_detect_grid_ = OccGrid2D::fromTristate(
          visited_map_->width(), visited_map_->height(),
          visited_map_->resolution(),
          visited_map_->originX(), visited_map_->originY(),
          visited_map_->data());
    } else {
      current_detect_grid_ = occ_grid_2d_;
    }
    const auto& detect_grid = *current_detect_grid_;
    const VisitedMap* visited_filter =
        par_.expl_detect_on_visited_map ? nullptr : visited_map_.get();
    auto clusters = frontier_detector_->detect(detect_grid, robot_xy,
                                               visited_filter);

    // ESDF-based clearance filter: drop frontiers whose centroid is closer
    // than `expl_min_obstacle_distance_m` to the nearest obstacle. The ESDF
    // gives a meter-accurate distance and matches what HGP/L-BFGS use for
    // collision avoidance, so frontiers we keep are guaranteed to have
    // breathing room from walls. Disabled when threshold <= 0 or no ESDF.
    if (esdf_grid_ && par_.expl_min_obstacle_distance_m > 0.0) {
      const double thresh = par_.expl_min_obstacle_distance_m;
      clusters.erase(
          std::remove_if(
              clusters.begin(), clusters.end(),
              [&](const FrontierCluster& c) {
                return esdf_grid_->queryDistance(c.centroid.x(),
                                                 c.centroid.y()) < thresh;
              }),
          clusters.end());
    }

    const auto active_peers = par_.expl_use_minpos
        ? peer_tracker_.getActivePeers(this->now().seconds(),
                                       par_.expl_peer_timeout_sec)
        : std::vector<PeerPose>{};
    frontier_manager_->update(clusters, detect_grid, robot_pose,
                              this->now().seconds(), active_peers);

    // Also retroactively invalidate existing records that drifted too close
    // to obstacles (e.g. via EMA centroid updates) or that were inserted
    // before the ESDF caught up. Without this, stale records that already
    // hugged a wall would never be cleared.
    if (esdf_grid_ && par_.expl_min_obstacle_distance_m > 0.0) {
      const double thresh = par_.expl_min_obstacle_distance_m;
      for (const auto& r : frontier_manager_->records()) {
        if (r.state != FrontierState::ACTIVE &&
            r.state != FrontierState::DORMANT) continue;
        if (esdf_grid_->queryDistance(r.centroid_xy.x(),
                                      r.centroid_xy.y()) < thresh) {
          frontier_manager_->markInvalidated(r.id, this->now().seconds());
        }
      }
    }

    // Re-validate persistent records against the detector's own frontier
    // definition. update()'s Step-d keeps a record ACTIVE as long as a single
    // UNKNOWN cell sits within verify_radius_cells of its centroid, which
    // disagrees with the detector's cluster_min_cells threshold and leaves
    // phantom frontiers on cells that were a border earlier but got explored
    // (sub-threshold remnants, or permanently-occluded specks that never
    // clear). Ask the detector "is this still a frontier here?" and retire
    // (VISITED) any ACTIVE/DORMANT record whose location no longer qualifies.
    // Uses the same grid + visited filter detection just ran on, so it is
    // exactly consistent with what the detector would (not) emit.
    if (par_.expl_reconcile_stale && frontier_detector_) {
      const int min_cells = frontier_detector_->params().cluster_min_cells;
      const double res = detect_grid.resolution();
      // Search out to the manager's merge radius: anything within matching
      // distance of a live frontier is "the same" frontier, so only records
      // with no qualifying frontier within that radius are retired.
      const int search_cells = std::max(1, static_cast<int>(std::ceil(
          frontier_manager_->params().merge_radius_m / std::max(1e-6, res))));
      for (const auto& r : frontier_manager_->records()) {
        if (r.state != FrontierState::ACTIVE &&
            r.state != FrontierState::DORMANT) continue;
        if (!frontier_detector_->isStillFrontier(
                detect_grid, r.centroid_xy, min_cells, search_cells,
                visited_filter)) {
          frontier_manager_->markVisited(r.id);
        }
      }
    }

    // Publish the persistent occupancy map ~1 Hz so RViz can layer it behind
    // the sliding occ_2d. Throttling matters because each publish copies the
    // full persistent buffer (~444 KB at the 100×100 m default; bigger if the
    // user enlarges expl_visited_map_width_m / height_m) and the map only
    // changes incrementally between frames. transient_local QoS guarantees
    // late-joining RViz still gets the latest snapshot.
    if (visited_map_ && par_.expl_publish_visited_map) {
      const double t_now = this->now().seconds();
      if (t_now - last_visited_publish_t_ >= 1.0) {
        publishVisitedMap();
        last_visited_publish_t_ = t_now;
      }
    }

    // Broadcast our visited map to peers so they can skip already-explored
    // frontiers. Same 1 Hz throttle as the local RViz publish.
    if (pub_peer_visited_map_ && visited_map_ && !visited_map_->empty()) {
      const double t_now = this->now().seconds();
      if (t_now - last_peer_visited_publish_t_ >= 1.0) {
        nav_msgs::msg::OccupancyGrid msg;
        msg.header.frame_id = ns_;
        msg.header.stamp    = this->now();
        msg.info.resolution = static_cast<float>(visited_map_->resolution());
        msg.info.width      = static_cast<unsigned>(visited_map_->width());
        msg.info.height     = static_cast<unsigned>(visited_map_->height());
        msg.info.origin.position.x = visited_map_->originX();
        msg.info.origin.position.y = visited_map_->originY();
        msg.info.origin.orientation.w = 1.0;
        const auto& v = visited_map_->data();
        msg.data.assign(v.begin(), v.end());
        pub_peer_visited_map_->publish(msg);
        last_peer_visited_publish_t_ = t_now;
      }
    }

    // Throttled diagnostic so it's obvious whether detection is finding
    // anything (and gives a hint about *why* if the answer is "no").
    {
      // Quick scan of the published grid for cell-state distribution.
      int n_unknown = 0, n_free = 0, n_occ = 0;
      const auto& occ = occ_grid_2d_->occupiedData();
      const auto& unk = occ_grid_2d_->unknownData();
      for (size_t i = 0; i < occ.size(); ++i) {
        if (unk[i])      ++n_unknown;
        else if (occ[i]) ++n_occ;
        else             ++n_free;
      }
      RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 2000,
          "[expl] grid %dx%d  free=%d  occ=%d  unknown=%d  fresh=%zu  db=%zu",
          occ_grid_2d_->width(), occ_grid_2d_->height(),
          n_free, n_occ, n_unknown, clusters.size(), frontier_manager_->size());
    }

    if (par_.expl_publish_markers) publishFrontierMarkers();
    publishFrontierData();

    // Drive the goal-selection loop immediately. This makes the robot start
    // moving as soon as the first frontier exists, rather than waiting for
    // the next 1 Hz explore-select tick.
    exploreSelectCallback();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Frontier exploration: pick the next goal from the global frontier DB
 *        and issue it through the same pathway as a manual term_goal. Skipped
 *        when a manual goal is active or our own previous goal is still being
 *        pursued.
 */
void MIGHTY_NODE::exploreSelectCallback() {
  if (!par_.expl_enabled) return;
  if (manual_goal_active_) return;
  if (!occ_grid_2d_ || !frontier_manager_) return;
  if (!current_detect_grid_) return;
  if (!state_initialized_) return;

  // Return-home handling. While returning home, suppress frontier selection so
  // the agent isn't yanked back out mid-return by a newly-detected frontier.
  // But once it actually REACHES home, if the latch came from the internal
  // "nothing left to explore" path, lift it and fall through to re-check
  // frontiers — so a frontier discovered during the return (or a leftover that
  // reappeared) gets explored instead of sitting drawn-but-unexplored forever.
  // An external /exploration/return_home command leaves the latch set (the
  // operator asked the agent to come home and stay parked).
  if (home_return_requested_) {
    if (!home_return_resume_on_arrival_) return;  // external command: stay parked
    state cur_home;
    mighty_ptr_->getState(cur_home);
    const double d_home = std::hypot(cur_home.pos.x() - exploration_start_pos_.x(),
                                     cur_home.pos.y() - exploration_start_pos_.y());
    if (d_home > par_.goal_radius) return;  // still en route home — keep the latch
    // Reached home: lift the latch and re-evaluate for remaining frontiers.
    home_return_requested_ = false;
    exploration_active_ = false;
    unreachable_consec_count_ = 0;
    RCLCPP_INFO(this->get_logger(),
                "Exploration: reached home — re-evaluating for remaining frontiers");
  }

  state cur;
  mighty_ptr_->getState(cur);
  Eigen::Vector3d robot_pose(cur.pos.x(), cur.pos.y(), cur.yaw);

  // If we still have an in-progress exploration goal that hasn't been
  // marked VISITED/INVALIDATED yet, leave it alone — unless we have stopped
  // making progress towards it.
  //
  // WATCHDOG. Every other way out of a pursued frontier depends on something
  // that can itself stop happening:
  //   * the pursuit timeout and the dwell/visited checks live in
  //     FrontierManager::update(), which only runs from occ2DCallback, so they
  //     stop firing if map updates stop;
  //   * goalReachedCheck() needs drone_status_ == GOAL_REACHED *and*
  //     checkReadyToReplan();
  //   * the HGP-unreachable path needs the replan loop to keep failing.
  // If none of them fires, this branch returns early forever and the robot sits
  // still — observed once in a 10-run campaign, stranded 26.8 m from home with
  // 99.3% coverage and no further goals issued. This check runs on the
  // select timer, so it is independent of all of the above.
  if (exploration_active_) {
    auto* r = frontier_manager_->find(current_explore_id_);
    if (r && (r->state == FrontierState::ACTIVE ||
              r->state == FrontierState::DORMANT)) {
      const double t_now = this->now().seconds();
      const Eigen::Vector2d pos_xy(cur.pos.x(), cur.pos.y());
      // Progress = the robot moved, or we are pursuing a different frontier
      // than when the timer was last reset.
      if (explore_stall_id_ != current_explore_id_ ||
          (pos_xy - explore_stall_pos_).norm() > par_.expl_stall_eps_m) {
        explore_stall_id_  = current_explore_id_;
        explore_stall_pos_ = pos_xy;
        explore_stall_t_   = t_now;
      } else if (par_.expl_stall_timeout_sec > 0.0 &&
                 t_now - explore_stall_t_ > par_.expl_stall_timeout_sec) {
        RCLCPP_WARN(this->get_logger(),
                    "Exploration watchdog: no progress towards frontier %lu for "
                    "%.0f s at (%.2f, %.2f) — invalidating so selection can move on",
                    static_cast<unsigned long>(current_explore_id_),
                    t_now - explore_stall_t_, pos_xy.x(), pos_xy.y());
        // markInvalidated() also counts a strike, so a spot that repeatedly
        // stalls escalates to a permanent keep-out instead of being re-picked.
        frontier_manager_->markInvalidated(current_explore_id_, t_now);
        exploration_active_ = false;
        unreachable_consec_count_ = 0;
        explore_stall_t_ = t_now;
      }
      if (exploration_active_) return;
    }
    // Otherwise (record gone or already terminal) fall through and pick a new one.
  }

  // Multi-agent timing guard: peer-shared maps can populate frontiers before
  // this agent's own state has settled past the (0,0,0) default, causing the
  // start capture below to grab origin instead of the real spawn pose. Only
  // applies when at least one peer has been heard — on solo hardware (e.g.
  // RR04 DLIO) the rover legitimately starts at origin and the guard would
  // otherwise deadlock the planner (goal never issued -> rover never moves
  // -> guard never lifts).
  if (!exploration_start_captured_ && par_.expl_use_minpos) {
    const auto peers = peer_tracker_.getActivePeers(
        this->now().seconds(), par_.expl_peer_timeout_sec);
    if (!peers.empty()) {
      const double dist_to_origin =
          std::sqrt(cur.pos.x() * cur.pos.x() + cur.pos.y() * cur.pos.y());
      if (dist_to_origin < 0.5) {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 2000,
          "Exploration deferred: robot pose (%.3f, %.3f) within 0.5m of origin "
          "with %zu peer(s) active — waiting for state to settle",
          cur.pos.x(), cur.pos.y(), peers.size());
        return;
      }
    }
  }

  std::optional<FrontierRecord> next;
  if (par_.expl_select_nearest) {
    // Strict "always go to the closest frontier" policy (utility-independent).
    // Commitment: keep pursuing the current frontier while it stays selectable
    // unless a different one is closer by more than select_commit_margin_m.
    // Without this, two near-equidistant frontiers make the robot flip-flop,
    // which reads as "not going to the closest one".
    next = frontier_manager_->selectNearest(robot_pose);
    if (next && exploration_active_) {
      const FrontierRecord* cur = frontier_manager_->find(current_explore_id_);
      if (cur && (cur->state == FrontierState::ACTIVE ||
                  cur->state == FrontierState::DORMANT)) {
        const Eigen::Vector2d rxy = robot_pose.head<2>();
        const double d_cur  = (rxy - cur->centroid_xy).norm();
        const double d_near = (rxy - next->centroid_xy).norm();
        if (d_cur <= d_near + par_.expl_select_commit_margin_m) {
          FrontierRecord keep = *cur;
          keep.cached_utility = -d_cur;
          next = keep;  // stick with the current goal
        }
      }
    }
  } else if (par_.expl_use_minpos) {
    auto peers = peer_tracker_.getActivePeers(
        this->now().seconds(), par_.expl_peer_timeout_sec);
    next = frontier_manager_->selectNextGoalMinPos(
        robot_pose, *current_detect_grid_, peers,
        par_.expl_min_frontier_dist_to_peers_m);
  } else {
    next = frontier_manager_->selectNextGoal(robot_pose, *current_detect_grid_);
  }
  if (!next) {
    // Having nothing to pursue is what triggers the return home — NOT whether a
    // goal happened to be active when selection came up empty.
    //
    // This used to be gated on exploration_active_, which stranded the robot.
    // The unreachable-frontier path (5 consecutive HGP failures) clears
    // exploration_active_, so on the next tick selection returns nothing, this
    // whole block was skipped, and the callback returned silently — every tick,
    // forever. Observed with 322 frontier records still in the database and the
    // detector reporting fresh=10 per cycle: the agent simply stopped acting,
    // 26.7 m from home, until the harness killed the run.
    //
    // exploration_start_captured_ keeps this from firing before exploration has
    // actually begun, when an empty selection just means the map is still
    // filling in.
    if (exploration_start_captured_) {
      // Grace period: don't conclude "exploration finished" on the first empty
      // tick. selectNearest can momentarily return nothing while a real frontier
      // is transiently INVALIDATED (one-off HGP failure / ESDF filter / dwell)
      // or between detection cycles. Require the empty condition to persist for
      // expl_home_grace_sec so the detector/reconcile can resurface any
      // remaining reachable frontier — the agent keeps exploring instead of
      // heading home prematurely.
      const double t_now_grace = this->now().seconds();
      if (no_frontier_since_t_ < 0.0) no_frontier_since_t_ = t_now_grace;
      if (par_.expl_home_grace_sec > 0.0 &&
          (t_now_grace - no_frontier_since_t_) < par_.expl_home_grace_sec) {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 1000,
            "Exploration: no selectable frontier — waiting %.1fs grace before home "
            "(%.1fs elapsed)", par_.expl_home_grace_sec,
            t_now_grace - no_frontier_since_t_);
        return;  // stay in exploration; wait for a frontier to (re)appear
      }
      // Already parked at the captured start with nothing selectable: we are
      // simply finished. Do NOT issue a return-home goal here — it would be
      // satisfied on the spot, the arrival check on the next tick would lift the
      // return-home latch, selection would come up empty again, and this branch
      // would re-enter. Measured at 102 give-up/resume cycles in a single run,
      // several per second, once the give-up path stopped being gated on
      // exploration_active_.
      const double d_home_now =
          (Eigen::Vector2d(cur.pos.x(), cur.pos.y())
           - exploration_start_pos_.head<2>()).norm();
      if (d_home_now <= par_.goal_radius) {
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 10000,
            "Exploration complete — parked at the start (%.2f, %.2f) with nothing "
            "left to explore", cur.pos.x(), cur.pos.y());
        exploration_active_ = false;
        return;
      }

      std::cout << "[mighty] No frontiers left. Robot now at ("
                << cur.pos.x() << ", " << cur.pos.y() << ", " << cur.pos.z()
                << "). Returning to captured start ("
                << exploration_start_pos_.x() << ", "
                << exploration_start_pos_.y() << ", "
                << exploration_start_pos_.z() << ")"
                << std::endl;
      RCLCPP_INFO(this->get_logger(),
                  "Exploration: nothing left to explore — returning to start at (%.2f, %.2f, %.2f)",
                  exploration_start_pos_.x(), exploration_start_pos_.y(), exploration_start_pos_.z());
      geometry_msgs::msg::PoseStamped home;
      home.header.frame_id    = par_.map_frame_id;
      home.header.stamp       = this->now();
      home.pose.position.x    = exploration_start_pos_.x();
      home.pose.position.y    = exploration_start_pos_.y();
      home.pose.position.z    = par_.expl_default_goal_z;
      home.pose.orientation.w = 1.0;
      if (par_.expl_external_selector) {
        // External selector owns term_goal — surface the return-home suggestion
        // through the same channel as frontier candidates so the external node
        // can decide. We still flip home_return_requested_ to stop suggesting
        // new frontiers from this loop.
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 5000,
          "Exploration: external_selector=true — not auto-issuing return-home goal");
      } else {
        terminalGoalCallbackImpl(home, /*from_user=*/false);
      }
      // Lock in the return-home: the gate at the top of exploreSelectCallback
      // now short-circuits, so a peer's visited_map merge or a fresh frontier
      // pop-up mid-transit can't yank the agent back out. Same flag the
      // /exploration/return_home topic uses — single source of truth.
      home_return_requested_ = true;
      // Internal give-up: once home is reached, re-evaluate so a frontier that
      // reappears/gets discovered during the return still gets explored.
      home_return_resume_on_arrival_ = true;
    }
    exploration_active_ = false;
    // Note: exploration_start_captured_ is intentionally NOT reset here.
    // While the robot is en route home, new frontiers may be discovered as the
    // map updates, causing exploreSelectCallback to auto-restart. Resetting
    // would let it re-capture the mid-transit pose as a new "start" and the
    // final return-home would land there instead of the true session start.
    // Only a manual user goal (true session boundary) resets the flag.
    return;
  }

  // A selectable frontier exists — reset the give-up grace timer.
  no_frontier_since_t_ = -1.0;

  geometry_msgs::msg::PoseStamped g;
  g.header.frame_id    = par_.map_frame_id;
  g.header.stamp       = this->now();
  g.pose.position.x    = next->centroid_xy.x();
  g.pose.position.y    = next->centroid_xy.y();
  g.pose.position.z    = par_.expl_default_goal_z;
  g.pose.orientation.w = 1.0;

  RCLCPP_INFO(this->get_logger(),
              "Exploration: -> frontier %lu at (%.2f, %.2f), state=%d, u=%.3f",
              static_cast<unsigned long>(next->id),
              next->centroid_xy.x(), next->centroid_xy.y(),
              static_cast<int>(next->state), next->cached_utility);

  if (par_.expl_external_selector) {
    // External selector owns term_goal. Skip the internal pursuit bookkeeping
    // entirely so we don't lock onto a frontier the external node might not
    // pick. Frontier markers and current_goal viz still publish elsewhere in
    // this function's lifecycle hooks.
    return;
  }

  terminalGoalCallbackImpl(g, /*from_user=*/false);

  if (!exploration_start_captured_) {
    exploration_start_pos_ = Eigen::Vector3d(cur.pos.x(), cur.pos.y(), cur.pos.z());
    exploration_start_captured_ = true;
    std::cout << "[mighty] Exploration start captured at ("
              << exploration_start_pos_.x() << ", "
              << exploration_start_pos_.y() << ", "
              << exploration_start_pos_.z() << ") — will return here when done"
              << std::endl;
    RCLCPP_INFO(this->get_logger(),
                "Exploration: starting from (%.2f, %.2f, %.2f) — will return here when done",
                exploration_start_pos_.x(), exploration_start_pos_.y(), exploration_start_pos_.z());
  }

  current_explore_id_       = next->id;
  exploration_active_       = true;
  unreachable_consec_count_ = 0;
  frontier_manager_->markSelected(
      next->id, Eigen::Vector2d(robot_pose.x(), robot_pose.y()),
      this->now().seconds());
  publishExplorationCurrentGoal(*next);
}

// ----------------------------------------------------------------------------

namespace {

std_msgs::msg::ColorRGBA makeColor(double r, double g, double b, double a) {
  std_msgs::msg::ColorRGBA c;
  c.r = static_cast<float>(r);
  c.g = static_cast<float>(g);
  c.b = static_cast<float>(b);
  c.a = static_cast<float>(a);
  return c;
}

std_msgs::msg::ColorRGBA colorForState(FrontierState s) {
  switch (s) {
    case FrontierState::ACTIVE:      return makeColor(0.0, 0.8, 1.0, 0.9);  // cyan
    case FrontierState::DORMANT:     return makeColor(0.6, 0.6, 0.6, 0.7);  // gray
    case FrontierState::VISITED:     return makeColor(0.0, 0.9, 0.0, 0.7);  // green
    case FrontierState::INVALIDATED: return makeColor(0.9, 0.0, 0.0, 0.7);  // red
  }
  return makeColor(1.0, 1.0, 1.0, 0.7);
}

}  // namespace

/**
 * @brief Publish a structured FrontierList for external consumers (e.g. the
 *        VLM goal selector). Mirrors the records the visualization shows, but
 *        with a stable schema and full per-frontier metadata. Only ACTIVE and
 *        DORMANT records are included — terminal states stay in the DB but
 *        aren't candidates for an external selector.
 */
void MIGHTY_NODE::publishFrontierData() {
  if (!pub_frontier_data_ || !frontier_manager_) return;

  dynus_interfaces::msg::FrontierList msg;
  msg.header.frame_id = par_.map_frame_id;
  msg.header.stamp    = this->now();

  for (const auto& r : frontier_manager_->records()) {
    if (r.state != FrontierState::ACTIVE && r.state != FrontierState::DORMANT) {
      continue;
    }
    dynus_interfaces::msg::Frontier f;
    f.id           = r.id;
    f.state        = (r.state == FrontierState::ACTIVE)
                       ? dynus_interfaces::msg::Frontier::STATE_ACTIVE
                       : dynus_interfaces::msg::Frontier::STATE_DORMANT;
    f.centroid.x   = r.centroid_xy.x();
    f.centroid.y   = r.centroid_xy.y();
    f.centroid.z   = par_.expl_default_goal_z;
    f.size_cells   = static_cast<uint64_t>(std::max(0, r.size_cells));
    const double cell_m = (occ_grid_2d_ ? occ_grid_2d_->resolution() : 0.0);
    f.size_m2      = static_cast<double>(f.size_cells) * cell_m * cell_m;
    f.aabb_min.x   = r.aabb_min.x();
    f.aabb_min.y   = r.aabb_min.y();
    f.aabb_max.x   = r.aabb_max.x();
    f.aabb_max.y   = r.aabb_max.y();
    f.cached_utility = r.cached_utility;
    msg.frontiers.push_back(f);
  }

  pub_frontier_data_->publish(msg);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish frontier visualization markers. We only show ACTIVE and
 *        DORMANT centroids plus a small `id=N` text label per centroid, plus
 *        the yellow robot→goal line. VISITED and INVALIDATED frontiers are
 *        kept in the database (so they still gate selection and DB eviction)
 *        but suppressed from RViz to keep the scene readable. One MarkerArray
 *        per cycle, prefixed with DELETEALL so stale markers don't accumulate.
 */
void MIGHTY_NODE::publishFrontierMarkers() {
  if (!pub_frontiers_ || !frontier_manager_ || !occ_grid_2d_) return;

  visualization_msgs::msg::MarkerArray arr;

  // DELETEALL prefix to wipe stale markers from previous cycles.
  {
    visualization_msgs::msg::Marker del;
    del.header.frame_id = par_.map_frame_id;
    del.header.stamp    = this->now();
    del.action          = visualization_msgs::msg::Marker::DELETEALL;
    arr.markers.push_back(del);
  }

  const auto& records = frontier_manager_->records();

  visualization_msgs::msg::Marker centroids;
  centroids.header.frame_id = par_.map_frame_id;
  centroids.header.stamp    = this->now();
  centroids.ns              = "frontier_centroids";
  centroids.id              = 0;
  centroids.type            = visualization_msgs::msg::Marker::SPHERE_LIST;
  centroids.action          = visualization_msgs::msg::Marker::ADD;
  centroids.scale.x = 0.3;
  centroids.scale.y = 0.3;
  centroids.scale.z = 0.3;
  centroids.pose.orientation.w = 1.0;

  int label_id = 0;
  for (const auto& r : records) {
    // Only ACTIVE and DORMANT are visualized — VISITED/INVALIDATED stay in
    // the DB but are hidden from RViz.
    if (r.state != FrontierState::ACTIVE && r.state != FrontierState::DORMANT) {
      continue;
    }

    // Centroid sphere, colored by state.
    geometry_msgs::msg::Point p;
    p.x = r.centroid_xy.x();
    p.y = r.centroid_xy.y();
    p.z = par_.expl_default_goal_z + 0.1;
    centroids.points.push_back(p);
    centroids.colors.push_back(colorForState(r.state));

    // Per-record text label — just the id, so you can track a specific
    // frontier across cycles. (Cluster size and cached utility are dropped.)
    visualization_msgs::msg::Marker label;
    label.header.frame_id = par_.map_frame_id;
    label.header.stamp    = this->now();
    label.ns              = "frontier_labels";
    label.id              = label_id++;
    label.type            = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
    label.action          = visualization_msgs::msg::Marker::ADD;
    label.pose.position.x = r.centroid_xy.x();
    label.pose.position.y = r.centroid_xy.y();
    label.pose.position.z = par_.expl_default_goal_z + 0.6;
    label.pose.orientation.w = 1.0;
    label.scale.z = 0.4;
    label.color   = makeColor(0.0, 0.0, 0.0, 1.0);
    {
      char buf[64];
      if (r.pursuit_deadline_t > 0.0 && r.pursuit_budget_sec > 0.0) {
        const double t_now = this->now().seconds();
        const double elapsed = r.pursuit_budget_sec
                             - std::max(0.0, r.pursuit_deadline_t - t_now);
        std::snprintf(buf, sizeof(buf), "id=%lu\n%.1fs/%.1fs",
                      static_cast<unsigned long>(r.id),
                      elapsed, r.pursuit_budget_sec);
      } else {
        std::snprintf(buf, sizeof(buf), "id=%lu",
                      static_cast<unsigned long>(r.id));
      }
      label.text = buf;
    }
    arr.markers.push_back(label);
  }
  if (!centroids.points.empty()) arr.markers.push_back(centroids);

  // Yellow line from robot to current exploration goal.
  if (exploration_active_ && state_initialized_) {
    state cur;
    mighty_ptr_->getState(cur);
    auto* r = frontier_manager_->find(current_explore_id_);
    if (r) {
      visualization_msgs::msg::Marker line;
      line.header.frame_id = par_.map_frame_id;
      line.header.stamp    = this->now();
      line.ns              = "frontier_goal_line";
      line.id              = 0;
      line.type            = visualization_msgs::msg::Marker::LINE_STRIP;
      line.action          = visualization_msgs::msg::Marker::ADD;
      line.scale.x         = 0.15;
      line.color           = makeColor(1.0, 1.0, 0.0, 0.9);
      line.pose.orientation.w = 1.0;
      geometry_msgs::msg::Point a;
      a.x = cur.pos.x();
      a.y = cur.pos.y();
      a.z = par_.expl_default_goal_z + 0.1;
      geometry_msgs::msg::Point b;
      b.x = r->centroid_xy.x();
      b.y = r->centroid_xy.y();
      b.z = par_.expl_default_goal_z + 0.1;
      line.points.push_back(a);
      line.points.push_back(b);
      arr.markers.push_back(line);
    }
  }

  // Yellow rectangle showing the user-configured exploration bounds, plus a
  // text label so it's obvious in RViz what the rectangle means.
  if (par_.expl_bounds_enabled) {
    const double z = par_.expl_default_goal_z + 0.05;
    const double x0 = par_.expl_bounds_min_x;
    const double x1 = par_.expl_bounds_max_x;
    const double y0 = par_.expl_bounds_min_y;
    const double y1 = par_.expl_bounds_max_y;

    visualization_msgs::msg::Marker rect;
    rect.header.frame_id = par_.map_frame_id;
    rect.header.stamp    = this->now();
    rect.ns              = "exploration_bounds";
    rect.id              = 0;
    rect.type            = visualization_msgs::msg::Marker::LINE_STRIP;
    rect.action          = visualization_msgs::msg::Marker::ADD;
    rect.scale.x         = 0.10;            // line thickness
    rect.color           = makeColor(1.0, 1.0, 0.0, 0.9);  // yellow
    rect.pose.orientation.w = 1.0;
    auto pt = [&](double x, double y) {
      geometry_msgs::msg::Point p; p.x = x; p.y = y; p.z = z; return p;
    };
    rect.points.push_back(pt(x0, y0));
    rect.points.push_back(pt(x1, y0));
    rect.points.push_back(pt(x1, y1));
    rect.points.push_back(pt(x0, y1));
    rect.points.push_back(pt(x0, y0));   // close the loop
    arr.markers.push_back(rect);

    visualization_msgs::msg::Marker label;
    label.header.frame_id = par_.map_frame_id;
    label.header.stamp    = this->now();
    label.ns              = "exploration_bounds_label";
    label.id              = 0;
    label.type            = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
    label.action          = visualization_msgs::msg::Marker::ADD;
    label.pose.position.x = 0.5 * (x0 + x1);
    label.pose.position.y = y1 + 0.5;     // sit just above the top edge
    label.pose.position.z = z + 0.5;
    label.pose.orientation.w = 1.0;
    label.scale.z = 0.6;
    label.color   = makeColor(1.0, 1.0, 0.0, 0.9);
    label.text    = "Exploration Area";
    arr.markers.push_back(label);
  }

  pub_frontiers_->publish(arr);
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::publishExplorationCurrentGoal(const FrontierRecord& r) {
  if (!pub_explore_current_goal_) return;
  geometry_msgs::msg::PoseStamped g;
  g.header.frame_id    = par_.map_frame_id;
  g.header.stamp       = this->now();
  g.pose.position.x    = r.centroid_xy.x();
  g.pose.position.y    = r.centroid_xy.y();
  g.pose.position.z    = par_.expl_default_goal_z;
  g.pose.orientation.w = 1.0;
  pub_explore_current_goal_->publish(g);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the persistent occupancy map as a nav_msgs/OccupancyGrid.
 *        The buffer already holds tristate values (-1 unknown / 0 free /
 *        100 occupied) matching the message encoding, so this is a flat
 *        copy. RViz subscribes to /exploration/visited_map and renders it
 *        behind the sliding occ_2d so revisited cells keep their last-known
 *        FREE/OCCUPIED color instead of flickering UNKNOWN.
 */
void MIGHTY_NODE::publishVisitedMap() {
  if (!pub_visited_map_ || !visited_map_ || visited_map_->empty()) return;

  nav_msgs::msg::OccupancyGrid msg;
  msg.header.frame_id = par_.map_frame_id;
  msg.header.stamp    = this->now();
  msg.info.resolution = static_cast<float>(visited_map_->resolution());
  msg.info.width      = static_cast<unsigned>(visited_map_->width());
  msg.info.height     = static_cast<unsigned>(visited_map_->height());
  msg.info.origin.position.x    = visited_map_->originX();
  msg.info.origin.position.y    = visited_map_->originY();
  // Match whatever z the live occ_2d layer is at (global_mapper uses
  // z_ground). Falls back to the exploration default goal z until the first
  // occ_2d arrives so a late RViz subscriber still sees a sensible plane.
  msg.info.origin.position.z    =
      occ2d_origin_z_.value_or(par_.expl_default_goal_z);
  msg.info.origin.orientation.w = 1.0;

  const auto& v = visited_map_->data();
  msg.data.assign(v.begin(), v.end());  // raw tristate copy
  pub_visited_map_->publish(msg);
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::publishHeatCloud() {
  if (!pub_heat_cloud_) return;

  auto map_util = mighty_ptr_->getMapUtil();
  if (!map_util) return;

  // Check if any heat source is enabled (terrain cost counts as heat in 2D mode)
  const bool has_heat =
      par_.dynamic_heat_enabled || par_.static_heat_enabled || par_.use_2d_planning;
  if (!par_.use_heat_map || !has_heat) return;

  // -------- Tunables --------
  const int stride = 1;              // 1 = every voxel, 2 = every 2 voxels, etc.
  const float heat_min = 0.001f;     // only publish voxels with heat >= this
  const size_t max_points = 200000;  // hard cap for safety
  // --------------------------

  const auto dim = map_util->getDim();  // Veci<3>
  const int nx = dim(0);
  const int ny = dim(1);
  const int nz = dim(2);

  // Use Eigen-aligned vector type
  vec_Vec3f pts;
  std::vector<float> intens;
  pts.reserve(50000);
  intens.reserve(50000);

  for (int x = 0; x < nx; x += stride) {
    for (int y = 0; y < ny; y += stride) {
      for (int z = 0; z < nz; z += stride) {
        const float h = map_util->getHeat(x, y, z);
        if (h < heat_min) continue;

        const Vec3f p = map_util->intToFloat(Veci<3>(x, y, z));
        pts.push_back(p);
        intens.push_back(h);

        if (pts.size() >= max_points) goto BUILD_MSG;
      }
    }
  }

BUILD_MSG:
  // Normalize intensities to [0, 1] so the color gradient is visible in RViz
  float max_intensity = 0.0f;
  for (const float v : intens) max_intensity = std::max(max_intensity, v);

  if (max_intensity > 0.0f) {
    const float inv_max = 1.0f / max_intensity;
    for (float& v : intens) v *= inv_max;
  }

  sensor_msgs::msg::PointCloud2 msg;
  msg.header.frame_id = par_.map_frame_id;
  msg.header.stamp = this->now();

  sensor_msgs::PointCloud2Modifier modifier(msg);
  modifier.setPointCloud2Fields(4, "x", 1, sensor_msgs::msg::PointField::FLOAT32, "y", 1,
                                sensor_msgs::msg::PointField::FLOAT32, "z", 1,
                                sensor_msgs::msg::PointField::FLOAT32, "intensity", 1,
                                sensor_msgs::msg::PointField::FLOAT32);
  modifier.resize(pts.size());

  sensor_msgs::PointCloud2Iterator<float> iter_x(msg, "x");
  sensor_msgs::PointCloud2Iterator<float> iter_y(msg, "y");
  sensor_msgs::PointCloud2Iterator<float> iter_z(msg, "z");
  sensor_msgs::PointCloud2Iterator<float> iter_i(msg, "intensity");

  for (size_t k = 0; k < pts.size(); ++k, ++iter_x, ++iter_y, ++iter_z, ++iter_i) {
    const auto& p = pts[k];
    *iter_x = static_cast<float>(p(0));
    *iter_y = static_cast<float>(p(1));
    *iter_z = static_cast<float>(p(2));
    *iter_i = intens[k];
  }

  pub_heat_cloud_->publish(msg);
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::publishGround2DOccupied() {
  if (!pub_ground_2d_occ_) return;

  auto map_util = mighty_ptr_->getMapUtil();
  if (!map_util) {
    RCLCPP_WARN_ONCE(this->get_logger(), "publishGround2DOccupied: map_util is null");
    return;
  }
  if (!map_util->has2DMap()) {
    RCLCPP_WARN_ONCE(this->get_logger(), "publishGround2DOccupied: has2DMap()=false");
    return;
  }
  RCLCPP_INFO_ONCE(this->get_logger(), "publishGround2DOccupied: publishing 2D occupied map");

  int dimX, dimY;
  map_util->get2DDimensions(dimX, dimY);
  const auto origin = map_util->getOrigin();
  const float res = static_cast<float>(map_util->getRes());
  const float goal_z = static_cast<float>(par_.default_goal_z);

  // Collect 2D occupied cells as 3D points projected to default_goal_z
  pcl::PointCloud<pcl::PointXYZI> cloud;
  for (int x = 0; x < dimX; ++x) {
    for (int y = 0; y < dimY; ++y) {
      if (map_util->get2DOccupancy(x, y) != 0) {  // occupied
        pcl::PointXYZI pt;
        pt.x = origin(0) + (x + 0.5f) * res;
        pt.y = origin(1) + (y + 0.5f) * res;
        pt.z = goal_z;
        pt.intensity = 1.0f;
        cloud.push_back(pt);
      }
    }
  }

  sensor_msgs::msg::PointCloud2 msg;
  pcl::toROSMsg(cloud, msg);
  msg.header.frame_id = par_.map_frame_id;
  msg.header.stamp = this->now();
  pub_ground_2d_occ_->publish(msg);
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::publishGround2DHeat() {
  if (!pub_ground_2d_heat_ || !par_.use_2d_planning) return;

  auto map_util = mighty_ptr_->getMapUtil();
  if (!map_util || !map_util->has2DMap()) return;

  int dimX, dimY;
  map_util->get2DDimensions(dimX, dimY);
  const auto origin = map_util->getOrigin();
  const float res = static_cast<float>(map_util->getRes());

  // Collect 2D heat (dynamic + static + terrain) as colored point cloud
  pcl::PointCloud<pcl::PointXYZI> cloud;
  for (int x = 0; x < dimX; ++x) {
    for (int y = 0; y < dimY; ++y) {
      const float h = map_util->getHeat2D(x, y);
      const float tc = map_util->getTerrainCost(x, y);
      const float total = std::max(h, tc);
      if (total > 0.001f) {
        pcl::PointXYZI pt;
        pt.x = origin(0) + (x + 0.5f) * res;
        pt.y = origin(1) + (y + 0.5f) * res;
        pt.z = 0.05f;  // slightly above ground for visibility
        pt.intensity = total;
        cloud.push_back(pt);
      }
    }
  }

  sensor_msgs::msg::PointCloud2 msg;
  pcl::toROSMsg(cloud, msg);
  msg.header.frame_id = par_.map_frame_id;
  msg.header.stamp = this->now();
  pub_ground_2d_heat_->publish(msg);
}

}  // namespace mighty

// ----------------------------------------------------------------------------

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);

  // Initialize multi-threaded executor
  rclcpp::executors::MultiThreadedExecutor executor;

  // add node to executor
  auto node = std::make_shared<mighty::MIGHTY_NODE>();
  executor.add_node(node);

  // spin
  executor.spin();

  rclcpp::shutdown();
  return 0;
}