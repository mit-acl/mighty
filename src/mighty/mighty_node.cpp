/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <mighty/mighty_node.hpp>

// ----------------------------------------------------------------------------

/**
 * @brief Constructor
 */
MIGHTY_NODE::MIGHTY_NODE() : Node("mighty_node")
{

  // Get id from ns
  ns_ = this->get_namespace();
  ns_ = ns_.substr(ns_.find_last_of("/") + 1);
  id_str_ = ns_.substr(ns_.size() - 2); // ns is like NX01, so we get the last two characters and convert to int
  id_ = std::stoi(id_str_);

  // Declare, set, and print parameters
  this->declareParameters();
  this->setParameters();
  this->printParameters();

  // Qos policy settings
  rclcpp::QoS critical_qos(rclcpp::KeepLast(10));
  critical_qos.reliable().durability_volatile();

  // Create callbackgroup
  this->cb_group_mu_1_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_2_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_3_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_4_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_5_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_6_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_7_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_8_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_mu_9_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_re_1_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_2_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_3_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_4_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_5_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_6_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_7_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_8_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_re_9_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);
  this->cb_group_map_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_replan_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
  this->cb_group_goal_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);

  // Options for callback group
  rclcpp::SubscriptionOptions options_re_1;
  options_re_1.callback_group = this->cb_group_re_1_;
  rclcpp::SubscriptionOptions options_re_2;
  options_re_2.callback_group = this->cb_group_re_2_;
  rclcpp::SubscriptionOptions options_map;
  options_map.callback_group = this->cb_group_map_;

  // Visulaization publishers
  pub_dynamic_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>("dynamic_occupied_grid", 10);                                              // visual level 2 (no longer used)
  pub_static_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("static_map_marker", 10);                                     // visual level 2
  pub_dynamic_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("dynamic_map_marker", 10);                                   // visual level 2
  pub_free_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("free_map_marker", 10);                                         // visual level 2
  pub_unknown_map_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("unknown_map_marker", 10);                                   // visual level 2
  pub_free_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>("free_grid", 10);                                                             // visual level 2 (no longer used)
  pub_unknown_map_ = this->create_publisher<sensor_msgs::msg::PointCloud2>("unknown_grid", 10);                                                       // visual level 2 (no longer used)
  pub_dgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("dgp_path_marker", 10);                                         // visual level 1
  pub_original_dgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("original_dgp_path_marker", 10);                       // visual level 1
  pub_free_dgp_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("free_dgp_path_marker", 10);                               // visual level 1
  pub_local_global_path_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("local_global_path_marker", 10);                       // visual level 1
  pub_local_global_path_after_push_marker_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("local_global_path_after_push_marker", 10); // visual level 1
  pub_poly_whole_ = this->create_publisher<decomp_ros_msgs::msg::PolyhedronArray>("poly_whole", 10);                                                  // visual level 1
  pub_poly_safe_ = this->create_publisher<decomp_ros_msgs::msg::PolyhedronArray>("poly_safe", 10);                                                    // visual level 1
  pub_traj_committed_colored_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("traj_committed_colored", 10);                           // visual level 1
  pub_traj_subopt_colored_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("traj_subopt_colored", 10);                                 // visual level 1
  pub_setpoint_ = this->create_publisher<geometry_msgs::msg::PointStamped>("setpoint_vis", 10);                                                       // visual level 1
  pub_actual_traj_ = this->create_publisher<visualization_msgs::msg::Marker>("actual_traj", 10);                                                      // visual level 1
  pub_fov_ = this->create_publisher<visualization_msgs::msg::Marker>("fov", 10);                                                                      // visual level 1
  pub_cp_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("cp", 10);                                                                   // visual level 1
  pub_static_push_points_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("static_push_points", 10);                                   // visual level 1
  pub_p_points_ = this->create_publisher<visualization_msgs::msg::MarkerArray>("p_points", 10);                                                       // visual level 1
  pub_point_A_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_A", 10);                                                             // visual level 1
  pub_point_G_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_G", 10);                                                             // visual level 1
  pub_point_E_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_E", 10);                                                             // visual level 1
  pub_point_G_term_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_G_term", 10);                                                   // visual level 1
  pub_current_state_ = this->create_publisher<geometry_msgs::msg::PointStamped>("point_current_state", 10);                                           // visual level 1
  pub_vel_text_ = this->create_publisher<visualization_msgs::msg::Marker>("vel_text", 10);                                                            // visual level 1

  // Debug publishers
  pub_yaw_output_ = this->create_publisher<dynus_interfaces::msg::YawOutput>("yaw_output", 10);

  // Essential publishers
  pub_own_traj_ = this->create_publisher<dynus_interfaces::msg::DynTraj>("/trajs", critical_qos);
  pub_goal_ = this->create_publisher<dynus_interfaces::msg::Goal>("goal", critical_qos);
  pub_goal_reached_ = this->create_publisher<std_msgs::msg::Empty>("goal_reached", critical_qos);

  // Subscribers
  sub_traj_ = this->create_subscription<dynus_interfaces::msg::DynTraj>("/trajs", critical_qos, std::bind(&MIGHTY_NODE::trajCallback, this, std::placeholders::_1), options_re_1);
  sub_predicted_traj_ = this->create_subscription<dynus_interfaces::msg::DynTraj>("predicted_trajs", critical_qos, std::bind(&MIGHTY_NODE::trajCallback, this, std::placeholders::_1), options_re_1);
  sub_state_ = this->create_subscription<dynus_interfaces::msg::State>("state", critical_qos, std::bind(&MIGHTY_NODE::stateCallback, this, std::placeholders::_1), options_re_1);
  sub_terminal_goal_ = this->create_subscription<geometry_msgs::msg::PoseStamped>("term_goal", critical_qos, std::bind(&MIGHTY_NODE::terminalGoalCallback, this, std::placeholders::_1));

  // Timer for callback
  timer_replanning_ = this->create_wall_timer(10ms, std::bind(&MIGHTY_NODE::replanCallback, this), this->cb_group_replan_);
  timer_goal_ = this->create_wall_timer(std::chrono::duration<double>(par_.dc), std::bind(&MIGHTY_NODE::publishGoal, this), this->cb_group_goal_);
  if (use_benchmark_)
    timer_goal_reached_check_ = this->create_wall_timer(100ms, std::bind(&MIGHTY_NODE::goalReachedCheckCallback, this), this->cb_group_re_3_);
  timer_cleanup_old_trajs_ = this->create_wall_timer(500ms, std::bind(&MIGHTY_NODE::cleanUpOldTrajsCallback, this), this->cb_group_mu_5_);
  if (par_.use_hardware)
    timer_initial_pose_ = this->create_wall_timer(100ms, std::bind(&MIGHTY_NODE::getInitialPoseHwCallback, this), this->cb_group_mu_9_);

  // Stop the timer for callback
  if (timer_replanning_)
    timer_replanning_->cancel();
  if (timer_goal_)
    timer_goal_->cancel();
  if (!use_benchmark_ && timer_goal_reached_check_)
    timer_goal_reached_check_->cancel();

  // Initialize the DYNUS object
  mighty_ptr_ = std::make_shared<MIGHTY>(par_);

  // Initialize the tf2 buffer and listener
  tf2_buffer_ = std::make_shared<tf2_ros::Buffer>(this->get_clock());
  tf2_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf2_buffer_);
  auto timer_interface = std::make_shared<tf2_ros::CreateTimerROS>(this->get_node_base_interface(), this->get_node_timers_interface());
  tf2_buffer_->setCreateTimerInterface(timer_interface);

  // Initialize the d435 depth frame ID and camera
  d435_depth_frame_id_ = ns_ + "/" + ns_ + "_d435_depth_optical_frame";
  lidar_frame_id_ = ns_ + "/NX01_livox";

  // Construct the FOV marker
  constructFOVMarker();

  // Initialize the initial pose topic name
  initial_pose_topic_ = ns_ + "/init_pose";

  // Synchronize the occupancy grid and unknown grid
  occup_grid_sub_.subscribe(this, "occupancy_grid", rmw_qos_profile_sensor_data, options_map);
  unknown_grid_sub_.subscribe(this, "unknown_grid", rmw_qos_profile_sensor_data, options_map);
  sync_.reset(new Sync(MySyncPolicy(10), occup_grid_sub_, unknown_grid_sub_));
  sync_->registerCallback(std::bind(&MIGHTY_NODE::mapCallback, this, std::placeholders::_1, std::placeholders::_2));
}

// ----------------------------------------------------------------------------

/**
 * @brief Destructor
 */
MIGHTY_NODE::~MIGHTY_NODE()
{
  // release the memory
  mighty_ptr_.reset();
}

// ----------------------------------------------------------------------------

/**
 * @brief Declare the parameters
 */
void MIGHTY_NODE::declareParameters()
{
  // UAV or Ground robot
  this->declare_parameter("vehicle_type", "uav");
  this->declare_parameter("provide_goal_in_global_frame", false);
  this->declare_parameter("use_hardware", false);

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
  this->declare_parameter("global_planner_huristic_weight", 1.0);
  this->declare_parameter("factor_dgp", 1.0);
  this->declare_parameter("inflation_dgp", 0.5);
  this->declare_parameter("x_min", -100.0);
  this->declare_parameter("x_max", 100.0);
  this->declare_parameter("y_min", -100.0);
  this->declare_parameter("y_max", 100.0);
  this->declare_parameter("z_min", 0.0);
  this->declare_parameter("z_max", 5.0);
  this->declare_parameter("dgp_timeout_duration_ms", 1000);
  this->declare_parameter("use_free_start", false);
  this->declare_parameter("free_start_factor", 1.0);
  this->declare_parameter("use_free_goal", false);
  this->declare_parameter("free_goal_factor", 1.0);
  this->declare_parameter("num_N", 3);
  this->declare_parameter("max_dist_vertexes", 5.0);
  this->declare_parameter("w_unknown", 1.0);
  this->declare_parameter("w_align", 60.0);
  this->declare_parameter("decay_len_cells", 20.0);
  this->declare_parameter("w_side", 0.2);

  // LOS post processing parameters
  this->declare_parameter("los_cells", 3);
  this->declare_parameter("min_len", 0.5);   // [m] minimum length between two waypoints after post processing
  this->declare_parameter("min_turn", 10.0); // [deg] minimum turn angle after post processing

  // Path push visualization parameters
  this->declare_parameter("use_state_update", true);
  this->declare_parameter("use_random_color_for_global_path", false);
  this->declare_parameter("use_path_push_for_visualization", false);

  // Decomposition parameters
  this->declare_parameter("local_box_size", std::vector<float>{2.0, 2.0, 2.0});
  this->declare_parameter("min_dist_from_agent_to_traj", 6.0);
  this->declare_parameter("use_shrinked_box", false);
  this->declare_parameter("shrinked_box_size", 0.2);

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

  // Initial guess parameters
  this->declare_parameter("use_multiple_initial_guesses", true);
  this->declare_parameter("num_perturbation_for_ig", 8);
  this->declare_parameter("r_max_for_ig", 1.0);

  // Optimization parameters
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

  // Debug flag
  this->declare_parameter("debug_verbose", false);
}

// ----------------------------------------------------------------------------

/**
 * @brief Set the parameters
 */
void MIGHTY_NODE::setParameters()
{
  // Set the parameters

  // Vehicle type (UAV, Wheeled Robit, or Quadruped)
  par_.vehicle_type = this->get_parameter("vehicle_type").as_string();
  par_.provide_goal_in_global_frame = this->get_parameter("provide_goal_in_global_frame").as_bool();
  par_.use_hardware = this->get_parameter("use_hardware").as_bool();

  // Flight mode
  par_.flight_mode = this->get_parameter("flight_mode").as_string();

  // Visual level
  par_.visual_level = this->get_parameter("visual_level").as_int();

  // Global Planner parameters
  file_path_ = this->get_parameter("file_path").as_string();
  use_benchmark_ = this->get_parameter("use_benchmark").as_bool();
  par_.global_planner = this->get_parameter("global_planner").as_string();
  par_.global_planner_verbose = this->get_parameter("global_planner_verbose").as_bool();
  par_.global_planner_huristic_weight = this->get_parameter("global_planner_huristic_weight").as_double();
  par_.factor_dgp = this->get_parameter("factor_dgp").as_double();
  par_.inflation_dgp = this->get_parameter("inflation_dgp").as_double();
  par_.x_min = this->get_parameter("x_min").as_double();
  par_.x_max = this->get_parameter("x_max").as_double();
  par_.y_min = this->get_parameter("y_min").as_double();
  par_.y_max = this->get_parameter("y_max").as_double();
  par_.z_min = this->get_parameter("z_min").as_double();
  par_.z_max = this->get_parameter("z_max").as_double();
  par_.dgp_timeout_duration_ms = this->get_parameter("dgp_timeout_duration_ms").as_int();
  par_.use_free_start = this->get_parameter("use_free_start").as_bool();
  par_.free_start_factor = this->get_parameter("free_start_factor").as_double();
  par_.use_free_goal = this->get_parameter("use_free_goal").as_bool();
  par_.free_goal_factor = this->get_parameter("free_goal_factor").as_double();
  par_.num_N = this->get_parameter("num_N").as_int();
  par_.max_dist_vertexes = this->get_parameter("max_dist_vertexes").as_double();
  par_.w_unknown = this->get_parameter("w_unknown").as_double();
  par_.w_align = this->get_parameter("w_align").as_double();
  par_.decay_len_cells = this->get_parameter("decay_len_cells").as_double();
  par_.w_side = this->get_parameter("w_side").as_double();

  // LOS post processing parameters
  par_.los_cells = this->get_parameter("los_cells").as_int();
  par_.min_len = this->get_parameter("min_len").as_double();
  par_.min_turn = this->get_parameter("min_turn").as_double();

  // Path push visualization parameters
  par_.use_state_update = this->get_parameter("use_state_update").as_bool();
  par_.use_random_color_for_global_path = this->get_parameter("use_random_color_for_global_path").as_bool();
  par_.use_path_push_for_visualization = this->get_parameter("use_path_push_for_visualization").as_bool();

  // Static obstacle push parameters

  // Decomposition parameters
  par_.local_box_size = this->get_parameter("local_box_size").as_double_array();
  par_.min_dist_from_agent_to_traj = this->get_parameter("min_dist_from_agent_to_traj").as_double();
  par_.use_shrinked_box = this->get_parameter("use_shrinked_box").as_bool();
  par_.shrinked_box_size = this->get_parameter("shrinked_box_size").as_double();

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

  // Initial guess parameters
  par_.use_multiple_initial_guesses = this->get_parameter("use_multiple_initial_guesses").as_bool();
  par_.num_perturbation_for_ig = this->get_parameter("num_perturbation_for_ig").as_int();
  par_.r_max_for_ig = this->get_parameter("r_max_for_ig").as_double();

  // Optimization parameters
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

  if (par_.default_goal_z <= par_.z_min)
  {
    RCLCPP_ERROR(this->get_logger(), "Default goal z is lower than the ground level");
  }

  if (par_.default_goal_z >= par_.z_max)
  {
    RCLCPP_ERROR(this->get_logger(), "Default goal z is higher than the max level");
  }

  // Debug flag
  par_.debug_verbose = this->get_parameter("debug_verbose").as_bool();
}

// ----------------------------------------------------------------------------

/**
 * @brief Print the parameters
 */
void MIGHTY_NODE::printParameters()
{
  // Print the parameters

  // Vehicle type (UAV, Wheeled Robit, or Quadruped)
  RCLCPP_INFO(this->get_logger(), "Vehicle Type: %d", par_.vehicle_type);
  RCLCPP_INFO(this->get_logger(), "Provide Goal in Global Frame: %d", par_.provide_goal_in_global_frame);
  RCLCPP_INFO(this->get_logger(), "Use Hardware: %d", par_.use_hardware);

  // Flight mode
  RCLCPP_INFO(this->get_logger(), "Flight Mode: %s", par_.flight_mode.c_str());

  // Visual
  RCLCPP_INFO(this->get_logger(), "Visual Level: %d", par_.visual_level);

  // DGP parameters
  RCLCPP_INFO(this->get_logger(), "File Path: %s", file_path_.c_str());
  RCLCPP_INFO(this->get_logger(), "Perform Benchmark?: %d", use_benchmark_);
  RCLCPP_INFO(this->get_logger(), "Initial Guess Planner: %s", par_.global_planner.c_str());
  RCLCPP_INFO(this->get_logger(), "DGP Planner Verbose: %d", par_.global_planner_verbose);
  RCLCPP_INFO(this->get_logger(), "Global Planner Huristic Weight: %f", par_.global_planner_huristic_weight);
  RCLCPP_INFO(this->get_logger(), "Factor DGP: %f", par_.factor_dgp);
  RCLCPP_INFO(this->get_logger(), "Inflation DGP: %f", par_.inflation_dgp);
  RCLCPP_INFO(this->get_logger(), "X Min: %f", par_.x_min);
  RCLCPP_INFO(this->get_logger(), "X Max: %f", par_.x_max);
  RCLCPP_INFO(this->get_logger(), "Y Min: %f", par_.y_min);
  RCLCPP_INFO(this->get_logger(), "Y Max: %f", par_.y_max);
  RCLCPP_INFO(this->get_logger(), "Z Ground: %f", par_.z_min);
  RCLCPP_INFO(this->get_logger(), "Z Max: %f", par_.z_max);
  RCLCPP_INFO(this->get_logger(), "DGP Timeout Duration: %d", par_.dgp_timeout_duration_ms);
  RCLCPP_INFO(this->get_logger(), "Use Free Start?: %d", par_.use_free_start);
  RCLCPP_INFO(this->get_logger(), "Free Start Factor: %f", par_.free_start_factor);
  RCLCPP_INFO(this->get_logger(), "Use Free Goal?: %d", par_.use_free_goal);
  RCLCPP_INFO(this->get_logger(), "Free Goal Factor: %f", par_.free_goal_factor);
  RCLCPP_INFO(this->get_logger(), "Num N: %d", par_.num_N);
  RCLCPP_INFO(this->get_logger(), "max_dist_vertexes: %f", par_.max_dist_vertexes);
  RCLCPP_INFO(this->get_logger(), "w_unknown: %f", par_.w_unknown);
  RCLCPP_INFO(this->get_logger(), "w_align: %f", par_.w_align);
  RCLCPP_INFO(this->get_logger(), "decay_len_cells: %f", par_.decay_len_cells);
  RCLCPP_INFO(this->get_logger(), "w_side: %f", par_.w_side);

  // LOS post processing parameters
  RCLCPP_INFO(this->get_logger(), "LOS Cells: %d", par_.los_cells);
  RCLCPP_INFO(this->get_logger(), "Min Len: %f", par_.min_len);
  RCLCPP_INFO(this->get_logger(), "Min Turn: %f", par_.min_turn);

  // Path push visualization parameters
  RCLCPP_INFO(this->get_logger(), "Use State Update?: %d", par_.use_state_update);
  RCLCPP_INFO(this->get_logger(), "Use Random Color for Global Path?: %d", par_.use_random_color_for_global_path);
  RCLCPP_INFO(this->get_logger(), "Use Path Push for Paper?: %d", par_.use_path_push_for_visualization);

  // Static obstacle push parameters

  RCLCPP_INFO(this->get_logger(), "Local Box Size: (%f, %f, %f)", par_.local_box_size[0], par_.local_box_size[1], par_.local_box_size[2]);
  RCLCPP_INFO(this->get_logger(), "Min Dist from Agent to Traj: %f", par_.min_dist_from_agent_to_traj);
  RCLCPP_INFO(this->get_logger(), "Use Shrinked Box: %d", par_.use_shrinked_box);
  RCLCPP_INFO(this->get_logger(), "Shrinked Box Size: %f", par_.shrinked_box_size);

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
  RCLCPP_INFO(this->get_logger(), "Comm Delay Inflation Alpha: %f", par_.comm_delay_inflation_alpha);
  RCLCPP_INFO(this->get_logger(), "Comm Delay Inflation Max: %f", par_.comm_delay_inflation_max);
  RCLCPP_INFO(this->get_logger(), "Comm Delay Filter Alpha: %f", par_.comm_delay_filter_alpha);

  // Simulation parameters
  RCLCPP_INFO(this->get_logger(), "D435 Depth Max: %f", par_.depth_camera_depth_max);
  RCLCPP_INFO(this->get_logger(), "FOV Visual Depth: %f", par_.fov_visual_depth);
  RCLCPP_INFO(this->get_logger(), "FOV Visual X Deg: %f", par_.fov_visual_x_deg);
  RCLCPP_INFO(this->get_logger(), "FOV Visual Y Deg: %f", par_.fov_visual_y_deg);

  // Initial guess parameters
  RCLCPP_INFO(this->get_logger(), "Use Multiple Initial Guesses: %d", par_.use_multiple_initial_guesses);
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
  RCLCPP_INFO(this->get_logger(), "Violation of Bodyrate Constr. Weight: %f", par_.dyn_constr_bodyrate_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Tilt Constr. Weight: %f", par_.dyn_constr_tilt_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Thrust Constr. Weight: %f", par_.dyn_constr_thrust_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Vel Constr. Weight: %f", par_.dyn_constr_vel_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Aeccel Constr. Weight: %f", par_.dyn_constr_acc_weight);
  RCLCPP_INFO(this->get_logger(), "Violation of Jerk Constr. Weight: %f", par_.dyn_constr_jerk_weight);
  RCLCPP_INFO(this->get_logger(), "Num Dynamic Obstacles Samples: %d", par_.num_dyn_obst_samples);
  RCLCPP_INFO(this->get_logger(), "Local Traj Co: %f", par_.planner_Co);
  RCLCPP_INFO(this->get_logger(), "Local Traj Cw: %f", par_.planner_Cw);
  RCLCPP_INFO(this->get_logger(), "Verbose Computation Time: %d", verbose_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Drone Bbox: (%f, %f, %f)", par_.drone_bbox[0], par_.drone_bbox[1], par_.drone_bbox[2]);
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
  RCLCPP_INFO(this->get_logger(), "Num Replanning Before Adapt: %d", par_.num_replanning_before_adapt);
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
void MIGHTY_NODE::cleanUpOldTrajsCallback()
{
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
void MIGHTY_NODE::trajCallback(const dynus_interfaces::msg::DynTraj::SharedPtr msg)
{

  // Filter out its own traj
  if (msg->id == id_)
    return;

  // Get current time
  double current_time = this->now().seconds();

  // Get dynTraj from the message
  auto traj = std::make_shared<dynTraj>();
  convertDynTrajMsg2DynTraj(*msg, traj, current_time);

  // Pass the dynTraj to mighty.cpp
  mighty_ptr_->addTraj(traj, current_time);
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function for the state of the agent
 * @param msg State message
 */
void MIGHTY_NODE::stateCallback(const dynus_interfaces::msg::State::SharedPtr msg)
{

  if (par_.use_state_update)
  {
    state current_state;
    current_state.setPos(msg->pos.x, msg->pos.y, msg->pos.z);
    current_state.setVel(msg->vel.x, msg->vel.y, msg->vel.z);
    current_state.setAccel(0.0, 0.0, 0.0);
    double roll, pitch, yaw;
    quaternion2Euler(msg->quat, roll, pitch, yaw);
    current_state.setYaw(yaw);
    mighty_ptr_->updateState(current_state);

    // publish the state
    publishCurrentState(current_state);

    // Publish the velocity in text
    if (par_.visual_level >= 1)
      publishVelocityInText(current_state.pos, current_state.vel.norm());
  }

  if (!state_initialized_)
  {

    // If we don't use state update, we need to initialize the state
    if (!par_.use_state_update)
    {
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

  if (par_.visual_level >= 1)
    publishActualTraj();
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function for replanning
 */
void MIGHTY_NODE::replanCallback()
{

  // Get the current time as double
  double current_time = this->now().seconds();

  // Set computation times to zero
  setComputationTimesToZero();

  // Replan (TODO: clean up)
  auto [replanning_result, dgp_result] = mighty_ptr_->replan(replanning_computation_time_, current_time);

  // Get computation time (used to find point A) - note this value is not updated in the replan function
  if (replanning_result)
  {
    // Get the replanning computation time
    replanning_computation_time_ = this->now().seconds() - current_time;
    if (par_.debug_verbose)
      printf("Total Replanning: %f ms\n", replanning_computation_time_ * 1000.0);
  }

  // To share trajectory with other agents
  if (replanning_result)
    publishOwnTraj();

  // For visualization of global path
  if (dgp_result && par_.visual_level >= 1)
    publishGlobalPath();

  // For visualization of free global path
  if (dgp_result && par_.visual_level >= 1)
    publishFreeGlobalPath();

  // For visualization of local_global_path and local_global_path_after_push_
  if (dgp_result && par_.visual_level >= 1)
    publishLocalGlobalPath();

  // For visualization of the local trajectory
  if (replanning_result && par_.visual_level >= 1)
    publishTraj();

  // For visualization of the safe corridor
  if (dgp_result && par_.visual_level >= 1)
    publishPoly();

  // For visualization of point G and point A
  if (replanning_result && par_.visual_level >= 1)
  {
    publishPointG();
    publishPointE();
    publishPointA();
  }

  // For visualization of control points
  if (replanning_result && par_.visual_level >= 1)
    publisCps();

  // For visualization of static push points and P points
  if (replanning_result && par_.visual_level >= 1)
  {
    mighty_ptr_->getStaticPushPoints(static_push_points_);
    publishStaticPushPoints();
  }

  // If verbose_computation_time_ or use_benchmark_ is true, we need to retrieve data from mighty_ptr_
  if (verbose_computation_time_ || use_benchmark_)
    retrieveData();

  // Verbose computation time to the terminal
  if (verbose_computation_time_)
    printComputationTime(replanning_result);

  // Record the data
  // if (replanning_result && use_benchmark_)
  if (use_benchmark_)
    recordData(replanning_result);

  // Usually this is done is goal callback but becuase we don't call that in push path test, we need to call it here
  if (par_.use_path_push_for_visualization)
    publishFOV();
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function for the terminal goal
 * @param msg Terminal goal message
 */
void MIGHTY_NODE::terminalGoalCallback(const geometry_msgs::msg::PoseStamped &msg)
{

  // Set the terminal goal
  state G_term;
  double goal_z;

  // If force_goal_z is true, set the goal_z to default_goal_z
  if (par_.force_goal_z)
    goal_z = par_.default_goal_z;
  else
    goal_z = msg.pose.position.z;

  // Check if the goal_z is within the limits
  if (goal_z < par_.z_min || goal_z > par_.z_max)
  {
    RCLCPP_ERROR(this->get_logger(), "Goal z is out of bounds: %f", goal_z);
    return;
  }

  // Set the terminal goal
  G_term.setPos(msg.pose.position.x, msg.pose.position.y, goal_z);

  // Update the terminal goal
  mighty_ptr_->setTerminalGoal(G_term);

  // Publish the term goal for visualization
  publishState(G_term, pub_point_G_term_);

  // Start replanning
  timer_replanning_->reset();

  // clear all the trajectories on we receive a new goal
  // clearMarkerActualTraj();
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish velocity in text
 * @param position Position to publish the text
 * @param velocity Velocity to publish
 */
void MIGHTY_NODE::publishVelocityInText(const Eigen::Vector3d &position, double velocity)
{

  // Set velocity's precision to 2 decimal points
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << velocity;

  // Make a string
  std::string text = oss.str() + "m/s";

  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = "map";
  marker.header.stamp = this->get_clock()->now();
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.ns = "velocity";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
  marker.scale.z = 1.0;
  marker.color.a = 1.0;
  marker.color.r = 1.0;
  marker.color.g = 1.0;
  marker.color.b = 1.0;
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
void MIGHTY_NODE::goalReachedCheckCallback()
{
  if (mighty_ptr_->goalReachedCheck())
  {
    logData();
    pub_goal_reached_->publish(std_msgs::msg::Empty());
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Callback function to get the initial pose from tf (for hardware use case)
 */
void MIGHTY_NODE::getInitialPoseHwCallback()
{
  // First find the transformation matrix from map to camera
  try
  {
    init_pose_transform_stamped_ = tf2_buffer_->lookupTransform("map", initial_pose_topic_, tf2::TimePointZero);

    // Print out the initial pose
    RCLCPP_INFO(this->get_logger(), "Initial pose received: (%f, %f, %f)", init_pose_transform_stamped_.transform.translation.x,
                init_pose_transform_stamped_.transform.translation.y, init_pose_transform_stamped_.transform.translation.z);

    // Push the initial pose to dynus
    mighty_ptr_->setInitialPose(init_pose_transform_stamped_);
  }
  catch (tf2::TransformException &ex)
  {
    RCLCPP_WARN(this->get_logger(), "Transform error: %s", ex.what());
    return;
  }

  // flag
  if (!initial_pose_received_)
  {
    initial_pose_received_ = true;
    timer_initial_pose_->cancel();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Convert dynTraj message to dynTraj
 * @param msg dynTraj message
 * @param traj dynTraj
 * @param current_time current time
 */
void MIGHTY_NODE::convertDynTrajMsg2DynTraj(const dynus_interfaces::msg::DynTraj &msg, std::shared_ptr<dynTraj> &traj, double current_time)
{

  // Inflate bbox using drone_bbox
  // We need to use the obstacle's bbox as well as ego drone's bbox
  traj->bbox << msg.bbox[0] / 2.0 + par_.drone_bbox[0] / 2.0, msg.bbox[1] / 2.0 + par_.drone_bbox[1] / 2.0, msg.bbox[2] / 2.0 + par_.drone_bbox[2] / 2.0;

  // Get id
  traj->id = msg.id;

  // Get pwp
  if (msg.mode == "pwp")
  {
    traj->pwp = mighty_utils::convertPwpMsg2Pwp(msg.pwp);
    traj->mode = dynTraj::Mode::Piecewise;
  }

  // Find quihtic coefficients from the given pwp
  if (msg.mode == "quintic")
  {
    traj->cx = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_x);
    traj->cy = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_y);
    traj->cz = mighty_utils::convertCoeffMsg2Coeff(msg.poly_coeffs_z);
    traj->poly_start_time = msg.poly_start_time;
    traj->poly_end_time = msg.poly_end_time;
    traj->mode = dynTraj::Mode::Quintic;
  }

  if (msg.mode == "analytic")
  {
    traj->mode = dynTraj::Mode::Analytic;
  }

  // Get covariances
  if (!msg.is_agent)
  {

    if (msg.ekf_cov_p.size() != 0)
    {
      traj->ekf_cov_p = mighty_utils::convertCovMsg2Cov(msg.ekf_cov_p); // ekf cov
    }

    if (msg.ekf_cov_q.size() != 0)
    {
      traj->ekf_cov_q = mighty_utils::convertCovMsg2Cov(msg.ekf_cov_q); // ekf cov
    }

    if (msg.poly_cov.size() != 0)
    {
      traj->poly_cov = mighty_utils::convertCovMsg2Cov(msg.poly_cov); // future traj cov
    }

    if (msg.function.size() == 3)
    {
      traj->traj_x = msg.function[0];
      traj->traj_y = msg.function[1];
      traj->traj_z = msg.function[2];
    }

    if (msg.velocity.size() == 3)
    {
      traj->traj_vx = msg.velocity[0];
      traj->traj_vy = msg.velocity[1];
      traj->traj_vz = msg.velocity[2];
    }

    if (msg.function.size() == 3 && msg.velocity.size() == 3)
    {
      if (traj->compileAnalytic())
      {
        // Change the mode only when we successfully compiled the analytic trajectory
        traj->mode = dynTraj::Mode::Analytic;
        // printf("Successfully compiled analytic traj id=%d\n", traj->id);
      }
      else
      {
        RCLCPP_ERROR(
            this->get_logger(),
            "Failed to compile analytic traj id=%d, falling back to zeros.",
            traj->id);
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
  if (traj->is_agent)
    traj->goal << msg.goal[0], msg.goal[1], msg.goal[2];

  // Get communication delay
  if (traj->is_agent && par_.use_comm_delay_inflation)
  {
    // Get the delay (current time - msg time)
    traj->communication_delay = this->now().seconds();
    -msg.header.stamp.sec + msg.header.stamp.nanosec * 1e-9;

    // Sanity check - if the delay is negative, set it to 0 - send warning message: it's probably due to the clock synchronization issue
    if (traj->communication_delay < 0)
    {
      traj->communication_delay = 0;
      RCLCPP_WARN(this->get_logger(), "Communication delay is negative. Setting it to 0.");
    }
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish control points
 */

void MIGHTY_NODE::publisCps()
{

  // Retrieve control points
  mighty_ptr_->retrieveCPs(cps_);

  // Create a marker array
  visualization_msgs::msg::MarkerArray marker_array;
  marker_array.markers.resize(cps_.size());

  // Loop through the control points (std::vector<Eigen::Matrix<double, 3, 4>>)
  for (int seg = 0; seg < cps_.size(); seg++)
  {

    // Create a marker
    visualization_msgs::msg::Marker marker;
    marker.header.frame_id = "map";
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
    if (seg % 3 == 0)
    {
      marker.color.r = 1.0;
      marker.color.g = 0.0;
      marker.color.b = 0.0;
    }
    else if (seg % 3 == 1)
    {
      marker.color.r = 0.0;
      marker.color.g = 1.0;
      marker.color.b = 0.0;
    }
    else
    {
      marker.color.r = 0.0;
      marker.color.g = 0.0;
      marker.color.b = 1.0;
    }

    auto cp = cps_[seg];

    // Loop through the control points for each segment
    for (int i = 0; i < cp.cols(); i++)
    {
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
void MIGHTY_NODE::publishStaticPushPoints()
{

  // Create a marker array
  visualization_msgs::msg::MarkerArray marker_array;
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = "map";
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
  for (int idx = 0; idx < static_push_points_.size(); idx++)
  {
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
void MIGHTY_NODE::setComputationTimesToZero()
{
  final_g_ = 0.0;
  global_planning_time_ = 0.0;
  dgp_static_jps_time_ = 0.0;
  dgp_check_path_time_ = 0.0;
  dgp_dynamic_astar_time_ = 0.0;
  dgp_recover_path_time_ = 0.0;
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
void MIGHTY_NODE::retrieveData()
{
  mighty_ptr_->retrieveData(final_g_,
                            global_planning_time_,
                            dgp_static_jps_time_,
                            dgp_check_path_time_,
                            dgp_dynamic_astar_time_,
                            dgp_recover_path_time_,
                            cvx_decomp_time_,
                            initial_guess_computation_time_,
                            local_traj_computation_time_,
                            safety_check_time_,
                            safe_paths_time_,
                            yaw_sequence_time_,
                            yaw_fitting_time_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Print the computation times
 */
void MIGHTY_NODE::printComputationTime(bool result)
{
  // Print the computation times
  RCLCPP_INFO(this->get_logger(), "Planner: %s", par_.global_planner.c_str());
  RCLCPP_INFO(this->get_logger(), "Result: %d", result);
  RCLCPP_INFO(this->get_logger(), "Cost (final node's g): %f", final_g_);
  RCLCPP_INFO(this->get_logger(), "Total replanning time [ms]: %f", replanning_computation_time_ * 1000.0);
  RCLCPP_INFO(this->get_logger(), "Global Planning Time [ms]: %f", global_planning_time_);
  RCLCPP_INFO(this->get_logger(), "CVX Decomposition Time [ms]: %f", cvx_decomp_time_);
  RCLCPP_INFO(this->get_logger(), "Initial Guess Time [ms]: %f", initial_guess_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Local Traj Time [ms]: %f", local_traj_computation_time_);
  RCLCPP_INFO(this->get_logger(), "Safe Paths Time [ms]: %f", safe_paths_time_);
  RCLCPP_INFO(this->get_logger(), "Safety Check Time [ms]: %f", safety_check_time_);
  RCLCPP_INFO(this->get_logger(), "Yaw Sequence Time [ms]: %f", yaw_sequence_time_);
  RCLCPP_INFO(this->get_logger(), "Yaw Fitting Time [ms]: %f", yaw_fitting_time_);
  if (par_.global_planner == "dgp")
  {
    RCLCPP_INFO(this->get_logger(), "Static JPS Time [ms]: %f", dgp_static_jps_time_);
    RCLCPP_INFO(this->get_logger(), "Check Path Time [ms]: %f", dgp_check_path_time_);
    RCLCPP_INFO(this->get_logger(), "Dynamic A* Time [ms]: %f", dgp_dynamic_astar_time_);
    RCLCPP_INFO(this->get_logger(), "Recover Path Time [ms]: %f", dgp_recover_path_time_);
  }
  RCLCPP_INFO(this->get_logger(), "------------------------");
}

// ----------------------------------------------------------------------------

/**
 * @brief Record the data
 * @param result result of the replanning
 */
void MIGHTY_NODE::recordData(bool result)
{

  // Record all the data into global_path_benchmark_
  std::tuple<bool, double, double, double, double, double, double, double, double, double, double, double, double, double, double> data;
  if (par_.global_planner == "dgp")
  {
    data = std::make_tuple(result, final_g_, replanning_computation_time_, global_planning_time_, cvx_decomp_time_, initial_guess_computation_time_, local_traj_computation_time_, safe_paths_time_, safety_check_time_, yaw_sequence_time_, yaw_fitting_time_, dgp_static_jps_time_, dgp_check_path_time_, dgp_dynamic_astar_time_, dgp_recover_path_time_);
  }
  else
  {
    data = std::make_tuple(result, final_g_, replanning_computation_time_, global_planning_time_, cvx_decomp_time_, initial_guess_computation_time_, local_traj_computation_time_, safe_paths_time_, safety_check_time_, yaw_sequence_time_, yaw_fitting_time_, 0.0, 0.0, 0.0, 0.0);
  }

  global_path_benchmark_.push_back(data);
}

// ----------------------------------------------------------------------------

/**
 * @brief Log the data to a csv file
 */
void MIGHTY_NODE::logData()
{

  // Loc the computation times to csv file
  // std::ofstream log_file(file_path_, std::ios_base::app); // Open the file in append mode
  std::ofstream log_file(file_path_); // Open the file in overwrite mode
  if (log_file.is_open())
  {
    if (par_.global_planner == "dgp")
    {
      // Header
      log_file << "Planner,Result,Cost (final node's g),Total replanning time [ms],Global Planning Time [ms],CVX Decomposition Time [ms],Initial Guess Time [ms],Local Traj Time [ms],Safe Paths Time [ms],Safety Check Time [ms],Yaw Sequence Time [ms],Yaw Fitting Time [ms],Static JPS Time [ms],Check Path Time [ms],Dynamic A* Time [ms],Recover Path Time [ms]\n";
    }
    else
    {
      // Header
      log_file << "Planner,Result,Cost (final node's g),Total replanning time [ms],Global Planning Time [ms],CVX Decomposition Time [ms],Initial Guess Time [ms],Local Traj Time [ms],Safe Paths Time [ms],Safety Check Time [ms],Yaw Sequence Time [ms],Yaw Fitting Time [ms]\n";
    }

    // Data
    for (const auto &row : global_path_benchmark_)
    {
      if (par_.global_planner == "dgp")
      {
        log_file << par_.global_planner << "," << std::get<0>(row) << "," << std::get<1>(row) << "," << std::get<2>(row) * 1000.0 << "," << std::get<3>(row) << "," << std::get<4>(row) << "," << std::get<5>(row) << "," << std::get<6>(row) << "," << std::get<7>(row) << "," << std::get<8>(row) << "," << std::get<9>(row) << "," << std::get<10>(row) << "," << std::get<11>(row) << "," << std::get<12>(row) << "," << std::get<13>(row) << std::get<14>(row) << "\n";
      }
      else
      {
        log_file << par_.global_planner << "," << std::get<0>(row) << "," << std::get<1>(row) << "," << std::get<2>(row) * 1000.0 << "," << std::get<3>(row) << "," << std::get<4>(row) << "," << std::get<5>(row) << "," << std::get<6>(row) << "," << std::get<7>(row) << "," << std::get<8>(row) << "," << std::get<9>(row) << "," << std::get<10>(row) << "\n";
      }
    }

    log_file.close();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the Point G (sub goal)
 */
void MIGHTY_NODE::publishPointG() const
{

  // get projected goal (G)
  state G;
  mighty_ptr_->getG(G);

  // Publish the goal for visualization
  publishState(G, pub_point_G_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the Point E (sub goal)
 */
void MIGHTY_NODE::publishPointE() const
{

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
void MIGHTY_NODE::publishPointA() const
{

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
void MIGHTY_NODE::publishCurrentState(const state &state) const
{
  // Publish the goal for visualization
  publishState(state, pub_current_state_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish state
 */
void MIGHTY_NODE::publishState(const state &data, const rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr &publisher) const
{
  geometry_msgs::msg::PointStamped p;
  p.header.frame_id = "map";
  p.header.stamp = this->now();
  p.point = eigen2point(data.pos);
  publisher->publish(p);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish its own trajectory for deconfliction
 */
void MIGHTY_NODE::publishOwnTraj()
{

  dynus_interfaces::msg::DynTraj msg;
  msg.header.stamp = this->now();
  msg.header.frame_id = "map";
  msg.bbox.push_back(par_.drone_bbox[0]);
  msg.bbox.push_back(par_.drone_bbox[1]);
  msg.bbox.push_back(par_.drone_bbox[2]);
  msg.id = id_;
  // msg.pwp = mighty_utils::convertPwp2PwpMsg(pwp_to_share_);
  msg.is_agent = true;

  // Get the terminal goal
  state G;
  mighty_ptr_->getG(G);
  msg.goal.push_back(G.pos(0));
  msg.goal.push_back(G.pos(1));
  msg.goal.push_back(G.pos(2));

  // Publish the trajectory
  pub_own_traj_->publish(msg);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the trajectory the agent actually followed for visualization
 */
void MIGHTY_NODE::publishActualTraj()
{
  // Initialize the previous point
  static geometry_msgs::msg::Point prev_p = pointOrigin();

  // Get the current state and position
  state current_state;
  mighty_ptr_->getState(current_state);
  Eigen::Vector3d current_pos = current_state.pos;

  if (current_pos.norm() < 1e-2)
    return; // because the state is not updated yet

  // If we use UAV, we can just use the state topic published by fake_sim, but if we use ground robot, since we use TF for state publisher, we cannot get velocity info from the state topic. So we will approximiate
  if (par_.vehicle_type != "uav")
  {

    // Initialize the previous position and time
    if (!publish_actual_traj_called_)
    {
      actual_traj_prev_pos_ = current_pos;
      actual_traj_prev_time_ = this->now().seconds();
      publish_actual_traj_called_ = true;
      return;
    }

    // Get the velocity
    current_state.vel = (current_pos - actual_traj_prev_pos_) / (this->now().seconds() - actual_traj_prev_time_);
  }

  // Set up the marker
  visualization_msgs::msg::Marker m;
  m.type = visualization_msgs::msg::Marker::ARROW;
  m.action = visualization_msgs::msg::Marker::ADD;
  m.id = actual_traj_id_;
  m.ns = "actual_traj_" + id_str_;
  m.color = getColorJet(current_state.vel.norm(), 0, par_.v_max); // note that par_.v_max is per axis
  m.scale.x = 0.15;
  m.scale.y = 0.0001;
  m.scale.z = 0.0001;
  m.header.stamp = this->now();
  m.header.frame_id = "map";

  // pose is actually not used in the marker, but if not RVIZ complains about the quaternion
  m.pose.position = pointOrigin();
  m.pose.orientation.x = 0.0;
  m.pose.orientation.y = 0.0;
  m.pose.orientation.z = 0.0;
  m.pose.orientation.w = 1.0;

  // Set the points
  geometry_msgs::msg::Point p;
  p = mighty_utils::convertEigen2Point(current_pos);
  m.points.push_back(prev_p);
  m.points.push_back(p);
  prev_p = p;

  // Return if the actual_traj_id_ is 0 - avoid publishing the first point which goes from the origin to the first point
  if (actual_traj_id_ == 0)
  {
    actual_traj_id_++;
    return;
  }
  actual_traj_id_++;

  // Publish the marker
  pub_actual_traj_->publish(m);
}

// ----------------------------------------------------------------------------

/**
 * @brief Clear the marker array
 */
void MIGHTY_NODE::clearMarkerActualTraj()
{
  visualization_msgs::msg::Marker m;
  m.type = visualization_msgs::msg::Marker::ARROW;
  m.action = visualization_msgs::msg::Marker::DELETEALL;
  m.id = 0;
  m.scale.x = 0.02;
  m.scale.y = 0.04;
  m.scale.z = 1;
  pub_actual_traj_->publish(m);
  actual_traj_id_ = 0;
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish goal (setpoint)
 */
void MIGHTY_NODE::publishGoal()
{

  // Initialize the goal
  state next_goal;

  // Get the next goal
  if (mighty_ptr_->getNextGoal(next_goal) && par_.use_state_update)
  {

    // Publish the goal (actual setpoint)
    dynus_interfaces::msg::Goal quadGoal;
    quadGoal.header.stamp = this->now();
    quadGoal.header.frame_id = "map";
    quadGoal.p = eigen2rosvector(next_goal.pos);
    quadGoal.v = eigen2rosvector(next_goal.vel);
    quadGoal.a = eigen2rosvector(next_goal.accel);
    quadGoal.j = eigen2rosvector(next_goal.jerk);
    quadGoal.yaw = next_goal.yaw;
    quadGoal.dyaw = next_goal.dyaw;
    pub_goal_->publish(quadGoal);

    // Publish the goal (setpoint) for visualization
    if (par_.visual_level >= 1)
      publishState(next_goal, pub_setpoint_);
  }

  // Publish FOV
  if (par_.visual_level >= 1)
    publishFOV();
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish Sefe Corridor Polyhedra
 */
void MIGHTY_NODE::publishPoly()
{

  // retrieve the polyhedra
  mighty_ptr_->retrievePolytopes(poly_whole_, poly_safe_);

  // For whole trajectory
  if (!poly_whole_.empty())
  {
    decomp_ros_msgs::msg::PolyhedronArray poly_whole_msg = DecompROS::polyhedron_array_to_ros(poly_whole_);
    poly_whole_msg.header.stamp = this->now();
    poly_whole_msg.header.frame_id = "map";
    poly_whole_msg.lifetime = rclcpp::Duration::from_seconds(1.0);
    pub_poly_whole_->publish(poly_whole_msg);
  }

  // For safe trajectory
  if (!poly_safe_.empty())
  {
    decomp_ros_msgs::msg::PolyhedronArray poly_safe_msg = DecompROS::polyhedron_array_to_ros(poly_safe_);
    poly_safe_msg.header.stamp = this->now();
    poly_safe_msg.header.frame_id = "map";
    poly_safe_msg.lifetime = rclcpp::Duration::from_seconds(1.0);
    pub_poly_safe_->publish(poly_safe_msg);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the trajectory
 */
void MIGHTY_NODE::publishTraj()
{
  auto now = this->now();

  // 1) DELETEALL on both topics
  {
    visualization_msgs::msg::MarkerArray clear_msg;
    visualization_msgs::msg::Marker clear_m;
    clear_m.header.frame_id = "map";
    clear_m.header.stamp = now;
    clear_m.action = visualization_msgs::msg::Marker::DELETEALL;
    clear_msg.markers.push_back(clear_m);

    pub_traj_committed_colored_->publish(clear_msg);
    pub_traj_subopt_colored_->publish(clear_msg);
  }

  // 2) Publish the committed (best) trajectory
  mighty_ptr_->retrieveGoalSetpoints(goal_setpoints_);
  {
    auto committed_ma = stateVector2ColoredMarkerArray(
        goal_setpoints_,
        /*type=*/1,
        par_.v_max,
        now);
    pub_traj_committed_colored_->publish(committed_ma);
  }

  // 3) Publish all sub-optimal trajectories
  if (par_.use_multiple_initial_guesses)
  {
    mighty_ptr_->retrieveListSubOptGoalSetpoints(list_subopt_goal_setpoints_);

    visualization_msgs::msg::MarkerArray subopt_ma;
    for (int i = 0; i < (int)list_subopt_goal_setpoints_.size(); ++i)
    {
      auto single = stateVector2ColoredMarkerArray(
          list_subopt_goal_setpoints_[i],
          /*type=*/i + 2,
          par_.v_max,
          now);
      // append all markers from this one:
      subopt_ma.markers.insert(
          subopt_ma.markers.end(),
          single.markers.begin(),
          single.markers.end());
    }
    pub_traj_subopt_colored_->publish(subopt_ma);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the global path (that can go through unknown space)
 */
void MIGHTY_NODE::publishGlobalPath()
{

  int global_path_color = RED;
  int original_global_path_color = ORANGE;

  // Generate random integer from 1 to 10 to generate random color
  if (par_.use_random_color_for_global_path)
    global_path_color = rand() % 10 + 1;

  // Get global_path
  vec_Vecf<3> global_path;
  mighty_ptr_->getGlobalPath(global_path);

  if (!global_path.empty())
  {
    // Publish global_path
    clearMarkerArray(dgp_path_marker_, pub_dgp_path_marker_);
    vectorOfVectors2MarkerArray(global_path, &dgp_path_marker_, color(global_path_color));
    pub_dgp_path_marker_->publish(dgp_path_marker_);
  }

  // Get the original global path
  vec_Vecf<3> original_global_path;
  mighty_ptr_->getOriginalGlobalPath(original_global_path);

  if (!original_global_path.empty())
  {
    // Publish original_global_path
    clearMarkerArray(original_dgp_path_marker_, pub_original_dgp_path_marker_);
    vectorOfVectors2MarkerArray(original_global_path, &original_dgp_path_marker_, color(original_global_path_color));
    pub_original_dgp_path_marker_->publish(original_dgp_path_marker_);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the free global path (that only goes through free space)
 */
void MIGHTY_NODE::publishFreeGlobalPath()
{

  // Get free_global_path
  vec_Vecf<3> free_global_path;
  mighty_ptr_->getFreeGlobalPath(free_global_path);

  if (free_global_path.empty())
    return;

  // Publish free_global_path
  clearMarkerArray(dgp_free_path_marker_, pub_free_dgp_path_marker_);
  vectorOfVectors2MarkerArray(free_global_path, &dgp_free_path_marker_, color(GREEN));
  pub_free_dgp_path_marker_->publish(dgp_free_path_marker_);
}

// ----------------------------------------------------------------------------

/**
 * @brief Publish the local_global_path and local_global_path_after_push_
 */
void MIGHTY_NODE::publishLocalGlobalPath()
{

  // Get the local global path and local global path after push
  vec_Vecf<3> local_global_path;
  vec_Vecf<3> local_global_path_after_push;
  mighty_ptr_->getLocalGlobalPath(local_global_path, local_global_path_after_push);

  if (!local_global_path.empty())
  {
    // Publish local_global_path
    clearMarkerArray(dgp_local_global_path_marker_, pub_local_global_path_marker_);
    vectorOfVectors2MarkerArray(local_global_path, &dgp_local_global_path_marker_, color(BLUE));
    pub_local_global_path_marker_->publish(dgp_local_global_path_marker_);
  }

  if (!local_global_path_after_push.empty())
  {
    // Publish local_global_path_after_push
    clearMarkerArray(dgp_local_global_path_after_push_marker_, pub_local_global_path_after_push_marker_);
    vectorOfVectors2MarkerArray(local_global_path_after_push, &dgp_local_global_path_after_push_marker_, color(ORANGE));
    pub_local_global_path_after_push_marker_->publish(dgp_local_global_path_after_push_marker_);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Create MarkerArray from vec_Vec3f
 */
void MIGHTY_NODE::createMarkerArrayFromVec_Vec3f(
    const vec_Vec3f &occupied_cells, const std_msgs::msg::ColorRGBA &color, int namespace_id, double scale, visualization_msgs::msg::MarkerArray *marker_array)
{

  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = "map";              // Set the appropriate frame
  marker.header.stamp = rclcpp::Clock().now(); // Use ROS2 clock
  marker.ns = "namespace_" + std::to_string(namespace_id);
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::CUBE_LIST; // Each point will be visualized as a cube
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.scale.x = par_.res;
  marker.scale.y = par_.res;
  marker.scale.z = par_.res;
  marker.color = color;

  for (const auto &cell : occupied_cells)
  {
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
void MIGHTY_NODE::clearMarkerArray(visualization_msgs::msg::MarkerArray &path_marker, rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr publisher)
{

  // If the marker array is empty, return
  if (path_marker.markers.size() == 0)
    return;

  // Clear the marker array
  int id_begin = path_marker.markers[0].id;

  for (int i = 0; i < path_marker.markers.size(); i++)
  {
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
void MIGHTY_NODE::constructFOVMarker()
{

  marker_fov_.header.stamp = this->now();
  marker_fov_.header.frame_id = d435_depth_frame_id_;
  marker_fov_.ns = "marker_fov";
  marker_fov_.id = marker_fov_id_++;
  marker_fov_.frame_locked = true;
  marker_fov_.type = marker_fov_.LINE_LIST;
  marker_fov_.action = marker_fov_.ADD;
  marker_fov_.pose = mighty_utils::identityGeometryMsgsPose();

  double delta_y = par_.fov_visual_depth * fabs(tan((par_.fov_visual_x_deg * M_PI / 180) / 2.0));
  double delta_z = par_.fov_visual_depth * fabs(tan((par_.fov_visual_y_deg * M_PI / 180) / 2.0));

  geometry_msgs::msg::Point v0 = eigen2point(Eigen::Vector3d(0.0, 0.0, 0.0));
  geometry_msgs::msg::Point v1 = eigen2point(Eigen::Vector3d(-delta_y, delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v2 = eigen2point(Eigen::Vector3d(delta_y, delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v3 = eigen2point(Eigen::Vector3d(delta_y, -delta_z, par_.fov_visual_depth));
  geometry_msgs::msg::Point v4 = eigen2point(Eigen::Vector3d(-delta_y, -delta_z, par_.fov_visual_depth));

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
void MIGHTY_NODE::publishFOV()
{
  marker_fov_.header.stamp = this->now();
  pub_fov_->publish(marker_fov_);
  return;
}

// ----------------------------------------------------------------------------

void MIGHTY_NODE::mapCallback(
    const sensor_msgs::msg::PointCloud2::ConstPtr &map_msg,
    const sensor_msgs::msg::PointCloud2::ConstPtr &unk_msg)
{
  // use PCL’s own Ptr (boost::shared_ptr)
  pcl::PointCloud<pcl::PointXYZ>::Ptr map_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*map_msg, *map_pc);

  pcl::PointCloud<pcl::PointXYZ>::Ptr unk_pc(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromROSMsg(*unk_msg, *unk_pc);

  mighty_ptr_->updateMap(map_pc, unk_pc);
}

// ----------------------------------------------------------------------------

int main(int argc, char **argv)
{

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