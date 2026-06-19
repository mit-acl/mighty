/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <algorithm>  // for std::clamp
#include <iostream>
#include <limits>     // for std::numeric_limits in dynTraj::getHorizon
#include <memory>
#include <vector>

#include <Eigen/Core>

#include <mighty/lbfgs_solver_utils.hpp>
#include <sim/exprtk.hpp>

#include "hgp/data_type.hpp"

struct StateDeriv {
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Vector3d accel;
  Eigen::Vector3d jerk;
};

struct polytope {
  Eigen::MatrixXd A;
  Eigen::MatrixXd b;
};

struct parameters {
  // Sim enviroment
  std::string sim_env;

  // UAV or Ground robot
  std::string vehicle_type;
  bool provide_goal_in_global_frame;
  bool use_hardware;
  std::string map_frame_id{"map"};
  bool share_traj{true};
  bool use_frame_alignment{false};
  int num_agents{10};

  // Simulated frame offset (for testing frame alignment in fake_sim)
  // Applied to published trajectory to simulate operating in a rotated frame
  double sim_frame_offset_qx{0.0};
  double sim_frame_offset_qy{0.0};
  double sim_frame_offset_qz{0.0};
  double sim_frame_offset_qw{1.0};

  // Formation flight (per-agent membership + desired pairwise offsets)
  bool use_formation{false};
  double formation_weight{0.0};
  std::vector<double> formation_self_offset{0.0, 0.0, 0.0};   // δ_i (length 3)
  std::vector<int64_t> formation_neighbor_ids;                // neighbor agent IDs
  std::vector<double> formation_neighbor_offsets;             // flat 3*N: δ_ij

  // Flight mode
  std::string flight_mode;

  // Visual level
  int visual_level;

  // Global planner parameters
  std::string global_planner;
  bool global_planner_verbose;
  double global_planner_heuristic_weight;
  double factor_hgp;
  double inflation_hgp;
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;
  double drone_radius;
  int hgp_timeout_duration_ms;
  int max_expand{10000};

  int max_num_expansion{10000};

  bool use_free_start;
  double free_start_factor;
  bool use_free_goal;
  double free_goal_factor;
  bool relocate_occupied_goal;

  // LOS post processing parameters
  int los_cells;
  double min_len;   // [m] minimum length between two waypoints after post processing
  double min_turn;  // [deg] minimum turn angle after post processing
  double heat_cutoff_ratio{
      0.5};  // [-] Cells with heat > ratio * Hmax are impassable (0=disabled, 1=all blocked)
  bool disable_all_smoothing{false};  // [-] Disable all path smoothing (use raw A* output)
  bool skip_path_smoothing{false};    // [-] Skip LoS shortcutting, use Laplacian smoothing instead
  int smooth_iterations{50};          // [-] Number of heat-aware Laplacian smoothing passes
  double smooth_alpha{0.3};           // [-] Smoothing relaxation factor (0=none, 1=aggressive)

  // Path push visualization parameters
  bool use_state_update;
  bool use_random_color_for_global_path;
  bool use_path_push_for_visualization;

  // Decomposition parameters
  std::vector<double> local_box_size;
  double min_dist_from_agent_to_traj;
  bool use_shrinked_box;
  double shrinked_box_size;

  // Map parameters
  double map_buffer;
  double center_shift_factor;
  double initial_wdx;
  double initial_wdy;
  double initial_wdz;
  double min_wdx;
  double min_wdy;
  double min_wdz;
  double res;

  // Heat map parameters (soft cost layer)
  bool use_heat_map{false};
  float heat_weight{1.0f};
  bool dynamic_heat_enabled{false};
  bool dynamic_as_occupied_current{true};
  bool dynamic_as_occupied_future{false};
  float heat_alpha0{1.0f};
  float heat_alpha1{2.0f};
  int heat_p{2};
  int heat_q{2};
  float heat_tau_ratio{0.5f};
  float heat_gamma{0.0f};
  float heat_Hmax{10.0f};
  float obst_max_vel{1.0f};
  float dyn_base_inflation_m{0.5f};
  float dyn_heat_tube_radius_m{0.5f};
  int heat_num_samples{15};
  double prediction_horizon{3.0};
  double prediction_mask_distance{2.0};  ///< [m] Drop predicted trajs if an agent traj is within this distance
  bool static_heat_enabled{false};
  float static_heat_alpha{2.0f};
  int static_heat_p{2};
  float static_heat_Hmax{50.0f};
  float static_heat_rmax_m{1.0f};
  float static_heat_default_radius_m{0.5f};
  bool static_heat_boundary_only{true};
  bool static_heat_apply_on_unknown{false};
  bool static_heat_exclude_dynamic{true};

  // Soft-cost obstacle parameters
  bool use_soft_cost_obstacles{false};
  float obstacle_soft_cost{100.0f};

  // HGP-specific parameters (from sando)
  double obst_position_error{0.0};
  bool inflate_unknown_boundary{false};
  bool sfc_use_unknown_as_obstacle{false};

  // Communication delay parameters
  bool use_comm_delay_inflation;
  double comm_delay_inflation_alpha;
  double comm_delay_inflation_max;
  double comm_delay_filter_alpha;

  // Simulation parameters
  double depth_camera_depth_max;
  double fov_visual_depth;
  double fov_visual_x_deg;
  double fov_visual_y_deg;
  bool use_sphere_sensing;
  double sphere_sensing_radius;

  // number of segments parameters
  double max_dist_vertexes;  // [m] Maximum distance between two consecutive vertexes
  double w_unknown;          // [-] Weight for the unknown cells in the global planner
  double w_align;            // strength of alignment penalty (cells)
  double decay_len_cells;    // e-folding distance from the start (cells)
  double w_side;             // side (handedness) tie-break strength (cells)

  // Initial guess parameters
  bool use_multiple_initial_guesses;  // [-] Use multiple initial guesses
  int num_perturbation_for_ig = 8;    // [-] Number of perturbations for the initial guess
  double r_max_for_ig = 1.0;          // [m] radius for the initial guess perturbation

  // Optimiztion parameters
  int num_N;
  double horizon;
  double dc;
  double v_max;
  double a_max;
  double j_max;
  bool closed_form_traj_verbose;
  double jerk_weight;
  double dynamic_weight;
  double time_weight;
  double pos_anchor_weight;
  double stat_weight;
  double dyn_constr_bodyrate_weight;
  double dyn_constr_tilt_weight;
  double dyn_constr_thrust_weight;
  double dyn_constr_vel_weight;
  double dyn_constr_acc_weight;
  double dyn_constr_jerk_weight;
  int num_dyn_obst_samples;  // Number of dynamic obstacle samples
  double planner_Co;         // for static obstacle avoidance
  double planner_Cw;
  std::vector<double> drone_bbox;
  double goal_radius;
  double goal_seen_radius;
  double init_turn_bf;      // initial turn buffer in degrees
  int integral_resolution;  // resolution for the integral in the optimization
  double hinge_mu;          // hinge mu for the optimization
  double omega_max;         // max body rate in rad/s
  double tilt_max_rad;      // max tilt in radians
  double f_min;             // min thrust in N
  double f_max;             // max thrust in N
  double mass;              // mass in kg
  double g;                 // gravity in m/s^2
  double fopt_threshold;    // threshold for the fopt to consider the optimization successful
  int spline_degree{5};     // [-] Hermite spline degree: 3 (cubic) or 5 (quintic)

  // L-BFGS parameters
  double f_dec_coeff;      // allow larger Armijo steps
  double cautious_factor;  // always accept BFGS update
  int past;                // number of past iterations to use
  int max_linesearch;      // fewer backtracking tries
  int max_iterations;      // allow more iterations
  double g_epsilon;        // gradient norm threshold for convergence
  double delta;            // stop once f-improvement is minimal

  // Dynamic obstacles parameters
  double traj_lifetime;

  // Dynamic k_value parameters
  int num_replanning_before_adapt;
  int default_k_value;
  double alpha_k_value_filtering;
  double k_value_factor;

  // Yaw-related parameters
  double alpha_filter_dyaw;
  double w_max;
  int yaw_spinning_threshold;
  double yaw_spinning_dyaw;

  // Simulation env parameters
  bool force_goal_z;
  double default_goal_z;

  // Debug flags
  bool debug_verbose;

  // Ground robot control parameters
  double ground_robot_kx{0.1};     // Forward/backward position gain
  double ground_robot_ky{0.1};     // Lateral position gain
  double ground_robot_kyaw{1.0};   // Yaw heading gain
  double ground_robot_eps{0.1};    // Distance threshold for goal achievement
  double ground_robot_v_max{1.0};  // Max linear velocity (m/s)
  double ground_robot_w_max{1.5};  // Max angular velocity (rad/s)
  double ground_robot_L_min{1.0};  // [m] Minimum lookahead distance for pure pursuit and point A

  // 2D ground robot planning parameters
  bool use_2d_planning{false};         // Master toggle for 2D ground planning
  double robot_height{0.5};            // [m] Robot height for obstacle column detection
  double obstacle_min_height{0.3};     // [m] Min height span in column to classify as obstacle
  bool use_column_any_occupied{true};  // [-] Any occupied voxel in column → 2D occupied
  double column_min_z{0.15};  // [m] Ignore occupied voxels below this z (filters ground floor)
  double terrain_cost_weight{1.0};       // [-] Multiplier for terrain gradient cost in A*
  std::string terrain_cost_mode{"max"};  // "max" or "avg" of neighbor height deltas
  double ground_slab_margin{0.3};  // [m] z margin for floor/ceiling virtual obstacles in corridors

  // ESDF-based obstacle avoidance (ground robot only)
  bool use_esdf_cost{false};            // [-] Enable ESDF distance cost in L-BFGS
  double esdf_weight{1e+3};            // [-] ESDF cost weight in objective
  double esdf_d_safe{1.0};             // [m] Safety distance threshold
  int esdf_truncation_distance{10};    // [voxels] Must match mapper config

  // Bend pre-alignment (ground robot only).
  // At each sharp HGP bend, back off the subgoal along the incoming direction
  // by corridor_backoff_m so the robot stops short of the inside corner,
  // turns in place to face the outgoing segment, then drives straight through
  // the bend. Works in narrow corridors because the robot body never
  // approaches the inside wall during the turn.
  //
  // Detection runs on the RAW A* path (not the resampled / smoothed path) so
  // corner angles are not blurred by the resample phase. A windowed detector
  // averages direction over corridor_detection_window_m on each side of every
  // candidate cell, which kills A* diagonal-grid noise while preserving real
  // 90° bends.
  bool corridor_hop_enabled{false};
  double corridor_corner_angle_deg{75.0};    // min direction change [deg] to count as a real bend
  double corridor_detection_window_m{0.8};   // arclength used to average incoming/outgoing direction
  double corridor_backoff_m{0.4};            // subgoal offset backward along the incoming segment
  double corridor_min_leg_m{0.3};            // min distance from A to the (backed-off) subgoal
  // Legacy ESDF gradient-ascent snap (kept for build compat — typically off).
  double corridor_clearance_threshold_m{0.7};
  double corridor_max_ascent_m{0.0};         // 0 disables the ascent
  double corridor_ascent_step_m{0.05};
  int    corridor_ascent_max_iters{0};       // 0 disables the ascent

  // Trajectory publishing parameters
  int trajectory_downsample_points{500};  // Number of points to downsample trajectory to
  double mpc_path_spacing{0.05};          // [m] Spacing between waypoints in MPC path

  // Frontier-based exploration (ground robot only).
  // Master toggle is `expl_enabled`. When enabled, mighty_node runs a frontier
  // detector + persistent global frontier database on the published 2D occupancy
  // grid (occ_2d_topic), and autonomously issues exploration goals via the same
  // pathway as a manual term_goal. A manual term_goal preempts exploration; the
  // robot resumes exploration after the manual goal is reached.
  bool   expl_enabled{false};
  double expl_select_rate_hz{1.0};
  double expl_default_goal_z{0.0};
  // Detector
  int    expl_cluster_min_cells{6};
  int    expl_border_margin_cells{2};
  int    expl_obstacle_clearance_cells{1};
  double expl_robot_snap_radius_m{1.0};
  // ESDF-based obstacle clearance for frontier centroids. Frontiers whose
  // centroid is within this distance of any obstacle (per the 2D ESDF) are
  // dropped and existing records are invalidated. Set to 0 (or no ESDF
  // available) to disable. Meters.
  double expl_min_obstacle_distance_m{0.0};
  // Optional axis-aligned exploration bounds (world frame). Frontier seeds
  // outside the box are dropped, so the robot only receives goals inside it.
  bool   expl_bounds_enabled{false};
  double expl_bounds_min_x{-50.0};
  double expl_bounds_max_x{ 50.0};
  double expl_bounds_min_y{-50.0};
  double expl_bounds_max_y{ 50.0};
  // Utility weights
  double expl_w_size{1.0};
  double expl_w_dist{2.0};
  double expl_w_info{1.0};
  double expl_w_revisit{0.5};
  double expl_w_heading{0.3};
  double expl_size_ref_m2{5.0};
  double expl_dist_ref_m{25.0};
  double expl_sensor_radius_m{5.0};
  double expl_goal_select_threshold{-1.0e9};
  // Manager / lifecycle
  double expl_merge_radius_m{1.0};
  double expl_centroid_ema_alpha{0.5};
  double expl_visit_radius_m{2.0};
  double expl_visit_dwell_sec{1.0};
  int    expl_verify_radius_cells{2};
  int    expl_max_frontiers{1000};
  int    expl_unreachable_consec_thresh{5};
  // Pursuit timeout — auto-invalidate a frontier we've been chasing too long.
  // Budget = max(min_sec, dist / v_ref * factor). Set factor <= 0 to disable.
  double expl_pursuit_timeout_factor{10.0};
  double expl_pursuit_timeout_v_ref{0.5};
  double expl_pursuit_timeout_min_sec{10.0};
  // Invalidation keep-out — drop fresh clusters that fall within radius_m of
  // any INVALIDATED record whose invalidation is still inside the cooldown
  // window. Set radius_m <= 0 to disable; cooldown_sec <= 0 = permanent.
  double expl_invalidation_keep_out_radius_m{1.5};
  double expl_invalidation_cooldown_sec{30.0};
  // Persistent visited bitmap (suppresses re-detection of revisited frontiers)
  double expl_visited_map_center_x{0.0};
  double expl_visited_map_center_y{0.0};
  double expl_visited_map_width_m{100.0};
  double expl_visited_map_height_m{100.0};
  double expl_visited_map_resolution_m{0.15};
  bool   expl_publish_visited_map{true};
  // When true, fuse persistent visited_map_ values into UNKNOWN cells of the
  // freshly received local OccupancyGrid before it is consumed by HGP / the
  // local optimizer / the frontier detector. This makes revisited regions
  // immediately come back with their last-known FREE/OCCUPIED state instead
  // of one-frame UNKNOWN flicker. Static-environment only — re-introduces
  // stale OCCUPIED for moving obstacles that have since left.
  bool   expl_fuse_persistent_into_local{true};
  // When true, run the frontier detector on the persistent visited_map_
  // (entire mission history) instead of the freshly received local sliding
  // window. Surfaces frontiers at the boundary of explored area no matter
  // where the robot currently is — fixes "revisited corridor doesn't show
  // its far-end frontier" symptom. Inherits the phantom-OCCUPIED caveat of
  // the persistent map for dynamic environments.
  bool   expl_detect_on_visited_map{true};
  // MinPos multi-robot exploration
  bool   expl_use_minpos{false};              // enable rank-based peer-aware allocation
  double expl_peer_timeout_sec{5.0};          // drop peer after this silence (seconds)
  double expl_peer_publish_rate_hz{5.0};      // throttle for pose broadcast (Hz)
  double expl_min_frontier_dist_to_peers_m{0.0};  // reject frontier candidates within this radius of any active peer; 0 disables
  double expl_peer_visit_radius_m{2.0};       // mark frontier VISITED when any active peer is within this radius (sticky); 0 disables
  // Visualization
  bool   expl_publish_markers{true};
  // When true, skip publishing exploration goals from exploreSelectCallback so
  // an external selector (e.g. a VLM node) can own term_goal. Frontier
  // detection, scoring, and marker publication continue as normal.
  bool   expl_external_selector{false};
};

struct BasisConverter {
  Eigen::Matrix<double, 4, 4> A_pos_mv_rest;
  Eigen::Matrix<double, 4, 4> A_pos_mv_rest_inv;
  Eigen::Matrix<double, 3, 3> A_vel_mv_rest;
  Eigen::Matrix<double, 3, 3> A_vel_mv_rest_inv;
  Eigen::Matrix<double, 2, 2> A_accel_mv_rest;
  Eigen::Matrix<double, 2, 2> A_accel_mv_rest_inv;
  Eigen::Matrix<double, 4, 4> A_pos_be_rest;
  Eigen::Matrix<double, 4, 4> A_pos_bs_seg0, A_pos_bs_seg1, A_pos_bs_rest, A_pos_bs_seg_last2,
      A_pos_bs_seg_last;

  Eigen::Matrix<double, 4, 4> M_pos_bs2mv_seg0, M_pos_bs2mv_seg1, M_pos_bs2mv_rest,
      M_pos_bs2mv_seg_last2, M_pos_bs2mv_seg_last;

  Eigen::Matrix<double, 4, 4> M_pos_bs2be_seg0, M_pos_bs2be_seg1, M_pos_bs2be_rest,
      M_pos_bs2be_seg_last2, M_pos_bs2be_seg_last;

  Eigen::Matrix<double, 3, 3> M_vel_bs2mv_seg0, M_vel_bs2mv_rest, M_vel_bs2mv_seg_last;
  Eigen::Matrix<double, 3, 3> M_vel_bs2be_seg0, M_vel_bs2be_rest, M_vel_bs2be_seg_last;

  BasisConverter() {
    // See matlab.
    // This is for t \in [0 1];

    //////MATRICES A FOR MINVO POSITION///////// (there is only one)
    A_pos_mv_rest << -3.4416308968564117698463178385282, 6.9895481477801393310755884158425,
        -4.4622887507045296828778191411402, 0.91437149978080234369315348885721,
        6.6792587327074839365081970754545, -11.845989901556746914934592496138,
        5.2523596690684613008670567069203, 0, -6.6792587327074839365081970754545,
        8.1917862965657040064115790301003, -1.5981560640774179482548333908198,
        0.085628500219197656306846511142794, 3.4416308968564117698463178385282,
        -3.3353445427890959784633650997421, 0.80808514571348655231020075007109,
        -0.0000000000000000084567769453869345852581318467855;

    //////INVERSE OF A_pos_mv_rest
    A_pos_mv_rest_inv = A_pos_mv_rest.inverse();

    //////MATRICES A FOR MINVO VELOCITY///////// (there is only one)
    A_vel_mv_rest << 1.4999999992328318931811281800037, -2.3660254034601951866889635311964,
        0.9330127021136816189983420599674, -2.9999999984656637863622563600074,
        2.9999999984656637863622563600074, 0, 1.4999999992328318931811281800037,
        -0.6339745950054685996732928288111, 0.066987297886318325490506708774774;

    //////INVERSE OF A_vel_mv_rest
    A_vel_mv_rest_inv = A_vel_mv_rest.inverse();

    //////MATRICES A FOR MINVO ACCELERATION///////// (there is only one)
    A_accel_mv_rest << -1.0, 1.0, 1.0, 0.0;

    /////INVERSE OF A_accel_mv_rest
    A_accel_mv_rest_inv = A_accel_mv_rest.inverse();

    //////MATRICES A FOR Bezier POSITION///////// (there is only one)
    A_pos_be_rest <<

        -1.0,
        3.0, -3.0, 1.0, 3.0, -6.0, 3.0, 0, -3.0, 3.0, 0, 0, 1.0, 0, 0, 0;

    //////MATRICES A FOR BSPLINE POSITION/////////
    A_pos_bs_seg0 <<

        -1.0000,
        3.0000, -3.0000, 1.0000, 1.7500, -4.5000, 3.0000, 0, -0.9167, 1.5000, 0, 0, 0.1667, 0, 0, 0;

    A_pos_bs_seg1 <<

        -0.2500,
        0.7500, -0.7500, 0.2500, 0.5833, -1.2500, 0.2500, 0.5833, -0.5000, 0.5000, 0.5000, 0.1667,
        0.1667, 0, 0, 0;

    A_pos_bs_rest <<

        -0.1667,
        0.5000, -0.5000, 0.1667, 0.5000, -1.0000, 0, 0.6667, -0.5000, 0.5000, 0.5000, 0.1667,
        0.1667, 0, 0, 0;

    A_pos_bs_seg_last2 <<

        -0.1667,
        0.5000, -0.5000, 0.1667, 0.5000, -1.0000, 0.0000, 0.6667, -0.5833, 0.5000, 0.5000, 0.1667,
        0.2500, 0, 0, 0;

    A_pos_bs_seg_last <<

        -0.1667,
        0.5000, -0.5000, 0.1667, 0.9167, -1.2500, -0.2500, 0.5833, -1.7500, 0.7500, 0.7500, 0.2500,
        1.0000, 0, 0, 0;

    //////BSPLINE to MINVO POSITION/////////

    M_pos_bs2mv_seg0 <<

        1.1023313949144333268037598827505,
        0.34205724556666972091534262290224, -0.092730934245582874453361910127569,
        -0.032032766697130621302846975595457, -0.049683556253749178166501110354147,
        0.65780347324677179710050722860615, 0.53053863760186903419935333658941,
        0.21181027098212013015654520131648, -0.047309044211162346038612724896666,
        0.015594436894155586093013710069499, 0.5051827557159349613158383363043,
        0.63650059656260427054519368539331, -0.0053387944495217444854096022766043,
        -0.015455155707597083292181849856206, 0.057009540927778303009976212933907,
        0.18372189915240558222286892942066;

    M_pos_bs2mv_seg1 <<

        0.27558284872860833170093997068761,
        0.085514311391667430228835655725561, -0.023182733561395718613340477531892,
        -0.0080081916742826553257117438988644, 0.6099042761975865811763242163579,
        0.63806904207840509091198555324809, 0.29959938009132258684985572472215,
        0.12252106674808682651445224109921, 0.11985166952332682033244282138185,
        0.29187180223752445806795208227413, 0.66657381254229419731416328431806,
        0.70176522577378930289881964199594, -0.0053387944495217444854096022766043,
        -0.015455155707597083292181849856206, 0.057009540927778303009976212933907,
        0.18372189915240558222286892942066;

    M_pos_bs2mv_rest <<

        0.18372189915240555446729331379174,
        0.057009540927778309948870116841135, -0.015455155707597117986651369392348,
        -0.0053387944495218164764338553140988, 0.70176522577378919187651717948029,
        0.66657381254229419731416328431806, 0.29187180223752384744528853843804,
        0.11985166952332582113172065874096, 0.11985166952332682033244282138185,
        0.29187180223752445806795208227413, 0.66657381254229419731416328431806,
        0.70176522577378930289881964199594, -0.0053387944495217444854096022766043,
        -0.015455155707597083292181849856206, 0.057009540927778303009976212933907,
        0.18372189915240558222286892942066;

    M_pos_bs2mv_seg_last2 <<

        0.18372189915240569324517139193631,
        0.057009540927778309948870116841135, -0.015455155707597145742226985021261,
        -0.0053387944495218164764338553140988, 0.70176522577378952494342456702725,
        0.66657381254229453038107067186502, 0.29187180223752412500104469472717,
        0.11985166952332593215402312125661, 0.1225210667480875342816304396365,
        0.29959938009132280889446064975346, 0.63806904207840497988968309073243,
        0.60990427619758624810941682881094, -0.0080081916742826154270717964323012,
        -0.023182733561395621468825822830695, 0.085514311391667444106623463540018,
        0.27558284872860833170093997068761;

    M_pos_bs2mv_seg_last <<

        0.18372189915240555446729331379174,
        0.057009540927778309948870116841135, -0.015455155707597117986651369392348,
        -0.0053387944495218164764338553140988, 0.63650059656260415952289122287766,
        0.5051827557159349613158383363043, 0.015594436894155294659469745965907,
        -0.047309044211162887272337229660479, 0.21181027098212068526805751389475,
        0.53053863760186914522165579910506, 0.65780347324677146403359984105919,
        -0.049683556253749622255710960416764, -0.032032766697130461708287185729205,
        -0.09273093424558248587530329132278, 0.34205724556666977642649385416007,
        1.1023313949144333268037598827505;

    //////BSPLINE to BEZIER POSITION/////////

    M_pos_bs2be_seg0 <<

        1.0000,
        0.0000, -0.0000, 0, 0, 1.0000, 0.5000, 0.2500, 0, -0.0000, 0.5000, 0.5833, 0, 0, 0, 0.1667;

    M_pos_bs2be_seg1 <<

        0.2500,
        0.0000, -0.0000, 0, 0.5833, 0.6667, 0.3333, 0.1667, 0.1667, 0.3333, 0.6667, 0.6667, 0, 0, 0,
        0.1667;

    M_pos_bs2be_rest <<

        0.1667,
        0.0000, 0, 0, 0.6667, 0.6667, 0.3333, 0.1667, 0.1667, 0.3333, 0.6667, 0.6667, 0, 0, 0,
        0.1667;

    M_pos_bs2be_seg_last2 <<

        0.1667,
        0, -0.0000, 0, 0.6667, 0.6667, 0.3333, 0.1667, 0.1667, 0.3333, 0.6667, 0.5833, 0, 0, 0,
        0.2500;

    M_pos_bs2be_seg_last <<

        0.1667,
        0.0000, 0, 0, 0.5833, 0.5000, 0, 0, 0.2500, 0.5000, 1.0000, 0, 0, 0, 0, 1.0000;

    /////BSPLINE to MINVO VELOCITY
    M_vel_bs2mv_seg0 <<

        1.077349059083916,
        0.1666702138890985, -0.07735049175615138, -0.03867488648729411, 0.7499977187062712,
        0.5386802643920123, -0.03867417280506149, 0.08333206631563977, 0.538670227146185;

    M_vel_bs2mv_rest <<

        0.538674529541958,
        0.08333510694454926, -0.03867524587807569, 0.4999996430546639, 0.8333328256508203,
        0.5000050185139366, -0.03867417280506149, 0.08333206631563977, 0.538670227146185;

    M_vel_bs2mv_seg_last <<

        0.538674529541958,
        0.08333510694454926, -0.03867524587807569, 0.5386738158597254, 0.7500007593351806,
        -0.03866520863224832, -0.07734834561012298, 0.1666641326312795, 1.07734045429237;

    /////BSPLINE to BEZIER VELOCITY
    M_vel_bs2be_seg0 <<

        1.0000,
        0, 0, 0, 1.0000, 0.5000, 0, 0, 0.5000;

    M_vel_bs2be_rest <<

        0.5000,
        0, 0, 0.5000, 1.0000, 0.5000, 0, 0, 0.5000;

    M_vel_bs2be_seg_last <<

        0.5000,
        0, 0, 0.5000, 1.0000, 0, 0, 0, 1.0000;
  }

  //////MATRIX A FOR MINVO POSITION/////////
  Eigen::Matrix<double, 4, 4> getArestMinvo() { return A_pos_mv_rest; }
  //////MATRIX A FOR Bezier POSITION/////////
  Eigen::Matrix<double, 4, 4> getArestBezier() { return A_pos_be_rest; }

  //////MATRIX A FOR BSPLINE POSITION/////////
  Eigen::Matrix<double, 4, 4> getArestBSpline() { return A_pos_bs_rest; }

  //////MATRICES A FOR MINVO POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getAMinvo(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> A_pos_mv;  // will have as many elements as num_pol
    for (int i = 0; i < num_pol; i++) {
      A_pos_mv.push_back(A_pos_mv_rest);
    }
    return A_pos_mv;
  }

  //////MATRICES A FOR Bezier POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getABezier(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> A_pos_be;  // will have as many elements as num_pol
    for (int i = 0; i < num_pol; i++) {
      A_pos_be.push_back(A_pos_be_rest);
    }
    return A_pos_be;
  }

  //////MATRICES A FOR BSPLINE POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getABSpline(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> A_pos_bs;  // will have as many elements as num_pol
    A_pos_bs.push_back(A_pos_bs_seg0);
    A_pos_bs.push_back(A_pos_bs_seg1);
    for (int i = 0; i < (num_pol - 4); i++) {
      A_pos_bs.push_back(A_pos_bs_rest);
    }
    A_pos_bs.push_back(A_pos_bs_seg_last2);
    A_pos_bs.push_back(A_pos_bs_seg_last);
    return A_pos_bs;
  }

  //////BSPLINE to MINVO POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getMinvoPosConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2mv;  // will have as many elements as num_pol
    M_pos_bs2mv.push_back(M_pos_bs2mv_seg0);
    M_pos_bs2mv.push_back(M_pos_bs2mv_seg1);
    for (int i = 0; i < (num_pol - 4); i++) {
      M_pos_bs2mv.push_back(M_pos_bs2mv_rest);
    }
    M_pos_bs2mv.push_back(M_pos_bs2mv_seg_last2);
    M_pos_bs2mv.push_back(M_pos_bs2mv_seg_last);
    return M_pos_bs2mv;
  }

  //////BEZIER to MINVO POSITION/////////
  //////Q_{MINVO} = M_{BEZIER2MINVO} * Q_{BEZIER}
  Eigen::Matrix<double, 4, 4> getMinvoPosConverterFromBezier() {
    // Compute the conversion matrix for one segment
    Eigen::Matrix<double, 4, 4> M_be2mv = A_pos_mv_rest_inv * A_pos_be_rest;

    return M_be2mv;
  }

  //////BSPLINE to BEZIER POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getBezierPosConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2be;  // will have as many elements as num_pol
    M_pos_bs2be.push_back(M_pos_bs2be_seg0);
    M_pos_bs2be.push_back(M_pos_bs2be_seg1);
    for (int i = 0; i < (num_pol - 4); i++) {
      M_pos_bs2be.push_back(M_pos_bs2be_rest);
    }
    M_pos_bs2be.push_back(M_pos_bs2be_seg_last2);
    M_pos_bs2be.push_back(M_pos_bs2be_seg_last);
    return M_pos_bs2be;
  }

  //////BSPLINE to BSPLINE POSITION/////////
  std::vector<Eigen::Matrix<double, 4, 4>> getBSplinePosConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2bs;  // will have as many elements as num_pol
    for (int i = 0; i < num_pol; i++) {
      M_pos_bs2bs.push_back(Eigen::Matrix<double, 4, 4>::Identity());
    }
    return M_pos_bs2bs;
  }

  //////BSPLINE to MINVO Velocity/////////
  std::vector<Eigen::Matrix<double, 3, 3>> getMinvoVelConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 3, 3>> M_vel_bs2mv;  // will have as many elements as num_pol
    M_vel_bs2mv.push_back(M_vel_bs2mv_seg0);
    for (int i = 0; i < (num_pol - 2 - 1); i++) {
      M_vel_bs2mv.push_back(M_vel_bs2mv_rest);
    }
    M_vel_bs2mv.push_back(M_vel_bs2mv_seg_last);
    return M_vel_bs2mv;
  }

  //////BSPLINE to BEZIER Velocity/////////
  std::vector<Eigen::Matrix<double, 3, 3>> getBezierVelConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 3, 3>> M_vel_bs2be;  // will have as many elements as segments
    M_vel_bs2be.push_back(M_vel_bs2be_seg0);
    for (int i = 0; i < (num_pol - 2 - 1); i++) {
      M_vel_bs2be.push_back(M_vel_bs2be_rest);
    }
    M_vel_bs2be.push_back(M_vel_bs2be_seg_last);
    return M_vel_bs2be;
  }

  //////BSPLINE to BSPLINE Velocity/////////
  std::vector<Eigen::Matrix<double, 3, 3>> getBSplineVelConverters(int num_pol) {
    std::vector<Eigen::Matrix<double, 3, 3>> M_vel_bs2bs;  // will have as many elements as num_pol
    for (int i = 0; i < num_pol; i++) {
      M_vel_bs2bs.push_back(Eigen::Matrix<double, 3, 3>::Identity());
    }
    return M_vel_bs2bs;
  }
};

struct PieceWisePol {
  // Interval 0: t\in[t0, t1)
  // Interval 1: t\in[t1, t2)
  // Interval 2: t\in[t2, t3)
  //...
  // Interval n-1: t\in[tn, tn+1)

  // n intervals in total

  // times has n+1 elements
  std::vector<double> times;  // [t0,t1,t2,...,tn+1]

  // coefficients has n elements
  // The coeffients are such that pol(t)=coeff_of_that_interval*[u^3 u^2 u 1]
  // with u=(t-t_min_that_interval)/(t_max_that_interval- t_min_that_interval)
  std::vector<Eigen::Matrix<double, 4, 1>> coeff_x;  // [a b c d]' of Int0 , [a b c d]' of Int1,...
  std::vector<Eigen::Matrix<double, 4, 1>> coeff_y;  // [a b c d]' of Int0 , [a b c d]' of Int1,...
  std::vector<Eigen::Matrix<double, 4, 1>> coeff_z;  // [a b c d]' of Int0 , [a b c d]' of Int1,...

  void clear() {
    times.clear();
    coeff_x.clear();
    coeff_y.clear();
    coeff_z.clear();
  }

  // Get the end time of the trajectory
  double getEndTime() const { return times.back(); }

  inline Eigen::Vector3d eval(double t) const {
    Eigen::Vector3d result;

    // return the last value of the polynomial in the last interval
    if (t >= times.back()) {
      Eigen::Matrix<double, 4, 1> tmp;
      // double u = 1;
      double u = times.back() - times[times.size() - 2];
      tmp << u * u * u, u * u, u, 1.0;
      result.x() = coeff_x.back().transpose() * tmp;
      result.y() = coeff_y.back().transpose() * tmp;
      result.z() = coeff_z.back().transpose() * tmp;
      return result;
    }

    // return the first value of the polynomial in the first interval
    if (t < times.front()) {
      Eigen::Matrix<double, 4, 1> tmp;
      double u = 0;
      tmp << u * u * u, u * u, u, 1.0;
      result.x() = coeff_x.front().transpose() * tmp;
      result.y() = coeff_y.front().transpose() * tmp;
      result.z() = coeff_z.front().transpose() * tmp;
      return result;
    }

    // Find the interval where t is
    //(times - 1) is the number of intervals
    for (int i = 0; i < (times.size() - 1); i++) {
      if (times[i] <= t && t < times[i + 1]) {
        // double u = (t - times[i]) / (times[i + 1] - times[i]);
        double u = t - times[i];

        Eigen::Matrix<double, 4, 1> tmp;
        tmp << u * u * u, u * u, u, 1.0;

        result.x() = coeff_x[i].transpose() * tmp;
        result.y() = coeff_y[i].transpose() * tmp;
        result.z() = coeff_z[i].transpose() * tmp;

        break;
      }
    }
    return result;
  }

  inline Eigen::Vector3d velocity(double t) const {
    Eigen::Vector3d vel;

    // Handle the case where t is after the last interval
    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      vel.x() = 3 * coeff_x.back()(0) * u * u + 2 * coeff_x.back()(1) * u + coeff_x.back()(2);
      vel.y() = 3 * coeff_y.back()(0) * u * u + 2 * coeff_y.back()(1) * u + coeff_y.back()(2);
      vel.z() = 3 * coeff_z.back()(0) * u * u + 2 * coeff_z.back()(1) * u + coeff_z.back()(2);
      return vel;
    }

    // Handle the case where t is before the first interval
    if (t < times.front()) {
      vel.x() = coeff_x.front()(2);
      vel.y() = coeff_y.front()(2);
      vel.z() = coeff_z.front()(2);
      return vel;
    }

    // Find the interval where t lies and calculate velocity
    for (int i = 0; i < (times.size() - 1); i++) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        vel.x() = 3 * coeff_x[i](0) * u * u + 2 * coeff_x[i](1) * u + coeff_x[i](2);
        vel.y() = 3 * coeff_y[i](0) * u * u + 2 * coeff_y[i](1) * u + coeff_y[i](2);
        vel.z() = 3 * coeff_z[i](0) * u * u + 2 * coeff_z[i](1) * u + coeff_z[i](2);
        break;
      }
    }

    return vel;
  }

  /// Evaluate acceleration (second derivative) at time t
  inline Eigen::Vector3d acceleration(double t) const {
    Eigen::Vector3d a{0, 0, 0};

    // 1) If t is after the last interval:
    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      // coeff_* .back()(0) is the cubic a, .back()(1) is the quadratic b
      a.x() = 6.0 * coeff_x.back()(0) * u + 2.0 * coeff_x.back()(1);
      a.y() = 6.0 * coeff_y.back()(0) * u + 2.0 * coeff_y.back()(1);
      a.z() = 6.0 * coeff_z.back()(0) * u + 2.0 * coeff_z.back()(1);
      return a;
    }

    // 2) If t is before the first interval:
    if (t < times.front()) {
      a.x() = 2.0 * coeff_x.front()(1);
      a.y() = 2.0 * coeff_y.front()(1);
      a.z() = 2.0 * coeff_z.front()(1);
      return a;
    }

    // 3) Otherwise, find the correct interval i
    for (size_t i = 0; i + 1 < times.size(); ++i) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        a.x() = 6.0 * coeff_x[i](0) * u + 2.0 * coeff_x[i](1);
        a.y() = 6.0 * coeff_y[i](0) * u + 2.0 * coeff_y[i](1);
        a.z() = 6.0 * coeff_z[i](0) * u + 2.0 * coeff_z[i](1);
        break;
      }
    }
    return a;
  }

  void print() const {
    std::cout << "coeff_x.size()= " << coeff_x.size() << std::endl;
    std::cout << "times.size()= " << times.size() << std::endl;
    std::cout << "Note that coeff_x.size() == times.size()-1" << std::endl;

    for (int i = 0; i < times.size(); i++) {
      printf("Time: %f\n", times[i]);
    }

    for (int i = 0; i < (times.size() - 1); i++) {
      std::cout << "From " << times[i] << " to " << times[i + 1] << std::endl;
      std::cout << "  Coeff_x= " << coeff_x[i].transpose() << std::endl;
      std::cout << "  Coeff_y= " << coeff_y[i].transpose() << std::endl;
      std::cout << "  Coeff_z= " << coeff_z[i].transpose() << std::endl;
    }
  }
};

struct PieceWiseQuinticPol {
  // breakpoints [t0, t1, …, tN]
  std::vector<double> times;

  // per‐segment coefficients [a b c d e f]ᵀ for a·u⁵ + b·u⁴ + c·u³ + d·u² + e·u + f
  std::vector<Eigen::Matrix<double, 6, 1>> coeff_x, coeff_y, coeff_z;

  /// Remove all intervals
  void clear() {
    times.clear();
    coeff_x.clear();
    coeff_y.clear();
    coeff_z.clear();
  }

  /// End time of the trajectory
  double getEndTime() const { return times.back(); }

  /// Position p = a·u⁵ + b·u⁴ + c·u³ + d·u² + e·u + f
  inline Eigen::Vector3d eval(double t) const {
    Eigen::Vector3d r{0, 0, 0};

    // (1) after last: clamp to end of last interval
    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      Eigen::Matrix<double, 6, 1> U;
      U << u * u * u * u * u, u * u * u * u, u * u * u, u * u, u, 1.0;
      r.x() = coeff_x.back().dot(U);
      r.y() = coeff_y.back().dot(U);
      r.z() = coeff_z.back().dot(U);
      return r;
    }

    // (2) before first: u=0
    if (t < times.front()) {
      Eigen::Matrix<double, 6, 1> U;
      U << 0, 0, 0, 0, 0, 1.0;
      r.x() = coeff_x.front().dot(U);
      r.y() = coeff_y.front().dot(U);
      r.z() = coeff_z.front().dot(U);
      return r;
    }

    // (3) somewhere in the middle
    for (size_t i = 0; i + 1 < times.size(); ++i) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        Eigen::Matrix<double, 6, 1> U;
        U << u * u * u * u * u, u * u * u * u, u * u * u, u * u, u, 1.0;
        r.x() = coeff_x[i].dot(U);
        r.y() = coeff_y[i].dot(U);
        r.z() = coeff_z[i].dot(U);
        break;
      }
    }
    return r;
  }

  /// Velocity p′ = 5a·u⁴ + 4b·u³ + 3c·u² + 2d·u + e
  inline Eigen::Vector3d velocity(double t) const {
    Eigen::Vector3d v{0, 0, 0};

    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      v.x() = 5 * coeff_x.back()(0) * u * u * u * u + 4 * coeff_x.back()(1) * u * u * u +
              3 * coeff_x.back()(2) * u * u + 2 * coeff_x.back()(3) * u + coeff_x.back()(4);
      v.y() = 5 * coeff_y.back()(0) * u * u * u * u + 4 * coeff_y.back()(1) * u * u * u +
              3 * coeff_y.back()(2) * u * u + 2 * coeff_y.back()(3) * u + coeff_y.back()(4);
      v.z() = 5 * coeff_z.back()(0) * u * u * u * u + 4 * coeff_z.back()(1) * u * u * u +
              3 * coeff_z.back()(2) * u * u + 2 * coeff_z.back()(3) * u + coeff_z.back()(4);
      return v;
    }

    if (t < times.front()) {
      v.x() = coeff_x.front()(4);
      v.y() = coeff_y.front()(4);
      v.z() = coeff_z.front()(4);
      return v;
    }

    for (size_t i = 0; i + 1 < times.size(); ++i) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        v.x() = 5 * coeff_x[i](0) * u * u * u * u + 4 * coeff_x[i](1) * u * u * u +
                3 * coeff_x[i](2) * u * u + 2 * coeff_x[i](3) * u + coeff_x[i](4);
        v.y() = 5 * coeff_y[i](0) * u * u * u * u + 4 * coeff_y[i](1) * u * u * u +
                3 * coeff_y[i](2) * u * u + 2 * coeff_y[i](3) * u + coeff_y[i](4);
        v.z() = 5 * coeff_z[i](0) * u * u * u * u + 4 * coeff_z[i](1) * u * u * u +
                3 * coeff_z[i](2) * u * u + 2 * coeff_z[i](3) * u + coeff_z[i](4);
        break;
      }
    }
    return v;
  }

  /// Acceleration p″ = 20a·u³ + 12b·u² + 6c·u + 2d
  inline Eigen::Vector3d acceleration(double t) const {
    Eigen::Vector3d a{0, 0, 0};

    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      a.x() = 20 * coeff_x.back()(0) * u * u * u + 12 * coeff_x.back()(1) * u * u +
              6 * coeff_x.back()(2) * u + 2 * coeff_x.back()(3);
      a.y() = 20 * coeff_y.back()(0) * u * u * u + 12 * coeff_y.back()(1) * u * u +
              6 * coeff_y.back()(2) * u + 2 * coeff_y.back()(3);
      a.z() = 20 * coeff_z.back()(0) * u * u * u + 12 * coeff_z.back()(1) * u * u +
              6 * coeff_z.back()(2) * u + 2 * coeff_z.back()(3);
      return a;
    }

    if (t < times.front()) {
      a.x() = 2 * coeff_x.front()(3);
      a.y() = 2 * coeff_y.front()(3);
      a.z() = 2 * coeff_z.front()(3);
      return a;
    }

    for (size_t i = 0; i + 1 < times.size(); ++i) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        a.x() = 20 * coeff_x[i](0) * u * u * u + 12 * coeff_x[i](1) * u * u +
                6 * coeff_x[i](2) * u + 2 * coeff_x[i](3);
        a.y() = 20 * coeff_y[i](0) * u * u * u + 12 * coeff_y[i](1) * u * u +
                6 * coeff_y[i](2) * u + 2 * coeff_y[i](3);
        a.z() = 20 * coeff_z[i](0) * u * u * u + 12 * coeff_z[i](1) * u * u +
                6 * coeff_z[i](2) * u + 2 * coeff_z[i](3);
        break;
      }
    }
    return a;
  }

  /// Jerk p‴ = 60a·u² + 24b·u + 6c
  inline Eigen::Vector3d jerk(double t) const {
    Eigen::Vector3d j{0, 0, 0};

    if (t >= times.back()) {
      double u = times.back() - times[times.size() - 2];
      j.x() = 60 * coeff_x.back()(0) * u * u + 24 * coeff_x.back()(1) * u + 6 * coeff_x.back()(2);
      j.y() = 60 * coeff_y.back()(0) * u * u + 24 * coeff_y.back()(1) * u + 6 * coeff_y.back()(2);
      j.z() = 60 * coeff_z.back()(0) * u * u + 24 * coeff_z.back()(1) * u + 6 * coeff_z.back()(2);
      return j;
    }

    if (t < times.front()) {
      j.x() = 6 * coeff_x.front()(2);
      j.y() = 6 * coeff_y.front()(2);
      j.z() = 6 * coeff_z.front()(2);
      return j;
    }

    for (size_t i = 0; i + 1 < times.size(); ++i) {
      if (times[i] <= t && t < times[i + 1]) {
        double u = t - times[i];
        j.x() = 60 * coeff_x[i](0) * u * u + 24 * coeff_x[i](1) * u + 6 * coeff_x[i](2);
        j.y() = 60 * coeff_y[i](0) * u * u + 24 * coeff_y[i](1) * u + 6 * coeff_y[i](2);
        j.z() = 60 * coeff_z[i](0) * u * u + 24 * coeff_z[i](1) * u + 6 * coeff_z[i](2);
        break;
      }
    }
    return j;
  }

  /// Print internal data
  void print() const {
    std::cout << "PieceWiseQuinticPol: " << times.size() - 1 << " segments\n";
    for (size_t i = 0; i + 1 < times.size(); ++i) {
      std::cout << " [" << times[i] << "," << times[i + 1] << "): "
                << "cx=" << coeff_x[i].transpose() << "  cy=" << coeff_y[i].transpose()
                << "  cz=" << coeff_z[i].transpose() << "\n";
    }
  }
};

struct dynTraj {
  /// Which representation to use
  enum class Mode { Piecewise, Quintic, Analytic } mode{Mode::Analytic};

  // --- piecewise quintic branch ---
  PieceWiseQuinticPol pwp;

  // --- single‐segment quintic branch ---
  double poly_start_time = 0.0;
  double poly_end_time = 0.0;
  Eigen::Matrix<double, 6, 1> cx, cy, cz;

  // --- analytic expression branch ---
  std::string traj_x, traj_y, traj_z;
  std::string traj_vx, traj_vy, traj_vz;  // optional (velocity expressions)
  double t_var{0.0};
  exprtk::symbol_table<double> symbol_table;
  exprtk::expression<double> expr_x, expr_y, expr_z;
  exprtk::expression<double> expr_vx, expr_vy, expr_vz;
  bool analytic_compiled{false};

  // shared metadata
  Eigen::Vector3d ekf_cov_p;
  Eigen::Vector3d ekf_cov_q;
  Eigen::Vector3d poly_cov;
  std::vector<Eigen::Matrix<double, 3, 6>> control_points;
  Eigen::Vector3d bbox;
  Eigen::Vector3d goal;
  bool is_agent = false;
  int id = -1;
  double time_received = 0.0;
  double tracking_utility = 0.0;
  double communication_delay = 0.0;

  dynTraj() = default;

  /// Quintic‐based constructor for quick tests
  dynTraj(const Eigen::Vector3d& x0, const Eigen::Vector3d& v0, const Eigen::Vector3d& a0,
          const Eigen::Vector3d& xf, const Eigen::Vector3d& vf, const Eigen::Vector3d& af,
          double poly_start_time_, double poly_end_time_)
      : mode(Mode::Quintic), poly_start_time(poly_start_time_), poly_end_time(poly_end_time_) {
    double duration = poly_end_time - poly_start_time;
    if (duration <= 0.0) {
      std::cerr << "Error: Invalid trajectory time range." << std::endl;
      cx.setZero();
      cy.setZero();
      cz.setZero();
    } else {
      lbfgs_solver_utils::fit_quintic(x0, v0, a0, xf, vf, af, duration, cx, cy, cz);
    }
  }

  /// Switch to a piecewise cubic representation
  inline void setPiecewise(const PieceWiseQuinticPol& poly) {
    mode = Mode::Piecewise;
    pwp = poly;
  }

  bool compileAnalytic() {
    symbol_table.clear();
    symbol_table.add_variable("t", t_var);
    symbol_table.add_constants();

    auto reg = [&](exprtk::expression<double>& e) { e.register_symbol_table(symbol_table); };
    reg(expr_x);
    reg(expr_y);
    reg(expr_z);
    reg(expr_vx);
    reg(expr_vy);
    reg(expr_vz);

    exprtk::parser<double> parser;

    auto compile_one = [&](const std::string& label, const std::string& src,
                           exprtk::expression<double>& expr) -> bool {
      if (src.empty()) {
        if (label == "traj_x" || label == "traj_y" || label == "traj_z") {
          std::cerr << "Missing required analytic expression " << label << "\n";
          return false;
        }
        // otherwise it was a velocity string → OK to skip
        return true;
      }

      if (!parser.compile(src, expr)) {
        std::ostringstream oss;
        oss << "ExprTk compile failure (" << label << "): '" << src << "' errors:";
        for (std::size_t i = 0; i < parser.error_count(); ++i) {
          auto e = parser.get_error(i);
          oss << " [pos " << e.token.position << " type " << exprtk::parser_error::to_str(e.mode)
              << " msg '" << e.diagnostic << "']";
        }
        std::cerr << oss.str() << std::endl;
        return false;
      }
      return true;
    };

    bool ok = true;
    ok &= compile_one("traj_x", traj_x, expr_x);
    ok &= compile_one("traj_y", traj_y, expr_y);
    ok &= compile_one("traj_z", traj_z, expr_z);
    ok &= compile_one("traj_vx", traj_vx, expr_vx);
    ok &= compile_one("traj_vy", traj_vy, expr_vy);
    ok &= compile_one("traj_vz", traj_vz, expr_vz);

    analytic_compiled = ok;
    return ok;
  }

  /// Evaluate position at time t
  inline Eigen::Vector3d eval(double t) const {
    switch (mode) {
      case Mode::Piecewise:
        return pwp.eval(t);
      case Mode::Quintic:
        return evalQuinticPos(t);
      case Mode::Analytic:
        return evalAnalyticPos(t);
    }
    return Eigen::Vector3d::Zero();
  }

  // p(τ) = a τ^5 + b τ^4 + c τ^3 + d τ^2 + e τ + f  with c = [a,b,c,d,e,f]
  static inline double poly5_abs(const Eigen::Matrix<double, 6, 1>& c, double tau) {
    double v = c(0);
    v = v * tau + c(1);
    v = v * tau + c(2);
    v = v * tau + c(3);
    v = v * tau + c(4);
    v = v * tau + c(5);
    return v;
  }

  // p'(τ) = 5 a τ^4 + 4 b τ^3 + 3 c τ^2 + 2 d τ + e
  static inline double dpoly5_abs(const Eigen::Matrix<double, 6, 1>& c, double tau) {
    double v = 5 * c(0);
    v = v * tau + 4 * c(1);
    v = v * tau + 3 * c(2);
    v = v * tau + 2 * c(3);
    v = v * tau + c(4);
    return v;
  }

  // p''(τ) = 20 a τ^3 + 12 b τ^2 + 6 c τ + 2 d
  static inline double ddpoly5_abs(const Eigen::Matrix<double, 6, 1>& c, double tau) {
    double v = 20 * c(0);
    v = v * tau + 12 * c(1);
    v = v * tau + 6 * c(2);
    v = v * tau + 2 * c(3);
    return v;
  }

  // --- normalized-time evaluation: u = (t - t0) / (tf - t0) clamped to [0,1] ---

  inline Eigen::Vector3d evalQuinticPos(double t) const {
    const double tau = t - poly_start_time;  // absolute time from segment start
    return {poly5_abs(cx, tau), poly5_abs(cy, tau), poly5_abs(cz, tau)};
  }

  inline Eigen::Vector3d evalAnalyticPos(double t) const {
    if (!analytic_compiled) {
      // this should never happen if steps 1+2 are correct
      std::cerr << "[dynTraj] evalAnalyticPos called but analytic_compiled==false\n";
      return Eigen::Vector3d::Zero();
    }

    const_cast<dynTraj*>(this)->t_var = t;
    return {expr_x.value(), expr_y.value(), expr_z.value()};
  }

  /// Evaluate velocity at time t
  inline Eigen::Vector3d velocity(double t) const {
    switch (mode) {
      case Mode::Piecewise:
        return pwp.velocity(t);
      case Mode::Quintic:
        return velocityQuintic(t);
      case Mode::Analytic:
        return velocityAnalytic(t);
    }
    return Eigen::Vector3d::Zero();
  }

  inline Eigen::Vector3d velocityQuintic(double t) const {
    const double tau = t - poly_start_time;
    // derivative w.r.t. absolute time (NO 1/duration factor)
    return {dpoly5_abs(cx, tau), dpoly5_abs(cy, tau), dpoly5_abs(cz, tau)};
  }

  inline Eigen::Vector3d velocityAnalytic(double t) const {
    if (!analytic_compiled) return Eigen::Vector3d::Zero();
    const_cast<dynTraj*>(this)->t_var = t;
    // If velocity expressions provided
    if (!traj_vx.empty() && !traj_vy.empty() && !traj_vz.empty())
      return {expr_vx.value(), expr_vy.value(), expr_vz.value()};

    // Fallback numerical diff (dt small):
    double dt = 1e-3;
    const_cast<dynTraj*>(this)->t_var = t;
    double x0 = expr_x.value(), y0 = expr_y.value(), z0 = expr_z.value();
    const_cast<dynTraj*>(this)->t_var = t + dt;
    double x1 = expr_x.value(), y1 = expr_y.value(), z1 = expr_z.value();
    return {(x1 - x0) / dt, (y1 - y0) / dt, (z1 - z0) / dt};
  }

  /// Evaluate acceleration at time t
  inline Eigen::Vector3d accel(double t) const {
    switch (mode) {
      case Mode::Piecewise:
        return pwp.acceleration(t);
      case Mode::Quintic:
        return accelQuintic(t);
      case Mode::Analytic:
        return accelAnalytic(t);
    }
    return Eigen::Vector3d::Zero();
  }

  inline Eigen::Vector3d accelQuintic(double t) const {
    const double tau = t - poly_start_time;
    return {ddpoly5_abs(cx, tau), ddpoly5_abs(cy, tau), ddpoly5_abs(cz, tau)};
  }

  inline Eigen::Vector3d accelAnalytic(double t) const {
    // If you add analytic second derivatives later, evaluate them here.
    // For now numeric second derivative:
    if (!analytic_compiled) return Eigen::Vector3d::Zero();
    double dt = 1e-3;
    Eigen::Vector3d v1 = velocity(t - dt);
    Eigen::Vector3d v2 = velocity(t + dt);
    return (v2 - v1) / (2 * dt);
  }

  /// Return [t_min, t_max] over which this trajectory is actually defined.
  ///
  /// `Piecewise` reports its breakpoint range — the only mode used by /trajs
  /// from other agents at the moment, since `publishOwnTraj()` sends mode="pwp".
  /// `Quintic` uses the explicit `poly_start_time`/`poly_end_time` fields.
  /// `Analytic` returns an "infinite" range — analytic expressions are valid
  /// for all t (or at least the caller has no better information than that).
  ///
  /// `has_horizon` is false only for the "uninitialized Quintic" case where
  /// poly_start_time == poly_end_time == 0. Callers can use that signal to
  /// fall back to "evaluate everywhere" behavior.
  inline void getHorizon(double& t_min, double& t_max, bool& has_horizon) const {
    switch (mode) {
      case Mode::Piecewise:
        if (!pwp.times.empty()) {
          t_min = pwp.times.front();
          t_max = pwp.times.back();
          has_horizon = (t_max > t_min);
          return;
        }
        break;
      case Mode::Quintic:
        if (poly_end_time > poly_start_time) {
          t_min = poly_start_time;
          t_max = poly_end_time;
          has_horizon = true;
          return;
        }
        break;
      case Mode::Analytic:
        // Analytic expressions are defined everywhere; report no clamp.
        t_min = -std::numeric_limits<double>::infinity();
        t_max = std::numeric_limits<double>::infinity();
        has_horizon = true;
        return;
    }
    t_min = 0.0;
    t_max = 0.0;
    has_horizon = false;
  }

  static const char* modeName(dynTraj::Mode m) {
    switch (m) {
      case dynTraj::Mode::Piecewise:
        return "Piecewise";
      case dynTraj::Mode::Quintic:
        return "Quintic";
      case dynTraj::Mode::Analytic:
        return "Analytic";
      default:
        return "Unknown";
    }
  }

  /// Print debug info
  inline void print() const {
    std::cout << "dynTraj id=" << id << " mode=" << modeName(mode) << "\n";

    if (mode == Mode::Piecewise) {
      pwp.print();
    } else if (mode == Mode::Quintic) {
      std::cout << "  cx=" << cx.transpose() << std::endl;
      std::cout << "  cy=" << cy.transpose() << std::endl;
      std::cout << "  cz=" << cz.transpose() << std::endl;
      std::cout << "  poly_start_time=" << poly_start_time << std::endl;
      std::cout << "  poly_end_time=" << poly_end_time << std::endl;
    } else if (mode == Mode::Analytic) {
      std::cout << "  traj_x='" << traj_x << "'\n";
      std::cout << "  traj_y='" << traj_y << "'\n";
      std::cout << "  traj_z='" << traj_z << "'\n";
      if (!traj_vx.empty() || !traj_vy.empty() || !traj_vz.empty()) {
        std::cout << "  traj_vx='" << traj_vx << "'\n";
        std::cout << "  traj_vy='" << traj_vy << "'\n";
        std::cout << "  traj_vz='" << traj_vz << "'\n";
      }
      std::cout << "  analytic_compiled=" << analytic_compiled << "\n";
    }
  }
};

struct state {
  // time stamp
  double t = 0.0;

  // pos, vel, accel, jerk, yaw, dyaw
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d vel = Eigen::Vector3d::Zero();
  Eigen::Vector3d accel = Eigen::Vector3d::Zero();
  Eigen::Vector3d jerk = Eigen::Vector3d::Zero();
  double yaw = 0.0;
  double dyaw = 0.0;

  // flag for tracking
  bool use_tracking_yaw = false;

  void setTimeStamp(const double data) { t = data; }

  void setPos(const double x, const double y, const double z) { pos << x, y, z; }
  void setVel(const double x, const double y, const double z) { vel << x, y, z; }
  void setAccel(const double x, const double y, const double z) { accel << x, y, z; }

  void setJerk(const double x, const double y, const double z) { jerk << x, y, z; }

  void setPos(const Eigen::Vector3d& data) { pos << data.x(), data.y(), data.z(); }

  void setVel(const Eigen::Vector3d& data) { vel << data.x(), data.y(), data.z(); }

  void setAccel(const Eigen::Vector3d& data) { accel << data.x(), data.y(), data.z(); }

  void setJerk(const Eigen::Vector3d& data) { jerk << data.x(), data.y(), data.z(); }

  void setState(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel,
                const Eigen::Vector3d& accel, const Eigen::Vector3d& jerk) {
    this->pos = pos;
    this->vel = vel;
    this->accel = accel;
    this->jerk = jerk;
  }

  void setYaw(const double data) { yaw = data; }

  void setDYaw(const double data) { dyaw = data; }

  void setZero() {
    pos = Eigen::Vector3d::Zero();
    vel = Eigen::Vector3d::Zero();
    accel = Eigen::Vector3d::Zero();
    jerk = Eigen::Vector3d::Zero();
    yaw = 0;
    dyaw = 0;
  }

  void printPos() { std::cout << "Pos= " << pos.transpose() << std::endl; }

  void print() {
    std::cout << "Time= " << t << std::endl;
    std::cout << "Pos= " << pos.transpose() << std::endl;
    std::cout << "Vel= " << vel.transpose() << std::endl;
    std::cout << "Accel= " << accel.transpose() << std::endl;
  }

  void printHorizontal() {
    std::cout << "Pos, Vel, Accel, Jerk= " << pos.transpose() << " " << vel.transpose() << " "
              << accel.transpose() << " " << jerk.transpose() << std::endl;
  }
};