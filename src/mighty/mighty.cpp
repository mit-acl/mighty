/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "mighty/mighty.hpp"

using namespace mighty;
using namespace termcolor;

typedef timer::Timer MyTimer;

// ----------------------------------------------------------------------------

/**
 * @brief Constructor for MIGHTY.
 * @param parameters par: Input configuration parameters.
 */
MIGHTY::MIGHTY(parameters par) : par_(par) {
  // Set up hgp_manager
  hgp_manager_.setParameters(par_);

  // Set up the planner parameters
  planner_params_.verbose = false;     // enable verbose output
  planner_params_.V_max = par_.v_max;  // max velocity
  planner_params_.A_max = par_.a_max;  // max acceleration
  planner_params_.J_max = par_.j_max;  // max jerk
  planner_params_.num_perturbation =
      par_.num_perturbation_for_ig;                // number of perturbations for initial guesses
  planner_params_.r_max = par_.r_max_for_ig;       // perturbation radius for initial guesses
  planner_params_.time_weight = par_.time_weight;  // weight for time cost
  planner_params_.pos_anchor_weight = par_.pos_anchor_weight;
  planner_params_.dyn_weight = par_.dynamic_weight;
  planner_params_.stat_weight = par_.stat_weight;
  planner_params_.jerk_weight = par_.jerk_weight;
  planner_params_.dyn_constr_vel_weight = par_.dyn_constr_vel_weight;
  planner_params_.dyn_constr_acc_weight = par_.dyn_constr_acc_weight;
  planner_params_.dyn_constr_jerk_weight = par_.dyn_constr_jerk_weight;
  planner_params_.dyn_constr_bodyrate_weight = par_.dyn_constr_bodyrate_weight;
  planner_params_.dyn_constr_tilt_weight = par_.dyn_constr_tilt_weight;
  planner_params_.dyn_constr_thrust_weight = par_.dyn_constr_thrust_weight;
  planner_params_.num_dyn_obst_samples =
      par_.num_dyn_obst_samples;         // Number of dynamic obstacle samples
  planner_params_.Co = par_.planner_Co;  // for static obstacle avoidance
  planner_params_.Cw = par_.planner_Cw;  // for dynamic obstacle avoidance
  planner_params_.BIG = 1e8;
  planner_params_.dc = par_.dc;  // discretization constant
  planner_params_.init_turn_bf = par_.init_turn_bf;
  planner_params_.is_2d_mode = par_.use_2d_planning && (par_.vehicle_type == "ground_robot");
  planner_params_.use_esdf_cost = par_.use_esdf_cost && (par_.vehicle_type == "ground_robot");
  planner_params_.esdf_weight = par_.esdf_weight;
  planner_params_.esdf_d_safe = par_.esdf_d_safe;
  planner_params_.spline_degree = par_.spline_degree;

  // Formation flight: copy weights and configured neighbor (id, offset) pairs
  planner_params_.formation_weight = par_.use_formation ? par_.formation_weight : 0.0;
  planner_params_.formation_neighbor_ids.clear();
  planner_params_.formation_offsets.clear();
  if (par_.use_formation) {
    for (size_t k = 0; k < par_.formation_neighbor_ids.size(); ++k) {
      planner_params_.formation_neighbor_ids.push_back(
          static_cast<int>(par_.formation_neighbor_ids[k]));
      planner_params_.formation_offsets.emplace_back(
          par_.formation_neighbor_offsets[3 * k + 0],
          par_.formation_neighbor_offsets[3 * k + 1],
          par_.formation_neighbor_offsets[3 * k + 2]);
    }
  }

  // Set up the L-BFGS parameters
  lbfgs_params_.mem_size = 256;
  lbfgs_params_.min_step = 1.0e-32;
  lbfgs_params_.f_dec_coeff = par_.f_dec_coeff;          // allow larger Armijo steps
  lbfgs_params_.cautious_factor = par_.cautious_factor;  // always accept BFGS update
  lbfgs_params_.past = par_.past;
  lbfgs_params_.max_linesearch = par_.max_linesearch;  // fewer backtracking tries
  lbfgs_params_.max_iterations = par_.max_iterations;  // allow more iterations
  lbfgs_params_.g_epsilon = par_.g_epsilon;
  lbfgs_params_.delta = par_.delta;  // stop once f-improvement is minimal

  // Set up unconstrained optimization solver for whole trajectory
  whole_traj_solver_ptr_ = std::make_shared<lbfgs::SolverLBFGS>();
  whole_traj_solver_ptr_->initializeSolver(planner_params_);

  // Set up basis converter
  BasisConverter basis_converter;
  A_rest_pos_basis_ = basis_converter.getArestMinvo();  // Use Minvo basis
  A_rest_pos_basis_inverse_ = A_rest_pos_basis_.inverse();

  // Parameters
  v_max_3d_ = Eigen::Vector3d(par_.v_max, par_.v_max, par_.v_max);
  v_max_ = par_.v_max;
  a_max_3d_ = Eigen::Vector3d(par_.a_max, par_.a_max, par_.a_max);
  j_max_3d_ = Eigen::Vector3d(par_.j_max, par_.j_max, par_.j_max);

  // Initialize the state
  changeDroneStatus(DroneStatus::GOAL_REACHED);

  // Initialize the map size
  wdx_ = par_.initial_wdx;
  wdy_ = par_.initial_wdy;
  wdz_ = par_.initial_wdz;

  // Map resolution
  map_res_ = par_.res;
}

// ----------------------------------------------------------------------------

/**
 * @brief Starts adaptive k-value.
 */
void MIGHTY::startAdaptKValue() {
  // Compute the average computation time
  for (int i = 0; i < store_computation_times_.size(); i++) {
    est_comp_time_ += store_computation_times_[i];
  }
  est_comp_time_ = est_comp_time_ / store_computation_times_.size();

  // Start k_value adaptation
  use_adapt_k_value_ = true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Computes the subgoal.
 * @param const state &A: starting state.
 * @param const state &G_term: goal state.
 * @return bool
 */
void MIGHTY::computeG(const state& A, const state& G_term, double horizon) {
  // Initialize the result
  state local_G;

  // Compute pos for G
  local_G.pos = mighty_utils::projectPointToSphere(A.pos, G_term.pos, horizon);

  // Compute yaw for G
  Eigen::Vector3d dir = (G_term.pos - local_G.pos).normalized();
  local_G.yaw = atan2(dir[1], dir[0]);

  // Set G
  setG(local_G);
}

// ----------------------------------------------------------------------------

void MIGHTY::computeG_corridorHop(const state& A, const state& G_term,
                                   vec_Vecf<3>& global_path,
                                   const vec_Vecf<3>& raw_global_path) {
  // Ground-robot-only "bend pre-alignment" subgoal selector.
  //
  // Detects sharp bends on the raw A* path using a windowed direction
  // average. At the first qualifying bend, back off the subgoal from the
  // corner cell by corridor_backoff_m along the windowed incoming direction
  // and set G.yaw to the windowed outgoing direction. Replace global_path
  // with raw_global_path[0..corner_idx] with the last point snapped to the
  // backed-off subgoal so L-BFGS plans only up to the subgoal.
  //
  // Windowed detection is key: the resampled `global_path` can blend a
  // single 90° cell across two waypoints, giving detected angles as low as
  // ~45° for a real right-angle corner. Running the detector on the raw
  // path with a ~0.8 m window averages out A* diagonal grid noise while
  // preserving the true bend angle.
  //
  // Hysteresis: once a subgoal is selected for a leg, cache it in
  // held_subgoal_* and keep reusing it across replans. Cleared on arrival
  // or on YAWING -> TRAVELING transition.
  if (global_path.size() < 2) return;

  // ----- Held-subgoal short circuit (hysteresis) -----
  if (held_subgoal_valid_) {
    const double dx = A.pos.x() - held_subgoal_pos_.x();
    const double dy = A.pos.y() - held_subgoal_pos_.y();
    const double dist_to_held = std::sqrt(dx * dx + dy * dy);
    if (dist_to_held > par_.goal_radius) {
      // Still en route to the held subgoal. Hold G; truncate global_path
      // at the waypoint nearest the held subgoal so L-BFGS stops there.
      state local_G;
      local_G.pos = held_subgoal_pos_;
      local_G.yaw = held_subgoal_yaw_;
      setG(local_G);

      size_t nearest_idx = global_path.size() - 1;
      double best_d2 = std::numeric_limits<double>::max();
      for (size_t i = 0; i < global_path.size(); ++i) {
        const double ddx = global_path[i].x() - held_subgoal_pos_.x();
        const double ddy = global_path[i].y() - held_subgoal_pos_.y();
        const double d2 = ddx * ddx + ddy * ddy;
        if (d2 < best_d2) {
          best_d2 = d2;
          nearest_idx = i;
        }
      }
      if (std::sqrt(best_d2) <= 2.0 * par_.corridor_backoff_m) {
        global_path[nearest_idx] = held_subgoal_pos_;
        global_path.resize(nearest_idx + 1);
        return;
      }
      held_subgoal_valid_ = false;  // stale, fall through to re-pick
    } else {
      // Arrived — let needReplan's TIP trigger fire. Hold G until YAWING
      // exits and clears this flag.
      state local_G;
      local_G.pos = held_subgoal_pos_;
      local_G.yaw = held_subgoal_yaw_;
      setG(local_G);
      return;
    }
  }

  // ----- Windowed bend detection on the raw A* path -----
  // Need at least a few cells to compute directions.
  const vec_Vecf<3>& raw = raw_global_path.size() >= 3 ? raw_global_path : global_path;
  if (raw.size() < 3) {
    // Degenerate — fall through to end-of-path.
    state local_G;
    local_G.pos = global_path.back();
    Eigen::Vector3d dir = G_term.pos - local_G.pos;
    local_G.yaw = (dir.norm() > 1e-6) ? std::atan2(dir.y(), dir.x()) : 0.0;
    held_subgoal_valid_ = false;
    setG(local_G);
    return;
  }

  const double cos_thresh =
      std::cos(par_.corridor_corner_angle_deg * M_PI / 180.0);
  const double window = par_.corridor_detection_window_m;
  const double min_leg_sq = par_.corridor_min_leg_m * par_.corridor_min_leg_m;
  const double backoff = par_.corridor_backoff_m;

  // Helper: compute the arclength distance between two indices along raw path.
  auto arclen_xy = [&](size_t lo, size_t hi) -> double {
    double d = 0.0;
    for (size_t k = lo; k + 1 <= hi; ++k) {
      const double ddx = raw[k + 1].x() - raw[k].x();
      const double ddy = raw[k + 1].y() - raw[k].y();
      d += std::sqrt(ddx * ddx + ddy * ddy);
    }
    return d;
  };

  size_t pick_idx = 0;
  bool found = false;
  Eigen::Vector2d incoming_unit(1.0, 0.0);
  Eigen::Vector2d outgoing_unit(1.0, 0.0);

  for (size_t i = 1; i + 1 < raw.size(); ++i) {
    // Walk back from i until cumulative arclength >= window.
    size_t i_back = i;
    double d_back = 0.0;
    while (i_back > 0 && d_back < window) {
      const double ddx = raw[i_back].x() - raw[i_back - 1].x();
      const double ddy = raw[i_back].y() - raw[i_back - 1].y();
      d_back += std::sqrt(ddx * ddx + ddy * ddy);
      --i_back;
    }

    // Walk forward from i until cumulative arclength >= window.
    size_t i_fwd = i;
    double d_fwd = 0.0;
    while (i_fwd + 1 < raw.size() && d_fwd < window) {
      const double ddx = raw[i_fwd + 1].x() - raw[i_fwd].x();
      const double ddy = raw[i_fwd + 1].y() - raw[i_fwd].y();
      d_fwd += std::sqrt(ddx * ddx + ddy * ddy);
      ++i_fwd;
    }

    // Need at least some arclength on each side to get a meaningful direction.
    if (d_back < 0.5 * window || d_fwd < 0.5 * window) continue;

    Eigen::Vector2d v1(raw[i].x() - raw[i_back].x(),
                       raw[i].y() - raw[i_back].y());
    Eigen::Vector2d v2(raw[i_fwd].x() - raw[i].x(),
                       raw[i_fwd].y() - raw[i].y());
    const double n1 = v1.norm();
    const double n2 = v2.norm();
    if (n1 < 1e-6 || n2 < 1e-6) continue;
    v1 /= n1;
    v2 /= n2;
    const double cos_turn = v1.dot(v2);
    if (cos_turn >= cos_thresh) continue;  // not sharp enough

    // Must be far enough ahead for a useful leg AFTER the backoff.
    const double dx_to_curr = raw[i].x() - A.pos.x();
    const double dy_to_curr = raw[i].y() - A.pos.y();
    if (dx_to_curr * dx_to_curr + dy_to_curr * dy_to_curr <
        min_leg_sq + backoff * backoff) {
      continue;
    }

    pick_idx = i;
    incoming_unit = v1;
    outgoing_unit = v2;
    found = true;
    break;
  }

  state local_G;

  if (found) {
    const Vecf<3>& p_corner = raw[pick_idx];
    local_G.pos = Eigen::Vector3d(p_corner.x() - backoff * incoming_unit.x(),
                                  p_corner.y() - backoff * incoming_unit.y(),
                                  p_corner.z());
    local_G.yaw = std::atan2(outgoing_unit.y(), outgoing_unit.x());

    // Replace global_path (the L-BFGS input) with raw[0..pick_idx] and snap
    // the last point to the backed-off subgoal. This makes the L-BFGS plan
    // stop exactly at the backed-off point.
    global_path.clear();
    global_path.reserve(pick_idx + 1);
    for (size_t k = 0; k <= pick_idx; ++k) {
      global_path.push_back(raw[k]);
    }
    if (!global_path.empty()) {
      global_path.back() = local_G.pos;
    }

    // Cache for hysteresis.
    held_subgoal_pos_ = local_G.pos;
    held_subgoal_yaw_ = local_G.yaw;
    held_subgoal_valid_ = true;
  } else {
    // No sharp bend — use the path end as G with yaw toward the terminal
    // goal (legacy horizon-projection fallback).
    local_G.pos = global_path.back();
    Eigen::Vector3d dir = G_term.pos - local_G.pos;
    local_G.yaw = (dir.norm() > 1e-6) ? std::atan2(dir.y(), dir.x()) : 0.0;
    held_subgoal_valid_ = false;
  }

  setG(local_G);
}

// ----------------------------------------------------------------------------

/**
 * @brief Checks if we need to replan.
 * @return bool
 */

bool MIGHTY::needReplan(const state& local_state, const state& local_G_term,
                        const state& last_plan_state) {
  // Compute the distance to the terminal goal
  double dist_to_term_G = (local_state.pos - local_G_term.pos).norm();
  double dist_from_last_plan_state_to_term_G = (last_plan_state.pos - local_G_term.pos).norm();

  if (dist_to_term_G < par_.goal_radius) {
    changeDroneStatus(DroneStatus::GOAL_REACHED);
    return false;
  }

  if (dist_to_term_G < par_.goal_seen_radius) {
    changeDroneStatus(
        DroneStatus::GOAL_SEEN);  // This triggers to use the hard final state constraint
  }

  if (drone_status_ == DroneStatus::GOAL_SEEN &&
      dist_from_last_plan_state_to_term_G < par_.goal_radius) {
    return false;
  }

  // Corridor-center subgoal arrival -> turn in place (ground robot only).
  // When the robot enters goal_radius of an intermediate corridor subgoal
  // (G_ != G_term), flip to YAWING with the target yaw already stored in
  // G_.yaw (set by computeG_corridorHop). The YAWING branch in the yaw
  // controller reads corridor_hop_yawing_ and uses G_.yaw as the target.
  // No replan happens while YAWING (see the gate below); replanning
  // resumes automatically after the yaw-error threshold at line ~1567.
  if (par_.vehicle_type == "ground_robot" && par_.corridor_hop_enabled &&
      drone_status_ == DroneStatus::TRAVELING) {
    state G_now;
    getG(G_now);
    const double dist_to_G = (local_state.pos - G_now.pos).norm();
    const double dist_G_to_term = (G_now.pos - local_G_term.pos).norm();
    const bool G_is_final = dist_G_to_term < 1e-3;
    if (!G_is_final && dist_to_G < par_.goal_radius) {
      corridor_hop_yawing_ = true;
      changeDroneStatus(DroneStatus::YAWING);
      return false;
    }
  }

  // Don't plan if drone is not traveling
  if (drone_status_ == DroneStatus::GOAL_REACHED || (drone_status_ == DroneStatus::YAWING))
    return false;

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::findAandAtime(state& A, double& A_time, double current_time,
                           double last_replaning_computation_time) {
  int plan_size;
  {
    std::lock_guard<std::mutex> lock(mtx_plan_);
    plan_size = plan_.size();
  }

  if (plan_size == 0) {
    std::cout << bold << red << "plan_size == 0" << reset << std::endl;
    return false;
  }

  if (par_.use_state_update) {
    // Change k_value dynamically
    // To get stable results, we will use a default value of k_value until we have enough
    // computation time
    if (!use_adapt_k_value_) {
      // Use default k_value
      k_value_ = std::max((int)plan_size - par_.default_k_value, 0);

      // Store computation times
      if (num_replanning_ != 1)  // Don't store the very first computation time (because we don't
                                 // have a previous computation time)
        store_computation_times_.push_back(last_replaning_computation_time);
    } else {
      // Computation time filtering
      est_comp_time_ = par_.alpha_k_value_filtering * last_replaning_computation_time +
                       (1 - par_.alpha_k_value_filtering) * est_comp_time_;

      // Get state number based on est_comp_time_ and dc
      k_value_ =
          std::max((int)plan_size - (int)(par_.k_value_factor * est_comp_time_ / par_.dc), 0);
    }

    // Check if k_value_ is valid
    if (plan_size - 1 - k_value_ < 0 || plan_size - 1 - k_value_ >= plan_size) {
      k_value_ =
          plan_size - 1;  // If k_value_ is larger than the plan size, we set it to the last state
    }

    // Get A
    {
      std::lock_guard<std::mutex> lock(mtx_plan_);
      A = plan_[plan_size - 1 - k_value_];
    }

    // Get A_time
    A_time = current_time +
             (plan_size - 1 - k_value_) *
                 par_.dc;  // time to A from current_pos is (plan_size - 1 - k_value_) * par_.dc;
  } else  // If we don't update state - this is for global planner benchmarking purposes
  {
    // Get state
    getState(A);
    A_time = current_time;
  }

  // Check if A is within the map (especially for z)
  if (A.pos[2] < par_.z_min || A.pos[2] > par_.z_max ||
      A.pos[0] < par_.x_min || A.pos[0] > par_.x_max ||
      A.pos[1] < par_.y_min || A.pos[1] > par_.y_max) {
    printf("A (%f, %f, %f) is out of the map\n", A.pos[0], A.pos[1], A.pos[2]);
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::checkIfPointOccupied(const Vec3f& point) {
  // Check if the point is free
  return hgp_manager_.checkIfPointOccupied(point);
}

// ----------------------------------------------------------------------------

bool MIGHTY::checkIfPointFree(const Vec3f& point) {
  // Check if the point is free
  return hgp_manager_.checkIfPointFree(point);
}

// ----------------------------------------------------------------------------

bool MIGHTY::getSafeCorridor(const vec_Vecf<3>& global_path, const state& A) {
  // Skip decomposition entirely when static obstacle weight is zero
  if (par_.stat_weight <= 0.0) {
    safe_corridor_polytopes_whole_.clear();
    poly_out_whole_.clear();
    poly_out_safe_.clear();
    cvx_decomp_time_ = 0.0;
    return true;
  }

  // Timer for computing the safe corridor
  MyTimer cvx_decomp_timer(true);

  // Debug
  if (par_.debug_verbose) std::cout << "Convex decomposition" << std::endl;

  // Set up ellipsoid decomposition utility
  EllipsoidDecomp3D ellip;
  ellip.set_local_bbox(
      Vec3f(par_.local_box_size[0], par_.local_box_size[1], par_.local_box_size[2]));
  ellip.set_z_min_and_max(par_.z_min, par_.z_max);
  ellip.set_inflate_distance(par_.drone_radius);

  // Get obstacle points for decomposition
  vec_Vec3f base_uo;
  if (par_.sfc_use_unknown_as_obstacle)
    hgp_manager_.getVecUnknownOccupied(base_uo);
  else
    hgp_manager_.getVecOccupied(base_uo);

  // Empty obstacle arrays (dynamic obstacles handled via heat map, not corridors)
  vec_Vecf<3> obst_pos, obst_bbox;
  // seg_end_times must have one positive entry per segment (path.size()-1)
  std::vector<double> seg_end_times(global_path.size() > 1 ? global_path.size() - 1 : 0, 1.0);

  // Check if the convex decomposition failed
  if (!hgp_manager_.cvxEllipsoidDecomp(ellip, global_path, base_uo, obst_pos, obst_bbox,
                                       seg_end_times, safe_corridor_polytopes_whole_,
                                       poly_out_whole_)) {
    std::cout << bold << red << "Convex decomposition failed" << reset << std::endl;
    poly_out_whole_.clear();
    poly_out_safe_.clear();
    return false;
  }

  // Generate safe corridor (treat unknown regions as occupied)
  {
    EllipsoidDecomp3D ellip_safe;
    ellip_safe.set_local_bbox(
        Vec3f(par_.local_box_size[0], par_.local_box_size[1], par_.local_box_size[2]));
    ellip_safe.set_z_min_and_max(par_.z_min, par_.z_max);
    ellip_safe.set_inflate_distance(par_.drone_radius);

    vec_Vec3f base_uo_safe;
    hgp_manager_.getVecUnknownOccupied(base_uo_safe);

    std::vector<LinearConstraint3D> safe_constraints;
    if (!hgp_manager_.cvxEllipsoidDecomp(ellip_safe, global_path, base_uo_safe, obst_pos, obst_bbox,
                                         seg_end_times, safe_constraints, poly_out_safe_)) {
      // Fallback: use whole corridor if safe decomposition fails
      poly_out_safe_ = poly_out_whole_;
    }
  }

  // Get computation time [ms]
  cvx_decomp_time_ = cvx_decomp_timer.getElapsedMicros() / 1000.0;

  return true;
}

// ----------------------------------------------------------------------------

void MIGHTY::computeMapSize(const Eigen::Vector3d& min_pos, const Eigen::Vector3d& max_pos) {
  // Get local_A
  state local_A;
  getA(local_A);

  // Increase the effective buffer size based on the number of HGP failures.
  double dynamic_buffer = par_.map_buffer;

  // Increase the effective buffer size based on velocity.
  double dynamic_buffer_x = dynamic_buffer;
  double dynamic_buffer_y = dynamic_buffer;
  double dynamic_buffer_z = dynamic_buffer;

  // Compute the distance to the terminal goal for each axis.
  double dist_x = std::abs(min_pos[0] - max_pos[0]);
  double dist_y = std::abs(min_pos[1] - max_pos[1]);
  double dist_z = std::abs(min_pos[2] - max_pos[2]);

  // Update the map size based on the min and max positions.
  wdx_ = std::max(dist_x + 2 * dynamic_buffer_x, par_.min_wdx);
  wdy_ = std::max(dist_y + 2 * dynamic_buffer_y, par_.min_wdy);
  wdz_ = std::max(dist_z + 2 * dynamic_buffer_z, par_.min_wdz);

  // Use the agent's position as the map center (no shifting toward goal)
  map_center_ = min_pos;
}

// ----------------------------------------------------------------------------

bool MIGHTY::checkPointWithinMap(const Eigen::Vector3d& point) const {
  // Check if the point is within the map boundaries for each axis
  return (std::abs(point[0] - map_center_[0]) <= wdx_ / 2.0) &&
         (std::abs(point[1] - map_center_[1]) <= wdy_ / 2.0) &&
         (std::abs(point[2] - map_center_[2]) <= wdz_ / 2.0);
}

// ----------------------------------------------------------------------------

void MIGHTY::getStaticPushPoints(vec_Vecf<3>& static_push_points) {
  static_push_points = static_push_points_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getLocalGlobalPath(vec_Vecf<3>& local_global_path,
                                vec_Vecf<3>& local_global_path_after_push) {
  local_global_path = local_global_path_;
  local_global_path_after_push = local_global_path_after_push_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getGlobalPath(vec_Vecf<3>& global_path) {
  std::lock_guard<std::mutex> lock(mtx_global_path_);
  global_path = global_path_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getOriginalGlobalPath(vec_Vecf<3>& original_global_path) {
  std::lock_guard<std::mutex> lock(mtx_original_global_path_);
  original_global_path = original_global_path_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getFreeGlobalPath(vec_Vecf<3>& free_global_path) {
  free_global_path = free_global_path_;
}

// ----------------------------------------------------------------------------

void MIGHTY::resetData() {
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

  poly_out_whole_.clear();
  poly_out_safe_.clear();
  goal_setpoints_.clear();
  pwp_to_share_.clear();
  optimal_yaw_sequence_.clear();
  yaw_control_points_.clear();
  yaw_knots_.clear();
  cps_.clear();
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveData(double& final_g, double& global_planning_time,
                          double& hgp_static_jps_time, double& hgp_check_path_time,
                          double& hgp_dynamic_astar_time, double& hgp_recover_path_time,
                          double& cvx_decomp_time, double& initial_guess_computation_time,
                          double& local_traj_computatoin_time, double& safety_check_time,
                          double& safe_paths_time, double& yaw_sequence_time,
                          double& yaw_fitting_time) {
  final_g = final_g_;
  global_planning_time = global_planning_time_;
  hgp_static_jps_time = hgp_static_jps_time_;
  hgp_check_path_time = hgp_check_path_time_;
  hgp_dynamic_astar_time = hgp_dynamic_astar_time_;
  hgp_recover_path_time = hgp_recover_path_time_;
  cvx_decomp_time = cvx_decomp_time_;
  initial_guess_computation_time = initial_guess_computation_time_;
  local_traj_computatoin_time = local_traj_computation_time_;
  safe_paths_time = safe_paths_time_;
  safety_check_time = safety_check_time_;
  yaw_sequence_time = yaw_sequence_time_;
  yaw_fitting_time = yaw_fitting_time_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrievePolytopes(vec_E<Polyhedron<3>>& poly_out_whole,
                               vec_E<Polyhedron<3>>& poly_out_safe) {
  poly_out_whole = poly_out_whole_;
  poly_out_safe = poly_out_safe_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveGoalSetpoints(std::vector<state>& goal_setpoints) {
  goal_setpoints = goal_setpoints_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveListSubOptGoalSetpoints(
    std::vector<std::vector<state>>& list_subopt_goal_setpoints) {
  list_subopt_goal_setpoints = list_subopt_goal_setpoints_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveCPs(std::vector<Eigen::Matrix<double, 3, 6>>& cps) { cps = cps_; }

// ----------------------------------------------------------------------------

/**
 * @brief Replans the trajectory.
 * @param double last_replaning_computation_time: Last replanning computation time.
 * @param double current_time: Current timestamp.
 */
std::tuple<bool, bool> MIGHTY::replan(double last_replaning_computation_time, double current_time) {
  /* -------------------- Housekeeping -------------------- */

  MyTimer timer_housekeeping(true);

  // Reset Data
  resetData();

  // Check if we need to replan
  if (!checkReadyToReplan()) return std::make_tuple(false, false);

  // Get states we need
  state local_state, local_G_term, last_plan_state;
  getState(local_state);
  getGterm(local_G_term);
  getLastPlanState(last_plan_state);

  // Re-run goal sanitization each replan: as the planning window slides
  // forward and previously-unknown cells become observed, an occupied goal
  // may only become detectable now. sanitizeTerminalGoal is a no-op when
  // the goal is already free or unknown, so this is cheap in the common case.
  if (par_.relocate_occupied_goal) {
    state goal_before = local_G_term;
    if (sanitizeTerminalGoal(local_G_term)) {
      if ((local_G_term.pos - goal_before.pos).norm() > 1e-6) {
        setGterm(local_G_term);  // persist the relocation for the next cycle
      }
    }
    // If sanitize returned false, leave the stored G_term_ unchanged and
    // let the planner try; the next sensor update may make a non-occupied
    // cell available.
  }

  // Check if we need to replan based on the distance to the terminal goal
  if (!needReplan(local_state, local_G_term, last_plan_state)) return std::make_tuple(false, false);

  if (par_.debug_verbose)
    printf("[TIMING] Housekeeping: %.2f ms\n", timer_housekeeping.getElapsedMicros() / 1000.0);

  /* -------------------- Global Planning -------------------- */

  MyTimer timer_global(true);
  vec_Vecf<3> global_path;
  if (!generateGlobalPath(global_path, current_time, last_replaning_computation_time)) {
    if (par_.debug_verbose)
      printf("[TIMING] Global Planning (FAILED): %.2f ms\n",
             timer_global.getElapsedMicros() / 1000.0);
    return std::make_tuple(false, false);
  }
  if (par_.debug_verbose)
    printf("[TIMING] Global Planning: %.2f ms\n", timer_global.getElapsedMicros() / 1000.0);

  /* -------------------- Local Trajectory Optimization -------------------- */

  MyTimer timer_local(true);
  if (!planLocalTrajectory(global_path)) {
    if (par_.debug_verbose)
      printf("[TIMING] Local Traj Opt (FAILED): %.2f ms\n",
             timer_local.getElapsedMicros() / 1000.0);
    return std::make_tuple(false, true);
  }

  /* -------------------- Append to Plan -------------------- */

  // For ground robots, we don't need to maintain plan_ since we use pure pursuit lookahead directly
  if (par_.vehicle_type != "uav") {
    if (par_.debug_verbose)
      std::cout << "Append to Plan: 0 ms (skipped for ground robot)" << std::endl;
  } else {
    MyTimer timer_append(true);
    if (!appendToPlan()) {
      if (par_.debug_verbose)
        std::cout << "Append to Plan: " << timer_append.getElapsedMicros() / 1000.0 << " ms"
                  << std::endl;
      return std::make_tuple(false, true);
    }
    if (par_.debug_verbose)
      std::cout << "Append to Plan: " << timer_append.getElapsedMicros() / 1000.0 << " ms"
                << std::endl;
  }

  /* -------------------- Final Housekeeping -------------------- */

  MyTimer timer_final(true);

  if (par_.debug_verbose)
    std::cout << bold << green << "Replanning succeeded" << reset << std::endl;

  // Reset the replanning failure count
  replanning_failure_count_ = 0;
  if (par_.debug_verbose)
    std::cout << "Final Housekeeping: " << timer_final.getElapsedMicros() / 1000.0 << " ms"
              << std::endl;

  return std::make_tuple(true, true);
}

// ----------------------------------------------------------------------------

bool MIGHTY::generateGlobalPath(vec_Vecf<3>& global_path, double current_time,
                                double last_replaning_computation_time) {
  // Get G and G_term
  state local_G, local_G_term;
  getG(local_G);
  getGterm(local_G_term);

  // Declare local variables
  state local_A;
  double A_time;

  // For ground robots, use pure pursuit lookahead point for velocity direction
  // For UAVs, compute point A from trajectory to account for momentum
  if (par_.vehicle_type != "uav") {
    // Ground robot: use current position but zero velocity
    // This lets the optimizer plan freely without velocity direction constraint
    state current_state;
    getState(current_state);

    local_A = current_state;  // Use current velocity as initial condition
    A_time = current_time;
  } else {
    // UAV: find A and A_time from trajectory
    if (!findAandAtime(local_A, A_time, current_time, last_replaning_computation_time)) {
      replanning_failure_count_++;
      return false;
    }
  }

  // Set A and A_time
  setA(local_A);
  setA_time(A_time);

  // Compute G (writes the freshly-projected subgoal into the member G_)
  computeG(local_A, local_G_term, par_.horizon);

  // Refresh local_G with the value computeG just wrote. The earlier getG() at
  // the top of this function returned the previous replan's subgoal, which
  // becomes stale the moment computeG updates G_. Without this re-read,
  // solveHGP below would plan toward the *previous* subgoal — and right after
  // a new term_goal is set (or a large G_term jump), the stale subgoal can
  // point in the opposite direction of the new one, producing a brief
  // backward HGP path.
  getG(local_G);

  // Set up the HGP planner (since updateVmax() needs to be called after setupHGPPlanner, we use
  // v_max_ from the last replan)
  // Pass ESDF and occ grids to HGP manager for ground robot A* planning
  if (par_.use_esdf_cost && par_.vehicle_type == "ground_robot") {
    if (occ_grid_2d_) hgp_manager_.setOccGrid2D(occ_grid_2d_);
    if (esdf_grid_) hgp_manager_.setEsdfGrid(esdf_grid_, par_.esdf_weight, par_.esdf_d_safe);
  }

  MyTimer timer_setup(true);
  hgp_manager_.setupHGPPlanner(
      par_.global_planner, par_.global_planner_verbose, map_res_, v_max_, par_.a_max, par_.j_max,
      par_.hgp_timeout_duration_ms, par_.max_num_expansion, par_.w_unknown, par_.w_align,
      par_.decay_len_cells, par_.w_side, par_.los_cells, par_.min_len, par_.min_turn);

  // For ground robot 2D mode: set max distance between path waypoints
  if (par_.use_2d_planning && par_.vehicle_type == "ground_robot") {
    hgp_manager_.setMaxDistVertexes2D(par_.max_dist_vertexes);
  }

  // Free start and goal if necessary
  if (par_.use_free_start) hgp_manager_.freeStart(local_A.pos, par_.free_start_factor);
  if (par_.use_free_goal) hgp_manager_.freeGoal(local_G.pos, par_.free_goal_factor);
  if (par_.debug_verbose)
    printf("[TIMING]   HGP setup: %.2f ms\n", timer_setup.getElapsedMicros() / 1000.0);

  // if using ground robot, fix z to default_goal_z
  if (par_.vehicle_type != "uav") {
    local_A.pos[2] = par_.default_goal_z;
    local_G.pos[2] = par_.default_goal_z;
  }

  // 1) Build a direction hint from the *previous* global path
  vec_Vecf<3> prev_global;
  getGlobalPath(prev_global);  // last successful global path

  Eigen::Vector3d dir_hint = (local_G.pos - local_A.pos).normalized();
  if (prev_global.size() >= 2) {
    Eigen::Vector3d s0 = prev_global[0];
    Eigen::Vector3d s1 = prev_global[1];
    Eigen::Vector3d seg = s1 - s0;
    if (seg.norm() > 1e-8) {
      dir_hint = seg.normalized();
    }
  } else {
    dir_hint = local_G.pos - local_A.pos;
    if (dir_hint.norm() > 1e-8) dir_hint.normalize();
  }

  // Keep ground robots planar
  if (par_.vehicle_type != "uav") dir_hint[2] = 0.0;

  Vec3f start_dir_hint(dir_hint.x(), dir_hint.y(), dir_hint.z());

  // Solve HGP
  MyTimer timer_solve(true);
  vec_Vecf<3> raw_global_path;
  if (!hgp_manager_.solveHGP(local_A.pos, start_dir_hint, local_G.pos, final_g_,
                             par_.global_planner_heuristic_weight, A_time, global_path,
                             raw_global_path)) {
    if (par_.debug_verbose)
      printf("[TIMING]   HGP solve (FAILED): %.2f ms\n", timer_solve.getElapsedMicros() / 1000.0);
    hgp_failure_count_++;
    replanning_failure_count_++;
    return false;
  }
  if (par_.debug_verbose)
    printf("[TIMING]   HGP solve: %.2f ms\n", timer_solve.getElapsedMicros() / 1000.0);

  // Bend pre-alignment (ground robot only). Detects sharp bends on the raw
  // A* path using a windowed direction average (robust to resample phase),
  // then replaces global_path with a truncated version ending at a backed-
  // off subgoal so L-BFGS stops short of each inside corner. UAVs unaffected.
  if (par_.vehicle_type == "ground_robot" && par_.corridor_hop_enabled) {
    computeG_corridorHop(local_A, local_G_term, global_path, raw_global_path);
    // Refresh local_G with the (possibly) new subgoal; local_G is read again
    // later in this function and fed into the safe-corridor / local-traj
    // pipeline in planLocalTrajectory.
    getG(local_G);
  }

  // use this for map resizing
  {
    std::lock_guard<std::mutex> lock(mtx_global_path_);
    global_path_ = global_path;
  }

  // For visualization — store the raw A* path (before smoothing/resampling)
  {
    std::lock_guard<std::mutex> lock(mtx_original_global_path_);
    original_global_path_ = raw_global_path;
    // Fix z for ground robots
    if (par_.vehicle_type != "uav") {
      for (auto& wp : original_global_path_) wp[2] = par_.default_goal_z;
    }
  }

  // Debug
  if (par_.debug_verbose) std::cout << "global_path.size(): " << global_path.size() << std::endl;

  // Get computation time
  hgp_manager_.getComputationTime(global_planning_time_, hgp_static_jps_time_, hgp_check_path_time_,
                                  hgp_dynamic_astar_time_, hgp_recover_path_time_);

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::planLocalTrajectory(vec_Vecf<3>& global_path) {
  // Get local_A, local_G and A_time
  state local_A, local_G;
  double A_time;
  getA(local_A);
  getG(local_G);
  getA_time(A_time);

  // If the global path's size is < 2 after trimming, we cannot proceed
  if (global_path.empty() || global_path.size() < 2) {
    // std::cout << bold << red << "Global path's size is < 2 after trimming" << reset << std::endl;
    replanning_failure_count_++;
    return false;
  }

  // Truncate global path to num_N + 1 segments (num_N pieces need num_N+1 waypoints)
  // This limits the local optimization horizon similar to dynus
  if (par_.num_N > 0 && global_path.size() > static_cast<size_t>(par_.num_N + 1)) {
    global_path.resize(par_.num_N + 1);
    if (par_.debug_verbose)
      std::cout << "[Local Traj] Truncated global path to " << global_path.size()
                << " points (num_N=" << par_.num_N << ")" << std::endl;
  }

  // If the global path has exactly 2 points, add a midpoint to make it 3 points
  if (global_path.size() == 2) {
    Eigen::Vector3d midpoint = (global_path[0] + global_path[1]) / 2.0;
    global_path.insert(global_path.begin() + 1, midpoint);
  }

  // convex decomposition
  MyTimer timer_sfc(true);
  if (!getSafeCorridor(global_path, local_A)) {
    if (par_.debug_verbose)
      printf("[TIMING]   Safe corridor (FAILED): %.2f ms\n", timer_sfc.getElapsedMicros() / 1000.0);
    replanning_failure_count_++;
    return false;
  }
  if (par_.debug_verbose)
    printf("[TIMING]   Safe corridor: %.2f ms\n", timer_sfc.getElapsedMicros() / 1000.0);

  // Initialize flag
  bool optimization_succeeded = false;

  if (par_.vehicle_type != "uav") {
    local_A.pos[2] = par_.default_goal_z;
  }

  MyTimer timer_opt(true);
  optimization_succeeded =
      generateLocalTrajectory(local_A, A_time, global_path, initial_guess_computation_time_,
                              local_traj_computation_time_, whole_traj_solver_ptr_);
  if (par_.debug_verbose)
    printf("[TIMING]   L-BFGS optimization: %.2f ms\n", timer_opt.getElapsedMicros() / 1000.0);

  if (par_.debug_verbose) {
    std::cout << "initial_guess_computation_time_: " << initial_guess_computation_time_ << " ms"
              << std::endl;
    std::cout << "Local trajectory optimization finished" << std::endl;
  }

  if (optimization_succeeded) {
    // For the optimal solution
    whole_traj_solver_ptr_->reconstructPVATCPopt(
        zopt_);  // First recover the final control points and times from z_opt
    whole_traj_solver_ptr_->getGoalSetpoints(goal_setpoints_);

    // print out the goal setpoints
    whole_traj_solver_ptr_->getPieceWisePol(pwp_to_share_);
    whole_traj_solver_ptr_->getControlPoints(cps_);  // Bezier control points

    // Get goal setpoints for suboptimal solutions for visualization
    if (par_.use_multiple_initial_guesses) {
      if (par_.debug_verbose)
        std::cout << "Size of list_z_subopt_: " << list_z_subopt_.size() << std::endl;

      // Initialize list_subopt_goal_setpoints_
      list_subopt_goal_setpoints_.clear();

      // Loop over list_z_subopt
      for (int idx = 0; idx < list_z_subopt_.size(); ++idx) {
        // Reconstruct control points and times for each suboptimal solution
        whole_traj_solver_ptr_->reconstructPVATCPopt(
            list_z_subopt_[idx]);  // First recover the final control points and times from z_opt
        std::vector<state> subopt_goal_setpoints;
        whole_traj_solver_ptr_->getGoalSetpoints(subopt_goal_setpoints);
        list_subopt_goal_setpoints_.push_back(subopt_goal_setpoints);
      }
    }
  } else {
    replanning_failure_count_++;
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------

void MIGHTY::getPiecewiseQuinticPol(PieceWiseQuinticPol& pwp) { pwp = pwp_to_share_; }

// ----------------------------------------------------------------------------

bool MIGHTY::generateLocalTrajectory(const state& local_A, double A_time, vec_Vec3f& global_path,
                                     double& initial_guess_computation_time,
                                     double& local_traj_computation_time,
                                     std::shared_ptr<lbfgs::SolverLBFGS>& whole_traj_solver_ptr) {
  if (par_.debug_verbose) std::cout << "Preparing solver for replan" << std::endl;

  auto local_trajs = getTrajs();  // RVO/NRVO elides copy

  // Get local_G
  state local_G;
  getG(local_G);

  state local_E;
  Vec3f mean_point;

  if (drone_status_ == DroneStatus::GOAL_REACHED || drone_status_ == DroneStatus::GOAL_SEEN) {
    local_E = local_G;
  } else {
    local_E.pos = global_path.back();
  }

  // if using ground robot, we fix the z
  if (par_.vehicle_type != "uav") {
    local_E.pos[2] = par_.default_goal_z;
  }

  // Pass ESDF grid to solver (ground robot only, thread-safe immutable snapshot)
  if (esdf_grid_) whole_traj_solver_ptr->setEsdfGrid(esdf_grid_);

  whole_traj_solver_ptr->prepareSolverForReplan(
      A_time, global_path, safe_corridor_polytopes_whole_, local_trajs, local_A, local_E,
      initial_guess_computation_time, par_.use_multiple_initial_guesses);

  // Skip the solver's pushed path — use the original HGP path instead
  // whole_traj_solver_ptr->getGlobalPath(global_path);

  {
    std::lock_guard<std::mutex> lock(mtx_global_path_);
    global_path_ = global_path;  // Update the global path
  }

  // update local_E
  local_E.pos = global_path.back();
  if (par_.vehicle_type != "uav") local_E.pos[2] = par_.default_goal_z;
  {
    std::lock_guard<std::mutex> lock(mtx_E_);
    E_ = local_E;  // Update the local_E
  }

  if (par_.debug_verbose) std::cout << "Solver prepared" << std::endl;

  // Get initial guesses
  auto list_z0 = whole_traj_solver_ptr_->getInitialGuesses();
  auto list_initial_guess_wps = whole_traj_solver_ptr_->getInitialGuessWaypoints();

  if (par_.debug_verbose) std::cout << "Initial guesses size: " << list_z0.size() << std::endl;

  // Update L-BFGS parameters.
  // lbfgs_params_.mem_size = static_cast<int>(list_z0[0].size());
  lbfgs_params_.mem_size = 256;

  // Prepare vectors for parallelization
  int status = -1;  // Initialize status
  int size_of_list_z0 = list_z0.size();

  std::vector<int> list_status(size_of_list_z0, -1);  // Initialize status for each thread
  std::vector<Eigen::VectorXd> list_zopt(size_of_list_z0);
  std::vector<double> list_fopt(size_of_list_z0, 0.0);
  std::vector<double> list_initial_guess_computation_time(size_of_list_z0);

  // Solve the optimiation problem.
  if (!par_.use_multiple_initial_guesses || size_of_list_z0 == 1) {
    auto t_start = std::chrono::high_resolution_clock::now();
    status = whole_traj_solver_ptr_->optimize(list_z0[0], zopt_, fopt_, lbfgs_params_);
    auto t_end = std::chrono::high_resolution_clock::now();
    local_traj_computation_time =
        std::chrono::duration<double, std::milli>(t_end - t_start).count();
    // std::cout << lbfgs::lbfgs_strerror(status) << std::endl;

    // status = 10;
    // zopt_ = list_z0[0]; // Initialize zopt_ with the correct size
    // fopt_ = 0.0; // Initialize fopt_ to zero
    list_initial_guess_wps_subopt_.clear();
    list_initial_guess_wps_subopt_.push_back(
        list_initial_guess_wps[0]);  // Store the initial guess waypoints for the optimal solution
  } else {
    // Parallelization approach (is there any faster way to do this?)
    std::vector<std::future<void>> futures;
    futures.reserve(size_of_list_z0);

    auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < size_of_list_z0; ++i) {
      futures.emplace_back(std::async(std::launch::async, [&, i]() {
        // make a fresh solver for thread-safety
        std::shared_ptr<lbfgs::SolverLBFGS> solver_ptr = std::make_shared<lbfgs::SolverLBFGS>();
        solver_ptr->initializeSolver(planner_params_);
        if (esdf_grid_) solver_ptr->setEsdfGrid(esdf_grid_);
        double initial_guess_computation_time = 0.0;
        solver_ptr->prepareSolverForReplan(
            A_time, global_path, safe_corridor_polytopes_whole_, local_trajs, local_A, local_E,
            initial_guess_computation_time);  // initial time t0 = 0.0

        // copy initial guess
        Eigen::VectorXd z0 = list_z0[i];
        Eigen::VectorXd zopt;
        double fopt;

        // run optimization
        int status = solver_ptr->optimize(z0, zopt, fopt, lbfgs_params_);
        list_status[i] = status;
        list_zopt[i] = zopt;
        list_fopt[i] = fopt;
        list_initial_guess_computation_time[i] = initial_guess_computation_time;
      }));
    }

    // wait for all to finish
    for (auto& f : futures) f.get();

    auto t_end = std::chrono::high_resolution_clock::now();
    local_traj_computation_time =
        std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // Find the best solution
    fopt_ = std::numeric_limits<double>::max();
    zopt_ = list_z0[0];  // Initialize zopt_ with the correct size
    int best_index = -1;
    for (size_t i = 0; i < size_of_list_z0; ++i) {
      // First check the status
      if (list_status[i] < 0) {
        std::cout << "list_status[" << i << "] = " << list_status[i] << ", skipping this solution"
                  << std::endl;
        std::cout << lbfgs::lbfgs_strerror(list_status[i]) << std::endl;
        continue;  // Skip this solution if it failed
      }

      if (list_fopt[i] < fopt_) {
        fopt_ = list_fopt[i];
        zopt_ = list_zopt[i];
        best_index = i;
      }
    }

    initial_guess_computation_time =
        list_initial_guess_computation_time[best_index];  // Get the initial guess computation time
                                                          // for the best solution

    // If no solution was found, return
    if (best_index < 0) {
      std::cout << bold << red << "No solution found in multiple initial guesses" << reset
                << std::endl;
      replanning_failure_count_++;
      initial_guess_computation_time = -100000.0;
      return false;
    }

    list_z_subopt_.clear();                  // Clear the suboptimal solutions list
    list_initial_guess_wps_subopt_.clear();  // Clear the suboptimal waypoints list

    // Grab sub optimal solutions for visualization (solutions that are valid (so no negative status
    // and no crazy big f) but no better than the best one)
    for (size_t i = 0; i < size_of_list_z0; ++i) {
      if (list_status[i] >= 0 && i != best_index && list_fopt[i] < planner_params_.BIG * 0.9) {
        // Add the solution to the trajs_ for visualization
        list_z_subopt_.push_back(list_zopt[i]);

        // Add the initial guess waypoints for visualization
        list_initial_guess_wps_subopt_.push_back(list_initial_guess_wps[i]);
      }
    }

    // Set the status to the best one found
    status = list_status[best_index];

  }  // End of parallelization

  if (par_.debug_verbose)
    std::cout << "Optimization status: " << status << ", fopt: " << fopt_
              << ", computation time: " << local_traj_computation_time << " ms" << std::endl;

  // If no solution is found, return.
  if (status < 0 || fopt_ > par_.fopt_threshold) {
    // do the same output in red with printf
    printf("\033[1;31mLocal Optimization Failed with status: %d, fopt: %.2f\033[0m\n", status,
           fopt_);
    return false;
  }

  // If the optimization succeeded, we can update the goal setpoints in green
  // printf("\033[1;32mLocal Optimization Succeeded with status: %d, fopt: %.2f\033[0m\n", status,
  // fopt_);
  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::appendToPlan() {
  if (par_.debug_verbose)
    std::cout << "goal_setpoints_.size(): " << goal_setpoints_.size() << std::endl;

  {
    std::lock_guard<std::mutex> lock(mtx_plan_);

    // get the size of the plan and plan_safe_paths
    int plan_size = plan_.size();

    // If the plan size is less than k_value_, which means we already passed point A, we cannot use
    // this plan
    if (plan_size < k_value_) {
      if (par_.debug_verbose)
        std::cout << bold << red << "(plan_size - k_value_) = " << (plan_size - k_value_) << " < 0"
                  << reset << std::endl;
      k_value_ = std::max(1, plan_size - 1);  // Decrease k_value_ to plan_size - 1 but at least 1
    } else  // If the plan size is greater than k_value_, which means we haven't passed point A yet,
            // we can use this plan
    {
      plan_.erase(plan_.end() - k_value_, plan_.end());
      plan_.insert(plan_.end(), goal_setpoints_.begin(), goal_setpoints_.end());
    }
  }

  // k_value adaptation initialization
  if (!got_enough_replanning_) {
    if (store_computation_times_.size() < par_.num_replanning_before_adapt) {
      num_replanning_++;
    } else {
      startAdaptKValue();
      got_enough_replanning_ = true;
    }
  }

  return true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the terminal goal state.
 * @param state &G_term: Output terminal goal state.
 */
void MIGHTY::getGterm(state& G_term) {
  std::lock_guard<std::mutex> lock(mtx_G_term_);
  G_term = G_term_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the terminal goal state.
 * @param state G_term: Terminal goal state to set.
 */
void MIGHTY::setGterm(const state& G_term) {
  std::lock_guard<std::mutex> lock(mtx_G_term_);
  G_term_ = G_term;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the subgoal.
 * @param state &G: Output subgoal.
 */
void MIGHTY::getG(state& G) {
  std::lock_guard<std::mutex> lock(mtx_G_);
  G = G_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets point E
 * @param state &G: Output point E
 */
void MIGHTY::getE(state& E) {
  std::lock_guard<std::mutex> lock(mtx_E_);
  E = E_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the subgoal.
 * @param state G: Subgoal to set.
 */
void MIGHTY::setG(const state& G) {
  std::lock_guard<std::mutex> lock(mtx_G_);
  G_ = G;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets A (starting point for global planning).
 * @param state &G: Output A.
 */
void MIGHTY::getA(state& A) {
  std::lock_guard<std::mutex> lock(mtx_A_);
  A = A_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets A (starting point for global planning).
 * @param state &G: Input A.
 */
void MIGHTY::setA(const state& A) {
  std::lock_guard<std::mutex> lock(mtx_A_);
  A_ = A;
}

// ----------------------------------------------------------------------------

void MIGHTY::getA_time(double& A_time) {
  std::lock_guard<std::mutex> lock(mtx_A_time_);
  A_time = A_time_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets A (starting point for global planning)'s time
 * @param state &G: Input A time
 */
void MIGHTY::setA_time(double A_time) {
  std::lock_guard<std::mutex> lock(mtx_A_time_);
  A_time_ = A_time;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the current state.
 * @param state &state: Output current state.
 */
void MIGHTY::getState(state& state) {
  std::lock_guard<std::mutex> lock(mtx_state_);
  state = state_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the last plan state
 * @param state &state: Output last plan state
 */
void MIGHTY::getLastPlanState(state& state) {
  // For ground robots, we don't maintain plan_, so just return the terminal goal
  if (par_.vehicle_type != "uav") {
    getGterm(state);
  } else {
    std::lock_guard<std::mutex> lock(mtx_plan_);
    state = plan_.back();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets trajs_ (returns copy of shared_ptrs, not deep copy of trajectory data)
 * @return std::vector<std::shared_ptr<dynTraj>>: Copy of trajs_ vector
 */
std::vector<std::shared_ptr<dynTraj>> MIGHTY::getTrajs() {
  std::lock_guard<std::mutex> lock(mtx_trajs_);
  return trajs_;  // RVO/NRVO will elide the copy in most cases
}

// ----------------------------------------------------------------------------

/**
 * @brief Cleans up old trajectories.
 * @param double current_time: Current timestamp.
 */
void MIGHTY::cleanUpOldTrajs(double current_time) {
  std::lock_guard<std::mutex> lock(mtx_trajs_);

  // remove_if moves all “expired” to the end, then erase() chops them off
  trajs_.erase(std::remove_if(trajs_.begin(), trajs_.end(),
                              [&](const std::shared_ptr<dynTraj>& t) {
                                return (current_time - t->time_received) > par_.traj_lifetime;
                              }),
               trajs_.end());
}

// ----------------------------------------------------------------------------

/**
 * @brief Adds or updates a trajectory.
 * @param dynTraj new_traj: New trajectory to add.
 * @param double current_time: Current timestamp.
 */
void MIGHTY::addTraj(std::shared_ptr<dynTraj> new_traj, double current_time) {
  if (new_traj->mode == dynTraj::Mode::Analytic && !new_traj->analytic_compiled) {
    printf("Dropping analytic traj id=%d because compile failed", new_traj->id);
    return;
  }

  if (new_traj->mode == dynTraj::Mode::Piecewise && new_traj->pwp.times.empty()) return;

  Eigen::Vector3d p = new_traj->eval(current_time);

  // If this is a predicted (non-agent) trajectory, mask it if an agent trajectory
  // is nearby — avoids double-counting when ACL mapper tracks another agent
  if (!new_traj->is_agent && par_.prediction_mask_distance > 0.0) {
    std::lock_guard<std::mutex> lock(mtx_trajs_);
    for (const auto& t : trajs_) {
      if (t->is_agent) {
        Eigen::Vector3d agent_pos = t->eval(current_time);
        if ((p - agent_pos).norm() < par_.prediction_mask_distance) {
          return;  // drop — an agent trajectory already covers this obstacle
        }
      }
    }
  }

  {
    std::lock_guard<std::mutex> lock(mtx_trajs_);
    auto it = std::find_if(trajs_.begin(), trajs_.end(), [&](const std::shared_ptr<dynTraj>& t) {
      return t->id == new_traj->id;
    });

    if (it != trajs_.end())
      *it = new_traj;
    else
      trajs_.push_back(new_traj);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Updates the current state.
 * @param state data: New state data.
 */
void MIGHTY::updateState(state data) {
  // If we are doing hardware and provide goal in global frame (e.g. vicon), we need to transform
  // the goal to the local frame

  if (par_.use_hardware && par_.provide_goal_in_global_frame) {
    // Apply transformation to position
    Eigen::Vector4d homo_pos(data.pos[0], data.pos[1], data.pos[2], 1.0);
    Eigen::Vector4d global_pos = init_pose_transform_ * homo_pos;
    data.pos = Eigen::Vector3d(global_pos[0], global_pos[1], global_pos[2]);

    // Apply rotation to velocity
    data.vel = init_pose_transform_rotation_ * data.vel;

    // Apply rotation to accel
    data.accel = init_pose_transform_rotation_ * data.accel;

    // Apply rotation to jerk
    data.jerk = init_pose_transform_rotation_ * data.jerk;

    // Apply yaw
    data.yaw += yaw_init_offset_;
  }

  {
    std::lock_guard<std::mutex> lock(mtx_state_);
    state_ = data;
  }

  if ((state_initialized_ == false || drone_status_ == DroneStatus::YAWING) &&
      (!par_.use_hardware || terminal_goal_initialized_)) {
    // create temporary state
    state tmp;
    tmp.pos = data.pos;
    tmp.yaw = data.yaw;
    previous_yaw_ = data.yaw;

    // Push the state to the plan
    {
      std::lock_guard<std::mutex> lock(mtx_plan_);
      plan_.clear();
      plan_.push_back(tmp);
    }

    // Update Point A
    setA(tmp);

    // Update Point G
    setG(tmp);

    // Update the flag
    state_initialized_ = true;
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Retrieves the next goal (setpoint) from the plan.
 * @param state &next_goal: Output next goal state.
 * @return bool
 */
bool MIGHTY::getNextGoal(state& next_goal) {
  // Check if the planner is initialized
  if (!checkReadyToReplan()) {
    return false;
  }

  // Pop the front of the plan
  next_goal.setZero();

  // Get local copy and potentially pop the plan
  std::deque<state> local_plan;
  {
    std::lock_guard<std::mutex> lock(mtx_plan_);
    local_plan = plan_;

    // If there's more than one goal setpoint, pop the front
    if (plan_.size() > 1) {
      plan_.pop_front();
    }
  }

  // Get the next goal
  next_goal = local_plan.front();

  if (par_.use_hardware && par_.provide_goal_in_global_frame) {
    // Apply transformation to position
    Eigen::Vector4d homo_pos(next_goal.pos[0], next_goal.pos[1], next_goal.pos[2], 1.0);
    Eigen::Vector4d global_pos = init_pose_transform_inv_ * homo_pos;

    // Apply transformation to velocity
    Eigen::Vector3d global_vel = init_pose_transform_rotation_inv_ * next_goal.vel;

    // Apply transformation to accel
    Eigen::Vector3d global_accel = init_pose_transform_rotation_inv_ * next_goal.accel;

    // Apply transformation to jerk
    Eigen::Vector3d global_jerk = init_pose_transform_rotation_inv_ * next_goal.jerk;

    next_goal.pos = Eigen::Vector3d(global_pos[0], global_pos[1], global_pos[2]);
    next_goal.vel = global_vel;
    next_goal.accel = global_accel;
    next_goal.jerk = global_jerk;
  }

  if (!(drone_status_ == DroneStatus::GOAL_REACHED)) {
    // Get the desired yaw
    // If the planner keeps failing, just keep spinning
    if (replanning_failure_count_ > par_.yaw_spinning_threshold) {
      next_goal.yaw = previous_yaw_ + par_.yaw_spinning_dyaw * par_.dc;
      next_goal.dyaw = par_.yaw_spinning_dyaw;
      previous_yaw_ = next_goal.yaw;
    } else {
      // If the local_plan is small just use the previous yaw with no dyaw
      if (local_plan.size() < 5) {
        next_goal.yaw = previous_yaw_;
        next_goal.dyaw = 0.0;
      } else {
        getDesiredYaw(next_goal);
      }
    }

    if (par_.use_hardware && par_.provide_goal_in_global_frame) {
      next_goal.yaw -= yaw_init_offset_;
    }

    next_goal.dyaw = std::clamp(next_goal.dyaw, -par_.w_max, par_.w_max);
  } else {
    next_goal.yaw = previous_yaw_;
    next_goal.dyaw = 0.0;
  }

  // Hardware: zero out yaw tracking (ground robots can't track yaw commands)
  if (par_.use_hardware) {
    next_goal.yaw = 0.0;
    next_goal.dyaw = 0.0;
  }

  return true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Computes the desired yaw for the next goal.
 * @param state &next_goal: Next goal state to update with desired yaw.
 */
void MIGHTY::getDesiredYaw(state& next_goal) {
  double diff = 0.0;
  double desired_yaw = 0.0;

  // Get state
  state local_state;
  getState(local_state);

  // Get G_term
  state G_term;
  {
    std::lock_guard<std::mutex> lock(mtx_G_term_);
    G_term = G_term_;
  }

  switch (drone_status_) {
    case DroneStatus::YAWING:
      if (corridor_hop_yawing_) {
        // Corridor-hop TIP: target is the heading toward the NEXT subgoal,
        // stored in G_.yaw by computeG_corridorHop.
        state G_now;
        getG(G_now);
        desired_yaw = G_now.yaw;
      } else {
        desired_yaw = atan2(G_term.pos[1] - next_goal.pos[1], G_term.pos[0] - next_goal.pos[0]);
      }
      diff = desired_yaw - local_state.yaw;
      // std::cout << "diff1= " << diff << std::endl;
      break;
    case DroneStatus::TRAVELING:
    case DroneStatus::GOAL_SEEN:
      desired_yaw =
          atan2(next_goal.pos[1] - local_state.pos.y(), next_goal.pos[0] - local_state.pos.x());
      diff = desired_yaw - local_state.yaw;
      next_goal.yaw = desired_yaw;
      break;
    case DroneStatus::GOAL_REACHED:
      next_goal.dyaw = 0.0;
      next_goal.yaw = previous_yaw_;
      return;
  }

  mighty_utils::angle_wrap(diff);
  if (fabs(diff) < 0.04 && drone_status_ == DroneStatus::YAWING) {
    // Leaving YAWING: clear the corridor-hop override and drop the held
    // subgoal so needReplan picks a fresh one for the next leg.
    corridor_hop_yawing_ = false;
    held_subgoal_valid_ = false;
    changeDroneStatus(DroneStatus::TRAVELING);
  }

  yaw(diff, next_goal);
}

// ----------------------------------------------------------------------------

void MIGHTY::yaw(double diff, state& next_goal) {
  saturate(diff, -par_.dc * par_.w_max, par_.dc * par_.w_max);
  dyaw_filtered_ = (1 - par_.alpha_filter_dyaw) * (copysign(1, diff) * par_.w_max) +
                   par_.alpha_filter_dyaw * dyaw_filtered_;
  next_goal.dyaw = dyaw_filtered_;
  next_goal.yaw = previous_yaw_ + dyaw_filtered_ * par_.dc;
  previous_yaw_ = next_goal.yaw;
}

// ----------------------------------------------------------------------------

/**
 * @brief If the goal is inside an occupied voxel, relocate it to the nearest
 *        free or unknown cell, then push outward along that direction by
 *        ||drone_bbox|| so the new goal sits clear of the obstacle. No-op if
 *        the goal is already non-occupied or if the map isn't ready yet.
 */
bool MIGHTY::sanitizeTerminalGoal(state& goal) {
  // Map not initialized yet (e.g. very first goal at startup): trust the user.
  if (!hgp_manager_.isMapInitialized()) return true;

  // Already free or unknown: nothing to do.
  if (!hgp_manager_.checkIfPointOccupied(goal.pos)) return true;

  // BFS for the closest free/unknown cell.
  Vec3f closest;
  hgp_manager_.findClosestNonOccupiedPoint(goal.pos, closest);

  // findClosestNonOccupiedPoint returns the input unchanged if it found nothing.
  if ((closest - goal.pos).norm() < 1e-6) {
    printf("[MIGHTY] sanitizeTerminalGoal: goal is occupied and no free/unknown "
           "cell found within BFS budget; dropping goal.\n");
    return false;
  }

  // Direction outward from the occupied goal toward free/unknown space.
  const Vec3f dir = (closest - goal.pos).normalized();

  // Required clearance from the original occupied point: ||drone_bbox||.
  // drone_bbox is a [x,y,z] full-extent vector; using its norm gives the
  // bounding-box diagonal as the safety margin.
  double clearance = 0.0;
  if (par_.drone_bbox.size() >= 3) {
    clearance = std::sqrt(par_.drone_bbox[0] * par_.drone_bbox[0] +
                          par_.drone_bbox[1] * par_.drone_bbox[1] +
                          par_.drone_bbox[2] * par_.drone_bbox[2]);
  }

  Vec3f new_pos = goal.pos + dir * clearance;

  // Safety loop: if pushing by `clearance` lands us back inside an obstacle
  // (thin wall, etc.), keep walking outward by one resolution at a time.
  // Cap iterations so we never spin forever.
  const double step = par_.res > 0.0 ? par_.res : 0.1;
  for (int i = 0; i < 20 && hgp_manager_.checkIfPointOccupied(new_pos); ++i) {
    new_pos += dir * step;
  }

  if (hgp_manager_.checkIfPointOccupied(new_pos)) {
    printf("[MIGHTY] sanitizeTerminalGoal: could not escape occupied region "
           "after %d steps; dropping goal.\n", 20);
    return false;
  }

  const Vec3f original = goal.pos;
  goal.pos = new_pos;
  printf("[MIGHTY] Goal relocated from (%.2f,%.2f,%.2f) to (%.2f,%.2f,%.2f) "
         "(was inside occupied voxel, clearance=%.2f m)\n",
         original.x(), original.y(), original.z(),
         new_pos.x(), new_pos.y(), new_pos.z(), clearance);
  return true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the terminal goal.
 * @param const state &term_goal: Desired terminal goal state.
 */
void MIGHTY::setTerminalGoal(const state& term_goal) {
  // Get the state
  state local_state;
  getState(local_state);

  // Set the terminal goal
  setGterm(term_goal);

  // Project the terminal goal to the sphere
  {
    std::lock_guard<std::mutex> lock(mtx_G_);
    G_.pos = mighty_utils::projectPointToSphere(local_state.pos, term_goal.pos, par_.horizon);
  }

  changeDroneStatus(DroneStatus::TRAVELING);

  if (!terminal_goal_initialized_) terminal_goal_initialized_ = true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Changes the drone's status (YAWING, TRAVELING, GOAL_SEEN, GOAL_REACHED).
 * @param int new_status: New status value.
 */
void MIGHTY::changeDroneStatus(int new_status) {
  if (new_status == drone_status_) return;

  std::cout << "Changing DroneStatus from ";

  switch (drone_status_) {
    case DroneStatus::YAWING:
      std::cout << bold << "status_=YAWING" << reset;
      break;
    case DroneStatus::TRAVELING:
      std::cout << bold << "status_=TRAVELING" << reset;
      break;
    case DroneStatus::GOAL_SEEN:
      std::cout << bold << "status_=GOAL_SEEN" << reset;
      break;
    case DroneStatus::GOAL_REACHED:
      std::cout << bold << "status_=GOAL_REACHED" << reset;
      break;
  }

  std::cout << " to ";

  switch (new_status) {
    case DroneStatus::YAWING:
      std::cout << bold << "status_=YAWING" << reset;
      break;
    case DroneStatus::TRAVELING:
      std::cout << bold << "status_=TRAVELING" << reset;
      break;
    case DroneStatus::GOAL_SEEN:
      std::cout << bold << "status_=GOAL_SEEN" << reset;
      break;
    case DroneStatus::GOAL_REACHED:
      std::cout << bold << "status_=GOAL_REACHED" << reset;
      break;
  }

  std::cout << std::endl;

  drone_status_ = new_status;
}

// ----------------------------------------------------------------------------

/**
 * @brief Checks if all necessary components are initialized.
 * @return bool
 */
bool MIGHTY::checkReadyToReplan() {
  bool is_ready = state_initialized_ && terminal_goal_initialized_ &&
                  hgp_manager_.isMapInitialized() &&
                  (!par_.use_hardware || (kdtree_map_initialized_
                                          // && kdtree_unk_initialized_
                                          ));

  return is_ready;
}

// ----------------------------------------------------------------------------

void MIGHTY::extractDynamicHeatData(const std::vector<std::shared_ptr<dynTraj>>& trajs,
                                    double current_time, const Eigen::Vector3d& agent_pos,
                                    vec_Vecf<3>& obst_pos, vec_Vecf<3>& obst_bbox,
                                    double& traj_max_time) {
  if (!par_.use_heat_map || !par_.dynamic_heat_enabled || trajs.empty()) return;

  // Filter obstacles within planning horizon
  std::vector<std::shared_ptr<dynTraj>> selected_trajs;
  selected_trajs.reserve(trajs.size());

  for (const auto& traj : trajs) {
    if (!traj) continue;
    Eigen::Vector3d pos = traj->eval(current_time);
    double dist = (pos - agent_pos).norm();
    if (dist > par_.horizon) continue;

    obst_pos.push_back(pos.cast<double>());
    Eigen::Vector3d bbox_half = traj->bbox / 2.0;
    obst_bbox.push_back(bbox_half.cast<double>());
    selected_trajs.push_back(traj);
  }

  // Compute prediction horizon: use the max of config prediction_horizon and
  // the actual trajectory duration (so agent trajectories are fully covered)
  double Th = par_.prediction_horizon;
  for (const auto& traj : selected_trajs) {
    if (traj->mode == dynTraj::Mode::Piecewise && traj->pwp.times.size() >= 2) {
      double remaining = traj->pwp.times.back() - std::max(traj->pwp.times.front(), current_time);
      Th = std::max(Th, std::max(0.0, remaining));
    }
  }
  traj_max_time = Th;

  // Sample predicted positions along each trajectory at future times from NOW
  if (par_.heat_num_samples > 0 && Th > 0.0 && !selected_trajs.empty()) {
    const int M = par_.heat_num_samples;

    std::vector<float> pred_times(M);
    for (int j = 0; j < M; ++j) {
      const double a = (M == 1) ? 0.0 : (double)j / (double)(M - 1);
      pred_times[j] = static_cast<float>(a * Th);
    }

    std::vector<vec_Vecf<3>> pred_samples(selected_trajs.size());
    for (size_t k = 0; k < selected_trajs.size(); ++k) {
      pred_samples[k].resize(M);
      for (int j = 0; j < M; ++j) {
        const double t_abs = current_time + (double)pred_times[j];
        Eigen::Vector3d pk = selected_trajs[k]->eval(t_abs);
        if (!std::isfinite(pk.x()) || !std::isfinite(pk.y()) || !std::isfinite(pk.z())) {
          pk = selected_trajs[k]->eval(current_time);
        }
        pred_samples[k][j] = pk;
      }
    }

    hgp_manager_.setDynamicPredictedSamples(pred_samples, pred_times);
  }
}

// ----------------------------------------------------------------------------

void MIGHTY::updateMap(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_map,
                       const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_unk) {
  // 1) Atomically store the incoming clouds
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_map_);
    pclptr_map_ = pclptr_map;
  }
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_unk_);
    pclptr_unk_ = pclptr_unk;
  }

  // Update the map size
  state local_state, local_G;
  getState(local_state);
  getG(local_G);
  computeMapSize(local_state.pos, local_G.pos);

  // 2) Extract obstacle data from dynamic trajectories for HGP heat map
  auto trajs = getTrajs();
  double current_time = local_state.t;
  vec_Vecf<3> obst_pos, obst_bbox;
  double traj_max_time = 0.0;

  extractDynamicHeatData(trajs, current_time, local_state.pos, obst_pos, obst_bbox, traj_max_time);

  // 3) Call HGP updateMap with extracted obstacle data
  hgp_manager_.updateMap(wdx_, wdy_, wdz_, map_center_, pclptr_map_, pclptr_unk_, obst_pos,
                         obst_bbox, traj_max_time);

  // 4) Known‐space KD‐tree
  if (pclptr_map_ && !pclptr_map_->points.empty()) {
    std::lock_guard<std::mutex> lk(mtx_kdtree_map_);
    kdtree_map_.setInputCloud(pclptr_map_);
    kdtree_map_initialized_ = true;
    hgp_manager_.updateVecOccupied(pclptr_to_vec(pclptr_map_));
  } else {
    RCLCPP_WARN(rclcpp::get_logger("mighty"),
                "updateMap: member pclptr_map_ was null or empty; skipping KD-tree update");
  }

  // 5) Unknown‐space KD‐tree
  if (pclptr_unk_ && !pclptr_unk_->points.empty()) {
    std::lock_guard<std::mutex> lk(mtx_kdtree_unk_);
    kdtree_unk_.setInputCloud(pclptr_unk_);
    kdtree_unk_initialized_ = true;
    hgp_manager_.updateVecUnknownOccupied(pclptr_to_vec(pclptr_unk_));
    hgp_manager_.insertVecOccupiedToVecUnknownOccupied();
  } else {
    RCLCPP_WARN(rclcpp::get_logger("mighty"),
                "updateMap: member pclptr_unk_ was null or empty; skipping KD‐tree update");
  }
}

// ----------------------------------------------------------------------------

void MIGHTY::updateOccupancyMap(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_map) {
  // 1) Atomically store the incoming clouds
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_map_);
    pclptr_map_ = pclptr_map;
  }

  // Update the map size
  state local_state, local_G;
  getState(local_state);
  getG(local_G);
  computeMapSize(local_state.pos, local_G.pos);

  // 2) Extract obstacle data from dynamic trajectories for HGP heat map
  auto trajs = getTrajs();
  double current_time = local_state.t;
  vec_Vecf<3> obst_pos, obst_bbox;
  double traj_max_time = 0.0;

  extractDynamicHeatData(trajs, current_time, local_state.pos, obst_pos, obst_bbox, traj_max_time);

  // No unknown cloud for occupancy-only update
  pcl::PointCloud<pcl::PointXYZ>::Ptr empty_unk(new pcl::PointCloud<pcl::PointXYZ>());

  hgp_manager_.updateMap(wdx_, wdy_, wdz_, map_center_, pclptr_map_, empty_unk, obst_pos, obst_bbox,
                         traj_max_time);

  // 3) Known‐space KD‐tree
  if (pclptr_map_ && !pclptr_map_->points.empty()) {
    std::lock_guard<std::mutex> lk(mtx_kdtree_map_);
    kdtree_map_.setInputCloud(pclptr_map_);
    kdtree_map_initialized_ = true;
    hgp_manager_.updateVecOccupied(pclptr_to_vec(pclptr_map_));
  } else {
    RCLCPP_WARN(rclcpp::get_logger("mighty"),
                "updateMap: member pclptr_map_ was null or empty; skipping KD-tree update");
  }
}

// ----------------------------------------------------------------------------

void MIGHTY::updateUnknownCloud(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr_unk) {
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_unk_);
    pclptr_unk_ = pclptr_unk;
  }
  if (pclptr_unk_ && !pclptr_unk_->points.empty()) {
    std::lock_guard<std::mutex> lk(mtx_kdtree_unk_);
    kdtree_unk_.setInputCloud(pclptr_unk_);
    hgp_manager_.updateVecUnknownOccupied(pclptr_to_vec(pclptr_unk_));
    hgp_manager_.insertVecOccupiedToVecUnknownOccupied();
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Set the initial pose.
 * @param const geometry_msgs::msg::TransformStamped &init_pose: Initial pose.
 */
void MIGHTY::setInitialPose(const geometry_msgs::msg::TransformStamped& init_pose) {
  init_pose_ = init_pose;

  // First compute transformation matrix from init_pose_ (geometry_msgs::msg::TransformStamped)
  Eigen::Matrix4d init_pose_transform = Eigen::Matrix4d::Identity();
  Eigen::Quaterniond init_pose_quat(
      init_pose_.transform.rotation.w, init_pose_.transform.rotation.x,
      init_pose_.transform.rotation.y, init_pose_.transform.rotation.z);
  Eigen::Vector3d init_pose_translation(init_pose_.transform.translation.x,
                                        init_pose_.transform.translation.y,
                                        init_pose_.transform.translation.z);
  init_pose_transform.block<3, 3>(0, 0) = init_pose_quat.toRotationMatrix();
  init_pose_transform.block<3, 1>(0, 3) = init_pose_translation;

  // Get initial pose
  init_pose_transform_ = init_pose_transform;
  init_pose_transform_rotation_ = init_pose_quat.toRotationMatrix();
  yaw_init_offset_ =
      std::atan2(init_pose_transform_rotation_(1, 0), init_pose_transform_rotation_(0, 0));

  std::cout << bold << green << "yaw_init_offset_: " << yaw_init_offset_ << reset << std::endl;

  // Get the inverse of init_pose_ (geometry_msgs::msg::TransformStamped)
  init_pose_transform_inv_ = init_pose_transform.inverse();
  init_pose_transform_rotation_inv_ = init_pose_quat.toRotationMatrix().inverse();
  // yaw_init_offset_ = std::atan2(init_pose_transform_rotation_inv_(1, 0),
  // init_pose_transform_rotation_inv_(0, 0));
}

// ----------------------------------------------------------------------------

// Apply the initial pose transformation to the pwp
void MIGHTY::applyInitiPoseTransform(PieceWisePol& pwp) {
  // Loop thru the intervals
  for (int i = 0; i < pwp.coeff_x.size(); i++) {
    // Loop thru a, b, c, and d
    for (int j = 0; j < 4; j++) {
      Eigen::Vector4d coeff;
      coeff[0] = pwp.coeff_x[i][j];
      coeff[1] = pwp.coeff_y[i][j];
      coeff[2] = pwp.coeff_z[i][j];
      coeff[3] = 1.0;

      // Apply multiplication
      coeff = init_pose_transform_ * coeff;

      // cout agent frame pose
      pwp.coeff_x[i][j] = coeff[0];
      pwp.coeff_y[i][j] = coeff[1];
      pwp.coeff_z[i][j] = coeff[2];
    }
  }
}

// ----------------------------------------------------------------------------

// Apply the inverse of initial pose transformation to the pwp
void MIGHTY::applyInitiPoseInverseTransform(PieceWisePol& pwp) {
  // Loop thru the intervals
  for (int i = 0; i < pwp.coeff_x.size(); i++) {
    // Loop thru a, b, c, and d
    for (int j = 0; j < 4; j++) {
      Eigen::Vector4d coeff;
      coeff[0] = pwp.coeff_x[i][j];
      coeff[1] = pwp.coeff_y[i][j];
      coeff[2] = pwp.coeff_z[i][j];
      coeff[3] = 1.0;

      // Apply multiplication
      coeff = init_pose_transform_inv_ * coeff;

      pwp.coeff_x[i][j] = coeff[0];
      pwp.coeff_y[i][j] = coeff[1];
      pwp.coeff_z[i][j] = coeff[2];
    }
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Checks if the goal is reached.
 * @return bool
 */
bool MIGHTY::goalReachedCheck() {
  if (checkReadyToReplan() && drone_status_ == DroneStatus::GOAL_REACHED) {
    return true;
  }
  return false;
}