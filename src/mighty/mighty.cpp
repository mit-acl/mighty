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
MIGHTY::MIGHTY(parameters par) : par_(par)
{

  // Set up dgp_manager
  dgp_manager_.setParameters(par_);

  // Set up the planner parameters (TODO: move to parameters)
  planner_params_.verbose = false;                                 // enable verbose output
  planner_params_.V_max = par_.v_max;                              // max velocity
  planner_params_.A_max = par_.a_max;                              // max acceleration
  planner_params_.J_max = par_.j_max;                              // max jerk
  planner_params_.num_perturbation = par_.num_perturbation_for_ig; // number of perturbations for initial guesses
  planner_params_.r_max = par_.r_max_for_ig;                       // perturbation radius for initial guesses
  planner_params_.time_weight = par_.time_weight;                  // weight for time cost
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
  planner_params_.num_dyn_obst_samples = par_.num_dyn_obst_samples; // Number of dynamic obstacle samples
  planner_params_.Co = par_.planner_Co;                             // for static obstacle avoidance
  planner_params_.Cw = par_.planner_Cw;                             // for dynamic obstacle avoidance
  planner_params_.BIG = 1e8;
  planner_params_.dc = par_.dc;                                             // descretiation constant
  planner_params_.init_turn_bf = par_.init_turn_bf;

  // Set up the L-BFGS parameters
  lbfgs_params_.mem_size = 256;
  lbfgs_params_.min_step = 1.0e-32;
  lbfgs_params_.f_dec_coeff = par_.f_dec_coeff;         // allow larger Armijo steps
  lbfgs_params_.cautious_factor = par_.cautious_factor; // always accept BFGS update
  lbfgs_params_.past = par_.past;
  lbfgs_params_.max_linesearch = par_.max_linesearch; // fewer backtracking tries
  lbfgs_params_.max_iterations = par_.max_iterations; // allow more iterations
  lbfgs_params_.g_epsilon = par_.g_epsilon;
  lbfgs_params_.delta = par_.delta; // stop once f-improvement is minimal

  // Set up unconstrained optimization solver for whole trajectory
  whole_traj_solver_ptr_ = std::make_shared<lbfgs::SolverLBFGS>();
  whole_traj_solver_ptr_->initializeSolver(planner_params_);

  // Set up basis converter
  BasisConverter basis_converter;
  A_rest_pos_basis_ = basis_converter.getArestMinvo(); // Use Minvo basis
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
void MIGHTY::startAdaptKValue()
{

  // Compute the average computation time
  for (int i = 0; i < store_computation_times_.size(); i++)
  {
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
void MIGHTY::computeG(const state &A, const state &G_term, double horizon)
{
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

/**
 * @brief Checks if we need to replan.
 * @return bool
 */

bool MIGHTY::needReplan(const state &local_state, const state &local_G_term, const state &last_plan_state)
{

  // Compute the distance to the terminal goal
  double dist_to_term_G = (local_state.pos - local_G_term.pos).norm();
  double dist_from_last_plan_state_to_term_G = (last_plan_state.pos - local_G_term.pos).norm();

  if (dist_to_term_G < par_.goal_radius)
  {
    changeDroneStatus(DroneStatus::GOAL_REACHED);
    return false;
  }

  if (dist_to_term_G < par_.goal_seen_radius)
  {
    changeDroneStatus(DroneStatus::GOAL_SEEN); // This triggers to use the hard final state constraint
  }

  if (drone_status_ == DroneStatus::GOAL_SEEN && dist_from_last_plan_state_to_term_G < par_.goal_radius)
  {
    return false;
  }

  // Don't plan if drone is not traveling
  if (drone_status_ == DroneStatus::GOAL_REACHED || (drone_status_ == DroneStatus::YAWING))
    return false;

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::findAandAtime(state &A, double &A_time, double current_time, double last_replaning_computation_time)
{

  mtx_plan_.lock();
  int plan_size = plan_.size();
  mtx_plan_.unlock();

  if (plan_size == 0)
  {
    std::cout << bold << red << "plan_size == 0" << reset << std::endl;
    return false;
  }

  if (par_.use_state_update)
  {
    // Change k_value dynamically
    // To get stable results, we will use a default value of k_value until we have enough computation time
    if (!use_adapt_k_value_)
    {
      // Use default k_value
      k_value_ = std::max((int)plan_size - par_.default_k_value, 0);

      // Store computation times
      if (num_replanning_ != 1) // Don't store the very first computation time (because we don't have a previous computation time)
        store_computation_times_.push_back(last_replaning_computation_time);
    }
    else
    {

      // Computation time filtering
      est_comp_time_ = par_.alpha_k_value_filtering * last_replaning_computation_time + (1 - par_.alpha_k_value_filtering) * est_comp_time_;

      // Get state number based on est_comp_time_ and dc
      k_value_ = std::max((int)plan_size - (int)(par_.k_value_factor * est_comp_time_ / par_.dc), 0);
    }

    // Check if k_value_ is valid
    if (plan_size - 1 - k_value_ < 0 || plan_size - 1 - k_value_ >= plan_size)
    {
      k_value_ = plan_size - 1; // If k_value_ is larger than the plan size, we set it to the last state
    }

    // Get A
    mtx_plan_.lock();
    A = plan_[plan_size - 1 - k_value_];
    mtx_plan_.unlock();

    // Get A_time
    A_time = current_time + (plan_size - 1 - k_value_) * par_.dc; // time to A from current_pos is (plan_size - 1 - k_value_) * par_.dc;
  }
  else // If we don't update state - this is for global planner benchmarking purposes
  {
    // Get state
    getState(A);
    A_time = current_time;
  }

  // Check if A is within the map (especially for z)
  if (A.pos[2] < par_.z_min || A.pos[2] > par_.z_max, A.pos[0] < par_.x_min || A.pos[0] > par_.x_max, A.pos[1] < par_.y_min || A.pos[1] > par_.y_max)
  {
    printf("A (%f, %f, %f) is out of the map\n", A.pos[0], A.pos[1], A.pos[2]);
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::checkIfPointOccupied(const Vec3f &point)
{
  // Check if the point is free
  return dgp_manager_.checkIfPointOccupied(point);
}

// ----------------------------------------------------------------------------

bool MIGHTY::checkIfPointFree(const Vec3f &point)
{
  // Check if the point is free
  return dgp_manager_.checkIfPointFree(point);
}

// ----------------------------------------------------------------------------

bool MIGHTY::getSafeCorridor(const vec_Vecf<3> &global_path, const state &A)
{

  // Timer for computing the safe corridor
  MyTimer cvx_decomp_timer(true);

  // Debug
  if (par_.debug_verbose)
    std::cout << "Convex decomposition" << std::endl;

  // Check if the convex decomposition failed
  if (!dgp_manager_.cvxEllipsoidDecomp(A, global_path, safe_corridor_polytopes_whole_, poly_out_whole_, false))
  {
    std::cout << bold << red << "Convex decomposition failed" << reset << std::endl;
    poly_out_whole_.clear();
    poly_out_safe_.clear();
    return false;
  }

  // Get computation time [ms]
  cvx_decomp_time_ = cvx_decomp_timer.getElapsedMicros() / 1000.0;

  return true;
}

// ----------------------------------------------------------------------------

void MIGHTY::computeMapSize(const Eigen::Vector3d &min_pos, const Eigen::Vector3d &max_pos)
{

  // Get local_A
  state local_A;
  getA(local_A);

  // Increase the effective buffer size based on the number of DGP failures.
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

  // Compute the base map center as the midpoint between the min and max positions.
  map_center_ = (min_pos + max_pos) / 2.0;

}

// ----------------------------------------------------------------------------

bool MIGHTY::checkPointWithinMap(const Eigen::Vector3d &point) const
{
  // Check if the point is within the map boundaries for each axis
  return (std::abs(point[0] - map_center_[0]) <= wdx_ / 2.0) && (std::abs(point[1] - map_center_[1]) <= wdy_ / 2.0) && (std::abs(point[2] - map_center_[2]) <= wdz_ / 2.0);
}

// ----------------------------------------------------------------------------

void MIGHTY::getStaticPushPoints(vec_Vecf<3> &static_push_points)
{
  static_push_points = static_push_points_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getLocalGlobalPath(vec_Vecf<3> &local_global_path, vec_Vecf<3> &local_global_path_after_push)
{
  local_global_path = local_global_path_;
  local_global_path_after_push = local_global_path_after_push_;
}

// ----------------------------------------------------------------------------

void MIGHTY::getGlobalPath(vec_Vecf<3> &global_path)
{
  mtx_global_path_.lock();
  global_path = global_path_;
  mtx_global_path_.unlock();
}

// ----------------------------------------------------------------------------

void MIGHTY::getOriginalGlobalPath(vec_Vecf<3> &original_global_path)
{
  mtx_original_global_path_.lock();
  original_global_path = original_global_path_;
  mtx_original_global_path_.unlock();
}

// ----------------------------------------------------------------------------

void MIGHTY::getFreeGlobalPath(vec_Vecf<3> &free_global_path)
{
  free_global_path = free_global_path_;
}

// ----------------------------------------------------------------------------

void MIGHTY::resetData()
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

  poly_out_whole_.clear();
  poly_out_safe_.clear();
  goal_setpoints_.clear();
  // pwp_to_share_.clear();
  optimal_yaw_sequence_.clear();
  yaw_control_points_.clear();
  yaw_knots_.clear();
  cps_.clear();
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveData(double &final_g,
                         double &global_planning_time,
                         double &dgp_static_jps_time,
                         double &dgp_check_path_time,
                         double &dgp_dynamic_astar_time,
                         double &dgp_recover_path_time,
                         double &cvx_decomp_time,
                         double &initial_guess_computation_time,
                         double &local_traj_computatoin_time,
                         double &safety_check_time,
                         double &safe_paths_time,
                         double &yaw_sequence_time,
                         double &yaw_fitting_time)
{
  final_g = final_g_;
  global_planning_time = global_planning_time_;
  dgp_static_jps_time = dgp_static_jps_time_;
  dgp_check_path_time = dgp_check_path_time_;
  dgp_dynamic_astar_time = dgp_dynamic_astar_time_;
  dgp_recover_path_time = dgp_recover_path_time_;
  cvx_decomp_time = cvx_decomp_time_;
  initial_guess_computation_time = initial_guess_computation_time_;
  local_traj_computatoin_time = local_traj_computation_time_;
  safe_paths_time = safe_paths_time_;
  safety_check_time = safety_check_time_;
  yaw_sequence_time = yaw_sequence_time_;
  yaw_fitting_time = yaw_fitting_time_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrievePolytopes(vec_E<Polyhedron<3>> &poly_out_whole, vec_E<Polyhedron<3>> &poly_out_safe)
{
  poly_out_whole = poly_out_whole_;
  poly_out_safe = poly_out_safe_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveGoalSetpoints(std::vector<state> &goal_setpoints)
{
  goal_setpoints = goal_setpoints_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveListSubOptGoalSetpoints(std::vector<std::vector<state>> &list_subopt_goal_setpoints)
{
  list_subopt_goal_setpoints = list_subopt_goal_setpoints_;
}

// ----------------------------------------------------------------------------

void MIGHTY::retrieveCPs(std::vector<Eigen::Matrix<double, 3, 6>> &cps)
{
  cps = cps_;
}

// ----------------------------------------------------------------------------

/**
 * @brief Replans the trajectory.
 * @param double last_replaning_computation_time: Last replanning computation time.
 * @param double current_time: Current timestamp.
 */
std::tuple<bool, bool> MIGHTY::replan(double last_replaning_computation_time, double current_time)
{

  /* -------------------- Housekeeping -------------------- */

  MyTimer timer_housekeeping(true);

  // Reset Data
  resetData();

  // Check if we need to replan
  if (!checkReadyToReplan())
  {
    std::cout << bold << red << "Planner is not ready to replan" << reset << std::endl;
    return std::make_tuple(false, false);
  }

  // Get states we need
  state local_state, local_G_term, last_plan_state;
  getState(local_state);
  getGterm(local_G_term);
  getLastPlanState(last_plan_state);

  // Check if we need to replan based on the distance to the terminal goal
  if (!needReplan(local_state, local_G_term, last_plan_state))
    return std::make_tuple(false, false);

  if (par_.debug_verbose)
    std::cout << "Housekeeping: " << timer_housekeeping.getElapsedMicros() / 1000.0 << " ms" << std::endl;

  /* -------------------- Global Planning -------------------- */

  MyTimer timer_global(true);
  vec_Vecf<3> global_path;
  if (!generateGlobalPath(global_path, current_time, last_replaning_computation_time))
  {
    if (par_.debug_verbose)
      std::cout << "Global Planning: " << timer_global.getElapsedMicros() / 1000.0 << " ms" << std::endl;
    return std::make_tuple(false, false);
  }
  if (par_.debug_verbose)
    std::cout << "Global Planning: " << timer_global.getElapsedMicros() / 1000.0 << " ms" << std::endl;

  /* -------------------- Local Trajectory Optimization -------------------- */

  MyTimer timer_local(true);
  if (!planLocalTrajectory(global_path))
  {
    if (par_.debug_verbose)
      std::cout << "Local Trajectory Optimization: " << timer_local.getElapsedMicros() / 1000.0 << " ms" << std::endl;
    return std::make_tuple(false, true);
  }
  if (par_.debug_verbose)
    std::cout << "Local Trajectory Optimization: " << timer_local.getElapsedMicros() / 1000.0 << " ms" << std::endl;

  /* -------------------- Append to Plan -------------------- */

  MyTimer timer_append(true);
  if (!appendToPlan())
  {
    if (par_.debug_verbose)
      std::cout << "Append to Plan: " << timer_append.getElapsedMicros() / 1000.0 << " ms" << std::endl;
    return std::make_tuple(false, true);
  }
  if (par_.debug_verbose)
    std::cout << "Append to Plan: " << timer_append.getElapsedMicros() / 1000.0 << " ms" << std::endl;

  /* -------------------- Final Housekeeping -------------------- */

  MyTimer timer_final(true);

  if (par_.debug_verbose)
    std::cout << bold << green << "Replanning succeeded" << reset << std::endl;

  // Reset the replanning failure count
  replanning_failure_count_ = 0;
  if (par_.debug_verbose)
    std::cout << "Final Housekeeping: " << timer_final.getElapsedMicros() / 1000.0 << " ms" << std::endl;

  return std::make_tuple(true, true);
}

// ----------------------------------------------------------------------------

bool MIGHTY::generateGlobalPath(vec_Vecf<3> &global_path, double current_time, double last_replaning_computation_time)
{

  // Get G and G_term
  state local_G, local_G_term;
  getG(local_G);
  getGterm(local_G_term);

  // Declare local variables
  state local_A;
  double A_time;

  // Find A and A_time
  if (!findAandAtime(local_A, A_time, current_time, last_replaning_computation_time))
  {
    replanning_failure_count_++;
    return false;
  }

  // Set A and A_time
  setA(local_A);
  setA_time(A_time);

  // Compute G
  computeG(local_A, local_G_term, par_.horizon);

  // Set up the DGP planner (since updateVmax() needs to be called after setupDGPPlanner, we use v_max_ from the last replan)
  dgp_manager_.setupDGPPlanner(par_.global_planner, par_.global_planner_verbose, map_res_, v_max_, par_.a_max, par_.j_max, par_.dgp_timeout_duration_ms, par_.w_unknown, par_.w_align, par_.decay_len_cells, par_.w_side, par_.los_cells, par_.min_len, par_.min_turn);

  // Free start and goal if necessary
  if (par_.use_free_start)
    dgp_manager_.freeStart(local_A.pos, par_.free_start_factor);
  if (par_.use_free_goal)
    dgp_manager_.freeGoal(local_G.pos, par_.free_goal_factor);

  // Debug
  if (par_.debug_verbose)
    std::cout << "Solving DGP" << std::endl;

  // if using ground robot, we fix the z
  if (par_.vehicle_type != "uav")
  {
    local_A.pos[2] = 1.0;
    local_G.pos[2] = 1.0;
  }

  // 1) Build a direction hint from the *previous* global path
  vec_Vecf<3> prev_global;
  getGlobalPath(prev_global); // last successful global path

  Eigen::Vector3d dir_hint = (local_G.pos - local_A.pos).normalized();
  if (prev_global.size() >= 2)
  {
    Eigen::Vector3d s0 = prev_global[0];
    Eigen::Vector3d s1 = prev_global[1];
    Eigen::Vector3d seg = s1 - s0;
    if (seg.norm() > 1e-8)
    {
      dir_hint = seg.normalized();
    }
  }
  else
  {
    dir_hint = local_G.pos - local_A.pos;
    if (dir_hint.norm() > 1e-8)
      dir_hint.normalize();
  }

  // Keep ground robots planar
  if (par_.vehicle_type != "uav")
    dir_hint[2] = 0.0;

  // 2) Use this as the "start_vel" argument (magnitude doesn't matter; we use the direction)
  Vec3f start_dir_hint(dir_hint.x(), dir_hint.y(), dir_hint.z());

  // Solve DGP
  // if (!dgp_manager_.solveDGP(local_A.pos, local_A.vel, local_G.pos, final_g_, par_.global_planner_huristic_weight, A_time, global_path))
  if (!dgp_manager_.solveDGP(local_A.pos, start_dir_hint, local_G.pos, final_g_, par_.global_planner_huristic_weight, A_time, global_path))
  {
    if (par_.debug_verbose)
      std::cout << bold << red << "DGP did not find a solution" << reset << std::endl;
    dgp_failure_count_++;
    replanning_failure_count_++;
    return false;
  }

  // use this for map resizing
  mtx_global_path_.lock();
  global_path_ = global_path;
  mtx_global_path_.unlock();

  // For visualization
  mtx_original_global_path_.lock();
  original_global_path_ = global_path;
  mtx_original_global_path_.unlock();

  // Debug
  if (par_.debug_verbose)
    std::cout << "global_path.size(): " << global_path.size() << std::endl;

  // Get computation time
  dgp_manager_.getComputationTime(global_planning_time_, dgp_static_jps_time_, dgp_check_path_time_, dgp_dynamic_astar_time_, dgp_recover_path_time_);

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::planLocalTrajectory(vec_Vecf<3> &global_path)
{

  // Get local_A, local_G and A_time
  state local_A, local_G;
  double A_time;
  getA(local_A);
  getG(local_G);
  getA_time(A_time);

  // If the global path's size is < 3 after trimming, we cannot proceed
  if (global_path.empty() || global_path.size() < 3)
  {
    std::cout << bold << red << "Global path's size is < 3 after trimming" << reset << std::endl;
    replanning_failure_count_++;
    return false;
  }

  // convex decomposition
  if (!getSafeCorridor(global_path, local_A))
  {
    replanning_failure_count_++;
    return false;
  }

  // Initialize flag
  bool optimization_succeeded = false;

  if (par_.vehicle_type != "uav")
  {
    local_A.pos[2] = 1.0;
  }

  optimization_succeeded = generateLocalTrajectory(
      local_A, A_time, global_path, initial_guess_computation_time_,
      local_traj_computation_time_, whole_traj_solver_ptr_);

  if (par_.debug_verbose)
  {
    std::cout << "initial_guess_computation_time_: " << initial_guess_computation_time_ << " ms" << std::endl;
    std::cout << "Local trajectory optimization finished" << std::endl;
  }

  if (optimization_succeeded)
  {

    // For the optimal solution
    whole_traj_solver_ptr_->reconstructPVATCPopt(zopt_); // First recover the final control points and times from z_opt
    whole_traj_solver_ptr_->getGoalSetpoints(goal_setpoints_);

    // print out the goal setpoints
    // whole_traj_solver_ptr_->getPieceWisePol(pwp_to_share_);
    whole_traj_solver_ptr_->getControlPoints(cps_); // Bezier control points

    // Get goal setpoints for suboptimal solutions for visualization
    if (par_.use_multiple_initial_guesses)
    {

      if (par_.debug_verbose)
        std::cout << "Size of list_z_subopt_: " << list_z_subopt_.size() << std::endl;

      // Initialize list_subopt_goal_setpoints_
      list_subopt_goal_setpoints_.clear();

      // Loop over list_z_subopt
      for (int idx = 0; idx < list_z_subopt_.size(); ++idx)
      {
        // Reconstruct control points and times for each suboptimal solution
        whole_traj_solver_ptr_->reconstructPVATCPopt(list_z_subopt_[idx]); // First recover the final control points and times from z_opt
        std::vector<state> subopt_goal_setpoints;
        whole_traj_solver_ptr_->getGoalSetpoints(subopt_goal_setpoints);
        list_subopt_goal_setpoints_.push_back(subopt_goal_setpoints);
      }
    }
  }
  else
  {
    replanning_failure_count_++;
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::generateLocalTrajectory(const state &local_A, double A_time,
                                    vec_Vec3f &global_path,
                                    double &initial_guess_computation_time,
                                    double &local_traj_computation_time,
                                    std::shared_ptr<lbfgs::SolverLBFGS> &whole_traj_solver_ptr)
{

  if (par_.debug_verbose)
    std::cout << "Preparing solver for replan" << std::endl;

  std::vector<std::shared_ptr<dynTraj>> local_trajs;
  getTrajs(local_trajs);

  // Get local_G
  state local_G;
  getG(local_G);

  state local_E;
  Vec3f mean_point;

  if (drone_status_ == DroneStatus::GOAL_REACHED || drone_status_ == DroneStatus::GOAL_SEEN)
  {
    local_E = local_G;
  }
  else
  {
    local_E.pos = global_path.back();
  }

  // if using ground robot, we fix the z
  if (par_.vehicle_type != "uav")
  {
    local_E.pos[2] = 1.0;
  }

  whole_traj_solver_ptr->prepareSolverForReplan(A_time, global_path, safe_corridor_polytopes_whole_, local_trajs, local_A, local_E, initial_guess_computation_time, par_.use_multiple_initial_guesses);

  // It's pushed in prepareSolverForReplan() so we get the pushed global path
  whole_traj_solver_ptr->getGlobalPath(global_path);

  mtx_global_path_.lock();
  global_path_ = global_path; // Update the global path
  mtx_global_path_.unlock();

  // update local_E
  local_E.pos = global_path.back();
  mtx_E_.lock();
  E_ = local_E; // Update the local_E
  mtx_E_.unlock();

  if (par_.debug_verbose)
    std::cout << "Solver prepared" << std::endl;

  // Get initial guesses
  auto list_z0 = whole_traj_solver_ptr_->getInitialGuesses();
  auto list_initial_guess_wps = whole_traj_solver_ptr_->getInitialGuessWaypoints();

  if (par_.debug_verbose)
    std::cout << "Initial guesses size: " << list_z0.size() << std::endl;

  // Update L-BFGS parameters.
  // lbfgs_params_.mem_size = static_cast<int>(list_z0[0].size());
  lbfgs_params_.mem_size = 256;

  // Prepare vectors for parallelization
  int status = -1; // Initialize status
  int size_of_list_z0 = list_z0.size();

  std::vector<int> list_status(size_of_list_z0, -1); // Initialize status for each thread
  std::vector<Eigen::VectorXd> list_zopt(size_of_list_z0);
  std::vector<double> list_fopt(size_of_list_z0, 0.0);
  std::vector<double> list_initial_guess_computation_time(size_of_list_z0);

  // Solve the optimiation problem.
  if (!par_.use_multiple_initial_guesses || size_of_list_z0 == 1)
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    status = whole_traj_solver_ptr_->optimize(list_z0[0], zopt_, fopt_, lbfgs_params_);
    auto t_end = std::chrono::high_resolution_clock::now();
    local_traj_computation_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    // std::cout << lbfgs::lbfgs_strerror(status) << std::endl;

    // status = 10;
    // zopt_ = list_z0[0]; // Initialize zopt_ with the correct size
    // fopt_ = 0.0; // Initialize fopt_ to zero
    list_initial_guess_wps_subopt_.clear();
    list_initial_guess_wps_subopt_.push_back(list_initial_guess_wps[0]); // Store the initial guess waypoints for the optimal solution
  }
  else
  {
    // Parallelization approach (is there any faster way to do this?)
    std::vector<std::future<void>> futures;
    futures.reserve(size_of_list_z0);

    auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < size_of_list_z0; ++i)
    {
      futures.emplace_back(std::async(std::launch::async,
                                      [&, i]()
                                      {
                                        // make a fresh solver for thread-safety
                                        std::shared_ptr<lbfgs::SolverLBFGS> solver_ptr = std::make_shared<lbfgs::SolverLBFGS>();
                                        solver_ptr->initializeSolver(planner_params_);
                                        double initial_guess_computation_time = 0.0;
                                        solver_ptr->prepareSolverForReplan(A_time, global_path, safe_corridor_polytopes_whole_, local_trajs, local_A, local_E, initial_guess_computation_time); // initial time t0 = 0.0

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
    for (auto &f : futures)
      f.get();

    auto t_end = std::chrono::high_resolution_clock::now();
    local_traj_computation_time = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    // Find the best solution
    fopt_ = std::numeric_limits<double>::max();
    zopt_ = list_z0[0]; // Initialize zopt_ with the correct size
    int best_index = -1;
    for (size_t i = 0; i < size_of_list_z0; ++i)
    {
      // First check the status
      if (list_status[i] < 0)
      {
        std::cout << "list_status[" << i << "] = " << list_status[i] << ", skipping this solution" << std::endl;
        std::cout << lbfgs::lbfgs_strerror(list_status[i]) << std::endl;
        continue; // Skip this solution if it failed
      }

      if (list_fopt[i] < fopt_)
      {
        fopt_ = list_fopt[i];
        zopt_ = list_zopt[i];
        best_index = i;
      }
    }

    initial_guess_computation_time = list_initial_guess_computation_time[best_index]; // Get the initial guess computation time for the best solution

    // If no solution was found, return
    if (best_index < 0)
    {
      std::cout << bold << red << "No solution found in multiple initial guesses" << reset << std::endl;
      replanning_failure_count_++;
      initial_guess_computation_time = -100000.0;
      return false;
    }

    list_z_subopt_.clear();                 // Clear the suboptimal solutions list
    list_initial_guess_wps_subopt_.clear(); // Clear the suboptimal waypoints list

    // Grab sub optimal solutions for visualization (solutions that are valid (so no negative status and no crazy big f) but no better than the best one)
    for (size_t i = 0; i < size_of_list_z0; ++i)
    {
      if (list_status[i] >= 0 && i != best_index && list_fopt[i] < planner_params_.BIG * 0.9)
      {
        // Add the solution to the trajs_ for visualization
        list_z_subopt_.push_back(list_zopt[i]);

        // Add the initial guess waypoints for visualization
        list_initial_guess_wps_subopt_.push_back(list_initial_guess_wps[i]);
      }
    }

    // Set the status to the best one found
    status = list_status[best_index];

  } // End of parallelization

  if (par_.debug_verbose)
    std::cout << "Optimization status: " << status << ", fopt: " << fopt_ << ", computation time: " << local_traj_computation_time << " ms" << std::endl;

  // If no solution is found, return.
  if (status < 0 || fopt_ > par_.fopt_threshold)
  {
    // do the same output in red with printf
    printf("\033[1;31mLocal Optimization Failed with status: %d, fopt: %.2f\033[0m\n", status, fopt_);
    return false;
  }

  // If the optimization succeeded, we can update the goal setpoints in green
  // printf("\033[1;32mLocal Optimization Succeeded with status: %d, fopt: %.2f\033[0m\n", status, fopt_);
  return true;
}

// ----------------------------------------------------------------------------

bool MIGHTY::appendToPlan()
{

  if (par_.debug_verbose)
    std::cout << "goal_setpoints_.size(): " << goal_setpoints_.size() << std::endl;

  // mutex lock
  mtx_plan_.lock();

  // get the size of the plan and plan_safe_paths
  int plan_size = plan_.size();

  // If the plan size is less than k_value_, which means we already passed point A, we cannot use this plan
  if (plan_size < k_value_)
  {
    if (par_.debug_verbose)
      std::cout << bold << red << "(plan_size - k_value_) = " << (plan_size - k_value_) << " < 0" << reset << std::endl;
    k_value_ = std::max(1, plan_size - 1); // Decrease k_value_ to plan_size - 1 but at least 1
  }
  else // If the plan size is greater than k_value_, which means we haven't passed point A yet, we can use this plan
  {
    plan_.erase(plan_.end() - k_value_, plan_.end());
    plan_.insert(plan_.end(), goal_setpoints_.begin(), goal_setpoints_.end());
  }

  // mutex unlock
  mtx_plan_.unlock();

  // k_value adaptation initialization
  if (!got_enough_replanning_)
  {
    if (store_computation_times_.size() < par_.num_replanning_before_adapt)
    {
      num_replanning_++;
    }
    else
    {
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
void MIGHTY::getGterm(state &G_term)
{
  mtx_G_term_.lock();
  G_term = G_term_;
  mtx_G_term_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the terminal goal state.
 * @param state G_term: Terminal goal state to set.
 */
void MIGHTY::setGterm(const state &G_term)
{
  mtx_G_term_.lock();
  G_term_ = G_term;
  mtx_G_term_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the subgoal.
 * @param state &G: Output subgoal.
 */
void MIGHTY::getG(state &G)
{
  mtx_G_.lock();
  G = G_;
  mtx_G_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets point E
 * @param state &G: Output point E
 */
void MIGHTY::getE(state &E)
{
  mtx_E_.lock();
  E = E_;
  mtx_E_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the subgoal.
 * @param state G: Subgoal to set.
 */
void MIGHTY::setG(const state &G)
{
  mtx_G_.lock();
  G_ = G;
  mtx_G_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets A (starting point for global planning).
 * @param state &G: Output A.
 */
void MIGHTY::getA(state &A)
{
  mtx_A_.lock();
  A = A_;
  mtx_A_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets A (starting point for global planning).
 * @param state &G: Input A.
 */
void MIGHTY::setA(const state &A)
{
  mtx_A_.lock();
  A_ = A;
  mtx_A_.unlock();
}

// ----------------------------------------------------------------------------

void MIGHTY::getA_time(double &A_time)
{
  mtx_A_time_.lock();
  A_time = A_time_;
  mtx_A_time_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets A (starting point for global planning)'s time
 * @param state &G: Input A time
 */
void MIGHTY::setA_time(double A_time)
{
  mtx_A_time_.lock();
  A_time_ = A_time;
  mtx_A_time_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the current state.
 * @param state &state: Output current state.
 */
void MIGHTY::getState(state &state)
{
  mtx_state_.lock();
  state = state_;
  mtx_state_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets the last plan state
 * @param state &state: Output last plan state
 */
void MIGHTY::getLastPlanState(state &state)
{
  mtx_plan_.lock();
  state = plan_.back();
  mtx_plan_.unlock();
}

// ----------------------------------------------------------------------------

/**
 * @brief Gets trajs_
 * @param std::vector<std::shared_ptr<dynTraj>> &trajs: Output trajs_
 */
void MIGHTY::getTrajs(std::vector<std::shared_ptr<dynTraj>> &out)
{
  std::lock_guard<std::mutex> lock(mtx_trajs_);
  out = trajs_; // copies shared_ptr only, not expressions
}

// ----------------------------------------------------------------------------

/**
 * @brief Cleans up old trajectories.
 * @param double current_time: Current timestamp.
 */
void MIGHTY::cleanUpOldTrajs(double current_time)
{
  std::lock_guard<std::mutex> lock(mtx_trajs_);

  // remove_if moves all “expired” to the end, then erase() chops them off
  trajs_.erase(
      std::remove_if(
          trajs_.begin(),
          trajs_.end(),
          [&](const std::shared_ptr<dynTraj> &t)
          {
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
void MIGHTY::addTraj(std::shared_ptr<dynTraj> new_traj, double current_time)
{

  if (new_traj->mode == dynTraj::Mode::Analytic && !new_traj->analytic_compiled)
  {
    printf("Dropping analytic traj id=%d because compile failed", new_traj->id);
    return;
  }

  // Evaluate once (no copies)
  Eigen::Vector3d p = new_traj->eval(current_time);
  if (!checkPointWithinMap(p))
    return;
  if ((p - state_.pos).norm() > par_.horizon)
    return;

  {
    std::lock_guard<std::mutex> lock(mtx_trajs_);
    auto it = std::find_if(trajs_.begin(), trajs_.end(),
                           [&](const std::shared_ptr<dynTraj> &t)
                           { return t->id == new_traj->id; });

    if (it != trajs_.end())
      *it = new_traj; // replace pointer
    else
      trajs_.push_back(new_traj);
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Updates the current state.
 * @param state data: New state data.
 */
void MIGHTY::updateState(state data)
{

  // If we are doing hardware and provide goal in global frame (e.g. vicon), we need to transform the goal to the local frame

  if (par_.use_hardware && par_.provide_goal_in_global_frame)
  {
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

  mtx_state_.lock();
  state_ = data;
  mtx_state_.unlock();

  if (state_initialized_ == false || drone_status_ == DroneStatus::YAWING)
  {

    // create temporary state
    state tmp;
    tmp.pos = data.pos;
    tmp.yaw = data.yaw;
    previous_yaw_ = data.yaw;

    // Push the state to the plan
    mtx_plan_.lock();
    plan_.clear();
    plan_.push_back(tmp);
    mtx_plan_.unlock();

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
bool MIGHTY::getNextGoal(state &next_goal)
{

  // Check if the planner is initialized
  if (!checkReadyToReplan())
  {
    return false;
  }

  // Pop the front of the plan
  next_goal.setZero();

  // If the plan is empty, return false
  mtx_plan_.lock(); // Lock the mutex
  auto local_plan = plan_;
  mtx_plan_.unlock(); // Unlock the mutex

  // Get the next goal
  next_goal = local_plan.front();

  // If there's more than one goal setpoint, pop the front
  if (local_plan.size() > 1)
  {
    mtx_plan_.lock();
    plan_.pop_front();
    mtx_plan_.unlock();
  }

  if (par_.use_hardware && par_.provide_goal_in_global_frame)
  {
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

  if (!(drone_status_ == DroneStatus::GOAL_REACHED))
  {
    // Get the desired yaw
    // If the planner keeps failing, just keep spinning
    if (replanning_failure_count_ > par_.yaw_spinning_threshold)
    {
      next_goal.yaw = previous_yaw_ + par_.yaw_spinning_dyaw * par_.dc;
      next_goal.dyaw = par_.yaw_spinning_dyaw;
      previous_yaw_ = next_goal.yaw;
    }
    else
    {
      // If the local_plan is small just use the previous yaw with no dyaw
      if (local_plan.size() < 5)
      {
        next_goal.yaw = previous_yaw_;
        next_goal.dyaw = 0.0;
      }
      else
      {
        getDesiredYaw(next_goal);
      }
    }

    if (par_.use_hardware && par_.provide_goal_in_global_frame)
    {
      next_goal.yaw -= yaw_init_offset_;
    }

    next_goal.dyaw = std::clamp(next_goal.dyaw, -par_.w_max, par_.w_max);
  }
  else
  {
    next_goal.yaw = previous_yaw_;
    next_goal.dyaw = 0.0;
  }

  return true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Computes the desired yaw for the next goal.
 * @param state &next_goal: Next goal state to update with desired yaw.
 */
void MIGHTY::getDesiredYaw(state &next_goal)
{

  double diff = 0.0;
  double desired_yaw = 0.0;

  // Get state
  state local_state;
  getState(local_state);

  // Get G_term
  mtx_G_term_.lock();
  state G_term = G_term_;
  mtx_G_term_.unlock();

  switch (drone_status_)
  {
  case DroneStatus::YAWING:
    desired_yaw = atan2(G_term.pos[1] - next_goal.pos[1], G_term.pos[0] - next_goal.pos[0]);
    diff = desired_yaw - local_state.yaw;
    // std::cout << "diff1= " << diff << std::endl;
    break;
  case DroneStatus::TRAVELING:
  case DroneStatus::GOAL_SEEN:
    desired_yaw = atan2(next_goal.pos[1] - local_state.pos.y(), next_goal.pos[0] - local_state.pos.x());
    diff = desired_yaw - local_state.yaw;
    next_goal.yaw = desired_yaw;
    break;
  case DroneStatus::GOAL_REACHED:
    next_goal.dyaw = 0.0;
    next_goal.yaw = previous_yaw_;
    return;
  }

  mighty_utils::angle_wrap(diff);
  if (fabs(diff) < 0.04 && drone_status_ == DroneStatus::YAWING)
  {
    changeDroneStatus(DroneStatus::TRAVELING);
  }

  yaw(diff, next_goal);
}

// ----------------------------------------------------------------------------

void MIGHTY::yaw(double diff, state &next_goal)
{
  saturate(diff, -par_.dc * par_.w_max, par_.dc * par_.w_max);
  dyaw_filtered_ = (1 - par_.alpha_filter_dyaw) * (copysign(1, diff) * par_.w_max) + par_.alpha_filter_dyaw * dyaw_filtered_;
  next_goal.dyaw = dyaw_filtered_;
  next_goal.yaw = previous_yaw_ + dyaw_filtered_ * par_.dc;
  previous_yaw_ = next_goal.yaw;
}

// ----------------------------------------------------------------------------

/**
 * @brief Sets the terminal goal.
 * @param const state &term_goal: Desired terminal goal state.
 */
void MIGHTY::setTerminalGoal(const state &term_goal)
{

  // Get the state
  state local_state;
  getState(local_state);

  // Set the terminal goal
  setGterm(term_goal);

  // Project the terminal goal to the sphere
  mtx_G_.lock();
  G_.pos = mighty_utils::projectPointToSphere(local_state.pos, term_goal.pos, par_.horizon);
  mtx_G_.unlock();

  changeDroneStatus(DroneStatus::TRAVELING);

  if (!terminal_goal_initialized_)
    terminal_goal_initialized_ = true;
}

// ----------------------------------------------------------------------------

/**
 * @brief Changes the drone's status (YAWING, TRAVELING, GOAL_SEEN, GOAL_REACHED).
 * @param int new_status: New status value.
 */
void MIGHTY::changeDroneStatus(int new_status)
{
  if (new_status == drone_status_)
    return;

  std::cout << "Changing DroneStatus from ";

  switch (drone_status_)
  {
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

  switch (new_status)
  {
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
bool MIGHTY::checkReadyToReplan()
{

  if (!kdtree_map_initialized_ ||
      // !kdtree_unk_initialized_ ||
      !state_initialized_ ||
      !terminal_goal_initialized_ ||
      !dgp_manager_.isMapInitialized())
  {
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------

void MIGHTY::updateMap(
    const pcl::PointCloud<pcl::PointXYZ>::ConstPtr &pclptr_map,
    const pcl::PointCloud<pcl::PointXYZ>::ConstPtr &pclptr_unk)
{
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

  // 2) map update (unlocked)
  dgp_manager_.updateMap(wdx_, wdy_, wdz_, map_center_, pclptr_map_);

  // 3) Known‐space KD‐tree
  if (pclptr_map_ && !pclptr_map_->points.empty())
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_map_);
    kdtree_map_.setInputCloud(pclptr_map_);
    kdtree_map_initialized_ = true;
    dgp_manager_.updateVecOccupied(pclptr_to_vec(pclptr_map_));
  }
  else
  {
    RCLCPP_WARN(
        rclcpp::get_logger("mighty"),
        "updateMap: member pclptr_map_ was null or empty; skipping KD‐tree update");
  }

  // 4) Unknown‐space KD‐tree
  if (pclptr_unk_ && !pclptr_unk_->points.empty())
  {
    std::lock_guard<std::mutex> lk(mtx_kdtree_unk_);
    kdtree_unk_.setInputCloud(pclptr_unk_);
    kdtree_unk_initialized_ = true;
    // merge known into unknown vector
    dgp_manager_.updateVecUnknownOccupied(pclptr_to_vec(pclptr_unk_));
    dgp_manager_.insertVecOccupiedToVecUnknownOccupied();
  }
  else
  {
    RCLCPP_WARN(
        rclcpp::get_logger("mighty"),
        "updateMap: member pclptr_unk_ was null or empty; skipping KD‐tree update");
  }
}

// ----------------------------------------------------------------------------

/**
 * @brief Set the initial pose.
 * @param const geometry_msgs::msg::TransformStamped &init_pose: Initial pose.
 */
void MIGHTY::setInitialPose(const geometry_msgs::msg::TransformStamped &init_pose)
{
  init_pose_ = init_pose;

  // First compute transformation matrix from init_pose_ (geometry_msgs::msg::TransformStamped)
  Eigen::Matrix4d init_pose_transform = Eigen::Matrix4d::Identity();
  Eigen::Quaterniond init_pose_quat(init_pose_.transform.rotation.w, init_pose_.transform.rotation.x, init_pose_.transform.rotation.y, init_pose_.transform.rotation.z);
  Eigen::Vector3d init_pose_translation(init_pose_.transform.translation.x, init_pose_.transform.translation.y, init_pose_.transform.translation.z);
  init_pose_transform.block<3, 3>(0, 0) = init_pose_quat.toRotationMatrix();
  init_pose_transform.block<3, 1>(0, 3) = init_pose_translation;

  // Get initial pose
  init_pose_transform_ = init_pose_transform;
  init_pose_transform_rotation_ = init_pose_quat.toRotationMatrix();
  yaw_init_offset_ = std::atan2(init_pose_transform_rotation_(1, 0),
                                init_pose_transform_rotation_(0, 0));

  std::cout << bold << green << "yaw_init_offset_: " << yaw_init_offset_ << reset << std::endl;

  // Get the inverse of init_pose_ (geometry_msgs::msg::TransformStamped)
  init_pose_transform_inv_ = init_pose_transform.inverse();
  init_pose_transform_rotation_inv_ = init_pose_quat.toRotationMatrix().inverse();
  // yaw_init_offset_ = std::atan2(init_pose_transform_rotation_inv_(1, 0),
  // init_pose_transform_rotation_inv_(0, 0));
}

// ----------------------------------------------------------------------------

// Apply the initial pose transformation to the pwp
void MIGHTY::applyInitiPoseTransform(PieceWisePol &pwp)
{
  // Loop thru the intervals
  for (int i = 0; i < pwp.coeff_x.size(); i++)
  {
    // Loop thru a, b, c, and d
    for (int j = 0; j < 4; j++)
    {
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
void MIGHTY::applyInitiPoseInverseTransform(PieceWisePol &pwp)
{
  // Loop thru the intervals
  for (int i = 0; i < pwp.coeff_x.size(); i++)
  {

    // Loop thru a, b, c, and d
    for (int j = 0; j < 4; j++)
    {
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
bool MIGHTY::goalReachedCheck()
{
  if (checkReadyToReplan() && drone_status_ == DroneStatus::GOAL_REACHED)
  {
    return true;
  }
  return false;
}