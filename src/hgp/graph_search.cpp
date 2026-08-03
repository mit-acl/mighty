/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "hgp/graph_search.hpp"

#include <cmath>

using namespace mighty;
using namespace termcolor;

typedef timer::Timer MyTimer;

// Help functions
static inline double clamp01(double x) { return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x); }
static inline double hypot2d(double x, double y) { return std::sqrt(x * x + y * y); }

GraphSearch::GraphSearch(const int8_t* cMap, const std::shared_ptr<mighty::VoxelMapUtil>& map_util,
                         int xDim, int yDim, int zDim, double eps, bool verbose,
                         std::string global_planner, double w_unknown, double w_align,
                         double decay_len_cells, double w_side)
    : cMap_(cMap),
      map_util_(map_util),
      xDim_(xDim),
      yDim_(yDim),
      zDim_(zDim),
      eps_(eps),
      verbose_(verbose),
      global_planner_(global_planner),
      w_unknown_(w_unknown),
      w_align_(w_align),
      decay_len_cells_(decay_len_cells),
      w_side_(w_side) {
  hm_.assign(xDim_ * yDim_ * zDim_, nullptr);
  seen_.assign(xDim_ * yDim_ * zDim_, false);

  // Set neighbors: 2D (8-connected) when zDim==1, else 3D (26-connected)
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      const int z_lo = (zDim_ == 1) ? 0 : -1;
      const int z_hi = (zDim_ == 1) ? 0 : 1;
      for (int z = z_lo; z <= z_hi; z++) {
        if (x == 0 && y == 0 && z == 0) continue;
        ns_.push_back(std::vector<int>{x, y, z});
      }
    }
  }

  jn3d_ = std::make_shared<JPS3DNeib>();
}

inline int GraphSearch::coordToId(int x, int y, int z) const {
  return x + y * xDim_ + z * xDim_ * yDim_;
}

inline bool GraphSearch::isFree(int x, int y, int z) const {
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ && z >= 0 && z < zDim_ &&
         cMap_[coordToId(x, y, z)] <= val_free_;
}

inline bool GraphSearch::isOccupied(int x, int y, int z) const {
  return x < 0 || x >= xDim_ || y < 0 || y >= yDim_ || z < 0 || z >= zDim_ ||
         cMap_[coordToId(x, y, z)] == val_occupied_;
}

inline bool GraphSearch::isUnknown(int x, int y, int z) const {
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ && z >= 0 && z < zDim_ &&
         cMap_[coordToId(x, y, z)] == val_unknown_;
}

inline double GraphSearch::getHeur(int x, int y, int z) const {
  return eps_ * std::sqrt((x - xGoal_) * (x - xGoal_) + (y - yGoal_) * (y - yGoal_) +
                          (z - zGoal_) * (z - zGoal_));
}

bool GraphSearch::plan(int xStart, int yStart, int zStart, int xGoal, int yGoal, int zGoal,
                       double initial_g, double& global_planning_time, double& hgp_static_jps_time,
                       double& hgp_check_path_time, double& hgp_dynamic_astar_time,
                       double& hgp_recover_path_time, double current_time, const Vec3f& start_vel,
                       int max_expand, int timeout_duration_ms) {
  pq_.clear();
  path_.clear();
  hm_.assign(xDim_ * yDim_ * zDim_, nullptr);
  seen_.assign(xDim_ * yDim_ * zDim_, false);

  start_vel_ = start_vel;

  // Set goal
  int goal_id = coordToId(xGoal, yGoal, zGoal);
  xGoal_ = xGoal;
  yGoal_ = yGoal;
  zGoal_ = zGoal;
  // Set start node

  int start_id = coordToId(xStart, yStart, zStart);
  StatePtr currNode_ptr =
      std::make_shared<State>(State(start_id, xStart, yStart, zStart, 0, 0, 0, 0));

  currNode_ptr->g = initial_g;
  currNode_ptr->h = getHeur(xStart, yStart, zStart);

  // Set current time
  current_time_ = current_time;

  // Start a timer
  MyTimer global_planning_timer(true);

  // Run the planner
  bool result = select_planner(currNode_ptr, max_expand, start_id, goal_id,
                               std::chrono::milliseconds(timeout_duration_ms));

  // Print the time taken [ms]
  global_planning_time_ = global_planning_timer.getElapsedMicros() / 1000.0;

  // Assign the time taken
  global_planning_time = global_planning_time_;
  hgp_static_jps_time = hgp_static_jps_time_;
  hgp_check_path_time = hgp_check_path_time_;
  hgp_dynamic_astar_time = hgp_dynamic_astar_time_;
  hgp_recover_path_time = hgp_recover_path_time_;

  return result;
}

bool GraphSearch::select_planner(StatePtr& currNode_ptr, int max_expand, int start_id, int goal_id,
                                 std::chrono::milliseconds timeout_duration) {
  // Make planner choice unambiguous:
  //  - sjps         : JPS successors (getJpsSucc)
  //  - astar        : A* grid successors (getSucc) with NO heat
  //  - astar_heat   : A* grid successors (getSucc) WITH heat
  use_jps_ = false;
  use_heat_ = false;

  if (global_planner_ == "sjps") {
    use_jps_ = true;
    return static_jps_plan(currNode_ptr, max_expand, start_id, goal_id, timeout_duration);
  } else if (global_planner_ == "sastar") {
    return static_jps_plan(currNode_ptr, max_expand, start_id, goal_id, timeout_duration);
  } else if (global_planner_ == "astar_heat") {
    use_heat_ = true;
    return static_jps_plan(currNode_ptr, max_expand, start_id, goal_id, timeout_duration);
  } else {
    printf("Unknown planner: %s\n", global_planner_.c_str());
    return false;
  }
}

bool GraphSearch::static_jps_plan(StatePtr& currNode_ptr, int max_expand, int start_id, int goal_id,
                                  std::chrono::milliseconds timeout_duration) {
  // Record the start time
  auto start_time = std::chrono::steady_clock::now();

  if (!currNode_ptr) {
    std::cerr << "Error: currNode_ptr is null!" << std::endl;
    return false;
  }

  // Insert start node
  currNode_ptr->heapkey = pq_.push(currNode_ptr);
  currNode_ptr->opened = true;
  hm_[currNode_ptr->id] = currNode_ptr;
  seen_[currNode_ptr->id] = true;

  int expand_iteration = 0;
  cMap_ = (map_util_->map_).data();

  // Track the best (closest-to-goal) node for partial path recovery
  StatePtr best_node = currNode_ptr;
  double best_h = currNode_ptr->h;

  if (verbose_) {
    std::cout << "[GraphSearch] planner=\"" << global_planner_
              << "\" use_jps=" << (use_jps_ ? 1 : 0) << " use_heat=" << (use_heat_ ? 1 : 0)
              << std::endl;
  }

  while (true) {
    expand_iteration++;

    // Check if timeout has occurred
    auto current_time = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time) >
        timeout_duration) {
      // std::cerr << "astar_heat: timeout after " << expand_iteration
      //           << " expansions, recovering partial path\n";
      path_ = recoverPath(best_node, start_id);
      return !path_.empty() && path_.size() > 1;
    }

    if (pq_.empty()) {
      // std::cerr << "astar_heat: priority queue empty after " << expand_iteration
      //           << " expansions, recovering partial path\n";
      path_ = recoverPath(best_node, start_id);
      return !path_.empty() && path_.size() > 1;
    }

    // get element with smallest cost
    currNode_ptr = pq_.top();
    pq_.pop();
    currNode_ptr->closed = true;  // Add to closed list

    // Update best node (closest to goal by heuristic)
    if (currNode_ptr->h < best_h) {
      best_h = currNode_ptr->h;
      best_node = currNode_ptr;
    }

    if (currNode_ptr->id == goal_id) {
      if (verbose_) printf("Goal Reached!\n");
      break;
    }

    std::vector<int> succ_ids;
    std::vector<double> succ_costs;

    if (use_jps_) {
      getJpsSucc(currNode_ptr, succ_ids, succ_costs);
    } else {
      // Heat is applied ONLY when global_planner_=="astar_heat" (use_heat_==true).
      getSucc(currNode_ptr, succ_ids, succ_costs, use_heat_);
    }

    // Process successors
    for (int s = 0; s < (int)succ_ids.size(); s++) {
      // see if we can improve the value of succstate
      StatePtr& child_ptr = atHM(succ_ids[s]);
      if (!child_ptr) {
        // (Optional) reconstruct or skip
        fprintf(stderr, "[WARN] child_ptr null at id=%d\n", succ_ids[s]);
        continue;
      }

      double tentative_gval = currNode_ptr->g + succ_costs[s];

      if (tentative_gval < child_ptr->g) {
        child_ptr->parentId = currNode_ptr->id;  // Assign new parent
        child_ptr->g = tentative_gval;           // Update gval

        // if currently in OPEN, update
        if (child_ptr->opened && !child_ptr->closed) {
          pq_.increase(child_ptr->heapkey);  // update heap
          child_ptr->dx = (child_ptr->x - currNode_ptr->x);
          child_ptr->dy = (child_ptr->y - currNode_ptr->y);
          child_ptr->dz = (child_ptr->z - currNode_ptr->z);
          if (child_ptr->dx != 0) child_ptr->dx /= std::abs(child_ptr->dx);
          if (child_ptr->dy != 0) child_ptr->dy /= std::abs(child_ptr->dy);
          if (child_ptr->dz != 0) child_ptr->dz /= std::abs(child_ptr->dz);
        }
        // if currently in CLOSED
        else if (child_ptr->opened && child_ptr->closed) {
        } else  // new node, add to heap
        {
          child_ptr->heapkey = pq_.push(child_ptr);
          child_ptr->opened = true;
        }
      }  //
    }    // Process successors

    if (max_expand > 0 && expand_iteration >= max_expand) {
      // std::cerr << "astar_heat: max_expand [" << max_expand
      //           << "] reached, recovering partial path\n";
      path_ = recoverPath(best_node, start_id);
      return !path_.empty() && path_.size() > 1;
    }
  }

  path_ = recoverPath(currNode_ptr, start_id);

  return true;
}

std::vector<StatePtr> GraphSearch::removeLinePts(const std::vector<StatePtr>& path) {
  if (path.size() < 3) return path;

  std::vector<StatePtr> new_path;
  new_path.push_back(path.front());
  for (unsigned int i = 1; i < path.size() - 1; i++) {
    Vecf<3> path_i_p_1 = Vecf<3>(path[i + 1]->x, path[i + 1]->y, path[i + 1]->z);
    Vecf<3> path_i = Vecf<3>(path[i]->x, path[i]->y, path[i]->z);
    Vecf<3> path_i_m_1 = Vecf<3>(path[i - 1]->x, path[i - 1]->y, path[i - 1]->z);
    Vecf<3> p = (path_i_p_1 - path_i) - (path_i - path_i_m_1);
    if (fabs(p(0)) + fabs(p(1)) + fabs(p(2)) > 1e-2) {
      new_path.push_back(path[i]);
    }
  }
  new_path.push_back(path.back());
  return new_path;
}

std::vector<StatePtr> GraphSearch::removeCornerPts(const std::vector<StatePtr>& path) {
  if (path.size() < 2) return path;

  // cut zigzag segment
  std::vector<StatePtr> optimized_path;
  Vecf<3> pose1 = Vecf<3>(path[0]->x, path[0]->y, path[0]->z);
  Vecf<3> pose2 = Vecf<3>(path[1]->x, path[1]->y, path[1]->z);
  Vecf<3> prev_pose = pose1;
  optimized_path.push_back(path[0]);
  decimal_t cost1, cost2, cost3;

  if (!map_util_->isBlocked(pose1, pose2))
    cost1 = (pose1 - pose2).norm();
  else
    cost1 = std::numeric_limits<decimal_t>::infinity();

  for (unsigned int i = 1; i < path.size() - 1; i++) {
    pose1 = Vecf<3>(path[i]->x, path[i]->y, path[i]->z);
    pose2 = Vecf<3>(path[i + 1]->x, path[i + 1]->y, path[i + 1]->z);
    if (!map_util_->isBlocked(pose1, pose2))
      cost2 = (pose1 - pose2).norm();
    else
      cost2 = std::numeric_limits<decimal_t>::infinity();

    if (!map_util_->isBlocked(prev_pose, pose2))
      cost3 = (prev_pose - pose2).norm();
    else
      cost3 = std::numeric_limits<decimal_t>::infinity();

    if (cost3 < cost1 + cost2)
      cost1 = cost3;
    else {
      optimized_path.push_back(path[i]);
      cost1 = (pose1 - pose2).norm();
      prev_pose = pose1;
    }
  }

  optimized_path.push_back(path.back());
  return optimized_path;
}

//// RIGHT NOW THE JPS GIVES BACK A VERY SPARSE PATH, SO TO CHECK COLLISION AGAINST DYNAMIC
/// OBSTACLES, WE NEED TO INTERPOLATE THE PATH
void GraphSearch::updateGValues() {
  // first reverse the path
  std::reverse(path_.begin(), path_.end());

  // first create a path with float values
  std::vector<Vecf<3>> path_float;
  for (int i = 0; i < path_.size(); i++) {
    path_float.push_back(map_util_->intToFloat(Veci<3>(path_[i]->x, path_[i]->y, path_[i]->z)));
  }

  // convert the path back to int and update g values
  std::vector<StatePtr> new_path_int;
  for (int i = 0; i < path_float.size(); i++) {
    Veci<3> temp =
        map_util_->floatToInt(Vecf<3>(path_float[i](0), path_float[i](1), path_float[i](2)));
    StatePtr temp_ptr = std::make_shared<State>(
        State(coordToId(temp(0), temp(1), temp(2)), temp(0), temp(1), temp(2), 0, 0,
              0));  // for now dx, dy, dz are 0

    // compute total path length and update g
    if (i > 0) {
      Vecf<3> p1 = path_float[i - 1];
      Vecf<3> p2 = path_float[i];
      temp_ptr->g = new_path_int[i - 1]->g + (p2 - p1).norm();
    } else {
      temp_ptr->g = path_[0]->g;
    }

    new_path_int.push_back(temp_ptr);
  }

  // update the path
  path_ = new_path_int;

  // reverse the path
  std::reverse(path_.begin(), path_.end());
}

/// Compute the node time
double GraphSearch::computeNodeTime(const Eigen::Vector3d& pf, const Eigen::Vector3d& p0,
                                    const Eigen::Vector3d& v0, StatePtr& currNode_ptr) {
  // Initialize final velocity vector (for each axis)
  Eigen::Vector3d vf = Eigen::Vector3d::Zero();

  // Absolute displacement between nodes
  Eigen::Vector3d dist = (pf - p0).cwiseAbs();

  // Will store travel time for each dimension; the overall node time is the max.
  std::vector<double> node_times;

  // Small tolerance for checking if v0 is nearly equal to v_max.
  double epsilon = 1e-6;

  // Iterate over the dimensions (x, y, z)
  for (int i = 0; i < 3; i++) {
    double a_max = a_max_3d_(i);   // maximum acceleration in this dimension
    double v_max = v_max_3d_(i);   // maximum velocity in this dimension (assumed defined)
    double delta = pf(i) - p0(i);  // signed displacement along axis i
    double node_time = 0;          // time to traverse this dimension

    // If delta is very small, we can assume no motion is needed.
    if (std::abs(delta) < epsilon) {
      node_time = 0;
      vf(i) = v0(i);
    } else if (delta > 0)  // Desired displacement is positive.
    {
      if (v0(i) >= 0)  // Already moving in the positive direction
      {
        // If already cruising (and no deceleration is needed)
        if (v0(i) >= (v_max - epsilon)) {
          vf(i) = v_max;
          node_time = dist(i) / v_max;
        } else {
          double v_candidate = std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i));
          if (v_candidate <= v_max) {
            vf(i) = v_candidate;
            node_time = (vf(i) - v0(i)) / a_max;
          } else {
            double d_accel = (v_max * v_max - v0(i) * v0(i)) / (2 * a_max);
            double t_accel = (v_max - v0(i)) / a_max;
            if (dist(i) > d_accel) {
              double d_cruise = dist(i) - d_accel;
              double t_cruise = d_cruise / v_max;
              node_time = t_accel + t_cruise;
            } else {
              // Fallback if numerical issues arise.
              node_time = (std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i)) - v0(i)) / a_max;
            }
            vf(i) = v_max;
          }
        }
      } else  // v0(i) < 0: initially moving in the wrong (or zero) direction.
      {
        double d_stop = (v0(i) * v0(i)) / (2 * a_max);
        double t_stop = (-v0(i)) / a_max;
        double d_accel = dist(i) + d_stop;
        double v_candidate = std::sqrt(2 * a_max * d_accel);
        if (v_candidate <= v_max) {
          double t_accel = std::sqrt(2 * d_accel / a_max);
          node_time = t_stop + t_accel;
          vf(i) = v_candidate;
        } else {
          double d_to_vmax = (v_max * v_max) / (2 * a_max);
          double t_accel = v_max / a_max;
          if (d_accel > d_to_vmax) {
            double d_cruise = d_accel - d_to_vmax;
            double t_cruise = d_cruise / v_max;
            node_time = t_stop + t_accel + t_cruise;
          } else {
            node_time = t_stop + t_accel;
          }
          vf(i) = v_max;
        }
      }
    } else  // delta <= 0: Desired displacement is negative.
    {
      if (v0(i) <= 0)  // Already moving in the negative direction.
      {
        if (std::fabs(v0(i)) >= (v_max - epsilon)) {
          vf(i) = -v_max;
          node_time = dist(i) / v_max;
        } else {
          double v_candidate = std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i));
          v_candidate = -v_candidate;  // assign the negative sign
          if (std::fabs(v_candidate) <= v_max) {
            vf(i) = v_candidate;
            node_time = (std::fabs(vf(i)) - std::fabs(v0(i))) / a_max;
          } else {
            double d_accel = ((v_max * v_max) - (v0(i) * v0(i))) / (2 * a_max);
            double t_accel = (v_max - std::fabs(v0(i))) / a_max;
            if (dist(i) > d_accel) {
              double d_cruise = dist(i) - d_accel;
              double t_cruise = d_cruise / v_max;
              node_time = t_accel + t_cruise;
            } else {
              node_time =
                  (std::fabs(std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i))) - std::fabs(v0(i))) /
                  a_max;
            }
            vf(i) = -v_max;
          }
        }
      } else  // v0(i) > 0: initially moving in the positive direction while desired displacement is
              // negative.
      {
        double d_stop = (v0(i) * v0(i)) / (2 * a_max);
        double t_stop = v0(i) / a_max;
        double d_accel = dist(i) + d_stop;
        double v_candidate = std::sqrt(2 * a_max * d_accel);
        v_candidate = -v_candidate;  // final velocity is negative.
        if (std::fabs(v_candidate) <= v_max) {
          double t_accel = std::sqrt(2 * d_accel / a_max);
          node_time = t_stop + t_accel;
          vf(i) = v_candidate;
        } else {
          double d_to_vmax = (v_max * v_max) / (2 * a_max);
          double t_accel = v_max / a_max;
          if (d_accel > d_to_vmax) {
            double d_cruise = d_accel - d_to_vmax;
            double t_cruise = d_cruise / v_max;
            node_time = t_stop + t_accel + t_cruise;
          } else {
            node_time = t_stop + t_accel;
          }
          vf(i) = -v_max;
        }
      }
    }

    // Save the computed time for this axis.
    node_times.push_back(node_time);
  }  // end for loop over dimensions

  // Update the current node's velocity with the computed final velocity vector.
  currNode_ptr->vel = vf;

  // Return the overall node time (the maximum among x, y, and z dimensions)
  double overall_time = *std::max_element(node_times.begin(), node_times.end());
  return overall_time;
}

void GraphSearch::setStartAndGoal(const Vecf<3>& start, const Vecf<3>& goal) {
  start_ = start;
  goal_ = goal;
}

void GraphSearch::setEsdfGrid(std::shared_ptr<const EsdfGrid2D> grid, double weight, double d_safe) {
  esdf_grid_ = grid;
  esdf_weight_astar_ = weight;
  esdf_d_safe_astar_ = d_safe;
}

void GraphSearch::setBounds(double max_values[3]) {
  v_max_ = max_values[0];
  v_max_3d_ << max_values[0], max_values[0], max_values[0];
  a_max_ = max_values[1];
  a_max_3d_ << max_values[1], max_values[1], max_values[1];
  j_max_ = max_values[2];
  j_max_3d_ << max_values[2], max_values[2], max_values[2];
}

std::vector<StatePtr> GraphSearch::recoverPath(StatePtr node, int start_id) {
  std::vector<StatePtr> path;
  path.push_back(node);
  while (node && node->id != start_id) {
    node = atHM(node->parentId);
    path.push_back(node);
  }

  return path;
}

void GraphSearch::getSucc(const StatePtr& curr, std::vector<int>& succ_ids,
                          std::vector<double>& succ_costs, bool use_heat) {
  succ_ids.clear();
  succ_costs.clear();

  // --- Build a forward reference once per node ---
  // Prefer start_vel_ when available; if it's anti-aligned with goal, flip it.
  Eigen::Vector2d vref(start_vel_(0), start_vel_(1));
  if (vref.squaredNorm() < 1e-9) {
    vref = Eigen::Vector2d(goal_(0) - start_(0), goal_(1) - start_(1));
  } else {
    Eigen::Vector2d vg(goal_(0) - start_(0), goal_(1) - start_(1));
    if (vg.squaredNorm() > 1e-9 && vref.dot(vg) < 0.0) vref = -vref;  // never bias away from goal
  }
  const double vref_n = vref.norm();

  for (const auto& d : ns_) {
    int new_x = curr->x + d[0];
    int new_y = curr->y + d[1];
    int new_z = curr->z + d[2];

    // Bounds check FIRST. Without this, when use_soft_cost_obstacles is on,
    // isOccupied() returns true for out-of-grid voxels but the soft-cost
    // branch lets the successor through. coordToId then aliases negative
    // coordinates to a valid id (e.g. (-1, 58, 9) → id 125694 = (114, 57, 9))
    // because the linear index x + y*xDim + z*xDim*yDim wraps around when x
    // is negative. The wraparound corrupts hm_ and produces backward paths
    // through the negative-x boundary. JPS doesn't hit this because jump()
    // uses isFree() which is bounds-correct.
    if (new_x < 0 || new_x >= xDim_ ||
        new_y < 0 || new_y >= yDim_ ||
        new_z < 0 || new_z >= zDim_) {
      continue;
    }

    // Occupancy check
    if (map_util_ && map_util_->has2DMap() && zDim_ == 1 && esdf_grid_) {
      // ESDF mode: only check the 2D ESDF-derived map (skip inflated 3D grid)
      if (map_util_->get2DOccupancy(new_x, new_y) != 0) {
        continue;
      }
    } else {
      // Standard mode: check 3D grid + 2D projection
      if (isOccupied(new_x, new_y, new_z)) {
        if (!map_util_ || !map_util_->useSoftCostObstacles()) {
          continue;
        }
      }
      if (map_util_ && map_util_->has2DMap() && map_util_->get2DOccupancy(new_x, new_y) != 0) {
        continue;
      }
    }

    // Hard heat cutoff: treat cells with heat > cutoff_ratio * Hmax as impassable
    if (use_heat && map_util_ && map_util_->staticHeatEnabled() && zDim_ == 1 &&
        map_util_->heat_cutoff_ratio_ > 0.0f) {
      float h = map_util_->getHeat2D(new_x, new_y);
      if (h > map_util_->heat_cutoff_ratio_ * map_util_->static_heat_Hmax_) {
        continue;
      }
    }

    int new_id = coordToId(new_x, new_y, new_z);
    if (new_id < 0 || new_id >= (int)hm_.size()) continue;

    if (!seen_[new_id]) {
      seen_[new_id] = true;
      hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, new_z, d[0], d[1], d[2]);
      if (new_id < 0 || new_id >= (int)hm_.size()) {
        fprintf(stderr, "[BUG] Creation failed id=%d\n", new_id);
        std::abort();
      }
      hm_[new_id]->h = getHeur(new_x, new_y, new_z);
    }

    // Base geometric step cost (1, √2, √3)
    double step_cost = std::sqrt(double(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));

    // -------- Dynamic heat-map cost (soft) --------
    // Static obstacles remain hard-blocked by isOccupied() above.
    // Heat is time-invariant and precomputed in map_util_ during readMap().
    if (use_heat && map_util_ &&
        (map_util_->dynamicHeatEnabled() || map_util_->staticHeatEnabled())) {
      const float w_heat = map_util_->getHeatWeight();
      if (w_heat > 0.0f) {
        // Use 2D heat for ground robots (zDim==1), 3D heat otherwise
        const float h = (zDim_ == 1) ? map_util_->getHeat2D(new_x, new_y)
                                     : map_util_->getHeat(new_x, new_y, new_z);
        step_cost += (double)(w_heat * h);
      }
    }

    // -------- Add soft obstacle cost if enabled --------
    if (use_heat && map_util_ && map_util_->useSoftCostObstacles()) {
      if (isOccupied(new_x, new_y, new_z))  // Check if this cell is occupied
      {
        const float w_heat = map_util_->getHeatWeight();
        const float soft_cost = map_util_->getObstacleSoftCost();
        step_cost += (double)(w_heat * soft_cost);
      }
    }

    // -------- ESDF distance cost (ground robot 2D only) --------
    // Disabled: ESDF avoidance is now handled via heat_2d_ built from buildMap2DFromEsdf().
    // Keeping this block would double-count obstacle proximity cost.
    if (false && esdf_grid_ && zDim_ == 1 && esdf_weight_astar_ > 0.0) {
      // Convert grid coords to world coords for ESDF query
      double wx, wy;
      if (map_util_) {
        Vecf<3> world_pos = map_util_->intToFloat(Veci<3>(new_x, new_y, 0));
        wx = world_pos(0);
        wy = world_pos(1);
      } else {
        wx = new_x;
        wy = new_y;
      }
      if (esdf_grid_->isInBounds(wx, wy)) {
        double d = esdf_grid_->queryDistance(wx, wy);
        if (d < esdf_d_safe_astar_) {
          double pen = esdf_d_safe_astar_ - d;
          step_cost += esdf_weight_astar_ * pen * pen;
        }
      }
    }

    // -------- Unknown cell penalty --------
    if (w_unknown_ > 0.0 && isUnknown(new_x, new_y, new_z)) {
      step_cost += w_unknown_;
    }

    succ_ids.push_back(new_id);
    succ_costs.push_back(step_cost);
  }
}

void GraphSearch::getJpsSucc(const StatePtr& curr, std::vector<int>& succ_ids,
                             std::vector<double>& succ_costs) {
  const int norm1 = std::abs(curr->dx) + std::abs(curr->dy) + std::abs(curr->dz);
  int num_neib = jn3d_->nsz[norm1][0];
  int num_fneib = jn3d_->nsz[norm1][1];
  int id = (curr->dx + 1) + 3 * (curr->dy + 1) + 9 * (curr->dz + 1);

  for (int dev = 0; dev < num_neib + num_fneib; ++dev) {
    int new_x, new_y, new_z;
    int dx, dy, dz;
    if (dev < num_neib) {
      dx = jn3d_->ns[id][0][dev];
      dy = jn3d_->ns[id][1][dev];
      dz = jn3d_->ns[id][2][dev];
      if (!jump(curr->x, curr->y, curr->z, dx, dy, dz, new_x, new_y, new_z)) continue;
    } else {
      int nx = curr->x + jn3d_->f1[id][0][dev - num_neib];
      int ny = curr->y + jn3d_->f1[id][1][dev - num_neib];
      int nz = curr->z + jn3d_->f1[id][2][dev - num_neib];
      if (isOccupied(nx, ny, nz)) {
        dx = jn3d_->f2[id][0][dev - num_neib];
        dy = jn3d_->f2[id][1][dev - num_neib];
        dz = jn3d_->f2[id][2][dev - num_neib];
        if (!jump(curr->x, curr->y, curr->z, dx, dy, dz, new_x, new_y, new_z)) continue;
      } else
        continue;
    }

    int new_id = coordToId(new_x, new_y, new_z);
    if (!seen_[new_id]) {
      seen_[new_id] = true;
      hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, new_z, dx, dy, dz);
      if (new_id < 0 || new_id >= (int)hm_.size()) {
        fprintf(stderr, "[BUG] Creation failed id=%d\n", new_id);
        std::abort();
      }
      hm_[new_id]->h = getHeur(new_x, new_y, new_z);
    }

    succ_ids.push_back(new_id);
    succ_costs.push_back(std::sqrt((new_x - curr->x) * (new_x - curr->x) +
                                   (new_y - curr->y) * (new_y - curr->y) +
                                   (new_z - curr->z) * (new_z - curr->z)));
  }
}

bool GraphSearch::jump(int x, int y, int z, int dx, int dy, int dz, int& new_x, int& new_y,
                       int& new_z) {
  new_x = x + dx;
  new_y = y + dy;
  new_z = z + dz;
  if (!isFree(new_x, new_y, new_z)) return false;

  // Prevent corner-cutting in 2D mode: for diagonal moves, check axis-aligned neighbors
  if (zDim_ == 1) {
    const int move_norm = std::abs(dx) + std::abs(dy);
    if (move_norm >= 2 && dx != 0 && dy != 0) {
      if (!isFree(x + dx, y, z) || !isFree(x, y + dy, z)) return false;
    }
  }

  if (new_x == xGoal_ && new_y == yGoal_ && new_z == zGoal_) return true;

  if (hasForced(new_x, new_y, new_z, dx, dy, dz)) return true;

  const int id = (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1);
  const int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
  int num_neib = jn3d_->nsz[norm1][0];
  for (int k = 0; k < num_neib - 1; ++k) {
    int new_new_x, new_new_y, new_new_z;
    if (jump(new_x, new_y, new_z, jn3d_->ns[id][0][k], jn3d_->ns[id][1][k], jn3d_->ns[id][2][k],
             new_new_x, new_new_y, new_new_z))
      return true;
  }

  return jump(new_x, new_y, new_z, dx, dy, dz, new_x, new_y, new_z);
}

inline bool GraphSearch::hasForced(int x, int y, int z, int dx, int dy, int dz) {
  int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
  int id = (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1);
  switch (norm1) {
    case 1:
      // 1-d move, check 8 neighbors
      for (int fn = 0; fn < 8; ++fn) {
        int nx = x + jn3d_->f1[id][0][fn];
        int ny = y + jn3d_->f1[id][1][fn];
        int nz = z + jn3d_->f1[id][2][fn];
        if (isOccupied(nx, ny, nz)) return true;
      }
      return false;
    case 2:
      // 2-d move, check 8 neighbors
      for (int fn = 0; fn < 8; ++fn) {
        int nx = x + jn3d_->f1[id][0][fn];
        int ny = y + jn3d_->f1[id][1][fn];
        int nz = z + jn3d_->f1[id][2][fn];
        if (isOccupied(nx, ny, nz)) return true;
      }
      return false;
    case 3:
      // 3-d move, check 6 neighbors
      for (int fn = 0; fn < 6; ++fn) {
        int nx = x + jn3d_->f1[id][0][fn];
        int ny = y + jn3d_->f1[id][1][fn];
        int nz = z + jn3d_->f1[id][2][fn];
        if (isOccupied(nx, ny, nz)) return true;
      }
      return false;
    default:
      return false;
  }
}

std::vector<StatePtr> GraphSearch::getPath() const { return path_; }

std::vector<StatePtr> GraphSearch::getOpenSet() const {
  std::vector<StatePtr> ss;
  for (const auto& it : hm_) {
    if (it && it->opened && !it->closed) ss.push_back(it);
  }
  return ss;
}

std::vector<StatePtr> GraphSearch::getCloseSet() const {
  std::vector<StatePtr> ss;
  for (const auto& it : hm_) {
    if (it && it->closed) ss.push_back(it);
  }
  return ss;
}

std::vector<StatePtr> GraphSearch::getAllSet() const {
  std::vector<StatePtr> ss;
  for (const auto& it : hm_) {
    if (it) ss.push_back(it);
  }
  return ss;
}

constexpr int JPS2DNeib::nsz[3][2];

JPS2DNeib::JPS2DNeib() {
  int id = 0;
  for (int dy = -1; dy <= 1; ++dy) {
    for (int dx = -1; dx <= 1; ++dx) {
      int norm1 = std::abs(dx) + std::abs(dy);
      for (int dev = 0; dev < nsz[norm1][0]; ++dev)
        Neib(dx, dy, norm1, dev, ns[id][0][dev], ns[id][1][dev]);
      for (int dev = 0; dev < nsz[norm1][1]; ++dev) {
        FNeib(dx, dy, norm1, dev, f1[id][0][dev], f1[id][1][dev], f2[id][0][dev], f2[id][1][dev]);
      }
      id++;
    }
  }
}

void JPS2DNeib::print() {
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      int id = (dx + 1) + 3 * (dy + 1);
      printf("[dx: %d, dy: %d]-->id: %d:\n", dx, dy, id);
      for (unsigned int i = 0; i < sizeof(f1[id][0]) / sizeof(f1[id][0][0]); i++)
        printf("                f1: [%d, %d]\n", f1[id][0][i], f1[id][1][i]);
    }
  }
}

void JPS2DNeib::Neib(int dx, int dy, int norm1, int dev, int& tx, int& ty) {
  switch (norm1) {
    case 0:
      switch (dev) {
        case 0:
          tx = 1;
          ty = 0;
          return;
        case 1:
          tx = -1;
          ty = 0;
          return;
        case 2:
          tx = 0;
          ty = 1;
          return;
        case 3:
          tx = 1;
          ty = 1;
          return;
        case 4:
          tx = -1;
          ty = 1;
          return;
        case 5:
          tx = 0;
          ty = -1;
          return;
        case 6:
          tx = 1;
          ty = -1;
          return;
        case 7:
          tx = -1;
          ty = -1;
          return;
      }
    case 1:
      tx = dx;
      ty = dy;
      return;
    case 2:
      switch (dev) {
        case 0:
          tx = dx;
          ty = 0;
          return;
        case 1:
          tx = 0;
          ty = dy;
          return;
        case 2:
          tx = dx;
          ty = dy;
          return;
      }
  }
}

void JPS2DNeib::FNeib(int dx, int dy, int norm1, int dev, int& fx, int& fy, int& nx, int& ny) {
  switch (norm1) {
    case 1:
      switch (dev) {
        case 0:
          fx = 0;
          fy = 1;
          break;
        case 1:
          fx = 0;
          fy = -1;
          break;
      }

      // switch order if different direction
      if (dx == 0) fx = fy, fy = 0;

      nx = dx + fx;
      ny = dy + fy;
      return;
    case 2:
      switch (dev) {
        case 0:
          fx = -dx;
          fy = 0;
          nx = -dx;
          ny = dy;
          return;
        case 1:
          fx = 0;
          fy = -dy;
          nx = dx;
          ny = -dy;
          return;
      }
  }
}

constexpr int JPS3DNeib::nsz[4][2];

JPS3DNeib::JPS3DNeib() {
  int id = 0;
  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
        for (int dev = 0; dev < nsz[norm1][0]; ++dev)
          Neib(dx, dy, dz, norm1, dev, ns[id][0][dev], ns[id][1][dev], ns[id][2][dev]);
        for (int dev = 0; dev < nsz[norm1][1]; ++dev) {
          FNeib(dx, dy, dz, norm1, dev, f1[id][0][dev], f1[id][1][dev], f1[id][2][dev],
                f2[id][0][dev], f2[id][1][dev], f2[id][2][dev]);
        }
        id++;
      }
    }
  }
}

void JPS3DNeib::Neib(int dx, int dy, int dz, int norm1, int dev, int& tx, int& ty, int& tz) {
  switch (norm1) {
    case 0:
      switch (dev) {
        case 0:
          tx = 1;
          ty = 0;
          tz = 0;
          return;
        case 1:
          tx = -1;
          ty = 0;
          tz = 0;
          return;
        case 2:
          tx = 0;
          ty = 1;
          tz = 0;
          return;
        case 3:
          tx = 1;
          ty = 1;
          tz = 0;
          return;
        case 4:
          tx = -1;
          ty = 1;
          tz = 0;
          return;
        case 5:
          tx = 0;
          ty = -1;
          tz = 0;
          return;
        case 6:
          tx = 1;
          ty = -1;
          tz = 0;
          return;
        case 7:
          tx = -1;
          ty = -1;
          tz = 0;
          return;
        case 8:
          tx = 0;
          ty = 0;
          tz = 1;
          return;
        case 9:
          tx = 1;
          ty = 0;
          tz = 1;
          return;
        case 10:
          tx = -1;
          ty = 0;
          tz = 1;
          return;
        case 11:
          tx = 0;
          ty = 1;
          tz = 1;
          return;
        case 12:
          tx = 1;
          ty = 1;
          tz = 1;
          return;
        case 13:
          tx = -1;
          ty = 1;
          tz = 1;
          return;
        case 14:
          tx = 0;
          ty = -1;
          tz = 1;
          return;
        case 15:
          tx = 1;
          ty = -1;
          tz = 1;
          return;
        case 16:
          tx = -1;
          ty = -1;
          tz = 1;
          return;
        case 17:
          tx = 0;
          ty = 0;
          tz = -1;
          return;
        case 18:
          tx = 1;
          ty = 0;
          tz = -1;
          return;
        case 19:
          tx = -1;
          ty = 0;
          tz = -1;
          return;
        case 20:
          tx = 0;
          ty = 1;
          tz = -1;
          return;
        case 21:
          tx = 1;
          ty = 1;
          tz = -1;
          return;
        case 22:
          tx = -1;
          ty = 1;
          tz = -1;
          return;
        case 23:
          tx = 0;
          ty = -1;
          tz = -1;
          return;
        case 24:
          tx = 1;
          ty = -1;
          tz = -1;
          return;
        case 25:
          tx = -1;
          ty = -1;
          tz = -1;
          return;
      }
    case 1:
      tx = dx;
      ty = dy;
      tz = dz;
      return;
    case 2:
      switch (dev) {
        case 0:
          if (dz == 0) {
            tx = 0;
            ty = dy;
            tz = 0;
            return;
          } else {
            tx = 0;
            ty = 0;
            tz = dz;
            return;
          }
        case 1:
          if (dx == 0) {
            tx = 0;
            ty = dy;
            tz = 0;
            return;
          } else {
            tx = dx;
            ty = 0;
            tz = 0;
            return;
          }
        case 2:
          tx = dx;
          ty = dy;
          tz = dz;
          return;
      }
    case 3:
      switch (dev) {
        case 0:
          tx = dx;
          ty = 0;
          tz = 0;
          return;
        case 1:
          tx = 0;
          ty = dy;
          tz = 0;
          return;
        case 2:
          tx = 0;
          ty = 0;
          tz = dz;
          return;
        case 3:
          tx = dx;
          ty = dy;
          tz = 0;
          return;
        case 4:
          tx = dx;
          ty = 0;
          tz = dz;
          return;
        case 5:
          tx = 0;
          ty = dy;
          tz = dz;
          return;
        case 6:
          tx = dx;
          ty = dy;
          tz = dz;
          return;
      }
  }
}

void JPS3DNeib::FNeib(int dx, int dy, int dz, int norm1, int dev, int& fx, int& fy, int& fz,
                      int& nx, int& ny, int& nz) {
  switch (norm1) {
    case 1:
      switch (dev) {
        case 0:
          fx = 0;
          fy = 1;
          fz = 0;
          break;
        case 1:
          fx = 0;
          fy = -1;
          fz = 0;
          break;
        case 2:
          fx = 1;
          fy = 0;
          fz = 0;
          break;
        case 3:
          fx = 1;
          fy = 1;
          fz = 0;
          break;
        case 4:
          fx = 1;
          fy = -1;
          fz = 0;
          break;
        case 5:
          fx = -1;
          fy = 0;
          fz = 0;
          break;
        case 6:
          fx = -1;
          fy = 1;
          fz = 0;
          break;
        case 7:
          fx = -1;
          fy = -1;
          fz = 0;
          break;
      }
      nx = fx;
      ny = fy;
      nz = dz;
      // switch order if different direction
      if (dx != 0) {
        fz = fx;
        fx = 0;
        nz = fz;
        nx = dx;
      }
      if (dy != 0) {
        fz = fy;
        fy = 0;
        nz = fz;
        ny = dy;
      }
      return;
    case 2:
      if (dx == 0) {
        switch (dev) {
          case 0:
            fx = 0;
            fy = 0;
            fz = -dz;
            nx = 0;
            ny = dy;
            nz = -dz;
            return;
          case 1:
            fx = 0;
            fy = -dy;
            fz = 0;
            nx = 0;
            ny = -dy;
            nz = dz;
            return;
          case 2:
            fx = 1;
            fy = 0;
            fz = 0;
            nx = 1;
            ny = dy;
            nz = dz;
            return;
          case 3:
            fx = -1;
            fy = 0;
            fz = 0;
            nx = -1;
            ny = dy;
            nz = dz;
            return;
          case 4:
            fx = 1;
            fy = 0;
            fz = -dz;
            nx = 1;
            ny = dy;
            nz = -dz;
            return;
          case 5:
            fx = 1;
            fy = -dy;
            fz = 0;
            nx = 1;
            ny = -dy;
            nz = dz;
            return;
          case 6:
            fx = -1;
            fy = 0;
            fz = -dz;
            nx = -1;
            ny = dy;
            nz = -dz;
            return;
          case 7:
            fx = -1;
            fy = -dy;
            fz = 0;
            nx = -1;
            ny = -dy;
            nz = dz;
            return;
          // Extras
          case 8:
            fx = 1;
            fy = 0;
            fz = 0;
            nx = 1;
            ny = dy;
            nz = 0;
            return;
          case 9:
            fx = 1;
            fy = 0;
            fz = 0;
            nx = 1;
            ny = 0;
            nz = dz;
            return;
          case 10:
            fx = -1;
            fy = 0;
            fz = 0;
            nx = -1;
            ny = dy;
            nz = 0;
            return;
          case 11:
            fx = -1;
            fy = 0;
            fz = 0;
            nx = -1;
            ny = 0;
            nz = dz;
            return;
        }
      } else if (dy == 0) {
        switch (dev) {
          case 0:
            fx = 0;
            fy = 0;
            fz = -dz;
            nx = dx;
            ny = 0;
            nz = -dz;
            return;
          case 1:
            fx = -dx;
            fy = 0;
            fz = 0;
            nx = -dx;
            ny = 0;
            nz = dz;
            return;
          case 2:
            fx = 0;
            fy = 1;
            fz = 0;
            nx = dx;
            ny = 1;
            nz = dz;
            return;
          case 3:
            fx = 0;
            fy = -1;
            fz = 0;
            nx = dx;
            ny = -1;
            nz = dz;
            return;
          case 4:
            fx = 0;
            fy = 1;
            fz = -dz;
            nx = dx;
            ny = 1;
            nz = -dz;
            return;
          case 5:
            fx = -dx;
            fy = 1;
            fz = 0;
            nx = -dx;
            ny = 1;
            nz = dz;
            return;
          case 6:
            fx = 0;
            fy = -1;
            fz = -dz;
            nx = dx;
            ny = -1;
            nz = -dz;
            return;
          case 7:
            fx = -dx;
            fy = -1;
            fz = 0;
            nx = -dx;
            ny = -1;
            nz = dz;
            return;
          // Extras
          case 8:
            fx = 0;
            fy = 1;
            fz = 0;
            nx = dx;
            ny = 1;
            nz = 0;
            return;
          case 9:
            fx = 0;
            fy = 1;
            fz = 0;
            nx = 0;
            ny = 1;
            nz = dz;
            return;
          case 10:
            fx = 0;
            fy = -1;
            fz = 0;
            nx = dx;
            ny = -1;
            nz = 0;
            return;
          case 11:
            fx = 0;
            fy = -1;
            fz = 0;
            nx = 0;
            ny = -1;
            nz = dz;
            return;
        }
      } else {  // dz==0
        switch (dev) {
          case 0:
            fx = 0;
            fy = -dy;
            fz = 0;
            nx = dx;
            ny = -dy;
            nz = 0;
            return;
          case 1:
            fx = -dx;
            fy = 0;
            fz = 0;
            nx = -dx;
            ny = dy;
            nz = 0;
            return;
          case 2:
            fx = 0;
            fy = 0;
            fz = 1;
            nx = dx;
            ny = dy;
            nz = 1;
            return;
          case 3:
            fx = 0;
            fy = 0;
            fz = -1;
            nx = dx;
            ny = dy;
            nz = -1;
            return;
          case 4:
            fx = 0;
            fy = -dy;
            fz = 1;
            nx = dx;
            ny = -dy;
            nz = 1;
            return;
          case 5:
            fx = -dx;
            fy = 0;
            fz = 1;
            nx = -dx;
            ny = dy;
            nz = 1;
            return;
          case 6:
            fx = 0;
            fy = -dy;
            fz = -1;
            nx = dx;
            ny = -dy;
            nz = -1;
            return;
          case 7:
            fx = -dx;
            fy = 0;
            fz = -1;
            nx = -dx;
            ny = dy;
            nz = -1;
            return;
          // Extras
          case 8:
            fx = 0;
            fy = 0;
            fz = 1;
            nx = dx;
            ny = 0;
            nz = 1;
            return;
          case 9:
            fx = 0;
            fy = 0;
            fz = 1;
            nx = 0;
            ny = dy;
            nz = 1;
            return;
          case 10:
            fx = 0;
            fy = 0;
            fz = -1;
            nx = dx;
            ny = 0;
            nz = -1;
            return;
          case 11:
            fx = 0;
            fy = 0;
            fz = -1;
            nx = 0;
            ny = dy;
            nz = -1;
            return;
        }
      }
      break;
    case 3:
      switch (dev) {
        case 0:
          fx = -dx;
          fy = 0;
          fz = 0;
          nx = -dx;
          ny = dy;
          nz = dz;
          return;
        case 1:
          fx = 0;
          fy = -dy;
          fz = 0;
          nx = dx;
          ny = -dy;
          nz = dz;
          return;
        case 2:
          fx = 0;
          fy = 0;
          fz = -dz;
          nx = dx;
          ny = dy;
          nz = -dz;
          return;
        // Need to check up to here for forced!
        case 3:
          fx = 0;
          fy = -dy;
          fz = -dz;
          nx = dx;
          ny = -dy;
          nz = -dz;
          return;
        case 4:
          fx = -dx;
          fy = 0;
          fz = -dz;
          nx = -dx;
          ny = dy;
          nz = -dz;
          return;
        case 5:
          fx = -dx;
          fy = -dy;
          fz = 0;
          nx = -dx;
          ny = -dy;
          nz = dz;
          return;
        // Extras
        case 6:
          fx = -dx;
          fy = 0;
          fz = 0;
          nx = -dx;
          ny = 0;
          nz = dz;
          return;
        case 7:
          fx = -dx;
          fy = 0;
          fz = 0;
          nx = -dx;
          ny = dy;
          nz = 0;
          return;
        case 8:
          fx = 0;
          fy = -dy;
          fz = 0;
          nx = 0;
          ny = -dy;
          nz = dz;
          return;
        case 9:
          fx = 0;
          fy = -dy;
          fz = 0;
          nx = dx;
          ny = -dy;
          nz = 0;
          return;
        case 10:
          fx = 0;
          fy = 0;
          fz = -dz;
          nx = 0;
          ny = dy;
          nz = -dz;
          return;
        case 11:
          fx = 0;
          fy = 0;
          fz = -dz;
          nx = dx;
          ny = 0;
          nz = -dz;
          return;
      }
  }
}