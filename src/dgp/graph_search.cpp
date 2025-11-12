#include "dgp/graph_search.hpp"
#include <cmath>

using namespace mighty;
using namespace termcolor;

typedef timer::Timer MyTimer;

// Help functions
static inline double clamp01(double x) { return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x); }
static inline double hypot2d(double x, double y) { return std::sqrt(x * x + y * y); }

GraphSearch::GraphSearch(const int *cMap, const std::shared_ptr<mighty::VoxelMapUtil> &map_util, int xDim, int yDim, int zDim, double eps, bool verbose, std::string global_planner, double w_unknown, double w_align, double decay_len_cells, double w_side) : cMap_(cMap), map_util_(map_util), xDim_(xDim), yDim_(yDim), zDim_(zDim), eps_(eps), verbose_(verbose), global_planner_(global_planner), w_unknown_(w_unknown), w_align_(w_align), decay_len_cells_(decay_len_cells), w_side_(w_side)
{
  hm_.assign(xDim_ * yDim_ * zDim_, nullptr);
  seen_.assign(xDim_ * yDim_ * zDim_, false);

  // Set 3D neighbors
  for (int x = -1; x <= 1; x++)
  {
    for (int y = -1; y <= 1; y++)
    {
      for (int z = -1; z <= 1; z++)
      {
        if (x == 0 && y == 0 && z == 0)
          continue;
        ns_.push_back(std::vector<int>{x, y, z});
      }
    }
  }

  // Set neighbors (not including the diagonal neighbors)
  // ns_.push_back(std::vector<int>{1, 0, 0});
  // ns_.push_back(std::vector<int>{-1, 0, 0});
  // ns_.push_back(std::vector<int>{0, 1, 0});
  // ns_.push_back(std::vector<int>{0, -1, 0});
  // ns_.push_back(std::vector<int>{0, 0, 1});
  // ns_.push_back(std::vector<int>{0, 0, -1});

  jn3d_ = std::make_shared<JPS3DNeib>();
}

inline int GraphSearch::coordToId(int x, int y, int z) const
{
  return x + y * xDim_ + z * xDim_ * yDim_;
}

inline bool GraphSearch::isFree(int x, int y, int z) const
{
  // return x > 0 && x < xDim_ && y > 0 && y < yDim_ && z > 0 && z < zDim_ && cMap_[coordToId(x, y, z)] == val_free_;
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ && z >= 0 && z < zDim_ && cMap_[coordToId(x, y, z)] <= val_free_;
}

inline bool GraphSearch::isOccupied(int x, int y, int z) const
{
  return x < 0 || x >= xDim_ || y < 0 || y >= yDim_ || z < 0 || z >= zDim_ || cMap_[coordToId(x, y, z)] == val_occupied_;
}

inline bool GraphSearch::isUnknown(int x, int y, int z) const
{
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ && z >= 0 && z < zDim_ && cMap_[coordToId(x, y, z)] == val_unknown_;
}

inline double GraphSearch::getHeur(int x, int y, int z) const
{
  return eps_ * std::sqrt((x - xGoal_) * (x - xGoal_) + (y - yGoal_) * (y - yGoal_) + (z - zGoal_) * (z - zGoal_));
}

bool GraphSearch::plan(int xStart, int yStart, int zStart, int xGoal, int yGoal, int zGoal, double initial_g, double &global_planning_time, double &dgp_static_jps_time, double &dgp_check_path_time, double &dgp_dynamic_astar_time, double &dgp_recover_path_time, double current_time, const Vec3f &start_vel, int max_expand, int timeout_duration_ms)
{
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
  StatePtr currNode_ptr = std::make_shared<State>(State(start_id, xStart, yStart, zStart, 0, 0, 0, 0));

  currNode_ptr->g = initial_g;
  currNode_ptr->h = getHeur(xStart, yStart, zStart);

  // Set current time
  current_time_ = current_time;

  // Start a timer
  MyTimer global_planning_timer(true);

  // Run the planner
  bool result = select_planner(currNode_ptr, max_expand, start_id, goal_id, std::chrono::milliseconds(timeout_duration_ms));

  // Print the time taken [ms]
  global_planning_time_ = global_planning_timer.getElapsedMicros() / 1000.0;

  // Assign the time taken
  global_planning_time = global_planning_time_;
  dgp_static_jps_time = dgp_static_jps_time_;
  dgp_check_path_time = dgp_check_path_time_;
  dgp_dynamic_astar_time = dgp_dynamic_astar_time_;
  dgp_recover_path_time = dgp_recover_path_time_;

  return result;
}

bool GraphSearch::select_planner(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration)
{

  if (global_planner_ == "sjps") // Static JPS planner
  {
    return static_jps_plan(currNode_ptr, max_expand, start_id, goal_id, timeout_duration);
  }
  else // Error
  {
    printf("Unknown planner: %s\n", global_planner_.c_str());
    return false;
  }
}

bool GraphSearch::static_jps_plan(StatePtr &currNode_ptr, int max_expand, int start_id, int goal_id, std::chrono::milliseconds timeout_duration)
{

  // Record the start time
  auto start_time = std::chrono::steady_clock::now();

  if (!currNode_ptr)
  {
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

  while (true)
  {
    expand_iteration++;

    // Check if timeout has occurred
    auto current_time = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time) > timeout_duration)
    {
      std::cerr << "Timeout occurred. Exiting safely.\n";
      return false;
    }

    if (pq_.empty())
    {
      std::cerr << "Error: Priority queue is empty!" << std::endl;
      return false;
    }

    // get element with smallest cost
    currNode_ptr = pq_.top();
    pq_.pop();
    currNode_ptr->closed = true; // Add to closed list

    if (currNode_ptr->id == goal_id)
    {
      if (verbose_)
        printf("Goal Reached!\n");
      break;
    }

    // // early stopping for planning horizon
    if (currNode_ptr->g > 70.0)
    {
      break;
    }

    std::vector<int> succ_ids;
    std::vector<double> succ_costs;

    // getJpsSucc(currNode_ptr, succ_ids, succ_costs);
    getSucc(currNode_ptr, succ_ids, succ_costs);

    // Process successors
    for (int s = 0; s < (int)succ_ids.size(); s++)
    {

      // see if we can improve the value of succstate
      StatePtr &child_ptr = atHM(succ_ids[s]);
      if (!child_ptr)
      {
        // (Optional) reconstruct or skip
        fprintf(stderr, "[WARN] child_ptr null at id=%d\n", succ_ids[s]);
        continue;
      }

      // Directional bias: only on the very first expansion
      // ----- Start-direction bias (first expansion only): penalize ONLY backward -----
      double penalty = 0.0;
      // if (currNode_ptr->id == start_id && start_vel_.squaredNorm() > 1e-8)
      // {
      //   Eigen::Vector3d pref = start_vel_.cast<double>();
      //   Eigen::Vector3d toGoal(xGoal_ - currNode_ptr->x,
      //                          yGoal_ - currNode_ptr->y,
      //                          zGoal_ - currNode_ptr->z);
      //   // If start_vel_ points >90° away from goal, flip it
      //   if (pref.dot(toGoal) < 0)
      //     pref = -pref;
      //   pref.normalize();

      //   Eigen::Vector3d step(static_cast<double>(child_ptr->x - currNode_ptr->x),
      //                        static_cast<double>(child_ptr->y - currNode_ptr->y),
      //                        static_cast<double>(child_ptr->z - currNode_ptr->z));
      //   double n = step.norm();
      //   if (n > 1e-9)
      //   {
      //     double cosang = (step / n).dot(pref);      // [-1, 1]
      //     double k_back = 1.0;                       // try 0.5–1.0
      //     penalty = k_back * std::max(0.0, -cosang); // penalize only backward
      //   }
      // }

      double tentative_gval = currNode_ptr->g + succ_costs[s] + penalty;

      if (tentative_gval < child_ptr->g)
      {
        child_ptr->parentId = currNode_ptr->id; // Assign new parent
        child_ptr->g = tentative_gval;          // Update gval

        // double fval = child_ptr->g + child_ptr->h;

        // if currently in OPEN, update
        if (child_ptr->opened && !child_ptr->closed)
        {
          pq_.increase(child_ptr->heapkey); // update heap
          child_ptr->dx = (child_ptr->x - currNode_ptr->x);
          child_ptr->dy = (child_ptr->y - currNode_ptr->y);
          child_ptr->dz = (child_ptr->z - currNode_ptr->z);
          if (child_ptr->dx != 0)
            child_ptr->dx /= std::abs(child_ptr->dx);
          if (child_ptr->dy != 0)
            child_ptr->dy /= std::abs(child_ptr->dy);
          if (child_ptr->dz != 0)
            child_ptr->dz /= std::abs(child_ptr->dz);
        }
        // if currently in CLOSED
        else if (child_ptr->opened && child_ptr->closed)
        {
          // printf("STATIC JPS ERROR!\n");
        }
        else // new node, add to heap
        {
          // printf("add to open set: %d, %d\n", child_ptr->x, child_ptr->y);
          child_ptr->heapkey = pq_.push(child_ptr);
          child_ptr->opened = true;
        }
      } //
    } // Process successors

    if (max_expand > 0 && expand_iteration >= max_expand)
    {
      if (verbose_)
      {
        printf("max_expandStep [%d] Reached\n\n", max_expand);
      }
      return false;
    }
  }

  if (verbose_)
  {
    // printf("goal g: %f, h: %f!\n", currNode_ptr->g, currNode_ptr->h);
    // printf("Expand [%d] nodes!\n", expand_iteration);
  }

  path_ = recoverPath(currNode_ptr, start_id);
  // updateGValues(); // this gives correct g values

  return true;
}

std::vector<StatePtr> GraphSearch::removeLinePts(const std::vector<StatePtr> &path)
{
  if (path.size() < 3)
    return path;

  std::vector<StatePtr> new_path;
  new_path.push_back(path.front());
  for (unsigned int i = 1; i < path.size() - 1; i++)
  {
    Vecf<3> path_i_p_1 = Vecf<3>(path[i + 1]->x, path[i + 1]->y, path[i + 1]->z);
    Vecf<3> path_i = Vecf<3>(path[i]->x, path[i]->y, path[i]->z);
    Vecf<3> path_i_m_1 = Vecf<3>(path[i - 1]->x, path[i - 1]->y, path[i - 1]->z);
    Vecf<3> p = (path_i_p_1 - path_i) - (path_i - path_i_m_1);
    if (fabs(p(0)) + fabs(p(1)) + fabs(p(2)) > 1e-2)
    {
      new_path.push_back(path[i]);
    }
  }
  new_path.push_back(path.back());
  return new_path;
}

std::vector<StatePtr> GraphSearch::removeCornerPts(const std::vector<StatePtr> &path)
{
  if (path.size() < 2)
    return path;

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

  for (unsigned int i = 1; i < path.size() - 1; i++)
  {
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
    else
    {
      optimized_path.push_back(path[i]);
      cost1 = (pose1 - pose2).norm();
      prev_pose = pose1;
    }
  }

  optimized_path.push_back(path.back());
  return optimized_path;
}

//// RIGHT NOW THE JPS GIVES BACK A VERY SPARSE PATH, SO TO CHECK COLLISION AGAINST DYNAMIC OBSTACLES, WE NEED TO INTERPOLATE THE PATH
void GraphSearch::updateGValues()
{

  // first reverse the path
  std::reverse(path_.begin(), path_.end());

  // first create a path with float values
  std::vector<Vecf<3>> path_float;
  for (int i = 0; i < path_.size(); i++)
  {
    path_float.push_back(map_util_->intToFloat(Veci<3>(path_[i]->x, path_[i]->y, path_[i]->z)));
  }

  // convert the path back to int and update g values
  std::vector<StatePtr> new_path_int;
  for (int i = 0; i < path_float.size(); i++)
  {
    Veci<3> temp = map_util_->floatToInt(Vecf<3>(path_float[i](0), path_float[i](1), path_float[i](2)));
    StatePtr temp_ptr = std::make_shared<State>(State(coordToId(temp(0), temp(1), temp(2)), temp(0), temp(1), temp(2), 0, 0, 0)); // for now dx, dy, dz are 0 (TODO: this could be problematic?)

    // compute total path length and update g
    if (i > 0)
    {
      Vecf<3> p1 = path_float[i - 1];
      Vecf<3> p2 = path_float[i];
      temp_ptr->g = new_path_int[i - 1]->g + (p2 - p1).norm();
    }
    else
    {
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
double GraphSearch::computeNodeTime(const Eigen::Vector3d &pf,
                                    const Eigen::Vector3d &p0,
                                    const Eigen::Vector3d &v0,
                                    StatePtr &currNode_ptr)
{
  // Initialize final velocity vector (for each axis)
  Eigen::Vector3d vf = Eigen::Vector3d::Zero();

  // Absolute displacement between nodes
  Eigen::Vector3d dist = (pf - p0).cwiseAbs();

  // Will store travel time for each dimension; the overall node time is the max.
  std::vector<double> node_times;

  // Small tolerance for checking if v0 is nearly equal to v_max.
  double epsilon = 1e-6;

  // Iterate over the dimensions (x, y, z)
  for (int i = 0; i < 3; i++)
  {
    double a_max = a_max_3d_(i);  // maximum acceleration in this dimension
    double v_max = v_max_3d_(i);  // maximum velocity in this dimension (assumed defined)
    double delta = pf(i) - p0(i); // signed displacement along axis i
    double node_time = 0;         // time to traverse this dimension

    // std::cout << "----- Dimension " << i << " -----" << std::endl;
    // std::cout << "delta = " << delta << ", dist = " << dist(i) << std::endl;

    // If delta is very small, we can assume no motion is needed.
    if (std::abs(delta) < epsilon)
    {
      // std::cout << "Delta is negligible, no motion needed." << std::endl;
      node_time = 0;
      vf(i) = v0(i);
    }
    else if (delta > 0) // Desired displacement is positive.
    {
      // std::cout << "Positive displacement branch." << std::endl;
      if (v0(i) >= 0) // Already moving in the positive direction
      {
        // std::cout << "v0(" << i << ") > 0: v0 = " << v0(i) << std::endl;
        // If already cruising (and no deceleration is needed)
        if (v0(i) >= (v_max - epsilon))
        {
          // std::cout << "Already cruising: v0(" << i << ") = " << v0(i) << " >= v_max - epsilon (" << v_max - epsilon << ")" << std::endl;
          vf(i) = v_max;
          node_time = dist(i) / v_max;
          // std::cout << "Computed node_time (cruise) = " << node_time << std::endl;
        }
        else
        {
          double v_candidate = std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i));
          // std::cout << "Computed v_candidate = " << v_candidate << std::endl;
          if (v_candidate <= v_max)
          {
            // std::cout << "v_candidate <= v_max" << std::endl;
            vf(i) = v_candidate;
            node_time = (vf(i) - v0(i)) / a_max;
            // std::cout << "Computed node_time (acceleration only) = " << node_time << std::endl;
          }
          else
          {
            // std::cout << "v_candidate > v_max, splitting into acceleration and cruise phases." << std::endl;
            double d_accel = (v_max * v_max - v0(i) * v0(i)) / (2 * a_max);
            double t_accel = (v_max - v0(i)) / a_max;
            // std::cout << "d_accel = " << d_accel << ", t_accel = " << t_accel << std::endl;
            if (dist(i) > d_accel)
            {
              double d_cruise = dist(i) - d_accel;
              double t_cruise = d_cruise / v_max;
              node_time = t_accel + t_cruise;
              // std::cout << "d_cruise = " << d_cruise << ", t_cruise = " << t_cruise << std::endl;
              // std::cout << "Total node_time = " << node_time << std::endl;
            }
            else
            {
              // Fallback if numerical issues arise.
              // std::cout << "Fallback branch: dist <= d_accel." << std::endl;
              node_time = (std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i)) - v0(i)) / a_max;
              // std::cout << "Computed fallback node_time = " << node_time << std::endl;
            }
            vf(i) = v_max;
          }
        }
      }
      else // v0(i) < 0: initially moving in the wrong (or zero) direction.
      {
        // std::cout << "v0(" << i << ") <= 0: initial velocity = " << v0(i) << " (wrong direction or stationary)" << std::endl;
        double d_stop = (v0(i) * v0(i)) / (2 * a_max);
        double t_stop = (-v0(i)) / a_max;
        // std::cout << "d_stop = " << d_stop << ", t_stop = " << t_stop << std::endl;
        double d_accel = dist(i) + d_stop;
        double v_candidate = std::sqrt(2 * a_max * d_accel);
        // std::cout << "Phase 2: d_accel = " << d_accel << ", v_candidate = " << v_candidate << std::endl;
        if (v_candidate <= v_max)
        {
          double t_accel = std::sqrt(2 * d_accel / a_max);
          node_time = t_stop + t_accel;
          vf(i) = v_candidate;
          // std::cout << "Computed node_time (decelerate then accelerate) = " << node_time << std::endl;
        }
        else
        {
          // std::cout << "v_candidate > v_max in deceleration branch, splitting phases." << std::endl;
          double d_to_vmax = (v_max * v_max) / (2 * a_max);
          double t_accel = v_max / a_max;
          // std::cout << "d_to_vmax = " << d_to_vmax << ", t_accel = " << t_accel << std::endl;
          if (d_accel > d_to_vmax)
          {
            double d_cruise = d_accel - d_to_vmax;
            double t_cruise = d_cruise / v_max;
            node_time = t_stop + t_accel + t_cruise;
            // std::cout << "d_cruise = " << d_cruise << ", t_cruise = " << t_cruise << std::endl;
          }
          else
          {
            node_time = t_stop + t_accel;
          }
          vf(i) = v_max;
          // std::cout << "Total node_time = " << node_time << std::endl;
        }
      }
    }
    else // delta <= 0: Desired displacement is negative.
    {
      // std::cout << "Negative displacement branch." << std::endl;
      if (v0(i) <= 0) // Already moving in the negative direction.
      {
        // std::cout << "v0(" << i << ") < 0: already moving negative, v0 = " << v0(i) << std::endl;
        if (std::fabs(v0(i)) >= (v_max - epsilon))
        {
          // std::cout << "Already cruising in the negative direction: |v0| >= v_max - epsilon" << std::endl;
          vf(i) = -v_max;
          node_time = dist(i) / v_max;
          // std::cout << "Computed node_time (cruise negative) = " << node_time << std::endl;
        }
        else
        {
          double v_candidate = std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i));
          v_candidate = -v_candidate; // assign the negative sign
          // std::cout << "Computed negative v_candidate = " << v_candidate << std::endl;
          if (std::fabs(v_candidate) <= v_max)
          {
            vf(i) = v_candidate;
            node_time = (std::fabs(vf(i)) - std::fabs(v0(i))) / a_max;
            // std::cout << "Computed node_time (negative acceleration only) = " << node_time << std::endl;
          }
          else
          {
            // std::cout << "v_candidate exceeds v_max in negative branch, splitting phases." << std::endl;
            double d_accel = ((v_max * v_max) - (v0(i) * v0(i))) / (2 * a_max);
            double t_accel = (v_max - std::fabs(v0(i))) / a_max;
            // std::cout << "d_accel = " << d_accel << ", t_accel = " << t_accel << std::endl;
            if (dist(i) > d_accel)
            {
              double d_cruise = dist(i) - d_accel;
              double t_cruise = d_cruise / v_max;
              node_time = t_accel + t_cruise;
              // std::cout << "d_cruise = " << d_cruise << ", t_cruise = " << t_cruise << std::endl;
            }
            else
            {
              node_time = (std::fabs(std::sqrt(v0(i) * v0(i) + 2 * a_max * dist(i))) - std::fabs(v0(i))) / a_max;
              // std::cout << "Computed fallback node_time (negative) = " << node_time << std::endl;
            }
            vf(i) = -v_max;
          }
        }
      }
      else // v0(i) > 0: initially moving in the positive direction while desired displacement is negative.
      {
        // std::cout << "v0(" << i << ") >= 0: Need to decelerate and reverse, v0 = " << v0(i) << std::endl;
        double d_stop = (v0(i) * v0(i)) / (2 * a_max);
        double t_stop = v0(i) / a_max;
        // std::cout << "d_stop = " << d_stop << ", t_stop = " << t_stop << std::endl;
        double d_accel = dist(i) + d_stop;
        double v_candidate = std::sqrt(2 * a_max * d_accel);
        v_candidate = -v_candidate; // final velocity is negative.
        // std::cout << "d_accel = " << d_accel << ", negative v_candidate = " << v_candidate << std::endl;
        if (std::fabs(v_candidate) <= v_max)
        {
          double t_accel = std::sqrt(2 * d_accel / a_max);
          node_time = t_stop + t_accel;
          vf(i) = v_candidate;
          // std::cout << "Computed node_time (decelerate then reverse) = " << node_time << std::endl;
        }
        else
        {
          // std::cout << "v_candidate exceeds v_max in reversal branch, splitting phases." << std::endl;
          double d_to_vmax = (v_max * v_max) / (2 * a_max);
          double t_accel = v_max / a_max;
          // std::cout << "d_to_vmax = " << d_to_vmax << ", t_accel = " << t_accel << std::endl;
          if (d_accel > d_to_vmax)
          {
            double d_cruise = d_accel - d_to_vmax;
            double t_cruise = d_cruise / v_max;
            node_time = t_stop + t_accel + t_cruise;
            // std::cout << "d_cruise = " << d_cruise << ", t_cruise = " << t_cruise << std::endl;
          }
          else
          {
            node_time = t_stop + t_accel;
          }
          vf(i) = -v_max;
          // std::cout << "Computed total node_time = " << node_time << std::endl;
        }
      }
    }

    // Save the computed time for this axis.
    node_times.push_back(node_time);
    // std::cout << "Final node_time for dimension " << i << " = " << node_time << std::endl;
    // std::cout << "Final velocity for dimension " << i << " = " << vf(i) << std::endl;
  } // end for loop over dimensions

  // Update the current node's velocity with the computed final velocity vector.
  currNode_ptr->vel = vf;

  // Return the overall node time (the maximum among x, y, and z dimensions)
  double overall_time = *std::max_element(node_times.begin(), node_times.end());
  // std::cout << "Overall node time: " << overall_time << std::endl;
  return overall_time;
}

void GraphSearch::setStartAndGoal(const Vecf<3> &start, const Vecf<3> &goal)
{
  start_ = start;
  goal_ = goal;
}

void GraphSearch::setBounds(double max_values[3])
{
  v_max_ = max_values[0];
  v_max_3d_ << max_values[0], max_values[0], max_values[0];
  a_max_ = max_values[1];
  a_max_3d_ << max_values[1], max_values[1], max_values[1];
  j_max_ = max_values[2];
  j_max_3d_ << max_values[2], max_values[2], max_values[2];
}

std::vector<StatePtr> GraphSearch::recoverPath(StatePtr node, int start_id)
{
  std::vector<StatePtr> path;
  path.push_back(node);
  while (node && node->id != start_id)
  {
    node = atHM(node->parentId);
    path.push_back(node);
  }

  return path;
}

void GraphSearch::getSucc(const StatePtr &curr, std::vector<int> &succ_ids, std::vector<double> &succ_costs)
{
  succ_ids.clear();
  succ_costs.clear();

  // --- Build a forward reference once per node ---
  // Prefer start_vel_ when available; if it's anti-aligned with goal, flip it.
  Eigen::Vector2d vref(start_vel_(0), start_vel_(1));
  if (vref.squaredNorm() < 1e-9)
  {
    vref = Eigen::Vector2d(goal_(0) - start_(0), goal_(1) - start_(1));
  }
  else
  {
    Eigen::Vector2d vg(goal_(0) - start_(0), goal_(1) - start_(1));
    if (vg.squaredNorm() > 1e-9 && vref.dot(vg) < 0.0)
      vref = -vref; // never bias away from goal
  }
  const double vref_n = vref.norm();

  for (const auto &d : ns_)
  {
    int new_x = curr->x + d[0];
    int new_y = curr->y + d[1];
    int new_z = curr->z + d[2];
    if (isOccupied(new_x, new_y, new_z))
      continue;

    int new_id = coordToId(new_x, new_y, new_z);
    if (new_id < 0 || new_id >= (int)hm_.size())
      continue;

    if (!seen_[new_id])
    {
      seen_[new_id] = true;
      hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, new_z, d[0], d[1], d[2]);
      if (new_id < 0 || new_id >= (int)hm_.size())
      {
        fprintf(stderr, "[BUG] Creation failed id=%d\n", new_id);
        std::abort();
      }
      hm_[new_id]->h = getHeur(new_x, new_y, new_z);
    }

    // Base geometric step cost (1, √2, √3)
    const double base = std::sqrt(double(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
    double step_cost = base;
    // const double dist_from_start = std::sqrt((new_x - start_(0)) * (new_x - start_(0)) + (new_y - start_(1)) * (new_y - start_(1)) + (new_z - start_(2)) * (new_z - start_(2)));
    // // -------- Alignment (XY) + optional side tie-breaker --------
    // if ((w_align_ > 0.0 || w_side_ > 0.0) && dist_from_start < 10.0)
    // {
    //   Eigen::Vector2d step2d(d[0], d[1]);
    //   const double step_n = step2d.norm();
    //   if (vref_n > 1e-9 && step_n > 1e-9)
    //   {
    //     // Decay by *pure* g if you track it; else use distance_score or path length in cells.
    //     // Here we assume curr->g holds pure accumulated step costs.
    //     const double decay = std::exp(-curr->g / std::max(1e-9, decay_len_cells_));

    //     const double cosang = std::clamp(step2d.dot(vref) / (step_n * vref_n), -1.0, 1.0);
    //     // Misalignment penalty (0 when perfectly aligned; π gives 2*w_align_*step_n*decay)
    //     step_cost += w_align_ * (1.0 - cosang) * decay * step_n;

    //     // Optional: left/right tie-breaker (small!); gated to forward-ish (cos>=0)
    //     if (w_side_ > 0.0)
    //     {
    //       const double crossz = (vref.x() * step2d.y() - vref.y() * step2d.x()) / (vref_n * step_n); // [-1,1]
    //       Eigen::Vector2d vg(goal_(0) - start_(0), goal_(1) - start_(1));
    //       const double sign_pref = (vref.x() * vg.y() - vref.y() * vg.x()) >= 0.0 ? 1.0 : -1.0;
    //       const double cospos = std::max(0.0, cosang);
    //       step_cost += -w_side_ * sign_pref * crossz * cospos * decay * step_n;
    //     }
    //   }
    // }

    // // -------- If unknown, add a penalty --------
    // if (isUnknown(new_x, new_y, new_z))
    // {
    //   step_cost += w_unknown_ * base;
    // }

    succ_ids.push_back(new_id);
    succ_costs.push_back(step_cost);

    // succ_ids.push_back(new_id);
    // succ_costs.push_back(std::sqrt((new_x - curr->x) * (new_x - curr->x) + (new_y - curr->y) * (new_y - curr->y) + (new_z - curr->z) * (new_z - curr->z)));
  }
}

void GraphSearch::getJpsSucc(const StatePtr &curr, std::vector<int> &succ_ids, std::vector<double> &succ_costs)
{

  const int norm1 = std::abs(curr->dx) + std::abs(curr->dy) + std::abs(curr->dz);
  int num_neib = jn3d_->nsz[norm1][0];
  int num_fneib = jn3d_->nsz[norm1][1];
  int id = (curr->dx + 1) + 3 * (curr->dy + 1) + 9 * (curr->dz + 1);

  for (int dev = 0; dev < num_neib + num_fneib; ++dev)
  {
    int new_x, new_y, new_z;
    int dx, dy, dz;
    if (dev < num_neib)
    {
      dx = jn3d_->ns[id][0][dev];
      dy = jn3d_->ns[id][1][dev];
      dz = jn3d_->ns[id][2][dev];
      if (!jump(curr->x, curr->y, curr->z,
                dx, dy, dz, new_x, new_y, new_z))
        continue;
    }
    else
    {
      int nx = curr->x + jn3d_->f1[id][0][dev - num_neib];
      int ny = curr->y + jn3d_->f1[id][1][dev - num_neib];
      int nz = curr->z + jn3d_->f1[id][2][dev - num_neib];
      if (isOccupied(nx, ny, nz))
      {
        dx = jn3d_->f2[id][0][dev - num_neib];
        dy = jn3d_->f2[id][1][dev - num_neib];
        dz = jn3d_->f2[id][2][dev - num_neib];
        if (!jump(curr->x, curr->y, curr->z,
                  dx, dy, dz, new_x, new_y, new_z))
          continue;
      }
      else
        continue;
    }

    int new_id = coordToId(new_x, new_y, new_z);
    if (!seen_[new_id])
    {
      seen_[new_id] = true;
      hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, new_z, dx, dy, dz);
      if (new_id < 0 || new_id >= (int)hm_.size())
      {
        fprintf(stderr, "[BUG] Creation failed id=%d\n", new_id);
        std::abort();
      }
      hm_[new_id]->h = getHeur(new_x, new_y, new_z);
    }

    // // Base geometric step cost:
    // double base = std::sqrt((new_x - curr->x) * (new_x - curr->x) +
    //                         (new_y - curr->y) * (new_y - curr->y) +
    //                         (new_z - curr->z) * (new_z - curr->z));

    // // --- Velocity alignment penalty (XY plane, decays with distance from start) ---
    // Eigen::Vector2d vref(start_vel_(0), start_vel_(1));
    // if (vref.squaredNorm() < 1e-9)
    // {
    //   // fallback: toward goal if we're nearly stopped
    //   vref = Eigen::Vector2d(goal_(0) - start_(0), goal_(1) - start_(1));
    // }
    // double vref_n = vref.norm();
    // Eigen::Vector2d step(new_x - curr->x, new_y - curr->y);
    // double step_n = step.norm();

    // double align_pen = 0.0, side_pen = 0.0;
    // if (vref_n > 1e-9 && step_n > 1e-9 && (w_align_ > 0.0 || w_side_ > 0.0))
    // {
    //   double cosang = clamp01(step.dot(vref) / (step_n * vref_n));
    //   double decay = std::exp(-curr->g / decay_cells_); // curr->g is in *cells*
    //   align_pen = w_align_ * (1.0 - cosang) * decay;    // 0 for perfectly aligned

    //   // Side tie-break: prefer consistently left or right relative to vref,
    //   // determined by geometry: sign( cross(vref, goal-start) )
    //   Eigen::Vector2d vg(goal_(0) - start_(0), goal_(1) - start_(1));
    //   double side_pref = (vref.x() * vg.y() - vref.y() * vg.x());
    //   double crossz = (vref.x() * step.y() - vref.y() * step.x()) / (vref_n * step_n); // [-1,1]
    //   double sign_pref = (side_pref >= 0.0) ? 1.0 : -1.0;                              // stable across replans
    //   side_pen = -w_side_ * sign_pref * crossz * decay;                                // cheaper on preferred side
    // }

    // // Push final cost (base + small penalties)
    // succ_ids.push_back(new_id);
    // succ_costs.push_back(base + align_pen + side_pen);

    succ_ids.push_back(new_id);
    succ_costs.push_back(std::sqrt((new_x - curr->x) * (new_x - curr->x) +
                                   (new_y - curr->y) * (new_y - curr->y) +
                                   (new_z - curr->z) * (new_z - curr->z)));
  }
}

bool GraphSearch::jump(int x, int y, int z, int dx, int dy, int dz, int &new_x, int &new_y, int &new_z)
{
  new_x = x + dx;
  new_y = y + dy;
  new_z = z + dz;
  if (!isFree(new_x, new_y, new_z))
    return false;

  if (new_x == xGoal_ && new_y == yGoal_ && new_z == zGoal_)
    return true;

  if (hasForced(new_x, new_y, new_z, dx, dy, dz))
    return true;

  const int id = (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1);
  const int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
  int num_neib = jn3d_->nsz[norm1][0];
  for (int k = 0; k < num_neib - 1; ++k)
  {
    int new_new_x, new_new_y, new_new_z;
    if (jump(new_x, new_y, new_z,
             jn3d_->ns[id][0][k], jn3d_->ns[id][1][k], jn3d_->ns[id][2][k],
             new_new_x, new_new_y, new_new_z))
      return true;
  }

  return jump(new_x, new_y, new_z, dx, dy, dz, new_x, new_y, new_z);
}

inline bool GraphSearch::hasForced(int x, int y, int z, int dx, int dy, int dz)
{
  int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
  int id = (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1);
  switch (norm1)
  {
  case 1:
    // 1-d move, check 8 neighbors
    for (int fn = 0; fn < 8; ++fn)
    {
      int nx = x + jn3d_->f1[id][0][fn];
      int ny = y + jn3d_->f1[id][1][fn];
      int nz = z + jn3d_->f1[id][2][fn];
      if (isOccupied(nx, ny, nz))
        return true;
    }
    return false;
  case 2:
    // 2-d move, check 8 neighbors
    for (int fn = 0; fn < 8; ++fn)
    {
      int nx = x + jn3d_->f1[id][0][fn];
      int ny = y + jn3d_->f1[id][1][fn];
      int nz = z + jn3d_->f1[id][2][fn];
      if (isOccupied(nx, ny, nz))
        return true;
    }
    return false;
  case 3:
    // 3-d move, check 6 neighbors
    for (int fn = 0; fn < 6; ++fn)
    {
      int nx = x + jn3d_->f1[id][0][fn];
      int ny = y + jn3d_->f1[id][1][fn];
      int nz = z + jn3d_->f1[id][2][fn];
      if (isOccupied(nx, ny, nz))
        return true;
    }
    return false;
  default:
    return false;
  }
}

std::vector<StatePtr> GraphSearch::getPath() const
{
  return path_;
}

std::vector<StatePtr> GraphSearch::getOpenSet() const
{
  std::vector<StatePtr> ss;
  for (const auto &it : hm_)
  {
    if (it && it->opened && !it->closed)
      ss.push_back(it);
  }
  return ss;
}

std::vector<StatePtr> GraphSearch::getCloseSet() const
{
  std::vector<StatePtr> ss;
  for (const auto &it : hm_)
  {
    if (it && it->closed)
      ss.push_back(it);
  }
  return ss;
}

std::vector<StatePtr> GraphSearch::getAllSet() const
{
  std::vector<StatePtr> ss;
  for (const auto &it : hm_)
  {
    if (it)
      ss.push_back(it);
  }
  return ss;
}

constexpr int JPS2DNeib::nsz[3][2];

JPS2DNeib::JPS2DNeib()
{
  int id = 0;
  for (int dy = -1; dy <= 1; ++dy)
  {
    for (int dx = -1; dx <= 1; ++dx)
    {
      int norm1 = std::abs(dx) + std::abs(dy);
      for (int dev = 0; dev < nsz[norm1][0]; ++dev)
        Neib(dx, dy, norm1, dev, ns[id][0][dev], ns[id][1][dev]);
      for (int dev = 0; dev < nsz[norm1][1]; ++dev)
      {
        FNeib(dx, dy, norm1, dev,
              f1[id][0][dev], f1[id][1][dev],
              f2[id][0][dev], f2[id][1][dev]);
      }
      id++;
    }
  }
}

void JPS2DNeib::print()
{
  for (int dx = -1; dx <= 1; dx++)
  {
    for (int dy = -1; dy <= 1; dy++)
    {
      int id = (dx + 1) + 3 * (dy + 1);
      printf("[dx: %d, dy: %d]-->id: %d:\n", dx, dy, id);
      for (unsigned int i = 0; i < sizeof(f1[id][0]) / sizeof(f1[id][0][0]); i++)
        printf("                f1: [%d, %d]\n", f1[id][0][i], f1[id][1][i]);
    }
  }
}

void JPS2DNeib::Neib(int dx, int dy, int norm1, int dev, int &tx, int &ty)
{
  switch (norm1)
  {
  case 0:
    switch (dev)
    {
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
    switch (dev)
    {
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

void JPS2DNeib::FNeib(int dx, int dy, int norm1, int dev,
                      int &fx, int &fy, int &nx, int &ny)
{
  switch (norm1)
  {
  case 1:
    switch (dev)
    {
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
    if (dx == 0)
      fx = fy, fy = 0;

    nx = dx + fx;
    ny = dy + fy;
    return;
  case 2:
    switch (dev)
    {
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

JPS3DNeib::JPS3DNeib()
{
  int id = 0;
  for (int dz = -1; dz <= 1; ++dz)
  {
    for (int dy = -1; dy <= 1; ++dy)
    {
      for (int dx = -1; dx <= 1; ++dx)
      {
        int norm1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
        for (int dev = 0; dev < nsz[norm1][0]; ++dev)
          Neib(dx, dy, dz, norm1, dev,
               ns[id][0][dev], ns[id][1][dev], ns[id][2][dev]);
        for (int dev = 0; dev < nsz[norm1][1]; ++dev)
        {
          FNeib(dx, dy, dz, norm1, dev,
                f1[id][0][dev], f1[id][1][dev], f1[id][2][dev],
                f2[id][0][dev], f2[id][1][dev], f2[id][2][dev]);
        }
        id++;
      }
    }
  }
}

void JPS3DNeib::Neib(int dx, int dy, int dz, int norm1, int dev,
                     int &tx, int &ty, int &tz)
{
  switch (norm1)
  {
  case 0:
    switch (dev)
    {
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
    switch (dev)
    {
    case 0:
      if (dz == 0)
      {
        tx = 0;
        ty = dy;
        tz = 0;
        return;
      }
      else
      {
        tx = 0;
        ty = 0;
        tz = dz;
        return;
      }
    case 1:
      if (dx == 0)
      {
        tx = 0;
        ty = dy;
        tz = 0;
        return;
      }
      else
      {
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
    switch (dev)
    {
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

void JPS3DNeib::FNeib(int dx, int dy, int dz, int norm1, int dev,
                      int &fx, int &fy, int &fz,
                      int &nx, int &ny, int &nz)
{
  switch (norm1)
  {
  case 1:
    switch (dev)
    {
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
    if (dx != 0)
    {
      fz = fx;
      fx = 0;
      nz = fz;
      nx = dx;
    }
    if (dy != 0)
    {
      fz = fy;
      fy = 0;
      nz = fz;
      ny = dy;
    }
    return;
  case 2:
    if (dx == 0)
    {
      switch (dev)
      {
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
    }
    else if (dy == 0)
    {
      switch (dev)
      {
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
    }
    else
    { // dz==0
      switch (dev)
      {
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
  case 3:
    switch (dev)
    {
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