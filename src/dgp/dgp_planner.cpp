/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <dgp/dgp_planner.hpp>

using namespace termcolor;

DGPPlanner::DGPPlanner(std::string global_planner, bool verbose, double v_max, double a_max, double j_max, int dgp_timeout_duration_ms, double w_unknown, double w_align, double decay_len_cells, double w_side, int los_cells, double min_len, double min_turn) : global_planner_(global_planner), planner_verbose_(verbose), v_max_(v_max), a_max_(a_max), j_max_(j_max), dgp_timeout_duration_ms_(dgp_timeout_duration_ms), w_unknown_(w_unknown), w_align_(w_align), decay_len_cells_(decay_len_cells), w_side_(w_side), los_cells_(los_cells), min_len_(min_len), min_turn_(min_turn)
{
  if (planner_verbose_)
    printf(ANSI_COLOR_CYAN "DGP PLANNER VERBOSE ON\n" ANSI_COLOR_RESET);
}

// --- In dgp_planner.cpp (definitions) ---
bool DGPPlanner::lineOfSightCapsule(const Vecf<3> &a,
                                    const Vecf<3> &b,
                                    int inflate_radius_cells) const
{
  // Sample along the segment and check a "capsule" of radius r around it.
  // r is in *cells* (so r=1.5 checks a segment thickened by 1.5 voxels).
  const double res = map_util_->getRes();
  const Vecf<3> d = b - a;
  const double L = d.norm();
  if (L < 1e-6)
    return true;

  const Vecf<3> u = d / L;
  const double step = std::max(0.25 * res, 0.05); // fine step
  const int radial = std::max(1, inflate_radius_cells);

  for (double s = 0.0; s <= L; s += step)
  {
    const Vecf<3> c = a + s * u;

    // Check the centerline first.
    if (map_util_->isOccupied(c))
      return false;

    // Then check a small cross-section around the center, approximating the capsule.
    for (int ix = -radial; ix <= radial; ++ix)
      for (int iy = -radial; iy <= radial; ++iy)
        for (int iz = -radial; iz <= radial; ++iz)
        {
          if (ix == 0 && iy == 0 && iz == 0)
            continue;
          Vecf<3> p = c + Vecf<3>(ix * res, iy * res, iz * res);
          // Skip far away corners outside the target radius
          if ((p - c).norm() > inflate_radius_cells * res)
            continue;
          if (map_util_->isOccupied(p))
            return false;
        }
  }
  return true;
}

vec_Vecf<3> DGPPlanner::shortCutByLoS(const vec_Vecf<3> &in,
                                      int inflate_radius_cells) const
{
  if (in.size() < 2)
    return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  size_t i = 0;
  out.push_back(in.front());

  while (i + 1 < in.size())
  {
    size_t best = i + 1;
    // Greedily jump as far as LoS (with inflation) allows
    for (size_t j = i + 2; j < in.size(); ++j)
    {
      if (lineOfSightCapsule(in[i], in[j], inflate_radius_cells))
        best = j;
      else
        break; // once blocked, further j will also be blocked in most maps
    }
    out.push_back(in[best]);
    i = best;
  }
  // Merge duplicates if any
  if (out.size() >= 2 && (out.end()[-1] - out.end()[-2]).norm() < 1e-9)
    out.pop_back();
  return out;
}

vec_Vecf<3> DGPPlanner::collapseShortEdges(const vec_Vecf<3> &in, double min_len) const
{
  if (in.size() < 2)
    return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  out.push_back(in.front());
  for (size_t i = 1; i + 1 < in.size(); ++i)
  {
    if ((in[i] - out.back()).norm() >= min_len)
      out.push_back(in[i]);
  }
  out.push_back(in.back());
  return out;
}

static inline double angleDeg(const Vecf<3> &u, const Vecf<3> &v)
{
  const double nu = u.norm(), nv = v.norm();
  if (nu < 1e-9 || nv < 1e-9)
    return 0.0;
  double c = std::max(-1.0, std::min(1.0, u.dot(v) / (nu * nv)));
  return std::acos(c) * 180.0 / M_PI;
}

vec_Vecf<3> DGPPlanner::angleSpacingFilter(const vec_Vecf<3> &in,
                                           double min_turn_deg,
                                           double min_seg_len) const
{
  if (in.size() < 3)
    return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  out.push_back(in.front());

  for (size_t i = 1; i + 1 < in.size(); ++i)
  {
    const Vecf<3> u = in[i] - out.back();
    const Vecf<3> v = in[i + 1] - in[i];
    const bool short_edge = (in[i] - out.back()).norm() < min_seg_len;
    const bool tiny_turn = angleDeg(u, v) < min_turn_deg;

    // Drop in[i] if it's both a tiny turn AND too close to its predecessor.
    if (!(short_edge && tiny_turn))
      out.push_back(in[i]);
  }
  out.push_back(in.back());
  return out;
}

void DGPPlanner::updateVmax(double v_max)
{
  v_max_ = v_max;
}

void DGPPlanner::setMapUtil(const std::shared_ptr<mighty::MapUtil<3>> &map_util)
{
  // Deep copy the map_util
  map_util_ = std::make_shared<mighty::MapUtil<3>>(*map_util);
  // map_util_ = map_util;
}

int DGPPlanner::status()
{
  return status_;
}

vec_Vecf<3> DGPPlanner::getPath()
{
  return path_;
}

vec_Vecf<3> DGPPlanner::getRawPath()
{
  return raw_path_;
}

vec_Vecf<3> DGPPlanner::removeCornerPts(const vec_Vecf<3> &path)
{
  if (path.size() < 2)
    return path;

  // cut zigzag segment
  vec_Vecf<3> optimized_path;
  Vecf<3> pose1 = path[0];
  Vecf<3> pose2 = path[1];
  Vecf<3> prev_pose = pose1;
  optimized_path.push_back(pose1);
  decimal_t cost1, cost2, cost3;

  if (!map_util_->isBlocked(pose1, pose2))
    cost1 = (pose1 - pose2).norm();
  else
    cost1 = std::numeric_limits<decimal_t>::infinity();

  for (unsigned int i = 1; i < path.size() - 1; i++)
  {
    pose1 = path[i];
    pose2 = path[i + 1];
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

vec_Vecf<3> DGPPlanner::removeLinePts(const vec_Vecf<3> &path)
{
  if (path.size() < 3)
    return path;

  vec_Vecf<3> new_path;
  new_path.push_back(path.front());
  for (unsigned int i = 1; i < path.size() - 1; i++)
  {
    Vecf<3> p = (path[i + 1] - path[i]) - (path[i] - path[i - 1]);
    if (fabs(p(0)) + fabs(p(1)) + fabs(p(2)) > 1e-2)
      new_path.push_back(path[i]);
  }
  new_path.push_back(path.back());
  return new_path;
}

vec_Vecf<3> DGPPlanner::getOpenSet() const
{
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getOpenSet();
  for (const auto &it : ss)
  {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

vec_Vecf<3> DGPPlanner::getCloseSet() const
{
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getCloseSet();
  for (const auto &it : ss)
  {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

vec_Vecf<3> DGPPlanner::getAllSet() const
{
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getAllSet();
  for (const auto &it : ss)
  {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

bool DGPPlanner::plan(const Vecf<3> &start, const Vecf<3> &start_vel, const Vecf<3> &goal, double &final_g, double current_time, decimal_t eps)
{

  if (map_util_->map_.size() == 0)
  {
    std::cout << "map size: " << map_util_->map_.size() << std::endl;
    printf(ANSI_COLOR_RED "need to set the map!\n" ANSI_COLOR_RESET);
    return false;
  }

  if (planner_verbose_)
  {
    std::cout << "Start: " << start.transpose() << std::endl;
    std::cout << "Goal:  " << goal.transpose() << std::endl;
    std::cout << "Epsilon:  " << eps << std::endl;
  }

  path_.clear();
  raw_path_.clear();
  status_ = 0;

  const Veci<3> start_int = map_util_->floatToInt(start);
  if (map_util_->isOutside(start_int) || map_util_->isOccupied(start_int))
  {
    if (planner_verbose_)
    {
      if (map_util_->isOccupied(start_int))
        printf(ANSI_COLOR_RED "start is occupied!\n" ANSI_COLOR_RESET);
      else if (map_util_->isUnknown(start_int))
        printf(ANSI_COLOR_RED "start is unknown!\n" ANSI_COLOR_RESET);
      else
      {
        printf(ANSI_COLOR_RED "start is outside!\n" ANSI_COLOR_RESET);
        std::cout << "start: " << start.transpose() << std::endl;
        std::cout << "start_int: " << start_int.transpose() << std::endl;
        std::cout << "Map origin: " << map_util_->getOrigin().transpose() << std::endl;
        std::cout << "Map dim: " << map_util_->getDim().transpose() << std::endl;
      }
    }
    status_ = 1;
    std::cout << bold << red << "Start is not free" << reset << std::endl;
    return false;
  }

  const Veci<3> goal_int = map_util_->floatToInt(goal);
  if (map_util_->isOutside(goal_int) || map_util_->isOccupied(goal_int))
  {
    if (planner_verbose_)
    {
      printf(ANSI_COLOR_RED "goal_int: %d %d %d\n" ANSI_COLOR_RESET, goal_int(0), goal_int(1), goal_int(2));
      printf(ANSI_COLOR_RED "is outside in x: %d\n" ANSI_COLOR_RESET, map_util_->isOutsideXYZ(goal_int, 0));
      printf(ANSI_COLOR_RED "is outside in y: %d\n" ANSI_COLOR_RESET, map_util_->isOutsideXYZ(goal_int, 1));
      printf(ANSI_COLOR_RED "is outside in z: %d\n" ANSI_COLOR_RESET, map_util_->isOutsideXYZ(goal_int, 2));
      std::cout << "Map origin: " << map_util_->getOrigin().transpose() << std::endl;
    }
    status_ = 2;
    std::cout << bold << red << "goal is not free!" << reset << std::endl;
    return false;
  }

  if ((map_util_->map_).empty())
  {
    if (planner_verbose_)
      printf(ANSI_COLOR_RED "need to set the map!\n" ANSI_COLOR_RESET);
    return -1;
  }

  const Veci<3> dim = map_util_->getDim();

  // compute initial g value (this is due to the fact that the actual initial position is not on the grid)
  double initial_g = (start - map_util_->intToFloat(start_int)).norm();

  // should we initialize the planner in constructor?
  graph_search_ = std::make_shared<mighty::GraphSearch>((map_util_->map_).data(), map_util_, dim(0), dim(1), dim(2), eps, planner_verbose_, global_planner_, w_unknown_);
  graph_search_->setStartAndGoal(start, goal);
  double max_values[3] = {v_max_, a_max_, j_max_};
  graph_search_->setBounds(max_values);

  // Run Initial guess module
  int max_expand = 10000;
  graph_search_->plan(start_int(0), start_int(1), start_int(2), goal_int(0), goal_int(1), goal_int(2), initial_g, global_planning_time_, dgp_static_jps_time_, dgp_check_path_time_, dgp_dynamic_astar_time_, dgp_recover_path_time_, current_time, start_vel, max_expand, dgp_timeout_duration_ms_);

  const auto path = graph_search_->getPath();

  if (path.size() < 1)
  {
    std::cout << ANSI_COLOR_RED "Cannot find a path from " << start.transpose() << " to " << goal.transpose() << " Abort!" ANSI_COLOR_RESET << std::endl;
    status_ = -1;
    return false;
  }

  // print the path
  // std::cout << "initial guess path found" << std::endl;
  // for (int i = 0; i < path.size(); i++)
  // {
  //   Vecf<3> p = map_util_->intToFloat(Veci<3>(path[i]->x, path[i]->y, path[i]->z));
  //   printf("path[%d]: %f %f %f and g: %f\n", i, p(0), p(1), p(2), path[i]->g);
  // }

  // get the final g value
  final_g = path.front()->g;

  //**** raw path, s --> g
  vec_Vecf<3> ps;
  for (const auto &it : path)
  {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }

  raw_path_ = ps;
  std::reverse(std::begin(raw_path_), std::end(raw_path_));

  // // Check if the path has really close points and if so remove them
  // if (path_.size() >= 2)
  // {
  //   vec_Vecf<3> new_path;
  //   new_path.push_back(path_.front());
  //   for (unsigned int i = 1; i < path_.size(); i++)
  //   {
  //     if ((path_[i] - new_path.back()).norm() > 0.3)
  //       new_path.push_back(path_[i]);
  //   }
  //   path_ = new_path;
  // }

  // 1) coarse LoS shortcut with inflation
  path_ = shortCutByLoS(raw_path_, los_cells_);

  // 2) collapse very short edges
  path_ = collapseShortEdges(path_, min_len_);

  // 3) a light angle/spacing filter
  path_ = angleSpacingFilter(path_, min_turn_, min_len_);

  // 4) (optional) one more LoS pass to knit long spans
  path_ = shortCutByLoS(path_, los_cells_);

  // // // 5) clean up path
  // cleanUpPath(path_);

  return true;
}

void DGPPlanner::cleanUpPath(vec_Vecf<3> &path)
{

  // Simplify the raw path
  path = removeLinePts(path);

  // If you activate removeCornerPts, the path will be too simplified for dynamic A*
  path = removeCornerPts(path);
  std::reverse(std::begin(path), std::end(path));
  path = removeCornerPts(path);
  std::reverse(std::begin(path), std::end(path));
}

double DGPPlanner::getInitialGuessPlanningTime()
{
  return global_planning_time_;
}

double DGPPlanner::getStaticJPSPlanningTime()
{
  return dgp_static_jps_time_;
}

double DGPPlanner::getCheckPathTime()
{
  return dgp_check_path_time_;
}

double DGPPlanner::getDynamicAstarTime()
{
  return dgp_dynamic_astar_time_;
}

double DGPPlanner::getRecoverPathTime()
{
  return dgp_recover_path_time_;
}
