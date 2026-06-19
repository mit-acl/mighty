/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <hgp/hgp_planner.hpp>

#include <fstream>
#include <iomanip>

using namespace termcolor;

// Per-call HGP debug log. Opened lazily on first use; writes voxel start/goal,
// whether the goal was clamped (outside the map), the resulting path length,
// and the path's first / last waypoint. Used to diagnose backward-pointing
// global paths around x=10.
namespace {
std::ofstream& hgp_debug_log() {
  static std::ofstream s;
  if (!s.is_open()) {
    s.open("/tmp/mighty_hgp_debug.log", std::ios::out | std::ios::trunc);
    if (s.is_open()) {
      s << std::fixed << std::setprecision(4);
      s << "# event | start_w | goal_w | start_int | goal_int_orig | goal_int_clamped | dim | origin | path_size | path_front_w | path_back_w\n";
      s.flush();
    }
  }
  return s;
}
}  // namespace

HGPPlanner::HGPPlanner(std::string global_planner, bool verbose, double v_max, double a_max,
                       double j_max, int hgp_timeout_duration_ms, double w_unknown, double w_align,
                       double decay_len_cells, double w_side, int los_cells, double min_len,
                       double min_turn)
    : global_planner_(global_planner),
      planner_verbose_(verbose),
      v_max_(v_max),
      a_max_(a_max),
      j_max_(j_max),
      hgp_timeout_duration_ms_(hgp_timeout_duration_ms),
      w_unknown_(w_unknown),
      w_align_(w_align),
      decay_len_cells_(decay_len_cells),
      w_side_(w_side),
      los_cells_(los_cells),
      min_len_(min_len),
      min_turn_(min_turn) {
  if (planner_verbose_) printf(ANSI_COLOR_CYAN "HGP PLANNER VERBOSE ON\n" ANSI_COLOR_RESET);
}

// --- In hgp_planner.cpp (definitions) ---
bool HGPPlanner::lineOfSightCapsule(const Vecf<3>& a, const Vecf<3>& b,
                                    int inflate_radius_cells) const {
  // Sample along the segment and check a "capsule" of radius r around it.
  // r is in *cells* (so r=1.5 checks a segment thickened by 1.5 voxels).
  const double res = map_util_->getRes();
  const Vecf<3> d = b - a;
  const double L = d.norm();
  if (L < 1e-6) return true;

  const Vecf<3> u = d / L;
  const double step = std::max(0.25 * res, 0.05);  // fine step
  const int radial = std::max(1, inflate_radius_cells);

  for (double s = 0.0; s <= L; s += step) {
    const Vecf<3> c = a + s * u;

    // Check the centerline first.
    if (map_util_->isOccupied(c)) return false;

    // Then check a small cross-section around the center, approximating the capsule.
    for (int ix = -radial; ix <= radial; ++ix)
      for (int iy = -radial; iy <= radial; ++iy)
        for (int iz = -radial; iz <= radial; ++iz) {
          if (ix == 0 && iy == 0 && iz == 0) continue;
          Vecf<3> p = c + Vecf<3>(ix * res, iy * res, iz * res);
          // Skip far away corners outside the target radius
          if ((p - c).norm() > inflate_radius_cells * res) continue;
          if (map_util_->isOccupied(p)) return false;
        }
  }
  return true;
}

vec_Vecf<3> HGPPlanner::shortCutByLoS(const vec_Vecf<3>& in, int inflate_radius_cells) const {
  if (in.size() < 2) return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  size_t i = 0;
  out.push_back(in.front());

  while (i + 1 < in.size()) {
    size_t best = i + 1;
    // Greedily jump as far as LoS (with inflation) allows
    for (size_t j = i + 2; j < in.size(); ++j) {
      if (lineOfSightCapsule(in[i], in[j], inflate_radius_cells))
        best = j;
      else
        break;  // once blocked, further j will also be blocked in most maps
    }
    out.push_back(in[best]);
    i = best;
  }
  // Merge duplicates if any
  if (out.size() >= 2 && (out.end()[-1] - out.end()[-2]).norm() < 1e-9) out.pop_back();
  return out;
}

vec_Vecf<3> HGPPlanner::collapseShortEdges(const vec_Vecf<3>& in, double min_len) const {
  if (in.size() < 2) return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  out.push_back(in.front());
  for (size_t i = 1; i + 1 < in.size(); ++i) {
    if ((in[i] - out.back()).norm() >= min_len) out.push_back(in[i]);
  }
  out.push_back(in.back());
  return out;
}

static inline double angleDeg(const Vecf<3>& u, const Vecf<3>& v) {
  const double nu = u.norm(), nv = v.norm();
  if (nu < 1e-9 || nv < 1e-9) return 0.0;
  double c = std::max(-1.0, std::min(1.0, u.dot(v) / (nu * nv)));
  return std::acos(c) * 180.0 / M_PI;
}

vec_Vecf<3> HGPPlanner::angleSpacingFilter(const vec_Vecf<3>& in, double min_turn_deg,
                                           double min_seg_len) const {
  if (in.size() < 3) return in;
  vec_Vecf<3> out;
  out.reserve(in.size());
  out.push_back(in.front());

  for (size_t i = 1; i + 1 < in.size(); ++i) {
    const Vecf<3> u = in[i] - out.back();
    const Vecf<3> v = in[i + 1] - in[i];
    const bool short_edge = (in[i] - out.back()).norm() < min_seg_len;
    const bool tiny_turn = angleDeg(u, v) < min_turn_deg;

    // Drop in[i] if it's both a tiny turn AND too close to its predecessor.
    if (!(short_edge && tiny_turn)) out.push_back(in[i]);
  }
  out.push_back(in.back());
  return out;
}

void HGPPlanner::updateVmax(double v_max) { v_max_ = v_max; }

void HGPPlanner::setMapUtil(const std::shared_ptr<mighty::MapUtil<3>>& map_util) {
  // Deep copy the map_util
  map_util_ = std::make_shared<mighty::MapUtil<3>>(*map_util);
}

int HGPPlanner::status() { return status_; }

vec_Vecf<3> HGPPlanner::getPath() { return path_; }

vec_Vecf<3> HGPPlanner::getRawPath() { return raw_path_; }

vec_Vecf<3> HGPPlanner::removeCornerPts(const vec_Vecf<3>& path) {
  if (path.size() < 2) return path;

  // Cost-aware shortcutting:
  // - Collision feasibility still enforced by isBlocked(a,b).
  // - "Cost" matches global planner intent: length + w_heat * integral(heat) along segment.
  // - Optional guardrail on peak heat to prevent "cutting through" hot zones when w_heat is
  // moderate.

  const bool heat_on =
      (map_util_ && (map_util_->dynamicHeatEnabled() || map_util_->staticHeatEnabled()) &&
       map_util_->getHeatWeight() > 0.0f);
  const double w_heat = heat_on ? (double)map_util_->getHeatWeight() : 0.0;
  const double res = map_util_ ? (double)map_util_->getRes() : 0.1;

  // Sampling step for segment cost integration (tradeoff: accuracy vs runtime)
  const double ds = std::max(0.5 * res, 0.05);  // ~half-voxel, but not too tiny

  // Guardrail: allow some peak-heat increase, but not large.
  // Increase if you want more aggressive smoothing; decrease if you want stronger avoidance.
  const double peak_heat_relax = 0.10;  // 10%

  auto segCostAndPeakHeat = [&](const Vecf<3>& a, const Vecf<3>& b,
                                double* peak_heat_out) -> double {
    if (peak_heat_out) *peak_heat_out = 0.0;

    if (!map_util_ || !lineOfSightCapsule(a, b, los_cells_))
      return std::numeric_limits<decimal_t>::infinity();

    const Vecf<3> d = b - a;
    const double L = (double)d.norm();
    if (L < 1e-9) return 0.0;

    // Base geometric term
    double cost = L;

    // Heat integral term (only when enabled)
    if (heat_on) {
      const Vecf<3> u = d / (decimal_t)L;
      double heat_int = 0.0;
      double peak_h = 0.0;

      // Sample along [0, L]
      for (double s = 0.0; s <= L; s += ds) {
        const Vecf<3> p = a + (decimal_t)s * u;
        const Veci<3> pi = map_util_->floatToInt(p);
        const float h = map_util_->getHeat(pi(0), pi(1), pi(2));
        peak_h = std::max(peak_h, (double)h);
        heat_int += (double)h * ds;
      }

      // include endpoint exactly (optional but cheap)
      {
        const Veci<3> bi = map_util_->floatToInt(b);
        const float h = map_util_->getHeat(bi(0), bi(1), bi(2));
        if ((double)h > peak_h) peak_h = (double)h;
      }

      if (peak_heat_out) *peak_heat_out = peak_h;

      // Add weighted heat integral
      cost += w_heat * heat_int;
    }

    return cost;
  };

  vec_Vecf<3> optimized_path;
  optimized_path.reserve(path.size());

  Vecf<3> prev_pose = path.front();
  optimized_path.push_back(prev_pose);

  // cost1 represents cost(prev_pose -> path[i]) for the "current kept edge" in the logic.
  // Initialize with edge (0->1).
  double peak1 = 0.0, peak2 = 0.0, peak3 = 0.0;
  decimal_t cost1 = (decimal_t)segCostAndPeakHeat(path[0], path[1], &peak1);

  for (unsigned int i = 1; i < path.size() - 1; i++) {
    const Vecf<3> pose1 = path[i];
    const Vecf<3> pose2 = path[i + 1];

    decimal_t cost2 = (decimal_t)segCostAndPeakHeat(pose1, pose2, &peak2);
    decimal_t cost3 = (decimal_t)segCostAndPeakHeat(prev_pose, pose2, &peak3);

    bool accept_shortcut = false;

    if (cost3 < cost1 + cost2) {
      accept_shortcut = true;

      // Peak-heat guardrail: don't accept a shortcut that significantly increases peak heat
      // relative to the two edges it replaces.
      if (heat_on) {
        const double ref_peak = std::max(peak1, peak2);
        const double allowed = ref_peak * (1.0 + peak_heat_relax);

        if (peak3 > allowed) accept_shortcut = false;
      }
    }

    if (accept_shortcut) {
      // Keep prev_pose, skip pose1, and update the "current edge cost" to prev_pose->pose2.
      cost1 = cost3;
      peak1 = peak3;
    } else {
      // Keep pose1.
      optimized_path.push_back(pose1);

      // Now the "current kept point" advances.
      prev_pose = pose1;

      // Update cost1 for the new kept edge pose1->pose2 (note: peak1 should track that edge too).
      cost1 = (decimal_t)segCostAndPeakHeat(pose1, pose2, &peak1);
    }
  }

  optimized_path.push_back(path.back());
  return optimized_path;
}

vec_Vecf<3> HGPPlanner::removeLinePts(const vec_Vecf<3>& path) {
  if (path.size() < 3) return path;

  vec_Vecf<3> new_path;
  new_path.push_back(path.front());
  for (unsigned int i = 1; i < path.size() - 1; i++) {
    Vecf<3> p = (path[i + 1] - path[i]) - (path[i] - path[i - 1]);
    if (fabs(p(0)) + fabs(p(1)) + fabs(p(2)) > 1e-2) new_path.push_back(path[i]);
  }
  new_path.push_back(path.back());
  return new_path;
}

vec_Vecf<3> HGPPlanner::getOpenSet() const {
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getOpenSet();
  for (const auto& it : ss) {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

vec_Vecf<3> HGPPlanner::getCloseSet() const {
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getCloseSet();
  for (const auto& it : ss) {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

vec_Vecf<3> HGPPlanner::getAllSet() const {
  vec_Vecf<3> ps;
  const auto ss = graph_search_->getAllSet();
  for (const auto& it : ss) {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    ps.push_back(map_util_->intToFloat(pn));
  }
  return ps;
}

bool HGPPlanner::plan(const Vecf<3>& start, const Vecf<3>& start_vel, const Vecf<3>& goal,
                      double& final_g, double current_time, decimal_t eps) {
  if (map_util_->map_.size() == 0) {
    std::cout << "map size: " << map_util_->map_.size() << std::endl;
    printf(ANSI_COLOR_RED "need to set the map!\n" ANSI_COLOR_RESET);
    return false;
  }

  if (planner_verbose_) {
    std::cout << "Start: " << start.transpose() << std::endl;
    std::cout << "Goal:  " << goal.transpose() << std::endl;
    std::cout << "Epsilon:  " << eps << std::endl;
  }

  path_.clear();
  raw_path_.clear();
  status_ = 0;

  Veci<3> start_int = map_util_->floatToInt(start);

  // In 2D mode, validate against 2D map instead of 3D map (ground points would block start/goal)
  if (is_2d_mode_ && map_util_->has2DMap()) {
    // Check 2D occupancy for start
    const int sx = start_int(0), sy = start_int(1);
    const Veci<3> dim = map_util_->getDim();
    if (sx < 0 || sx >= dim(0) || sy < 0 || sy >= dim(1)) {
      std::cout << bold << red << "Start is outside 2D map" << reset << std::endl;
      status_ = 1;
      return false;
    }
    if (!map_util_->useSoftCostObstacles() && map_util_->get2DOccupancy(sx, sy) != 0) {
      std::cout << bold << red << "Start is occupied in 2D map" << reset << std::endl;
      status_ = 1;
      return false;
    }
  } else {
    if (map_util_->isOutside(start_int) ||
        (!map_util_->useSoftCostObstacles() && map_util_->isOccupied(start_int))) {
      if (planner_verbose_) {
        if (map_util_->isOccupied(start_int))
          printf(ANSI_COLOR_RED "start is occupied!\n" ANSI_COLOR_RESET);
        else if (map_util_->isUnknown(start_int))
          printf(ANSI_COLOR_RED "start is unknown!\n" ANSI_COLOR_RESET);
        else {
          printf(ANSI_COLOR_RED "start is outside!\n" ANSI_COLOR_RESET);
          std::cout << "start: " << start.transpose() << std::endl;
          std::cout << "start_int: " << start_int.transpose() << std::endl;
          std::cout << "Map origin: " << map_util_->getOrigin().transpose() << std::endl;
          std::cout << "Map dim: " << map_util_->getDim().transpose() << std::endl;
        }
      }
      status_ = 1;
      return false;
    }
  }

  Veci<3> goal_int = map_util_->floatToInt(goal);
  const Veci<3> dim = map_util_->getDim();

  // Save original (pre-clamp) goal voxel for debug logging.
  const Veci<3> goal_int_orig = goal_int;

  // If goal is outside the map, clamp to nearest boundary cell
  if (map_util_->isOutside(goal_int)) {
    for (int i = 0; i < 3; ++i) goal_int(i) = std::clamp(goal_int(i), 0, dim(i) - 1);
  }

  // In non-soft-cost mode, reject if goal cell is occupied
  if (is_2d_mode_ && map_util_->has2DMap()) {
    const int gx = goal_int(0), gy = goal_int(1);
    if (!map_util_->useSoftCostObstacles() && map_util_->get2DOccupancy(gx, gy) != 0) {
      std::cout << bold << red << "Goal is occupied in 2D map" << reset << std::endl;
      status_ = 2;
      return false;
    }
  } else if (!map_util_->useSoftCostObstacles() && map_util_->isOccupied(goal_int)) {
    std::cout << bold << red << "goal is occupied!"
              << " goal=" << goal.transpose() << " goal_int=" << goal_int.transpose() << reset
              << std::endl;
    status_ = 2;
    return false;
  }

  if ((map_util_->map_).empty()) {
    if (planner_verbose_) printf(ANSI_COLOR_RED "need to set the map!\n" ANSI_COLOR_RESET);
    return -1;
  }

  // compute initial g value (this is due to the fact that the actual initial position is not on the
  // grid)
  double initial_g = (start - map_util_->intToFloat(start_int)).norm();

  // In 2D mode, use the projected 2D map with zDim=1
  const int zDim_for_search = is_2d_mode_ ? 1 : dim(2);
  const int8_t* cMap_for_search = nullptr;
  if (is_2d_mode_ && map_util_->has2DMap()) {
    cMap_for_search = map_util_->get2DMapData();
  } else {
    cMap_for_search = (map_util_->map_).data();
  }

  graph_search_ = std::make_shared<mighty::GraphSearch>(cMap_for_search, map_util_, dim(0), dim(1),
                                                        zDim_for_search, eps, planner_verbose_,
                                                        global_planner_, w_unknown_);
  graph_search_->setStartAndGoal(start, goal);
  double max_values[3] = {v_max_, a_max_, j_max_};
  graph_search_->setBounds(max_values);

  // Pass ESDF grid for ground robot A* cost
  if (esdf_grid_ && is_2d_mode_) {
    graph_search_->setEsdfGrid(esdf_grid_, esdf_weight_astar_, esdf_d_safe_astar_);
  }

  // In 2D mode, force start/goal z to 0
  const int sz = is_2d_mode_ ? 0 : start_int(2);
  const int gz = is_2d_mode_ ? 0 : goal_int(2);

  // Run global plan module
  int max_expand = max_expand_;
  graph_search_->plan(start_int(0), start_int(1), sz, goal_int(0), goal_int(1), gz, initial_g,
                      global_planning_time_, hgp_static_jps_time_, hgp_check_path_time_,
                      hgp_dynamic_astar_time_, hgp_recover_path_time_, current_time, start_vel,
                      max_expand, hgp_timeout_duration_ms_);

  const auto path = graph_search_->getPath();

  if (path.size() < 1) {
    std::cout << ANSI_COLOR_RED "Cannot find a path from " << start.transpose() << " to "
              << goal.transpose() << " Abort!" ANSI_COLOR_RESET << std::endl;
    status_ = -1;
    return false;
  }

  // get the final g value
  final_g = path.front()->g;

  //**** raw path, s --> g
  vec_Vecf<3> ps;
  for (const auto& it : path) {
    Veci<3> pn;
    pn << it->x, it->y, it->z;
    Vecf<3> pt = map_util_->intToFloat(pn);
    if (is_2d_mode_) {
      // In 2D mode, project path to z=0 reference plane
      pt(2) = 0.0;
    }
    ps.push_back(pt);
  }

  raw_path_ = ps;
  std::reverse(std::begin(raw_path_), std::end(raw_path_));

  if (disable_all_smoothing_) {
    // No smoothing — use raw A* path directly
    path_ = raw_path_;
  } else if (skip_path_smoothing_) {
    // Heat-aware Laplacian smoothing (moves waypoints, doesn't remove them)
    path_ = smoothPathHeatAware(raw_path_, smooth_iterations_, smooth_alpha_);
  } else if (is_2d_mode_) {
    // 2D mode: sample points from raw A* path such that consecutive waypoints
    // never exceed max_dist_vertexes_2d_. No smoothing or LoS.
    if (raw_path_.size() >= 2 && max_dist_vertexes_2d_ > 0.0) {
      path_.clear();
      path_.push_back(raw_path_.front());
      double accum_dist = 0.0;
      for (size_t i = 1; i < raw_path_.size(); ++i) {
        accum_dist += (raw_path_[i] - raw_path_[i - 1]).norm();
        if (accum_dist >= max_dist_vertexes_2d_) {
          path_.push_back(raw_path_[i]);
          accum_dist = 0.0;
        }
      }
      // Always include the last point
      if ((path_.back() - raw_path_.back()).norm() > 1e-6) {
        path_.push_back(raw_path_.back());
      }
    } else {
      path_ = raw_path_;
    }
  } else if (global_planner_ == "sjps" || global_planner_ == "sastar") {
    // 1) coarse LoS shortcut with inflation
    path_ = shortCutByLoS(raw_path_, los_cells_);

    // 2) collapse very short edges
    path_ = collapseShortEdges(path_, min_len_);

    // 3) a light angle/spacing filter
    path_ = angleSpacingFilter(path_, min_turn_, min_len_);

    // 4) (optional) one more LoS pass to knit long spans
    path_ = shortCutByLoS(path_, los_cells_);

    // 5) clean up path
    cleanUpPath(path_);
  } else {
    auto tmp = raw_path_;
    cleanUpPath(tmp);
    path_ = tmp;
  }

  // Corridor-center corner snap (ground robot only). No-op if snap_enabled_
  // is false (UAV flow is never touched). Runs after all other smoothing so
  // the input already has collinear/corner redundancy removed.
  snapCornersToClearance(path_);

  // Debug log: dump full HGP plan call (input start/goal voxels, clamping,
  // raw_path size, processed path size, path endpoints, map metadata).
  if (auto& s = hgp_debug_log(); s.is_open()) {
    auto fmtv = [](const Vecf<3>& v) {
      std::ostringstream o;
      o << std::fixed << std::setprecision(3)
        << v.x() << "," << v.y() << "," << v.z();
      return o.str();
    };
    auto fmti = [](const Veci<3>& v) {
      std::ostringstream o;
      o << v.x() << "," << v.y() << "," << v.z();
      return o.str();
    };
    s << "HGP"
      << " start_w=" << fmtv(start)
      << " goal_w=" << fmtv(goal)
      << " start_int=" << fmti(start_int)
      << " goal_int_orig=" << fmti(goal_int_orig)
      << " goal_int_clamped=" << fmti(goal_int)
      << " dim=" << fmti(dim)
      << " origin=" << fmtv(map_util_->getOrigin())
      << " raw_n=" << raw_path_.size()
      << " path_n=" << path_.size();
    if (!path_.empty()) {
      s << " front=" << fmtv(path_.front())
        << " back=" << fmtv(path_.back());
    }
    s << "\n";
    s.flush();
  }

  return true;
}

void HGPPlanner::cleanUpPath(vec_Vecf<3>& path) {
  // 1) Remove perfectly collinear points (purely geometric, safe)
  path = removeLinePts(path);

  // 2) Cost-aware corner shortcutting (respects heat cost + peak-heat guardrail)
  path = removeCornerPts(path);
  std::reverse(std::begin(path), std::end(path));
  path = removeCornerPts(path);
  std::reverse(std::begin(path), std::end(path));
}

void HGPPlanner::snapCornersToClearance(vec_Vecf<3>& path) const {
  // Ground-robot post-processing: push each sharp-corner waypoint toward the
  // local ESDF max-clearance point via gradient ascent. Open-field corners
  // (already ≥ snap_clearance_threshold_m_ from any obstacle) are left alone.
  if (!snap_enabled_ || !esdf_grid_) return;
  if (path.size() < 3) return;

  const double cos_thresh = std::cos(snap_corner_angle_rad_);

  for (size_t i = 1; i + 1 < path.size(); ++i) {
    const Vecf<3>& p_prev = path[i - 1];
    const Vecf<3>& p_curr = path[i];
    const Vecf<3>& p_next = path[i + 1];

    // Direction-change check: compare unit vectors of the two incident edges
    // in the xy plane (path is 2D for ground robot). If the turn is shallower
    // than snap_corner_angle_rad_ it's not a corner — skip.
    Eigen::Vector2d v1(p_curr.x() - p_prev.x(), p_curr.y() - p_prev.y());
    Eigen::Vector2d v2(p_next.x() - p_curr.x(), p_next.y() - p_curr.y());
    const double n1 = v1.norm();
    const double n2 = v2.norm();
    if (n1 < 1e-6 || n2 < 1e-6) continue;
    const double cos_turn = (v1 / n1).dot(v2 / n2);
    if (cos_turn >= cos_thresh) continue;  // nearly straight, not a corner

    // Open-field guard: already ≥ threshold from obstacles, no push needed.
    if (!esdf_grid_->isInBounds(p_curr.x(), p_curr.y())) continue;
    const double d0 = esdf_grid_->queryDistance(p_curr.x(), p_curr.y());
    if (d0 >= snap_clearance_threshold_m_) continue;

    // Gradient ascent from the corner toward the local max-clearance point.
    // Cap by max_ascent_m and stop once clearance meets the threshold.
    double px = p_curr.x();
    double py = p_curr.y();
    double total_ascent = 0.0;

    for (int iter = 0; iter < snap_ascent_max_iters_; ++iter) {
      if (!esdf_grid_->isInBounds(px, py)) break;
      const Eigen::Vector2d g = esdf_grid_->queryGradient(px, py);
      const double g_norm = g.norm();
      if (g_norm < 1e-3) break;  // local max or flat region

      // Cap the step so we never exceed the total ascent budget.
      double step_len = snap_ascent_step_m_;
      if (total_ascent + step_len > snap_max_ascent_m_) {
        step_len = snap_max_ascent_m_ - total_ascent;
        if (step_len <= 0.0) break;
      }
      const Eigen::Vector2d step = (g / g_norm) * step_len;
      px += step.x();
      py += step.y();
      total_ascent += step_len;

      const double d = esdf_grid_->queryDistance(px, py);
      if (d >= snap_clearance_threshold_m_) break;  // reached corridor center
      if (total_ascent >= snap_max_ascent_m_) break;
    }

    // Commit the snapped xy; preserve z (operating in the 2D plane).
    path[i].x() = px;
    path[i].y() = py;
  }
}

vec_Vecf<3> HGPPlanner::smoothPathHeatAware(const vec_Vecf<3>& path, int iterations,
                                            double alpha) const {
  if (path.size() < 3 || !map_util_) return path;

  const bool heat_on = (map_util_->dynamicHeatEnabled() || map_util_->staticHeatEnabled()) &&
                       map_util_->getHeatWeight() > 0.0f;
  const bool has_2d = map_util_->has2DMap();

  vec_Vecf<3> smoothed = path;

  for (int iter = 0; iter < iterations; iter++) {
    for (size_t i = 1; i < smoothed.size() - 1; i++) {
      Vecf<3> mid = (smoothed[i - 1] + smoothed[i + 1]) / 2.0;
      Vecf<3> candidate = smoothed[i] + alpha * (mid - smoothed[i]);

      // Reject if occupied in 3D
      if (map_util_->isOccupied(candidate)) continue;

      // Reject if occupied in 2D projection
      if (has_2d) {
        Veci<3> ci = map_util_->floatToInt(candidate);
        if (map_util_->get2DOccupancy(ci(0), ci(1)) != 0) continue;
      }

      // Reject if heat increases too much (use 2D heat in 2D mode)
      if (heat_on) {
        Veci<3> oi = map_util_->floatToInt(smoothed[i]);
        Veci<3> ci = map_util_->floatToInt(candidate);
        float h_old = is_2d_mode_ ? map_util_->getHeat2D(oi(0), oi(1))
                                  : map_util_->getHeat(oi(0), oi(1), oi(2));
        float h_new = is_2d_mode_ ? map_util_->getHeat2D(ci(0), ci(1))
                                  : map_util_->getHeat(ci(0), ci(1), ci(2));
        if (h_new > h_old * 1.1 + 0.01) continue;
      }

      smoothed[i] = candidate;
    }
  }
  return smoothed;
}

double HGPPlanner::getInitialGuessPlanningTime() { return global_planning_time_; }

double HGPPlanner::getStaticJPSPlanningTime() { return hgp_static_jps_time_; }

double HGPPlanner::getCheckPathTime() { return hgp_check_path_time_; }

double HGPPlanner::getDynamicAstarTime() { return hgp_dynamic_astar_time_; }

double HGPPlanner::getRecoverPathTime() { return hgp_recover_path_time_; }
