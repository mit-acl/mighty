// Header-only geometry helpers for global-path post-processing.
//
// Kept header-only (and dependency-light: Eigen + hgp/data_type only) so the
// logic is unit-testable without pulling in ROS / the full planner. See
// docs/plans/persistent_global_map_plan.md (Phase 2).

#pragma once

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

#include "hgp/data_type.hpp"

namespace mighty_utils {

/**
 * @brief Truncate a global path at the planning horizon so a local optimizer
 *        only ever sees the next `horizon` meters of a (stable) route.
 *
 * Walks `path` from the end and finds the last vertex within `horizon` of
 * `center`, then interpolates the exact point where the path crosses the sphere
 * of radius `horizon` centered at `center`. `path` is resized to end at that
 * crossing (kept vertices [0..idx] + the crossing point), and the crossing is
 * written to `subgoal_out`.
 *
 * Matches the "last exit" semantics of getLastIntersectionWithSphere so a route
 * that briefly dips back inside the horizon still commits as far along as
 * possible. When the whole path lies within `horizon` (the goal is reachable
 * inside one horizon), `path` is left unchanged and `subgoal_out = path.back()`.
 *
 * @param path        In/out global path (world frame). Unchanged if it already
 *                    ends within the horizon.
 * @param center      Sphere center — normally the local planning start A.pos.
 * @param horizon     Sphere radius (meters). Must be > 0 to truncate.
 * @param subgoal_out Output: the horizon-boundary subgoal (== path.back() when
 *                    no truncation happens).
 * @return false if `path` has fewer than 2 vertices (nothing to truncate);
 *         true otherwise.
 */
inline bool truncateGlobalPathAtHorizon(vec_Vecf<3>& path,
                                        const Eigen::Vector3d& center,
                                        double horizon,
                                        Eigen::Vector3d& subgoal_out) {
  if (path.size() < 2) return false;

  // Last vertex inside the horizon sphere (scan from the end).
  int index = -1;
  for (int i = static_cast<int>(path.size()) - 1; i >= 0; --i) {
    if ((path[i] - center).norm() < horizon) {
      index = i;
      break;
    }
  }

  // Whole path already within the horizon, or the last vertex is the boundary:
  // keep the path intact and use its end as the subgoal (goal within reach).
  if (index < 0 || index == static_cast<int>(path.size()) - 1) {
    subgoal_out = path.back();
    return true;
  }

  // Interpolate the crossing on segment [A=path[index] (inside), B=path[index+1]
  // (outside)]: solve ||A + t*(B-A) - center||^2 = horizon^2 for t in [0,1],
  // taking the exit (larger) root.
  const Eigen::Vector3d A = path[index];
  const Eigen::Vector3d B = path[index + 1];
  const Eigen::Vector3d d = B - A;
  const Eigen::Vector3d f = A - center;
  const double a = d.dot(d);
  const double b = 2.0 * f.dot(d);
  const double c = f.dot(f) - horizon * horizon;
  double t = 1.0;
  if (a > 1e-12) {
    const double disc = b * b - 4.0 * a * c;
    if (disc >= 0.0) {
      t = (-b + std::sqrt(disc)) / (2.0 * a);
      t = std::max(0.0, std::min(1.0, t));
    }
  }
  subgoal_out = A + t * d;

  // Keep [0..index], then append the boundary crossing as the new endpoint.
  path.resize(index + 1);
  path.push_back(subgoal_out);
  return true;
}

/**
 * @brief Re-anchor a cached route to a new start pose.
 *
 * Finds the vertex of `path` nearest to `start`, drops everything before it,
 * and prepends `start` — so the returned path begins at the robot's current
 * position and continues along the committed route. Used by the global-path
 * reuse gate as the robot advances along a cached route.
 *
 * @param path          In/out route (world frame). Replaced with the re-anchored
 *                      route [start, nearest_vertex, ... , end].
 * @param start         New start pose (robot A.pos).
 * @param deviation_out Output: distance from `start` to the nearest route vertex
 *                      (an approximation of how far the robot has strayed from
 *                      the route; caller uses it to reject a too-far reuse).
 * @return false if `path` is empty or the re-anchored route has < 2 vertices.
 */
inline bool reanchorPathToStart(vec_Vecf<3>& path, const Eigen::Vector3d& start,
                                double& deviation_out) {
  if (path.empty()) return false;

  int best = 0;
  double best_d = (path[0] - start).norm();
  for (int i = 1; i < static_cast<int>(path.size()); ++i) {
    const double d = (path[i] - start).norm();
    if (d < best_d) {
      best_d = d;
      best = i;
    }
  }
  deviation_out = best_d;

  vec_Vecf<3> out;
  out.reserve(path.size() - best + 1);
  out.push_back(start);
  for (int i = best; i < static_cast<int>(path.size()); ++i) out.push_back(path[i]);
  path = std::move(out);
  return path.size() >= 2;
}

}  // namespace mighty_utils
