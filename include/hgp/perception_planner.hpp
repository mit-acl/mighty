// /* ----------------------------------------------------------------------------
//  * Perception-aware lattice A* (C++ port of the Python prototype).
//  *
//  * State  = (ix, iy, ith, moved) on a 2D grid with n_headings discrete headings.
//  *          `moved` = 1 iff the state was reached by a forward primitive (so the
//  *          pose ~L behind along the incoming heading was really occupied by the
//  *          robot). `moved` = 0 after an in-place turn or at the start.
//  *
//  * Belief = tristate grid (OccGrid2D): FREE / OCC / UNK. Planning is optimistic:
//  *          UNK is traversable iff the observation-coverage invariant holds.
//  *
//  * Sensor = bearing-dependent range window (Mid360FOV for the tilted rover, or a
//  *          simple annulus+FOV), occluded by known-OCC cells (UNK is optimistically
//  *          transparent).
//  *
//  * Hard invariant (per forward primitive p -> p'): every swept cell is known FREE,
//  * or UNK and predicted-visible from the primitive start pose p, or UNK and
//  * predicted-visible from a virtual previous pose p (-) k*step*dir(theta) (only if
//  * moved == 1). Any other case makes the edge INFEASIBLE (not merely expensive).
//  *
//  * This is a faithful transfer of the prototype; the belief source is the existing
//  * tristate OccGrid2D rather than a bespoke BeliefMap. Kept ROS-free and Eigen-free
//  * so it is buildable/testable in isolation; the HGP stack converts its pose output
//  * to vec_Vecf<3> at the wiring boundary.
//  * -------------------------------------------------------------------------- */
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "mighty/occ_grid_2d.hpp"

namespace hgp {

// ---------------------------------------------------------------------------
// Sensor models
// ---------------------------------------------------------------------------

/** @brief Bearing-dependent range window [r_lo, r_hi] about the heading. */
class SensorModel {
 public:
  virtual ~SensorModel() = default;
  /** @return true and sets (r_lo, r_hi) if `rel_bearing` (pose-relative, rad) is
   *  inside the FOV; false otherwise. */
  virtual bool rangeBounds(double rel_bearing, double& r_lo, double& r_hi) const = 0;
  virtual double rMax() const = 0;
  virtual double fovDeg() const = 0;
};

/** @brief Simple annulus [r_min, r_max] x azimuth wedge (prototype SensorModel). */
class AnnulusSensor : public SensorModel {
 public:
  AnnulusSensor(double r_min = 1.0, double r_max = 5.0, double fov_deg = 360.0)
      : r_min_(r_min), r_max_(r_max), fov_deg_(fov_deg) {}
  bool rangeBounds(double rel_bearing, double& r_lo, double& r_hi) const override;
  double rMax() const override { return r_max_; }
  double fovDeg() const override { return fov_deg_; }

 private:
  double r_min_, r_max_, fov_deg_;
};

/** @brief Tilted Livox Mid-360 ground-rover FOV (prototype Mid360FOV).
 *
 *  Near boundary = ground-start horn of the pitched elevation band; far boundary
 *  = elliptical trust envelope ((x-xc)/(F-xc))^2 + (y/W)^2 = 1 in the pose frame.
 *  Occlusion is handled separately by Visibility. Constants default to the values
 *  in the prototype and can be overridden from config. */
class Mid360FOV : public SensorModel {
 public:
  Mid360FOV(double tilt_deg = 20.0, double e_lo = -7.0, double e_hi = 52.0,
            double ground_start_fwd = 1.10, double F = 2.30, double W = 0.315,
            double xc = 0.40, double bin_deg = 0.25);
  bool rangeBounds(double rel_bearing, double& r_lo, double& r_hi) const override;
  double rMax() const override { return r_max_; }
  double fovDeg() const override { return fov_deg_; }

 private:
  double bin_deg_;
  int n_ = 0;
  std::vector<double> r_near_;  // +inf where the ground horn is absent
  std::vector<double> r_far_;
  double r_max_ = 0.0;
  double fov_deg_ = 0.0;
};

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

struct PerceptionParams {
  double res = 0.125;          ///< grid resolution [m] (should match belief res)
  int n_headings = 12;         ///< heading bins (12 -> 30 deg)
  double prim_len = 1.25;      ///< forward primitive length [m] (>= sensor r_min!)
  double robot_radius = 0.25;  ///< disc footprint for collision inflation [m]
  double turn_cost = 0.35;     ///< cost of one in-place heading step [m-equiv]
  double w_unknown = 0.20;     ///< soft cost per metre of unknown traversed
  double goal_tol = 0.55;      ///< goal position tolerance [m]
  bool use_coverage_rule = true;   ///< false -> naive plan-through-unknown baseline
  int back_projection_steps = 5;   ///< ladder depth for the virtual-previous-pose credit
  double back_projection_step = 0.25;  ///< ladder spacing [m]
  // Resource guards (mirror MIGHTY's grid A*): on hitting either limit, or on the
  // open set emptying before the goal, the search returns the partial path to the
  // best (closest-to-goal) node reached rather than failing outright.
  int max_expand = 10000;          ///< max A* expansions; <= 0 means unlimited
  int timeout_ms = 1000;           ///< wall-clock planning budget [ms]; <= 0 means none
};

// ---------------------------------------------------------------------------
// Motion primitives on the heading lattice
// ---------------------------------------------------------------------------

struct Primitive {
  enum Kind { FWD, ARC, TURN };
  Kind kind = FWD;
  int dth = 0;                    ///< heading change in bins
  int end_dx = 0, end_dy = 0;     ///< integer end cell offset
  std::vector<std::array<int, 2>> sweep;  ///< centerline swept cell offsets (dx, dy)
  double cost = 0.0;
  double end_dx_m = 0.0, end_dy_m = 0.0, dtheta = 0.0;  ///< exact end pose (viz)
};

/** @return per start-heading list of primitives (straight, gentle arcs, turns). */
std::vector<std::vector<Primitive>> buildPrimitives(const PerceptionParams& P);

// ---------------------------------------------------------------------------
// Visibility (predicted observation) with memoized ray casting
// ---------------------------------------------------------------------------

/** @brief Full-width (four 32-bit ints) ray-cache key -- no bit-packing, so
 *  negative/out-of-bounds coordinates and indices > 65535 can never alias. */
struct RayKey {
  int ix0, iy0, ix1, iy1;
  bool operator==(const RayKey& o) const {
    return ix0 == o.ix0 && iy0 == o.iy0 && ix1 == o.ix1 && iy1 == o.iy1;
  }
};
struct RayKeyHash {
  std::size_t operator()(const RayKey& k) const {
    std::size_t h = 1469598103934665603ull;  // FNV-1a
    auto mix = [&h](int v) {
      h = (h ^ static_cast<std::size_t>(static_cast<uint32_t>(v))) * 1099511628211ull;
    };
    mix(k.ix0);
    mix(k.iy0);
    mix(k.ix1);
    mix(k.iy1);
    return h;
  }
};

class Visibility {
 public:
  Visibility(const OccGrid2D& belief, const SensorModel& sensor)
      : m_(belief), s_(sensor) {}
  /** Predicted-visible test for cell (cix,ciy) from metric pose (px,py,pth). */
  bool visible(double px, double py, double pth, int cix, int ciy);

 private:
  bool rayClear(int ix0, int iy0, int ix1, int iy1);
  const OccGrid2D& m_;
  const SensorModel& s_;
  std::unordered_map<RayKey, char, RayKeyHash> ray_cache_;  // 0/1 = clear cached false/true
};

// ---------------------------------------------------------------------------
// Planner entry
// ---------------------------------------------------------------------------

/** @brief Why the search stopped -- lets the caller log/monitor (e.g. distinguish
 *  "hit the expansion cap" from "genuinely unreachable"). */
enum class StopReason { GOAL, MAX_EXPAND, TIMEOUT, EXHAUSTED, NO_PROGRESS };
const char* stopReasonStr(StopReason r);

struct PerceptionPlanResult {
  bool ok = false;
  bool partial = false;   ///< true if `states` is a best-node partial path (goal not reached)
  StopReason stop_reason = StopReason::EXHAUSTED;
  std::vector<std::array<double, 3>> states;  ///< dense pose path (x, y, theta)
  double cost = 0.0;
  int expanded = 0;
  int blind_unknown_entries = 0;              ///< unknown centerline cells the coverage audit
                                              ///< could not credit (ladder-based; diagnostic only)
  std::vector<std::array<double, 2>> blind_cells;  ///< their world coords (audit)
};

/** @brief Perception-aware lattice A*.
 *
 *  @param belief        tristate belief (FREE/OCC/UNK).
 *  @param sensor        sensor model used for the coverage invariant.
 *  @param P             planner parameters.
 *  @param start_x/y/theta  start pose (world metres / rad).
 *  @param goal_x/y      goal position (world metres).
 *  @param audit_sensor  optional sensor used only to audit realised coverage of
 *                       the returned path (defaults to `sensor`). */
PerceptionPlanResult planPerceptionAware(const OccGrid2D& belief, const SensorModel& sensor,
                                         const PerceptionParams& P, double start_x, double start_y,
                                         double start_theta, double goal_x, double goal_y,
                                         const SensorModel* audit_sensor = nullptr);

}  // namespace hgp
