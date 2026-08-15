// /* ----------------------------------------------------------------------------
//  * Perception-aware lattice A* -- C++ port of the Python prototype.
//  * See include/hgp/perception_planner.hpp for the model/invariant description.
//  * -------------------------------------------------------------------------- */
#include "hgp/perception_planner.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

namespace hgp {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kInf = std::numeric_limits<double>::infinity();

inline double wrapPi(double a) {
  // wrap to (-pi, pi]
  double r = std::fmod(a + kPi, 2.0 * kPi);
  if (r < 0.0) r += 2.0 * kPi;
  return r - kPi;
}
}  // namespace

const char* stopReasonStr(StopReason r) {
  switch (r) {
    case StopReason::GOAL: return "GOAL";
    case StopReason::MAX_EXPAND: return "MAX_EXPAND";
    case StopReason::TIMEOUT: return "TIMEOUT";
    case StopReason::EXHAUSTED: return "EXHAUSTED";
    case StopReason::NO_PROGRESS: return "NO_PROGRESS";
  }
  return "?";
}

// ===========================================================================
// AnnulusSensor
// ===========================================================================
bool AnnulusSensor::rangeBounds(double rel_bearing, double& r_lo, double& r_hi) const {
  if (fov_deg_ < 359.9) {
    const double half = (fov_deg_ * kPi / 180.0) / 2.0;
    if (std::abs(wrapPi(rel_bearing)) > half) return false;
  }
  r_lo = r_min_;
  r_hi = r_max_;
  return true;
}

// ===========================================================================
// Mid360FOV -- faithful port of the prototype precompute
// ===========================================================================
Mid360FOV::Mid360FOV(double tilt_deg, double e_lo, double e_hi, double ground_start_fwd, double F,
                     double W, double xc, double bin_deg)
    : bin_deg_(bin_deg) {
  const double tau = tilt_deg * kPi / 180.0;
  const double h = ground_start_fwd * std::tan((-e_lo) * kPi / 180.0 + tau);

  // Azimuth bins: edges = arange(-180, 180 + bin_deg, bin_deg) -> n = round(360/bin_deg).
  n_ = static_cast<int>(std::llround(360.0 / bin_deg_));
  r_near_.assign(n_, kInf);
  r_far_.assign(n_, kInf);

  // --- near boundary: min ground-hit range per azimuth bin over the elevation band ---
  const double phi_lo = -180.0, phi_hi = 180.0, phi_step = 0.1;  // deg, [lo, hi)
  const double ele_step = 0.2;                                   // deg, [e_lo, e_hi]
  const int n_ele = static_cast<int>(std::floor((e_hi - e_lo) / ele_step + 1e-9)) + 1;
  const int n_phi = static_cast<int>(std::llround((phi_hi - phi_lo) / phi_step));  // 3600

  for (int ie = 0; ie < n_ele; ++ie) {
    const double ele = (e_lo + ie * ele_step) * kPi / 180.0;
    const double ce = std::cos(ele), se = std::sin(ele);
    for (int ip = 0; ip < n_phi; ++ip) {
      const double phi = (phi_lo + ip * phi_step) * kPi / 180.0;
      const double cp = std::cos(phi), sp = std::sin(phi);
      const double dxw = ce * cp * std::cos(tau) + se * std::sin(tau);
      const double dyw = ce * sp;
      const double dzw = -ce * cp * std::sin(tau) + se * std::cos(tau);
      if (dzw >= -1e-9) continue;  // ray does not hit the ground plane below
      const double t = h / std::max(-dzw, 1e-12);
      const double gx = t * dxw, gy = t * dyw;
      const double brg_deg = std::atan2(gy, gx) * 180.0 / kPi;  // [-180, 180]
      const double r = std::hypot(gx, gy);
      // digitize(brg, edges) - 1, clipped to [0, n_-1]
      int idx = static_cast<int>(std::floor((brg_deg - phi_lo) / bin_deg_));
      if (idx < 0) idx = 0;
      if (idx > n_ - 1) idx = n_ - 1;
      if (r < r_near_[idx]) r_near_[idx] = r;
    }
  }

  // --- far boundary: elliptical trust envelope solved per bin-center bearing ---
  const double a = F - xc, b = W;
  const double A = 1.0 / (a * a), B = 1.0 / (b * b), C = -2.0 * xc / (a * a),
               D = xc * xc / (a * a) - 1.0;
  for (int i = 0; i < n_; ++i) {
    const double bc = (phi_lo + (i + 0.5) * bin_deg_) * kPi / 180.0;  // bin center, rad
    const double cb = std::cos(bc), sb = std::sin(bc);
    const double qa = A * cb * cb + B * sb * sb;
    const double qb = C * cb;
    const double disc = qb * qb - 4.0 * qa * D;
    if (qa > 1e-12 && disc >= 0.0) {
      r_far_[i] = (-qb + std::sqrt(disc)) / (2.0 * qa);
    } else {
      r_far_[i] = kInf;
    }
  }

  // --- derived: r_max and fov over bins where both boundaries are valid ---
  double rmax = 0.0, absmax_deg = 0.0;
  for (int i = 0; i < n_; ++i) {
    const bool okb = std::isfinite(r_near_[i]) && std::isfinite(r_far_[i]) && (r_far_[i] > r_near_[i]);
    if (!okb) continue;
    rmax = std::max(rmax, r_far_[i]);
    const double bc_deg = phi_lo + (i + 0.5) * bin_deg_;
    absmax_deg = std::max(absmax_deg, std::abs(bc_deg));
  }
  r_max_ = rmax;
  fov_deg_ = 2.0 * absmax_deg;
}

bool Mid360FOV::rangeBounds(double rel_bearing, double& r_lo, double& r_hi) const {
  const double b_deg = wrapPi(rel_bearing) * 180.0 / kPi;  // [-180, 180)
  int i = static_cast<int>((b_deg + 180.0) / bin_deg_);
  if (i < 0) i = 0;
  if (i > n_ - 1) i = n_ - 1;
  const double lo = r_near_[i], hi = r_far_[i];
  if (!std::isfinite(lo) || hi <= lo) return false;
  r_lo = lo;
  r_hi = hi;
  return true;
}

// ===========================================================================
// Motion primitives
// ===========================================================================
std::vector<std::vector<Primitive>> buildPrimitives(const PerceptionParams& P) {
  const int NH = P.n_headings;
  const double dth_rad = 2.0 * kPi / NH;
  std::vector<std::vector<Primitive>> prims(NH);
  const int n_samp = std::max(8, static_cast<int>(std::llround(P.prim_len / (P.res * 0.4))));

  for (int hh = 0; hh < NH; ++hh) {
    const double th0 = hh * dth_rad;

    for (int dth : {0, -1, 1}) {
      double ex, ey;
      if (dth == 0) {
        ex = P.prim_len * std::cos(th0);
        ey = P.prim_len * std::sin(th0);
      } else {
        const double dpsi = dth * dth_rad;
        const double R = P.prim_len / std::abs(dpsi);
        const double sign = dth > 0 ? 1.0 : -1.0;
        const double cx = -R * std::sin(th0) * sign;
        const double cy = R * std::cos(th0) * sign;
        const double psi = sign * P.prim_len / R;
        ex = cx + R * std::sin(th0 + psi) * sign;
        ey = cy - R * std::cos(th0 + psi) * sign;
      }
      const int end_dx = static_cast<int>(std::floor(ex / P.res + 0.5));
      const int end_dy = static_cast<int>(std::floor(ey / P.res + 0.5));
      const double Px = end_dx * P.res, Py = end_dy * P.res;
      const double c = std::hypot(Px, Py);
      const double alpha = wrapPi(std::atan2(Py, Px) - th0);

      Primitive pr;
      pr.kind = (dth == 0) ? Primitive::FWD : Primitive::ARC;
      pr.dth = dth;
      pr.end_dx = end_dx;
      pr.end_dy = end_dy;

      std::vector<std::array<int, 2>> sweep;
      std::unordered_map<uint64_t, char> seen;
      auto push_cell = [&](int cxi, int cyi) {
        const uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(cxi)) << 32) |
                             static_cast<uint32_t>(cyi);
        if (seen.find(key) == seen.end()) {
          seen.emplace(key, 1);
          sweep.push_back({cxi, cyi});
        }
      };

      double length;
      if (std::abs(alpha) < 1e-9) {
        length = c;
        for (int k = 1; k <= n_samp; ++k) {
          const double u = static_cast<double>(k) / n_samp;
          const double x = u * Px, y = u * Py;
          push_cell(static_cast<int>(std::floor(x / P.res + 0.5)),
                    static_cast<int>(std::floor(y / P.res + 0.5)));
        }
      } else {
        const double Rt = c / (2.0 * std::sin(alpha));  // signed radius
        length = std::abs(c * alpha / std::sin(alpha));
        for (int k = 1; k <= n_samp; ++k) {
          const double u = static_cast<double>(k) / n_samp;
          const double th_u = th0 + 2.0 * alpha * u;
          const double x = Rt * (std::sin(th_u) - std::sin(th0));
          const double y = -Rt * (std::cos(th_u) - std::cos(th0));
          push_cell(static_cast<int>(std::floor(x / P.res + 0.5)),
                    static_cast<int>(std::floor(y / P.res + 0.5)));
        }
      }

      pr.sweep = std::move(sweep);
      pr.cost = length + 5e-4 * std::abs(dth);  // epsilon tie-break per heading change
      pr.end_dx_m = Px;
      pr.end_dy_m = Py;
      pr.dtheta = dth * dth_rad;
      prims[hh].push_back(std::move(pr));
    }

    for (int dth : {-1, 1}) {
      Primitive pr;
      pr.kind = Primitive::TURN;
      pr.dth = dth;
      pr.cost = P.turn_cost;
      pr.dtheta = dth * dth_rad;
      prims[hh].push_back(std::move(pr));
    }
  }
  return prims;
}

// ===========================================================================
// Visibility
// ===========================================================================
bool Visibility::rayClear(int ix0, int iy0, int ix1, int iy1) {
  const RayKey key{ix0, iy0, ix1, iy1};
  auto it = ray_cache_.find(key);
  if (it != ray_cache_.end()) return it->second != 0;

  // Amanatides & Woo grid DDA over the segment between cell centers: visits every
  // cell the ray actually crosses (no gaps, unlike fractional sampling). Cell i
  // spans [i, i+1) with center i+0.5. Both endpoints (observer cell and target
  // cell) are excluded; UNK is transparent (optimism), OOB counts as occupied.
  bool ok = true;
  if (!(ix0 == ix1 && iy0 == iy1)) {
    const double x0 = ix0 + 0.5, y0 = iy0 + 0.5;
    const double x1 = ix1 + 0.5, y1 = iy1 + 0.5;
    const double dx = x1 - x0, dy = y1 - y0;

    const int stepx = (dx > 0) - (dx < 0);
    const int stepy = (dy > 0) - (dy < 0);
    // t (in [0,1]) to reach the first cell boundary, and to advance one full cell.
    double tMaxX = kInf, tDeltaX = kInf;
    if (stepx != 0) {
      const double bx = (stepx > 0) ? (ix0 + 1) : ix0;  // next x grid line
      tMaxX = (bx - x0) / dx;
      tDeltaX = std::abs(1.0 / dx);
    }
    double tMaxY = kInf, tDeltaY = kInf;
    if (stepy != 0) {
      const double by = (stepy > 0) ? (iy0 + 1) : iy0;
      tMaxY = (by - y0) / dy;
      tDeltaY = std::abs(1.0 / dy);
    }

    int cx = ix0, cy = iy0;
    while (true) {
      if (tMaxX < tMaxY) {
        cx += stepx;
        tMaxX += tDeltaX;
      } else {
        cy += stepy;
        tMaxY += tDeltaY;
      }
      if (cx == ix1 && cy == iy1) break;         // reached target cell (endpoint)
      if (tMaxX > 1.0 && tMaxY > 1.0) break;      // passed the segment end
      if (m_.isOccupied(cx, cy)) {                // interior cell blocks the ray
        ok = false;
        break;
      }
    }
  }
  ray_cache_.emplace(key, ok ? 1 : 0);
  return ok;
}

bool Visibility::visible(double px, double py, double pth, int cix, int ciy) {
  double cx, cy;
  m_.gridToWorld(cix, ciy, cx, cy);
  const double dx = cx - px, dy = cy - py;
  const double d = std::hypot(dx, dy);
  double r_lo, r_hi;
  if (!s_.rangeBounds(std::atan2(dy, dx) - pth, r_lo, r_hi)) return false;
  if (d < r_lo || d > r_hi) return false;
  int pix, piy;
  m_.worldToGrid(px, py, pix, piy);
  return rayClear(pix, piy, cix, ciy);
}

// ===========================================================================
// A* search
// ===========================================================================
namespace {

// Open-list entry: min-heap on f, then g (matches the prototype's (f, g, state)).
struct OpenNode {
  double f;
  double g;
  int64_t key;
  bool operator>(const OpenNode& o) const {
    if (f != o.f) return f > o.f;
    return g > o.g;
  }
};

}  // namespace

PerceptionPlanResult planPerceptionAware(const OccGrid2D& belief, const SensorModel& sensor,
                                         const PerceptionParams& P, double start_x, double start_y,
                                         double start_theta, double goal_x, double goal_y,
                                         const SensorModel* audit_sensor) {
  PerceptionPlanResult res;
  const int NH = P.n_headings;
  const double dth_rad = 2.0 * kPi / NH;
  const int W = belief.width();
  const int H = belief.height();
  if (W <= 0 || H <= 0) return res;

  const auto prims = buildPrimitives(P);
  Visibility vis(belief, sensor);
  Visibility vis_audit(belief, audit_sensor ? *audit_sensor : sensor);

  // --- footprint inflation: cells blocked for the robot center (disc) + borders ---
  const int rad = static_cast<int>(std::ceil(P.robot_radius / P.res));
  std::vector<char> blocked(static_cast<size_t>(W) * H, 0);
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      if (!belief.isOccupied(x, y)) continue;
      for (int dy = -rad; dy <= rad; ++dy) {
        for (int dx = -rad; dx <= rad; ++dx) {
          if (dx * dx + dy * dy > rad * rad) continue;
          const int nx = x + dx, ny = y + dy;
          if (nx < 0 || nx >= W || ny < 0 || ny >= H) continue;
          blocked[static_cast<size_t>(ny) * W + nx] = 1;
        }
      }
    }
  }
  // Borders blocked (footprint would leave the map).
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (x <= rad || x >= W - rad - 1 || y <= rad || y >= H - rad - 1)
        blocked[static_cast<size_t>(y) * W + x] = 1;

  auto inBounds = [&](int ix, int iy) { return ix >= 0 && ix < W && iy >= 0 && iy < H; };
  auto isBlocked = [&](int ix, int iy) {
    return !inBounds(ix, iy) || blocked[static_cast<size_t>(iy) * W + ix] != 0;
  };

  int six, siy, gix, giy;
  belief.worldToGrid(start_x, start_y, six, siy);
  belief.worldToGrid(goal_x, goal_y, gix, giy);
  (void)gix;
  (void)giy;
  int sih = static_cast<int>(std::llround(start_theta / dth_rad)) % NH;
  if (sih < 0) sih += NH;

  auto keyOf = [&](int ix, int iy, int ih, int mv) -> int64_t {
    return (((static_cast<int64_t>(iy) * W + ix) * NH + ih) * 2) + mv;
  };
  auto hFn = [&](int ix, int iy) {
    double wx, wy;
    belief.gridToWorld(ix, iy, wx, wy);
    return std::hypot(goal_x - wx, goal_y - wy);
  };

  struct Came {
    int64_t prev;
    const Primitive* prim;  // nullptr at the start
  };
  std::unordered_map<int64_t, double> g;
  std::unordered_map<int64_t, Came> came;
  std::unordered_map<int64_t, std::array<int, 4>> unpack;  // key -> (ix,iy,ih,mv)

  const int64_t start_key = keyOf(six, siy, sih, 0);
  g[start_key] = 0.0;
  came[start_key] = {-1, nullptr};
  unpack[start_key] = {six, siy, sih, 0};

  // Best (closest-to-goal) node seen, for partial-path recovery on a resource
  // limit or an exhausted open set -- mirrors MIGHTY's grid A* best_node behavior.
  const double start_h = hFn(six, siy);
  int64_t best_key = start_key;
  double best_h = start_h;

  std::priority_queue<OpenNode, std::vector<OpenNode>, std::greater<OpenNode>> open;
  open.push({start_h, 0.0, start_key});

  int64_t goal_key = -1;
  int expanded = 0;
  StopReason stop = StopReason::EXHAUSTED;  // default: open set emptied
  const auto t_start = std::chrono::steady_clock::now();

  while (!open.empty()) {
    const OpenNode top = open.top();
    open.pop();
    const int64_t st = top.key;
    auto git = g.find(st);
    if (git == g.end() || top.g > git->second + 1e-9) continue;
    const auto s = unpack[st];
    const int ix = s[0], iy = s[1], ih = s[2], mv = s[3];
    ++expanded;

    // Resource guards (match MIGHTY's grid A*): on hitting either limit, stop and
    // recover the best partial path found so far (handled after the loop).
    if (P.max_expand > 0 && expanded >= P.max_expand) {
      stop = StopReason::MAX_EXPAND;
      break;
    }
    if (P.timeout_ms > 0 &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - t_start)
                .count() >= P.timeout_ms) {
      stop = StopReason::TIMEOUT;
      break;
    }

    double wx, wy;
    belief.gridToWorld(ix, iy, wx, wy);
    const double h_cur = std::hypot(goal_x - wx, goal_y - wy);
    if (h_cur < best_h) {  // track closest-to-goal node for partial recovery
      best_h = h_cur;
      best_key = st;
    }
    if (h_cur <= P.goal_tol) {
      goal_key = st;
      stop = StopReason::GOAL;
      break;
    }

    const double px = wx, py = wy;
    const double pth = ih * dth_rad;
    const double gc = top.g;

    for (const Primitive& pr : prims[ih]) {
      if (pr.kind == Primitive::TURN) {
        const int nih = ((ih + pr.dth) % NH + NH) % NH;
        const int64_t nst = keyOf(ix, iy, nih, 0);
        const double ng = gc + pr.cost;
        auto it = g.find(nst);
        if (it == g.end() || ng < it->second - 1e-9) {
          g[nst] = ng;
          came[nst] = {st, &pr};
          unpack[nst] = {ix, iy, nih, 0};
          open.push({ng + hFn(ix, iy), ng, nst});
        }
        continue;
      }

      // forward / arc primitive: collision + coverage over the CENTERLINE swept
      // cells (matches the Python prototype invariant). OCC is handled for the
      // whole disc footprint by the robot_radius map inflation
      // (centerline-vs-inflated == footprint-vs-raw occupied); each swept UNKNOWN
      // cell must be predicted-observable before traversal.
      bool feasible = true;
      int unk_cells = 0;
      for (const auto& off : pr.sweep) {
        const int cix = ix + off[0], ciy = iy + off[1];
        if (isBlocked(cix, ciy)) {
          feasible = false;
          break;
        }
        if (!belief.isUnknown(cix, ciy)) continue;
        ++unk_cells;
        if (!P.use_coverage_rule) continue;
        bool ok = vis.visible(px, py, pth, cix, ciy);
        if (!ok && mv) {
          // Back-projected virtual poses along the incoming heading. Armed after
          // any primitive that reached this state with mv==1 (see mv_after), i.e.
          // whenever the incoming lattice step is fine enough (<= 45 deg) for the
          // straight-history assumption behind this ladder to hold.
          for (int k = 1; k <= P.back_projection_steps; ++k) {
            const double dlt = P.back_projection_step * k;
            const double bx = px - dlt * std::cos(pth);
            const double by = py - dlt * std::sin(pth);
            if (vis.visible(bx, by, pth, cix, ciy)) {
              ok = true;
              break;
            }
          }
        }
        if (!ok) {
          feasible = false;
          break;
        }
      }
      if (!feasible) continue;

      const int nih = ((ih + pr.dth) % NH + NH) % NH;
      // Faithful to the Python prototype: the back-projection ladder assumes a
      // ~straight incoming history along the current heading, which holds within one
      // heading bin. Straights always arm it; arcs arm it too as long as the lattice
      // step is <= 45 deg (coarser bins would break the straight-history assumption).
      const int mv_after = (pr.dth == 0 || dth_rad <= kPi / 4.0 + 1e-9) ? 1 : 0;
      const int nix = ix + pr.end_dx, niy = iy + pr.end_dy;
      if (!inBounds(nix, niy)) continue;
      const int64_t nst = keyOf(nix, niy, nih, mv_after);
      const double ng = gc + pr.cost + P.w_unknown * unk_cells * P.res;
      auto it = g.find(nst);
      if (it == g.end() || ng < it->second - 1e-9) {
        g[nst] = ng;
        came[nst] = {st, &pr};
        unpack[nst] = {nix, niy, nih, mv_after};
        open.push({ng + hFn(nix, niy), ng, nst});
      }
    }
  }

  res.expanded = expanded;

  // Terminal node: the goal if reached, else the best (closest-to-goal) node for a
  // partial path (MIGHTY best_node behavior). If nothing improved on the start
  // (open set exhausted / limit hit with zero progress), report failure so the
  // caller can fall back to the guarded grid search.
  int64_t terminal = goal_key;
  bool partial = false;
  if (terminal < 0) {
    if (best_key == start_key) {
      res.stop_reason = StopReason::NO_PROGRESS;
      res.ok = false;
      return res;
    }
    terminal = best_key;
    partial = true;
  }
  res.stop_reason = (goal_key >= 0) ? StopReason::GOAL : stop;  // MAX_EXPAND/TIMEOUT/EXHAUSTED

  // --- reconstruct pose path from the terminal node ---
  std::vector<int64_t> chain;
  for (int64_t st = terminal; st != -1; st = came[st].prev) chain.push_back(st);
  std::reverse(chain.begin(), chain.end());

  for (int64_t st : chain) {
    const auto s = unpack[st];
    double wx, wy;
    belief.gridToWorld(s[0], s[1], wx, wy);
    res.states.push_back({wx, wy, s[2] * dth_rad});
  }
  res.cost = g[terminal];
  res.ok = true;
  res.partial = partial;

  // --- coverage audit (faithful to the Python prototype) ----------------------
  // Diagnostic re-check of the returned path: walk each forward/arc primitive and,
  // for every unknown CENTERLINE swept cell, confirm it was predicted-visible from
  // the primitive's start node pose -- or, if that node was reached by motion
  // (mv == 1), from one of the back-projected virtual poses along the incoming
  // heading (the same ladder the planner uses). Uses vis_audit (the audit sensor if
  // one was supplied, else the planner's sensor). Reported, never used to reject.
  int blind = 0;
  for (size_t i = 1; i < chain.size(); ++i) {
    const Came& cm = came[chain[i]];
    if (cm.prim == nullptr || cm.prim->kind == Primitive::TURN) continue;
    const auto sprev = unpack[chain[i - 1]];
    const int ix = sprev[0], iy = sprev[1], ih = sprev[2], mv = sprev[3];
    double px, py;
    belief.gridToWorld(ix, iy, px, py);
    const double pth = ih * dth_rad;
    for (const auto& off : cm.prim->sweep) {
      const int cix = ix + off[0], ciy = iy + off[1];
      if (!belief.isUnknown(cix, ciy)) continue;
      bool ok = vis_audit.visible(px, py, pth, cix, ciy);
      if (!ok && mv) {
        for (int k = 1; k <= P.back_projection_steps; ++k) {
          const double dlt = P.back_projection_step * k;
          const double bx = px - dlt * std::cos(pth);
          const double by = py - dlt * std::sin(pth);
          if (vis_audit.visible(bx, by, pth, cix, ciy)) {
            ok = true;
            break;
          }
        }
      }
      if (!ok) {
        ++blind;
        double bwx, bwy;
        belief.gridToWorld(cix, ciy, bwx, bwy);
        res.blind_cells.push_back({bwx, bwy});
      }
    }
  }
  res.blind_unknown_entries = blind;
  return res;
}

}  // namespace hgp
