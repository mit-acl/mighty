/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <mighty/mighty_solver.hpp>
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

using namespace lbfgs;

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// GCOPTER τ→T and dT/dτ (exactly matches Python)
inline double tau_to_T(double tau)
{
    if (tau >= 0.0)
        return 0.5 * tau * tau + tau + 1.0;
    const double den = 0.5 * tau * tau - tau + 1.0;
    return 1.0 / den;
}
inline double dT_dtau(double tau)
{
    if (tau >= 0.0)
        return tau + 1.0;
    const double den = 0.5 * tau * tau - tau + 1.0;
    return (1.0 - tau) / (den * den);
}
inline double T_to_tau(double T)
{
    if (T >= 1.0)
    {
        double x = std::max(2.0 * T - 1.0, 0.0);
        return std::sqrt(x) - 1.0;
    }
    double x = std::max(2.0 / T - 1.0, 0.0);
    return 1.0 - std::sqrt(x);
}

// Build per-knot T̄ = [ T0, 0.5(T0+T1), ..., 0.5(T_{M-2}+T_{M-1}), T_{M-1} ]
inline void build_Tbar(const std::vector<double> &T, std::vector<double> &Tbar)
{
    const int M = (int)T.size();
    Tbar.resize(M + 1);
    Tbar[0] = T[0];
    for (int i = 1; i < M; ++i)
        Tbar[i] = 0.5 * (T[i - 1] + T[i]);
    Tbar[M] = T[M - 1];
}

// Distribute gTbar[k] (knot derivatives) onto segment derivatives gTextra[s]:
//   dTbar_0/dT0 = 1
//   dTbar_i/dT_{i-1} = 1/2, dTbar_i/dT_i = 1/2  for i=1..M-1
//   dTbar_M/dT_{M-1} = 1
inline void distribute_gTbar_to_segments(const std::vector<double> &gTbar,
                                         std::vector<double> &gTextra)
{
    const int M = (int)gTextra.size();
    // safety: expect gTbar.size() == M+1
    if ((int)gTbar.size() != M + 1)
        return;
    // clear
    std::fill(gTextra.begin(), gTextra.end(), 0.0);
    // endpoints
    gTextra[0] += gTbar[0];
    gTextra[M - 1] += gTbar[M];
    // interiors
    for (int i = 1; i < M; ++i)
    {
        const double v = 0.5 * gTbar[i];
        gTextra[i - 1] += v; // from Tbar_i wrt T_{i-1}
        gTextra[i] += v;     // from Tbar_i wrt T_i
    }
}

// --- GCOPTER cubic hinge (C²) and its derivative ---
// ψ(x; μ) = 0                         , x ≤ 0
//        = (μ - x/2) (x/μ)^3          , 0 < x < μ
//        = x - μ/2                    , x ≥ μ
inline double smoothed_l1(double x, double mu)
{
    if (x <= 0.0)
        return 0.0;
    if (x >= mu)
        return x - 0.5 * mu;
    const double inv_mu = 1.0 / mu;
    const double xd = x * inv_mu;
    const double xd3 = xd * xd * xd;
    return (mu - 0.5 * x) * xd3;
}

// dψ/dx = (x^2 / μ^3) (3μ − 2x)
inline double smoothed_l1_prime(double x, double mu)
{
    if (x <= 0.0)
        return 0.0;
    if (x >= mu)
        return 1.0;
    const double inv_mu = 1.0 / mu;
    const double x2 = x * x;
    return (x2 * (3.0 * mu - 2.0 * x)) * (inv_mu * inv_mu * inv_mu);
}

inline void safe_acos_with_grad(double c_raw, double &theta, double &dtheta_dc)
{
    // clamp to avoid NaNs
    const double c = (c_raw < -1.0 ? -1.0 : (c_raw > 1.0 ? 1.0 : c_raw));
    theta = std::acos(c);

    // if clamped (or extremely close), treat the derivative as zero;
    // this matches the behavior of the forward pass under clamp.
    const double s2 = 1.0 - c * c;
    if (s2 <= 1e-14)
    {
        dtheta_dc = 0.0;
    }
    else
    {
        dtheta_dc = -1.0 / std::sqrt(s2);
    }
}

// --- helper: accumulate CP gradients into (P,V,A) and ∂J/∂T for one segment ---
static inline void pushCPGradToHermite(
    const Eigen::Vector3d gCP[6], // in: ∂J/∂CP_j, j=0..5
    double Ts,                    // in: segment duration
    const Eigen::Vector3d &V0,    // in: V[s]
    const Eigen::Vector3d &A0,    // in: A[s]
    const Eigen::Vector3d &V1,    // in: V[s+1]
    const Eigen::Vector3d &A1,    // in: A[s+1]
    Eigen::Vector3d &gP0,         // out+= ∂J/∂P[s]
    Eigen::Vector3d &gV0,         // out+= ∂J/∂V[s]
    Eigen::Vector3d &gA0,         // out+= ∂J/∂A[s]
    Eigen::Vector3d &gP1,         // out+= ∂J/∂P[s+1]
    Eigen::Vector3d &gV1,         // out+= ∂J/∂V[s+1]
    Eigen::Vector3d &gA1,         // out+= ∂J/∂A[s+1]
    double &gT                    // out+= ∂J/∂T[s] (from CP(T) only)
)
{
    const double T2 = Ts * Ts;

    // CP mapping used in reconstruct():
    // c0=P0
    // c1=P0 + (T/5)V0
    // c2=P0 + (2T/5)V0 + (T^2/20)A0
    // c3=P1 - (2T/5)V1 + (T^2/20)A1
    // c4=P1 - (T/5)V1
    // c5=P1

    // (i) push to endpoints
    gP0 += gCP[0] + gCP[1] + gCP[2];
    gP1 += gCP[3] + gCP[4] + gCP[5];

    // (ii) push to V,A at both ends
    gV0 += (Ts / 5.0) * gCP[1] + (2.0 * Ts / 5.0) * gCP[2];
    gA0 += (T2 / 20.0) * gCP[2];

    gV1 += (-2.0 * Ts / 5.0) * gCP[3] + (-Ts / 5.0) * gCP[4];
    gA1 += (T2 / 20.0) * gCP[3];

    // (iii) explicit ∂J/∂T via the CP(T) dependence
    gT += gCP[1].dot(V0 / 5.0);
    gT += gCP[2].dot((2.0 / 5.0) * V0 + (Ts / 10.0) * A0);
    gT += gCP[3].dot((-2.0 / 5.0) * V1 + (Ts / 10.0) * A1);
    gT += gCP[4].dot((-1.0 / 5.0) * V1);
}

// --- helper: closed-form jerk contribution (cost + ∂ wrt CP + ∂ wrt T) ---
static inline void jerkClosedFormForSegment(
    const std::array<Eigen::Vector3d, 6> &CP, double Ts,
    double &J,              // out+= segment jerk cost
    Eigen::Vector3d gCP[6], // out+= ∂J/∂CP_j
    double &gT,             // out+= ∂J/∂T[s] (1/T^5 part only)
    double jerk_weight)
{
    // finite differences of Bezier control points (third forward differences)
    const Eigen::Vector3d d30 = CP[3] - 3.0 * CP[2] + 3.0 * CP[1] - CP[0];
    const Eigen::Vector3d d31 = CP[4] - 3.0 * CP[3] + 3.0 * CP[2] - CP[1];
    const Eigen::Vector3d d32 = CP[5] - 3.0 * CP[4] + 3.0 * CP[3] - CP[2];

    const double invT5 = 1.0 / (Ts * Ts * Ts * Ts * Ts);
    const double C = 3600.0 * invT5;

    const double S = d30.squaredNorm() + d31.squaredNorm() + d32.squaredNorm();
    J += C * S;

    auto add = [&](int i, const Eigen::Vector3d &v)
    { gCP[i] += jerk_weight * 2.0 * C * v; };

    // backprop through d3* = linear combos of CP
    // m=0
    add(0, -d30);
    add(1, 3.0 * d30);
    add(2, -3.0 * d30);
    add(3, d30);
    // m=1
    add(1, -d31);
    add(2, 3.0 * d31);
    add(3, -3.0 * d31);
    add(4, d31);
    // m=2
    add(2, -d32);
    add(3, 3.0 * d32);
    add(4, -3.0 * d32);
    add(5, d32);

    // explicit derivative of 1/T^5 factor
    gT += jerk_weight * -5.0 * (3600.0 / (Ts * Ts * Ts * Ts * Ts * Ts)) * S;
}

inline void bernstein5(double u, double B[6])
{
    const double um = 1.0 - u;
    const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, u5 = u4 * u;
    const double m2 = um * um, m3 = m2 * um, m4 = m3 * um, m5 = m4 * um;
    B[0] = m5;
    B[1] = 5.0 * u * m4;
    B[2] = 10.0 * u2 * m3;
    B[3] = 10.0 * u3 * m2;
    B[4] = 5.0 * u4 * um;
    B[5] = u5;
}

inline void bernstein4(double u, double b[5])
{
    const double um = 1.0 - u;
    const double u2 = u * u, u3 = u2 * u, u4 = u3 * u;
    const double m2 = um * um, m3 = m2 * um, m4 = m3 * um;
    b[0] = m4;
    b[1] = 4.0 * u * m3;
    b[2] = 6.0 * u2 * m2;
    b[3] = 4.0 * u3 * um;
    b[4] = u4;
}

// Cubic Bernstein basis B^3_k(s), k=0..3
inline void bernstein3(double s, double B[4])
{
    // (optional) clamp for numeric safety
    if (s < 0.0)
        s = 0.0;
    if (s > 1.0)
        s = 1.0;

    const double t = 1.0 - s;
    const double t2 = t * t;
    const double s2 = s * s;

    B[0] = t * t2;       // (1 - s)^3
    B[1] = 3.0 * s * t2; // 3 s (1 - s)^2
    B[2] = 3.0 * s2 * t; // 3 s^2 (1 - s)
    B[3] = s * s2;       // s^3
}

// Quadratic Bernstein basis B^2_k(s), k=0..2
inline void bernstein2(double s, double B[3])
{
    // (optional) clamp for numeric safety
    if (s < 0.0)
        s = 0.0;
    if (s > 1.0)
        s = 1.0;

    const double t = 1.0 - s;

    B[0] = t * t;       // (1 - s)^2
    B[1] = 2.0 * s * t; // 2 s (1 - s)
    B[2] = s * s;       // s^2
}

// Precompute Bernstein bases for all sample nodes j=0..kappa.
// Safe to call repeatedly; it only rebuilds when kappa changes.
void SolverMIGHTY::ensureBernsteinCache(int kappa) const
{
    if (kappa <= 0)
        return;
    if (bern_.kappa == kappa)
        return;

    bern_.kappa = kappa;
    const int N = kappa + 1;
    bern_.B5.resize(N);
    bern_.B4.resize(N);
    bern_.B3.resize(N);
    bern_.B2.resize(N);
    bern_.W.resize(N);

    for (int j = 0; j < N; ++j)
    {
        const double tau = static_cast<double>(j) / static_cast<double>(kappa);
        bernstein5(tau, bern_.B5[j].data());
        bernstein4(tau, bern_.B4[j].data());
        bernstein3(tau, bern_.B3[j].data());
        bernstein2(tau, bern_.B2[j].data());
        bern_.W[j] = (j == 0 || j == kappa) ? 0.5 : 1.0; // trapezoid weight
    }
}

// Numerically stable central-difference directional derivative
double SolverMIGHTY::centralDiff(const VecXd &z,
                                const VecXd &d,
                                double eps_base)
{
    // Scale step per IEEE advice
    double eps = eps_base / std::max(1.0, d.norm());
    VecXd zp = z + eps * d;
    VecXd zm = z - eps * d;
    return (evaluateObjective(zp) - evaluateObjective(zm)) / (2.0 * eps);
}

void SolverMIGHTY::checkGradDirectional(const VecXd &z0, int num_dirs, double eps, unsigned seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> N(0.0, 1.0);

    VecXd g(z0.size());
    (void)evaluateObjectiveAndGradient(z0, g);

    double max_rel_err = 0.0;
    for (int k = 0; k < num_dirs; ++k)
    {
        VecXd d(z0.size());
        for (int i = 0; i < d.size(); ++i)
            d[i] = N(rng);

        // No projection needed; all variables are free.
        double n = d.norm();
        if (n == 0.0)
        {
            --k;
            continue;
        } // extremely unlikely
        d /= n; // makes eps a consistent step

        const double fd = centralDiff(z0, d, eps);
        const double ad = g.dot(d);
        const double denom = std::max(1e-12, std::abs(fd) + std::abs(ad));
        const double rel = std::abs(fd - ad) / denom;

        max_rel_err = std::max(max_rel_err, rel);
        if (rel > 1e-3)
            printf("\033[1;31m [dir %d] fd=%.6f ad=%.6f rel_err=%.6f \033[0m\n", k, fd, ad, rel);
    }
}

void SolverMIGHTY::checkGradCoordinates(const VecXd &z0, int max_coords, double eps, unsigned seed)
{
    VecXd g(z0.size());
    (void)evaluateObjectiveAndGradient(z0, g);

    // All coordinates are free now.
    std::vector<int> idx(z0.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::mt19937 rng(seed);
    std::shuffle(idx.begin(), idx.end(), rng);
    if ((int)idx.size() > max_coords)
        idx.resize(max_coords);

    double worst = 0.0;
    for (int i : idx)
    {
        VecXd ei = VecXd::Zero(z0.size());
        ei[i] = 1.0; // unit basis direction (no projection)

        const double g_fd = centralDiff(z0, ei, eps);
        const double g_ad = g[i];
        const double abserr = std::abs(g_fd - g_ad);
        const double rel = abserr / std::max(1e-12, std::abs(g_fd) + std::abs(g_ad));
        worst = std::max(worst, rel);

        if (rel > 1e-3)
            printf("\033[1;31m [idx %d] g_fd=%.6f g_ad=%.6f abs_err=%.6f rel_err=%.6f \033[0m\n",
                   i, g_fd, g_ad, abserr, rel);
    }
}

void SolverMIGHTY::seamBackward(int seam,
                               const Eigen::VectorXd &q,
                               const Eigen::Vector3d &gp,
                               Eigen::VectorXd &gq) const
{
    const auto &OB = vPolys_OB_[2 * seam + 1];
    const int k = (int)OB.cols();

    // ∂J/∂w = B^T gp, where B = OB.rightCols(k-1)
    const Eigen::VectorXd gw = OB.rightCols(k - 1).transpose() * gp;

    const double n = std::max(q.norm(), 1e-12);
    const Eigen::VectorXd u = q / n;

    // ∂J/∂u: [2 u_head ⊙ gw; 0]
    Eigen::VectorXd gu = Eigen::VectorXd::Zero(k);
    gu.head(k - 1).array() = 2.0 * u.head(k - 1).array() * gw.array();

    // gq += ((I - u u^T)/n) * gu  ==>  gq += (gu - u*(u.dot(gu))) / n
    gq += (gu - u * (u.dot(gu))) / n;
}

// -----------------------------------------------------------------------------

SolverMIGHTY::SolverMIGHTY()
{
    // Constructor implementation
}

// -----------------------------------------------------------------------------

SolverMIGHTY::~SolverMIGHTY()
{
    // Destructor implementation
    // Clean up resources if necessary
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::initializeSolver(const planner_params_t &params, const Eigen::VectorXd &physical_params)
{
    // Initialize the solver with parameters that won't change throughout the mission
    V_max_ = params.V_max;                                 // Max velocity
    time_weight_ = params.time_weight;                     // Weight for time in the objective
    stat_weight_ = params.stat_weight;                     // Weight for static avoidance in the objective
    jerk_weight_ = params.jerk_weight;                     // Weight for jerk in the objective
    dyn_constr_vel_weight_ = params.dyn_constr_vel_weight; // Weight for dynamic velocity constraints
    Co_ = params.Co;                                             // Clearance distance for static obstacle avoidance
    BIG_ = params.BIG;                                           // A large constant for static constraints
    dc_ = params.dc;                                             // Discretization constant
    V_min_ = 0.0;                                                // Minimum speed
    turn_buf_ = params.init_turn_bf * M_PI / 180;                // 15° in radians
    turn_span_ = M_PI - turn_buf_;                               // over which we ramp down

    integral_resolution_ = params.integral_resolution; // e.g., 30
    hinge_mu_ = params.hinge_mu;                       // e.g., 1e-2
    omege_max_ = params.omega_max;                     // e.g., 6.0
    tilt_max_rad_ = params.tilt_max_rad;               // e.g., 35° in rad
    f_min_ = params.f_min;
    f_max_ = params.f_max;
    mass_ = params.mass;
    g_ = params.g;

    dyn_constr_bodyrate_weight_ = params.dyn_constr_bodyrate_weight;
    dyn_constr_tilt_weight_ = params.dyn_constr_tilt_weight;
    dyn_constr_thrust_weight_ = params.dyn_constr_thrust_weight;

    // flat map
    flatmap_.reset(physical_params(0), physical_params(1), physical_params(2), physical_params(3), physical_params(4), physical_params(5));
}

// -----------------------------------------------------------------------------

// sign function
template <typename T>
int sgn(T val)
{
return (T(0) < val) - (val < T(0));
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::prepareSolverForReplan(double t0,
                                         const vec_Vec3f &global_wps,
                                         const std::vector<LinearConstraint3D> &safe_corridor,
                                         const Eigen::VectorXd &initial_xi,
                                         const Eigen::VectorXd &init_times,
                                         const state &initial_state,
                                         const state &goal_state,
                                         double &initial_guess_computation_time,
                                         const PolyhedraV &vPolys,
                                         bool use_multiple_initial_guesses)
{

    // ---- initialize variables ----
    t0_ = t0;
    x0_ = initial_state.pos;
    v0_ = initial_state.vel;
    a0_ = initial_state.accel;
    const Vec3 goal_pos = goal_state.pos; // we’ll set xf_ after we build waypoints
    vf_ = goal_state.vel;
    af_ = goal_state.accel;

    // ---- static corridor planes ----
    setStaticConstraints(safe_corridor);

    // ---- corridor OBs BEFORE xi decode ----
    vPolys_OB_ = vPolys; // copy corridor polyhedra

    // Derive M_ from corridor: expect size = 2*M_ - 1
    const int Ncorr = static_cast<int>(vPolys_OB_.size());
    M_ = (Ncorr + 1) / 2;

    // Preallocate
    P_s_.resize(M_ + 1);
    V_s_.resize(M_ + 1);
    A_s_.resize(M_ + 1);
    CP_s_.resize(M_);
    T_s_.resize(M_);
    gP_s_.resize(M_ + 1);
    gV_s_.resize(M_ + 1);
    gA_s_.resize(M_ + 1);
    gT_s_.resize(M_);

    // Layout: sizes/offsets for seam q-blocks
    seamSizes_.resize(M_ - 1);
    int sum_k = 0;
    for (int s = 0; s < M_ - 1; ++s)
    {
        const int k = static_cast<int>(vPolys_OB_[2 * s + 1].cols());
        seamSizes_(s) = k;
        sum_k += k;
    }
    seamOffsets_.resize(M_ - 1);

    // xi (seams) unchanged — after you computed sum_k
    offXi_ = 0;
    offPM_ = offXi_ + sum_k; // sentinel = end of xi (no P0/PM slots)

    // Only INTERIOR v̂, â: i = 1..M_-1  →  count = M_-1 (can be zero)
    const int nVhatFree = std::max(0, M_ - 1);
    const int nAhatFree = nVhatFree;

    offVhat_ = offPM_;
    offAhat_ = offVhat_ + 3 * nVhatFree;
    offTau_ = offAhat_ + 3 * nAhatFree;

    K_ = offTau_ + M_; // τ unchanged (one per segment)

    for (int s = 0, acc = offXi_; s < M_ - 1; ++s)
    {
        seamOffsets_[s] = acc;
        acc += seamSizes_(s);
    }
    corridor_q_active_ = true;

    // ---- now it’s safe to DECODE xi into internal waypoints ----
    global_wps_.clear();
    std::vector<Eigen::VectorXd> q_seams;
    global_wps_.resize(M_ + 1);
    global_wps_[0] = x0_;
    global_wps_[M_] = goal_pos;

    Eigen::VectorXd r;
    int overlaps = M_ - 1; // number of corridor overlaps (seams)
    for (int i = 0, j = 0, k; i < overlaps; i++, j += k)
    {
        k = vPolys_OB_[2 * i + 1].cols();
        Eigen::Map<const Eigen::VectorXd> q(initial_xi.data() + j, k);
        r = q.normalized().head(k - 1);
        global_wps_[i + 1] = vPolys_OB_[2 * i + 1].rightCols(k - 1) * r.cwiseProduct(r) + vPolys_OB_[2 * i + 1].col(0);
    }

    // Set xf_ **after** building the waypoint list
    xf_ = global_wps_.back();

    // ---- cheap segment times, then V/A ----
    std::vector<double> T;
    std::vector<Vec3> P(global_wps_.begin(), global_wps_.end());
    std::vector<Vec3> V, A;
    findInitialGuess(T, V, A);

    // ---- pack z0 (uses xi blocks directly) ----
    Eigen::VectorXd z0;
    packDecisionVariables(initial_xi, P, V, A, T, z0);

    // (debug) reconstruct & check
    std::vector<Vec3> Pchk, Vchk, Achk;
    std::vector<std::array<Vec3, 6>> CP;
    std::vector<double> Tchk;
    reconstruct(z0, Pchk, Vchk, Achk, CP, Tchk);

    for (int i = 1; i < M_; ++i)
    {
        if ((Pchk[i] - P[i]).norm() > 1e-4)
            std::cerr << "[init] seam " << (i - 1) << " decode error = "
                      << (Pchk[i] - P[i]).norm() << "\n";

        if ((Vchk[i] - V[i]).norm() > 1e-4)
            std::cerr << "[init] seam " << (i - 1) << " V decode error = "
                      << (Vchk[i] - V[i]).norm() << "\n";

        if ((Achk[i] - A[i]).norm() > 1e-4)
            std::cerr << "[init] seam " << (i - 1) << " A decode error = "
                      << (Achk[i] - A[i]).norm() << "\n";

        if (std::abs(Tchk[i - 1] - T[i - 1]) > 1e-4)
            std::cerr << "[init] segment " << (i - 1) << " T decode error = "
                      << std::abs(Tchk[i - 1] - T[i - 1]) << "\n";
    }

    list_z0_.clear();
    list_z0_.push_back(z0);

    // wherever you have a valid z (e.g., right before starting L-BFGS)
    // checkGradDirectional(z0, /*num_dirs=*/8, /*eps=*/1e-5, /*seed=*/42);
    // checkGradCoordinates(z0, /*max_coords=*/256, /*eps=*/1e-5, /*seed=*/43);

}

// -----------------------------------------------------------------------------
// Extract q-blocks (one per seam). Works for either "full xi" (endpoints+seams)
// or "seams-only xi". Returns false if the length doesn't match expectations.
bool SolverMIGHTY::extractSeamQFromInitialXi(
    const Eigen::VectorXd &initial_xi,
    std::vector<Eigen::VectorXd> &q_seams) const
{
    // Preconditions
    if (M_ < 2 || vPolys_OB_.empty())
    {
        std::cerr << "[xi] precondition fail: M_=" << M_
                  << " vPolys_OB_.size()=" << vPolys_OB_.size() << "\n";
        return false;
    }
    const int expected = 2 * M_ - 1;
    if ((int)vPolys_OB_.size() != expected)
    {
        std::cerr << "[xi] corridor OB list size " << vPolys_OB_.size()
                  << " != 2*M_-1 (" << expected << ")\n";
        return false;
    }

    auto k_of_idx = [&](int idx)
    { return (int)vPolys_OB_[idx].cols(); };

    // Compute expected lengths
    int seams_len = 0;
    for (int s = 0; s < M_ - 1; ++s)
        seams_len += k_of_idx(2 * s + 1);

    int full_len = 0;
    for (int i = 0; i <= M_; ++i)
    {
        const int l = (i == 0) ? 0 : (i == M_ ? 2 * M_ - 2 : 2 * i - 1);
        full_len += k_of_idx(l);
    }

    q_seams.clear();
    q_seams.resize(M_ - 1);

    if (initial_xi.size() == seams_len)
    {
        int acc = 0;
        for (int s = 0; s < M_ - 1; ++s)
        {
            const int k = k_of_idx(2 * s + 1);
            q_seams[s] = initial_xi.segment(acc, k);
            acc += k;
        }
        return true;
    }
    if (initial_xi.size() == full_len)
    {
        int j = 0;
        for (int i = 0; i <= M_; ++i)
        {
            const int l = (i == 0) ? 0 : (i == M_ ? 2 * M_ - 2 : 2 * i - 1);
            const int k = k_of_idx(l);
            if (i >= 1 && i <= M_ - 1)
            {
                const int s = i - 1;
                q_seams[s] = initial_xi.segment(j, k);
            }
            j += k;
        }
        return true;
    }

    return false;
}

// -----------------------------------------------------------------------------

// Layout assumed (as you set offsets already):
// z = [ P0(3),  q_seams (variable sizes),  PM(3),  Vhat(3*(M+1)),  Ahat(3*(M+1)),  tau(M) ]
void SolverMIGHTY::packDecisionVariables(
    const Eigen::VectorXd &initial_xi, // <- NEW
    const std::vector<Vec3> &P,
    const std::vector<Vec3> &V,
    const std::vector<Vec3> &A,
    const std::vector<double> &T,
    VecXd &z) const
{
    z.setZero(K_);

    // 2) seams: copy xi blocks directly into z
    if (corridor_q_active_)
    {
        std::vector<Eigen::VectorXd> q_seams;
        const bool ok = extractSeamQFromInitialXi(initial_xi, q_seams);
        if (!ok)
        {
            // Fallback to encoding from P[i]
            std::cerr << "[pack] initial_xi size mismatch, re-encoding seams from P\n";
        }
        else
        {
            for (int s = 0; s < M_ - 1; ++s)
            {
                const int off = seamOffsets_[s];
                const int k = seamSizes_(s);
                // size check (defensive)
                if (q_seams[s].size() == k)
                    z.segment(off, k) = q_seams[s];
                else
                {
                    std::cerr << "[pack] seam " << s << " xi size " << q_seams[s].size()
                              << " != expected k=" << k << ", re-encoding from P\n";
                }
            }
        }
    }
    else
    {
        std::cerr << "[pack] warning: packing seams without corridor_q_active_\n";
    }

    // 3) v̂, â  (with T̄)
    // 3) v-blocks and a-blocks (INTERIOR only: i = 1..M_-1)
    std::vector<double> Tbar; // still OK to build, we just won't use it when unscaled
    build_Tbar(T, Tbar);

    for (int i = 1; i < M_; ++i)
    {
        const int off_v = vhatOffset(i);
        const int off_a = ahatOffset(i);
        if (scale_derivs_)
        {
            const double Tb = Tbar[i];
            z.segment<3>(off_v) = Tb * V[i];
            z.segment<3>(off_a) = (Tb * Tb) * A[i];
        }
        else
        {
            // Unscaled benchmark: store raw V and A
            z.segment<3>(off_v) = V[i];
            z.segment<3>(off_a) = A[i];
        }
    }

    // 4) τ
    for (int s = 0; s < M_; ++s)
        z[offTau_ + s] = T_to_tau(T[s]);
}

// -----------------------------------------------------------------------------

// >>> scatter ∂J/∂P into z: endpoints go to Cartesian slots; seams backprop to q-blocks
void SolverMIGHTY::scatterPosGrads(const std::vector<Vec3> &gP,
                                  const VecXd &z, VecXd &grad) const
{
    for (int i = 1; i < M_; ++i)
    {
        const int seam = i - 1, k = seamSizes_(seam), off = seamOffsets_[seam];
        Eigen::VectorXd gq = Eigen::VectorXd::Zero(k);
        seamBackward(seam, z.segment(off, k), gP[i], gq);
        grad.segment(off, k) += gq;
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::reconstruct(
    const VecXd &z,
    std::vector<Vec3> &P,
    std::vector<Vec3> &V,
    std::vector<Vec3> &A,
    std::vector<std::array<Vec3, 6>> &CP,
    std::vector<double> &T) const
{
    const int knotCount = M_ + 1;
    P.resize(knotCount);
    V.resize(knotCount);
    A.resize(knotCount);
    T.resize(M_);
    CP.resize(M_);

    // 1) τ -> T
    for (int s = 0; s < M_; ++s)
        T[s] = tau_to_T(z[offTau_ + s]);

    // 2) T̄
    std::vector<double> Tbar(knotCount);
    build_Tbar(T, Tbar);

    // 3) decode P (endpoints + seams)
    Eigen::VectorXd xi = z.segment(offXi_, offPM_ - offXi_).eval();
    for (int i = 0, j = 0, k; i < M_ - 1; ++i, j += k)
    {
        k = vPolys_OB_[2 * i + 1].cols();
        Eigen::Map<const Eigen::VectorXd> q(xi.data() + j, k);
        Eigen::VectorXd r = q.normalized().head(k - 1);
        P[i + 1] = vPolys_OB_[2 * i + 1].rightCols(k - 1) * r.cwiseProduct(r) + vPolys_OB_[2 * i + 1].col(0);
    }

    // Endpoints: fixed
    P.front() = x0_;
    P.back() = xf_;

    // v̂/â → V/A
    // endpoints fixed:
    V[0] = v0_;
    A[0] = a0_;
    V[M_] = vf_;
    A[M_] = af_;

    // 4) v̂, â -> V, A
    // 4) derivative blocks -> V, A (INTERIOR only)
    for (int i = 1; i < M_; ++i)
    {
        const int voff = vhatOffset(i);
        const int aoff = ahatOffset(i);
        if (scale_derivs_)
        {
            const double Tb = std::max(1e-12, Tbar[i]);
            V[i] = z.segment<3>(voff) / Tb;
            A[i] = z.segment<3>(aoff) / (Tb * Tb);
        }
        else
        {
            // Unscaled benchmark: raw values
            V[i] = z.segment<3>(voff);
            A[i] = z.segment<3>(aoff);
        }
    }

    // 6) build CP from (P,V,A,T)
    for (int s = 0; s < M_; ++s)
    {
        const Vec3 &p0 = P[s], &v0s = V[s], &a0s = A[s];
        const Vec3 &p1 = P[s + 1], &v1s = V[s + 1], &a1s = A[s + 1];
        const double Ts = T[s], T2 = Ts * Ts;
        auto &c = CP[s];
        c[0] = p0;
        c[1] = p0 + (Ts / 5.0) * v0s;
        c[2] = p0 + (2.0 * Ts / 5.0) * v0s + (T2 / 20.0) * a0s;
        c[3] = p1 - (2.0 * Ts / 5.0) * v1s + (T2 / 20.0) * a1s;
        c[4] = p1 - (Ts / 5.0) * v1s;
        c[5] = p1;
    }

    auto insideIntersection = [&](int seam, const Eigen::Vector3d &p)
    {
        const auto &A1 = A_stat_[seam];
        const auto &b1 = b_stat_[seam];
        const auto &A2 = A_stat_[seam + 1];
        const auto &b2 = b_stat_[seam + 1];

        const double eps = 1e-8; // inward slack
        // Inside means A p >= b (+eps)
        for (int r = 0; r < A1.rows(); ++r)
            if (A1.row(r).dot(p) - b1[r] < eps) // <-- flipped sign
                return false;
        for (int r = 0; r < A2.rows(); ++r)
            if (A2.row(r).dot(p) - b2[r] < eps) // <-- flipped sign
                return false;
        return true;
    };

    // 7) (optional) sanity: seam points are inside intersection
    for (int i = 1; i < M_; ++i)
    {
        if (!insideIntersection(i - 1, P[i]))
        {
            // std::cerr << "[sanity] seam " << (i - 1) << " waypoint violates intersection\n";
        }
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::findInitialGuess(
    std::vector<double> &T,
    std::vector<Vec3> &V,
    std::vector<Vec3> &A)
{
    // --- 1) clear outputs and reserve ---
    T.clear();
    T.reserve(M_);
    V.clear();
    V.reserve(M_ + 1);
    A.clear();
    A.reserve(M_ + 1);

    // --- 3) precompute segment directions ---
    std::vector<Vec3> dirs(M_);
    for (int i = 0; i < M_; ++i)
    {
        Vec3 Δ = global_wps_[i + 1] - global_wps_[i];
        double d = Δ.norm();
        dirs[i] = (d > 0 ? (Δ / d).eval() : Vec3::UnitX());
    }

    // --- 4) per‐waypoint speeds and velocity vectors ---
    std::vector<double> s(M_ + 1);

    s[0] = v0_.norm();
    s[M_] = vf_.norm();

    for (int i = 1; i < M_; ++i)
    {
        double ca = std::clamp(dirs[i - 1].dot(dirs[i]), -1.0, 1.0);
        double angle = std::acos(ca); // 0..π

        double factor;
        if (angle <= turn_buf_)
        {
            factor = 1.0; // within buffer, no slow-down
        }
        else
        {
            // linearly ramp from 1 at turn_buf_ → 0 at π
            factor = 1.0 - (angle - turn_buf_) / turn_span_;
            factor = std::clamp(factor, 0.0, 1.0);
        }

        s[i] = std::max(V_min_, factor * V_max_);
    }

    // build velocity vectors at each waypoint
    V[0] = v0_;
    for (int i = 1; i < M_; ++i)
    {
        Vec3 bis = dirs[i - 1] + dirs[i];
        double n = bis.norm();
        if (n > 0)
            bis /= n;
        else
            bis = dirs[i];
        V[i] = bis * s[i];
    }
    V[M_] = vf_;

    // if second to last waypoint is close to the last, scale down its velocity
    // double d_last = (global_wps_[M_] - global_wps_[M_ - 1]).norm();
    // if (d_last < 2.0)
    // {
    //     double f = d_last / 2.0;
    //     f = std::clamp(f, 0.1, 1.0);
    //     V[M_ - 1] *= f;
    // }

    // --- 5) allocate time per segment using average vector velocity ---
    for (int i = 0; i < M_; ++i)
    {
        double dist = (global_wps_[i + 1] - global_wps_[i]).norm();
        Vec3 v_avg = 0.5 * (V[i] + V[i + 1]);
        double speed = std::max(V_min_, v_avg.norm());
        T.push_back(dist / speed);
    }

    // --- 6) generate initial A ---
    solveMinJerkAccOnlyClosedForm(global_wps_, T, V, a0_, af_, A);

}

// -----------------------------------------------------------------------------

//  Closed-form quintic control‐points for one segment
std::array<Vec3, 6> SolverMIGHTY::computeQuinticCP(
    const Vec3 &P0, const Vec3 &V0, const Vec3 &A0,
    const Vec3 &P1, const Vec3 &V1, const Vec3 &A1,
    double T)
{
    double T2 = T * T, T3 = T2 * T, T4 = T3 * T, T5 = T4 * T;
    Vec3 c0 = P0;
    Vec3 c1 = V0;
    Vec3 c2 = 0.5 * A0;
    Vec3 c3 = (-20.0 * P0 + 20.0 * P1 - (8.0 * V1 + 12.0 * V0) * T - (3.0 * A0 - A1) * T2) / (2.0 * T3);
    Vec3 c4 = (30.0 * P0 - 30.0 * P1 + (14.0 * V1 + 16.0 * V0) * T + (3.0 * A0 - 2.0 * A1) * T2) / (2.0 * T4);
    Vec3 c5 = (-12.0 * P0 + 12.0 * P1 - (6.0 * V1 + 6.0 * V0) * T - (A0 - A1) * T2) / (2.0 * T5);

    std::array<Vec3, 6> CP;
    for (int i = 0; i < 6; ++i)
    {
        double t = T * (double(i) / 5.0);
        double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
        CP[i] = c0 + c1 * t + c2 * t2 + c3 * t3 + c4 * t4 + c5 * t5;
    }
    return CP;
}

// -----------------------------------------------------------------------------

inline Eigen::Matrix4d SolverMIGHTY::K_r(double T)
{
    double t2 = T * T, t3 = t2 * T;
    double i2 = 1.0 / t2, i3 = 1.0 / t3;
    Eigen::Matrix4d K;
    K << 192 * i3, 36 * i2, 168 * i3, -24 * i2,
        36 * i2, 9.0 / T, 24 * i2, -3.0 / T,
        168 * i3, 24 * i2, 192 * i3, -36 * i2,
        -24 * i2, -3.0 / T, -36 * i2, 9.0 / T;
    return K;
}

// -----------------------------------------------------------------------------

inline std::array<double, 4> SolverMIGHTY::k_r(double T, double P0, double P1)
{
    double k0 = (-66960 * P0 * T + 66960 * P1 * T) / (2 * T * T * T * T * T) + (-24480 * P0 * T + 24480 * P1 * T) / (2 * T * T * T * T * T) + (4320 * P0 * T - 4320 * P1 * T) / (2 * T * T * T * T * T) + (25920 * P0 * T - 25920 * P1 * T) / (2 * T * T * T * T * T) + (61920 * P0 * T - 61920 * P1 * T) / (2 * T * T * T * T * T);

    double k1 = (-11880 * P0 * T * T + 11880 * P1 * T * T) / (2 * T * T * T * T * T) + (-5400 * P0 * T * T + 5400 * P1 * T * T) / (2 * T * T * T * T * T) + (1080 * P0 * T * T - 1080 * P1 * T * T) / (2 * T * T * T * T * T) + (4320 * P0 * T * T - 4320 * P1 * T * T) / (2 * T * T * T * T * T) + (12000 * P0 * T * T - 12000 * P1 * T * T) / (2 * T * T * T * T * T);

    double k2 = (-62640 * P0 * T + 62640 * P1 * T) / (2 * T * T * T * T * T) + (-18720 * P0 * T + 18720 * P1 * T) / (2 * T * T * T * T * T) + (2880 * P0 * T - 2880 * P1 * T) / (2 * T * T * T * T * T) + (25920 * P0 * T - 25920 * P1 * T) / (2 * T * T * T * T * T) + (53280 * P0 * T - 53280 * P1 * T) / (2 * T * T * T * T * T);

    double k3 = (-7680 * P0 * T * T + 7680 * P1 * T * T) / (2 * T * T * T * T * T) + (-4320 * P0 * T * T + 4320 * P1 * T * T) / (2 * T * T * T * T * T) + (-360 * P0 * T * T + 360 * P1 * T * T) / (2 * T * T * T * T * T) + (2520 * P0 * T * T - 2520 * P1 * T * T) / (2 * T * T * T * T * T) + (9720 * P0 * T * T - 9720 * P1 * T * T) / (2 * T * T * T * T * T);
    return {k0, k1, k2, k3};
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::assemble_H_b(
    const std::vector<Vec3> &wps,
    const std::vector<double> &T,
    const Vec3 &v0, const Vec3 &a0,
    const Vec3 &vf, const Vec3 &af,
    Eigen::MatrixXd &H,
    Eigen::MatrixXd &b)
{
    int M = (int)T.size();
    int N = M - 1; // # interior junctions
    int D = 3;     // 3d

    H = Eigen::MatrixXd::Zero(2 * N, 2 * N);
    b = Eigen::MatrixXd::Zero(2 * N, D);

    // for each segment r = 0..M-1
    for (int r = 0; r < M; ++r)
    {
        // 1) Hessian block
        Eigen::Matrix4d Kr = K_r(T[r]);

        // 2) constant-term block (4×3)
        Eigen::Matrix<double, 4, 3> kr;
        for (int dim = 0; dim < 3; ++dim)
        {
            auto c = k_r(T[r],
                         wps[r](dim),
                         wps[r + 1](dim));
            for (int i = 0; i < 4; ++i)
                kr(i, dim) = c[i];
        }

        // 3) left half → junction r-1
        if (r >= 1)
        {
            int i = r - 1;
            H.block<2, 2>(2 * i, 2 * i) += Kr.block<2, 2>(0, 0);
            b.block<2, 3>(2 * i, 0) += kr.block<2, 3>(0, 0);
        }

        // 4) right half → junction r
        if (r <= M - 2)
        {
            int j = r;
            H.block<2, 2>(2 * j, 2 * j) += Kr.block<2, 2>(2, 2);
            b.block<2, 3>(2 * j, 0) += kr.block<2, 3>(2, 0);
        }

        // 5) coupling off-diagonals (H only)
        if (r >= 1 && r <= M - 2)
        {
            int i = r - 1, j = r;
            H.block<2, 2>(2 * i, 2 * j) += Kr.block<2, 2>(0, 2);
            H.block<2, 2>(2 * j, 2 * i) += Kr.block<2, 2>(2, 0);
        }
    }

    // 6) endpoint cross-terms into b
    {
        // segment 0: x0 (known) → x1 interior
        Eigen::Matrix4d K0 = K_r(T[0]);
        Eigen::Matrix<double, 2, 3> v0a0;
        v0a0.row(0) = v0.transpose();
        v0a0.row(1) = a0.transpose();
        // b[0:2] += K0[2:4,0:2]*[v0;a0]
        b.block<2, 3>(0, 0) += K0.block<2, 2>(2, 0) * v0a0;
    }
    {
        // segment M-1: x_{M} known → x_{M-1} interior
        Eigen::Matrix4d KM = K_r(T.back());
        Eigen::Matrix<double, 2, 3> vfaf;
        vfaf.row(0) = vf.transpose();
        vfaf.row(1) = af.transpose();
        int i = 2 * (M - 2);
        // b[i:i+2] += KM[0:2,2:4]*[vf;af]
        b.block<2, 3>(i, 0) += KM.block<2, 2>(0, 2) * vfaf;
    }
}

// -----------------------------------------------------------------------------

//  Build and solve the banded min-jerk system:
void SolverMIGHTY::solveMinJerkVelAcc(
    const std::vector<Vec3> &wps, // size M+1
    const std::vector<double> &T, // size M
    const Vec3 &v0, const Vec3 &a0,
    const Vec3 &vf, const Vec3 &af,
    std::vector<Vec3> &V, // out size M+1
    std::vector<Vec3> &A  // out size M+1
)
{

    int M = (int)T.size();

    Eigen::MatrixXd H, b;
    assemble_H_b(wps, T, v0, a0, vf, af, H, b);

    // Solve H x = -b   →   x is (2N×3)
    Eigen::MatrixXd X = H.ldlt().solve(-b);

    // Unpack into V,A
    V.clear();
    A.clear();
    V.resize(M + 1);
    A.resize(M + 1);
    V[0] = v0;
    A[0] = a0;
    V[M] = vf;
    A[M] = af;
    for (int i = 1; i < M; ++i)
    {
        // row 2*(i-1) of X is the velocity at junction i
        V[i] = X.row(2 * (i - 1)).transpose(); // Vec3 = (3×1)
        // row 2*(i-1)+1 of X is the acceleration at junction i
        A[i] = X.row(2 * (i - 1) + 1).transpose();
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::assemble_H_b_acc_only(
    const std::vector<Vec3> &wps,         // size M+1
    const std::vector<double> &T,         // size M
    const std::vector<Vec3> &V,           // size M+1 (ALL velocities fixed)
    const Vec3 &a0, const Vec3 &af,       // fixed endpoint accelerations
    Eigen::MatrixXd &H,                   // out: (N×N)
    Eigen::MatrixXd &b                    // out: (N×3)
)
{
    const int M = static_cast<int>(T.size());
    const int N = M - 1;
    if (N <= 0) {
        H.resize(0, 0);
        b.resize(0, 3);
        return;
    }

    H = Eigen::MatrixXd::Zero(N, N);
    b = Eigen::MatrixXd::Zero(N, 3);

    for (int r = 0; r < M; ++r)
    {
        const Eigen::Matrix4d Kr = K_r(T[r]);  // 4×4 over [V_r, A_r, V_{r+1}, A_{r+1}]

        // k_r as a (4×3) block, column-wise per dimension
        Eigen::Matrix<double, 4, 3> kr;
        for (int dim = 0; dim < 3; ++dim) {
            const auto c = k_r(T[r], wps[r](dim), wps[r+1](dim)); // {k0,k1,k2,k3}
            kr(0, dim) = c[0];
            kr(1, dim) = c[1];
            kr(2, dim) = c[2];
            kr(3, dim) = c[3];
        }

        if (r == 0)
        {
            // Only A1 is unknown → row uses A_{r+1} (index 3 in Kr)
            const int j = 0; // maps to A1
            H(j, j) += Kr(3, 3);
            // b_j += Kr[3,0]*V_r + Kr[3,2]*V_{r+1} + Kr[3,1]*a0 + k_r[3]
            b.row(j) += (Kr(3, 0) * V[r].transpose()
                      +  Kr(3, 2) * V[r+1].transpose()
                      +  Kr(3, 1) * a0.transpose()
                      +  kr.row(3));
        }
        else if (r == M - 1)
        {
            // Only A_{M-1} is unknown → row uses A_r (index 1 in Kr)
            const int i = N - 1; // maps to A_{M-1}
            H(i, i) += Kr(1, 1);
            // b_i += Kr[1,0]*V_r + Kr[1,2]*V_{r+1} + Kr[1,3]*af + k_r[1]
            b.row(i) += (Kr(1, 0) * V[r].transpose()
                      +  Kr(1, 2) * V[r+1].transpose()
                      +  Kr(1, 3) * af.transpose()
                      +  kr.row(1));
        }
        else
        {
            // Both A_r and A_{r+1} are unknown
            const int i = r - 1; // maps to A_r
            const int j = r;     // maps to A_{r+1}

            // Quadratic terms: (A_r, A_{r+1})
            H(i, i) += Kr(1, 1);
            H(j, j) += Kr(3, 3);
            H(i, j) += Kr(1, 3);
            H(j, i) += Kr(3, 1);

            // Linear terms (velocities + constants)
            // b_i += Kr[1,0]*V_r + Kr[1,2]*V_{r+1} + k_r[1]
            b.row(i) += (Kr(1, 0) * V[r].transpose()
                      +  Kr(1, 2) * V[r+1].transpose()
                      +  kr.row(1));
            // b_j += Kr[3,0]*V_r + Kr[3,2]*V_{r+1} + k_r[3]
            b.row(j) += (Kr(3, 0) * V[r].transpose()
                      +  Kr(3, 2) * V[r+1].transpose()
                      +  kr.row(3));
        }
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::solveMinJerkAccOnlyClosedForm(
    const std::vector<Vec3> &wps,       // size M+1
    const std::vector<double> &T,       // size M
    const std::vector<Vec3> &V,         // size M+1 (fixed)
    const Vec3 &a0, const Vec3 &af,     // fixed endpoint accelerations
    std::vector<Vec3> &A_out            // out: size M+1
)
{
    const int M = static_cast<int>(T.size());
    const int N = M - 1;

    A_out.clear();
    A_out.resize(M + 1);
    A_out.front() = a0;
    A_out.back()  = af;

    if (N <= 0) return;  // no interior unknowns

    Eigen::MatrixXd H, b;            // (N×N), (N×3)
    assemble_H_b_acc_only(wps, T, V, a0, af, H, b);

    // Solve H * X = -b, where X is (N×3) stacked [A1; ...; A_{M-1}]
    Eigen::MatrixXd X = H.ldlt().solve(-b);

    for (int i = 1; i <= M - 1; ++i) {
        A_out[i] = X.row(i - 1).transpose(); // Vec3
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::getGoalSetpoints(std::vector<state> &goal_setpoints)
{
    // 1) total time & sample count
    double total_time = std::accumulate(T_opt_.begin(), T_opt_.end(), 0.0);
    // BUG #1: truncating here can undershoot the end of the trajectory.
    // It’s usually better to do:
    int N = static_cast<int>(std::ceil(total_time / dc_));
    N = std::max(N, 2); // at least two points

    // 2) resize output
    goal_setpoints.resize(N);

    // 3) precompute timestamps
    std::vector<double> timestamps(N);
    for (int i = 0; i < N; ++i)
    {
        timestamps[i] = (i + 1) * dc_;
    }

    // 4) cumulative segment‑end times
    int M = static_cast<int>(T_opt_.size());
    std::vector<double> ends(M);
    if (M > 0)
    {
        ends[0] = T_opt_[0];
        for (int s = 1; s < M; ++s)
        {
            ends[s] = ends[s - 1] + T_opt_[s];
        }
    }

    // 5) parallel fill of goal_setpoints
    #pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        double ti = timestamps[i];

        // find which segment we're in
        auto it = std::upper_bound(ends.begin(), ends.end(), ti);
        int seg = static_cast<int>(it - ends.begin());
        seg = std::min(seg, std::max(0, M - 1));

        double t0 = (seg > 0 ? ends[seg - 1] : 0.0);
        double tau = (T_opt_[seg] > 0.0)
                         ? (ti - t0) / T_opt_[seg]
                         : 0.0;

        StateDeriv D = evalStateDeriv(seg, tau);

        state st;
        st.setTimeStamp(t0_ + ti);
        st.setPos(D.pos.x(), D.pos.y(), D.pos.z());
        st.setVel(D.vel.x(), D.vel.y(), D.vel.z());
        st.setAccel(D.accel.x(), D.accel.y(), D.accel.z());
        st.setJerk(D.jerk.x(), D.jerk.y(), D.jerk.z());
        goal_setpoints[i] = st;
    }

    // 6) zero out final derivatives explicitly
    {
        auto &last = goal_setpoints.back();
        // BUG #2: .transpose() is unnecessary if vel/accel are Vector3d
        last.vel = Eigen::Vector3d::Zero();
        last.accel = Eigen::Vector3d::Zero();
        last.jerk = Eigen::Vector3d::Zero();
    }
}

// -----------------------------------------------------------------------------

inline StateDeriv SolverMIGHTY::evalStateDeriv(int s, double tau) const
{
    // endpoints and duration
    const auto &P0 = P_opt_[s];
    const auto &V0 = V_opt_[s];
    const auto &A0 = A_opt_[s];
    const auto &P1 = P_opt_[s + 1];
    const auto &V1 = V_opt_[s + 1];
    const auto &A1 = A_opt_[s + 1];
    double T = T_opt_[s];
    double t2 = tau * tau, t3 = t2 * tau, t4 = t3 * tau, t5 = t4 * tau;

    // Quintic Hermite basis H_i(τ)
    double h0 = 1 - 10 * t3 + 15 * t4 - 6 * t5;
    double h1 = tau - 6 * t3 + 8 * t4 - 3 * t5;
    double h2 = 0.5 * (t2 - 3 * t3 + 3 * t4 - t5);
    double h3 = 10 * t3 - 15 * t4 + 6 * t5;
    double h4 = -4 * t3 + 7 * t4 - 3 * t5;
    double h5 = 0.5 * (t3 - 2 * t4 + t5);

    // position
    Eigen::Vector3d p =
        P0 * h0 + V0 * (h1 * T) + A0 * (h2 * T * T) + P1 * h3 + V1 * (h4 * T) + A1 * (h5 * T * T);

    // first derivatives of H_i(τ)
    double dh0 = -30 * t2 + 60 * t3 - 30 * t4;
    double dh1 = 1 - 18 * t2 + 32 * t3 - 15 * t4;
    double dh2 = 0.5 * (2 * tau - 9 * t2 + 12 * t3 - 5 * t4);
    double dh3 = 30 * t2 - 60 * t3 + 30 * t4;
    double dh4 = -12 * t2 + 28 * t3 - 15 * t4;
    double dh5 = 0.5 * (3 * t2 - 8 * t3 + 5 * t4);

    // velocity = (1/T) * d/dτ
    Eigen::Vector3d v = (P0 * dh0 + V0 * (dh1 * T) + A0 * (dh2 * T * T) + P1 * dh3 + V1 * (dh4 * T) + A1 * (dh5 * T * T)) / T;

    // second derivatives of H_i(τ)
    double d2h0 = -60 * tau + 180 * t2 - 120 * t3;
    double d2h1 = -36 * tau + 96 * t2 - 60 * t3;
    double d2h2 = 0.5 * (2 - 18 * tau + 36 * t2 - 20 * t3);
    double d2h3 = 60 * tau - 180 * t2 + 120 * t3;
    double d2h4 = -24 * tau + 84 * t2 - 60 * t3;
    double d2h5 = 0.5 * (6 * tau - 24 * t2 + 20 * t3);

    // acceleration = (1/T^2) * d²/dτ²
    Eigen::Vector3d a = (P0 * d2h0 + V0 * (d2h1 * T) + A0 * (d2h2 * T * T) + P1 * d2h3 + V1 * (d2h4 * T) + A1 * (d2h5 * T * T)) / (T * T);

    // third derivatives of H_i(τ)
    double d3h0 = -60 + 360 * tau - 360 * t2;
    double d3h1 = -36 + 192 * tau - 180 * t2;
    double d3h2 = 0.5 * (-18 + 72 * tau - 60 * t2);
    double d3h3 = 60 - 360 * tau + 360 * t2;
    double d3h4 = -24 + 168 * tau - 180 * t2;
    double d3h5 = 0.5 * (6 - 48 * tau + 60 * t2);

    // jerk = (1/T^3) * d³/dτ³
    Eigen::Vector3d j = (P0 * d3h0 + V0 * (d3h1 * T) + A0 * (d3h2 * T * T) + P1 * d3h3 + V1 * (d3h4 * T) + A1 * (d3h5 * T * T)) / (T * T * T);

    return {p, v, a, j};
}

double SolverMIGHTY::evaluateObjectiveAndGradient(
    const Eigen::VectorXd &z,
    Eigen::VectorXd &grad)
{
    // ---- Reconstruct once (shared by f & g) ----
    reconstruct_inplace(z);
    auto &P = P_s_;
    auto &V = V_s_;
    auto &A = A_s_;
    auto &CP = CP_s_;
    auto &T = T_s_;

    const int M = static_cast<int>(T.size());
    if (M == 0)
    {
        grad.setZero(K_);
        return 0.0;
    }
    const int knots = M + 1;

    // =========================================================================
    // Precompute shape bases (Bernstein) once for the chosen integral resolution
    // =========================================================================
    const int kappa = (integral_resolution_ > 0 ? integral_resolution_ : 30);
    ensureBernsteinCache(kappa); // fills bern_.B5, B4, B3, B2, W

    // =========================================================================
    // Cost accumulators (unweighted)
    // =========================================================================
    double J_time = 0.0; // Σ T
    double J_jerk = 0.0; // closed-form per segment
    double J_stat = 0.0; // corridor clearance
    double J_vel = 0.0;  // velocity limit
    double J_om = 0.0;   // body-rate limit
    double J_tilt = 0.0; // tilt limit
    double J_thr = 0.0;  // thrust ring

    for (double Ts : T)
        J_time += Ts;

    const double mu = (hinge_mu_ > 0.0 ? hinge_mu_ : 1e-2);
    const double Vmax2 = V_max_ * V_max_;
    const double Om2 = omege_max_ * omege_max_;

    // GCOPTER thrust ring parameters
    const double f_mean = 0.5 * (f_min_ + f_max_);
    const double f_radi = 0.5 * std::abs(f_max_ - f_min_);
    const double f_radi2 = f_radi * f_radi;

    // =========================================================================
    // Gradient accumulators in (P,V,A,T) space
    // =========================================================================
    std::vector<Vec3> gP(knots, Vec3::Zero());
    std::vector<Vec3> gV(knots, Vec3::Zero());
    std::vector<Vec3> gA(knots, Vec3::Zero());
    std::vector<double> gT(M, 0.0);

    Eigen::setNbThreads(1);

    // =========================================================================
    // Loop over segments
    // =========================================================================
    for (int s = 0; s < M; ++s)
    {
        const double Ts = T[s];
        const double invT = 1.0 / (Ts + 1e-16);
        const double invT2 = invT * invT;
        const double invT3 = invT2 * invT;
        const double invT4 = invT2 * invT2;

        // ---- (1) JERK closed-form (cost + grads wrt CP and T) ----
        Eigen::Vector3d gCP_jerk[6] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(),
                                       Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
        double gT_jerk = 0.0;
        jerkClosedFormForSegment(CP[s], Ts, J_jerk, gCP_jerk, gT_jerk, jerk_weight_);

        // ---- collect per-sample contributions (corridor + limits) ----
        Eigen::Vector3d gCP_samp[6] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(),
                                       Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
        double gT_samp = 0.0;

        // Precompute Bezier finite differences (shape-space)
        Vec3 D1[5], D2[4], D3[3];
        for (int j = 0; j < 5; ++j)
            D1[j] = CP[s][j + 1] - CP[s][j];
        for (int j = 0; j < 4; ++j)
            D2[j] = CP[s][j + 2] - 2.0 * CP[s][j + 1] + CP[s][j];
        for (int j = 0; j < 3; ++j)
            D3[j] = CP[s][j + 3] - 3.0 * CP[s][j + 2] + 3.0 * CP[s][j + 1] - CP[s][j];

        // Pack control points & finite diffs once per segment (3×N)
        Eigen::Matrix<double, 3, 6> Cseg;
        for (int c = 0; c < 6; ++c)
            Cseg.col(c) = CP[s][c];
        Eigen::Matrix<double, 3, 5> D1m;
        for (int c = 0; c < 5; ++c)
            D1m.col(c) = D1[c];
        Eigen::Matrix<double, 3, 4> D2m;
        for (int c = 0; c < 4; ++c)
            D2m.col(c) = D2[c];
        Eigen::Matrix<double, 3, 3> D3m;
        for (int c = 0; c < 3; ++c)
            D3m.col(c) = D3[c];

        // Static planes for this segment
        const auto &Aseg = A_stat_[s]; // (H x 3)
        const auto &bseg = b_stat_[s]; // (H)
        const bool has_planes = (Aseg.rows() > 0);

        // Sampling (trapezoid in scaled s)
        for (int j = 0; j <= kappa; ++j)
        {
            const double wj = bern_.W[j]; // 0.5 at ends, else 1
            const double tau = (kappa > 0 ? (double)j / (double)kappa : 0.0);
            const double dt = Ts / (double)kappa;
            const double wseg = wj * dt; // time weight per sample

            // Cached Bernstein basis values
            const auto &B5 = bern_.B5[j];
            const auto &B4 = bern_.B4[j];
            const auto &B3 = bern_.B3[j];
            const auto &B2 = bern_.B2[j];

            // Evaluate x(s), d1(s), d2(s), d3(s) in shape-space
            const Eigen::Map<const Eigen::Matrix<double, 6, 1>> b5(B5.data());
            const Eigen::Map<const Eigen::Matrix<double, 5, 1>> b4(B4.data());
            const Eigen::Map<const Eigen::Matrix<double, 4, 1>> b3(B3.data());
            const Eigen::Map<const Eigen::Matrix<double, 3, 1>> b2v(B2.data());

            // shape-space evals (no 5/20/60 here)
            Vec3 x = Cseg * b5;   // 3×6 · 6×1
            Vec3 d1s = D1m * b4;  // 3×5 · 5×1
            Vec3 d2s = D2m * b3;  // 3×4 · 4×1
            Vec3 d3s = D3m * b2v; // 3×3 · 3×1

            // Convert to physical derivatives
            const Vec3 v = (5.0 * invT) * d1s;
            const Vec3 acc = (20.0 * invT2) * d2s;
            const Vec3 jrk = (60.0 * invT3) * d3s;

            // ---- Static corridor clearance: penalize Co_ - (A x - b) ----
            if (has_planes)
            {
                for (int h = 0; h < Aseg.rows(); ++h)
                {
                    const double gval = Aseg.row(h).dot(x) - bseg[h]; // A x - b
                    const double viol = Co_ - gval;
                    const double phi = smoothed_l1(viol, mu);
                    J_stat += phi * wseg;

                    if (viol > 0.0)
                    {
                        const double dphi = smoothed_l1_prime(viol, mu);
                        const Vec3 gx = -(dphi)*Aseg.row(h).transpose() * (stat_weight_ * wseg);
                        for (int k = 0; k < 6; ++k)
                            gCP_samp[k] += B5[k] * gx;

                        // dt-only path: ∂(wseg)/∂T = wj/kappa
                        gT_samp += stat_weight_ * (wj / (double)kappa) * phi;
                    }
                }
            }

            // ---- Velocity limit (skip absolute endpoints) ----
            const bool is_global_start = (s == 0 && j == 0);
            const bool is_global_finish = (s == M - 1 && j == kappa);
            if (!is_global_start && !is_global_finish)
            {
                const double yv = v.squaredNorm() - Vmax2;
                const double phi = smoothed_l1(yv, mu);
                J_vel += phi * wseg;

                if (yv > 0.0)
                {
                    const double dphi = smoothed_l1_prime(yv, mu);
                    const Vec3 gv = (2.0 * dphi) * v * (dyn_constr_vel_weight_ * wseg); // ∂J/∂v

                    // v = (5/T) d1  ⇒ backprop to CP
                    const Vec3 g1 = gv * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP_samp[k] -= G1k;
                        gCP_samp[k + 1] += G1k;
                    }

                    // dt-only + T-scaling
                    gT_samp += dyn_constr_vel_weight_ * (wj / (double)kappa) * phi;
                    gT_samp += -5.0 * gv.dot(d1s) * invT2; // ∂v/∂T = -(5/T^2) d1
                }
            }

            // ---- Flatness: thr, quat, omg from (v, a, j) ----
            double thr;
            Eigen::Vector4d quat;
            Eigen::Vector3d omg;
            flatmap_.forward(v, acc, jrk, 0.0, 0.0, thr, quat, omg);

            // ---- Body-rate: y = ||ω||^2 - Ω^2 (C¹ hinge) ----
            {
                const double y = omg.squaredNorm() - Om2;
                const double phi = smoothed_l1(y, mu);
                J_om += phi * wseg;

                if (y > 0.0)
                {
                    const double dphi = smoothed_l1_prime(y, mu);

                    const Vec3 pos_grad = Vec3::Zero();
                    const Vec3 vel_grad = Vec3::Zero();
                    const double thr_grad = 0.0;
                    Eigen::Vector4d quat_grad = Eigen::Vector4d::Zero();
                    Vec3 omg_grad = (2.0 * dphi) * omg * (dyn_constr_bodyrate_weight_ * wseg);

                    Vec3 totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;
                    flatmap_.backward(pos_grad, vel_grad, thr_grad, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer, totalPsi, totalPsiD);

                    // Map (v,a,j) grads onto CP via d1s,d2s,d3s
                    const Vec3 g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP_samp[k] -= G1k;
                        gCP_samp[k + 1] += G1k;
                    }
                    const Vec3 g2 = totalAcc * invT2;
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3[k]) * g2;
                        gCP_samp[k] += G2k;
                        gCP_samp[k + 1] -= 2.0 * G2k;
                        gCP_samp[k + 2] += G2k;
                    }
                    const Vec3 g3 = totalJer * invT3;
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3;
                        gCP_samp[k] -= G3k;
                        gCP_samp[k + 1] += 3.0 * G3k;
                        gCP_samp[k + 2] -= 3.0 * G3k;
                        gCP_samp[k + 3] += G3k;
                    }

                    // dt-only + T-scale
                    gT_samp += dyn_constr_bodyrate_weight_ * (wj / (double)kappa) * phi;
                    gT_samp += -5.0 * totalVel.dot(d1s) * invT2;
                    gT_samp += -40.0 * totalAcc.dot(d2s) * invT3;
                    gT_samp += -180.0 * totalJer.dot(d3s) * invT4;
                }
            }

            // ---- Tilt: θ = acos(R33) with R33 = q0^2 - qx^2 - qy^2 + qz^2 ----
            {
                const double q0 = quat(0), qx = quat(1), qy = quat(2), qz = quat(3);
                const double R33 = q0 * q0 - qx * qx - qy * qy + qz * qz;

                double theta, dtheta_dR33;
                safe_acos_with_grad(R33, theta, dtheta_dR33);

                const double y = theta - tilt_max_rad_;
                const double phi = smoothed_l1(y, mu);
                J_tilt += phi * wseg;

                if (y > 0.0)
                {
                    const double dphi = smoothed_l1_prime(y, mu);

                    // ∂R33/∂q = [2q0, -2qx, -2qy, 2qz]
                    Eigen::Vector4d dR33_dq;
                    dR33_dq << 2.0 * q0, -2.0 * qx, -2.0 * qy, 2.0 * qz;

                    // raw ∂J/∂q
                    Eigen::Vector4d quat_grad = (dphi * dtheta_dR33) * dR33_dq;

                    // project to tangent of unit quaternion (guards renormalization in forward)
                    const double qn2 = quat.squaredNorm();
                    if (qn2 > 1e-16)
                        quat_grad -= (quat.dot(quat_grad) / qn2) * quat;

                    quat_grad *= (dyn_constr_tilt_weight_ * wseg);

                    // Backprop through flatness
                    const Vec3 pos_grad = Vec3::Zero(), vel_grad = Vec3::Zero(), omg_grad = Vec3::Zero();
                    const double thr_grad = 0.0;
                    Vec3 totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;
                    flatmap_.backward(pos_grad, vel_grad, thr_grad, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer, totalPsi, totalPsiD);

                    // Map (v,a,j) grads onto CP
                    const Vec3 g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP_samp[k] -= G1k;
                        gCP_samp[k + 1] += G1k;
                    }
                    const Vec3 g2 = totalAcc * invT2;
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3[k]) * g2;
                        gCP_samp[k] += G2k;
                        gCP_samp[k + 1] -= 2.0 * G2k;
                        gCP_samp[k + 2] += G2k;
                    }
                    const Vec3 g3 = totalJer * invT3;
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3;
                        gCP_samp[k] -= G3k;
                        gCP_samp[k + 1] += 3.0 * G3k;
                        gCP_samp[k + 2] -= 3.0 * G3k;
                        gCP_samp[k + 3] += G3k;
                    }

                    // T paths
                    gT_samp += dyn_constr_tilt_weight_ * (wj / (double)kappa) * phi;
                    gT_samp += -5.0 * totalVel.dot(d1s) * invT2;
                    gT_samp += -40.0 * totalAcc.dot(d2s) * invT3;
                    gT_samp += -180.0 * totalJer.dot(d3s) * invT4;
                }
            }

            // ---- Thrust ring: ((thr - f_mean)^2 - f_radi^2) ----
            {
                const double df = thr - f_mean;
                const double y = df * df - f_radi2;
                const double phi = smoothed_l1(y, mu);
                J_thr += phi * wseg;

                if (y > 0.0)
                {
                    const double dphi = smoothed_l1_prime(y, mu);
                    const double dJ_dthr = 2.0 * df * dphi * (dyn_constr_thrust_weight_ * wseg);

                    const Vec3 pos_grad = Vec3::Zero(), vel_grad = Vec3::Zero(), omg_grad = Vec3::Zero();
                    const Eigen::Vector4d quat_grad = Eigen::Vector4d::Zero();
                    Vec3 totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;

                    flatmap_.backward(pos_grad, vel_grad, dJ_dthr, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer, totalPsi, totalPsiD);

                    // Map (v,a,j) grads onto CP
                    const Vec3 g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP_samp[k] -= G1k;
                        gCP_samp[k + 1] += G1k;
                    }
                    const Vec3 g2 = totalAcc * invT2;
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3[k]) * g2;
                        gCP_samp[k] += G2k;
                        gCP_samp[k + 1] -= 2.0 * G2k;
                        gCP_samp[k + 2] += G2k;
                    }
                    const Vec3 g3 = totalJer * invT3;
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3;
                        gCP_samp[k] -= G3k;
                        gCP_samp[k + 1] += 3.0 * G3k;
                        gCP_samp[k + 2] -= 3.0 * G3k;
                        gCP_samp[k + 3] += G3k;
                    }

                    // T paths
                    gT_samp += dyn_constr_thrust_weight_ * (wj / (double)kappa) * phi;
                    gT_samp += -5.0 * totalVel.dot(d1s) * invT2;
                    gT_samp += -40.0 * totalAcc.dot(d2s) * invT3;
                    gT_samp += -180.0 * totalJer.dot(d3s) * invT4;
                }
            }
        } // end samples

        // ---- (2) Push both jerk+sample CP grads into (P,V,A) and T ----
        Eigen::Vector3d gCP_all[6];
        for (int i = 0; i < 6; ++i)
            gCP_all[i] = gCP_jerk[i] + gCP_samp[i];

        pushCPGradToHermite(
            gCP_all, Ts,
            V[s], A[s], V[s + 1], A[s + 1],
            gP[s], gV[s], gA[s], gP[s + 1], gV[s + 1], gA[s + 1],
            gT[s] // += CP(T) path
        );
        gT[s] += gT_jerk + gT_samp; // add factor paths

    } // segments

    // =========================================================================
    // Scatter once into 'grad' and do T̄ coupling + τ chain + time term
    // =========================================================================
    grad.setZero(K_);

    // positions
    scatterPosGrads(gP, z, grad);

    // derivatives blocks (INTERIOR only)
    if (scale_derivs_)
    {
        std::vector<double> Tbar;
        build_Tbar(T, Tbar);

        for (int i = 1; i < M_; ++i)
        {
            const double Tb = std::max(1e-12, Tbar[i]);
            const int voff = vhatOffset(i);
            const int aoff = ahatOffset(i);
            grad.segment<3>(voff) += gV[i] / Tb;
            grad.segment<3>(aoff) += gA[i] / (Tb * Tb);
        }

        // T̄ coupling
        std::vector<double> gTbar(M_ + 1, 0.0);
        for (int i = 1; i < M_; ++i)
        {
            const double Tb = std::max(1e-12, Tbar[i]);
            const int voff = vhatOffset(i);
            const int aoff = ahatOffset(i);
            Eigen::Map<const Vec3> vhat_i(z.data() + voff);
            Eigen::Map<const Vec3> ahat_i(z.data() + aoff);
            gTbar[i] += -gV[i].dot(vhat_i) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_i) / (Tb * Tb * Tb);
        }

        std::vector<double> gTextra(M, 0.0);
        distribute_gTbar_to_segments(gTbar, gTextra);
        for (int s = 0; s < M_; ++s)
        {
            gT[s] += gTextra[s]; // <- add only in scaled mode
        }
    }
    else
    {
        // Unscaled benchmark: no T̄ path at all
        for (int i = 1; i < M_; ++i)
        {
            const int voff = vhatOffset(i);
            const int aoff = ahatOffset(i);
            grad.segment<3>(voff) += gV[i];
            grad.segment<3>(aoff) += gA[i];
        }
    }

    // chain ∂J/∂T_s → ∂J/∂τ_s (common to both modes)
    for (int s = 0; s < M_; ++s)
        grad[offTau_ + s] += gT[s] * dT_dtau(z[offTau_ + s]);

    // time regularization (common)
    for (int s = 0; s < M_; ++s)
        grad[offTau_ + s] += time_weight_ * dT_dtau(z[offTau_ + s]);

    // =========================================================================
    // Final weighted objective (matches evaluateObjective)
    // =========================================================================
    const double f =
        time_weight_ * J_time +
        jerk_weight_ * J_jerk +
        stat_weight_ * J_stat +
        dyn_constr_vel_weight_ * J_vel +
        dyn_constr_bodyrate_weight_ * J_om +
        dyn_constr_tilt_weight_ * J_tilt +
        dyn_constr_thrust_weight_ * J_thr;

    return f;
}

double SolverMIGHTY::evaluateObjective(const VecXd &z)
{
    // Reconstruct
    std::vector<Vec3> P, V, A;
    std::vector<std::array<Vec3, 6>> CP;
    std::vector<double> T;
    reconstruct(z, P, V, A, CP, T);
    const int M = static_cast<int>(T.size());
    if (M == 0)
        return 0.0;

    // ---- 1) time ----
    double J_time = 0.0;
    for (double Ts : T)
        J_time += Ts;

    // ---- 2) jerk (closed form per segment) ----
    double J_jerk = 0.0;
    for (int s = 0; s < M; ++s)
    {
        const double Ts = T[s];
        const double C = 3600.0 / std::pow(Ts, 5);
        const Vec3 d30 = CP[s][3] - 3.0 * CP[s][2] + 3.0 * CP[s][1] - CP[s][0];
        const Vec3 d31 = CP[s][4] - 3.0 * CP[s][3] + 3.0 * CP[s][2] - CP[s][1];
        const Vec3 d32 = CP[s][5] - 3.0 * CP[s][4] + 3.0 * CP[s][3] - CP[s][2];
        J_jerk += C * (d30.squaredNorm() + d31.squaredNorm() + d32.squaredNorm());
    }

    // ---- 3) sampled terms in s ∈ [0,1] with dt = T_s / kappa (trapezoid) ----
    const int kappa = (integral_resolution_ > 0 ? integral_resolution_ : 30);
    if (kappa <= 0)
        return time_weight_ * J_time + jerk_weight_ * J_jerk;

    const double mu = (hinge_mu_ > 0.0 ? hinge_mu_ : 1e-2);
    const double Vmax2 = V_max_ * V_max_;
    const double Om2 = omege_max_ * omege_max_;
    const double m = mass_;
    const double g = g_;
    const double eps = 1e-12;
    const Vec3 e3(0.0, 0.0, 1.0);

    // GCOPTER thrust ring parameters
    const double f_mean = 0.5 * (f_min_ + f_max_);
    const double f_radi = 0.5 * std::abs(f_max_ - f_min_);
    const double f_radi2 = f_radi * f_radi;

    // new accumulators
    double J_stat = 0.0, J_vel = 0.0, J_om = 0.0, J_tilt = 0.0, J_thr = 0.0, J_om_smooth = 0.0, J_tilt_bias = 0.0, J_Tfloor = 0.0;

    for (int s = 0; s < M; ++s)
    {
        const double Ts = T[s];
        const double dt = Ts / static_cast<double>(kappa);

        // Precompute finite differences of CP (Bezier derivatives in s-space)
        Vec3 D1[5], D2[4], D3[3];
        for (int j = 0; j < 5; ++j)
            D1[j] = CP[s][j + 1] - CP[s][j];
        for (int j = 0; j < 4; ++j)
            D2[j] = CP[s][j + 2] - 2.0 * CP[s][j + 1] + CP[s][j];
        for (int j = 0; j < 3; ++j)
            D3[j] = CP[s][j + 3] - 3.0 * CP[s][j + 2] + 3.0 * CP[s][j + 1] - CP[s][j];

        // Static corridor planes (can be empty)
        const auto &Aseg = A_stat_[s]; // (H x 3)
        const auto &bseg = b_stat_[s]; // (H)
        const bool has_planes = (Aseg.rows() > 0);

        for (int j = 0; j <= kappa; ++j)
        {
            const double wj = (j == 0 || j == kappa) ? 0.5 : 1.0;
            const double tau = static_cast<double>(j) / static_cast<double>(kappa);

            double B5[6], B4[5], B3b[4], B2[3];
            bernstein5(tau, B5);
            bernstein4(tau, B4);
            bernstein3(tau, B3b);
            bernstein2(tau, B2);

            // x(s), dx/ds, d2x/ds2, d3x/ds3
            Vec3 x = Vec3::Zero(), dxs = Vec3::Zero(), d2s = Vec3::Zero(), d3s = Vec3::Zero();
            for (int k = 0; k < 6; ++k)
                x += B5[k] * CP[s][k];
            for (int k = 0; k < 5; ++k)
                dxs += B4[k] * D1[k];
            for (int k = 0; k < 4; ++k)
                d2s += B3b[k] * D2[k];
            for (int k = 0; k < 3; ++k)
                d3s += B2[k] * D3[k];

            // convert to t-derivatives: v = (5/T)*dxs, a = (20/T^2)*d2s, j = (60/T^3)*d3s
            const double invT = 1.0 / (Ts + 1e-16);
            const Vec3 v = (5.0 * invT) * dxs;
            const Vec3 a = (20.0 * invT * invT) * d2s;
            const Vec3 jrk = (60.0 * invT * invT * invT) * d3s;

            const double wseg = wj * dt;

            // --- Static corridor: penalize when A x - b < Co_  (same as Python)
            if (has_planes)
            {
                for (int h = 0; h < Aseg.rows(); ++h)
                {
                    const double gval = Aseg.row(h).dot(x) - bseg[h];
                    const double viol = Co_ - gval;         // <-- key: use Co_
                    J_stat += smoothed_l1(viol, mu) * wseg; // smoothed_l1 zeros negatives
                }
            }

            // --- Velocity: y = ||v||^2 - Vmax^2
            const bool is_global_start = (s == 0 && j == 0);
            const bool is_global_finish = (s == M - 1 && j == kappa);

            if (!is_global_start && !is_global_finish)
            {
                const double yv = v.squaredNorm() - Vmax2;
                J_vel += smoothed_l1(yv, mu) * wseg;
            }

            // Get thr, quat, omg
            double thr;
            Eigen::Vector4d quat;
            Eigen::Vector3d omg;
            flatmap_.forward(v, a, jrk, 0.0, 0.0, thr, quat, omg);

            // --- Body-rate from (a, j) with yaẇ=0
            const double viol_omg = omg.squaredNorm() - Om2;
            J_om += smoothed_l1(viol_omg, mu) * wseg;

            // --- Tilt: y = cosθ_max - cosθ, cosθ = e3·b3
            const double cos_theta = quat(0) * quat(0) - quat(1) * quat(1) - quat(2) * quat(2) + quat(3) * quat(3);
            const double viol_theta = acos(cos_theta) - tilt_max_rad_;
            J_tilt += smoothed_l1(viol_theta, mu) * wseg;

            // --- Thrust ring (GCOPTER): y = ((f - f_mean)^2 - f_radi^2)
            const double viol_thr = (thr - f_mean) * (thr - f_mean) - f_radi2;
            J_thr += smoothed_l1(viol_thr, mu) * wseg;
        }
    }

    // final weighted sum (add the 4 new terms)
    return time_weight_ * J_time + jerk_weight_ * J_jerk + stat_weight_ * J_stat + dyn_constr_vel_weight_ * J_vel + dyn_constr_bodyrate_weight_ * J_om + dyn_constr_tilt_weight_ * J_tilt + dyn_constr_thrust_weight_ * J_thr;
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::computeAnalyticalGrad(const Eigen::VectorXd &z, Eigen::VectorXd &grad)
{
    // Reconstruct (used by helpers)
    std::vector<Vec3> P, V, A;
    std::vector<std::array<Vec3, 6>> CP;
    std::vector<double> T;
    reconstruct(z, P, V, A, CP, T);

    grad.setZero(K_);

    // 1) time term: ∂J_time/∂τ_s = w_time * dT/dτ_s
    for (int s = 0; s < M_; ++s)
    {
        const double dT = dT_dtau(z[offTau_ + s]);
        grad[offTau_ + s] += time_weight_ * dT;
    }

    // 2) jerk piece
    {
        Eigen::VectorXd g(K_);
        g.setZero();
        dJ_jerk_dz(z, P, V, A, CP, T, g);
        grad += jerk_weight_ * g;
    }

    // 3) sampled static + dynamic limits (vel, body-rate, tilt, thrust)
    {
        Eigen::VectorXd g(K_);
        g.setZero();
        dJ_limits_and_static_dz(z, P, V, A, CP, T, g); // thrust update is inside this helper
        grad += g;                                     // weights applied per-term inside the helper
    }
}

// -----------------------------------------------------------------------------

int SolverMIGHTY::progressCallback(void *instance,
                                  const Eigen::VectorXd &z,
                                  const Eigen::VectorXd &g,
                                  const double f,
                                  const double /*step*/,
                                  const int /*k*/,
                                  const int /*ls*/)
{
    auto *self = static_cast<SolverMIGHTY *>(instance);

    // Remember the latest iterate and objective (so we can return it if we abort)
    self->last_z_ = z;
    self->last_f_ = f;

    if (self->have_deadline_)
    {
        const auto now = std::chrono::steady_clock::now();
        if (now >= self->opt_deadline_)
        {
            self->timed_out_ = true;
            return 1; // non-zero => abort now
        }
    }

    return 0; // continue
}

// -----------------------------------------------------------------------------

double SolverMIGHTY::evalObjGradCallback(
    void *instance,
    const Eigen::VectorXd &x,
    Eigen::VectorXd &g)
{
    return static_cast<SolverMIGHTY *>(instance)->evaluateObjectiveAndGradient(x, g);
}

// -----------------------------------------------------------------------------

int SolverMIGHTY::optimize(
    const Eigen::VectorXd &z0,
    Eigen::VectorXd &z_opt,
    double &f_opt,
    const lbfgs::lbfgs_parameter_t &param) const
{
    // copy initial guess
    Eigen::VectorXd z = z0;

    int status = lbfgs::lbfgs_optimize(
        z,                                 // in/out decision vector
        f_opt,                             // out final cost
        &SolverMIGHTY::evalObjGradCallback, // objective+gradient callback
        /*stepbound*/ nullptr,
        // /*progress*/ &SolverMIGHTY::progressCallback,
        /*progress*/ nullptr,
        const_cast<SolverMIGHTY *>(this),
        param);

    // copy result
    z_opt = z;
    return status;
}

// -----------------------------------------------------------------------------

int SolverMIGHTY::optimize(const Eigen::VectorXd &z0,
                          Eigen::VectorXd &z_opt,
                          double &f_opt,
                          const lbfgs::lbfgs_parameter_t &param,
                          std::chrono::milliseconds wall_budget) const
{
    // copy initial guess
    Eigen::VectorXd z = z0;

    // Set up deadline and reset bookkeeping
    opt_start_ = std::chrono::steady_clock::now();
    opt_deadline_ = opt_start_ + wall_budget;
    have_deadline_ = (wall_budget.count() > 0);
    timed_out_ = false;

    last_z_ = z0;
    last_f_ = std::numeric_limits<double>::infinity();

    // IMPORTANT: enable the progress callback (not nullptr)
    int status = lbfgs::lbfgs_optimize(
        z,                                 // in/out decision vector
        f_opt,                             // out final cost
        &SolverMIGHTY::evalObjGradCallback,
        /*stepbound*/ nullptr,
        /*progress*/ &SolverMIGHTY::progressCallback, // <- enable it
        const_cast<SolverMIGHTY *>(this),
        param);

    // If we aborted due to time, hand back the latest/best iterate we kept.
    if (timed_out_)
    {
        z_opt = last_z_;
        // If the library didn't pass us f at every step, f_opt may be stale.
        // Use last_f_ when it is finite; otherwise keep f_opt.
        if (std::isfinite(last_f_))
            f_opt = last_f_;
        return status;
    }

    // Normal termination
    z_opt = z;
    return status;
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::setStaticConstraints(
    const std::vector<LinearConstraint3D> &cons)
{
    A_stat_.clear();
    b_stat_.clear();
    A_stat_.reserve(cons.size());
    b_stat_.reserve(cons.size());

    for (auto const &lc : cons)
    {
        // lc.A() is Eigen::Matrix<decimal_t,Dynamic,3>
        // lc.b() is Eigen::Matrix<decimal_t,Dynamic,1>
        // cast to double if decimal_t isn’t already double:
        Eigen::Matrix<double, Eigen::Dynamic, 3> A =
            lc.A().template cast<double>();
        Eigen::VectorXd b =
            lc.b().template cast<double>();

        A_stat_.push_back(std::move(A));
        b_stat_.push_back(std::move(b));
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::dJ_jerk_dz(const VecXd &z,
                             const std::vector<Vec3> &P,
                             const std::vector<Vec3> &V,
                             const std::vector<Vec3> &A,
                             const std::vector<std::array<Vec3, 6>> &CP,
                             const std::vector<double> &T,
                             VecXd &grad) const
{
    const int M = static_cast<int>(T.size());
    const int knots = M + 1;

    // --- accumulate in (P,V,A,T) space ---
    std::vector<Vec3> gP(knots, Vec3::Zero());
    std::vector<Vec3> gV(knots, Vec3::Zero());
    std::vector<Vec3> gA(knots, Vec3::Zero());
    std::vector<double> gT(M, 0.0);

    for (int s = 0; s < M; ++s)
    {
        const double Ts = T[s];
        const double T2 = Ts * Ts;
        const double C = 3600.0 / std::pow(Ts, 5);

        const Vec3 d30 = CP[s][3] - 3.0 * CP[s][2] + 3.0 * CP[s][1] - CP[s][0];
        const Vec3 d31 = CP[s][4] - 3.0 * CP[s][3] + 3.0 * CP[s][2] - CP[s][1];
        const Vec3 d32 = CP[s][5] - 3.0 * CP[s][4] + 3.0 * CP[s][3] - CP[s][2];

        Vec3 dJdCP[6] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(),
                         Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
        auto add = [&](int i, const Vec3 &v)
        { dJdCP[i] += 2.0 * C * v; };
        // m = 0
        add(0, -d30);
        add(1, 3.0 * d30);
        add(2, -3.0 * d30);
        add(3, d30);
        // m = 1
        add(1, -d31);
        add(2, 3.0 * d31);
        add(3, -3.0 * d31);
        add(4, d31);
        // m = 2
        add(2, -d32);
        add(3, 3.0 * d32);
        add(4, -3.0 * d32);
        add(5, d32);

        // push to (P,V,A)
        gP[s] += dJdCP[0] + dJdCP[1] + dJdCP[2];
        gV[s] += (Ts / 5.0) * dJdCP[1] + (2.0 * Ts / 5.0) * dJdCP[2];
        gA[s] += (T2 / 20.0) * dJdCP[2];
        gP[s + 1] += dJdCP[3] + dJdCP[4] + dJdCP[5];
        gV[s + 1] += (-2.0 * Ts / 5.0) * dJdCP[3] + (-Ts / 5.0) * dJdCP[4];
        gA[s + 1] += (T2 / 20.0) * dJdCP[3];

        // explicit ∂/∂T from CP(T)
        gT[s] += dJdCP[1].dot(V[s] / 5.0);
        gT[s] += dJdCP[2].dot((2.0 / 5.0) * V[s] + (Ts / 10.0) * A[s]);
        gT[s] += dJdCP[3].dot((-2.0 / 5.0) * V[s + 1] + (Ts / 10.0) * A[s + 1]);
        gT[s] += dJdCP[4].dot((-1.0 / 5.0) * V[s + 1]);

        // explicit ∂/∂T of 1/T^5
        const double S = d30.squaredNorm() + d31.squaredNorm() + d32.squaredNorm();
        gT[s] += (-5.0) * (3600.0 / std::pow(Ts, 6)) * S;
    }

    // --- scatter (P,V,A) once ---
    scatterPosGrads(gP, z, grad);

    std::vector<double> Tbar;
    build_Tbar(T, Tbar);
    for (int i = 1; i < M_; ++i)
    {
        const double Tb = std::max(1e-12, Tbar[i]);
        grad.segment<3>(offVhat_ + 3 * i) += gV[i] / Tb;
        grad.segment<3>(offAhat_ + 3 * i) += gA[i] / (Tb * Tb);
    }

    // Keep T̄ coupling interior-only here as well:
    std::vector<double> gTbar(M + 1, 0.0);
    for (int i = 1; i < M_; ++i)
    {
        const double Tb = std::max(1e-12, Tbar[i]);
        const Eigen::Vector3d vhat = getVhat(z, i);
        const Eigen::Vector3d ahat = getAhat(z, i);
        gTbar[i] += -gV[i].dot(vhat) / (Tb * Tb) - 2.0 * gA[i].dot(ahat) / (Tb * Tb * Tb);
    }

    // --- add T̄→{T_s} to segment ∂J/∂T and chain to τ ---
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
    {
        gT[s] += gTextra[s];
        grad[offTau_ + s] += gT[s] * dT_dtau(z[offTau_ + s]);
    }
}

// -----------------------------------------------------------------------------

void SolverMIGHTY::dJ_limits_and_static_dz(const VecXd &z,
                                          const std::vector<Vec3> &P,
                                          const std::vector<Vec3> &V,
                                          const std::vector<Vec3> &A,
                                          const std::vector<std::array<Vec3, 6>> &CP,
                                          const std::vector<double> &T,
                                          VecXd &grad)
{
    const int M = static_cast<int>(T.size());
    const int knots = M + 1;

    // accumulators in (P,V,A,T) space
    std::vector<Vec3> gP(knots, Vec3::Zero());
    std::vector<Vec3> gV(knots, Vec3::Zero());
    std::vector<Vec3> gA(knots, Vec3::Zero());
    std::vector<double> gT(M, 0.0);

    const int kappa = (integral_resolution_ > 0 ? integral_resolution_ : 30);
    if (kappa <= 0)
        return;

    const double mu = (hinge_mu_ > 0.0 ? hinge_mu_ : 1e-2);
    const double Vmax2 = V_max_ * V_max_;
    const double Om2 = omege_max_ * omege_max_;
    const double cos_tilt_max = std::cos(tilt_max_rad_);
    const double m = mass_, g = g_;
    const Vec3 e3(0.0, 0.0, 1.0);

    // GCOPTER thrust ring parameters (mean + radius)
    const double f_mean = 0.5 * (f_min_ + f_max_);
    const double f_radi = 0.5 * std::abs(f_max_ - f_min_);
    const double f_radi2 = f_radi * f_radi;

    for (int s = 0; s < M; ++s)
    {
        const double Ts = T[s];
        const double dt = Ts / static_cast<double>(kappa);
        const double invT = 1.0 / (Ts + 1e-16);

        // Precompute Bézier differences (shape-space, no 5/20/60 factors here)
        Vec3 D1[5], D2[4], D3[3];
        for (int j = 0; j < 5; ++j)
            D1[j] = CP[s][j + 1] - CP[s][j];
        for (int j = 0; j < 4; ++j)
            D2[j] = CP[s][j + 2] - 2.0 * CP[s][j + 1] + CP[s][j];
        for (int j = 0; j < 3; ++j)
            D3[j] = CP[s][j + 3] - 3.0 * CP[s][j + 2] + 3.0 * CP[s][j + 1] - CP[s][j];

        // Static planes (can be empty)
        const auto &Aseg = A_stat_[s];
        const auto &bseg = b_stat_[s];
        const bool has_planes = (Aseg.rows() > 0);

        // per-segment CP gradient bucket
        Vec3 gCP[6] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(),
                       Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};

        for (int j = 0; j <= kappa; ++j)
        {
            const double wj = (j == 0 || j == kappa) ? 0.5 : 1.0;
            const double tau = (kappa == 0 ? 0.0 : double(j) / double(kappa));

            double B5[6], B4[5], B3b[4], B2[3];
            bernstein5(tau, B5);
            bernstein4(tau, B4);
            bernstein3(tau, B3b);
            bernstein2(tau, B2);

            // shape-space derivatives (no factors)
            Vec3 x = Vec3::Zero(), d1s = Vec3::Zero(), d2s = Vec3::Zero(), d3s = Vec3::Zero();
            for (int k = 0; k < 6; ++k)
                x += B5[k] * CP[s][k];
            for (int k = 0; k < 5; ++k)
                d1s += B4[k] * D1[k];
            for (int k = 0; k < 4; ++k)
                d2s += B3b[k] * D2[k];
            for (int k = 0; k < 3; ++k)
                d3s += B2[k] * D3[k];

            // physical derivatives
            const Vec3 v = (5.0 * invT) * d1s;                  // v = (5/T) d1
            const Vec3 acc = (20.0 * invT * invT) * d2s;        // a = (20/T^2) d2
            const Vec3 jrk = (60.0 * invT * invT * invT) * d3s; // j = (60/T^3) d3

            const double wseg = wj * dt;

            // --- Static corridor: Ax - b >= 0  -> penalize (-g)_+ ---
            if (has_planes)
            {
                for (int h = 0; h < Aseg.rows(); ++h)
                {
                    const double gval = Aseg.row(h).dot(x) - bseg[h];
                    if (gval < Co_) // match Python: penalize when g < Co
                    {
                        const double viol = Co_ - gval; // >0 inside the clearance band
                        const double phi_p = smoothed_l1_prime(viol, mu);

                        // ∂J/∂x = -φ'(viol) * Aᵀ * (w_stat * wseg)
                        const Vec3 gx = -(phi_p)*Aseg.row(h).transpose() * (stat_weight_ * wseg);

                        // x = Σ B5 * CP  ⇒  ∂J/∂CP_j += B5_j * gx
                        for (int k = 0; k < 6; ++k)
                            gCP[k] += B5[k] * gx;

                        // dt-only path: ∂(wseg)/∂T = wj / kappa
                        gT[s] += stat_weight_ * (wj / (double)kappa) * smoothed_l1(viol, mu);
                    }
                }
            }

            // --- Velocity: y = ||v||^2 - Vmax^2 ---
            const bool is_global_start = (s == 0 && j == 0);
            const bool is_global_finish = (s == M - 1 && j == kappa);

            if (!is_global_start && !is_global_finish)
            {
                const double yv = v.squaredNorm() - Vmax2;
                if (yv > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yv, mu);
                    const Vec3 gv = (2.0 * phi_p) * v * (dyn_constr_vel_weight_ * wseg); // ∂J/∂v

                    // v=(5/T)d1 ⇒ backprop to CP
                    const Vec3 g1 = gv * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP[k] -= G1k;
                        gCP[k + 1] += G1k;
                    }

                    // dt-only + T-scale path (for this term only)
                    gT[s] += dyn_constr_vel_weight_ * (wj / (double)kappa) * smoothed_l1(yv, mu);
                    gT[s] += -5.0 * gv.dot(d1s) / (Ts * Ts); // ∂v/∂T = -(5/T^2)d1
                }
            }

            // Common pieces for a/j
            double thr;
            Eigen::Vector4d quat;
            Eigen::Vector3d omg;
            flatmap_.forward(v, acc, jrk, 0.0, 0.0, thr, quat, omg);

            // --- Body-rate via flatness map: y = ||ω||^2 - Ω^2  -------------------------
            // ---- Body-rate: y = ||ω||^2 - Ω^2 (use C¹ hinge to match objective/fused) ----
            {
                const double y = omg.squaredNorm() - Om2;
                const double phi = smoothed_l1(y, mu);
                if (y > 0.0)
                {
                    const double dphi = smoothed_l1_prime(y, mu);

                    const Vec3 pos_grad = Vec3::Zero();
                    const Vec3 vel_grad = Vec3::Zero();
                    const double thr_grad = 0.0;
                    Eigen::Vector4d quat_grad = Eigen::Vector4d::Zero();
                    Vec3 omg_grad = (2.0 * dphi) * omg * (dyn_constr_bodyrate_weight_ * wseg);

                    Vec3 totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;
                    flatmap_.backward(pos_grad, vel_grad, thr_grad, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer,
                                      totalPsi, totalPsiD);

                    // v=(5/T)d1, a=(20/T^2)d2, j=(60/T^3)d3
                    const Vec3 g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP[k] -= G1k;
                        gCP[k + 1] += G1k;
                    }
                    const Vec3 g2_d2s = totalAcc / (Ts * Ts);
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * g2_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }
                    const Vec3 g3_d3s = totalJer / (Ts * Ts * Ts);
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3_d3s;
                        gCP[k] -= G3k;
                        gCP[k + 1] += 3.0 * G3k;
                        gCP[k + 2] -= 3.0 * G3k;
                        gCP[k + 3] += G3k;
                    }

                    // dt-only + T-scale paths
                    gT[s] += dyn_constr_bodyrate_weight_ * (wj / (double)kappa) * phi;
                    gT[s] += -5.0 * totalVel.dot(d1s) / (Ts * Ts);
                    gT[s] += -40.0 * totalAcc.dot(d2s) / (Ts * Ts * Ts);
                    gT[s] += -180.0 * totalJer.dot(d3s) / (Ts * Ts * Ts * Ts);
                }
            }

            // --- Tilt via quaternion: θ = acos(cosθ),  cosθ = 1 - 2(qx^2 + qy^2) --------
            {
                // components of the quaternion returned by flatmap_.forward
                const double q0 = quat(0), qx = quat(1), qy = quat(2), qz = quat(3);

                // R33 = q0^2 - qx^2 - qy^2 + qz^2 (exact for unit quaternions)
                const double R33 = q0 * q0 - qx * qx - qy * qy + qz * qz;

                double theta, dtheta_dR33;
                safe_acos_with_grad(R33, theta, dtheta_dR33); // sets derivative to 0 at clamp

                const double y = theta - tilt_max_rad_;
                const double phi = smoothed_l1(y, mu);

                if (y > 0.0) // active hinge
                {
                    const double dphi = smoothed_l1_prime(y, mu); // dJ/dθ (without weights)

                    // ∂R33/∂q = [2 q0, -2 qx, -2 qy,  2 qz]
                    Eigen::Vector4d dR33_dq;
                    dR33_dq << 2.0 * q0, -2.0 * qx, -2.0 * qy, 2.0 * qz;

                    // raw ∂J/∂q from chain θ(R33(q))
                    Eigen::Vector4d quat_grad = (dphi * dtheta_dR33) * dR33_dq;

                    // OPTIONAL but recommended: project to tangent space of unit quaternion
                    // to match forward's normalization (invariant to scaling/sign flips)
                    const double qnorm2 = quat.squaredNorm();
                    if (qnorm2 > 1e-16)
                        quat_grad -= (quat.dot(quat_grad) / qnorm2) * quat;

                    // weight + time quadrature
                    quat_grad *= (dyn_constr_tilt_weight_ * wseg);

                    // Backprop through the flatness map to (v,a,j)
                    const Eigen::Vector3d pos_grad = Eigen::Vector3d::Zero();
                    const Eigen::Vector3d vel_grad = Eigen::Vector3d::Zero();
                    const Eigen::Vector3d omg_grad = Eigen::Vector3d::Zero();
                    const double thr_grad = 0.0;

                    Eigen::Vector3d totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;
                    flatmap_.backward(pos_grad, vel_grad, thr_grad, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer, totalPsi, totalPsiD);

                    // push to CP: v=(5/T)d1, a=(20/T^2)d2, j=(60/T^3)d3
                    const Eigen::Vector3d g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Eigen::Vector3d G1k = (5.0 * B4[k]) * g1;
                        gCP[k] -= G1k;
                        gCP[k + 1] += G1k;
                    }
                    const Eigen::Vector3d g2_d2s = totalAcc / (Ts * Ts);
                    for (int k = 0; k < 4; ++k)
                    {
                        const Eigen::Vector3d G2k = (20.0 * B3b[k]) * g2_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }
                    const Eigen::Vector3d g3_d3s = totalJer / (Ts * Ts * Ts);
                    for (int k = 0; k < 3; ++k)
                    {
                        const Eigen::Vector3d G3k = (60.0 * B2[k]) * g3_d3s;
                        gCP[k] -= G3k;
                        gCP[k + 1] += 3.0 * G3k;
                        gCP[k + 2] -= 3.0 * G3k;
                        gCP[k + 3] += G3k;
                    }

                    // T-derivatives (scale paths) — keep exactly as you had:
                    gT[s] += dyn_constr_tilt_weight_ * (wj / (double)kappa) * phi; // dt-only
                    gT[s] += -5.0 * totalVel.dot(d1s) / (Ts * Ts);                 // v = 5/T d1
                    gT[s] += -40.0 * totalAcc.dot(d2s) / (Ts * Ts * Ts);           // a = 20/T^2 d2
                    gT[s] += -180.0 * totalJer.dot(d3s) / (Ts * Ts * Ts * Ts);     // j = 60/T^3 d3
                }
            }

            // --- Thrust ring via flatness map: ((thr - f_mean)^2 - f_radi^2) -------------
            // ---- Thrust ring: ((thr - f_mean)^2 - f_radi^2) (use C¹ hinge) ----
            {
                const double df = thr - f_mean;
                const double y = df * df - f_radi2;
                const double phi = smoothed_l1(y, mu);

                if (y > 0.0)
                {
                    const double dphi = smoothed_l1_prime(y, mu);
                    const double thr_grad = 2.0 * df * dphi * (dyn_constr_thrust_weight_ * wseg);

                    const Vec3 pos_grad = Vec3::Zero(), vel_grad = Vec3::Zero(), omg_grad = Vec3::Zero();
                    const Eigen::Vector4d quat_grad = Eigen::Vector4d::Zero();

                    Vec3 totalPos, totalVel, totalAcc, totalJer;
                    double totalPsi, totalPsiD;
                    flatmap_.backward(pos_grad, vel_grad, thr_grad, quat_grad, omg_grad,
                                      totalPos, totalVel, totalAcc, totalJer,
                                      totalPsi, totalPsiD);

                    const Vec3 g1 = totalVel * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP[k] -= G1k;
                        gCP[k + 1] += G1k;
                    }
                    const Vec3 g2_d2s = totalAcc / (Ts * Ts);
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * g2_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }
                    const Vec3 g3_d3s = totalJer / (Ts * Ts * Ts);
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3_d3s;
                        gCP[k] -= G3k;
                        gCP[k + 1] += 3.0 * G3k;
                        gCP[k + 2] -= 3.0 * G3k;
                        gCP[k + 3] += G3k;
                    }

                    gT[s] += dyn_constr_thrust_weight_ * (wj / (double)kappa) * phi;
                    gT[s] += -5.0 * totalVel.dot(d1s) / (Ts * Ts);
                    gT[s] += -40.0 * totalAcc.dot(d2s) / std::pow(Ts, 3);
                    gT[s] += -180.0 * totalJer.dot(d3s) / std::pow(Ts, 4);
                }
            }

        } // end samples

        // push gCP -> (P,V,A) for this segment and explicit CP(T) terms
        const double Ts2 = Ts * Ts;
        gP[s] += gCP[0] + gCP[1] + gCP[2];
        gV[s] += (Ts / 5.0) * gCP[1] + (2.0 * Ts / 5.0) * gCP[2];
        gA[s] += (Ts2 / 20.0) * gCP[2];

        gP[s + 1] += gCP[3] + gCP[4] + gCP[5];
        gV[s + 1] += (-2.0 * Ts / 5.0) * gCP[3] + (-Ts / 5.0) * gCP[4];
        gA[s + 1] += (Ts2 / 20.0) * gCP[3];

        // ∂J/∂T via CP(T) dependence
        gT[s] += gCP[1].dot(V[s] / 5.0);
        gT[s] += gCP[2].dot((2.0 / 5.0) * V[s] + (Ts / 10.0) * A[s]);
        gT[s] += gCP[3].dot((-2.0 / 5.0) * V[s + 1] + (Ts / 10.0) * A[s + 1]);
        gT[s] += gCP[4].dot((-1.0 / 5.0) * V[s + 1]);
    }

    // --- scatter (P,V,A) once ---
    scatterPosGrads(gP, z, grad);

    // 2) explicit interior-only scatter to v̂/â  (replaces accumVhatGrad/accumAhatGrad)
    std::vector<double> Tbar;
    build_Tbar(T, Tbar);

    // OLD (buggy): for (int i = 1; i <= M_; ++i) { accumVhatGrad(i, gV[i], grad); accumAhatGrad(i, gA[i], grad); }
    // NEW (correct):
    // v̂/â (interior only)
    for (int i = 1; i < M_; ++i)
    {
        const double Tb = std::max(1e-12, Tbar[i]);
        const int voff = vhatOffset(i);
        const int aoff = ahatOffset(i);
        grad.segment<3>(voff) += gV[i] / Tb;
        grad.segment<3>(aoff) += gA[i] / (Tb * Tb);
    }

    std::vector<double> gTbar(M_ + 1, 0.0);
    for (int i = 1; i < M_; ++i)
    {
        const double Tb = std::max(1e-12, Tbar[i]);
        const int voff = vhatOffset(i);
        const int aoff = ahatOffset(i);
        Eigen::Map<const Vec3> vhat_i(z.data() + voff);
        Eigen::Map<const Vec3> ahat_i(z.data() + aoff);
        gTbar[i] += -gV[i].dot(vhat_i) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_i) / (Tb * Tb * Tb);
    }

    // --- add T̄→{T_s} to segment ∂J/∂T and chain to τ ---
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
    {
        gT[s] += gTextra[s];
        grad[offTau_ + s] += gT[s] * dT_dtau(z[offTau_ + s]);
    }
}