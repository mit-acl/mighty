/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <mighty/lbfgs_solver.hpp>
#include <chrono>
#include <random>

using namespace lbfgs;

// H-polyhedron: rows are [n_x n_y n_z d] meaning n^T x <= d
using PolyhedronH = Eigen::Matrix<double, Eigen::Dynamic, 4>;
using PolyhedronV = Eigen::Matrix3Xd;
using PolyhedraH = std::vector<PolyhedronH>;
using PolyhedraV = std::vector<PolyhedronV>;

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

// Enumerate vertices of H-polyhedron by all 3-plane intersections.
inline bool enumerateVertices(const PolyhedronH &H, PolyhedronV &V, double tol = 1e-8)
{
    const int m = H.rows();
    if (m < 4)
        return false;

    const Eigen::MatrixXd A = H.leftCols<3>();
    const Eigen::VectorXd b = H.col(3);

    std::vector<Eigen::Vector3d> verts;
    verts.reserve(m * m);

    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j)
            for (int k = j + 1; k < m; ++k)
            {
                Eigen::Matrix3d A3;
                A3.row(0) = A.row(i);
                A3.row(1) = A.row(j);
                A3.row(2) = A.row(k);
                if (std::abs(A3.determinant()) < 1e-10)
                    continue;

                Eigen::Vector3d b3(b(i), b(j), b(k));
                Eigen::Vector3d x = A3.colPivHouseholderQr().solve(b3);

                // Feasibility: A x <= b (+tol)
                if (((A * x).array() <= (b.array() + tol)).all())
                {
                    bool unique = true;
                    for (auto &p : verts)
                        if ((p - x).norm() < 1e-7)
                        {
                            unique = false;
                            break;
                        }
                    if (unique)
                        verts.push_back(x);
                }
            }

    if (verts.size() < 4)
        return false; // likely unbounded/degenerate
    V.resize(3, (int)verts.size());
    for (int i = 0; i < (int)verts.size(); ++i)
        V.col(i) = verts[i];
    return true;
}

inline int zcp_offset_for_knot(int i, int M)
{
    return (i <= 0 || i >= M) ? -1 : 9 * (i - 1);
}

// === GCOPTER-style corridor "shortcut" ===
//
// A simple, effective pass that removes a corridor cell if its two neighbors
// still intersect after removing it: i.e., if  H[i-1] ∩ H[i+1]  is non-empty,
// then H[i] is unnecessary for connectivity and can be erased.
//
// This is the same spirit as GCOPTER's sfc_shortcut: aggressively remove
// overlapped/contained polytopes while preserving a chain of pairwise
// intersections.

inline bool hasNonEmptyIntersection(const PolyhedronH &Ha,
                                    const PolyhedronH &Hb,
                                    double tol = 1e-8)
{
    // Stack the halfspaces and test feasibility by enumerating vertices.
    PolyhedronH H(Ha.rows() + Hb.rows(), 4);
    H.topRows(Ha.rows()) = Ha;
    H.bottomRows(Hb.rows()) = Hb;

    PolyhedronV V;
    return enumerateVertices(H, V, tol); // true ⇒ bounded intersection with ≥4 verts
}

// Keep signature compatible with existing calls, but allow (optional) A_stat_/b_stat_ sync.
inline void sfc_shortcut(
    PolyhedraH &H,
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3>> *A_stat = nullptr,
    std::vector<Eigen::VectorXd> *b_stat = nullptr,
    int max_passes = 4,
    double tol = 1e-8)
{
    if (H.size() < 3)
        return;

    // We only synchronize A_stat_/b_stat_ if both are provided and their sizes match H.
    const bool sync_AB = (A_stat && b_stat &&
                          A_stat->size() == H.size() &&
                          b_stat->size() == H.size());

    bool changed = true;
    int passes = 0;

    while (changed && passes < max_passes)
    {
        changed = false;
        // Walk triples (i-1, i, i+1) and remove i if (i-1)∩(i+1) ≠ ∅
        for (size_t i = 1; i + 1 < H.size(); /* increment inside */)
        {
            if (hasNonEmptyIntersection(H[i - 1], H[i + 1], tol))
            {
                // erase the corridor cell
                H.erase(H.begin() + static_cast<long>(i));

                // erase the matching per-segment plane set (keep arrays aligned)
                if (sync_AB)
                {
                    A_stat->erase(A_stat->begin() + static_cast<long>(i));
                    b_stat->erase(b_stat->begin() + static_cast<long>(i));
                }

                changed = true;
                // don't ++i; re-check new triple at this position
            }
            else
            {
                ++i;
            }
        }
        ++passes;
    }

    // Defensive checks when syncing
    if (sync_AB)
    {
        // By construction, A_stat_/b_stat_ should track H size 1:1 (one cell ↔ one segment).
        if (A_stat->size() != H.size() || b_stat->size() != H.size())
        {
            std::cerr << "[sfc_shortcut] size mismatch after shortcut: "
                      << "H=" << H.size()
                      << " A_stat=" << A_stat->size()
                      << " b_stat=" << b_stat->size() << std::endl;
        }
    }
    else if (A_stat || b_stat)
    {
        // Warn if someone passed only one pointer or already mismatched sizes.
        std::cerr << "[sfc_shortcut] A_stat_/b_stat_ not synchronized: "
                  << "either nullptr or size did not match H at entry.\n";
    }
}

// Convert LinearConstraint3D → normalized H-polytopes (n_x n_y n_z d) per segment
inline void linearConstraintsToHPolys(
    const std::vector<LinearConstraint3D> &cons,
    PolyhedraH &hPolys)
{
    hPolys.clear();
    hPolys.reserve(cons.size());
    for (const auto &lc : cons)
    {
        // A x <= b ; rows are planes
        Eigen::Matrix<double, Eigen::Dynamic, 3> A =
            lc.A().template cast<double>();
        Eigen::VectorXd b = lc.b().template cast<double>();

        PolyhedronH H(A.rows(), 4);
        H.leftCols<3>() = A;
        H.col(3) = b;

        // Normalize rows (GCOPTER does this for better numerics)
        Eigen::ArrayXd nrm = H.leftCols<3>().rowwise().norm();
        for (int r = 0; r < H.rows(); ++r)
            if (nrm(r) > 1e-12)
                H.row(r) /= nrm(r);

        hPolys.emplace_back(std::move(H));
    }
}

// Convert a sequence of H-polytopes into V-representation lists:
// vPs = [poly0, inter01, poly1, inter12, ..., polyN]
inline bool processCorridor(const PolyhedraH &hPs, PolyhedraV &vPs)
{
    const int N = (int)hPs.size();
    if (N == 0)
        return false;

    vPs.clear();
    vPs.reserve(2 * N - 1);

    PolyhedronV V, VI;
    for (int i = 0; i < N - 1; ++i)
    {
        if (!enumerateVertices(hPs[i], V))
            return false;
        // OB form: col0 = origin, others = (vertex - origin)
        PolyhedronV O;
        O.resize(3, V.cols());
        O.col(0) = V.col(0);
        O.rightCols(V.cols() - 1) = V.rightCols(V.cols() - 1).colwise() - V.col(0);
        vPs.push_back(O);

        // Intersection of consecutive polytopes -> stack H rows then enumerate
        PolyhedronH HI(hPs[i].rows() + hPs[i + 1].rows(), 4);
        HI.topRows(hPs[i].rows()) = hPs[i];
        HI.bottomRows(hPs[i + 1].rows()) = hPs[i + 1];
        if (!enumerateVertices(HI, VI))
            return false;
        PolyhedronV OI;
        OI.resize(3, VI.cols());
        OI.col(0) = VI.col(0);
        OI.rightCols(VI.cols() - 1) = VI.rightCols(VI.cols() - 1).colwise() - VI.col(0);
        vPs.push_back(OI);
    }
    // last poly
    if (!enumerateVertices(hPs.back(), V))
        return false;
    PolyhedronV O;
    O.resize(3, V.cols());
    O.col(0) = V.col(0);
    O.rightCols(V.cols() - 1) = V.rightCols(V.cols() - 1).colwise() - V.col(0);
    vPs.push_back(O);
    return true;
}

// Shortest-path cost over intersection points (smoothed L2)
inline double costDistance(void *ptr, const Eigen::VectorXd &xi, Eigen::VectorXd &gradXi)
{
    void **data = (void **)ptr;
    const double &eps = *((const double *)(data[0]));
    const Eigen::Vector3d &p0 = *((const Eigen::Vector3d *)(data[1]));
    const Eigen::Vector3d &pf = *((const Eigen::Vector3d *)(data[2]));
    const PolyhedraV &vPolys = *((PolyhedraV *)(data[3]));

    const int overlaps = (int)vPolys.size() / 2; // # intersection polytopes
    gradXi.setZero(xi.size());

    Eigen::Matrix3Xd gradP = Eigen::Matrix3Xd::Zero(3, overlaps);
    Eigen::Vector3d a, b, d;
    double cost = 0.0;

    // forward pass: accumulate smooth distances and waypoint grads
    for (int i = 0, j = 0, k = 0; i <= overlaps; ++i, j += k)
    {
        a = (i == 0) ? p0 : b;
        if (i < overlaps)
        {
            const auto &OB = vPolys[2 * i + 1]; // intersection poly in OB form
            k = (int)OB.cols();
            Eigen::Map<const Eigen::VectorXd> q(xi.data() + j, k);

            // r are "square-root barycentric" weights; last dim is slack to enforce unit norm
            Eigen::VectorXd r = q.normalized().head(k - 1);
            b = OB.rightCols(k - 1) * r.cwiseProduct(r) + OB.col(0);
        }
        else
            b = pf;

        d = b - a;
        const double sm = std::sqrt(d.squaredNorm() + eps);
        cost += sm;
        if (i < overlaps)
            gradP.col(i) += d / sm;
        if (i > 0)
            gradP.col(i - 1) -= d / sm;
    }

    // backward pass: dcost/dxi via chain rule on r(q)
    for (int i = 0, j = 0, k; i < overlaps; ++i, j += k)
    {
        const auto &OB = vPolys[2 * i + 1];
        k = (int)OB.cols();
        Eigen::Map<const Eigen::VectorXd> q(xi.data() + j, k);
        Eigen::Map<Eigen::VectorXd> g(gradXi.data() + j, k);

        const double sq = q.squaredNorm();
        const double inv = 1.0 / std::sqrt(sq);
        const Eigen::VectorXd uq = q * inv;

        g.head(k - 1) = (OB.rightCols(k - 1).transpose() * gradP.col(i)).array() * uq.head(k - 1).array() * 2.0;
        g(k - 1) = 0.0;                 // slack dim
        g = (g - uq * uq.dot(g)) * inv; // project to sphere tangent space

        // soft penalty to keep ||q|| >= 1 (same trick as GCOPTER)
        const double viol = sq - 1.0;
        if (viol > 0.0)
        {
            const double c = viol * viol * viol;
            const double dc = 3.0 * viol * viol;
            g += dc * 2.0 * q;
            cost += c;
        }
    }

    return cost;
}

// Solve for shortest path waypoints: returns [p0, w1, w2, ..., pf]
inline void getShortestPath(const Eigen::Vector3d &p0,
                            const Eigen::Vector3d &pf,
                            const PolyhedraV &vPolys,
                            double smooth_eps,
                            Eigen::Matrix3Xd &path)
{
    const int overlaps = (int)vPolys.size() / 2;
    // Decision: for each intersection poly, one vector q \in R^{k} (k = #verts of poly in OB form)
    Eigen::VectorXi sizes(overlaps);
    for (int i = 0; i < overlaps; ++i)
        sizes(i) = vPolys[2 * i + 1].cols();

    Eigen::VectorXd xi(sizes.sum());
    // init q blocks with uniform weights on sphere
    for (int i = 0, offset = 0; i < overlaps; ++i)
    {
        const int k = sizes(i);
        xi.segment(offset, k).setConstant(std::sqrt(1.0 / k));
        offset += k;
    }

    void *ptrs[4];
    ptrs[0] = (void *)(&smooth_eps);
    ptrs[1] = (void *)(&p0);
    ptrs[2] = (void *)(&pf);
    ptrs[3] = (void *)(&vPolys);

    double minDist = 0.0;
    lbfgs::lbfgs_parameter_t sp_params;
    sp_params.past = 3;
    sp_params.delta = 1.0e-3;
    sp_params.g_epsilon = 1.0e-5;

    lbfgs::lbfgs_optimize(xi, minDist, &costDistance, nullptr, nullptr, ptrs, sp_params);

    // decode xi -> points
    path.resize(3, overlaps + 2);
    path.col(0) = p0;
    path.col(overlaps + 1) = pf;

    for (int i = 0, off = 0; i < overlaps; ++i)
    {
        const auto &OB = vPolys[2 * i + 1];
        const int k = (int)OB.cols();
        Eigen::Map<const Eigen::VectorXd> q(xi.data() + off, k);
        Eigen::VectorXd r = q.normalized().head(k - 1);
        path.col(i + 1) = OB.rightCols(k - 1) * r.cwiseProduct(r) + OB.col(0);
        off += k;
    }
}

// ====== Helpers ======

// Compute geometric median of a 3xK vertex matrix (Weiszfeld).
inline Eigen::Vector3d geometricMedian(const Eigen::Matrix3Xd &V)
{
    const int n = (int)V.cols();
    Eigen::Vector3d x = V.rowwise().mean();
    const int kMax = 50;
    const double eps = 1e-8;
    for (int k = 0; k < kMax; ++k)
    {
        Eigen::Vector3d num = Eigen::Vector3d::Zero();
        double den = 0.0;
        bool hit = false;
        for (int i = 0; i < n; ++i)
        {
            Eigen::Vector3d d = x - V.col(i);
            double w = d.norm();
            if (w < 1e-12)
            {
                x = V.col(i);
                hit = true;
                break;
            }
            w = 1.0 / w;
            num += w * V.col(i);
            den += w;
        }
        if (hit)
            break;
        Eigen::Vector3d xnew = num / std::max(den, 1e-16);
        if ((xnew - x).norm() < eps)
        {
            x = xnew;
            break;
        }
        x = xnew;
    }
    return x;
}

// Build actual vertex list from OB form (col0=V0, cols>=1 are offsets).
inline Eigen::Matrix3Xd OB_to_vertices(const Eigen::Matrix3Xd &OB)
{
    const int k = (int)OB.cols();
    Eigen::Matrix3Xd V(3, k);
    V.col(0) = OB.col(0);
    if (k > 1)
        V.rightCols(k - 1) = OB.rightCols(k - 1).colwise() + OB.col(0);
    return V;
}

// Compute a robust interior point for a *cell poly* in V/OB form (even indices of vPolys).
inline Eigen::Vector3d cellInteriorPointFromV(const Eigen::Matrix3Xd &cellOB)
{
    Eigen::Matrix3Xd V = OB_to_vertices(cellOB);
    return geometricMedian(V); // centroid(V) also OK but less robust
}

// Given H=[A|b] with ||row_i(A)||=1, compute the per-cell safe margin from point c:
// m_safe = tau * min_i (b_i - a_i^T c). Clamp to >=0.
inline double safeMarginFromPoint(const Eigen::MatrixXd &H, const Eigen::Vector3d &c, double tau)
{
    const Eigen::MatrixXd A = H.leftCols<3>();
    const Eigen::VectorXd b = H.col(3);
    Eigen::VectorXd slack = b - (A * c);
    double min_slack = slack.minCoeff();
    return std::max(0.0, tau * min_slack);
}

// Apply per-cell uniform inset: b <- b - m_i
template <typename Derived>
inline void insetCellH(Eigen::MatrixBase<Derived> &H, double m)
{
    if (m <= 0.0)
        return;
    H.derived().col(3).array() -= m; // b <- b - m  (assuming H is [A|b] and rows normalized)
}

inline void accumulate_cp_to_pvaT(
    int s,
    const std::array<Vec3, 6> &gCPs,
    const std::vector<Vec3> &V,
    const std::vector<Vec3> &A,
    const std::vector<double> &T,
    std::vector<Vec3> &gP,
    std::vector<Vec3> &gV,
    std::vector<Vec3> &gA,
    std::vector<double> &gT_cp)
{
    const double Ts = T[s];
    const double Ts2 = Ts * Ts;

    // left endpoint (s)
    gP[s] += gCPs[0] + gCPs[1] + gCPs[2];
    gV[s] += (Ts / 5.0) * gCPs[1] + (2.0 * Ts / 5.0) * gCPs[2];
    gA[s] += (Ts2 / 20.0) * gCPs[2];

    // right endpoint (s+1)
    gP[s + 1] += gCPs[3] + gCPs[4] + gCPs[5];
    gV[s + 1] += (-2.0 * Ts / 5.0) * gCPs[3] + (-Ts / 5.0) * gCPs[4];
    gA[s + 1] += (Ts2 / 20.0) * gCPs[3];

    // ∂/∂T via CP(T) dependence
    gT_cp[s] += gCPs[1].dot(V[s] / 5.0);
    gT_cp[s] += gCPs[2].dot((2.0 / 5.0) * V[s] + (Ts / 10.0) * A[s]);
    gT_cp[s] += gCPs[3].dot((-2.0 / 5.0) * V[s + 1] + (Ts / 10.0) * A[s + 1]);
    gT_cp[s] += gCPs[4].dot((-1.0 / 5.0) * V[s + 1]);
}

// ---- end: helpers ported from mighty_solver.cpp ----

// -----------------------------------------------------------------------------

SolverLBFGS::SolverLBFGS()
{
    // Constructor implementation
}

// -----------------------------------------------------------------------------

SolverLBFGS::~SolverLBFGS()
{
    // Destructor implementation
    // Clean up resources if necessary
}

// -----------------------------------------------------------------------------

void SolverLBFGS::initializeSolver(const planner_params_t &params)
{
    // Initialize the solver with parameters that won't change throughout the mission
    verbose_ = params.verbose;                               // Verbosity level
    V_max_ = params.V_max;                                   // Max velocity
    A_max_ = params.A_max;                                   // Max acceleration
    J_max_ = params.J_max;                                   // Max jerk
    num_perturbation_ = params.num_perturbation;             // Number of perturbations for initial guesses
    r_max_ = params.r_max;                                   // Perturbation radius for initial guesses
    time_weight_ = params.time_weight;                       // Weight for time in the objective
    pos_anchor_weight_ = params.pos_anchor_weight;           // Weight for position anchors in the objective
    dyn_weight_ = params.dyn_weight;                         // Weight for dynamic avoidance in the objective
    stat_weight_ = params.stat_weight;                       // Weight for static avoidance in the objective
    jerk_weight_ = params.jerk_weight;                       // Weight for jerk in the objective
    dyn_constr_vel_weight_ = params.dyn_constr_vel_weight;   // Weight for dynamic velocity constraints
    dyn_constr_acc_weight_ = params.dyn_constr_acc_weight;   // Weight for dynamic acceleration constraints
    dyn_constr_jerk_weight_ = params.dyn_constr_jerk_weight; // Weight for dynamic jerk constraints
    dyn_constr_bodyrate_weight_ = params.dyn_constr_bodyrate_weight;
    dyn_constr_tilt_weight_ = params.dyn_constr_tilt_weight;
    dyn_constr_thrust_weight_ = params.dyn_constr_thrust_weight;
    num_dyn_obst_samples_ = params.num_dyn_obst_samples;         // Number of dynamic obstacle samples
    Co_ = params.Co;                                             // Clearance distance for static obstacle avoidance
    Cw_ = params.Cw;                                             // Clearance distance for dynamic obstacle avoidance
    BIG_ = params.BIG;                                           // A large constant for static constraints
    dc_ = params.dc;                                             // Discretization constant
    V_min_ = 0.0;                                                // Minimum speed
    turn_buf_ = params.init_turn_bf * M_PI / 180;                // 15° in radians
    turn_span_ = M_PI - turn_buf_;                               // over which we ramp down
    cos_thresh_ = std::cos(turn_buf_ * M_PI / 180.0);
    Cw2_ = Cw_ * Cw_;
    integral_resolution_ = params.integral_resolution; // e.g., 30
    hinge_mu_ = params.hinge_mu;                       // e.g., 1e-2
    omege_max_ = params.omega_max;                     // e.g., 6.0
    tilt_max_rad_ = params.tilt_max_rad;               // e.g., 35° in rad
    f_min_ = params.f_min;
    f_max_ = params.f_max;
    mass_ = params.mass;
    g_ = params.g;
}

// -----------------------------------------------------------------------------

// Numerically stable central-difference directional derivative
double SolverLBFGS::centralDiff(const VecXd &z,
                                const VecXd &d,
                                double eps_base)
{
    // Scale step per IEEE advice
    double eps = eps_base / std::max(1.0, d.norm());
    VecXd zp = z + eps * d;
    VecXd zm = z - eps * d;
    // return (evaluateObjective(zp) - evaluateObjective(zm)) / (2.0 * eps);
    double fp, fm;
    VecXd g(z.size());
    fp = evaluateObjectiveAndGradientFused(zp, g);
    fm = evaluateObjectiveAndGradientFused(zm, g);
    return (fp - fm) / (2.0 * eps);
}

void SolverLBFGS::checkGradDirectional(const VecXd &z0, int num_dirs, double eps, unsigned seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> N(0.0, 1.0);

    VecXd g(z0.size());
    // (void)evaluateObjectiveAndGradient(z0, g);
    (void)evaluateObjectiveAndGradientFused(z0, g);

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

void SolverLBFGS::checkGradCoordinates(const VecXd &z0, int max_coords, double eps, unsigned seed)
{
    VecXd g(z0.size());
    // (void)evaluateObjectiveAndGradient(z0, g);
    (void)evaluateObjectiveAndGradientFused(z0, g);

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

// Build corridor-shortest polyline and overwrite global_wps_.
// Optionally apply sfc_shortcut before processCorridor.
bool SolverLBFGS::replaceGlobalPathWithCorridorShortest(bool do_shortcut,
                                                        double smooth_eps)
{
    // Ensure we have per-segment planes (A_stat_, b_stat_) already set
    if (A_stat_.empty() || A_stat_.size() != b_stat_.size())
        return false;

    // 1) H-polys from (A_stat_, b_stat_) with GCOPTER row normalization
    PolyhedraH hPolys;
    hPolys.reserve(A_stat_.size());
    for (size_t i = 0; i < A_stat_.size(); ++i)
    {
        const auto &A = A_stat_[i];
        const auto &b = b_stat_[i];
        if (A.rows() == 0 || b.size() == 0)
            return false;

        PolyhedronH H(A.rows(), 4);
        H.leftCols<3>() = A;
        H.col(3) = b;

        Eigen::ArrayXd nrm = H.leftCols<3>().rowwise().norm();
        for (int r = 0; r < H.rows(); ++r)
            if (nrm(r) > 1e-12)
                H.row(r) /= nrm(r);

        hPolys.emplace_back(std::move(H));
    }

    // 2) (Optional) GCOPTER-style shortcut of the corridor
    if (do_shortcut)
        sfc_shortcut(hPolys, &A_stat_, &b_stat_);

    // 3) H→V with intersections: [poly0, inter01, poly1, ...]
    PolyhedraV vPolys;
    if (!processCorridor(hPolys, vPolys)) // uses your existing helper
        return false;

    // 3.2) Get the center points in the overlapping polytopes (odd indices of vPolys)

    const int num_odd_cells = (int)vPolys.size() / 2; // odd indices are intersections
    std::vector<Eigen::Vector3d> centers;
    centers.reserve(num_odd_cells);

    for (int ci = 0; ci < num_odd_cells; ++ci)
    {
        const int vIdx = 2 * ci + 1; // odd indices are intersections
        const auto &cellOB = vPolys[vIdx];
        centers.push_back(cellInteriorPointFromV(cellOB));
    }

    P_anchor_.clear();
    P_anchor_.reserve(centers.size() + 2);
    P_anchor_.push_back(x0_);
    for (const auto &c : centers)
        P_anchor_.push_back(c);
    P_anchor_.push_back(xf_);

    // 4) Solve for the shortest polyline p0→...→pf through intersections
    Eigen::Matrix3Xd path;
    getShortestPath(x0_, xf_, vPolys, smooth_eps, path); // existing helper

    // We'll build global_wps_ from 'path' while skipping close seams.
    global_wps_.clear();
    global_wps_.reserve(path.cols());

    auto colAsVec3 = [&](int c)
    { return Vec3(path(0, c), path(1, c), path(2, c)); };

    // push start
    global_wps_.push_back(colAsVec3(0));

    // i runs over interior seams [1 .. path.cols()-2]
    // If seam i is removed, we also erase corridor cell i.
    for (int i = 1; i + 1 < path.cols(); /* increment inside */)
    {
        Vec3 a = colAsVec3(i);
        // Keep seam 'i'
        global_wps_.push_back(a);
        ++i;
    }

    // push final
    global_wps_.push_back(colAsVec3(path.cols() - 1));

    // --- after pruning, recompute sizes that depend on #segments ----------
    M_ = static_cast<int>(global_wps_.size()) - 1;

    // One corridor plane set per segment:
    if (M_ != static_cast<int>(A_stat_.size()) || M_ != static_cast<int>(b_stat_.size()))
    {
        std::cerr << "[sfc] size mismatch after seam pruning: "
                  << "M_=" << M_ << " A_stat=" << A_stat_.size()
                  << " b_stat=" << b_stat_.size() << std::endl;
        return false;
    }

    // Use the interior‑knot layout used elsewhere in this file:
    K_cp_ = 9 * std::max(0, M_ - 1);
    K_sig_ = M_;
    K_ = K_cp_ + K_sig_;

    return true;
}

// -----------------------------------------------------------------------------

void SolverLBFGS::prepareSolverForReplan(double t0,
                                         const vec_Vec3f &global_wps,
                                         const std::vector<LinearConstraint3D> &safe_corridor,
                                         const std::vector<std::shared_ptr<dynTraj>> &obstacles,
                                         const state &initial_state,
                                         const state &goal_state,
                                         double &initial_guess_computation_time,
                                         bool use_for_safe_path,
                                         bool use_multiple_initial_guesses)
{
    // Initialize the solver with parameters that changes for replanning
    t0_ = t0; // Current time in the trajectory

    // Global waypoints
    global_wps_.clear();
    for (const auto &wp : global_wps)
    {
        global_wps_.emplace_back(wp.x(), wp.y(), wp.z());
    }

    // Number of segments
    M_ = global_wps_.size() - 1; // Number of segments is one less than the number of waypoints

    // free control points indices
    K_cp_ = 9 * (M_ + 1); // Number of decision variables (CPs + slack times)
    K_sig_ = M_;          // Number of slack times (K_sig_ = M_)
    K_ = K_cp_ + K_sig_;  // Total number of decision variables

    // Set initial and goal states
    x0_ = initial_state.pos;
    v0_ = initial_state.vel;
    a0_ = initial_state.accel;

    // Copy the dynamic obstacles
    obstacles_.clear();
    obstacles_ = obstacles;

    // Create A_stat_ and b_stat_ from safe corridor
    if (use_for_safe_path)
    {
        setStaticConstraintsForSafePath(safe_corridor);
    }
    else
    {
        setStaticConstraints(safe_corridor);
    }

    // Push the waypoints by static corridor
    // pushWaypointsByStaticCorridor(global_wps_);

    xf_ = global_wps_.back();
    vf_ = goal_state.vel;
    af_ = goal_state.accel;

    // ---- REPLACE HEURISTIC ROUTE WITH SHORTEST-PATH THROUGH THE CORRIDOR ----
    if (!use_for_safe_path)
    {

        if (!replaceGlobalPathWithCorridorShortest(/*do_shortcut=*/false, /*smooth_eps=*/1e-6))
        {
            if (verbose_)
                std::cout << "[init] Corridor-shortest init failed; "
                             "falling back to original route.\n";
        }
    }
    // ------------------------------------------------------------------------------

    // Single initial guess: use min-jerk helper
    initial_guess_wps_.clear();
    initial_guess_wps_.push_back(global_wps_);

    list_z0_.clear();

    // Allocate containers
    Eigen::VectorXd z0;
    std::vector<Vec3> P, V, A;
    std::vector<double> T;

    // Time can be computed by deviding the distance by the max velocity
    for (int i = 0; i < M_; ++i)
    {
        double dist = (global_wps_[i + 1] - global_wps_[i]).norm();
        double t = dist / V_max_;
        T.push_back(t);
    }

    for (int i = 0; i < M_ + 1; ++i)
    {
        P.push_back(global_wps_[i]);
    }

    // // Save this P as anchors
    // if (pos_anchor_weight_ > 0.0) {
    //     P_anchor_ = P;
    // }

    // Initialize velocities and accelerations
    auto start = std::chrono::high_resolution_clock::now();
    findInitialGuess(T, V, A);
    auto end = std::chrono::high_resolution_clock::now();
    // Measure time in microseconds and convert to milliseconds
    initial_guess_computation_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

    // Generate the initial guess in z form
    packDecisionVariables(P, V, A, T, z0);

    // Store the results
    list_z0_.push_back(z0);

    // Sanity check
    sanityCheck();

    // Check cost objective and gradient
    // VecXd g_test1(z0.size());
    // double f1 = evaluateObjectiveAndGradientFused(z0, g_test1);
    // VecXd g_test2(z0.size());
    // double f2 = evaluateObjectiveAndGradient(z0, g_test2);

    // // Check the difference
    // double err_g = (g_test1 - g_test2).norm();
    // double err_f = std::abs(f1 - f2);

    // if (err_g > 1e-5 || err_f > 1e-5)
    // {
    //     // print in red
    //     std::cout << "\033[1;31m |f1-f2|=" << err_f << " ||g1-g2||=" << err_g << "\033[0m\n";
    // }

    // wherever you have a valid z (e.g., right before starting L-BFGS)
    // checkGradDirectional(z0, /*num_dirs=*/8, /*eps=*/1e-5, /*seed=*/42);
    // checkGradCoordinates(z0, /*max_coords=*/256, /*eps=*/1e-5, /*seed=*/43);
}

// -----------------------------------------------------------------------------

void SolverLBFGS::getGlobalPath(vec_Vecf<3> &global_path)
{
    global_path.clear();
    global_path.reserve(global_wps_.size());
    for (auto const &wp : global_wps_)
    {
        // wp is Eigen::Vector3d; Vec3f is Eigen::Matrix<double,3,1>
        global_path.emplace_back(wp.x(), wp.y(), wp.z());
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::packDecisionVariables(
    const std::vector<Vec3> &P,
    const std::vector<Vec3> &V,
    const std::vector<Vec3> &A,
    const std::vector<double> &T,
    VecXd &z) const
{
    // CP block (interior only) + time block (all segments)
    z.resize(K_);

    // per-knot T̄ for v̂, â scaling
    std::vector<double> Tbar(M_ + 1);
    build_Tbar(T, Tbar);

    // ---- pack interior knots only: i = 1..M_-1 ----
    for (int i = 1; i < M_; ++i)
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = Tbar[i];
        const double Tb2 = Tb * Tb;

        // positions
        z(base + 0) = P[i].x();
        z(base + 1) = P[i].y();
        z(base + 2) = P[i].z();

        // v̂ = T̄ * V
        const Vec3 vhat = Tb * V[i];
        z(base + 3) = vhat.x();
        z(base + 4) = vhat.y();
        z(base + 5) = vhat.z();

        // â = T̄^2 * A
        const Vec3 ahat = Tb2 * A[i];
        z(base + 6) = ahat.x();
        z(base + 7) = ahat.y();
        z(base + 8) = ahat.z();
    }

    // ---- pack τ for time (all segments) ----
    for (int s = 0; s < M_; ++s)
        z[K_cp_ + s] = T_to_tau(T[s]);
}

// -----------------------------------------------------------------------------

void SolverLBFGS::reconstruct(
    const VecXd &z,
    std::vector<Vec3> &P,
    std::vector<Vec3> &V,
    std::vector<Vec3> &A,
    std::vector<std::array<Vec3, 6>> &CP,
    std::vector<double> &T) const
{
    const int knots = M_ + 1;
    P.resize(knots);
    V.resize(knots);
    A.resize(knots);
    T.resize(M_);
    CP.resize(M_);

    // 1) τ → T for all segments
    for (int s = 0; s < M_; ++s)
        T[s] = tau_to_T(z[K_cp_ + s]);

    // 2) T̄ for v̂/â decode of interior knots
    std::vector<double> Tbar(knots);
    build_Tbar(T, Tbar);

    // 3) set endpoints from boundary conditions
    P[0] = x0_;
    V[0] = v0_;
    A[0] = a0_;

    P[M_] = xf_;
    V[M_] = vf_;
    A[M_] = af_;

    // 4) decode interior knots from z
    for (int i = 1; i < M_; ++i)
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = std::max(1e-12, Tbar[i]);

        // position
        P[i] = Vec3(z[base + 0], z[base + 1], z[base + 2]);

        // v̂ → V, â → A
        const Vec3 vhat(z[base + 3], z[base + 4], z[base + 5]);
        const Vec3 ahat(z[base + 6], z[base + 7], z[base + 8]);
        V[i] = vhat / Tb;
        A[i] = ahat / (Tb * Tb);
    }

    // 5) control points (unchanged)
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
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Build H_acc (N×N) and b_acc (N×3) for the "accelerations only" min-jerk problem.
// Unknown vector is [A1, A2, ..., A_{M-1}], each Ai is 3D.
// Uses your K_r(T) and k_r(T,P0,P1).
// -----------------------------------------------------------------------------
void SolverLBFGS::assemble_H_b_acc_only(
    const std::vector<Vec3> &wps,   // size M+1
    const std::vector<double> &T,   // size M
    const std::vector<Vec3> &V,     // size M+1 (ALL velocities fixed)
    const Vec3 &a0, const Vec3 &af, // fixed endpoint accelerations
    Eigen::MatrixXd &H,             // out: (N×N)
    Eigen::MatrixXd &b              // out: (N×3)
)
{
    const int M = static_cast<int>(T.size());
    const int N = M - 1;
    if (N <= 0)
    {
        H.resize(0, 0);
        b.resize(0, 3);
        return;
    }

    H = Eigen::MatrixXd::Zero(N, N);
    b = Eigen::MatrixXd::Zero(N, 3);

    for (int r = 0; r < M; ++r)
    {
        const Eigen::Matrix4d Kr = K_r(T[r]); // 4×4 over [V_r, A_r, V_{r+1}, A_{r+1}]

        // k_r as a (4×3) block, column-wise per dimension
        Eigen::Matrix<double, 4, 3> kr;
        for (int dim = 0; dim < 3; ++dim)
        {
            const auto c = k_r(T[r], wps[r](dim), wps[r + 1](dim)); // {k0,k1,k2,k3}
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
            b.row(j) += (Kr(3, 0) * V[r].transpose() + Kr(3, 2) * V[r + 1].transpose() + Kr(3, 1) * a0.transpose() + kr.row(3));
        }
        else if (r == M - 1)
        {
            // Only A_{M-1} is unknown → row uses A_r (index 1 in Kr)
            const int i = N - 1; // maps to A_{M-1}
            H(i, i) += Kr(1, 1);
            // b_i += Kr[1,0]*V_r + Kr[1,2]*V_{r+1} + Kr[1,3]*af + k_r[1]
            b.row(i) += (Kr(1, 0) * V[r].transpose() + Kr(1, 2) * V[r + 1].transpose() + Kr(1, 3) * af.transpose() + kr.row(1));
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
            b.row(i) += (Kr(1, 0) * V[r].transpose() + Kr(1, 2) * V[r + 1].transpose() + kr.row(1));
            // b_j += Kr[3,0]*V_r + Kr[3,2]*V_{r+1} + k_r[3]
            b.row(j) += (Kr(3, 0) * V[r].transpose() + Kr(3, 2) * V[r + 1].transpose() + kr.row(3));
        }
    }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Solve H_acc * a = -b_acc for interior accelerations,
// then pack full A with A[0]=a0 and A[M]=af.
// -----------------------------------------------------------------------------
void SolverLBFGS::solveMinJerkAccOnlyClosedForm(
    const std::vector<Vec3> &wps,   // size M+1
    const std::vector<double> &T,   // size M
    const std::vector<Vec3> &V,     // size M+1 (fixed)
    const Vec3 &a0, const Vec3 &af, // fixed endpoint accelerations
    std::vector<Vec3> &A_out        // out: size M+1
)
{
    const int M = static_cast<int>(T.size());
    const int N = M - 1;

    A_out.clear();
    A_out.resize(M + 1);
    A_out.front() = a0;
    A_out.back() = af;

    if (N <= 0)
        return; // no interior unknowns

    Eigen::MatrixXd H, b; // (N×N), (N×3)
    assemble_H_b_acc_only(wps, T, V, a0, af, H, b);

    // Solve H * X = -b, where X is (N×3) stacked [A1; ...; A_{M-1}]
    Eigen::MatrixXd X = H.ldlt().solve(-b);

    for (int i = 1; i <= M - 1; ++i)
    {
        A_out[i] = X.row(i - 1).transpose(); // Vec3
    }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

void SolverLBFGS::findInitialGuess(
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

    // V[0] = v0_;
    // for (int i = 1; i < M_; ++i)
    // {
    //     V[i] = dirs[i] * s[i];
    // }
    // V[M_] = vf_;

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

    // --- 6) generate initial V, A via your min-jerk helper ---
    // solveMinJerkVelAcc(global_wps_, T, v0_, a0_, vf_, af_, V, A);
}

// -----------------------------------------------------------------------------

//  Closed-form quintic control‐points for one segment
std::array<Vec3, 6> SolverLBFGS::computeQuinticCP(
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

inline Eigen::Matrix4d SolverLBFGS::K_r(double T)
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

inline std::array<double, 4> SolverLBFGS::k_r(double T, double P0, double P1)
{
    // Example placeholder — you must substitute your real formulas:
    double k0 = (-66960 * P0 * T + 66960 * P1 * T) / (2 * T * T * T * T * T) + (-24480 * P0 * T + 24480 * P1 * T) / (2 * T * T * T * T * T) + (4320 * P0 * T - 4320 * P1 * T) / (2 * T * T * T * T * T) + (25920 * P0 * T - 25920 * P1 * T) / (2 * T * T * T * T * T) + (61920 * P0 * T - 61920 * P1 * T) / (2 * T * T * T * T * T);

    double k1 = (-11880 * P0 * T * T + 11880 * P1 * T * T) / (2 * T * T * T * T * T) + (-5400 * P0 * T * T + 5400 * P1 * T * T) / (2 * T * T * T * T * T) + (1080 * P0 * T * T - 1080 * P1 * T * T) / (2 * T * T * T * T * T) + (4320 * P0 * T * T - 4320 * P1 * T * T) / (2 * T * T * T * T * T) + (12000 * P0 * T * T - 12000 * P1 * T * T) / (2 * T * T * T * T * T);

    double k2 = (-62640 * P0 * T + 62640 * P1 * T) / (2 * T * T * T * T * T) + (-18720 * P0 * T + 18720 * P1 * T) / (2 * T * T * T * T * T) + (2880 * P0 * T - 2880 * P1 * T) / (2 * T * T * T * T * T) + (25920 * P0 * T - 25920 * P1 * T) / (2 * T * T * T * T * T) + (53280 * P0 * T - 53280 * P1 * T) / (2 * T * T * T * T * T);

    double k3 = (-7680 * P0 * T * T + 7680 * P1 * T * T) / (2 * T * T * T * T * T) + (-4320 * P0 * T * T + 4320 * P1 * T * T) / (2 * T * T * T * T * T) + (-360 * P0 * T * T + 360 * P1 * T * T) / (2 * T * T * T * T * T) + (2520 * P0 * T * T - 2520 * P1 * T * T) / (2 * T * T * T * T * T) + (9720 * P0 * T * T - 9720 * P1 * T * T) / (2 * T * T * T * T * T);
    return {k0, k1, k2, k3};
}

// -----------------------------------------------------------------------------

void SolverLBFGS::assemble_H_b(
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
void SolverLBFGS::solveMinJerkVelAcc(
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

void SolverLBFGS::getGoalSetpoints(std::vector<state> &goal_setpoints)
{
    // 1) total time & sample count
    double total_time = std::accumulate(T_opt_.begin(), T_opt_.end(), 0.0);
    // BUG #1: truncating here can undershoot the end of the trajectory.
    // It’s usually better to do:
    //   int N = static_cast<int>(std::ceil(total_time / dc_));
    int N = static_cast<int>(total_time / dc_);
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

void SolverLBFGS::getPieceWisePol(PieceWiseQuinticPol &pwp)
{
    // 1) reset
    pwp.clear();

    // 2) build the time‐breaks
    pwp.times.push_back(t0_);
    double t_acc = 0.0;
    for (int i = 0; i < M_; ++i)
    {
        t_acc += T_opt_[i];
        pwp.times.push_back(t0_ + t_acc);
    }

    // 3) for each segment, convert Hermite → power basis:
    for (int s = 0; s < M_; ++s)
    {
        double T = T_opt_[s];
        double T2 = T * T, T3 = T2 * T, T4 = T3 * T, T5 = T4 * T;

        // Hermite endpoints
        const auto &P0 = P_opt_[s], &P1 = P_opt_[s + 1];
        const auto &V0 = V_opt_[s], &V1 = V_opt_[s + 1];
        const auto &A0 = A_opt_[s], &A1 = A_opt_[s + 1];

        // helper deltas
        Eigen::Vector3d C0 = P1 - (P0 + V0 * T + 0.5 * A0 * T2);
        Eigen::Vector3d C1 = V1 - (V0 + A0 * T);
        Eigen::Vector3d C2 = A1 - A0;

        // build coeff vector [a0, a1, a2, a3, a4, a5]^T per axis
        Eigen::Matrix<double, 6, 1> cx, cy, cz;

        // a0, a1, a2
        cx(0) = P0.x();
        cx(1) = V0.x();
        cx(2) = 0.5 * A0.x();
        // a3, a4, a5
        cx(3) = (10 * C0.x() - 4 * C1.x() + 0.5 * C2.x()) / T3;
        cx(4) = (-15 * C0.x() + 7 * C1.x() - C2.x()) / T4;
        cx(5) = (6 * C0.x() - 3 * C1.x() + 0.5 * C2.x()) / T5;

        // same for y
        cy(0) = P0.y();
        cy(1) = V0.y();
        cy(2) = 0.5 * A0.y();
        cy(3) = (10 * C0.y() - 4 * C1.y() + 0.5 * C2.y()) / T3;
        cy(4) = (-15 * C0.y() + 7 * C1.y() - C2.y()) / T4;
        cy(5) = (6 * C0.y() - 3 * C1.y() + 0.5 * C2.y()) / T5;

        // …and z
        cz(0) = P0.z();
        cz(1) = V0.z();
        cz(2) = 0.5 * A0.z();
        cz(3) = (10 * C0.z() - 4 * C1.z() + 0.5 * C2.z()) / T3;
        cz(4) = (-15 * C0.z() + 7 * C1.z() - C2.z()) / T4;
        cz(5) = (6 * C0.z() - 3 * C1.z() + 0.5 * C2.z()) / T5;

        // store it
        pwp.coeff_x.push_back(cx);
        pwp.coeff_y.push_back(cy);
        pwp.coeff_z.push_back(cz);
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::getControlPoints(std::vector<Eigen::Matrix<double, 3, 6>> &cps)
{
    // Resize the control points vector
    cps.resize(M_);

    // Fill the control points
    for (int i = 0; i < M_; i++)
    {
        Eigen::Matrix<double, 3, 6> cp_i;
        for (int j = 0; j < 6; j++)
        {
            cp_i(0, j) = CP_opt_[i][j][0];
            cp_i(1, j) = CP_opt_[i][j][1];
            cp_i(2, j) = CP_opt_[i][j][2];
        }
        cps[i] = cp_i;
    }
}

// -----------------------------------------------------------------------------

inline StateDeriv SolverLBFGS::evalStateDeriv(int s, double tau) const
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

// -----------------------------------------------------------------------------

// degree-5 Bézier: p(u) = Σ_{j=0..5} B5_j(u) * CP[j]
inline Eigen::Vector3d evalBezierQuintic(const std::array<Eigen::Vector3d, 6> &cp, double u)
{
    u = std::min(std::max(u, 0.0), 1.0);
    const double omu = 1.0 - u;
    const double omu2 = omu * omu;
    const double omu3 = omu2 * omu;
    const double omu4 = omu3 * omu;
    const double omu5 = omu4 * omu;
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double u4 = u3 * u;
    const double u5 = u4 * u;

    const double b0 = omu5;
    const double b1 = 5.0 * u * omu4;
    const double b2 = 10.0 * u2 * omu3;
    const double b3 = 10.0 * u3 * omu2;
    const double b4 = 5.0 * u4 * omu;
    const double b5 = u5;

    return b0 * cp[0] + b1 * cp[1] + b2 * cp[2] + b3 * cp[3] + b4 * cp[4] + b5 * cp[5];
}

// Sample robot positions at absolute times t_i = t0 + i*dt, i=1..N
// Uses a single pass over segments (monotone times) for speed.
inline void sampleRobotPositionsUniform(
    const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
    const std::vector<double> &T,
    double t0, int N,
    std::vector<double> &t_samples,
    std::vector<Eigen::Vector3d> &p_samples)
{
    const int M = static_cast<int>(T.size());
    t_samples.resize(N);
    p_samples.resize(N);

    // absolute edges: [t0, t0+T0, t0+T0+T1, ...]
    std::vector<double> edges(M + 1, t0);
    for (int i = 1; i <= M; ++i)
        edges[i] = edges[i - 1] + T[i - 1];

    const double total_T = edges.back() - t0;
    const double dt = (N > 0) ? (total_T / N) : 0.0;

    int s = 0; // current segment
    for (int i = 0; i < N; ++i)
    {
        const double ti = t0 + (i + 1) * dt; // [t0+dt, ..., t0+N*dt]
        t_samples[i] = ti;

        while (s + 1 < (int)edges.size() && ti > edges[s + 1])
            ++s;
        const double Ts = T[std::min(s, M - 1)];
        double u = 0.0;
        if (Ts > 0.0)
        {
            u = (ti - edges[s]) / Ts;
            if (u < 0.0)
                u = 0.0;
            if (u > 1.0)
                u = 1.0;
        }
        p_samples[i] = evalBezierQuintic(CP[std::min(s, M - 1)], u);
    }
}

// -----------------------------------------------------------------------------

double SolverLBFGS::evaluateObjective(const VecXd &z) const
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
    const double Amax2 = A_max_ * A_max_;
    const double Jmax2 = J_max_ * J_max_;
    const double Om2 = omege_max_ * omege_max_;
    const double cos_tilt_max = std::cos(tilt_max_rad_);
    const double m = mass_;
    const double g = g_;
    const double eps = 1e-12;
    const Vec3 e3(0.0, 0.0, 1.0);

    // GCOPTER thrust ring parameters
    const double f_mean = 0.5 * (f_min_ + f_max_);
    const double f_radi = 0.5 * std::abs(f_max_ - f_min_);
    const double f_radi2 = f_radi * f_radi;

    // new accumulators
    double J_stat = 0.0, J_vel = 0.0, J_acc = 0.0, J_jmax = 0.0, J_om = 0.0, J_tilt = 0.0, J_thr = 0.0, J_dyn = 0.0;

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

            // Constraints are: Ax < b (Ax - b < 0)
            // With margin Co_, we penalize: Ax - b + Co_ < 0 -> violation = Co_ + (Ax - b)
            if (has_planes)
            {
                for (int h = 0; h < Aseg.rows(); ++h)
                {
                    const double gval = Aseg.row(h).dot(x) - bseg[h];
                    const double viol = gval + Co_;         // <-- key: use Co_
                    J_stat += smoothed_l1(viol, mu) * wseg; // smoothed_l1 zeros negatives
                }
            }

            // --- Velocity: y = ||v||^2 - Vmax^2
            const double yv = v.squaredNorm() - Vmax2;
            J_vel += smoothed_l1(yv, mu) * wseg;

            // --- Acceleration: y = ||a||^2 - Amax^2
            const double ya = a.squaredNorm() - Amax2;
            J_acc += smoothed_l1(ya, mu) * wseg;

            // --- jerk max: y = ||j||^2 - Jmax^2 ---
            const double yj = jrk.squaredNorm() - Jmax2;
            J_jmax += smoothed_l1(yj, mu) * wseg; // add hinge only

            // --- Body-rate from (a, j) with yaẇ=0
            const Vec3 n = a + g * e3;
            const double r = n.norm() + eps; // match Python placement of eps
            const Vec3 b3 = n / r;
            const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - b3 * b3.transpose();
            const Vec3 b3dot = (P * jrk) / r;
            const double yom = b3dot.squaredNorm() - Om2;
            J_om += smoothed_l1(yom, mu) * wseg;

            // --- Tilt: y = cosθ_max - cosθ, cosθ = e3·b3
            const double cos_t = b3.z();
            const double yt = cos_tilt_max - cos_t;
            J_tilt += smoothed_l1(yt, mu) * wseg;

            // --- Thrust ring (GCOPTER): y = ((f - f_mean)^2 - f_radi^2)
            const double f = m * r;
            const double df = f - f_mean;
            const double yth = df * df - f_radi2;
            J_thr += smoothed_l1(yth, mu) * wseg;
        }
    }

    // --- 4) dynamic‐obstacle cost with uniform time sampling (right Riemann)
    //    times: [t0+dt, t0+2dt, ..., t0+N*dt],   dt = total_T / N
    // if (dyn_weight_ > 0.0 && !obstacles_.empty())
    // {
    //     int N = (num_dyn_obst_samples_ > 0 ? num_dyn_obst_samples_ : 10); // set this member or replace as needed
    //     N = std::max(N, 1);
    //     std::vector<double> t_samples;
    //     std::vector<Vec3> p_samples;
    //     sampleRobotPositionsUniform(CP, T, t0_, N, t_samples, p_samples);

    //     const double total_T = t_abs_.back() - t0_;
    //     const double dt = total_T / N;

    //     for (const auto &obs : obstacles_)
    //     {
    //         for (int i = 0; i < N; ++i)
    //         {
    //             const double t = t_samples[i];
    //             const Vec3 &pi = p_samples[i];
    //             const Vec3 ki = obs->eval(t);

    //             const Vec3 diff = pi - ki;
    //             const double d2 = diff.squaredNorm();
    //             const double viol = Cw2_ - d2;

    //             // std::cout << "t: " << t << ", pos: " << pi.transpose() << ", obs: " << ki.transpose() << ", d2: " << d2 << ", viol: " << viol << std::endl;

    //             if (viol > 0.0)
    //                 J_dyn += (viol * viol * viol) * dt; // Riemann scaling
    //         }
    //     }
    // }

    // final weighted sum (add the 4 new terms)
    return time_weight_ * J_time + jerk_weight_ * J_jerk + stat_weight_ * J_stat + dyn_constr_vel_weight_ * J_vel + dyn_constr_acc_weight_ * J_acc + dyn_constr_jerk_weight_ * J_jmax + dyn_constr_bodyrate_weight_ * J_om + dyn_constr_tilt_weight_ * J_tilt + dyn_constr_thrust_weight_ * J_thr + dyn_weight_ * J_dyn;
}

// -----------------------------------------------------------------------------

void SolverLBFGS::computeAnalyticalGrad(const Eigen::VectorXd &z, Eigen::VectorXd &grad) const
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
        const double dT = dT_dtau(z[K_cp_ + s]);
        grad[K_cp_ + s] += time_weight_ * dT;
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

    // 4) dynamic‐obstacle cost gradient
    if (dyn_weight_ > 0.0 && !obstacles_.empty())
    {
        Eigen::VectorXd g_dyn(K_);
        g_dyn.setZero();
        dJ_dyn_dz(z, P, V, A, CP, T, g_dyn);
        grad += dyn_weight_ * g_dyn;
    }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

double SolverLBFGS::evaluateObjectiveAndGradient(
    const Eigen::VectorXd &z,
    Eigen::VectorXd &g) const
{

    // auto t_objective_start = std::chrono::high_resolution_clock::now();

    // 1) compute the objective alone
    double f = evaluateObjective(z);

    // auto t_objective_end = std::chrono::high_resolution_clock::now();
    // auto duration_objective = std::chrono::duration_cast<std::chrono::microseconds>(t_objective_end - t_objective_start);
    // std::cout << "Objective evaluation took: "
    //           << duration_objective.count() / 1000.0 << " ms." << std::endl;

    // auto t_gradient_start = std::chrono::high_resolution_clock::now();

    // 2) compute analytic gradient into g
    g.resize(z.size());
    computeAnalyticalGrad(z, g);

    // auto t_gradient_end = std::chrono::high_resolution_clock::now();
    // auto duration_gradient = std::chrono::duration_cast<std::chrono::microseconds>(t_gradient_end - t_gradient_start);
    // std::cout << "Gradient evaluation took: "
    //           << duration_gradient.count() / 1000.0 << " ms." << std::endl;

    return f;
}

// -----------------------------------------------------------------------------

double SolverLBFGS::evaluateObjectiveAndGradientFused(const Eigen::VectorXd &z, Eigen::VectorXd &grad)
{
    // -------------------------------------------------------------------------
    // Reconstruct once (P, V, A, CP, T) from z  [matches LBFGS layout / scaling]
    // -------------------------------------------------------------------------
    std::vector<Vec3> P, V, A;
    std::vector<std::array<Vec3, 6>> CP;
    std::vector<double> T;
    reconstruct(z, P, V, A, CP, T); // uses τ→T, T̄ decode for interior v̂/â to V/A
    const int M = static_cast<int>(T.size());
    const int knots = M + 1;
    grad.setZero(K_);

    if (M == 0)
        return 0.0;

    // -------------------------------------------------------------------------
    // Objective accumulators (UNWEIGHTED); weights applied in final sum
    // -------------------------------------------------------------------------
    double J_dyn = 0.0;
    double J_time = 0.0;
    double J_panchor = 0.0;
    double J_jerk = 0.0;
    double J_stat = 0.0;
    double J_vel = 0.0;
    double J_acc = 0.0;
    double J_jmax = 0.0;
    double J_om = 0.0;
    double J_tilt = 0.0;
    double J_thr = 0.0;

    // Grad accumulators in (P,V,A,T) space
    std::vector<Vec3> gP(knots, Vec3::Zero());
    std::vector<Vec3> gV(knots, Vec3::Zero());
    std::vector<Vec3> gA(knots, Vec3::Zero());
    std::vector<double> gT(M, 0.0);

    // -------------------------------------------------------------------------
    // Constants / sampling setup
    // -------------------------------------------------------------------------
    const int kappa = (integral_resolution_ > 0 ? integral_resolution_ : 30);
    const double mu = (hinge_mu_ > 0.0 ? hinge_mu_ : 1e-2);
    const double Vmax2 = V_max_ * V_max_;
    const double Amax2 = A_max_ * A_max_;
    const double Jmax2 = J_max_ * J_max_;
    const double Om2 = omege_max_ * omege_max_;
    const double cos_tilt_max = std::cos(tilt_max_rad_);
    const double m = mass_, g = g_;
    const Vec3 e3(0.0, 0.0, 1.0);

    // GCOPTER thrust ring parameters (mean + radius)
    const double f_mean = 0.5 * (f_min_ + f_max_);
    const double f_radi = 0.5 * std::abs(f_max_ - f_min_);
    const double f_radi2 = f_radi * f_radi;

    // Precompute Bernstein basis for j=0..kappa (shared by all segments)
    std::vector<std::array<double, 6>> B5(kappa + 1);
    std::vector<std::array<double, 5>> B4(kappa + 1);
    std::vector<std::array<double, 4>> B3(kappa + 1);
    std::vector<std::array<double, 3>> B2(kappa + 1);
    std::vector<double> wj(kappa + 1, 1.0);
    if (kappa > 0)
    {
        for (int j = 0; j <= kappa; ++j)
        {
            const double tau = static_cast<double>(j) / static_cast<double>(kappa);
            bernstein5(tau, B5[j].data());
            bernstein4(tau, B4[j].data());
            bernstein3(tau, B3[j].data());
            bernstein2(tau, B2[j].data());
            wj[j] = (j == 0 || j == kappa) ? 0.5 : 1.0;
        }
    }

    // -------------------------------------------------------------------------
    // Segment loop
    // -------------------------------------------------------------------------
    for (int s = 0; s < M; ++s)
    {

        const double Ts = T[s];
        const double invT = 1.0 / (Ts + 1e-16);
        const double invT2 = invT * invT;
        const double invT3 = invT2 * invT;
        const double invT4 = invT2 * invT2;

        // ---- (1) Closed-form jerk cost and its gradient wrt CP and T
        //      J_jerk_s = (3600 / T^5) * (||d30||^2 + ||d31||^2 + ||d32||^2)
        const Vec3 d30 = CP[s][3] - 3.0 * CP[s][2] + 3.0 * CP[s][1] - CP[s][0];
        const Vec3 d31 = CP[s][4] - 3.0 * CP[s][3] + 3.0 * CP[s][2] - CP[s][1];
        const Vec3 d32 = CP[s][5] - 3.0 * CP[s][4] + 3.0 * CP[s][3] - CP[s][2];
        const double S = d30.squaredNorm() + d31.squaredNorm() + d32.squaredNorm();
        const double Cj = 3600.0 * std::pow(invT, 5);

        J_jerk += Cj * S;

        // ∂J_jerk/∂CP (vector form): 2 * Cj * [ -d30, (3d30-d31), (-3d30+3d31-d32), (d30-3d31+3d32), (d31-3d32), d32 ]
        std::array<Vec3, 6> gCP_jerk;
        gCP_jerk[0] = -2.0 * Cj * d30;
        gCP_jerk[1] = 2.0 * Cj * (3.0 * d30 - d31);
        gCP_jerk[2] = 2.0 * Cj * (-3.0 * d30 + 3.0 * d31 - d32);
        gCP_jerk[3] = 2.0 * Cj * (d30 - 3.0 * d31 + 3.0 * d32);
        gCP_jerk[4] = 2.0 * Cj * (d31 - 3.0 * d32);
        gCP_jerk[5] = 2.0 * Cj * d32;

        // ∂J_jerk/∂T (factor path): d/dT[(3600 T^{-5})] * S = (-5)*3600*T^{-6} * S
        double gT_jerk = (-5.0) * 3600.0 * std::pow(invT, 6) * S;

        // ---- (2) Sampled terms (static corridor + dynamic limits), accumulate CP/T grads
        std::array<Vec3, 6> gCP_samp; // accumulates all sampled constraints in this segment
        for (int k = 0; k < 6; ++k)
            gCP_samp[k].setZero();
        double gT_samp = 0.0;

        const auto &Aseg = A_stat_[s];
        const auto &bseg = b_stat_[s];
        const bool has_planes = (Aseg.rows() > 0);

        if (kappa > 0)
        {
            const double dt = Ts / static_cast<double>(kappa);

            // Precompute Bezier finite differences in shape-space (independent of tau)
            Vec3 D1[5], D2s[4], D3s[3];
            for (int j = 0; j < 5; ++j)
                D1[j] = CP[s][j + 1] - CP[s][j];
            for (int j = 0; j < 4; ++j)
                D2s[j] = CP[s][j + 2] - 2.0 * CP[s][j + 1] + CP[s][j];
            for (int j = 0; j < 3; ++j)
                D3s[j] = CP[s][j + 3] - 3.0 * CP[s][j + 2] + 3.0 * CP[s][j + 1] - CP[s][j];

            for (int j = 0; j <= kappa; ++j)
            {
                // Basis (cached)
                const auto &b5 = B5[j];
                const auto &b4 = B4[j];
                const auto &b3 = B3[j];
                const auto &b2 = B2[j];
                const double wseg = wj[j] * dt;

                // x, d1s, d2s, d3s in shape-space
                Vec3 x = Vec3::Zero();
                Vec3 d1 = Vec3::Zero();
                Vec3 d2 = Vec3::Zero();
                Vec3 d3 = Vec3::Zero();

                for (int k = 0; k < 6; ++k)
                    x += b5[k] * CP[s][k];
                for (int k = 0; k < 5; ++k)
                    d1 += b4[k] * D1[k];
                for (int k = 0; k < 4; ++k)
                    d2 += b3[k] * D2s[k];
                for (int k = 0; k < 3; ++k)
                    d3 += b2[k] * D3s[k];

                // Convert to physical derivatives
                const Vec3 v = 5.0 * invT * d1;
                const Vec3 acc = 20.0 * invT2 * d2;
                const Vec3 jrk = 60.0 * invT3 * d3;

                // ---------- Static corridor: viol = Co_ + (A x - b)
                if (has_planes)
                {
                    for (int h = 0; h < Aseg.rows(); ++h)
                    {
                        const double gval = Aseg.row(h).dot(x) - bseg[h];
                        const double viol = gval + Co_;
                        J_stat += smoothed_l1(viol, mu) * wseg;

                        if (viol > 0.0)
                        {
                            const double phi_p = smoothed_l1_prime(viol, mu);
                            const Vec3 gx = (phi_p)*Aseg.row(h).transpose() * (stat_weight_ * wseg);
                            for (int k = 0; k < 6; ++k)
                                gCP_samp[k] += b5[k] * gx;

                            // dt-only path
                            gT_samp += stat_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(viol, mu);
                        }
                    }
                }

                // ---------- Velocity: y = ||v||^2 - Vmax^2
                if (dyn_constr_vel_weight_ > 0.0)
                {
                    const double yv = v.squaredNorm() - Vmax2;
                    J_vel += smoothed_l1(yv, mu) * wseg;

                    if (yv > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(yv, mu);
                        const Vec3 gv = (2.0 * phi_p) * v * (dyn_constr_vel_weight_ * wseg); // ∂J/∂v
                        // v = (5/T) d1  ⇒  ∂J/∂d1 = (5/T) gv
                        const Vec3 g1 = gv * invT;
                        for (int k = 0; k < 5; ++k)
                        {
                            const Vec3 G1k = (5.0 * b4[k]) * g1;
                            gCP_samp[k] -= G1k;
                            gCP_samp[k + 1] += G1k;
                        }
                        // dt-only + time-scaling path for v
                        gT_samp += dyn_constr_vel_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(yv, mu);
                        gT_samp += -5.0 * gv.dot(d1) * invT2; // ∂v/∂T = -(5/T^2) d1
                    }
                }

                // ---------- Acceleration: y = ||a||^2 - Amax^2  (present in LBFGS)
                if (dyn_constr_acc_weight_ > 0.0)
                {
                    const double ya = acc.squaredNorm() - Amax2;
                    J_acc += smoothed_l1(ya, mu) * wseg;

                    if (ya > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(ya, mu);
                        const Vec3 ga = (2.0 * phi_p) * acc * (dyn_constr_acc_weight_ * wseg); // ∂J/∂a
                        // a = (20/T^2) d2 ⇒ ∂J/∂d2 = (20/T^2) ga  ⇒  push via 2nd-diff stencil
                        const Vec3 g2 = ga * invT2; // we multiply by 20 in the stencil below
                        for (int k = 0; k < 4; ++k)
                        {
                            const Vec3 G2k = (20.0 * b3[k]) * g2;
                            gCP_samp[k] += G2k;
                            gCP_samp[k + 1] -= 2.0 * G2k;
                            gCP_samp[k + 2] += G2k;
                        }
                        // dt-only + time-scaling path for a
                        gT_samp += dyn_constr_acc_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(ya, mu);
                        gT_samp += -40.0 * ga.dot(d2) * invT3; // ∂a/∂T = -(40/T^3) d2
                    }
                }

                // --- jerk max constraint y = ||j||^2 − Jmax^2
                if (dyn_constr_jerk_weight_ > 0.0)
                {
                    const double yj = jrk.squaredNorm() - Jmax2;
                    J_jmax += smoothed_l1(yj, mu) * wseg;

                    if (yj > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(yj, mu);
                        const Vec3 gj_phys = (2.0 * phi_p) * jrk * (dyn_constr_jerk_weight_ * wseg); // ∂J/∂j

                        // j = (60/T^3) d3  ⇒  ∂J/∂d3 = (60/T^3) * ∂J/∂j
                        const Vec3 g3 = gj_phys * invT3; // 60 is applied in the stencil below

                        // push to CPs via 3rd‑derivative Bézier stencil
                        for (int k = 0; k < 3; ++k)
                        {
                            const Vec3 G3k = (60.0 * b2[k]) * g3;
                            gCP_samp[k] -= G3k;
                            gCP_samp[k + 1] += 3.0 * G3k;
                            gCP_samp[k + 2] -= 3.0 * G3k;
                            gCP_samp[k + 3] += G3k;
                        }

                        // dt‑only term + T‑scaling path
                        gT_samp += dyn_constr_jerk_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(yj, mu);
                        // ∂j/∂T = -(180/T^4) d3
                        gT_samp += -180.0 * gj_phys.dot(d3) * invT4;
                    }
                }

                // ---------- Body-rate: b3dot = (P j)/r,  r=||a+ge3||
                if (dyn_constr_bodyrate_weight_ > 0.0)
                {
                    const Vec3 n = acc + g * e3;
                    const double r = std::sqrt(n.squaredNorm() + 1e-24);
                    const Vec3 b3_rate = n / r;
                    const Eigen::Matrix3d Pmat = Eigen::Matrix3d::Identity() - b3_rate * b3_rate.transpose();
                    const Vec3 b3dot = (Pmat * jrk) / r;
                    const double yom = b3dot.squaredNorm() - Om2;
                    J_om += smoothed_l1(yom, mu) * wseg;

                    if (yom > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(yom, mu);
                        const Vec3 gb = (2.0 * phi_p) * b3dot * (dyn_constr_bodyrate_weight_ * wseg); // ∂J/∂b3dot

                        const Vec3 Pgb = Pmat * gb;
                        const Vec3 Pj = Pmat * jrk;
                        const double alpha = Pj.dot(gb);
                        const double beta = b3_rate.dot(jrk);
                        const double gamma = b3_rate.dot(gb);

                        // ∂J/∂j, ∂J/∂a (physical)
                        const Vec3 g3_phys = Pgb / r;
                        const Vec3 ga_phys = -(beta * Pgb + alpha * b3_rate + gamma * Pj) / (r * r);

                        // Map to shape derivatives and push to CP
                        const Vec3 g2 = ga_phys * invT2; // 20 in stencil
                        const Vec3 g3 = g3_phys * invT3; // 60 in stencil
                        for (int k = 0; k < 4; ++k)
                        {
                            const Vec3 G2k = (20.0 * b3[k]) * g2;
                            gCP_samp[k] += G2k;
                            gCP_samp[k + 1] -= 2.0 * G2k;
                            gCP_samp[k + 2] += G2k;
                        }
                        for (int k = 0; k < 3; ++k)
                        {
                            const Vec3 G3k = (60.0 * b2[k]) * g3;
                            gCP_samp[k] -= G3k;
                            gCP_samp[k + 1] += 3.0 * G3k;
                            gCP_samp[k + 2] -= 3.0 * G3k;
                            gCP_samp[k + 3] += G3k;
                        }

                        // dt-only + time-scaling paths for a and j
                        gT_samp += dyn_constr_bodyrate_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(yom, mu);
                        gT_samp += -40.0 * ga_phys.dot(d2) * invT3;  // a path
                        gT_samp += -180.0 * g3_phys.dot(d3) * invT4; // j path
                    }
                }

                // ---------- Tilt: y = cosθ_max - (e3·b3)
                if (dyn_constr_tilt_weight_ > 0.0)
                {
                    const Vec3 n = acc + g * e3;
                    const double r = std::sqrt(n.squaredNorm() + 1e-24);
                    const Vec3 b3_rate = n / r;
                    const double cos_t = b3_rate.z();
                    const double yt = cos_tilt_max - cos_t;
                    J_tilt += smoothed_l1(yt, mu) * wseg;

                    if (yt > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(yt, mu);
                        const Vec3 dcos_da = (e3 - cos_t * b3_rate) / r;
                        const Vec3 ga_phys = -(phi_p)*dcos_da * (dyn_constr_tilt_weight_ * wseg);

                        const Vec3 g2 = ga_phys * invT2;
                        for (int k = 0; k < 4; ++k)
                        {
                            const Vec3 G2k = (20.0 * b3[k]) * g2;
                            gCP_samp[k] += G2k;
                            gCP_samp[k + 1] -= 2.0 * G2k;
                            gCP_samp[k + 2] += G2k;
                        }

                        gT_samp += dyn_constr_tilt_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(yt, mu);
                        gT_samp += -40.0 * ga_phys.dot(d2) * invT3;
                    }
                }

                // ---------- Thrust ring: φ((f - f_mean)^2 - f_radi^2), f = m * ||a+ge3||
                if (dyn_constr_thrust_weight_ > 0.0)
                {
                    const Vec3 n = acc + g * e3;
                    const double r = std::sqrt(n.squaredNorm() + 1e-24);
                    const double f = m * r;
                    const double df = f - f_mean;
                    const double yth = df * df - f_radi2;
                    J_thr += smoothed_l1(yth, mu) * wseg;

                    if (yth > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(yth, mu);
                        const double dphi_df = 2.0 * df * phi_p;
                        const Vec3 ga_phys = (dyn_constr_thrust_weight_ * wseg) * (dphi_df * m) * (n / r);

                        const Vec3 g2 = ga_phys * invT2;
                        for (int k = 0; k < 4; ++k)
                        {
                            const Vec3 G2k = (20.0 * b3[k]) * g2;
                            gCP_samp[k] += G2k;
                            gCP_samp[k + 1] -= 2.0 * G2k;
                            gCP_samp[k + 2] += G2k;
                        }

                        gT_samp += dyn_constr_thrust_weight_ * (wj[j] / static_cast<double>(kappa)) * smoothed_l1(yth, mu);
                        gT_samp += -40.0 * ga_phys.dot(d2) * invT3;
                    }
                }

            } // samples
        } // if kappa > 0

        // ---- (3) Push jerk+sample CP grads into (P,V,A) and T (via CP(T))
        // Combine jerk (weighted) + sampled CP gradients
        std::array<Vec3, 6> gCP_all;
        for (int i = 0; i < 6; ++i)
            gCP_all[i] = jerk_weight_ * gCP_jerk[i] + gCP_samp[i];

        accumulate_cp_to_pvaT(
            s, gCP_all,
            V, A, T,
            gP, gV, gA, gT);                       // adds the CP(T) path into gT[s]
        gT[s] += jerk_weight_ * gT_jerk + gT_samp; // add the factor/time-scaling paths

        // ---- (4) Time part of objective
        J_time += Ts;
    } // segments

    // ---- Dynamic obstacles (fused): cost + gradient (stable midpoint sampling) ----
    if (dyn_weight_ > 0.0 && !obstacles_.empty())
    {
        const int N = std::max(1, (num_dyn_obst_samples_ > 0 ? num_dyn_obst_samples_ : 10));
        const double total_T = std::accumulate(T.begin(), T.end(), 0.0);

        if (total_T > 0.0)
        {
            const double dt = total_T / static_cast<double>(N);

            // Absolute segment edges: [t0_, t0_+T0, t0_+T0+T1, ...]
            std::vector<double> edges(M + 1, t0_);
            for (int i = 1; i <= M; ++i)
                edges[i] = edges[i - 1] + T[i - 1];

            std::vector<std::array<Vec3, 6>> gCP_dyn(M);
            for (int s = 0; s < M; ++s)
                for (int j = 0; j < 6; ++j)
                    gCP_dyn[s][j].setZero();
            std::vector<double> gT_dyn(M, 0.0);

            constexpr double kEpsT = 1e-9;   // keep time strictly inside horizon
            constexpr double kEpsTau = 1e-9; // keep tau strictly inside (0,1)

            for (const auto &obs : obstacles_)
            {
                if (!obs)
                    continue;

                for (int i = 0; i < N; ++i)
                {
                    // Midpoint in (0,1) → avoids endpoints by construction
                    const double u = (i + 0.5) / static_cast<double>(N);
                    double t_abs = t0_ + u * total_T;
                    t_abs = std::min(std::max(t_abs, edges.front() + kEpsT), edges.back() - kEpsT);

                    // Locate segment and local tau (strictly inside)
                    int s = int(std::upper_bound(edges.begin(), edges.end(), t_abs) - edges.begin()) - 1;
                    s = std::clamp(s, 0, M - 1);

                    const double Ts = T[s];
                    if (Ts <= 0.0)
                        continue;

                    double tau = (t_abs - edges[s]) / Ts;
                    tau = std::min(std::max(tau, kEpsTau), 1.0 - kEpsTau);

                    // Bernstein basis and derivative wrt tau
                    double B[6], dB[6];
                    evalBernstein5(tau, B, dB);

                    // Robot position
                    Vec3 p = Vec3::Zero();
                    for (int j = 0; j < 6; ++j)
                        p += B[j] * CP[s][j];

                    // Obstacle position & velocity
                    const Vec3 k = obs->eval(t_abs);
                    const Vec3 v_obs = obs->velocity(t_abs);

                    // Hinge h = Cw^2 − ||p − k||^2
                    const Vec3 diff = p - k;
                    const double d2 = diff.squaredNorm();
                    const double h = Cw2_ - d2;
                    if (h <= 0.0)
                        continue;

                    const double h2 = h * h;
                    const double h3 = h2 * h;

                    // ∂/∂p (h^3) = -6 h^2 (p - k)
                    const double factor = -6.0 * h2;

                    // (1) CP path: ∂J/∂CP = (∂J/∂p)(∂p/∂CP) dt
                    for (int j = 0; j < 6; ++j)
                        gCP_dyn[s][j] += (factor * dt * B[j]) * diff;

                    // (2) τ(T) path: ∂J/∂T_r via τ
                    Vec3 dp_dtau = Vec3::Zero();
                    for (int j = 0; j < 6; ++j)
                        dp_dtau += dB[j] * CP[s][j];
                    const double term_tau = diff.dot(dp_dtau);

                    for (int r = 0; r < M; ++r)
                    {
                        const double is_lt = (r < s) ? 1.0 : 0.0;
                        const double is_eq = (r == s) ? 1.0 : 0.0;
                        const double dtau_dTr = (u - is_lt - tau * is_eq) / Ts; // correct piecewise formula
                        gT_dyn[r] += (factor * dt) * term_tau * dtau_dTr;
                    }

                    // (3) Obstacle time path: t_abs = t0 + u*total_T ⇒ ∂t_abs/∂T_r = u
                    const double inner_obs = -diff.dot(v_obs); // = (∂J/∂k)·k̇ with the sign folded
                    for (int r = 0; r < M; ++r)
                        gT_dyn[r] += (factor * dt) * inner_obs * u;

                    // (4) dt scaling: dt = total_T/N ⇒ ∂dt/∂T_r = 1/N
                    for (int r = 0; r < M; ++r)
                        gT_dyn[r] += h3 / static_cast<double>(N);

                    // accumulate objective
                    J_dyn += h3 * dt;
                }
            }

            // Push CP(T) chain into (gP,gV,gA) and merge explicit time-path pieces
            for (int s = 0; s < M; ++s)
            {
                std::array<Vec3, 6> tmp;
                for (int j = 0; j < 6; ++j)
                    tmp[j] = dyn_weight_ * gCP_dyn[s][j];
                accumulate_cp_to_pvaT(s, tmp, V, A, T, gP, gV, gA, gT); // also adds CP→T pieces
                gT[s] += dyn_weight_ * gT_dyn[s];
            }
        }
    }

    // --- NEW: prox-to-initial P regularizer ------------------------------------
    // Penalize squared distance from the initial P (first reconstruct).
    // Only adds to gP (no dependency on V, A, or T).
    if (pos_anchor_weight_ > 0.0)
    {
        // If your endpoints are fixed (common), restrict to interior knots so that
        // the objective matches the variables that receive gradients.
        const int i0 = 1;
        const int i1 = knots - 1;
        for (int i = i0; i < i1; ++i)
        {
            const Vec3 diff = P[i] - P_anchor_[i];
            J_panchor += diff.squaredNorm();            // unweighted objective accumulator
            gP[i] += (2.0 * pos_anchor_weight_) * diff; // gradient scaled by weight (∂/∂P_i  ||P_i-P0_i||^2)
        }

        // (Optional) If you prefer averaging to make the scale independent of number of knots:
        // if (i1 > i0) J_panchor /= double(i1 - i0);
    }

    // -------------------------------------------------------------------------
    // Scatter once into z, add T̄ coupling (interior only), then chain ∂T→∂τ
    // -------------------------------------------------------------------------
    // (a) scatter (P,V,A) -> z=[p, v̂, â, τ] with T̄-decoding
    std::vector<double> Tbar(knots);
    build_Tbar(T, Tbar);

    auto vhat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 3], z[base + 4], z[base + 5]);
    };
    auto ahat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 6], z[base + 7], z[base + 8]);
    };

    std::vector<double> gTbar(knots, 0.0);
    for (int i = 1; i < M_; ++i) // interior knots only
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = std::max(1e-12, Tbar[i]);

        grad.segment<3>(base + 0) += gP[i];
        grad.segment<3>(base + 3) += (1.0 / Tb) * gV[i];
        grad.segment<3>(base + 6) += (1.0 / (Tb * Tb)) * gA[i];

        // T̄ coupling (decode) — only where v̂, â exist
        gTbar[i] += -gV[i].dot(vhat_from_z(i)) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_from_z(i)) / (Tb * Tb * Tb);
    }

    // (b) add the T̄→{T_s} contributions to segment ∂J/∂T
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
        gT[s] += gTextra[s];

    // (c) chain to τ via dT/dτ and add time term ∂(time)/∂τ
    for (int s = 0; s < M; ++s)
    {
        const double dTdtau = dT_dtau(z[K_cp_ + s]);
        grad[K_cp_ + s] += gT[s] * dTdtau;
        grad[K_cp_ + s] += time_weight_ * dTdtau; // time term gradient
    }

    // -------------------------------------------------------------------------
    // Final weighted objective (mirrors lbfgs evaluateObjective)
    // -------------------------------------------------------------------------
    const double f =
        dyn_weight_ * J_dyn +
        time_weight_ * J_time +
        pos_anchor_weight_ * J_panchor +
        jerk_weight_ * J_jerk +
        stat_weight_ * J_stat +
        dyn_constr_vel_weight_ * J_vel +
        dyn_constr_acc_weight_ * J_acc +
        dyn_constr_jerk_weight_ * J_jmax +
        dyn_constr_bodyrate_weight_ * J_om +
        dyn_constr_tilt_weight_ * J_tilt +
        dyn_constr_thrust_weight_ * J_thr;

    return f;
}

// -----------------------------------------------------------------------------

int SolverLBFGS::progressCallback(
    void *instance,
    const Eigen::VectorXd &z,
    const Eigen::VectorXd &g,
    const double f,
    const double step,
    const int k,
    const int ls)
{
    // std::cout << "iter=" << k << " f=" << f << " z=" << z.transpose() << " |g|=" << g.norm() << std::endl;
    // std::cout << "iter=" << k << " f=" << f << " |g|=" << g.norm() << std::endl;
    // printf("iter=%d f=%.6f |g|=%.6f\n", k, f, g.norm());
    return 0; // return non‐zero to abort optimization early
}

// -----------------------------------------------------------------------------

double SolverLBFGS::evalObjGradCallback(
    void *instance,
    const Eigen::VectorXd &x,
    Eigen::VectorXd &g)
{
    // dispatch into the instance
    return static_cast<SolverLBFGS *>(instance)->evaluateObjectiveAndGradientFused(x, g);
    // return static_cast<SolverLBFGS *>(instance)->evaluateObjectiveAndGradient(x, g);
}

// -----------------------------------------------------------------------------

int SolverLBFGS::optimize(
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
        &SolverLBFGS::evalObjGradCallback, // objective+gradient callback
        /*stepbound*/ nullptr,
        /*progress*/ &SolverLBFGS::progressCallback,
        // /*progress*/ nullptr,
        const_cast<SolverLBFGS *>(this),
        param);

    // copy result
    z_opt = z;
    return status;
}

// -----------------------------------------------------------------------------

// void SolverLBFGS::setStaticConstraints(
//     const std::vector<LinearConstraint3D> &cons)
// {
//     A_stat_.clear();
//     b_stat_.clear();
//     A_stat_.reserve(cons.size());
//     b_stat_.reserve(cons.size());

//     for (auto const &lc : cons)
//     {
//         // lc.A() is Eigen::Matrix<decimal_t,Dynamic,3>
//         // lc.b() is Eigen::Matrix<decimal_t,Dynamic,1>
//         // cast to double if decimal_t isn’t already double:
//         Eigen::Matrix<double, Eigen::Dynamic, 3> A =
//             lc.A().template cast<double>();
//         Eigen::VectorXd b =
//             lc.b().template cast<double>();

//         A_stat_.push_back(std::move(A));
//         b_stat_.push_back(std::move(b));
//     }
// }

void SolverLBFGS::setStaticConstraints(const std::vector<LinearConstraint3D> &cons)
{
    A_stat_.clear();
    b_stat_.clear();
    A_stat_.reserve(cons.size());
    b_stat_.reserve(cons.size());

    for (auto const &lc : cons)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 3> A = lc.A().template cast<double>();
        Eigen::VectorXd b = lc.b().template cast<double>();
        for (int r = 0; r < A.rows(); ++r)
        {
            const double n = A.row(r).norm();
            if (n > 1e-12)
            {
                A.row(r) /= n;
                b[r] /= n;
            }
        }
        A_stat_.push_back(std::move(A));
        b_stat_.push_back(std::move(b));
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::setStaticConstraintsForSafePath(const std::vector<LinearConstraint3D> &cons)
{

    // safe path should only have one static constraint set 
    // but the path should have two segments (or 3 points) -> so that we can compute trajectory (if it's only two points, the initial and final conditions determine the trajectory completely)

    assert(cons.size() == 1 && "Only one static constraint set is supported for safe path.");

    A_stat_.clear();
    b_stat_.clear();

    for (auto const &lc : cons)
    {
        Eigen::Matrix<double, Eigen::Dynamic, 3> A = lc.A().template cast<double>();
        Eigen::VectorXd b = lc.b().template cast<double>();
        for (int r = 0; r < A.rows(); ++r)
        {
            const double n = A.row(r).norm();
            if (n > 1e-12)
            {
                A.row(r) /= n;
                b[r] /= n;
            }
        }
        A_stat_.push_back(std::move(A));
        A_stat_.push_back(std::move(A));
        b_stat_.push_back(std::move(b));
        b_stat_.push_back(std::move(b));
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::reconstructPVATCPopt(const Eigen::VectorXd &z)
{

    reconstruct(z, P_opt_, V_opt_, A_opt_, CP_opt_, T_opt_);

    // std::cout << "Reconstructed P_opt size: " << P_opt_.size() << std::endl;
    // for (const auto &p : P_opt_)
    // {
    //     std::cout << p.transpose() << std::endl;
    // }
    // std::cout << "Reconstructed V_opt size: " << V_opt_.size() << std::endl;
    // for (const auto &v : V_opt_)
    // {
    //     std::cout << v.transpose() << std::endl;
    // }
    // std::cout << "Reconstructed A_opt size: " << A_opt_.size() << std::endl;
    // for (const auto &a : A_opt_)
    // {
    //     std::cout << a.transpose() << std::endl;
    // }
    // std::cout << "Reconstructed T_opt size: " << T_opt_.size() << std::endl;
    // for (const auto &t : T_opt_)
    // {
    //     std::cout << t << " ";
    // }
}

//------------------------------------------------------------------------------

void SolverLBFGS::dJ_dyn_dz(const VecXd &z,
                            const std::vector<Vec3> &P,
                            const std::vector<Vec3> &V,
                            const std::vector<Vec3> &A,
                            const std::vector<std::array<Vec3, 6>> &CP,
                            const std::vector<double> &T,
                            VecXd &grad) const
{
    const int M = static_cast<int>(T.size());
    if (M == 0 || obstacles_.empty())
        return;

    const int N = std::max(1, num_dyn_obst_samples_);
    const double Cw2 = Cw2_;

    // --- absolute edges and uniform time step ---
    std::vector<double> edges(M + 1);
    edges[0] = t0_;
    for (int i = 1; i <= M; ++i)
        edges[i] = edges[i - 1] + T[i - 1];

    const double total_T = edges.back() - edges.front();
    const double invN = 1.0 / static_cast<double>(N);
    const double dt = total_T * invN;

    const int knots = M + 1;
    std::vector<Vec3> gP(knots, Vec3::Zero());
    std::vector<Vec3> gV(knots, Vec3::Zero());
    std::vector<Vec3> gA(knots, Vec3::Zero());
    std::vector<double> gT(M, 0.0);          // ∂J/∂T (segment times)
    std::vector<double> gT_cp(M, 0.0);       // CP(T) path
    std::vector<std::array<Vec3, 6>> gCP(M); // ∂J/∂CP
    for (int s = 0; s < M; ++s)
        for (int j = 0; j < 6; ++j)
            gCP[s][j].setZero();

    // --- loop obstacles and right-Riemann samples i=1..N ---
    for (const auto &obs : obstacles_)
    {
        for (int i = 1; i <= N; ++i)
        {
            const double i_over_N = static_cast<double>(i) * invN;
            const double t_abs = edges.front() + i_over_N * total_T;

            // find active segment s
            int s = int(std::upper_bound(edges.begin(), edges.end(), t_abs) - edges.begin()) - 1;
            if (s < 0)
                s = 0;
            if (s >= M)
                s = M - 1;

            const double Ts = T[s];
            if (Ts <= 0.0)
                continue;

            // local parameter u in this segment
            double u = (t_abs - edges[s]) / Ts;
            if (u < 0.0)
                u = 0.0;
            else if (u > 1.0)
                u = 1.0;

            double B[6], b4[5];
            bernstein5(u, B);
            bernstein4(u, b4);

            // position at the sample
            Vec3 p = Vec3::Zero();
            for (int j = 0; j < 6; ++j)
                p += B[j] * CP[s][j];

            const Vec3 k = obs->eval(t_abs);
            const Vec3 v_obs = obs->velocity(t_abs);

            const Vec3 diff = p - k;
            const double d2 = diff.squaredNorm();
            if (d2 >= Cw2)
                continue; // outside hinge support

            const double h = (Cw2 - d2);
            const double fac = -6.0 * h * h; // 3*h^2 * (-2*diff)

            // ---- ∂J/∂CP via p(u) ----
            for (int j = 0; j < 6; ++j)
                gCP[s][j] += (fac * dt * B[j]) * diff;

            // ---- u-chain to times: du/dT_r = (i/N - 1_{r<s} - u 1_{r=s}) / Ts ----
            Vec3 dp_du = Vec3::Zero();
            for (int j = 0; j < 5; ++j)
                dp_du += b4[j] * (CP[s][j + 1] - CP[s][j]);
            dp_du *= 5.0;

            const double term_u = diff.dot(dp_du); // enters with fac*dt*term_u*du_dTr
            for (int r = 0; r < M; ++r)
            {
                const double is_lt = (r < s) ? 1.0 : 0.0;
                const double is_eq = (r == s) ? 1.0 : 0.0;
                const double du_dTr = (i_over_N - is_lt - u * is_eq) / Ts;
                gT[r] += (fac * dt) * term_u * du_dTr;
            }

            // ---- obstacle motion term: ∂t_i/∂T_r = i/N,  ∂h/∂t = +2 diff·k̇(t) ----
            const double inner_obs = -diff.dot(v_obs); // so fac*inner_obs = +6 h^2 (diff·k̇)
            for (int r = 0; r < M; ++r)
                gT[r] += (fac * dt) * inner_obs * i_over_N;

            // ---- Riemann weight dt = total_T/N → ∂dt/∂T_r = 1/N ----
            const double hinge3 = h * h * h;
            for (int r = 0; r < M; ++r)
                gT[r] += hinge3 * invN;
        }
    }

    // ---- push CP grads → (P,V,A,T) via CP(T) linkage and add to gT ----
    for (int s = 0; s < M; ++s)
        accumulate_cp_to_pvaT(s, gCP[s], V, A, T, gP, gV, gA, gT_cp);
    for (int s = 0; s < M; ++s)
        gT[s] += gT_cp[s];

    // ---- scatter into z with interior‑only & T̄‑decode, and add T̄‑coupling back to gT ----
    std::vector<double> Tbar(knots);
    build_Tbar(T, Tbar);

    auto vhat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 3], z[base + 4], z[base + 5]);
    };
    auto ahat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 6], z[base + 7], z[base + 8]);
    };

    std::vector<double> gTbar(knots, 0.0);
    for (int i = 1; i < M_; ++i) // interior knots only
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = std::max(1e-12, Tbar[i]);

        grad.segment<3>(base + 0) += gP[i];
        grad.segment<3>(base + 3) += (1.0 / Tb) * gV[i];
        grad.segment<3>(base + 6) += (1.0 / (Tb * Tb)) * gA[i];

        // T̄ decode coupling
        gTbar[i] += -gV[i].dot(vhat_from_z(i)) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_from_z(i)) / (Tb * Tb * Tb);
    }

    // distribute T̄ coupling to segment times and add into gT
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
        gT[s] += gTextra[s];

    // ---- final chain to τ ----
    for (int s = 0; s < M; ++s)
    {
        const double dTdtau = dT_dtau(z[K_cp_ + s]);
        grad[K_cp_ + s] += gT[s] * dTdtau;
    }
}

//------------------------------------------------------------------------------
int SolverLBFGS::findSegment(double ti,
                             const std::vector<double> &T) const
{
    int M = int(T.size());
    static std::vector<double> ends;
    ends.resize(M);
    ends[0] = T[0];
    for (int i = 1; i < M; ++i)
        ends[i] = ends[i - 1] + T[i];
    auto it = std::upper_bound(ends.begin(), ends.end(), ti);
    int s = int(it - ends.begin());
    return (s < M ? s : M - 1);
}

//------------------------------------------------------------------------------

// Fast and stable version
void SolverLBFGS::evalBernstein5(double tau, double B[6], double dB[6]) const
{
    // clamp into [0,1] to avoid tiny overshoots
    if (tau <= 0.0)
        tau = 0.0;
    else if (tau >= 1.0)
        tau = 1.0;

    double u = tau;
    double v = 1.0 - tau;

    // compute all powers with just multiplies
    double u2 = u * u;
    double u3 = u2 * u;
    double u4 = u3 * u;
    double u5 = u4 * u;

    double v2 = v * v;
    double v3 = v2 * v;
    double v4 = v3 * v;
    double v5 = v4 * v;

    // degree‐5 Bernstein basis
    B[0] = v5;
    B[1] = 5.0 * u * v4;
    B[2] = 10.0 * u2 * v3;
    B[3] = 10.0 * u3 * v2;
    B[4] = 5.0 * u4 * v;
    B[5] = u5;

    // degree‐4 Bernstein used for the derivative formula
    //    B4[k] = C(4,k) u^k v^(4−k)
    double B4_0 = v4;
    double B4_1 = 4 * u * v3;
    double B4_2 = 6 * u2 * v2;
    double B4_3 = 4 * u3 * v;
    double B4_4 = u4;

    // dB[j] = 5 * ( B4[j−1] − B4[j] ), with B4[−1]=B4[5]=0
    dB[0] = -5.0 * B4_0;
    dB[1] = 5.0 * (B4_0 - B4_1);
    dB[2] = 5.0 * (B4_1 - B4_2);
    dB[3] = 5.0 * (B4_2 - B4_3);
    dB[4] = 5.0 * (B4_3 - B4_4);
    dB[5] = 5.0 * B4_4;
}

//------------------------------------------------------------------------------

Eigen::Vector3d SolverLBFGS::evalSample(
    const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
    const std::vector<double> &T,
    int sample_i,
    int num_samples,
    int &seg,
    double &tau,
    double &t0) const
{
    // 1) absolute time of this sample
    double total_time = std::accumulate(T.begin(), T.end(), 0.0);
    double alpha = double(sample_i) / double(num_samples - 1);
    double ti = total_time * alpha;

    // 2) find segment index
    seg = findSegment(ti, T);

    // 3) compute t0 and τ
    t0 = (seg > 0)
             ? std::accumulate(T.begin(), T.begin() + seg, 0.0)
             : 0.0;
    double Ti = T[seg];
    tau = (ti - t0) / Ti;

    // 4) Bernstein basis at τ
    double B[6], dB[6];
    evalBernstein5(tau, B, dB);

    // 5) form position
    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    for (int j = 0; j < 6; ++j)
    {
        p += CP[seg][j] * B[j];
    }
    return p;
}

//------------------------------------------------------------------------------

Eigen::Vector3d SolverLBFGS::evalObs(
    const Eigen::Matrix<double, 6, 1> &cx,
    const Eigen::Matrix<double, 6, 1> &cy,
    const Eigen::Matrix<double, 6, 1> &cz,
    double t) const
{
    double tp = 1.0;
    double x = 0, y = 0, zv = 0;
    for (int j = 0; j < 6; ++j)
    {
        x += cx[j] * tp;
        y += cy[j] * tp;
        zv += cz[j] * tp;
        tp *= t;
    }
    return {x, y, zv};
}

// -----------------------------------------------------------------------------

void SolverLBFGS::dJ_jerk_dz(const VecXd &z,
                             const std::vector<Vec3> &P,
                             const std::vector<Vec3> &V,
                             const std::vector<Vec3> &A,
                             const std::vector<std::array<Vec3, 6>> &CP,
                             const std::vector<double> &T,
                             VecXd &grad) const
{
    const int M = static_cast<int>(T.size());
    const int knots = M + 1;

    // accumulators in (P,V,A,T) space
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

        // dJ/dCP via ∂/∂(Δ3)^T(Δ3) chain
        Vec3 dJdCP[6] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(),
                         Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
        auto add = [&](int i, const Vec3 &v)
        { dJdCP[i] += 2.0 * C * v; };
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

        // CP -> (P,V,A) and explicit T terms
        gP[s] += dJdCP[0] + dJdCP[1] + dJdCP[2];
        gV[s] += (Ts / 5.0) * dJdCP[1] + (2.0 * Ts / 5.0) * dJdCP[2];
        gA[s] += (T2 / 20.0) * dJdCP[2];

        gP[s + 1] += dJdCP[3] + dJdCP[4] + dJdCP[5];
        gV[s + 1] += (-2.0 * Ts / 5.0) * dJdCP[3] + (-Ts / 5.0) * dJdCP[4];
        gA[s + 1] += (T2 / 20.0) * dJdCP[3];

        gT[s] += dJdCP[1].dot(V[s] / 5.0);
        gT[s] += dJdCP[2].dot((2.0 / 5.0) * V[s] + (Ts / 10.0) * A[s]);
        gT[s] += dJdCP[3].dot((-2.0 / 5.0) * V[s + 1] + (Ts / 10.0) * A[s + 1]);
        gT[s] += dJdCP[4].dot((-1.0 / 5.0) * V[s + 1]);

        // explicit ∂/∂T of the 1/T^5 factor
        const double S = d30.squaredNorm() + d31.squaredNorm() + d32.squaredNorm();
        gT[s] += (-5.0) * (3600.0 / std::pow(Ts, 6)) * S;
    }

    // (a) scatter (P,V,A) -> z=[p,v̂,â,τ] with the T̄-decoding
    std::vector<double> Tbar;
    build_Tbar(T, Tbar);

    // read v̂, â from z only for interior knots
    auto vhat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 3], z[base + 4], z[base + 5]);
    };
    auto ahat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 6], z[base + 7], z[base + 8]);
    };

    std::vector<double> gTbar(knots, 0.0);
    for (int i = 1; i < M_; ++i) // interior knots only
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = std::max(1e-12, Tbar[i]);

        grad.segment<3>(base + 0) += gP[i];
        grad.segment<3>(base + 3) += (1.0 / Tb) * gV[i];
        grad.segment<3>(base + 6) += (1.0 / (Tb * Tb)) * gA[i];

        // decode → time chain only where v̂, â exist (interior)
        gTbar[i] += -gV[i].dot(vhat_from_z(i)) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_from_z(i)) / (Tb * Tb * Tb);
    }

    // (b) add the T̄→{T_s} contributions to segment ∂J/∂T
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
        gT[s] += gTextra[s];

    // (c) chain to τ via dT/dτ (GCOPTER map) and accumulate into z
    for (int s = 0; s < M; ++s)
    {
        const double dTdtau = dT_dtau(z[K_cp_ + s]);
        grad[K_cp_ + s] += gT[s] * dTdtau;
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::dJ_limits_and_static_dz(const VecXd &z,
                                          const std::vector<Vec3> &P,
                                          const std::vector<Vec3> &V,
                                          const std::vector<Vec3> &A,
                                          const std::vector<std::array<Vec3, 6>> &CP,
                                          const std::vector<double> &T,
                                          VecXd &grad) const
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
    const double Amax2 = A_max_ * A_max_;
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

            // --- Static corridor: Ax - b < 0
            // -> penalize when viol = Co + (Ax - b) > 0
            if (has_planes)
            {
                for (int h = 0; h < Aseg.rows(); ++h)
                {
                    const double gval = Aseg.row(h).dot(x) - bseg[h];
                    const double viol = gval + Co_; // >0 inside the clearance band
                    if (viol > 0.0)
                    {
                        const double phi_p = smoothed_l1_prime(viol, mu);

                        // ∂J/∂x = -φ'(viol) * Aᵀ * (w_stat * wseg)
                        const Vec3 gx = (phi_p)*Aseg.row(h).transpose() * (stat_weight_ * wseg);

                        // x = Σ B5 * CP  ⇒  ∂J/∂CP_j += B5_j * gx
                        for (int k = 0; k < 6; ++k)
                            gCP[k] += B5[k] * gx;

                        // dt-only path: ∂(wseg)/∂T = wj / kappa
                        gT[s] += stat_weight_ * (wj / (double)kappa) * smoothed_l1(viol, mu);
                    }
                }
            }

            // --- Velocity: y = ||v||^2 - Vmax^2 ---
            {
                const double yv = v.squaredNorm() - Vmax2;
                if (yv > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yv, mu);
                    const Vec3 gv = (2.0 * phi_p) * v * (dyn_constr_vel_weight_ * wseg); // ∂J/∂v

                    // v = (5/T) d1  ⇒  ∂J/∂d1 = (5/T) gv  (we do: g1 = gv * (1/T), then multiply by 5 in stencil)
                    const Vec3 g1 = gv * invT;
                    for (int k = 0; k < 5; ++k)
                    {
                        const Vec3 G1k = (5.0 * B4[k]) * g1;
                        gCP[k] -= G1k;
                        gCP[k + 1] += G1k;
                    }

                    // dt-only + T-scale path (v)
                    gT[s] += dyn_constr_vel_weight_ * (wj / (double)kappa) * smoothed_l1(yv, mu);
                    gT[s] += -5.0 * gv.dot(d1s) / (Ts * Ts); // ∂v/∂T = -(5/T^2)d1
                }
            }

            // --- Acceleration: y = ||a||^2 - Amax^2 ---
            {
                // y_a = ||a||^2 - A_max^2  with  a = (20/T^2) d2s
                const double ya = acc.squaredNorm() - Amax2;
                if (ya > 0.0)
                {
                    // ∂J/∂a = 2 a * φ'(y_a) * w
                    const double phi_p = smoothed_l1_prime(ya, mu);
                    const Vec3 ga_phys = (2.0 * phi_p) * acc * (dyn_constr_acc_weight_ * wseg);

                    // Chain to shape-space second derivative:
                    // a = (20/T^2) d2s  ⇒  ∂J/∂d2s = (20/T^2) * ∂J/∂a  = ga_phys / T^2  (then multiply by 20 in stencil)
                    const Vec3 g2_d2s = ga_phys / (Ts * Ts);

                    // Push to CPs: d2s = Σ B3b[k] * (CP_{k+2} - 2 CP_{k+1} + CP_k)
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * g2_d2s; // 20 from a=(20/T^2)d2s
                        gCP[k] += G2k;
                        gCP[k + 1] += -2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }

                    // dt-only path: contribution from the trapezoidal weight
                    gT[s] += dyn_constr_acc_weight_ * (wj / static_cast<double>(kappa)) * smoothed_l1(ya, mu);

                    // T-scale path: a(T) = (20/T^2) d2s  ⇒  ∂a/∂T = -(40/T^3) d2s
                    gT[s] += -40.0 * ga_phys.dot(d2s) / std::pow(Ts, 3);
                }
            }

            // --- Jerk max: y = ||j||^2 - Jmax^2 ---
            {
                const double Jmax2 = J_max_ * J_max_; // define once per function if not already
                const double yj = jrk.squaredNorm() - Jmax2;
                if (yj > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yj, mu);
                    const Vec3 gj_phys = (2.0 * phi_p) * jrk * (dyn_constr_jerk_weight_ * wseg); // ∂J/∂j

                    // j = (60/T^3) d3s  ⇒  ∂J/∂d3s = (60/T^3) * ∂J/∂j
                    const Vec3 g3_d3s = gj_phys / (Ts * Ts * Ts); // 60 is applied in the stencil

                    // push to control points via 3rd‑derivative stencil
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3_d3s;
                        gCP[k] -= G3k;
                        gCP[k + 1] += 3.0 * G3k;
                        gCP[k + 2] -= 3.0 * G3k;
                        gCP[k + 3] += G3k;
                    }

                    // dt‑only contribution (trapezoid weight) + T‑scaling path:
                    gT[s] += dyn_constr_jerk_weight_ * (wj / static_cast<double>(kappa)) * smoothed_l1(yj, mu);
                    // ∂j/∂T = -(180/T^4) d3s  ⇒  (∂J/∂T)_j = (∂J/∂j)·∂j/∂T
                    gT[s] += -180.0 * gj_phys.dot(d3s) / std::pow(Ts, 4);
                }
            }

            // Common pieces for a/j
            const Vec3 n = acc + g * e3;
            const double r = std::sqrt(n.squaredNorm() + 1e-24);
            const Vec3 b3 = n / r;
            const Eigen::Matrix3d Pmat = Eigen::Matrix3d::Identity() - b3 * b3.transpose();

            // --- Body-rate: y = ||b3dot||^2 - Ω^2, with b3dot = (P j)/r ---
            {
                const Vec3 b3dot = (Pmat * jrk) / r;
                const double yom = b3dot.squaredNorm() - Om2;
                if (yom > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yom, mu);
                    const Vec3 gb = (2.0 * phi_p) * b3dot * (dyn_constr_bodyrate_weight_ * wseg); // ∂J/∂b3dot

                    const Vec3 Pgb = Pmat * gb;
                    const Vec3 Pj = Pmat * jrk;
                    const double alpha = Pj.dot(gb);
                    const double beta = b3.dot(jrk);
                    const double gamma = b3.dot(gb);

                    // ∂J/∂j, ∂J/∂a (physical)
                    const Vec3 g3_phys = Pgb / r;
                    const Vec3 ga_phys = -(beta * Pgb + alpha * b3 + gamma * Pj) / (r * r);

                    // Convert to shape-space derivatives and push to CP:
                    const Vec3 g2_d2s = ga_phys / (Ts * Ts);
                    const Vec3 g3_d3s = g3_phys / (Ts * Ts * Ts);

                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * g2_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }
                    for (int k = 0; k < 3; ++k)
                    {
                        const Vec3 G3k = (60.0 * B2[k]) * g3_d3s;
                        gCP[k] -= G3k;
                        gCP[k + 1] += 3.0 * G3k;
                        gCP[k + 2] -= 3.0 * G3k;
                        gCP[k + 3] += G3k;
                    }

                    // dt-only + T-scale paths for a and j
                    gT[s] += dyn_constr_bodyrate_weight_ * (wj / (double)kappa) * smoothed_l1(yom, mu);
                    gT[s] += -40.0 * ga_phys.dot(d2s) / std::pow(Ts, 3);  // a = (20/T^2)d2
                    gT[s] += -180.0 * g3_phys.dot(d3s) / std::pow(Ts, 4); // j = (60/T^3)d3
                }
            }

            // --- Tilt: y = cosθ_max - cosθ,  cosθ = e3·b3 ---
            {
                const double cos_t = b3.z();
                const double yt = cos_tilt_max - cos_t;
                if (yt > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yt, mu);
                    // ∂cos/∂a = (e3 - cosθ b3)/r
                    const Vec3 dcos_da = (e3 - cos_t * b3) / r;
                    const Vec3 ga_phys = -(phi_p)*dcos_da * (dyn_constr_tilt_weight_ * wseg); // ∂J/∂a

                    const Vec3 ga_d2s = ga_phys / (Ts * Ts);
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * ga_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }

                    gT[s] += dyn_constr_tilt_weight_ * (wj / (double)kappa) * smoothed_l1(yt, mu);
                    gT[s] += -40.0 * ga_phys.dot(d2s) / std::pow(Ts, 3);
                }
            }

            // --- Thrust (GCOPTER ring): φ( (f - f_mean)^2 - f_radi^2 ) ---
            {
                const double f = m * r;
                const double df = f - f_mean;
                const double yth = df * df - f_radi2;
                if (yth > 0.0)
                {
                    const double phi_p = smoothed_l1_prime(yth, mu);
                    // ∂g/∂f = 2(f - f_mean),  ∂f/∂a = m b3
                    const double dphi_df = 2.0 * df * phi_p;
                    const Vec3 ga_phys = (dyn_constr_thrust_weight_ * wseg) * (dphi_df * m) * b3;

                    const Vec3 ga_d2s = ga_phys / (Ts * Ts);
                    for (int k = 0; k < 4; ++k)
                    {
                        const Vec3 G2k = (20.0 * B3b[k]) * ga_d2s;
                        gCP[k] += G2k;
                        gCP[k + 1] -= 2.0 * G2k;
                        gCP[k + 2] += G2k;
                    }

                    gT[s] += dyn_constr_thrust_weight_ * (wj / (double)kappa) * smoothed_l1(yth, mu);
                    gT[s] += -40.0 * ga_phys.dot(d2s) / std::pow(Ts, 3);
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

    // (a) scatter (P,V,A) -> z=[p, v̂, â, τ] with T̄-decoding
    std::vector<double> Tbar;
    build_Tbar(T, Tbar);

    // read v̂, â from z only for interior knots
    auto vhat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 3], z[base + 4], z[base + 5]);
    };
    auto ahat_from_z = [&](int i) -> Vec3
    {
        const int base = zcp_offset_for_knot(i, M_);
        return (base < 0) ? Vec3::Zero() : Vec3(z[base + 6], z[base + 7], z[base + 8]);
    };

    std::vector<double> gTbar(knots, 0.0);
    for (int i = 1; i < M_; ++i) // interior knots only
    {
        const int base = zcp_offset_for_knot(i, M_);
        const double Tb = std::max(1e-12, Tbar[i]);

        grad.segment<3>(base + 0) += gP[i];
        grad.segment<3>(base + 3) += (1.0 / Tb) * gV[i];
        grad.segment<3>(base + 6) += (1.0 / (Tb * Tb)) * gA[i];

        // decode → time chain only where v̂, â exist (interior)
        gTbar[i] += -gV[i].dot(vhat_from_z(i)) / (Tb * Tb) - 2.0 * gA[i].dot(ahat_from_z(i)) / (Tb * Tb * Tb);
    }

    // (b) distribute decode extras to segments and add to gT
    std::vector<double> gTextra(M, 0.0);
    distribute_gTbar_to_segments(gTbar, gTextra);
    for (int s = 0; s < M; ++s)
        gT[s] += gTextra[s];

    // (c) chain to τ via dT/dτ and accumulate into z
    for (int s = 0; s < M; ++s)
    {
        const double dTdtau = dT_dtau(z[K_cp_ + s]);
        grad[K_cp_ + s] += gT[s] * dTdtau;
    }
}

// -----------------------------------------------------------------------------

void SolverLBFGS::sanityCheck() const
{

    // If fewer than 2 segments, that's not valid
    if (M_ < 2)
    {
        std::cout << "\033[31mError: M_ must be at least 2.\033[0m" << std::endl;
    }

    if (K_ <= 0)
    {
        std::cout << "\033[31mError: K_ must be positive.\033[0m" << std::endl;
    }

    // Check if the number of segments matches (the number of global waypoints - 1)
    if (M_ != global_wps_.size() - 1)
    {
        // red color
        std::cout << "\033[31mError: Number of segments does not match the number of global waypoints - 1.\033[0m" << std::endl;
        std::cout << "M_ = " << M_ << ", global_wps_.size() = " << global_wps_.size() << std::endl;
    }
}