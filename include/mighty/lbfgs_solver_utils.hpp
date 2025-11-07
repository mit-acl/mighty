/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef LBFGS_SOLVER_UTILS_HPP
#define LBFGS_SOLVER_UTILS_HPP

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <tuple>   // for std::tie
#include <numeric> // for std::accumulate

namespace lbfgs_solver_utils
{

    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using Mat3xN = std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>>;
    using MatXd = Eigen::MatrixXd;
    using VecXd = Eigen::VectorXd;
    using Plane = std::pair<Eigen::RowVector3d, double>;
    using Set = std::vector<Plane>;
    using Sets = std::vector<Set>;
    // A block of parallel half‐spaces: A_block is (n×3), b_block is (n)
    using PlaneBlock = std::pair<Eigen::Matrix<double, Eigen::Dynamic, 3>,
                                 Eigen::VectorXd>;
    // one block per segment‐pair:
    using ConstraintBlocks = std::vector<std::vector<PlaneBlock>>;

    /**
     * @brief Given an absolute time ti and segment durations T,
     *        compute the segment index, the segment start time t0,
     *        and the normalized local time tau = (ti - t0)/T[seg].
     *
     * @param ti    Absolute time along the trajectory
     * @param T     Vector of segment durations (size M)
     * @param[out] seg Index of the segment to which ti belongs (0-based)
     * @param[out] tau Normalized time within that segment ∈ [0,1]
     */
    void get_segment_and_tau(
        double ti,
        const std::vector<double> &T,
        int &seg,
        double &tau);

    // --------------------------------------------------------------

    /**
     * @brief Compute the binomial coefficient "n choose k" robustly.
     * @param n  The total number of items.
     * @param k  The number of items to choose.
     * @return   The value of C(n, k).
     */
    static double binomial(int n, int k);

    // --------------------------------------------------------------

    /**
     * @brief Compute degree-n Bernstein basis and its derivative at parameter τ.
     * @param n    Degree of the Bernstein polynomial.
     * @param tau  Parameter in [0,1].
     * @param[out] B   Output vector of length n+1, where
     *                 B[j] = C(n, j) · τ^j · (1−τ)^(n−j).
     * @param[out] dB  Output vector of length n+1, the derivative dB[j]/dτ.
     */
    void bernstein_basis_and_derivative(
        int n,
        double tau,
        std::vector<double> &B,
        std::vector<double> &dB);

    // --------------------------------------------------------------

    /**
     * @brief Bernstein basis polynomial of degree n, index j, at parameter tau.
     * @param n     Degree of the Bernstein polynomial.
     * @param j     Basis index (0 ≤ j ≤ n).
     * @param tau   Parameter in [0,1].
     * @return      Value of B_{n,j}(tau).
     */
    double bernstein(int n, int j, double tau);

    // --------------------------------------------------------------

    /**
     * @brief Evaluate a composite quintic Bézier trajectory and its first two derivatives.
     *
     * Given M segments each with 6 control‐points in R³ and their durations T,
     * sample positions, velocities, and accelerations at num_samples equally
     * spaced times over [0, sum(T)].
     *
     * @param CP          vector of M arrays of 6 control‐points (each Eigen::Vector3d)
     * @param T           segment durations (size M)
     * @param num_samples number of samples along the whole trajectory
     * @param[out] pts    output positions (size = num_samples)
     * @param[out] vels   output velocities (size = num_samples)
     * @param[out] accs   output accelerations (size = num_samples)
     */
    void eval_traj_and_derivs(
        const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
        const std::vector<double> &T,
        int num_samples,
        std::vector<Eigen::Vector3d> &pts,
        std::vector<Eigen::Vector3d> &vels,
        std::vector<Eigen::Vector3d> &accs);

    // --------------------------------------------------------------

    /**
     * @brief Numerically stable SoftPlus activation: log(1+exp(x)).
     * @param x input
     * @return softplus(x)
     */
    double softplus(double x);

    // -----------------------------------------------------------------------------

    /**
     * @brief Sample positions along a composite quintic Bézier in ℝ³.
     *
     * Given M segments each with 6 control‐points in CP, and segment durations T,
     * uniformly sample num_samples points over the total trajectory time.
     *
     * @param CP          Control points: vector of M segments, each an array of 6 Vec3
     * @param T           Segment durations (size M)
     * @param num_samples Number of samples to take along [0, sum(T)]
     * @param[out] pts    Output positions (will be resized to num_samples)
     */
    void eval_traj(
        const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
        const std::vector<double> &T,
        int num_samples,
        std::vector<Eigen::Vector3d> &pts);

    // -----------------------------------------------------------------------------

    /**
     * @brief Fit a single‐segment quintic polynomial in ℝ³:
     *        p(t)=∑_{i=0}^5 c_i t^i
     *        s.t. p(0)=x0, p′(0)=v0, p″(0)=a0,
     *             p(T)=xf, p′(T)=vf, p″(T)=af.
     *
     * @param x0  Initial position
     * @param v0  Initial velocity
     * @param a0  Initial acceleration
     * @param xf  Final position
     * @param vf  Final velocity
     * @param af  Final acceleration
     * @param T   Duration of the segment
     * @param[out] cx  Coefficients c₀…c₅ for the x‐component
     * @param[out] cy  Coefficients c₀…c₅ for the y‐component
     * @param[out] cz  Coefficients c₀…c₅ for the z‐component
     */
    void fit_quintic(
        const Eigen::Vector3d &x0,
        const Eigen::Vector3d &v0,
        const Eigen::Vector3d &a0,
        const Eigen::Vector3d &xf,
        const Eigen::Vector3d &vf,
        const Eigen::Vector3d &af,
        double T,
        Eigen::Matrix<double, 6, 1> &cx,
        Eigen::Matrix<double, 6, 1> &cy,
        Eigen::Matrix<double, 6, 1> &cz);

    // -----------------------------------------------------------------------------

    /**
     * @brief Evaluate a degree-5 polynomial obstacle trajectory in ℝ³.
     *
     * Given a vector of times `t` and three length-6 coefficient vectors
     * `cx, cy, cz`, computes for each time the point
     *   p(t[i]) = (∑ cx[j] t[i]^j, ∑ cy[j] t[i]^j, ∑ cz[j] t[i]^j).
     *
     * @param t    Input time samples (size N)
     * @param cx   Polynomial coefficients for x(t) (6×1)
     * @param cy   Polynomial coefficients for y(t) (6×1)
     * @param cz   Polynomial coefficients for z(t) (6×1)
     * @param[out] pts  Output positions (size N)
     */
    void f_obs_poly(
        const std::vector<double> &t,
        const Eigen::Matrix<double, 6, 1> &cx,
        const Eigen::Matrix<double, 6, 1> &cy,
        const Eigen::Matrix<double, 6, 1> &cz,
        std::vector<Eigen::Vector3d> &pts);

    // -----------------------------------------------------------------------------

} // namespace lbfgs_solver_utils

#endif // DYNUS_LBFGS_SOLVER_UTILS_HPP
