/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <mighty/lbfgs_solver_utils.hpp>

namespace lbfgs_solver_utils
{

    // --------------------------------------------------------------

    void get_segment_and_tau(
        double ti,
        const std::vector<double> &T,
        int &seg,
        double &tau)
    {
        int M = static_cast<int>(T.size());
        // build cumulative end times
        std::vector<double> ends(M);
        if (M > 0)
        {
            ends[0] = T[0];
            for (int i = 1; i < M; ++i)
            {
                ends[i] = ends[i - 1] + T[i];
            }
        }

        // find first end time > ti
        auto it = std::upper_bound(ends.begin(), ends.end(), ti);
        seg = static_cast<int>(it - ends.begin());
        if (seg >= M)
            seg = M - 1; // clamp to last segment

        // segment start time
        double t0 = (seg > 0 ? ends[seg - 1] : 0.0);

        // normalized local time
        if (seg >= 0 && seg < M && T[seg] > 0.0)
        {
            tau = (ti - t0) / T[seg];
        }
        else
        {
            tau = 0.0;
        }
    }

    // --------------------------------------------------------------

    static double binomial(int n, int k)
    {
        if (k < 0 || k > n)
            return 0.0;
        k = std::min(k, n - k);
        double c = 1.0;
        for (int i = 1; i <= k; ++i)
        {
            c *= double(n - k + i) / double(i);
        }
        return c;
    }

    // --------------------------------------------------------------

    void bernstein_basis_and_derivative(
        int n,
        double tau,
        std::vector<double> &B,
        std::vector<double> &dB)
    {
        // Resize outputs
        B.resize(n + 1);
        dB.resize(n + 1);

        // Compute degree-n Bernstein basis
        for (int j = 0; j <= n; ++j)
        {
            double coeff = binomial(n, j);
            B[j] = coeff * std::pow(tau, j) * std::pow(1.0 - tau, n - j);
        }

        if (n == 0)
        {
            dB[0] = 0.0;
            return;
        }

        // Compute degree-(n−1) basis for recurrence
        std::vector<double> Bm1(n);
        for (int j = 0; j < n; ++j)
        {
            double coeff = binomial(n - 1, j);
            Bm1[j] = coeff * std::pow(tau, j) * std::pow(1.0 - tau, (n - 1) - j);
        }

        // Use recurrence: dB_j = n [ B_{n-1,j-1} − B_{n-1,j} ]
        for (int j = 0; j <= n; ++j)
        {
            double v = 0.0;
            if (j > 0)
                v += n * Bm1[j - 1];
            if (j < n)
                v -= n * Bm1[j];
            dB[j] = v;
        }
    }

    // --------------------------------------------------------------

    double bernstein(int n, int j, double tau)
    {
        return binomial(n, j) * std::pow(tau, j) * std::pow(1.0 - tau, n - j);
    }

    // --------------------------------------------------------------

    void eval_traj_and_derivs(
        const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
        const std::vector<double> &T,
        int num_samples,
        std::vector<Eigen::Vector3d> &pts,
        std::vector<Eigen::Vector3d> &vels,
        std::vector<Eigen::Vector3d> &accs)
    {
        int M = static_cast<int>(CP.size());
        // 1) total trajectory time
        double total_time = std::accumulate(T.begin(), T.end(), 0.0);

        // 2) absolute sample times
        std::vector<double> t_abs(num_samples);
        for (int i = 0; i < num_samples; ++i)
        {
            t_abs[i] = total_time * double(i) / double(num_samples - 1);
        }

        // 3) cumulative segment end-times
        std::vector<double> ends(M);
        ends[0] = T[0];
        for (int i = 1; i < M; ++i)
        {
            ends[i] = ends[i - 1] + T[i];
        }

        // 4) prepare outputs
        pts.resize(num_samples);
        vels.resize(num_samples);
        accs.resize(num_samples);

        // 5) sample each time
        for (int k = 0; k < num_samples; ++k)
        {
            double ti = t_abs[k];
            // find segment index
            int seg = int(std::upper_bound(ends.begin(), ends.end(), ti) - ends.begin());
            if (seg >= M)
                seg = M - 1;

            double t0 = (seg > 0 ? ends[seg - 1] : 0.0);
            double Ti = T[seg];
            double tau = (ti - t0) / Ti;

            // grab the six control‐points for this segment
            const auto &P = CP[seg];

            // position
            Eigen::Vector3d p = Eigen::Vector3d::Zero();
            for (int j = 0; j < 6; ++j)
            {
                p += P[j] * bernstein(5, j, tau);
            }
            pts[k] = p;

            // velocity
            Eigen::Vector3d v = Eigen::Vector3d::Zero();
            for (int j = 0; j < 5; ++j)
            {
                v += (P[j + 1] - P[j]) * bernstein(4, j, tau);
            }
            vels[k] = (5.0 / Ti) * v;

            // acceleration
            Eigen::Vector3d a = Eigen::Vector3d::Zero();
            for (int j = 0; j < 4; ++j)
            {
                a += (P[j + 2] - 2.0 * P[j + 1] + P[j]) * bernstein(3, j, tau);
            }
            accs[k] = (20.0 / (Ti * Ti)) * a;
        }
    }

    // --------------------------------------------------------------

    double softplus(double x)
    {
        if (x > 0.0)
        {
            // x + log(1 + exp(-x)) is stable when x is large
            return x + std::log1p(std::exp(-x));
        }
        else
        {
            // log(1 + exp(x)) is stable when x is negative
            return std::log1p(std::exp(x));
        }
    }

    void eval_traj(
        const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
        const std::vector<double> &T,
        int num_samples,
        std::vector<Eigen::Vector3d> &pts)
    {
        const int M = static_cast<int>(CP.size());
        // total trajectory time
        double total_time = std::accumulate(T.begin(), T.end(), 0.0);

        // absolute sample times
        std::vector<double> t_abs(num_samples);
        for (int i = 0; i < num_samples; ++i)
        {
            t_abs[i] = (i == num_samples - 1)
                           ? total_time
                           : total_time * double(i) / double(num_samples - 1);
        }

        pts.resize(num_samples);

        for (int i = 0; i < num_samples; ++i)
        {
            double ti = t_abs[i];

            // 1) find segment index and normalized time tau
            int seg; double tau;
            get_segment_and_tau(ti, T, seg, tau);

            // 2) compute Bernstein basis of degree 5 at tau
            std::vector<double> B(6), dB(6);
            bernstein_basis_and_derivative(5, tau, B, dB);

            // 3) evaluate position via p = sum_{j=0..5} B[j] * CP[seg][j]
            Eigen::Vector3d p = Eigen::Vector3d::Zero();
            for (int j = 0; j <= 5; ++j)
            {
                p += B[j] * CP[seg][j];
            }
            pts[i] = p;
        }
    }

    // --------------------------------------------------------------

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
        Eigen::Matrix<double, 6, 1> &cz)
    {
        // Build the 6×6 system A · c = b for each dimension
        Eigen::Matrix<double, 6, 6> A = Eigen::Matrix<double, 6, 6>::Zero();
        Eigen::Matrix<double, 6, 1> bx = Eigen::Matrix<double, 6, 1>::Zero();
        Eigen::Matrix<double, 6, 1> by = Eigen::Matrix<double, 6, 1>::Zero();
        Eigen::Matrix<double, 6, 1> bz = Eigen::Matrix<double, 6, 1>::Zero();

        // 1) Boundary at t=0: p(0)=x0, p'(0)=v0, p''(0)=a0
        A(0, 0) = 1.0;
        bx(0) = x0.x();
        by(0) = x0.y();
        bz(0) = x0.z();

        A(1, 1) = 1.0;
        bx(1) = v0.x();
        by(1) = v0.y();
        bz(1) = v0.z();

        A(2, 2) = 2.0;
        bx(2) = a0.x();
        by(2) = a0.y();
        bz(2) = a0.z();

        // 2) Position at t=T: p(T)=xf
        for (int i = 0; i < 6; ++i)
        {
            A(3, i) = std::pow(T, i);
        }
        bx(3) = xf.x();
        by(3) = xf.y();
        bz(3) = xf.z();

        // 3) Velocity at t=T: p'(T)=vf
        for (int i = 1; i < 6; ++i)
        {
            A(4, i) = i * std::pow(T, i - 1);
        }
        bx(4) = vf.x();
        by(4) = vf.y();
        bz(4) = vf.z();

        // 4) Acceleration at t=T: p''(T)=af
        for (int i = 2; i < 6; ++i)
        {
            A(5, i) = i * (i - 1) * std::pow(T, i - 2);
        }
        bx(5) = af.x();
        by(5) = af.y();
        bz(5) = af.z();

        // Solve A c = b
        // You can pick whichever decomposition you like; here we use full pivot LU:
        cx = A.fullPivLu().solve(bx);
        cy = A.fullPivLu().solve(by);
        cz = A.fullPivLu().solve(bz);
    }

    // --------------------------------------------------------------

    void f_obs_poly(
        const std::vector<double> &t,
        const Eigen::Matrix<double, 6, 1> &cx,
        const Eigen::Matrix<double, 6, 1> &cy,
        const Eigen::Matrix<double, 6, 1> &cz,
        std::vector<Eigen::Vector3d> &pts)
    {
        const int N = static_cast<int>(t.size());
        pts.resize(N);

        for (int i = 0; i < N; ++i)
        {
            double ti = t[i];
            double ti_pow = 1.0;

            // accumulate monomials t^0 … t^5
            Eigen::Vector3d p(0.0, 0.0, 0.0);
            for (int j = 0; j < 6; ++j)
            {
                if (j > 0)
                    ti_pow *= ti;
                p.x() += cx[j] * ti_pow;
                p.y() += cy[j] * ti_pow;
                p.z() += cz[j] * ti_pow;
            }
            pts[i] = p;
        }
    }

}