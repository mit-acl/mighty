/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef LBFGS_SOLVER_HPP
#define LBFGS_SOLVER_HPP

// #define EIGEN_DONT_ALIGN
// #define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

//
#include <Eigen/Dense>
#include <vector>
#include <mighty/lbfgs_solver_utils.hpp>
#include <numeric>
#include <cmath> // for std::log, std::pow
#include <mighty/lbfgs.hpp>
#include <omp.h>
#include <set>

// ROS
#include <mighty/mighty_type.hpp>
#include <decomp_util/ellipsoid_decomp.h>
#include <decomp_util/seed_decomp.h>

namespace lbfgs
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

    struct planner_params_t
    {
        bool verbose = false;     // verbosity level
        double V_max = 2.0;       // max velocity
        double A_max = 10.0;      // max acceleration
        double J_max = 100.0;     // max jerk
        int num_perturbation = 8; // number of perturbations for initial guesses
        double r_max = 1.0;       // perturbation radius for initial guesses
        double time_weight = 1.0;
        double pos_anchor_weight = 1.0;
        double dyn_weight = 1.0;
        double stat_weight = 1.0;
        double jerk_weight = 10.0;
        double dyn_constr_bodyrate_weight = 10.0; // weight for body rate constraints
        double dyn_constr_tilt_weight = 10.0;     // weight for tilt constraints
        double dyn_constr_thrust_weight = 10.0;   // weight for
        double dyn_constr_vel_weight = 10.0;
        double dyn_constr_acc_weight = 10.0;
        double dyn_constr_jerk_weight = 0.0;
        int num_dyn_obst_samples = 10;         // Number of dynamic obstacle samples
        double Co = 0.5;                       // for static obstacle avoidance
        double Cw = 1.0;                       // for dynamic obstacle avoidance
        double BIG = 1e8;                      // a large constant for static constraints
        double dc = 0.01;                      // descretiation constant
        double init_turn_bf = 15.0;            // initial turn buffer in degrees
        int integral_resolution = 30;          // resolution for the integral in the optimization
        double hinge_mu = 1e-2;                // hinge mu for the optimization
        double omega_max = 1.0;                // max body rate in rad/s
        double tilt_max_rad = 0.6;             // max tilt in radians (e.g., 35° in rad)
        double f_min = 0.0;                    // min thrust in N
        double f_max = 20.0;                   // max thrust in N
        double mass = 1.0;                     // mass in kg
        double g = 9.81;                       // gravity in m/s^2
    };

    class SolverLBFGS
    {

    public:
        SolverLBFGS();
        ~SolverLBFGS();

        // -----------------------------------------------------------------------------

        // ------------------------------
        // Functions
        // ------------------------------

        void pushWaypointsByStaticCorridor(std::vector<Vec3> &wps);
        void getGlobalPath(vec_Vecf<3> &global_path);

        // -----------------------------------------------------------------------------
        /**
         * @brief Build a collection of initial‐guess vectors z₀ (and the
         *        corresponding control‐points & times) by sampling
         *        systematic perturbations of a “global” waypoint list.
         *
         * @param global_wps       Base waypoints (size M+1)
         * @param joint_indices    Which indices to perturb
         * @param r_max            Perturbation radius
         * @param constraint_sets  For each joint, list of (A_row,b) enforcing A·pt>b
         * @param[out] list_z0     Packed [freeCPs; sigmas] for each waypoint set
         * @param[out] list_cp0    Full reconstructed CP arrays for each z₀
         * @param[out] list_t0     Reconstructed segment‐time vectors for each z₀
         */
        void buildInitialGuesses(const ConstraintBlocks &constraint_sets);

        // -----------------------------------------------------------------------------

        /**
         * @brief Systematically perturb joint waypoints in 3D by sampling
         *        eight equally spaced directions in the bisector plane.
         *
         * @param constraint_sets  For each joint, list of (A_row,b) pairs enforcing A_row·pt > b
         * @param[out] samples     Output array of waypoint sets (each size K)
         */
        std::vector<std::vector<Vec3>>
        sample_systematic_perturbed_waypoints(
            const ConstraintBlocks &constraint_sets) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Evaluate the full trajectory objective (time + dynamic‐avoidance +
         *        static barrier + jerk + dynamic‐constraint penalties).
         *
         * @param z  Decision vector [p₀,v₀,a₀ … p_M,v_M,a_M; σ₀…σ_{M-1}]
         * @return   Scalar value of the objective at z
         */
        double evaluateObjective(const VecXd &z) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Compute the analytic gradient ∇J(z).
         * @param z      Decision variable vector [free CPs; slack times] (size K).
         * @param[out] grad  Gradient vector (size K), overwritten on output.
         */
        void computeAnalyticalGrad(
            const Eigen::VectorXd &z,
            Eigen::VectorXd &grad) const;

        // -----------------------------------------------------------------------------

        double evaluateObjectiveAndGradientFused(const Eigen::VectorXd& z, Eigen::VectorXd& grad);

        // -----------------------------------------------------------------------------

        /**
         * @brief Callback for L-BFGS: evaluate cost & analytic gradient in one shot.
         *
         * @param z Current decision variables
         * @param g Output gradient
         * @return cost = J(z)
         */
        double evaluateObjectiveAndGradient(const Eigen::VectorXd &z,
                                            Eigen::VectorXd &g) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief [Debug] Progress callback for L-BFGS, called at each iteration.
         *
         * @param instance  Pointer to the SolverLBFGS instance (unused).
         * @param x         Current decision variables.
         * @param g         Current gradient.
         * @param fx        Current objective value.
         * @param step      Step size (unused).
         * @param k         Iteration number.
         * @param ls        Line search status (unused).
         * @return          0 to continue, non-zero to stop.
         */
        static int progressCallback(
            void *instance,
            const Eigen::VectorXd &z,
            const Eigen::VectorXd &g,
            const double f,
            const double step,
            const int k,
            const int ls);

        // -----------------------------------------------------------------------------

        /**
         * @brief Optimize the trajectory using L-BFGS.
         *
         * @param z0      Initial guess for decision variables (size K).
         * @param z_opt   Output optimized decision variables (size K).
         * @param f_opt   Output final objective value.
         * @param param   L-BFGS parameters.
         * @return       Status code: 0 on success, negative on error.
         */
        int optimize(
            const Eigen::VectorXd &z0,
            Eigen::VectorXd &z_opt,
            double &f_opt,
            const lbfgs::lbfgs_parameter_t &param) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Reconstruct knot states and Hermite control points from decision vector z.
         *
         * @param z    Decision vector of length K_cp_ + K_sig_:
         *             [p₀,v₀,a₀, …, p_M,v_M,a_M, σ₀, …, σ_{M-1}]
         * @param P    (output) Knot positions, size M_+1
         * @param V    (output) Knot velocities, size M_+1
         * @param A    (output) Knot accelerations, size M_+1
         * @param CP   (output) Hermite control points, size M_×6
         * @param T    (output) Segment durations, size M_
         */
        void reconstruct(
            const VecXd &z,
            std::vector<Vec3> &P,
            std::vector<Vec3> &V,
            std::vector<Vec3> &A,
            std::vector<std::array<Vec3, 6>> &CP,
            std::vector<double> &T) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Pack knot states and durations into the decision vector z.
         *
         * Inverse of reconstruct(...).  Given
         *   P: size M_+1 of knot positions,
         *   V: size M_+1 of knot velocities,
         *   A: size M_+1 of knot accelerations,
         *   T: size M_   of segment durations,
         * this writes
         *   z = [p₀,v₀,a₀, …, p_M,v_M,a_M, σ₀…σ_{M-1}]
         * where σ_s = log(T[s]).
         *
         * @param P   Knot positions, length M_+1
         * @param V   Knot velocities, length M_+1
         * @param A   Knot accelerations, length M_+1
         * @param T   Segment durations, length M_
         * @param z   (output) decision vector, length = 9*(M_+1) + M_
         */
        void packDecisionVariables(
            const std::vector<Vec3> &P,
            const std::vector<Vec3> &V,
            const std::vector<Vec3> &A,
            const std::vector<double> &T,
            VecXd &z) const;

        // -----------------------------------------------------------------------------

        /// Overwrite the static‐constraint arrays from a list of LinearConstraint3D

        void setStaticConstraints(const std::vector<LinearConstraint3D> &cons);

        void setStaticConstraintsForSafePath(const std::vector<LinearConstraint3D> &cons);

        // -----------------------------------------------------------------------------
        // ROS-related functions
        // -----------------------------------------------------------------------------

        /**
         * @brief Sanity check function to verify the solver's setup.
         */
        void sanityCheck() const;

        /**
         * @brief initialize the solver with default parameters.
         */
        void initializeSolver(const planner_params_t &params);

        bool replaceGlobalPathWithCorridorShortest(bool do_shortcut, double smooth_eps);

        /**
         *
         */
        void prepareSolverForReplan(double t0,
                                    const vec_Vec3f &global_wps,
                                    const std::vector<LinearConstraint3D> &safe_corridor,
                                    const std::vector<std::shared_ptr<dynTraj>> &obstacles,
                                    const state &initial_state,
                                    const state &goal_state,
                                    double &initial_guess_computation_time,
                                    bool use_for_safe_path = false,
                                    bool use_multiple_initial_guesses = true);

        /**
         * @brief Set the initial guess waypoints for the trajectory.
         */
        static std::array<Vec3, 6> computeQuinticCP(
            const Vec3 &P0, const Vec3 &V0, const Vec3 &A0,
            const Vec3 &P1, const Vec3 &V1, const Vec3 &A1,
            double T);
        static void solveMinJerkVelAcc(
            const std::vector<Vec3> &wps, // size M+1
            const std::vector<double> &T, // size M
            const Vec3 &v0, const Vec3 &a0,
            const Vec3 &vf, const Vec3 &af,
            std::vector<Vec3> &V,  // out size M+1
            std::vector<Vec3> &A); // out size M+1

        void assemble_H_b_acc_only(
            const std::vector<Vec3> &wps,   // size M+1
            const std::vector<double> &T,   // size M
            const std::vector<Vec3> &V,     // size M+1 (ALL velocities fixed)
            const Vec3 &a0, const Vec3 &af, // fixed endpoint accelerations
            Eigen::MatrixXd &H,             // out: (N×N)
            Eigen::MatrixXd &b              // out: (N×3)
        );

        void solveMinJerkAccOnlyClosedForm(
            const std::vector<Vec3> &wps,   // size M+1
            const std::vector<double> &T,   // size M
            const std::vector<Vec3> &V,     // size M+1 (fixed)
            const Vec3 &a0, const Vec3 &af, // fixed endpoint accelerations
            std::vector<Vec3> &A_out        // out: size M+1
        );

        void findInitialGuess(std::vector<double> &T,
                              std::vector<Vec3> &V,
                              std::vector<Vec3> &A);

        /// Elementwise softplus on a Vec3
        static Vec3 softplus_vec(const Vec3 &v)
        {
            return v.unaryExpr([](double x)
                               { return std::log1p(std::exp(x)); });
        }

        /// Scalar softplus
        static double softplus(double x)
        {
            return std::log1p(std::exp(x));
        }

        // Helper functions for closed-form min-jerk
        static inline Eigen::Matrix4d K_r(double T);
        static inline std::array<double, 4> k_r(double T, double P0, double P1);

        /**
         * Assemble the block‐tridiagonal Hessian H and RHS b for the closed-form min-jerk.
         *
         * Inputs:
         *   wps  – M+1 waypoints (Vec3)
         *   T    – M segment times
         *   v0,a0,vf,af – endpoint velocity/accel
         *
         * Outputs:
         *   H (2N×2N), b (2N×3), with N = M−1
         */
        static void assemble_H_b(
            const std::vector<Vec3> &wps,
            const std::vector<double> &T,
            const Vec3 &v0, const Vec3 &a0,
            const Vec3 &vf, const Vec3 &af,
            Eigen::MatrixXd &H,
            Eigen::MatrixXd &b);

        /**
         * @brief Get the goal setpoints
         */
        void getGoalSetpoints(std::vector<state> &goal_setpoints);

        /**
         * @brief
         */
        void getPieceWisePol(PieceWiseQuinticPol &pwp);

        /**
         *
         */
        void getControlPoints(std::vector<Eigen::Matrix<double, 3, 6>> &cps);

        /**
         *
         */
        void reconstructPVATCPopt(const Eigen::VectorXd &z);

        /**
         *
         */
        inline StateDeriv evalStateDeriv(int segment, double tau) const;
        inline double getPosDouble(int segment, double tau, int axis) const;
        inline double getVelDouble(int segment, double tau, int axis) const;
        inline double getAccelDouble(int segment, double tau, int axis) const;
        inline double getJerkDouble(int segment, double tau, int axis) const;

        // -----------------------------------------------------------------------------

        // -------------------------------
        // Getters
        // -------------------------------

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the number of decision variables (free CPs + slack times).
         * @return Number of decision variables K = 3F + M.
         */
        int getNumDecisionVariables() const
        {
            return K_;
        }

        // -----------------------------------------------------------------------------

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the number of segments in the trajectory.
         * @return Number of segments M = global_wps_.size() - 1.
         */
        int getNumSegments() const
        {
            return static_cast<int>(global_wps_.size()) - 1;
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Get list_z0, the initial guesses for decision variables.
         * @return List of initial guesses for decision variables (size N).
         */
        const std::vector<Eigen::VectorXd> &getInitialGuesses() const
        {
            return list_z0_;
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the initial guess waypoints.
         * @return List of initial guess waypoints (size N).
         */
        const std::vector<std::vector<Vec3>> &getInitialGuessWaypoints() const
        {
            return initial_guess_wps_;
        }

    private:

        // obstacles
        std::vector<std::shared_ptr<dynTraj>> obstacles_;


        /**
         * @brief Static wrapper for evaluateObjectiveAndGradient.
         * @return The objective value at x; also fills g.
         */
        static double evalObjGradCallback(
            void *instance,
            const Eigen::VectorXd &x,
            Eigen::VectorXd &g);

        // -----------------------------------------------------------------------------

        /**
         * @brief Analytic gradient of the dynamic‐obstacle cost at knot times.
         *
         * Cost: J_dyn = ∑_obs ∑_{i=1..M-1} [max( Cw² – ‖E·(P_i–K_i)‖², 0)]³
         *
         * @param z     Decision vector [p₀,v₀,a₀,…,p_M,v_M,a_M, σ₀…σ_{M-1}]
         * @param P     Knot positions (size M+1).
         * @param T     Segment durations (size M).
         * @param grad  (output) gradient ∇_z J_dyn, same size as z
         */
        void dJ_dyn_dz(const VecXd &z,
                       const std::vector<Vec3> &P,
                       const std::vector<Vec3> &V,
                       const std::vector<Vec3> &A,
                       const std::vector<std::array<Vec3, 6>> &CP,
                       const std::vector<double> &T,
                       VecXd &grad) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Analytical gradient of the static‐barrier cost
         *   J_stat = ∑_{s=0..M-1} ∑_{j∈relevant(s)} −log( A_stat[s]·P[s,j] − b_stat[s] )
         *   w.r.t. decision vector z = [free CPs; slack times].
         *
         * @param CP         Bezier control points, M segments each with 6×3 points.
         * @param P     Knot positions (size M+1).
         * @param T     Segment durations (size M).
         * @param grad  (output) gradient ∇_z J_stat, same size as z
         */
        void dJ_stat_dz(const VecXd &z,
                        const std::vector<Vec3> &P,
                        const std::vector<Vec3> &V,
                        const std::vector<Vec3> &A,
                        const std::vector<std::array<Vec3, 6>> &CP,
                        const std::vector<double> &T,
                        VecXd &grad) const;

        // -----------------------------------------------------------------------------

        /**
         * @brief Analytical gradient of the jerk penalty
         *        J_jerk = ∑ₛ (3600/Tₛ⁵) ∑_{m=0..2} ‖Δ³P[s,m]‖²
         * w.r.t. z = [free CPs; slack times].
         */
        void dJ_jerk_dz(const VecXd &z,
                        const std::vector<Vec3> &P,
                        const std::vector<Vec3> &V,
                        const std::vector<Vec3> &A,
                        const std::vector<std::array<Vec3, 6>> &CP,
                        const std::vector<double> &T,
                        VecXd &grad) const;

        // -----------------------------------------------------------------------------

        void dJ_limits_and_static_dz(const VecXd &z,
                                     const std::vector<Vec3> &P,
                                     const std::vector<Vec3> &V,
                                     const std::vector<Vec3> &A,
                                     const std::vector<std::array<Vec3, 6>> &CP,
                                     const std::vector<double> &T,
                                     VecXd &grad) const;

        // -----------------------------------------------------------------------------

        double centralDiff(const VecXd &z,
                           const VecXd &d,
                           double eps_base);

        void checkGradDirectional(const VecXd &z0,
                                  int num_dirs /*=8*/,
                                  double eps /*=1e-6*/,
                                  unsigned seed /*=0*/);

        void checkGradCoordinates(const VecXd &z0, int max_coords, double eps, unsigned seed);

    private:
        /// Find which segment index s contains time t_i
        int findSegment(double ti, const std::vector<double> &T) const;

        /// Evaluate Bernstein basis (degree=5) and derivative at tau
        void evalBernstein5(double tau, double B[6], double dB[6]) const;

        /**
         * @brief Evaluate one trajectory sample u_i and its segment/tau.
         *
         * @param CP           full control–point array [M][6]
         * @param T            segment durations (size M)
         * @param sample_i     index i in [0..num_samples-1]
         * @param num_samples total number of samples
         * @param[out] seg    segment index containing sample i
         * @param[out] tau    normalized time in [0,1] within that segment
         * @param[out] t0     absolute start time of that segment
         * @return            u_i position in ℝ³
         */
        Eigen::Vector3d evalSample(
            const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
            const std::vector<double> &T,
            int sample_i,
            int num_samples,
            int &seg,
            double &tau,
            double &t0) const;

        /// Evaluate obstacle position at time t
        Eigen::Vector3d evalObs(const Eigen::Matrix<double, 6, 1> &cx,
                                const Eigen::Matrix<double, 6, 1> &cy,
                                const Eigen::Matrix<double, 6, 1> &cz,
                                double t) const;

        /**
         * Build a small Jacobian ∂u_i/∂z for one sample:
         *   - only the free‐CP columns in segment s (3 per cp)
         *   - plus the two slack columns σ_s and σ_{s+1}
         */
        // in lbfgs_solver.hpp, under your class SolverLBFGS
        void buildSparseDuDzSample(
            const Eigen::VectorXd &z,
            const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
            const std::vector<double> &T,
            int sample_i,
            int num_samples,
            int seg,
            double tau,
            const std::vector<std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>>> &dPdz,
            const Eigen::MatrixXd &dTdz,
            std::vector<int> &cols,
            Eigen::MatrixXd &Jsmall) const;

        /**
         * Build a small Jacobian ∂k_i/∂z for one sample:
         *   - only the M slack columns (but will be the same α factor for each)
         */
        void buildSparseDkDzSample(
            const Eigen::VectorXd &z,
            const std::vector<double> &T,
            int sample_i,
            int num_samples,
            const Eigen::Matrix<double, 6, 1> &cx,
            const Eigen::Matrix<double, 6, 1> &cy,
            const Eigen::Matrix<double, 6, 1> &cz,
            std::vector<int> &cols,
            Eigen::MatrixXd &Ksmall) const;

        void buildSparseDkDzSample(
            const Eigen::VectorXd &z,
            const std::vector<double> &T,
            int sample_i,
            int num_samples,
            const PieceWisePol &pw,
            std::vector<int> &cols,
            Eigen::MatrixXd &Ksmall) const;

        // -----------------------------------------------------------------------------
        // Helper functions for ROS
        // -----------------------------------------------------------------------------

        inline int nCk5(int n, int k)
        {
            // Pascal’s triangle up to row 5
            static constexpr int C[6][6] = {
                {1, 0, 0, 0, 0, 0},
                {1, 1, 0, 0, 0, 0},
                {1, 2, 1, 0, 0, 0},
                {1, 3, 3, 1, 0, 0},
                {1, 4, 6, 4, 1, 0},
                {1, 5, 10, 10, 5, 1}};
            return C[n][k];
        }

    protected:
        // ------------------------------
        // Parameters
        // ------------------------------

        // Parameters
        int M_;               // Number of segments (this now is global_wps_.size() - 1)
        int K_;               // Number of decision variables (CPs + slack times)
        int K_cp_;            // Number of control points
        int K_sig_;           // Number of slack times (K_sig_ = M_)
        bool verbose_{false}; // Verbosity flag

        // Ego-agent boundary conditions
        Vec3 x0_;
        Vec3 v0_;
        Vec3 a0_;
        Vec3 xf_;
        Vec3 vf_;
        Vec3 af_;

        // static constraints
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3>> A_stat_; // each row is a plane normal
        std::vector<Eigen::VectorXd> b_stat_;                          // each entry is the corresponding offset

        // Global waypoints
        std::vector<Vec3> global_wps_;

        // Dynamic constraints
        double V_max_;
        double A_max_;
        double J_max_;

        // Perturbation parameters
        int num_perturbation_;
        double r_max_;

        // sampling & smoothing
        int integral_resolution_ = 30;
        double hinge_mu_ = 1e-2;

        // vehicle/dynamics limits and constants
        double omege_max_ = 6.0;                    // rad/s
        double tilt_max_rad_ = M_PI * 35.0 / 180.0; // 35 deg
        double f_min_ = 0.0;                        // N
        double f_max_ = 20.0;                       // N
        double mass_ = 1.0;                         // kg
        double g_ = 9.81;                           // m/s^2

        // Position anchoring

        // Optimization weights and settings
        double time_weight_;
        double pos_anchor_weight_;
        double dyn_weight_;
        double stat_weight_;
        double jerk_weight_;
        double dyn_constr_vel_weight_;
        double dyn_constr_acc_weight_;
        double dyn_constr_jerk_weight_;
        double dyn_constr_bodyrate_weight_ = 1.0;
        double dyn_constr_tilt_weight_ = 1.0;
        double dyn_constr_thrust_weight_ = 1.0;
        int num_dyn_obst_samples_; // Number of dynamic obstacle samples
        double Co_; // for static obstacle avoidance
        double Cw_;
        double V_min_;      // used for initial guess
        double turn_buf_;   // used for initial guess
        double turn_span_;  // used for initial guess
        double cos_thresh_; // used for initial guess

        // Position anchoring
        std::vector<Vec3> P_anchor_;

        // Constants
        double Cw2_;

        // Optimization variables
        mutable std::vector<double> t_abs_;

        // Lists for initial guesses
        std::vector<std::vector<Vec3>> initial_guess_wps_;
        std::vector<Eigen::VectorXd> list_z0_;

        // Precomputed binomial coefficients for Bernstein basis
        const int C5_[6] = {1, 5, 10, 10, 5, 1}; // currently not used
        const int C4_[5] = {1, 4, 6, 4, 1};      // currently not used

        // Large value
        double BIG_;
        // Small value
        double eps_ = 1e-6;

        // ROS parameters
        std::vector<state> goal_setpoints_; // vector of goal setpoints
        double dc_;                         // small value
        std::vector<std::array<Vec3, 6>> CP_opt_;
        std::vector<double> T_opt_;
        std::vector<Vec3> P_opt_, V_opt_, A_opt_;
        double t0_{0.0}; // start time of the trajectory

    }; // class SolverLBFGS

} // namespace lbfgs

#endif