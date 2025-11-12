#pragma once
/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#define EIGEN_DONT_ALIGN
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include <Eigen/Dense>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <gcopter/flatness.hpp>
#include <gcopter/lbfgs.hpp>
#include <mighty/mighty_type.hpp>
#include <decomp_util/ellipsoid_decomp.h>
#include <decomp_util/seed_decomp.h>
#include <decomp_util/ellipsoid_decomp.h>
#include <decomp_util/seed_decomp.h>

namespace lbfgs
{

    using Vec3 = Eigen::Vector3d;
    using MatXd = Eigen::MatrixXd;
    using VecXd = Eigen::VectorXd;
    using Plane = std::pair<Eigen::RowVector3d, double>;
    using Set = std::vector<Plane>;
    using Sets = std::vector<Set>;
    using PlaneBlock = std::pair<Eigen::Matrix<double, Eigen::Dynamic, 3>, Eigen::VectorXd>;
    using ConstraintBlocks = std::vector<std::vector<PlaneBlock>>;
    using PolyhedronH = Eigen::Matrix<double, Eigen::Dynamic, 4>;
    using PolyhedronV = Eigen::Matrix3Xd;
    using PolyhedraH = std::vector<PolyhedronH>;
    using PolyhedraV = std::vector<PolyhedronV>;

    /**
     * @brief Static parameters and penalty weights for the MIGHTY L-BFGS solver.
     *
     * Values are read once via initializeSolver() and then treated as constants
     * during optimization/re-planning.
     */
    struct planner_params_t
    {
        // Vehicle / limits
        double V_max{5.0};        ///< Max allowed speed [m/s]
        double omega_max{6.0};    ///< Max body-rate magnitude [rad/s]
        double tilt_max_rad{0.6}; ///< Max tilt angle [rad]
        double f_min{0.0};        ///< Min thrust [N]
        double f_max{20.0};       ///< Max thrust [N]
        double mass{1.0};         ///< Vehicle mass [kg]
        double g{9.81};           ///< Gravity [m/s^2]

        // Objective weights
        double time_weight{1.0};                ///< ΣT weight
        double jerk_weight{1.0};                ///< ∫||j||^2 weight (closed form)
        double stat_weight{1.0};                ///< Static corridor clearance weight
        double dyn_constr_vel_weight{1.0};      ///< Velocity-limit hinge weight
        double dyn_constr_acc_weight{0.0};      ///< (Reserved) Acceleration-limit weight
        double dyn_constr_bodyrate_weight{1.0}; ///< Body-rate hinge weight
        double dyn_constr_tilt_weight{1.0};     ///< Tilt hinge weight
        double dyn_constr_thrust_weight{1.0};   ///< Thrust-ring hinge weight

        // Corridor / smoothing
        double Co{0.0};              ///< Clearance distance used inside corridor hinge
        double BIG{1e9};             ///< Large constant for internal constraints (if needed)
        double dc{0.01};             ///< Export sampling step for setpoints [s]
        int integral_resolution{30}; ///< # of samples per segment for time integrals
        double hinge_mu{1e-2};       ///< C¹/C² hinge smoothing parameter

        // Turn buffer used to reduce waypoint speeds at sharp corners (degrees)
        double init_turn_bf{15.0};
    };

    /**
     * @brief Piecewise quintic (power-basis) representation of the optimized path.
     *
     * One set of 6 coefficients per axis and per segment, plus absolute breakpoints.
     */
    struct PieceWiseQuinticPol
    {
        std::vector<double> times;                        ///< Absolute breakpoints t_0..t_M (size M+1)
        std::vector<Eigen::Matrix<double, 6, 1>> coeff_x; ///< x(t) coefficients per segment
        std::vector<Eigen::Matrix<double, 6, 1>> coeff_y; ///< y(t) coefficients per segment
        std::vector<Eigen::Matrix<double, 6, 1>> coeff_z; ///< z(t) coefficients per segment

        /// Clear all buffers.
        inline void clear()
        {
            times.clear();
            coeff_x.clear();
            coeff_y.clear();
            coeff_z.clear();
        }
    };

    /**
     * @brief L‑BFGS trajectory optimizer that parameterizes a multi-segment
     *        path with quintic Hermite end-states and equivalent Bézier control points.
     *
     * Key features:
     *  - Uses corridor planes for static avoidance (sampled hinge).
     *  - Enforces/penalizes speed, body-rate, tilt, and thrust-ring via flatness map.
     *  - Supports fused objective+gradient for efficiency.
     */
    class SolverMIGHTY
    {
    public:
        /// Default constructor.
        SolverMIGHTY();

        /// Destructor.
        ~SolverMIGHTY();

        /**
         * @brief Initialize mission‑wide constants and penalty weights.
         *
         * @param params            All scalar weights/limits for the objective and constraints.
         * @param physical_params   Vector of physical parameters forwarded to the flatness map:
         *                          [ mass, g, horizDrag, vertDrag, parasDrag, speedEps ].
         */
        void initializeSolver(const planner_params_t &params,
                              const Eigen::VectorXd &physical_params);

        /**
         * @brief Prepare a replan: ingest corridor & waypoints, decode the seam vector @p initial_xi,
         *        build internal offsets and an initial guess @c z0 (stored in getInitialGuesses()).
         *
         * The expected V‑polytope list has size 2*M−1 with SFC cells at even indices and seam cells at odd indices.
         * Seams are parameterized by a unit vector @f$q@f$ whose squared components define barycentric weights.
         *
         * @param t0          Absolute time stamp for the start of the trajectory.
         * @param global_wps  Start→goal waypoints (M+1 knots).
         * @param safe_corridor Static corridor planes per segment (A x ≥ b).
         * @param initial_xi  Packed seam variables q from upstream (either seams‑only or "full").
         * @param init_times  Upstream initial segment times (not required for packing τ; kept for completeness).
         * @param initial_state  Start boundary state (pos/vel/acc).
         * @param goal_state     Goal boundary state (pos/vel/acc).
         * @param initial_guess_computation_time  Filled with the initial guess wall time [ms].
         * @param vPolys       Corridor V‑polytopes: [SFC_0, seam_0, SFC_1, seam_1, ..., SFC_M].
         * @param use_multiple_initial_guesses  If true, may seed multiple z0 candidates.
         */
        void prepareSolverForReplan(double t0,
                                    const vec_Vec3f &global_wps,
                                    const std::vector<LinearConstraint3D> &safe_corridor,
                                    const Eigen::VectorXd &initial_xi,
                                    const Eigen::VectorXd &init_times,
                                    const state &initial_state,
                                    const state &goal_state,
                                    double &initial_guess_computation_time,
                                    const PolyhedraV &vPolys,
                                    bool use_multiple_initial_guesses);

        /**
         * @brief Extract per‑seam q‑blocks from a packed @p initial_xi vector.
         *
         * Accepts either the "seams‑only" packing or the "full" packing (endpoints + seams).
         * The geometry (number of columns in each seam polytope) determines expected sizes.
         *
         * @param initial_xi Input packed xi.
         * @param q_seams    Output vector of seam blocks @f$q_s@f$ (size M−1).
         * @return true if @p initial_xi length matches the corridor expectations; false otherwise.
         */
        bool extractSeamQFromInitialXi(const Eigen::VectorXd &initial_xi,
                                       std::vector<Eigen::VectorXd> &q_seams) const;

        /**
         * @brief Pack decision vector in the internal layout
         * @code
         * z = [ xi(seams),  v̂_i (i=1..M-1),  â_i (i=1..M-1),  τ_s (s=0..M-1) ]
         * @endcode
         * Endpoints (P0, PM, V0, A0, VM, AM) are fixed and not part of @p z.
         */
        void packDecisionVariables(const Eigen::VectorXd &initial_xi,
                                   const std::vector<Vec3> &P,
                                   const std::vector<Vec3> &V,
                                   const std::vector<Vec3> &A,
                                   const std::vector<double> &T,
                                   VecXd &z) const;

        /**
         * @brief Decode @p z into Hermite end-states and Bézier control points.
         *
         * @param z  Decision vector.
         * @param P  Output knot positions (size M+1).
         * @param V  Output knot velocities (size M+1).
         * @param A  Output knot accelerations (size M+1).
         * @param CP Output per‑segment 6 Bézier control points.
         * @param T  Output per‑segment durations.
         */
        void reconstruct(const VecXd &z,
                         std::vector<Vec3> &P,
                         std::vector<Vec3> &V,
                         std::vector<Vec3> &A,
                         std::vector<std::array<Vec3, 6>> &CP,
                         std::vector<double> &T) const;

        /**
         * @brief Build a quick, direction‑aware initial guess for (T, V, A) over @c global_wps_.
         *
         * Speeds at interior waypoints are reduced based on a turn buffer; accelerations are
         * obtained by solving a small banded min‑jerk system in closed form.
         */
        void findInitialGuess(std::vector<double> &T,
                              std::vector<Vec3> &V,
                              std::vector<Vec3> &A);

        /**
         * @brief Closed‑form generation of quintic Bézier control points for a single segment,
         *        given Hermite boundary values.
         */
        std::array<Vec3, 6> computeQuinticCP(const Vec3 &P0, const Vec3 &V0, const Vec3 &A0,
                                             const Vec3 &P1, const Vec3 &V1, const Vec3 &A1,
                                             double T);

        /**
         * @brief Solve the coupled min‑jerk system with interior velocities and accelerations free.
         *
         * @param wps  Knots (positions).
         * @param T    Segment durations.
         * @param v0,a0  Fixed start velocity/acceleration.
         * @param vf,af  Fixed final velocity/acceleration.
         * @param V,A  Outputs at all knots.
         */
        void solveMinJerkVelAcc(const std::vector<Vec3> &wps,
                                const std::vector<double> &T,
                                const Vec3 &v0, const Vec3 &a0,
                                const Vec3 &vf, const Vec3 &af,
                                std::vector<Vec3> &V,
                                std::vector<Vec3> &A);

        /**
         * @brief Solve the reduced min‑jerk system with interior accelerations free
         *        and all velocities fixed.
         *
         * @param wps Knots (positions).
         * @param T   Segment durations.
         * @param V   Fixed velocities at all knots.
         * @param a0,af Fixed endpoint accelerations.
         * @param A_out Output accelerations at all knots.
         */
        void solveMinJerkAccOnlyClosedForm(const std::vector<Vec3> &wps,
                                           const std::vector<double> &T,
                                           const std::vector<Vec3> &V,
                                           const Vec3 &a0, const Vec3 &af,
                                           std::vector<Vec3> &A_out);

        /**
         * @brief Fused evaluation of objective and analytical gradient.
         *
         * Reconstructs once, then accumulates time, jerk (closed form) and
         * sampled constraints (corridor, velocity, body-rate, tilt, thrust).
         *
         * @param z     Decision vector.
         * @param grad  Output gradient (same size as @p z).
         * @return Objective value @f$f(z)@f$.
         */
        double evaluateObjectiveAndGradient(const Eigen::VectorXd &z,
                                            Eigen::VectorXd &grad);

        /**
         * @brief Objective-only evaluation (slower than the fused path).
         *
         * @param z Decision vector.
         * @return Objective value.
         */
        double evaluateObjective(const VecXd &z);

        /**
         * @brief Standalone analytical gradient (alternative path to the fused one).
         *
         * @param z    Decision vector.
         * @param grad Output gradient.
         */
        void computeAnalyticalGrad(const Eigen::VectorXd &z,
                                   Eigen::VectorXd &grad);

        /**
         * @brief Run L‑BFGS from an initial guess.
         *
         * @param z0     Initial decision vector.
         * @param z_opt  Output: optimized vector.
         * @param f_opt  Output: final objective.
         * @param param  L‑BFGS parameters (mem size, LS, tolerances, etc.).
         * @return Library termination status.
         */
        int optimize(const Eigen::VectorXd &z0,
                     Eigen::VectorXd &z_opt,
                     double &f_opt,
                     const lbfgs_parameter_t &param) const;

        /**
         * @brief Run L‑BFGS with a wall‑clock budget; if time elapses, returns the latest iterate.
         *
         * @param z0     Initial decision vector.
         * @param z_opt  Output: optimized or latest vector.
         * @param f_opt  Output: final or latest objective.
         * @param param  L‑BFGS parameters.
         * @param wall_budget Time budget.
         * @return Library termination status.
         */
        int optimize(const Eigen::VectorXd &z0,
                     Eigen::VectorXd &z_opt,
                     double &f_opt,
                     const lbfgs_parameter_t &param,
                     std::chrono::milliseconds wall_budget) const;

        /**
         * @brief Replace static corridor planes (A x ≥ b) used by the sampling penalties.
         */
        void setStaticConstraints(const std::vector<LinearConstraint3D> &cons);

        /**
         * @brief Return uniformly-spaced setpoints over the current @c T_opt_/P_opt_/V_opt_/A_opt_.
         *
         * The output is sampled in steps of @c dc_ (see planner_params_t::dc).
         */
        void getGoalSetpoints(std::vector<state> &goal_setpoints);

        /**
         * @brief Convert the optimized Hermite form to piecewise quintic power basis.
         */
        void getPieceWisePol(PieceWiseQuinticPol &pwp);

        /**
         * @brief Evaluate position/velocity/acceleration/jerk at segment @p s and τ∈[0,1].
         */
        StateDeriv evalStateDeriv(int s, double tau) const;

        /**
         * @brief Central-difference directional derivative for gradient checking.
         *
         * @param z        Base point.
         * @param d        Direction (will be normalized inside the routine).
         * @param eps_base Base step size before scaling.
         * @return (f(z+εd)−f(z−εd)) / (2ε).
         */
        double centralDiff(const VecXd &z, const VecXd &d, double eps_base);

        /**
         * @brief Random directional gradient check.
         *
         * @param z0       Point at which to check.
         * @param num_dirs Number of random directions.
         * @param eps      Finite-difference step.
         * @param seed     RNG seed.
         */
        void checkGradDirectional(const VecXd &z0, int num_dirs, double eps, unsigned seed);

        /**
         * @brief Random coordinate-wise gradient check over up to @p max_coords entries.
         */
        void checkGradCoordinates(const VecXd &z0, int max_coords, double eps, unsigned seed);

        /**
         * @brief Access the prepared initial guesses (usually a single @c z0).
         */
        inline const std::vector<Eigen::VectorXd> &getInitialGuesses() const { return list_z0_; }

    private:
        // ---- L-BFGS glue (static callbacks) ----
        static int progressCallback(void *instance,
                                    const Eigen::VectorXd &z,
                                    const Eigen::VectorXd &g,
                                    double f,
                                    double step,
                                    int k,
                                    int ls);

        static double evalObjGradCallback(void *instance,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g);

        // ---- Bernstein cache (mutable: filled in const contexts) ----
        struct BernsteinCache
        {
            int kappa{-1};
            std::vector<std::array<double, 6>> B5;
            std::vector<std::array<double, 5>> B4;
            std::vector<std::array<double, 4>> B3;
            std::vector<std::array<double, 3>> B2;
            std::vector<double> W;
        };
        mutable BernsteinCache bern_; // mutable so const methods can use it
        void ensureBernsteinCache(int kappa) const;

        // Timing and last-iterate bookkeeping (marked mutable so optimize can remain const)
        mutable std::chrono::steady_clock::time_point opt_start_;
        mutable std::chrono::steady_clock::time_point opt_deadline_;
        mutable bool have_deadline_ = false;
        mutable bool timed_out_ = false;

        mutable Eigen::VectorXd last_z_;
        mutable double last_f_ = std::numeric_limits<double>::infinity();

        // ---- Internal helpers (implemented in .cpp) ----
        inline void reconstruct_inplace(const VecXd &z) const
        {
            reconstruct(z, P_s_, V_s_, A_s_, CP_s_, T_s_);
        }
        void scatterPosGrads(const std::vector<Vec3> &gP, const VecXd &z, VecXd &grad) const;
        void seamBackward(int seam, const Eigen::VectorXd &q, const Eigen::Vector3d &gp, Eigen::VectorXd &gq) const;

        Eigen::Matrix4d K_r(double T);
        std::array<double, 4> k_r(double T, double P0, double P1);

        void assemble_H_b(const std::vector<Vec3> &wps,
                          const std::vector<double> &T,
                          const Vec3 &v0, const Vec3 &a0,
                          const Vec3 &vf, const Vec3 &af,
                          Eigen::MatrixXd &H,
                          Eigen::MatrixXd &b);

        void assemble_H_b_acc_only(const std::vector<Vec3> &wps,
                                   const std::vector<double> &T,
                                   const std::vector<Vec3> &V,
                                   const Vec3 &a0, const Vec3 &af,
                                   Eigen::MatrixXd &H,
                                   Eigen::MatrixXd &b);

        void dJ_jerk_dz(const VecXd &z,
                        const std::vector<Vec3> &P,
                        const std::vector<Vec3> &V,
                        const std::vector<Vec3> &A,
                        const std::vector<std::array<Vec3, 6>> &CP,
                        const std::vector<double> &T,
                        VecXd &grad) const;

        void dJ_limits_and_static_dz(const VecXd &z,
                                     const std::vector<Vec3> &P,
                                     const std::vector<Vec3> &V,
                                     const std::vector<Vec3> &A,
                                     const std::vector<std::array<Vec3, 6>> &CP,
                                     const std::vector<double> &T,
                                     VecXd &grad);

        // Offsets for interior derivative blocks
        inline int vhatOffset(int i) const { return offVhat_ + 3 * (i - 1); }
        inline int ahatOffset(int i) const { return offAhat_ + 3 * (i - 1); }
        inline Eigen::Vector3d getVhat(const VecXd &z, int i) const { return z.segment<3>(vhatOffset(i)); }
        inline Eigen::Vector3d getAhat(const VecXd &z, int i) const { return z.segment<3>(ahatOffset(i)); }

    private:
        // ------------------------------
        // Parameters
        // ------------------------------

        // Parameters
        int M_;               // Number of segments (this now is global_wps_.size() - 1)
        int K_;               // Number of decision variables (CPs + slack times)

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

        // sampling & smoothing
        int integral_resolution_ = 30;
        double hinge_mu_ = 1e-2;
        
        // vehicle/dynamics limits and constants
        double V_max_;
        double omege_max_ = 6.0;                    // rad/s
        double tilt_max_rad_ = M_PI * 35.0 / 180.0; // 35 deg
        double f_min_ = 0.0;                        // N
        double f_max_ = 20.0;                       // N
        double mass_ = 1.0;                         // kg
        double g_ = 9.81;                           // m/s^2

        // parameters
        flatness::FlatnessMap flatmap_;

        // Optimization weights and settings
        double time_weight_;
        double dyn_weight_;
        double stat_weight_;
        double jerk_weight_;
        double dyn_constr_vel_weight_;
        double dyn_constr_bodyrate_weight_ = 1.0;
        double dyn_constr_tilt_weight_ = 1.0;
        double dyn_constr_thrust_weight_ = 1.0;
        double Co_;        // for static obstacle avoidance
        double V_min_;     // used for initial guess
        double turn_buf_;  // used for initial guess
        double turn_span_; // used for initial guess

        // Lists for initial guesses
        std::vector<Eigen::VectorXd> list_z0_;

        // Large value
        double BIG_;
        // Small value
        double eps_ = 1e-6;

        // ROS parameters
        double dc_; // small value
        std::vector<double> T_opt_;
        std::vector<Vec3> P_opt_, V_opt_, A_opt_;
        double t0_{0.0}; // start time of the trajectory

        // >>> add to SolverLBFGS class (private:)
        PolyhedraV vPolys_OB_;         // [poly0, inter01, poly1, inter12, ..., polyM]
        Eigen::VectorXi seamSizes_;    // k_i = vPolys_OB_[2*i-1].cols() for i=1..M-1 (intersections)
        std::vector<int> seamOffsets_; // offsets in z for each q-block

        // z layout offsets (computed per problem)
        int offP0_ = 0;   // 3 doubles
        int offXi_ = 0;   // start of all q blocks
        int offPM_ = 0;   // 3 doubles
        int offVhat_ = 0; // 3*(M+1)
        int offAhat_ = 0; // 3*(M+1)
        int offTau_ = 0;  // M

        // --- layout state ---
        bool corridor_q_active_ = true; // true when we use [P0 | q's | PM | v̂ | â | τ]

        // pre-allocate for vectors
        // scratch reused by every f/g eval
        mutable std::vector<Vec3> P_s_, V_s_, A_s_;
        mutable std::vector<std::array<Vec3, 6>> CP_s_;
        mutable std::vector<double> T_s_;
        mutable std::vector<Vec3> gP_s_, gV_s_, gA_s_;
        mutable std::vector<double> gT_s_;

        // in class SolverLBFGS (private:)
        bool scale_derivs_ = true; // default = current behavior (v̂ = T̄ V, â = T̄² A)
    };

} // namespace lbfgs
