/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

 #include "misc/visualizer.hpp"
 #include "gcopter/flatness.hpp"
 #include "gcopter/sfc_gen.hpp"
 #include "gcopter/voxel_map.hpp"
 #include "gcopter/gcopter.hpp"
 #include "gcopter/trajectory.hpp"
 
 #include <mighty/mighty_solver.hpp>
 
 #include <rclcpp/rclcpp.hpp>
 #include <geometry_msgs/msg/pose_stamped.hpp>
 #include <sensor_msgs/msg/point_cloud2.hpp>
 #include <std_msgs/msg/float64.hpp>
 
 #include <Eigen/Dense>
 
 #include <algorithm>
 #include <array>
 #include <chrono>
 #include <cmath>
 #include <iostream>
 #include <memory>
 #include <numeric>
 #include <string>
 #include <vector>
 
 // Short-hands for types used by MIGHTY
 using lbfgs::planner_params_t;
 using lbfgs::SolverMIGHTY;
 using lbfgs::Vec3;
 
 using PolyhedronV = Eigen::Matrix3Xd;        // vertices of one convex polyhedron (V-rep)
 using PolyhedraV = std::vector<PolyhedronV>; // list of polyhedra (SFC V-rep)
 
 // ============================================================================
 // Helpers (math / conversions)
 // ============================================================================
 
 /**
  * @brief Evaluate a 5th-degree Bezier segment and its first three derivatives.
  *
  * Control points: CP[0..5] in R^3.
  * Parameterization: s = u in [0,1], physical time t scaled by segment duration T.
  *
  * @param CP  Six 3D control points (P0..P5).
  * @param T   Segment duration (seconds).
  * @param u   Normalized parameter in [0,1].
  * @param p   Out: position at u.
  * @param v   Out: velocity at u (dt units).
  * @param a   Out: acceleration at u (dt^2 units).
  * @param j   Out: jerk at u (dt^3 units).
  */
 static inline void evalBezier5_PVAJ(const std::array<Eigen::Vector3d, 6> &CP,
                                     double T,
                                     double u,
                                     Eigen::Vector3d &p,
                                     Eigen::Vector3d &v,
                                     Eigen::Vector3d &a,
                                     Eigen::Vector3d &j)
 {
     const double om = 1.0 - u;
 
     // Position basis (degree-5)
     const double b0 = std::pow(om, 5);
     const double b1 = 5.0 * u * std::pow(om, 4);
     const double b2 = 10.0 * u * u * std::pow(om, 3);
     const double b3 = 10.0 * std::pow(u, 3) * om * om;
     const double b4 = 5.0 * std::pow(u, 4) * om;
     const double b5 = std::pow(u, 5);
     p = b0 * CP[0] + b1 * CP[1] + b2 * CP[2] + b3 * CP[3] + b4 * CP[4] + b5 * CP[5];
 
     // Forward differences for derivative bases
     std::array<Eigen::Vector3d, 5> D1;
     for (int i = 0; i < 5; ++i) D1[i] = CP[i + 1] - CP[i];
     std::array<Eigen::Vector3d, 4> D2;
     for (int i = 0; i < 4; ++i) D2[i] = D1[i + 1] - D1[i];
     std::array<Eigen::Vector3d, 3> D3;
     for (int i = 0; i < 3; ++i) D3[i] = D2[i + 1] - D2[i];
 
     const double invT  = 1.0 / std::max(T, 1.0e-9);
     const double invT2 = invT * invT;
     const double invT3 = invT2 * invT;
 
     // Velocity (degree-4 basis)
     const double c0 = std::pow(om, 4);
     const double c1 = 4.0 * u * std::pow(om, 3);
     const double c2 = 6.0 * u * u * om * om;
     const double c3 = 4.0 * std::pow(u, 3) * om;
     const double c4 = std::pow(u, 4);
     v = 5.0 * invT * (c0 * D1[0] + c1 * D1[1] + c2 * D1[2] + c3 * D1[3] + c4 * D1[4]);
 
     // Acceleration (degree-3 basis)
     const double d0 = std::pow(om, 3);
     const double d1 = 3.0 * u * om * om;
     const double d2 = 3.0 * u * u * om;
     const double d3 = std::pow(u, 3);
     a = 20.0 * invT2 * (d0 * D2[0] + d1 * D2[1] + d2 * D2[2] + d3 * D2[3]);
 
     // Jerk (degree-2 basis)
     const double e0 = om * om;
     const double e1 = 2.0 * u * om;
     const double e2 = u * u;
     j = 60.0 * invT3 * (e0 * D3[0] + e1 * D3[1] + e2 * D3[2]);
 }
 
 /**
  * @brief Convert GCOPTER half-space representation (A x + d ≤ 0) to
  *        the MIGHTY format (−A x ≤ d).
  */
 static std::vector<LinearConstraint3D>
 toLinearConstraints(const std::vector<Eigen::MatrixX4d> &hPolys)
 {
     std::vector<LinearConstraint3D> out;
     out.reserve(hPolys.size());
 
     for (const auto &H : hPolys)
     {
         if (H.rows() == 0 || H.cols() != 4)
         {
             out.emplace_back(); // keep an empty placeholder
             continue;
         }
 
         const int m = static_cast<int>(H.rows());
         Eigen::MatrixXd A(m, 3);
         Eigen::VectorXd b(m);
         A = -H.leftCols<3>(); // -A
         b = H.col(3);         // d
         out.emplace_back(A, b);
     }
     return out;
 }
 
 /**
  * @brief Minimal wrapper around the MIGHTY solver to produce {CP, T}.
  *
  * This function prepares the L-BFGS parameters and invokes the solver with
  * the provided initial guess (xi, segment times), corridor constraints,
  * and physical parameters. Returns control points and times on success.
  */
 struct MightyOut
 {
     std::vector<std::array<Eigen::Vector3d, 6>> CP; // per segment control points
     std::vector<double> T;                          // per segment durations
     double obj{0.0};                                // final objective value
     double wall_ms{0.0};                            // wall-clock time (ms)
 };
 
 static MightyOut runMighty(const vec_Vecf<3> &global_wps,
                            const PolyhedraV &vPolytopes,
                            const Eigen::VectorXd &initial_xi,
                            const Eigen::VectorXd &init_times,
                            const std::vector<LinearConstraint3D> &safe_corridor,
                            const state &initial_state,
                            const state &final_state,
                            const planner_params_t &params,
                            const Eigen::VectorXd &physical_params,
                            int mem_size,
                            int max_linesearch,
                            int past,
                            double min_step,
                            int max_iterations,
                            double g_epsilon)
 {
     MightyOut out;
 
     auto solver = std::make_shared<SolverMIGHTY>();
     solver->initializeSolver(params, physical_params);
 
     double ig_ms = 0.0;                               // initial-guess timing (unused here)
     const double t0 = 0.0;
 
     // Prepare for (re)planning from t0
     solver->prepareSolverForReplan(t0,
                                    global_wps,
                                    safe_corridor,
                                    initial_xi,
                                    init_times,
                                    initial_state,
                                    final_state,
                                    ig_ms,
                                    vPolytopes,
                                    /*use_multiple_initial_guesses*/ false);
 
     // Use the first initial guess (we only asked for one)
     const auto &z0_list = solver->getInitialGuesses();
     if (z0_list.empty())
         return out;
 
     Eigen::VectorXd z0 = z0_list.front(), zopt;
     double fopt = 0.0;
 
     // L-BFGS parameters (tuned externally via Config)
     lbfgs::lbfgs_parameter_t lb{};
     lb.mem_size = mem_size;
     lb.max_linesearch = max_linesearch;
     lb.past = past;
     lb.min_step = min_step;
     lb.max_iterations = max_iterations;
     lb.g_epsilon = g_epsilon;
     lb.delta = 0.0; // relative tol handled internally in gradient code
 
     // Optimize
     const auto t_start = std::chrono::high_resolution_clock::now();
     solver->optimize(z0, zopt, fopt, lb);
     const auto t_end = std::chrono::high_resolution_clock::now();
 
     out.wall_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() * 1e-3;
     out.obj = fopt;
 
     // Recover Bezier control points and segment times
     {
         std::vector<Vec3> P, V, A;
         std::vector<std::array<Vec3, 6>> CPv;
         std::vector<double> Tv;
         solver->reconstruct(zopt, P, V, A, CPv, Tv);
 
         out.CP.resize(CPv.size());
         for (size_t s = 0; s < CPv.size(); ++s)
             for (int j = 0; j < 6; ++j)
                 out.CP[s][j] = Eigen::Vector3d(CPv[s][j].x(), CPv[s][j].y(), CPv[s][j].z());
         out.T = Tv;
     }
     return out;
 }
 
 // ============================================================================
 // ROS 2 Node Configuration
 // ============================================================================
 
 /**
  * @brief Runtime configuration (loaded via ROS 2 parameters with defaults).
  *
  * All values here are read only once at startup (node constructor).
  * Tuning can be done via parameter YAMLs/launch files.
  */
 struct Config
 {
     // ---- Topics / map pre-processing ----
     std::string mapTopic;     ///< PointCloud2 input for occupancy
     std::string targetTopic;  ///< PoseStamped goal selection
     double dilateRadius;      ///< Obstacle dilation (m)
     double voxelWidth;        ///< Voxel size (m)
     std::vector<double> mapBound; ///< [xmin, xmax, ymin, ymax, zmin, zmax]
     double timeoutRRT;        ///< RRT timeout (s)
 
     // ---- Physical / dynamic limits ----
     double maxVelMag;     ///< Max speed (m/s)
     double maxBdrMag;     ///< Max body rate (rad/s)
     double maxTiltAngle;  ///< Max tilt (rad)
     double minThrust;     ///< Min thrust (N)   (for ring)
     double maxThrust;     ///< Max thrust (N)   (for ring)
     double vehicleMass;   ///< Mass (kg)
     double gravAcc;       ///< Gravity (m/s^2)
     double horizDrag;     ///< Aerodynamic horizontal drag
     double vertDrag;      ///< Aerodynamic vertical drag
     double parasDrag;     ///< Parasitic drag
     double speedEps;      ///< Small speed epsilon
 
     // ---- GCOPTER (for initial guess only) ----
     double weightT;            ///< Time weight
     std::vector<double> chiVec;///< Penalties for constraints
     double smoothingEps;       ///< Hinge smoothing
     int    integralIntervs;    ///< Sampling resolution
     double relCostTol;         ///< Relative cost tolerance
 
     // ---- MIGHTY weights (objective & soft limits) ----
     double mightyJerkWeight;   ///< Jerk penalty
     double mightyWeightT;      ///< Time term weight (w_T)
     double mightyPosWeight;    ///< Corridor clearance weight
     double mightyVelWeight;    ///< Velocity limit weight
     double mightyOmegaWeight;  ///< Body-rate limit weight
     double mightyThetaWeight;  ///< Tilt limit weight
     double mightyThrustWeight; ///< Thrust ring weight
     double mightyInitTurnBF;   ///< Initial turn buffer (degrees)
 
     // ---- Safe Flight Corridor (SFC) ----
     double sfc_progress; ///< Sliding window progress (m)
     double sfc_range;    ///< SFC inflation range (m)
 
     // ---- L-BFGS tuning ----
     int    lbfgs_mem_size;
     int    lbfgs_max_linesearch;
     int    lbfgs_past;
     double lbfgs_min_step;
     int    lbfgs_max_iterations;
     double lbfgs_g_epsilon;
 
     explicit Config(rclcpp::Node &node)
     {
         // Declare parameters with sensible defaults
         node.declare_parameter("MapTopic", std::string("/voxel_map"));
         node.declare_parameter("TargetTopic", std::string("/goal_pose"));
         node.declare_parameter("DilateRadius", 0.5);
         node.declare_parameter("VoxelWidth", 0.25);
         node.declare_parameter("MapBound", std::vector<double>{-25.0, 25.0, -25.0, 25.0, 0.0, 5.0});
         node.declare_parameter("TimeoutRRT", 0.02);
 
         node.declare_parameter("MaxVelMag", 4.0);
         node.declare_parameter("MaxBdrMag", 2.1);
         node.declare_parameter("MaxTiltAngle", 1.05);
         node.declare_parameter("MinThrust", 2.0);
         node.declare_parameter("MaxThrust", 12.0);
         node.declare_parameter("VehicleMass", 0.61);
         node.declare_parameter("GravAcc", 9.8);
         node.declare_parameter("HorizDrag", 0.70);
         node.declare_parameter("VertDrag", 0.80);
         node.declare_parameter("ParasDrag", 0.01);
         node.declare_parameter("SpeedEps", 1.0e-4);
 
         node.declare_parameter("WeightT", 20.0);
         node.declare_parameter("ChiVec", std::vector<double>{1.0e4, 1.0e4, 1.0e4, 1.0e4, 1.0e5});
         node.declare_parameter("SmoothingEps", 1.0e-2);
         node.declare_parameter("IntegralIntervs", 16);
         node.declare_parameter("RelCostTol", 1.0e-5);
 
         // MIGHTY
         node.declare_parameter("MIGHTYJerkWeight", 0.1);
         node.declare_parameter("MIGHTYWeightT", 10.0);
         node.declare_parameter("MIGHTYPosWeight", 1.0e4);
         node.declare_parameter("MIGHTYVelWeight", 1.0e4);
         node.declare_parameter("MIGHTYOmegaWeight", 1.0e4);
         node.declare_parameter("MIGHTYThetaWeight", 1.0e4);
         node.declare_parameter("MIGHTYThrustWeight", 1.0e5);
         node.declare_parameter("MIGHTYInitTurnBF", 15.0);
 
         // SFC
         node.declare_parameter("SFC_PROGRESS", 3.0);
         node.declare_parameter("SFC_RANGE", 7.0);
 
         // L-BFGS
         node.declare_parameter("LBFGS_Mem", 50);
         node.declare_parameter("LBFGS_MaxLineSearch", 20);
         node.declare_parameter("LBFGS_Past", 3);
         node.declare_parameter("LBFGS_MinStep", 1.0e-32);
         node.declare_parameter("LBFGS_MaxIters", 1000);
         node.declare_parameter("LBFGS_gEps", 1.0e-5);
 
         // Fetch parameters
         node.get_parameter("MapTopic", mapTopic);
         node.get_parameter("TargetTopic", targetTopic);
         node.get_parameter("DilateRadius", dilateRadius);
         node.get_parameter("VoxelWidth", voxelWidth);
         node.get_parameter("MapBound", mapBound);
         node.get_parameter("TimeoutRRT", timeoutRRT);
 
         node.get_parameter("MaxVelMag", maxVelMag);
         node.get_parameter("MaxBdrMag", maxBdrMag);
         node.get_parameter("MaxTiltAngle", maxTiltAngle);
         node.get_parameter("MinThrust", minThrust);
         node.get_parameter("MaxThrust", maxThrust);
         node.get_parameter("VehicleMass", vehicleMass);
         node.get_parameter("GravAcc", gravAcc);
         node.get_parameter("HorizDrag", horizDrag);
         node.get_parameter("VertDrag", vertDrag);
         node.get_parameter("ParasDrag", parasDrag);
         node.get_parameter("SpeedEps", speedEps);
 
         node.get_parameter("WeightT", weightT);
         node.get_parameter("ChiVec", chiVec);
         node.get_parameter("SmoothingEps", smoothingEps);
         node.get_parameter("IntegralIntervs", integralIntervs);
         node.get_parameter("RelCostTol", relCostTol);
 
         node.get_parameter("MIGHTYJerkWeight", mightyJerkWeight);
         node.get_parameter("MIGHTYWeightT", mightyWeightT);
         node.get_parameter("MIGHTYPosWeight", mightyPosWeight);
         node.get_parameter("MIGHTYVelWeight", mightyVelWeight);
         node.get_parameter("MIGHTYOmegaWeight", mightyOmegaWeight);
         node.get_parameter("MIGHTYThetaWeight", mightyThetaWeight);
         node.get_parameter("MIGHTYThrustWeight", mightyThrustWeight);
         node.get_parameter("MIGHTYInitTurnBF", mightyInitTurnBF);
 
         node.get_parameter("SFC_PROGRESS", sfc_progress);
         node.get_parameter("SFC_RANGE", sfc_range);
 
         node.get_parameter("LBFGS_Mem", lbfgs_mem_size);
         node.get_parameter("LBFGS_MaxLineSearch", lbfgs_max_linesearch);
         node.get_parameter("LBFGS_Past", lbfgs_past);
         node.get_parameter("LBFGS_MinStep", lbfgs_min_step);
         node.get_parameter("LBFGS_MaxIters", lbfgs_max_iterations);
         node.get_parameter("LBFGS_gEps", lbfgs_g_epsilon);
     }
 };
 
 // ============================================================================
 // ROS 2 Node (GlobalPlanner)
 // ============================================================================
 
 /**
  * @brief ROS2 node that builds an SFC from a voxel map and optimizes a
  *        dynamically-feasible trajectory with MIGHTY.
  */
 class GlobalPlanner : public rclcpp::Node
 {
 private:
     Config config;
 
     // Subscriptions
     rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr mapSub;
     rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr targetSub;
 
     // Visualization/output & occupancy map
     Visualizer visualizer;
     voxel_map::VoxelMap voxelMap;
 
     bool mapInitialized{false};
     std::vector<Eigen::Vector3d> startGoal; // [start, goal] once both are selected
 
     // MIGHTY result cache for animation/telemetry
     std::vector<std::array<Eigen::Vector3d, 6>> mightyCP_;
     std::vector<double> mightyT_;
     std::vector<double> mightyEdges_; // cumulative times for segment lookups
     double mightyStamp_{0.0};
 
     [[maybe_unused]] static inline double clamp(double x, double lo, double hi)
     {
         return std::max(lo, std::min(x, hi));
     }
 
 public:
     GlobalPlanner()
         : Node("mighty_node")
         , config(*this)
         , visualizer(*this)
     {
         // Build voxel map from bounds and resolution
         const auto &mb = config.mapBound;
         Eigen::Vector3i xyz(
             static_cast<int>((mb[1] - mb[0]) / config.voxelWidth),
             static_cast<int>((mb[3] - mb[2]) / config.voxelWidth),
             static_cast<int>((mb[5] - mb[4]) / config.voxelWidth));
         Eigen::Vector3d offset(mb[0], mb[2], mb[4]);
         voxelMap = voxel_map::VoxelMap(xyz, offset, config.voxelWidth);
 
         // Subscriptions
         mapSub = this->create_subscription<sensor_msgs::msg::PointCloud2>(
             config.mapTopic, 1,
             std::bind(&GlobalPlanner::mapCallBack, this, std::placeholders::_1));
 
         targetSub = this->create_subscription<geometry_msgs::msg::PoseStamped>(
             config.targetTopic, 1,
             std::bind(&GlobalPlanner::targetCallBack, this, std::placeholders::_1));
 
         RCLCPP_INFO(this->get_logger(), "MIGHTY global planner ready.");
     }
 
     /**
      * @brief One-shot voxel ingestion and dilation.
      *        The first received cloud is used; later clouds are ignored.
      */
     void mapCallBack(const sensor_msgs::msg::PointCloud2::SharedPtr msg)
     {
         if (mapInitialized) return;
 
         const size_t total  = msg->data.size() / msg->point_step;
         const float *fdata  = reinterpret_cast<const float *>(&msg->data[0]);
         const size_t stride = msg->point_step / sizeof(float);
 
         for (size_t i = 0; i < total; ++i)
         {
             const size_t cur = stride * i;
             const float x = fdata[cur + 0];
             const float y = fdata[cur + 1];
             const float z = fdata[cur + 2];
 
             if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z))
                 voxelMap.setOccupied(Eigen::Vector3d(x, y, z));
         }
 
         // Inflate obstacles by the given radius (in voxel units)
         voxelMap.dilate(std::ceil(config.dilateRadius / voxelMap.getScale()));
         mapInitialized = true;
 
         RCLCPP_INFO(this->get_logger(), "Voxel map ingested & dilated.");
     }
 
     /**
      * @brief Collects two clicks (start, goal) from PoseStamped messages.
      *        z is mapped into the [zmin+radius, zmax-radius] band using orientation.z as a slider.
      */
     void targetCallBack(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
     {
         if (!mapInitialized) return;
 
         if (startGoal.size() >= 2) startGoal.clear();
 
         // Map orientation.z ∈ [-1,1] to a feasible altitude band
         const double zGoal = config.mapBound[4] + config.dilateRadius +
                              std::fabs(msg->pose.orientation.z) *
                                  (config.mapBound[5] - config.mapBound[4] - 2.0 * config.dilateRadius);
 
         const Eigen::Vector3d goal(msg->pose.position.x,
                                    msg->pose.position.y,
                                    zGoal);
 
         if (voxelMap.query(goal) == 0)
         {
             visualizer.visualizeStartGoal(goal, 0.5, startGoal.size());
             startGoal.emplace_back(goal);
         }
         else
         {
             RCLCPP_WARN(this->get_logger(), "Infeasible Position Selected!");
         }
 
         plan();
     }
 
     /**
      * @brief Full pipeline: RRT → SFC → GCOPTER initial guess → MIGHTY → visualize.
      */
     void plan()
     {
         if (startGoal.size() != 2) return;
 
         // ---------------------------------------------------------------------
         // 1) Global route via RRT (coarse, collision-free waypoints)
         // ---------------------------------------------------------------------
         std::vector<Eigen::Vector3d> route;
         sfc_gen::planPath<voxel_map::VoxelMap>(
             startGoal[0], startGoal[1],
             voxelMap.getOrigin(), voxelMap.getCorner(),
             &voxelMap, config.timeoutRRT, route);
 
         if (route.size() <= 1)
         {
             RCLCPP_WARN(this->get_logger(), "Global route failed.");
             return;
         }
 
         // ---------------------------------------------------------------------
         // 2) Safe Flight Corridor (SFC) from route & surface points
         // ---------------------------------------------------------------------
         std::vector<Eigen::MatrixX4d> hPolys; // H-rep (Ax + d ≤ 0)
         std::vector<Eigen::Vector3d> pc;      // surface points
         voxelMap.getSurf(pc);
 
         sfc_gen::convexCover(route, pc,
                              voxelMap.getOrigin(), voxelMap.getCorner(),
                              config.sfc_progress, config.sfc_range, hPolys);
         sfc_gen::shortCut(hPolys);  // remove redundancy / tighten
         visualizer.visualizePolytope(hPolys);
 
         // ---------------------------------------------------------------------
         // 3) Boundary states for GCOPTER (position + zero v/a)
         // ---------------------------------------------------------------------
         Eigen::Matrix3d iniState, finState;
         iniState << route.front(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero();
         finState << route.back(),  Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero();
 
         // ---------------------------------------------------------------------
         // 4) GCOPTER setup (for initial guess + V-polytopes only)
         // ---------------------------------------------------------------------
         gcopter::GCOPTER_PolytopeSFC gc;
 
         Eigen::VectorXd magnitudeBounds(5), penaltyWeights(5), physicalParams(6);
         magnitudeBounds << config.maxVelMag, config.maxBdrMag, config.maxTiltAngle,
                             config.minThrust, config.maxThrust;
         penaltyWeights << config.chiVec[0], config.chiVec[1], config.chiVec[2],
                            config.chiVec[3], config.chiVec[4];
         physicalParams << config.vehicleMass, config.gravAcc, config.horizDrag,
                            config.vertDrag, config.parasDrag, config.speedEps;
 
         Trajectory<5> dummy_traj; // internal pipeline exercise (not visualized)
         dummy_traj.clear();
 
         if (!gc.setup(config.weightT, iniState, finState, hPolys,
                       INFINITY, config.smoothingEps, config.integralIntervs,
                       magnitudeBounds, penaltyWeights, physicalParams))
         {
             RCLCPP_ERROR(this->get_logger(), "GCOPTER setup failed.");
             return;
         }
         if (std::isinf(gc.optimize(dummy_traj, config.relCostTol)))
         {
             RCLCPP_ERROR(this->get_logger(), "GCOPTER failed to produce initial guess.");
             return;
         }
 
         // Initial guess points/times and seam variables
         Eigen::Matrix3Xd init_points;
         Eigen::VectorXd init_times, initial_xi;
         gc.getInitialGuess(init_points, init_times, initial_xi);
 
         // Corridor as V-polytopes (for MIGHTY seams/constraints)
         PolyhedraV vPolys;
         gc.getVPolytopes(vPolys);
 
         // ---------------------------------------------------------------------
         // 5) Build MIGHTY parameter pack
         // ---------------------------------------------------------------------
         planner_params_t mighty{};
         mighty.V_max                         = config.maxVelMag;
         mighty.time_weight                   = config.mightyWeightT;
         mighty.jerk_weight                   = config.mightyJerkWeight;
         mighty.stat_weight                   = config.mightyPosWeight;
         mighty.dyn_constr_vel_weight         = config.mightyVelWeight;
         mighty.dyn_constr_bodyrate_weight    = config.mightyOmegaWeight;
         mighty.dyn_constr_tilt_weight        = config.mightyThetaWeight;
         mighty.dyn_constr_thrust_weight      = config.mightyThrustWeight;
         mighty.init_turn_bf                  = config.mightyInitTurnBF; // degrees
         mighty.Co                            = 0.0;     // corridor clearance target
         mighty.BIG                           = 1e9;     // big-M for constraints
         mighty.dc                            = 0.01;    // sampling step for outputs
         mighty.integral_resolution           = config.integralIntervs;
         mighty.hinge_mu                      = config.smoothingEps;
         mighty.omega_max                     = config.maxBdrMag;
         mighty.tilt_max_rad                  = config.maxTiltAngle;
         mighty.f_min                         = config.minThrust;
         mighty.f_max                         = config.maxThrust;
         mighty.mass                          = config.vehicleMass;
         mighty.g                             = config.gravAcc;
 
         // ---------------------------------------------------------------------
         // 6) Endpoint states for MIGHTY (pos with zero vel/acc)
         // ---------------------------------------------------------------------
         state init_state, goal_state;
         init_state.pos   = Vec3(route.front().x(), route.front().y(), route.front().z());
         init_state.vel   = Vec3::Zero();
         init_state.accel = Vec3::Zero();
         goal_state.pos   = Vec3(route.back().x(),  route.back().y(),  route.back().z());
         goal_state.vel   = Vec3::Zero();
         goal_state.accel = Vec3::Zero();
 
         // ---------------------------------------------------------------------
         // 7) Waypoints for MIGHTY’s initial guess: start + (GCOPTER interior) + goal
         // ---------------------------------------------------------------------
         vec_Vecf<3> wps;
         wps.emplace_back(init_state.pos.x(), init_state.pos.y(), init_state.pos.z());
         for (int c = 0; c < init_points.cols(); ++c)
             wps.emplace_back(init_points(0, c), init_points(1, c), init_points(2, c));
         wps.emplace_back(goal_state.pos.x(), goal_state.pos.y(), goal_state.pos.z());
 
         // ---------------------------------------------------------------------
         // 8) Static corridor planes (linear constraints per segment)
         // ---------------------------------------------------------------------
         const auto lin_constraints = toLinearConstraints(hPolys);
 
         // ---------------------------------------------------------------------
         // 9) Physical params for flatness map (mass/drag/gravity/speedEps)
         // ---------------------------------------------------------------------
         Eigen::VectorXd phys(6);
         phys << config.vehicleMass, config.gravAcc, config.horizDrag,
                 config.vertDrag,    config.parasDrag, config.speedEps;
 
         // ---------------------------------------------------------------------
         // 10) Run MIGHTY optimization (L-BFGS)
         // ---------------------------------------------------------------------
         const MightyOut out = runMighty(wps, vPolys, initial_xi, init_times, lin_constraints,
                                         init_state, goal_state, mighty, phys,
                                         config.lbfgs_mem_size, config.lbfgs_max_linesearch,
                                         config.lbfgs_past, config.lbfgs_min_step,
                                         config.lbfgs_max_iterations, config.lbfgs_g_epsilon);
 
         if (out.CP.empty() || out.T.empty())
         {
             RCLCPP_ERROR(this->get_logger(), "MIGHTY produced no trajectory.");
             return;
         }
 
         // ---------------------------------------------------------------------
         // 11) Visualize & cache for live telemetry/animation
         // ---------------------------------------------------------------------
         mightyCP_ = out.CP;
         mightyT_  = out.T;
 
         // Build cumulative edges for segment lookup
         mightyEdges_.assign(mightyT_.size() + 1, 0.0);
         for (size_t s = 0; s < mightyT_.size(); ++s)
             mightyEdges_[s + 1] = mightyEdges_[s] + mightyT_[s];
 
         // Stamp “now” for animation base time
         mightyStamp_ = this->now().seconds();
 
         // Draw the full Bezier path (green)
         visualizer.visualizeBezier(mightyCP_, mightyT_,
                                    /*ns=*/"mighty",
                                    /*width=*/0.3f,
                                    /*samples=*/120,
                                    /*r=*/0.0f, /*g=*/1.0f, /*b=*/0.0f, /*a=*/1.0f,
                                    /*frame_id=*/"odom",
                                    /*speed_scale=*/config.maxVelMag);
     }
 
     /**
      * @brief Called continuously in the main loop to publish dynamic telemetry
      *        (speed, thrust, tilt, body-rate) and draw a moving sphere.
      */
     void process()
     {
         if (mightyCP_.empty() || mightyT_.empty())
             return;
 
         const double Ttot = mightyEdges_.empty()
                                 ? std::accumulate(mightyT_.begin(), mightyT_.end(), 0.0)
                                 : mightyEdges_.back();
 
         const double t = this->now().seconds() - mightyStamp_;
         if (t < 0.0 || t > Ttot)
             return;
 
         // Locate the active segment for time t
         const int seg = std::clamp(
             static_cast<int>(std::upper_bound(mightyEdges_.begin(), mightyEdges_.end(), t) -
                              mightyEdges_.begin()) - 1,
             0,
             static_cast<int>(mightyCP_.size()) - 1);
 
         const double Ts = std::max(1.0e-9, mightyT_[seg]);
         const double u  = std::clamp((t - mightyEdges_[seg]) / Ts, 0.0, 1.0);
 
         // Evaluate P/V/A/J on current segment
         Eigen::Vector3d p, v, a, j;
         evalBezier5_PVAJ(mightyCP_[seg], Ts, u, p, v, a, j);
 
         // Flatness → thrust, attitude(quat), body-rates
         double            thr   = 0.0;
         Eigen::Vector4d   quat  = Eigen::Vector4d::Zero(); // [w, x, y, z]
         Eigen::Vector3d   omg   = Eigen::Vector3d::Zero();
 
         Eigen::VectorXd phys(6);
         phys << config.vehicleMass, config.gravAcc, config.horizDrag,
                 config.vertDrag,    config.parasDrag, config.speedEps;
 
         flatness::FlatnessMap flatmap;
         flatmap.reset(phys(0), phys(1), phys(2), phys(3), phys(4), phys(5));
         flatmap.forward(v, a, j, /*yaw*/ 0.0, /*yawdot*/ 0.0, thr, quat, omg);
 
         const double speed        = v.norm();
         const double bodyratemag  = omg.norm();
         const double tiltangle    = std::acos(std::max(-1.0, std::min(1.0,
                                    1.0 - 2.0 * (quat(1) * quat(1) + quat(2) * quat(2)))));
 
         // Publish telemetry
         std_msgs::msg::Float64 speedMsg, thrMsg, tiltMsg, bdrMsg;
         speedMsg.data = speed;
         thrMsg.data   = thr;
         tiltMsg.data  = tiltangle;
         bdrMsg.data   = bodyratemag;
 
         visualizer.speedPub->publish(speedMsg);
         visualizer.thrPub->publish(thrMsg);
         visualizer.tiltPub->publish(tiltMsg);
         visualizer.bdrPub->publish(bdrMsg);
 
         // Moving dot at current position
         visualizer.visualizeSphere(p, config.dilateRadius);
     }
 };
 
 // ============================================================================
 // main()
 // ============================================================================
 
 int main(int argc, char **argv)
 {
     rclcpp::init(argc, argv);
     auto node = std::make_shared<GlobalPlanner>();
 
     // High-rate spin to keep telemetry smooth
     rclcpp::Rate rate(1000);
     while (rclcpp::ok())
     {
         node->process();
         rclcpp::spin_some(node);
         rate.sleep();
     }
     rclcpp::shutdown();
     return 0;
 }
 