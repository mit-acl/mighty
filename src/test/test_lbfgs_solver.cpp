/*
command to run: cbd && ./build/dynus/test_lbfgs_solver && cat results.txt
and in test_unconstrained_opt_3d.ipynb, you can visualize the results saved in results.txt.
*/

// test_lbfgs_solver.cpp
#include <mighty/lbfgs_solver.hpp>
#include <mighty/lbfgs.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <future>
#include <vector>
#include <chrono>

#include <decomp_util/ellipsoid_decomp.h>
#include <decomp_util/seed_decomp.h>

/// @brief Build a toy static‐constraint set of 3 segments, each with one plane
static std::vector<LinearConstraint3D> makeSampleStaticConstraints()
{
    std::vector<LinearConstraint3D> constraints;
    constraints.reserve(3);

    // Segment 0: A0 = [ -0.5, 0.0, 0.1 ], b0 = [ -0.5 ]
    {
        Eigen::Matrix<double, 1, 3> A0;
        A0 << -0.5, 0.0, 0.1;
        Eigen::VectorXd b0(1);
        b0 << -0.5;
        constraints.emplace_back(A0, b0);
    }

    // Segment 1: A1 = [ -0.5, 0.5, 0.1 ], b1 = [ 1.0 ]
    {
        Eigen::Matrix<double, 1, 3> A1;
        A1 << -0.5, 0.3, 0.01;
        Eigen::VectorXd b1(1);
        b1 << -0.1;
        constraints.emplace_back(A1, b1);
    }

    // Segment 2: A2 = [ 0.0, 0.5, 0.1 ], b2 = [ 2.0 ]
    {
        Eigen::Matrix<double, 1, 3> A2;
        A2 << 0.0, 0.5, 0.1;
        Eigen::VectorXd b2(1);
        b2 << 2.0;
        constraints.emplace_back(A2, b2);
    }

    return constraints;
}

/// @brief Build a topy dynamic obstacles
static std::vector<std::shared_ptr<dynTraj>> makeSampleDynamicObstacles()
{
    std::vector<std::shared_ptr<dynTraj>> obstacles;
    obstacles.reserve(2);

    // Obstacle 0
    {
        dynTraj d(
            Eigen::Vector3d{1.0, 3.0, 1.0},  // x0
            Eigen::Vector3d{-0.3, 0.3, 0.0}, // v0
            Eigen::Vector3d{0.0, 0.0, 0.0},  // a0
            Eigen::Vector3d{5.0, 5.0, 1.0},  // xf
            Eigen::Vector3d{0.0, 0.0, 0.0},  // vf
            Eigen::Vector3d{0.0, 0.0, 0.0},  // af
            0.0,                             // poly_start_time
            8.0                              // poly_end_time
        );
        d.id = 10;
        d.is_agent = false;

        // Create a shared pointer to the dynTraj object
        std::shared_ptr<dynTraj> d_ptr = std::make_shared<dynTraj>(std::move(d));
        // Add the shared pointer to the obstacles vector
        obstacles.push_back(d_ptr);
    }

    // Obstacle 1
    {
        dynTraj d(
            Eigen::Vector3d{0.0, 0.5, 0.0},
            Eigen::Vector3d{0.5, 0.0, 0.0},
            Eigen::Vector3d{0.0, 0.0, 0.0},
            Eigen::Vector3d{6.0, 5.0, -0.5},
            Eigen::Vector3d{0.0, 0.0, 0.0},
            Eigen::Vector3d{0.0, 0.0, 0.0},
            0.0, // poly_start_time
            8.0  // poly_end_time
        );
        d.id = 11;
        d.is_agent = false;

        // Create a shared pointer to the dynTraj object
        std::shared_ptr<dynTraj> d_ptr = std::make_shared<dynTraj>(std::move(d));
        // Add the shared pointer to the obstacles vector
        obstacles.push_back(d_ptr);
    }

    return obstacles;
}

int main()
{
    using namespace lbfgs;

    // 0) Set up the problem data

    // Constant parameters
    planner_params_t planner_params;
    planner_params.verbose = false;      // enable verbose output
    planner_params.V_max = 2.0;          // max velocity
    planner_params.A_max = 10.0;         // max acceleration
    planner_params.J_max = 100.0;        // max jerk
    planner_params.num_perturbation = 8; // number of perturbations for initial guesses
    planner_params.r_max = 1.0;          // perturbation radius for initial guesses
    planner_params.dyn_weight = 0.0;
    planner_params.time_weight = 1.0;
    planner_params.pos_anchor_weight = 1.0;
    planner_params.jerk_weight = 1.0; // weight for jerk in the objective function
    planner_params.stat_weight = 1.0;
    planner_params.dyn_constr_vel_weight = 1.0;
    planner_params.dyn_constr_acc_weight = 1.0;
    planner_params.dyn_constr_jerk_weight = 1.0;
    planner_params.dyn_constr_bodyrate_weight = 1.0; // weight for body rate constraints
    planner_params.dyn_constr_tilt_weight = 1.0;     // weight for tilt constraints
    planner_params.dyn_constr_thrust_weight = 1.0;   // weight for thrust constraints
    planner_params.Co = 0.2; // for static obstacle avoidance
    planner_params.Cw = 1.0; // for dynamic obstacle avoidance
    planner_params.BIG = 1e8;
    planner_params.dc = 0.01;                // descretiation constant
    planner_params.integral_resolution = 30; // resolution for the integral in the optimization
    planner_params.hinge_mu = 1e-2;          // hinge mu for the optimization
    planner_params.omega_max = 0.5;          // max body rate in rad
    planner_params.tilt_max_rad = 0.174533;  // max tilt in radians (e.g., 35° in rad)
    planner_params.f_min = 0.0;              // min thrust in N
    planner_params.f_max = 10.0;             // max thrust in N
    planner_params.mass = 1.0;               // mass in kg
    planner_params.g = 9.81;                 // gravity in m/s^2

    // Global waypoints
    vec_Vecf<3> global_wps = {
        Vecf<3>{0.0, 0.0, 0.0}, // Start point
        Vecf<3>{0.0, 2.0, 0.0}, // Intermediate point
        Vecf<3>{0.0, 4.0, 0.0}, // Another intermediate point
        Vecf<3>{5.0, 5.0, 0.0}  // End point
    };
    // Static constraints for 3 segments, each with one plane
    std::vector<LinearConstraint3D> safe_corridor = makeSampleStaticConstraints();
    // Dynamic obstacles
    std::vector<std::shared_ptr<dynTraj>> obstacles = makeSampleDynamicObstacles();
    // Initial and final conditions
    state initial_state, final_state;
    initial_state.setPos(0.0, 0.0, 0.0);
    initial_state.setVel(0.0, 0.0, 0.0);
    initial_state.setAccel(0.0, 0.0, 0.0);
    final_state.setPos(5.0, 5.0, 0.0);
    final_state.setVel(0.0, 0.0, 0.0);
    final_state.setAccel(0.0, 0.0, 0.0);

    // 1) Construct your trajectory‐optimization solver (holds problem data)
    std::shared_ptr<SolverLBFGS> prototype_solver_ptr = std::make_shared<SolverLBFGS>();
    prototype_solver_ptr->initializeSolver(planner_params);
    double t0 = 0.0; // initial time
    double initial_guess_computation_time;
    prototype_solver_ptr->prepareSolverForReplan(t0, global_wps, safe_corridor, obstacles, initial_state, final_state, initial_guess_computation_time, false);

    // CHANGE HERE USING POYTHON RESULTS
    static const double z0_data[39] = {
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 0, -0.050662662931596956, 4.87257505985749, 0, -0.094808548455064326,
        -0.51372833982035027, 0, 0, 4, 0, 1.2695909971048445, 20.878459542727704, 0,
        32.007205070390974, -109.61237615412993, 0, 5, 5, 0, 0, 0,
        0, 0, 0, 0, 0.73205080756887719, 0.34183291757628842, 3.3740550707847063};

    double py_cost = 38962.797717491572;
    Eigen::VectorXd py_gradient(prototype_solver_ptr->getNumDecisionVariables());
    py_gradient << 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, -216.314747808777, 23427.219173763511, -7.028109920136, -102.470196836984, 10833.165241272422, 2.849014573552, -7.143723000717, 1194.308615690353, -0.740841338578, 244.820220206900, -24685.640050017631, -1.293869191514, -31.575452670201, 3039.270650175118, 1.202832900445, 2.406450934129, -128.335797118355, 0.126886298017, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, -25746.114496695474, -80891.995004942903, -34939.713362082628;

    // Implementation
    Eigen::Map<const Eigen::VectorXd> z0(z0_data, 39);
    // define list_z0 for single initial guess
    std::vector<Eigen::VectorXd> list_z0 = {z0};
    size_t N = list_z0.size();

    // 3) Prepare storage for results
    std::vector<Eigen::VectorXd> list_zopt(N);
    std::vector<double> list_fopt(N);
    std::vector<std::future<void>> futures;
    futures.reserve(N);

    // 4) Set up L-BFGS parameters
    lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = static_cast<int>(list_z0[0].size());
    lbfgs_params.f_dec_coeff = 1e-1; // allow larger Armijo steps
    lbfgs_params.past = 3;
    lbfgs_params.max_linesearch = 64; // fewer backtracking tries
    lbfgs_params.max_iterations = 500;  // allow more iterations
    // lbfgs_params.cautious_factor = 1e-4;
    lbfgs_params.g_epsilon = 1e-5;       // stop when ‖g‖ isn’t super small
    lbfgs_params.delta = 1e-5;      // stop once f-improvement is minimal
    bool optimize_multiple = false; // set to false for single optimization

    // Test 1: test cost & gradient evaluation
    std::cout << "Testing cost and gradient evaluation...\n";
    Eigen::VectorXd g;
    double f = prototype_solver_ptr->evaluateObjectiveAndGradient(list_z0[0], g);

    // This is the value from Python, so we can compare:
    // Compare the difference between C++ and Python results
    if (std::abs(f - py_cost) > 1e-6)
    {
        std::cout << "\033[31m"; // red text
    }
    else
    {
        std::cout << "\033[32m"; // green text
    }
    std::cout << "Cost difference: " << std::abs(f - py_cost) << "\033[0m\n"; // reset color

    // Compare the norm of the gradient
    double grad_norm = g.norm();
    double py_grad_norm = py_gradient.norm();
    if (std::abs(grad_norm - py_grad_norm) > 1e-6)
    {
        std::cout << "\033[31m"; // red text
        std::cout << "c++ grad: " << g.transpose() << "\n";
        std::cout << "py grad: " << py_gradient.transpose() << "\n";
    }
    else
    {
        std::cout << "\033[32m"; // green text
    }
    std::cout << "Gradient norm difference: " << std::abs(grad_norm - py_grad_norm) << "\033[0m\n"; // reset color}

    // 5) Launch one async optimization per initial guess
    auto t_start = std::chrono::high_resolution_clock::now();

    // Just single optimization z0 = list_z0[0];
    if (!optimize_multiple)
    {
        Eigen::VectorXd zopt;
        double fopt;

        int status = prototype_solver_ptr->optimize(list_z0[0], zopt, fopt, lbfgs_params);
        list_zopt[0] = zopt;
        list_fopt[0] = fopt;
    }
    else
    {
        // Parallelization approach (is there any faster way to do this?)
        for (size_t i = 0; i < N; ++i)
        {
            futures.emplace_back(std::async(std::launch::async,
                                            [&, i]()
                                            {
                                                // make a fresh solver for thread‐safety
                                                std::shared_ptr<SolverLBFGS> solver_ptr = std::make_shared<SolverLBFGS>();
                                                solver_ptr->initializeSolver(planner_params);
                                                solver_ptr->prepareSolverForReplan(0.0, global_wps, safe_corridor, obstacles, initial_state, final_state, initial_guess_computation_time); // initial time t0 = 0.0

                                                // copy initial guess
                                                Eigen::VectorXd z0 = list_z0[i];
                                                Eigen::VectorXd zopt;
                                                double fopt;

                                                // run optimization
                                                int status = solver_ptr->optimize(z0, zopt, fopt, lbfgs_params);
                                                list_zopt[i] = zopt;
                                                list_fopt[i] = fopt;
                                            }));
        }
        // wait for all to finish
        for (auto &f : futures)
            f.get();
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    // 6) ROS-related debugging
    if (!optimize_multiple)
    {
        std::cout << "ROS-related debugging:\n";

        // First recover the final control points and times from z_opt
        prototype_solver_ptr->reconstructPVATCPopt(list_zopt[0]);

        // Initialize the size of the goal setpoints

        std::vector<state> goal_setpoints;
        PieceWisePol pwp_to_share;
        std::vector<Eigen::Matrix<double, 3, 6>> cps;
        prototype_solver_ptr->getGoalSetpoints(goal_setpoints);
        // prototype_solver_ptr->getPieceWisePol(pwp_to_share);
        prototype_solver_ptr->getControlPoints(cps); // Bezier control points

        // Print velocity norm history
        // std::cout << "Vel Norm:\n";
        // for (const auto &sp : goal_setpoints)
        // {
        //     std::cout << "Time: " << sp.t << ", Vel Norm: " << sp.vel.norm() << "\n";
        // }

        // // Print the goal setpoints
        // std::cout << "Number of goal setpoints: " << goal_setpoints.size() << "\n";
        // std::cout << "Goal setpoints:\n";
        // for (const auto &sp : goal_setpoints)
        // {
        //     std::cout << "Time: " << sp.t << ", Pos: " << sp.pos.transpose()
        //               << ", Vel: " << sp.vel.transpose()
        //               << ", Accel: " << sp.accel.transpose()
        //               << ", Jerk: " << sp.jerk.transpose()
        //               << ", Yaw: " << sp.yaw
        //               << ", DYaw: " << sp.dyaw
        //               << "\n";
        // }

        // // Print the piecewise polynomial
        // std::cout << "Piecewise polynomial:\n";
        // pwp_to_share.print();

        // // Print the control points
        // std::cout << "Control points:\n";
        // for (const auto &cp : cps)
        // {
        //     std::cout << cp.transpose() << "\n";
        // }
    }

    // 7) Dump results to a Python‐readable file
    std::ofstream out("results.txt");
    out << "# Auto-generated by test_lbfgs_solver.cpp\n";
    out << "# total optimization time: " << elapsed_ms.count() << " ms\n\n";

    // z0_list
    out << "z0_list = [\n";
    for (size_t i = 0; i < N; ++i)
    {
        out << "    np.array([";
        for (int j = 0; j < list_z0[i].size(); ++j)
        {
            out << std::fixed << std::setprecision(8)
                << list_z0[i][j];
            if (j + 1 < list_z0[i].size())
                out << ", ";
        }
        out << "], dtype=float)";
        if (i + 1 < N)
            out << ",";
        out << "\n";
    }
    out << "]\n\n";

    // zopt_list
    out << "zopt_list = [\n";
    for (size_t i = 0; i < N; ++i)
    {
        out << "    np.array([";
        for (int j = 0; j < list_zopt[i].size(); ++j)
        {
            out << std::fixed << std::setprecision(8)
                << list_zopt[i][j];
            if (j + 1 < list_zopt[i].size())
                out << ", ";
        }
        out << "], dtype=float)";
        if (i + 1 < N)
            out << ",";
        out << "\n";
    }
    out << "]\n\n";

    // Optionally save fopt_list too
    out << "fopt_list = [";
    for (size_t i = 0; i < N; ++i)
    {
        out << std::fixed << std::setprecision(8)
            << list_fopt[i];
        if (i + 1 < N)
            out << ", ";
    }
    out << "]\n";

    out.close();

    std::cout << "Done. Wrote " << N << " results to results.txt"
              << " (total " << elapsed_ms.count() << " ms).\n";

    return 0;
}
