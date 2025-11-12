
/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef INITIAL_GUESS_HPP
#define INITIAL_GUESS_HPP

#include <vector>
#include <queue>
#include <Eigen/Dense>
#include <mighty/mighty_type.hpp>
#include <decomp_rviz_plugins/data_ros_utils.hpp>

class InitialGuessGenerator {
public:
    struct Node {
        Eigen::Vector3d position;
        int index;
        double cost;
        Node* parent;

        Node(const Eigen::Vector3d& pos, int idx, double c, Node* p)
            : position(pos), index(idx), cost(c), parent(p) {}

        bool operator<(const Node& other) const {
            return cost > other.cost; // Min-heap based on cost
        }
    };

    InitialGuessGenerator(int numControlPoints, double maxVelocity, double maxAcceleration);

    // Set up polytopes for trajectory generation
    void setPolytopes(const std::vector<LinearConstraint3D>& polytopes);

    // Set initial conditions (position, velocity, acceleration)
    void setInitialConditions(const Eigen::Vector3d& initialPos,
                              const Eigen::Vector3d& initialVel,
                              const Eigen::Vector3d& initialAcc);

    // Generate control points for the initial guess
    bool generateControlPoints(std::vector<Eigen::Vector3d>& controlPoints);

private:

    // Parameters

    // Variables
    Eigen::Vector3d q0_, q1_, q2_, q_final_m_2_, q_final_m_1_, q_final_;
    double closest_dist_so_far_ = std::numeric_limits<double>::max(); // stores the closest node found
    Node* closest_result_so_far_ptr_ = nullptr;                          // stores the closest node found          

    // Queues/Maps
    std::unordered_map<Eigen::Vector3i, bool, matrix_hash<Eigen::Vector3i>> map_open_list_; // Save the discrete points that are in the open list


    int numControlPoints_;
    double maxVelocity_;
    double maxAcceleration_;
    std::vector<LinearConstraint3D> polytopes_;

    Eigen::Vector3d initialPosition_;
    Eigen::Vector3d initialVelocity_;
    Eigen::Vector3d initialAcceleration_;

    // Helper function to check if a point is inside a polytope
    bool isInsidePolytope(const Eigen::Vector3d& point, const LinearConstraint3D& polytope);

    // Graph search utility
    void expandNode(const Node& currentNode, std::priority_queue<Node>& openSet);
};

#endif // INITIAL_GUESS_HPP
