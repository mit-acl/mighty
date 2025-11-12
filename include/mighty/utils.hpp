/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef MIGHTY_UTILS_HPP
#define MIGHTY_UTILS_HPP

#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/color_rgba.hpp"
#include <geometry_msgs/msg/vector3.hpp>
#include <geometry_msgs/msg/point.hpp>

#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <mighty/mighty_type.hpp>
#include <dgp/utils.hpp>
#include <deque>

#include "dynus_interfaces/msg/pwp_traj.hpp"
#include "dynus_interfaces/msg/coeff_poly3.hpp"
#include "dynus_interfaces/msg/dyn_traj.hpp"

namespace mighty_utils
{

    // Define colors
    static constexpr int red_normal = 1;
    static constexpr int red_trans = 2;
    static constexpr int red_trans_trans = 3;
    static constexpr int green_normal = 4;
    static constexpr int blue_normal = 5;
    static constexpr int blue_trans = 6;
    static constexpr int blue_trans_trans = 7;
    static constexpr int blue_light = 8;
    static constexpr int yellow_normal = 9;
    static constexpr int orange_trans = 10;
    static constexpr int black_trans = 11;
    static constexpr int teal_normal = 12;
    static constexpr int green_trans_trans = 13;

    // Convert a PieceWisePol Message to a PieceWisePol
    PieceWisePol convertPwpMsg2Pwp(const dynus_interfaces::msg::PWPTraj& pwp_msg);

    // Convert a PieceWisePol to a PieceWisePol Message
    dynus_interfaces::msg::PWPTraj convertPwp2PwpMsg(const PieceWisePol& pwp);

    // Convert a PieceWisePol to a Colored Marker Array
    visualization_msgs::msg::MarkerArray convertPwp2ColoredMarkerArray(PieceWisePol& pwp, int samples);

    // Convert a PieceWisePol to a Marker Array
    geometry_msgs::msg::Point convertEigen2Point(Eigen::Vector3d vector);

    // Function to convert std::vector<float> to Eigen::Vector3d
    Eigen::Vector3d convertCovMsg2Cov(const std::vector<float>& msg_cov);

    // Function to convert Eigen::Vector3d to std::vector<float>
    std::vector<float> convertCov2CovMsg(const Eigen::Vector3d& cov);

    // Function to convert std::vector<float> to Eigen::Matrix<double, 6, 1>
    Eigen::Matrix<double, 6, 1> convertCoeffMsg2Coeff(const std::vector<float>& msg_coeff);

    // Provide a color based on the id
    std_msgs::msg::ColorRGBA getColor(int id);

    // Convert coefficients to control points
    std::vector<Eigen::Matrix<double, 3, 4>> convertCoefficients2ControlPoints(const PieceWisePol& pwp, const Eigen::Matrix<double, 4, 4>& A_rest_pos_basis_inverse);

    // Get the time of the trajectory
    double getMinTimeDoubleIntegrator1D(const double p0, const double v0, const double pf, const double vf,
                                    const double v_max, const double a_max);
    double getMinTimeDoubleIntegrator3D(const Eigen::Vector3d& p0, const Eigen::Vector3d& v0, const Eigen::Vector3d& pf,
                                    const Eigen::Vector3d& vf, const Eigen::Vector3d& v_max,
                                    const Eigen::Vector3d& a_max);

    //## From Wikipedia - http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    void quaternion2Euler(tf2::Quaternion q, double& roll, double& pitch, double& yaw);
    void quaternion2Euler(Eigen::Quaterniond q, double& roll, double& pitch, double& yaw);
    void quaternion2Euler(geometry_msgs::msg::Quaternion q, double& roll, double& pitch, double& yaw);

    // saturate function
    void saturate(double& var, double min, double max);

    // angle wrap function
    void angle_wrap(double& diff);

    // project a point to a box
    Eigen::Vector3d projectPointToBox(Eigen::Vector3d& P1, Eigen::Vector3d& P2, double wdx, double wdy, double wdz);

    // project a point to a sphere
    Eigen::Vector3d projectPointToSphere(const Eigen::Vector3d &P1, const Eigen::Vector3d &P2, double radius);

    // convert visualization_msgs::msg::MarkerArray to vec_Vecf<3>
    void convertMarkerArray2Vec_Vec_Vecf3(const visualization_msgs::msg::MarkerArray& marker_array, std::vector<vec_Vecf<3>>& vec, std::vector<double>& scale_vec);

    // create more vertexes in path
    void createMoreVertexes(vec_Vecf<3>& path, double d);

    // euclidean distance
    double euclideanDistance(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

    // get 4x4 transformation matrix from geometry_msgs::msg::TransformStamped
    Eigen::Matrix4d transformStampedToMatrix(const geometry_msgs::msg::TransformStamped& transform_stamped);

    // estimate velocities of the path
    void findVelocitiesInPath(const vec_Vecf<3>& path, vec_Vecf<3>& velocities, const state& A, const Eigen::Vector3d& v_max_3d, bool verbose);

    // get travel times of the path
    std::vector<double> getTravelTimes(const vec_Vecf<3>& path, const state& A, bool debug_verbose, const Eigen::Vector3d& v_max_3d, const Eigen::Vector3d& a_max_3d);

    // identity geometry_msgs::msg::Pose
    geometry_msgs::msg::Pose identityGeometryMsgsPose();

    // sign function
    template <typename T>
    int sgn(T val)
    {
    return (T(0) < val) - (val < T(0));
    }

}  // namespace mighty_utils

#endif
