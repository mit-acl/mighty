/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef DGP_UTILS_HPP
#define DGP_UTILS_HPP
#include <iostream>
#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/color_rgba.hpp>
#include <geometry_msgs/msg/vector3.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <dgp/data_utils.hpp>
#include <dgp/termcolor.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include "mighty/mighty_type.hpp"
#include <deque>

#define RED 1
#define RED_TRANS 2
#define RED_TRANS_TRANS 3
#define GREEN 4
#define BLUE 5
#define BLUE_TRANS 6
#define BLUE_TRANS_TRANS 7
#define BLUE_LIGHT 8
#define YELLOW 9
#define ORANGE_TRANS 10
#define BLACK_TRANS 11
#define ORANGE 12
#define GREEN_TRANS_TRANS 13

#define STATE 0
#define INPUT 1

#define WHOLE_TRAJ 0
#define RESCUE_PATH 1

#define OCCUPIED_SPACE 1
#define UNKOWN_AND_OCCUPIED_SPACE 2

void printStateDeque(std::deque<state> &data);

void printStateVector(std::vector<state> &data);

void vectorOfVectors2MarkerArray(vec_Vecf<3> traj, visualization_msgs::msg::MarkerArray *m_array, std_msgs::msg::ColorRGBA color,
                                 int type = visualization_msgs::msg::Marker::ARROW,
                                 std::vector<double> radii = std::vector<double>());

std_msgs::msg::ColorRGBA getColorJet(double v, double vmin, double vmax);

std_msgs::msg::ColorRGBA color(int id);

// ## From Wikipedia - http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
void quaternion2Euler(tf2::Quaternion q, double &roll, double &pitch, double &yaw);

void quaternion2Euler(Eigen::Quaterniond q, double &roll, double &pitch, double &yaw);

void quaternion2Euler(geometry_msgs::msg::Quaternion q, double &roll, double &pitch, double &yaw);

void saturate(double &var, double min, double max);

visualization_msgs::msg::Marker getMarkerSphere(double scale, int my_color);

double angleBetVectors(const Eigen::Vector3d &a, const Eigen::Vector3d &b);

// returns the points around B sampled in the sphere with radius r and center center.
std::vector<Eigen::Vector3d> samplePointsSphere(Eigen::Vector3d &B, double r, Eigen::Vector3d &center);

void printElementsOfJPS(vec_Vecf<3> &path);

// returns the points around B sampled in the sphere with radius r and center center, and sampled intelligently with
// the given path
// last_index_inside_sphere is the the index of the last point that is inside the sphere (should be provided as a
// parameter to this function)
// B is the first intersection of JPS with the sphere
std::vector<Eigen::Vector3d> samplePointsSphereWithJPS(Eigen::Vector3d &B, double r, Eigen::Vector3d &center_sent,
                                                       vec_Vecf<3> &path_sent, int last_index_inside_sphere);

void angle_wrap(double &diff);

pcl::PointXYZ eigenPoint2pclPoint(Eigen::Vector3d &p);

vec_Vec3f pclptr_to_vec(const pcl::KdTreeFLANN<pcl::PointXYZ>::PointCloudConstPtr ptr_cloud);

vec_Vec3f pclptr_to_vec(const pcl::KdTreeFLANN<pcl::PointXYZ>::PointCloudConstPtr ptr_cloud1,
                        const pcl::KdTreeFLANN<pcl::PointXYZ>::PointCloudConstPtr ptr_cloud2);

float solvePolyOrder2(Eigen::Vector3f coeff);

// coeff is from highest degree to lowest degree. Returns the smallest positive real solution. Returns -1 if a
// root is imaginary or if it's negative

geometry_msgs::msg::Point pointOrigin();

Eigen::Vector3d vec2eigen(geometry_msgs::msg::Vector3 vector);

geometry_msgs::msg::Vector3 eigen2rosvector(Eigen::Vector3d vector);

geometry_msgs::msg::Point eigen2point(Eigen::Vector3d vector);

geometry_msgs::msg::Vector3 vectorNull();

geometry_msgs::msg::Vector3 vectorUniform(double a);

template <typename T>
using vec_E = std::vector<T, Eigen::aligned_allocator<T>>;

template <int N>
using Vecf = Eigen::Matrix<decimal_t, N, 1>; // Be CAREFUL, because this is with doubles!

template <int N>
using vec_Vecf = vec_E<Vecf<N>>;

// returns 1 if there is an intersection between the segment P1-P2 and the plane given by coeff=[A B C D]
// (Ax+By+Cz+D==0)  returns 0 if there is no intersection.
// The intersection point is saved in "intersection"
bool getIntersectionWithPlane(const Eigen::Vector3d &P1, const Eigen::Vector3d &P2, const Eigen::Vector4d &coeff,
                              Eigen::Vector3d &intersection);

double normJPS(vec_Vecf<3> &path, int index_start);

// Crop the end of a JPS path by a given distance
void reduceJPSbyDistance(vec_Vecf<3> &path, double d);

// given 2 points (A inside and B outside the sphere) it computes the intersection of the lines between
// that 2 points and the sphere
Eigen::Vector3d getIntersectionWithSphere(Eigen::Vector3d &A, Eigen::Vector3d &B, double r, Eigen::Vector3d &center);

// Given a path (starting inside the sphere and finishing outside of it) expressed by a vector of 3D-vectors (points),
// it returns its first intersection with a sphere of radius=r and center=center
// the center is added as the first point of the path to ensure that the first element of the path is inside the sphere
// (to avoid issues with the first point of JPS2)
Eigen::Vector3d getFirstIntersectionWithSphere(vec_Vecf<3> &path, double r, Eigen::Vector3d &center,
                                               int *last_index_inside_sphere = NULL,
                                               bool *noPointsOutsideSphere = NULL);

// Given a path (starting inside the sphere and finishing outside of it) expressed by a vector of 3D-vectors (points),
// it returns its first intersection with a sphere of radius=r and center=center
Eigen::Vector3d getLastIntersectionWithSphere(vec_Vecf<3> path, double r, Eigen::Vector3d center);

double getDistancePath(vec_Vecf<3> &path);

// Same as the previous one, but also returns dist = the distance form the last intersection to the goal (following
// the path)
Eigen::Vector3d getLastIntersectionWithSphere(vec_Vecf<3> path, double r, Eigen::Vector3d center, double *Jdist);

// returns the point placed between two concentric spheres with radii ra, rb, and center=center
// If the path goes out from the 1st sphere, and then enters again, these points are also considered!
// I.e: it returns all the points between the first point that goes out from the sphere and the last point that is
// inside Sb
vec_Vecf<3> getPointsBw2Spheres(vec_Vecf<3> path, double ra, double rb, Eigen::Vector3d center);

vec_Vecf<3> copyJPS(vec_Vecf<3> path);

// Overload to be able to print a std::vector
template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
{
  if (!v.empty())
  {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

visualization_msgs::msg::MarkerArray stateVector2ColoredMarkerArray(const std::vector<state> &data, int type,
                                                                    double max_value,
                                                                    const rclcpp::Time &stamp);

void deleteVertexes(vec_Vecf<3> &JPS_path, int max_value);

#endif
