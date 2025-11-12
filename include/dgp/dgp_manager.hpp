/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

// ROS2
#include <rclcpp/rclcpp.hpp>

// Convex Decomposition includes
#include <decomp_rviz_plugins/data_ros_utils.hpp>
#include <decomp_util/ellipsoid_decomp.h>
#include <decomp_util/seed_decomp.h>

// Map includes
#include <dgp/read_map.hpp>

// DGP includes
#include <dgp/data_utils.hpp>
#include <dgp/dgp_planner.hpp>
#include <dgp/utils.hpp>

// Other includes
#include <Eigen/Dense>
#include <mutex>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include "timer.hpp"

// prefix
using namespace std::chrono_literals;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

class DGPManager
{
public:

    DGPManager();

    void setParameters(const parameters &par);
    void setupDGPPlanner(const std::string &global_planner, bool global_planner_verbose, double res, double v_max, double a_max, double j_max, int timeout_duration_ms, double w_unknown, double w_align, double decay_len_cells, double w_side, int los_cells = 3, double min_len = 0.5, double min_turn = 10.0);
    void getOccupiedCells(vec_Vecf<3> &occupied_cells);
    void getFreeCells(vec_Vecf<3> &free_cells);
    void getOccupiedCellsForCvxDecomp(vec_Vecf<3> &occupied_cells, const vec_Vecf<3> &path, bool use_for_safe_path);
    void getDynamicOccupiedCellsForVis(vec_Vecf<3> &occupied_cells, vec_Vecf<3> &free_cells, vec_Vecf<3> &unknown_cells, double current_time);
    void updateMap(double wdx, double wdy, double wdz, const Vec3f &center_map, const pcl::PointCloud<pcl::PointXYZ>::ConstPtr& pclptr);
    void freeStart(Vec3f &start_sent, double factor);
    void freeGoal(Vec3f &goal_sent, double factor);
    bool checkIfPointOccupied(const Vec3f &point);
    bool solveDGP(const Vec3f &start_sent, const Vec3f &start_vel, const Vec3f &goal_sent, double &final_g, double weight, double current_time, vec_Vecf<3> &path);
    bool checkIfPathInFree(const vec_Vecf<3> &path, vec_Vecf<3> &free_path);
    void getComputationTime(double &global_planning_time, double &dgp_static_jps_time, double &dgp_check_path_time, double &dgp_dynamic_astar_time, double &dgp_recover_path_time);
    bool cvxEllipsoidDecomp(const state &A, const vec_Vecf<3> &path, std::vector<LinearConstraint3D> &l_constraints, vec_E<Polyhedron<3>> &poly_out, bool use_for_safe_path = false);
    bool checkIfPointFree(const Vec3f &point) const;
    void updateReadMapUtil();
    void pushPathIntoFreeSpace(const vec_Vecf<3> &path, vec_Vecf<3> &free_path);
    void findClosestFreePoint(const Vec3f &point, Vec3f &closest_free_point);
    int countUnknownCells() const;
    int getTotalNumCells() const;
    void updateVmax(double v_max);
    void cleanUpPath(vec_Vecf<3> &path);
    bool isMapInitialized() const;
    bool checkIfPointHasNonFreeNeighbour(const Vec3f& point) const;
    void getVecOccupied(vec_Vec3f& vec_o);
    void updateVecOccupied(const vec_Vec3f& vec_o);
    void getVecUnknownOccupied(vec_Vec3f& vec_uo);
    void updateVecUnknownOccupied(const vec_Vec3f& vec_uo);
    void insertVecOccupiedToVecUnknownOccupied();

    std::shared_ptr<mighty::VoxelMapUtil> map_util_;
    std::shared_ptr<mighty::VoxelMapUtil> map_util_for_planning_;
    std::unique_ptr<DGPPlanner> planner_ptr_;
    
private:

    vec_Vec3f vec_o_;   // Vector that contains the occupied points
    vec_Vec3f vec_uo_;  // Vector that contains the unkown and occupied points

    // Mutex
    std::mutex mtx_map_util_;
    std::mutex mtx_vec_o_;  
    std::mutex mtx_vec_uo_;

    // Convex decomposition
    EllipsoidDecomp3D ellip_decomp_util_;
    std::vector<float> local_box_size_;

    // Parameters
    parameters par_;
    double weight_;
    double res_, drone_radius_;
    double v_max_, a_max_, j_max_;
    Eigen::Vector3d v_max_3d_, a_max_3d_, j_max_3d_;
    double max_dist_vertexes_ = 2.0;
    bool use_shrinked_box_ = false;
    double shrinked_box_size_ = 0.0;

    // Constants
    const Eigen::Vector3d unitX_ = Eigen::Vector3d::UnitX();
    const Eigen::Vector3d unitY_ = Eigen::Vector3d::UnitY();
    const Eigen::Vector3d unitZ_ = Eigen::Vector3d::UnitZ();

    // Flats
    bool map_initialized_ = false;
};
