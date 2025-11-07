/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "rclcpp/rclcpp.hpp"
#include "mighty/mighty_type.hpp"
#include "mighty/mighty.hpp"
#include "dgp/utils.hpp"
#include <dynus_interfaces/msg/dyn_traj.hpp>
#include <dynus_interfaces/msg/dyn_traj_array.hpp>
#include "dynus_interfaces/msg/state.hpp"
#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/yaw_output.hpp"
#include "dynus_interfaces/msg/pn_adaptation.hpp"
#include <mighty/utils.hpp>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include <tf2_eigen/tf2_eigen.h>
#include <Eigen/Dense>

#include <visualization_msgs/msg/marker_array.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <std_msgs/msg/empty.hpp>
#include <std_msgs/msg/string.hpp>

#include <pcl_conversions/pcl_conversions.h>
#include <fstream>                                         // for flie operations
#include "ament_index_cpp/get_package_share_directory.hpp" // for getting the package path
#include <vector>
#include <unordered_map>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations" // pcl::SAC_SAMPLE_SIZE is protected since PCL 1.8.0 // NOLINT
#include <pcl/sample_consensus/model_types.h>
#pragma GCC diagnostic pop
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/passthrough.h>
#include "visualization_msgs/msg/marker_array.hpp"
#include "nav_msgs/msg/occupancy_grid.hpp"
#include "std_msgs/msg/color_rgba.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"
#include "std_srvs/srv/empty.hpp"
#include "pcl_conversions/pcl_conversions.h"
#include "pcl_ros/transforms.hpp"
#include "tf2_ros/buffer.h"
#include "tf2_ros/create_timer_ros.h"
#include "tf2_ros/message_filter.h"
#include "tf2_ros/transform_listener.h"
#include "message_filters/subscriber.h"
#include <ament_index_cpp/get_package_share_directory.hpp>
#include <algorithm>
#include <memory>
#include <string>
#include <execution>

// prefix
using namespace std::chrono_literals;

// Define the synchronization policy
typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::msg::PointCloud2, sensor_msgs::msg::PointCloud2> MySyncPolicy;
typedef message_filters::Synchronizer<MySyncPolicy> Sync;

// ROS2 component (run-time composition)
// common practice is to use namespace for components
namespace mighty
{

    using PCLPoint = pcl::PointXYZ;
    using PCLPointCloud = pcl::PointCloud<pcl::PointXYZ>;

    class MIGHTY_NODE : public rclcpp::Node
    {

    public:
        MIGHTY_NODE();  
        ~MIGHTY_NODE(); 

    private:

        // Callbacks
        void replanCallback();
        void trajCallback(const dynus_interfaces::msg::DynTraj::SharedPtr msg);
        void stateCallback(const dynus_interfaces::msg::State::SharedPtr msg);
        void terminalGoalCallback(const geometry_msgs::msg::PoseStamped &msg);
        void mapCallback(const sensor_msgs::msg::PointCloud2::ConstPtr &pcl2ptr_map_ros, const sensor_msgs::msg::PointCloud2::ConstPtr &pcl2ptr_unk_ros);
        void goalReachedCheckCallback();
        void convertDynTrajMsg2DynTraj(const dynus_interfaces::msg::DynTraj &msg, std::shared_ptr<dynTraj> &traj, double current_time);
        void cleanUpOldTrajsCallback();
        void getInitialPoseHwCallback();

        // Others
        void declareParameters();                                                                       
        void setParameters();                                                                           
        void printParameters();                                                                         
        void createMarkerArrayFromVec_Vec3f(const vec_Vec3f &occupied_cells, const std_msgs::msg::ColorRGBA &color, int namespace_id, double scale, visualization_msgs::msg::MarkerArray *marker_array);
        void clearMarkerArray(visualization_msgs::msg::MarkerArray &path_marker, rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr publisher); 
        void clearMarkerActualTraj();                                                                                                                           
        void runSim();                                                                                                                                          
        void printComputationTime(bool result);                                                                                                              
        void recordData(bool result);                                                                                                                          
        void logData();                                                                                                                         
        void setComputationTimesToZero();
        void constructFOVMarker();
        void retrieveData();
        void retrieveDataForVisualizationAndTopics();

        // Functions to publish
        void publishGlobalPath();
        void publishFreeGlobalPath();
        void publishPoly();
        void publishTraj();
        void publishOwnTraj();
        void publishActualTraj();
        void publishGoal();         // Publish the goal (trajectory points)
        void publishPointG() const; // Publish the point G (projected terminal goal)
        void publishPointE() const; // Publish the point E (Local trajectory's goal)
        void publishPointA() const; // Publish the point A (starting point)
        void publishCurrentState(const state &state) const;
        void publishState(const state &data, const rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr &publisher) const;
        void publishFOV();
        void publisCps();
        void publishStaticPushPoints();
        void publishLocalGlobalPath();
        void publishVelocityInText(const Eigen::Vector3d &position, double velocity);

        // Timers for callback
        rclcpp::TimerBase::SharedPtr timer_replanning_;
        rclcpp::TimerBase::SharedPtr timer_goal_;
        rclcpp::TimerBase::SharedPtr timer_update_tmap_;
        rclcpp::TimerBase::SharedPtr timer_goal_reached_check_;
        rclcpp::TimerBase::SharedPtr timer_cleanup_old_trajs_;
        rclcpp::TimerBase::SharedPtr timer_trajs_update_for_tmap_;
        rclcpp::TimerBase::SharedPtr timer_check_plan_safety_;
        rclcpp::TimerBase::SharedPtr timer_initial_pose_;

        // Callback groups
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_1_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_2_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_3_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_4_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_5_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_6_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_7_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_8_;
        rclcpp::CallbackGroup::SharedPtr cb_group_mu_9_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_1_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_2_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_3_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_4_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_5_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_6_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_7_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_8_;
        rclcpp::CallbackGroup::SharedPtr cb_group_re_9_;
        rclcpp::CallbackGroup::SharedPtr cb_group_map_;
        rclcpp::CallbackGroup::SharedPtr cb_group_replan_;
        rclcpp::CallbackGroup::SharedPtr cb_group_goal_;

        // Publishers
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_dgp_path_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_original_dgp_path_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_free_dgp_path_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_local_global_path_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_local_global_path_after_push_marker_;
        rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_dynamic_map_;
        rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_free_map_;
        rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr pub_unknown_map_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_static_map_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_dynamic_map_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_free_map_marker_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_unknown_map_marker_;
        rclcpp::Publisher<decomp_ros_msgs::msg::PolyhedronArray>::SharedPtr pub_poly_whole_;
        rclcpp::Publisher<decomp_ros_msgs::msg::PolyhedronArray>::SharedPtr pub_poly_safe_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_traj_committed_colored_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_traj_subopt_colored_;
        rclcpp::Publisher<dynus_interfaces::msg::DynTraj>::SharedPtr pub_own_traj_;
        rclcpp::Publisher<dynus_interfaces::msg::Goal>::SharedPtr pub_goal_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_G_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_E_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_G_term_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_point_A_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_current_state_;
        rclcpp::Publisher<std_msgs::msg::Empty>::SharedPtr pub_goal_reached_;
        rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_setpoint_;
        rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_actual_traj_;
        rclcpp::Publisher<dynus_interfaces::msg::YawOutput>::SharedPtr pub_yaw_output_;
        rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_fov_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_cp_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_static_push_points_;
        rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_p_points_;
        rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_vel_text_;

        // Subscribers
        rclcpp::Subscription<dynus_interfaces::msg::DynTraj>::SharedPtr sub_traj_;
        rclcpp::Subscription<dynus_interfaces::msg::DynTraj>::SharedPtr sub_predicted_traj_;
        rclcpp::Subscription<dynus_interfaces::msg::State>::SharedPtr sub_state_;
        rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr sub_terminal_goal_;

        // Time synchronizer
        message_filters::Subscriber<sensor_msgs::msg::PointCloud2> occup_grid_sub_;
        message_filters::Subscriber<sensor_msgs::msg::PointCloud2> unknown_grid_sub_;
        std::shared_ptr<Sync> sync_;

        // Visualization
        visualization_msgs::msg::MarkerArray dgp_path_marker_;
        visualization_msgs::msg::MarkerArray original_dgp_path_marker_;
        visualization_msgs::msg::MarkerArray dgp_free_path_marker_;
        visualization_msgs::msg::MarkerArray dgp_local_global_path_marker_;
        visualization_msgs::msg::MarkerArray dgp_local_global_path_after_push_marker_;
        visualization_msgs::msg::MarkerArray traj_committed_colored_;
        visualization_msgs::msg::MarkerArray traj_subopt_colored_;
        visualization_msgs::msg::Marker marker_fov_;

        // Mutex
        std::mutex cloud_callback_mutex_; // Mutex for cloud callback

        // Parameters
        int id_;
        std::string ns_;
        std::string id_str_;
        parameters par_;
        double final_g_ = 0.0;                  // only for debugging
        bool verbose_computation_time_ = false; // only for debugging
        int marker_fov_id_ = 0;

        // DYNUS pointer
        std::shared_ptr<MIGHTY> mighty_ptr_;

        // Global Path Benchmarking
        std::string file_path_;        // only for benchmarking
        std::vector<std::tuple<bool,   // Result
                               double, // Cost
                               double, // Total replanning time
                               double, // Global planning time
                               double, // CVX decomposition time
                               double, // Initial guess computation time
                               double, // Local trajectory time
                               double, // Safe paths time
                               double, // Safety Check time
                               double, // Yaw sequence time
                               double, // Yaw fitting time
                               double, // Static JPS time in DGP
                               double, // Check path time in DGP
                               double, // Dynamic A* time in DGP
                               double  // Recover path time in DGP
                               >>
            global_path_benchmark_; // only for benchmarking

        // debug
        double replan_last_time_called_ = 0.0;

        // Record computation time
        double global_planning_time_ = 0.0;
        double dgp_static_jps_time_ = 0.0;
        double dgp_check_path_time_ = 0.0;
        double dgp_dynamic_astar_time_ = 0.0;
        double dgp_recover_path_time_ = 0.0;
        double cvx_decomp_time_ = 0.0;
        double initial_guess_computation_time_ = 0.0; // Time for computing initial guess
        double local_traj_computation_time_ = 0.0;
        double safe_paths_time_ = 0.0;
        double safety_check_time_ = 0.0;
        double yaw_sequence_time_ = 0.0;
        double yaw_fitting_time_ = 0.0;
        double replanning_computation_time_ = 0.0;
        double current_time_for_debug_ = 0.0;

        // Visualization
        vec_E<Polyhedron<3>> poly_whole_;
        vec_E<Polyhedron<3>> poly_safe_;
        std::vector<state> goal_setpoints_;
        std::vector<std::vector<state>> list_subopt_goal_setpoints_;
        int actual_traj_id_ = 0;

        // Yaw Optimization Debugging
        std::vector<double> optimal_yaw_sequence_;
        std::vector<double> yaw_control_points_;
        std::vector<double> yaw_knots_;

        // Local trajectory debugging
        std::vector<Eigen::Matrix<double, 3, 6>> cps_;

        // Static push points and p points
        vec_Vecf<3> static_push_points_;
        vec_Vecf<3> p_points_;
        int static_push_points_id_ = 0;
        int p_points_id_ = 0;

        // Trajectory sharing
        // PieceWisePol pwp_to_share_; // Piecewise polynomial

        // Flags
        bool state_initialized_ = false;             // State initialized
        bool replan_timer_started_ = false;          // Replan timer started
        bool publish_actual_traj_called_ = false;    // Publish actual trajectory called
        bool use_benchmark_ = false;                 // Use benchmark

        // Actual trajectory variables for ground robots
        Eigen::Vector3d actual_traj_prev_pos_;
        double actual_traj_prev_time_;

        int last_subopt_count_{0}; // tracks how many subopt strips we drew last time

        // D435 parameters
        std::string d435_depth_frame_id_;
        std::string lidar_frame_id_;
        std::string d435_camera_info_topic_;

        // TF2 buffer and listener
        std::shared_ptr<tf2_ros::Buffer> tf2_buffer_;
        std::shared_ptr<tf2_ros::TransformListener> tf2_listener_;

        // initial pose (for hardware)
        bool initial_pose_received_ = false;
        std::string initial_pose_topic_;
        geometry_msgs::msg::TransformStamped init_pose_transform_stamped_;

        // Timer to make sure we don't sample point cloud too often
        rclcpp::Time last_lidar_callback_time_;
        rclcpp::Time last_depth_camera_callback_time_;

    };

} // namespace mighty