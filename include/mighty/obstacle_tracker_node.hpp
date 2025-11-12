#ifndef OBSTACLE_TRACKER_NODE_HPP
#define OBSTACLE_TRACKER_NODE_HPP

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/filter.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/common/transforms.h>
#include <Eigen/Dense>
#include <visualization_msgs/msg/marker_array.hpp>
#include <dynus_interfaces/msg/dyn_traj.hpp>
#include <mighty/mighty_type.hpp>
#include <mighty/utils.hpp>
#include <tf2_eigen/tf2_eigen.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>


// EKF Parameters for 3D
struct EKFState {

    // Core variables
    Eigen::VectorXd x;  // state vector [x, y, z, vx, vy, vz, ax, ay, az]
    Eigen::MatrixXd P;  // state covariance matrix
    Eigen::MatrixXd Q;  // process noise covariance matrix
    Eigen::MatrixXd R;  // measurement noise covariance matrix

    // For visualization
    double time_updated = 0.0;
    Eigen::Vector3d bbox;
    int id;
    std_msgs::msg::ColorRGBA color;

    EKFState() {} // Constructor for Cluster struct
    EKFState(int state_size, Eigen::MatrixXd Q, Eigen::MatrixXd R, double time_updated, Eigen::Vector3d bbox, int id) {
        x = Eigen::VectorXd::Zero(state_size);
        P = Eigen::MatrixXd::Identity(state_size, state_size);
        this->Q = Q;
        this->R = R;
        this->time_updated = time_updated;
        this->bbox = bbox;
        this->id = id;
        setColor();
    }

    void setColor() {
        this->color.r = static_cast<float>(rand()) / RAND_MAX;  // Random red
        this->color.g = static_cast<float>(rand()) / RAND_MAX;  // Random green
        this->color.b = static_cast<float>(rand()) / RAND_MAX;  // Random blue
        // If you want to set a specific color
        // this->color.r = 0.0 / 255.0;
        // this->color.g = 0.0 / 255.0;
        // this->color.b = 255.0 / 255.0;
        this->color.a = 0.4;  // Opacity
    }
};

// Cluster struct
struct Cluster {
    EKFState ekf_state;
    Eigen::Vector3d centroid;
    Cluster() {} // Constructor
    void setEKFStateAndCentroid(EKFState ekf_state, Eigen::Vector3d centroid) 
    {
        this->ekf_state = ekf_state;
        this->centroid = centroid;
    }
};

class ObstacleTrackerNode : public rclcpp::Node {
public:
    ObstacleTrackerNode();

private:

    // Parameters
    int visual_level_;
    bool use_adaptive_kf_;
    double adaptive_kf_alpha_;
    double adaptive_kf_dt_;
    double cluster_tolerance_;
    int min_cluster_size_;
    int max_cluster_size_;
    double prediction_horizon_;
    double prediction_dt_;
    double time_to_delete_old_obstacles_;
    double cluster_bbox_cutoff_size_;
    bool use_life_time_for_box_visualization_;
    double box_visualization_duration_;
    double mighty_map_res_;
    double velocity_threshold_;
    double acceleration_threshold_;
    bool use_hardware_;

    // Subscriber and publisher
    rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr sub_pointcloud_;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_markers_;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_bboxes_;
    rclcpp::Publisher<dynus_interfaces::msg::DynTraj>::SharedPtr pub_predicted_traj_;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_unc_sphere_;

    // EKF states for multiple objects
    std::vector<EKFState> ekf_states_;  // Vector of EKF states for multiple objects

    // frame id
    std::string frame_id_ = "map";

    // id 
    int marker_id_ = 0;
    int ekf_state_id_ = 0;

    // TF2 buffer and listener
    std::shared_ptr<tf2_ros::Buffer> tf2_buffer_;
    std::shared_ptr<tf2_ros::TransformListener> tf2_listener_;

    // Prediction parameters
    double degree_for_pwp_ = 3;
    double degree_for_poly_ = 5;

    // Functions
    void declareAndsetParameters();
    void pointcloudCallback(const sensor_msgs::msg::PointCloud2::SharedPtr msg);
    void calculateAverageQandR(Eigen::MatrixXd &Q_avg, Eigen::MatrixXd &R_avg);
    void deleteOldEKFstates();
    void publishPredictions(const std::vector<Cluster> &clusters);
    void publishBoxes(const std::vector<Cluster>& clusters);
    void getCentroidsAndSizesOfClusters(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, const std::vector<pcl::PointIndices> &cluster_indices, std::vector<Eigen::Vector3d> &cluster_centroids, std::vector<Eigen::Vector3d> &cluster_sizes);
    Eigen::VectorXd polyfit(const std::vector<double>& t, const std::vector<double>& y, int degree);
    double calculateVariance(const std::vector<double>& t, const std::vector<double>& y, const Eigen::VectorXd& beta, int degree);
    void filterStaticObstacles();

};

#endif // OBSTACLE_TRACKER_NODE_HPP
