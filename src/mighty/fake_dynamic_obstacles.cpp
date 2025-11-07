/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <rclcpp/rclcpp.hpp>
#include <random>
#include <cmath>
#include <vector>
#include <std_msgs/msg/color_rgba.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <geometry_msgs/msg/vector3.hpp>
#include <geometry_msgs/msg/pose.hpp>
#include <geometry_msgs/msg/quaternion.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include "dynus_interfaces/msg/dyn_traj.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include <gazebo_msgs/msg/model_state.hpp>
#include <std_msgs/msg/header.hpp>
#include <chrono>

using namespace std::chrono_literals;

class DynCorridor : public rclcpp::Node 
{
public:
    DynCorridor(int total_num_obs, bool gazebo)
        : Node("dynamic_obstacles"), total_num_obs(total_num_obs), gazebo(gazebo) 
    {


        // Set up timer to publish transforms periodically
        timer_ = this->create_wall_timer(std::chrono::milliseconds(10), std::bind(&DynCorridor::pubTF, this));
    }

    // Functions from Python code, adjusted to C++
    void getTrajectoryPosMeshBBox(double x, double y, double z, double scale, double offset, double slower) 
    {

        // Get string representation of positions
        wave_in_z(x_string, y_string, z_string, x, y, z, scale, offset, slower);

        // Get mesh path
        std::string mesh_path = "package://mighty/share/dynus/meshes/obstacles/model2.dae";

        // Get bounding box
        std::vector<double> bbox = {0.4, 0.4, 4.0};

    }

    void wave_in_z(std::string& x_string, std::string& y_string, std::string& z_string, double x, double y, double z, double scale, double offset, double slower) 
    {

        t_string="t/" + std::to_string(slower) + "+";
        x_string=std::to_string(x);
        y_string=std::to_string(y);
        z_string=std::to_string(scale)+"*(-sin( "+t_string +std::to_string(offset)+"))" + "+" + std::to_string(z);  

    }

    visualization_msgs::msg::Marker generateMarker(const std::string &mesh, const std::vector<double> &bbox, int i) {
        visualization_msgs::msg::Marker marker;
        marker.id = i;
        marker.ns = "mesh";
        marker.header.frame_id = "world";
        marker.type = visualization_msgs::msg::Marker::MESH_RESOURCE;
        marker.action = visualization_msgs::msg::Marker::ADD;
        marker.pose.position.x = 0.0;
        marker.pose.position.y = 0.0;
        marker.pose.position.z = 0.0;
        marker.pose.orientation.x = 0.0;
        marker.pose.orientation.y = 0.0;
        marker.pose.orientation.z = 0.0;
        marker.pose.orientation.w = 1.0;
        marker.lifetime = rclcpp::Duration(0.0);
        marker.mesh_use_embedded_materials = true;
        marker.mesh_resource = mesh;
        marker.scale.x = bbox[0];
        marker.scale.y = bbox[1];
        marker.scale.z = bbox[2];
        return marker;
    }

    void sendGazeboState()
    {
        auto request = std::make_shared<gazebo_msgs::srv::SetEntityState::Request>();
        request->state.name = ns_;

        request->state.pose.position.x = state_.pos.x;
        request->state.pose.position.y = state_.pos.y;
        request->state.pose.position.z = state_.pos.z;
        request->state.pose.orientation.x = state_.quat.x;
        request->state.pose.orientation.y = state_.quat.y;
        request->state.pose.orientation.z = state_.quat.z;
        request->state.pose.orientation.w = state_.quat.w;

        auto future = gazebo_client_->async_send_request(request);
        
        try
        {
            auto result = future.get();
            // RCLCPP_INFO(this->get_logger(), "Gazebo state sent");
        }
        catch(const std::exception& e)
        {
            // RCLCPP_ERROR(this->get_logger(), "Failed to call service");
        }
        
    }

    void pubTF()
    {
        static tf2_ros::TransformBroadcaster br(this);
        geometry_msgs::msg::TransformStamped t;

        t.header.stamp = this->get_clock()->now();
        t.header.frame_id = "map";
        t.child_frame_id = ns_ + "/base_link";

        t.transform.translation.x = state_.pos.x;
        t.transform.translation.y = state_.pos.y;
        t.transform.translation.z = state_.pos.z;

        t.transform.rotation.x = state_.quat.x;
        t.transform.rotation.y = state_.quat.y;
        t.transform.rotation.z = state_.quat.z;
        t.transform.rotation.w = state_.quat.w;

        br.sendTransform(t);
    }

    double eval(const std::string &expression, double t) {
        // Placeholder for the actual evaluation of the expression
        return t;  // Dummy implementation
    }

private:
    rclcpp::Publisher<panther_msgs::msg::DynTraj>::SharedPtr pub_traj_;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr pub_shapes_dynamic_mesh;
    std::vector<panther_msgs::msg::DynTraj> all_dyn_traj;
    visualization_msgs::msg::MarkerArray marker_array;
    rclcpp::TimerBase::SharedPtr timer_;
    int total_num_obs;
    int num_of_dyn_objects;
    int num_of_stat_objects;
    double x_min, x_max, y_min, y_max, z_min, z_max;
    std::vector<double> scale;
    std::vector<std::string> available_meshes_static = {"package://panther/meshes/ConcreteDamage01b/model3.dae", "package://panther/meshes/ConcreteDamage01b/model2.dae"};
    std::vector<std::string> available_meshes_dynamic = {"package://panther/meshes/ConcreteDamage01b/model4.dae"};
    std::vector<double> bbox_dynamic = {0.4, 0.4, 0.4};
    std::vector<double> bbox_static_vert = {0.4, 0.4, 4.0};
    bool gazebo;
    std::mt19937 random_generator;
};

int main(int argc, char **argv) 
{
    rclcpp::init(argc, argv);

    // Parse command line arguments
    if (argc > 1) {
        total_num_obs = std::stoi(argv[1]);
    }

    if (argc > 2) {
        std::string gazebo_str = argv[2];
        gazebo = (gazebo_str == "True" || gazebo_str == "true");
    }

    auto corridor_node = std::make_shared<DynCorridor>(total_num_obs, gazebo);
    rclcpp::executors::MultiThreadedExecutor executor;
    executor.add_node(corridor_node);
    executor.spin();

    rclcpp::shutdown();
    return 0;
}
