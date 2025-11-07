/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "rclcpp/rclcpp.hpp"
#include "rclcpp/executors/multi_threaded_executor.hpp"
#include "rclcpp/callback_group.hpp"
#include <ament_index_cpp/get_package_share_directory.hpp>
#include "geometry_msgs/msg/transform_stamped.hpp"
#include "gazebo_msgs/srv/set_entity_state.hpp"
#include "gazebo_msgs/msg/entity_state.hpp"
#include "tf2_ros/transform_broadcaster.h"
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include "visualization_msgs/msg/marker.hpp"
#include "tf2/LinearMath/Quaternion.h"
#include <math.h>
#include <chrono>
#include <thread>
#include <Eigen/StdVector>

using namespace std::chrono_literals;

class FakeSim : public rclcpp::Node
{

public:

    // Constructor
    FakeSim() : Node("fake_sim"), br_(this, rclcpp::QoS(10).reliable().durability_volatile()) // Note: adding rclcpp::QoS(10).reliable().durability_volatile() to tf broadcaster helped a lot with the visualization lag.
    {

        RCLCPP_INFO(this->get_logger(), "Initializing FakeSim...");

        // Initialize callback groups
        cb_group_me_1_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
        cb_group_re_1_ = this->create_callback_group(rclcpp::CallbackGroupType::Reentrant);

        // Declare parameters
        this->declare_parameter<std::vector<double>>("start_pos", {0.0, 0.0, 4.0});
        this->declare_parameter<double>("start_yaw", -1.57);
        this->declare_parameter<bool>("send_state_to_gazebo", true);
        this->declare_parameter<double>("default_goal_z", 0.3);
        this->declare_parameter<int>("visual_level", 0);

        // Get parameters
        auto start_pos = this->get_parameter("start_pos").as_double_array();
        double yaw = this->get_parameter("start_yaw").as_double();
        send_state_to_gazebo_ = this->get_parameter("send_state_to_gazebo").as_bool();
        default_goal_z_ = this->get_parameter("default_goal_z").as_double();
        int visual_level = this->get_parameter("visual_level").as_int();

        // Print parameters
        RCLCPP_INFO(this->get_logger(), "Start position: %f, %f, %f", start_pos[0], start_pos[1], start_pos[2]);
        RCLCPP_INFO(this->get_logger(), "Start yaw: %f", yaw);
        RCLCPP_INFO(this->get_logger(), "Send state to Gazebo: %d", send_state_to_gazebo_);
        RCLCPP_INFO(this->get_logger(), "Default goal z: %f", default_goal_z_);
        RCLCPP_INFO(this->get_logger(), "Visual level: %d", visual_level);

        // Initialize state
        state_ = dynus_interfaces::msg::State();
        state_.header.frame_id = "map";
        state_.pos.x = start_pos[0];
        state_.pos.y = start_pos[1];
        state_.pos.z = start_pos[2];
        double pitch = 0.0, roll = 0.0;
        tf2::Quaternion quat;
        quat.setRPY(roll, pitch, yaw);
        state_.quat.x = quat.x();
        state_.quat.y = quat.y();
        state_.quat.z = quat.z();
        state_.quat.w = quat.w();

        // Get namespace
        ns_ = get_namespace();
        if (!ns_.empty() && ns_[0] == '/')
            ns_ = ns_.substr(1); // Create a substring starting from the second character

        // Publishers
        pub_state_ = this->create_publisher<dynus_interfaces::msg::State>("state", rclcpp::QoS(10).reliable().durability_volatile());
        pub_marker_drone_ = this->create_publisher<visualization_msgs::msg::Marker>("drone_marker", rclcpp::QoS(10).reliable().durability_volatile());

        // Subscribers
        sub_goal_ = this->create_subscription<dynus_interfaces::msg::Goal>("goal", 10, std::bind(&FakeSim::goalCallback, this, std::placeholders::_1));

        // Timer to simulate TF broadcast
        // Not sure buf calling TF periodically suffers from latency
        timer_ = this->create_wall_timer(10ms, std::bind(&FakeSim::pubCallback, this), cb_group_me_1_);

        // Initialize the tf2 buffer and listener
        tf2_buffer_ = std::make_shared<tf2_ros::Buffer>(this->get_clock());
        tf2_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf2_buffer_);
        target_frame_ = ns_ + "/base_link";

        // Initialize t_
        t_.header.frame_id = "map";
        t_.child_frame_id = target_frame_;

        // Gazebo service client
        gazebo_client_ = this->create_client<gazebo_msgs::srv::SetEntityState>("/plug/set_entity_state");
        while (!gazebo_client_->wait_for_service(10s))
        {
            RCLCPP_INFO(this->get_logger(), "Gazebo service not available, waiting again...");
        }

        // Delay before sending the initial state to Gazebo
        // threads_.push_back(std::thread(&FakeSim::sendGazeboState, this));
        std::thread(&FakeSim::sendGazeboState, this).detach();

        // Flag to publish drone marker
        if (visual_level > 0)
        {
            publish_marker_drone_ = true;
        }
        else
        {
            publish_marker_drone_ = false;
        }

        // Publish the initial state
        std::this_thread::sleep_for(5s);
        state_.header.stamp = this->get_clock()->now();
        pub_state_->publish(state_);

        // Package path
        package_path_ = ament_index_cpp::get_package_share_directory("mighty");

        RCLCPP_INFO(this->get_logger(), "Package path: %s", package_path_.c_str());
        RCLCPP_INFO(this->get_logger(), "FakeSim initialized");

    } // end of constructor

private:
    // std::vector<std::thread> threads_; // thread to handle service calls
    std::string package_path_;
    std::string ns_;
    rclcpp::CallbackGroup::SharedPtr cb_group_me_1_;
    rclcpp::CallbackGroup::SharedPtr cb_group_re_1_;
    rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr pub_state_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_marker_drone_;
    rclcpp::Subscription<dynus_interfaces::msg::Goal>::SharedPtr sub_goal_;
    rclcpp::Client<gazebo_msgs::srv::SetEntityState>::SharedPtr gazebo_client_;
    rclcpp::TimerBase::SharedPtr timer_;
    dynus_interfaces::msg::State state_;
    bool publish_marker_drone_;
    bool send_state_to_gazebo_ = true;
    std::shared_ptr<tf2_ros::Buffer> tf2_buffer_;
    std::shared_ptr<tf2_ros::TransformListener> tf2_listener_;
    std::string target_frame_;
    double default_goal_z_;
    int drone_marker_id_ = 1;
    tf2_ros::TransformBroadcaster br_;
    geometry_msgs::msg::TransformStamped t_;

    // This is for ground robot
    void updateStateFromTF()
    {
        try
        {
            // Lookup the transform from "map" to the namespace-specific base link
            geometry_msgs::msg::TransformStamped transform = tf2_buffer_->lookupTransform("map", target_frame_, tf2::TimePointZero);

            // Update state based on the transform
            state_.header.stamp = this->get_clock()->now();
            state_.pos.x = transform.transform.translation.x;
            state_.pos.y = transform.transform.translation.y;

            // Ground is detected as obstacles, so DYNUS is planning at default_goal_z_ altitude, so we set the z position to default_goal_z_ here
            // This doesn't affect ground robot controller in convert_goal_to_cmd_vel.cpp since it doesn't use z position
            state_.pos.z = default_goal_z_; // instead of state_.pos.z = transform.transform.translation.z;

            state_.quat.x = transform.transform.rotation.x;
            state_.quat.y = transform.transform.rotation.y;
            state_.quat.z = transform.transform.rotation.z;
            state_.quat.w = transform.transform.rotation.w;

            // Publish the updated state
            pub_state_->publish(state_);
        }
        catch (const tf2::TransformException &ex)
        {
            RCLCPP_WARN(this->get_logger(), "Could not transform: %s", ex.what());
        }
    }

    void goalCallback(const dynus_interfaces::msg::Goal::SharedPtr data)
    {
        // Hopf fibration approach
        Eigen::Vector3d thrust;
        thrust << data->a.x, data->a.y, data->a.z + 9.81;
        Eigen::Vector3d thrust_normaized = thrust.normalized();

        double a, b, c;
        a = thrust_normaized.x();
        b = thrust_normaized.y();
        c = thrust_normaized.z();

        tf2::Quaternion qabc;
        tf2::Quaternion qpsi;

        double tmp = 1 / std::sqrt(2 * (1 + c));

        qabc.setValue(
            -b * tmp,
            a * tmp,
            0.0,
            tmp * (1 + c));

        qpsi.setValue(
            0.0,
            0.0,
            sin(data->yaw / 2),
            cos(data->yaw / 2));

        tf2::Quaternion w_q_b = qabc * qpsi;

        // Publish state if not using ground robot (if using ground robot, we publish state using tf with updateStateFromTF)
        state_.header.stamp = this->get_clock()->now();
        state_.pos = data->p;
        state_.vel = data->v;
        state_.quat.w = w_q_b.w();
        state_.quat.x = w_q_b.x();
        state_.quat.y = w_q_b.y();
        state_.quat.z = w_q_b.z();

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
        catch (const std::exception &e)
        {
            // RCLCPP_ERROR(this->get_logger(), "Failed to call service");
        }
    }

    void getTransformStamped()
    {
        t_.header.stamp = this->get_clock()->now();

        t_.transform.translation.x = state_.pos.x;
        t_.transform.translation.y = state_.pos.y;
        t_.transform.translation.z = state_.pos.z;

        t_.transform.rotation.x = state_.quat.x;
        t_.transform.rotation.y = state_.quat.y;
        t_.transform.rotation.z = state_.quat.z;
        t_.transform.rotation.w = state_.quat.w;
    }

    void pubCallback()
    {
        // Publish the transform
        getTransformStamped();
        br_.sendTransform(t_);

        // Publish drone marker
        if (publish_marker_drone_)
            pub_marker_drone_->publish(getDroneMarker());

        // Send the state to Gazebo
        if (send_state_to_gazebo_)
        {
            // threads_.push_back(std::thread(&FakeSim::sendGazeboState, this));
            std::thread(&FakeSim::sendGazeboState, this).detach();
        }

        // Publish the state
        state_.header.stamp = this->get_clock()->now();
        pub_state_->publish(state_);
    }

    visualization_msgs::msg::Marker getDroneMarker()
    {
        visualization_msgs::msg::Marker marker;
        marker.id = drone_marker_id_;
        marker.ns = std::string("mesh_") + this->get_namespace();
        marker.header.frame_id = "map";
        marker.header.stamp = this->get_clock()->now();
        marker.type = marker.MESH_RESOURCE;
        marker.action = marker.ADD;

        marker.pose.position.x = state_.pos.x;
        marker.pose.position.y = state_.pos.y;
        marker.pose.position.z = state_.pos.z;
        marker.pose.orientation.x = state_.quat.x;
        marker.pose.orientation.y = state_.quat.y;
        marker.pose.orientation.z = state_.quat.z;
        marker.pose.orientation.w = state_.quat.w;

        marker.mesh_use_embedded_materials = true;
        marker.mesh_resource = "package://mighty/meshes/quadrotor/quadrotor.dae";
        std::cout << "Mesh resource: " << marker.mesh_resource << std::endl;
        marker.scale.x = 0.75;
        marker.scale.y = 0.75;
        marker.scale.z = 0.75;

        return marker;
    }
};

int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<FakeSim>();
    rclcpp::executors::MultiThreadedExecutor executor;
    executor.add_node(node);
    executor.spin();

    rclcpp::shutdown();
    return 0;
}
