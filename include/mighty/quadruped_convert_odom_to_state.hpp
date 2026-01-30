
#include "rclcpp/rclcpp.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "std_msgs/msg/header.hpp"
#include "geometry_msgs/msg/vector3.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "mighty/mighty.hpp" // TODO<investigate>: To make dynus_interfaces::msg::State work, we need to include mighty/mighty.hpp, but I'm not sure why. Investigate this.
#include <dynus_interfaces/msg/state.hpp>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/msg/transform_stamped.hpp>

#ifndef ODOMETRY_TO_STATE_NODE_HPP
#define ODOMETRY_TO_STATE_NODE_HPP

class OdometryToStateNode : public rclcpp::Node
{
public:
    OdometryToStateNode();
    void timer_callback();

private:
    void callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg);

    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odometry_sub_;
    rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr state_publisher_;

    std::shared_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
    geometry_msgs::msg::TransformStamped tf_msg_;

    // Set up a timer to publish the TF transform at 100 Hz
    rclcpp::TimerBase::SharedPtr timer_;

    // mutex for the tf_msg_
    std::mutex tf_msg_mutex_;
};

#endif // ODOMETRY_TO_STATE_NODE_HPP