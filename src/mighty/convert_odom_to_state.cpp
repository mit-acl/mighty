#include <mighty/convert_odom_to_state.hpp>

OdometryToStateNode::OdometryToStateNode() : Node("odometry_to_state_node")
{
    // Subscriber for nav_msgs/Odometry
    odometry_sub_ = this->create_subscription<nav_msgs::msg::Odometry>(
        "odom", 10,
        std::bind(&OdometryToStateNode::callback, this, std::placeholders::_1));

    // Publisher for dynus_interfaces/State
    state_publisher_ = this->create_publisher<dynus_interfaces::msg::State>("state", 10);
}

void OdometryToStateNode::callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg)
{
    // Construct the State message
    auto state_msg = dynus_interfaces::msg::State();
    state_msg.header.stamp = odom_msg->header.stamp;
    state_msg.header.frame_id = odom_msg->header.frame_id;

    // Set position from Odometry
    state_msg.pos.x = odom_msg->pose.pose.position.x;
    state_msg.pos.y = odom_msg->pose.pose.position.y;
    state_msg.pos.z = odom_msg->pose.pose.position.z;

    // Set velocity from Odometry
    state_msg.vel.x = odom_msg->twist.twist.linear.x;
    state_msg.vel.y = odom_msg->twist.twist.linear.y;
    state_msg.vel.z = odom_msg->twist.twist.linear.z;

    // Set orientation from Odometry
    state_msg.quat.x = odom_msg->pose.pose.orientation.x;
    state_msg.quat.y = odom_msg->pose.pose.orientation.y;
    state_msg.quat.z = odom_msg->pose.pose.orientation.z;
    state_msg.quat.w = odom_msg->pose.pose.orientation.w;

    // Publish the State message
    state_publisher_->publish(state_msg);
}

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<OdometryToStateNode>();

    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
