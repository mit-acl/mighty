/* ----------------------------------------------------------------------------
 * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/point_cloud2.hpp"

// This node listens to /velodyne_points (which are stamped with simulation time)
// and republishes them with the header.stamp set to ROS2 system time.
class VelodyneTimeCorrector : public rclcpp::Node
{
public:
  VelodyneTimeCorrector()
  : Node("velodyne_time_corrector")
  {
    RCLCPP_INFO(this->get_logger(), "VelodyneTimeCorrector node started.");

    // Publisher: republish the PointCloud2 message with ROS2 system time stamp.
    publisher_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(
      "/velodyne_points_ros_time", 10);

    // Subscription: listen to the original velodyne points.
    subscription_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(
      "/velodyne_points", 10,
      std::bind(&VelodyneTimeCorrector::pointCloudCallback, this, std::placeholders::_1));
  }

private:
  void pointCloudCallback(const sensor_msgs::msg::PointCloud2::SharedPtr msg)
  {
    // Make a copy of the incoming message so we can update the timestamp.
    auto corrected_msg = std::make_shared<sensor_msgs::msg::PointCloud2>(*msg);

    // Use a ROS2 time stamp.
    corrected_msg->header.stamp = this->now();

    // Optionally, you can also change the frame_id or any other header information.
    // corrected_msg->header.frame_id = "your_preferred_frame";

    publisher_->publish(*corrected_msg);
  }

  rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr publisher_;
  rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr subscription_;
};

int main(int argc, char * argv[])
{
  rclcpp::init(argc, argv);
  auto node = std::make_shared<VelodyneTimeCorrector>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
