#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>

#include "maps.hpp"

void optimizeMap(mocka::Maps::BasicInfo& in)
{
  std::vector<int> temp;

  pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

  cloud->width  = in.cloud->width;
  cloud->height = in.cloud->height;
  cloud->points.resize(cloud->width * cloud->height);

  for (size_t i = 0; i < cloud->width; i++)
  {
    cloud->points[i].x = in.cloud->points[i].x;
    cloud->points[i].y = in.cloud->points[i].y;
    cloud->points[i].z = in.cloud->points[i].z;
  }

  kdtree.setInputCloud(cloud);
  double radius = 1.75 / in.scale;  // 1.75 is rounded up sqrt(3)

  for (size_t i = 0; i < cloud->width; i++)
  {
    std::vector<int> pointIdxRadiusSearch;
    std::vector<float> pointRadiusSquaredDistance;

    if (kdtree.radiusSearch(cloud->points[i], radius, pointIdxRadiusSearch,
                            pointRadiusSquaredDistance) >= 27)
    {
      temp.push_back(i);
    }
  }

  // Erase points in reverse order to maintain valid indices
  for (auto it = temp.rbegin(); it != temp.rend(); ++it)
  {
    in.cloud->points.erase(in.cloud->points.begin() + *it);
  }

  in.cloud->width -= static_cast<int>(temp.size());

  pcl::toROSMsg(*in.cloud, *in.output);
  in.output->header.frame_id = "odom";

  RCLCPP_INFO(rclcpp::get_logger("maps"), "finish: number of points after optimization %d", in.cloud->width);
}

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);

  auto node = rclcpp::Node::make_shared("mockamap");

  auto pcl_pub = node->create_publisher<sensor_msgs::msg::PointCloud2>("mock_map", 1);

  auto cloud = std::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
  auto output = std::make_shared<sensor_msgs::msg::PointCloud2>();

  int seed;
  int sizeX, sizeY, sizeZ;
  double scale;
  double update_freq;
  int type;

  node->declare_parameter("seed", 4546);
  node->declare_parameter("update_freq", 1.0);
  node->declare_parameter("resolution", 0.38);
  node->declare_parameter("x_length", 100);
  node->declare_parameter("y_length", 100);
  node->declare_parameter("z_length", 10);
  node->declare_parameter("type", 1);

  node->get_parameter("seed", seed);
  node->get_parameter("update_freq", update_freq);
  node->get_parameter("resolution", scale);
  node->get_parameter("x_length", sizeX);
  node->get_parameter("y_length", sizeY);
  node->get_parameter("z_length", sizeZ);
  node->get_parameter("type", type);

  scale = 1 / scale;
  sizeX = static_cast<int>(sizeX * scale);
  sizeY = static_cast<int>(sizeY * scale);
  sizeZ = static_cast<int>(sizeZ * scale);

  mocka::Maps::BasicInfo info;
  // Note: You need to update Maps::BasicInfo to hold rclcpp::Node::SharedPtr instead of ros::NodeHandle*
  info.node = node; 
  info.sizeX      = sizeX;
  info.sizeY      = sizeY;
  info.sizeZ      = sizeZ;
  info.seed       = seed;
  info.scale      = scale;
  info.output     = std::shared_ptr<sensor_msgs::msg::PointCloud2>(output);
  info.cloud      = cloud;

  mocka::Maps map;
  map.setInfo(info);
  map.generate(type);

  // optimizeMap(info);

  rclcpp::Rate loop_rate(update_freq);

  while (rclcpp::ok())
  {
    pcl_pub->publish(*output);
    rclcpp::spin_some(node);
    loop_rate.sleep();
  }

  rclcpp::shutdown();
  return 0;
}
