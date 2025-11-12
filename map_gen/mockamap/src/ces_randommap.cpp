#include <iostream>
#include <vector>
#include <deque>
#include <cmath>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/kdtree.h>

#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geometry_msgs/msg/vector3.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>

#include <rclcpp/rclcpp.hpp>

#include <Eigen/Eigen>
#include <Eigen/SVD>

#include <random>
#include <chrono>

#include "visualization_msgs/msg/marker.hpp"
#include "visualization_msgs/msg/marker_array.hpp"

// Your project headers
#include "maps.hpp"

using namespace std;
using namespace mocka;

#if MAP_OR_WORLD
const std::string kFrameIdNs_ = "/map";
#else
const std::string kFrameIdNs_ = "/world";
#endif

// Globals — should ideally be class members
pcl::search::KdTree<pcl::PointXYZ> kdtreeLocalMap;
std::vector<int> pointIdxRadiusSearch;
std::vector<float> pointRadiusSquaredDistance;

rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr local_map_pub;
rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr local_map_inflate_pub;
rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr global_map_pub;

rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr map_sub;
rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub;

std::deque<nav_msgs::msg::Odometry> odom_queue;
std::vector<double> state;
const size_t odom_queue_size = 200;
nav_msgs::msg::Odometry odom;

double z_limit;
double sensing_rate;
double sensing_range;

bool map_ok = false;
bool has_odom = false;

sensor_msgs::msg::PointCloud2 globalMap_pcd;
sensor_msgs::msg::PointCloud2 localMap_pcd;
pcl::PointCloud<pcl::PointXYZ> cloudMap;
rclcpp::Time begin_time{0};

typedef Eigen::Vector3d ObsPos;
typedef Eigen::Vector3d ObsSize; // x, y, height --- z
typedef std::pair<ObsPos, ObsPos> Obstacle;
std::vector<Obstacle> obstacle_list;

void fixedMapGenerate()
{
  double resolution = 0.1; // You need to initialize this properly

  cloudMap.points.clear();

  obstacle_list.push_back(make_pair(ObsPos(-7.0, 1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(-1.0, 1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(10.0, 1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(16.0, 1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(-4.0, -1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(13.0, -1.0, 0.0), ObsSize(1.0, 3.0, 5.0)));

  obstacle_list.push_back(make_pair(ObsPos(5.0, 2.5, 0.0), ObsSize(30.0, 1.0, 5.0)));
  obstacle_list.push_back(make_pair(ObsPos(5.0, -2.5, 0.0), ObsSize(30.0, 1.0, 5.0)));

  int num_total_obs = obstacle_list.size();
  pcl::PointXYZ pt_insert;

  for (int i = 0; i < num_total_obs; i++)
  {
    double x = obstacle_list[i].first[0];
    double y = obstacle_list[i].first[1];
    double z = obstacle_list[i].first[2];
    double lx = obstacle_list[i].second[0];
    double ly = obstacle_list[i].second[1];
    double lz = obstacle_list[i].second[2];

    int num_mesh_x = std::ceil(lx / resolution);
    int num_mesh_y = std::ceil(ly / resolution);
    int num_mesh_z = std::ceil(lz / resolution);

    int left_x = -num_mesh_x / 2;
    int right_x = num_mesh_x / 2;
    int left_y = -num_mesh_y / 2;
    int right_y = num_mesh_y / 2;
    int left_z = 0;
    int right_z = num_mesh_z;

    for (int r = left_x; r < right_x; r++)
      for (int s = left_y; s < right_y; s++)
        for (int t = left_z; t < right_z; t++)
        {
          if ((r - left_x) * (r - right_x + 1) * (s - left_y) * (s - right_y + 1) * (t - left_z) * (t - right_z + 1) == 0)
          {
            pt_insert.x = x + r * resolution;
            pt_insert.y = y + s * resolution;
            pt_insert.z = z + t * resolution;
            cloudMap.points.push_back(pt_insert);
          }
        }
  }

  cloudMap.width = cloudMap.points.size();
  cloudMap.height = 1;
  cloudMap.is_dense = true;

  RCLCPP_WARN(rclcpp::get_logger("mocka"), "Finished generate random map");
  std::cout << cloudMap.size() << std::endl;
  kdtreeLocalMap.setInputCloud(cloudMap.makeShared());
  map_ok = true;
}

void rcvOdometryCallback(const nav_msgs::msg::Odometry::SharedPtr odom_msg)
{
  if (odom_msg->child_frame_id == "X" || odom_msg->child_frame_id == "O")
    return;

  odom = *odom_msg;
  has_odom = true;

  state = {
    odom.pose.pose.position.x,
    odom.pose.pose.position.y,
    odom.pose.pose.position.z,
    odom.twist.twist.linear.x,
    odom.twist.twist.linear.y,
    odom.twist.twist.linear.z,
    0.0,
    0.0,
    0.0};

  odom_queue.push_back(*odom_msg);
  while (odom_queue.size() > odom_queue_size)
    odom_queue.pop_front();
}

int frequency_division_global = 40;

void publishAllPoints()
{
  if (!map_ok)
    return;

  rclcpp::Time now = rclcpp::Clock().now();

  if ((now - begin_time).seconds() > 7.0)
    return;

  frequency_division_global--;
  if (frequency_division_global == 0)
  {
    pcl::toROSMsg(cloudMap, globalMap_pcd);
    globalMap_pcd.header.frame_id = kFrameIdNs_;
    global_map_pub->publish(globalMap_pcd);
    frequency_division_global = 40;
    RCLCPP_ERROR(rclcpp::get_logger("mocka"), "[SERVER] Publish one global map");
  }
}

void pubSensedPoints()
{
  if (!map_ok || !has_odom)
    return;

  rclcpp::Time time_bef_sensing = rclcpp::Clock().now();

  pcl::PointCloud<pcl::PointXYZ> localMap;

  pcl::PointXYZ searchPoint(state[0], state[1], state[2]);
  pointIdxRadiusSearch.clear();
  pointRadiusSquaredDistance.clear();

  pcl::PointXYZ ptInInflation;

  if (kdtreeLocalMap.radiusSearch(searchPoint, sensing_range, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
  {
    for (size_t i = 0; i < pointIdxRadiusSearch.size(); ++i)
    {
      ptInInflation = cloudMap.points[pointIdxRadiusSearch[i]];
      localMap.points.push_back(ptInInflation);
    }
  }

  pcl::PointXYZ pt_fix;
  pt_fix.x = state[0];
  pt_fix.y = state[1];
  pt_fix.z = 0.0;
  localMap.points.push_back(pt_fix);

  localMap.width = localMap.points.size();
  localMap.height = 1;
  localMap.is_dense = true;

  pcl::toROSMsg(localMap, localMap_pcd);

  localMap_pcd.header.frame_id = kFrameIdNs_;
  local_map_pub->publish(localMap_pcd);

  rclcpp::Time time_aft_sensing = rclcpp::Clock().now();

  if ((time_aft_sensing - begin_time).seconds() > 5.0)
    return;

  frequency_division_global--;
  if (frequency_division_global == 0)
  {
    pcl::toROSMsg(cloudMap, globalMap_pcd);
    globalMap_pcd.header.frame_id = kFrameIdNs_;
    global_map_pub->publish(globalMap_pcd);
    frequency_division_global = 40;
    RCLCPP_INFO(rclcpp::get_logger("mocka"), "[SERVER] Publish one global map");
  }
}
