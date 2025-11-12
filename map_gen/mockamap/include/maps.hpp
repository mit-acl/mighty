#ifndef MAPS_HPP
#define MAPS_HPP

#include <rclcpp/rclcpp.hpp>  // ROS2 core
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <sensor_msgs/msg/point_cloud2.hpp>  // ROS2 message

namespace mocka {

class Maps {
public:
  struct BasicInfo {
    // Use rclcpp::Node::SharedPtr or reference instead of ros::NodeHandle*
    rclcpp::Node::SharedPtr node;  

    int sizeX;
    int sizeY;
    int sizeZ;
    int seed;
    double scale;

    // Use smart pointers for ROS2 messages
    sensor_msgs::msg::PointCloud2::SharedPtr output;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud;
  };

  BasicInfo getInfo() const;
  void setInfo(const BasicInfo &value);

public:
  Maps();

  void generate(int type);

private:
  BasicInfo info;

  void pcl2ros();

  void perlin3D();
  void maze2D();
  void randomMapGenerate();
  void Maze3DGen();
  void recursiveDivision(int xl, int xh, int yl, int yh, Eigen::MatrixXi &maze);
  void recursizeDivisionMaze(Eigen::MatrixXi &maze);
  void optimizeMap();
};

class MazePoint {
private:
  pcl::PointXYZ point;
  double dist1;
  double dist2;
  int point1;
  int point2;
  bool isdoor;

public:
  pcl::PointXYZ getPoint();
  int getPoint1();
  int getPoint2();
  double getDist1();
  double getDist2();
  void setPoint(pcl::PointXYZ p);
  void setPoint1(int p);
  void setPoint2(int p);
  void setDist1(double set);
  void setDist2(double set);
};

} // namespace mocka

#endif // MAPS_HPP
