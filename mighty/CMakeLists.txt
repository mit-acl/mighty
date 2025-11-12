cmake_minimum_required(VERSION 3.5)
project(mighty)

# Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
add_compile_options(-Wall -O3 -fPIC)

# Dependencies
find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(ompl REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(nav_msgs REQUIRED)

# Include directories
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS}
)

# Library
add_library(mighty src/mighty.cpp)

target_link_libraries(mighty
  ${OMPL_LIBRARIES}
)

ament_target_dependencies(mighty
  rclcpp
  std_msgs
  geometry_msgs
  sensor_msgs
  visualization_msgs
)

install(DIRECTORY
  launch
  config
  DESTINATION share/${PROJECT_NAME}
)

# Executable
add_executable(mighty_node src/mighty.cpp src/mighty_solver.cpp)
target_link_libraries(mighty_node mighty ${OMPL_LIBRARIES})
ament_target_dependencies(mighty_node
  rclcpp
  visualization_msgs
  nav_msgs
  geometry_msgs
  std_msgs
  sensor_msgs
)

# Install targets
install(TARGETS
  mighty
  mighty_node
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY include/
  DESTINATION include/
)

# Export package
ament_export_include_directories(include)
ament_export_libraries(mighty)
ament_export_dependencies(
  rclcpp
  std_msgs
  geometry_msgs
  sensor_msgs
  visualization_msgs
  ompl
  Eigen3
)

ament_package()
