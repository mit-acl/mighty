cmake_minimum_required(VERSION 3.8)
project(mighty)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Default build type" FORCE)
  message(STATUS "mighty: defaulting CMAKE_BUILD_TYPE to Release")
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# avoid Found-related error (ref: https://stackoverflow.com/questions/18327197/opencv-is-considered-to-be-not-found)
set(decomp_rviz_plugins_FOUND 1)

# Option to disable Gazebo-dependent targets (e.g. on Jetson/ARM64)
option(BUILD_SIMULATION "Build Gazebo simulation targets" ON)

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(ament_cmake_python REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclpy REQUIRED)
find_package(rclcpp_components REQUIRED)
find_package(example_interfaces REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(decomp_util REQUIRED)
find_package(decomp_ros_msgs REQUIRED)
find_package(decomp_rviz_plugins REQUIRED)
find_package(decomp_test_node REQUIRED)
find_package(PCL REQUIRED)
find_package(pcl_conversions REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(dynus_interfaces REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(tf2_eigen REQUIRED)
find_package(std_srvs REQUIRED)
find_package(pcl_ros REQUIRED)
find_package(nlohmann_json REQUIRED)

find_package(gazebo_msgs REQUIRED)
if(BUILD_SIMULATION)
  find_package(gazebo REQUIRED)
  find_package(gazebo_ros REQUIRED)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# include directories
include_directories(${EIGEN3_INCLUDE_DIR} ${PCL_INCLUDE_DIRS} ${DECOMP_UTIL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS} ${BOOST_LIBRARY_DIRS})
if(BUILD_SIMULATION)
  include_directories(${GAZEBO_INCLUDE_DIRS})
  link_directories(${GAZEBO_LIBRARY_DIRS})
endif()
add_definitions(${PCL_DEFINITIONS})

# add executables
add_executable(mighty src/mighty/mighty_node.cpp src/mighty/mighty.cpp src/hgp/hgp_manager.cpp src/hgp/utils.cpp src/hgp/hgp_planner.cpp src/hgp/graph_search.cpp src/mighty/utils.cpp src/mighty/lbfgs_solver.cpp src/mighty/lbfgs_solver_utils.cpp src/mighty/frontier_detector.cpp src/mighty/frontier_manager.cpp)
target_link_libraries(mighty ${PCL_LIBRARIES})
ament_target_dependencies(mighty rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 decomp_ros_msgs sensor_msgs std_msgs nav_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs tf2 tf2_ros tf2_eigen std_srvs pcl_ros)

# Add pure pursuit controller for ground robot
add_executable(pure_pursuit src/mighty/pure_pursuit.cpp)
target_link_libraries(pure_pursuit)
ament_target_dependencies(pure_pursuit rclcpp dynus_interfaces geometry_msgs visualization_msgs)

add_executable(trajectory_tracker src/mighty/trajectory_tracker.cpp)
target_link_libraries(trajectory_tracker)
ament_target_dependencies(trajectory_tracker rclcpp dynus_interfaces geometry_msgs visualization_msgs)

# Add odom to state converter
add_executable(convert_odom_to_state src/mighty/convert_odom_to_state.cpp)
target_link_libraries(convert_odom_to_state)
ament_target_dependencies(convert_odom_to_state rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs tf2 tf2_ros tf2_eigen)

# Add vicon to state converter
add_executable(convert_vicon_to_state src/mighty/convert_vicon_to_state.cpp)
target_link_libraries(convert_vicon_to_state)
ament_target_dependencies(convert_vicon_to_state rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 decomp_ros_msgs sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs tf2 tf2_ros tf2_eigen)

# Add the executable for the lbfgs solver test file
add_executable(test_lbfgs_solver src/test/test_lbfgs_solver.cpp src/mighty/lbfgs_solver.cpp src/mighty/lbfgs_solver_utils.cpp)
ament_target_dependencies(test_lbfgs_solver rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 decomp_ros_msgs sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs)

# Add the fake_sim executable (uses gazebo_msgs only, no Gazebo library needed)
add_executable(fake_sim src/mighty/fake_sim.cpp)
target_link_libraries(fake_sim)
ament_target_dependencies(fake_sim rclcpp rclcpp_components visualization_msgs dynus_interfaces sensor_msgs gazebo_msgs sensor_msgs tf2 tf2_ros nav_msgs)

# Simulation-only targets (require Gazebo, not available on ARM64/Jetson)
if(BUILD_SIMULATION)
  add_executable(convert_velodyne_to_ros_time src/mighty/convert_velodyne_to_ros_time.cpp)
  target_link_libraries(convert_velodyne_to_ros_time)
  ament_target_dependencies(convert_velodyne_to_ros_time rclcpp rclcpp_components visualization_msgs dynus_interfaces sensor_msgs gazebo_msgs sensor_msgs tf2 tf2_ros nav_msgs)

  add_library(move_model SHARED src/sim/move_model.cpp)
  target_link_libraries(move_model ${GAZEBO_LIBRARIES})
  ament_target_dependencies(move_model rclcpp tf2 tf2_ros geometry_msgs)

  add_library(imu_plugin SHARED src/sim/gazebo_ros_imu_sensor.cpp)
  target_link_libraries(imu_plugin ${GAZEBO_LIBRARIES})
  ament_target_dependencies(imu_plugin sensor_msgs gazebo_ros)

  add_executable(dynamic_forest_node src/sim/dynamic_forest_node.cpp)
  ament_target_dependencies(dynamic_forest_node
    rclcpp
    geometry_msgs
    gazebo_msgs
    visualization_msgs
    std_msgs
    tf2_ros
    dynus_interfaces
  )
endif()

# install for launch and config files
install(DIRECTORY
launch
config
rviz
worlds
urdf
meshes
DESTINATION share/${PROJECT_NAME}
)

# install for executables (always built)
install(TARGETS
  mighty
  fake_sim
  pure_pursuit
  trajectory_tracker
  convert_odom_to_state
  convert_vicon_to_state
  DESTINATION lib/${PROJECT_NAME}
)

# install for simulation executables and libraries
if(BUILD_SIMULATION)
  install(TARGETS
    convert_velodyne_to_ros_time
    dynamic_forest_node
    DESTINATION lib/${PROJECT_NAME}
  )
  install(TARGETS
    move_model
    imu_plugin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )
endif()

# install Python modules
ament_python_install_package(${PROJECT_NAME})

# install Python executables
install(PROGRAMS
  scripts/goal_sender.py
  scripts/goal_monitor_node.py
  scripts/generate_random_forest.py
  benchmarking/unc_benchmark_node.py
  scripts/repub_rviz_2Dgoal.py
  scripts/frame_align_publisher.py
  scripts/formation_viz_node.py
  DESTINATION lib/${PROJECT_NAME}
)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(test_occ_grid_2d_unknown
    src/test/test_occ_grid_2d_unknown.cpp)
  ament_target_dependencies(test_occ_grid_2d_unknown
    rclcpp Eigen3 nav_msgs)

  ament_add_gtest(test_frontier_detector
    src/test/test_frontier_detector.cpp
    src/mighty/frontier_detector.cpp)
  ament_target_dependencies(test_frontier_detector
    rclcpp Eigen3 nav_msgs)

  ament_add_gtest(test_frontier_manager
    src/test/test_frontier_manager.cpp
    src/mighty/frontier_detector.cpp
    src/mighty/frontier_manager.cpp)
  ament_target_dependencies(test_frontier_manager
    rclcpp Eigen3 nav_msgs)

  ament_add_gtest(test_gradient_check
    src/test/test_gradient_check.cpp
    src/mighty/lbfgs_solver.cpp
    src/mighty/lbfgs_solver_utils.cpp)
  ament_target_dependencies(test_gradient_check
    rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3
    decomp_ros_msgs sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node
    tf2_geometry_msgs)

  ament_add_gtest(test_formation_cost
    src/test/test_formation_cost.cpp
    src/mighty/lbfgs_solver.cpp
    src/mighty/lbfgs_solver_utils.cpp)
  ament_target_dependencies(test_formation_cost
    rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3
    decomp_ros_msgs sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node
    tf2_geometry_msgs)
endif()

ament_package()
