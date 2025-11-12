cmake_minimum_required(VERSION 3.8)
project(mighty)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# avoid Found-related error (ref: https://stackoverflow.com/questions/18327197/opencv-is-considered-to-be-not-found)
set(decomp_rviz_plugins_FOUND 1)

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(ament_cmake_python REQUIRED)
find_package(gazebo REQUIRED) # Note: if you put this line later on, it gives an error
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
find_package(gazebo_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(gazebo_ros REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(tf2_eigen REQUIRED)
find_package(std_srvs REQUIRED)
find_package(pcl_ros REQUIRED)
find_package(nlohmann_json REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# include directories
include_directories(${EIGEN3_INCLUDE_DIR} ${PCL_INCLUDE_DIRS} ${DECOMP_UTIL_INCLUDE_DIRS} ${GAZEBO_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS} ${BOOST_LIBRARY_DIRS} ${GAZEBO_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# add executables
add_executable(mighty src/mighty/mighty_node.cpp src/mighty/mighty.cpp src/dgp/dgp_manager.cpp src/dgp/utils.cpp src/dgp/dgp_planner.cpp src/dgp/graph_search.cpp src/mighty/utils.cpp src/mighty/lbfgs_solver.cpp src/mighty/lbfgs_solver_utils.cpp)
target_link_libraries(mighty ${PCL_LIBRARIES})
ament_target_dependencies(mighty rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 decomp_ros_msgs sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs tf2 tf2_ros tf2_eigen std_srvs pcl_ros)

# Add the fake_sim executable
add_executable(fake_sim src/mighty/fake_sim.cpp)
target_link_libraries(fake_sim)
ament_target_dependencies(fake_sim rclcpp rclcpp_components visualization_msgs dynus_interfaces sensor_msgs gazebo_msgs sensor_msgs tf2 tf2_ros nav_msgs)

# Add the obstracle tracker
add_executable(obstacle_tracker_node src/mighty/obstacle_tracker_node.cpp src/dgp/utils.cpp src/mighty/utils.cpp)
target_link_libraries(obstacle_tracker_node ${PCL_LIBRARIES})
ament_target_dependencies(obstacle_tracker_node rclcpp rclcpp_components visualization_msgs Eigen3 dynus_interfaces std_msgs gazebo_msgs sensor_msgs tf2 tf2_ros pcl_conversions tf2_geometry_msgs tf2_eigen)

# Add goal to cmd_vel converter
add_executable(convert_goal_to_cmd_vel src/mighty/convert_goal_to_cmd_vel.cpp)
target_link_libraries(convert_goal_to_cmd_vel)
ament_target_dependencies(convert_goal_to_cmd_vel rclcpp example_interfaces dynus_interfaces visualization_msgs Eigen3 sensor_msgs std_msgs decomp_rviz_plugins decomp_test_node tf2_geometry_msgs tf2 tf2_ros tf2_eigen)

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

# Add the convert velodyne executable
add_executable(convert_velodyne_to_ros_time src/mighty/convert_velodyne_to_ros_time.cpp)
target_link_libraries(convert_velodyne_to_ros_time)
ament_target_dependencies(convert_velodyne_to_ros_time rclcpp rclcpp_components visualization_msgs dynus_interfaces sensor_msgs gazebo_msgs sensor_msgs tf2 tf2_ros nav_msgs)

# For dynamic obstacles in Gazebo
add_library(move_model SHARED src/sim/move_model.cpp)
target_link_libraries(move_model ${GAZEBO_LIBRARIES})
ament_target_dependencies(move_model rclcpp tf2 tf2_ros geometry_msgs)

# For gazebo plugins (IMU)
add_library(imu_plugin SHARED src/sim/gazebo_ros_imu_sensor.cpp)
target_link_libraries(imu_plugin ${GAZEBO_LIBRARIES})
ament_target_dependencies(imu_plugin sensor_msgs gazebo_ros)

# For dynamic forest node
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

# install for executables
install(TARGETS
  mighty
  fake_sim
  obstacle_tracker_node
  convert_goal_to_cmd_vel
  convert_odom_to_state
  convert_vicon_to_state
  convert_velodyne_to_ros_time
  dynamic_forest_node
  DESTINATION lib/${PROJECT_NAME}
)

# install for libraries
install(TARGETS
  move_model
  imu_plugin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

# install Python modules
ament_python_install_package(${PROJECT_NAME})

# install Python executables
install(PROGRAMS
  scripts/goal_sender.py
  scripts/goal_monitor_node.py
  scripts/generate_random_forest.py
  benchmarking/dynamic_obstacle_benchmark_node.py
  benchmarking/unc_benchmark_node.py
  DESTINATION lib/${PROJECT_NAME}
)

ament_package()
