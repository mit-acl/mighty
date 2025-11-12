cmake_minimum_required(VERSION 3.5)
project(mockamap)

# Require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
add_compile_options(-Wall -O3 -fPIC)

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(pcl_conversions REQUIRED)
find_package(pcl_ros REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)

include_directories(
  include
)

# Gather source files
file(GLOB_RECURSE LIB_SRCS "src/*.cpp")      # For library sources, put your lib code in src/lib/
file(GLOB_RECURSE NODE_SRCS "src/*.cpp")    # For node sources, put your executable code in src/node/

# Create library target
add_library(${PROJECT_NAME}_lib src/maps.cpp)

ament_target_dependencies(${PROJECT_NAME}_lib
  rclcpp
  pcl_ros
  pcl_conversions
  nav_msgs
  sensor_msgs
  geometry_msgs
  visualization_msgs
)

target_include_directories(${PROJECT_NAME}_lib PUBLIC include)

# Create executable target
add_executable(${PROJECT_NAME}_node ${NODE_SRCS})

# Link executable with library
target_link_libraries(${PROJECT_NAME}_node ${PROJECT_NAME}_lib)

ament_target_dependencies(${PROJECT_NAME}_node
  rclcpp
  pcl_ros
  pcl_conversions
  nav_msgs
  sensor_msgs
  geometry_msgs
  visualization_msgs
)

install(TARGETS
  ${PROJECT_NAME}_node
  ${PROJECT_NAME}_lib
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY include/
  DESTINATION include/
)

ament_package()
