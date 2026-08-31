# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import os
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration, PythonExpression
from launch_ros.actions import Node
import launch_ros.actions
import launch_ros.descriptions
from ament_index_python.packages import get_package_share_directory
from launch.conditions import IfCondition, UnlessCondition

def generate_launch_description():

    # LaunchConfigurations
    init_x = LaunchConfiguration('init_x_', default=0.0)
    init_y = LaunchConfiguration('init_y_', default=0.0)
    init_z = LaunchConfiguration('init_z_', default=0.0)
    obj_num = LaunchConfiguration('obj_num', default=1)
    map_size_x_ = LaunchConfiguration('map_size_x_', default=12.5)
    map_size_y_ = LaunchConfiguration('map_size_y_', default=12.5)
    map_size_z_ = LaunchConfiguration('map_size_z_', default=6.0)
    c_num = LaunchConfiguration('c_num', default=0)
    p_num = LaunchConfiguration('p_num', default=30)
    min_dist = LaunchConfiguration('min_dist', default=2.0)
    odometry_topic = LaunchConfiguration('odometry_topic', default='visual_slam/odom')

    # DeclareLaunchArguments
    init_x_arg = DeclareLaunchArgument('init_x_', default_value=init_x, description='Initial X position')
    init_y_arg = DeclareLaunchArgument('init_y_', default_value=init_y, description='Initial Y position')
    init_z_arg = DeclareLaunchArgument('init_z_', default_value=init_z, description='Initial Z position')
    obj_num_arg = DeclareLaunchArgument('obj_num', default_value=obj_num, description='Number of objects')
    map_size_x_arg = DeclareLaunchArgument('map_size_x_', default_value=map_size_x_, description='Map size X')
    map_size_y_arg = DeclareLaunchArgument('map_size_y_', default_value=map_size_y_, description='Map size Y')
    map_size_z_arg = DeclareLaunchArgument('map_size_z_', default_value=map_size_z_, description='Map size Z')
    c_num_arg = DeclareLaunchArgument('c_num', default_value=c_num, description='Circle number')
    p_num_arg = DeclareLaunchArgument('p_num', default_value=p_num, description='Polygon number')
    min_dist_arg = DeclareLaunchArgument('min_dist', default_value=min_dist, description='Minimum distance')
    odometry_topic_arg = DeclareLaunchArgument('odometry_topic', default_value=odometry_topic, description='Odometry topic')
    
    use_mockamap = LaunchConfiguration('use_mockamap', default=False) # map_generator or mockamap
    use_mockamap_arg = DeclareLaunchArgument('use_mockamap', default_value=use_mockamap, description='Choose map type, map_generator or mockamap')

    # Obstacle-free worlds (e.g. the ground-robot position-exchange demo):
    # launch_map_generator:=false drops the random forest, and
    # use_empty_world_cloud:=true starts a constant SUB-FLOOR point grid on
    # /map_generator/global_cloud instead. The planner cannot run without a
    # continuously arriving sensor cloud (checkReadyToReplan requires an
    # initialized map, and pcl_render publishes nothing without a global map
    # — an EMPTY cloud would leave its FLANN index unbuilt). Points at
    # z = -3 m keep pcl_render fed while staying below every agent's z_min,
    # so the voxel map sees a world with no obstacles at all.
    launch_map_generator = LaunchConfiguration('launch_map_generator', default='true')
    launch_map_generator_arg = DeclareLaunchArgument(
        'launch_map_generator', default_value='true',
        description='Start the random_forest obstacle generator (false = obstacle-free world)')
    use_empty_world_cloud = LaunchConfiguration('use_empty_world_cloud', default='false')
    use_empty_world_cloud_arg = DeclareLaunchArgument(
        'use_empty_world_cloud', default_value='false',
        description='Publish the sub-floor global cloud that feeds pcl_render in obstacle-free worlds')
    use_dynamic = LaunchConfiguration('use_dynamic', default=True)  
    use_dynamic_arg = DeclareLaunchArgument('use_dynamic', default_value=use_dynamic, description='Use Drone Simulation Considering Dynamics or Not')
    default_rviz_config = os.path.join(get_package_share_directory('mighty'), 'rviz', 'mighty.rviz')
    rviz_config = LaunchConfiguration('rviz_config', default=default_rviz_config)
    rviz_config_arg = DeclareLaunchArgument('rviz_config', default_value=default_rviz_config, description='Path to RViz config file')

    # Node Definitions
    random_forest_node = Node(
        package='map_generator',
        executable='random_forest',
        name='random_forest',
        output='screen',
        remappings=[
            ('odometry', odometry_topic)
        ],
        parameters=[
            {'map/x_size': map_size_x_},
            {'map/y_size': map_size_y_},
            {'map/z_size': map_size_z_},
            {'map/resolution': 0.1},
            {'ObstacleShape/seed': 0},
            {'map/obs_num': p_num},
            {'ObstacleShape/lower_rad': 0.8},
            {'ObstacleShape/upper_rad': 0.8},
            {'ObstacleShape/lower_hei': 3.0},
            {'ObstacleShape/upper_hei': 3.0},
            {'map/circle_num': c_num},
            {'ObstacleShape/radius_l': 0.8},
            {'ObstacleShape/radius_h': 0.8},
            {'ObstacleShape/z_l': 1.0},
            {'ObstacleShape/z_h': 1.0},
            {'ObstacleShape/theta': 0.5},
            {'sensing/radius': 10.0},
            {'sensing/rate': 50.0},
            {'min_distance': min_dist}
        ],
        # Runs only when a map is wanted at all AND mockamap wasn't chosen —
        # preserves the historical use_mockamap semantics.
        condition=IfCondition(PythonExpression([
            "'", launch_map_generator, "'.lower() == 'true' and '",
            use_mockamap, "'.lower() != 'true'"]))
    )

    subfloor_cloud_node = Node(
        package='mighty',
        executable='publish_subfloor_cloud.py',
        name='subfloor_cloud',
        output='screen',
        condition=IfCondition(use_empty_world_cloud)
    )

    # use_rviz:=false runs the fake_sim scenarios fully headless, which is what
    # the automated campaign (scripts/exploration_test.py) needs: RViz on a
    # shared display just contends for the GPU during long unattended runs.
    use_rviz = LaunchConfiguration('use_rviz', default='true')
    use_rviz_arg = DeclareLaunchArgument(
        'use_rviz', default_value='true',
        description='Launch RViz (set false for headless / unattended runs)')

    rviz_node = launch_ros.actions.Node(
            package='rviz2', executable='rviz2', output='screen',
            arguments=['--display-config', rviz_config],
            condition=IfCondition(use_rviz))

    # Create LaunchDescription
    ld = LaunchDescription()

    # Add LaunchArguments
    ld.add_action(init_x_arg)
    ld.add_action(init_y_arg)
    ld.add_action(init_z_arg)
    ld.add_action(obj_num_arg)
    ld.add_action(map_size_x_arg)
    ld.add_action(map_size_y_arg)
    ld.add_action(map_size_z_arg)
    ld.add_action(c_num_arg)
    ld.add_action(p_num_arg)
    ld.add_action(min_dist_arg)
    ld.add_action(odometry_topic_arg)
    ld.add_action(use_mockamap_arg)
    ld.add_action(use_dynamic_arg)
    ld.add_action(rviz_config_arg)
    ld.add_action(use_rviz_arg)
    ld.add_action(launch_map_generator_arg)
    ld.add_action(use_empty_world_cloud_arg)

    ld.add_action(random_forest_node)
    ld.add_action(subfloor_cloud_node)
    ld.add_action(rviz_node)

    return ld