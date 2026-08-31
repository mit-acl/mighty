# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import os
from launch import LaunchDescription
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
from launch.actions import IncludeLaunchDescription, DeclareLaunchArgument, OpaqueFunction
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import PathJoinSubstitution, LaunchConfiguration
from launch_ros.substitutions import FindPackageShare

def convert_str_to_bool(str):
    return True if (str == 'true' or str == 'True' or str == 1 or str == '1') else False

def generate_launch_description():
    
    # Declare a new argument "env" to choose the world file.
    env_arg = DeclareLaunchArgument(
        'env',
        default_value='easy_forest',
        description='Environment name to determine the Gazebo world file'
    )

    # Remove the previous world_path argument and use env instead.
    use_rviz_arg = DeclareLaunchArgument(
        'use_rviz', default_value='true', description='Flag to enable or disable RViz'
    )
    use_gazebo_gui_arg = DeclareLaunchArgument(
        'use_gazebo_gui', default_value='false', description='Flag to enable or disable gazebo gui'
    )
    use_dyn_obs_arg = DeclareLaunchArgument(
        'use_dyn_obs', default_value='false', description='Flag to enable or disable dynamic obstacles'
    )
    use_ground_robot_arg = DeclareLaunchArgument(
        'use_ground_robot', default_value='false', description='Use ground robot (affects RViz config)'
    )
    rviz_config_arg = DeclareLaunchArgument(
        'rviz_config', default_value='',
        description='RViz config filename to use instead of the auto-selected one '
                    '(bare name, resolved against mighty/rviz). Multi-robot ground '
                    'runs pass multi_mighty_sim_ground_robot.rviz.'
    )

    # Dynamic obstacle spawn region
    num_dyn_obstacles_arg = DeclareLaunchArgument('num_dyn_obstacles', default_value='1', description='Number of dynamic obstacles')
    dyn_x_min_arg = DeclareLaunchArgument('dyn_x_min', default_value='5.0', description='Dynamic obstacle spawn x min')
    dyn_x_max_arg = DeclareLaunchArgument('dyn_x_max', default_value='105.0', description='Dynamic obstacle spawn x max')
    dyn_y_min_arg = DeclareLaunchArgument('dyn_y_min', default_value='-5.0', description='Dynamic obstacle spawn y min')
    dyn_y_max_arg = DeclareLaunchArgument('dyn_y_max', default_value='5.0', description='Dynamic obstacle spawn y max')
    dyn_z_min_arg = DeclareLaunchArgument('dyn_z_min', default_value='1.0', description='Dynamic obstacle spawn z min')
    dyn_z_max_arg = DeclareLaunchArgument('dyn_z_max', default_value='5.0', description='Dynamic obstacle spawn z max')
    dyn_scale_z_min_arg = DeclareLaunchArgument('dyn_scale_z_min', default_value='2.0', description='Dynamic obstacle z scale min')
    dyn_scale_z_max_arg = DeclareLaunchArgument('dyn_scale_z_max', default_value='4.0', description='Dynamic obstacle z scale max')

    # benchmark name
    benchmark_name_arg = DeclareLaunchArgument('benchmark_name', default_value='benchmark_name', description='Benchmark name')

    # Opaque function to launch nodes
    def launch_setup(context, *args, **kwargs):
        
        # Get the environment value from the 'env' launch argument.
        env_value = LaunchConfiguration('env').perform(context)
        # Map environment names to corresponding Gazebo world file names.
        world_mapping = {
            'high_res_forest': 'big_forest_high_res.world',
            'static_uncertainty_test2': 'static_uncertainty_test2.world',
            'static_uncertainty_test3': 'static_uncertainty_test3.world',
            'static_uncertainty_test4': 'static_uncertainty_test4.world',
            'office_faster': 'office.world',
            'office': 'office.world',
            'cave_start': 'simple_tunnel.world',
            'cave_vertical': 'simple_tunnel.world',
            'cave_person': 'simple_tunnel.world',
            'forest3': 'forest3.world',
            'yaw_benchmark': 'forest3.world',
            'global_planner': 'forest3.world',
            'multiagent_performance': 'forest3.world',
            'path_push': 'forest3.world',
            'ACL_office': 'ACL_office.world',
            'ground_robot': 'ACL_office.world',
            'ground_robot_forest': 'ground_robot_forest.world',
            'multiagent_testing': 'empty.world',
            'empty_wo_ground': 'empty_wo_ground.world',
            'empty': 'empty.world',
            'hospital': 'hospital.world',
            'easy_forest': 'easy_forest.world',
            'easy_high_forest': 'easy_high_forest.world',
            'medium_forest': 'medium_forest.world',
            'hard_forest': 'hard_forest.world',
            'dynamic_forest': 'dynamic_forest.world',
        }

        # Choose the world file based on the provided environment.
        world_file = world_mapping.get(env_value, 'easy_forest.world')
        world_path = PathJoinSubstitution([FindPackageShare('mighty'), 'worlds', world_file])

        use_rviz = convert_str_to_bool(LaunchConfiguration('use_rviz').perform(context))
        use_dyn_obs = convert_str_to_bool(LaunchConfiguration('use_dyn_obs').perform(context))
        use_gazebo_gui = LaunchConfiguration('use_gazebo_gui').perform(context)
        use_ground_robot = convert_str_to_bool(LaunchConfiguration('use_ground_robot').perform(context))

        # Create a rviz node - prefer source rviz config over installed one.
        # rviz_config overrides the auto-selection: multi-robot ground runs pass
        # multi_mighty_sim_ground_robot.rviz, whose per-namespace display groups
        # for NX02+ are enabled. The single-robot config leaves them off on
        # purpose — a solo run should not open with three empty topic trees.
        rviz_config_override = LaunchConfiguration('rviz_config').perform(context)
        if rviz_config_override:
            rviz_config_filename = rviz_config_override
        else:
            rviz_config_filename = 'mighty_sim_ground_robot.rviz' if use_ground_robot else 'mighty.rviz'
        # Derive workspace root from installed package path (install/mighty/share/mighty -> workspace root)
        install_share_dir = get_package_share_directory('mighty')
        ws_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(install_share_dir))))
        src_rviz_dir = os.path.join(ws_root, 'src', 'mighty', 'rviz')
        rviz_config_file = os.path.join(src_rviz_dir, rviz_config_filename)
        # Fallback to installed config if source path doesn't exist
        if not os.path.exists(rviz_config_file):
            rviz_config_file = os.path.join(install_share_dir, 'rviz', rviz_config_filename)
        # Fallback to default mighty.rviz if ground robot config doesn't exist
        if use_ground_robot and not os.path.exists(rviz_config_file):
            rviz_config_file = os.path.join(src_rviz_dir, 'mighty.rviz')
            if not os.path.exists(rviz_config_file):
                rviz_config_file = os.path.join(install_share_dir, 'rviz', 'mighty.rviz')

        rviz_node = Node(
                    package='rviz2',
                    executable='rviz2',
                    name='rviz2',
                    output='log',
                    emulate_tty=True,
                    arguments=['-d', rviz_config_file, '--ros-args', '--log-level', 'error'],
                    parameters=[{'use_sim_time': False}]
                )

        # Include Gazebo launch file
        gazebo_launch = IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                PathJoinSubstitution([FindPackageShare('gazebo_ros'), 'launch', 'gazebo.launch.py'])
            ),
            launch_arguments={'world': world_path, 'use_sim_time': 'false', 'gui': use_gazebo_gui, 'enable_gpu': 'true'}.items()
        )

        # Dynamic obstacles (optional)
        num_dyn_obstacles = LaunchConfiguration('num_dyn_obstacles').perform(context)
        dyn_x_min = LaunchConfiguration('dyn_x_min').perform(context)
        dyn_x_max = LaunchConfiguration('dyn_x_max').perform(context)
        dyn_y_min = LaunchConfiguration('dyn_y_min').perform(context)
        dyn_y_max = LaunchConfiguration('dyn_y_max').perform(context)
        dyn_z_min = LaunchConfiguration('dyn_z_min').perform(context)
        dyn_z_max = LaunchConfiguration('dyn_z_max').perform(context)
        dyn_scale_z_min = LaunchConfiguration('dyn_scale_z_min').perform(context)
        dyn_scale_z_max = LaunchConfiguration('dyn_scale_z_max').perform(context)
        dynamic_obstacles_launch = IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                PathJoinSubstitution([FindPackageShare('mighty'), 'launch', 'dyn_obstacles.launch.py'])
            ),
            launch_arguments={
                "num_obstacles": num_dyn_obstacles,
                "x_min": dyn_x_min,
                "x_max": dyn_x_max,
                "y_min": dyn_y_min,
                "y_max": dyn_y_max,
                "z_min": dyn_z_min,
                "z_max": dyn_z_max,
                "scale_z_min": dyn_scale_z_min,
                "scale_z_max": dyn_scale_z_max,
                "publish_rate_hz": "50.0",
                "seed": "0",
                "launch_forest_node": "true",
                "forest_start_delay": "2.0",
                "spawn_interval": "1.0",
            }.items()
        )

        # Return launch description
        nodes_to_start = [gazebo_launch]
        nodes_to_start.append(rviz_node) if use_rviz else None
        nodes_to_start.append(dynamic_obstacles_launch) if use_dyn_obs else None

        return nodes_to_start

    return LaunchDescription([
        env_arg,
        use_rviz_arg,
        use_gazebo_gui_arg,
        use_dyn_obs_arg,
        use_ground_robot_arg,
        rviz_config_arg,
        num_dyn_obstacles_arg,
        dyn_x_min_arg,
        dyn_x_max_arg,
        dyn_y_min_arg,
        dyn_y_max_arg,
        dyn_z_min_arg,
        dyn_z_max_arg,
        dyn_scale_z_min_arg,
        dyn_scale_z_max_arg,
        benchmark_name_arg,
        OpaqueFunction(function=launch_setup)
    ])
