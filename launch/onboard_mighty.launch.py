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
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import Command, PathJoinSubstitution, LaunchConfiguration
from launch_ros.substitutions import FindPackageShare
from launch_ros.parameter_descriptions import ParameterValue
import yaml
from math import radians

# Robot type constants
QUADROTOR = 'quadrotor'
RED_ROVER = 'red_rover'
STAR_ROBOT = 'star_robot'

def convert_str_to_bool(str):
    return True if (str == 'true' or str == 'True' or str == 1 or str == '1') else False

def generate_launch_description():

    # Declare launch arguments

    # initial position and yaw of the quadrotor
    x_arg = DeclareLaunchArgument('x', default_value='20.0', description='Initial x position of the quadrotor')
    y_arg = DeclareLaunchArgument('y', default_value='9.0', description='Initial y position of the quadrotor')
    z_arg = DeclareLaunchArgument('z', default_value='2.0', description='Initial z position of the quadrotor')
    yaw_arg = DeclareLaunchArgument('yaw', default_value='180', description='Initial yaw angle of the quadrotor')
    namespace_arg = DeclareLaunchArgument('namespace', default_value='NX01', description='Namespace of the nodes') # namespace
    data_file_arg = DeclareLaunchArgument('data_file', default_value='', description='File name to store data')
    global_planner_arg = DeclareLaunchArgument('global_planner', default_value='sjps', description='Global planner to use') # global planner
    use_benchmark_arg = DeclareLaunchArgument('use_benchmark', default_value='false', description='Flag to indicate whether to use the global planner benchmark') # global planner benchmark
    use_hardware_arg = DeclareLaunchArgument('use_hardware', default_value='false', description='Flag to indicate whether to use hardware or simulation') # flag to indicte if this is hardware or simulation
    publish_odom_arg  = DeclareLaunchArgument('publish_odom', default_value='true')
    odom_topic_arg    = DeclareLaunchArgument('odom_topic', default_value='visual_slam/odom')
    odom_frame_id_arg = DeclareLaunchArgument('odom_frame_id', default_value='map')
    sim_env_arg = DeclareLaunchArgument('sim_env', default_value='', description='Simulation environment: gazebo or fake_sim (empty = use mighty.yaml default)')
    config_file_arg = DeclareLaunchArgument('config_file', default_value='', description='Override planner config yaml (bare filename in mighty/config, or absolute path; empty = auto-select mighty.yaml / mighty_ground_robot.yaml)')
    use_ground_robot_arg = DeclareLaunchArgument('use_ground_robot', default_value='false', description='Enable ground robot mode (spawns p3at, uses cmd_vel control)')
    use_onboard_localization_arg = DeclareLaunchArgument('use_onboard_localization', default_value='false', description='Use onboard localization (DLIO) vs Vicon')
    depth_camera_name_arg = DeclareLaunchArgument('depth_camera_name', default_value='d435', description='Depth camera name for topic remapping')
    robot_type_arg = DeclareLaunchArgument('robot_type', default_value='quadrotor', description='Robot type: quadrotor, red_rover, star_robot')
    num_agents_arg = DeclareLaunchArgument('num_agents', default_value='10', description='Number of agents (for frame alignment subscriptions)')
    external_selector_arg = DeclareLaunchArgument(
        'external_selector', default_value='false',
        description='Set exploration.external_selector=true so mighty_node yields term_goal control (used by --use-vlm)')
    map_frame_id_arg = DeclareLaunchArgument('map_frame_id', default_value='',
        description='Override map frame ID (empty = auto from use_hardware)')
    use_frame_alignment_arg = DeclareLaunchArgument('use_frame_alignment', default_value='',
        description='Override use_frame_alignment (empty = use config default)')
    sim_frame_offset_qz_arg = DeclareLaunchArgument('sim_frame_offset_qz', default_value='',
        description='Simulated frame offset qz (empty = use config default)')
    sim_frame_offset_qw_arg = DeclareLaunchArgument('sim_frame_offset_qw', default_value='',
        description='Simulated frame offset qw (empty = use config default)')
    twist_topic_arg = DeclareLaunchArgument('twist_topic', default_value='twist',
        description='Twist topic for Vicon/mocap state converter (e.g., mocap/twist)')

    # Formation flight (empty string = use config default)
    use_formation_arg = DeclareLaunchArgument('use_formation', default_value='',
        description='Override use_formation (empty = use config default)')
    formation_weight_arg = DeclareLaunchArgument('formation_weight', default_value='',
        description='Override formation_weight (empty = use config default)')
    formation_self_offset_arg = DeclareLaunchArgument('formation_self_offset', default_value='',
        description='Comma-separated per-agent offset "dx,dy,dz" (empty = use config default)')
    formation_neighbor_ids_arg = DeclareLaunchArgument('formation_neighbor_ids', default_value='',
        description='Comma-separated neighbor IDs "1,2,3" (empty = use config default)')
    formation_neighbor_offsets_arg = DeclareLaunchArgument('formation_neighbor_offsets', default_value='',
        description='Comma-separated flat 3*N offsets "dx1,dy1,dz1,dx2,dy2,dz2,..." (empty = use config default)')

    # Need to be the same as simulartor.launch.py
    map_size_x_arg = DeclareLaunchArgument('map_size_x', default_value='20.0')
    map_size_y_arg = DeclareLaunchArgument('map_size_y', default_value='20.0')
    map_size_z_arg = DeclareLaunchArgument('map_size_z', default_value='6.0')
    odometry_topic_arg = DeclareLaunchArgument('odometry_topic', default_value='visual_slam/odom')

    # Opaque function to launch nodes
    def launch_setup(context, *args, **kwargs):

        x = LaunchConfiguration('x').perform(context)
        y = LaunchConfiguration('y').perform(context)
        z = LaunchConfiguration('z').perform(context)
        yaw = LaunchConfiguration('yaw').perform(context)
        namespace = LaunchConfiguration('namespace').perform(context)
        data_file = LaunchConfiguration('data_file').perform(context)
        global_planner = LaunchConfiguration('global_planner').perform(context)
        use_benchmark = convert_str_to_bool(LaunchConfiguration('use_benchmark').perform(context))
        use_hardware = convert_str_to_bool(LaunchConfiguration('use_hardware').perform(context))
        publish_odom = convert_str_to_bool(LaunchConfiguration('publish_odom').perform(context))
        odom_topic = LaunchConfiguration('odom_topic').perform(context)
        odom_frame_id = LaunchConfiguration('odom_frame_id').perform(context)
        base_frame_id = namespace + '/base_link'
        map_size_x = float(LaunchConfiguration('map_size_x').perform(context))
        map_size_y = float(LaunchConfiguration('map_size_y').perform(context))
        map_size_z = float(LaunchConfiguration('map_size_z').perform(context))
        odometry_topic = LaunchConfiguration('odometry_topic').perform(context)
        sim_env = LaunchConfiguration('sim_env').perform(context)
        use_ground_robot = convert_str_to_bool(LaunchConfiguration('use_ground_robot').perform(context))
        use_onboard_localization = convert_str_to_bool(LaunchConfiguration('use_onboard_localization').perform(context))
        depth_camera_name = LaunchConfiguration('depth_camera_name').perform(context)
        robot_type = LaunchConfiguration('robot_type').perform(context)
        num_agents = int(LaunchConfiguration('num_agents').perform(context))
        map_frame_id_override = LaunchConfiguration('map_frame_id').perform(context)
        use_frame_alignment_str = LaunchConfiguration('use_frame_alignment').perform(context)
        sim_frame_offset_qz_str = LaunchConfiguration('sim_frame_offset_qz').perform(context)
        sim_frame_offset_qw_str = LaunchConfiguration('sim_frame_offset_qw').perform(context)

        # The path to the urdf file - select based on robot type
        urdf_filename = 'p3at.urdf.xacro' if use_ground_robot else 'quadrotor.urdf.xacro'
        urdf_path=PathJoinSubstitution([FindPackageShare('mighty'), 'urdf', urdf_filename])
        # Planner config: explicit override wins (bare filename resolved against
        # mighty/config, or an absolute path used as-is); else auto-select by robot type.
        config_file_override = LaunchConfiguration('config_file').perform(context)
        if config_file_override:
            parameters_path = config_file_override if os.path.isabs(config_file_override) \
                else os.path.join(get_package_share_directory('mighty'), 'config', config_file_override)
        else:
            config_filename = 'mighty_ground_robot.yaml' if use_ground_robot else 'mighty.yaml'
            parameters_path = os.path.join(get_package_share_directory('mighty'), 'config', config_filename)

        # Get the dict of parameters from the yaml file
        with open(parameters_path, 'r') as file:
            parameters = yaml.safe_load(file)

        # Extract specific node parameters
        parameters = parameters['mighty_node']['ros__parameters']

        # Override sim_env if provided via launch argument
        if sim_env:
            parameters['sim_env'] = sim_env

        # Override vehicle_type if using ground robot
        if use_ground_robot:
            parameters['vehicle_type'] = 'ground_robot'

        # Override with HW config if using hardware
        if use_hardware:
            if robot_type in [RED_ROVER, STAR_ROBOT]:
                hw_config_filename = 'hw_mighty_ground_robot.yaml'
            else:  # quadrotor
                hw_config_filename = 'hw_mighty.yaml'
            hw_parameters_path = os.path.join(get_package_share_directory('mighty'), 'config', hw_config_filename)
            with open(hw_parameters_path, 'r') as f:
                hw_params = yaml.safe_load(f)['mighty_node']['ros__parameters']
            parameters.update(hw_params)

        # Update parameters for benchmarking
        parameters['file_path'] = data_file
        parameters['use_benchmark'] = bool(use_benchmark)
        if use_benchmark:
            parameters['global_planner'] = global_planner
   
        # Map frame id: hardware uses per-agent map frame, simulation uses global "map"
        map_frame_id = map_frame_id_override if map_frame_id_override else (f'{namespace}/map' if use_hardware else 'map')
        parameters['map_frame_id'] = map_frame_id
        parameters['num_agents'] = num_agents

        external_selector_str = LaunchConfiguration('external_selector').perform(context)
        if external_selector_str:
            parameters['exploration.external_selector'] = convert_str_to_bool(external_selector_str)
        if use_frame_alignment_str:
            parameters['use_frame_alignment'] = convert_str_to_bool(use_frame_alignment_str)
        if sim_frame_offset_qz_str:
            parameters['sim_frame_offset_qz'] = float(sim_frame_offset_qz_str)
        if sim_frame_offset_qw_str:
            parameters['sim_frame_offset_qw'] = float(sim_frame_offset_qw_str)

        # Formation flight overrides (empty launch arg means "use config default")
        use_formation_str = LaunchConfiguration('use_formation').perform(context)
        formation_weight_str = LaunchConfiguration('formation_weight').perform(context)
        formation_self_offset_str = LaunchConfiguration('formation_self_offset').perform(context)
        formation_neighbor_ids_str = LaunchConfiguration('formation_neighbor_ids').perform(context)
        formation_neighbor_offsets_str = LaunchConfiguration('formation_neighbor_offsets').perform(context)
        if use_formation_str:
            parameters['use_formation'] = convert_str_to_bool(use_formation_str)
        if formation_weight_str:
            parameters['formation_weight'] = float(formation_weight_str)
        if formation_self_offset_str:
            parameters['formation_self_offset'] = [float(v) for v in formation_self_offset_str.split(',')]
        if formation_neighbor_ids_str:
            # rclcpp's IntegerArray parameter requires int64, and empty strings
            # must become an empty list (not [0]).
            parameters['formation_neighbor_ids'] = [int(v) for v in formation_neighbor_ids_str.split(',') if v != '']
        if formation_neighbor_offsets_str:
            parameters['formation_neighbor_offsets'] = [float(v) for v in formation_neighbor_offsets_str.split(',') if v != '']

        # rclpy can't infer the type of empty lists in dict-form parameters
        # (raises "got '()' of type 'tuple'"). Drop them so the C++ node falls
        # back to its declared default.
        parameters = {k: v for k, v in parameters.items() if not (isinstance(v, list) and len(v) == 0)}

        # Lidar topic remapping for hardware vs simulation
        lidar_point_cloud_topic = 'livox/lidar' if use_hardware else 'mid360_PointCloud2'

        # Create a Dynus node
        mighty_node = Node(
                    package='mighty',
                    executable='mighty',
                    name='mighty_node',
                    namespace=namespace,
                    output='screen',
                    emulate_tty=True,
                    parameters=[parameters],
                    remappings=[('lidar_cloud_in', lidar_point_cloud_topic),
                                ('depth_camera_cloud_in', f'{depth_camera_name}/depth/color/points')],
                    arguments=['--ros-args', '--log-level', 'error'],
                    # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
        )

        # Robot state publisher node
        robot_state_publisher_node = Node(
            package='robot_state_publisher',
            executable='robot_state_publisher',
            name='robot_state_publisher',
            output='screen',
            namespace=namespace,
            parameters=[{
                'robot_description': ParameterValue(Command(['xacro ', urdf_path, ' namespace:=', namespace, ' d435_range_max_depth:=', str(parameters['depth_camera_depth_max'])]), value_type=str),
                'use_sim_time': False,
                'frame_prefix': namespace + '/',
            }],
            arguments=['--ros-args', '--log-level', 'error']
        )

        # Spawn entity node for Gazebo
        # Get the start position and yaw from the parameters
        yaw = str(radians(float(yaw)))
        spawn_entity_node = Node(
            package='gazebo_ros',
            executable='spawn_entity.py',
            name='spawn_entity',
            namespace=namespace,
            parameters=[{
                'use_sim_time': False,
            }],
            arguments=['-topic', 'robot_description', '-entity', namespace, '-x', x, '-y', y, '-z', z, '-Y', yaw, '--ros-args', '--log-level', 'error'],
        )
        
        # Convert pose and twist (from Vicon/mocap) to state
        twist_topic = LaunchConfiguration('twist_topic').perform(context)
        pose_twist_to_state_node = Node(
            package='mighty',
            executable='convert_vicon_to_state',
            name='convert_vicon_to_state',
            namespace=namespace,
            remappings=[
                ('world', 'world'),        # Remap incoming PoseStamped topic
                ('twist', twist_topic),     # Remap incoming TwistStamped topic (e.g., mocap/twist)
                ('state', 'state')          # Remap outgoing State topic
            ],
            emulate_tty=True,
            output='screen',
            arguments=['--ros-args', '--log-level', 'error'],
        )

        # MPC controller for ground robot (subscribes to mpc_waypoints, publishes cmd_vel).
        # Only ground robots (sim) or red_rover/star_robot (HW) use it. Guard the package
        # lookup so UAV sims don't require the 'mpc' package to be installed (e.g. Docker).
        need_mpc = (use_ground_robot and not use_hardware) or (use_hardware and robot_type in [RED_ROVER, STAR_ROBOT])
        mpc_node = None
        if need_mpc:
            mpc_params_filename = 'mpc.yaml' if use_hardware else 'mpc_sim.yaml'
            mpc_params_file = os.path.join(get_package_share_directory('mpc'), 'config', mpc_params_filename)
            mpc_cmd_vel_topic = 'cmd_vel_auto' if use_hardware else 'cmd_vel'
            mpc_overrides = {'cmd_vel_topic': mpc_cmd_vel_topic}
            # In sim, fake_sim publishes a single global "map" frame (not per-agent).
            # Prefix '/' tells mpc_node to treat the frame as absolute (skip namespace prefix).
            if not use_hardware:
                mpc_overrides['tracking_frame'] = '/' + map_frame_id
            # On HW with mocap, pose comes from the Vicon "world" topic (not DLIO)
            if use_hardware and not use_onboard_localization:
                mpc_overrides['pose_topic'] = 'world'
            mpc_node = Node(
                package='mpc',
                executable='mpc_node',
                name='mpc',
                namespace=namespace,
                output='screen',
                parameters=[mpc_params_file, mpc_overrides],
            )

        # Create a fake sim node
        fake_sim_node = Node(
                    package='mighty',
                    executable='fake_sim',
                    name='fake_sim',
                    namespace=namespace,
                    emulate_tty=True,
                    parameters=[{"start_pos": [float(x), float(y), float(z)],
                                 "start_yaw": float(yaw),
                                 "send_state_to_gazebo": parameters['sim_env'] == 'gazebo',
                                 "publish_tf": True,
                                 "publish_state": True,
                                 "use_ground_robot": use_ground_robot,
                                 "publish_odom": publish_odom,
                                 "odom_topic": odom_topic,
                                 "odom_frame_id": odom_frame_id,
                                 "base_frame_id": base_frame_id,
                                 "map_frame_id": map_frame_id,
                                 "visual_level": 1}],
                    output='screen',
        )

        camera_file = os.path.join( 
        get_package_share_directory('local_sensing'), 
        'config', 
        'camera.yaml' 
        )

        pcl_render_node = Node(
            package='local_sensing',
            executable='pcl_render_node',
            namespace=namespace,
            name='pcl_render_node',
            output='screen',
            parameters=[
                {'sensing_horizon': 5.0},
                {'sensing_rate': 30.0},
                {'estimation_rate': 30.0},
                {'map/x_size': map_size_x},
                {'map/y_size': map_size_y},
                {'map/z_size': map_size_z},
                {'use_sphere_sensing': bool(parameters.get('use_sphere_sensing', False))},
                {'sphere_sensing_radius': float(parameters.get('sphere_sensing_radius', 5.0))},
                camera_file
            ],
            remappings=[
                ('global_map', '/map_generator/global_cloud'),
                ('odometry', odometry_topic),
                ('depth', 'pcl_render_node/depth')
            ],
        )

        # HW: Odom to state (DLIO remapping)
        hw_odom_to_state_node = Node(
            package='mighty', executable='convert_odom_to_state',
            name='convert_odom_to_state', namespace=namespace,
            remappings=[('odom', 'dlio/odom_node/odom'), ('state', 'state')],
            output='screen', emulate_tty=True)

        # HW: Static TF (map->odom identity, for robots using external localization)
        static_tf_node = Node(
            package='tf2_ros', executable='static_transform_publisher',
            name='static_tf_map_to_odom', output='screen',
            arguments=['0','0','0','0','0','0','1', f'{namespace}/map', f'{namespace}/odom'])
        
        # Return launch description
        nodes_to_start = [mighty_node]
        if use_hardware:
            if use_onboard_localization:
                if robot_type == QUADROTOR:
                    nodes_to_start.append(hw_odom_to_state_node)
                elif robot_type in [STAR_ROBOT, RED_ROVER]:
                    nodes_to_start.extend([hw_odom_to_state_node, mpc_node]) #static_tf_node
            else:
                nodes_to_start.append(pose_twist_to_state_node)  # Vicon
                if robot_type in [STAR_ROBOT, RED_ROVER]:
                    nodes_to_start.append(mpc_node)
        else:
            # === EXISTING SIM CODE — COMPLETELY UNCHANGED ===
            nodes_to_start.append(pose_twist_to_state_node) if use_hardware else None
            nodes_to_start.append(fake_sim_node) if not use_hardware else None
            nodes_to_start.append(robot_state_publisher_node) if parameters['sim_env'] == 'gazebo' else None
            nodes_to_start.append(spawn_entity_node) if parameters['sim_env'] == 'gazebo' else None
            if use_ground_robot:
                nodes_to_start.append(mpc_node)
            nodes_to_start.append(pcl_render_node) if parameters['sim_env'] == 'fake_sim' else None

        return nodes_to_start

    # Create launch description
    return LaunchDescription([
        x_arg,
        y_arg,
        z_arg,
        yaw_arg,
        namespace_arg,
        data_file_arg,
        global_planner_arg,
        use_benchmark_arg,
        use_hardware_arg,
        publish_odom_arg,
        odom_topic_arg,
        odom_frame_id_arg,
        map_size_x_arg,
        map_size_y_arg,
        map_size_z_arg,
        odometry_topic_arg,
        sim_env_arg,
        config_file_arg,
        use_ground_robot_arg,
        use_onboard_localization_arg,
        depth_camera_name_arg,
        robot_type_arg,
        num_agents_arg,
        external_selector_arg,
        map_frame_id_arg,
        use_frame_alignment_arg,
        sim_frame_offset_qz_arg,
        sim_frame_offset_qw_arg,
        twist_topic_arg,
        use_formation_arg,
        formation_weight_arg,
        formation_self_offset_arg,
        formation_neighbor_ids_arg,
        formation_neighbor_offsets_arg,
        OpaqueFunction(function=launch_setup)
    ])