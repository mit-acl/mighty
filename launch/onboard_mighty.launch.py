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
    use_obstacle_tracker_arg = DeclareLaunchArgument('use_obstacle_tracker', default_value='false', description='Flag to indicate whether to start the obstacle tracker node') # flag to indicate whether to start the obstacle tracker node
    data_file_arg = DeclareLaunchArgument('data_file', default_value='/media/kkondo/T7/dynus/tro_paper/global_planner_benchmarking/dgp.csv', description='File name to store data') # file name to store data
    global_planner_arg = DeclareLaunchArgument('global_planner', default_value='sjps', description='Global planner to use') # global planner
    use_benchmark_arg = DeclareLaunchArgument('use_benchmark', default_value='false', description='Flag to indicate whether to use the global planner benchmark') # global planner benchmark
    use_hardware_arg = DeclareLaunchArgument('use_hardware', default_value='false', description='Flag to indicate whether to use hardware or simulation') # flag to indicte if this is hardware or simulation
    use_onboard_localization_arg = DeclareLaunchArgument('use_onboard_localization', default_value='false', description='Flag to indicate whether to use t265 or vicon for localization') # flag to indicate whether to use t265 (odom) or vicon (pose & twist) for localization
    depth_camera_name_arg = DeclareLaunchArgument('depth_camera_name', default_value='d435', description='Depth camera name') # depth camera name

    # Opaque function to launch nodes
    def launch_setup(context, *args, **kwargs):

        x = LaunchConfiguration('x').perform(context)
        y = LaunchConfiguration('y').perform(context)
        z = LaunchConfiguration('z').perform(context)
        yaw = LaunchConfiguration('yaw').perform(context)
        namespace = LaunchConfiguration('namespace').perform(context)
        use_obstacle_tracker = convert_str_to_bool(LaunchConfiguration('use_obstacle_tracker').perform(context))
        data_file = LaunchConfiguration('data_file').perform(context)
        global_planner = LaunchConfiguration('global_planner').perform(context)
        use_benchmark = convert_str_to_bool(LaunchConfiguration('use_benchmark').perform(context))
        use_hardware = convert_str_to_bool(LaunchConfiguration('use_hardware').perform(context))
        use_onboard_localization = convert_str_to_bool(LaunchConfiguration('use_onboard_localization').perform(context))
        depth_camera_name = LaunchConfiguration('depth_camera_name').perform(context)

        # The path to the urdf file
        urdf_path=PathJoinSubstitution([FindPackageShare('mighty'), 'urdf', 'quadrotor.urdf.xacro'])
        parameters_path=os.path.join(get_package_share_directory('mighty'), 'config', 'mighty.yaml')

        # Get the dict of parameters from the yaml file
        with open(parameters_path, 'r') as file:
            parameters = yaml.safe_load(file)

        # Extract specific node parameters
        parameters = parameters['mighty_node']['ros__parameters']
    
        # Update parameters for benchmarking
        parameters['file_path'] = data_file
        parameters['use_benchmark'] = bool(use_benchmark)
        if use_benchmark:
            parameters['global_planner'] = global_planner
        lidar_point_could_topic = 'livox/lidar' if use_hardware else 'mid360_PointCloud2'
   
        # Create a Dynus node
        mighty_node = Node(
                    package='mighty',
                    executable='mighty',
                    name='mighty_node',
                    namespace=namespace,
                    output='screen',
                    emulate_tty=True,
                    parameters=[parameters],
                    remappings=[('lidar_cloud_in', lidar_point_could_topic),
                                ('depth_camera_cloud_in', f'{depth_camera_name}/depth/color/points')],
                    # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
                    # arguments=['--ros-args', '--log-level', 'error']
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
        
        # Create an obstacle tracker node
        obstacle_tracker_node = Node(
            package='mighty',
            executable='obstacle_tracker_node',
            namespace=namespace,
            name='obstacle_tracker_node',
            emulate_tty=True,
            parameters=[parameters],
            # prefix='xterm -e gdb -ex run --args', # gdb debugging
            output='screen',
            # remappings=[('point_cloud', lidar_point_could_topic)],
            remappings=[('point_cloud', f'{depth_camera_name}/depth/color/points')],
        )

        # Convert odom (from T265) to state
        odom_to_state_node = Node(
            package='mighty',
            executable='convert_odom_to_state',
            name='convert_odom_to_state',
            namespace=namespace,
            remappings=[
                ('odom', 'dlio/odom_node/odom'),  # Remap incoming Odometry topic
                ('state', 'state')  # Remap outgoing State topic
            ],
            emulate_tty=True,
            output='screen',
            # prefix='xterm -e gdb -ex run --args', # gdb debugging
            # arguments=['--ros-args', '--log-level', 'error']
        )

        # Convert pose and twist (from Vicon) to state
        pose_twist_to_state_node = Node(
            package='mighty',
            executable='convert_vicon_to_state',
            name='convert_vicon_to_state',
            namespace=namespace,
            remappings=[
                ('world', 'world'),  # Remap incoming PoseStamped topic
                ('twist', 'twist'),  # Remap incoming TwistStamped topic
                ('state', 'state')   # Remap outgoing State topic
            ],
            emulate_tty=True,
            output='screen',
            # prefix='xterm -e gdb -ex run --args', # gdb debugging
            # arguments=['--ros-args', '--log-level', 'error']
        )

        # When using ground robot, we don't need to send the exact state to gazebo - the state will be taken care of by wheel controllers
        # send_state_to_gazebo = False if use_ground_robot else True
        # Create a fake sim node
        fake_sim_node = Node(
                    package='mighty',
                    executable='fake_sim',
                    name='fake_sim',
                    namespace=namespace,
                    emulate_tty=True,
                    parameters=[{"start_pos": [float(x), float(y), float(z)], 
                                 "start_yaw": float(yaw),
                                 "send_state_to_gazebo": True,
                                 "visual_level": parameters['visual_level']}],
                    # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
        )

        # Return launch description
        if use_hardware and use_onboard_localization:
            nodes_to_start = [mighty_node, odom_to_state_node] # use T265 for localization
        elif use_hardware and not use_onboard_localization:
            nodes_to_start = [mighty_node, pose_twist_to_state_node] # use Vicon for localization
        else:
            nodes_to_start = [mighty_node, robot_state_publisher_node, spawn_entity_node, fake_sim_node] # simulation

        nodes_to_start.append(obstacle_tracker_node) if use_obstacle_tracker else None

        return nodes_to_start

    # Create launch description
    return LaunchDescription([
        x_arg,
        y_arg,
        z_arg,
        yaw_arg,
        namespace_arg,
        use_obstacle_tracker_arg,
        data_file_arg,
        global_planner_arg,
        use_benchmark_arg,
        use_hardware_arg,
        use_onboard_localization_arg,
        depth_camera_name_arg,
        OpaqueFunction(function=launch_setup)
    ])