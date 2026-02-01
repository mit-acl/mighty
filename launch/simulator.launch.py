import os
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from launch_ros.actions import ComposableNodeContainer
import launch_ros.actions
import launch_ros.descriptions
from ament_index_python.packages import get_package_share_directory
from launch.conditions import IfCondition, UnlessCondition, LaunchConfigurationEquals, LaunchConfigurationNotEquals
# from launch.substitutions import EqualsSubstitution, NotEqualsSubstitution
from launch.substitutions import PythonExpression

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
    
    # 地图属性以及是否使用动力学仿真
    use_mockamap = LaunchConfiguration('use_mockamap', default=False) # map_generator or mockamap 
    use_mockamap_arg = DeclareLaunchArgument('use_mockamap', default_value=use_mockamap, description='Choose map type, map_generator or mockamap')
    use_dynamic = LaunchConfiguration('use_dynamic', default=True)  
    use_dynamic_arg = DeclareLaunchArgument('use_dynamic', default_value=use_dynamic, description='Use Drone Simulation Considering Dynamics or Not')
    rviz_config_path = os.path.join(get_package_share_directory('mighty'), 'rviz', 'mighty.rviz')

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
        condition = UnlessCondition(use_mockamap)
    )

    # mockamap_node = Node(
    #     package='mockamap',
    #     executable='mockamap_node',
    #     name='mockamap_node',
    #     output='screen',
    #     remappings=[
    #         ('/mock_map', '/map_generator/global_cloud')
    #     ],
    #     parameters=[
    #         {'seed': 127},
    #         {'update_freq': 0.5},
    #         {'resolution': 0.1},
    #         {'x_length': PythonExpression(['int(', map_size_x_, ')'])},
    #         {'y_length': PythonExpression(['int(', map_size_y_, ')'])},
    #         {'z_length': PythonExpression(['int(', map_size_z_, ')'])},
    #         {'type': 1},
    #         {'complexity': 0.05},
    #         {'fill': 0.12},
    #         {'fractal': 1},
    #         {'attenuation': 0.1}
    #     ],
    #     condition = IfCondition(use_mockamap)
    # )

    rviz_node = launch_ros.actions.Node(
            package='rviz2', executable='rviz2', output='screen',
            arguments=['--display-config', rviz_config_path])

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

    ld.add_action(random_forest_node)
    # ld.add_action(mockamap_node)
    ld.add_action(rviz_node)

    return ld