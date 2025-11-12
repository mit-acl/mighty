from launch import LaunchDescription
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare
from launch.substitutions import PathJoinSubstitution

def generate_launch_description():
    mighty_config_path = PathJoinSubstitution(
        [FindPackageShare('mighty'), 'config', 'mighty.yaml']
    )

    rviz_config_path = PathJoinSubstitution(
        [FindPackageShare('mighty'), 'config', 'mighty.rviz']
    )

    return LaunchDescription([
        Node(
            package='rviz2',
            executable='rviz2',
            name='rviz',
            output='screen',
            arguments=['-d', rviz_config_path],
        ),

        Node(
            package='mockamap',
            executable='mockamap_node',
            name='mockamap_node',
            output='screen',
            parameters=[{
                'seed': 1024,
                'update_freq': 1.0,
                'resolution': 0.25,
                'x_length': 50,
                'y_length': 50,
                'z_length': 5,
                'type': 1,
                'complexity': 0.025,
                'fill': 0.3,
                'fractal': 1,
                'attenuation': 0.1,
            }],
            remappings=[('/mock_map', '/voxel_map')],
        ),

        Node(
            package='mighty',
            executable='mighty_node',
            name='mighty_node',
            output='screen',
            parameters=[mighty_config_path],
        ),
    ])
