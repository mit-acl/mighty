from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node

def generate_launch_description():
    # keep your tolerance arg
    goal_tol_arg = DeclareLaunchArgument(
        'goal_tolerance',
        default_value='0.5',
        description='Distance tolerance to consider a goal reached'
    )

    # list out the four namespaces you want
    namespaces = ['NX01']

    # for each namespace, create one Node
    nodes = []
    for ns in namespaces:
        nodes.append(
            Node(
                package='mighty',
                executable='goal_monitor_node.py',
                namespace=ns,
                name='goal_monitor_node',  # this will live under /<ns>/goal_monitor_node
                output='screen',
                parameters=[{
                    'goal_tolerance': LaunchConfiguration('goal_tolerance')
                }]
            )
        )

    return LaunchDescription([
        goal_tol_arg,
        *nodes
    ])
