from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def discover_namespaces():
    """Query the ROS graph and return the namespaces of agents that publish a
    'state' topic (e.g. /NX01/state -> 'NX01')."""
    import rclpy
    from rclpy.node import Node as RclpyNode

    rclpy.init()
    tmp_node = RclpyNode('goal_monitor_namespace_discovery')
    try:
        # Let discovery settle so topics from already-running agents show up.
        end = tmp_node.get_clock().now() + rclpy.duration.Duration(seconds=2.0)
        namespaces = set()
        while tmp_node.get_clock().now() < end:
            rclpy.spin_once(tmp_node, timeout_sec=0.2)
            for name, _types in tmp_node.get_topic_names_and_types():
                # match '/<ns>/state' (single-level namespace)
                parts = name.split('/')
                if len(parts) == 3 and parts[0] == '' and parts[2] == 'state':
                    namespaces.add(parts[1])
        return sorted(namespaces)
    finally:
        tmp_node.destroy_node()
        rclpy.shutdown()


def launch_setup(context, *args, **kwargs):
    namespaces = discover_namespaces()

    if not namespaces:
        print('[goal_monitor.launch.py] WARNING: no /<ns>/state topics found; '
              'no goal_monitor_node launched.')

    nodes = []
    for ns in namespaces:
        nodes.append(
            Node(
                package='mighty',
                executable='goal_monitor_node.py',
                namespace=ns,
                name='goal_monitor_node',  # lives under /<ns>/goal_monitor_node
                output='screen',
                parameters=[{
                    'goal_tolerance': LaunchConfiguration('goal_tolerance')
                }]
            )
        )
    return nodes


def generate_launch_description():
    # tolerance arg
    goal_tol_arg = DeclareLaunchArgument(
        'goal_tolerance',
        default_value='1.0',
        description='Distance tolerance to consider a goal reached'
    )

    return LaunchDescription([
        goal_tol_arg,
        OpaqueFunction(function=launch_setup),
    ])
