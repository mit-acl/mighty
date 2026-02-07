from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node

def generate_launch_description():
    # tolerance arg
    goal_tol_arg = DeclareLaunchArgument(
        'goal_tolerance',
        default_value='0.5',
        description='Distance tolerance to consider a goal reached'
    )
    agent_prefix_arg = DeclareLaunchArgument(
        'agent_prefix',
        default_value='NX',
        description='Agent namespace prefix: NX for drones, RR for ground robots'
    )
    use_hardware_arg = DeclareLaunchArgument(
        'use_hardware',
        default_value='false',
        description='Use hardware mode (affects goal frame_id)'
    )

    def launch_setup(context):
        prefix = LaunchConfiguration('agent_prefix').perform(context)
        use_hardware = LaunchConfiguration('use_hardware').perform(context)
        goal_tolerance = LaunchConfiguration('goal_tolerance').perform(context)

        namespaces = [f'{prefix}{i:02d}' for i in range(1, 11)]

        nodes = []
        for ns in namespaces:
            nodes.append(
                Node(
                    package='mighty',
                    executable='goal_monitor_node.py',
                    namespace=ns,
                    name='goal_monitor_node',
                    output='screen',
                    parameters=[{
                        'goal_tolerance': float(goal_tolerance),
                        'use_hardware': use_hardware.lower() in ('true', '1'),
                    }]
                )
            )
        return nodes

    return LaunchDescription([
        goal_tol_arg,
        agent_prefix_arg,
        use_hardware_arg,
        OpaqueFunction(function=launch_setup),
    ])
