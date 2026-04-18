# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

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
    num_agents_arg = DeclareLaunchArgument(
        'num_agents',
        default_value='10',
        description='Number of agents (must match the actual number launched)'
    )
    radius_arg = DeclareLaunchArgument(
        'radius',
        default_value='10.0',
        description='Circle formation radius'
    )
    use_ground_robot_arg = DeclareLaunchArgument(
        'use_ground_robot',
        default_value='false',
        description='Ground robot mode (sets goal z to 0.2 instead of 1.0)'
    )
    angle_offset_arg = DeclareLaunchArgument(
        'angle_offset',
        default_value='0.0',
        description='Angular offset in radians for circle positions (e.g. 0.7854 for 45 deg)'
    )

    def launch_setup(context):
        prefix = LaunchConfiguration('agent_prefix').perform(context)
        use_hardware = LaunchConfiguration('use_hardware').perform(context)
        goal_tolerance = LaunchConfiguration('goal_tolerance').perform(context)
        num_agents = int(LaunchConfiguration('num_agents').perform(context))
        radius = float(LaunchConfiguration('radius').perform(context))
        use_ground_robot = LaunchConfiguration('use_ground_robot').perform(context)
        angle_offset = float(LaunchConfiguration('angle_offset').perform(context))

        namespaces = [f'{prefix}{i:02d}' for i in range(1, num_agents + 1)]

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
                        'use_ground_robot': use_ground_robot.lower() in ('true', '1'),
                        'num_agents': num_agents,
                        'radius': radius,
                        'angle_offset': angle_offset,
                    }]
                )
            )
        return nodes

    return LaunchDescription([
        goal_tol_arg,
        agent_prefix_arg,
        use_hardware_arg,
        use_ground_robot_arg,
        num_agents_arg,
        radius_arg,
        angle_offset_arg,
        OpaqueFunction(function=launch_setup),
    ])
