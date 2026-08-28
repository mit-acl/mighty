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
    mode_arg = DeclareLaunchArgument(
        'mode',
        default_value='legacy',
        description='legacy: one goal_monitor_node per agent (endless swap); '
                    'exchange/random: single central coordinator with synchronized rounds'
    )
    num_cycles_arg = DeclareLaunchArgument(
        'num_cycles',
        default_value='0',
        description='exchange/random: number of cycles (exchange: out + back = 1 cycle; '
                    'random: 1 round = 1 cycle); <= 0 runs forever'
    )
    clockwise_arg = DeclareLaunchArgument(
        'clockwise',
        default_value='true',
        description='exchange/random: number agents clockwise along the circle (JR '
                    'hardware convention); false matches the CCW sim spawn convention'
    )
    seed_arg = DeclareLaunchArgument(
        'seed',
        default_value='-1',
        description='random mode RNG seed; < 0 for nondeterministic'
    )
    state_source_arg = DeclareLaunchArgument(
        'state_source',
        default_value='state',
        description='exchange/random: where agent positions come from — state: '
                    'dynus_interfaces/State on /<ns>/state; world: '
                    'geometry_msgs/PoseStamped on /<ns>/world (raw mocap)'
    )
    list_agents_arg = DeclareLaunchArgument(
        'list_agents',
        default_value='',
        description='exchange/random: explicit comma-separated agent list (e.g. '
                    '"JR01,JR02,RR03") overriding agent_prefix/num_agents; '
                    'list order = clockwise slot order, first agent at angle_offset'
    )

    def launch_setup(context):
        prefix = LaunchConfiguration('agent_prefix').perform(context)
        use_hardware = LaunchConfiguration('use_hardware').perform(context)
        goal_tolerance = LaunchConfiguration('goal_tolerance').perform(context)
        num_agents = int(LaunchConfiguration('num_agents').perform(context))
        radius = float(LaunchConfiguration('radius').perform(context))
        use_ground_robot = LaunchConfiguration('use_ground_robot').perform(context)
        angle_offset = float(LaunchConfiguration('angle_offset').perform(context))
        mode = LaunchConfiguration('mode').perform(context).lower()

        # Central coordinator: one node driving all agents in synchronized rounds
        if mode in ('exchange', 'random'):
            return [
                Node(
                    package='mighty',
                    executable='goal_monitor_central_node.py',
                    name='goal_monitor_central_node',
                    output='screen',
                    parameters=[{
                        'mode': mode,
                        'agent_prefix': prefix,
                        'num_agents': num_agents,
                        'radius': radius,
                        'angle_offset': angle_offset,
                        'clockwise': LaunchConfiguration('clockwise').perform(context).lower() in ('true', '1'),
                        'num_cycles': int(LaunchConfiguration('num_cycles').perform(context)),
                        'goal_tolerance': float(goal_tolerance),
                        'use_hardware': use_hardware.lower() in ('true', '1'),
                        'use_ground_robot': use_ground_robot.lower() in ('true', '1'),
                        'seed': int(LaunchConfiguration('seed').perform(context)),
                        'state_source': LaunchConfiguration('state_source').perform(context).lower(),
                        'list_agents': LaunchConfiguration('list_agents').perform(context),
                    }]
                )
            ]

        if mode != 'legacy':
            raise RuntimeError(
                f"goal_monitor.launch.py: invalid mode '{mode}' "
                f"(expected legacy, exchange, or random)")

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
        mode_arg,
        num_cycles_arg,
        clockwise_arg,
        seed_arg,
        state_source_arg,
        list_agents_arg,
        OpaqueFunction(function=launch_setup),
    ])
