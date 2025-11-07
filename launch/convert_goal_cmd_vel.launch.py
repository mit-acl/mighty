#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration

def generate_launch_description():
    # Namespace argument
    namespace_arg = DeclareLaunchArgument(
        'namespace',
        default_value='quadrotor',
        description='Namespace of the nodes'
    )

    # Parameters for the GoalToCmdVel node
    param_x = DeclareLaunchArgument(
        'x', default_value='0.0', description='Initial x position'
    )
    param_y = DeclareLaunchArgument(
        'y', default_value='0.0', description='Initial y position'
    )
    param_z = DeclareLaunchArgument(
        'z', default_value='0.0', description='Initial z position'
    )

    def launch_setup(context, *args, **kwargs):
        namespace = LaunchConfiguration('namespace').perform(context)

        # GoalToCmdVel Node
        node = Node(
            package='mighty',  # Replace with your actual package name
            executable='convert_goal_to_cmd_vel',  # The name of your compiled executable
            name='convert_goal_to_cmd_vel',
            namespace=namespace,
            output='screen',
            parameters=[
                {'x': LaunchConfiguration('x')},
                {'y': LaunchConfiguration('y')},
                {'z': LaunchConfiguration('z')}
            ],
            remappings=[
                ('goal', 'goal'),  # Incoming goal topic
                ('state', 'state'),  # Incoming state topic
                ('cmd_vel', 'cmd_vel_auto')  # Outgoing command velocity topic
            ]
        )

        return [node]

    return LaunchDescription([
        namespace_arg,
        param_x,
        param_y,
        param_z,
        OpaqueFunction(function=launch_setup)
    ])
