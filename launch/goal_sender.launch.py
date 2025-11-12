#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import os
import yaml
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration


def generate_launch_description():

    # Parameters
    parameters_path = os.path.join(
        get_package_share_directory('mighty'),
        'config',
        'mighty.yaml'
        )

    # Get the dict of parameters from the yaml file
    with open(parameters_path, 'r') as file:
        parameters = yaml.safe_load(file)

    # Extract specific node parameters
    parameters = parameters['mighty_node']['ros__parameters']
    default_goal_z = parameters['default_goal_z']

    return LaunchDescription([

        # Declare the launch arguments
        DeclareLaunchArgument(
            'list_agents',
            default_value="['NX01', 'NX02', 'NX03', 'NX04', 'NX05']",
            # default_value="['NX01', 'NX02']",
            description='List of agent names'
        ),
        DeclareLaunchArgument(
            'list_goals',
            # default_value="['[50.0, 20.0, 0.0]', '[50.0, 10.0, 0.0]', '[50.0, 0.0, 0.0]', '[50.0, -10.0, 0.0]', '[50.0, -20.0, 0.0]']",
            default_value="['[50.0, -30.0, 0.0]', '[50.0, -15.0, 0.0]', '[50.0, 0.0, 0.0]', '[50.0, 15.0, 0.0]', '[50.0, 30.0, 0.0]']",
            # default_value="['[35.0, 20.0, 0.0]', '[35.0, 10.0, 0.0]', '[35.0, 0.0, 0.0]', '[35.0, -10.0, 0.0]', '[35.0, -20.0, 0.0]']",
            # default_value="['[5.0, 0.0, 2.0]', '[-5.0, 0.0, 2.0]']",
            description='List of goal coordinates'
        ),
        DeclareLaunchArgument(
            'default_goal_z',
            default_value=str(default_goal_z),
            description='Default Z value for the goal'
        ),

        # Launch the GoalSender node
        Node(
            package='mighty',
            executable='goal_sender.py',
            name='goal_sender',
            output='screen',
            parameters=[{
                'list_agents': LaunchConfiguration('list_agents'),
                'list_goals': LaunchConfiguration('list_goals'),
                'default_goal_z': LaunchConfiguration('default_goal_z'),
            }]
        ),
    ])
