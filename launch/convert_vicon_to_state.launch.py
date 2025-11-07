#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import launch
from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():


    return LaunchDescription([
        Node(
            package='mighty',
            executable='convert_vicon_to_state', 
            name='convert_vicon_to_state',
            remappings=[
                ('world', 'world'),  # Remap incoming PoseStamped topic
                ('twist', 'twist'),  # Remap incoming TwistStamped topic
                ('state', 'state')   # Remap outgoing State topic
            ],
            output='screen'
        )
    ])
