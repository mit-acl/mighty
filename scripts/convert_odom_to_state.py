#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry
from std_msgs.msg import Header
from geometry_msgs.msg import Vector3, Quaternion
from dynus_interfaces.msg import State

class OdometryToStateNode(Node):
    def __init__(self):
        super().__init__('odometry_to_state_node')

        # Subscriber for nav_msgs/Odometry
        self.odometry_sub = self.create_subscription(
            Odometry, 'odom', self.callback, 10)

        # Publisher for dynus_interfaces/State
        self.state_publisher = self.create_publisher(State, 'state', 10)

    def callback(self, odom_msg):
        # Construct the State message
        state_msg = State()
        state_msg.header = Header()
        state_msg.header.stamp = odom_msg.header.stamp
        state_msg.header.frame_id = odom_msg.header.frame_id

        # Set position from Odometry
        state_msg.pos = Vector3(
            x=odom_msg.pose.pose.position.x,
            y=odom_msg.pose.pose.position.y,
            z=odom_msg.pose.pose.position.z
        )

        # Set velocity from Odometry
        state_msg.vel = Vector3(
            x=odom_msg.twist.twist.linear.x,
            y=odom_msg.twist.twist.linear.y,
            z=odom_msg.twist.twist.linear.z
        )

        # Set orientation from Odometry
        state_msg.quat = Quaternion(
            x=odom_msg.pose.pose.orientation.x,
            y=odom_msg.pose.pose.orientation.y,
            z=odom_msg.pose.pose.orientation.z,
            w=odom_msg.pose.pose.orientation.w
        )

        # Publish the State message
        self.state_publisher.publish(state_msg)
        self.get_logger().info(f'Published State: {state_msg}')

def main(args=None):
    rclpy.init(args=args)
    node = OdometryToStateNode()

    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()
