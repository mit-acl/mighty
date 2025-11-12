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
from geometry_msgs.msg import PoseStamped, TwistStamped, Vector3, Quaternion
from std_msgs.msg import Header
from dynus_interfaces.msg import State
from message_filters import Subscriber, ApproximateTimeSynchronizer

class PoseTwistToStateNode(Node):
    def __init__(self):
        super().__init__('pose_twist_to_state_node')

        # Create subscribers using message_filters
        self.pose_sub = Subscriber(self, PoseStamped, 'world')
        self.twist_sub = Subscriber(self, TwistStamped, 'twist')

        # Synchronize pose and twist topics
        self.sync = ApproximateTimeSynchronizer(
            [self.pose_sub, self.twist_sub], queue_size=10, slop=0.1)
        self.sync.registerCallback(self.callback)

        # Publisher
        self.state_publisher = self.create_publisher(State, 'state', 10)

    def callback(self, pose_msg, twist_msg):
        # Construct the State message
        state_msg = State()
        state_msg.header = Header()
        state_msg.header.stamp = self.get_clock().now().to_msg()
        state_msg.header.frame_id = pose_msg.header.frame_id

        # Set position from PoseStamped
        state_msg.pos = Vector3(
            x=pose_msg.pose.position.x,
            y=pose_msg.pose.position.y,
            z=pose_msg.pose.position.z
        )

        # Set velocity from TwistStamped
        state_msg.vel = Vector3(
            x=twist_msg.twist.linear.x,
            y=twist_msg.twist.linear.y,
            z=twist_msg.twist.linear.z
        )

        # Set orientation from PoseStamped
        state_msg.quat = Quaternion(
            x=pose_msg.pose.orientation.x,
            y=pose_msg.pose.orientation.y,
            z=pose_msg.pose.orientation.z,
            w=pose_msg.pose.orientation.w
        )

        # Publish the State message
        self.state_publisher.publish(state_msg)
        self.get_logger().info(f'Published State: {state_msg}')

def main(args=None):
    rclpy.init(args=args)
    node = PoseTwistToStateNode()

    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()
