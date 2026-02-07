#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""
Publishes a ground-truth TransformStamped at 10 Hz on
  /lidar_registration/frame_align/{ego_name}/{other_name}

Parameters:
  ego_name   (str)   – e.g. "RR01"
  other_name (str)   – e.g. "RR02"
  tx, ty, tz (float) – translation components
  qx, qy, qz, qw (float) – quaternion components
"""

import rclpy
from rclpy.node import Node
from geometry_msgs.msg import TransformStamped


class FrameAlignPublisher(Node):

    def __init__(self):
        super().__init__('frame_align_publisher')

        # Declare parameters
        self.declare_parameter('ego_name', 'RR01')
        self.declare_parameter('other_name', 'RR02')
        self.declare_parameter('tx', 0.0)
        self.declare_parameter('ty', 0.0)
        self.declare_parameter('tz', 0.0)
        self.declare_parameter('qx', 0.0)
        self.declare_parameter('qy', 0.0)
        self.declare_parameter('qz', 0.0)
        self.declare_parameter('qw', 1.0)

        # Get parameters
        ego = str(self.get_parameter('ego_name').value)
        other = str(self.get_parameter('other_name').value)
        self.tx = float(self.get_parameter('tx').value)
        self.ty = float(self.get_parameter('ty').value)
        self.tz = float(self.get_parameter('tz').value)
        self.qx = float(self.get_parameter('qx').value)
        self.qy = float(self.get_parameter('qy').value)
        self.qz = float(self.get_parameter('qz').value)
        self.qw = float(self.get_parameter('qw').value)

        # Build the topic name
        topic = f'/lidar_registration/frame_align/{ego}/{other}'
        self.get_logger().info(f'Publishing frame alignment on {topic}')

        # Build the static parts of the message
        self.msg = TransformStamped()
        self.msg.header.frame_id = f'{ego}/map'
        self.msg.child_frame_id = f'{other}/map'
        self.msg.transform.translation.x = self.tx
        self.msg.transform.translation.y = self.ty
        self.msg.transform.translation.z = self.tz
        self.msg.transform.rotation.x = self.qx
        self.msg.transform.rotation.y = self.qy
        self.msg.transform.rotation.z = self.qz
        self.msg.transform.rotation.w = self.qw

        # Publisher + timer at 10 Hz
        self.pub = self.create_publisher(TransformStamped, topic, 10)
        self.timer = self.create_timer(0.1, self.timer_callback)

    def timer_callback(self):
        self.msg.header.stamp = self.get_clock().now().to_msg()
        self.pub.publish(self.msg)


def main(args=None):
    rclpy.init(args=args)
    node = FrameAlignPublisher()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
