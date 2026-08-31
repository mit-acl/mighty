#!/usr/bin/env python3
# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */
"""Constant SUB-FLOOR global cloud for obstacle-free fake_sim worlds.

The planner cannot run on literally nothing: MIGHTY::checkReadyToReplan()
requires an initialized map even in simulation, the map only initializes when
a sensor_point_cloud arrives, and pcl_render publishes nothing unless its
radiusSearch over the global cloud returns points (an EMPTY global cloud
would leave the FLANN index unbuilt). So an "empty world" still needs a
global cloud — just one the voxel map cannot see.

This node publishes a flat grid of points at z = -3 m. Geometry that makes it
work: |z| = 3 m is inside the agents' 5 m sensing sphere (horizontal reach
sqrt(5^2 - 3^2) = 4 m >> the 2 m grid spacing), so every agent's pcl_render
always finds points and keeps feeding the planner — while z = -3 m is below
z_min, so map_util drops every point and the planned world is obstacle-free.

Republished at 1 Hz so late-starting agents' pcl_render nodes (which latch
the first global map they hear) are all served.
"""

import struct

import rclpy
from rclpy.node import Node
from rclpy.qos import QoSProfile, ReliabilityPolicy, HistoryPolicy
from sensor_msgs.msg import PointCloud2, PointField
from std_msgs.msg import Header


class SubfloorCloudPublisher(Node):
    def __init__(self):
        super().__init__('subfloor_cloud')
        self.declare_parameter('extent_m', 14.0)
        self.declare_parameter('spacing_m', 2.0)
        self.declare_parameter('z_m', -3.0)
        self.declare_parameter('frame_id', 'map')

        extent = self.get_parameter('extent_m').value
        spacing = self.get_parameter('spacing_m').value
        z = self.get_parameter('z_m').value
        self.frame_id = self.get_parameter('frame_id').value

        pts = []
        x = -extent
        while x <= extent + 1e-6:
            y = -extent
            while y <= extent + 1e-6:
                pts.append((x, y, z))
                y += spacing
            x += spacing
        self.msg = self._make_cloud(pts)
        self.get_logger().info(
            f'Sub-floor cloud: {len(pts)} points at z={z} m over ±{extent} m '
            f'(feeds pcl_render; invisible to the voxel map via z_min)')

        qos = QoSProfile(depth=1, history=HistoryPolicy.KEEP_LAST,
                         reliability=ReliabilityPolicy.RELIABLE)
        self.pub = self.create_publisher(PointCloud2, '/map_generator/global_cloud', qos)
        self.timer = self.create_timer(1.0, self._tick)

    def _make_cloud(self, pts):
        msg = PointCloud2()
        msg.header = Header(frame_id=self.frame_id)
        msg.height = 1
        msg.width = len(pts)
        msg.fields = [
            PointField(name='x', offset=0, datatype=PointField.FLOAT32, count=1),
            PointField(name='y', offset=4, datatype=PointField.FLOAT32, count=1),
            PointField(name='z', offset=8, datatype=PointField.FLOAT32, count=1),
        ]
        msg.is_bigendian = False
        msg.point_step = 12
        msg.row_step = msg.point_step * len(pts)
        msg.is_dense = True
        msg.data = b''.join(struct.pack('<fff', *p) for p in pts)
        return msg

    def _tick(self):
        self.msg.header.stamp = self.get_clock().now().to_msg()
        self.pub.publish(self.msg)


def main():
    rclpy.init()
    node = SubfloorCloudPublisher()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()


if __name__ == '__main__':
    main()
