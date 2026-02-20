#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2024, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import rclpy
from rclpy.node import Node
from dynus_interfaces.msg import State
from geometry_msgs.msg import PoseStamped, Vector3
from std_msgs.msg import Header
import math
import os

# republish clicked 2D goal pose in rviz as dynus State type
class RepubGoalNode(Node):
    def __init__(self):
        super().__init__('repub_goal_node')

        # Get rover name
        self.namespace = self.get_namespace().strip('/')
        self.get_logger().info(f"Namespace: {self.namespace}")
        vehtype = os.getenv("VEHTYPE", default = None)
        vehnum = os.getenv("VEHNUM", default = None)
        self.rover_name = vehtype + vehnum

        # Publishers and Subscribers
        self.rviz_sub = self.create_subscription(PoseStamped, f'{self.rover_name}/term_goal_rviz', self.rviz_goal_cb, 10)
        # self.term_goal_pub = self.create_publisher(PoseStamped, f'{self.rover_name}/term_goal', 10)
        self.term_goal_pub = self.create_publisher(State, f'{self.rover_name}/term_goal', 10)

        self.get_logger().info("Rviz 2D republisher node initialized.")

    def rviz_goal_cb(self, msg: PoseStamped):

        # Get goal position
        goal_pos = msg.pose.position

        # init PoseStamped message
        # term_goal = PoseStamped()
        term_goal = State()
        term_goal.header = Header()
        term_goal.header.stamp = self.get_clock().now().to_msg()
        term_goal.header.frame_id = "world"

        # # set goal position
        # term_goal.pose.position.x = goal_pos.x
        # term_goal.pose.position.y = goal_pos.y
        # term_goal.pose.position.z = goal_pos.z

        # set goal position
        term_goal.pos.x = goal_pos.x
        term_goal.pos.y = goal_pos.y
        term_goal.pos.z = goal_pos.z
        
        # goal velocity of zero
        term_goal.vel.x = 0.
        term_goal.vel.y = 0.
        term_goal.vel.z = 0.

        # ignoring orientation
        term_goal.quat.x = 0.0
        term_goal.quat.y = 0.0
        term_goal.quat.z = 0.0
        term_goal.quat.w = 1.0

        # publish the term goal
        self.term_goal_pub.publish(term_goal)
        self.get_logger().info(f"Published term goal: [{goal_pos.x}, {goal_pos.y}, {goal_pos.z}]")

def main(args=None):
    rclpy.init(args=args)
    node = RepubGoalNode()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
