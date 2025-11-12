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
from dynus_interfaces.msg import State
from geometry_msgs.msg import PoseStamped, Vector3
from std_msgs.msg import Header
import math

class GoalMonitorNode(Node):
    def __init__(self):
        super().__init__('goal_monitor_node')

        # Get namespace
        self.namespace = self.get_namespace().strip('/')
        self.get_logger().info(f"Namespace: {self.namespace}")

        # Parameters
        self.declare_parameter('goal_tolerance', 1.0)  # Distance tolerance to consider goal reached
        self.goal_tolerance = self.get_parameter('goal_tolerance').value
        self.distance_check_frequency = 1.0  # Frequency to check the distance to the goal
        self.current_goal_index = 0

        # Define goal points (x, y, z) in the world frame
        if self.namespace == 'NX01':
            self.goal_points = [[105.0, 15.0, 3.0]]
            # for y = 10, …, -10, append [0,y,3] then [105,y,3]
            for y in range(10, -15, -5):
                self.goal_points.append([-5.0, float(y), 3.0])
                self.goal_points.append([105.0, float(y), 3.0])
        elif self.namespace == 'NX02':
            self.goal_points = [[13.0, -3.0, 3.0], [-3.0, 1.5, 3.0]]
        elif self.namespace == 'NX03':
            self.goal_points = [[13.0, 3.0, 3.0], [-3.0, -1.5, 3.0]]
        elif self.namespace == 'NX04':
            self.goal_points = [[15.0, 1.5, 3.0], [0.0, -3.0, 3.0]]
        else:
            self.get_logger().error(f"Unknown namespace: {self.namespace}. No goal points defined.")
            self.goal_points = [[0.0, 0.0, 0.0]]  # Default goal point if namespace is unknown

        # Publishers and Subscribers
        self.state_sub = self.create_subscription(State, 'state', self.state_callback, 10)
        self.term_goal_pub = self.create_publisher(PoseStamped, 'term_goal', 10)

        # Timer to check the distance to the current goal
        self.goal_timer = self.create_timer(self.distance_check_frequency, self.distance_check_callback)

        # Timer to publish the current goal periodically
        self.term_goal_timer = self.create_timer(1.0, self.publish_term_goal)

        # Data to store
        self.current_position = Vector3()

        self.get_logger().info("Goal Monitor Node initialized.")

    def state_callback(self, msg: State):

        """Callback for monitoring the drone's position."""
        self.current_position = msg.pos

    def distance_check_callback(self):

        # Get the current goal point
        goal_x, goal_y, goal_z = self.goal_points[self.current_goal_index]

        # Compute the Euclidean distance to the current goal
        distance = math.sqrt(
            (self.current_position.x - goal_x) ** 2 +
            (self.current_position.y - goal_y) ** 2 +
            (self.current_position.z - goal_z) ** 2
        )

        self.get_logger().info(f"Distance to goal: {distance:.2f}")

        # Check if the drone has reached the current goal and next goal is not out of bounds
        if distance < self.goal_tolerance and self.current_goal_index < len(self.goal_points) - 1:
            self.get_logger().info(f"Goal {self.current_goal_index} reached!")
            self.current_goal_index = self.current_goal_index + 1

    def publish_term_goal(self):
        """Publishes the current goal as a PoseStamped message."""
        goal_x, goal_y, goal_z = self.goal_points[self.current_goal_index]

        # Create PoseStamped message
        term_goal = PoseStamped()
        term_goal.header = Header()
        term_goal.header.stamp = self.get_clock().now().to_msg()
        term_goal.header.frame_id = "world"

        term_goal.pose.position.x = goal_x
        term_goal.pose.position.y = goal_y
        term_goal.pose.position.z = goal_z

        term_goal.pose.orientation.x = 0.0
        term_goal.pose.orientation.y = 0.0
        term_goal.pose.orientation.z = 0.0
        term_goal.pose.orientation.w = 1.0  # Identity quaternion

        # Publish the term goal
        self.term_goal_pub.publish(term_goal)
        self.get_logger().info(f"Published term goal: [{goal_x}, {goal_y}, {goal_z}]")

def main(args=None):
    rclpy.init(args=args)
    node = GoalMonitorNode()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
