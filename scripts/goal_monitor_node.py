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


def circle_position(agent_index, num_agents, radius, z=1.0, angle_offset=0.0):
    """Compute own and opposite positions on a circle for swap goals.

    Agent i sits at angle = 2*pi*(i-1)/N + angle_offset on a circle of the
    given radius; its swap target is the diametrically opposite point
    (angle + pi).

    Returns:
        (own_x, own_y, z, opp_x, opp_y, z)
    """
    angle = 2.0 * math.pi * (agent_index - 1) / num_agents + angle_offset
    own_x = round(radius * math.cos(angle), 3)
    own_y = round(radius * math.sin(angle), 3)
    opp_x = round(radius * math.cos(angle + math.pi), 3)
    opp_y = round(radius * math.sin(angle + math.pi), 3)
    return own_x, own_y, z, opp_x, opp_y, z


class GoalMonitorNode(Node):
    def __init__(self):
        super().__init__('goal_monitor_node')

        # Get namespace
        self.namespace = self.get_namespace().strip('/')
        self.get_logger().info(f"Namespace: {self.namespace}")

        # Parameters
        self.declare_parameter('goal_tolerance', 1.0)
        self.goal_tolerance = self.get_parameter('goal_tolerance').value
        self.declare_parameter('use_hardware', False)
        self.use_hardware = self.get_parameter('use_hardware').value
        self.declare_parameter('use_ground_robot', False)
        self.use_ground_robot = self.get_parameter('use_ground_robot').value
        self.declare_parameter('odom_type', 'dlio')  # "dlio" or "mocap"
        odom_type = self.get_parameter('odom_type').value
        self.declare_parameter('goal_type', 1)  # 1 or 2 (only for mocap mode)
        goal_type = self.get_parameter('goal_type').value
        self.declare_parameter('num_agents', 10)
        num_agents = self.get_parameter('num_agents').value
        self.declare_parameter('radius', 10.0)
        radius = self.get_parameter('radius').value
        self.declare_parameter('angle_offset', 0.0)
        angle_offset = self.get_parameter('angle_offset').value
        self.distance_check_frequency = 1.0
        self.current_goal_index = 0

        z = 0.2 if self.use_ground_robot else 1.0

        # Hardware ground robot: fixed goal pairs based on odom_type
        if self.use_hardware and self.use_ground_robot:
            if odom_type == 'mocap' or odom_type == 'dlio_in_mocap':
                if goal_type == 1:
                    self.goal_points = [[3.5, 3.5, z], [-3.5, -3.5, z]]
                else:
                    self.goal_points = [[-3.5, 3.5, z], [3.5, -3.5, z]]
                self.get_logger().info(
                    f"HW ground robot mocap goals (type {goal_type}): "
                    f"{self.goal_points[0]} <-> {self.goal_points[1]}")
            else:  # dlio
                self.goal_points = [[0.0, 0.0, z], [16.0, 0.0, z]]
                self.get_logger().info(
                    f"HW ground robot DLIO goals: "
                    f"{self.goal_points[0]} <-> {self.goal_points[1]}")

        # Simulation: circle formation swap
        elif self.namespace.startswith('NX'):
            agent_index = int(self.namespace[2:])  # NX01 -> 1
            own_x, own_y, _, opp_x, opp_y, _ = circle_position(agent_index, num_agents, radius, z=z, angle_offset=angle_offset)
            self.goal_points = [[opp_x, opp_y, z], [own_x, own_y, z]]
            self.get_logger().info(
                f"Circle swap goals (N={num_agents}, R={radius}): "
                f"start ({own_x},{own_y}) <-> opposite ({opp_x},{opp_y})")

        elif self.namespace.startswith('RR'):
            agent_index = int(self.namespace[2:])  # RR01 -> 1
            own_x, own_y, _, opp_x, opp_y, _ = circle_position(agent_index, num_agents, radius, z=z, angle_offset=angle_offset)
            self.goal_points = [[opp_x, opp_y, z], [own_x, own_y, z]]
            self.get_logger().info(
                f"Circle swap goals (N={num_agents}, R={radius}): "
                f"start ({own_x},{own_y}) <-> opposite ({opp_x},{opp_y})")

        else:
            self.get_logger().error(f"Unknown namespace: {self.namespace}. No goal points defined.")
            self.goal_points = [[0.0, 0.0, 0.0]]

        # No need to repeat — we wrap around forever

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

        # Compute distance to goal (2D for ground robots, 3D for UAVs)
        if self.use_ground_robot:
            distance = math.sqrt(
                (self.current_position.x - goal_x) ** 2 +
                (self.current_position.y - goal_y) ** 2
            )
        else:
            distance = math.sqrt(
                (self.current_position.x - goal_x) ** 2 +
                (self.current_position.y - goal_y) ** 2 +
                (self.current_position.z - goal_z) ** 2
            )

        self.get_logger().info(f"Distance to goal: {distance:.2f}")

        # Check if the drone has reached the current goal, then wrap around
        if distance < self.goal_tolerance:
            self.get_logger().info(f"Goal {self.current_goal_index} reached!")
            self.current_goal_index = (self.current_goal_index + 1) % len(self.goal_points)

    def publish_term_goal(self):
        """Publishes the current goal as a PoseStamped message."""
        goal_x, goal_y, goal_z = self.goal_points[self.current_goal_index]

        # Create PoseStamped message
        term_goal = PoseStamped()
        term_goal.header = Header()
        term_goal.header.stamp = self.get_clock().now().to_msg()
        term_goal.header.frame_id = f'{self.namespace}/map' if self.use_hardware else 'map'

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
