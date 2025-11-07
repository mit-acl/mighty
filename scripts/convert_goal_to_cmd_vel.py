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
from geometry_msgs.msg import Twist
from dynus_interfaces.msg import Goal, State
from tf_transformations import euler_from_quaternion
import numpy as np
from numpy import linalg as LA
import math

class GoalToCmdVel(Node):

    def __init__(self):

        super().__init__('goal_to_cmd_vel')

        # Initialize state
        self.state = State()
        self.state.pos.x = self.declare_parameter('x', 0.0).value
        self.state.pos.y = self.declare_parameter('y', 0.0).value
        self.state.pos.z = self.declare_parameter('z', 0.0).value
        self.state.quat.x = 0.0
        self.state.quat.y = 0.0
        self.state.quat.z = 0.0
        self.state.quat.w = 1.0
        self.current_yaw = 0.0

        # Initialize goal
        self.goal = Goal()
        self.goal.p.x = 0.0
        self.goal.p.y = 0.0
        self.goal.p.z = 0.0
        self.goal.v.x = 0.0
        self.goal.v.y = 0.0
        self.goal.v.z = 0.0
        self.goal.a.x = 0.0
        self.goal.a.y = 0.0
        self.goal.a.z = 0.0

        # Publishers
        self.pub_cmd_vel = self.create_publisher(Twist, 'cmd_vel_auto', 10)

        # Subscribers
        self.sub_goal = self.create_subscription(Goal, 'goal', self.goal_cb, 10)
        self.sub_state = self.create_subscription(State, 'state', self.state_cb, 10)

        # Timers
        self.create_timer(0.1, self.cmd_vel_cb)

        # Control parameters
        self.kv = 2.0
        self.kdist = 3.5
        self.kw = 5.0
        self.kyaw = 5.0
        self.kalpha = 1.5

        # Flags
        self.state_initialized = False
        self.goal_initialized = False

    def state_cb(self, msg):
        
        self.state = msg
        (yaw, _, _) = euler_from_quaternion(( msg.quat.x, msg.quat.y, msg.quat.z, msg.quat.w ), "szyx")
        self.current_yaw = yaw
        self.state_initialized = True

    def goal_cb(self, msg):
        self.goal = msg
        self.goal_initialized = True

    def cmd_vel_cb(self):

        if not self.state_initialized or not self.goal_initialized:
            return

        # Initialize the twist message
        twist = Twist()

        # Get the goal position, velocity, and acceleration
        x = self.goal.p.x
        y = self.goal.p.y
        xd = self.goal.v.x
        yd = self.goal.v.y
        xd2 = self.goal.a.x
        yd2 = self.goal.a.y

        # Calculate the desired velocity and acceleration
        v_desired = math.sqrt(xd**2 + yd**2)
        alpha = self.current_yaw - math.atan2(y - self.state.pos.y, x - self.state.pos.x)
        alpha = self.wrap_pi(alpha)
        forward = 1 if -math.pi / 2.0 < alpha <= math.pi / 2.0 else -1

        dist_error = forward * math.sqrt((x - self.state.pos.x)**2 + (y - self.state.pos.y)**2)

        if abs(dist_error) < 0.03:
            alpha = 0

        vel_norm = LA.norm([self.goal.v.x, self.goal.v.y, self.goal.v.z])

        if abs(dist_error) < 0.10 and vel_norm < 0.05:
            yaw_error = self.wrap_pi(self.current_yaw - self.goal.yaw)
            twist.linear.x = 0.0
            twist.angular.z = -self.kyaw * yaw_error
        else:
            numerator = xd * yd2 - yd * xd2
            denominator = xd**2 + yd**2
            w_desired = numerator / denominator if denominator > 0.01 else 0.0

            desired_yaw = math.atan2(yd, xd)
            yaw_error = self.wrap_pi(self.current_yaw - desired_yaw)

            twist.linear.x = self.kv * v_desired + self.kdist * dist_error
            twist.angular.z = self.kw * w_desired - self.kyaw * yaw_error - self.kalpha * alpha

        self.pub_cmd_vel.publish(twist)

    def wrap_pi(self, x):
        x = (x + math.pi) % (2 * math.pi)
        return x - math.pi if x >= 0 else x + math.pi

def main(args=None):
    rclpy.init(args=args)
    node = GoalToCmdVel()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
