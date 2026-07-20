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
from geometry_msgs.msg import PoseStamped
from std_msgs.msg import Header
import math
import numpy as np
from scipy.spatial.transform import Rotation


def stamp_to_sec(stamp):
    """Convert a builtin_interfaces/Time stamp to float seconds."""
    return stamp.sec + stamp.nanosec * 1e-9


class GoalMonitorNode(Node):
    def __init__(self):
        super().__init__('goal_monitor_node')

        # Get namespace
        self.namespace = self.get_namespace().strip('/')
        self.get_logger().info(f"Namespace: {self.namespace}")

        # Parameters
        self.declare_parameter('goal_tolerance', 1.0)  # Distance tolerance to consider goal reached
        self.goal_tolerance = self.get_parameter('goal_tolerance').value

        # --- Initial world<-ego offset estimation ---------------------------------
        # The onboard state estimator (T265/DLIO) reports the pose in an "ego" frame
        # that starts at (0, 0, 0), whereas the goals are expressed in the global
        # "world" frame (w.r.t. the mocap origin). The mocap (/<ns>/world) is only
        # available at the very beginning (it covers one side of the space and drops
        # out mid-flight), so we use it to lock a fixed rigid transform T_world_ego
        # once at startup and then apply it to the onboard state for the rest of the
        # flight. This lets us compare the (offset-corrected) onboard state against
        # the world-frame goals.
        self.declare_parameter('num_init_samples', 10)   # mocap/state pairs to average for the offset
        self.declare_parameter('init_sync_slop', 0.15)   # [s] max stamp mismatch to pair mocap with state
        self.declare_parameter('init_timeout', 5.0)      # [s] after this, pair regardless of stamp skew
        self.num_init_samples = int(self.get_parameter('num_init_samples').value)
        self.init_sync_slop = float(self.get_parameter('init_sync_slop').value)
        self.init_timeout = float(self.get_parameter('init_timeout').value)

        self.distance_check_frequency = 1.0  # Frequency to check the distance to the goal
        self.current_goal_index = 0

        # Define goal points (x, y, z) in the world frame
        # Agents are on a circle of radius 10.0 at z=3.0 and swap with their opposite partner.
        if self.namespace == 'PX01':
            self.goal_points = [[14.0,  5.0, 1.0], [ -4.0,  -2.0, 1.5]]
        elif self.namespace == 'PX02':
            self.goal_points = [[14.0,  -3.0, 1.0], [ -4.0,  2.0, 1.5]]
        # elif self.namespace == 'PX03':
            # self.goal_points = [[14.0,  0.0, 2.0], [ -4.0,  0.0, 2.0]]
        elif self.namespace == 'PX04':
            self.goal_points = [[12.5,  3.5, 1.0], [ -2.2,  -4.0, 1.5]]
        elif self.namespace == 'PX05':
            self.goal_points = [[12.5,  -1.5, 1.0], [ -2.2,  4.0, 1.5]]
        elif self.namespace == 'PX06':
            # opposite of NX02
            self.goal_points = [[4.0,  0.0, 2.0], [ -4.0,  0.0, 2.0]]
        elif self.namespace == 'PX07':
            # opposite of NX03
            self.goal_points = [[4.0,  0.0, 2.0], [ -4.0,  0.0, 2.0]]
        elif self.namespace == 'PX08':
            # opposite of NX04
            self.goal_points = [[4.0,  0.0, 2.0], [ -4.0,  0.0, 2.0]]
        else:
            self.get_logger().error(f"Unknown namespace: {self.namespace}. No goal points defined.")
            self.goal_points = [[4.0,  0.0, 2.0], [ -4.0,  0.0, 2.0]]

        # repeat pattern
        num_iterations = 2
        self.goal_points = self.goal_points * num_iterations

        # Publishers and Subscribers
        self.state_sub = self.create_subscription(State, 'state', self.state_callback, 10)
        # Mocap pose in the global/world frame; only used to lock the initial offset.
        self.world_sub = self.create_subscription(PoseStamped, 'world', self.world_callback, 10)
        self.term_goal_pub = self.create_publisher(PoseStamped, 'term_goal', 10)

        # Timer to check the distance to the current goal
        self.goal_timer = self.create_timer(self.distance_check_frequency, self.distance_check_callback)

        # Timer to publish the current goal periodically
        self.term_goal_timer = self.create_timer(1.0, self.publish_term_goal)

        # Data to store
        self.latest_state = None        # most recent onboard State (ego frame)

        # Locked world<-ego transform (applied to the onboard state once initialized)
        self.R_world_ego = None         # scipy Rotation mapping ego -> world
        self.t_world_ego = None         # np.array (3,): translation ego -> world
        self.offset_locked = False
        self._init_quats = []           # accumulated candidate offset quaternions (xyzw)
        self._init_trans = []           # accumulated candidate offset translations
        self._init_start_time = None    # time of first usable mocap+state pair [s]

        self.get_logger().info("Goal Monitor Node initialized.")

    # ---------------------------------------------------------------------------
    # Callbacks
    # ---------------------------------------------------------------------------

    def state_callback(self, msg: State):
        """Store the latest onboard state (ego frame, starts at 0, 0, 0)."""
        self.latest_state = msg

    def world_callback(self, msg: PoseStamped):
        """Lock the initial world<-ego offset from the mocap pose.

        Runs only until the offset is locked; the mocap is neither required nor
        expected once the drone flies out of coverage.
        """
        if self.offset_locked:
            return
        if self.latest_state is None:
            self.get_logger().info(
                "Received /world but no onboard /state yet; waiting to lock offset...",
                throttle_duration_sec=2.0)
            return

        now = self.get_clock().now().nanoseconds * 1e-9
        if self._init_start_time is None:
            self._init_start_time = now

        # Pair this mocap pose with the latest onboard state. Require the two
        # stamps to be close (so the correspondence is simultaneous and the
        # transform is correct); relax the requirement after init_timeout so we
        # are guaranteed to eventually lock even if the clocks are skewed.
        dt = abs(stamp_to_sec(msg.header.stamp) - stamp_to_sec(self.latest_state.header.stamp))
        slop = self.init_sync_slop if (now - self._init_start_time) < self.init_timeout else float('inf')
        if dt > slop:
            self.get_logger().warning(
                f"/world and /state stamps differ by {dt:.3f}s (> {self.init_sync_slop:.3f}s); "
                "waiting for a synchronized pair to lock the offset...",
                throttle_duration_sec=2.0)
            return

        # Body pose in the world frame (mocap) and in the ego frame (onboard state).
        q_wb = np.array([msg.pose.orientation.x, msg.pose.orientation.y,
                         msg.pose.orientation.z, msg.pose.orientation.w])
        q_eb = np.array([self.latest_state.quat.x, self.latest_state.quat.y,
                         self.latest_state.quat.z, self.latest_state.quat.w])
        if np.linalg.norm(q_wb) < 1e-6 or np.linalg.norm(q_eb) < 1e-6:
            # Degenerate quaternion (estimator not initialized yet) -> skip sample.
            return

        R_wb = Rotation.from_quat(q_wb)
        R_eb = Rotation.from_quat(q_eb)
        t_wb = np.array([msg.pose.position.x, msg.pose.position.y, msg.pose.position.z])
        t_eb = np.array([self.latest_state.pos.x, self.latest_state.pos.y, self.latest_state.pos.z])

        # T_world_ego = T_world_body * (T_ego_body)^{-1}
        R_we = R_wb * R_eb.inv()
        t_we = t_wb - R_we.apply(t_eb)

        self._init_quats.append(R_we.as_quat())
        self._init_trans.append(t_we)
        self.get_logger().info(
            f"Collected initial-offset sample {len(self._init_trans)}/{self.num_init_samples}",
            throttle_duration_sec=0.5)

        if len(self._init_trans) >= self.num_init_samples:
            self._finalize_offset()

    def _finalize_offset(self):
        """Average the collected candidate offsets and lock the transform."""
        self.R_world_ego = Rotation.from_quat(np.array(self._init_quats)).mean()
        self.t_world_ego = np.mean(np.array(self._init_trans), axis=0)
        self.offset_locked = True
        self._init_quats = []
        self._init_trans = []

        yaw_deg = math.degrees(self.R_world_ego.as_euler('xyz')[2])
        t = self.t_world_ego
        self.get_logger().info(
            f"Locked world<-ego offset from {self.num_init_samples} samples: "
            f"t=[{t[0]:.3f}, {t[1]:.3f}, {t[2]:.3f}], yaw={yaw_deg:.1f} deg. "
            "Using offset-corrected onboard state from now on.")

    # ---------------------------------------------------------------------------
    # Helpers
    # ---------------------------------------------------------------------------

    def get_estimated_world_position(self):
        """Onboard state position expressed in the world frame.

        Returns None until both an onboard state and the initial offset are
        available.
        """
        if self.latest_state is None or not self.offset_locked:
            return None
        p_ego = np.array([self.latest_state.pos.x,
                          self.latest_state.pos.y,
                          self.latest_state.pos.z])
        return self.R_world_ego.apply(p_ego) + self.t_world_ego

    # ---------------------------------------------------------------------------
    # Timers
    # ---------------------------------------------------------------------------

    def distance_check_callback(self):

        # Current position in the world frame (onboard state + locked offset)
        pos_world = self.get_estimated_world_position()
        if pos_world is None:
            if self.latest_state is None:
                self.get_logger().info("Waiting for onboard /state...", throttle_duration_sec=2.0)
            else:
                self.get_logger().info(
                    "Waiting for initial mocap /world to lock the world<-ego offset...",
                    throttle_duration_sec=2.0)
            return

        # Get the current goal point (world frame)
        goal_x, goal_y, goal_z = self.goal_points[self.current_goal_index]

        # Compute the Euclidean distance to the current goal
        distance = math.sqrt(
            (pos_world[0] - goal_x) ** 2 +
            (pos_world[1] - goal_y) ** 2 +
            (pos_world[2] - goal_z) ** 2
        )

        self.get_logger().info(
            f"[world] pos=[{pos_world[0]:.2f}, {pos_world[1]:.2f}, {pos_world[2]:.2f}] "
            f"goal=[{goal_x}, {goal_y}, {goal_z}] dist={distance:.2f}")

        # Check if the drone has reached the current goal and next goal is not out of bounds
        if distance < self.goal_tolerance and self.current_goal_index < len(self.goal_points) - 1:
            self.get_logger().info(f"Goal {self.current_goal_index} reached!")
            self.current_goal_index = self.current_goal_index + 1

    def publish_term_goal(self):
        """Publishes the current goal as a PoseStamped message (world frame)."""
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
