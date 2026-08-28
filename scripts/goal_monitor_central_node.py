#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""Central multi-agent goal coordinator with synchronized rounds.

Unlike goal_monitor_node.py (one independent node per agent, endless swap),
this single node watches every agent's state and only advances to the next
round of goals once ALL agents are within goal_tolerance of their current
goals. Two modes:

  exchange: every agent goes to the point diametrically opposite its start
            slot, then back home; out + back = 1 cycle, repeated num_cycles
            times.
  random:   every round, agents are re-assigned the N circle slots by a
            random derangement (no agent keeps its slot, no two agents share
            a target), repeated num_cycles rounds.

Launched via goal_monitor.launch.py with mode:=exchange or mode:=random.
"""

import math
import random

import yaml

import rclpy
from rclpy.node import Node
from rclpy.qos import QoSProfile, HistoryPolicy, ReliabilityPolicy, DurabilityPolicy
from rclpy.task import Future

from dynus_interfaces.msg import State
from geometry_msgs.msg import PoseStamped


def slot_angle(slot, num_agents, angle_offset, clockwise):
    """Angle of circle slot `slot` (0-based; slot k is agent k+1's start).

    Clockwise numbering matches the JR rover placement convention (agent 1 at
    angle_offset, agent 2 at angle_offset - 2*pi/N, ...); counterclockwise
    matches the sim spawn convention used by run_sim.py and circle_position().
    """
    sign = -1.0 if clockwise else 1.0
    return angle_offset + sign * 2.0 * math.pi * slot / num_agents


def slot_point(slot, num_agents, radius, angle_offset, clockwise):
    angle = slot_angle(slot, num_agents, angle_offset, clockwise)
    return round(radius * math.cos(angle), 3), round(radius * math.sin(angle), 3)


def antipode_point(slot, num_agents, radius, angle_offset, clockwise):
    angle = slot_angle(slot, num_agents, angle_offset, clockwise) + math.pi
    return round(radius * math.cos(angle), 3), round(radius * math.sin(angle), 3)


def random_derangement(n, rng):
    """Random permutation of range(n) with no fixed points (requires n >= 2)."""
    perm = list(range(n))
    while True:
        rng.shuffle(perm)
        if all(perm[k] != k for k in range(n)):
            return perm


class GoalMonitorCentralNode(Node):
    def __init__(self):
        super().__init__('goal_monitor_central_node')

        # Parameters
        self.declare_parameter('mode', 'exchange')
        self.mode = self.get_parameter('mode').value
        self.declare_parameter('agent_prefix', 'NX')
        prefix = self.get_parameter('agent_prefix').value
        self.declare_parameter('num_agents', 10)
        self.num_agents = self.get_parameter('num_agents').value
        self.declare_parameter('radius', 10.0)
        self.radius = self.get_parameter('radius').value
        self.declare_parameter('angle_offset', 0.0)
        self.angle_offset = self.get_parameter('angle_offset').value
        self.declare_parameter('clockwise', True)
        self.clockwise = self.get_parameter('clockwise').value
        self.declare_parameter('num_cycles', 0)
        self.num_cycles = self.get_parameter('num_cycles').value
        self.declare_parameter('goal_tolerance', 0.5)
        self.goal_tolerance = self.get_parameter('goal_tolerance').value
        self.declare_parameter('use_hardware', False)
        self.use_hardware = self.get_parameter('use_hardware').value
        self.declare_parameter('use_ground_robot', False)
        self.use_ground_robot = self.get_parameter('use_ground_robot').value
        self.declare_parameter('seed', -1)
        seed = self.get_parameter('seed').value
        # 'state': dynus_interfaces/State on /<ns>/state (sim, or hardware with
        # convert_*_to_state running); 'world': geometry_msgs/PoseStamped on
        # /<ns>/world (raw mocap, same topic convert_vicon_to_state consumes)
        self.declare_parameter('state_source', 'state')
        self.state_source = self.get_parameter('state_source').value
        # Explicit agent list for mixed-prefix fleets, overriding
        # agent_prefix/num_agents; accepts "JR01,JR02,RR03" or the
        # goal_sender.launch.py-style "['JR01', 'JR02', 'RR03']". List order
        # defines the circle slot order (first agent at angle_offset).
        self.declare_parameter('list_agents', '')
        list_agents_str = self.get_parameter('list_agents').value.strip()

        if self.mode not in ('exchange', 'random'):
            raise ValueError(f"Invalid mode '{self.mode}' (expected 'exchange' or 'random')")
        if self.state_source not in ('state', 'world'):
            raise ValueError(
                f"Invalid state_source '{self.state_source}' (expected 'state' or 'world')")

        if list_agents_str:
            if list_agents_str.startswith('['):
                self.namespaces = [str(a).strip() for a in yaml.safe_load(list_agents_str)]
            else:
                self.namespaces = [a.strip() for a in list_agents_str.split(',') if a.strip()]
            self.num_agents = len(self.namespaces)
        else:
            self.namespaces = [f'{prefix}{i:02d}' for i in range(1, self.num_agents + 1)]

        if self.mode == 'random' and self.num_agents < 2:
            raise ValueError("random mode needs at least 2 agents (no derangement exists)")
        half_min_chord = self.radius * math.sin(math.pi / self.num_agents)
        if self.goal_tolerance >= half_min_chord:
            self.get_logger().warn(
                f"goal_tolerance ({self.goal_tolerance:.2f} m) >= half the minimum "
                f"inter-slot distance ({half_min_chord:.2f} m); arrivals may trigger "
                f"instantly or ambiguously")

        self.z = 0.2 if self.use_ground_robot else 1.0

        # latest (x, y, z) per agent; None until the first state message, so an
        # agent that has not reported yet can never count as arrived
        self.positions = {ns: None for ns in self.namespaces}
        self.current_slot = {ns: k for k, ns in enumerate(self.namespaces)}
        self.targets = {}
        self.arrived = {ns: False for ns in self.namespaces}
        self.round_active = False
        self.leg = 0  # exchange only: even = outbound to antipode, odd = back home
        self.rounds_started = 0
        self.cycles_completed = 0
        self.rng = random.Random(seed if seed >= 0 else None)
        self.done_future = Future()

        # QoS: freshest sample per agent (see formation_viz_node.py) — BEST_EFFORT
        # depth 1 is compatible with the RELIABLE state publishers and never
        # drains a stale FIFO backlog
        state_qos = QoSProfile(
            depth=1,
            history=HistoryPolicy.KEEP_LAST,
            reliability=ReliabilityPolicy.BEST_EFFORT,
            durability=DurabilityPolicy.VOLATILE,
        )
        if self.state_source == 'world':
            self.state_subs = [
                self.create_subscription(
                    PoseStamped, f'/{ns}/world', self._make_world_cb(ns), state_qos)
                for ns in self.namespaces
            ]
        else:
            self.state_subs = [
                self.create_subscription(
                    State, f'/{ns}/state', self._make_state_cb(ns), state_qos)
                for ns in self.namespaces
            ]
        # term_goal must be RELIABLE to match the planner's critical_qos
        # subscription; the 1 Hz republish covers late-starting planners
        self.goal_pubs = {
            ns: self.create_publisher(PoseStamped, f'/{ns}/term_goal', 10)
            for ns in self.namespaces
        }

        self.control_timer = self.create_timer(1.0, self.control_loop)

        cycles_str = str(self.num_cycles) if self.num_cycles > 0 else 'infinite'
        self.get_logger().info(
            f"Central goal monitor: mode={self.mode}, agents={self.namespaces}, "
            f"radius={self.radius}, angle_offset={self.angle_offset:.4f}, "
            f"clockwise={self.clockwise}, cycles={cycles_str}, "
            f"tolerance={self.goal_tolerance}, state_source={self.state_source}")
        geom = (self.num_agents, self.radius, self.angle_offset, self.clockwise)
        for ns in self.namespaces:
            hx, hy = slot_point(self.current_slot[ns], *geom)
            self.get_logger().info(f"  {ns}: home slot ({hx}, {hy})")

    def _make_state_cb(self, ns):
        def cb(msg: State):
            self.positions[ns] = (msg.pos.x, msg.pos.y, msg.pos.z)
        return cb

    def _make_world_cb(self, ns):
        def cb(msg: PoseStamped):
            p = msg.pose.position
            self.positions[ns] = (p.x, p.y, p.z)
        return cb

    def control_loop(self):
        if self.done_future.done():
            return

        # Gate: publish nothing until every agent has reported its state once
        missing = [ns for ns in self.namespaces if self.positions[ns] is None]
        if missing:
            self.get_logger().info(
                f"Waiting for first state from: {', '.join(missing)}",
                throttle_duration_sec=5.0)
            return

        if not self.round_active:
            self.start_round()
            return

        self.publish_goals()

        for ns in self.namespaces:
            # latched: drift while waiting for slower agents cannot un-arrive
            if not self.arrived[ns] and self.distance_to_target(ns) < self.goal_tolerance:
                self.arrived[ns] = True
                self.get_logger().info(
                    f"{ns} reached its goal "
                    f"({sum(self.arrived.values())}/{self.num_agents})")

        if all(self.arrived.values()):
            self.complete_round()

    def distance_to_target(self, ns):
        x, y, z = self.positions[ns]
        gx, gy, gz = self.targets[ns]
        if self.use_ground_robot:
            return math.hypot(x - gx, y - gy)
        return math.sqrt((x - gx) ** 2 + (y - gy) ** 2 + (z - gz) ** 2)

    def start_round(self):
        geom = (self.num_agents, self.radius, self.angle_offset, self.clockwise)
        if self.mode == 'exchange':
            outbound = self.leg % 2 == 0
            point_fn = antipode_point if outbound else slot_point
            for ns in self.namespaces:
                x, y = point_fn(self.current_slot[ns], *geom)
                self.targets[ns] = (x, y, self.z)
            leg_name = 'outbound to antipode' if outbound else 'return home'
        else:  # random
            perm = random_derangement(self.num_agents, self.rng)
            for ns in self.namespaces:
                self.current_slot[ns] = perm[self.current_slot[ns]]
                x, y = slot_point(self.current_slot[ns], *geom)
                self.targets[ns] = (x, y, self.z)
            leg_name = 'random derangement'

        self.arrived = {ns: False for ns in self.namespaces}
        self.round_active = True
        self.rounds_started += 1
        goals_str = ', '.join(
            f"{ns}->({self.targets[ns][0]}, {self.targets[ns][1]})" for ns in self.namespaces)
        self.get_logger().info(f"Round {self.rounds_started} [{leg_name}]: {goals_str}")
        self.publish_goals()

    def complete_round(self):
        self.round_active = False
        if self.mode == 'exchange':
            self.leg += 1
            if self.leg % 2 == 0:
                self.cycles_completed += 1
                self.get_logger().info(
                    f"Exchange cycle {self.cycles_completed} complete (out and back)")
        else:
            self.cycles_completed += 1
            self.get_logger().info(f"Random round {self.cycles_completed} complete")

        if self.num_cycles > 0 and self.cycles_completed >= self.num_cycles:
            self.finish()

    def finish(self):
        self.get_logger().info(
            f"All {self.num_cycles} {self.mode} cycles complete "
            f"after {self.rounds_started} rounds. Shutting down.")
        self.control_timer.cancel()
        self.done_future.set_result(True)

    def publish_goals(self):
        now = self.get_clock().now().to_msg()
        for ns in self.namespaces:
            gx, gy, gz = self.targets[ns]
            term_goal = PoseStamped()
            term_goal.header.stamp = now
            if self.state_source == 'world':
                # mocap: goals live in the shared mocap frame (same convention
                # as repub_rviz_2Dgoal.py)
                term_goal.header.frame_id = 'world'
            else:
                term_goal.header.frame_id = f'{ns}/map' if self.use_hardware else 'map'
            term_goal.pose.position.x = float(gx)
            term_goal.pose.position.y = float(gy)
            term_goal.pose.position.z = float(gz)
            term_goal.pose.orientation.w = 1.0  # planner ignores orientation
            self.goal_pubs[ns].publish(term_goal)


def main(args=None):
    rclpy.init(args=args)
    node = GoalMonitorCentralNode()
    try:
        # completes (and the process exits, ending the launch) when all cycles
        # are done; num_cycles <= 0 runs until Ctrl-C
        rclpy.spin_until_future_complete(node, node.done_future)
    except KeyboardInterrupt:
        pass
    node.destroy_node()
    rclpy.try_shutdown()


if __name__ == '__main__':
    main()
