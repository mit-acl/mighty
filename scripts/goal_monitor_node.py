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
from std_msgs.msg import Header, String
import math


def circle_position(agent_index, num_agents, radius, z=1.0, angle_offset=0.0,
                    crossing_offset_m=0.0):
    """Compute own and opposite positions on a circle for swap goals.

    Agent i sits at angle = 2*pi*(i-1)/N + angle_offset on a circle of the
    given radius; its swap target is the diametrically opposite point
    (angle + pi).

    crossing_offset_m > 0 turns the centre crossing into a ROUNDABOUT: every
    goal is shifted tangentially by the same handedness (counter-clockwise),
    so instead of N chords intersecting at one point — where every conflict
    is head-on and unresolvable by lateral avoidance on non-holonomic robots
    — the flows curve around the centre in the same rotational sense.
    Head-on pairs become offset passes 2*offset apart. Both endpoints of the
    swap are shifted, so the return leg flows the same way. 0 preserves the
    exact historical geometry.

    Returns:
        (own_x, own_y, z, opp_x, opp_y, z)
    """
    angle = 2.0 * math.pi * (agent_index - 1) / num_agents + angle_offset
    own_x = radius * math.cos(angle)
    own_y = radius * math.sin(angle)
    opp_x = radius * math.cos(angle + math.pi)
    opp_y = radius * math.sin(angle + math.pi)
    if crossing_offset_m != 0.0:
        # Tangential (CCW) unit vector at each endpoint; shifting a GOAL by
        # the tangent at that goal biases every approach to pass the centre
        # on the same side.
        opp_x += crossing_offset_m * -math.sin(angle + math.pi)
        opp_y += crossing_offset_m * math.cos(angle + math.pi)
        own_x += crossing_offset_m * -math.sin(angle)
        own_y += crossing_offset_m * math.cos(angle)
    return (round(own_x, 3), round(own_y, 3), z,
            round(opp_x, 3), round(opp_y, 3), z)


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
        # Per-agent start stagger: agent i holds its FIRST goal until
        # i * start_stagger_sec after node start. With N agents swapping
        # antipodes, releasing all goals in the same second sends everyone
        # through the circle centre simultaneously — measured with 10 ground
        # robots, that crowd defeats every avoidance setting (non-holonomic
        # robots can't strafe out of a scrum). Staggering keeps only a few
        # robots in the crossing zone at a time, which pairwise avoidance
        # handles. Default 0.0 = historical behavior (UAV modes unchanged).
        self.declare_parameter('start_stagger_sec', 0.0)
        start_stagger_sec = self.get_parameter('start_stagger_sec').value
        # Roundabout offset for the centre crossing (see circle_position).
        self.declare_parameter('crossing_offset_m', 0.0)
        crossing_offset_m = self.get_parameter('crossing_offset_m').value
        # Goal pattern:
        #   'swap'        — classic antipodal exchange (historical behavior).
        #   'random_slot' — on every arrival, pick the next goal uniformly
        #                   from the N evenly-spaced start slots (excluding
        #                   the slot just reached). Traffic spreads naturally
        #                   in space and time instead of one simultaneous
        #                   centre event: legs vary from adjacent-slot hops
        #                   to full crossings. Seeded per-namespace so runs
        #                   are reproducible; two robots MAY briefly target
        #                   the same slot (goal_tolerance absorbs it — both
        #                   arrive and immediately draw new goals).
        self.declare_parameter('goal_pattern', 'swap')
        self.goal_pattern = self.get_parameter('goal_pattern').value
        self.declare_parameter('random_seed', 0)
        random_seed = self.get_parameter('random_seed').value
        # random_slot only: redraw after this long without progress (see
        # distance_check_callback for the deadlock chain this breaks).
        self.declare_parameter('redraw_stagnation_sec', 20.0)
        self.redraw_stagnation_sec = self.get_parameter('redraw_stagnation_sec').value
        self._stagnation_best = None
        self._stagnation_t = self.get_clock().now()
        self.distance_check_frequency = 1.0
        self.current_goal_index = 0
        self.stagger_agent_index = 0
        self.node_start_time = self.get_clock().now()

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
            self.stagger_agent_index = agent_index - 1   # NX01 starts immediately
            own_x, own_y, _, opp_x, opp_y, _ = circle_position(agent_index, num_agents, radius, z=z, angle_offset=angle_offset, crossing_offset_m=crossing_offset_m)
            if self.goal_pattern == 'random_slot':
                # All N start slots (agent i's own point for i = 1..N).
                self.slots = []
                for i in range(1, num_agents + 1):
                    sx, sy, _, _, _, _ = circle_position(
                        i, num_agents, radius, z=z, angle_offset=angle_offset)
                    self.slots.append([sx, sy, z])
                # Stable per-namespace seed: Python's hash() is randomized
                # per process (PYTHONHASHSEED), which silently made every run
                # a different draw sequence — impossible to reproduce a
                # failing run. Sum-of-ords is stable across processes.
                ns_key = sum(ord(c) * (31 ** i) for i, c in enumerate(self.namespace))
                self.rng = __import__('random').Random(ns_key ^ int(random_seed))
                # Peer slot claims (global topic, all monitors pub+sub): a
                # draw EXCLUDES slots other robots currently target, killing
                # same-slot contention at the source. DEFAULT OFF, measured
                # both ways on a 48-core host: claims removed the 20 s
                # contention waits as designed, but those waits were
                # accidental load-shedding — simultaneous-mover density
                # rose and fleet-wide close passes got WORSE (below-1m time
                # per 900 s: 13-28 s without claims vs 59-238 s with).
                # Stale claims (5 s) expire so a crashed monitor never
                # fences off a slot forever.
                self.declare_parameter('use_slot_claims', False)
                self.use_slot_claims = self.get_parameter('use_slot_claims').value
                self.peer_claims = {}
                self.current_slot = agent_index - 1
                first = self._draw_next_slot()
                self.goal_points = [self.slots[first]]
                self.current_slot = first
                self.get_logger().info(
                    f"Random-slot goals (N={num_agents}, R={radius}): "
                    f"start slot {agent_index - 1} -> first goal slot {first}")
            else:
                self.goal_points = [[opp_x, opp_y, z], [own_x, own_y, z]]
                self.get_logger().info(
                    f"Circle swap goals (N={num_agents}, R={radius}): "
                    f"start ({own_x},{own_y}) <-> opposite ({opp_x},{opp_y})")

        elif self.namespace.startswith('RR'):
            agent_index = int(self.namespace[2:])  # RR01 -> 1
            own_x, own_y, _, opp_x, opp_y, _ = circle_position(agent_index, num_agents, radius, z=z, angle_offset=angle_offset, crossing_offset_m=crossing_offset_m)
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

        # Slot-claim exchange (random_slot only; global topic, all monitors).
        # The initial draw in __init__ ran before this publisher existed;
        # publish_term_goal's 1 Hz tick re-broadcasts the current claim, so
        # peers converge within a second of startup.
        if (getattr(self, 'goal_pattern', 'swap') == 'random_slot'
                and hasattr(self, 'slots')
                and getattr(self, 'use_slot_claims', False)):
            self.claim_pub = self.create_publisher(String, '/swap/slot_claims', 10)
            self.claim_sub = self.create_subscription(
                String, '/swap/slot_claims', self._claim_callback, 20)

        self.hold_until_sec = self.stagger_agent_index * start_stagger_sec

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

    def _draw_next_slot(self):
        """Uniform draw over the slot indices, excluding the current one and
        slots freshly claimed by peers (see peer_claims). Falls back to
        exclude-current-only if peers have claimed everything else — the
        draw must never dead-end."""
        now_s = self.get_clock().now().nanoseconds * 1e-9
        claimed = set()
        if getattr(self, 'use_slot_claims', False):
            claimed = {idx for ns, (idx, t) in getattr(self, 'peer_claims', {}).items()
                       if now_s - t < 5.0}
        choices = [i for i in range(len(self.slots))
                   if i != self.current_slot and i not in claimed]
        if not choices:
            choices = [i for i in range(len(self.slots)) if i != self.current_slot]
        pick = self.rng.choice(choices)
        self._publish_claim(pick)
        return pick

    def _publish_claim(self, slot_idx):
        if getattr(self, 'claim_pub', None) is not None:
            self.claim_pub.publish(String(data=f'{self.namespace}:{slot_idx}'))

    def _claim_callback(self, msg: String):
        try:
            ns, idx = msg.data.rsplit(':', 1)
        except ValueError:
            return
        if ns == self.namespace:
            return
        now_s = self.get_clock().now().nanoseconds * 1e-9
        self.peer_claims[ns] = (int(idx), now_s)

    def _still_staggered(self):
        """True while this agent's start-stagger hold is active."""
        if self.hold_until_sec <= 0.0:
            return False
        elapsed = (self.get_clock().now() - self.node_start_time).nanoseconds * 1e-9
        return elapsed < self.hold_until_sec

    def distance_check_callback(self):
        if self._still_staggered():
            return

        # Get the current goal point
        goal_x, goal_y, goal_z = self.goal_points[self.current_goal_index]

        # random_slot self-healing: a drawn slot can be REJECTED by the
        # planner ("goal in occupied space") when a peer happens to be parked
        # on it — and a robot whose goal is rejected parks in place, blocking
        # ITS slot for others: a deadlock chain (measured: NX07 pinned at one
        # slot for 800 s across reproducible runs while the planner rejected
        # its goal once per second). The monitor cannot see rejections, but
        # it can see their symptom — zero progress — and the fix is local:
        # after redraw_stagnation_sec without the goal distance improving,
        # draw a DIFFERENT slot. Bounded pause instead of permanent deadlock,
        # and moving again frees this robot's own slot for everyone else.
        if self.goal_pattern == 'random_slot':
            d_now = math.sqrt((self.current_position.x - goal_x) ** 2 +
                              (self.current_position.y - goal_y) ** 2)
            if self._stagnation_best is None or d_now < self._stagnation_best - 0.3:
                self._stagnation_best = d_now
                self._stagnation_t = self.get_clock().now()
            else:
                held = (self.get_clock().now() - self._stagnation_t).nanoseconds * 1e-9
                if held > self.redraw_stagnation_sec:
                    nxt = self._draw_next_slot()
                    self.get_logger().warn(
                        f"Random-slot: no progress toward slot {self.current_slot} "
                        f"for {held:.0f}s (peer parked on it?) — redrawing to slot {nxt}")
                    self.current_slot = nxt
                    self.goal_points = [self.slots[nxt]]
                    self.current_goal_index = 0
                    self._stagnation_best = None
                    return

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

        # Check if the drone has reached the current goal, then advance:
        # swap wraps around its two fixed points; random_slot draws fresh.
        if distance < self.goal_tolerance:
            self.get_logger().info(f"Goal {self.current_goal_index} reached!")
            if getattr(self, 'goal_pattern', 'swap') == 'random_slot':
                nxt = self._draw_next_slot()
                self.current_slot = nxt
                self.goal_points = [self.slots[nxt]]
                self.current_goal_index = 0
                self.get_logger().info(f"Random-slot: next goal slot {nxt}")
            else:
                self.current_goal_index = (self.current_goal_index + 1) % len(self.goal_points)

    def publish_term_goal(self):
        """Publishes the current goal as a PoseStamped message."""
        if self._still_staggered():
            return
        # Re-broadcast the slot claim at this timer's 1 Hz: covers the
        # initial draw (which predates the publisher) and lets peers'
        # 5 s claim expiry work as a liveness check, not a churn source.
        if getattr(self, 'goal_pattern', 'swap') == 'random_slot' and hasattr(self, 'slots'):
            self._publish_claim(self.current_slot)
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
