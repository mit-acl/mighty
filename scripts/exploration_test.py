#!/usr/bin/env python3
# /* ----------------------------------------------------------------------------
#  * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */
"""Automated regression harness for the MIGHTY simulations.

Launches a scenario through ``run_sim.py`` (headless), watches it over the ROS
graph, decides pass/fail against explicit gates, and tears everything down. Built
for unattended campaigns, so the guiding rule is **a run can never hang**: three
independent stoppers bound every run (per-run wall-clock cap, early abort the
moment a stuck predicate trips, and a campaign-wide deadline).

Usage (the workspace overlay must be sourced — the monitor subscribes to
``dynus_interfaces`` messages)::

    source ~/code/mighty_ws/install/setup.bash

    # one run
    python3 src/mighty/scripts/exploration_test.py run --scenario ground-exploration

    # ten runs plus summary.md / summary.json
    python3 src/mighty/scripts/exploration_test.py campaign \
        --scenario ground-exploration --runs 10 --out /tmp/campaign

    # list what can be run
    python3 src/mighty/scripts/exploration_test.py scenarios

Artifacts per run: ``run<NN>.json`` (metrics + verdict) and ``run<NN>/pane*.log``
(complete pane output, streamed via ``tmux pipe-pane`` so nothing is lost to the
scrollback limit).
"""

import argparse
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
WS_DIR = REPO_DIR.parent.parent          # .../mighty_ws
DEFAULT_SETUP_BASH = WS_DIR / 'install' / 'setup.bash'
TMUX_SESSION = 'mighty_sim'

# The planner colourises status names, so log matching must strip ANSI first.
ANSI_RE = re.compile(r'\x1b\[[0-9;]*m')

# Log patterns that mean a node fell over. Checked against every pane.
CRASH_PATTERNS = [
    'process has died',
    'Traceback (most recent call last)',
    'terminate called',
    'what():',
    'Segmentation fault',
    'core dumped',
]

# ---------------------------------------------------------------------------
# Scenario registry
# ---------------------------------------------------------------------------
# Each scenario declares how to launch it and how to judge it. `kind` selects
# the monitor's state machine:
#   'exploration' — explore until no frontiers remain, then return to start
#   'goal'        — fly/drive to a fixed goal once
#   'swap'        — continuous position swapping; time-boxed, scored on progress

SCENARIOS = {
    'ground-exploration': dict(
        kind='exploration',
        mode='exploration-singleagent-ground',
        # log_level=info is required, not cosmetic: mighty_node defaults to
        # --log-level error, which suppresses every RCLCPP_INFO/WARN. That hides
        # the exploration watchdog's WARN, each "Exploration: -> frontier ..."
        # selection, and the throttled "[expl] grid ... fresh=N db=N" counters —
        # i.e. exactly the evidence needed to tell a wedged state machine from
        # selection churn.
        extra_args=['--no-rviz', '--log-level', 'info'],
        namespaces=['NX01'],
        # Observed completions span 669-1320 s, and returning home is not
        # one-way (the planner re-checks for leftover frontiers on arrival), so
        # a healthy run can be long. Sized well clear of that so the cap catches
        # hangs rather than manufacturing failures out of slow-but-correct runs.
        cap_s=2400.0,
        startup_cap_s=240.0,
        home=(0.0, 2.0),
        home_tol_m=0.5,                  # config goal_radius
        home_confirm_s=90.0,             # grace for the planner to declare arrival
        bounds=(-7.0, 22.2, -27.0, 8.0),  # config exploration.bounds
        min_coverage=0.95,               # reference run reached 0.986
        esdf_truncation_voxels=10,       # config esdf_truncation_distance
        # Reported, not gated — see the collision gate in evaluate(). Values
        # below this are worth eyeballing against the map.
        clearance_warn_m=0.25,
        max_tilt_deg=20.0,
        gazebo=True,
        description='Single ground robot, autonomous frontier exploration '
                    '(ACL_office), must return to start when done.',
    ),

    # --- UAV scenarios advertised in README.md (Native Linux section) --------
    # fake_sim scenarios have no Gazebo and no 2D ESDF, so clearance/tilt gates
    # do not apply; they are judged on goal progress, separation and altitude.
    'uav-interactive': dict(
        kind='goal',
        mode='interactive',
        extra_args=['--no-rviz'],
        namespaces=['NX01'],
        cap_s=300.0,
        startup_cap_s=180.0,
        # RViz's "2D Goal Pose" is unavailable headless, so the harness
        # publishes the same goal the click would produce.
        # The fake_sim random forest is ~12.5 m across (simulator.launch.py
        # map_size_x_/y_), so a goal beyond that would fly through empty space
        # and prove nothing. z is rewritten to default_goal_z regardless.
        publish_goal=(10.0, 0.0, 1.0),
        goal_tol_m=1.5,
        min_altitude_m=0.2,
        gazebo=False,
        description='Single UAV, fake_sim, goal injected the way RViz 2D Goal '
                    'Pose would (README: --mode interactive).',
    ),
    'uav-multiagent-10': dict(
        kind='swap',
        mode='multiagent',
        extra_args=['--no-rviz'],
        namespaces=[f'NX{i:02d}' for i in range(1, 11)],
        cap_s=480.0,               # continuous scenario: time-boxed, not "done"
        startup_cap_s=240.0,
        min_swap_legs=1,           # every agent must complete >=1 crossing
        # drone_bbox is [0.5, 0.5, 0.5], so bounding boxes overlap below 0.5 m
        # centre-to-centre. Gate on actual overlap; the observed minimum is
        # reported either way, so tight-but-legal passes show up without failing.
        min_separation_m=0.5,
        min_altitude_m=0.2,
        gazebo=False,
        description='10 UAVs swapping antipodal positions on a circle, '
                    'time-boxed (README: --mode multiagent).',
    ),
    'uav-multiagent-5': dict(
        kind='swap',
        mode='multiagent',
        extra_args=['--no-rviz', '--num-agents', '5'],
        namespaces=[f'NX{i:02d}' for i in range(1, 6)],
        cap_s=480.0,
        startup_cap_s=240.0,
        min_swap_legs=1,
        # drone_bbox is [0.5, 0.5, 0.5], so bounding boxes overlap below 0.5 m
        # centre-to-centre. Gate on actual overlap; the observed minimum is
        # reported either way, so tight-but-legal passes show up without failing.
        min_separation_m=0.5,
        min_altitude_m=0.2,
        gazebo=False,
        description='5 UAVs swapping antipodal positions '
                    '(README: --mode multiagent --num-agents 5).',
    ),
    'uav-gazebo-hard-forest': dict(
        kind='goal',
        mode='gazebo',
        extra_args=['--no-rviz'],
        namespaces=['NX01'],
        cap_s=900.0,
        startup_cap_s=300.0,
        goal=(105.0, 0.0, 3.0),    # run_sim.py --goal default (Docker uses 305)
        goal_tol_m=2.0,
        min_altitude_m=0.2,
        gazebo=True,
        description='Single UAV through hard_forest to the default goal '
                    '(README: --mode gazebo).',
    ),
    'uav-gazebo-goal': dict(
        kind='goal',
        mode='gazebo',
        extra_args=['--no-rviz', '--goal', '100', '0', '3'],
        namespaces=['NX01'],
        cap_s=900.0,
        startup_cap_s=300.0,
        goal=(100.0, 0.0, 3.0),
        goal_tol_m=2.0,
        min_altitude_m=0.2,
        gazebo=True,
        description='Single UAV to an explicit goal '
                    '(README: --mode gazebo --goal 100 0 3).',
    ),
    'uav-gazebo-easy-forest': dict(
        kind='goal',
        mode='gazebo',
        extra_args=['--no-rviz', '--env', 'easy_forest'],
        namespaces=['NX01'],
        cap_s=900.0,
        startup_cap_s=300.0,
        goal=(105.0, 0.0, 3.0),
        goal_tol_m=2.0,
        min_altitude_m=0.2,
        gazebo=True,
        description='Single UAV through easy_forest '
                    '(README: --mode gazebo --env easy_forest).',
    ),
}


def scenario(name):
    if name not in SCENARIOS:
        sys.exit(f"unknown scenario '{name}'. Known: {', '.join(sorted(SCENARIOS))}")
    return SCENARIOS[name]


# ---------------------------------------------------------------------------
# Process / tmux plumbing
# ---------------------------------------------------------------------------

def _run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL, **kw)


def _pgrep(pattern):
    """PIDs whose full command line matches `pattern` (regex)."""
    out = subprocess.run(['pgrep', '-f', pattern], stdout=subprocess.PIPE,
                         text=True)
    return [int(p) for p in out.stdout.split() if p.strip()]


# Only ever touch processes belonging to THIS workspace, plus the nodes this
# harness launched. Deliberately not using scripts/kill_ros_processes.py, which
# matches any process containing "ros" and would kill unrelated work.
# Matched against full command lines (pgrep -f). Note the planner's executable is
# named `mighty`, not `mighty_node`, and mpc_node lives in a separate package —
# match the installed paths so nothing survives into the next run and poisons it.
TEARDOWN_PATTERNS = [
    r'gz(server|client)',
    r'install/mighty/lib/mighty/',        # mighty, fake_sim, convert_odom_to_state
    r'install/mpc/lib/mpc/mpc_node',
    r'install/global_mapper_ros/lib/',
    r'ros2 launch (mighty|global_mapper_ros)',
    r'rviz2',
]


def teardown(verbose=True, timeout_s=45.0):
    """Kill the sim and *verify* it is gone. Escalates SIGTERM -> SIGKILL."""
    _run(['tmux', 'kill-session', '-t', TMUX_SESSION])
    deadline = time.time() + timeout_s
    escalated = False
    while True:
        alive = []
        for pat in TEARDOWN_PATTERNS:
            alive += _pgrep(pat)
        alive = sorted(set(alive) - {os.getpid()})
        if not alive:
            return True
        sig = signal.SIGKILL if (escalated or time.time() > deadline - 20) else signal.SIGTERM
        for pid in alive:
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
        if sig == signal.SIGKILL:
            escalated = True
        if time.time() > deadline:
            if verbose:
                print(f"[teardown] WARNING: {len(alive)} process(es) survived: {alive}")
            return False
        time.sleep(2.0)


def launch(sc, yaml_path, setup_bash, ros_domain_id, log_dir):
    """Emit the tmuxp YAML via run_sim.py, launch detached, start pane logging."""
    cmd = [sys.executable, str(SCRIPT_DIR / 'run_sim.py'),
           '--mode', sc['mode'], '-s', str(setup_bash),
           '--ros-domain-id', str(ros_domain_id),
           *sc.get('extra_args', []),
           '--emit-yaml', str(yaml_path)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    if r.returncode != 0 or not yaml_path.exists():
        raise RuntimeError(f"run_sim.py --emit-yaml failed:\n{r.stdout}")

    env = os.environ.copy()
    env['SETUP_BASH'] = str(setup_bash)
    env['ROS_DOMAIN_ID'] = str(ros_domain_id)
    r = subprocess.run(['tmuxp', 'load', '-d', str(yaml_path)], env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    if r.returncode != 0:
        raise RuntimeError(f"tmuxp load failed:\n{r.stdout}")

    # Stream every pane to disk. tmux scrollback is capped (the mighty pane holds
    # only a couple hundred lines), so capture-pane alone loses the very lines
    # that decide pass/fail. pipe-pane tees from now on, losing nothing.
    log_dir.mkdir(parents=True, exist_ok=True)
    panes = subprocess.run(['tmux', 'list-panes', '-t', f'{TMUX_SESSION}:main',
                            '-F', '#{pane_index}'],
                           stdout=subprocess.PIPE, text=True).stdout.split()
    for idx in panes:
        target = f'{TMUX_SESSION}:main.{idx}'
        out = log_dir / f'pane{idx}.log'
        _run(['tmux', 'pipe-pane', '-o', '-t', target, f'cat >> {out}'])
    return [log_dir / f'pane{i}.log' for i in panes]


# Filters Gazebo's contact stream down to real collisions at the source. The raw
# stream is ~15k lines/s, far too much to hand to Python, but awk reduces it to
# one line per foreign contact — normally zero. A contact is a collision when
# exactly one side is the robot and the other is not the floor; robot-vs-robot
# (wheel against its own chassis) and anything-vs-ground_plane are just driving.
CONTACT_AWK = r'''
/collision1:/ { split($2, a, "::"); m1 = a[1]; next }
/collision2:/ {
  split($2, b, "::"); m2 = b[1];
  if ((m1 == NS) != (m2 == NS)) {
    other = (m1 == NS) ? m2 : m1;
    if (other != "ground_plane") { print systime() " " other; fflush(); }
  }
}
'''

# Contacts arrive once per physics step while touching, so raw counts are
# meaningless as severity. Group them into events separated by this gap: a 2 s
# graze that the robot backs out of is a very different result from a wedge.
CONTACT_EVENT_GAP_S = 3.0


def start_contact_watch(ns, out_file):
    """Watch every Gazebo contact for the whole run (not a sample).

    Restarts itself if gzserver isn't up yet or the stream drops, so it can be
    started before the simulation finishes launching.
    """
    script = (
        f'while true; do '
        f'  gz topic -e /gazebo/default/physics/contacts 2>/dev/null '
        f"    | awk -F'\"' -v NS={ns} '{CONTACT_AWK}' >> {out_file}; "
        f'  sleep 2; '
        f'done')
    return subprocess.Popen(['bash', '-c', script],
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL,
                            preexec_fn=os.setsid)


def stop_contact_watch(proc):
    if proc is None:
        return
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def read_logs(log_paths):
    text = []
    for p in log_paths:
        try:
            text.append(p.read_text(errors='replace'))
        except FileNotFoundError:
            pass
    return '\n'.join(text)


# ---------------------------------------------------------------------------
# Monitor
# ---------------------------------------------------------------------------

def import_ros():
    try:
        import rclpy                                             # noqa: F401
        from dynus_interfaces.msg import State                    # noqa: F401
        from nav_msgs.msg import OccupancyGrid                    # noqa: F401
        from geometry_msgs.msg import PoseStamped                 # noqa: F401
        from gazebo_msgs.msg import ModelStates                   # noqa: F401
    except ImportError as e:
        sys.exit(
            f"ROS 2 imports failed ({e}).\n"
            f"Source the workspace overlay first:\n"
            f"  source {DEFAULT_SETUP_BASH}")


class Monitor:
    """Samples the ROS graph and decides when (and why) a run is over.

    Terminal outcomes:
        PASS, FAIL_<gate>, STUCK_NO_PROGRESS, STUCK_NO_COVERAGE,
        STUCK_OSCILLATING, TIMEOUT, LAUNCH_FAILED
    """

    # --- Stuck thresholds -------------------------------------------------
    # Tuned to be SPECIFIC rather than sensitive: a false positive sends the
    # campaign chasing a phantom bug, while a slow catch costs at most one
    # wall-clock cap. Every predicate therefore needs corroborating evidence
    # that the robot is not simply travelling or turning in place, and the
    # per-run cap is the universal backstop.
    TICK_S = 0.25                      # metric sampling period

    # Wedged: commanded but not translating at all.
    NO_PROGRESS_WINDOW_S = 90.0
    NO_PROGRESS_EPS_M = 0.05

    # Exploration stalled: the only trustworthy signal that a *frontier
    # explorer* is stuck is that it has stopped learning about the world.
    # Motion-based predicates alone raise false alarms, because greedy
    # nearest-frontier selection legitimately makes long round trips: finish a
    # dead-end corridor, then drive all the way back because the next-nearest
    # frontier is behind you. That looks identical to flip-flopping if you only
    # watch net displacement, and an earlier revision of this harness killed a
    # healthy run at 90% coverage for exactly that reason.
    NO_COVERAGE_WINDOW_S = 240.0
    NO_COVERAGE_EPS = 0.001            # 0.1% of the bounds area

    # Once coverage has already met the scenario's gate, stagnation is not a
    # coverage failure — the run has done its job. NO_PROGRESS and the
    # wall-clock cap still bound it.
    # A stalled run is reported as OSCILLATING rather than NO_COVERAGE when the
    # evidence shows it driving hard and getting nowhere. A retarget only counts
    # when the goal JUMPS: the manager's EMA centroid update
    # (centroid_ema_alpha) slides the published goal continuously, so
    # position-change alone measures jitter, not target switching.
    OSC_WINDOW_S = 90.0
    OSC_MIN_JUMP_M = 2.0
    OSC_SWITCHES = 6
    OSC_NET_DISP_M = 1.5
    OSC_MIN_PATH_M = 6.0

    # The robot must have travelled at least this far from its start before a
    # home-pointing term_goal counts as "returning" rather than "just spawned".
    RETURN_MIN_EXCURSION_M = 5.0

    def __init__(self, sc, ns, log_paths):
        import rclpy
        from rclpy.node import Node
        from rclpy.qos import (QoSProfile, DurabilityPolicy, ReliabilityPolicy,
                               HistoryPolicy, qos_profile_sensor_data)
        from dynus_interfaces.msg import State
        from nav_msgs.msg import OccupancyGrid
        from geometry_msgs.msg import PoseStamped
        from gazebo_msgs.msg import ModelStates

        self.sc = sc
        self.ns = ns                     # primary agent (metrics/phase machine)
        self.all_ns = sc['namespaces']
        self.kind = sc['kind']
        self.log_paths = log_paths

        self.node = Node('exploration_test_monitor')
        tl_qos = QoSProfile(depth=1, history=HistoryPolicy.KEEP_LAST,
                            reliability=ReliabilityPolicy.RELIABLE,
                            durability=DurabilityPolicy.TRANSIENT_LOCAL)

        # BEST_EFFORT + VOLATILE everywhere it is safe: a best-effort subscriber
        # matches reliable AND best-effort publishers, so the monitor can't
        # silently see nothing because a node picked a stricter profile than we
        # guessed. Dropping the odd sample is harmless for monitoring.
        lax_qos = QoSProfile(depth=10, history=HistoryPolicy.KEEP_LAST,
                             reliability=ReliabilityPolicy.BEST_EFFORT,
                             durability=DurabilityPolicy.VOLATILE)

        self.state = None
        # Every agent's pose, for separation / altitude / per-agent progress.
        self.peer_pos = {}
        for a in self.all_ns:
            self.node.create_subscription(
                State, f'/{a}/state',
                lambda m, aa=a: self._on_any_state(aa, m), lax_qos)
            # goal_monitor republishes term_goal each time an agent arrives, so
            # counting distinct term_goals counts completed swap legs.
            self.node.create_subscription(
                PoseStamped, f'/{a}/term_goal',
                lambda m, aa=a: self._on_term_goal(aa, m), lax_qos)

        if self.kind == 'exploration':
            self.node.create_subscription(
                OccupancyGrid, f'/{ns}/exploration/visited_map',
                self._on_visited, tl_qos)
            self.node.create_subscription(
                OccupancyGrid, f'/{ns}/esdf_2d_topic', self._on_esdf,
                qos_profile_sensor_data)
            self.node.create_subscription(
                PoseStamped, f'/{ns}/exploration/current_goal', self._on_goal,
                lax_qos)
        if sc.get('gazebo'):
            self.node.create_subscription(
                ModelStates, '/plug/model_states_plug', self._on_models, lax_qos)

        # Headless equivalent of clicking "2D Goal Pose" in RViz. Only the
        # scenarios that declare publish_goal need it; --mode gazebo already has
        # goal_sender.py issuing the goal.
        self._goal_pub = None
        if sc.get('publish_goal'):
            self._goal_pub = self.node.create_publisher(
                PoseStamped, f'/{ns}/term_goal', 10)

        # --- accumulated metrics ---
        self.t0 = time.time()
        self.t_first_state = None
        self.phase = 'WAIT_START'        # WAIT_START -> EXPLORING -> RETURNING
        self.coverage = 0.0
        self.coverage_peak = 0.0
        self.min_clearance = math.inf
        self.min_clearance_pos = None
        self.clearance_samples = []
        self.max_tilt_deg = 0.0
        self.obstacle_contacts = 0
        self.contact_models = []
        self.contact_events = []
        self.contact_file = None
        self.min_altitude = math.inf
        self.path_length_m = 0.0
        self.goal_switches = 0
        self.samples = 0
        self.goal_history = []           # [[t_rel, x, y], ...] on each retarget
        self.pos_track = []              # [[t_rel, x, y], ...] every 2 s
        self._t_last_track = 0.0
        self.swap_legs = {}              # ns -> completed swap legs
        self._last_term_goal = {}
        self.min_separation_m = math.inf
        self.goal_reached = False
        self._goal_published = False

        self._esdf = None                # (arr, W, H, res, ox, oy, dmax)
        self._bounds_mask = None
        self._last_pos = None
        self._last_goal = None
        self._t_last_tick = 0.0
        # Trailing (t, x, y, cumulative_path) ring for windowed net-displacement
        # and windowed path-length queries. Decimated to TICK_S so it stays small.
        self._track = deque()
        self._prog_ref = (None, 0.0)     # (pos, t) for no-progress window
        self._cov_ref = (0.0, 0.0)       # (coverage, t) for no-coverage window
        self._osc_events = []            # timestamps of REAL retargets
        # Diagnostics captured at the moment a stuck predicate trips, so a
        # false positive can be audited from the run JSON instead of guessed at.
        self.stuck_evidence = {}
        self._log_cache = ''
        self._t_last_log = 0.0
        self.saw_no_frontiers = False
        self.n_giveups = 0
        self._giveups_handled = 0
        self._reached_home_this_return = False
        self.saw_arrival_confirmed = False
        self._t_arrived = None
        self.max_dist_home = 0.0

    # -- callbacks ---------------------------------------------------------
    def _on_any_state(self, ns, msg):
        self.peer_pos[ns] = (msg.pos.x, msg.pos.y, msg.pos.z)
        self.min_altitude = min(self.min_altitude, msg.pos.z)
        if ns == self.ns:
            self.state = msg
            if self.t_first_state is None:
                self.t_first_state = time.time()

    def _on_term_goal(self, ns, msg):
        p = (msg.pose.position.x, msg.pose.position.y, msg.pose.position.z)
        prev = self._last_term_goal.get(ns)
        if prev is None or math.dist(p, prev) > 1.0:
            if prev is not None:
                # A fresh term_goal means goal_monitor saw this agent arrive.
                self.swap_legs[ns] = self.swap_legs.get(ns, 0) + 1
            self._last_term_goal[ns] = p

    def _on_goal(self, msg):
        p = (msg.pose.position.x, msg.pose.position.y)
        if self._last_goal is None or \
                math.dist(p, self._last_goal) >= self.OSC_MIN_JUMP_M:
            if self._last_goal is not None:
                self.goal_switches += 1
                self._osc_events.append(time.time())
            # Retarget history, for reconstructing oscillation geometry offline.
            self.goal_history.append(
                [round(time.time() - self.t0, 1), round(p[0], 2), round(p[1], 2)])
        self._last_goal = p

    def _on_visited(self, msg):
        arr = np.frombuffer(bytes(msg.data), dtype=np.int8)
        W, H = msg.info.width, msg.info.height
        if arr.size != W * H:
            return
        if self._bounds_mask is None or self._bounds_mask.shape != (H, W):
            res = msg.info.resolution
            ox, oy = msg.info.origin.position.x, msg.info.origin.position.y
            xs = ox + (np.arange(W) + 0.5) * res
            ys = oy + (np.arange(H) + 0.5) * res
            x0, x1, y0, y1 = self.sc['bounds']
            mx = (xs >= x0) & (xs <= x1)
            my = (ys >= y0) & (ys <= y1)
            self._bounds_mask = np.outer(my, mx)
        inb = arr.reshape(H, W)[self._bounds_mask]
        if inb.size:
            self.coverage = float((inb >= 0).sum()) / float(inb.size)
            self.coverage_peak = max(self.coverage_peak, self.coverage)

    def _on_esdf(self, msg):
        arr = np.frombuffer(bytes(msg.data), dtype=np.int8).astype(np.float32)
        W, H = msg.info.width, msg.info.height
        if arr.size != W * H:
            return
        res = msg.info.resolution
        dmax = self.sc['esdf_truncation_voxels'] * res
        # Mapper encodes value = 100 * (1 - d / dmax); mirror EsdfGrid2D exactly.
        dist = np.where(arr <= 0, dmax,
                        np.where(arr >= 100, 0.0, dmax * (1.0 - arr / 100.0)))
        self._esdf = (dist.reshape(H, W), W, H, res,
                      msg.info.origin.position.x + 0.5 * res,
                      msg.info.origin.position.y + 0.5 * res, dmax)

    def _on_models(self, msg):
        try:
            i = msg.name.index(self.ns)
        except ValueError:
            return
        q = msg.pose[i].orientation
        # roll/pitch from quaternion
        sinr = 2.0 * (q.w * q.x + q.y * q.z)
        cosr = 1.0 - 2.0 * (q.x * q.x + q.y * q.y)
        roll = math.atan2(sinr, cosr)
        sinp = max(-1.0, min(1.0, 2.0 * (q.w * q.y - q.z * q.x)))
        pitch = math.asin(sinp)
        tilt = math.degrees(max(abs(roll), abs(pitch)))
        self.max_tilt_deg = max(self.max_tilt_deg, tilt)

    # -- helpers -----------------------------------------------------------
    def _clearance_at(self, x, y):
        if self._esdf is None:
            return None
        dist, W, H, res, ox, oy, dmax = self._esdf
        ix = int(math.floor((x - ox) / res))
        iy = int(math.floor((y - oy) / res))
        if ix < 0 or ix >= W or iy < 0 or iy >= H:
            return dmax                      # outside window: no penalty
        return float(dist[iy, ix])

    # --- Gazebo contact oracle -------------------------------------------
    # Collisions are read from Gazebo's own contact stream, filtered by awk in a
    # side process (see start_contact_watch) that runs for the ENTIRE run. This
    # is ground truth, and it is the only evidence here that can establish
    # "didn't hit anything" — no perception-derived metric can, as this campaign
    # demonstrated: the ESDF at the robot's centre reported 0.00 m on a run where
    # Gazebo recorded zero contacts.
    def _read_contacts(self):
        if self.contact_file is None:
            return
        try:
            raw = Path(self.contact_file).read_text().splitlines()
        except (FileNotFoundError, OSError):
            return
        stamps = []
        for line in raw:
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            try:
                stamps.append((int(parts[0]), parts[1].strip()))
            except ValueError:
                continue
        self.obstacle_contacts = len(stamps)
        for _, n in stamps:
            if n not in self.contact_models:
                self.contact_models.append(n)

        # Group into events and locate each one, so a brief recoverable graze is
        # distinguishable from being wedged against a wall.
        events, cur = [], None
        for t, n in sorted(stamps):
            if cur and t - cur['end'] <= CONTACT_EVENT_GAP_S:
                cur['end'] = t
            else:
                if cur:
                    events.append(cur)
                cur = dict(start=t, end=t, model=n)
        if cur:
            events.append(cur)
        for e in events:
            e['duration_s'] = e['end'] - e['start'] + 1
            e['t_into_run_s'] = round(e['start'] - self.t0, 1)
            near = [p for p in self.pos_track
                    if abs(p[0] - e['t_into_run_s']) <= 3]
            e['pos'] = near[len(near) // 2][1:] if near else None
            del e['start'], e['end']
        self.contact_events = events

    def settle_and_refresh(self, seconds=6.0):
        """Keep spinning briefly after a terminal state, then re-read the logs.

        Terminal detection fires the instant the robot is inside the home
        tolerance, which is *before* the controller prints its own confirmation.
        Without this settle the finished_cleanly gate loses a race it should win.
        """
        import rclpy
        t_end = time.time() + seconds
        while time.time() < t_end:
            rclpy.spin_once(self.node, timeout_sec=0.2)
        self._t_last_log = 0.0
        self._refresh_logs()
        self._read_contacts()

    def _window(self, now, span_s):
        """(net_displacement_m, path_length_m) over the trailing `span_s`."""
        if not self._track:
            return 0.0, 0.0
        cutoff = now - span_s
        oldest = None
        for t, x, y, cum in self._track:
            if t >= cutoff:
                oldest = (t, x, y, cum)
                break
        if oldest is None:
            oldest = self._track[-1]
        t_o, x_o, y_o, cum_o = oldest
        _, x_n, y_n, cum_n = self._track[-1]
        return math.dist((x_o, y_o), (x_n, y_n)), cum_n - cum_o

    def _refresh_logs(self):
        now = time.time()
        if now - self._t_last_log < 5.0:
            return
        self._t_last_log = now
        # Strip ANSI colour codes: the planner wraps status names in them, so
        # naive substring matching on "to status_=X" would never match.
        self._log_cache = ANSI_RE.sub('', read_logs(self.log_paths))
        self.n_giveups = self._log_cache.count('No frontiers left')
        idx = self._log_cache.find('No frontiers left')
        if idx >= 0:
            self.saw_no_frontiers = True
            # Arrival must be confirmed AFTER the decision to go home, and it
            # must be a transition INTO GOAL_REACHED. Matching the bare token
            # would also match the *source* side of "from GOAL_REACHED to
            # TRAVELING", which is how every run starts and how each
            # intermediate subgoal completes — both false positives.
            tail = self._log_cache[idx:]
            self.saw_arrival_confirmed = ('MPC: goal reached' in tail
                                          or 'to status_=GOAL_REACHED' in tail)
        # For fixed-goal scenarios, the planner declaring arrival is
        # authoritative even if relocate_occupied_goal moved the goal.
        if self.kind == 'goal':
            self.saw_arrival_confirmed = (
                'to status_=GOAL_REACHED' in self._log_cache
                or 'MPC: goal reached' in self._log_cache)

    # -- main loop ---------------------------------------------------------
    def spin(self, deadline_epoch=None):
        import rclpy
        cap = self.sc['cap_s']
        startup_cap = self.sc['startup_cap_s']

        while True:
            rclpy.spin_once(self.node, timeout_sec=0.25)
            now = time.time()
            elapsed = now - self.t0
            self._refresh_logs()

            # Stopper 3: campaign-wide deadline.
            if deadline_epoch and now >= deadline_epoch:
                return 'CAMPAIGN_DEADLINE'

            # Never wait forever for the sim to come up.
            if self.state is None:
                if elapsed > startup_cap:
                    return 'LAUNCH_FAILED'
                continue

            # Decimate metric updates: /state arrives at >100 Hz, so sampling
            # every message would inflate path length with jitter and grow the
            # trailing track without bound.
            if now - self._t_last_tick < self.TICK_S:
                continue
            self._t_last_tick = now

            pos = (self.state.pos.x, self.state.pos.y, self.state.pos.z)
            self.samples += 1
            self.min_altitude = min(self.min_altitude, pos[2])
            c = self._clearance_at(pos[0], pos[1])
            if c is not None:
                self.clearance_samples.append(c)
                if c < self.min_clearance:
                    self.min_clearance = c
                    # Where the closest approach happened, so it can be checked
                    # against the map instead of trusted blindly.
                    self.min_clearance_pos = [round(pos[0], 2), round(pos[1], 2)]
            if self._last_pos is not None:
                self.path_length_m += math.dist(self._last_pos[:2], pos[:2])
            self._last_pos = pos
            self._track.append((now, pos[0], pos[1], self.path_length_m))
            while self._track and now - self._track[0][0] > 2 * self.OSC_WINDOW_S:
                self._track.popleft()
            if now - self._t_last_track >= 2.0:
                self._t_last_track = now
                self.pos_track.append([round(now - self.t0, 1),
                                       round(pos[0], 2), round(pos[1], 2)])

            # Track inter-agent separation (the real safety property of the
            # antipodal swap scenario: everyone crosses through the centre).
            if len(self.peer_pos) > 1:
                pts = list(self.peer_pos.values())
                for i in range(len(pts)):
                    for j in range(i + 1, len(pts)):
                        self.min_separation_m = min(
                            self.min_separation_m, math.dist(pts[i], pts[j]))

            if self.phase == 'WAIT_START':
                self.phase = 'EXPLORING'
                self._prog_ref = (pos, now)
                self._cov_ref = (self.coverage, now)
                # Inject the goal once the agent is actually alive, so it isn't
                # published into the void before the planner subscribes.
                if self._goal_pub is not None and not self._goal_published:
                    from geometry_msgs.msg import PoseStamped as _PS
                    g = _PS()
                    g.header.frame_id = 'map'
                    g.header.stamp = self.node.get_clock().now().to_msg()
                    gx, gy, gz = self.sc['publish_goal']
                    g.pose.position.x, g.pose.position.y, g.pose.position.z = gx, gy, gz
                    g.pose.orientation.w = 1.0
                    for _ in range(5):           # best-effort: repeat a few times
                        self._goal_pub.publish(g)
                        rclpy.spin_once(self.node, timeout_sec=0.1)
                    self._goal_published = True

            # Stopper 1: per-run wall-clock cap. For the continuous swap
            # scenario the cap IS the intended end of the run, not a failure.
            if elapsed > cap:
                if self.kind == 'swap':
                    return 'TIME_BOX_ELAPSED'
                return 'TIMEOUT'

            # --- kind-specific success ---------------------------------------
            if self.kind == 'goal':
                gx, gy, _ = self.sc.get('goal') or self.sc['publish_goal']
                # Horizontal distance only: mighty.yaml sets force_goal_z with
                # default_goal_z=1.0, so the planner rewrites the commanded
                # altitude regardless of what the goal asked for (even
                # `--goal 105 0 3` flies at 1.0 m). Comparing in 3D would make
                # this gate unreachable by construction. Altitude is checked
                # separately by the stayed_airborne gate.
                if (math.dist(pos[:2], (gx, gy)) <= self.sc['goal_tol_m']
                        or self.saw_arrival_confirmed):
                    self.goal_reached = True
                    return 'DONE'

            # Phase is derived from the CURRENT term_goal every tick, not latched
            # from a log line. Returning home is not a one-way trip: on arrival
            # the planner lifts its return-home latch and re-checks for leftover
            # frontiers (home_return_resume_on_arrival_), so a run can legitimately
            # go explore -> home -> explore -> home. A latched phase mis-handles
            # that; an earlier revision kept an arrival timer running while the
            # robot was 26 m away and wrongly reported the arrival unconfirmed.
            if self.kind == 'exploration':
                d_home_now = math.dist(pos[:2], self.sc['home'])
                self.max_dist_home = max(self.max_dist_home, d_home_now)
                # Driven by the planner's own give-up message, because there is
                # no ROS-side signal to use: exploration goals are issued by a
                # direct internal call to terminalGoalCallbackImpl(), NOT
                # published on term_goal, so term_goal never carries either the
                # frontier goals or the home goal.
                #
                # Counting occurrences (rather than latching on the first) is
                # what makes explore -> home -> explore -> home work: the planner
                # lifts its return-home latch on arrival and re-checks for
                # leftover frontiers, so a single run can give up several times.
                new_phase = self.phase
                if self.n_giveups > self._giveups_handled:
                    self._giveups_handled = self.n_giveups
                    new_phase = 'RETURNING'
                    self._reached_home_this_return = False
                elif self.phase == 'RETURNING':
                    if d_home_now <= self.sc['home_tol_m']:
                        self._reached_home_this_return = True
                    elif (self._reached_home_this_return
                          and d_home_now > self.RETURN_MIN_EXCURSION_M):
                        # Arrived, then set off again: the planner lifted its
                        # latch and found more to explore. Only *this* counts as
                        # resuming. Testing "far from home" alone would fire on
                        # the tick after the give-up, while the robot is still
                        # 20 m out and has not returned yet — which flipped the
                        # phase straight back and left the robot's eventual
                        # parking at home looking like a stall.
                        new_phase = 'EXPLORING'
                if new_phase != self.phase:
                    self.phase = new_phase
                    self._prog_ref = (pos, now)
                    self._cov_ref = (self.coverage, now)
                    self._t_arrived = None    # a fresh approach starts fresh

            # Success: parked at the captured start after a clean finish. Being
            # inside the tolerance is not enough on its own — the planner can
            # still be TRAVELING. Wait (bounded) for it to declare arrival,
            # otherwise the harness tears the sim down mid-manoeuvre and the
            # "did it finish cleanly?" question is decided by a race.
            if self.phase == 'RETURNING':
                hx, hy = self.sc['home']
                if math.dist((pos[0], pos[1]), (hx, hy)) <= self.sc['home_tol_m']:
                    if self._t_arrived is None:
                        self._t_arrived = now
                    if self.saw_arrival_confirmed:
                        return 'DONE'
                    if now - self._t_arrived > self.sc['home_confirm_s']:
                        # Physically home but never declared it. Real finding,
                        # distinct from both success and failure-to-return.
                        return 'DONE_UNCONFIRMED'
                else:
                    # Still manoeuvring toward home: the confirmation clock only
                    # runs while actually inside the tolerance.
                    self._t_arrived = None

            # Stopper 2: early abort on stuck. Each predicate records the
            # evidence that tripped it so a false positive is auditable.
            parked_at_home = (self.kind == 'exploration'
                              and self.phase == 'RETURNING'
                              and math.dist(pos[:2], self.sc['home'])
                                  <= self.sc['home_tol_m'])
            ref_pos, ref_t = self._prog_ref
            if parked_at_home:
                # Standing still at the captured start is success, not a stall.
                self._prog_ref = (pos, now)
            elif math.dist(ref_pos[:2], pos[:2]) > self.NO_PROGRESS_EPS_M:
                self._prog_ref = (pos, now)
            elif now - ref_t > self.NO_PROGRESS_WINDOW_S:
                self.stuck_evidence = dict(
                    window_s=round(now - ref_t, 1),
                    displacement_m=round(math.dist(ref_pos[:2], pos[:2]), 3),
                    pos=[round(v, 2) for v in pos], phase=self.phase)
                return 'STUCK_NO_PROGRESS'

            if self.kind == 'exploration' and self.phase == 'EXPLORING':
                self._osc_events = [t for t in self._osc_events
                                    if now - t <= self.OSC_WINDOW_S]
                cov_ref, cov_t = self._cov_ref
                if self.coverage - cov_ref > self.NO_COVERAGE_EPS:
                    self._cov_ref = (self.coverage, now)
                elif (now - cov_t > self.NO_COVERAGE_WINDOW_S and
                      self.coverage < self.sc['min_coverage']):
                    # Learned nothing new for four minutes and hasn't yet met the
                    # coverage gate: genuinely stalled. Label it by the evidence.
                    net, path = self._window(now, self.OSC_WINDOW_S)
                    self.stuck_evidence = dict(
                        stagnant_s=round(now - cov_t, 1),
                        coverage=round(self.coverage, 4),
                        retargets_in_window=len(self._osc_events),
                        window_net_disp_m=round(net, 2),
                        window_path_m=round(path, 2))
                    if (len(self._osc_events) >= self.OSC_SWITCHES and
                            net < self.OSC_NET_DISP_M and
                            path > self.OSC_MIN_PATH_M):
                        return 'STUCK_OSCILLATING'
                    return 'STUCK_NO_COVERAGE'

    def destroy(self):
        try:
            self.node.destroy_node()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Gate evaluation
# ---------------------------------------------------------------------------

def evaluate(sc, mon, outcome, duration_s, logs):
    """Turn raw metrics into per-gate pass/fail plus an overall verdict."""
    gates = {}
    crashes = sorted({p for p in CRASH_PATTERNS if p in logs})

    # A continuous swap scenario has no "done" — elapsing the time box is the
    # intended end, so it counts as reaching a terminal state.
    ok_terminal = ('DONE', 'DONE_UNCONFIRMED', 'TIME_BOX_ELAPSED')
    gates['reached_terminal_state'] = (outcome in ok_terminal)
    gates['no_stuck'] = not outcome.startswith('STUCK')
    gates['no_timeout'] = outcome not in ('TIMEOUT', 'CAMPAIGN_DEADLINE')
    gates['launched'] = outcome != 'LAUNCH_FAILED'
    gates['clean_logs'] = not crashes

    home_d = None
    if sc['kind'] == 'exploration':
        gates['coverage'] = mon.coverage_peak >= sc['min_coverage']
        gates['finished_cleanly'] = (mon.saw_no_frontiers
                                     and mon.saw_arrival_confirmed)
        home_d = (math.dist(mon._last_pos[:2], sc['home'])
                  if mon._last_pos else float('inf'))
        gates['returned_home'] = home_d <= sc['home_tol_m']
        gates['not_tilted'] = mon.max_tilt_deg < sc['max_tilt_deg']
        # Collision gate uses Gazebo's own contact stream, i.e. ground truth.
        # The ESDF-derived clearance is reported but NOT gated: sampled at the
        # robot's centre it under-reports true distance (measured 0.42 m where
        # the nearest occupied cell was 0.645 m), because the mapper's 2D ESDF
        # is built from an inflated/height-filtered projection. It is a useful
        # trend, not a threshold.
        gates['no_obstacle_contact'] = (mon.obstacle_contacts == 0)

    elif sc['kind'] == 'goal':
        gates['goal_reached'] = mon.goal_reached
        gates['stayed_airborne'] = (mon.min_altitude >= sc['min_altitude_m']
                                    if math.isfinite(mon.min_altitude) else False)

    elif sc['kind'] == 'swap':
        # Every agent must have completed at least min_swap_legs crossings, and
        # all agents must have reported state at all (catches partial launches).
        legs = {a: mon.swap_legs.get(a, 0) for a in sc['namespaces']}
        gates['all_agents_alive'] = len(mon.peer_pos) == len(sc['namespaces'])
        gates['all_agents_swapped'] = all(v >= sc['min_swap_legs']
                                          for v in legs.values())
        gates['separation'] = (mon.min_separation_m >= sc['min_separation_m']
                               if math.isfinite(mon.min_separation_m) else False)
        gates['stayed_airborne'] = (mon.min_altitude >= sc['min_altitude_m']
                                    if math.isfinite(mon.min_altitude) else False)

    metrics = dict(
        outcome=outcome,
        duration_s=round(duration_s, 1),
        coverage_known_frac=round(mon.coverage_peak, 4),
        min_clearance_m=(round(mon.min_clearance, 3)
                         if math.isfinite(mon.min_clearance) else None),
        min_clearance_pos=mon.min_clearance_pos,
        clearance_p01_m=(round(float(np.percentile(mon.clearance_samples, 1)), 3)
                         if mon.clearance_samples else None),
        clearance_p05_m=(round(float(np.percentile(mon.clearance_samples, 5)), 3)
                         if mon.clearance_samples else None),
        clearance_warn=(math.isfinite(mon.min_clearance) and
                        mon.min_clearance < sc.get('clearance_warn_m', 0.0)),
        max_tilt_deg=round(mon.max_tilt_deg, 2),
        min_altitude_m=(round(mon.min_altitude, 3)
                        if math.isfinite(mon.min_altitude) else None),
        final_dist_to_home_m=(round(home_d, 3)
                              if home_d is not None and math.isfinite(home_d)
                              else None),
        path_length_m=round(mon.path_length_m, 1),
        goal_switches=mon.goal_switches,
        saw_no_frontiers=mon.saw_no_frontiers,
        saw_arrival_confirmed=mon.saw_arrival_confirmed,
        crash_patterns=crashes,
        metric_ticks=mon.samples,
        stuck_evidence=mon.stuck_evidence,
        goal_reached=mon.goal_reached,
        obstacle_contacts=mon.obstacle_contacts,
        contact_models=mon.contact_models,
        contact_events=mon.contact_events,
        min_separation_m=(round(mon.min_separation_m, 3)
                          if math.isfinite(mon.min_separation_m) else None),
        swap_legs={a: mon.swap_legs.get(a, 0) for a in sc['namespaces']}
                  if sc['kind'] == 'swap' else None,
        agents_seen=len(mon.peer_pos),
    )
    return gates, metrics, all(gates.values())


# ---------------------------------------------------------------------------
# Run / campaign drivers
# ---------------------------------------------------------------------------

def do_run(sc, name, out_dir, setup_bash, ros_domain_id, run_idx,
           deadline_epoch=None, attach=False, log_file=None):
    """Run and judge one scenario.

    attach=True monitors a simulation somebody else started and leaves teardown
    to them — used for Docker validation, where the sim runs inside a container
    with its own tmux server that the host cannot reach. In that mode the pane
    logs come from `log_file` (the container's stdout) instead of tmux.
    """
    import rclpy

    run_dir = out_dir / f'run{run_idx:02d}'
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = run_dir / 'sim.yaml'
    started = time.time()
    outcome, mon = None, None
    contact_watch, contact_file = None, None
    try:
        if attach:
            print(f"[run {run_idx:02d}] attaching to a running {name} ...",
                  flush=True)
            log_paths = [Path(log_file)] if log_file else []
        else:
            print(f"[run {run_idx:02d}] tearing down any previous sim ...",
                  flush=True)
            teardown()
            print(f"[run {run_idx:02d}] launching {name} ...", flush=True)
            log_paths = launch(sc, yaml_path, setup_bash, ros_domain_id, run_dir)
        if sc.get('gazebo'):
            contact_file = run_dir / 'contacts.log'
            contact_file.touch()
            contact_watch = start_contact_watch(sc['namespaces'][0], contact_file)
        if not rclpy.ok():
            rclpy.init()
        mon = Monitor(sc, sc['namespaces'][0], log_paths)
        mon.contact_file = contact_file
        outcome = mon.spin(deadline_epoch=deadline_epoch)
        # Let trailing controller/planner log lines land before judging logs.
        mon.settle_and_refresh()
    except Exception as e:                      # launch failure, etc.
        outcome = f'HARNESS_ERROR: {e}'
        log_paths = sorted(run_dir.glob('pane*.log'))
    finally:
        stop_contact_watch(contact_watch)
    duration = time.time() - started

    logs = read_logs(log_paths)
    if mon is not None and not outcome.startswith('HARNESS_ERROR'):
        gates, metrics, passed = evaluate(sc, mon, outcome, duration, logs)
        mon.destroy()
    else:
        gates = {'launched': False}
        metrics = dict(outcome=outcome, duration_s=round(duration, 1))
        passed = False

    print(f"[run {run_idx:02d}] {'PASS' if passed else 'FAIL'} "
          f"outcome={metrics['outcome']} "
          f"coverage={metrics.get('coverage_known_frac')} "
          f"clearance={metrics.get('min_clearance_m')} "
          f"{metrics['duration_s']}s", flush=True)
    failed = [g for g, ok in gates.items() if not ok]
    if failed:
        print(f"[run {run_idx:02d}] failed gates: {', '.join(failed)}", flush=True)

    record = dict(scenario=name, run=run_idx, passed=passed,
                  gates=gates, metrics=metrics,
                  started_epoch=round(started, 1))
    if mon is not None:
        # Kept out of `metrics` so summary tables stay readable; these are the
        # traces used to reconstruct a failure offline.
        record['traces'] = dict(goal_history=mon.goal_history,
                                pos_track=mon.pos_track)
    (out_dir / f'run{run_idx:02d}.json').write_text(json.dumps(record, indent=2))

    if not attach:            # in attach mode the caller owns the sim's lifetime
        teardown()
    return record


def write_summary(out_dir, name, records, deadline_note=''):
    n = len(records)
    npass = sum(1 for r in records if r['passed'])
    kind = SCENARIOS.get(name, {}).get('kind', 'exploration')

    def vals(key):
        return [r['metrics'][key] for r in records
                if r['metrics'].get(key) is not None]

    lines = [f'# Campaign summary — `{name}`', '',
             f'**{npass}/{n} passed.**', '']
    if deadline_note:
        lines += [deadline_note, '']

    dur = vals('duration_s')
    agg = [f'| duration mean / max | {np.mean(dur):.0f} s / {max(dur):.0f} s |'] if dur else []
    if kind == 'exploration':
        cov, clr = vals('coverage_known_frac'), vals('min_clearance_m')
        if cov:
            agg.insert(0, f'| coverage mean / min | {np.mean(cov):.4f} / {min(cov):.4f} |')
        if clr:
            agg.append(f'| min clearance over all runs | {min(clr):.3f} m |')
    else:
        sep = vals('min_separation_m')
        if sep:
            agg.append(f'| min inter-agent separation | {min(sep):.3f} m |')
    if agg:
        lines += ['| Aggregate | Value |', '|---|---|'] + agg + ['']

    # Per-kind columns: only show what the scenario actually measures.
    if kind == 'exploration':
        cols = [('Coverage', 'coverage_known_frac'), ('Clearance (m)', 'min_clearance_m'),
                ('Tilt (deg)', 'max_tilt_deg'), ('Home dist (m)', 'final_dist_to_home_m')]
    elif kind == 'goal':
        cols = [('Goal reached', 'goal_reached'), ('Min alt (m)', 'min_altitude_m'),
                ('Path (m)', 'path_length_m')]
    else:
        cols = [('Agents', 'agents_seen'), ('Swap legs', 'swap_legs'),
                ('Min sep (m)', 'min_separation_m')]

    header = ['Run', 'Verdict', 'Outcome'] + [c[0] for c in cols] + ['Dur (s)', 'Failed gates']
    lines += ['| ' + ' | '.join(header) + ' |',
              '|' + '---|' * len(header)]
    for r in records:
        m = r['metrics']
        bad = ', '.join(g for g, ok in r['gates'].items() if not ok) or '—'
        row = [f"{r['run']:02d}", 'PASS' if r['passed'] else 'FAIL', str(m.get('outcome'))]
        for _, key in cols:
            v = m.get(key)
            if isinstance(v, dict):          # swap_legs: show min..max
                v = f"{min(v.values())}..{max(v.values())}" if v else '—'
            row.append(str(v))
        row += [str(m.get('duration_s')), bad]
        lines.append('| ' + ' | '.join(row) + ' |')
    lines.append('')

    (out_dir / 'summary.md').write_text('\n'.join(lines))
    (out_dir / 'summary.json').write_text(json.dumps(
        dict(scenario=name, kind=kind, passed=npass, total=n, records=records),
        indent=2))
    print(f'--- {name}: {npass}/{n} passed (summary at {out_dir / "summary.md"})')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    for cmd in ('run', 'campaign'):
        p = sub.add_parser(cmd)
        p.add_argument('--scenario', default='ground-exploration')
        p.add_argument('--out', default='/tmp/mighty_campaign', type=Path)
        p.add_argument('-s', '--setup-bash', default=DEFAULT_SETUP_BASH, type=Path)
        p.add_argument('--ros-domain-id', type=int, default=30)
        if cmd == 'campaign':
            p.add_argument('--runs', type=int, default=10)
            p.add_argument('--max-hours', type=float, default=None,
                           help='Absolute campaign deadline; remaining runs are '
                                'reported as not-run rather than assumed passing.')
        else:
            p.add_argument('--index', type=int, default=1)
            p.add_argument('--attach', action='store_true',
                           help='Monitor a simulation started elsewhere (e.g. in '
                                'a Docker container) instead of launching one; '
                                'teardown is left to the caller.')
            p.add_argument('--log-file', default=None,
                           help='With --attach: read this file for planner log '
                                'markers instead of tmux panes.')
    sub.add_parser('scenarios')

    args = ap.parse_args()

    if args.cmd == 'scenarios':
        for k, v in sorted(SCENARIOS.items()):
            print(f"{k:24s} {v['description']}")
        return

    import_ros()
    os.environ['ROS_DOMAIN_ID'] = str(args.ros_domain_id)
    sc = scenario(args.scenario)
    args.out.mkdir(parents=True, exist_ok=True)

    if args.cmd == 'run':
        do_run(sc, args.scenario, args.out, args.setup_bash,
               args.ros_domain_id, args.index,
               attach=args.attach, log_file=args.log_file)
        return

    deadline = (time.time() + args.max_hours * 3600.0) if args.max_hours else None
    records, note = [], ''
    for i in range(1, args.runs + 1):
        if deadline and time.time() >= deadline:
            note = (f'> Campaign deadline hit after {len(records)} run(s); '
                    f'runs {i}-{args.runs} were **not run** (not passed).')
            print(f"[campaign] {note}", flush=True)
            break
        records.append(do_run(sc, args.scenario, args.out, args.setup_bash,
                             args.ros_domain_id, i, deadline_epoch=deadline))
        write_summary(args.out, args.scenario, records, note)
    write_summary(args.out, args.scenario, records, note)


if __name__ == '__main__':
    main()
