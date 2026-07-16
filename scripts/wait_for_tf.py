#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""Block until a set of TF transforms are available, then exit 0.

A startup gate for launch/tmux panes: run this BEFORE a node that depends on
transforms which may not be published yet. It fixes the rmw_zenoh /tf_static
race where a consumer that subscribes before a one-shot latched static
transform is published never receives it -- a fresh late-joining listener
(this script, and then the gated node) queries and receives it reliably.

Usage:
    wait_for_tf.py TARGET1 SOURCE1 [TARGET2 SOURCE2 ...] [--timeout S] [--poll S]

Each TARGET SOURCE pair is checked with tf2 can_transform(target, source, latest).
Exits 0 once ALL pairs are available. On timeout it logs a warning and STILL
exits 0, so the launch degrades (the node prints its own tf warnings) rather
than deadlocking the whole stack.
"""

import argparse
import sys

import rclpy
from rclpy.duration import Duration
from rclpy.node import Node
from rclpy.time import Time
from tf2_ros import Buffer, TransformListener


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('frames', nargs='+',
                    help='flat list of TARGET SOURCE pairs (even count)')
    ap.add_argument('--timeout', type=float, default=60.0,
                    help='seconds to wait before proceeding anyway (default: 60)')
    ap.add_argument('--poll', type=float, default=0.5,
                    help='seconds between availability checks (default: 0.5)')
    args = ap.parse_args()

    if len(args.frames) % 2 != 0:
        print('wait_for_tf: frames must be TARGET SOURCE pairs (even count)',
              file=sys.stderr)
        sys.exit(2)
    pairs = list(zip(args.frames[0::2], args.frames[1::2]))

    rclpy.init()
    # Global namespace on purpose: the TransformListener then uses the global
    # /tf and /tf_static (where the static_transform_publishers publish), exactly
    # like a plain `ros2 run tf2_ros tf2_echo`.
    node = Node('wait_for_tf')
    buf = Buffer()
    TransformListener(buf, node)
    node.get_logger().info(
        f'waiting for transforms {pairs} (timeout {args.timeout}s)')

    start = node.get_clock().now()
    timeout = Duration(seconds=args.timeout)
    remaining = list(pairs)
    try:
        while rclpy.ok():
            rclpy.spin_once(node, timeout_sec=args.poll)
            remaining = [(t, s) for (t, s) in remaining
                         if not buf.can_transform(t, s, Time())]
            if not remaining:
                node.get_logger().info('all transforms available; proceeding')
                return
            if node.get_clock().now() - start > timeout:
                node.get_logger().warn(
                    f'TIMEOUT after {args.timeout}s; still missing {remaining}; '
                    f'proceeding anyway (degraded)')
                return
    finally:
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()


if __name__ == '__main__':
    main()
