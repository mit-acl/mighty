#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""
Remote MIGHTY Launcher

Opens a tmux session with one pane per agent, each running mighty via SSH.

Usage:
    python3 scripts/run_mighty_remote.py RR03 RR04 RR05
"""

import argparse
import subprocess
import sys


SESSION_NAME = 'mighty_remote'


def main():
    parser = argparse.ArgumentParser(
        description='Launch mighty on remote agents via SSH in tmux panes',
    )
    parser.add_argument(
        'agents',
        nargs='+',
        help='Agent names to SSH into (e.g., RR03 RR04 RR05)',
    )
    args = parser.parse_args()

    # Kill any existing session
    subprocess.run(
        ['tmux', 'kill-session', '-t', SESSION_NAME],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Create a new detached session with the first agent
    first = args.agents[0]
    subprocess.run(
        ['tmux', 'new-session', '-d', '-s', SESSION_NAME, f'ssh {first} "mighty"'],
        check=True,
    )

    # Split a new pane for each remaining agent
    for agent in args.agents[1:]:
        subprocess.run(
            ['tmux', 'split-window', '-t', SESSION_NAME, f'ssh {agent} "mighty"'],
            check=True,
        )
        # Re-tile after each split so panes stay evenly distributed
        subprocess.run(
            ['tmux', 'select-layout', '-t', SESSION_NAME, 'tiled'],
            check=True,
        )

    # Attach to the session
    subprocess.run(['tmux', 'attach-session', '-t', SESSION_NAME])


if __name__ == '__main__':
    main()
