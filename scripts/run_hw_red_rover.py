#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

"""
Hardware Red Rover Launcher

Launches the MIGHTY planner on a Red Rover ground robot with DLIO, mocap, or
mocap-seeded-DLIO localization.

Usage:
    # DLIO localization (starts at origin)
    python3 scripts/run_hw_red_rover.py

    # Mocap localization (state from /world topic)
    python3 scripts/run_hw_red_rover.py --odom-type mocap

    # DLIO seeded with mocap pose: DLIO subscribes to /<ns>/world once at
    # startup and anchors its odom frame to that pose, so its outputs (and the
    # downstream map) live in mocap world coords from the first scan.
    python3 scripts/run_hw_red_rover.py --odom-type dlio_in_mocap

    # Mocap with diagonal goal type 2
    python3 scripts/run_hw_red_rover.py --odom-type mocap --goal-type 2

    # Preview generated YAML without launching
    python3 scripts/run_hw_red_rover.py --odom-type mocap --dry-run
"""

import argparse
import os
import subprocess
import sys
import tempfile
from datetime import datetime
import yaml
from pathlib import Path


# Paths
MIGHTY_WS = Path('/home/swarm/code/mighty_ws')
SETUP_BASH = MIGHTY_WS / 'install' / 'setup.bash'
DECOMP_SETUP_BASH = Path('/home/swarm/code/decomp_ws/install/setup.bash')


def generate_yaml(odom_type: str, rover_name: str, goal_type: int,
                  two_d_only: bool = False) -> str:
    """Generate tmuxp YAML for hardware red rover."""

    source_ws = f'source {SETUP_BASH}'
    timestamp = datetime.now().strftime('%Y%m%d%H%M')

    # Mighty launch: mocap vs DLIO differ in use_onboard_localization and twist_topic.
    # `dlio_in_mocap` is a DLIO setup with one knob flipped (the DLIO seed),
    # so the mighty side is identical to plain `dlio`.
    if odom_type == 'mocap':
        mighty_cmd = (
            f'ros2 launch mighty onboard_mighty.launch.py'
            f' x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:={rover_name}'
            f' use_hardware:=true'
            f' use_onboard_localization:=false robot_type:=red_rover'
            f' depth_camera_name:=d455 twist_topic:=mocap/twist'
        )
    else:  # 'dlio' or 'dlio_in_mocap'
        mighty_cmd = (
            f'ros2 launch mighty onboard_mighty.launch.py'
            f' x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:={rover_name}'
            f' use_hardware:=true'
            f' use_onboard_localization:=true robot_type:=red_rover'
            f' depth_camera_name:=d455'
        )

    # DLIO pane: replaced by static TFs when using mocap; in `dlio_in_mocap`
    # mode DLIO runs as usual but is told to seed its initial pose from the
    # mocap PoseStamped on `<ns>/world`, so its odom frame is anchored to the
    # mocap world from the first scan.
    if odom_type == 'mocap':
        dlio_cmd = (
            f'ros2 run tf2_ros static_transform_publisher'
            f' --frame-id {rover_name} --child-frame-id {rover_name}/base_link'
            f' & ros2 run tf2_ros static_transform_publisher'
            f' --frame-id world --child-frame-id {rover_name}/map'
        )
    elif odom_type == 'dlio_in_mocap':
        dlio_cmd = (
            f'ros2 launch direct_lidar_inertial_odometry dlio.launch.py'
            f' namespace:={rover_name} initial_pose_topic:=world'
            f' two_d_only:={"true" if two_d_only else "false"}'
        )
    else:
        dlio_cmd = (
            f'ros2 launch direct_lidar_inertial_odometry dlio.launch.py'
            f' namespace:={rover_name}'
            f' two_d_only:={"true" if two_d_only else "false"}'
        )

    # Global mapper: mocap uses pose_stamped on "world", DLIO (and dlio_in_mocap)
    # use dlio/odom_node/pose. Same mapper config in both DLIO modes — the
    # difference is solely whether DLIO's frame is anchored to mocap or to (0,0,0).
    if odom_type == 'mocap':
        mapper_cmd = (
            f'ros2 launch global_mapper_ros global_mapper_node.launch.py'
            f' hardware:=true ground_robot:=true quad:={rover_name}'
            f' depth_pointcloud_topic:=livox/lidar'
            f' pose_topic:=world pose_type:=pose_stamped'
            f' use_obstacle_tracker:=false'
        )
    else:  # 'dlio' or 'dlio_in_mocap'
        mapper_cmd = (
            f'ros2 launch global_mapper_ros global_mapper_node.launch.py'
            f' hardware:=true ground_robot:=true quad:={rover_name}'
            f' depth_pointcloud_topic:=livox/lidar'
            f' pose_topic:=dlio/odom_node/pose'
            f' use_obstacle_tracker:=false'
        )

    # Goal monitor: odom_type and goal_type are now baked in directly
    goal_monitor_cmd = (
        f'ros2 run mighty goal_monitor_node.py --ros-args'
        f' -r __ns:=/{rover_name}'
        f' -p use_hardware:=true -p use_ground_robot:=true'
        f' -p odom_type:={odom_type} -p goal_type:={goal_type}'
        f' -p goal_tolerance:=1.0 -p num_agents:=1 -p radius:=10.0'
    )

    panes = [
        # Onboard mighty
        {
            'shell_command': [
                source_ws,
                f'source {DECOMP_SETUP_BASH}',
                mighty_cmd,
            ]
        },
        # Livox LiDAR
        {
            'shell_command': [
                source_ws,
                f'ros2 launch livox_ros_driver2 run_MID360_launch.py namespace:={rover_name}',
            ]
        },
        # DLIO (dlio mode) or static TFs (mocap mode)
        {
            'shell_command': [
                source_ws,
                'sleep 5',
                dlio_cmd,
            ]
        },
        # Republish rviz 2D goal to term goal
        {
            'shell_command': [
                source_ws,
                'ros2 run mighty repub_rviz_2Dgoal.py',
            ]
        },
        # Static TF (lidar tilt)
        {
            'shell_command': [
                source_ws,
                'sleep 5',
                f'ros2 run tf2_ros static_transform_publisher 0 0 0 0 0.3490659 0 {rover_name}/base_link {rover_name}/lidar',
            ]
        },
        # Static TF (odom to map)
        {
            'shell_command': [
                source_ws,
                'sleep 5',
                f'ros2 run tf2_ros static_transform_publisher 0 0 0 0 0 0 {rover_name}/map {rover_name}/odom',
            ]
        },
        # Global Mapper
        {
            'shell_command': [
                source_ws,
                mapper_cmd,
            ]
        },
        # MPC controller is spawned by onboard_mighty.launch.py (hw + red_rover branch)
        # using hw_mighty_ground_robot.yaml + mpc.yaml, publishing to cmd_vel_auto.
        # Goal Monitor
        # {
        #     'shell_command': [
        #         source_ws,
        #         goal_monitor_cmd,
        #     ]
        # },
        # Bag recorder
        # {
        #     'shell_command': [
        #         source_ws,
        #         f'echo "odom_type={odom_type}  rover={rover_name}  goal_type={goal_type}"',
        #         f'python3 {MIGHTY_WS / "src" / "mighty" / "scripts" / "hw_bag_record.py"}'
        #         f' --bag_path /home/swarm/data/multi_mighty'
        #         f' --agents {rover_name}'
        #         f' --bag_name {rover_name}_{timestamp}',
        #     ]
        # },
    ]

    yaml_content = {
        'session_name': 'hw_mighty',
        'windows': [{
            'window_name': 'main',
            'layout': 'tiled',
            'shell_command_before': [
                f'source /opt/ros/humble/setup.bash',
            ],
            'panes': panes,
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def main():
    parser = argparse.ArgumentParser(
        description='Hardware Red Rover Launcher',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        '--odom-type', '-o',
        choices=['dlio', 'mocap', 'dlio_in_mocap'],
        default='dlio',
        help='Localization source (default: dlio). dlio_in_mocap = DLIO seeds '
             'its initial pose from one mocap PoseStamped on /<ns>/world, so '
             'its odom frame is anchored to mocap world from the first scan.',
    )

    parser.add_argument(
        '--goal-type', '-g',
        type=int,
        choices=[1, 2],
        default=1,
        help='Mocap goal pattern: 1 = (4,4)<->(-4,-4), 2 = (-4,4)<->(4,-4) (default: 1)',
    )

    parser.add_argument(
        '--two-d-only',
        action='store_true',
        help='Tell DLIO to overwrite its published z with a constant (pinned to '
             'the mocap seed z when --odom-type=dlio_in_mocap, else averaged '
             'from the first samples). For flat 2D environments.',
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print generated YAML without launching',
    )

    args = parser.parse_args()

    # Read rover name from environment (set in .bashrc)
    rover_name = os.environ.get('ROVER_NAME')
    if not rover_name:
        print('[ERROR] ROVER_NAME not set. Check your .bashrc (VEHTYPE/VEHNUM).', file=sys.stderr)
        sys.exit(1)

    yaml_content = generate_yaml(args.odom_type, rover_name, args.goal_type,
                                 two_d_only=args.two_d_only)

    print(f'[INFO] Odom type: {args.odom_type}')
    print(f'[INFO] Rover: {rover_name}')
    if args.odom_type == 'mocap':
        print(f'[INFO] Goal type: {args.goal_type}')
    if args.two_d_only and args.odom_type != 'mocap':
        print('[INFO] DLIO 2D-only output: ON (z pinned to constant)')

    if args.dry_run:
        print('\n[DRY RUN] Generated YAML:')
        print('-' * 60)
        print(yaml_content)
        print('-' * 60)
        return

    # Kill any existing session with this name
    subprocess.run(['tmux', 'kill-session', '-t', 'hw_mighty'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Write temp YAML and launch
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        temp_path = f.name

    try:
        print('[INFO] Launching hardware red rover...')
        subprocess.run(['tmuxp', 'load', '-d', temp_path], check=True)

        # Attach to the session
        subprocess.run(['tmux', 'attach-session', '-t', 'hw_mighty'])
    except subprocess.CalledProcessError as e:
        print(f'[ERROR] Failed to launch: {e}', file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print('[ERROR] tmuxp not found. Install with: pip install tmuxp', file=sys.stderr)
        sys.exit(1)
    finally:
        os.unlink(temp_path)


if __name__ == '__main__':
    main()
