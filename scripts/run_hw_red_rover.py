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

Launches the MIGHTY planner on a Red Rover ground robot with either DLIO or
mocap localization.

Usage:
    # DLIO localization (default)
    python3 scripts/run_hw_red_rover.py

    # Mocap localization
    python3 scripts/run_hw_red_rover.py --odom-type mocap

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
MPC_CONFIG = MIGHTY_WS / 'src' / 'mpc' / 'config' / 'mpc.yaml'


def generate_yaml(odom_type: str, rover_name: str, goal_type: int) -> str:
    """Generate tmuxp YAML for hardware red rover."""

    source_ws = f'source {SETUP_BASH}'
    timestamp = datetime.now().strftime('%Y%m%d%H%M')

    # Mighty launch: mocap vs DLIO differ in use_onboard_localization and twist_topic
    if odom_type == 'mocap':
        mighty_cmd = (
            f'ros2 launch mighty onboard_mighty.launch.py'
            f' x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:={rover_name}'
            f' use_hardware:=true use_obstacle_tracker:=false'
            f' use_onboard_localization:=false robot_type:=red_rover'
            f' depth_camera_name:=d455 twist_topic:=mocap/twist'
        )
    else:
        mighty_cmd = (
            f'ros2 launch mighty onboard_mighty.launch.py'
            f' x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:={rover_name}'
            f' use_hardware:=true use_obstacle_tracker:=false'
            f' use_onboard_localization:=true robot_type:=red_rover'
            f' depth_camera_name:=d455'
        )

    # DLIO pane: replaced by static TFs when using mocap
    if odom_type == 'mocap':
        dlio_cmd = (
            f'ros2 run tf2_ros static_transform_publisher'
            f' --frame-id {rover_name} --child-frame-id {rover_name}/base_link'
            f' & ros2 run tf2_ros static_transform_publisher'
            f' --frame-id world --child-frame-id {rover_name}/map'
        )
    else:
        dlio_cmd = f'ros2 launch direct_lidar_inertial_odometry dlio.launch.py namespace:={rover_name}'

    # Global mapper: mocap uses pose_stamped on "world" topic, DLIO uses dlio/odom_node/pose.
    # `ground_robot:=true` plus `hardware:=true` picks hw_ground_robot.yaml from the launch default.
    if odom_type == 'mocap':
        mapper_cmd = (
            f'ros2 launch global_mapper_ros global_mapper_node.launch.py'
            f' hardware:=true ground_robot:=true quad:={rover_name}'
            f' depth_pointcloud_topic:=livox/lidar'
            f' pose_topic:=world pose_type:=pose_stamped'
        )
    else:
        mapper_cmd = (
            f'ros2 launch global_mapper_ros global_mapper_node.launch.py'
            f' hardware:=true ground_robot:=true quad:={rover_name}'
            f' depth_pointcloud_topic:=livox/lidar'
            f' pose_topic:=dlio/odom_node/pose'
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
        # Global Mapper
        {
            'shell_command': [
                source_ws,
                mapper_cmd,
            ]
        },
        # MPC Controller (mocap: override pose_topic from dlio/odom_node/pose to world)
        # {
        #     'shell_command': [
        #         source_ws,
        #         f'ros2 launch mpc mpc.launch.py namespace:={rover_name} hardware:=true params_file:={MPC_CONFIG}'
        #         if odom_type != 'mocap' else
        #         f'ros2 run mpc mpc_node --ros-args -r __ns:=/{rover_name}'
        #         f' --params-file {MPC_CONFIG}'
        #         f' -p cmd_vel_topic:=cmd_vel_auto'
        #         f' -p pose_topic:=world',
        #     ]
        # },
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
        choices=['dlio', 'mocap'],
        default='dlio',
        help='Localization source (default: dlio)',
    )

    parser.add_argument(
        '--goal-type', '-g',
        type=int,
        choices=[1, 2],
        default=1,
        help='Mocap goal pattern: 1 = (4,4)<->(-4,-4), 2 = (-4,4)<->(4,-4) (default: 1)',
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

    yaml_content = generate_yaml(args.odom_type, rover_name, args.goal_type)

    print(f'[INFO] Odom type: {args.odom_type}')
    print(f'[INFO] Rover: {rover_name}')
    if args.odom_type == 'mocap':
        print(f'[INFO] Goal type: {args.goal_type}')

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
