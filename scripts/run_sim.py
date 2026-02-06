#!/usr/bin/env python3
"""
MIGHTY Simulation Launcher

This script provides a unified interface to launch MIGHTY simulations in two modes:
1. Multi-agent simulation with fake sensing (fake_sim)
2. Single-agent simulation with Gazebo and ACL mapper (gazebo)

Usage:
    # Multi-agent fake simulation (10 agents in a circle) - auto-detects workspace
    python3 scripts/run_sim.py --mode multiagent

    # Single-agent UAV Gazebo simulation with default goal
    python3 scripts/run_sim.py --mode gazebo

    # Single-agent ground robot simulation (Pioneer 3-AT)
    python3 scripts/run_sim.py --mode gazebo --ground-robot

    # Ground robot with custom goal and environment
    python3 scripts/run_sim.py --mode gazebo --ground-robot --env easy_forest --goal 50 30 1

    # UAV with custom goal
    python3 scripts/run_sim.py --mode gazebo --goal 100 50 3

    # Custom number of agents for multiagent mode
    python3 scripts/run_sim.py --mode multiagent --num-agents 5

    # Custom environment for Gazebo mode
    python3 scripts/run_sim.py --mode gazebo --env easy_forest

    # Explicitly specify setup.bash if auto-detection fails
    python3 scripts/run_sim.py --mode gazebo --setup-bash /path/to/install/setup.bash
"""

import argparse
import math
import os
import subprocess
import sys
import tempfile
import yaml
from pathlib import Path


def find_setup_bash(args_setup_bash: str = None) -> Path:
    """Find setup.bash path. Auto-detects workspace if not specified."""
    if args_setup_bash:
        # User provided explicit path
        path = Path(args_setup_bash)
        if path.exists():
            return path
        print(f"[ERROR] Specified setup.bash not found: {args_setup_bash}", file=sys.stderr)
        sys.exit(1)

    # Auto-detect: try to find workspace root
    script_path = Path(__file__).resolve()
    # Assume script is in: <workspace>/src/mighty/scripts/run_sim.py
    workspace_root = script_path.parent.parent.parent.parent
    setup_bash = workspace_root / "install" / "setup.bash"

    if setup_bash.exists():
        print(f"[INFO] Auto-detected setup.bash at: {setup_bash}")
        return setup_bash

    print("[ERROR] Could not auto-detect setup.bash. Please specify with --setup-bash", file=sys.stderr)
    print(f"  Searched at: {setup_bash}", file=sys.stderr)
    print("  Example: python3 run_sim.py --mode gazebo --setup-bash /path/to/install/setup.bash", file=sys.stderr)
    sys.exit(1)




def generate_multiagent_positions(num_agents: int, radius: float = 10.0, z: float = 1.0):
    """Generate agent positions in a circle formation."""
    agents = []
    for i in range(num_agents):
        angle = 2 * math.pi * i / num_agents
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        # Yaw points toward center (opposite of position angle)
        yaw_deg = math.degrees(angle + math.pi)
        # Normalize to [-180, 180]
        if yaw_deg > 180:
            yaw_deg -= 360
        agents.append({
            'namespace': f'NX{i+1:02d}',
            'x': round(x, 3),
            'y': round(y, 3),
            'z': z,
            'yaw': round(yaw_deg, 1)
        })
    return agents


def generate_multiagent_yaml(setup_bash: Path, agents: list, sim_env: str, ros_domain_id: int = 7) -> str:
    """Generate YAML for multi-agent fake simulation."""
    panes = []

    # Base station (simulator)
    panes.append({
        'shell_command': [
            'ros2 launch mighty simulator.launch.py'
        ]
    })

    # Agent panes
    for agent in agents:
        panes.append({
            'shell_command': [
                'sleep 10',
                f"ros2 launch mighty onboard_mighty.launch.py namespace:={agent['namespace']} "
                f"x:={agent['x']} y:={agent['y']} z:={agent['z']} yaw:={agent['yaw']} sim_env:={sim_env}"
            ]
        })

    # Goal monitor
    panes.append({
        'shell_command': [
            'sleep 20',
            'ros2 launch mighty goal_monitor.launch.py'
        ]
    })

    yaml_content = {
        'session_name': 'mighty_sim',
        'windows': [{
            'window_name': 'main',
            'layout': 'tiled',
            'shell_command_before': [
                f'''if [ -z "$SETUP_BASH" ] || [ ! -f "$SETUP_BASH" ]; then
  echo "[ERROR] SETUP_BASH is missing or invalid: $SETUP_BASH" >&2
  exit 1
fi
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_gazebo_yaml(setup_bash: Path, goal: tuple, sim_env: str,
                         env: str = 'hard_forest',
                         start_pos: tuple = (0, 0, 3.0), start_yaw: float = 1.57,
                         ros_domain_id: int = 7, use_rviz: bool = True,
                         use_gazebo_gui: bool = False, use_ground_robot: bool = False) -> str:
    """Generate YAML for single-agent Gazebo simulation."""
    goal_x, goal_y, goal_z = goal
    start_x, start_y, start_z = start_pos

    panes = [
        # Base station with Gazebo
        {
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                f'ros2 launch mighty base_mighty.launch.py use_dyn_obs:=false '
                f'use_gazebo_gui:={str(use_gazebo_gui).lower()} use_rviz:={str(use_rviz).lower()} '
                f'env:={env} use_ground_robot:={str(use_ground_robot).lower()}'
            ]
        },
        # Ground robot odom-to-state converter (only for ground robot)
        # Converts /NX01/odom (from Gazebo diff_drive) to /NX01/state (for mapper and planner)
        {
            'shell_command': [
                'sleep 10',
                'ros2 run mighty convert_odom_to_state --ros-args -r __ns:=/NX01 -r odom:=odom -r state:=state'
            ] if use_ground_robot else ['echo "Skipping convert_odom_to_state (UAV mode)"']
        },
        # ACL mapper
        {
            'shell_command': [
                'sleep 10',
                f'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                f'param_file:={"global_mapper_ground_robot.yaml" if use_ground_robot else "global_mapper.yaml"}'
            ]
        },
        # Onboard agent NX01
        {
            'shell_command': [
                'sleep 10',
                f'ros2 launch mighty onboard_mighty.launch.py x:={start_x} y:={start_y} z:={start_z} yaw:={start_yaw} '
                f'sim_env:={sim_env} use_ground_robot:={str(use_ground_robot).lower()}'
            ]
        },
        # Goal sender
        {
            'shell_command': [
                'sleep 20',
                f"ros2 launch mighty goal_sender.launch.py list_agents:=\"['NX01']\" list_goals:=\"['[{goal_x}, {goal_y}, {goal_z}]']\""
            ]
        }
    ]

    yaml_content = {
        'session_name': 'mighty_sim',
        'windows': [{
            'window_name': 'main',
            'layout': 'tiled',
            'shell_command_before': [
                f'''if [ -z "$SETUP_BASH" ] || [ ! -f "$SETUP_BASH" ]; then
  echo "[ERROR] SETUP_BASH is missing or invalid: $SETUP_BASH" >&2
  exit 1
fi
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def main():
    parser = argparse.ArgumentParser(
        description='MIGHTY Simulation Launcher',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--mode', '-m',
        choices=['multiagent', 'gazebo'],
        required=True,
        help='Simulation mode: multiagent (fake sensing) or gazebo (single agent with ACL mapper)'
    )

    parser.add_argument(
        '--setup-bash', '-s',
        type=str,
        required=False,
        default=None,
        help='Path to setup.bash (required)'
    )

    parser.add_argument(
        '--goal', '-g',
        type=float,
        nargs=3,
        metavar=('X', 'Y', 'Z'),
        default=[105.0, 0.0, 3.0],
        help='Goal position for gazebo mode (default: 105.0 0.0 3.0)'
    )

    parser.add_argument(
        '--start', '-p',
        type=float,
        nargs=3,
        metavar=('X', 'Y', 'Z'),
        default=[0.0, 0.0, 3.0],
        help='Start position for gazebo mode (default: 0.0 0.0 3.0)'
    )

    parser.add_argument(
        '--start-yaw',
        type=float,
        default=1.57,
        help='Start yaw in radians for gazebo mode (default: 1.57)'
    )

    parser.add_argument(
        '--num-agents', '-n',
        type=int,
        default=10,
        help='Number of agents for multiagent mode (default: 10)'
    )

    parser.add_argument(
        '--radius', '-r',
        type=float,
        default=10.0,
        help='Circle radius for multiagent formation (default: 10.0)'
    )

    parser.add_argument(
        '--env', '-e',
        type=str,
        default='hard_forest',
        help='Gazebo environment (default: hard_forest)'
    )

    parser.add_argument(
        '--ros-domain-id',
        type=int,
        default=7,
        help='ROS_DOMAIN_ID (default: 7)'
    )

    parser.add_argument(
        '--rviz',
        action='store_true',
        default=True,
        help='Enable RViz (default: True)'
    )

    parser.add_argument(
        '--no-rviz',
        action='store_true',
        help='Disable RViz'
    )

    parser.add_argument(
        '--gazebo-gui',
        action='store_true',
        help='Enable Gazebo GUI (default: False)'
    )

    parser.add_argument(
        '--ground-robot',
        action='store_true',
        help='Use ground robot (Pioneer 3-AT) instead of UAV'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print the generated YAML without launching'
    )

    args = parser.parse_args()

    # Find setup.bash path
    setup_bash = find_setup_bash(args.setup_bash)
    print(f"[INFO] Using setup.bash: {setup_bash}")

    # Determine sim_env and generate YAML
    if args.mode == 'multiagent':
        sim_env = 'fake_sim'
        agents = generate_multiagent_positions(args.num_agents, args.radius)
        yaml_content = generate_multiagent_yaml(setup_bash, agents, sim_env, args.ros_domain_id)
        print(f"[INFO] Mode: Multi-agent simulation with {args.num_agents} agents (sim_env={sim_env})")
    else:  # gazebo
        sim_env = 'gazebo'
        use_rviz = args.rviz and not args.no_rviz

        # Determine if using ground robot
        use_ground_robot = args.ground_robot

        # Map environment names to world files
        env_to_world_mapping = {
            'ACL_office': 'ACL_office',
            'easy_forest': 'easy_forest',
            'hard_forest': 'hard_forest',
        }
        world_name = env_to_world_mapping.get(args.env, args.env)

        # Adjust start position z for ground robot (ground level vs flying)
        start_pos = list(args.start)
        if use_ground_robot and start_pos[2] == 3.0:  # Only adjust if using default z
            start_pos[2] = 0.0  # Ground robot base_link at z=0 (wheels at ground)
        start_pos = tuple(start_pos)

        yaml_content = generate_gazebo_yaml(
            setup_bash,
            goal=tuple(args.goal),
            sim_env=sim_env,
            env=world_name,
            start_pos=start_pos,
            start_yaw=args.start_yaw,
            ros_domain_id=args.ros_domain_id,
            use_rviz=use_rviz,
            use_gazebo_gui=args.gazebo_gui,
            use_ground_robot=use_ground_robot
        )
        print(f"[INFO] Mode: Single-agent Gazebo simulation (sim_env={sim_env})")
        print(f"[INFO] Environment: {args.env} (world: {world_name})")
        if use_ground_robot:
            print(f"[INFO] Vehicle: Ground robot (Pioneer 3-AT)")
        print(f"[INFO] Start: ({start_pos[0]}, {start_pos[1]}, {start_pos[2]})")
        print(f"[INFO] Goal: ({args.goal[0]}, {args.goal[1]}, {args.goal[2]})")

    if args.dry_run:
        print("\n[DRY RUN] Generated YAML:")
        print("-" * 60)
        print(yaml_content)
        print("-" * 60)
        return

    # Write temporary YAML file and launch with tmuxp
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        temp_yaml_path = f.name

    try:
        print(f"[INFO] Launching simulation...")
        env = os.environ.copy()
        env['SETUP_BASH'] = str(setup_bash)
        subprocess.run(['tmuxp', 'load', temp_yaml_path], env=env, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to launch simulation: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("[ERROR] tmuxp not found. Install with: pip install tmuxp", file=sys.stderr)
        sys.exit(1)
    finally:
        # Clean up temp file
        os.unlink(temp_yaml_path)


if __name__ == '__main__':
    main()
