#!/usr/bin/env python3

# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

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
import re
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


def find_workspace_root() -> Path:
    """Find workspace root relative to this script."""
    # Script is in: <workspace>/src/mighty/scripts/run_sim.py
    return Path(__file__).resolve().parent.parent.parent.parent


def exploration_enabled_in_yaml(config_path: Path) -> bool:
    """Return True if `exploration.enabled: true` appears in the given YAML.
    Tolerant to the dotted-key form used by ROS 2 parameter files. Used to
    auto-skip goal_sender when frontier-based exploration would drive the
    robot itself.
    """
    if not config_path.exists():
        return False
    try:
        with open(config_path, 'r') as f:
            for raw in f:
                line = raw.split('#', 1)[0].strip()
                if not line:
                    continue
                # Match `exploration.enabled: true` (with optional whitespace).
                if line.startswith('exploration.enabled'):
                    _, _, val = line.partition(':')
                    return val.strip().lower() == 'true'
    except OSError:
        pass
    return False


def find_rviz_config(name: str = 'mighty.rviz') -> Path:
    """Find the RViz config in the source tree (relative to this script)."""
    script_path = Path(__file__).resolve()
    # Script is in: <package>/scripts/run_sim.py, rviz is in: <package>/rviz/<name>
    return script_path.parent.parent / 'rviz' / name


def generate_multiagent_positions(num_agents: int, radius: float = 10.0, z: float = 1.0, prefix: str = 'NX', angle_offset: float = 0.0):
    """Generate agent positions in a circle formation."""
    agents = []
    for i in range(num_agents):
        angle = 2 * math.pi * i / num_agents + angle_offset
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        # Yaw points toward center (opposite of position angle)
        yaw_deg = math.degrees(angle + math.pi)
        # Normalize to [-180, 180]
        if yaw_deg > 180:
            yaw_deg -= 360
        agents.append({
            'namespace': f'{prefix}{i+1:02d}',
            'x': round(x, 3),
            'y': round(y, 3),
            'z': z,
            'yaw': round(yaw_deg, 1)
        })
    return agents


def generate_multiagent_yaml(setup_bash: Path, agents: list, sim_env: str, ros_domain_id: int = 20, radius: float = 10.0, no_goal: bool = False, rviz_config: Path = None, use_ground_robot: bool = False, agent_prefix: str = 'NX', config_file: Path = None) -> str:
    """Generate YAML for multi-agent fake simulation."""
    panes = []

    # Base station (simulator)
    sim_cmd = 'ros2 launch mighty simulator.launch.py'
    if rviz_config:
        sim_cmd += f' rviz_config:={rviz_config}'
    panes.append({
        'shell_command': [sim_cmd]
    })

    # Agent panes
    ground_robot_flag = f' use_ground_robot:={str(use_ground_robot).lower()}'
    config_flag = f' config_file:={config_file}' if config_file else ''
    for agent in agents:
        panes.append({
            'shell_command': [
                'sleep 10',
                f"ros2 launch mighty onboard_mighty.launch.py namespace:={agent['namespace']} "
                f"x:={agent['x']} y:={agent['y']} z:={agent['z']} yaw:={agent['yaw']} sim_env:={sim_env}"
                f"{ground_robot_flag}{config_flag}"
            ]
        })

    # Goal monitor
    if not no_goal:
        num_agents = len(agents)
        panes.append({
            'shell_command': [
                'sleep 20',
                f'ros2 launch mighty goal_monitor.launch.py num_agents:={num_agents} radius:={radius} agent_prefix:={agent_prefix} use_ground_robot:={str(use_ground_robot).lower()}'
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_interactive_yaml(setup_bash: Path, ros_domain_id: int = 20, rviz_config: Path = None) -> str:
    """Generate YAML for single-agent interactive simulation (click goals in RViz)."""
    sim_cmd = 'ros2 launch mighty simulator.launch.py'
    if rviz_config:
        sim_cmd += f' rviz_config:={rviz_config}'
    panes = [
        # Base station (random forest map + RViz)
        {
            'shell_command': [sim_cmd]
        },
        # Single agent NX01 at center
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch mighty onboard_mighty.launch.py namespace:=NX01 '
                'x:=0.0 y:=0.0 z:=1.0 yaw:=0.0 sim_env:=fake_sim'
            ]
        },
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_multiagent_ground_yaml(setup_bash: Path, agents: list, radius: float,
                                    ros_domain_id: int = 20) -> str:
    """Generate YAML for multi-agent ground robot simulation in Gazebo with MPC."""
    panes = []

    # Base station: Gazebo with ground_robot_forest world + RViz (no dynamic obstacles)
    panes.append({
        'shell_command': [
            'source /usr/share/gazebo/setup.bash',
            'ros2 launch mighty base_mighty.launch.py '
            'use_gazebo_gui:=false use_rviz:=true env:=ground_robot_forest use_ground_robot:=true'
        ]
    })

    # Per-agent: odom converter + ACL mapper + mighty (MPC) + MPC controller
    for i, agent in enumerate(agents):
        ns = agent['namespace']

        # Odom-to-state converter
        panes.append({
            'shell_command': [
                'sleep 10',
                f'ros2 run mighty convert_odom_to_state --ros-args -r __ns:=/{ns} -r odom:=odom -r state:=state'
            ]
        })

        # ACL mapper (obstacle tracker disabled for static env)
        panes.append({
            'shell_command': [
                'sleep 10',
                f'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                f'use_obstacle_tracker:=false param_file:=sim_ground_robot.yaml quad:={ns}'
            ]
        })

        # Mighty planner with trajectory tracker
        panes.append({
            'shell_command': [
                'sleep 12',
                f"ros2 launch mighty onboard_mighty.launch.py namespace:={ns} "
                f"x:={agent['x']} y:={agent['y']} z:={agent['z']} yaw:={agent['yaw']} "
                f"sim_env:=gazebo use_ground_robot:=true "
                f"num_agents:={len(agents)}"
            ]
        })

    # Goal monitor (swap pattern)
    num_agents = len(agents)
    panes.append({
        'shell_command': [
            'sleep 25',
            f'ros2 launch mighty goal_monitor.launch.py num_agents:={num_agents} '
            f'radius:={radius} agent_prefix:=NX goal_tolerance:=1.0 use_ground_robot:=true'
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_exploration_multiagent_ground_yaml(
        setup_bash: Path, agents: list, ros_domain_id: int = 20,
        rviz_config: Path = None, sim_env: str = 'fake_sim',
        env: str = 'ACL_office', use_vlm: bool = False,
        use_follow: bool = False) -> str:
    """Generate YAML for multi-agent ground robot exploration.

    Each agent runs:  onboard_mighty (ground robot, exploration enabled)
                    + global_mapper  (2D occ/ESDF for frontier detection)
                    + convert_odom_to_state (Gazebo only)
    No goal monitor — frontier-based exploration is self-driven.
    MinPos + visited-map sharing coordinate the agents.
    """
    panes = []
    use_gazebo = (sim_env == 'gazebo')

    # Base station
    if use_gazebo:
        panes.append({
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                f'ros2 launch mighty base_mighty.launch.py '
                f'use_gazebo_gui:=false use_rviz:=true env:={env} use_ground_robot:=true'
            ]
        })
    else:
        sim_cmd = 'ros2 launch mighty simulator.launch.py'
        if rviz_config:
            sim_cmd += f' rviz_config:={rviz_config}'
        panes.append({
            'shell_command': [sim_cmd]
        })

    # Per-agent nodes
    for i, agent in enumerate(agents):
        ns = agent['namespace']
        delay = 10 + i * 2  # stagger startup to avoid resource spikes

        # Gazebo: odom-to-state converter (Gazebo publishes odom, mighty needs state)
        if use_gazebo:
            panes.append({
                'shell_command': [
                    f'sleep {delay}',
                    f'ros2 run mighty convert_odom_to_state '
                    f'--ros-args -r __ns:=/{ns} -r odom:=odom -r state:=state'
                ]
            })

        # ACL mapper (provides occ_2d, esdf_2d for frontier detection)
        gazebo_flag = ' use_gazebo:=true' if use_gazebo else ' hardware:=false'
        panes.append({
            'shell_command': [
                f'sleep {delay}',
                f'ros2 launch global_mapper_ros global_mapper_node.launch.py'
                f'{gazebo_flag} ground_robot:=true '
                f'param_file:=sim_ground_robot.yaml quad:={ns} '
                # NOTE: global_mapper_node.launch.py defaults
                # use_obstacle_tracker:=true and OVERRIDES the yaml's
                # obstacle_tracker.enabled / use_temporal_grid values.
                # Forcing it false here so the yaml's `enabled: false`
                # actually takes effect — no /tracked_obstacles publish.
                f'use_obstacle_tracker:=false'
            ]
        })

        # Mighty planner (ground robot, exploration + MinPos enabled via config)
        external_selector_arg = ' external_selector:=true' if use_vlm else ''
        panes.append({
            'shell_command': [
                f'sleep {delay + 2}',
                f"ros2 launch mighty onboard_mighty.launch.py namespace:={ns} "
                f"x:={agent['x']} y:={agent['y']} z:={agent['z']} yaw:={agent['yaw']} "
                f"sim_env:={sim_env} use_ground_robot:=true "
                f"num_agents:={len(agents)}{external_selector_arg}"
            ]
        })

    if use_vlm:
        # VLM selector + chat CLI. Single-agent only (the VLM node defaults to
        # NX01 and is not multi-agent-aware in this iteration).
        ns = agents[0]['namespace']
        enable_follow = 'true' if use_follow else 'false'
        panes.append({
            'shell_command': [
                f'sleep {15 + len(agents) * 2}',
                f'ros2 launch vlm_goal_selector vlm_selector.launch.py '
                f'robot_namespace:={ns} enable_follow:={enable_follow}',
            ]
        })
        panes.append({
            'shell_command': [
                f'sleep {17 + len(agents) * 2}',
                'ros2 run vlm_goal_selector chat_cli',
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_swap_multiagent_ground_yaml(
        setup_bash: Path, agents: list, radius: float, angle_offset: float,
        ros_domain_id: int = 20, env: str = 'ACL_office') -> str:
    """Generate YAML for multi-agent ground robot position swapping in Gazebo.

    Each agent swaps with its diametrically opposite peer on a circle.
    No exploration — uses goal_monitor for continuous swap behavior.
    """
    panes = []

    # Base station: Gazebo world + RViz
    panes.append({
        'shell_command': [
            'source /usr/share/gazebo/setup.bash',
            f'ros2 launch mighty base_mighty.launch.py '
            f'use_gazebo_gui:=false use_rviz:=true env:={env} use_ground_robot:=true'
        ]
    })

    # Per-agent: odom converter + ACL mapper + mighty
    for i, agent in enumerate(agents):
        ns = agent['namespace']
        delay = 10 + i * 2

        panes.append({
            'shell_command': [
                f'sleep {delay}',
                f'ros2 run mighty convert_odom_to_state '
                f'--ros-args -r __ns:=/{ns} -r odom:=odom -r state:=state'
            ]
        })

        panes.append({
            'shell_command': [
                f'sleep {delay}',
                f'ros2 launch global_mapper_ros global_mapper_node.launch.py '
                f'use_gazebo:=true use_obstacle_tracker:=false '
                f'param_file:=sim_ground_robot.yaml quad:={ns}'
            ]
        })

        panes.append({
            'shell_command': [
                f'sleep {delay + 2}',
                f"ros2 launch mighty onboard_mighty.launch.py namespace:={ns} "
                f"x:={agent['x']} y:={agent['y']} z:={agent['z']} yaw:={agent['yaw']} "
                f"sim_env:=gazebo use_ground_robot:=true "
                f"num_agents:={len(agents)}"
            ]
        })

    # Goal monitor for position swapping
    num_agents = len(agents)
    panes.append({
        'shell_command': [
            'sleep 25',
            f'ros2 launch mighty goal_monitor.launch.py num_agents:={num_agents} '
            f'radius:={radius} angle_offset:={angle_offset} '
            f'agent_prefix:=NX goal_tolerance:=1.0 use_ground_robot:=true'
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_dyn_test_yaml(setup_bash: Path, ros_domain_id: int = 7) -> str:
    """Generate YAML for dynamic obstacle test: one drone + one dyn obstacle in Gazebo."""
    panes = [
        # Base station with Gazebo + 1 dynamic obstacle
        {
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                'ros2 launch mighty base_mighty.launch.py use_dyn_obs:=true '
                'num_dyn_obstacles:=1 dyn_x_min:=3.0 dyn_x_max:=3.0 dyn_y_min:=0.0 dyn_y_max:=0.0 '
                'use_gazebo_gui:=false use_rviz:=true env:=empty'
            ]
        },
        # ACL mapper (with obstacle tracker enabled for dynamic obstacle test)
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                'use_obstacle_tracker:=true param_file:=sim_uav.yaml'
            ]
        },
        # Onboard agent NX01 — stationary, no goal sent
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch mighty onboard_mighty.launch.py x:=0.0 y:=0.0 z:=3.0 yaw:=0.0 sim_env:=gazebo'
            ]
        },
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_dyn_test_ground_yaml(setup_bash: Path, ros_domain_id: int = 7) -> str:
    """Generate YAML for ground robot + dynamic obstacle test in Gazebo."""
    panes = [
        # Base station with Gazebo (static obstacles only)
        {
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                'ros2 launch mighty base_mighty.launch.py '
                'use_gazebo_gui:=false use_rviz:=true env:=ground_robot_forest use_ground_robot:=true'
            ]
        },
        # Odom-to-state converter for ground robot
        {
            'shell_command': [
                'sleep 10',
                'ros2 run mighty convert_odom_to_state --ros-args -r __ns:=/NX01 -r odom:=odom -r state:=state'
            ]
        },
        # ACL mapper (obstacle tracker disabled for static environment)
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                'use_obstacle_tracker:=false param_file:=sim_ground_robot.yaml'
            ]
        },
        # Onboard agent NX01 — ground robot, no goal sent
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch mighty onboard_mighty.launch.py x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 '
                'sim_env:=gazebo use_ground_robot:=true'
            ]
        },
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_dyn_test_ground_mpc_yaml(setup_bash: Path, ros_domain_id: int = 7) -> str:
    """Generate YAML for ground robot + MPC + static obstacle test in Gazebo."""
    panes = [
        # Base station with Gazebo + 1 dynamic obstacle
        {
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                'ros2 launch mighty base_mighty.launch.py use_dyn_obs:=true '
                'num_dyn_obstacles:=1 dyn_x_min:=3.0 dyn_x_max:=3.0 dyn_y_min:=0.0 dyn_y_max:=0.0 '
                'dyn_z_min:=0.3 dyn_z_max:=0.3 dyn_scale_z_min:=0.0 dyn_scale_z_max:=0.0 '
                'use_gazebo_gui:=false use_rviz:=true env:=ground_robot_forest use_ground_robot:=true'
            ]
        },
        # Odom-to-state converter for ground robot
        {
            'shell_command': [
                'sleep 10',
                'ros2 run mighty convert_odom_to_state --ros-args -r __ns:=/NX01 -r odom:=odom -r state:=state'
            ]
        },
        # ACL mapper (obstacle tracker enabled for dynamic obstacle)
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                'use_obstacle_tracker:=true param_file:=sim_ground_robot.yaml'
            ]
        },
        # Onboard agent NX01 — ground robot (onboard_mighty.launch.py spawns the MPC controller)
        {
            'shell_command': [
                'sleep 10',
                'ros2 launch mighty onboard_mighty.launch.py x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 '
                'sim_env:=gazebo use_ground_robot:=true'
            ]
        },
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
        }]
    }

    return yaml.dump(yaml_content, default_flow_style=False, sort_keys=False)


def generate_gazebo_yaml(setup_bash: Path, goal: tuple, sim_env: str,
                         env: str = 'hard_forest',
                         start_pos: tuple = (0, 0, 3.0), start_yaw: float = 1.57,
                         ros_domain_id: int = 7, use_rviz: bool = True,
                         use_gazebo_gui: bool = False, use_ground_robot: bool = False,
                         no_goal: bool = False) -> str:
    """Generate YAML for single-agent Gazebo simulation."""
    goal_x, goal_y, goal_z = goal
    start_x, start_y, start_z = start_pos

    panes = [
        # Base station with Gazebo
        {
            'shell_command': [
                'source /usr/share/gazebo/setup.bash',
                f'ros2 launch mighty base_mighty.launch.py '
                f'use_gazebo_gui:={str(use_gazebo_gui).lower()} use_rviz:={str(use_rviz).lower()} '
                f'env:={env} use_ground_robot:={str(use_ground_robot).lower()}'
            ]
        },
        # Ground robot odom-to-state converter (only for ground robot)
        # Converts /NX01/odom (from Gazebo diff_drive) to /NX01/state (for mapper and planner)
        {
            'shell_command': [
                'sleep 3',
                'ros2 run mighty convert_odom_to_state --ros-args -r __ns:=/NX01 -r odom:=odom -r state:=state'
            ] if use_ground_robot else ['echo "Skipping convert_odom_to_state (UAV mode)"']
        },
        # ACL mapper
        {
            'shell_command': [
                'sleep 3',
                f'ros2 launch global_mapper_ros global_mapper_node.launch.py use_gazebo:=true '
                f'use_obstacle_tracker:=false '
                f'param_file:={"sim_ground_robot.yaml" if use_ground_robot else "sim_uav.yaml"}'
            ]
        },
        # Onboard agent NX01
        {
            'shell_command': [
                'sleep 3',
                f'ros2 launch mighty onboard_mighty.launch.py x:={start_x} y:={start_y} z:={start_z} yaw:={start_yaw} '
                f'sim_env:={sim_env} use_ground_robot:={str(use_ground_robot).lower()}'
            ]
        },
    ]

    if not no_goal:
        # Goal sender
        panes.append({
            'shell_command': [
                'sleep 20',
                f"ros2 launch mighty goal_sender.launch.py list_agents:=\"['NX01']\" list_goals:=\"['[{goal_x}, {goal_y}, {goal_z}]']\""
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
unset AMENT_PREFIX_PATH COLCON_PREFIX_PATH CMAKE_PREFIX_PATH PYTHONPATH LD_LIBRARY_PATH
. "$SETUP_BASH"''',
                f'export ROS_DOMAIN_ID={ros_domain_id}'
            ],
            'panes': panes + [{'shell_command': ['# free pane — type commands here, e.g. `ros2 topic echo /NX01/term_goal`']}]
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
        choices=['multiagent', 'multiagent-ground', 'exploration-singleagent-ground', 'exploration-multiagent-ground', 'swap-multiagent-ground', 'gazebo', 'interactive', 'dyn-test', 'dyn-test-ground', 'dyn-test-ground-mpc'],
        required=True,
        help='Simulation mode: multiagent, exploration-singleagent-ground, exploration-multiagent-ground, swap-multiagent-ground, gazebo, interactive, dyn-test, dyn-test-ground'
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
        default=20,
        help='ROS_DOMAIN_ID (default: 20)'
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
        '--no-goal',
        action='store_true',
        help='Do not auto-publish a terminal goal (lets you send goals manually)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print the generated YAML without launching'
    )

    parser.add_argument(
        '--use-vlm',
        action='store_true',
        help='Use the VLM frontier selector (vlm_goal_selector). Only valid with exploration-singleagent-ground.'
    )

    parser.add_argument(
        '--follow',
        action='store_true',
        help='Spawn the YOLOv8 person_tracker_node alongside the VLM selector. '
             'Off by default because the tracker\'s CPU load (~70%% of one '
             'core for YOLO inference) contends with mpc_node and causes yaw '
             'wobble. Enable when you actually want to use follow-mode.'
    )

    args = parser.parse_args()

    # Find setup.bash path and rviz config
    setup_bash = find_setup_bash(args.setup_bash)
    rviz_config = find_rviz_config()
    print(f"[INFO] Using setup.bash: {setup_bash}")
    print(f"[INFO] Using rviz config: {rviz_config}")

    # Determine sim_env and generate YAML
    if args.mode == 'dyn-test':
        yaml_content = generate_dyn_test_yaml(setup_bash, args.ros_domain_id)
        print(f"[INFO] Mode: Dynamic obstacle test (1 drone + 1 dyn obstacle in Gazebo)")
        print(f"[INFO] Drone at (0, 0, 3.0) facing +x, no goal — observe heat map in RViz")
    elif args.mode == 'dyn-test-ground':
        yaml_content = generate_dyn_test_ground_yaml(setup_bash, args.ros_domain_id)
        print(f"[INFO] Mode: Ground robot + dynamic obstacle test in Gazebo")
        print(f"[INFO] Ground robot at (0, 0, 0) facing +x, obstacle at ~(3, 0, 0.3)")
    elif args.mode == 'dyn-test-ground-mpc':
        yaml_content = generate_dyn_test_ground_mpc_yaml(setup_bash, args.ros_domain_id)
        print(f"[INFO] Mode: Ground robot + MPC + dynamic obstacle test in Gazebo")
        print(f"[INFO] Ground robot at (0, 0, 0), MPC controller, obstacle at ~(3, 0, 0.3)")
        print(f"[INFO] Use '2D Goal Pose' in RViz to send goals")
    elif args.mode == 'interactive':
        yaml_content = generate_interactive_yaml(setup_bash, args.ros_domain_id, rviz_config=rviz_config)
        print(f"[INFO] Mode: Interactive single-agent simulation (sim_env=fake_sim)")
        print(f"[INFO] Agent NX01 at (0, 0, 1.0) — use '2D Goal Pose' in RViz to send goals")
    elif args.mode == 'exploration-singleagent-ground':
        # Single ground robot, frontier-based exploration in Gazebo. Spawn at
        # y=2 (off origin) so the same multi-agent origin guard in
        # exploreSelectCallback doesn't apply — and num_agents=1 disables the
        # peer-coupled paths anyway.
        agents = [{
            'namespace': 'NX01',
            'x': 0.0,
            'y': 2.0,
            'z': 0.0,
            'yaw': 0.0,
        }]
        sim_env = 'gazebo'
        env = args.env if args.env != 'hard_forest' else 'ACL_office'
        yaml_content = generate_exploration_multiagent_ground_yaml(
            setup_bash, agents, args.ros_domain_id, rviz_config=rviz_config,
            sim_env=sim_env, env=env, use_vlm=args.use_vlm,
            use_follow=args.follow)
        print(f"[INFO] Mode: Single-agent ground robot exploration (Gazebo)")
        print(f"[INFO] Environment: {env}")
        print(f"[INFO]   {agents[0]['namespace']}: ({agents[0]['x']}, {agents[0]['y']}, {agents[0]['z']}) yaw={agents[0]['yaw']}")
        if args.use_vlm:
            print(f"[INFO] VLM frontier selector ENABLED — mighty_node yields term_goal to vlm_goal_selector")
            print(f"[INFO] Chat CLI launches in its own pane. Type instructions for the robot there.")
        else:
            print(f"[INFO] Exploration is self-driven — no goal needed")
    elif args.mode == 'exploration-multiagent-ground':
        num = args.num_agents if args.num_agents != 10 else 3
        # Arrange agents in a line at y=2, x-axis spaced 5m apart. The y offset
        # keeps every agent off (0,0,0) so the multi-agent origin guard in
        # exploreSelectCallback (which defers exploration when within 0.5m of
        # origin while peers are active) doesn't deadlock the middle agent.
        spacing = 5.0
        agents = []
        for i in range(num):
            x = -spacing * (num - 1) / 2.0 + spacing * i
            agents.append({
                'namespace': f'NX{i+1:02d}',
                'x': round(x, 3),
                'y': 2.0,
                'z': 0.0,
                'yaw': 0.0,
            })
        # Default to Gazebo + ACL_office; --env overrides the world
        sim_env = 'gazebo'
        env = args.env if args.env != 'hard_forest' else 'ACL_office'
        yaml_content = generate_exploration_multiagent_ground_yaml(
            setup_bash, agents, args.ros_domain_id, rviz_config=rviz_config,
            sim_env=sim_env, env=env)
        print(f"[INFO] Mode: Multi-agent ground robot exploration (Gazebo + MinPos) with {num} agents")
        print(f"[INFO] Environment: {env}")
        for a in agents:
            print(f"[INFO]   {a['namespace']}: ({a['x']}, {a['y']}, {a['z']}) yaw={a['yaw']}")
        print(f"[INFO] Exploration is self-driven — no goal needed")
    elif args.mode == 'swap-multiagent-ground':
        num = args.num_agents if args.num_agents != 10 else 4
        radius = math.sqrt(32)  # corners of 8x8 square → radius = sqrt(4²+4²)
        angle_offset = math.pi / 4  # 45° so agents land on (4,4), (-4,4), (-4,-4), (4,-4)
        agents = generate_multiagent_positions(num, radius, z=0.0, angle_offset=angle_offset)
        env = args.env if args.env != 'hard_forest' else 'ACL_office'
        yaml_content = generate_swap_multiagent_ground_yaml(
            setup_bash, agents, radius, angle_offset, args.ros_domain_id, env=env)
        print(f"[INFO] Mode: Multi-agent ground robot position swap (Gazebo) with {num} agents")
        print(f"[INFO] Environment: {env}")
        for a in agents:
            print(f"[INFO]   {a['namespace']}: ({a['x']}, {a['y']}, {a['z']}) yaw={a['yaw']}")
        print(f"[INFO] Agents swap to diametrically opposite positions")
    elif args.mode == 'multiagent-ground':
        num = args.num_agents if args.num_agents != 10 else 4
        radius = args.radius if args.radius != 10.0 else 12.0
        agents = generate_multiagent_positions(num, radius, z=0.0, prefix='NX')
        yaml_content = generate_multiagent_ground_yaml(setup_bash, agents, radius,
                                                       args.ros_domain_id)
        print(f"[INFO] Mode: Multi-agent ground robot (Gazebo + MPC) with {num} agents (radius={radius})")
        for a in agents:
            print(f"[INFO]   {a['namespace']}: ({a['x']}, {a['y']}, {a['z']}) yaw={a['yaw']}")
    elif args.mode == 'multiagent':
        sim_env = 'fake_sim'
        # Multi-agent uses a dedicated RViz config with per-agent NX01..NX10 display groups
        rviz_config = find_rviz_config('multi_mighty.rviz')
        # Multi-agent uses a dedicated planner config (config/multi_mighty.yaml)
        config_file = Path(__file__).resolve().parent.parent / 'config' / 'multi_mighty.yaml'
        if not config_file.exists():
            print(f"[WARN] {config_file} not found; falling back to default planner config (mighty.yaml)")
            config_file = None
        agents = generate_multiagent_positions(args.num_agents, args.radius)
        yaml_content = generate_multiagent_yaml(setup_bash, agents, sim_env, args.ros_domain_id, args.radius, no_goal=args.no_goal, rviz_config=rviz_config, config_file=config_file)
        print(f"[INFO] Mode: Multi-agent simulation with {args.num_agents} agents (sim_env={sim_env})")
        print(f"[INFO] Using multi-agent rviz config: {rviz_config}")
        if config_file:
            print(f"[INFO] Using multi-agent planner config: {config_file}")
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

        # If frontier-based exploration is enabled in the ground-robot config,
        # auto-skip the goal_sender pane — the exploration loop will issue
        # goals and a manual goal_sender publication would just preempt it.
        no_goal = args.no_goal
        if use_ground_robot and not no_goal:
            cfg_path = (Path(__file__).resolve().parent.parent
                        / 'config' / 'mighty_ground_robot.yaml')
            if exploration_enabled_in_yaml(cfg_path):
                no_goal = True
                print(f"[INFO] Exploration enabled in {cfg_path.name} — "
                      f"skipping goal_sender (frontier loop drives the robot)")

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
            use_ground_robot=use_ground_robot,
            no_goal=no_goal
        )
        print(f"[INFO] Mode: Single-agent Gazebo simulation (sim_env={sim_env})")
        print(f"[INFO] Environment: {args.env} (world: {world_name})")
        if use_ground_robot:
            print(f"[INFO] Vehicle: Ground robot (Pioneer 3-AT)")
        print(f"[INFO] Start: ({start_pos[0]}, {start_pos[1]}, {start_pos[2]})")
        if not no_goal:
            print(f"[INFO] Goal: ({args.goal[0]}, {args.goal[1]}, {args.goal[2]})")

    if args.dry_run:
        print("\n[DRY RUN] Generated YAML:")
        print("-" * 60)
        print(yaml_content)
        print("-" * 60)
        return

    # Kill ONLY the mighty_sim tmux session — other tmux sessions on the same
    # server (the user's own work, etc.) are left alone. If mighty_sim is the
    # only session, tmux will exit on its own, but unrelated sessions persist.
    subprocess.run(['tmux', 'kill-session', '-t', 'mighty_sim'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Kill ONLY gzserver/gzclient processes that came from a prior run of THIS
    # workspace's launch (i.e. command line references the mighty_ws install
    # path). Bare `killall gzserver gzclient` would kill any unrelated Gazebo
    # the user might have running for another project.
    workspace_marker = str(setup_bash.resolve().parent.parent)  # .../mighty_ws
    subprocess.run(
        ['pkill', '-f', f'gz(server|client).*{re.escape(workspace_marker)}'],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    # Write temporary YAML file and launch with tmuxp
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        temp_yaml_path = f.name

    try:
        print(f"[INFO] Launching simulation...")
        env = os.environ.copy()
        env['SETUP_BASH'] = str(setup_bash)
        # Foreground attach — tmuxp creates the session and attaches the
        # current terminal so the user sees panes immediately and doesn't
        # need a separate `tmux attach`.
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
