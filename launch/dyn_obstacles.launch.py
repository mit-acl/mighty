#!/usr/bin/env python3
"""
Single Source of Motion launch:
- Spawns static obstacle models (no internal motion plugin; no traj_x/y/z expressions).
- After spawning (or a fixed delay), starts dynamic_forest_node which:
    * Publishes /trajs (DynTraj for planners)
    * Updates Gazebo model poses (physical motion)
"""

import os
from launch import LaunchDescription
from launch.actions import (
    DeclareLaunchArgument,
    OpaqueFunction,
    TimerAction
)
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
from launch_ros.parameter_descriptions import ParameterValue
from launch.substitutions import Command
import json

# Global variable
_OBSTACLES_JSON_STORAGE = {}

# ---------- Helpers ----------

def _as_bool(context, name, default=False):
    raw = LaunchConfiguration(name).perform(context)
    if raw is None:
        return default
    s = raw.strip().lower()
    if s == '':
        return default
    return s in ['true', '1', 'yes', 'on', 't']

def _as(context, name, cast, default):
    """
    Safe retrieval & casting of a LaunchConfiguration.
    If the argument is absent or empty, returns default; otherwise casts.
    """
    raw = LaunchConfiguration(name).perform(context)
    if raw is None:
        print(f"[dyn_obstacles] {name}=None -> default {default}")
        return default
    s = raw.strip()
    if s == '':
        print(f"[dyn_obstacles] {name}=<empty> -> default {default}")
        return default
    try:
        return cast(s)
    except ValueError:
        raise RuntimeError(f"[dyn_obstacles] Argument '{name}' expected {cast.__name__}, got '{s}'")

# Optional: flip to True for one‑time verbose dump
DEBUG_DYN_OBS = False


# ---------- Spawn Static Obstacles ----------

def trefoil_expr_with_vel(x0, y0, z0, sx, sy, sz, offset, slower):
    """
    Return position and velocity expression strings for trefoil knot.

    Args:
        x0,y0,z0 : centers
        sx,sy,sz : scale factors
        offset   : phase offset
        slower   : time scaling (bigger => slower)

    Returns:
        (x_str, y_str, z_str, vx_str, vy_str, vz_str)
    """
    tt = f"t/{slower}+{offset}"

    # Position
    x_str = f"{sx/6.0}*(sin({tt})+2*sin(2*{tt}))+{x0}"
    y_str = f"{sy/5.0}*(cos({tt})-2*cos(2*{tt}))+{y0}"
    z_str = f"{(sz/2.0)}*(-sin(3*{tt}))+{z0}"

    inv_slow = f"(1/{slower})"  # explicit to keep expression exact

    # Velocity (derivatives)
    vx_str = f"{sx/6.0}*{inv_slow}*(cos({tt})+4*cos(2*{tt}))"
    vy_str = f"{sy/5.0}*{inv_slow}*(-sin({tt})+4*sin(2*{tt}))"
    vz_str = f"-{3*sz/2.0}*{inv_slow}*cos(3*{tt})"

    return x_str, y_str, z_str, vx_str, vy_str, vz_str

def _spawn_static_block(context):
    import random

    num_obstacles   = _as(context, 'num_obstacles', int,   30)
    x_min           = _as(context, 'x_min', float, 5.0)
    x_max           = _as(context, 'x_max', float, 105.0)
    y_min           = _as(context, 'y_min', float, -10.0)
    y_max           = _as(context, 'y_max', float, 10.0)
    z_min           = _as(context, 'z_min', float, 0.0)
    z_max           = _as(context, 'z_max', float, 6.0)
    slower_min      = _as(context, 'slower_min', float, 10.0)
    slower_max      = _as(context, 'slower_max', float, 12.0)
    spawn_interval  = _as(context, 'spawn_interval', float, 10.0)
    seed            = _as(context, 'seed', int, 0)
    use_sim_time    = _as_bool(context, 'use_sim_time', False)
    urdf_xacro      = LaunchConfiguration('urdf_xacro').perform(context) or 'dyn_obstacle1.urdf.xacro'

    scale_range = [[2.0, 5.0], [5.0, 10.0], [2.0, 4.0]]

    offset_range = [0.0, 3.0] # [offset_min, offset_max]
    slower_range = [slower_min, slower_max]

    size = 1.0  # Default size for the obstacle

    # set seed
    random.seed(seed)

    if DEBUG_DYN_OBS:
        print("[dyn_obstacles][spawn] num_obstacles:", num_obstacles,
              "seed:", seed, "spawn_interval:", spawn_interval)

    urdf_path = os.path.join(get_package_share_directory('mighty'), 'urdf', urdf_xacro)

    actions = []
    obstacles_meta = []
    for i in range(num_obstacles):
        entity = f"obstacle_{i}"

        x = x_min + (x_max - x_min) * random.random()
        y = y_min + (y_max - y_min) * random.random()
        z = z_min + (z_max - z_min) * random.random()

        sx = scale_range[0][0] + (scale_range[0][1] - scale_range[0][0]) * random.random()
        sy = scale_range[1][0] + (scale_range[1][1] - scale_range[1][0]) * random.random()
        sz = scale_range[2][0] + (scale_range[2][1] - scale_range[2][0]) * random.random()
        offset = random.uniform(*offset_range)
        slower = random.uniform(*slower_range)

        x_str, y_str, z_str, vx_str, vy_str, vz_str = trefoil_expr_with_vel(
            x0=x, y0=y, z0=z, sx=sx, sy=sy, sz=sz,
            offset=offset, slower=slower
        )

        # for Gazebo
        robot_description = ParameterValue(
            Command([
                'xacro ', urdf_path,
                ' traj_x:=', x_str,
                ' traj_y:=', y_str,
                ' traj_z:=', z_str,
                ' size:=', str(size),
                ' namespace:=', entity
            ]),
            value_type=str
        )

        rsp = Node(
            package='robot_state_publisher',
            executable='robot_state_publisher',
            name=f'{entity}_rsp',
            output='screen',
            parameters=[{
                'robot_description': robot_description,
                'use_sim_time': use_sim_time,
                'frame_prefix': entity + '/'
            }],
            remappings=[('/robot_description', f'/{entity}/robot_description')]
        )

        spawn = Node(
            package='gazebo_ros',
            executable='spawn_entity.py',
            name=f'{entity}_spawn',
            output='screen',
            arguments=[
                '-entity', entity,
                '-topic', f'/{entity}/robot_description',
                '-x', str(x), '-y', str(y), '-z', str(z)
            ]
        )

        actions.append(
            TimerAction(
                period=i * spawn_interval,
                actions=[rsp, 
                         TimerAction(period=0.3, actions=[spawn])]
            )
        )

        # Pass parameters to dynamic_forest_node (where we publish DynTraj)
        obstacles_meta.append({
            "name": entity,
            "x0": x, "y0": y, "z0": z,
            "scale_x": sx, "scale_y": sy, "scale_z": sz,
            "offset": offset, "slower": slower,
            "traj_x": x_str, "traj_y": y_str, "traj_z": z_str,
            "traj_vx": vx_str, "traj_vy": vy_str, "traj_vz": vz_str,
            "size": size,
        })



    # Save the metadata to a JSON file
    # get the home directory of the user
    print("Saving obstacles metadata to context.locals['obstacles_json']")
    _OBSTACLES_JSON_STORAGE['obstacles_json'] = json.dumps(obstacles_meta)

    return actions

# ---------- Launch dynamic_forest_node After Spawns ----------

def _maybe_launch_forest_node(context):
    if not _as_bool(context, 'launch_forest_node', True):
        return []
    
    print("maybe_launch_forest_node: launching dynamic_forest_node")

    delay            = _as(context, 'forest_start_delay', float, 1.0)
    total_num_obs    = _as(context, 'num_obstacles', int, 10)
    dynamic_ratio    = _as(context, 'dynamic_ratio', float, 0.5)
    publish_rate_hz  = _as(context, 'publish_rate_hz', float, 50.0)
    seed             = _as(context, 'seed', int, 0)
    publish_markers  = _as_bool(context, 'publish_markers', True)
    publish_tf       = _as_bool(context, 'publish_tf', True)

    if DEBUG_DYN_OBS:
        print("[dyn_obstacles][forest_node] delay:", delay,
              "total_num_obs:", total_num_obs,
              "ratio:", dynamic_ratio,
              "rate:", publish_rate_hz)

    obstacles_json = _OBSTACLES_JSON_STORAGE.get('obstacles_json', '[]')

    params = {
        'obstacles_json': obstacles_json,
        'use_external_obstacles_json': True,
        'total_num_obs': total_num_obs,
        'dynamic_ratio': dynamic_ratio,
        'publish_rate_hz': publish_rate_hz,
        'seed': seed,
        'publish_markers': publish_markers,
        'publish_tf': publish_tf,
        'use_spawn_origins': False
    }

    forest_node = Node(
        package='mighty',
        executable='dynamic_forest_node',
        name='dynamic_forest_trajs',
        output='screen',
        parameters=[params],
        # prefix='xterm -e gdb -q -ex run --args', # gdb debugging
    )

    return [TimerAction(period=delay, actions=[forest_node])]


# ---------- Main Launch Description ----------

def generate_launch_description():
    args = [
        DeclareLaunchArgument('num_obstacles', default_value='100'),
        DeclareLaunchArgument('seed', default_value='0'),
        DeclareLaunchArgument('dynamic_ratio', default_value='0.5'),
        DeclareLaunchArgument('publish_rate_hz', default_value='50.0'),
        DeclareLaunchArgument('spawn_interval', default_value='0.1'),
        DeclareLaunchArgument('forest_start_delay', default_value='3.0'),
        DeclareLaunchArgument('launch_forest_node', default_value='true'),
        DeclareLaunchArgument('publish_markers', default_value='true'),
        DeclareLaunchArgument('publish_tf', default_value='true'),
        DeclareLaunchArgument('use_sim_time', default_value='false'),

        # Spatial ranges
        DeclareLaunchArgument('x_min', default_value='5.0'),
        DeclareLaunchArgument('x_max', default_value='105.0'),
        DeclareLaunchArgument('y_min', default_value='-5.0'),
        DeclareLaunchArgument('y_max', default_value='5.0'),
        DeclareLaunchArgument('z_min', default_value='1.0'),
        DeclareLaunchArgument('z_max', default_value='5.0'),

        # Trajectory params
        DeclareLaunchArgument('slower_min', default_value='4.0'),
        DeclareLaunchArgument('slower_max', default_value='6.0'),
        DeclareLaunchArgument('scale_global', default_value='1.0'),

        # URDF
        DeclareLaunchArgument('urdf_xacro', default_value='dyn_obstacle1.urdf.xacro'),
    ]

    ld = LaunchDescription(args)
    ld.add_action(OpaqueFunction(function=_spawn_static_block))
    ld.add_action(OpaqueFunction(function=_maybe_launch_forest_node))
    return ld

if __name__ == '__main__':
    generate_launch_description()
