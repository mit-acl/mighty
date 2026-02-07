import os
from launch import LaunchDescription
from launch.actions import IncludeLaunchDescription, TimerAction
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory


def generate_launch_description():

    mighty_dir = get_package_share_directory('mighty')

    # ── 1. Simulator (map generator + RViz) ─────────────────────────────
    simulator = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(mighty_dir, 'launch', 'simulator.launch.py')
        ),
    )

    # ── 2. Agent RR01 ───────────────────────────────────────────────────
    rr01 = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(mighty_dir, 'launch', 'onboard_mighty.launch.py')
        ),
        launch_arguments={
            'namespace': 'RR01',
            'x': '0.0',
            'y': '0.0',
            'z': '1.0',
            'yaw': '0',
            'sim_env': 'fake_sim',
            'map_frame_id': 'RR01/map',
            'num_agents': '2',
            'use_frame_alignment': 'true',
        }.items(),
    )

    # ── 3. Agent RR02 (delayed 5 s so RR01 is up first) ─────────────────
    rr02 = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(mighty_dir, 'launch', 'onboard_mighty.launch.py')
        ),
        launch_arguments={
            'namespace': 'RR02',
            'x': '0.0',
            'y': '0.0',
            'z': '1.0',
            'yaw': '0',
            'sim_env': 'fake_sim',
            'map_frame_id': 'RR02/map',
            'num_agents': '2',
            'use_frame_alignment': 'true',
        }.items(),
    )

    # ── 4. Frame-alignment publishers (ground-truth transforms) ──────────
    #  T^{RR01/map}_{RR02/map}: R=Rz(90°), t=[5,0,0]
    fa_rr01_rr02 = Node(
        package='mighty',
        executable='frame_align_publisher.py',
        name='fa_rr01_rr02',
        output='screen',
        parameters=[{
            'ego_name': 'RR01',
            'other_name': 'RR02',
            'tx': 5.0,
            'ty': 0.0,
            'tz': 0.0,
            'qx': 0.0,
            'qy': 0.0,
            'qz': 0.7071,
            'qw': 0.7071,
        }],
    )

    #  T^{RR02/map}_{RR01/map}: R=Rz(-90°), t=[0,5,0]
    fa_rr02_rr01 = Node(
        package='mighty',
        executable='frame_align_publisher.py',
        name='fa_rr02_rr01',
        output='screen',
        parameters=[{
            'ego_name': 'RR02',
            'other_name': 'RR01',
            'tx': 0.0,
            'ty': 5.0,
            'tz': 0.0,
            'qx': 0.0,
            'qy': 0.0,
            'qz': -0.7071,
            'qw': 0.7071,
        }],
    )

    # ── Static TFs so RViz can use "map" as fixed frame ────────────────
    #  map → RR01/map  (identity — RR01's frame IS the world frame)
    static_tf_map_rr01 = Node(
        package='tf2_ros', executable='static_transform_publisher',
        name='static_tf_map_rr01', output='screen',
        arguments=['0','0','0','0','0','0','1', 'map', 'RR01/map'],
    )
    #  map → RR02/map  (RR02's origin is at [5,0,0] rotated 90° yaw in world)
    static_tf_map_rr02 = Node(
        package='tf2_ros', executable='static_transform_publisher',
        name='static_tf_map_rr02', output='screen',
        arguments=['5','0','0','0','0','0.7071','0.7071', 'map', 'RR02/map'],
    )

    # ── 5. Goal sender (use_hardware=true so goals use {agent}/map) ──────
    goal_sender = Node(
        package='mighty',
        executable='goal_sender.py',
        name='goal_sender',
        output='screen',
        parameters=[{
            'list_agents': ['RR01', 'RR02'],
            'list_goals': ['[5.0, 5.0]', '[-5.0, -5.0]'],
            'default_goal_z': 1.0,
            'use_hardware': True,
        }],
    )

    return LaunchDescription([
        # Simulator + static TFs first (so RViz has "map" frame immediately)
        simulator,
        static_tf_map_rr01,
        static_tf_map_rr02,
        # RR01 after 3 s
        TimerAction(period=3.0, actions=[rr01]),
        # RR02 after 8 s
        TimerAction(period=8.0, actions=[rr02]),
        # Frame-alignment publishers after 10 s
        TimerAction(period=10.0, actions=[fa_rr01_rr02, fa_rr02_rr01]),
        # Goals after 15 s
        TimerAction(period=15.0, actions=[goal_sender]),
    ])
