# /* ----------------------------------------------------------------------------
#  * Standalone diagonal-grid test for the HGP anti-corner-cutting logic.
#  *
#  * Brings up ONLY the MIGHTY planner (ground-robot / 2D mode), RViz, and the
#  * diag_grid_test publisher -- no Gazebo, no fake_sim, no mapper -- so the
#  * planned path reacts purely to the synthetic diagonal wall. The robot does not
#  * drive (nothing integrates the trajectory); we only watch the planned
#  * `hgp_path_marker` update as it routes through / around the wall.
#  *
#  * Usage:
#  *   ros2 launch mighty diag_test.launch.py
#  *
#  * Then in RViz set Fixed Frame = map, add the OccupancyGrid (/NX01/occ_2d_topic)
#  * and the MarkerArray (/NX01/hgp_path_marker). Rebuild first so the C++ fix is
#  * compiled in (colcon build --packages-select mighty). Compare branch `main`
#  * (path cuts through the wall) vs `mad_demo` (path goes around).
#  * -------------------------------------------------------------------------- */
import os
import yaml

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def launch_setup(context, *args, **kwargs):
    namespace = LaunchConfiguration("namespace").perform(context)
    use_rviz = LaunchConfiguration("use_rviz").perform(context).lower() in ("true", "1")

    share = get_package_share_directory("mighty")

    # Load the ground-robot planner config and force the 2D / ground-robot path.
    cfg_path = os.path.join(share, "config", "mighty_ground_robot.yaml")
    with open(cfg_path) as f:
        params = yaml.safe_load(f)["mighty_node"]["ros__parameters"]
    params["vehicle_type"] = "ground_robot"   # zDim_ == 1 -> diagonal check active
    params["map_frame_id"] = "map"            # single fixed frame for the test
    params["sim_env"] = "fake_sim"            # avoid gazebo/hardware code paths

    mighty_node = Node(
        package="mighty",
        executable="mighty",
        name="mighty_node",
        namespace=namespace,
        output="screen",
        emulate_tty=True,
        parameters=[params],
        # We don't feed a point cloud; point the sensor input at an unused topic
        # so the ground-robot 2D map comes solely from our occ_2d_topic.
        remappings=[("lidar_cloud_in", "unused_lidar_cloud")],
    )

    diag_pub = Node(
        package="mighty",
        executable="diag_grid_test.py",
        name="diag_grid_test",
        namespace=namespace,
        output="screen",
        parameters=[{
            "frame_id": "map",
            "base_frame_id": "base_link",  # -> <ns>/base_link
            "publish_state": True,          # no fake_sim here, so we own /state
        }],
    )

    rviz_cfg = os.path.join(share, "rviz", "mighty_sim_ground_robot.rviz")
    rviz_node = Node(
        package="rviz2",
        executable="rviz2",
        name="rviz2",
        arguments=(["-d", rviz_cfg] if os.path.exists(rviz_cfg) else []),
        output="screen",
        condition=None,
    )

    # Mid-360 coverage overlay: draws the elliptical far envelope + near/blind ring
    # at the robot pose (follows /state). Publishes MarkerArray on <ns>/mid360_fov.
    fov_viz = Node(
        package="mighty",
        executable="mid360_fov_viz.py",
        name="mid360_fov_viz",
        namespace=namespace,
        output="screen",
        parameters=[{"frame_id": "map", "use_state": True}],
    )

    nodes = [mighty_node, diag_pub, fov_viz]
    if use_rviz:
        nodes.append(rviz_node)
    return nodes


def generate_launch_description():
    return LaunchDescription([
        DeclareLaunchArgument("namespace", default_value="NX01"),
        DeclareLaunchArgument("use_rviz", default_value="true"),
        OpaqueFunction(function=launch_setup),
    ])
