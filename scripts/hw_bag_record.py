# /* ----------------------------------------------------------------------------
#  * Copyright 2025, Kota Kondo, Aerospace Controls Laboratory
#  * Massachusetts Institute of Technology
#  * All Rights Reserved
#  * Authors: Kota Kondo, et al.
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */

import os
import subprocess
import argparse

def record_ros2_bag(bag_name, bag_path, agents, topics=None):

    # Per-agent topics — this list mirrors what the RR04 rviz config displays
    # (rviz/mighty_RR04.rviz) plus a few essentials for offline analysis
    # (state, cmd_vel, dlio, lookahead).
    base_topics = [
        # --- Global paths / subgoals / trajectories ---
        "/original_hgp_path_marker",
        "/hgp_path_marker",
        "/traj_committed_colored",
        "/traj_subopt_colored",
        "/hover_avoidance_viz",
        "/point_A",
        "/point_E",
        "/point_G",
        "/point_G_term",
        "/term_goal",
        "/goal",
        # --- Raw sensors ---
        "/mid360_PointCloud2",
        "/livox/lidar",
        "/d435/depth/color/points",
        "/d435/color/image_raw",
        "/d435/color/image_raw/compressed",
        "/d435/color/camera_info",
        # --- 3D voxel grids ---
        "/free_grid",
        "/unknown_grid",
        "/occupancy_grid",
        "/dynamic_grid",
        "/heat_cloud",
        # --- 2D projections / ESDF / exploration ---
        "/occ_2d_topic",
        "/occ_2d_topic_updates",
        "/esdf_2d_topic",
        "/esdf_2d_topic_updates",
        "/ground_2d_occupied",
        "/ground_2d_heat",
        "/exploration/visited_map",
        "/exploration/visited_map_updates",
        "/exploration/frontiers",
        "/exploration/current_goal",
        # --- Vehicle state & viz ---
        "/state",
        "/drone_marker",
        "/drone_marker_array",
        "/vel_text",
        "/fov",
        "/tracked_obstacles",
        "/cluster_bounding_boxes",
        "/uncertainty_spheres",
        # --- Control (cmd_vel) ---
        "/cmd_vel",
        "/cmd_vel_auto",
        # --- DLIO odometry ---
        "/dlio/odom_node/pose",
        "/dlio/odom_node/odom",
    ]

    # Static / non-agent-namespaced topics
    static_topics = [
        "/clicked_point",
        "/clock",
        "/initialpose",
        "/parameter_events",
        "/rosout",
        "/tf",
        "/tf_static",
        "/map_generator/global_cloud",
        "/robot_description",
    ]

    # Generate topics for all agents
    all_topics = []
    for agent in agents:
        for topic in base_topics:
            all_topics.append(f"/{agent}{topic}")

    # Add static topics (non-agent-specific topics)
    all_topics.extend(static_topics)

    # Use provided topics if specified, otherwise default to generated topics
    if topics is None:
        topics = all_topics
    
    # Build the ros2 bag record command
    command = ["ros2", "bag", "record", "-o", os.path.join(bag_path, bag_name)] + topics
    
    # Execute the command
    try:
        subprocess.run(command, check=True)
        print(f"Recording started for bag: {bag_name} at path: {bag_path}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to start recording: {e}")

# Main function
if __name__ == "__main__":

    # Create arguments for record_ros2_bag
    parser = argparse.ArgumentParser(description="Record a ROS2 bag")
    parser.add_argument("--bag_number", type=int, help="Bag number to record")
    parser.add_argument("--bag_name", type=str, help="Custom bag name (overrides --bag_number)", default=None)
    parser.add_argument("--bag_path", type=str, help="Path to save the bag",
                        default="/home/swarm/data/mighty_office_exploration")
    parser.add_argument("--agents", nargs="+", help="List of agents to record",
                        default=["RR04"])
    args = parser.parse_args()

    # Customize bag name and path as needed
    if args.bag_name:
        bag_name = args.bag_name
    else:
        # bag_name = "num_" + str(args.bag_number)
        bag_name = "hw_" + str(args.bag_number)
    bag_path = args.bag_path

    # Convert string list (eg. "['NX01', 'NX02']") to list (eg. ['NX01', 'NX02'])
    if args.agents:
        args.agents = args.agents[0].replace("[", "").replace("]", "").replace("'", "").split(", ")

    print("Bag name:", bag_name)
    print("Bag path:", bag_path)
    print("Agents:", args.agents)

    # List of agents
    agents = args.agents

    record_ros2_bag(bag_name, bag_path, agents)