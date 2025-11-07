#!/usr/bin/env python3
import csv
import glob
import os
import sys
import math

# ROS 2 bag API
from rosbag2_py import SequentialReader, StorageOptions, ConverterOptions
import rclpy.serialization
from rosidl_runtime_py.utilities import get_message

# Helper function to deserialize a message.
def deserialize_message(serialized_msg, msg_type):
    return rclpy.serialization.deserialize_message(serialized_msg, msg_type)

def read_ros2_bag(bag_path, topic_name, message_type_str):
    """
    Reads a ROS 2 bag (folder) and extracts messages from a specified topic.
    
    Args:
        bag_path (str): Path to the ROS 2 bag folder.
        topic_name (str): The name of the topic to extract data from.
        message_type_str (str): The message type as a string (e.g., "my_package/msg/State").
    
    Returns:
        list: A list of messages from the specified topic.
    """
    storage_options = StorageOptions(uri=bag_path, storage_id="sqlite3")
    converter_options = ConverterOptions(
        input_serialization_format="cdr", output_serialization_format="cdr"
    )
    reader = SequentialReader()
    reader.open(storage_options, converter_options)
    
    all_topics = reader.get_all_topics_and_types()
    topic_names = [topic.name for topic in all_topics]
    
    if topic_name not in topic_names:
        print("Topic '{}' not found in bag '{}'.".format(topic_name, bag_path))
        return []
    
    # Dynamically resolve the message type.
    msg_type = get_message(message_type_str)
    
    data = []
    while reader.has_next():
        topic, serialized_msg, _ = reader.read_next()
        if topic == topic_name:
            msg = deserialize_message(serialized_msg, msg_type)
            data.append(msg)
    return data

def load_forest_parameters(csv_path):
    """Load forest parameters from CSV and return a list of cylinder dictionaries."""
    cylinders = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cyl = {
                "id": row["id"],
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
                "radius": float(row["radius"]),
                "height": float(row["height"])
            }
            cylinders.append(cyl)
    return cylinders

def check_collision(agent_pos, cylinder, drone_radius=0.1):
    """
    Check if agent_pos (tuple of x,y,z) collides with the given cylinder.
    Cylinder is defined by (x, y, z, radius, height), where z is the cylinder's base
    (here spawned with z = height/2, so it spans z=0 to z=height).
    If a collision occurs, return a tuple (True, cylinder_id, penetration_distance),
    where penetration_distance = (drone_radius + cylinder["radius"]) - horizontal distance.
    Otherwise, return (False, None, None).
    """
    dx = agent_pos[0] - cylinder["x"]
    dy = agent_pos[1] - cylinder["y"]
    horizontal_dist = math.sqrt(dx*dx + dy*dy)
    allowed_dist = cylinder["radius"] + drone_radius
    if horizontal_dist < allowed_dist:
        if 0 <= agent_pos[2] <= cylinder["height"]:
            penetration = allowed_dist - horizontal_dist
            print("\033[91mCollision with cylinder {} at position ({:.2f}, {:.2f}, {:.2f}), penetration {:.8f} m\033[0m".format(
                cylinder["id"], agent_pos[0], agent_pos[1], agent_pos[2], penetration))
            return (True, cylinder["id"], penetration)
    return (False, None, None)

def process_ros2_bag(bag_path, cylinders, drone_radius, topic, message_type_str):
    """
    Process a ROS 2 bag folder by reading messages from the specified topic.
    For each message, check for collision with any cylinder.
    Returns a tuple (collision_found, collided_obstacle, penetration).
    """
    messages = read_ros2_bag(bag_path, topic, message_type_str)
    for msg in messages:
        # Assume msg has a field 'pos' (geometry_msgs/Point)
        agent_pos = (msg.p.x, msg.p.y, msg.p.z)
        for cyl in cylinders:
            collision, coll_id, pen = check_collision(agent_pos, cyl, drone_radius)
            if collision:
                return (True, coll_id, pen)
    return (False, None, None)

def main():
    if len(sys.argv) < 6:
        print("Usage: {} <forest_csv_path> <bag_folder> <drone_radius> <pos_topic_name> <message_type_str>".format(sys.argv[0]))
        print("Example: python3 collision_checker.py /home/kkondo/data/easy_forest_obstacle_parameters.csv /media/kkondo/T7/dynus/tro_paper/global_planner_benchmarking/bags/dgp 0.1 /NX01/goal dynus_interfaces/msg/Goal")
        sys.exit(1)
    
    forest_csv = sys.argv[1]
    bag_folder = sys.argv[2]
    drone_radius = float(sys.argv[3])
    topic = sys.argv[4]
    message_type_str = sys.argv[5]
    
    cylinders = load_forest_parameters(forest_csv)
    
    # add num_0 to num9 to at the end of the bag paths
    # bag_paths = [os.path.join(bag_path, f"num_{i}/num_{i}") for i in range(10) for bag_path in glob.glob(bag_folder)]
    bag_paths = [os.path.join(bag_path, f"num_{i}") for i in range(10) for bag_path in glob.glob(bag_folder)]
    bag_paths = sorted(bag_paths)
    
    results = []
    for bag_path in bag_paths:
        print("Processing bag: {}".format(bag_path))
        collision, collided_obstacle, penetration = process_ros2_bag(bag_path, cylinders, drone_radius, topic, message_type_str)
        if collision:
            status = "collision"
            # print red
            print("\033[91mCollision detected in bag: {}\033[0m".format(bag_path))
        else:
            status = "no collision"
            collided_obstacle = ""
            penetration = ""
            # print green
            print("\033[92mNo collision detected in bag: {}\033[0m".format(bag_path))
        results.append((os.path.basename(bag_path), status, collided_obstacle, penetration))
    
    # Write results to collision_check.csv.
    output_csv = f"{bag_folder}/collision_check.csv"
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bag_folder", "collision_status", "collided_obstacle", "penetration"])
        for bag_name, status, coll, pen in results:
            writer.writerow([bag_name, status, coll, pen])
    
    print("Collision check complete. Results written to {}".format(output_csv))

if __name__ == "__main__":
    main()
