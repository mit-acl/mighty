#!/usr/bin/env python3
"""
Run Gazebo first:
  ros2 launch gazebo_ros gazebo.launch.py

Then:
  ros2 run mighty generate_random_forest.py --ros-args -p difficulty:=hard -p min_clearance:=2.0 -p shape_mode:=mixed -p box_probability:=0.5

Add this to .world for ROS2 bridge:
  <plugin name='gazebo_ros_state' filename='libgazebo_ros_state.so'>
    <ros>
      <namespace>/plug</namespace>
      <argument>model_states:=model_states_plug</argument>
      <argument>link_states:=link_states_plug</argument>
    </ros>
    <update_rate>100.0</update_rate>
  </plugin>

Usage:

# all boxes
ros2 run mighty generate_random_forest.py --ros-args -p shape_mode:=box

# mixed (70% boxes), random yaw, 2.0 m clearance
ros2 run mighty generate_random_forest.py --ros-args -p shape_mode:=mixed -p box_probability:=0.7 -p min_clearance:=2.0
  
Parameters (new ones starred ★):
- seed (int)
- difficulty (easy, medium, hard, dynamic)   # "dynamic" leaves x∈[35,65] free of static obstacles
- min_clearance (float, m)                   # gap between obstacle *footprints* (conservative via bounding circle)
- max_place_tries (int)
- shrink_after_ratio (float in [0,1])
- shrink_rate (float in (0,1))
- shrink_min_radius (float)                  # applies to cylinders; <=0 → use kMinRadius
★ shape_mode: "cylinder" | "box" | "mixed"
★ box_probability (float in [0,1])          # probability of box when shape_mode="mixed"
★ box_width_min/max, box_depth_min/max (m)  # footprint size ranges for boxes
★ yaw_random (bool)                          # random yaw for boxes
"""
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Pose
from gazebo_msgs.srv import SpawnEntity, DeleteEntity
from visualization_msgs.msg import Marker, MarkerArray
import random
import math
import time
import csv
import os

# ---------- helpers ----------
def yaw_to_quat(yaw: float):
    s = math.sin(0.5 * yaw)
    c = math.cos(0.5 * yaw)
    return 0.0, 0.0, s, c

class GazeboPrimitiveSpawner:
    def __init__(self, node: Node):
        self.node = node
        self.client = node.create_client(SpawnEntity, 'spawn_entity')
        while not self.client.wait_for_service(timeout_sec=1.0):
            self.node.get_logger().info('Waiting for spawn_entity service...')

    def spawn_cylinder(self, name: str, frame_id: str, x: float, y: float, z: float,
                       radius: float, height: float, mass: float = 2.0):
        self._spawn_primitive(name, False, frame_id, x, y, z,
                              0.0, 0.0, 0.0, 1.0,
                              radius, height, depth=1.0, mass=mass)

    def spawn_box(self, name: str, frame_id: str, x: float, y: float, z: float,
                  width: float, depth: float, height: float,
                  yaw: float = 0.0, mass: float = 2.0):
        qx, qy, qz, qw = yaw_to_quat(yaw)
        self._spawn_primitive(name, True, frame_id, x, y, z,
                              qx, qy, qz, qw,
                              width, height, depth, mass)

    def _spawn_primitive(self, name: str, do_cube: bool, frame_id: str,
                         x: float, y: float, z: float,
                         qx: float, qy: float, qz: float, qw: float,
                         width_or_radius: float, height: float, depth: float, mass: float):
        pose = Pose()
        pose.position.x = x; pose.position.y = y; pose.position.z = z
        pose.orientation.x = qx; pose.orientation.y = qy; pose.orientation.z = qz; pose.orientation.w = qw

        if do_cube:
            geometry_str = f"<box><size>{width_or_radius} {depth} {height}</size></box>"
            # Note: Gazebo box size order is X Y Z; we pass (width, depth, height)
            mass12 = mass / 12.0
            xx = mass12 * (depth**2 + height**2)
            yy = mass12 * (width_or_radius**2 + height**2)
            zz = mass12 * (width_or_radius**2 + depth**2)
        else:
            geometry_str = f"<cylinder><length>{height}</length><radius>{width_or_radius}</radius></cylinder>"
            mass12 = mass / 12.0
            xx = mass12 * (3 * width_or_radius**2 + height**2)
            yy = mass12 * (3 * width_or_radius**2 + height**2)
            zz = 0.5 * mass * width_or_radius**2

        sdf = (
            f"<?xml version='1.0'?>"
            f"<sdf version='1.4'>"
            f"<model name='{name}'>"
            f"<static>1</static>"
            f"<link name='link'>"
            f"<inertial><mass>{mass}</mass><inertia>"
            f"<ixx>{xx}</ixx><ixy>0.0</ixy><ixz>0.0</ixz>"
            f"<iyy>{yy}</iyy><iyz>0.0</iyz><izz>{zz}</izz>"
            f"</inertia></inertial>"
            f"<collision name='collision'>"
            f"<geometry>{geometry_str}</geometry>"
            f"<surface><contact><collide_without_contact>true</collide_without_contact></contact></surface>"
            f"</collision>"
            f"<visual name='visual'>"
            f"<geometry>{geometry_str}</geometry>"
            f"<material><script>"
            f"<uri>file://media/materials/scripts/gazebo.material</uri>"
            f"<name>Gazebo/Blue</name>"
            f"</script></material>"
            f"</visual>"
            f"<gravity>0</gravity>"
            f"</link>"
            f"</model>"
            f"</sdf>"
        )

        request = SpawnEntity.Request()
        request.name = name
        request.xml = sdf
        request.robot_namespace = "obstacle_spawner"
        request.initial_pose = pose
        request.reference_frame = frame_id

        future = self.client.call_async(request)
        rclpy.spin_until_future_complete(self.node, future)

class ForestSpawner(Node):
    def __init__(self):
        super().__init__('create_random_forest')

        # ---------------- Params ----------------
        self.declare_parameter("seed", 0)
        self.declare_parameter("difficulty", "dynamic")
        self.declare_parameter("min_clearance", 2.0)
        self.declare_parameter("max_place_tries", 2000)
        self.declare_parameter("shrink_after_ratio", 0.5)
        self.declare_parameter("shrink_rate", 0.5)
        self.declare_parameter("shrink_min_radius", -1.0)  # <= 0 → use kMinRadius

        # shape/mix params
        self.declare_parameter("shape_mode", "cylinder")     # "cylinder" | "box" | "mixed"
        self.declare_parameter("box_probability", 0.5)       # used when shape_mode="mixed"
        self.declare_parameter("box_width_min", 1.0)
        self.declare_parameter("box_width_max", 3.0)
        self.declare_parameter("box_depth_min", 10.0)
        self.declare_parameter("box_depth_max", 16.0)
        self.declare_parameter("yaw_random", True)

        seed = int(self.get_parameter("seed").value)
        difficulty = str(self.get_parameter("difficulty").value)
        self.min_clearance = max(0.0, float(self.get_parameter("min_clearance").value))
        self.max_place_tries = int(self.get_parameter("max_place_tries").value)
        self.shrink_after_ratio = min(max(float(self.get_parameter("shrink_after_ratio").value), 0.0), 1.0)
        self.shrink_rate = min(max(float(self.get_parameter("shrink_rate").value), 0.05), 0.999)
        self.param_shrink_min_radius = float(self.get_parameter("shrink_min_radius").value)

        self.shape_mode = str(self.get_parameter("shape_mode").value).lower()
        self.box_probability = float(self.get_parameter("box_probability").value)
        self.box_probability = min(max(self.box_probability, 0.0), 1.0)
        self.box_w_min = float(self.get_parameter("box_width_min").value)
        self.box_w_max = float(self.get_parameter("box_width_max").value)
        self.box_d_min = float(self.get_parameter("box_depth_min").value)
        self.box_d_max = float(self.get_parameter("box_depth_max").value)
        self.yaw_random = bool(self.get_parameter("yaw_random").value)

        # --- Difficulty & density ---
        self.dynamic = False
        if difficulty == "easy":
            self.density = 0.05
            CSV_PATH = "/home/kkondo/data/easy_forest_obstacle_parameters.csv"
        elif difficulty == "medium":
            self.density = 0.1
            CSV_PATH = "/home/kkondo/data/medium_forest_obstacle_parameters.csv"
        elif difficulty == "hard":
            self.density = 0.4
            CSV_PATH = "/home/kkondo/data/hard_forest_obstacle_parameters.csv"
        elif difficulty == "dynamic":
            self.density = 0.05
            self.dynamic = True
            self.dynamic_min_x = 35.0
            self.dynamic_max_x = 65.0
            CSV_PATH = "/home/kkondo/data/dynamic_forest_obstacle_parameters.csv"
        else:
            self.get_logger().warn(f"Unknown difficulty '{difficulty}', defaulting to medium.")
            self.density = 0.1
            CSV_PATH = "/home/kkondo/data/medium_forest_obstacle_parameters.csv"

        # Map bounds
        self.min_x = 3.0
        self.max_x = 303.0
        self.min_y = -20.0
        self.max_y =  20.0
        self.size_z = 5.0

        usable_area = (self.max_x - self.min_x) * (self.max_y - self.min_y)
        target_area = self.density * usable_area

        self.get_logger().info(
            f"Seed: {seed}, Difficulty: {difficulty}, min_clearance: {self.min_clearance} m, "
            f"max_place_tries: {self.max_place_tries}, shrink_after_ratio: {self.shrink_after_ratio}, "
            f"shrink_rate: {self.shrink_rate}, shape_mode: {self.shape_mode}, box_prob: {self.box_probability}"
            + (f", dynamic-free zone: [{self.dynamic_min_x},{self.dynamic_max_x}]" if self.dynamic else "")
        )
        self.get_logger().info(
            f"Usable area: {usable_area:.2f} m^2 | Target density: {self.density:.3f} → Target area: {target_area:.2f} m^2"
        )

        # Publisher & RNG
        self.marker_pub = self.create_publisher(MarkerArray, "forest", 10)
        random.seed(seed)

        # Spawner + ground
        self.spawner = GazeboPrimitiveSpawner(self)
        self.spawn_ground_plane()

        # Footprint/height ranges
        self.kMinHeight = 1.0 if difficulty in ("easy", "dynamic") else 6.0
        self.kMaxHeight = 5.0 if difficulty in ("easy", "dynamic") else 6.0
        self.kMinRadius = 1.0
        self.kMaxRadius = 1.5

        # Shrink thresholds
        self.shrink_min_radius = self.param_shrink_min_radius if self.param_shrink_min_radius > 0.0 else self.kMinRadius
        self.shrink_threshold = max(1, int(self.max_place_tries * self.shrink_after_ratio))

        occupied_area = 0.0
        marker_array = MarkerArray()
        n_objs = 0
        # store (x, y, eff_radius) for clearance checks
        placed = []

        def valid_with_clearance(x, y, eff_r):
            for (xj, yj, rj) in placed:
                if math.hypot(x - xj, y - yj) < (eff_r + rj + self.min_clearance):
                    return False
            return True

        # CSV
        os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
        with open(CSV_PATH, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["id", "shape", "x", "y", "z", "radius", "width", "depth", "height", "yaw"])

            while occupied_area < target_area:
                height = random.uniform(self.kMinHeight, self.kMaxHeight)

                # --- choose shape ---
                shape = "cylinder"
                if self.shape_mode == "box":
                    shape = "box"
                elif self.shape_mode == "mixed":
                    shape = "box" if (random.random() < self.box_probability) else "cylinder"

                # --- sample initial footprint ---
                if shape == "cylinder":
                    radius = random.uniform(self.kMinRadius, self.kMaxRadius)
                    width = depth = None
                    yaw = 0.0
                else:
                    width = random.uniform(self.box_w_min, self.box_w_max)
                    depth = random.uniform(self.box_d_min, self.box_d_max)
                    radius = None
                    yaw = (random.uniform(-math.pi, math.pi) if self.yaw_random else 0.0)

                # --- place with retries + shrinking ---
                placed_ok = False
                shrunk_once = False

                for attempt in range(1, self.max_place_tries + 1):
                    # Shrink if struggling
                    if attempt >= self.shrink_threshold:
                        if shape == "cylinder" and radius is not None and radius > self.shrink_min_radius + 1e-9:
                            new_r = max(radius * self.shrink_rate, self.shrink_min_radius)
                            if new_r < radius - 1e-9:
                                if not shrunk_once:
                                    self.get_logger().info(
                                        f"Shrinking cylinder r {radius:.3f}→{new_r:.3f} (attempt {attempt}/{self.max_place_tries})"
                                    )
                                    shrunk_once = True
                                radius = new_r
                        elif shape == "box" and width is not None and depth is not None:
                            min_w = max(0.05, self.box_w_min)
                            min_d = max(0.05, self.box_d_min)
                            new_w = max(width * self.shrink_rate, min_w)
                            new_d = max(depth * self.shrink_rate, min_d)
                            if (new_w < width - 1e-9) or (new_d < depth - 1e-9):
                                if not shrunk_once:
                                    self.get_logger().info(
                                        f"Shrinking box (w,d) ({width:.3f},{depth:.3f})→({new_w:.3f},{new_d:.3f}) "
                                        f"(attempt {attempt}/{self.max_place_tries})"
                                    )
                                    shrunk_once = True
                                width, depth = new_w, new_d

                    # effective bounding radius for boundary/clearance
                    # eff_r = radius if shape == "cylinder" else 0.5 * math.hypot(width, depth)
                    eff_r = radius if shape == "cylinder" else width

                    # if footprint too large to fit, keep shrinking/trying
                    low_x, high_x = self.min_x + eff_r, self.max_x - eff_r
                    low_y, high_y = self.min_y + eff_r, self.max_y - eff_r
                    if (low_x >= high_x) or (low_y >= high_y):
                        continue

                    pos_x = random.uniform(low_x, high_x)
                    pos_y = random.uniform(low_y, high_y)

                    # Respect dynamic free band in x
                    if self.dynamic:
                        if not (pos_x + eff_r <= self.dynamic_min_x or pos_x - eff_r >= self.dynamic_max_x):
                            continue

                    if valid_with_clearance(pos_x, pos_y, eff_r):
                        placed_ok = True
                        break

                if not placed_ok:
                    self.get_logger().warn(
                        "Could not place more obstacles without violating clearance; stopping early. "
                        "Consider reducing density, increasing tries, or allowing smaller sizes."
                    )
                    break

                # ---- spawn into Gazebo ----
                pos_z = height / 2.0
                if shape == "cylinder":
                    self.spawner.spawn_cylinder(str(n_objs), "world", pos_x, pos_y, pos_z, radius, height, mass=2.0)
                    area = math.pi * radius * radius
                    placed.append((pos_x, pos_y, radius))
                else:
                    self.spawner.spawn_box(str(n_objs), "world", pos_x, pos_y, pos_z, width, depth, height, yaw=yaw, mass=2.0)
                    area = width * depth
                    eff_r = 0.5 * math.hypot(width, depth)
                    placed.append((pos_x, pos_y, eff_r))

                occupied_area += area

                # ---- CSV + RViz marker ----
                if shape == "cylinder":
                    csv_writer.writerow([n_objs, "cylinder", pos_x, pos_y, pos_z, radius, 0.0, 0.0, height, 0.0])
                    marker = Marker()
                    marker.header.frame_id = "world"
                    marker.header.stamp = self.get_clock().now().to_msg()
                    marker.ns = "obstacles"; marker.id = n_objs
                    marker.type = Marker.CYLINDER; marker.action = Marker.ADD
                    marker.pose.position.x = pos_x; marker.pose.position.y = pos_y; marker.pose.position.z = pos_z
                    marker.pose.orientation.w = 1.0
                    marker.scale.x = 2 * radius; marker.scale.y = 2 * radius; marker.scale.z = height
                else:
                    qx, qy, qz, qw = yaw_to_quat(yaw)
                    csv_writer.writerow([n_objs, "box", pos_x, pos_y, pos_z, 0.0, width, depth, height, yaw])
                    marker = Marker()
                    marker.header.frame_id = "world"
                    marker.header.stamp = self.get_clock().now().to_msg()
                    marker.ns = "obstacles"; marker.id = n_objs
                    marker.type = Marker.CUBE; marker.action = Marker.ADD
                    marker.pose.position.x = pos_x; marker.pose.position.y = pos_y; marker.pose.position.z = pos_z
                    marker.pose.orientation.x = qx; marker.pose.orientation.y = qy
                    marker.pose.orientation.z = qz; marker.pose.orientation.w = qw
                    marker.scale.x = width; marker.scale.y = depth; marker.scale.z = height

                marker.color.r = 0.0; marker.color.g = 1.0; marker.color.b = 0.0; marker.color.a = 1.0
                marker_array.markers.append(marker)

                n_objs += 1

                # ---- progress log ----
                current_density = occupied_area / usable_area
                completion_pct = min(occupied_area / target_area * 100.0, 100.0)
                self.get_logger().info(
                    f"[{n_objs:04d}] {shape:8s}  density={current_density:.4f}  "
                    f"({completion_pct:.1f}% of target)  occupied={occupied_area:.2f}/{target_area:.2f} m^2"
                )

        self.get_logger().info(f"Final density: {(occupied_area/usable_area):.4f} "
                               f"({min(occupied_area/target_area*100.0,100.0):.1f}% of target)")
        self.get_logger().info(f"Total occupied area = {occupied_area:.2f} m^2, Number of objects: {n_objs}")
        time.sleep(3)
        self.get_logger().info("Publishing marker array for visualization.")
        self.marker_pub.publish(marker_array)

    def delete_model(self, model_name: str):
        delete_client = self.create_client(DeleteEntity, 'delete_entity')
        while not delete_client.wait_for_service(timeout_sec=1.0):
            self.get_logger().info('Waiting for delete_entity service...')
        req = DeleteEntity.Request(); req.name = model_name
        future = delete_client.call_async(req)
        rclpy.spin_until_future_complete(self, future)

    def spawn_ground_plane(self):
        margin = 5.0
        size_x = (self.max_x - self.min_x) + 2 * margin
        size_y = (self.max_y - self.min_y) + 2 * margin
        center_x = (self.min_x + self.max_x) / 2.0
        center_y = (self.min_y + self.max_y) / 2.0

        self.delete_model("ground_plane")

        sdf = f"""<?xml version="1.0" ?>
        <sdf version="1.5">
          <model name='ground_plane_map'>
            <static>1</static>
            <link name='link'>
              <collision name='collision'>
                <geometry><plane><normal>0 0 1</normal><size>{size_x} {size_y}</size></plane></geometry>
                <surface>
                  <contact><collide_without_contact>true</collide_without_contact></contact>
                  <friction><ode><mu>100</mu><mu2>50</mu2></ode></friction>
                </surface>
              </collision>
              <visual name='visual'>
                <cast_shadows>0</cast_shadows>
                <geometry><plane><normal>0 0 1</normal><size>{size_x} {size_y}</size></plane></geometry>
                <material><script>
                  <uri>file://media/materials/scripts/gazebo.material</uri>
                  <name>Gazebo/Grey</name>
                </script></material>
              </visual>
            </link>
          </model>
        </sdf>
        """
        pose = Pose()
        pose.position.x = center_x; pose.position.y = center_y; pose.position.z = 0.0
        pose.orientation.w = 1.0

        req = SpawnEntity.Request()
        req.name = "ground_plane_map"; req.xml = sdf
        req.robot_namespace = "ground_plane_map"
        req.initial_pose = pose; req.reference_frame = "world"

        future = self.spawner.client.call_async(req)
        rclpy.spin_until_future_complete(self, future)

def main(args=None):
    rclpy.init(args=args)
    node = ForestSpawner()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()