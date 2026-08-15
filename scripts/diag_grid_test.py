#!/usr/bin/env python3
# /* ----------------------------------------------------------------------------
#  * Diagonal-grid test harness for MIGHTY's HGP anti-corner-cutting logic
#  * (src/hgp/graph_search.cpp, getSucc()).
#  * -------------------------------------------------------------------------- */
"""
Publishes, in the `map` frame, everything the ground-robot 2D planner needs to
compute a path WITHOUT the full sim (no Gazebo, no fake_sim, no real mapper), so
the diagonal corner-cutting behavior can be seen in isolation in RViz:

  occ_2d_topic  nav_msgs/OccupancyGrid    A diagonal STAIRCASE wall of occupied
                                          cells. In ground-robot mode MIGHTY
                                          prefers this source (buildMap2DFromOcc2D
                                          -- binary, NO inflation), so the wall
                                          stays exactly one cell wide and the
                                          diagonal "pinholes" between consecutive
                                          occupied cells are preserved. Those
                                          pinholes are what the new check must
                                          refuse to cross.
  state         dynus_interfaces/State    Fixed start pose on one side of the wall
                                          (published only if publish_state:=true;
                                          set false if fake_sim already owns it).
  term_goal     geometry_msgs/PoseStamped Goal on the far side, so the straight
                                          path wants to cut across the wall.
  TF map->base_frame                      Static transform at the start pose.

Topics are RELATIVE, so launch this node INTO the planner's namespace (e.g.
NX01) and they resolve to /NX01/occ_2d_topic, /NX01/state, /NX01/term_goal.

What to watch in RViz (display the OccupancyGrid + the `hgp_path_marker`
MarkerArray, fixed frame = map):
  * WITHOUT the fix (branch `main`): the path zig-zags diagonally THROUGH the
    staircase corners -- i.e. straight across the wall.
  * WITH the fix (branch `mad_demo`): each diagonal whose two orthogonally-
    adjacent side cells are occupied is rejected, so the path detours around an
    END of the wall.

Occupancy encoding matches the mapper's occ_2d_topic: 0 = free, 100 = occupied,
-1 = unknown (MIGHTY treats unknown as traversable).
"""
import math

import rclpy
from rclpy.node import Node
from rclpy.qos import QoSProfile, QoSDurabilityPolicy, QoSReliabilityPolicy

from nav_msgs.msg import OccupancyGrid
from geometry_msgs.msg import PoseStamped, TransformStamped
from tf2_ros import StaticTransformBroadcaster

try:
    from dynus_interfaces.msg import State
    HAVE_STATE = True
except Exception:  # pragma: no cover - only when dynus_interfaces isn't sourced
    HAVE_STATE = False


class DiagGridTest(Node):
    def __init__(self):
        super().__init__("diag_grid_test")

        # --- parameters (override from the launch or CLI) ---
        self.declare_parameter("frame_id", "map")
        self.declare_parameter("base_frame_id", "base_link")  # relative -> <ns>/base_link
        self.declare_parameter("resolution", 0.15)            # match mighty_map_res
        self.declare_parameter("width_cells", 100)
        self.declare_parameter("height_cells", 100)
        self.declare_parameter("origin_x", -7.5)
        self.declare_parameter("origin_y", -7.5)
        # Diagonal wall drawn as a 1-cell-wide line between these two world points.
        self.declare_parameter("wall_x0", -3.0)
        self.declare_parameter("wall_y0", -3.0)
        self.declare_parameter("wall_x1", 3.0)
        self.declare_parameter("wall_y1", 3.0)
        # Start / goal placed on opposite sides of the diagonal so the straight
        # path crosses it (default: below-right start, above-left goal).
        self.declare_parameter("start_x", 2.5)
        self.declare_parameter("start_y", -2.5)
        self.declare_parameter("goal_x", -2.5)
        self.declare_parameter("goal_y", 2.5)
        self.declare_parameter("goal_z", 0.2)
        self.declare_parameter("publish_state", True)

        gp = self.get_parameter
        self.frame_id = gp("frame_id").value
        self.base_frame_id = gp("base_frame_id").value
        self.res = float(gp("resolution").value)
        self.w = int(gp("width_cells").value)
        self.h = int(gp("height_cells").value)
        self.ox = float(gp("origin_x").value)
        self.oy = float(gp("origin_y").value)
        self.start = (float(gp("start_x").value), float(gp("start_y").value))
        self.goal = (float(gp("goal_x").value), float(gp("goal_y").value), float(gp("goal_z").value))
        self.publish_state = bool(gp("publish_state").value)

        self.grid_msg = self._build_grid(
            float(gp("wall_x0").value), float(gp("wall_y0").value),
            float(gp("wall_x1").value), float(gp("wall_y1").value),
        )

        # Latched QoS (RELIABLE + TRANSIENT_LOCAL) is compatible with any
        # subscriber and delivers the map to late joiners; we also republish.
        latched = QoSProfile(depth=1)
        latched.durability = QoSDurabilityPolicy.TRANSIENT_LOCAL
        latched.reliability = QoSReliabilityPolicy.RELIABLE

        self.pub_occ = self.create_publisher(OccupancyGrid, "occ_2d_topic", latched)
        self.pub_goal = self.create_publisher(PoseStamped, "term_goal", latched)
        if self.publish_state:
            if not HAVE_STATE:
                raise RuntimeError(
                    "publish_state:=true but dynus_interfaces/State is not importable -- "
                    "source your mighty_ws install first."
                )
            self.pub_state = self.create_publisher(State, "state", 10)

        self.tf_static = StaticTransformBroadcaster(self)
        self._publish_static_tf()

        # Publish once immediately, then keep republishing (occ/goal latched, but
        # re-sending is cheap and covers any subscriber that starts later).
        self._publish_all()
        self.create_timer(0.5, self._publish_all)          # occ + goal @ 2 Hz
        if self.publish_state:
            self.create_timer(0.02, self._publish_state)    # state @ 50 Hz

        occ = sum(1 for v in self.grid_msg.data if v == 100)
        self.get_logger().info(
            f"diag_grid_test: {self.w}x{self.h} @ {self.res}m, {occ} occupied (diagonal wall), "
            f"start={self.start} goal={self.goal[:2]} frame={self.frame_id} "
            f"publish_state={self.publish_state}"
        )

    def _world_to_cell(self, wx, wy):
        cx = int(math.floor((wx - self.ox) / self.res))
        cy = int(math.floor((wy - self.oy) / self.res))
        return cx, cy

    def _build_grid(self, x0, y0, x1, y1):
        msg = OccupancyGrid()
        msg.header.frame_id = self.frame_id
        msg.info.resolution = self.res
        msg.info.width = self.w
        msg.info.height = self.h
        msg.info.origin.position.x = self.ox
        msg.info.origin.position.y = self.oy
        msg.info.origin.position.z = 0.0
        msg.info.origin.orientation.w = 1.0
        data = [0] * (self.w * self.h)  # 0 = free

        # Rasterize the diagonal as a 1-cell-wide line of occupied cells. Sampling
        # finely and marking the containing cell yields the classic staircase, so
        # consecutive occupied cells touch only at their corners -> the diagonal
        # pinholes the corner-cut check is meant to close.
        c0 = self._world_to_cell(x0, y0)
        c1 = self._world_to_cell(x1, y1)
        steps = max(abs(c1[0] - c0[0]), abs(c1[1] - c0[1])) * 4 + 1
        for i in range(steps + 1):
            t = i / steps
            wx = x0 + t * (x1 - x0)
            wy = y0 + t * (y1 - y0)
            cx, cy = self._world_to_cell(wx, wy)
            if 0 <= cx < self.w and 0 <= cy < self.h:
                data[cy * self.w + cx] = 100  # occupied
        msg.data = data
        return msg

    def _stamp(self):
        return self.get_clock().now().to_msg()

    def _publish_static_tf(self):
        tf = TransformStamped()
        tf.header.stamp = self._stamp()
        tf.header.frame_id = self.frame_id
        tf.child_frame_id = self.base_frame_id
        tf.transform.translation.x = self.start[0]
        tf.transform.translation.y = self.start[1]
        tf.transform.translation.z = 0.0
        tf.transform.rotation.w = 1.0
        self.tf_static.sendTransform(tf)

    def _publish_all(self):
        self.grid_msg.header.stamp = self._stamp()
        self.pub_occ.publish(self.grid_msg)

        goal = PoseStamped()
        goal.header.stamp = self._stamp()
        goal.header.frame_id = self.frame_id
        goal.pose.position.x = self.goal[0]
        goal.pose.position.y = self.goal[1]
        goal.pose.position.z = self.goal[2]
        goal.pose.orientation.w = 1.0
        self.pub_goal.publish(goal)

    def _publish_state(self):
        st = State()
        st.header.stamp = self._stamp()
        st.header.frame_id = self.frame_id
        st.pos.x = self.start[0]
        st.pos.y = self.start[1]
        st.pos.z = 0.0
        st.quat.w = 1.0
        self.pub_state.publish(st)


def main():
    rclpy.init()
    node = DiagGridTest()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.try_shutdown()


if __name__ == "__main__":
    main()
