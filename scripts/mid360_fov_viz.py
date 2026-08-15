#!/usr/bin/env python3
# /* ----------------------------------------------------------------------------
#  * RViz visualization of the tilted Livox Mid-360 ground-coverage envelope used
#  * by the perception-aware planner (hgp/perception_planner.hpp Mid360FOV).
#  *
#  * Draws, at the robot's current pose (from dynus_interfaces/State, else a static
#  * pose), a MarkerArray on `mid360_fov`:
#  *   - FAR envelope   : green ring = elliptical trust boundary (~2.30 m).
#  *   - NEAR/BLIND ring : red ring  = ground-start horn (~1.1 m); inside it is the
#  *                       sensor blind spot.
#  *   - covered sector : translucent green fill between near and far (what the
#  *                       robot can actually observe this instant).
#  *   - blind sector   : translucent red fill from the robot out to the near ring.
#  *
#  * Reuses the exact prototype Mid360FOV precompute (verified to match the C++:
#  * r_max=2.30 m, fov=29.75 deg). Standalone (no dependency on the C++ planner),
#  * so it works both in diag_test.launch.py and in the live ground-robot sim.
#  * -------------------------------------------------------------------------- */
import math

import numpy as np
import rclpy
from rclpy.node import Node
from rclpy.qos import QoSProfile, QoSDurabilityPolicy, QoSReliabilityPolicy

from geometry_msgs.msg import Point
from std_msgs.msg import ColorRGBA
from visualization_msgs.msg import Marker, MarkerArray

try:
    from dynus_interfaces.msg import State
    HAVE_STATE = True
except Exception:  # pragma: no cover
    HAVE_STATE = False


def compute_mid360_envelope(tilt_deg=20.0, e_lo=-7.0, e_hi=52.0, ground_start_fwd=1.10,
                            F=2.30, W=0.315, xc=0.40, bin_deg=0.25):
    """Return (bearings_rad, r_near, r_far) over the azimuth bins where both the
    ground horn and the ellipse are valid (the FOV wedge). Faithful to the
    prototype / C++ Mid360FOV."""
    tau = math.radians(tilt_deg)
    h = ground_start_fwd * math.tan(math.radians(-e_lo) + tau)
    phi = np.radians(np.arange(-180.0, 180.0, 0.1))
    ele = np.radians(np.arange(e_lo, e_hi + 1e-9, 0.2))
    PHI, ELE = np.meshgrid(phi, ele)
    dxw = np.cos(ELE) * np.cos(PHI) * math.cos(tau) + np.sin(ELE) * math.sin(tau)
    dyw = np.cos(ELE) * np.sin(PHI)
    dzw = -np.cos(ELE) * np.cos(PHI) * math.sin(tau) + np.sin(ELE) * math.cos(tau)
    hit = dzw < -1e-9
    t = np.where(hit, h / np.maximum(-dzw, 1e-12), np.nan)
    gx, gy = (t * dxw)[hit], (t * dyw)[hit]
    brg = np.degrees(np.arctan2(gy, gx))
    r = np.hypot(gx, gy)
    edges = np.arange(-180.0, 180.0 + bin_deg, bin_deg)
    idx = np.clip(np.digitize(brg, edges) - 1, 0, len(edges) - 2)
    rn = np.full(len(edges) - 1, np.inf)
    np.minimum.at(rn, idx, r)
    bc = np.radians(0.5 * (edges[:-1] + edges[1:]))
    a, b = F - xc, W
    A, B, C, D = 1 / a**2, 1 / b**2, -2 * xc / a**2, xc**2 / a**2 - 1
    qa = A * np.cos(bc) ** 2 + B * np.sin(bc) ** 2
    qb = C * np.cos(bc)
    rf = (-qb + np.sqrt(qb * qb - 4 * qa * D)) / (2 * qa)
    okb = np.isfinite(rn) & (rf > rn)
    order = np.argsort(bc[okb])
    return bc[okb][order], rn[okb][order], rf[okb][order]


class Mid360FovViz(Node):
    def __init__(self):
        super().__init__("mid360_fov_viz")
        self.declare_parameter("frame_id", "map")
        self.declare_parameter("z", 0.05)                # lay the fans just above ground
        self.declare_parameter("use_state", True)        # follow the robot via /state
        self.declare_parameter("static_x", 0.0)
        self.declare_parameter("static_y", 0.0)
        self.declare_parameter("static_theta", 0.0)
        # Mid-360 constants (defaults match the C++/prototype).
        self.declare_parameter("tilt_deg", 20.0)
        self.declare_parameter("e_lo", -7.0)
        self.declare_parameter("e_hi", 52.0)
        self.declare_parameter("ground_start_fwd", 1.10)
        self.declare_parameter("F", 2.30)
        self.declare_parameter("W", 0.315)
        self.declare_parameter("xc", 0.40)

        gp = self.get_parameter
        self.frame_id = gp("frame_id").value
        self.z = float(gp("z").value)
        self.use_state = bool(gp("use_state").value) and HAVE_STATE
        self.pose = (float(gp("static_x").value), float(gp("static_y").value),
                     float(gp("static_theta").value))

        self.bc, self.rn, self.rf = compute_mid360_envelope(
            float(gp("tilt_deg").value), float(gp("e_lo").value), float(gp("e_hi").value),
            float(gp("ground_start_fwd").value), float(gp("F").value), float(gp("W").value),
            float(gp("xc").value))
        self.get_logger().info(
            f"Mid360 FOV: {len(self.bc)} bins, r_near~{self.rn.min():.2f}m, "
            f"r_far(max)={self.rf.max():.2f}m, fov={2*math.degrees(abs(self.bc).max()):.1f}deg")

        latched = QoSProfile(depth=1)
        latched.durability = QoSDurabilityPolicy.TRANSIENT_LOCAL
        latched.reliability = QoSReliabilityPolicy.RELIABLE
        self.pub = self.create_publisher(MarkerArray, "mid360_fov", latched)

        if self.use_state:
            self.create_subscription(State, "state", self._on_state, 10)
        self.create_timer(0.1, self._publish)  # 10 Hz

    def _on_state(self, msg):
        yaw = 2.0 * math.atan2(msg.quat.z, msg.quat.w)  # planar yaw from quaternion
        self.pose = (msg.pos.x, msg.pos.y, yaw)

    def _pt(self, x, y, px, py, cth, sth):
        p = Point()
        p.x = px + x * cth - y * sth
        p.y = py + x * sth + y * cth
        p.z = self.z
        return p

    def _publish(self):
        px, py, th = self.pose
        cth, sth = math.cos(th), math.sin(th)
        near = [(self.rn[i] * math.cos(self.bc[i]), self.rn[i] * math.sin(self.bc[i]))
                for i in range(len(self.bc))]
        far = [(self.rf[i] * math.cos(self.bc[i]), self.rf[i] * math.sin(self.bc[i]))
               for i in range(len(self.bc))]

        arr = MarkerArray()
        stamp = self.get_clock().now().to_msg()

        def base(mid, mtype, r, g, b, a, scale):
            m = Marker()
            m.header.frame_id = self.frame_id
            m.header.stamp = stamp
            m.ns = "mid360_fov"
            m.id = mid
            m.type = mtype
            m.action = Marker.ADD
            m.pose.orientation.w = 1.0
            m.scale.x = m.scale.y = m.scale.z = scale
            m.color = ColorRGBA(r=r, g=g, b=b, a=a)
            return m

        # covered sector (green fill) between near and far arcs
        cov = base(0, Marker.TRIANGLE_LIST, 0.1, 0.8, 0.2, 0.25, 1.0)
        for i in range(len(self.bc) - 1):
            n0, n1, f0, f1 = near[i], near[i + 1], far[i], far[i + 1]
            for (ax, ay), (bx, by), (cx, cy) in ((n0, f0, f1), (n0, f1, n1)):
                cov.points += [self._pt(ax, ay, px, py, cth, sth),
                               self._pt(bx, by, px, py, cth, sth),
                               self._pt(cx, cy, px, py, cth, sth)]
        arr.markers.append(cov)

        # blind sector (red fill) from robot origin out to the near arc
        blind = base(1, Marker.TRIANGLE_LIST, 0.9, 0.1, 0.1, 0.30, 1.0)
        for i in range(len(self.bc) - 1):
            n0, n1 = near[i], near[i + 1]
            for (ax, ay), (bx, by), (cx, cy) in (((0.0, 0.0), n0, n1),):
                blind.points += [self._pt(ax, ay, px, py, cth, sth),
                                 self._pt(bx, by, px, py, cth, sth),
                                 self._pt(cx, cy, px, py, cth, sth)]
        arr.markers.append(blind)

        # far envelope ring (solid green line)
        far_line = base(2, Marker.LINE_STRIP, 0.1, 0.9, 0.2, 0.9, 0.04)
        far_line.points = [self._pt(x, y, px, py, cth, sth) for (x, y) in far]
        arr.markers.append(far_line)

        # near / blind boundary ring (solid red line)
        near_line = base(3, Marker.LINE_STRIP, 0.95, 0.15, 0.15, 0.95, 0.04)
        near_line.points = [self._pt(x, y, px, py, cth, sth) for (x, y) in near]
        arr.markers.append(near_line)

        self.pub.publish(arr)


def main():
    rclpy.init()
    node = Mid360FovViz()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.try_shutdown()


if __name__ == "__main__":
    main()
