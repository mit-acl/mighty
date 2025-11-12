#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "gcopter/trajectory.hpp"
#include "gcopter/quickhull.hpp"
#include "gcopter/geo_utils.hpp"

#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64.hpp>
#include <geometry_msgs/msg/point.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <visualization_msgs/msg/marker_array.hpp>

#include <chrono>
#include <memory>
#include <cmath>
#include <iostream>

// Visualizer for the planner
class Visualizer
{
private:
    rclcpp::Node &node_; // Node reference

    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr routePub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr wayPointsPub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr trajectoryPub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr meshPub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr edgePub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr spherePub;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr text_pub_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr mighty_waypoints_pub_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr mighty_traj_pub_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr mighty_text_pub_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr mighty_sphere_pub_;

public:
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr speedPub;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr thrPub;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr tiltPub;
    rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr bdrPub;

public:
    Visualizer(rclcpp::Node &node) : node_(node)
    {
        routePub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/route", 10);
        wayPointsPub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/waypoints", 10);
        trajectoryPub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/trajectory", 10);
        meshPub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/mesh", 10);
        edgePub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/edge", 10);
        spherePub = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/spheres", 10);
        speedPub = node_.create_publisher<std_msgs::msg::Float64>("/visualizer/speed", 10);
        thrPub = node_.create_publisher<std_msgs::msg::Float64>("/visualizer/total_thrust", 10);
        tiltPub = node_.create_publisher<std_msgs::msg::Float64>("/visualizer/tilt_angle", 10);
        bdrPub = node_.create_publisher<std_msgs::msg::Float64>("/visualizer/body_rate", 10);
        text_pub_ = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/vel_text", 10);
        mighty_waypoints_pub_ = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/mighty/waypoints", 10);
        mighty_traj_pub_ = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/mighty/trajectory", 10);
        mighty_text_pub_ = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/mighty/vel_text", 10);
        mighty_sphere_pub_ = node_.create_publisher<visualization_msgs::msg::Marker>("/visualizer/mighty/spheres", 10);
    }

    template <int D>
    inline void visualize(const Trajectory<D> &traj,
                          const std::vector<Eigen::Vector3d> &route)
    {
        using namespace visualization_msgs::msg;

        Marker routeMarker, wayPointsMarker, trajMarker;

        routeMarker.id = 0;
        routeMarker.type = Marker::LINE_LIST;
        routeMarker.header.stamp = node_.now();
        routeMarker.header.frame_id = "odom";
        routeMarker.pose.orientation.w = 1.0;
        routeMarker.action = Marker::ADD;
        routeMarker.ns = "route";
        routeMarker.color.r = 1.0;
        routeMarker.color.g = 0.0;
        routeMarker.color.b = 0.0;
        routeMarker.color.a = 1.0;
        routeMarker.scale.x = 0.1;

        wayPointsMarker = routeMarker;
        wayPointsMarker.id = -1;
        wayPointsMarker.type = Marker::SPHERE_LIST;
        wayPointsMarker.ns = "waypoints";
        wayPointsMarker.scale.x = wayPointsMarker.scale.y = wayPointsMarker.scale.z = 0.35;

        trajMarker = routeMarker;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.0;
        trajMarker.color.g = 0.5;
        trajMarker.color.b = 1.0;
        trajMarker.scale.x = 0.3;

        if (!route.empty())
        {
            bool first = true;
            Eigen::Vector3d last;
            for (const auto &pt : route)
            {
                if (first)
                {
                    last = pt;
                    first = false;
                    continue;
                }
                geometry_msgs::msg::Point p;
                p.x = last(0);
                p.y = last(1);
                p.z = last(2);
                routeMarker.points.push_back(p);
                p.x = pt(0);
                p.y = pt(1);
                p.z = pt(2);
                routeMarker.points.push_back(p);
                last = pt;
            }
            routePub->publish(routeMarker);
        }

        if (traj.getPieceNum() > 0)
        {
            Eigen::MatrixXd wps = traj.getPositions();
            for (int i = 0; i < wps.cols(); ++i)
            {
                geometry_msgs::msg::Point p;
                p.x = wps(0, i);
                p.y = wps(1, i);
                p.z = wps(2, i);
                wayPointsMarker.points.push_back(p);
            }
            wayPointsPub->publish(wayPointsMarker);
        }

        if (traj.getPieceNum() > 0)
        {
            double T = 0.01;
            Eigen::Vector3d lastX = traj.getPos(0.0);
            for (double t = T; t < traj.getTotalDuration(); t += T)
            {
                Eigen::Vector3d X = traj.getPos(t);
                geometry_msgs::msg::Point p1, p2;
                p1.x = lastX(0);
                p1.y = lastX(1);
                p1.z = lastX(2);
                p2.x = X(0);
                p2.y = X(1);
                p2.z = X(2);
                trajMarker.points.push_back(p1);
                trajMarker.points.push_back(p2);
                lastX = X;
            }
            trajectoryPub->publish(trajMarker);
        }
    }

    inline void visualizePolytope(const std::vector<Eigen::MatrixX4d> &hPolys)
    {
        Eigen::Matrix3Xd mesh(3, 0), curTris(3, 0), oldTris(3, 0);

        for (const auto &hPoly : hPolys)
        {
            oldTris = mesh;
            Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
            geo_utils::enumerateVs(hPoly, vPoly);

            quickhull::QuickHull<double> qh;
            const auto hull = qh.getConvexHull(vPoly.data(), vPoly.cols(), false, true);
            const auto &idxBuffer = hull.getIndexBuffer();
            int hNum = idxBuffer.size() / 3;

            curTris.resize(3, hNum * 3);
            for (int i = 0; i < hNum * 3; ++i)
                curTris.col(i) = vPoly.col(idxBuffer[i]);

            mesh.resize(3, oldTris.cols() + curTris.cols());
            mesh.leftCols(oldTris.cols()) = oldTris;
            mesh.rightCols(curTris.cols()) = curTris;
        }

        using namespace visualization_msgs::msg;

        Marker meshMarker, edgeMarker;
        meshMarker.id = 0;
        meshMarker.header.stamp = node_.now();
        meshMarker.header.frame_id = "odom";
        meshMarker.pose.orientation.w = 1.0;
        meshMarker.action = Marker::ADD;
        meshMarker.type = Marker::TRIANGLE_LIST;
        meshMarker.ns = "mesh";
        meshMarker.color.r = 0.0;
        meshMarker.color.g = 0.0;
        meshMarker.color.b = 1.0;
        meshMarker.color.a = 0.15;
        meshMarker.scale.x = 1.0;
        meshMarker.scale.y = 1.0;
        meshMarker.scale.z = 1.0;

        edgeMarker = meshMarker;
        edgeMarker.type = Marker::LINE_LIST;
        edgeMarker.ns = "edge";
        edgeMarker.color.g = 1.0;
        edgeMarker.color.a = 1.0;
        edgeMarker.scale.x = 0.02;

        for (int i = 0; i < mesh.cols(); ++i)
        {
            geometry_msgs::msg::Point p;
            p.x = mesh(0, i);
            p.y = mesh(1, i);
            p.z = mesh(2, i);

            meshMarker.points.push_back(p);
        }

        for (int i = 0; i < mesh.cols() / 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                geometry_msgs::msg::Point p1, p2;
                p1.x = mesh(0, 3 * i + j);
                p1.y = mesh(1, 3 * i + j);
                p1.z = mesh(2, 3 * i + j);
                p2.x = mesh(0, 3 * i + (j + 1) % 3);
                p2.y = mesh(1, 3 * i + (j + 1) % 3);
                p2.z = mesh(2, 3 * i + (j + 1) % 3);
                edgeMarker.points.push_back(p1);
                edgeMarker.points.push_back(p2);
            }
        }

        meshPub->publish(meshMarker);
        edgePub->publish(edgeMarker);
    }

    inline void visualizeSphere(const Eigen::Vector3d &center, const double &radius)
    {
        using namespace visualization_msgs::msg;

        Marker sphereMarker;
        sphereMarker.id = 0;
        sphereMarker.type = Marker::SPHERE_LIST;
        sphereMarker.header.stamp = node_.now();
        sphereMarker.header.frame_id = "odom";
        sphereMarker.pose.orientation.w = 1.0;
        sphereMarker.action = Marker::ADD;
        sphereMarker.ns = "spheres";
        sphereMarker.color.b = 1.0;
        sphereMarker.color.a = 1.0;
        sphereMarker.scale.x = sphereMarker.scale.y = sphereMarker.scale.z = radius * 2.0;

        geometry_msgs::msg::Point p;
        p.x = center(0);
        p.y = center(1);
        p.z = center(2);
        sphereMarker.points.push_back(p);

        spherePub->publish(sphereMarker);
    }

    inline void visualizeStartGoal(const Eigen::Vector3d &center, const double &radius, const int sg)
    {
        using namespace visualization_msgs::msg;

        Marker marker;
        marker.id = sg;
        marker.type = Marker::SPHERE_LIST;
        marker.header.stamp = node_.now();
        marker.header.frame_id = "odom";
        marker.pose.orientation.w = 1.0;
        marker.action = Marker::ADD;
        marker.ns = "StartGoal";
        marker.color.r = 1.0;
        marker.color.a = 1.0;
        marker.scale.x = marker.scale.y = marker.scale.z = radius * 2.0;

        geometry_msgs::msg::Point p;
        p.x = center(0);
        p.y = center(1);
        p.z = center(2);
        marker.points.push_back(p);

        spherePub->publish(marker);
    }

    std_msgs::msg::ColorRGBA getColorJet(double v, double vmin, double vmax)
    {
        std_msgs::msg::ColorRGBA c;
        c.r = 1;
        c.g = 1;
        c.b = 1;
        c.a = 1;
        // white
        double dv;

        if (v < vmin)
            v = vmin;
        if (v > vmax)
            v = vmax;
        dv = vmax - vmin;

        if (v < (vmin + 0.25 * dv))
        {
            c.r = 0;
            c.g = 4 * (v - vmin) / dv;
        }
        else if (v < (vmin + 0.5 * dv))
        {
            c.r = 0;
            c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
        }
        else if (v < (vmin + 0.75 * dv))
        {
            c.r = 4 * (v - vmin - 0.5 * dv) / dv;
            c.b = 0;
        }
        else
        {
            c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
            c.b = 0;
        }

        return (c);
    }

    // Sample and publish a degree-5 Bezier trajectory (CP[seg][0..5], durations T[seg])
    // Colors are per-point using getColorJet(speed, 0, v_max_axis_or_auto).
    inline void visualizeBezier(
        const std::vector<std::array<Eigen::Vector3d, 6>> &CP,
        const std::vector<double> &T,
        const std::string &ns = "mighty", // use a different namespace than GCOPTER
        double width = 0.06,
        int samples_per_seg = 120,
        float r = 0.0f, float g = 1.0f, float b = 0.0f, float a = 1.0f,
        const std::string &frame_id = "odom", // set to your RViz fixed frame if different
        double v_max_axis = -1.0              // pass par_.v_max here; if <=0, auto-scale to sampled max
    )
    {
        if (!mighty_traj_pub_ || CP.empty() || T.empty())
            return;

        // Bernstein for degree-5 and degree-4
        auto B5 = [](double u)
        {
            const double w = 1.0 - u;
            return std::array<double, 6>{
                w * w * w * w * w,
                5 * u * w * w * w * w,
                10 * u * u * w * w * w,
                10 * u * u * u * w * w,
                5 * u * u * u * u * w,
                u * u * u * u * u};
        };
        auto B4 = [](double u)
        {
            const double w = 1.0 - u;
            return std::array<double, 5>{
                w * w * w * w,
                4 * u * w * w * w,
                6 * u * u * w * w,
                4 * u * u * u * w,
                u * u * u * u};
        };

        // First pass (optional): find max speed if v_max_axis not provided.
        double vmax_used = v_max_axis;
        if (vmax_used <= 0.0)
        {
            double vmax_obs = 1e-9;
            for (size_t s = 0; s < CP.size(); ++s)
            {
                const int N = std::max(4, samples_per_seg);
                const double Ts = (s < T.size()) ? std::max(T[s], 1e-9) : 1.0;
                for (int k = 0; k <= N; ++k)
                {
                    if (s > 0 && k == 0)
                        continue; // avoid duplicate seam point
                    const double u = static_cast<double>(k) / N;

                    // dp/du for degree-5 Bezier: 5 * sum_{i=0..4} (P_{i+1} - P_i) * B4_i(u)
                    const auto b4 = B4(u);
                    Eigen::Vector3d dpdu =
                        5.0 * ((CP[s][1] - CP[s][0]) * b4[0] + (CP[s][2] - CP[s][1]) * b4[1] + (CP[s][3] - CP[s][2]) * b4[2] + (CP[s][4] - CP[s][3]) * b4[3] + (CP[s][5] - CP[s][4]) * b4[4]);

                    const double speed = (dpdu / Ts).norm();
                    vmax_obs = std::max(vmax_obs, speed);
                }
            }
            vmax_used = vmax_obs;
        }

        visualization_msgs::msg::Marker mk;
        mk.header.stamp = rclcpp::Clock().now();
        mk.header.frame_id = frame_id;
        mk.ns = ns;
        mk.id = 1; // fixed id per-namespace so a new publish replaces the old line
        mk.type = visualization_msgs::msg::Marker::LINE_STRIP;
        mk.action = visualization_msgs::msg::Marker::ADD;
        mk.pose.orientation.w = 1.0;
        mk.scale.x = width; // use the provided width
        mk.color.r = r;     // fallback (RViz uses per-point if sizes match)
        mk.color.g = g;
        mk.color.b = b;
        mk.color.a = a;
        mk.lifetime = rclcpp::Duration(0, 0); // persistent

        mk.points.clear();
        mk.colors.clear();
        mk.points.reserve(CP.size() * (samples_per_seg + 1));
        mk.colors.reserve(CP.size() * (samples_per_seg + 1));

        for (size_t s = 0; s < CP.size(); ++s)
        {
            const int N = std::max(4, samples_per_seg);
            const double Ts = (s < T.size()) ? std::max(T[s], 1e-9) : 1.0;

            for (int k = 0; k <= N; ++k)
            {
                if (s > 0 && k == 0)
                    continue; // avoid duplicating seam point

                const double u = static_cast<double>(k) / N;

                const auto b5 = B5(u);
                Eigen::Vector3d p =
                    b5[0] * CP[s][0] + b5[1] * CP[s][1] + b5[2] * CP[s][2] +
                    b5[3] * CP[s][3] + b5[4] * CP[s][4] + b5[5] * CP[s][5];

                const auto b4 = B4(u);
                Eigen::Vector3d dpdu =
                    5.0 * ((CP[s][1] - CP[s][0]) * b4[0] + (CP[s][2] - CP[s][1]) * b4[1] + (CP[s][3] - CP[s][2]) * b4[2] + (CP[s][4] - CP[s][3]) * b4[3] + (CP[s][5] - CP[s][4]) * b4[4]);
                const double speed = (dpdu / Ts).norm();

                geometry_msgs::msg::Point q;
                q.x = p.x();
                q.y = p.y();
                q.z = p.z();
                mk.points.emplace_back(q);

                // Map speed -> color (warm = fast)
                std_msgs::msg::ColorRGBA c = getColorJet(speed, 0.0, vmax_used);
                mk.colors.emplace_back(c);
            }
        }

        // Sanity: ensure RViz uses per-point colors
        if (mk.colors.size() != mk.points.size())
        {
            mk.colors.clear(); // fall back to uniform color if something went off
        }

        mighty_traj_pub_->publish(mk);
    }

    inline void visualizeText(const Eigen::Vector3d &p,
                              const std::string &text,
                              double scale_z,
                              float r, float g, float b, float a,
                              const std::string &frame_id,
                              const std::string &ns,
                              int id,
                              double ttl_sec,
                              const std::string &planner = "gcopter")
    {
        visualization_msgs::msg::Marker m;
        m.header.stamp = rclcpp::Clock().now();
        m.header.frame_id = frame_id;
        m.ns = ns;
        m.id = id;
        m.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
        m.action = visualization_msgs::msg::Marker::ADD;

        // Slight offset upward so it doesn't overlap the sphere
        m.pose.position.x = p.x();
        m.pose.position.y = p.y();
        m.pose.position.z = p.z() + 1.0;
        m.pose.orientation.w = 1.0;

        // Only scale.z is used for text size
        m.scale.z = scale_z;

        m.color.r = r;
        m.color.g = g;
        m.color.b = b;
        m.color.a = a;

        m.text = text;

        // Short lifetime so publishing updates the same ID each tick
        m.lifetime = rclcpp::Duration::from_seconds(ttl_sec);

        if (planner == "gcopter")
        {
            text_pub_->publish(m);
        }
        else if (planner == "mighty")
        {
            mighty_text_pub_->publish(m);
        }
        else
        {
            RCLCPP_WARN(node_.get_logger(), "Unknown planner '%s' for text visualization", planner.c_str());
        }
    }

    inline void visualizePoints(const std::vector<Eigen::Vector3d> &pts,
                                float radius,
                                float r, float g, float b, float a,
                                const std::string &frame_id,
                                const std::string &ns,
                                double ttl_sec)
    {
        visualization_msgs::msg::Marker m;
        m.header.frame_id = frame_id;
        m.header.stamp = node_.now(); // assuming Visualizer stores node_
        m.ns = ns;
        m.id = 0;
        m.type = visualization_msgs::msg::Marker::SPHERE_LIST;
        m.action = visualization_msgs::msg::Marker::ADD;
        m.scale.x = m.scale.y = m.scale.z = radius;
        m.color.r = r;
        m.color.g = g;
        m.color.b = b;
        m.color.a = a;
        m.lifetime = rclcpp::Duration::from_seconds(ttl_sec);
        m.pose.orientation.w = 1.0;

        m.points.reserve(pts.size());
        for (const auto &P : pts)
        {
            geometry_msgs::msg::Point q;
            q.x = P.x();
            q.y = P.y();
            q.z = P.z();
            m.points.push_back(q);
        }
        mighty_waypoints_pub_->publish(m);
    }

    inline void visualizeSphereColor(const Eigen::Vector3d &p,
                                     float radius,
                                     float r, float g, float b, float a,
                                     const std::string &frame_id,
                                     const std::string &ns,
                                     int id,
                                     double ttl_sec)
    {
        visualization_msgs::msg::Marker m;
        m.header.frame_id = frame_id;
        m.header.stamp = node_.now();
        m.ns = ns;
        m.id = id;
        m.type = visualization_msgs::msg::Marker::SPHERE;
        m.action = visualization_msgs::msg::Marker::ADD;
        m.scale.x = m.scale.y = m.scale.z = 2.0 * radius;
        m.color.r = r;
        m.color.g = g;
        m.color.b = b;
        m.color.a = a;
        m.scale.x = m.scale.y = m.scale.z = radius * 2.0;
        m.lifetime = rclcpp::Duration::from_seconds(ttl_sec);
        m.pose.position.x = p.x();
        m.pose.position.y = p.y();
        m.pose.position.z = p.z();
        m.pose.orientation.w = 1.0;
        mighty_sphere_pub_->publish(m);
    }
};

#endif // VISUALIZER_HPP