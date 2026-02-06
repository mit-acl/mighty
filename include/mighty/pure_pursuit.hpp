#ifndef PURE_PURSUIT_HPP
#define PURE_PURSUIT_HPP

#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/point_stamped.hpp>
#include <visualization_msgs/msg/marker.hpp>
#include <dynus_interfaces/msg/state.hpp>
#include <dynus_interfaces/msg/trajectory.hpp>
#include <cmath>
#include <vector>

class PurePursuit : public rclcpp::Node
{
public:
    PurePursuit();

private:
    void trajectoryCallback(const dynus_interfaces::msg::Trajectory::SharedPtr msg);
    void stateCallback(const dynus_interfaces::msg::State::SharedPtr msg);
    void controlCallback();

    size_t findClosestWaypointIndex();
    size_t findLookaheadWaypointIndex(size_t start_idx, double lookahead_dist);
    double computeDynamicLookahead(double reference_velocity);
    double wrapPi(double angle);

    // Publishers and Subscribers
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr pub_cmd_vel_;
    rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr pub_lookahead_;
    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr pub_lookahead_marker_;
    rclcpp::Subscription<dynus_interfaces::msg::Trajectory>::SharedPtr sub_trajectory_;
    rclcpp::Subscription<dynus_interfaces::msg::State>::SharedPtr sub_state_;
    rclcpp::TimerBase::SharedPtr timer_;

    // State
    dynus_interfaces::msg::State current_state_;
    dynus_interfaces::msg::Trajectory trajectory_;
    bool state_initialized_;
    bool trajectory_initialized_;

    // Parameters
    double L_min_;        // Minimum lookahead distance (m)
    double k_v_;          // Velocity-dependent lookahead factor (s)
    double control_rate_; // Hz
    double max_velocity_; // Maximum commanded velocity (m/s)
    double max_angular_velocity_; // Maximum angular velocity (rad/s)
    double stopping_radius_;  // Stop when within this distance of goal (m)
    double adaptive_lookahead_distance_; // Start reducing lookahead when within this distance of goal (m)
    double turn_in_place_threshold_; // Heading error threshold for turn-in-place (rad)
    double slow_down_threshold_;     // Heading error threshold for speed reduction (rad)
    double w_smoothing_alpha_;       // Angular velocity smoothing factor (0 = no smoothing)
    double prev_w_command_ = 0.0;
};

#endif // PURE_PURSUIT_HPP
