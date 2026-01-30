
#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/twist.hpp"
#include <mighty/mighty.hpp>
#include "dynus_interfaces/msg/goal.hpp"
#include "dynus_interfaces/msg/state.hpp"
#include "tf2/LinearMath/Quaternion.h"
#include "tf2/LinearMath/Matrix3x3.h"
#include <cmath>
#include <Eigen/Dense>

#ifndef GOAL_TO_CMD_VEL_HPP
#define GOAL_TO_CMD_VEL_HPP

class GoalToCmdVel : public rclcpp::Node
{
public:
    GoalToCmdVel();

private:
    void stateCallback(const dynus_interfaces::msg::State::SharedPtr msg);
    void goalCallback(const dynus_interfaces::msg::Goal::SharedPtr msg);
    void cmdVelCallback();
    double wrapPi(double x);

    // State and Goal
    dynus_interfaces::msg::State state_;
    dynus_interfaces::msg::Goal goal_;
    double current_roll_;
    double current_pitch_;
    double current_yaw_;

    // Publishers and Subscribers
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr pub_cmd_vel_;
    rclcpp::Subscription<dynus_interfaces::msg::Goal>::SharedPtr sub_goal_;
    rclcpp::Subscription<dynus_interfaces::msg::State>::SharedPtr sub_state_;

    // Timers
    rclcpp::TimerBase::SharedPtr timer_;

    // PID gains (load these from parameters as needed)
    double kp_x_{1.0}, ki_x_{0.0}, kd_x_{0.0};
    double kp_y_{1.0}, ki_y_{0.0}, kd_y_{0.0};
    double kp_yaw_{0.2}, ki_yaw_{0.1}, kd_yaw_{0.0};

    // PID state for x and y control
    double error_x_integral_{0.0}, error_y_integral_{0.0}, yaw_error_integral_{0.0};
    double previous_error_x_{0.0}, previous_error_y_{0.0}, previous_yaw_error_{0.0};

    // Time of previous update
    rclcpp::Time last_time_{this->now()};

    // Flags
    bool state_initialized_;
    bool goal_initialized_;

    // Filter
    double alpha_filter_ = 0.9;
    double prev_angular_z_ = 0.0;
};

#endif // GOAL_TO_CMD_VEL_HPP
