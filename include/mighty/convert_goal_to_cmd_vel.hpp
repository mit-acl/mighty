
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

    // Control parameters
    double kv_;
    double kdist_;
    double kw_;
    double kyaw_;
    double kalpha_;
    double kx_;
    double ky_;
    double eps_;

    // Flags
    bool state_initialized_;
    bool goal_initialized_;
};

#endif // GOAL_TO_CMD_VEL_HPP
