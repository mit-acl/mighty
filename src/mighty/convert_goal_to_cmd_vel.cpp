#include <mighty/convert_goal_to_cmd_vel.hpp>

GoalToCmdVel::GoalToCmdVel() : Node("goal_to_cmd_vel"),
                               current_yaw_(0.0),
                               kx_(1.0),
                               ky_(1.0),
                               kyaw_(1.0),
                               eps_(1e-2),
                               state_initialized_(false),
                               goal_initialized_(false)
{
    // Initialize state
    this->declare_parameter("x", 0.0);
    this->declare_parameter("y", 0.0);
    this->declare_parameter("z", 0.0);
    this->declare_parameter("yaw", 0.0);
    this->declare_parameter("cmd_vel_topic_name", "cmd_vel");
    this->declare_parameter("ground_robot_kx", 1.0);
    this->declare_parameter("ground_robot_ky", 1.0);
    this->declare_parameter("ground_robot_kyaw", 1.0);
    this->declare_parameter("ground_robot_eps", 1e-2);

    state_.pos.x = this->get_parameter("x").as_double();
    state_.pos.y = this->get_parameter("y").as_double();
    state_.pos.z = this->get_parameter("z").as_double();
    double yaw = this->get_parameter("yaw").as_double();

    // Convert yaw, pitch, roll to quaternion
    double pitch = 0.0, roll = 0.0;
    tf2::Quaternion quat;
    quat.setRPY(roll, pitch, yaw);

    state_.quat.x = quat.x();
    state_.quat.y = quat.y();
    state_.quat.z = quat.z();
    state_.quat.w = quat.w();

    // Topic names
    std::string cmd_vel_topic_name = this->get_parameter("cmd_vel_topic_name").as_string();

    // Gain parameters
    kx_ = this->get_parameter("ground_robot_kx").as_double();
    ky_ = this->get_parameter("ground_robot_ky").as_double();
    kyaw_ = this->get_parameter("ground_robot_kyaw").as_double();
    eps_ = this->get_parameter("ground_robot_eps").as_double();

    std::cout << "kx_: " << kx_ << std::endl;
    std::cout << "ky_: " << ky_ << std::endl;
    std::cout << "kyaw_: " << kyaw_ << std::endl;
    std::cout << "eps_: " << eps_ << std::endl;

    // Initialize goal
    goal_.p.x = 0.0;
    goal_.p.y = 0.0;
    goal_.p.z = 0.0;
    goal_.v.x = 0.0;
    goal_.v.y = 0.0;
    goal_.v.z = 0.0;
    goal_.a.x = 0.0;
    goal_.a.y = 0.0;
    goal_.a.z = 0.0;

    // Publishers and Subscribers
    pub_cmd_vel_ = this->create_publisher<geometry_msgs::msg::Twist>(cmd_vel_topic_name, 10);
    sub_goal_ = this->create_subscription<dynus_interfaces::msg::Goal>(
        "goal", 10, std::bind(&GoalToCmdVel::goalCallback, this, std::placeholders::_1));
    sub_state_ = this->create_subscription<dynus_interfaces::msg::State>(
        "state", 10, std::bind(&GoalToCmdVel::stateCallback, this, std::placeholders::_1));

    // Timers
    timer_ = this->create_wall_timer(
        std::chrono::milliseconds(10),
        std::bind(&GoalToCmdVel::cmdVelCallback, this));
}

void GoalToCmdVel::stateCallback(const dynus_interfaces::msg::State::SharedPtr msg)
{
    state_ = *msg;

    tf2::Quaternion quat(state_.quat.x, state_.quat.y, state_.quat.z, state_.quat.w);
    tf2::Matrix3x3(quat).getRPY(current_roll_, current_pitch_, current_yaw_);

    state_initialized_ = true;
}

void GoalToCmdVel::goalCallback(const dynus_interfaces::msg::Goal::SharedPtr msg)
{
    goal_ = *msg;
    goal_initialized_ = true;
}

void GoalToCmdVel::cmdVelCallback()
{
    if (!state_initialized_ || !goal_initialized_) return;

    geometry_msgs::msg::Twist twist;

    double x_desired = goal_.p.x;
    double y_desired = goal_.p.y;
    double xd_desired = goal_.v.x;
    double yd_desired = goal_.v.y;
    double yaw_desired = goal_.yaw;
    double yawd_desired = goal_.dyaw;

    double v_desired = std::sqrt(xd_desired * xd_desired + yd_desired * yd_desired);

    // Compute errors in desired frame
    double ex = std::cos(yaw_desired) * (state_.pos.x - x_desired) + std::sin(yaw_desired) * (state_.pos.y - y_desired);
    double ey = -std::sin(yaw_desired) * (state_.pos.x - x_desired) + std::cos(yaw_desired) * (state_.pos.y - y_desired);
    double eyaw = current_yaw_ - yaw_desired;

    // Compute control commands
    double v_command = v_desired * std::cos(eyaw) - kx_ * ex;
    double yawd_command = yawd_desired - v_desired * (ky_ * ey + std::sin(eyaw)) / std::sqrt(ey*ey + eps_*eps_) - kyaw_ * eyaw;

    // Clip control commands. This helped improve robustness; modify as needed
    // v_command = std::clamp(v_command, 0.0, 2.0);
    // yawd_command = std::clamp(yawd_command, -0.5, 0.5);

    twist.linear.x = v_command;
    twist.angular.z = yawd_command;

    pub_cmd_vel_->publish(twist);
}

double GoalToCmdVel::wrapPi(double x)
{
    x = std::fmod(x + M_PI, 2 * M_PI);
    return x - M_PI < 0 ? x - M_PI : x - M_PI;
}

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<GoalToCmdVel>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}