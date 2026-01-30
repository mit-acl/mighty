#include <mighty/quadruped_convert_goal_to_cmd_vel.hpp>

GoalToCmdVel::GoalToCmdVel() : Node("goal_to_cmd_vel"),
                                current_yaw_(0.0),
                                kp_x_(1.0),
                                ki_x_(0.0),
                                kd_x_(0.0),
                                kp_y_(1.0),
                                ki_y_(0.0),
                                kd_y_(0.0),
                                state_initialized_(false),
                                goal_initialized_(false)
{
    // Initialize state
    this->declare_parameter("x", 0.0);
    this->declare_parameter("y", 0.0);
    this->declare_parameter("z", 0.0);
    this->declare_parameter("yaw", 0.0);
    this->declare_parameter("cmd_vel_topic_name", "cmd_vel");
    this->declare_parameter("qudruped_kp_x", 1.0);
    this->declare_parameter("qudruped_ki_x", 0.0);
    this->declare_parameter("qudruped_kd_x", 0.0);
    this->declare_parameter("qudruped_kp_y", 1.0);
    this->declare_parameter("qudruped_ki_y", 0.0);
    this->declare_parameter("qudruped_kd_y", 0.0);

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
    kp_x_ = this->get_parameter("qudruped_kp_x").as_double();
    ki_x_ = this->get_parameter("qudruped_ki_x").as_double();
    kd_x_ = this->get_parameter("qudruped_kd_x").as_double();
    kp_y_ = this->get_parameter("qudruped_kp_y").as_double();
    ki_y_ = this->get_parameter("qudruped_ki_y").as_double();
    kd_y_ = this->get_parameter("qudruped_kd_y").as_double();

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

    if (!state_initialized_) state_initialized_ = true;
}

void GoalToCmdVel::goalCallback(const dynus_interfaces::msg::Goal::SharedPtr msg)
{
    goal_ = *msg;
    if (!goal_initialized_) goal_initialized_ = true;
}

void GoalToCmdVel::cmdVelCallback()
{
    if (!state_initialized_ || !goal_initialized_) return;

    // Get current time and compute time difference (dt)
    rclcpp::Time current_time = this->now();
    double dt = (current_time - last_time_).seconds();
    if (dt < 1e-3) { dt = 1e-3; }  // safeguard against too-small dt

    // --- Compute Errors in World Frame ---
    double error_x_world = goal_.p.x - state_.pos.x;
    double error_y_world = goal_.p.y - state_.pos.y;

    // --- Transform Errors into the Robot's Body Frame ---
    // This is useful if the cmd_vel is interpreted in the robot's frame.
    double cos_yaw = std::cos(current_yaw_);
    double sin_yaw = std::sin(current_yaw_);
    double error_x = cos_yaw * error_x_world + sin_yaw * error_y_world;
    double error_y = -sin_yaw * error_x_world + cos_yaw * error_y_world;

    // --- PID Control Computation ---
    // Update the integral terms
    error_x_integral_ += error_x * dt;
    error_y_integral_ += error_y * dt;

    // Compute derivative terms
    double derivative_x = (error_x - previous_error_x_) / dt;
    double derivative_y = (error_y - previous_error_y_) / dt;

    // Compute control outputs for x and y
    double control_x = kp_x_ * error_x + ki_x_ * error_x_integral_ + kd_x_ * derivative_x;
    double control_y = kp_y_ * error_y + ki_y_ * error_y_integral_ + kd_y_ * derivative_y;

    // Save current errors and time for next iteration
    previous_error_x_ = error_x;
    previous_error_y_ = error_y;
    last_time_ = current_time;

    // --- Compose and Publish the Command Velocity ---
    geometry_msgs::msg::Twist twist;
    twist.linear.x = control_x;
    twist.linear.y = control_y;
    twist.angular.z = goal_.dyaw;

    pub_cmd_vel_->publish(twist);
}


double GoalToCmdVel::wrapPi(double x)
{
    while (x > M_PI) x -= 2 * M_PI;
    while (x < -M_PI) x += 2 * M_PI;
    return x;
}

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<GoalToCmdVel>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
