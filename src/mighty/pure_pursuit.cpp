#include <mighty/pure_pursuit.hpp>

PurePursuit::PurePursuit()
    : Node("pure_pursuit"),
      state_initialized_(false),
      trajectory_initialized_(false)
{
    // Declare parameters
    this->declare_parameter("L_min", 0.5);
    this->declare_parameter("k_v", 0.5);
    this->declare_parameter("control_rate", 50.0);
    this->declare_parameter("max_velocity", 1.5);
    this->declare_parameter("max_angular_velocity", 3.0);  // Max turn rate (rad/s)
    this->declare_parameter("stopping_radius", 0.1);
    this->declare_parameter("adaptive_lookahead_distance", 2.0);
    this->declare_parameter("turn_in_place_threshold_deg", 60.0);  // degrees
    this->declare_parameter("slow_down_threshold_deg", 30.0);      // degrees
    this->declare_parameter("w_smoothing_alpha", 0.3);
    this->declare_parameter("use_hardware", false);

    // Get parameters
    L_min_ = this->get_parameter("L_min").as_double();
    k_v_ = this->get_parameter("k_v").as_double();
    control_rate_ = this->get_parameter("control_rate").as_double();
    max_velocity_ = this->get_parameter("max_velocity").as_double();
    max_angular_velocity_ = this->get_parameter("max_angular_velocity").as_double();
    stopping_radius_ = this->get_parameter("stopping_radius").as_double();
    adaptive_lookahead_distance_ = this->get_parameter("adaptive_lookahead_distance").as_double();
    turn_in_place_threshold_ = this->get_parameter("turn_in_place_threshold_deg").as_double() * M_PI / 180.0;
    slow_down_threshold_ = this->get_parameter("slow_down_threshold_deg").as_double() * M_PI / 180.0;
    w_smoothing_alpha_ = this->get_parameter("w_smoothing_alpha").as_double();
    use_hardware_ = this->get_parameter("use_hardware").as_bool();

    // Publishers and Subscribers
    std::string cmd_vel_string = use_hardware_ ? "cmd_vel_auto" : "cmd_vel";
    pub_cmd_vel_ = this->create_publisher<geometry_msgs::msg::Twist>(cmd_vel_string, 10);
    pub_lookahead_ = this->create_publisher<geometry_msgs::msg::PointStamped>("lookahead_point", 10);
    pub_lookahead_marker_ = this->create_publisher<visualization_msgs::msg::Marker>("lookahead_marker", 10);

    sub_trajectory_ = this->create_subscription<dynus_interfaces::msg::Trajectory>(
        "trajectory", 10,
        std::bind(&PurePursuit::trajectoryCallback, this, std::placeholders::_1));

    sub_state_ = this->create_subscription<dynus_interfaces::msg::State>(
        "state", 10,
        std::bind(&PurePursuit::stateCallback, this, std::placeholders::_1));

    // Control timer
    timer_ = this->create_wall_timer(
        std::chrono::milliseconds(static_cast<int>(1000.0 / control_rate_)),
        std::bind(&PurePursuit::controlCallback, this));
}

void PurePursuit::trajectoryCallback(const dynus_interfaces::msg::Trajectory::SharedPtr msg)
{
    trajectory_ = *msg;
    trajectory_initialized_ = true;
}

void PurePursuit::stateCallback(const dynus_interfaces::msg::State::SharedPtr msg)
{
    current_state_ = *msg;
    state_initialized_ = true;
}

size_t PurePursuit::findClosestWaypointIndex()
{
    if (trajectory_.goals.empty()) return 0;

    size_t closest_idx = 0;
    double min_dist_sq = std::numeric_limits<double>::max();

    for (size_t i = 0; i < trajectory_.goals.size(); ++i)
    {
        double dx = trajectory_.goals[i].p.x - current_state_.pos.x;
        double dy = trajectory_.goals[i].p.y - current_state_.pos.y;
        double dist_sq = dx * dx + dy * dy;

        if (dist_sq < min_dist_sq)
        {
            min_dist_sq = dist_sq;
            closest_idx = i;
        }
    }

    return closest_idx;
}

size_t PurePursuit::findLookaheadWaypointIndex(size_t start_idx, double lookahead_dist)
{
    if (trajectory_.goals.empty()) return 0;

    size_t lookahead_idx = start_idx;
    double accumulated_dist = 0.0;

    for (size_t i = start_idx; i < trajectory_.goals.size() - 1; ++i)
    {
        double dx = trajectory_.goals[i + 1].p.x - trajectory_.goals[i].p.x;
        double dy = trajectory_.goals[i + 1].p.y - trajectory_.goals[i].p.y;
        accumulated_dist += std::sqrt(dx * dx + dy * dy);

        if (accumulated_dist >= lookahead_dist)
        {
            lookahead_idx = i + 1;
            break;
        }
        lookahead_idx = i + 1;
    }

    return lookahead_idx;
}

double PurePursuit::computeDynamicLookahead(double reference_velocity)
{
    return L_min_ + k_v_ * std::abs(reference_velocity);
}

double PurePursuit::wrapPi(double angle)
{
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

void PurePursuit::controlCallback()
{
    if (!state_initialized_ || !trajectory_initialized_ || trajectory_.goals.empty())
    {
        return;
    }

    // Find closest waypoint on trajectory
    size_t closest_idx = findClosestWaypointIndex();
    const auto& closest_waypoint = trajectory_.goals[closest_idx];

    // Check if we've reached the end of the trajectory
    const auto& last_waypoint = trajectory_.goals.back();
    double dx_to_goal = last_waypoint.p.x - current_state_.pos.x;
    double dy_to_goal = last_waypoint.p.y - current_state_.pos.y;
    double dist_to_goal = std::sqrt(dx_to_goal * dx_to_goal + dy_to_goal * dy_to_goal);

    // Stop if we're close to the goal
    if (dist_to_goal < stopping_radius_)
    {
        geometry_msgs::msg::Twist twist;
        twist.linear.x = 0.0;
        twist.angular.z = 0.0;
        pub_cmd_vel_->publish(twist);
        prev_w_command_ = 0.0;
        return;
    }

    // Get initial velocity estimate from closest waypoint for lookahead calculation
    double v_closest = std::sqrt(closest_waypoint.v.x * closest_waypoint.v.x +
                                 closest_waypoint.v.y * closest_waypoint.v.y);

    // Use a minimum velocity for lookahead calculation to avoid tiny lookahead distances
    double v_for_lookahead = std::max(v_closest, 0.5);

    // Compute dynamic lookahead distance: L = L_min + k_v * v
    double lookahead_dist = computeDynamicLookahead(v_for_lookahead);

    // Adaptive lookahead: reduce lookahead distance when approaching goal for graceful stopping
    if (dist_to_goal < adaptive_lookahead_distance_)
    {
        double reduction_factor = dist_to_goal / adaptive_lookahead_distance_;
        lookahead_dist *= reduction_factor;
        lookahead_dist = std::max(lookahead_dist, 0.1);
    }

    // Find lookahead waypoint
    size_t lookahead_idx = findLookaheadWaypointIndex(closest_idx, lookahead_dist);
    const auto& lookahead_waypoint = trajectory_.goals[lookahead_idx];

    // Publish lookahead point for planner
    geometry_msgs::msg::PointStamped lookahead_msg;
    lookahead_msg.header.stamp = this->now();
    lookahead_msg.header.frame_id = "world";
    lookahead_msg.point.x = lookahead_waypoint.p.x;
    lookahead_msg.point.y = lookahead_waypoint.p.y;
    lookahead_msg.point.z = lookahead_waypoint.p.z;
    pub_lookahead_->publish(lookahead_msg);

    // Publish lookahead marker for RViz visualization
    visualization_msgs::msg::Marker marker;
    marker.header.stamp = this->now();
    marker.header.frame_id = "map";
    marker.ns = "lookahead";
    marker.id = 0;
    marker.type = visualization_msgs::msg::Marker::SPHERE;
    marker.action = visualization_msgs::msg::Marker::ADD;
    marker.pose.position.x = lookahead_waypoint.p.x;
    marker.pose.position.y = lookahead_waypoint.p.y;
    marker.pose.position.z = lookahead_waypoint.p.z;
    marker.pose.orientation.w = 1.0;
    marker.scale.x = 0.25;
    marker.scale.y = 0.25;
    marker.scale.z = 0.25;
    marker.color.r = 0.0;
    marker.color.g = 0.8;
    marker.color.b = 1.0;
    marker.color.a = 0.9;
    marker.lifetime = rclcpp::Duration::from_seconds(0.1);
    pub_lookahead_marker_->publish(marker);

    // ============================================================
    // CONTROL LAW
    // ============================================================

    // Get current yaw from quaternion
    double current_yaw = std::atan2(2.0 * (current_state_.quat.w * current_state_.quat.z +
                                            current_state_.quat.x * current_state_.quat.y),
                                    1.0 - 2.0 * (current_state_.quat.y * current_state_.quat.y +
                                                  current_state_.quat.z * current_state_.quat.z));

    // Compute heading error to lookahead point
    double dx = lookahead_waypoint.p.x - current_state_.pos.x;
    double dy = lookahead_waypoint.p.y - current_state_.pos.y;
    double heading_to_lookahead = std::atan2(dy, dx);
    double alpha = wrapPi(heading_to_lookahead - current_yaw);
    double abs_alpha = std::abs(alpha);

    // ---- Turn-in-place for large heading errors ----
    if (abs_alpha > turn_in_place_threshold_)
    {
        // Proportional turn rate based on heading error
        double w_turn = (alpha / abs_alpha) * max_angular_velocity_ * 0.9;

        // Smooth the angular command
        w_turn = w_smoothing_alpha_ * prev_w_command_ + (1.0 - w_smoothing_alpha_) * w_turn;
        prev_w_command_ = w_turn;

        geometry_msgs::msg::Twist twist;
        twist.linear.x = 0.0;
        twist.angular.z = w_turn;
        pub_cmd_vel_->publish(twist);
        return;
    }

    // ---- Speed reduction based on heading error ----
    double speed_reduction = 1.0;
    if (abs_alpha > slow_down_threshold_)
    {
        // Linear ramp from 1.0 at slow_down_threshold down to 0.0 at turn_in_place_threshold
        speed_reduction = std::max(0.0, 1.0 - (abs_alpha - slow_down_threshold_) /
                                                (turn_in_place_threshold_ - slow_down_threshold_));
    }

    // Also slow down near the goal for smooth stopping
    double goal_speed_factor = 1.0;
    if (dist_to_goal < adaptive_lookahead_distance_)
    {
        goal_speed_factor = std::max(0.1, dist_to_goal / adaptive_lookahead_distance_);
    }

    // Get reference velocity from lookahead waypoint
    double v_ref = std::sqrt(lookahead_waypoint.v.x * lookahead_waypoint.v.x +
                             lookahead_waypoint.v.y * lookahead_waypoint.v.y);
    v_ref = std::min(v_ref, max_velocity_);

    // Apply speed reductions
    double v_command = v_ref * speed_reduction * goal_speed_factor;

    // Ensure minimum velocity when not turning in place (avoid stalling)
    if (v_command < 0.05 && abs_alpha < turn_in_place_threshold_)
    {
        v_command = 0.05;
    }

    // ---- Pure pursuit curvature using lookahead_dist (not dist_to_lookahead) ----
    double curvature = (2.0 * std::sin(alpha)) / std::max(lookahead_dist, 0.1);

    // Angular velocity command: pure pursuit + heading correction term
    double w_pursuit = v_command * curvature;

    // Proportional heading correction for stronger steering at all speeds
    double k_yaw = 4.0;
    double w_heading = k_yaw * alpha;

    // Sum both terms for aggressive steering
    double w_command = w_pursuit + w_heading;

    // Clamp angular velocity
    w_command = std::clamp(w_command, -max_angular_velocity_, max_angular_velocity_);

    // ---- Smooth angular velocity to reduce jitter ----
    w_command = w_smoothing_alpha_ * prev_w_command_ + (1.0 - w_smoothing_alpha_) * w_command;
    prev_w_command_ = w_command;

    // Publish cmd_vel
    geometry_msgs::msg::Twist twist;
    twist.linear.x = v_command;
    twist.angular.z = w_command;
    pub_cmd_vel_->publish(twist);
}

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<PurePursuit>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
