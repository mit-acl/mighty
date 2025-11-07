#include "mighty/convert_vicon_to_state.hpp"

PoseTwistToStateNode::PoseTwistToStateNode() : Node("pose_twist_to_state_node")
{

    // define quality of service: all messages that you want to receive must have the same
    rmw_qos_profile_t custom_qos_profile = rmw_qos_profile_default;
    custom_qos_profile.depth = 1;
    custom_qos_profile.reliability = rmw_qos_reliability_policy_t::RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT;
    custom_qos_profile.history = rmw_qos_history_policy_t::RMW_QOS_POLICY_HISTORY_KEEP_LAST;
    custom_qos_profile.durability = rmw_qos_durability_policy_t::RMW_QOS_POLICY_DURABILITY_TRANSIENT_LOCAL;

    sub_pose_.subscribe(this, "world", custom_qos_profile);
    sub_twist_.subscribe(this, "twist", custom_qos_profile);

    // Time synchronizer setup
    pose_twist_sync_.reset(new Sync(MySyncPolicy(10), sub_pose_, sub_twist_));
    pose_twist_sync_->registerCallback(std::bind(&PoseTwistToStateNode::callback, this, std::placeholders::_1, std::placeholders::_2));

    // Publisher setup
    pub_state_ = this->create_publisher<dynus_interfaces::msg::State>("state", 10);
}

void PoseTwistToStateNode::callback(const std::shared_ptr<const geometry_msgs::msg::PoseStamped>& pose_msg,
                                    const std::shared_ptr<const geometry_msgs::msg::TwistStamped>& twist_msg)
{
    // Construct the State message
    auto state_msg = dynus_interfaces::msg::State();
    state_msg.header.stamp = this->get_clock()->now();
    state_msg.header.frame_id = pose_msg->header.frame_id;

    // Set position from PoseStamped
    state_msg.pos.x = pose_msg->pose.position.x;
    state_msg.pos.y = pose_msg->pose.position.y;
    state_msg.pos.z = pose_msg->pose.position.z;

    // Set velocity from TwistStamped
    state_msg.vel.x = twist_msg->twist.linear.x;
    state_msg.vel.y = twist_msg->twist.linear.y;
    state_msg.vel.z = twist_msg->twist.linear.z;

    // Set orientation from PoseStamped
    state_msg.quat.x = pose_msg->pose.orientation.x;
    state_msg.quat.y = pose_msg->pose.orientation.y;
    state_msg.quat.z = pose_msg->pose.orientation.z;
    state_msg.quat.w = pose_msg->pose.orientation.w;

    // Publish the State message
    pub_state_->publish(state_msg);
    RCLCPP_INFO(this->get_logger(), "Published State");
}

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<PoseTwistToStateNode>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
