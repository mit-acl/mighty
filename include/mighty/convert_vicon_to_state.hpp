#pragma once

#include <mighty/mighty.hpp>

#include "rclcpp/rclcpp.hpp"

#include "dynus_interfaces/msg/state.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "geometry_msgs/msg/twist_stamped.hpp"
#include "geometry_msgs/msg/vector3.hpp"
#include "message_filters/subscriber.h"
#include "message_filters/sync_policies/approximate_time.h"
#include "message_filters/sync_policies/exact_time.h"
#include "message_filters/synchronizer.h"
#include "std_msgs/msg/header.hpp"

// Define the synchronization policy
typedef message_filters::sync_policies::ApproximateTime<geometry_msgs::msg::PoseStamped,
                                                        geometry_msgs::msg::TwistStamped>
    MySyncPolicy;
// typedef message_filters::sync_policies::ExactTime<geometry_msgs::msg::PoseStamped,
// geometry_msgs::msg::TwistStamped> MySyncPolicy;
typedef message_filters::Synchronizer<MySyncPolicy> Sync;

/** @brief ROS 2 node that converts synchronized PoseStamped + TwistStamped to dynus_interfaces/State.
 *
 *  Subscribes to Vicon pose and twist topics via approximate time synchronization
 *  and republishes the fused data as a State message for the MIGHTY planner.
 */
class PoseTwistToStateNode : public rclcpp::Node {
 public:
  /** @brief Construct the node, set up synchronized subscribers and publisher. */
  PoseTwistToStateNode();

 private:
  void callback(const std::shared_ptr<const geometry_msgs::msg::PoseStamped>& pose_msg,
                const std::shared_ptr<const geometry_msgs::msg::TwistStamped>& twist_msg);

  /** @brief Return the z to publish, pinned to a constant when two_d_only_ is set.
   *  @param raw_z The z straight from the mocap PoseStamped.
   */
  double effectivePublishZ(double raw_z);

  // 2D-only z pinning. Mocap z is the noisiest axis for a ground robot that
  // never leaves the floor, and that noise propagates into the planner state.
  // When two_d_only_ is set, the published z is the average of the first
  // two_d_init_samples_ samples and then held constant for the rest of the run
  // (two_d_init_samples_ <= 0 locks on the very first sample). Mirrors DLIO's
  // output/twoDOnly. No mutex: the sync callback is the only writer and the
  // node spins single-threaded.
  bool   two_d_only_{false};
  int    two_d_init_samples_{50};
  bool   two_d_z_locked_{false};
  double two_d_z_const_{0.0};
  double two_d_z_sum_{0.0};
  int    two_d_z_count_{0};

  // Synchronized subscription
  message_filters::Subscriber<geometry_msgs::msg::PoseStamped> sub_pose_;
  message_filters::Subscriber<geometry_msgs::msg::TwistStamped> sub_twist_;

  // Time synchronizer
  std::shared_ptr<Sync> pose_twist_sync_;

  // Publisher
  rclcpp::Publisher<dynus_interfaces::msg::State>::SharedPtr pub_state_;
};
