// Thread-safe tracker for peer robot positions used by MinPos frontier
// allocation.  Each robot publishes its position on a shared ROS 2 topic;
// the subscriber callback feeds updatePeer(), and the explore-select timer
// reads getActivePeers().  The internal mutex serializes these two threads.

#pragma once

#include <Eigen/Core>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

struct PeerPose {
  Eigen::Vector2d position{Eigen::Vector2d::Zero()};
};

class PeerTracker {
 public:
  void updatePeer(const std::string& id, double x, double y, double t) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto& p = peers_[id];
    p.x = x;
    p.y = y;
    p.timestamp = t;
  }

  std::vector<PeerPose> getActivePeers(double t_now, double timeout_sec) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<PeerPose> out;
    for (const auto& [id, p] : peers_) {
      if (t_now - p.timestamp <= timeout_sec) {
        out.push_back({Eigen::Vector2d(p.x, p.y)});
      }
    }
    return out;
  }

  size_t size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return peers_.size();
  }

 private:
  struct Entry {
    double x = 0.0;
    double y = 0.0;
    double timestamp = 0.0;
  };
  std::unordered_map<std::string, Entry> peers_;
  mutable std::mutex mutex_;
};
