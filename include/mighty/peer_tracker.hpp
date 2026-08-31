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

// A peer's currently-pursued frontier. Unlike PeerPose the id is kept: the
// claim tie-break needs to know WHO is contesting, not just where.
struct PeerClaim {
  std::string id;  // peer namespace, e.g. "NX02"
  Eigen::Vector2d position{Eigen::Vector2d::Zero()};
};

// Contested-claim tie-break: when two robots claim (frontiers within
// claim_block_radius of) the same spot, the lexicographically smaller
// namespace keeps it and the other yields. `<` on distinct strings is a
// strict total order, so exactly one side wins — both-yield and both-keep
// are impossible. Equal ids never contest: a robot's own messages are
// self-filtered at subscription (header.frame_id == ns_), so its own claim
// never enters its tracker.
inline bool winsClaimTie(const std::string& my_id, const std::string& peer_id) {
  return my_id < peer_id;
}

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

  // Identified positions of peers heard within the timeout. The yield
  // governor needs both WHO and WHERE: whether to slow for a nearby peer is
  // decided by the same namespace tie-break that settles claim contests, so
  // exactly one robot of any pair yields. Returns PeerClaim as a plain
  // (id, position) carrier.
  std::vector<PeerClaim> getActivePeersWithIds(double t_now,
                                               double timeout_sec) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<PeerClaim> out;
    for (const auto& [id, p] : peers_) {
      if (t_now - p.timestamp <= timeout_sec) {
        out.push_back({id, Eigen::Vector2d(p.x, p.y)});
      }
    }
    return out;
  }

  // Ids of peers heard within the timeout. getActivePeers deliberately drops
  // ids (MinPos rank only needs positions); the centralized-viz leader
  // election needs them — leader = lexicographic min of {own ns, active ids}.
  std::vector<std::string> getActivePeerIds(double t_now,
                                            double timeout_sec) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::string> out;
    for (const auto& [id, p] : peers_) {
      if (t_now - p.timestamp <= timeout_sec) {
        out.push_back(id);
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

// Thread-safe tracker for peer frontier claims, the intent half of MinPos
// coordination (PeerTracker above is the position half). Same shape and the
// same read-side staleness contract (<= timeout). One deliberate deviation:
// callers should feed updateClaim() the ARRIVAL time (receiver clock), not
// the message header stamp. Claims gate hard selection decisions — a peer
// with a skewed clock would otherwise have its claims either instantly stale
// (frontier stolen) or immortal (frontier blocked forever). Pose staleness
// tolerates skew because positions only bias ranking; claims must not.
class ClaimTracker {
 public:
  void updateClaim(const std::string& id, double x, double y, double t) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto& c = claims_[id];
    c.x = x;
    c.y = y;
    c.timestamp = t;
  }

  // A retarget overwrites the same key, so "re-claim IS the release": there
  // is no explicit un-claim message, and a robot that stops publishing
  // (finished, crashed, went home) ages out after timeout_sec.
  std::vector<PeerClaim> getActiveClaims(double t_now,
                                         double timeout_sec) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<PeerClaim> out;
    for (const auto& [id, c] : claims_) {
      if (t_now - c.timestamp <= timeout_sec) {
        out.push_back({id, Eigen::Vector2d(c.x, c.y)});
      }
    }
    return out;
  }

  size_t size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return claims_.size();
  }

 private:
  struct Entry {
    double x = 0.0;
    double y = 0.0;
    double timestamp = 0.0;
  };
  std::unordered_map<std::string, Entry> claims_;
  mutable std::mutex mutex_;
};
