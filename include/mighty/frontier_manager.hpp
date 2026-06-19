// Persistent global frontier database that survives the sliding map.
//
// The mapper publishes a robot-centered window that slides as the robot
// moves; cells that fall off the window are wiped to UNKNOWN. The detector
// can only find frontiers that are currently inside the window. This manager
// keeps a world-frame database of all frontiers seen so far, classifies them
// into ACTIVE/DORMANT/VISITED/INVALIDATED, and exposes a goal-selection API
// that ranks frontiers by an additive utility function.
//
// Threading: not thread-safe internally. All calls must come from the same
// callback group (in MIGHTY, occ2DCallback and the explore-select timer).

#pragma once

#include <Eigen/Core>
#include <cstdint>
#include <optional>
#include <vector>

#include "mighty/frontier_detector.hpp"
#include "mighty/occ_grid_2d.hpp"
#include "mighty/peer_tracker.hpp"

enum class FrontierState {
  ACTIVE,        // currently observable inside the local map
  DORMANT,       // outside the current local map, still pending
  VISITED,       // robot already explored it
  INVALIDATED,   // turned out to be unreachable or occupied
};

struct FrontierRecord {
  uint64_t        id            = 0;
  Eigen::Vector2d centroid_xy   = Eigen::Vector2d::Zero();   // world frame
  int             size_cells    = 0;
  double          first_seen_t  = 0.0;
  double          last_seen_t   = 0.0;
  FrontierState   state         = FrontierState::ACTIVE;
  int             visit_count   = 0;
  double          dwell_time_sec = 0.0;
  double          cached_utility = 0.0;
  Eigen::Vector2d aabb_min      = Eigen::Vector2d::Zero();
  Eigen::Vector2d aabb_max      = Eigen::Vector2d::Zero();
  // Pursuit deadline: absolute time (seconds) after which this record is
  // auto-invalidated if still ACTIVE/DORMANT. <=0 means "not being pursued".
  // Set by markSelected(), cleared on INVALIDATED/VISITED transition.
  double          pursuit_deadline_t = -1.0;
  // Total pursuit budget (seconds) allocated by markSelected(). Stored so
  // callers can display elapsed / total (e.g. "5.0s / 30.0s" in RViz).
  double          pursuit_budget_sec = 0.0;
  // Wall-clock time (seconds) when this record entered INVALIDATED. -1 if
  // never invalidated. Used by update() to suppress fresh clusters from
  // re-spawning a brand-new ACTIVE record right next to a frontier we just
  // gave up on. See invalidation_keep_out_radius_m / cooldown_sec params.
  double          invalidated_at_t = -1.0;
};

struct FrontierManagerParams {
  // Matching / lifecycle
  double merge_radius_m            = 1.0;
  double centroid_ema_alpha        = 0.5;
  double visit_radius_m            = 2.0;
  double visit_dwell_sec           = 1.0;
  int    verify_radius_cells       = 2;
  int    max_frontiers             = 1000;

  // Ranking weights (additive). All weights >= 0; set a weight to 0 to disable
  // the corresponding term.
  double w_size     = 1.0;
  double w_dist     = 2.0;
  double w_info     = 1.0;
  double w_revisit  = 0.5;
  double w_heading  = 0.3;

  // Normalizers
  double size_ref_m2     = 5.0;
  double dist_ref_m      = 25.0;
  double sensor_radius_m = 5.0;

  double goal_select_threshold = -1.0e9;  // -inf: always pick something

  // Pursuit timeout. When a frontier is selected as the current exploration
  // goal, it is given a deadline of
  //     max(pursuit_timeout_min_sec,
  //         dist / pursuit_timeout_v_ref * pursuit_timeout_factor)
  // seconds. If the deadline elapses while the record is still ACTIVE or
  // DORMANT, it is auto-INVALIDATED and the selector moves on. Set
  // pursuit_timeout_factor <= 0 to disable the feature entirely.
  double pursuit_timeout_factor  = 10.0;
  double pursuit_timeout_v_ref   = 0.5;   // reference velocity (m/s)
  double pursuit_timeout_min_sec = 10.0;  // floor, regardless of distance

  // Spatial keep-out around INVALIDATED records. A fresh cluster whose
  // centroid is within this radius of any INVALIDATED record (and inside
  // the cooldown window) is dropped instead of spawning a new ACTIVE record
  // — so HGP-unreachable / wall-hugger / pursuit-timeout decisions actually
  // stick instead of being immediately undone by the next detection cycle.
  // Set <= 0 to disable.
  double invalidation_keep_out_radius_m = 1.5;
  // Cooldown window (s). After this many seconds the keep-out lifts and the
  // area can be re-explored. Set <= 0 for permanent suppression (until the
  // record is evicted). Default 30s — enough that the agent picks a different
  // frontier first, short enough that transient HGP failures don't blacklist.
  double invalidation_cooldown_sec      = 30.0;

  // Peer-presence visit suppression: when an active peer's pose comes within
  // peer_visit_radius_m of an ACTIVE/DORMANT frontier centroid, flip that
  // record to VISITED immediately (sticky — the record stays VISITED for
  // the rest of the run). Robust to map-sharing failure under frame drift
  // because it relies only on per-peer poses (single-point sync), not on
  // grid alignment. Set <= 0 to disable.
  double peer_visit_radius_m            = 2.0;
};

class FrontierManager {
 public:
  explicit FrontierManager(FrontierManagerParams p) : params_(p) {}

  /** @brief Update the DB from a fresh detection batch.
   *  Match → EMA-update existing records; insert new ones; classify
   *  in-window-but-not-detected as VISITED/INVALIDATED; mark out-of-window
   *  records as DORMANT; apply robot-proximity dwell visit check.
   *
   *  @param fresh        Latest detector output.
   *  @param current_grid The grid the detector ran on (for in-bounds tests
   *                      and verify-cell queries).
   *  @param robot_pose   (x, y, yaw) in world frame.
   *  @param t_now        Current ROS time, seconds.
   */
  void update(const std::vector<FrontierCluster>& fresh,
              const OccGrid2D& current_grid,
              const Eigen::Vector3d& robot_pose,
              double t_now,
              const std::vector<PeerPose>& peers = {});

  /** @brief Pick the next exploration goal.
   *  Two-tier sort: ACTIVE first, then DORMANT. Skips VISITED/INVALIDATED.
   *  Returns nullopt if nothing scores above goal_select_threshold.
   */
  std::optional<FrontierRecord> selectNextGoal(
      const Eigen::Vector3d& robot_pose,
      const OccGrid2D& current_grid) const;

  /** @brief MinPos variant: pick the frontier where this robot has the lowest
   *  rank (fewest peers closer to it). When @p peers is empty, degenerates to
   *  nearest-frontier with utility tiebreak — identical to single-robot.
   */
  std::optional<FrontierRecord> selectNextGoalMinPos(
      const Eigen::Vector3d& robot_pose,
      const OccGrid2D& current_grid,
      const std::vector<PeerPose>& peers,
      double min_dist_to_peers_m = 0.0) const;

  void markVisited(uint64_t id);
  void markInvalidated(uint64_t id, double t_now);

  /** @brief Mark a frontier as the current pursuit goal and arm its timeout.
   *  Called by the exploration select tick right after picking a frontier.
   *  No-op if pursuit_timeout_factor <= 0.
   */
  void markSelected(uint64_t id, const Eigen::Vector2d& robot_xy, double t_now);

  const FrontierRecord* find(uint64_t id) const;
  const std::vector<FrontierRecord>& records() const { return records_; }
  const FrontierManagerParams& params() const { return params_; }

  // Test helpers
  size_t size() const { return records_.size(); }
  void clear() { records_.clear(); next_id_ = 0; last_update_t_ = 0.0; }

 private:
  double computeUtility(const FrontierRecord& r,
                        const Eigen::Vector3d& robot_pose,
                        const OccGrid2D& current_grid) const;
  bool   isInsideMap(const FrontierRecord& r, const OccGrid2D& grid) const;
  void   evictIfOverCap();

  FrontierManagerParams params_;
  std::vector<FrontierRecord> records_;
  uint64_t next_id_ = 0;
  double last_update_t_ = 0.0;
};
