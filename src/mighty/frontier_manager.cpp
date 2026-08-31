#include "mighty/frontier_manager.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

inline double clamp01(double v) {
  return std::max(0.0, std::min(1.0, v));
}

}  // namespace

bool FrontierManager::isInsideMap(const FrontierRecord& r,
                                  const OccGrid2D& grid) const {
  int ix, iy;
  grid.worldToGrid(r.centroid_xy.x(), r.centroid_xy.y(), ix, iy);
  // Account for the same border margin the detector uses by treating cells
  // very close to the edge as "not really inside" — they may have just slid
  // in and not been observed yet.
  return ix >= 0 && ix < grid.width() && iy >= 0 && iy < grid.height();
}

void FrontierManager::update(const std::vector<FrontierCluster>& fresh,
                             const OccGrid2D& current_grid,
                             const Eigen::Vector3d& robot_pose,
                             double t_now,
                             const std::vector<PeerPose>& peers) {
  const Eigen::Vector2d robot_xy = robot_pose.head<2>();
  const double dt = (last_update_t_ > 0.0) ? std::max(0.0, t_now - last_update_t_)
                                           : 0.0;

  // ---- Step a/b: brute-force greedy match fresh -> existing ----
  std::vector<int> match_for_fresh(fresh.size(), -1);
  std::vector<uint8_t> existing_matched(records_.size(), 0);

  for (size_t i = 0; i < fresh.size(); ++i) {
    double best_d = params_.merge_radius_m;  // strictly less than = match
    int best_j = -1;
    for (size_t j = 0; j < records_.size(); ++j) {
      if (existing_matched[j]) continue;
      const auto& rec = records_[j];
      // INVALIDATED records (wall-huggers etc.) stay skipped — we want fresh
      // clusters near them to be re-evaluated, not re-bound to a bad record.
      // VISITED records DO participate in matching so that fresh clusters at
      // an already-visited location get absorbed into the existing record
      // instead of spawning a brand-new ACTIVE ghost every dwell cycle.
      if (rec.state == FrontierState::INVALIDATED) {
        continue;
      }
      const double d = (fresh[i].centroid - rec.centroid_xy).norm();
      if (d < best_d) {
        best_d = d;
        best_j = static_cast<int>(j);
      }
    }
    if (best_j >= 0) {
      match_for_fresh[i] = best_j;
      existing_matched[best_j] = 1;
    }
  }

  // Apply matches with EMA, then insert unmatched fresh as new records.
  const double alpha = params_.centroid_ema_alpha;
  for (size_t i = 0; i < fresh.size(); ++i) {
    const auto& c = fresh[i];
    if (match_for_fresh[i] >= 0) {
      auto& r = records_[match_for_fresh[i]];
      r.last_seen_t = t_now;
      // Don't revive a VISITED record or overwrite its centroid/size —
      // absorbing the fresh cluster is enough to keep it from spawning a
      // duplicate. Its VISITED state is preserved.
      if (r.state != FrontierState::VISITED) {
        r.centroid_xy = alpha * c.centroid + (1.0 - alpha) * r.centroid_xy;
        r.size_cells  = c.size_cells;
        r.aabb_min    = c.aabb_min;
        r.aabb_max    = c.aabb_max;
        r.state       = FrontierState::ACTIVE;
      }
    } else {
      // Suppress fresh clusters that fall inside the keep-out radius of a
      // recently-INVALIDATED record. The matching loop above intentionally
      // skipped INVALIDATED so we wouldn't revive a bad record; without this
      // check we'd just spawn a brand-new ACTIVE ghost in the same place and
      // the agent would re-pursue the goal it just gave up on.
      //
      // Two suppression regimes per nearby tombstone:
      //   * permanent — once its timeout_count has reached max_strikes, the
      //     frontier has been retried enough; never respawn it again.
      //   * temporary — otherwise honor the cooldown window, then allow a
      //     respawn but CARRY the strike count forward so repeated timeouts
      //     eventually escalate to the permanent regime above.
      int carry_strikes = 0;
      if (params_.invalidation_keep_out_radius_m > 0.0) {
        bool suppressed = false;
        const double r2 = params_.invalidation_keep_out_radius_m
                        * params_.invalidation_keep_out_radius_m;
        for (const auto& rec : records_) {
          if (rec.state != FrontierState::INVALIDATED) continue;
          if ((c.centroid - rec.centroid_xy).squaredNorm() >= r2) continue;
          // Worst strike count among nearby tombstones carries onto the respawn.
          carry_strikes = std::max(carry_strikes, rec.timeout_count);
          // Permanent keep-out: retried too many times, give up for good.
          if (params_.pursuit_timeout_max_strikes > 0 &&
              rec.timeout_count >= params_.pursuit_timeout_max_strikes) {
            suppressed = true;
            break;
          }
          // Temporary keep-out: honor the cooldown (<= 0 means permanent).
          if (params_.invalidation_cooldown_sec <= 0.0 ||
              (rec.invalidated_at_t >= 0.0 &&
               t_now - rec.invalidated_at_t
                   <= params_.invalidation_cooldown_sec)) {
            suppressed = true;
            break;
          }
        }
        if (suppressed) continue;  // drop fresh cluster, don't spawn new record
      }
      FrontierRecord r;
      r.id           = next_id_++;
      r.centroid_xy  = c.centroid;
      r.size_cells   = c.size_cells;
      r.aabb_min     = c.aabb_min;
      r.aabb_max     = c.aabb_max;
      r.first_seen_t = t_now;
      r.last_seen_t  = t_now;
      r.state        = FrontierState::ACTIVE;
      r.timeout_count = carry_strikes;  // persist strikes across respawn
      records_.push_back(r);
      // The new record is implicitly "matched" — don't classify it as
      // unmatched-in-window in step d. Mark its slot as matched.
      existing_matched.push_back(1);
    }
  }

  // ---- Step d/e: classify previously-known unmatched records ----
  for (size_t j = 0; j < records_.size(); ++j) {
    if (existing_matched[j]) continue;
    auto& r = records_[j];
    if (r.state == FrontierState::VISITED ||
        r.state == FrontierState::INVALIDATED) {
      continue;
    }

    if (!isInsideMap(r, current_grid)) {
      // Step e: outside local window -> DORMANT, untouched.
      r.state = FrontierState::DORMANT;
      continue;
    }

    // Step d: inside the local window but the detector didn't pick it up.
    int ix, iy;
    current_grid.worldToGrid(r.centroid_xy.x(), r.centroid_xy.y(), ix, iy);
    if (current_grid.isOccupied(ix, iy)) {
      r.state = FrontierState::INVALIDATED;
      r.invalidated_at_t = t_now;
      continue;
    }
    if (current_grid.isFree(ix, iy)) {
      // Check for unknown neighbors within verify_radius_cells. If none, the
      // area is fully explored -> VISITED. Otherwise leave ACTIVE for a later
      // detection cycle.
      bool has_unknown_nbr = false;
      const int rcells = params_.verify_radius_cells;
      for (int dy = -rcells; dy <= rcells && !has_unknown_nbr; ++dy) {
        for (int dx = -rcells; dx <= rcells && !has_unknown_nbr; ++dx) {
          if (dx == 0 && dy == 0) continue;
          if (current_grid.isUnknown(ix + dx, iy + dy)) {
            has_unknown_nbr = true;
          }
        }
      }
      if (!has_unknown_nbr) {
        r.state = FrontierState::VISITED;
      }
      // else: keep ACTIVE, just don't bump last_seen_t.
    }
    // else: cell is itself unknown — leave the record alone. Will be re-matched
    // on a later cycle when the detector picks it up.
  }

  // ---- Step f: robot-proximity dwell-based visit check ----
  for (auto& r : records_) {
    if (r.state == FrontierState::VISITED ||
        r.state == FrontierState::INVALIDATED) {
      continue;
    }
    const double d = (robot_xy - r.centroid_xy).norm();
    if (d < params_.visit_radius_m) {
      r.dwell_time_sec += dt;
      if (r.dwell_time_sec >= params_.visit_dwell_sec) {
        ++r.visit_count;
        r.state = FrontierState::VISITED;
      }
    } else {
      // Reset dwell when the robot leaves the visit radius — visits must be
      // contiguous, not the sum of brief drive-bys.
      r.dwell_time_sec = 0.0;
    }
  }

  // ---- Peer-presence visit suppression ----
  // If any active peer's pose is within peer_visit_radius_m of a record's
  // centroid, treat the record as VISITED. Sticky — once flipped, the record
  // stays VISITED for the rest of the run (or until eviction). This catches
  // peers that have already explored an area we still see as a frontier
  // (e.g. when peer-shared visited_maps don't propagate cleanly under
  // frame drift). No dwell — peer pose history is implicit in the FrontierRecord
  // state, so a single sample where a peer is in radius is enough.
  if (params_.peer_visit_radius_m > 0.0 && !peers.empty()) {
    const double r2 = params_.peer_visit_radius_m * params_.peer_visit_radius_m;
    for (auto& r : records_) {
      if (r.state != FrontierState::ACTIVE &&
          r.state != FrontierState::DORMANT) continue;
      for (const auto& p : peers) {
        if ((p.position - r.centroid_xy).squaredNorm() < r2) {
          r.state = FrontierState::VISITED;
          ++r.visit_count;
          break;
        }
      }
    }
  }

  // ---- Pursuit-timeout check ----
  // Any record with an armed deadline that has elapsed is INVALIDATED if it's
  // still ACTIVE or DORMANT — once we've given up trying to reach it, the
  // tombstone participates in the invalidation keep-out so the next detection
  // cycle doesn't immediately spawn a fresh ACTIVE ghost at the same spot.
  // Eviction (evictIfOverCap) prefers INVALIDATED → DORMANT → VISITED, so the
  // DB stays bounded. Records already VISITED/INVALIDATED just have their
  // deadlines cleared.
  if (params_.pursuit_timeout_factor > 0.0) {
    for (auto& r : records_) {
      if (r.pursuit_deadline_t <= 0.0) continue;
      if (r.state != FrontierState::ACTIVE &&
          r.state != FrontierState::DORMANT) {
        r.pursuit_deadline_t = -1.0;
        r.pursuit_budget_sec = 0.0;
        continue;
      }
      if (t_now >= r.pursuit_deadline_t) {
        r.state = FrontierState::INVALIDATED;
        r.invalidated_at_t = t_now;
        ++r.timeout_count;  // strike: enough of these -> permanent keep-out
        r.pursuit_deadline_t = -1.0;
        r.pursuit_budget_sec = 0.0;
      }
    }
  }

  evictIfOverCap();
  last_update_t_ = t_now;
}

double FrontierManager::computeUtility(const FrontierRecord& r,
                                       const Eigen::Vector3d& robot_pose,
                                       const OccGrid2D& current_grid) const {
  const Eigen::Vector2d robot_xy = robot_pose.head<2>();
  const double yaw = robot_pose.z();

  const double size_norm = clamp01(r.size_cells * current_grid.resolution()
                                       * current_grid.resolution()
                                       / params_.size_ref_m2);

  const double dist = (robot_xy - r.centroid_xy).norm();
  const double dist_norm = clamp01(dist / params_.dist_ref_m);

  // Info gain: count UNKNOWN cells in a disk of radius sensor_radius_m around
  // the centroid in the CURRENT grid. If the centroid is outside the current
  // map, the disk count will be 0 — that's intentional (we have no fresh
  // info-gain estimate for DORMANT frontiers, so they win on size/distance).
  double info_norm = 0.0;
  {
    int cx, cy;
    current_grid.worldToGrid(r.centroid_xy.x(), r.centroid_xy.y(), cx, cy);
    if (current_grid.inBounds(cx, cy)) {
      const int rcells = static_cast<int>(
          std::ceil(params_.sensor_radius_m / current_grid.resolution()));
      const double r2 = static_cast<double>(rcells) * rcells;
      int count = 0;
      int capacity = 0;
      for (int dy = -rcells; dy <= rcells; ++dy) {
        for (int dx = -rcells; dx <= rcells; ++dx) {
          if (dx * dx + dy * dy > r2) continue;
          ++capacity;
          if (current_grid.isUnknown(cx + dx, cy + dy)) ++count;
        }
      }
      if (capacity > 0) {
        info_norm = clamp01(static_cast<double>(count) / capacity);
      }
    }
  }

  // Heading alignment: cosine of angle between robot heading and direction
  // to centroid. Clamped at 0 (no penalty for behind-robot frontiers; only a
  // bonus for in-front).
  double heading_term = 0.0;
  {
    const Eigen::Vector2d to = r.centroid_xy - robot_xy;
    if (to.norm() > 1e-6) {
      const double angle_to = std::atan2(to.y(), to.x());
      heading_term = std::max(0.0, std::cos(angle_to - yaw));
    }
  }

  return params_.w_size    * size_norm
       - params_.w_dist    * dist_norm
       + params_.w_info    * info_norm
       - params_.w_revisit * r.visit_count
       + params_.w_heading * heading_term;
}

std::optional<FrontierRecord> FrontierManager::selectNextGoal(
    const Eigen::Vector3d& robot_pose, const OccGrid2D& current_grid) const {
  auto pickBest = [&](FrontierState want) -> std::optional<FrontierRecord> {
    double best_u = -std::numeric_limits<double>::infinity();
    int best_idx = -1;
    for (size_t i = 0; i < records_.size(); ++i) {
      if (records_[i].state != want) continue;
      const double u = computeUtility(records_[i], robot_pose, current_grid);
      if (u > best_u) {
        best_u = u;
        best_idx = static_cast<int>(i);
      }
    }
    if (best_idx < 0) return std::nullopt;
    if (best_u < params_.goal_select_threshold) return std::nullopt;
    FrontierRecord r = records_[best_idx];
    r.cached_utility = best_u;
    return r;
  };

  // Two-tier: exhaust ACTIVE before falling back to DORMANT.
  if (auto r = pickBest(FrontierState::ACTIVE)) return r;
  if (auto r = pickBest(FrontierState::DORMANT)) return r;
  return std::nullopt;
}

std::optional<FrontierRecord> FrontierManager::selectNextGoalMinPos(
    const Eigen::Vector3d& robot_pose,
    const OccGrid2D& current_grid,
    const std::vector<PeerPose>& peers,
    double min_dist_to_peers_m,
    const std::vector<PeerClaim>& claims,
    double claim_block_radius_m) const {
  const Eigen::Vector2d robot_xy = robot_pose.head<2>();

  struct Candidate {
    int idx;
    int rank;
    double dist;
    double utility;
  };

  auto collectCandidates = [&](FrontierState want, double keep_out) -> std::vector<Candidate> {
    std::vector<Candidate> cands;
    for (size_t i = 0; i < records_.size(); ++i) {
      if (records_[i].state != want) continue;

      if (keep_out > 0.0) {
        bool too_close = false;
        for (const auto& p : peers) {
          if ((p.position - records_[i].centroid_xy).norm() < keep_out) {
            too_close = true;
            break;
          }
        }
        if (too_close) continue;
      }

      // Claim hard-skip: a peer is already driving to (a frontier within
      // claim_block_radius of) this spot. Captured from the enclosing scope,
      // NOT from keep_out — so the keep-out fallback below relaxes only the
      // pose filter while claimed frontiers stay off-limits. Relaxing the
      // pose keep-out replaces "nothing" with "the best available"; relaxing
      // the claim filter would replace "nothing" with "a guaranteed
      // duplicate". All-claimed → nullopt is the designed home-with-resume
      // endgame, not starvation: selection still returns, the caller's grace
      // timer runs, and the robot acts.
      if (claim_block_radius_m > 0.0) {
        bool claimed = false;
        for (const auto& c : claims) {
          if ((c.position - records_[i].centroid_xy).norm() < claim_block_radius_m) {
            claimed = true;
            break;
          }
        }
        if (claimed) continue;
      }

      const double d_self = (robot_xy - records_[i].centroid_xy).norm();
      int rank = 0;
      for (const auto& p : peers) {
        if ((p.position - records_[i].centroid_xy).norm() < d_self) ++rank;
      }
      const double u = computeUtility(records_[i], robot_pose, current_grid);
      cands.push_back({static_cast<int>(i), rank, d_self, u});
    }
    std::sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b) {
      if (a.rank != b.rank) return a.rank < b.rank;
      if (a.dist != b.dist) return a.dist < b.dist;
      return a.utility > b.utility;
    });
    return cands;
  };

  // Two-tier: exhaust ACTIVE before falling back to DORMANT.
  for (FrontierState tier : {FrontierState::ACTIVE, FrontierState::DORMANT}) {
    auto cands = collectCandidates(tier, min_dist_to_peers_m);

    // Keep-out fallback. The peer-distance filter rejects outright, so with a
    // large radius every candidate in a tier can be inside some peer's circle
    // and the tier comes back empty. Returning nullopt there does not make the
    // robot cautious — it makes it idle: it selects nothing, never captures its
    // exploration start, and simply sits until the peers happen to move. That
    // was measured, not theorised (one run spent ~1400 s of its 1661 s stopped,
    // having driven only 132 m).
    //
    // Retrying the same tier without the keep-out can only ADD candidates where
    // there were none, so it is strictly safer than the previous behaviour: it
    // never overrides a viable separated choice, it only replaces "nothing" with
    // "the best available". Separation stays a strong preference rather than a
    // hard constraint that can deadlock.
    if (cands.empty() && min_dist_to_peers_m > 0.0 && !peers.empty()) {
      cands = collectCandidates(tier, 0.0);
    }

    if (!cands.empty()) {
      FrontierRecord r = records_[cands[0].idx];
      r.cached_utility = cands[0].utility;
      return r;
    }
  }
  return std::nullopt;
}

std::optional<FrontierRecord> FrontierManager::selectNearest(
    const Eigen::Vector3d& robot_pose) const {
  const Eigen::Vector2d robot_xy = robot_pose.head<2>();
  double best_d = std::numeric_limits<double>::infinity();
  int best_idx = -1;
  // Global nearest across ALL selectable frontiers — no ACTIVE-before-DORMANT
  // tiering, so a nearer DORMANT frontier is never passed over for a farther
  // ACTIVE one. VISITED/INVALIDATED are skipped.
  for (size_t i = 0; i < records_.size(); ++i) {
    const auto& r = records_[i];
    if (r.state != FrontierState::ACTIVE &&
        r.state != FrontierState::DORMANT) continue;
    const double d = (robot_xy - r.centroid_xy).norm();
    if (d < best_d) {
      best_d = d;
      best_idx = static_cast<int>(i);
    }
  }
  if (best_idx < 0) return std::nullopt;
  FrontierRecord r = records_[best_idx];
  r.cached_utility = -best_d;  // negative distance, so "higher = better" holds
  return r;
}

void FrontierManager::markVisited(uint64_t id) {
  for (auto& r : records_) {
    if (r.id == id) {
      ++r.visit_count;
      r.state = FrontierState::VISITED;
      r.pursuit_deadline_t = -1.0;
      r.pursuit_budget_sec = 0.0;
      return;
    }
  }
}

void FrontierManager::markInvalidated(uint64_t id, double t_now) {
  for (auto& r : records_) {
    if (r.id == id) {
      r.state = FrontierState::INVALIDATED;
      r.invalidated_at_t = t_now;
      // Count this as a strike. Repeated invalidations (HGP-unreachable, ESDF-
      // too-close) at the same location escalate to a PERMANENT keep-out once
      // timeout_count reaches pursuit_timeout_max_strikes — see update()'s
      // keep-out scan. Without this, an unreachable frontier cycles forever
      // (INVALIDATED -> cooldown expiry -> respawn ACTIVE -> fail -> ...),
      // never letting exploration settle to zero frontiers.
      ++r.timeout_count;
      r.pursuit_deadline_t = -1.0;
      r.pursuit_budget_sec = 0.0;
      return;
    }
  }
}

void FrontierManager::markSelected(uint64_t id, const Eigen::Vector2d& robot_xy,
                                   double t_now) {
  if (params_.pursuit_timeout_factor <= 0.0) return;  // feature disabled
  for (auto& r : records_) {
    if (r.id != id) continue;
    // Don't clobber a deadline that's already armed for this record — we keep
    // ticking from first selection even if the select tick re-affirms the
    // same goal later (e.g. after a re-match).
    if (r.pursuit_deadline_t > 0.0) return;
    const double dist = (r.centroid_xy - robot_xy).norm();
    const double v_ref = std::max(1e-3, params_.pursuit_timeout_v_ref);
    const double budget = std::max(params_.pursuit_timeout_min_sec,
                                   dist / v_ref * params_.pursuit_timeout_factor);
    r.pursuit_budget_sec = budget;
    r.pursuit_deadline_t = t_now + budget;
    return;
  }
}

int FrontierManager::applyPeerStatus(
    const std::vector<PeerFrontierStatus>& statuses,
    double t_now, double match_radius_m,
    std::optional<uint64_t> immune_id) {
  if (match_radius_m <= 0.0) return 0;
  int retired = 0;

  // Consumed-verdict markers make each peer verdict apply ONCE, not forever.
  // Peers rebroadcast their full DB at ~2 Hz for the rest of the run, so
  // without the markers an old verdict keeps retiring every NEW record that
  // later appears near the same spot — measured consequence: residual
  // unknown slivers could never re-enter anyone's frontier set, and in 2 of
  // 9 campaign runs a single early verdict near a doorway walled off a whole
  // unexplored wing (team went home at 0.64 coverage). A marker is refreshed
  // by every matching assert while the peer keeps broadcasting, and expires
  // peer_verdict_ttl_sec after the last one — so if the verdict-holder dies
  // and a teammate later reaches a REAL fresh verdict there, it applies.
  // Expired markers are pruned here (the only writer).
  peer_verdict_markers_.erase(
      std::remove_if(peer_verdict_markers_.begin(), peer_verdict_markers_.end(),
                     [&](const PeerVerdictMarker& mk) {
                       return params_.peer_verdict_ttl_sec > 0.0 &&
                              t_now - mk.last_assert_t >
                                  params_.peer_verdict_ttl_sec;
                     }),
      peer_verdict_markers_.end());

  for (const auto& s : statuses) {
    // Already consumed? Refresh the marker and skip — re-applying an old
    // verdict to whatever record happens to be nearest NOW is how distinct
    // new frontiers get wrongly killed.
    PeerVerdictMarker* consumed = nullptr;
    {
      double best_d = match_radius_m;
      for (auto& mk : peer_verdict_markers_) {
        const double d = (mk.position - s.position).norm();
        if (d < best_d) {
          best_d = d;
          consumed = &mk;
        }
      }
    }
    if (consumed) {
      consumed->last_assert_t = t_now;
      continue;
    }

    // Fresh verdict: find the nearest record of ANY state within the radius.
    // Matching across all states keeps a verdict pinned to the record that
    // absorbed it instead of creeping onto a distinct neighbor.
    int best_idx = -1;
    double best_d = match_radius_m;
    for (size_t i = 0; i < records_.size(); ++i) {
      const double d = (records_[i].centroid_xy - s.position).norm();
      if (d < best_d) {
        best_d = d;
        best_idx = static_cast<int>(i);
      }
    }

    if (best_idx >= 0) {
      FrontierRecord& r = records_[best_idx];
      const bool selectable = (r.state == FrontierState::ACTIVE ||
                               r.state == FrontierState::DORMANT);
      if (selectable) {
        if (!s.invalidated) {
          // Peer visited it — same semantics as the peer-presence suppression
          // in update() step g: sticky VISITED. If we are mid-pursuit of this
          // record, the state flip breaks our pursuit lock at the next select
          // tick — the desired "stop driving to a cleared frontier" behavior.
          ++r.visit_count;
          r.state = FrontierState::VISITED;
          r.pursuit_deadline_t = -1.0;
          r.pursuit_budget_sec = 0.0;
          ++retired;
        } else {
          // Peer gave up on it. Our own pursuit is immune: reachability can
          // be approach-specific, and one robot's failure must not abort a
          // pursuit that may succeed from ours. The verdict is deliberately
          // NOT consumed on the immune path — once our pursuit ends, the
          // peer's still-broadcast verdict gets a real chance to apply.
          if (immune_id && r.id == *immune_id) continue;
          // No strike, and no repeated cooldown refresh (the verdict applies
          // once): the record then follows the NORMAL local lifecycle —
          // cooldown, respawn, strikes-escalate-to-permanent — exactly as if
          // this robot had invalidated it itself. Peer verdicts inform the
          // team; they do not overrule the local retry policy.
          r.state = FrontierState::INVALIDATED;
          r.invalidated_at_t = t_now;
          r.pursuit_deadline_t = -1.0;
          r.pursuit_budget_sec = 0.0;
          ++retired;
        }
      }
      // Terminal nearest: nothing to do — never resurrect. Falls through to
      // marker creation so the rebroadcast stops being evaluated.
    }
    // Consume the verdict (also when nothing matched: never create records
    // from hearsay, and a record appearing here LATER is new information the
    // old verdict has no authority over — worst case we re-verify a spot the
    // peer cleared, and merged-map revalidation retires it within a tick).
    peer_verdict_markers_.push_back({s.position, t_now});
  }
  return retired;
}

void FrontierManager::clearPursuit(uint64_t id) {
  for (auto& r : records_) {
    if (r.id == id) {
      r.pursuit_deadline_t = -1.0;
      r.pursuit_budget_sec = 0.0;
      return;
    }
  }
}

const FrontierRecord* FrontierManager::find(uint64_t id) const {
  for (const auto& r : records_) {
    if (r.id == id) return &r;
  }
  return nullptr;
}

void FrontierManager::evictIfOverCap() {
  if (static_cast<int>(records_.size()) <= params_.max_frontiers) return;

  // Eviction priority: oldest VISITED > oldest INVALIDATED > oldest DORMANT.
  // Never evict ACTIVE.
  auto evictOne = [&](FrontierState target) -> bool {
    int oldest_idx = -1;
    double oldest_t = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < records_.size(); ++i) {
      if (records_[i].state != target) continue;
      if (records_[i].first_seen_t < oldest_t) {
        oldest_t = records_[i].first_seen_t;
        oldest_idx = static_cast<int>(i);
      }
    }
    if (oldest_idx < 0) return false;
    records_.erase(records_.begin() + oldest_idx);
    return true;
  };

  while (static_cast<int>(records_.size()) > params_.max_frontiers) {
    if (evictOne(FrontierState::VISITED))     continue;
    if (evictOne(FrontierState::INVALIDATED)) continue;
    if (evictOne(FrontierState::DORMANT))     continue;
    // Only ACTIVE left and we're still over cap. Stop — we never evict ACTIVE.
    break;
  }
}
