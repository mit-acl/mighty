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
      if (params_.invalidation_keep_out_radius_m > 0.0) {
        bool suppressed = false;
        const double r2 = params_.invalidation_keep_out_radius_m
                        * params_.invalidation_keep_out_radius_m;
        for (const auto& rec : records_) {
          if (rec.state != FrontierState::INVALIDATED) continue;
          if (params_.invalidation_cooldown_sec > 0.0) {
            if (rec.invalidated_at_t < 0.0) continue;
            if (t_now - rec.invalidated_at_t
                > params_.invalidation_cooldown_sec) continue;
          }
          if ((c.centroid - rec.centroid_xy).squaredNorm() < r2) {
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
    double min_dist_to_peers_m) const {
  const Eigen::Vector2d robot_xy = robot_pose.head<2>();

  struct Candidate {
    int idx;
    int rank;
    double dist;
    double utility;
  };

  auto collectCandidates = [&](FrontierState want) -> std::vector<Candidate> {
    std::vector<Candidate> cands;
    for (size_t i = 0; i < records_.size(); ++i) {
      if (records_[i].state != want) continue;

      if (min_dist_to_peers_m > 0.0) {
        bool too_close = false;
        for (const auto& p : peers) {
          if ((p.position - records_[i].centroid_xy).norm() < min_dist_to_peers_m) {
            too_close = true;
            break;
          }
        }
        if (too_close) continue;
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
    auto cands = collectCandidates(tier);
    if (!cands.empty()) {
      FrontierRecord r = records_[cands[0].idx];
      r.cached_utility = cands[0].utility;
      return r;
    }
  }
  return std::nullopt;
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
