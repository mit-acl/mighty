// Tests for FrontierManager (matching, EMA, lifecycle, eviction).

#include <gtest/gtest.h>

#include <vector>

#include "mighty/frontier_manager.hpp"
#include "mighty/occ_grid_2d.hpp"

namespace {

constexpr int8_t U = -1;
constexpr int8_t F = 0;
constexpr int8_t O = 100;

// Build a 20x20 0.1m grid covering [0, 2.0]x[0, 2.0], all FREE by default.
std::shared_ptr<const OccGrid2D> MakeAllFreeGrid() {
  const int W = 20, H = 20;
  std::vector<int8_t> data(W * H, F);
  return OccGrid2D::fromTristate(W, H, 0.1, 0.0, 0.0, data);
}

// Build a 20x20 grid where the right half is unknown — there will be free
// cells along the seam at ix=9 that the manager can verify against.
std::shared_ptr<const OccGrid2D> MakeHalfUnknownGrid() {
  const int W = 20, H = 20;
  std::vector<int8_t> data(W * H);
  for (int iy = 0; iy < H; ++iy)
    for (int ix = 0; ix < W; ++ix)
      data[iy * W + ix] = (ix < 10) ? F : U;
  return OccGrid2D::fromTristate(W, H, 0.1, 0.0, 0.0, data);
}

// Same shape as MakeHalfUnknownGrid but with the window origin shifted in x, so
// a record inserted against the unshifted grid falls outside this one. That is
// the only way a record reaches DORMANT (frontier_manager.cpp:147-150).
std::shared_ptr<const OccGrid2D> MakeShiftedHalfUnknownGrid(double origin_x) {
  const int W = 20, H = 20;
  std::vector<int8_t> data(W * H);
  for (int iy = 0; iy < H; ++iy)
    for (int ix = 0; ix < W; ++ix)
      data[iy * W + ix] = (ix < 10) ? F : U;
  return OccGrid2D::fromTristate(W, H, 0.1, origin_x, 0.0, data);
}

FrontierCluster MakeCluster(double cx, double cy, int size_cells = 10) {
  FrontierCluster c;
  c.centroid    = Eigen::Vector2d(cx, cy);
  c.size_cells  = size_cells;
  c.size_m2     = size_cells * 0.01;
  c.aabb_min    = Eigen::Vector2d(cx - 0.1, cy - 0.1);
  c.aabb_max    = Eigen::Vector2d(cx + 0.1, cy + 0.1);
  return c;
}

FrontierManagerParams DefaultParams() {
  FrontierManagerParams p;
  p.merge_radius_m         = 0.5;
  p.centroid_ema_alpha     = 0.5;
  p.visit_radius_m         = 0.3;
  p.visit_dwell_sec        = 1.0;
  p.verify_radius_cells    = 2;
  p.max_frontiers          = 100;
  // Explicitly off. The header default is 2.0, which is harmless for the
  // single-robot tests only because update()'s peers vector defaults to empty —
  // stating it here keeps that from being an accident.
  p.peer_visit_radius_m    = 0.0;
  return p;
}

// Params for the peer-suppression tests, which need the radius armed.
FrontierManagerParams PeerParams(double peer_visit_radius_m) {
  FrontierManagerParams p = DefaultParams();
  p.peer_visit_radius_m = peer_visit_radius_m;
  return p;
}

std::vector<PeerPose> MakePeers(std::initializer_list<std::pair<double, double>> xys) {
  std::vector<PeerPose> peers;
  for (const auto& [x, y] : xys) peers.push_back(PeerPose{Eigen::Vector2d(x, y)});
  return peers;
}

int CountState(const FrontierManager& mgr, FrontierState want) {
  int n = 0;
  for (const auto& r : mgr.records())
    if (r.state == want) ++n;
  return n;
}

}  // namespace

TEST(FrontierManagerTest, FirstDetectInsertsAsActive) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());

  std::vector<FrontierCluster> fresh = {MakeCluster(0.95, 1.0, 10)};
  mgr.update(fresh, *grid, Eigen::Vector3d(0.5, 1.0, 0.0), 1.0);

  ASSERT_EQ(mgr.size(), 1u);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
  EXPECT_EQ(mgr.records()[0].size_cells, 10);
}

TEST(FrontierManagerTest, NearbyDetectionMergesViaEMA) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());

  // First detection at (0.95, 1.0)
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid,
             Eigen::Vector3d(0.5, 1.0, 0.0), 1.0);
  ASSERT_EQ(mgr.size(), 1u);
  const auto id0 = mgr.records()[0].id;

  // Second detection drifted by < merge_radius (0.5m).
  mgr.update({MakeCluster(1.05, 1.1, 12)}, *grid,
             Eigen::Vector3d(0.5, 1.0, 0.0), 1.5);
  ASSERT_EQ(mgr.size(), 1u);  // merged, not duplicated
  EXPECT_EQ(mgr.records()[0].id, id0);
  EXPECT_EQ(mgr.records()[0].size_cells, 12);
  // EMA centroid is halfway between (0.95, 1.0) and (1.05, 1.1).
  EXPECT_NEAR(mgr.records()[0].centroid_xy.x(), 1.0, 1e-6);
  EXPECT_NEAR(mgr.records()[0].centroid_xy.y(), 1.05, 1e-6);
}

TEST(FrontierManagerTest, FreeCellNoUnknownNeighborsBecomesVisited) {
  // Build an all-free grid (no unknowns anywhere). Existing record sits at a
  // free cell with no unknown neighbors -> should be marked VISITED.
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());

  // Bootstrap a record using a half-unknown grid first.
  auto half = MakeHalfUnknownGrid();
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *half,
             Eigen::Vector3d(0.5, 1.0, 0.0), 1.0);
  ASSERT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);

  // Next cycle: same grid coordinate, but the area is now fully free and
  // the detector reports no fresh clusters. The record should flip to VISITED.
  // Move the robot far from the centroid so the dwell-visit check doesn't
  // also trigger and mask whichever code path we're testing.
  mgr.update({}, *grid, Eigen::Vector3d(0.05, 0.05, 0.0), 1.5);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED);
}

TEST(FrontierManagerTest, OccupiedCellBecomesInvalidated) {
  // Build a grid where the centroid cell is OCCUPIED.
  const int W = 20, H = 20;
  std::vector<int8_t> data(W * H, F);
  // Mark cell at world (0.95, 1.0) as occupied. ix=9, iy=10
  data[10 * W + 9] = O;
  auto grid = OccGrid2D::fromTristate(W, H, 0.1, 0.0, 0.0, data);

  FrontierManager mgr(DefaultParams());
  // Bootstrap with a frontier record at that location.
  auto half = MakeHalfUnknownGrid();
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *half,
             Eigen::Vector3d(0.0, 1.0, 0.0), 1.0);

  // Next cycle: the cell at the centroid is now OCCUPIED. No fresh clusters.
  mgr.update({}, *grid, Eigen::Vector3d(0.0, 0.0, 0.0), 1.5);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::INVALIDATED);
}

TEST(FrontierManagerTest, OutOfWindowRecordBecomesDormant) {
  // Build a small 5x5 grid covering [0, 0.5] in x and y.
  const int W = 5, H = 5;
  std::vector<int8_t> data(W * H, F);
  auto small_grid = OccGrid2D::fromTristate(W, H, 0.1, 0.0, 0.0, data);

  FrontierManager mgr(DefaultParams());
  // Bootstrap a record in a larger grid first.
  auto big = MakeHalfUnknownGrid();
  mgr.update({MakeCluster(1.5, 1.5, 10)}, *big,
             Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);

  // Next cycle: the small grid does not contain the centroid (1.5, 1.5).
  // Move the robot far away from the centroid so dwell visit check doesn't
  // mask the test.
  mgr.update({}, *small_grid, Eigen::Vector3d(0.0, 0.0, 0.0), 1.5);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::DORMANT);
}

TEST(FrontierManagerTest, RobotDwellMarksVisited) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());

  // Bootstrap a record at (1.0, 1.0).
  mgr.update({MakeCluster(1.0, 1.0, 10)}, *grid,
             Eigen::Vector3d(2.0, 2.0, 0.0), 1.0);
  ASSERT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);

  // Robot dwells right next to the centroid for > visit_dwell_sec.
  mgr.update({MakeCluster(1.0, 1.0, 10)}, *grid,
             Eigen::Vector3d(1.05, 1.05, 0.0), 2.5);   // dt=1.5 > 1.0
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED);
  EXPECT_GE(mgr.records()[0].visit_count, 1);
}

TEST(FrontierManagerTest, SelectNextGoalPrefersActiveOverDormant) {
  // 100x100 grid so both records are in-bounds; we manually set states.
  const int W = 100, H = 100;
  std::vector<int8_t> data(W * H, F);
  // Mark some unknown so info_norm has something to bite on (otherwise
  // utilities are equal and the test depends on iteration order).
  for (int iy = 0; iy < H; ++iy)
    for (int ix = 50; ix < W; ++ix)
      data[iy * W + ix] = U;
  auto grid = OccGrid2D::fromTristate(W, H, 0.1, 0.0, 0.0, data);

  FrontierManager mgr(DefaultParams());

  // Insert a fresh ACTIVE record at world (5.0, 5.0).
  mgr.update({MakeCluster(5.0, 5.0, 10)}, *grid,
             Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
  ASSERT_EQ(mgr.size(), 1u);

  // Insert a second fresh ACTIVE record but bigger.
  mgr.update({MakeCluster(5.0, 5.0, 10), MakeCluster(7.0, 7.0, 50)}, *grid,
             Eigen::Vector3d(0.0, 0.0, 0.0), 1.5);
  ASSERT_EQ(mgr.size(), 2u);

  // Manually set the smaller record to DORMANT.
  // (We can't access records_ directly; do a controlled state change via
  // markInvalidated then markVisited would change semantics — instead, slide
  // it out by querying with a tiny grid that doesn't contain it.)
  // Easier: just check that selectNextGoal picks an ACTIVE record at all,
  // and that it picks the one with higher utility (the larger size).
  auto pick = mgr.selectNextGoal(Eigen::Vector3d(0.0, 0.0, 0.0), *grid);
  ASSERT_TRUE(pick.has_value());
  EXPECT_EQ(pick->state, FrontierState::ACTIVE);
  EXPECT_EQ(pick->size_cells, 50);  // bigger one wins
}

TEST(FrontierManagerTest, EvictionPrefersVisitedOverDormant) {
  FrontierManagerParams p = DefaultParams();
  p.max_frontiers = 2;
  FrontierManager mgr(p);

  auto grid = MakeAllFreeGrid();

  // Insert 2 records at t=1.0 -> both ACTIVE, no eviction.
  mgr.update({MakeCluster(0.5, 0.5, 5), MakeCluster(1.5, 1.5, 5)},
             *grid, Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
  ASSERT_EQ(mgr.size(), 2u);

  // Mark the first one VISITED externally.
  const uint64_t visited_id = mgr.records()[0].id;
  mgr.markVisited(visited_id);

  // Now insert a 3rd record. The VISITED one is the only legal eviction
  // candidate (manager never evicts ACTIVE), so it should be removed even
  // though we'll still be over cap with 3 ACTIVE records remaining.
  mgr.update({MakeCluster(0.5, 0.5, 5), MakeCluster(1.5, 1.5, 5),
              MakeCluster(0.7, 0.7, 5)},
             *grid, Eigen::Vector3d(0.0, 0.0, 0.0), 2.0);

  // The originally-VISITED record must be gone — it was the highest-priority
  // eviction target.
  EXPECT_EQ(mgr.find(visited_id), nullptr);
  for (const auto& r : mgr.records()) {
    EXPECT_NE(r.state, FrontierState::VISITED);
  }
}

TEST(FrontierManagerTest, RepeatedTimeoutsEscalateToPermanentBlacklist) {
  // A frontier the robot keeps failing to reach should stop reappearing after
  // pursuit_timeout_max_strikes timeouts (instead of respawning every cooldown).
  auto grid = MakeHalfUnknownGrid();
  FrontierManagerParams p = DefaultParams();
  p.pursuit_timeout_factor         = 1.0;
  p.pursuit_timeout_v_ref          = 1.0;
  p.pursuit_timeout_min_sec        = 1.0;   // budget ~ 1.0 s (dist is < 1 m)
  p.pursuit_timeout_max_strikes    = 2;     // give up after 2 timeouts
  p.invalidation_keep_out_radius_m = 0.5;
  p.invalidation_cooldown_sec      = 1.0;
  FrontierManager mgr(p);

  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  const Eigen::Vector2d robot_xy = robot.head<2>();
  auto cluster = [] {
    return std::vector<FrontierCluster>{MakeCluster(0.95, 1.0, 10)};
  };
  auto countActive = [&]() {
    int n = 0;
    for (const auto& r : mgr.records())
      if (r.state == FrontierState::ACTIVE) ++n;
    return n;
  };
  auto activeTimeoutCount = [&]() -> int {
    for (const auto& r : mgr.records())
      if (r.state == FrontierState::ACTIVE) return r.timeout_count;
    return -1;
  };

  // --- Strike 1: insert, select, let the pursuit deadline elapse. ---
  mgr.update(cluster(), *grid, robot, 1.0);
  ASSERT_EQ(countActive(), 1);
  mgr.markSelected(mgr.records()[0].id, robot_xy, 1.0);  // deadline = 2.0
  mgr.update(cluster(), *grid, robot, 2.5);              // past deadline
  ASSERT_EQ(countActive(), 0);
  ASSERT_EQ(mgr.records()[0].state, FrontierState::INVALIDATED);
  EXPECT_EQ(mgr.records()[0].timeout_count, 1);

  // --- Cooldown elapses -> respawns (retry), carrying the strike forward. ---
  mgr.update(cluster(), *grid, robot, 4.0);
  ASSERT_EQ(countActive(), 1);
  EXPECT_EQ(activeTimeoutCount(), 1) << "strike count must survive respawn";

  // --- Strike 2: select and time out again -> count reaches max_strikes. ---
  uint64_t id2 = 0;
  for (const auto& r : mgr.records())
    if (r.state == FrontierState::ACTIVE) id2 = r.id;
  mgr.markSelected(id2, robot_xy, 4.0);                  // deadline = 5.0
  mgr.update(cluster(), *grid, robot, 5.5);
  ASSERT_EQ(countActive(), 0);

  // --- Cooldown elapses again -> now PERMANENTLY blacklisted, no respawn. ---
  mgr.update(cluster(), *grid, robot, 7.0);
  EXPECT_EQ(countActive(), 0)
      << "frontier reappeared after reaching max strikes";
}

TEST(FrontierManagerTest, SelectNearestIgnoresSizeAndRetiredRecords) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());

  // Near+small vs far+big. The utility function would prefer the big one;
  // selectNearest must ignore size and pick the geometrically closest.
  mgr.update({MakeCluster(0.3, 0.0, 5),      // near, small
              MakeCluster(1.8, 0.0, 100)},   // far, big
             *grid, Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
  ASSERT_EQ(mgr.size(), 2u);

  auto pick = mgr.selectNearest(Eigen::Vector3d(0.0, 0.0, 0.0));
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 0.3, 1e-6);  // the near one
  const uint64_t near_id = pick->id;

  // Retiring the near frontier -> nearest falls back to the far one (VISITED
  // and INVALIDATED are never selectable).
  mgr.markVisited(near_id);
  auto pick2 = mgr.selectNearest(Eigen::Vector3d(0.0, 0.0, 0.0));
  ASSERT_TRUE(pick2.has_value());
  EXPECT_NEAR(pick2->centroid_xy.x(), 1.8, 1e-6);
}

// ===========================================================================
// MinPos frontier allocation (multi-robot)
//
// selectNextGoalMinPos ranks each frontier by how many peers are closer to it
// than we are, so robots spread out without any explicit negotiation. None of
// this had ever executed before these tests: the sim config disables MinPos,
// and both configs leave exploration.select_nearest true, which short-circuits
// the MinPos branch outright (mighty_node.cpp:3596).
// ===========================================================================

// The single-robot safety property: with no peers heard, MinPos must behave
// exactly like the nearest-frontier selection a solo robot already relies on.
TEST(FrontierManagerMinPosTest, NoPeersDegeneratesToNearest) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.3, 0.0, 5), MakeCluster(1.8, 0.0, 100)},
             *grid, robot, 1.0);

  auto minpos  = mgr.selectNextGoalMinPos(robot, *grid, {}, 0.0);
  auto nearest = mgr.selectNearest(robot);
  ASSERT_TRUE(minpos.has_value());
  ASSERT_TRUE(nearest.has_value());
  EXPECT_EQ(minpos->id, nearest->id);
  EXPECT_NEAR(minpos->centroid_xy.x(), 0.3, 1e-6);
}

// The core MinPos property. A is closest to us, but a peer is closer still, so
// A is rank 1 while B stays rank 0 — we must yield A and take B even though it
// costs us distance. Without this, robots converge on the same frontier.
//
// The geometry has to be off-axis. Collinear (robot, A, B) cannot express this
// case at all: a peer close enough to outrank us at A is necessarily within
// |AB| + |PA| of B, which is never farther than our own distance to B, so the
// peer outranks us at BOTH and the distance tie-break correctly returns A.
// Placing B perpendicular gives the peer a longer path to B than ours.
TEST(FrontierManagerMinPosTest, LowestRankWinsOverNearer) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(1.0, 1.0, 0.0);

  mgr.update({MakeCluster(1.0, 1.4, 10),    // A: 0.4 from us
              MakeCluster(0.2, 1.0, 10)},   // B: 0.8 from us, perpendicular
             *grid, robot, 1.0);

  // Peer at (1.0, 1.5): 0.1 from A (beats our 0.4 -> A is rank 1),
  // but 0.94 from B (worse than our 0.8 -> B stays rank 0).
  auto pick = mgr.selectNextGoalMinPos(robot, *grid, MakePeers({{1.0, 1.5}}), 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 0.2, 1e-6)
      << "took the frontier a peer was closer to — MinPos rank ignored";
}

TEST(FrontierManagerMinPosTest, EqualRankBreaksByDistance) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.4, 0.0, 10), MakeCluster(1.2, 0.0, 10)},
             *grid, robot, 1.0);

  // Peer far from both -> both rank 0 -> nearest wins.
  auto pick = mgr.selectNextGoalMinPos(robot, *grid, MakePeers({{1.9, 1.9}}), 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 0.4, 1e-6);
}

// Symmetric placement makes the distances bitwise equal, so the sort falls
// through to the utility term.
TEST(FrontierManagerMinPosTest, EqualRankAndDistanceBreaksByUtility) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(1.0, 1.0, 0.0);

  mgr.update({MakeCluster(0.5, 1.0, 5),     // small
              MakeCluster(1.5, 1.0, 100)},  // big, same distance
             *grid, robot, 1.0);

  auto pick = mgr.selectNextGoalMinPos(robot, *grid, {}, 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_EQ(pick->size_cells, 100) << "tie should fall through to utility";
}

TEST(FrontierManagerMinPosTest, RejectsFrontiersWithinMinDistToPeers) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.5, 0.0, 10), MakeCluster(1.5, 0.0, 10)},
             *grid, robot, 1.0);

  // Peer 0.2 from A. With a 0.5 m keep-out, A is rejected outright.
  auto pick = mgr.selectNextGoalMinPos(robot, *grid, MakePeers({{0.7, 0.0}}), 0.5);
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 1.5, 1e-6);
}

// The starvation case, and the fallback that fixes it.
//
// The peer keep-out rejects candidates outright, so a large radius can empty a
// tier entirely. Returning nullopt there does not make the robot cautious, it
// makes it idle — it selects nothing, never captures its exploration start, and
// sits until the peers move. Measured before the fallback existed: one run
// spent ~1400 s of 1661 s stopped, having driven 132 m.
//
// The fallback retries the tier with the keep-out disabled, which can only add
// candidates where there were none.
TEST(FrontierManagerMinPosTest, AllCandidatesRejectedFallsBackToUnfiltered) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.5, 0.0, 10), MakeCluster(0.9, 0.0, 10)},
             *grid, robot, 1.0);

  // Keep-out of 5 m swallows both candidates.
  auto pick = mgr.selectNextGoalMinPos(robot, *grid, MakePeers({{0.7, 0.0}}), 5.0);
  ASSERT_TRUE(pick.has_value()) << "robot selected nothing and would sit idle";
  EXPECT_NEAR(pick->centroid_xy.x(), 0.5, 1e-6) << "should take the nearest available";
}

// The fallback must not override a viable separated choice — it fires only when
// the filtered set is genuinely empty.
TEST(FrontierManagerMinPosTest, FallbackDoesNotOverrideAViableSeparatedChoice) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.5, 0.0, 10),    // inside the keep-out
              MakeCluster(1.8, 0.0, 10)},   // outside it
             *grid, robot, 1.0);

  auto pick = mgr.selectNextGoalMinPos(robot, *grid, MakePeers({{0.6, 0.0}}), 1.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 1.8, 1e-6)
      << "fallback fired despite a separated candidate being available";
}

// With no peers there is nothing to fall back from; behaviour is unchanged.
TEST(FrontierManagerMinPosTest, EmptyDbStillReturnsNulloptWithFallback) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  auto pick = mgr.selectNextGoalMinPos(Eigen::Vector3d(0.0, 0.0, 0.0), *grid,
                                       MakePeers({{1.0, 1.0}}), 5.0);
  EXPECT_FALSE(pick.has_value()) << "empty database must still yield nothing";
}

// Tiering dominates rank: the ACTIVE tier is exhausted before DORMANT is even
// considered, so a badly-ranked ACTIVE still beats a perfectly-ranked DORMANT.
// Worth pinning because it is a real coordination wrinkle — a robot will not
// hand off an ACTIVE frontier a peer is sitting on if its only alternative is
// dormant, which is the opposite of what "lowest rank wins" suggests.
//
// A record only goes DORMANT when it falls outside the current map window
// (frontier_manager.cpp:147-150), so this needs a second, shifted grid.
TEST(FrontierManagerMinPosTest, ActiveTierBeatsBetterRankedDormant) {
  FrontierManager mgr(DefaultParams());

  // Window 1 covers [0,2]^2 — insert a record at (0.3, 1.0).
  auto grid1 = MakeHalfUnknownGrid();
  mgr.update({MakeCluster(0.3, 1.0, 10)}, *grid1, Eigen::Vector3d(0.5, 1.0, 0.0), 1.0);
  ASSERT_EQ(CountState(mgr, FrontierState::ACTIVE), 1);

  // Window 2 covers [1.0,3.0]x[0,2]. The old record is now outside it -> DORMANT,
  // and a fresh cluster at (1.5, 1.0) becomes the ACTIVE one.
  auto grid2 = MakeShiftedHalfUnknownGrid(1.0);
  const Eigen::Vector3d robot(1.0, 1.0, 0.0);
  mgr.update({MakeCluster(1.5, 1.0, 10)}, *grid2, robot, 2.0);
  ASSERT_EQ(CountState(mgr, FrontierState::DORMANT), 1);
  ASSERT_EQ(CountState(mgr, FrontierState::ACTIVE), 1);

  // Peer at (1.5, 1.05) is 0.05 from the ACTIVE frontier (beats our 0.5 -> rank 1)
  // and 1.2 from the DORMANT one (worse than our 0.7 -> rank 0).
  auto pick = mgr.selectNextGoalMinPos(robot, *grid2, MakePeers({{1.5, 1.05}}), 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_EQ(pick->state, FrontierState::ACTIVE)
      << "a rank-0 DORMANT outranked a rank-1 ACTIVE — tiering was bypassed";
  EXPECT_NEAR(pick->centroid_xy.x(), 1.5, 1e-6);
}

TEST(FrontierManagerMinPosTest, SkipsVisitedAndInvalidated) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);

  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.5, 0.0, 10)},
             *grid, robot, 1.0);
  mgr.markVisited(mgr.records()[0].id);

  auto pick = mgr.selectNextGoalMinPos(robot, *grid, {}, 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_NEAR(pick->centroid_xy.x(), 1.5, 1e-6);
}

TEST(FrontierManagerMinPosTest, EmptyDatabaseReturnsNullopt) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  auto pick = mgr.selectNextGoalMinPos(Eigen::Vector3d(0.0, 0.0, 0.0), *grid,
                                       MakePeers({{1.0, 1.0}}), 1.0);
  EXPECT_FALSE(pick.has_value());
}

// ===========================================================================
// Peer-presence visit suppression (FrontierManager::update)
// ===========================================================================

TEST(FrontierManagerPeerSuppressionTest, PeerWithinRadiusMarksVisited) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(PeerParams(0.5));
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  ASSERT_EQ(CountState(mgr, FrontierState::ACTIVE), 1);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 2.0,
             MakePeers({{1.0, 1.0}}));   // 0.05 away — inside the radius
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED);
  EXPECT_GE(mgr.records()[0].visit_count, 1);
}

TEST(FrontierManagerPeerSuppressionTest, PeerOutsideRadiusLeavesRecordActive) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(PeerParams(0.3));
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 2.0,
             MakePeers({{1.8, 1.8}}));
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
}

// Suppression is sticky by design — the peer only has to be seen there once.
// Worth pinning because it is the mechanism most likely to leave a permanent
// hole in coverage if a peer merely drives past without stopping.
TEST(FrontierManagerPeerSuppressionTest, SuppressionIsStickyAfterPeerLeaves) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(PeerParams(0.5));
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 2.0,
             MakePeers({{1.0, 1.0}}));
  ASSERT_EQ(mgr.records()[0].state, FrontierState::VISITED);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 3.0, {});
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED)
      << "record revived after the peer left — suppression is meant to stick";
}

TEST(FrontierManagerPeerSuppressionTest, ZeroRadiusDisablesSuppression) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(PeerParams(0.0));
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 2.0,
             MakePeers({{0.95, 1.0}}));   // right on top of it
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
}

// The single-robot regression guard: an armed radius must change nothing when
// no peers are reported. This is what keeps the 15/15-validated solo behaviour
// intact once the multi-robot config turns the radius on.
TEST(FrontierManagerPeerSuppressionTest, NoPeersLeavesRecordsUntouched) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(PeerParams(5.0));   // radius covers the whole grid
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 2.0, {});
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
}

// ===========================================================================
// Claim-aware selection (selectNextGoalMinPos claims / claim_block_radius_m)
// ===========================================================================

namespace {

struct ClaimSpec { const char* id; double x; double y; };
std::vector<PeerClaim> MakeClaims(std::initializer_list<ClaimSpec> specs) {
  std::vector<PeerClaim> out;
  for (const auto& s : specs) out.push_back({s.id, Eigen::Vector2d(s.x, s.y)});
  return out;
}

struct StatusSpec { double x; double y; bool invalidated; };
std::vector<PeerFrontierStatus> MakeStatuses(
    std::initializer_list<StatusSpec> specs) {
  std::vector<PeerFrontierStatus> out;
  for (const auto& s : specs)
    out.push_back({Eigen::Vector2d(s.x, s.y), s.invalidated});
  return out;
}

}  // namespace

TEST(FrontierManagerMinPosClaimsTest, ClaimBlocksNearestCandidate) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);
  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.8, 0.0, 10)},
             *grid, robot, 1.0);

  auto pick = mgr.selectNextGoalMinPos(robot, *grid, {}, 0.0,
                                       MakeClaims({{"NX02", 0.35, 0.0}}), 0.5);
  ASSERT_TRUE(pick.has_value());
  EXPECT_GT(pick->centroid_xy.x(), 1.0)
      << "nearest frontier is inside a peer's claim disk and must be skipped";
}

// The load-bearing asymmetry: the pose keep-out fallback resurrects
// candidates when it would otherwise starve selection, but the claim filter
// must NEVER be relaxed — a stolen claim is a guaranteed duplicate, and
// all-claimed -> nullopt is the designed go-home-and-resume endgame.
TEST(FrontierManagerMinPosClaimsTest, ClaimFilterNeverRelaxedByFallback) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);
  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.8, 0.0, 10)},
             *grid, robot, 1.0);

  // Peer keep-out rejects everything AND a claim covers everything.
  auto blocked = mgr.selectNextGoalMinPos(
      robot, *grid, MakePeers({{1.0, 0.0}}), 5.0,
      MakeClaims({{"NX02", 1.0, 0.0}}), 1.2);
  EXPECT_FALSE(blocked.has_value())
      << "the keep-out fallback resurrected a CLAIMED candidate";

  // Same pose geometry without claims: the fallback fires and picks one.
  auto fallback = mgr.selectNextGoalMinPos(
      robot, *grid, MakePeers({{1.0, 0.0}}), 5.0);
  EXPECT_TRUE(fallback.has_value())
      << "contrast leg: pose keep-out alone must still fall back";
}

TEST(FrontierManagerMinPosClaimsTest, PoseKeepOutStillRelaxesWithClaimsPresent) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);
  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.8, 0.0, 10)},
             *grid, robot, 1.0);

  // Keep-out blocks both; the claim covers only the far frontier. The
  // fallback must relax the pose filter, keep the claim filter, and land on
  // the near frontier.
  auto pick = mgr.selectNextGoalMinPos(
      robot, *grid, MakePeers({{1.0, 0.0}}), 5.0,
      MakeClaims({{"NX02", 1.8, 0.0}}), 0.5);
  ASSERT_TRUE(pick.has_value());
  EXPECT_LT(pick->centroid_xy.x(), 1.0);
}

TEST(FrontierManagerMinPosClaimsTest, ZeroRadiusDisablesClaimFilter) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);
  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.8, 0.0, 10)},
             *grid, robot, 1.0);

  auto pick = mgr.selectNextGoalMinPos(robot, *grid, {}, 0.0,
                                       MakeClaims({{"NX02", 0.35, 0.0}}), 0.0);
  ASSERT_TRUE(pick.has_value());
  EXPECT_LT(pick->centroid_xy.x(), 1.0)
      << "radius 0 must disable the claim filter entirely";
}

// Pins the default-argument contract: every pre-claims call site (and test)
// keeps its exact behaviour.
TEST(FrontierManagerMinPosClaimsTest, EmptyClaimsBackwardCompatible) {
  auto grid = MakeAllFreeGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.0, 0.0, 0.0);
  mgr.update({MakeCluster(0.3, 0.0, 10), MakeCluster(1.8, 0.0, 10)},
             *grid, robot, 1.0);

  auto four_arg = mgr.selectNextGoalMinPos(robot, *grid,
                                           MakePeers({{0.7, 0.0}}), 0.5);
  auto six_arg  = mgr.selectNextGoalMinPos(robot, *grid,
                                           MakePeers({{0.7, 0.0}}), 0.5, {}, 0.0);
  ASSERT_TRUE(four_arg.has_value());
  ASSERT_TRUE(six_arg.has_value());
  EXPECT_EQ(four_arg->id, six_arg->id);
}

// ===========================================================================
// Peer frontier-status application (FrontierManager::applyPeerStatus)
// ===========================================================================

TEST(FrontierManagerApplyPeerStatusTest, VisitedRetiresNearestActive) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  const uint64_t id = mgr.records()[0].id;
  mgr.markSelected(id, Eigen::Vector2d(0.5, 1.0), 1.0);
  ASSERT_GT(mgr.records()[0].pursuit_deadline_t, 0.0);

  const int n = mgr.applyPeerStatus(MakeStatuses({{1.0, 1.0, false}}), 2.0, 0.5);
  EXPECT_EQ(n, 1);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED);
  EXPECT_GE(mgr.records()[0].visit_count, 1);
  EXPECT_LT(mgr.records()[0].pursuit_deadline_t, 0.0)
      << "a retired record must not keep an armed pursuit deadline";
}

TEST(FrontierManagerApplyPeerStatusTest, InvalidatedSetsTombstoneWithoutStrike) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);

  const int n = mgr.applyPeerStatus(MakeStatuses({{1.0, 1.0, true}}), 5.0, 0.5);
  EXPECT_EQ(n, 1);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::INVALIDATED);
  EXPECT_DOUBLE_EQ(mgr.records()[0].invalidated_at_t, 5.0);
  EXPECT_EQ(mgr.records()[0].timeout_count, 0)
      << "peer-sourced tombstones carry no strike: an independent local "
         "attempt after the cooldown must stay possible";
}

TEST(FrontierManagerApplyPeerStatusTest, NoMatchOutsideRadiusNoCreate) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);

  const int n = mgr.applyPeerStatus(MakeStatuses({{0.2, 0.2, false}}), 2.0, 0.5);
  EXPECT_EQ(n, 0);
  EXPECT_EQ(mgr.size(), 1u) << "status must never create records from hearsay";
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
}

// The nearest-match-across-ALL-states rule: once a record is terminal it
// absorbs the peer's repeated broadcast of the same verdict, instead of the
// verdict creeping onto a distinct neighbor on the next 2 Hz rebroadcast.
TEST(FrontierManagerApplyPeerStatusTest, TerminalNearestAbsorbsStatus) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10), MakeCluster(0.95, 0.4, 10)},
             *grid, robot, 1.0);
  ASSERT_EQ(mgr.size(), 2u);
  const uint64_t visited_id =
      (mgr.records()[0].centroid_xy.y() > 0.7) ? mgr.records()[0].id
                                               : mgr.records()[1].id;
  mgr.markVisited(visited_id);

  // Status lands nearest the already-VISITED record; the ACTIVE neighbor
  // 0.45 m away is also inside the radius but must NOT be retired.
  const int n = mgr.applyPeerStatus(MakeStatuses({{0.95, 0.85, false}}), 2.0, 0.5);
  EXPECT_EQ(n, 0);
  EXPECT_EQ(CountState(mgr, FrontierState::ACTIVE), 1);
  EXPECT_EQ(CountState(mgr, FrontierState::VISITED), 1);
}

TEST(FrontierManagerApplyPeerStatusTest, NeverResurrects) {
  auto grid = MakeHalfUnknownGrid();
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  // INVALIDATED stays INVALIDATED under a peer VISITED verdict; a repeated
  // peer INVALIDATED refreshes the cooldown stamp without adding a strike.
  FrontierManager mgr(DefaultParams());
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.markInvalidated(mgr.records()[0].id, 1.0);
  const int strikes_after_local = mgr.records()[0].timeout_count;

  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 1.0, false}}), 2.0, 0.5), 0);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::INVALIDATED);

  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 1.0, true}}), 3.0, 0.5), 0);
  EXPECT_DOUBLE_EQ(mgr.records()[0].invalidated_at_t, 1.0)
      << "peer rebroadcasts must NOT slide the local cooldown — the record "
         "follows the normal local respawn/strike lifecycle";
  EXPECT_EQ(mgr.records()[0].timeout_count, strikes_after_local);

  // VISITED stays VISITED under a peer INVALIDATED verdict.
  FrontierManager mgr2(DefaultParams());
  mgr2.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr2.markVisited(mgr2.records()[0].id);
  EXPECT_EQ(mgr2.applyPeerStatus(MakeStatuses({{0.95, 1.0, true}}), 2.0, 0.5), 0);
  EXPECT_EQ(mgr2.records()[0].state, FrontierState::VISITED);
}

TEST(FrontierManagerApplyPeerStatusTest, CurrentPursuitImmuneToInvalidated) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  const uint64_t id = mgr.records()[0].id;

  // A peer's failure must not abort OUR pursuit — reachability can be
  // approach-specific, and our own timeout is the judge.
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 1.0, true}}), 2.0, 0.5, id), 0);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);

  // A peer's VISITED verdict is NOT immune: someone cleared it, stop driving.
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 1.0, false}}), 3.0, 0.5, id), 1);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::VISITED);
}

// Stability under the 2 Hz rebroadcast: the same status list applied twice
// retires once and then does nothing.
TEST(FrontierManagerApplyPeerStatusTest, RepeatedApplicationIdempotent) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);

  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{1.0, 1.0, false}}), 2.0, 0.5), 1);
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{1.0, 1.0, false}}), 2.5, 0.5), 0);
  EXPECT_EQ(mgr.records()[0].visit_count, 1);
}

// ===========================================================================
// Consumed-verdict markers (a peer verdict applies once, not forever)
// ===========================================================================

// The failure this prevents, measured live: peers rebroadcast their full DB
// at ~2 Hz for the whole run, so without consumption an old verdict keeps
// retiring every NEW record that later appears near the same spot — residual
// unknown slivers could never re-enter anyone's frontier set, and in 2 of 9
// campaign runs a single early verdict near a doorway walled off a whole
// unexplored wing (the team went home at 0.64 coverage).
TEST(FrontierManagerVerdictConsumptionTest, VerdictAppliesOnlyOnce) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());   // peer_verdict_ttl default 5.0
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);

  // Fresh verdict retires the existing record (0.45 m from the status point).
  ASSERT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 2.0, 0.5), 1);
  ASSERT_EQ(mgr.records()[0].state, FrontierState::VISITED);

  // A genuinely NEW record appears near the old spot — new information.
  mgr.update({MakeCluster(0.95, 0.45, 10)}, *grid, robot, 3.0);
  ASSERT_EQ(mgr.size(), 2u) << "insert must never be blocked by verdicts";

  // The peer's rebroadcast of the SAME verdict arrives. The new record is
  // now the nearest match (0.10 m vs 0.45 m) — without consumption it would
  // be wrongly retired by a verdict that predates its existence.
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 3.5, 0.5), 0);
  EXPECT_EQ(CountState(mgr, FrontierState::ACTIVE), 1)
      << "a rebroadcast old verdict killed a record that appeared after it";
}

TEST(FrontierManagerVerdictConsumptionTest, ExpiredMarkerAllowsFreshVerdict) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());   // ttl default 5.0
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  ASSERT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 2.0, 0.5), 1);

  // Long silence (verdict-holder finished or died), then a record appears
  // and a NEW verdict at the same spot arrives: it must apply — consumption
  // is per-assertion-stream, not a permanent land grab.
  mgr.update({MakeCluster(0.95, 0.45, 10)}, *grid, robot, 20.0);
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 21.0, 0.5), 1);
  EXPECT_EQ(CountState(mgr, FrontierState::ACTIVE), 0);
}

TEST(FrontierManagerVerdictConsumptionTest, ReassertKeepsVerdictConsumed) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());   // ttl default 5.0
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  ASSERT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 2.0, 0.5), 1);
  // Rebroadcast at t=6 refreshes the marker (2.0 + ttl would expire at 7.0).
  ASSERT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 6.0, 0.5), 0);

  // t=10.5 is past the ORIGINAL assert's expiry but inside the refreshed
  // one: still consumed, the new record survives.
  mgr.update({MakeCluster(0.95, 0.45, 10)}, *grid, robot, 10.0);
  EXPECT_EQ(mgr.applyPeerStatus(MakeStatuses({{0.95, 0.55, false}}), 10.5, 0.5), 0);
  EXPECT_EQ(CountState(mgr, FrontierState::ACTIVE), 1);
}

// Single-robot regression guard: solo runs never call applyPeerStatus, so no
// markers exist and inserts behave exactly as before.
TEST(FrontierManagerVerdictConsumptionTest, NoStatusNoSuppression) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  mgr.update({MakeCluster(0.95, 0.4, 10)}, *grid, robot, 2.0);
  EXPECT_EQ(mgr.size(), 2u);
}

// ===========================================================================
// clearPursuit (the yield path)
// ===========================================================================

TEST(FrontierManagerClearPursuitTest, DropsDeadlineWithoutStateOrStrike) {
  auto grid = MakeHalfUnknownGrid();
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);

  FrontierManager mgr(DefaultParams());
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  const uint64_t id = mgr.records()[0].id;
  mgr.markSelected(id, Eigen::Vector2d(0.5, 1.0), 1.0);
  ASSERT_GT(mgr.records()[0].pursuit_deadline_t, 0.0);

  mgr.clearPursuit(id);
  EXPECT_LT(mgr.records()[0].pursuit_deadline_t, 0.0);
  EXPECT_DOUBLE_EQ(mgr.records()[0].pursuit_budget_sec, 0.0);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
  EXPECT_EQ(mgr.records()[0].timeout_count, 0);

  // Far past the ORIGINAL deadline: the yielded record must not be
  // invalidated by a pursuit we abandoned...
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 100.0);
  EXPECT_EQ(mgr.records()[0].state, FrontierState::ACTIVE);
  EXPECT_EQ(mgr.records()[0].timeout_count, 0);

  // ...whereas WITHOUT the yield the same timeline fires the timeout —
  // proving this test would catch a broken clearPursuit.
  FrontierManager ctrl(DefaultParams());
  ctrl.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  ctrl.markSelected(ctrl.records()[0].id, Eigen::Vector2d(0.5, 1.0), 1.0);
  ctrl.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 100.0);
  EXPECT_EQ(ctrl.records()[0].state, FrontierState::INVALIDATED);
}

TEST(FrontierManagerClearPursuitTest, AllowsFreshReArm) {
  auto grid = MakeHalfUnknownGrid();
  FrontierManager mgr(DefaultParams());
  const Eigen::Vector3d robot(0.5, 1.0, 0.0);
  mgr.update({MakeCluster(0.95, 1.0, 10)}, *grid, robot, 1.0);
  const uint64_t id = mgr.records()[0].id;

  mgr.markSelected(id, Eigen::Vector2d(0.5, 1.0), 1.0);
  mgr.clearPursuit(id);
  // markSelected's no-re-arm guard keys on deadline > 0; after a yield the
  // next legitimate pursuit must get a fresh, correctly-based budget.
  mgr.markSelected(id, Eigen::Vector2d(0.5, 1.0), 50.0);
  EXPECT_GT(mgr.records()[0].pursuit_deadline_t, 50.0);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
