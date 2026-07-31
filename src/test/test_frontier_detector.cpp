// Tests for FrontierDetector (Wavefront Frontier Detection).
//
// All tests build a synthetic OccGrid2D directly via fromTristate — no ROS, no
// rclcpp, no real planner.

#include <gtest/gtest.h>

#include <vector>

#include "mighty/frontier_detector.hpp"
#include "mighty/occ_grid_2d.hpp"

namespace {

constexpr int8_t U = -1;   // unknown
constexpr int8_t F = 0;    // free
constexpr int8_t O = 100;  // occupied

// Build a 10x10 OccGrid2D from a flat tristate buffer (row-major, iy=0 first).
std::shared_ptr<const OccGrid2D> MakeGrid(int w, int h,
                                          const std::vector<int8_t>& data,
                                          double res = 0.1) {
  return OccGrid2D::fromTristate(w, h, res, 0.0, 0.0, data);
}

// Detector with cluster_min_cells lowered for tiny test grids.
FrontierDetector MakeDetector(int min_cells = 2, int border_margin = 0) {
  FrontierDetectorParams p;
  p.cluster_min_cells   = min_cells;
  p.border_margin_cells = border_margin;
  p.robot_snap_radius_m = 5.0;  // generous snap so we never bail
  return FrontierDetector(p);
}

}  // namespace

TEST(FrontierDetectorTest, EmptyAllFreeGridYieldsNoFrontiers) {
  // 10x10, all free. No unknown anywhere -> no frontiers.
  const std::vector<int8_t> data(100, F);
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();
  auto out = det.detect(*grid, Eigen::Vector2d(0.5, 0.5));
  EXPECT_EQ(out.size(), 0u);
}

TEST(FrontierDetectorTest, HalfFreeHalfUnknownYieldsOneCluster) {
  // 10x10. Left 5 columns FREE, right 5 UNKNOWN.
  // Frontier seam should run along ix=4 (the rightmost free column).
  std::vector<int8_t> data(100);
  for (int iy = 0; iy < 10; ++iy) {
    for (int ix = 0; ix < 10; ++ix) {
      data[iy * 10 + ix] = (ix < 5) ? F : U;
    }
  }
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();
  auto out = det.detect(*grid, Eigen::Vector2d(0.15, 0.5));
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].size_cells, 10);  // 10 cells along iy=0..9 at ix=4
  // Centroid x should be near (4 + 0.5) * 0.1 = 0.45m
  EXPECT_NEAR(out[0].centroid.x(), 0.45, 1e-6);
}

TEST(FrontierDetectorTest, TwoDisconnectedFrontiersYieldTwoClusters) {
  // 11x11. Layout:
  //   columns 0..4 FREE, column 5 OCCUPIED (wall), columns 6..10 FREE
  //   row 0..3 FREE everywhere (left+right halves), row 4 OCCUPIED wall,
  //   row 5..10 alternating to make UNKNOWN above the wall.
  // Simpler: build two free pockets each surrounded by unknown on one side.
  //
  // 11x5 grid:
  //   row 0:  F F F F F O F F F F F
  //   row 1:  F F F F F O F F F F F
  //   row 2:  U U U U F O F U U U U
  //   row 3:  U U U U F O F U U U U
  //   row 4:  U U U U U U U U U U U
  //
  // Two free pockets connected only via row 0/1 around the wall — but the
  // wall blocks horizontal travel. Make the wall shorter so the BFS can
  // reach both halves: wall only on row 2 onwards.
  //
  // Actually, for two truly DISCONNECTED frontier clusters reachable from a
  // single robot pose, we need both pockets connected to the robot through
  // free space. Use a U-shaped corridor.
  //
  // 11x5 with two unknown blobs flanking a free corridor:
  //   row 0: F F F F F F F F F F F     <- free corridor (robot here)
  //   row 1: F F F F F F F F F F F
  //   row 2: U U U F F F F F U U U     <- frontier seams at iy=2 around blobs
  //   row 3: U U U F F F F F U U U
  //   row 4: U U U F F F F F U U U
  // Two unknown blobs: (0..2, 2..4) and (8..10, 2..4). Each adjacent to free
  // cells along iy=1 and iy=2 boundary. The frontier seeds are the free cells
  // adjacent to unknown.
  const int W = 11, H = 5;
  std::vector<int8_t> data(W * H, F);
  auto setCell = [&](int ix, int iy, int8_t v) { data[iy * W + ix] = v; };
  for (int iy = 2; iy <= 4; ++iy) {
    for (int ix = 0; ix <= 2; ++ix)  setCell(ix, iy, U);
    for (int ix = 8; ix <= 10; ++ix) setCell(ix, iy, U);
  }
  auto grid = MakeGrid(W, H, data);
  auto det = MakeDetector(/*min_cells=*/2, /*border_margin=*/0);
  auto out = det.detect(*grid, Eigen::Vector2d(0.55, 0.05));
  EXPECT_EQ(out.size(), 2u);
}

TEST(FrontierDetectorTest, MinClusterSizeFiltersTinyClusters) {
  // 10x10. Most free, but a single-cell unknown island. The frontier (the
  // free cell adjacent to that unknown) is one cluster of 1 cell — should be
  // filtered if cluster_min_cells > 1.
  std::vector<int8_t> data(100, F);
  data[5 * 10 + 5] = U;   // single unknown cell at (5,5)
  // Surround it with free; the 8 neighbors are frontier cells (one cluster).
  auto grid = MakeGrid(10, 10, data);

  // With min_cells=2, the cluster of 8 should pass.
  {
    auto det = MakeDetector(/*min_cells=*/2);
    auto out = det.detect(*grid, Eigen::Vector2d(0.05, 0.05));
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].size_cells, 8);
  }
  // With min_cells=20, the cluster should be filtered.
  {
    auto det = MakeDetector(/*min_cells=*/20);
    auto out = det.detect(*grid, Eigen::Vector2d(0.05, 0.05));
    EXPECT_EQ(out.size(), 0u);
  }
}

TEST(FrontierDetectorTest, BorderMarginExcludesEdgeFrontiers) {
  // 10x10. Inner 8x8 is free; outer ring is unknown (simulating fell-off
  // window). With border_margin=0, the inner-ring free cells along ix=1 or
  // iy=1 etc. become frontiers. With border_margin=2, all of those should
  // be excluded.
  std::vector<int8_t> data(100, U);
  for (int iy = 1; iy <= 8; ++iy)
    for (int ix = 1; ix <= 8; ++ix)
      data[iy * 10 + ix] = F;

  auto grid = MakeGrid(10, 10, data);

  // border_margin=0 -> seeds along inner ring -> non-empty.
  {
    FrontierDetectorParams p;
    p.cluster_min_cells   = 2;
    p.border_margin_cells = 0;
    p.robot_snap_radius_m = 5.0;
    FrontierDetector det(p);
    auto out = det.detect(*grid, Eigen::Vector2d(0.45, 0.45));
    EXPECT_GT(out.size(), 0u);
  }

  // border_margin=2 -> excludes the inner-ring; the only inner cells left are
  // (2..7, 2..7), none of which touch unknown (their neighbors are all free).
  // So 0 frontiers.
  {
    FrontierDetectorParams p;
    p.cluster_min_cells   = 2;
    p.border_margin_cells = 2;
    p.robot_snap_radius_m = 5.0;
    FrontierDetector det(p);
    auto out = det.detect(*grid, Eigen::Vector2d(0.45, 0.45));
    EXPECT_EQ(out.size(), 0u);
  }
}

TEST(FrontierDetectorTest, RobotOnOccupiedCellSnapsAndDetects) {
  // 10x10. Left 5 columns free, right 5 unknown, with the very center
  // occupied so the robot pose lands on it.
  std::vector<int8_t> data(100);
  for (int iy = 0; iy < 10; ++iy)
    for (int ix = 0; ix < 10; ++ix)
      data[iy * 10 + ix] = (ix < 5) ? F : U;
  data[5 * 10 + 2] = O;  // an obstacle in the free area near robot

  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();
  // Place robot directly on the occupied cell -> snap should find a free
  // neighbor and detection should still find the right-half frontier.
  auto out = det.detect(*grid, Eigen::Vector2d(0.25, 0.55));
  EXPECT_GT(out.size(), 0u);
}

TEST(FrontierDetectorTest, ObstacleClearanceFiltersWallHuggingFrontiers) {
  // 11x11 grid. Left half FREE, right half UNKNOWN. Plus a single OCCUPIED
  // cell at (5, 5) — right on the would-be frontier seam.
  //
  // With obstacle_clearance_cells=0, the frontier seam is the column ix=4
  // (10 cells along iy=0..9). Some of those cells have an occupied 8-neighbor.
  //
  // With obstacle_clearance_cells=1, every cell within 1 of (5,5) is
  // disqualified, so the frontier shrinks (or disappears for very small
  // sizes).
  const int W = 11, H = 11;
  std::vector<int8_t> data(W * H);
  for (int iy = 0; iy < H; ++iy)
    for (int ix = 0; ix < W; ++ix)
      data[iy * W + ix] = (ix < 5) ? F : U;
  data[5 * W + 5] = O;  // wall right on the seam at (5, 5)
  auto grid = MakeGrid(W, H, data);

  // No clearance -> the seam survives intact.
  {
    FrontierDetectorParams p;
    p.cluster_min_cells       = 2;
    p.border_margin_cells     = 0;
    p.obstacle_clearance_cells = 0;
    p.robot_snap_radius_m     = 5.0;
    FrontierDetector det(p);
    auto out = det.detect(*grid, Eigen::Vector2d(0.15, 0.55));
    ASSERT_EQ(out.size(), 1u);
    const int n_no_clearance = out[0].size_cells;
    EXPECT_GT(n_no_clearance, 0);
  }

  // Clearance=1 -> all seam cells within 1 of the wall are dropped.
  {
    FrontierDetectorParams p;
    p.cluster_min_cells       = 2;
    p.border_margin_cells     = 0;
    p.obstacle_clearance_cells = 1;
    p.robot_snap_radius_m     = 5.0;
    FrontierDetector det(p);
    auto out = det.detect(*grid, Eigen::Vector2d(0.15, 0.55));
    // The seam is 11 cells (iy=0..10 at ix=4); 3 of them (iy=4,5,6) are
    // within clearance=1 of the wall at (5,5) and get filtered. The
    // remainder is 8 cells, possibly split into 2 sub-clusters.
    int total_cells = 0;
    for (const auto& c : out) total_cells += c.size_cells;
    EXPECT_EQ(total_cells, 8);
    // No surviving cell may have an *in-bounds* occupied 8-neighbor.
    // (OOB doesn't count — the detector treats out-of-window as "outside",
    // not "occupied", consistent with the unknown-neighbor handling.)
    for (const auto& c : out) {
      for (const auto& cell : c.cells) {
        int ix, iy;
        grid->worldToGrid(cell.x(), cell.y(), ix, iy);
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            const int nx = ix + dx;
            const int ny = iy + dy;
            if (!grid->inBounds(nx, ny)) continue;
            EXPECT_FALSE(grid->isOccupied(nx, ny))
                << "frontier cell (" << ix << "," << iy
                << ") has occupied neighbor at (" << nx << "," << ny << ")";
          }
        }
      }
    }
  }
}

TEST(FrontierDetectorTest, FreeIslandDisconnectedFromUnknownYieldsNothing) {
  // 11x11. Center 3x3 is a free island, surrounded by occupied wall, with
  // unknown outside the wall. Robot is inside the island. The BFS can reach
  // none of the unknown cells (wall blocks it), so 0 frontiers.
  const int W = 11, H = 11;
  std::vector<int8_t> data(W * H, U);
  // wall ring at offset 4..6
  for (int iy = 3; iy <= 7; ++iy) {
    for (int ix = 3; ix <= 7; ++ix) {
      if (iy == 3 || iy == 7 || ix == 3 || ix == 7)
        data[iy * W + ix] = O;
      else
        data[iy * W + ix] = F;
    }
  }
  auto grid = MakeGrid(W, H, data);
  auto det = MakeDetector();
  auto out = det.detect(*grid, Eigen::Vector2d(0.55, 0.55));
  EXPECT_EQ(out.size(), 0u);
}

// ---- isStillFrontier: record re-validation against the current grid --------

TEST(FrontierDetectorTest, IsStillFrontierTrueOnLiveFrontier) {
  // Left 5 cols FREE, right 5 UNKNOWN -> a 10-cell frontier seam at ix=4.
  std::vector<int8_t> data(100);
  for (int iy = 0; iy < 10; ++iy)
    for (int ix = 0; ix < 10; ++ix)
      data[iy * 10 + ix] = (ix < 5) ? F : U;
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();

  double wx, wy;
  grid->gridToWorld(4, 5, wx, wy);  // right on the seam
  EXPECT_TRUE(det.isStillFrontier(*grid, Eigen::Vector2d(wx, wy),
                                  /*min_cells=*/5, /*search_radius_cells=*/2));
}

TEST(FrontierDetectorTest, IsStillFrontierFalseAfterAreaExplored) {
  // The former seam location is now fully observed FREE -> no frontier here,
  // even though it is a known/free cell (the exact phantom we want to retire).
  const std::vector<int8_t> data(100, F);
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();

  double wx, wy;
  grid->gridToWorld(4, 5, wx, wy);
  EXPECT_FALSE(det.isStillFrontier(*grid, Eigen::Vector2d(wx, wy),
                                   /*min_cells=*/5, /*search_radius_cells=*/2));
}

TEST(FrontierDetectorTest, IsStillFrontierFalseForSubThresholdRemnant) {
  // A lone UNKNOWN speck leaves an 8-cell frontier ring around it. With
  // min_cells=15 that is below threshold, so the record must be retired -- this
  // is the detector/manager asymmetry the reconcile pass exists to close.
  std::vector<int8_t> data(100, F);
  data[5 * 10 + 5] = U;  // single unexplored speck (e.g. an occlusion shadow)
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();

  double wx, wy;
  grid->gridToWorld(4, 5, wx, wy);
  EXPECT_FALSE(det.isStillFrontier(*grid, Eigen::Vector2d(wx, wy),
                                   /*min_cells=*/15, /*search_radius_cells=*/2));
}

TEST(FrontierDetectorTest, IsStillFrontierToleratesCentroidDrift) {
  // Live seam at ix=4, but query a centroid nudged 2 cells into free space
  // (EMA drift). A radius-1 search can't reach the seam; radius-2 can.
  std::vector<int8_t> data(100);
  for (int iy = 0; iy < 10; ++iy)
    for (int ix = 0; ix < 10; ++ix)
      data[iy * 10 + ix] = (ix < 5) ? F : U;
  auto grid = MakeGrid(10, 10, data);
  auto det = MakeDetector();

  double wx, wy;
  grid->gridToWorld(2, 5, wx, wy);  // 2 cells left of the seam
  EXPECT_FALSE(det.isStillFrontier(*grid, Eigen::Vector2d(wx, wy),
                                   /*min_cells=*/5, /*search_radius_cells=*/1));
  EXPECT_TRUE(det.isStillFrontier(*grid, Eigen::Vector2d(wx, wy),
                                  /*min_cells=*/5, /*search_radius_cells=*/2));
}

// --- unknown_bridge_cells: sensor near-field seam ---------------------------
//
// Reproduces the startup deadlock seen on the ground robot: the mapper clears a
// small cylinder around the robot, but the first ray-cleared cells sit farther
// out, so the robot's free pocket is ringed by UNKNOWN and not 8-connected to
// the free space it can actually see. A free-only BFS reaches the pocket alone,
// the only cluster it reports sits on top of the robot, the dwell check retires
// it VISITED within a second, and the robot never gets a goal.
namespace {

// 20x10: robot pocket at ix 1..3, a `seam` wide UNKNOWN gap, then a large free
// region at ix 6..14. Everything else UNKNOWN.
std::shared_ptr<const OccGrid2D> MakeSeamGrid() {
  std::vector<int8_t> data(200, U);
  for (int iy = 3; iy <= 6; ++iy)
    for (int ix = 1; ix <= 3; ++ix) data[iy * 20 + ix] = F;   // robot pocket
  for (int iy = 1; iy <= 8; ++iy)
    for (int ix = 6; ix <= 14; ++ix) data[iy * 20 + ix] = F;  // main region
  return OccGrid2D::fromTristate(20, 10, 0.1, 0.0, 0.0, data);
}

FrontierDetector MakeBridgeDetector(int bridge) {
  FrontierDetectorParams p;
  p.cluster_min_cells    = 2;
  p.border_margin_cells  = 0;
  p.robot_snap_radius_m  = 5.0;
  p.unknown_bridge_cells = bridge;
  return FrontierDetector(p);
}

}  // namespace

TEST(FrontierDetectorTest, NearFieldSeamStrandsRobotWithoutBridging) {
  auto grid = MakeSeamGrid();
  auto det = MakeBridgeDetector(0);
  // Robot sits in the pocket at (ix=2, iy=4).
  auto out = det.detect(*grid, Eigen::Vector2d(0.25, 0.45));
  // Only the pocket is reachable: one cluster, and it is centred on the robot.
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].size_cells, 10);  // 3x4 pocket minus its 2 interior cells
  EXPECT_LT(out[0].centroid.x(), 0.4);
}

TEST(FrontierDetectorTest, BridgingMustCoverTheFullSeamWidth) {
  auto grid = MakeSeamGrid();
  // Seam is 2 cells wide (ix=4,5); a 1-cell budget cannot cross it.
  auto out = MakeBridgeDetector(1).detect(*grid, Eigen::Vector2d(0.25, 0.45));
  ASSERT_EQ(out.size(), 1u);
  EXPECT_LT(out[0].centroid.x(), 0.4);
}

TEST(FrontierDetectorTest, BridgingReconnectsRobotToMainFreeRegion) {
  auto grid = MakeSeamGrid();
  auto out = MakeBridgeDetector(2).detect(*grid, Eigen::Vector2d(0.25, 0.45));
  // Both rings are now found, and they stay distinct clusters (Pass B only
  // links frontier cells, which are free by definition). Sorted size-descending.
  ASSERT_EQ(out.size(), 2u);
  EXPECT_EQ(out[0].size_cells, 30);  // boundary ring of the 9x8 main region
  EXPECT_GT(out[0].centroid.x(), 0.6);
  EXPECT_EQ(out[1].size_cells, 10);  // the pocket
}

TEST(FrontierDetectorTest, BridgingDoesNotLeakThroughOccupiedCells) {
  // Same layout, but the seam is a wall instead of unobserved space. Even a
  // generous budget must not cross it.
  std::vector<int8_t> data(200, U);
  for (int iy = 3; iy <= 6; ++iy)
    for (int ix = 1; ix <= 3; ++ix) data[iy * 20 + ix] = F;
  for (int iy = 0; iy < 10; ++iy)
    for (int ix = 4; ix <= 5; ++ix) data[iy * 20 + ix] = O;   // wall
  for (int iy = 1; iy <= 8; ++iy)
    for (int ix = 6; ix <= 14; ++ix) data[iy * 20 + ix] = F;
  auto grid = OccGrid2D::fromTristate(20, 10, 0.1, 0.0, 0.0, data);

  auto out = MakeBridgeDetector(8).detect(*grid, Eigen::Vector2d(0.25, 0.45));
  // Pocket only. obstacle_clearance_cells=1 drops the column hugging the wall,
  // so the reported cluster is smaller than the unobstructed case.
  ASSERT_EQ(out.size(), 1u);
  EXPECT_LT(out[0].centroid.x(), 0.4);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
