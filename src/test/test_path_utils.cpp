// Tests for mighty_utils::truncateGlobalPathAtHorizon (header-only, no ROS).

#include <gtest/gtest.h>

#include <vector>

#include "mighty/path_utils.hpp"

namespace {

vec_Vecf<3> MakePath(const std::vector<Eigen::Vector3d>& pts) {
  vec_Vecf<3> p;
  for (const auto& v : pts) p.push_back(v);
  return p;
}

}  // namespace

TEST(PathUtilsTest, GoalBeyondHorizonTruncatesAtBoundary) {
  // Straight path 0 -> 10 along +x; horizon 3 -> subgoal at (3,0,0).
  auto path = MakePath({{0, 0, 0}, {10, 0, 0}});
  Eigen::Vector3d subgoal;
  ASSERT_TRUE(mighty_utils::truncateGlobalPathAtHorizon(path, Eigen::Vector3d(0, 0, 0), 3.0, subgoal));
  EXPECT_NEAR(subgoal.x(), 3.0, 1e-6);
  EXPECT_NEAR(subgoal.y(), 0.0, 1e-6);
  // Path now ends exactly at the subgoal.
  EXPECT_NEAR((path.back() - subgoal).norm(), 0.0, 1e-6);
  EXPECT_NEAR(path.back().x(), 3.0, 1e-6);
}

TEST(PathUtilsTest, GoalWithinHorizonKeepsPath) {
  // Whole path within horizon -> unchanged, subgoal == goal.
  auto path = MakePath({{0, 0, 0}, {2, 0, 0}});
  Eigen::Vector3d subgoal;
  ASSERT_TRUE(mighty_utils::truncateGlobalPathAtHorizon(path, Eigen::Vector3d(0, 0, 0), 5.0, subgoal));
  EXPECT_NEAR(subgoal.x(), 2.0, 1e-6);
  EXPECT_EQ(path.size(), 2u);
}

TEST(PathUtilsTest, ShortPathReturnsFalse) {
  auto path = MakePath({{0, 0, 0}});
  Eigen::Vector3d subgoal;
  EXPECT_FALSE(mighty_utils::truncateGlobalPathAtHorizon(path, Eigen::Vector3d(0, 0, 0), 3.0, subgoal));
}

TEST(PathUtilsTest, CrossingInLaterSegment) {
  // 0 -> 2 -> 6 along +x; horizon 3 crosses in the 2nd segment at x=3.
  auto path = MakePath({{0, 0, 0}, {2, 0, 0}, {6, 0, 0}});
  Eigen::Vector3d subgoal;
  ASSERT_TRUE(mighty_utils::truncateGlobalPathAtHorizon(path, Eigen::Vector3d(0, 0, 0), 3.0, subgoal));
  EXPECT_NEAR(subgoal.x(), 3.0, 1e-6);
  EXPECT_EQ(path.size(), 3u);  // kept [0, 2] + subgoal(3)
}

TEST(PathUtilsTest, DiagonalCrossingInterpolatesCorrectly) {
  // Diagonal 0 -> (10,10,0); horizon 5 -> crossing at distance 5 from origin,
  // i.e. (5/sqrt2, 5/sqrt2) = (3.5355, 3.5355).
  auto path = MakePath({{0, 0, 0}, {10, 10, 0}});
  Eigen::Vector3d subgoal;
  ASSERT_TRUE(mighty_utils::truncateGlobalPathAtHorizon(path, Eigen::Vector3d(0, 0, 0), 5.0, subgoal));
  EXPECT_NEAR(subgoal.norm(), 5.0, 1e-6);
  EXPECT_NEAR(subgoal.x(), 5.0 / std::sqrt(2.0), 1e-6);
  EXPECT_NEAR(subgoal.y(), 5.0 / std::sqrt(2.0), 1e-6);
}

// ---- reanchorPathToStart -----------------------------------------------------

TEST(PathUtilsTest, ReanchorDropsPassedVerticesAndPrependsStart) {
  // Route 0->2->4->6 along +x; robot now at x=3 (nearest vertex is x=2 or x=4,
  // both dist 1 -> first found (x=2, index 1) wins by strict-less-than).
  auto path = MakePath({{0, 0, 0}, {2, 0, 0}, {4, 0, 0}, {6, 0, 0}});
  double deviation = -1.0;
  ASSERT_TRUE(mighty_utils::reanchorPathToStart(path, Eigen::Vector3d(3, 0, 0), deviation));
  EXPECT_NEAR(deviation, 1.0, 1e-6);
  // Starts at the robot; keeps from the nearest vertex onward.
  EXPECT_NEAR(path.front().x(), 3.0, 1e-6);
  EXPECT_NEAR(path.back().x(), 6.0, 1e-6);
  EXPECT_EQ(path.size(), 4u);  // [3, 2, 4, 6]
}

TEST(PathUtilsTest, ReanchorDeviationReflectsLateralOffset) {
  // Route along +x; robot offset laterally by 0.5 in y near x=2.
  auto path = MakePath({{0, 0, 0}, {2, 0, 0}, {4, 0, 0}});
  double deviation = -1.0;
  ASSERT_TRUE(mighty_utils::reanchorPathToStart(path, Eigen::Vector3d(2, 0.5, 0), deviation));
  EXPECT_NEAR(deviation, 0.5, 1e-6);
  EXPECT_NEAR(path.front().y(), 0.5, 1e-6);
}

TEST(PathUtilsTest, ReanchorEmptyReturnsFalse) {
  vec_Vecf<3> path;
  double deviation = -1.0;
  EXPECT_FALSE(mighty_utils::reanchorPathToStart(path, Eigen::Vector3d(0, 0, 0), deviation));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
