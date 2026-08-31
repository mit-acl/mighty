// Tests for VisitedMap — the persistent per-robot map, and specifically
// mergeFrom(), which is how one robot learns from another's exploration.
//
// mergeFrom is the receiving half of /exploration/visited_maps. Until
// multi-robot exploration was enabled it had never executed: the sim config
// left exploration.minpos.enabled false, so the publisher and subscriber that
// drive it were never even constructed. These tests pin its contract before it
// carries real data.

#include <gtest/gtest.h>

#include <vector>

#include "mighty/occ_grid_2d.hpp"
#include "mighty/visited_map.hpp"

namespace {

constexpr int8_t U = VisitedMap::kUnknown;
constexpr int8_t F = VisitedMap::kFree;
constexpr int8_t O = VisitedMap::kOccupied;

// 2x2 m map centred on the origin at 0.1 m -> 20x20 cells spanning [-1, 1]^2.
VisitedMap MakeMap() { return VisitedMap(0.0, 0.0, 2.0, 2.0, 0.1); }

// A peer buffer of w*h cells, all `fill`, whose lower-left corner is at
// (origin_x, origin_y).
std::vector<int8_t> MakePeerBuffer(int w, int h, int8_t fill) {
  return std::vector<int8_t>(static_cast<size_t>(w) * h, fill);
}

}  // namespace

TEST(VisitedMapTest, MergeFillsUnknownCells) {
  VisitedMap m = MakeMap();
  ASSERT_EQ(m.getStateWorld(0.05, 0.05), U);

  auto peer = MakePeerBuffer(4, 4, F);
  m.mergeFrom(peer.data(), 4, 4, 0.0, 0.0, 0.1);  // covers [0, 0.4]^2

  EXPECT_EQ(m.getStateWorld(0.05, 0.05), F);
  EXPECT_EQ(m.getStateWorld(0.35, 0.35), F);
  EXPECT_EQ(m.getStateWorld(0.55, 0.55), U) << "merged outside the peer's extent";
}

// The safety property: a peer must never be able to erase what this robot saw
// with its own sensor — in either direction. Peer-reported free space must not
// delete my wall, and peer-reported wall must not delete my free space.
TEST(VisitedMapTest, MergeNeverOverwritesOwnObservations) {
  VisitedMap m = MakeMap();
  m.setStateWorld(0.05, 0.05, O);   // I saw a wall here
  m.setStateWorld(0.15, 0.05, F);   // I saw free space here

  auto peer_free = MakePeerBuffer(4, 4, F);
  m.mergeFrom(peer_free.data(), 4, 4, 0.0, 0.0, 0.1);
  EXPECT_EQ(m.getStateWorld(0.05, 0.05), O) << "peer erased my wall";

  auto peer_occ = MakePeerBuffer(4, 4, O);
  m.mergeFrom(peer_occ.data(), 4, 4, 0.0, 0.0, 0.1);
  EXPECT_EQ(m.getStateWorld(0.15, 0.05), F) << "peer erased my free space";
}

TEST(VisitedMapTest, MergeSkipsPeerUnknownCells) {
  VisitedMap m = MakeMap();
  auto peer = MakePeerBuffer(4, 4, U);
  m.mergeFrom(peer.data(), 4, 4, 0.0, 0.0, 0.1);
  EXPECT_EQ(m.getStateWorld(0.05, 0.05), U);
}

// Peers do not share an origin, so the merge has to go through world
// coordinates rather than cell indices. Offset by a non-integer number of
// cells so an index-space shortcut would land in the wrong place.
TEST(VisitedMapTest, MergeWithDifferentOriginUsesWorldCoords) {
  VisitedMap m = MakeMap();
  auto peer = MakePeerBuffer(2, 2, O);
  // Peer cell (0,0) centre = -0.53 + 0.05 = -0.48 -> local cell containing -0.48.
  m.mergeFrom(peer.data(), 2, 2, -0.53, -0.53, 0.1);

  EXPECT_EQ(m.getStateWorld(-0.48, -0.48), O);
  EXPECT_EQ(m.getStateWorld(0.5, 0.5), U) << "merge landed far from the peer extent";
}

// Documents a real limitation rather than asserting a fix. mergeFrom iterates
// the PEER grid and writes one local cell per peer cell, so a coarser peer map
// leaves holes: a 0.2 m peer cell touches only one of the four 0.1 m local
// cells it geometrically covers. Harmless today because every robot runs the
// same config and therefore the same resolution — but the signature accepts a
// resolution argument, so the behaviour should be discovered here rather than
// in the field.
TEST(VisitedMapTest, MergeFromCoarserPeerLeavesGaps) {
  VisitedMap m = MakeMap();
  auto peer = MakePeerBuffer(2, 2, F);
  m.mergeFrom(peer.data(), 2, 2, 0.0, 0.0, 0.2);  // 0.2 m cells over [0, 0.4]^2

  // Peer cell (0,0) centre is (0.1, 0.1) -> exactly one local cell is written.
  EXPECT_EQ(m.getStateWorld(0.1, 0.1), F);
  // The rest of the area that peer cell covers stays unknown.
  EXPECT_EQ(m.getStateWorld(0.05, 0.05), U)
      << "coarse-peer merge unexpectedly filled the whole covered area";
}

TEST(VisitedMapTest, MergeOutOfBoundsIsDroppedWithoutCrashing) {
  VisitedMap m = MakeMap();
  auto peer = MakePeerBuffer(4, 4, F);
  m.mergeFrom(peer.data(), 4, 4, 50.0, 50.0, 0.1);   // entirely outside
  for (int8_t v : m.data()) EXPECT_EQ(v, U);
}

TEST(VisitedMapTest, MergeIntoDefaultConstructedMapIsNoOp) {
  VisitedMap m;   // empty
  ASSERT_TRUE(m.empty());
  auto peer = MakePeerBuffer(4, 4, F);
  m.mergeFrom(peer.data(), 4, 4, 0.0, 0.0, 0.1);
  EXPECT_TRUE(m.empty());
  EXPECT_EQ(m.getStateWorld(0.05, 0.05), U);
}

TEST(VisitedMapTest, AbsorbRecordsFreeAndOccupiedAndSkipsUnknown) {
  VisitedMap m = MakeMap();
  // 3x1 local grid over [0, 0.3] x [0, 0.1]: free, occupied, unknown.
  std::vector<int8_t> cells = {0, 100, -1};
  auto grid = OccGrid2D::fromTristate(3, 1, 0.1, 0.0, 0.0, cells);
  m.absorb(*grid);

  EXPECT_EQ(m.getStateWorld(0.05, 0.05), F);
  EXPECT_EQ(m.getStateWorld(0.15, 0.05), O);
  EXPECT_EQ(m.getStateWorld(0.25, 0.05), U) << "absorb recorded an unobserved cell";
}

// End-to-end precedence: what I absorbed from my own sensor survives a
// subsequent peer merge over the same area.
TEST(VisitedMapTest, AbsorbThenMergeKeepsOwnObservations) {
  VisitedMap m = MakeMap();
  std::vector<int8_t> cells = {100, 100};
  auto grid = OccGrid2D::fromTristate(2, 1, 0.1, 0.0, 0.0, cells);
  m.absorb(*grid);
  ASSERT_EQ(m.getStateWorld(0.05, 0.05), O);

  auto peer = MakePeerBuffer(4, 4, F);
  m.mergeFrom(peer.data(), 4, 4, 0.0, 0.0, 0.1);

  EXPECT_EQ(m.getStateWorld(0.05, 0.05), O);
  EXPECT_EQ(m.getStateWorld(0.15, 0.05), O);
  EXPECT_EQ(m.getStateWorld(0.25, 0.05), F) << "peer data should fill what I never saw";
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
