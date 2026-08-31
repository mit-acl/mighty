// Tests for PeerTracker — the receiving side of /exploration/peer_poses, which
// feeds MinPos frontier allocation.
//
// This is the only genuinely concurrent class in multi-robot exploration: the
// peer-pose subscriber writes from a reentrant callback group while the
// explore-select timer reads. Both the staleness contract and the mutex are
// pinned here.

#include <gtest/gtest.h>

#include <atomic>
#include <thread>
#include <vector>

#include "mighty/peer_tracker.hpp"

TEST(PeerTrackerTest, UpdateThenGetActiveReturnsPeer) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 2.0, 10.0);

  auto peers = t.getActivePeers(10.5, 5.0);
  ASSERT_EQ(peers.size(), 1u);
  EXPECT_DOUBLE_EQ(peers[0].position.x(), 1.0);
  EXPECT_DOUBLE_EQ(peers[0].position.y(), 2.0);
}

TEST(PeerTrackerTest, StalePeerIsDropped) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 2.0, 10.0);
  EXPECT_TRUE(t.getActivePeers(20.0, 5.0).empty());
}

// The boundary is inclusive (`t_now - timestamp <= timeout_sec`). Pinned so a
// later refactor to a strict `<` is caught: at 5 Hz publishing and a 5 s
// timeout the exact-boundary sample is a real sample, not a corner case.
TEST(PeerTrackerTest, TimeoutBoundaryIsInclusive) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 2.0, 10.0);
  EXPECT_EQ(t.getActivePeers(15.0, 5.0).size(), 1u) << "exact timeout excluded";
  EXPECT_TRUE(t.getActivePeers(15.0001, 5.0).empty());
}

TEST(PeerTrackerTest, SameIdUpdatesInPlace) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 2.0, 10.0);
  t.updatePeer("NX02", 3.0, 4.0, 11.0);

  EXPECT_EQ(t.size(), 1u) << "same peer stored twice";
  auto peers = t.getActivePeers(11.0, 5.0);
  ASSERT_EQ(peers.size(), 1u);
  EXPECT_DOUBLE_EQ(peers[0].position.x(), 3.0);
}

TEST(PeerTrackerTest, MultiplePeersAllReturned) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 0.0, 10.0);
  t.updatePeer("NX03", 2.0, 0.0, 10.0);
  EXPECT_EQ(t.getActivePeers(10.1, 5.0).size(), 2u);
}

TEST(PeerTrackerTest, MixedFreshAndStaleReturnsOnlyFresh) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 0.0, 10.0);   // stale
  t.updatePeer("NX03", 2.0, 0.0, 20.0);   // fresh
  auto peers = t.getActivePeers(21.0, 5.0);
  ASSERT_EQ(peers.size(), 1u);
  EXPECT_DOUBLE_EQ(peers[0].position.x(), 2.0);
}

TEST(PeerTrackerTest, EmptyTrackerReturnsEmpty) {
  PeerTracker t;
  EXPECT_TRUE(t.getActivePeers(1.0, 5.0).empty());
  EXPECT_EQ(t.size(), 0u);
}

// Stale peers are filtered on read but never erased, so size() keeps counting
// them. Harmless for a fixed robot team (the map is keyed by namespace, so it
// is bounded by the fleet size) but worth stating explicitly, since size() and
// getActivePeers() disagreeing would otherwise look like a bug.
TEST(PeerTrackerTest, SizeCountsStalePeersToo) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 0.0, 10.0);
  ASSERT_TRUE(t.getActivePeers(100.0, 5.0).empty());
  EXPECT_EQ(t.size(), 1u);
}

// One writer, one reader, as in the node. Every observed position must be a
// coherent (x, y) pair — never x from one update and y from another.
TEST(PeerTrackerTest, ConcurrentUpdateAndReadStaysCoherent) {
  PeerTracker t;
  std::atomic<bool> stop{false};
  std::atomic<int> incoherent{0};

  std::thread writer([&] {
    for (int i = 0; i < 20000; ++i) {
      // Invariant the reader checks: y == 10 * x.
      const double x = static_cast<double>(i % 100);
      t.updatePeer("NX02", x, 10.0 * x, 1000.0);
    }
    stop = true;
  });

  std::thread reader([&] {
    while (!stop) {
      for (const auto& p : t.getActivePeers(1000.0, 5.0)) {
        if (p.position.y() != 10.0 * p.position.x()) ++incoherent;
      }
    }
  });

  writer.join();
  reader.join();
  EXPECT_EQ(incoherent.load(), 0) << "torn read — peer x and y came from different updates";
}

// ===========================================================================
// getActivePeerIds — the centralized-viz leader election's input
// ===========================================================================

TEST(PeerTrackerTest, ActivePeerIdsMatchActivePeers) {
  PeerTracker t;
  t.updatePeer("NX02", 1.0, 0.0, 10.0);   // stale at read time
  t.updatePeer("NX03", 2.0, 0.0, 20.0);   // fresh

  auto ids = t.getActivePeerIds(21.0, 5.0);
  ASSERT_EQ(ids.size(), 1u);
  EXPECT_EQ(ids[0], "NX03");
  EXPECT_EQ(ids.size(), t.getActivePeers(21.0, 5.0).size())
      << "the two active views disagree about who is alive";
}

// ===========================================================================
// ClaimTracker — the receiving side of /exploration/peer_claims
// ===========================================================================

TEST(ClaimTrackerTest, UpdateThenGetRetainsId) {
  ClaimTracker t;
  t.updateClaim("NX02", 1.0, 2.0, 10.0);

  auto claims = t.getActiveClaims(10.5, 5.0);
  ASSERT_EQ(claims.size(), 1u);
  EXPECT_EQ(claims[0].id, "NX02") << "tie-break needs to know WHO is claiming";
  EXPECT_DOUBLE_EQ(claims[0].position.x(), 1.0);
  EXPECT_DOUBLE_EQ(claims[0].position.y(), 2.0);
}

TEST(ClaimTrackerTest, StaleClaimDropped) {
  ClaimTracker t;
  t.updateClaim("NX02", 1.0, 2.0, 10.0);
  EXPECT_TRUE(t.getActiveClaims(20.0, 5.0).empty())
      << "a crashed or finished claimant must free its frontier via staleness";
}

// Same inclusive boundary contract as PeerTracker.
TEST(ClaimTrackerTest, TimeoutBoundaryIsInclusive) {
  ClaimTracker t;
  t.updateClaim("NX02", 1.0, 2.0, 10.0);
  EXPECT_EQ(t.getActiveClaims(15.0, 5.0).size(), 1u);
  EXPECT_TRUE(t.getActiveClaims(15.0001, 5.0).empty());
}

// Re-claim IS the release: there is no explicit un-claim message. A retarget
// overwrites the same key, so at every peer the claimant holds exactly one
// claim — the current one.
TEST(ClaimTrackerTest, SameIdUpsertsInPlace) {
  ClaimTracker t;
  t.updateClaim("NX02", 1.0, 2.0, 10.0);
  t.updateClaim("NX02", 5.0, 6.0, 11.0);

  EXPECT_EQ(t.size(), 1u);
  auto claims = t.getActiveClaims(11.0, 5.0);
  ASSERT_EQ(claims.size(), 1u);
  EXPECT_DOUBLE_EQ(claims[0].position.x(), 5.0);
}

// Same one-writer-one-reader shape as the node (subscriber on a reentrant
// group, selection tick reading).
TEST(ClaimTrackerTest, ConcurrentUpdateAndReadStaysCoherent) {
  ClaimTracker t;
  std::atomic<bool> stop{false};
  std::atomic<int> incoherent{0};

  std::thread writer([&] {
    for (int i = 0; i < 20000; ++i) {
      const double x = static_cast<double>(i % 100);
      t.updateClaim("NX02", x, 10.0 * x, 1000.0);
    }
    stop = true;
  });

  std::thread reader([&] {
    while (!stop) {
      for (const auto& c : t.getActiveClaims(1000.0, 5.0)) {
        if (c.position.y() != 10.0 * c.position.x()) ++incoherent;
      }
    }
  });

  writer.join();
  reader.join();
  EXPECT_EQ(incoherent.load(), 0);
}

// ===========================================================================
// winsClaimTie — the contested-claim resolution rule
// ===========================================================================

// `<` on strings is a strict total order: for distinct ids exactly one side
// wins, so both-yield (deadlock) and both-keep (duplicate) are impossible.
TEST(ClaimTieBreakTest, StrictTotalOrder) {
  EXPECT_TRUE(winsClaimTie("NX01", "NX02"));
  EXPECT_FALSE(winsClaimTie("NX02", "NX01"));
  EXPECT_NE(winsClaimTie("NX01", "NX02"), winsClaimTie("NX02", "NX01"));

  // Equal ids never contest in practice (self-filter at subscription), but
  // the function must not declare a win either way.
  EXPECT_FALSE(winsClaimTie("NX01", "NX01"));

  // Lexicographic, not numeric — fine for the NX0x fleet, and pinned so a
  // future 10+-robot naming scheme revisits this deliberately.
  EXPECT_TRUE(winsClaimTie("NX02", "NX10"));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
