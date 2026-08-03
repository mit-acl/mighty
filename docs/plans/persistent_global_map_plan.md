# Implementation Plan: Persistent Global Map + True-Goal Planning (Ground Robot)

**Goal:** Make the ground-robot global A\* planner plan over a **persistent 2D map**
(built from the existing `visited_map_`) toward the **true terminal goal**
(`G_term`), then derive the local subgoal by **truncating** the resulting global
path at arc-length = `horizon`. This eliminates the two per-tick variation sources
that cause the left/right obstacle wobble: (1) the horizon-**projected** goal that
moves as the robot moves, and (2) the **sliding** local map that re-hides cells.

**Scope:** Ground robot 2D only (`vehicle_type == "ground_robot"`), gated behind a
new `use_persistent_global_map` flag. UAV behavior is untouched.

**Architecture parallel (Nav2):** global planner plans to the true goal over a
persistent global costmap; the local controller tracks a windowed slice. Here the
"window" becomes a horizon-truncation of a *stable* global route instead of a
straight-line goal projection.

---

## Phase 0 — Allowed APIs (verified against source; do NOT invent)

All facts below were read directly from the repo. Cite these; do not assume other
signatures exist.

### Map types
- `OccGrid2D` is **immutable**, shared as `std::shared_ptr<const OccGrid2D>`, and
  exposes **bit-vectors, not an int8 buffer**.
  - `static std::shared_ptr<const OccGrid2D> OccGrid2D::fromTristate(int width, int height, double resolution, double origin_x, double origin_y, const std::vector<int8_t>& data)` — `include/mighty/occ_grid_2d.hpp:51-70`. Encoding: `occupied = (data[i] >= 100)`, `unknown = (data[i] < 0)`.
  - `static ... fromOccupancyGrid(const nav_msgs::msg::OccupancyGrid&)` — `occ_grid_2d.hpp:30-48`.
  - Accessors: `width() :167`, `height() :168`, `originX() :169`, `originY() :170`, `resolution() :171`, `occupiedData() :173`, `unknownData() :174`, `isOccupied/isFree/isUnknown/inBounds/worldToGrid/gridToWorld`.
  - **No `data()` accessor** on OccGrid2D (that method belongs to `VisitedMap`).
- `VisitedMap` (`include/mighty/visited_map.hpp`): `absorb(const OccGrid2D&) :74`, `data() -> const std::vector<int8_t>& :147`, `width() :144`, `height() :145`, `resolution() :143`, `originX() :141`, `originY() :142`, `empty() :146`. Constants `kUnknown=-1, kFree=0, kOccupied=100 :45-47`. `absorb()` skips unknown local cells and overwrites observed cells with their most-recent FREE/OCCUPIED state (no raytrace clearing beyond re-observation).

### Map ingestion into the planner (already wired; reuse unchanged)
- `MIGHTY::setOccGrid2D(std::shared_ptr<const OccGrid2D>)` — `include/mighty/mighty.hpp:398` (stores into the **core's** `occ_grid_2d_` member `mighty.hpp:459`, separate from the node's member).
- Node passes the grid to the core at `mighty_node.cpp:3257`: `mighty_ptr_->setOccGrid2D(occ_grid_2d_);`
- Core → manager at `mighty.cpp:856`: `hgp_manager_.setOccGrid2D(occ_grid_2d_);` (during `setupHGPPlanner`).
- Manager → A* map at `hgp_manager.cpp:1214-1216`: `map_util_->buildMap2DFromOcc2D(*occ_grid_2d_, ...)`.
- **Unknown is traversable:** `buildMap2DFromOcc2D` (`map_util.hpp:2006-2075`) starts every cell FREE and only marks a cell occupied if a source cell is `>=100`; unknown source cells stay free. `get2DOccupancy` returns `0` for free+unknown, `100` for occupied and out-of-bounds. `w_unknown_` in `graph_search.cpp:676` keys off the **3D** `cMap_`, not the 2D map — so swapping the 2D source preserves plan-through-unknown behavior exactly.

### Copy-ready persistent-grid build (mirror this exactly)
`mighty_node.cpp:3282-3286` (currently builds the exploration detect grid):
```cpp
current_detect_grid_ = OccGrid2D::fromTristate(
    visited_map_->width(), visited_map_->height(),
    visited_map_->resolution(),
    visited_map_->originX(), visited_map_->originY(),
    visited_map_->data());
```
`absorb` call: `mighty_node.cpp:3274` — `if (visited_map_) visited_map_->absorb(*occ_grid_2d_);` (runs inside the `expl_enabled && state_initialized_` block starting `:3263`).

### Global path generation & goal projection
- `MIGHTY::computeG(const state& A, const state& G_term, double horizon)` — `mighty.cpp:135-148`. Sets `local_G.pos = mighty_utils::projectPointToSphere(A.pos, G_term.pos, horizon)` and `setG(local_G)`. **This is the projection to replace.**
- `MIGHTY::generateGlobalPath(...)` — `mighty.cpp:807-956`. Key lines:
  - `:841` `computeG(local_A, local_G_term, par_.horizon);` (writes projected subgoal into `G_`)
  - `:850` `getG(local_G);`
  - `:908-910` `hgp_manager_.solveHGP(local_A.pos, start_dir_hint, local_G.pos, final_g_, ...)` — **goal arg is `local_G.pos`; change to the true goal.**
  - `:924-930` corridor-hop block: `computeG_corridorHop(local_A, local_G_term, global_path, raw_global_path); getG(local_G);`
  - `:933-936` store `global_path_` under `mtx_global_path_`.
- `HGPManager::solveHGP(const Vec3f& start_sent, const Vec3f& start_vel, const Vec3f& goal_sent, double& final_g, double weight, double current_time, vec_Vecf<3>& path, vec_Vecf<3>& raw_path)` — `hgp_manager.cpp:259` / `hgp_manager.hpp:155`.
- Goal accessors/members: `getGterm :118`, `setGterm :123`, `getG :128`, `setG :138`, `getA :143`, `setA :148`, `getGlobalPath :289`; members `G_ :507`, `G_term_ :511`, `global_path_ :581`.

### Path truncation utilities (reuse; do NOT reimplement geometry)
- `mighty_utils::projectPointToSphere(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, double radius)` — `src/mighty/utils.cpp:567` / `include/mighty/utils.hpp:193`. Straight-line only (does not follow a path).
- `getLastIntersectionWithSphere(vec_Vecf<3> path, double r, Eigen::Vector3d center)` — `src/hgp/utils.cpp:697-716` / `include/hgp/utils.hpp:269`. Path-aware; scans from path END backward to first point inside radius `r`, interpolates. **Caveats:** takes `path` by value; assumes `path[0]` inside the sphere; returns `path.back()` if the whole path is inside; a second overload with `double* Jdist` (`utils.cpp:728-757`) **prints an ERROR** when the whole path is inside — prefer the no-`Jdist` overload or guard.
- `resamplePathUniform(vec_Vecf<3>& path, double spacing)` — `utils.cpp:686` / `utils.hpp:221`. In-place uniform arc-length resample.

### Copy-ready truncation pattern (mirror this for horizon truncation)
`MIGHTY::computeG_corridorHop` found-branch, `mighty.cpp:305-337`:
```cpp
const Vecf<3>& p_corner = raw[pick_idx];                 // :306
local_G.pos = ...backed-off point...;                    // :307-309
local_G.yaw = std::atan2(outgoing_unit.y(), outgoing_unit.x());  // :310
global_path.clear();                                     // :315
global_path.reserve(pick_idx + 1);                       // :316
for (size_t k = 0; k <= pick_idx; ++k) global_path.push_back(raw[k]); // :317-319
if (!global_path.empty()) global_path.back() = local_G.pos;          // :320-322
setG(local_G);                                           // :337
```

### Optimizer endpoint (why truncation is sufficient)
`generateLocalTrajectory` sets the L-BFGS endpoint from the **global path end**, not `local_G`:
- `mighty.cpp:1078-1082`: `if (GOAL_REACHED || GOAL_SEEN) local_E = local_G; else local_E.pos = global_path.back();`
- `mighty.cpp:1105`: `local_E.pos = global_path.back();`
So once `global_path` is truncated at the horizon, the optimizer re-targets automatically. `setG(...)` must still be updated (used for yaw, the GOAL_SEEN override, and corridor-hop arrival detection in `needReplan` `:375-387`).

### Third projection site (keep consistent)
`MIGHTY::setTerminalGoal` — `mighty.cpp:1789-1790` seeds `G_.pos = projectPointToSphere(local_state.pos, term_goal.pos, par_.horizon)` under `mtx_G_` when a new goal arrives.

### Param / config / test conventions
- **Param 3-site pattern:** declare at `mighty_node.cpp:732` (`this->declare_parameter("corridor_hop_enabled", false);`); read at `:1096` (`par_.corridor_hop_enabled = this->get_parameter(...).as_bool();`); struct field at `mighty_type.hpp:295` (`bool corridor_hop_enabled{false};`). Double example: declare `:733`, read `:1097`, field `mighty_type.hpp:296`.
- **Gating idiom:** `if (par_.vehicle_type == "ground_robot" && par_.corridor_hop_enabled)` — `mighty_node.cpp:2394`.
- **Config:** `config/mighty_ground_robot.yaml`; `global_planner` at line 89, `horizon` at 95, corridor block 286-304.
- **Struct insert point:** `mighty_type.hpp` after line 304 (end of corridor block) or near `horizon :192`.
- **Test registration:** `CMakeLists.txt:187-198` (`ament_add_gtest(<name> <test.cpp> <prod.cpp...>)` + `ament_target_dependencies(<name> rclcpp Eigen3 nav_msgs)`); insert before `endif()` at `:217`. A helper test that uses `mighty_utils` must add `src/mighty/utils.cpp` to the gtest sources (no existing test links it). Closest pure-helper test to copy: `src/test/test_occ_grid_2d_unknown.cpp`.

### Anti-patterns to avoid
- Do **not** call `OccGrid2D::data()` — it does not exist; use `visited_map_->data()` fed into `OccGrid2D::fromTristate`.
- Do **not** try to mutate an `OccGrid2D` in place — rebuild via `fromTristate`.
- Do **not** clobber the **node's** `occ_grid_2d_` member (exploration/frontier/absorb use it) — feed the persistent grid to the **core** via `mighty_ptr_->setOccGrid2D(...)`.
- Do **not** reimplement path/sphere geometry — use `getLastIntersectionWithSphere` / `projectPointToSphere`.
- Do **not** use the `getLastIntersectionWithSphere` overload with `Jdist` in the hot loop (it prints ERROR when the whole path is inside the sphere).
- Do **not** touch the UAV path — gate every change on `vehicle_type == "ground_robot" && use_persistent_global_map`.

---

## Phase 1 — Plan over the persistent map (goal unchanged)

**Objective:** Feed the ground-robot global planner an `OccGrid2D` built from
`visited_map_` instead of the sliding `occ_2d`. Keep planning to the projected
`local_G` for now, to isolate "does planning over the persistent map work?".

### What to implement
1. **Add param `use_persistent_global_map` (bool, default `false`).**
   - Declare: copy `mighty_node.cpp:732` idiom → `this->declare_parameter("use_persistent_global_map", false);`
   - Read: copy `mighty_node.cpp:1096` idiom → `par_.use_persistent_global_map = this->get_parameter("use_persistent_global_map").as_bool();`
   - Struct field: copy `mighty_type.hpp:295` idiom → add `bool use_persistent_global_map{false};` after line 304.
   - Config: add to `config/mighty_ground_robot.yaml` near `global_planner` (line 89) — `use_persistent_global_map: true` with a comment.

2. **In `occ2DCallback` (`mighty_node.cpp`), feed the core a persistent grid.**
   - The node's `occ_grid_2d_` member stays the sliding grid (exploration needs it).
   - Ensure `absorb` runs **before** building the planning grid. Currently `absorb` is at `:3274` inside the `expl_enabled` block; the persistent-planning path must not depend on `expl_enabled`. Restructure so that when `use_persistent_global_map && vehicle_type=="ground_robot"`:
     - call `visited_map_->absorb(*occ_grid_2d_)` (guard `visited_map_ && !visited_map_->empty()` after first absorb),
     - build `auto planning_grid = OccGrid2D::fromTristate(visited_map_->width(), visited_map_->height(), visited_map_->resolution(), visited_map_->originX(), visited_map_->originY(), visited_map_->data());` (copy from `mighty_node.cpp:3282-3286`),
     - call `mighty_ptr_->setOccGrid2D(planning_grid);`
   - Else (current behavior): `mighty_ptr_->setOccGrid2D(occ_grid_2d_);` (the existing `:3257` line).
   - **Note the ordering:** if exploration is enabled, `absorb` also runs at `:3274`; make the planning-grid build reuse the same absorbed map (build it after `:3274`, or hoist absorb up and skip the duplicate). Do not absorb twice per callback.

### Documentation references
- Persistent grid build: `mighty_node.cpp:3282-3286` (copy verbatim).
- `absorb`: `mighty_node.cpp:3274`; `VisitedMap::absorb` `visited_map.hpp:74`.
- setOccGrid2D chain: `mighty.hpp:398`, `mighty.cpp:856`, `hgp_manager.cpp:1214-1216`.
- Param pattern: `mighty_node.cpp:732/1096`, `mighty_type.hpp:295`.
- Gating idiom: `mighty_node.cpp:2394`.

### Verification checklist
- [ ] `colcon build --packages-select mighty --cmake-args -DCMAKE_BUILD_TYPE=Release` succeeds.
- [ ] `grep -n "use_persistent_global_map" src/mighty/mighty_node.cpp include/mighty/mighty_type.hpp config/mighty_ground_robot.yaml` shows all 4 sites.
- [ ] With the flag **off**, behavior is byte-identical to today (the `else` branch calls the original `:3257` line).
- [ ] With the flag **on**, run the ground-robot exploration sim; in RViz confirm the A\* path is planned over the accumulated persistent extent (not just the ~horizon window) — e.g. the robot re-enters a previously mapped area and the plan uses old geometry.
- [ ] No double-absorb per callback (add a temporary log or assert if unsure).

### Anti-pattern guards
- Do not replace the node's `occ_grid_2d_` member.
- Do not gate the planning-grid build on `expl_enabled` (planning must work without exploration).
- Do not call `OccGrid2D::data()`.

---

## Phase 2 — Plan to the true goal + horizon path-truncation

**Objective:** A\* plans to `G_term`; the local subgoal (`local_G` / `global_path.back()`)
is derived by truncating the returned global path at `horizon`. Requires Phase 1.

### What to implement
1. **Extract a pure, unit-testable helper** in `mighty_utils` (so we can test without ROS):
   - Signature (add to `include/mighty/utils.hpp`, define in `src/mighty/utils.cpp`):
     `bool truncateGlobalPathAtHorizon(vec_Vecf<3>& path, const Eigen::Vector3d& start, double horizon, Eigen::Vector3d& subgoal_out);`
   - Behavior: if `path.size() < 2` return false. Compute `subgoal_out = getLastIntersectionWithSphere(path, horizon, start)` (`hgp/utils.cpp:697`, no-`Jdist` overload). Find the last path index inside the sphere, rebuild `path` as `[0..idx]` then overwrite `path.back() = subgoal_out` — **mirror `mighty.cpp:315-322`**. If the whole path is within `horizon`, `subgoal_out == path.back()` (goal reached within horizon) — leave path intact, return true.
   - Reuse geometry from `getLastIntersectionWithSphere`; do not hand-roll intersection math.

2. **Modify `generateGlobalPath` (`mighty.cpp:807-956`), gated on `use_persistent_global_map && vehicle_type=="ground_robot"`:**
   - Keep the `computeG(...)` call at `:841` for the **UAV / flag-off** path only. For the new path, you still need a `local_G` for `setupHGPPlanner`/`freeGoal`; set it to `local_G_term` (the true goal) before solve — pin z to `default_goal_z` as done at `:878-881`.
   - Change the `solveHGP` goal arg at `:908-910` from `local_G.pos` to `local_G_term.pos`.
   - **After** `solveHGP` succeeds and **before** the `mtx_global_path_` store (`:933-936`):
     - If `corridor_hop_enabled`: leave the existing `computeG_corridorHop(...)` block (`:924-930`) to own truncation + `setG` (it already truncates to a backed-off subgoal). Do **not** also horizon-truncate.
     - Else: call `truncateGlobalPathAtHorizon(global_path, local_A.pos, par_.horizon, subgoal)`, build a `state local_G` with `pos = subgoal`, `yaw` = heading along the last kept segment (copy the yaw idiom from `mighty.cpp:310` or `:333`), and `setG(local_G)`.
   - Leave the direction-hint block (`:883-903`) as-is (harmless; still biases the first segment).

3. **Keep the seed consistent:** in `setTerminalGoal` (`mighty.cpp:1789-1790`), when `use_persistent_global_map && ground_robot`, seed `G_.pos = term_goal.pos` (the true goal) instead of the projection, so the first replan before a path exists is consistent. Guard with the flag; leave UAV/off path untouched.

### Documentation references
- Truncation pattern to mirror: `mighty.cpp:315-322` and commit `:337`.
- `getLastIntersectionWithSphere`: `hgp/utils.cpp:697-716`, decl `hgp/utils.hpp:269`.
- solveHGP goal arg: `mighty.cpp:908-910`; signature `hgp_manager.cpp:259`.
- Optimizer endpoint (auto-retarget): `mighty.cpp:1078-1082`, `:1105`.
- corridor-hop ownership of truncation: `mighty.cpp:924-930`, `:305-337`.
- Seed site: `mighty.cpp:1789-1790`.

### Verification checklist
- [ ] Build succeeds.
- [ ] New unit test `test_path_truncation` passes (see below) — register per `CMakeLists.txt:187-198`, adding `src/mighty/utils.cpp` (and any deps of `getLastIntersectionWithSphere`, e.g. `src/hgp/utils.cpp`) to the gtest sources.
- [ ] Test cases: (a) goal beyond horizon → path truncated, `subgoal` on the sphere, `path.back()==subgoal`; (b) goal within horizon → path unchanged, `subgoal==goal`; (c) `path.size()<2` → returns false.
- [ ] In the running ground-robot sim: log the chosen global path length and first/last points each replan; confirm the route is **stable frame-to-frame** near an obstacle (no left/right homotopy flip) and the robot makes steady forward progress.
- [ ] `grep -n "solveHGP" src/mighty/mighty.cpp` — confirm the ground-robot/persistent path passes `local_G_term.pos`.
- [ ] Flag **off** ⇒ identical to pre-Phase-2 (still `computeG` projection + `local_G.pos` goal).

### Anti-pattern guards
- Do not remove `computeG` / the projection for the UAV or flag-off paths.
- Do not double-truncate when `corridor_hop_enabled` is on.
- Do not forget `setG(...)` after truncation (yaw / GOAL_SEEN / needReplan depend on `G_`).
- Do not use the `Jdist` overload of `getLastIntersectionWithSphere` in the hot path.

---

## Phase 3 — (Optional) Reuse-until-invalid gate

**Objective:** With map+goal now stable, A\* already returns a consistent route each
tick, so this is a pure **compute** optimization (skip the re-solve), not a
correctness fix. Implement only if per-tick A\* over the persistent grid is too slow
(measure first — see Phase 4).

### What to implement (sketch; expand into full tasks only if pursued)
- Cache `global_path_` + the `G_term` it was planned for + a timestamp.
- In `replan()` (`mighty.cpp:708`), between `needReplan` (`:742`) and `generateGlobalPath` (`:751`): if the cached path is valid, reuse it (re-anchor to current `A`, truncate at horizon) and skip A\*. Validity = goal unchanged AND collision-free vs current 2D map (`get2DOccupancy` along waypoints) AND robot near the path AND age < timeout. These are the Nav2 `GoalUpdated` / `IsPathValid` / `PathExpiringTimer` triggers.
- Add params `reuse_global_path` (bool), `global_path_max_age_sec` (double), `global_path_max_deviation_m` (double).

### Verification checklist
- [ ] A\* solve count per second drops materially with reuse on (add a counter/log).
- [ ] A newly-observed obstacle on the committed path forces a fresh solve (drive the robot at a wall that appears; confirm replan).
- [ ] No behavioral regression vs Phase 2 with reuse off.

---

## Phase 4 — Final Verification

1. **Build:** `cd ~/code/mighty_ws && colcon build --packages-select mighty --cmake-args -DCMAKE_BUILD_TYPE=Release` — clean (only pre-existing `binomial` warnings).
2. **Tests:** run `./build/mighty/test_path_truncation`, `./build/mighty/test_frontier_detector`, `./build/mighty/test_frontier_manager`, `./build/mighty/test_occ_grid_2d_unknown` — all green.
3. **Anti-pattern greps:**
   - `grep -rn "OccGrid2D::data\|->data()" src/mighty/mighty_node.cpp` — ensure `data()` is only called on `visited_map_`, never on an `OccGrid2D`.
   - `grep -n "getLastIntersectionWithSphere" src/mighty` — ensure the `Jdist` overload isn't used in `generateGlobalPath`.
   - `grep -n "use_persistent_global_map" src/ include/ config/` — all sites present; every new logic branch gated by it + `vehicle_type=="ground_robot"`.
4. **Compute budget (decides Phase 3):** with the flag on, log `generateGlobalPath` timing (`par_.debug_verbose` prints `[TIMING] Global Planning`). Confirm it fits the replan period; if not, do Phase 3.
5. **Sim acceptance (user-run, cannot be self-verified):**
   - Ground-robot exploration; robot no longer flips left/right around obstacles and makes steady progress.
   - Plan visibly uses the persistent map when re-entering explored area.
   - Flag off ⇒ old behavior.
   - Stale-obstacle caveat noted: `visited_map_` has no raytrace clearing beyond re-observation — acceptable for static forests; revisit if dynamic obstacles are added.

---

## Known caveats / open decisions (carry forward)

- **Fusion detail (Phase 1):** pure `visited_map_` vs `visited_map_` overlaid with the current `occ_2d` occupied cells. `absorb()` already keeps them fresh near the robot, so pure `visited_map_` should suffice; overlay only if newly-appeared obstacles lag by a frame.
- **Map extent:** anything beyond the 100×100 m `visited_map_` isn't planned over. Bump `exploration.visited_map.width_m/height_m` for larger missions (~445 KB per 100 m²).
- **Stale obstacles:** no raytrace clearing (Nav2 has it); fine for static exploration.
- **Struct/config home:** `use_persistent_global_map` can live near `global_planner` (planner setting) or the corridor block (ground-robot behavior) — pick one and be consistent.
