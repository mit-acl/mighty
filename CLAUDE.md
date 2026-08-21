# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MIGHTY (Hermite Spline-based Efficient Trajectory Planning) is a ROS 2 (Humble) `ament_cmake` package for
real-time trajectory planning with obstacle avoidance. Originally a UAV planner (RA-L 2026 paper), it now
also drives **ground robots** (Pioneer 3-AT in sim, "Red Rover" on hardware) including a **frontier-based
autonomous exploration** stack, single- and multi-agent.

This package lives at `<ws>/src/mighty` inside a colcon workspace (`~/code/mighty_ws`). A second workspace,
`~/code/decomp_ws`, holds DecompROS2 and **must be sourced before building or running** — `decomp_util` is a
hard build dependency.

## Build

```bash
cd ~/code/mighty_ws
source /opt/ros/humble/setup.bash
source ~/code/decomp_ws/install/setup.bash

colcon build --packages-select mighty --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

- `CMAKE_BUILD_TYPE` defaults to `Release` in `CMakeLists.txt` if unset. Always build Release — the planner
  runs a 10 ms replan timer and is unusable in Debug.
- `-DBUILD_SIMULATION=OFF` drops every Gazebo-linked target (`move_model`, `imu_plugin`,
  `convert_velodyne_to_ros_time`, `dynamic_forest_node`). Required on ARM64/Jetson, where Gazebo Classic is
  unavailable. `mighty`, `fake_sim`, and the state converters always build.
- `./setup.sh [-j N] [--jetson]` does full first-time provisioning (ROS 2, apt deps, `vcs import` from
  `mighty.repos`, builds `decomp_ws`/`livox_ws`/`mighty_ws`, patches `~/.bashrc`). Auto-detects Jetson via
  `/etc/nv_tegra_release`. Safe to re-run.
- `mighty.repos` pins every dependency to an exact commit. When a dependency change is needed, update the
  pin there — do not rely on floating branches.

## Tests

GTest suites are registered under `if(BUILD_TESTING)` in `CMakeLists.txt`:
`test_occ_grid_2d_unknown`, `test_frontier_detector`, `test_frontier_manager`, `test_gradient_check`,
`test_formation_cost`.

```bash
cd ~/code/mighty_ws

# All tests for the package
colcon test --packages-select mighty
colcon test-result --verbose            # show failures

# A single suite
colcon test --packages-select mighty --ctest-args -R test_frontier_manager

# Fastest inner loop: run the gtest binary directly (supports --gtest_filter)
./build/mighty/test_frontier_manager --gtest_filter='*Dormant*'
```

`src/test/test_lbfgs_solver.cpp` is **not** a gtest — it is a standalone `main()` benchmark/driver installed
as an executable: `ros2 run mighty test_lbfgs_solver`.

The tested units are deliberately ROS-free (`OccGrid2D`, `VisitedMap`, `FrontierDetector`,
`FrontierManager`, the L-BFGS cost/gradient). Keep new logic in that shape when it needs coverage — anything
pulled into `mighty_node.cpp` becomes effectively untestable.

## Running

`scripts/run_sim.py` generates a **tmuxp** YAML session on the fly and launches it (session name
`mighty_sim`). `--dry-run` prints the YAML without launching — use it to see exactly which nodes/launch
files a mode wires up.

```bash
cd ~/code/mighty_ws
python3 src/mighty/scripts/run_sim.py --mode <MODE> [options]
```

Modes (`--mode`/`-m`):

| Mode | What it runs |
|---|---|
| `interactive` | Single UAV, `fake_sim`, no auto goal — send goals with RViz "2D Goal Pose" |
| `gazebo` | Single agent in Gazebo + ACL mapper; `--ground-robot` swaps the quadrotor for a P3-AT |
| `multiagent` | N agents (`--num-agents`, default 10) on a circle (`--radius`) with `fake_sim` |
| `multiagent-ground` | Same, ground robots |
| `exploration-singleagent-ground` | Frontier exploration, one ground robot in Gazebo (default env `ACL_office`) |
| `exploration-multiagent-ground` | Frontier exploration, N ground robots in a line at y=2 |
| `swap-multiagent-ground` | Goal-swap benchmark via `goal_monitor` (no exploration) |
| `dyn-test`, `dyn-test-ground`, `dyn-test-ground-mpc` | Single agent + dynamic obstacle, for heat-map/tracking debugging |

Common options: `--goal X Y Z`, `--start X Y Z`, `--start-yaw`, `--env <world>` (from `worlds/`),
`--ros-domain-id` (default 20), `--no-rviz`, `--gazebo-gui`, `--no-goal`. `--setup-bash/-s` is optional —
the script auto-detects `<ws>/install/setup.bash` from its own path.

Exploration modes intentionally spawn agents at `y=2.0`, off the origin: `exploreSelectCallback` defers
exploration for any agent within 0.5 m of the origin while peers are active, and spawning at (0,0) deadlocks
it. Preserve that offset when adding exploration modes.

Hardware (ground robot):

```bash
python3 src/mighty/scripts/run_hw_red_rover.py --odom-type {dlio,mocap,dlio_in_mocap} [--goal-type N] [--dry-run]
```

`ROVER_NAME` must be set in the environment (from `VEHTYPE`/`VEHNUM` in `~/.bashrc`).

Other helpers: `scripts/kill_ros_processes.py` (nuke stray ROS processes),
`scripts/bag_record.py` / `hw_bag_record.py`, `benchmarking/unc_benchmark_node.py`, plus analysis notebooks
in `benchmarking/`. Docker equivalents live in `docker/Makefile` (`make build`, `run-interactive`,
`run-gazebo`, `run-multiagent`, `run-mac-*`, `shell`).

## Architecture

### Core planning pipeline (`src/mighty/`, `include/mighty/`)

- **`mighty_node.hpp/cpp`** (~3.9 kLOC) — the ROS layer, and the largest file in the repo. Declares ~246
  parameters, owns all pubs/subs/timers, the frontier-exploration state machine, peer/map sharing, and
  visualization. `replanCallback` fires every 10 ms; a `dc`-rate timer publishes `goal` setpoints.
- **`mighty.hpp/cpp`** — `MIGHTY`, the planner core. Owns the `DroneStatus` state machine
  (`YAWING → TRAVELING → GOAL_SEEN → GOAL_REACHED`) and the per-replan sequence:
  global path (HGP A*) → safe flight corridor (ellipsoid decomposition) → local Hermite-spline trajectory
  (L-BFGS) → `pwp` piecewise polynomial + Bezier control points for sharing.
- **`lbfgs_solver.hpp/cpp`** + **`lbfgs_solver_utils.cpp`** (~3.7 kLOC) — the optimizer: cost terms
  (corridor, ESDF clearance, dynamic-obstacle, formation, jerk), analytic gradients, multiple initial
  guesses solved in parallel (OpenMP). `lbfgs.hpp` is the vendored line-search core (adapted from
  ZJU-FAST-Lab/GCOPTER — see README acknowledgment). Gradients are covered by `test_gradient_check`; if you
  touch a cost term, add a finite-difference case there.
- **`mighty_type.hpp`** — the `parameters` struct. Single source of truth for every tunable.
- **`occ_grid_2d.hpp`** / **`esdf_grid_2d.hpp`** — immutable snapshots built from the ACL mapper's
  `occ_2d_topic` / `esdf_2d_topic` `OccupancyGrid`s. Constructed once per message, shared as
  `shared_ptr<const>` so planner threads never lock. `OccGrid2D` is tristate (UNKNOWN/FREE/OCCUPIED) and
  precomputes a BFS distance heat map for A*; `EsdfGrid2D` pre-converts the mapper's linear int8 encoding to
  float meters for bilinear queries during optimization. **Ground robot only.**

### Global planner (`src/hgp/`, `include/hgp/`)

Formerly "DGP" — any lingering `dgp` reference is stale.

- **`hgp_manager.hpp/cpp`** — owns the 3D voxel map (`map_util.hpp`), obstacle inflation, and safe-flight-
  corridor generation via `decomp_util` ellipsoid/seed decomposition. Thread-safe wrapper around the map.
- **`hgp_planner.hpp/cpp`** — A*/JPS global search with line-of-sight shortcutting, collinear/corner
  removal, and heat-aware Laplacian smoothing. 3D for UAVs, 2D for ground robots (`astar_heat`, `sjps`, …).
- **`graph_search.hpp/cpp`** — the search kernel; consumes `OccGrid2D::occupiedData()` in 2D mode, so
  **unknown cells are traversable** for global planning ("plan through unknown").

### Frontier exploration (ground robot only)

Enabled by `expl_enabled`; ~80 `expl_*` parameters in `mighty_type.hpp`. Pipeline:

1. `occ2DCallback` builds a fresh `OccGrid2D` and absorbs it into **`VisitedMap`** — a fixed-size,
   world-frame, persistent tristate grid. Without it, cells that slide out of the mapper's local window
   re-appear as UNKNOWN and are re-detected as frontiers forever. Detection runs on the visited map when
   `expl_detect_on_visited_map` is set.
2. **`FrontierDetector`** (stateless) — two-pass BFS over the reachable known-free region: outer BFS tags
   free cells adjacent to unknown, inner BFS clusters them 8-connected. Cost is O(reachable cells).
3. **`FrontierManager`** — persistent records keyed by id, with lifecycle
   `ACTIVE → DORMANT → VISITED / INVALIDATED`, EMA centroid merging, and an additive utility ranking
   (size, distance, info gain, revisit penalty, heading alignment — the `expl_w_*` weights).
4. `exploreSelectCallback` (timer at `expl_select_rate_hz`) picks the top frontier and injects it as the
   terminal goal, with pursuit timeouts, keep-out zones after invalidation, and home-return handling.
5. **`PeerTracker`** + `/exploration/peer_poses` and `/exploration/visited_maps` implement multi-agent
   coordination: peers' visited maps are fused in, frontiers near an active peer are rejected or marked
   VISITED, and `expl_use_minpos` enables rank-based allocation.

Threading contract: `FrontierDetector`/`FrontierManager`/`VisitedMap` are **not** internally synchronized.
All access must stay in the `occ2DCallback` + explore-select callback group. `PeerTracker` is the exception
(mutex-guarded).

### Simulation and hardware glue

- **`fake_sim.cpp`** — Gazebo-free simulator: integrates the published trajectory and emits state/odom. Also
  built on Jetson, so it is the way to smoke-test the planner without Gazebo.
- **`convert_odom_to_state.cpp`** / **`convert_vicon_to_state.cpp`** — odom/mocap → `dynus_interfaces/State`.
- **`src/sim/`** — `move_model.cpp` (Gazebo plugin moving dynamic obstacles), `gazebo_ros_imu_sensor.cpp`,
  `dynamic_forest_node.cpp` (simulated dynamic obstacles).
- Python nodes installed via `install(PROGRAMS …)`: `goal_sender.py`, `goal_monitor_node.py`,
  `generate_random_forest.py`, `repub_rviz_2Dgoal.py`, `frame_align_publisher.py`, `formation_viz_node.py`.

### Key topics (all namespaced per agent, e.g. `/NX01/…`)

Subscribes: `state`, `term_goal`, `/swarm_goal`, `sensor_point_cloud`, `occupancy_grid`, `unknown_grid`,
`occ_2d_topic`, `esdf_2d_topic`, `/trajs` (peer trajectories), `predicted_trajs`,
`/exploration/{peer_poses,visited_maps,return_home}`.
Publishes: `goal` (setpoint at `dc` rate), `trajectory`, `/trajs` (own traj, gated by `share_traj`),
`goal_reached`, `mpc_waypoints`, `yaw_output`, plus many markers/point clouds gated by `visual_level`.

### Data flow

1. `stateCallback` (100 Hz) and the map callbacks update `MIGHTY_NODE`'s snapshots.
2. `HGPManager` maintains the voxel map and computes convex free-space corridors.
3. `MIGHTY::replan` produces a Hermite spline inside the corridors via L-BFGS.
4. `goal` setpoints stream to the controller (`fake_sim`, Gazebo, MPC, or hardware).

## Conventions and gotchas

- **Adding a parameter touches five places**: the field in `mighty_type.hpp::parameters`, a
  `declare_parameter` and a matching `get_parameter` in `mighty_node.cpp` (~246 of each, kept in the same
  order), and the YAML files it applies to. A missing YAML entry silently falls back to the declared
  default.
- **Config layering** (`launch/onboard_mighty.launch.py`): base config is `mighty_ground_robot.yaml` when
  `use_ground_robot`, else `mighty.yaml`; when `use_hardware`, `hw_mighty_ground_robot.yaml` /
  `hw_mighty.yaml` is loaded and `dict.update()`-ed **on top**. So the hardware YAMLs are sparse overrides,
  not full configs — and a hardware-only value must go in the hw file, not the base.
- YAMLs are tiered by comment banners: `[CONFIGURE]` (per-vehicle/environment), `[TUNE]`
  (performance/compute trade-offs), `[INTERNAL]` (algorithmic defaults). Keep new parameters in the right
  tier with a unit tag (`[m]`, `[m/s]`, `[s]`, `[-]`).
- Launch files prefer the **source** `rviz/` directory over the installed one, so RViz config edits take
  effect without a rebuild. Editing `config/`, `launch/`, `worlds/`, or `urdf/` does require a rebuild
  (they are installed to `share/`).
- `map_frame_id` is `map` in simulation but `<namespace>/map` on hardware.
- C++ style is enforced by `.clang-format` (Google, 100 cols, left-aligned pointers, include groups
  ordered system → C++ stdlib → Eigen → PCL → ROS → project). Match it; there is no linter in the build.
- Every source file carries the ACL copyright header — copy it into new files.
- `include/sim/exprtk.hpp` (39 kLOC) and `include/hgp/termcolor.hpp` are vendored third-party; do not edit
  or reformat.
- `include/mighty/initial_guess.hpp`, `pure_pursuit.hpp`, and `trajectory_tracker.hpp` are not included by
  any target — dead headers, not part of the build.
