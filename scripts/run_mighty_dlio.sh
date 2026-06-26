#!/usr/bin/env bash
# =============================================================================
# run_mighty_dlio.sh   (RR06)
#
# Bring up the MIGHTY hardware autonomy stack on RR06 using PURE DLIO odometry
# (no mocap). Pure DLIO seeds at origin (0,0) and needs no mocap, unlike
# dlio_in_mocap (the stock `mighty` alias), which needs /RR06/world publishing
# at startup to seed DLIO.
#
# Differences from the stock `mighty` alias:
#   1. --odom-type dlio   -> pure DLIO. Seeds at origin; needs NO mocap.
#   2. `colcon build` ENABLED -> rebuilds on each launch (rr06 runs from source).
#                                Comment it out below for a config-only relaunch.
#
# The launcher (run_hw_red_rover.py) internally does:
#     tmux  kill-session -t hw_mighty      # tear down any existing stack
#     tmuxp load -d <generated.yaml>       # start the 7-pane session DETACHED
#     tmux  attach-session -t hw_mighty    # attach (no-op without a TTY)
# so running this cleanly replaces an already-running stack.
# =============================================================================

WS=/home/swarm/code/mighty_ws

# --- ROS 2 / Zenoh environment (RR06) ----------------------------------------
# In an interactive `swarm` login shell these are already exported by ~/.bashrc.
# Repeated here so the script ALSO works non-interactively (ssh "bash ...", cron),
# where ~/.bashrc is NOT sourced.
export VEHTYPE=RR
export VEHNUM=06
export ROVER_NAME="${VEHTYPE}${VEHNUM}"            # -> RR06
export ROS_DOMAIN_ID=22                            # RR06 on domain 22
export RMW_IMPLEMENTATION=rmw_zenoh_cpp
export ZENOH_ROUTER_CONFIG_URI=/home/swarm/config/zenoh_config_isolated.json
# tmuxp is system-installed at /usr/bin on RR06 (already on PATH).

# --- Source ROS 2 + the mighty workspace overlay -----------------------------
source /opt/ros/humble/setup.bash
source "${WS}/install/setup.bash"

# --- Rebuild (ENABLED) -------------------------------------------------------
# rr06 rebuilds on each launch (matches the stock `mighty` alias' build step).
# Only needed if you changed C++/Python source (NOT just YAML config) -- comment
# out the two lines below for a config-only relaunch.
#
# Build from the workspace root so colcon discovers the --packages-select
# packages and writes build/install/log into ${WS}, regardless of the
# caller's CWD (makes this script safe to run from any directory).
cd "${WS}"
colcon build --packages-select \
    mighty mpc global_mapper global_mapper_ros direct_lidar_inertial_odometry \
    --cmake-args -DCMAKE_BUILD_TYPE=Release
source "${WS}/install/setup.bash"

# --- Launch ------------------------------------------------------------------
# Runs in the FOREGROUND and attaches to the tmux session (like the `mighty`
# alias). Detach from tmux with Ctrl-b d; the stack keeps running.
cd "${WS}"
python3 src/mighty/scripts/run_hw_red_rover.py --odom-type dlio --two-d-only

# -----------------------------------------------------------------------------
# To launch DETACHED instead (survives your shell closing):
#   nohup setsid bash /home/swarm/run_mighty_dlio.sh </dev/null \
#         >/tmp/mighty_relaunch.log 2>&1 &
#   tmux attach -t hw_mighty            (Ctrl-b d to detach)
# -----------------------------------------------------------------------------
#
# CONFIG STATE (RR06, as of 2026-06-26):
#   * Frontier exploration is ON: hw_mighty_ground_robot.yaml exploration.enabled=true
#     -> auto-explore mode. Set exploration.enabled=false for manual-goal only.
#   * Geofence is OFF: exploration.bounds.enabled=false (src + install copies).
#     Flip to true (with the bounds.min/max_x/y) to constrain exploration to a box.
#   * global_mapper z_ground = -0.35 in BOTH src and install copies of
#     acl-mapping/global_mapper_ros/cfg/hw_ground_robot.yaml (low-obstacle fix APPLIED,
#     matches RR08); min_occupied_neighbors = 1. z_ground is init-only -> relaunch to load.
#
# LOW-OBSTACLE FIX (APPLIED on RR06, parity with RR08):
#   The red-rover frame convention has DLIO anchor base_link at z=0 with base_link->lidar
#   translation [0,0,0], so the lidar is modeled at map z~=0 while it sits ~0.5 m above the
#   floor. The 3D occupancy PUBLISH gate in global_mapper_ros.cc
#   (`if (xyz[2] <= z_ground) continue;`) drops every occupied voxel at/below z_ground; the
#   same gate feeds mighty's 2D A* grid (global_mapper.cc:570). With the old z_ground=0.0,
#   low obstacles (shorter than ~0.67 m above the floor) were never mapped and the rover
#   could drive into them. Fix: z_ground -> -0.35 (keeps min_occupied_neighbors=1).
#   NOTE: colcon build is ENABLED above, so src copies the cfg over install each launch --
#   keep the src value at -0.35 or a rebuild will clobber the fix back to 0.0.
