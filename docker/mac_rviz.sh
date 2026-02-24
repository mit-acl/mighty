#!/bin/bash
# Launch native RViz2 on Mac, connecting to MIGHTY running in Docker via Zenoh
#
# Prerequisites (one-time setup):
#   micromamba create -n ros_rviz -c conda-forge -c robostack-staging -c robostack-humble \
#     ros-humble-rviz2 ros-humble-rmw-cyclonedds-cpp ros-humble-tf2-tools -y
#
# Usage:
#   1. Start Docker:  make run  (in docker/ directory)
#   2. Run this script: ./mac_rviz.sh
#   3. RViz2 opens with topics from Docker

set -e

# Activate RoboStack environment
export MAMBA_EXE="${HOME}/.local/bin/micromamba"
export MAMBA_ROOT_PREFIX="${HOME}/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2>/dev/null)"
micromamba activate ros_rviz

# Must use CycloneDDS (FastDDS crashes on macOS ARM64)
export RMW_IMPLEMENTATION=rmw_cyclonedds_cpp
export ROS_DOMAIN_ID=7
export ROS_DISTRO=humble

# Start Zenoh bridge connecting to Docker (force IPv4 to avoid connecting to self on IPv6)
echo "[Mac] Starting Zenoh bridge (connecting to Docker on 127.0.0.1:7447)..."
"${HOME}/.local/bin/zenoh-bridge-ros2dds" -e tcp/127.0.0.1:7448 --no-multicast-scouting -l "" &
ZENOH_PID=$!
sleep 3

echo "[Mac] Zenoh bridge running (PID: $ZENOH_PID)"
echo "[Mac] Starting RViz2..."

# Launch RViz2 (use MIGHTY's rviz config if available)
RVIZ_CONFIG="$(dirname "$0")/../rviz/mighty.rviz"
if [ -f "$RVIZ_CONFIG" ]; then
    rviz2 -d "$RVIZ_CONFIG"
else
    rviz2
fi

# Cleanup
kill $ZENOH_PID 2>/dev/null
echo "[Mac] Done."
