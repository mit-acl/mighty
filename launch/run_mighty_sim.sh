#!/usr/bin/env bash
set -euo pipefail
if [ $# -lt 1 ]; then
  echo "Usage: $0 /path/to/ros_ws/install/setup.bash" >&2
  exit 1
fi
SETUP_BASH="$1" tmuxp load src/mighty/launch/mighty_sim.yaml
