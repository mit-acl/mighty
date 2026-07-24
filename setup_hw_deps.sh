#!/usr/bin/env bash
# setup_hw_deps.sh — clone/pin the dependencies the CONTAINERIZED hardware stack
# (docker/mighty_hw.sh build) consumes, at the versions pinned in mighty.repos
# (single source of truth — this script has no version table of its own).
#
# Why not `vcs import src < mighty.repos` / setup.sh:
#   - a blanket vcs import drops DecompROS2 / livox_ros_driver2 / Livox-SDK2 into
#     mighty_ws/src, where colcon would try (and fail) to build the livox driver.
#     This script routes each repo to the directory the hw Dockerfile expects.
#   - setup.sh targets the v0.0.5 sim release, builds everything natively on the
#     host (the hw container builds in-image), and edits ~/.bashrc.
#
# Routing (CODE_DIR = parent of mighty_ws, normally /home/swarm/code):
#   DecompROS2                  -> decomp_ws/src        (compiled in-image)
#   Livox-SDK2                  -> CODE_DIR             (compiled in-image)
#   livox_ros_driver2           -> SKIPPED              (container bakes the
#       PREBUILT livox_ws/install; the pinned source does not compile against
#       the pinned SDK — copy install/ from a donor rover, see check below)
#   sim-only repos              -> SKIPPED for hw (dockerignored, BUILD_SIMULATION=OFF)
#   everything else             -> mighty_ws/src        (compiled in-image)
#
# Idempotent: existing repos are fetched and checked out to the pin (repos with
# local modifications are refused, never touched). If a remote branch tip equals
# the pin, a local tracking branch is attached instead of a detached HEAD.
# Run as a user with GitLab+GitHub SSH access (fleet convention: blakete).
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CODE_DIR="${CODE_DIR:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
MIGHTY_SRC="$CODE_DIR/mighty_ws/src"
DECOMP_SRC="$CODE_DIR/decomp_ws/src"
REPOS_FILE="$SCRIPT_DIR/mighty.repos"

SIM_ONLY="gazebo_ros_pkgs livox_laser_simulation_ros2 realsense_gazebo_plugin uav_simulator"

echo "CODE_DIR:   $CODE_DIR"
echo "REPOS_FILE: $REPOS_FILE"
mkdir -p "$MIGHTY_SRC" "$DECOMP_SRC"

# name url version, one repo per line, straight from the manifest
parse_repos() {
  python3 - "$REPOS_FILE" <<'PYEOF'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
for name, r in d["repositories"].items():
    print(name, r["url"], r["version"])
PYEOF
}

parse_repos | while read -r name url pin; do
  case "$name" in
    livox_ros_driver2)
      echo ""; echo "=== $name: SKIP (prebuilt livox_ws/install is baked; never rebuilt)"
      continue ;;
    DecompROS2)  dest="$DECOMP_SRC/$name" ;;
    Livox-SDK2)  dest="$CODE_DIR/$name" ;;
    *)
      if grep -qw "$name" <<<"$SIM_ONLY"; then
        echo ""; echo "=== $name: SKIP (sim-only; dockerignored for the hw image)"
        continue
      fi
      dest="$MIGHTY_SRC/$name" ;;
  esac

  echo ""
  echo "=== $name -> ${pin:0:9}  ($dest)"
  if [ -d "$dest/.git" ]; then
    if [ -n "$(git -C "$dest" status --porcelain)" ]; then
      echo "    SKIP: local modifications present — resolve by hand"
      continue
    fi
    git -C "$dest" fetch origin
  else
    git clone "$url" "$dest"
  fi

  # Prefer a tracking branch when some origin branch tip == the pin
  branch="$(git -C "$dest" branch -r --points-at "$pin" --format='%(refname:short)' 2>/dev/null \
            | grep -v HEAD | head -1 | sed 's#^origin/##')"
  if [ -n "$branch" ]; then
    git -C "$dest" checkout -B "$branch" "$pin"
    git -C "$dest" branch --set-upstream-to="origin/$branch" "$branch" 2>/dev/null || true
    echo "    OK: $branch @ $(git -C "$dest" rev-parse --short HEAD)"
  else
    git -C "$dest" checkout --detach "$pin"
    echo "    OK: detached @ $(git -C "$dest" rev-parse --short HEAD)"
  fi
done

echo ""
echo "=== livox prebuilt check ==="
LIVOX_NODE="$CODE_DIR/livox_ws/install/livox_ros_driver2/lib/livox_ros_driver2/livox_ros_driver2_node"
if [ -e "$LIVOX_NODE" ]; then
  echo "    OK: prebuilt livox_ros_driver2 install present"
else
  echo "    MISSING: $CODE_DIR/livox_ws/install"
  echo "    Copy it from a donor rover, e.g.:"
  echo "      rsync -a swarm@rr08.local:/home/swarm/code/livox_ws/install/ $CODE_DIR/livox_ws/install/"
fi

echo ""
echo "Done. Remaining manual steps for a new rover:"
echo "  1. apply any not-yet-committed fleet diffs (see rr08 acl-mapping/dlio)"
echo "  2. edit docker/compose.hw.yaml identity (ROVER_NAME/VEHTYPE/VEHNUM)"
echo "  3. docker/mighty_hw.sh build && docker/mighty_hw.sh start"
