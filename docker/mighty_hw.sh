#!/usr/bin/env bash
# mighty_hw.sh — SINGLE entrypoint for the CONTAINERIZED MIGHTY hardware stack,
# host-tmux variant (sibling of the rover drive/sensors services): the hw-mighty
# container just idles (`sleep infinity`, see compose.hw.yaml) and every node
# runs as a `docker exec` whose controlling pane lives in a HOST tmux session
# 'hw_mighty', built right here in bash — no tmuxp, no run_hw_red_rover.py
# (that script is the NATIVE `mighty` alias path and stays untouched).
#
#   mighty_hw.sh start [--odom-type dlio|dlio_in_mocap|mocap] [--two-d-only|--no-two-d]
#   mighty_hw.sh attach            # tmux attach (Ctrl-b d detaches; stack keeps running)
#   mighty_hw.sh stop              # kill session, then compose down
#   mighty_hw.sh status            # container + pane status
#   mighty_hw.sh logs [...]        # container PID-1 output (idle loop; panes hold the real logs)
#   mighty_hw.sh rebuild [...]     # stop + compose build + start (start flags forwarded)
#
# Panes by --odom-type (default dlio; two_d_only defaults ON, --no-two-d disables):
#   dlio           Onboard MIGHTY | RViz 2D goal | Global Mapper | DLIO | seed pose | tf map->odom
#                  Replicates what dlio.service ran: DLIO anchored by a constant
#                  seed-pose spoof on /<ns>/world at z=${DLIO_SEED_Z:-0.4}.
#   dlio_in_mocap  same minus the seed pose pane — a REAL mocap publishes
#                  /<ns>/world, and DLIO anchors to the FIRST pose it receives,
#                  so a running spoof would race it and win.
#   mocap          no DLIO at all: mighty flips to use_onboard_localization:=false
#                  + twist from mocap/twist, the mapper consumes raw livox/lidar
#                  with pose from /<ns>/world, and two extra static TFs bridge
#                  the mocap frames (<ns> -> <ns>/base_link, world -> <ns>/map).
#
# PREREQS in EVERY mode (this stack launches NO sensors and NO router):
#   drive.service    hosts the zenoh router on :7447
#   sensors.service  livox + D455 + the base_link->lidar static TF
# dlio.service must be DISABLED wherever this stack runs — its nodes and the
# DLIO panes here are the same nodes (double-publish). Rollback path: re-enable
# dlio.service and use the native `mighty` alias.
#
# NOTE: DLIO runs a ~3 s stationary IMU calibration whenever its pane (re)starts
# — keep the rover still. Identity (ROBOT_NAME, RMW, optional DLIO_SEED_Z) comes
# from /etc/rover/rover.env via compose env_file; nothing in this directory is
# rover-specific.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER=hw-mighty
SESSION=hw_mighty
COMPOSE=(docker compose -f "${SCRIPT_DIR}/compose.hw.yaml")

# In-container command preludes. Pane commands are sent single-quoted, so $VARs
# in them expand INSIDE the container (env_file provides ROBOT_NAME etc.) — and
# therefore a pane command must never contain a single quote. Every pane
# spin-waits on the router first; the guard shares the pane's history line so
# Ctrl-C -> Up -> Enter re-runs it too.
SETUP='source /opt/ros/humble/setup.bash && source /home/swarm/code/mighty_ws/install/setup.bash'
DECOMP='source /home/swarm/code/decomp_ws/install/setup.bash'
WAIT_ROUTER='until (echo >/dev/tcp/127.0.0.1/7447) 2>/dev/null; do echo "waiting for zenoh router (drive.service)..."; sleep 2; done'

dx() { echo "docker exec -it ${CONTAINER} bash -c '${SETUP} && ${WAIT_ROUTER} && $1'"; }

attach() { exec tmux attach -t "${SESSION}"; }

usage() {
    echo "usage: $0 {start [--odom-type dlio|dlio_in_mocap|mocap] [--two-d-only|--no-two-d] | attach | stop | status | logs | rebuild [start flags]}" >&2
    exit 2
}

start() {
    local odom_type=dlio two_d=true
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --odom-type)  odom_type="${2:?--odom-type needs a value}"; shift ;;
            --two-d-only) two_d=true ;;
            --no-two-d)   two_d=false ;;
            *) echo "[mighty_hw] unknown option: $1" >&2; usage ;;
        esac
        shift
    done
    case "${odom_type}" in
        dlio|dlio_in_mocap|mocap) ;;
        *) echo "[mighty_hw] bad --odom-type '${odom_type}'" >&2; usage ;;
    esac

    # Fail fast on the hard prereq, warn on the soft one — the pane-side
    # spin-waits still guard every node, this is just early readable feedback.
    if ! (echo >/dev/tcp/127.0.0.1/7447) 2>/dev/null; then
        echo "[mighty_hw] no zenoh router on 127.0.0.1:7447 — start drive.service first" >&2
        exit 1
    fi
    if ! docker ps --format '{{.Names}}' | grep -cx sensors >/dev/null; then
        echo "[mighty_hw] WARNING: no 'sensors' container (sensors.service down?) — livox/DLIO/mapper will sit idle" >&2
    fi

    "${COMPOSE[@]}" up -d
    echo "[mighty_hw] waiting for the ${CONTAINER} container..."
    local i
    for (( i = 0; i < 30; i += 2 )); do
        [[ -n "$(docker ps -q --filter "name=^${CONTAINER}$")" ]] && break
        sleep 2
    done
    if [[ -z "$(docker ps -q --filter "name=^${CONTAINER}$")" ]]; then
        echo "[mighty_hw] container did not come up; last logs:" >&2
        docker logs --tail 40 "${CONTAINER}" >&2 || true
        exit 1
    fi
    local robot_name
    robot_name="$(docker exec "${CONTAINER}" printenv ROBOT_NAME)"

    # ---- per-mode node commands (mind the single-quote rule above) ----------
    # Consumers gate on the TFs they need via wait_for_tf.py (the /tf_static
    # startup-race fix): the static publishers below latch one transient-local
    # sample, and a subscriber that forms too early would never receive it.
    local mighty_cmd mapper_cmd
    if [[ "${odom_type}" == mocap ]]; then
        mighty_cmd="${DECOMP}"' && ros2 run mighty wait_for_tf.py $ROBOT_NAME/map $ROBOT_NAME/odom && ros2 launch mighty onboard_mighty.launch.py x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:=$ROBOT_NAME use_hardware:=true use_onboard_localization:=false robot_type:=red_rover depth_camera_name:=d455 twist_topic:=mocap/twist'
        mapper_cmd='ros2 run mighty wait_for_tf.py $ROBOT_NAME/map $ROBOT_NAME/odom $ROBOT_NAME/base_link $ROBOT_NAME/lidar world $ROBOT_NAME/map && ros2 launch global_mapper_ros global_mapper_node.launch.py hardware:=true ground_robot:=true quad:=$ROBOT_NAME depth_pointcloud_topic:=livox/lidar pose_topic:=world pose_type:=pose_stamped use_obstacle_tracker:=false'
    else  # dlio | dlio_in_mocap — identical mighty/mapper; only the seed differs
        mighty_cmd="${DECOMP}"' && ros2 run mighty wait_for_tf.py $ROBOT_NAME/map $ROBOT_NAME/odom && ros2 launch mighty onboard_mighty.launch.py x:=0.0 y:=0.0 z:=0.0 yaw:=0.0 namespace:=$ROBOT_NAME use_hardware:=true use_onboard_localization:=true robot_type:=red_rover depth_camera_name:=d455'
        mapper_cmd='ros2 run mighty wait_for_tf.py $ROBOT_NAME/map $ROBOT_NAME/odom $ROBOT_NAME/base_link $ROBOT_NAME/lidar && ros2 launch global_mapper_ros global_mapper_node.launch.py hardware:=true ground_robot:=true quad:=$ROBOT_NAME depth_pointcloud_topic:=dlio/odom_node/pointcloud/deskewed use_sim_time:=false pose_topic:=dlio/odom_node/pose use_obstacle_tracker:=false'
    fi
    local dlio_cmd='ros2 launch direct_lidar_inertial_odometry dlio.launch.py namespace:=$ROBOT_NAME initial_pose_topic:=world two_d_only:='"${two_d}"
    local seed_cmd='ros2 topic pub -r 2 /$ROBOT_NAME/world geometry_msgs/msg/PoseStamped "{header: {frame_id: world}, pose: {position: {z: ${DLIO_SEED_Z:-0.4}}, orientation: {w: 1}}}"'
    local tf_map_odom='ros2 run tf2_ros static_transform_publisher 0 0 0 0 0 0 $ROBOT_NAME/map $ROBOT_NAME/odom'
    local tf_mocap_base='ros2 run tf2_ros static_transform_publisher --frame-id $ROBOT_NAME --child-frame-id $ROBOT_NAME/base_link'
    local tf_world_map='ros2 run tf2_ros static_transform_publisher --frame-id world --child-frame-id $ROBOT_NAME/map'

    local -a titles cmds
    titles=('Onboard MIGHTY' 'RViz 2D goal' 'Global Mapper')
    cmds=(
        "$(dx "${mighty_cmd}")"
        "$(dx 'ros2 run mighty repub_rviz_2Dgoal.py')"
        "$(dx "${mapper_cmd}")"
    )
    case "${odom_type}" in
        dlio)
            titles+=('DLIO (seed-anchored)' 'seed pose' 'tf map->odom')
            cmds+=("$(dx "${dlio_cmd}")" "$(dx "${seed_cmd}")" "$(dx "${tf_map_odom}")")
            ;;
        dlio_in_mocap)
            titles+=('DLIO (mocap-seeded)' 'tf map->odom')
            cmds+=("$(dx "${dlio_cmd}")" "$(dx "${tf_map_odom}")")
            ;;
        mocap)
            titles+=('tf mocap->base_link' 'tf world->map' 'tf map->odom')
            cmds+=("$(dx "${tf_mocap_base}")" "$(dx "${tf_world_map}")" "$(dx "${tf_map_odom}")")
            ;;
    esac

    # ---- build the HOST session (drop any stale one first). -x/-y: a detached
    # session started from a non-tty needs an explicit size; an attaching
    # client resizes it anyway.
    tmux kill-session -t "${SESSION}" 2>/dev/null || true
    tmux new-session -d -s "${SESSION}" -n main -x 220 -y 56
    tmux set-option -w -t "${SESSION}:main" pane-border-status top
    tmux set-option -w -t "${SESSION}:main" pane-border-format ' #{pane_index}: #{pane_title} '
    for (( i = 1; i < ${#cmds[@]}; i++ )); do
        tmux split-window -t "${SESSION}:main"
        tmux select-layout -t "${SESSION}:main" tiled
    done
    local -a panes
    mapfile -t panes < <(tmux list-panes -t "${SESSION}:main" -F '#{pane_id}')
    for i in "${!cmds[@]}"; do
        tmux select-pane -t "${panes[$i]}" -T "${titles[$i]}"
        tmux send-keys -t "${panes[$i]}" "${cmds[$i]}" C-m
    done
    tmux select-pane -t "${panes[0]}"

    echo "[mighty_hw] ${SESSION} session up (rover: ${robot_name}, odom: ${odom_type}, two_d_only: ${two_d})"
    if [[ -t 0 && -t 1 ]]; then
        attach
    else
        echo "[mighty_hw] no tty — attach with: tmux attach -t ${SESSION}"
    fi
}

cmd="${1:-start}"
[[ $# -gt 0 ]] && shift

case "${cmd}" in
    start)
        start "$@"
        ;;
    attach)
        attach
        ;;
    stop)
        # Session first: closing the panes HUPs the docker-exec'd nodes before
        # the container itself goes away.
        tmux kill-session -t "${SESSION}" 2>/dev/null || true
        exec "${COMPOSE[@]}" down
        ;;
    status)
        docker ps --filter "name=^${CONTAINER}$" --format 'container: {{.Names}}  {{.Status}}' | grep . \
            || { echo "container: not running"; exit 1; }
        tmux list-panes -t "${SESSION}" \
            -F 'pane #{pane_index}  #{pane_title}  (#{pane_current_command})' 2>/dev/null \
            || echo "no ${SESSION} tmux session on the host"
        ;;
    logs)
        exec docker logs "$@" "${CONTAINER}"
        ;;
    rebuild)
        "$0" stop || true
        "${COMPOSE[@]}" build
        exec "$0" start "$@"
        ;;
    *)
        usage
        ;;
esac
