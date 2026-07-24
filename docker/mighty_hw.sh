#!/usr/bin/env bash
# mighty_hw.sh — manual driver for the CONTAINERIZED MIGHTY hardware stack,
# HOST-tmux variant (mirrors rover-drive's drive_host_tmux.sh): the hw-mighty
# container just idles (`sleep infinity`, see compose.hw.yaml) and the 7-pane
# hw_mighty tmux session lives on the HOST — each pane docker-execs its node
# into the container. Opt-in: the native `mighty` / `kill_mighty` aliases are
# untouched and remain the rollback path. PREREQ: the host Zenoh router
# (`zenoh_route`) must be running.
#
#   mighty_hw.sh start [--odom-type dlio|mocap|dlio_in_mocap] [--two-d-only|--no-two-d] [--no-sensors]
#   mighty_hw.sh attach     # tmux attach -t hw_mighty (Ctrl-b d detaches; stack keeps running)
#   mighty_hw.sh stop       # kill host session + compose down
#   mighty_hw.sh status     # container + pane status
#   mighty_hw.sh logs       # container PID-1 output (idle loop; panes hold the real logs)
#
# start = docker compose up -d, then build the HOST session with
# run_hw_red_rover.py --docker-exec (run from the host checkout) and attach.
# The ODOM_TYPE/TWO_D_ONLY/NO_SENSORS defaults resolve through compose's
# `environment:` block and are read back from the container, so compose stays
# the single source of the defaults. Per-pane Ctrl-C / Up / Enter restarts work
# exactly as before.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER=hw-mighty
SESSION=hw_mighty
COMPOSE=(-f "${SCRIPT_DIR}/compose.hw.yaml")
LAUNCHER="${SCRIPT_DIR}/../scripts/run_hw_red_rover.py"

attach() {
    exec tmux attach -t "${SESSION}"
}

cmd="${1:-start}"
[[ $# -gt 0 ]] && shift

case "${cmd}" in
    start)
        while [[ $# -gt 0 ]]; do
            case "$1" in
                --odom-type)   export ODOM_TYPE="$2"; shift ;;
                --two-d-only)  export TWO_D_ONLY=1 ;;
                --no-two-d)    export TWO_D_ONLY= ;;
                --no-sensors)  export NO_SENSORS=1 ;;
                *) echo "unknown option: $1" >&2; exit 2 ;;
            esac
            shift
        done
        docker compose "${COMPOSE[@]}" up -d
        echo "[mighty_hw] waiting for the ${CONTAINER} container..."
        for (( i = 0; i < 30; i += 2 )); do
            [[ -n "$(docker ps -q --filter "name=^${CONTAINER}$")" ]] && break
            sleep 2
        done
        if [[ -z "$(docker ps -q --filter "name=^${CONTAINER}$")" ]]; then
            echo "[mighty_hw] container did not come up; last logs:" >&2
            docker logs --tail 40 "${CONTAINER}" >&2 || true
            exit 1
        fi
        # Identity + run flags come from the container env (compose resolved them).
        ROVER_NAME="$(docker exec "${CONTAINER}" printenv ROVER_NAME)"
        ODOM_TYPE="$(docker exec "${CONTAINER}" printenv ODOM_TYPE || true)"
        TWO_D_ONLY="$(docker exec "${CONTAINER}" printenv TWO_D_ONLY || true)"
        NO_SENSORS="$(docker exec "${CONTAINER}" printenv NO_SENSORS || true)"
        export ROVER_NAME
        echo "[mighty_hw] building HOST ${SESSION} session (rover: ${ROVER_NAME})"
        # The launcher loads the session on the host tmux server and attaches.
        exec python3 "${LAUNCHER}" --docker-exec "${CONTAINER}" \
            --odom-type "${ODOM_TYPE:-dlio_in_mocap}" \
            ${TWO_D_ONLY:+--two-d-only} ${NO_SENSORS:+--no-sensors}
        ;;
    attach)
        attach
        ;;
    stop)
        # Session first: closing the panes HUPs the docker-exec'd nodes before
        # the container itself goes away (same ordering as drive_host_tmux.sh).
        tmux kill-session -t "${SESSION}" 2>/dev/null || true
        exec docker compose "${COMPOSE[@]}" down
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
    *)
        echo "usage: $0 {start [--odom-type T] [--two-d-only|--no-two-d] [--no-sensors] | attach | stop | status | logs}" >&2
        exit 2
        ;;
esac
