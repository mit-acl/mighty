#!/usr/bin/env bash
# mighty_hw.sh — manual driver for the CONTAINERIZED MIGHTY hardware stack.
# Opt-in for now: the native `mighty` / `kill_mighty` aliases are untouched and
# remain the rollback path.
#
#   mighty_hw.sh start [--dev] [--odom-type dlio|mocap|dlio_in_mocap] [--no-two-d]
#   mighty_hw.sh attach     # tmux attach (Ctrl-b d detaches; stack keeps running)
#   mighty_hw.sh stop       # compose down — tears the whole stack down
#   mighty_hw.sh status     # container + pane status
#   mighty_hw.sh logs       # PID-1 wrapper output (session load messages)
#
# start = docker compose up -d, wait for the session, then attach. The 7-pane
# hw_mighty tmux session runs INSIDE the container; per-pane Ctrl-C / Up /
# Enter restarts work exactly as native. --dev mounts the host src/mighty +
# src/mpc and rebuilds them in-container first (see compose.hw.dev.yaml).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER=hw-mighty
SESSION=hw_mighty
BASE=(-f "${SCRIPT_DIR}/compose.hw.yaml")
DEVF=(-f "${SCRIPT_DIR}/compose.hw.yaml" -f "${SCRIPT_DIR}/compose.hw.dev.yaml")

attach() {
    exec docker exec -it "${CONTAINER}" tmux attach -t "${SESSION}"
}

cmd="${1:-start}"
[[ $# -gt 0 ]] && shift

case "${cmd}" in
    start)
        files=("${BASE[@]}")
        dev=0
        while [[ $# -gt 0 ]]; do
            case "$1" in
                --dev)         files=("${DEVF[@]}"); dev=1 ;;
                --odom-type)   export ODOM_TYPE="$2"; shift ;;
                --two-d-only)  export TWO_D_ONLY=1 ;;
                --no-two-d)    export TWO_D_ONLY= ;;
                *) echo "unknown option: $1" >&2; exit 2 ;;
            esac
            shift
        done
        docker compose "${files[@]}" up -d
        # Wait for the tmux session (dev mode colcon-builds first, so allow long).
        limit=$(( dev == 1 ? 300 : 60 ))
        echo "[mighty_hw] waiting for the ${SESSION} session (up to ${limit}s)..."
        for (( i = 0; i < limit; i += 2 )); do
            if [[ -z "$(docker ps -q --filter "name=^${CONTAINER}$")" ]]; then
                echo "[mighty_hw] container exited; last logs:" >&2
                docker logs --tail 80 "${CONTAINER}" >&2 || true
                exit 1
            fi
            if docker exec "${CONTAINER}" tmux has-session -t "${SESSION}" 2>/dev/null; then
                attach
            fi
            sleep 2
        done
        echo "[mighty_hw] session did not come up within ${limit}s; container logs:" >&2
        docker logs --tail 80 "${CONTAINER}" >&2
        exit 1
        ;;
    attach)
        attach
        ;;
    stop)
        exec docker compose "${BASE[@]}" down
        ;;
    status)
        docker ps --filter "name=^${CONTAINER}$" --format 'container: {{.Names}}  {{.Status}}' | grep . \
            || { echo "container: not running"; exit 1; }
        docker exec "${CONTAINER}" tmux list-panes -t "${SESSION}" \
            -F 'pane #{pane_index}  #{pane_title}  (#{pane_current_command})' 2>/dev/null \
            || echo "no ${SESSION} tmux session in the container"
        ;;
    logs)
        exec docker logs "$@" "${CONTAINER}"
        ;;
    *)
        echo "usage: $0 {start [--dev] [--odom-type T] [--two-d-only|--no-two-d] | attach | stop | status | logs}" >&2
        exit 2
        ;;
esac
