#!/bin/bash

source ~/.bashrc
source /home/kkondo/code/mighty_ws/install/setup.bash
source /home/kkondo/code/decomp_ws/install/setup.bash
source /usr/share/gazebo/setup.sh

# If USE_XPRA is set, start Xpra server for remote window forwarding
if [ "${USE_XPRA}" = "true" ]; then
    unset LIBGL_ALWAYS_INDIRECT
    export LIBGL_ALWAYS_SOFTWARE=1
    export GALLIUM_DRIVER=llvmpipe
    xpra start :100 \
        --bind-tcp=0.0.0.0:8080 \
        --html=on \
        --encoding=jpeg \
        --quality=85 \
        --speed=70 \
        --min-speed=50
    sleep 5
    export DISPLAY=:100
    echo "[Docker] Xpra server running on port 8080"
    echo "[Docker] Open in browser: http://localhost:8080"
fi

# If USE_FOXGLOVE is set, start the foxglove bridge
if [ "${USE_FOXGLOVE}" = "true" ]; then
    export ROS_DOMAIN_ID=7
    ros2 launch foxglove_bridge foxglove_bridge_launch.xml port:=8765 &
    sleep 2
    echo "[Docker] Foxglove bridge running on ws://localhost:8765"
fi

cd /home/kkondo/code/mighty_ws

# Default values
MODE="${MODE:-multiagent}"
GOAL_X="${GOAL_X:-305.0}"
GOAL_Y="${GOAL_Y:-0.0}"
GOAL_Z="${GOAL_Z:-3.0}"
NUM_AGENTS="${NUM_AGENTS:-10}"
ENV="${ENV:-hard_forest}"
GROUND_ROBOT="${GROUND_ROBOT:-false}"

# Build arguments
ARGS="--mode $MODE -s /home/kkondo/code/mighty_ws/install/setup.bash"

# If using Foxglove, skip RViz inside Docker (Xpra keeps RViz, it just forwards the window)
if [ "${USE_FOXGLOVE}" = "true" ]; then
    ARGS="$ARGS --no-rviz"
fi

if [ "$MODE" = "gazebo" ]; then
    ARGS="$ARGS --goal $GOAL_X $GOAL_Y $GOAL_Z --env $ENV"
    # Single-agent ground robot (Pioneer 3-AT) instead of the UAV
    if [ "$GROUND_ROBOT" = "true" ]; then
        ARGS="$ARGS --ground-robot"
    fi
elif [ "$MODE" = "interactive" ]; then
    : # interactive mode needs no extra args
elif [ "$MODE" = "exploration-singleagent-ground" ]; then
    # Single-agent ground robot autonomous frontier exploration.
    # run_sim.py maps the default env (hard_forest) to ACL_office for this mode.
    ARGS="$ARGS --env $ENV"
else
    ARGS="$ARGS --num-agents $NUM_AGENTS"
fi

echo "[Docker] Running: python3 src/mighty/scripts/run_sim.py $ARGS"
python3 src/mighty/scripts/run_sim.py $ARGS

# Keep the container alive so Xpra/tmux sessions remain accessible
if [ "${USE_XPRA}" = "true" ]; then
    echo "[Docker] Simulation launched. Container will stay alive for Xpra."
    echo "[Docker] Press Ctrl+C to stop."
    tail -f /dev/null
fi
