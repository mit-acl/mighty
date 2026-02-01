#!/bin/bash

source ~/.bashrc
source /home/kkondo/code/mighty_ws/install/setup.bash
source /home/kkondo/code/decomp_ws/install/setup.bash
source /usr/share/gazebo/setup.sh

cd /home/kkondo/code/mighty_ws

# Default values
MODE="${MODE:-multiagent}"
GOAL_X="${GOAL_X:-305.0}"
GOAL_Y="${GOAL_Y:-0.0}"
GOAL_Z="${GOAL_Z:-3.0}"
NUM_AGENTS="${NUM_AGENTS:-10}"
ENV="${ENV:-hard_forest}"

# Build arguments
ARGS="--mode $MODE -s /home/kkondo/code/mighty_ws/install/setup.bash"

if [ "$MODE" = "gazebo" ]; then
    ARGS="$ARGS --goal $GOAL_X $GOAL_Y $GOAL_Z --env $ENV"
else
    ARGS="$ARGS --num-agents $NUM_AGENTS"
fi

echo "[Docker] Running: python3 src/mighty/scripts/run_sim.py $ARGS"
python3 src/mighty/scripts/run_sim.py $ARGS
