#!/bin/bash

source ~/.bashrc
source /home/kkondo/code/mighty_ws/install/setup.bash
source /usr/share/gazebo/setup.sh

# (0) sim in static environment
python3 /home/kkondo/code/mighty_ws/src/mighty/launch/run_single_sim.py --env easy_forest

# (1) sim in dynamic environment
# python3 run_single_sim.py --env dynamic_forest