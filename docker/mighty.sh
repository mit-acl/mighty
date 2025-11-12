#!/bin/bash

source ~/.bashrc
source /home/kkondo/code/mighty_ws/install/setup.bash
source /home/kkondo/code/decomp_ws/install/setup.bash
source /usr/share/gazebo/setup.sh

# (0) sim in static environment
tmuxp load /home/kkondo/code/mighty_ws/src/mighty/launch/docker_mighty_sim.yaml
