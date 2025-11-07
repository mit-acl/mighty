#!/bin/bash

source ~/.bashrc
source /home/kkondo/code/mighty_ws/install/setup.bash
source /usr/share/gazebo/setup.sh

# (0) sim in static environment
./src/mighty/launch/run_mighty_sim.sh /home/kkondo/code/mighty_ws/install/setup.bash
