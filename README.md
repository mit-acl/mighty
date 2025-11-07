# MIGHTY: Quintic Hermite Spline-based Efficient Trajectory Planning #

### **Submitted to the IEEE Robotics and Automation Letters (RA-L)**

| **Trajectory** | **Forest** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_complex_benchmarks.gif" width="350" height="193" alt="Complex Benchmarks"></a> | <a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hard_forest.gif" width="350" height="193" alt="Static Forest"></a> |

| **Dynamic Obstacles** | **Long Flight** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_dynamic_sim.gif" width="350" height="193" alt="Dynamic Obstacles"></a> | <a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hw_long_flight.gif" width="350" height="193" alt="Hardware Long Flight"></a>

| **Fast Flight 1** | **Fast Flight 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hw_fast_flight_1.gif" width="350" height="193" alt="Hardware Fast Flight 1"></a> | <a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hw_fast_flight_2.gif" width="350" height="193" alt="Hardware Fast Flight 2"></a>

| **Dynamic Env 1** | **Dynamic Env 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hw_dynamic_1.gif" width="350" height="193" alt="Hardware Dynamic Env 1"></a> | <a target="_blank" href="https://youtu.be/SI8YbMS-wyw"><img src="./imgs/mighty_gifs_hw_dynamic_2.gif" width="350" height="193" alt="Hardware Dynamic Env 2"></a>

## General Setup

MIGHTY has been tested on both Docker and native installations on Ubuntu 22.04 with ROS 2 Humble.

### Use Docker (Recommended)

1. **Install Docker:**  
   Follow the [official Docker installation guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/).

2. **Clone the Repository and Navigate to the Docker Folder:**
   ```bash
   mkdir -p ~/code/ws/src
   cd ~/code/ws/src
   git clone https://github.com/kotakondo/dynus.git
   cd src/dynus/docker
   git checkout lbfgs
   ```
3. **Generate a Personal Access Token (PAT):**
    - Follow [this guide](https://www.geeksforgeeks.org/how-to-generate-personal-access-token-in-github/) to generate your PAT.
    - When generating your PAT, ensure:
      - **Repository access:** Set to **All repositories**
      - In the **Repository permissions** menu, find the **Contents** row and select **Access > Read and Write**

4. **BUILD:**
    - Navigate to the docker folder in your dynus repo (eg. `cd ~/code/dynus_ws/src/dynus/docker/`) and run this
    - ```bash
      make build GIT_ACCESS_TOKEN=(copy and paste your_token_here)
      ```
5. **Run Simulation**
    - Run the following command to start the simulation:
      ```bash
      make run-sim
      ```
6. **(Optional) If you want to run the sim in dynamic environment:**
    - comment out `python3 src/dynus/launch/run_single_sim.py --env easy_forest` in dynus.sh
    - uncomment `python3 src/dynus/launch/run_single_sim.py --env dynamic_forest` in dynus.sh
    - run 
      ```bash
      make run-sim
      ```

<details>
  <summary><b>Useful Docker Commands</b></summary>

  - **Remove all caches:**  
    ```bash
    docker builder prune
    ```

  - **Remove all containers:**  
    ```bash
    docker rm $(docker ps -a -q)
    ```

  - **Remove all images:**  
    ```bash
    docker rmi $(docker images -q)
    ```

</details>

5. **Run the Docker Container:**
   ```bash
   make run-sim
   ```
   This command will start the simulation in the Docker container. If you encounter any issues, please refer to the [Error Handling](#error-handling) section.
   

### Native Installation

1. **Clone the Repository and Navigate to the Workspace Folder:**
   ```bash
   mkdir -p ~/code/ws
   cd ~/code/ws
   git clone https://github.com/kotakondo/dynus.git
   cd dynus
   ```

2. **Run the Setup Script:**
   ```bash
   ./setup.sh
   ```
   This script will first install ROS 2 Humble, then DYNUS and its dependencies. Please note that this script modifies your `~/.bashrc` file.

<details>
  <summary><b>Error Handling</b></summary>

  - **Error 1: When running `colcon build`**
    ```
    By not providing "Finddecomp_util.cmake" in CMAKE_MODULE_PATH, this project
    has asked CMake to find a package configuration file provided by
    "decomp_util", but CMake did not find one.
    ```
    - **Solution:**  
      Source `install/setup.bash` and build again.

  - **Error 2: When running Python script**  
    ```
    No module named 'rclpy._rclpy_pybind11'
    ```
    - **Solution:**  
      Deactivate the Conda environment if you are using Conda.

  - **Runtime Error 3:**  
    ```
    "Spawn status: Entity pushed to spawn queue, but spawn service timed out
    waiting for entity to appear in simulation under the name [quadrotor]"
    ```
    - **Solution:**  
      Go to the world file you are using and make sure `sim_time` is set to `0`.

  - **Error 4: When building dlio package**  
    ```
    fatal error: numpy/ndarrayobject.h: No such file or directory
    ```
    - **Solution:**  
      Just cleaned up the workspace and rebuilt it.

  - **Error 5: Gazebo runtime error**
    ```
    [gzserver-1] gzserver: /usr/include/boost/smart_ptr/shared_ptr.hpp:728: 
    typename boost::detail::sp_member_access<T>::type boost::shared_ptr<T>::operator->() const 
    [with T = gazebo::rendering::Scene; typename boost::detail::sp_member_access<T>::type = gazebo::rendering::Scene*]: 
    Assertion `px != 0' failed.
    ```
    - **Solution:**  
      Source Gazebo's setup file:  
      ```
      source /usr/share/gazebo/setup.bash
      ```

  - **Error 6: Missing package configuration file**  
    ```
    Could not find a package configuration file provided by "diagnostic_updater" 
    with any of the following names: when building realsense-ros.
    ```
    - **Solution:**  
      Install the missing package:  
      ```
      sudo apt-get install ros-humble-diagnostic-updater
      ```
</details>

### Run the Simulation

Change the /path/to/install/ to your actual install path (eg. `./src/dynus/launch/run_mighty.sh /home/kkondo/code/dynus_ws/install/setup.bash`)
```bash
./src/dynus/launch/run_mighty.sh /path/to/install/setup.bash
```

### Notes

<details>
  <summary><b>ROS2 Multiagent Hardware Setup</b></summary>

  **Debugging Useful Commands:**
  - `ros2 multicast send`
  - `ros2 multicast receive`
  - `ros2 run demo_nodes_cpp talker`
  - `ros2 run demo_nodes_cpp listener`

  **Commands to Run:**
  - `sudo ufw disable`
  - `sudo ufw allow in proto udp to 224.0.0.0/4`  
    *(if you disable ufw, you don't need this)*
  - `sudo ufw allow in proto udp from 224.0.0.0/4`  
    *(if you disable ufw, you don't need this)*
  - `sudo ufw allow in proto udp from 192.168.1.0/24`  
    *(if you disable ufw, you don't need this)*

  **Things to Check:**
  - `export RMW_IMPLEMENTATION=rmw_fastrtps_cpp`
  - `export ROS_LOCALHOST_ONLY=0`
  - Ensure `ROS_DOMAIN_ID` is set to the same value for all agents.
  - *(For ground station)* You may need to unplug other network interfaces to encourage the computer to use the correct network.
</details>

<details>
  <summary><b>How to Monitor Compute</b></summary>

  - **Record Performance Data:**
    ```bash
    perf record -g -- tmuxp load src/dynus/launch/default_sim.yaml
    ```
  - **Generate Performance Script:**
    ```bash
    perf script >> script.txt
    ```
  - **Generate Flame Graph:**
    ```bash
    cat /home/kkondo/code/dynus_ws/script.txt | ./stackcollapse-perf.pl | ./flamegraph.pl > flame.html
    ```
  - **View Flame Graph:**  
    Open `flame.html` in your browser (e.g., type the path in your browser).
</details>

<details>
  <summary><b>Gazebo</b></summary>

  - **Hospital Models & World Files:**
    - [OpenRobotics Hospital Models](https://app.gazebosim.org/OpenRobotics/fuel/collections/Hospital)
    - [Useful World Files for Gazebo and ROS 2 Simulations](https://automaticaddison.com/useful-world-files-for-gazebo-and-ros-2-simulations/)
</details>

<details>
  <summary><b>For Octomap</b></summary>

  - **Installation Commands:**
    ```bash
    sudo apt-get install ros-humble-octomap
    sudo apt-get install ros-humble-octomap-msgs
    sudo apt-get install ros-humble-octomap-ros
    sudo apt-get install ros-humble-octomap-rviz-plugins
    sudo apt install ros-humble-rviz-common
    ```
</details>

<details>
  <summary><b>If You Use Conda for Dynus</b></summary>

  - Use **Python 3.10 or below** *(to avoid pybind errors)*.
  - Install GCC 12.1.0 to prevent gcc errors:
    ```bash
    conda install -c conda-forge gcc=12.1.0
    ```
  - Install required Python packages:
    ```bash
    pip install numpy lxml pyyaml empy==3.3.4 catkin_pkg lark
    ```
</details>

<details>
  <summary><b>How to Build Docker Using Private Repo with PAT</b></summary>

  - **Generate a Personal Access Token (PAT):**
    - Follow [this guide](https://www.geeksforgeeks.org/how-to-generate-personal-access-token-in-github/) to generate your PAT.
    - When generating your PAT, ensure:
      - **Repository access:** Set to **All repositories**
      - In the **Repository permissions** menu, find the **Contents** row and select **Access > Read and Write**
</details>

 <summary><b>Docker-related Issues</b></summary>
  - Rviz GUI not showing up? Run this in the terminal:
    ```bash
    xhost +
    ```
<details>
  <summary><b>Object Detection Module</b></summary>

  - **Install Dependencies:**
    ```bash
    pip install opencv-python
    pip install ultralytics
    ```
  - **Download YOLO Model:**
    - Use `yolo11n.pt` (for example):
      ```bash
      wget https://github.com/ultralytics/assets/releases/download/v8.3.0/yolo11n.pt
      ```
  - **Install Specific Numpy Version:**
    ```bash
    pip install numpy==1.26.4
    ```
</details>

<details>
  <summary><b>New Conda Environment for YOLO</b></summary>

  - **Create and Configure Environment:**
    ```bash
    conda install -n yolo -c conda-forge libstdcxx-ng
    ```
</details>

<details>
  <summary><b>When you get Found ROS process to kill: defunct </b></summary>
  - **Solution:**
    - Run the following command to kill the defunct process:
      ```bash
      ps -ef | grep defunct
      ```
    - And kill the process like this:
      ```bash
      kill -9 <PID> <PPID>
      ```
</details>

<details>
  <summary><b>Goal Commands</b></summary>

  - **Office Space:**
    ```bash
    ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 36.0, y: 42.0, z: 2.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
    ```
  - **Forest3:**
    ```bash
    ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 50.0, y: 0.0, z: 3.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
    ```
  - **Path Push Visualization:**
    ```bash
    ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 8.0, y: 0.0, z: 3.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
    ```
  - **Ground Robot:**
    ```bash
    ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 50.0, y: 0.0, z: 0.5}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
    ```
</details>

<details>
  <summary><b>Bag Record Commands</b></summary>

  - **Path Push Visualization Bag:**
    ```bash
    ros2 bag record /NX01/dgp_path_marker /tf /tf_static /rosout /NX01/point_G /NX01/cluster_bounding_boxes /NX01/tracked_obstacles /NX01/fov /NX01/uncertainty_spheres -o test7
    ```
  - **Ground Robot Visualization Bag:**
    ```bash
    ros2 bag record /NX01/d435/color/image_raw /NX01/fov /NX01/point_G /NX01/actual_traj /NX01/mid360_PointCloud2 /tf /tf_static -o test0
    ```
</details>

<details>
  <summary><b>(ACL Specific) Running a Large Scale Simulation on Lambda Machines</b></summary>

  **Setup:**
  - Connect a monitor to the machine.
  - Set the display:
    ```bash
    export DISPLAY=:1
    ```
  - Run:
    ```bash
    xhost +
    ```
  - *If you see "Depth Camera not found" in `base_mighty.launch..py`, try rebuilding the docker image.*

  **Error Troubleshooting:**
  - **Error:**
    ```
    docker: Error response from daemon: could not select device driver "" with capabilities: [[gpu]].
    make: *** [Makefile:11: run-sim] Error 125
    ```
  - **Solution:**
    - Follow [this StackOverflow thread](https://stackoverflow.com/questions/75118992/docker-error-response-from-daemon-could-not-select-device-driver-with-capab).
    - Use `run-sim-no-gpu` in `dynus/docker/Makefile` *(currently Lambda machines' GPU is not working—should be fixed soon)*.
    - Note: Realsense won’t work without a display. (SSH with `-X` may not work; a physical monitor is required.)
    - Once the simulation is running, the bag file is saved inside Docker. To copy it to your local machine:
      ```bash
      docker cp <container_id>:/path/to/bag/file /path/to/local/machine
      ```
      *(Find `<container_id>` by running `docker ps`.)*

  **Goal Sender Examples:**
  - **Forest3 Goal Sender (Straight Path):**
    ```bash
    ros2 launch dynus goal_sender.launch.py list_agents:="['NX01', 'NX02', 'NX03', 'NX04', 'NX05']" list_goals:="['[50.0, 20.0, 0.0]', '[50.0, 10.0, 0.0]', '[50.0, 0.0, 0.0]', '[50.0, -10.0, 0.0]', '[50.0, -20.0, 0.0]']" default_goal_z:=2.0
    ```
  - **Forest3 Goal Sender (Cross Path):**
    ```bash
    ros2 launch dynus goal_sender.launch.py list_agents:="['NX01', 'NX02', 'NX03', 'NX04', 'NX05']" list_goals:="['[50.0, -20.0, 0.0]', '[50.0, -10.0, 0.0]', '[50.0, 0.0, 0.0]', '[50.0, 10.0, 0.0]', '[50.0, 20.0, 0.0]']" default_goal_z:=2.0
    ```
  - **Big Forest High Res Goal Sender (Straight Path):**
    ```bash
    ros2 launch dynus goal_sender.launch.py list_agents:="['NX01', 'NX02', 'NX03', 'NX04', 'NX05']" list_goals:="['[35.0, 20.0, 0.0]', '[35.0, 10.0, 0.0]', '[35.0, 0.0, 0.0]', '[35.0, -10.0, 0.0]', '[35.0, -20.0, 0.0]']" default_goal_z:=3.0
    ```
  - **Big Forest High Res Goal Sender (Cross Path):**
    ```bash
    ros2 launch dynus goal_sender.launch.py list_agents:="['NX01', 'NX02', 'NX03', 'NX04', 'NX05']" list_goals:="['[35.0, -20.0, 0.0]', '[35.0, -10.0, 0.0]', '[35.0, 0.0, 0.0]', '[35.0, 10.0, 0.0]', '[35.0, 20.0, 0.0]']" default_goal_z:=3.0
    ```
  - **Quadruped Forest3 Goal Sender:**
    ```bash
    ros2 launch dynus goal_sender.launch.py list_agents:="['NX01']" list_goals:="['[45.0, 0.0, 0.0]']" default_goal_z:=1.0
    ```
</details>
