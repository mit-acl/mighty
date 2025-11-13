# MIGHTY: Hermite Spline-based Efficient Trajectory Planning #

### **Submitted to the IEEE Robotics and Automation Letters (RA-L)**

| **Trajectory** | **Forest** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_complex_benchmarks.gif" width="350" height="193" alt="Complex Benchmarks"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hard_forest.gif" width="350" height="193" alt="Static Forest"></a> |

| **Dynamic Obstacles** | **Long Flight** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_dynamic_sim.gif" width="350" height="193" alt="Dynamic Obstacles"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_long_flight.gif" width="350" height="193" alt="Hardware Long Flight"></a>

| **Fast Flight 1** | **Fast Flight 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_fast_flight_1.gif" width="350" height="193" alt="Hardware Fast Flight 1"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_fast_flight_2.gif" width="350" height="193" alt="Hardware Fast Flight 2"></a>

| **Dynamic Env 1** | **Dynamic Env 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_dynamic_1.gif" width="350" height="193" alt="Hardware Dynamic Env 1"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_dynamic_2.gif" width="350" height="193" alt="Hardware Dynamic Env 2"></a>

## Paper

## Video

The full video is available [https://youtu.be/Pvb-VPUdLvg](https://youtu.be/Pvb-VPUdLvg).

## Interactive Demo

If you are interested in an interactive demo of MIGHTY, please switch to the `interactive_demo` branch [https://github.com/mit-acl/mighty/tree/interactive_demo] by running:

```bash
git checkout interactive_demo
```
and follow the setup instructions in the README of that branch.

## Fork

Since you might want to use interactive demos, when you fork this repository, please make sure to also include the `interactive_demo` branch by unselecting the "Copy the main branch only" option.

## Setup

MIGHTY has been tested on both Docker and native installations on Ubuntu 22.04 with ROS 2 Humble.

### Use Docker (Recommended)

1. **Install Docker:**  
   Follow the [official Docker installation guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/).

2. **Clone the Repository and Navigate to the Docker Folder:**
   ```bash
   mkdir -p ~/code/ws/src
   cd ~/code/ws/src
   git clone https://github.com/mit-acl/mighty.git
   cd src/mighty/docker
   ```

3. **BUILD:**
    - Navigate to the docker folder in your mighty repo (eg. `cd ~/code/ws/src/mighty/docker/`) and run this
      ```bash
      make build
      ```

4. **Run Simulation**
    - Run the following command to start the simulation:
      ```bash
      make run
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

### Native Installation

1. **Clone the Repository and Navigate to the Workspace Folder:**
   ```bash
   mkdir -p ~/code/ws
   cd ~/code/ws
   git clone https://github.com/mit-acl/mighty.git
   cd mighty
   ```

2. **Run the Setup Script:**
   ```bash
   ./setup.sh
   ```
   This script will first install ROS 2 Humble, then MIGHTY and its dependencies. Please note that this script modifies your `~/.bashrc` file.

 3. **Run the Simulation**
    Run the simulation. You might need to change the path to `setup.bash` to its absolute path (eg. `/home/kkondo/code/ws/install/setup.bash`).
    ```bash
    cd ~/code/mighty_ws && ./src/mighty/launch/run_mighty_sim.sh ~/code/mighty_ws/install/setup.bash
    ```

### Notes

<details>
  <summary><b>Goal Command Example</b></summary>

  - ```bash
    ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 305.0, y: 0.0, z: 3.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
    ```
</details>
