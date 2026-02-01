# MIGHTY: Hermite Spline-based Efficient Trajectory Planning

If you like this project, please consider starring ⭐ the repo!

**Submitted to the IEEE Robotics and Automation Letters (RA-L)**

| **Trajectory** | **Forest** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_complex_benchmarks.gif" width="360" height="240" alt="Complex Benchmarks"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hard_forest.gif" width="360" height="240" alt="Static Forest"></a> |

| **Dynamic Obstacles** | **Long Flight** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_dynamic_sim.gif" width="360" height="240" alt="Dynamic Obstacles"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_long_flight.gif" width="360" height="240" alt="Hardware Long Flight"></a>

| **Fast Flight 1** | **Fast Flight 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_fast_flight_1.gif" width="360" height="240" alt="Hardware Fast Flight 1"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_fast_flight_2.gif" width="360" height="240" alt="Hardware Fast Flight 2"></a>

| **Dynamic Env 1** | **Dynamic Env 2** |
| ------------------------- | ------------------------- |
<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_dynamic_1.gif" width="360" height="240" alt="Hardware Dynamic Env 1"></a> | <a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_gifs_hw_dynamic_2.gif" width="360" height="240" alt="Hardware Dynamic Env 2"></a>

## Paper

MIGHTY: Hermite Spline-based Efficient Trajectory Planning is available at [https://arxiv.org/abs/2511.10822](https://arxiv.org/abs/2511.10822)!

```bibtex
@article{kondo2025mighty,
      title={MIGHTY: Hermite Spline-based Efficient Trajectory Planning},
      author={Kota Kondo and Yuwei Wu and Vijay Kumar and Jonathan P. How},
      year={2025},
      eprint={2511.10822},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      url={https://arxiv.org/abs/2511.10822},
}
```

## Video

The full video is available at [https://youtu.be/Pvb-VPUdLvg](https://youtu.be/Pvb-VPUdLvg).

## Interactive Demo

For an interactive demo of MIGHTY, please switch to the `interactive_demo` branch:
```bash
git checkout interactive_demo
```
and follow the setup instructions in the README of that branch.

**Note:** When forking this repository, unselect "Copy the main branch only" to include the `interactive_demo` branch.

---

## Installation

MIGHTY has been tested on Ubuntu 22.04 with ROS 2 Humble.

### Docker Installation (Recommended)

**1. Install Docker**

Follow the [official Docker installation guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/).

**2. Clone the Repository**

```bash
mkdir -p ~/code/ws/src
cd ~/code/ws/src
git clone https://github.com/mit-acl/mighty.git
cd mighty/docker
```

**3. Build the Docker Image**

```bash
make build

# Or build without cache (useful when dependencies change)
make build-no-cache
```

**4. Run Simulations**

```bash
# Multi-agent simulation (default: 10 agents)
make run

# Multi-agent with custom number of agents
make run-multiagent NUM_AGENTS=5

# Single-agent Gazebo simulation
make run-gazebo

# Gazebo with custom goal
make run-gazebo GOAL_X=100 GOAL_Y=50 GOAL_Z=3

# Gazebo with different environment (default: hard_forest)
make run-gazebo ENV=easy_forest

# Interactive shell (for debugging)
make shell
```

<details>
  <summary><b>Docker Make Targets Reference</b></summary>

  | Target | Description | Options |
  |--------|-------------|---------|
  | `make build` | Build the Docker image | - |
  | `make build-no-cache` | Build without cache (forces fresh build) | - |
  | `make run` | Run multiagent simulation (default: 10 agents) | - |
  | `make run-multiagent` | Run multiagent with custom agent count | `NUM_AGENTS=N` (default: 10) |
  | `make run-gazebo` | Run single-agent Gazebo simulation | `GOAL_X`, `GOAL_Y`, `GOAL_Z` (default: 305, 0, 3), `ENV` (default: hard_forest) |
  | `make shell` | Open interactive shell for debugging | - |

</details>

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

---

### Native Installation

**1. Clone the Repository**

```bash
mkdir -p ~/code
cd ~/code
git clone https://github.com/mit-acl/mighty.git mighty_ws/src/mighty
cd mighty_ws/src/mighty
```

**2. Run the Setup Script**

```bash
./setup.sh
```

This automated script will:
- ✅ Install ROS 2 Humble (if not already installed)
- ✅ Install all system dependencies
- ✅ Import all repositories from `mighty.repos` at tested commits
- ✅ Build DecompROS2, Livox-SDK2, and livox_ros_driver2
- ✅ Build MIGHTY and all ROS dependencies
- ✅ Configure your `~/.bashrc` for future use

**Notes:**
- You'll be prompted for sudo password once at the start
- Safe to re-run if something fails (skips already-installed components)
- After completion, run `source ~/.bashrc` to use MIGHTY immediately

**3. Run Simulations**

Use the unified simulation launcher script `run_sim.py`:

```bash
cd ~/code/mighty_ws

# Multi-agent simulation (default: 10 agents)
python3 src/mighty/scripts/run_sim.py --mode multiagent --setup-bash ~/code/mighty_ws/install/setup.bash

# Multi-agent with custom number of agents
python3 src/mighty/scripts/run_sim.py --mode multiagent -s ~/code/mighty_ws/install/setup.bash --num-agents 5

# Single-agent Gazebo simulation
python3 src/mighty/scripts/run_sim.py --mode gazebo -s ~/code/mighty_ws/install/setup.bash

# Gazebo with custom goal position
python3 src/mighty/scripts/run_sim.py --mode gazebo -s ~/code/mighty_ws/install/setup.bash --goal 100 0 3

# Gazebo with different environment
python3 src/mighty/scripts/run_sim.py --mode gazebo -s ~/code/mighty_ws/install/setup.bash --env easy_forest

# Gazebo with GUI enabled
python3 src/mighty/scripts/run_sim.py --mode gazebo -s ~/code/mighty_ws/install/setup.bash --gazebo-gui
```

<details>
  <summary><b>All Simulation Options</b></summary>

  ```
  --mode, -m          Required. 'multiagent' or 'gazebo'
  --setup-bash, -s    Required. Path to setup.bash
  --goal, -g          Goal position X Y Z for gazebo mode (default: 305 0 3)
  --start, -p         Start position X Y Z for gazebo mode (default: 0 0 3)
  --start-yaw         Start yaw in radians (default: 1.57)
  --num-agents, -n    Number of agents for multiagent mode (default: 10)
  --radius, -r        Circle radius for multiagent formation (default: 10)
  --env, -e           Gazebo environment (default: hard_forest)
  --ros-domain-id     ROS_DOMAIN_ID (default: 7)
  --no-rviz           Disable RViz
  --gazebo-gui        Enable Gazebo GUI
  --dry-run           Print generated YAML without launching
  ```
</details>

---

## Additional Information

<details>
  <summary><b>Dependencies</b></summary>

  All dependencies are version-controlled in `mighty.repos`:
  - **ROS 2 packages**: acl-mapping, dynus_interfaces, gazebo_ros_pkgs, livox_laser_simulation_ros2, realsense_gazebo_plugin, uav_simulator
  - **DecompROS2**: Decomposition utilities (requires decomp_util to build first)
  - **Livox-SDK2**: Livox LiDAR SDK (non-ROS binary)
  - **livox_ros_driver2**: Livox ROS 2 driver (uses custom build.sh)
</details>

<details>
  <summary><b>Paper Benchmarking</b></summary>

  The simple, complex, jerk weight sweep, and reference position/velocity benchmarking is available at:
  ```bash
  https://github.com/kotakondo/GCOPTER
  ```
</details>

<details>
  <summary><b>Bag Recording</b></summary>

  ```bash
  python3 src/mighty/scripts/bag_record.py --bag_number 3
  ```
</details>

<details>
  <summary><b>Goal Command Example</b></summary>

  ```bash
  ros2 topic pub /NX01/term_goal geometry_msgs/msg/PoseStamped "{header: {stamp: {sec: 0, nanosec: 0}, frame_id: 'map'}, pose: {position: {x: 305.0, y: 0.0, z: 3.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}}" --once
  ```
</details>
