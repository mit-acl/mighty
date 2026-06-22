# MIGHTY: Hermite Spline-based Efficient Trajectory Planning

If you like this project, please consider starring ⭐ the repo!

**Accepted to the IEEE Robotics and Automation Letters (RA-L)**

<video src="https://github.com/user-attachments/assets/a5127ce3-6662-4b5f-8ca6-2f84f38fddf8" width="100%" autoplay loop muted playsinline controls></video>

**Multi-Agent Trajectory Planning** — multiple aerial agents arranged on a circle swap to diametrically opposite positions, avoiding one another in real time:

<video src="https://github.com/user-attachments/assets/25ebadb2-050c-4762-81df-80c962d652d5" width="100%" autoplay loop muted playsinline controls></video>

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
@ARTICLE{kondo2026mighty,
  author={Kondo, Kota and Wu, Yuwei and Kumar, Vijay and How, Jonathan P.},
  journal={IEEE Robotics and Automation Letters}, 
  title={MIGHTY: Hermite Spline-based Efficient Trajectory Planning}, 
  year={2026},
  volume={},
  number={},
  pages={1-8},
  keywords={Anisotropic;Central Processing Unit;Filters;Radio access networks;Regional area networks;Location awareness;Mobile communication;Communication systems;High frequency;Indoor environment;Aerial Systems: Perception and Autonomy;Motion and Path Planning;Trajectory Optimization;Hermite Splines;Unmanned Aerial Vehicles},
  doi={10.1109/LRA.2026.3681187}
}
```

## Video

The full video is available at [https://youtu.be/Pvb-VPUdLvg](https://youtu.be/Pvb-VPUdLvg).

---

## Installation

MIGHTY has been tested on Ubuntu 22.04 with ROS 2 Humble. Four installation methods are available:

| Method | Platform | Notes |
|--------|----------|-------|
| [Docker (Linux)](#docker-installation-linux) | Linux | Uses Docker Engine (apt install, **not** Docker Desktop) |
| [Docker (Mac)](#docker-installation-mac) | macOS (Apple Silicon / Intel) | Uses Docker Desktop; Xpra for browser-based visualization |
| [Native (Linux)](#native-installation-linux) | Ubuntu 22.04 | Best for development and hardware deployment |
| [Jetson](#jetson-setup) | NVIDIA Jetson (Orin Nano, etc.) | ARM64 build; skips Gazebo, limits parallelism for low RAM |

---

### Docker Installation (Linux)

**1. Install Docker Engine**

Install Docker Engine via apt — do **not** use Docker Desktop on Linux as it may cause issues.

Follow the [official guide (Install using the apt repository)](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository).

**2. Clone the Repository**

```bash
git clone --depth 1 https://github.com/mit-acl/mighty.git
cd mighty/docker
```

**3. Build the Docker Image**

```bash
make build

# Or build without cache (useful when dependencies change)
make build-no-cache
```

**4. Run Simulations**

> **GPU errors?** On Linux, GPU support (NVIDIA) is enabled by default. If you get GPU-related errors (e.g., no NVIDIA driver or runtime), disable it with `GPU=false`:
> ```bash
> make run-interactive GPU=false
> ```

```bash
# Single-agent interactive simulation (click goals in RViz2 with "2D Goal Pose")
make run-interactive

# Multi-agent aerial simulation (default: 10 agents swapping positions on a circle)
make run-multiagent

# Multi-agent with a custom number of agents
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

In `run-interactive` mode, send goals from the RViz2 toolbar:
1. Click **"2D Goal Pose"**
2. Click and drag on the map to set the goal position and orientation
3. The agent will plan and navigate to the goal

<details>
  <summary><b>Docker Make Targets Reference</b></summary>

  | Target | Description | Options |
  |--------|-------------|---------|
  | `make build` | Build the Docker image | - |
  | `make build-no-cache` | Build without cache (forces fresh build) | - |
  | `make run-interactive` | Run single agent with manual goal (RViz2 2D Goal Pose) | - |
  | `make run-multiagent` | Run multi-agent aerial simulation (agents swap positions on a circle) | `NUM_AGENTS` (default: 10) |
  | `make run-gazebo` | Run single-agent Gazebo simulation | `GOAL_X`, `GOAL_Y`, `GOAL_Z` (default: 305, 0, 3), `ENV` (default: hard_forest) |
  | `make run-mac` | Run multi-agent aerial simulation on Mac (Xpra, browser at localhost:8080) | `NUM_AGENTS` (default: 10) |
  | `make run-mac-interactive` | Run single agent on Mac with manual goal (Xpra, browser at localhost:8080) | - |
  | `make run-mac-gazebo` | Run Gazebo on Mac (Xpra) | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` |
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

### Docker Installation (Mac)

MIGHTY runs on macOS via Docker with [Xpra](https://xpra.org/) for browser-based visualization. Xpra is installed inside the Docker image automatically — you do **not** need to install Xpra, X11, or XQuartz on your Mac.

#### Prerequisites

**1. Install Docker Desktop**

Download and install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/).

**2. Configure Docker Desktop**

Open Docker Desktop and go to **Settings** (gear icon). Two settings must be changed:

- **General** > Under "Virtual Machine Options", enable **"Use Rosetta for x86_64/amd64 emulation on Apple Silicon"** (Apple Silicon Macs only — significantly speeds up the amd64 build and runtime)
- **Resources** > Set **Memory** to at least **24 GB** (32 GB recommended). The default 4 GB is not enough and will cause out-of-memory failures during the build.

#### Building

```bash
git clone --depth 1 https://github.com/mit-acl/mighty.git
cd mighty/docker
make build
```

> The first build takes a while since it cross-compiles for amd64. Subsequent builds use Docker layer caching and are much faster. Use `make build-no-cache` to force a fresh build.

#### Running Simulations

All `run-mac*` targets use Xpra to stream GUI windows (RViz2, etc.) to your browser.

```bash
cd mighty/docker

# Multi-agent aerial simulation (default: 10 agents swapping positions on a circle)
make run-mac

# Multi-agent with a custom number of agents
make run-mac NUM_AGENTS=5

# Single agent with manual goal control (use RViz2's "2D Goal Pose" tool)
make run-mac-interactive

# Single-agent Gazebo simulation
make run-mac-gazebo
```

#### Visualization

After running any `run-mac*` command, open your browser and go to:

**http://localhost:8080**

RViz2 and other GUI windows will appear in the browser. You can interact with them just like native windows — pan, zoom, and use the RViz2 toolbar.

To send a goal manually (in `run-mac-interactive` mode):
1. In the RViz2 toolbar, click **"2D Goal Pose"**
2. Click and drag on the map to set the goal position and orientation
3. The agent will plan and navigate to the goal

<details>
  <summary><b>Mac Make Targets Reference</b></summary>

  | Target | Description | Options |
  |--------|-------------|---------|
  | `make run-mac` | Multi-agent aerial simulation (agents swap positions on a circle) | `NUM_AGENTS` (default: 10) |
  | `make run-mac-interactive` | Single agent with manual 2D Goal Pose | - |
  | `make run-mac-gazebo` | Single-agent Gazebo simulation | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` |

</details>

<details>
  <summary><b>Troubleshooting</b></summary>

  - **Build fails with out-of-memory error**: Increase Docker Desktop memory allocation (Settings > Resources > Memory). 24 GB minimum, 32 GB recommended.
  - **Build is very slow**: Make sure Rosetta emulation is enabled in Docker Desktop settings (General > Virtual Machine Options).
  - **Browser shows nothing at localhost:8080**: Wait 10-20 seconds after launching — Xpra takes a moment to start. Check the terminal for `[Docker] Xpra server running on port 8080`.
  - **Port 8080 already in use**: Another service is using port 8080. Stop it first, or modify the `-p 8080:8080` in the Makefile to use a different port (e.g., `-p 9090:8080`).

</details>

---

### Native Installation (Linux)

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
- Install ROS 2 Humble (if not already installed)
- Install all system dependencies
- Import all repositories from `mighty.repos` at tested commits
- Build DecompROS2, Livox-SDK2, and livox_ros_driver2
- Build MIGHTY and all ROS dependencies
- Configure your `~/.bashrc` for future use

**Notes:**
- You'll be prompted for sudo password once at the start
- Safe to re-run if something fails (skips already-installed components)
- The script appends ROS 2 sourcing, workspace sourcing, and environment variables to `~/.bashrc`
- After completion, run `source ~/.bashrc` to use MIGHTY immediately

**3. Run Simulations**

Use the unified simulation launcher script `run_sim.py`:

```bash
cd ~/code/mighty_ws

# Single-agent interactive simulation (click goals in RViz2 with "2D Goal Pose")
python3 src/mighty/scripts/run_sim.py --mode interactive --setup-bash ~/code/mighty_ws/install/setup.bash

# Multi-agent aerial simulation (default: 10 agents swapping positions on a circle)
python3 src/mighty/scripts/run_sim.py --mode multiagent -s ~/code/mighty_ws/install/setup.bash

# Multi-agent with a custom number of agents
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
  --mode, -m          Required. 'interactive', 'multiagent', or 'gazebo'
  --setup-bash, -s    Required. Path to setup.bash
  --num-agents, -n    Number of agents for multiagent mode (default: 10)
  --radius            Circle radius for the multiagent formation (default: 10.0)
  --goal, -g          Goal position X Y Z for gazebo mode (default: 305 0 3)
  --start, -p         Start position X Y Z for gazebo mode (default: 0 0 3)
  --start-yaw         Start yaw in radians (default: 1.57)
  --env, -e           Gazebo environment (default: hard_forest)
  --ros-domain-id     ROS_DOMAIN_ID (default: 20)
  --no-rviz           Disable RViz
  --gazebo-gui        Enable Gazebo GUI
  --dry-run           Print generated YAML without launching
  ```
</details>

---

### Jetson Setup

For NVIDIA Jetson boards (Orin Nano, Orin NX, AGX Orin, etc.) running JetPack with Ubuntu 22.04.

Gazebo Classic is not available on ARM64, so the setup script automatically detects the Jetson platform and skips simulation-only packages. The planner and all hardware nodes build normally.

**1. Add Swap Space (Recommended)**

Jetson boards have limited RAM (e.g., 8 GB on Orin Nano). Compilation of large C++ files can trigger the OOM killer without extra swap:

```bash
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Make it persistent across reboots
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

**2. Clone and Run Setup**

```bash
mkdir -p ~/code
cd ~/code
git clone https://github.com/mit-acl/mighty.git mighty_ws/src/mighty
cd mighty_ws/src/mighty
./setup.sh
```

The script auto-detects the Jetson via `/etc/nv_tegra_release` and will:
- Skip Gazebo packages (`gazebo_dev`, `gazebo_ros`, `gazebo_plugins`, `realsense_gazebo_plugin`, `ros2_livox_simulation`)
- Build with `-DBUILD_SIMULATION=OFF` (disables Gazebo-linked targets in the `mighty` package)
- Limit compiler parallelism (`-j1`) to avoid OOM on memory-constrained boards
- Use the correct `aarch64-linux-gnu` library path

You can also pass `--jetson` explicitly or combine with `-j`:

```bash
./setup.sh --jetson -j 2  # if you have enough RAM for 2 parallel compile jobs
```

**3. Test the Planner (without Gazebo)**

`fake_sim` is always built, even on Jetson. It provides simulated odometry and a lightweight environment without Gazebo so you can verify the planner works:

```bash
source ~/.bashrc
cd ~/code/mighty_ws
python3 src/mighty/scripts/run_sim.py --mode interactive --setup-bash ~/code/mighty_ws/install/setup.bash
```

Use RViz2's **"2D Goal Pose"** tool to send goals and confirm the planner is running correctly.

<details>
  <summary><b>What is skipped on Jetson</b></summary>

  The following packages are excluded from the build:

  | Package | Reason |
  |---------|--------|
  | `gazebo_dev` | Requires `libgazebo11-dev` (x86 only) |
  | `gazebo_ros` | Depends on `gazebo_dev` |
  | `gazebo_plugins` | Depends on `gazebo_dev` |
  | `gazebo_ros_pkgs` | Meta-package for all Gazebo ROS packages |
  | `realsense_gazebo_plugin` | Gazebo plugin for RealSense simulation |
  | `ros2_livox_simulation` | Gazebo plugin for Livox LiDAR simulation |

  Within the `mighty` package, these simulation-only targets are also disabled:
  - `move_model` (Gazebo model mover plugin)
  - `imu_plugin` (Gazebo IMU sensor plugin)
  - `convert_velodyne_to_ros_time` (Gazebo timestamp converter)
  - `dynamic_forest_node` (simulated dynamic obstacles)

  All hardware nodes (`mighty`, `fake_sim`, `convert_odom_to_state`, `convert_vicon_to_state`) build normally.

</details>

---

## Acknowledgments

The L-BFGS solver implementation in this repository (`src/mighty/lbfgs_solver.cpp`, `include/mighty/lbfgs_solver.hpp`) is adapted from [ZJU-FAST-Lab/GCOPTER](https://github.com/ZJU-FAST-Lab/GCOPTER/blob/main/gcopter/include/gcopter/lbfgs.hpp). We thank the authors for making their implementation publicly available.

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
