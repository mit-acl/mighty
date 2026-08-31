# MIGHTY: Hermite Spline-based Efficient Trajectory Planning

If you like this project, please consider starring ⭐ the repo!

**Accepted to the IEEE Robotics and Automation Letters (RA-L)**

<video src="https://github.com/user-attachments/assets/a5127ce3-6662-4b5f-8ca6-2f84f38fddf8" width="100%" autoplay loop muted playsinline controls></video>

**Multi-Agent Trajectory Planning** — multiple aerial agents arranged on a circle swap to diametrically opposite positions, avoiding one another in real time:

<video src="https://github.com/user-attachments/assets/25ebadb2-050c-4762-81df-80c962d652d5" width="100%" autoplay loop muted playsinline controls></video>

**Ground Robot Deployment** — a single Pioneer 3-AT on hardware and in simulation, both autonomously exploring an unknown space and navigating to a commanded goal:

<video src="https://github.com/user-attachments/assets/15bed70e-f922-44c7-9df2-51693ad33bb7" width="100%" autoplay loop muted playsinline controls></video>

**Multi-Robot Ground Exploration** — three Pioneer 3-ATs explore an unknown office together. Each drives to the frontier it is best placed for — the one with the fewest teammates closer to it — and shares its map as it goes, so they divide the space instead of re-covering each other's ground:

<video src="https://github.com/user-attachments/assets/REPLACE-WITH-ASSET-ID" width="100%" autoplay loop muted playsinline controls></video>

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

## Running the Simulations

Every scenario, and the one command that starts it. Native commands assume the
workspace is built and run from `~/code/mighty_ws`; Docker commands run from
`mighty/docker`. See [Installation](#installation) first.

| Scenario | Native | Docker |
|----------|--------|--------|
| **Ground robot, autonomous exploration** (Pioneer 3-AT explores `ACL_office` and returns to its start) | `python3 src/mighty/scripts/run_sim.py --mode exploration-singleagent-ground -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-exploration` |
| **Ground robots, multi-robot exploration** (3 Pioneer 3-ATs divide `ACL_office` via MinPos) | `python3 src/mighty/scripts/run_sim.py --mode exploration-multiagent-ground -s ~/code/mighty_ws/install/setup.bash` | `make run-multiagent-ground-exploration` |
| Ground robot in a forest world | `python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-robot` |
| Single UAV, goals clicked in RViz | `python3 src/mighty/scripts/run_sim.py --mode interactive -s ~/code/mighty_ws/install/setup.bash` | `make run-interactive` |
| Multi-agent UAVs swapping on a circle | `python3 src/mighty/scripts/run_sim.py --mode multiagent -s ~/code/mighty_ws/install/setup.bash` | `make run-multiagent` |
| **Ground robots swapping on a circle, no Gazebo** (10 Pioneer 3-ATs exchange antipodal positions in RViz; full planner→MPC→unicycle chain; ≥ 1 m separation floor) | `python3 src/mighty/scripts/run_sim.py --mode multiagent-ground-fake -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-swap-fake` |
| **Ground robots swapping on a circle, Gazebo + sensorless models** (same exchange, real P3AT models rendered in RViz; `--goal-pattern random_slot` for continuous traffic) | `python3 src/mighty/scripts/run_sim.py --mode multiagent-ground-gazebo -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-swap` |
| Single UAV through a Gazebo forest | `python3 src/mighty/scripts/run_sim.py --mode gazebo -s ~/code/mighty_ws/install/setup.bash` | `make run-gazebo` |

Useful flags (native): `--no-rviz` for headless, `--env <world>` to pick a world,
`--num-agents N`, `--goal X Y Z`, `--gazebo-gui`, and `--log-level info` to see
the planner's own diagnostics. `python3 src/mighty/scripts/run_sim.py --help`
lists them all, and the "All Simulation Options" block under
[Native Installation](#native-installation-linux) documents each one. Docker
equivalents are environment variables — `ENV=`, `NUM_AGENTS=`, `GOAL_X/Y/Z=`,
`NO_RVIZ=true`, `GPU=false`.

> **Ground robot users:** the ground-robot scenarios need the
> [`mpc`](https://github.com/kotakondo/mpc) path-tracking controller. It is
> listed in `mighty.repos`, so `setup.sh` and the Docker build fetch it
> automatically — but if you assembled the workspace by hand and see
> `PackageNotFoundError: mpc`, that is the missing piece. UAV scenarios do not
> need it.

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

`make build` clones MIGHTY at the pinned release tag (`MIGHTY_VERSION` in
`docker/Dockerfile`) and imports every dependency at the commit pinned in
`mighty.repos`, so the image is reproducible and independent of your working
tree. Override the tag with `make build MIGHTY_VERSION=v0.0.7`.

<details>
  <summary><b>Testing local changes in Docker</b></summary>

  Two mechanisms, and the difference matters:

  | | What it picks up | Cost |
  |---|---|---|
  | `DEV=1` on any `run-*` target | Bind-mounts your local `config/`, `launch/`, `rviz/` and `src/` over the image's. **C++ changes are NOT picked up** — the binaries in the image were compiled from the tag. | free |
  | `make build-local` | Rebuilds the image from your **working tree**, recompiling the mighty package. | minutes |

  ```bash
  # Iterate on parameters / launch files / RViz configs — no rebuild
  make run-ground-exploration DEV=1

  # Test uncommitted C++ against the image
  make build-local
  make run-ground-exploration

# Multi-robot: 3 ground robots dividing the map via MinPos
make run-multiagent-ground-exploration
  ```

  `build-local` reuses every cached dependency layer and recompiles only the
  mighty package, so it takes minutes rather than the hours a full `make build`
  needs. It requires a prior `make build` to populate those layers. Note that
  `mighty.repos` still comes from the pinned tag: if your change needs a *new
  dependency version*, bump the pin, push a tag, and use
  `make build MIGHTY_VERSION=<tag>` instead.
</details>

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
  | `make run-gazebo` | Run single-agent Gazebo simulation | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` — all optional; unset means the same defaults as the native command (hard_forest, goal 305 0 3) |
  | `make run-ground-swap` | 10 ground robots exchanging positions in Gazebo (sensorless P3AT models, ≥ 1 m separation floor) | `NUM_AGENTS`, `GOAL_PATTERN` (`random_slot`), `GOAL_STAGGER`, `RADIUS`, `MPC_RATE`, `MPC_HORIZON` — all optional, native defaults when unset |
  | `make run-ground-swap-fake` | Same exchange without Gazebo (colored boxes in RViz — lightest variant) | same as `run-ground-swap` |
  | `make run-mac` | Run multi-agent aerial simulation on Mac (Xpra, browser at localhost:8080) | `NUM_AGENTS` (default: 10) |
  | `make run-mac-ground-swap` | Ground swap on Mac (Xpra, browser at localhost:8080) | `NUM_AGENTS`, `GOAL_PATTERN`, `GOAL_STAGGER`, `RADIUS` |
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
- Install the Python packages the ground-robot controller needs
  (`casadi`, `do-mpc`) — see the note below
- Import all repositories from `mighty.repos` at tested commits
- Build DecompROS2, Livox-SDK2, and livox_ros_driver2
- Build MIGHTY and all ROS dependencies
- Configure your `~/.bashrc` for future use

> **If you installed before this note existed**, or you set the workspace up by
> hand, install the two Python dependencies of the `mpc` ground-robot
> controller explicitly:
>
> ```bash
> pip install casadi==3.6.7 do-mpc==5.1.1
> ```
>
> They are imported by `mpc_node.py` but not declared in `mpc/package.xml`
> (which lists only ROS packages), so `rosdep` will not install them and the
> build still succeeds — `mpc` is pure Python and nothing is compiled. Without
> them the ground robot simply never moves: `mpc_node` exits immediately with
> `ModuleNotFoundError`, while every other node starts normally, so the
> terminal looks healthy. If your robot sits motionless in the exploration sim,
> check this first.

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
  --mode, -m          Required. One of:
                        interactive                     single UAV, goals clicked in RViz
                        multiagent                      N UAVs swapping on a circle (fake_sim)
                        gazebo                          single UAV (or --ground-robot) in Gazebo
                        exploration-singleagent-ground  ground robot, autonomous exploration
                        exploration-multiagent-ground   multiple ground robots, MinPos exploration
                        multiagent-ground               N ground robots with MPC
                        swap-multiagent-ground          ground robots swapping positions
                        dyn-test, dyn-test-ground, dyn-test-ground-mpc
                                                        dynamic-obstacle test scenarios
  --setup-bash, -s    Required. Path to setup.bash
  --num-agents, -n    Number of agents for multiagent modes (default: 10)
  --radius, -r        Circle radius for the multiagent formation (default: 10.0)
  --goal, -g          Goal position X Y Z for gazebo mode (default: depends on --env —
                      305 for hard_forest, 105 for easy_forest/forest, 30 for the
                      small worlds; z drops to 0 for --ground-robot)
  --start, -p         Start position X Y Z for gazebo mode (default: 0 0 3)
  --start-yaw         Start yaw in radians (default: 1.57)
  --env, -e           Gazebo environment (default: hard_forest)
  --ros-domain-id     ROS_DOMAIN_ID (default: 30)
  --ground-robot      Use the Pioneer 3-AT ground robot instead of the UAV
  --no-rviz           Disable RViz (headless; required for unattended runs)
  --gazebo-gui        Enable Gazebo GUI
  --no-goal           Do not auto-publish a terminal goal
  --use-vlm           Use the VLM frontier selector (exploration-singleagent-ground only)
  --follow            Also spawn the YOLOv8 person tracker (needs --use-vlm)
  --dry-run           Print generated YAML without launching
  --emit-yaml PATH    Write the generated YAML to PATH and exit (used by the test harness)
  ```

  > The goal default is chosen from the world's actual extent, so the agent
  > crosses it rather than stopping partway. `hard_forest` runs to x=301.6
  > (hence 305, matching the Docker entrypoint), while `easy_forest` ends at
  > x=100. Pass `--goal X Y Z` to override.
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

## Ground Robot Simulation

Besides the aerial vehicle, MIGHTY drives wheeled ground robots (Pioneer 3-AT)
in Gazebo, alone or as a coordinated team:

| Scenario | What it does | Goal |
|----------|--------------|------|
| **Autonomous exploration** | One robot explores an unknown `ACL_office` map using frontier detection. | None — self-driven |
| **Multi-robot exploration** | Three robots explore `ACL_office` together, dividing frontiers between them via MinPos. | None — self-driven |
| **Forest navigation** | One robot crosses a forest world (e.g. `hard_forest`) to a commanded position. | `--goal X Y Z` |
| **Position exchange, no Gazebo** | Ten robots on a 10 m circle continuously swap antipodal positions in an obstacle-free RViz world (`--mode multiagent-ground-fake`). | None — goal_monitor drives |

> **How the no-Gazebo mode works.** Each robot runs the full real chain —
> planner → `mpc_node` → `fake_sim`'s unicycle integrator — with RViz as the
> only display (per-robot colored boxes, planned and driven trajectories).
> There is no mapper and no obstacle world; the planner is fed by a constant
> *sub-floor* point grid at z = −3 m (`scripts/publish_subfloor_cloud.py`),
> which keeps its map machinery alive while staying below `z_min`, so the
> planned world is genuinely empty. All goals are released together by
> default — the classic simultaneous crossing (`--goal-stagger 5` releases
> them sequentially instead), and `config/swap_mighty_ground_robot.yaml`
> documents every parameter choice with the measurements behind it.
> Both swap modes enforce a **guaranteed ≥ 1 m pairwise separation floor**
> — see the Gazebo variant below for the mechanism and the measured
> numbers; `scripts/separation_probe.py` is the verifier.

> **The Gazebo variant** (`--mode multiagent-ground-gazebo`) runs the same
> exchange with Gazebo in the loop, so RViz shows the actual **P3AT robot
> model** (RobotModel displays fed by each robot's `robot_description`)
> instead of colored boxes. What makes 10 Gazebo robots affordable:
> `use_lidar:=false use_camera:=false` strip every sensor from the spawned
> model — the mid360 alone costs 36,000 rays at 10 Hz per robot — and
> `map_source:=sensor_cloud` feeds each planner from the same sub-floor
> cloud as the no-Gazebo mode. Measured with all ten robots driving:
> **gzserver ~10% of one core** (three lidar-equipped robots cost
> 150–200%). Motion remains the teleport chain every Gazebo ground mode
> uses — fake_sim integrates `cmd_vel` and teleports the model — so robots
> are kinematic puppets. The harness scenario (`ground-swap-10-gazebo`)
> additionally gates on every namespace existing as a Gazebo entity, closing
> the documented spawn-flake trap where a robot drives around invisibly.
>
> **The ≥ 1 m separation floor — how it works.** Soft avoidance cannot
> floor a ten-way crossing: with all path-shaping costs healthy, minima of
> 0.02–0.27 m were measured at the shared centre. The guarantee is
> geometric instead, two mechanisms working as a pair
> (`peer_standoff_m` / `peer_static_disc_m`, enabled only by the swap
> overlay — every other config keeps them off):
>
> 1. **Standoff rings** (reference truncation, `publishMpcPath`): the
>    junior of each ID pair cuts its MPC reference at the first waypoint
>    whose time-matched distance to the senior's shared trajectory *or* to
>    the senior's live position would dip below 1.2 m (checking both worlds
>    matters: a yielding robot's advertised plan diverges from where it
>    actually is). The senior symmetrically stands off moving juniors at a
>    slightly smaller ring — in any head-on the junior stops first, which
>    preserves the priority order, and a robot already inside a ring may
>    only move *outward*, so near-misses cannot lock in.
> 2. **Parked robots become real obstacles** (`occupancyMapCallback`): a
>    peer whose live position sits still for ~0.25 s gets a *filled* disc
>    of synthetic obstacle points stamped into every other robot's sensor
>    cloud at the corridor-true radius (disc + robot radius), so A*, the
>    safe-corridor decomposition, and the heat layer all wall it off and
>    every replanned detour is born corridor-feasible. Detours therefore
>    clear the rings, which is what lets every hold self-release — no
>    timeouts anywhere. (Nearby peers get a body-size disc instead, so a
>    robot's own A* start is never swallowed by a wall.)
>
> Validated on a 48-core host (40 Hz MPC): on the shipping build, the
> antipodal campaign measured run minima of **1.183 / 1.215 m with zero
> time below 1 m** (three earlier 900 s runs on the pre-hardening build:
> 1.155–1.187 m, also zero), all ten robots completing legs.
> **`random_slot` is hardened, not guaranteed**: its continuous goal churn
> creates crawl-approach geometries the rings cover less perfectly, and
> the final 3 × 900 s campaign measured **zero or one brief sub-1 m
> episode per run** (3.4–14.3 s, minima 0.08–0.53 m) — down from 16–21
> episodes per run before the hardening. The residual event class and the
> next lever (reserved crossings / speed-scaled rings) are documented
> future work. The guarantee has an honest
> price: every close encounter now costs a park → wall → detour cycle, so
> fleet throughput at ten simultaneous robots is roughly **2–3× slower**
> than the old floorless free-for-all, and lower-priority (higher-ID)
> robots yield more often — that asymmetry is what makes the pattern
> deadlock-free. Verify any run yourself with
> `scripts/separation_probe.py --agents 10 --threshold 1.0`.
>
> **Two goal patterns.** The default is the classic simultaneous antipodal
> exchange — maximum drama, one big crossing wave. `--goal-pattern
> random_slot` makes each robot draw its next goal from the ten
> evenly-spaced circle slots on every arrival instead: traffic spreads
> organically in space and time, and fleet throughput is far higher than
> the antipodal wave because encounters are pairwise instead of ten-way
> (measured 129–160 slot arrivals per 900 s run with every robot
> progressing; see the separation note above for its floor status). Seeded per-robot (stable across runs; `random_seed` on
> `goal_monitor.launch.py`), so runs are reproducible. Two mechanics worth
> knowing: slots double as parking spots, so a goal drawn onto a
> momentarily-occupied slot is rejected by the planner ("goal in occupied
> space") — the monitor self-heals by **redrawing after 20 s without
> progress** (`redraw_stagnation_sec`; redraws are normal traffic
> management, not failure). Without that redraw, a rejected robot parks
> and blocks *its* slot for everyone — a measured deadlock chain (one
> robot pinned 800 s, reproducibly, before the fix). That waiting-out-a-
> redraw state is also why the planner keeps broadcasting on `/trajs` even
> when replanning fails: receivers prune agents silent for
> `traj_lifetime`, and a pruned robot is invisible to every peer's
> avoidance — the keepalive re-shares the ended trajectory, whose
> evaluation is the parked position.
>
> **Where visual "wobble" comes from — measured, not guessed.** A robot
> driving alone in this mode tracks perfectly: **zero yaw-rate reversals**
> over full crossings, 2 cm mean cross-track. Under dense traffic the
> residual weave is the avoidance system working at saturation, not a
> tracking error — an 11-cell tuning matrix over every path-shaping knob
> moved nothing outside noise. With the separation floor the worst churn is
> replaced by queue-and-detour behaviour by construction; for even calmer
> traffic use `--num-agents 5`, `--goal-stagger 2`, or a larger `--radius`.
>
> Standing up ten robots exposed and fixed three latent fleet-scale bugs
> that benefit every mode: (1) `fake_sim`'s Gazebo teleport spawned a
> blocking thread per 10 ms tick — once gzserver fell behind, the pile-up
> aborted the node; it is now a single in-flight async request that degrades
> gracefully. (2) The planner's replan attempt rate was hardcoded at 100 Hz;
> ten planners consumed ~8.5 cores and **starved the MPC control loops**
> (30 Hz design degraded to 16 Hz; the yaw command sat at its limit 33 % of
> driving time vs 1.3 % solo — visible as constant overshoot). The fix is a
> CPU budget: `replan_rate_hz: 20` in the swap overlay plus
> `mpc_control_rate_hz:=20 mpc_n_horizon:=10` launch overrides (exposed as
> `--mpc-rate` / `--mpc-horizon`; measured after: yaw saturation 2.9 %).
> (3) The planner's default `MultiThreadedExecutor` sizes itself to the
> machine's core count, and every callback thread lazily grows its own
> OpenMP pool — on a 48-core host that is **2,315 threads per planner**
> (25k threads and load 64 for the fleet). The executor is now capped at 8
> callback threads; on many-core hosts also set `OMP_NUM_THREADS=2
> OMP_WAIT_POLICY=PASSIVE` per robot. Every default users see is otherwise
> unchanged — solo and small-fleet modes keep 100 Hz replanning and the
> stock MPC settings.
>
> **Running the fleet on a big remote machine, watching locally.** The
> Gazebo swap is sensorless, so its DDS traffic (TF, markers, small clouds)
> is light enough to stream across a LAN: run the sim headless on the
> many-core box (`--no-rviz`, in Docker with `--network=host`) and start
> only RViz on your desk with the same `ROS_DOMAIN_ID`. On a 48-core host
> the whole 10-robot fleet runs at ~30 % load with MPC at
> `--mpc-rate 40 --mpc-horizon 20`. Two gotchas: pick a **non-default
> domain ID** so the sim cannot cross-talk with anything else on the LAN,
> and if you probe topics on the sim host itself, do it *inside* the sim's
> container — FastDDS's shared-memory transport silently drops data between
> sibling containers.

> **Exploring vs. driving to a goal.** The ground-robot scenarios differ only in
> who picks the goal:
>
> - `--mode exploration-singleagent-ground` / `make run-ground-exploration` —
>   the frontier loop drives, no goal needed.
> - `--mode exploration-multiagent-ground` / `make run-multiagent-ground-exploration` —
>   same, with three robots splitting the work (see Multi-Robot Exploration below).
> - `--mode gazebo --ground-robot` / `make run-ground-robot` — crosses the world
>   to the `--goal` / `GOAL_*` you give it, exactly like the UAV Gazebo
>   scenario. Exploration is switched off at launch
>   (`exploration_enabled:=false`), so it is honoured even though
>   `config/mighty_ground_robot.yaml` enables exploration.
>
> Both accept `--env` / `ENV`, so either behaviour is available in either world.
> The goal's z is forced to 0 for a ground robot when left at the UAV default.

### What a healthy exploration run looks like

The robot spawns at `(0, 2)` in `ACL_office`, explores until no reachable
frontiers remain, then drives back to its spawn point and parks.

| Stage | What you should see |
|-------|---------------------|
| ~10–30 s after launch | `Exploration start captured at (0, 2, 0)` in the planner pane; the robot starts moving |
| While exploring | Grey/blue frontier spheres with `id=N` labels, a yellow `Exploration Area` rectangle, and a yellow line from the robot to its current target |
| Finish | `No frontiers left. Robot now at (...). Returning to captured start (0, 2, 0)` followed by `MPC: goal reached` |

A full run maps ~99% of the exploration bounds (~1000 m² of free space). Expect
**9–20 minutes** — the robot uses greedy nearest-frontier selection, so it
sometimes finishes a dead end and drives all the way back across the map. That
is inefficient but normal, not a stall.

### Troubleshooting

**The robot never moves and `mpc_node` is missing from the node list.** Check
the launch terminal for:

```
[mpc_node-5] ModuleNotFoundError: No module named 'casadi'
[ERROR] [mpc_node-5]: process has died [exit code 1]
```

The `mpc` controller imports `casadi` and `do_mpc`, which are not declared in
its `package.xml`, so neither `rosdep` nor the build installs them. Everything
else starts cleanly, so the only symptom is a stationary robot. Fix with:

```bash
pip install casadi==3.6.7 do-mpc==5.1.1
```

`setup.sh` and the Docker image now install these automatically; this affects
workspaces created before that change.

**No frontier markers appear and the robot never moves.** The frontier detector
flood-fills outward from the robot over *known-free* cells. Sensors have a
near-field blind zone, so the mapper's small cleared cylinder around the robot
can be ringed by UNKNOWN and disconnected from everything the robot can see —
the flood fill then finds a single "frontier" on top of the robot, the dwell
check retires it within a second, and exploration ends before it starts.
`exploration.detector.unknown_bridge_cells` (default `4` for the ground robot)
lets the search step across that many UNKNOWN cells. Raise it if your sensor's
blind zone is wider; OCCUPIED cells always block, so it never bridges walls.

**The robot explores but never returns home.** Return-home is triggered by
frontier exhaustion. If `exploration.bounds.*` is larger than the reachable
space, unreachable frontiers keep it exploring — shrink the bounds to the area
you actually want covered.

**Seeing the planner's diagnostics.** `onboard_mighty.launch.py` runs
`mighty_node` at `--log-level error` by default, which keeps the terminal quiet
but also hides the planner's own INFO output. The most useful line for anything
frontier-related is throttled to 1 Hz and looks like:

```
[expl] grid 133x133  free=6521  occ=209  unknown=10959  fresh=12  db=27
```

`fresh` is how many frontier clusters the detector just found and `db` is how
many the persistent manager is tracking — `fresh=0 db=0` means detection is
finding nothing, which is a different problem from selection or navigation.
Enable it with:

```bash
ros2 launch mighty onboard_mighty.launch.py log_level:=info ...
```

Note that `std::cout` output (`No frontiers left`, `Changing DroneStatus`) is
never suppressed, so its presence does not tell you the log level is high enough.

### Multi-Robot Exploration

```bash
python3 src/mighty/scripts/run_sim.py --mode exploration-multiagent-ground \
  -s ~/code/mighty_ws/install/setup.bash            # native, 3 robots
make run-multiagent-ground-exploration               # Docker
make run-multiagent-ground-exploration NUM_AGENTS=4  # more robots
```

Three robots spawn on a line at `y = 2` (x = −5, 0, +5) and explore
`ACL_office` together. Each returns to **its own** spawn point when the team
runs out of frontiers.

**How they divide the work.** Every robot broadcasts its pose, and each frontier
is ranked by how many teammates are closer to it than you are. You take the
frontier where your rank is lowest — so the robot best placed for a given area
goes there, and the others look elsewhere. No negotiation, no central planner,
no assigned sectors: the allocation falls out of everyone applying the same rule
to the same information. Robots also exchange their persistent maps, so area one
robot has already seen stops generating frontiers for the others.

Positions alone still leave a race: two robots near-equidistant from a frontier
both compute rank 0, both commit, and both drive there. So intent and verdicts
travel too. Each robot **claims** the frontier it is pursuing; selection skips
anything inside a teammate's claim radius, and when two robots claim the same
spot before hearing each other, the lexicographically smaller namespace keeps it
and the other re-targets within about a second. Each robot also broadcasts its
frontier **status** list — including *visited* and *unreachable* verdicts — so
one robot finishing (or failing) a frontier retires it for the whole team in
under half a second, instead of a teammate driving across the office to a spot
that was already cleared.

Five global topics carry this (all agents publish and subscribe):

| Topic | Purpose |
|-------|---------|
| `/exploration/peer_poses` | pose broadcast that MinPos ranks against (5 Hz per robot) |
| `/exploration/peer_claims` | currently-pursued frontier; contested claims resolve by namespace (5 Hz while pursuing) |
| `/exploration/frontier_status` | full frontier DB incl. VISITED/INVALIDATED — one robot's verdict retires the record everywhere (2 Hz per robot) |
| `/exploration/visited_maps` | persistent map exchange (1 Hz per robot) |
| `/exploration/return_home` | operator command: everyone stop and go home |

RViz reads one more: `/exploration/frontier_markers`, the **single** team-wide
frontier layer (state-colored centroids plus an orange claim ring per robot),
published by the lexicographically smallest live namespace so three robots
don't paint three overlapping copies of the same frontiers.

`/trajs` deconfliction — the same mechanism the multi-UAV scenario uses — is
already active, which is why the robots avoid each other without any of the
above.

**Configuration.** All robots load `config/mighty_ground_robot.yaml` with
`config/multi_mighty_ground_robot.yaml` layered on top. That overlay is
deliberately tiny:

| Key | Solo | Team | Why |
|-----|------|------|-----|
| `exploration.minpos.enabled` | `false` | `true` | creates the peer-pose and map-sharing machinery |
| `exploration.select_nearest` | `true` | `false` | selection is `if (select_nearest) … else if (use_minpos) …`, so leaving this true makes MinPos unreachable |
| `exploration.minpos.min_frontier_dist_to_peers_m` | `3.0` | `4.0` | keeps the robots' work areas apart — the upstream reason they stop meeting in corridors |
| `exploration.minpos.claim_block_radius_m` | `4.0` | `4.0` | never target a spot a teammate is already driving to; the claim filter is never relaxed, unlike the keep-out above |
| `exploration.minpos.share_frontier_status` | `true` | `true` | visited/unreachable verdicts propagate team-wide (inert without `minpos.enabled`) |
| `exploration.visualization.centralized` | `false` | `true` | one global frontier-marker topic instead of per-namespace ones |
| `planner_Cw` | `1.0` | `1.4` | avoidance envelope; 1.4 + 0.6 m of robot = the ~2 m corridors exactly, so passing stays geometrically feasible (2.0 shoved trajectories into walls; 1.0 left no margin — see tuning note) |
| `dynamic_weight` | `1e+3` | `1e+3` | measured: raising it bought nothing and doubled the violence of close passes |
| `peer_yield_radius_m` / `peer_yield_speed_factor` | `0` (off) | `3.0` / `0.3` | head-on yield rule: tie-break loser crawls, winner proceeds — the asymmetry that resolves pinch encounters |
| `agent_dyn_horizon_sec` | `0` (off) | `8.0` | only react to peers' predicted trajectories 8 s out — see the smoothness note below |
| `jerk_weight` | `1e+3` | `1e+4` | matches the hardware ground config; stiffer spline against peer-induced bending |
| `heat_alpha1` / `dyn_heat_tube_radius_m` | `1.5` / `3.0` | `3.0` / `4.0` | stronger, wider heat tubes along peer trajectories, so the global route detours around teammates before the local planner has to |
| `mpc_path_publish_rate_hz` | `0` (every replan) | `10.0` | stop re-anchoring the MPC's reference ~85×/s — the yaw-hunting fix, see below |

**Why `agent_dyn_horizon_sec` exists.** Peers replan about once a second, so
their predicted positions more than a few seconds out are fiction — but the
local avoidance term used to sample the entire ~40 s trajectory against them,
bending the tail for phantom crossings on every replan, which accumulated
into visibly wavy driving *even with all robots far apart*. With the 8 s
horizon (a smooth fade, gradient-checked), the planned paths measure
0.36 rad/m of turning — the same as a single robot's 0.35 — and 8 s is still
4× the avoidance envelope for two closing robots, so nothing that can
physically happen before both robots replan goes unavoided.

**Why `mpc_path_publish_rate_hz` exists.** Every successful replan used to
re-publish the MPC's reference path — ~85×/s, three times per 30 Hz control
cycle — and each copy *starts at the robot's current pose*, so the
controller's target micro-rotated every ~12 ms. Solo, that tracks cleanly.
With three robots sharing one host, timing jitter turns each micro-rotation
into a hard yaw correction (the sim controller carries `w_max: 1.5` for
corner-clearance reasons, see `mpc_sim.yaml`): measured, the yaw command hit
its limit 7.4% of driving time vs 1.8% solo — visible as serpentine driving
and overshooting turns. Capping the reference at 10 Hz restored the solo yaw
profile exactly (limit-hits 2.3%, command RMS 0.382 vs solo's 0.381) and cut
driven-path turning by another 22%; worst-case reference staleness at
0.5 m/s is ~5 cm of travel. The final smoothness stack was validated by a
5-run campaign: 5/5 complete, coverage ≥ 0.99, zero duplicated-targeting
seconds, zero robot-robot contacts, and the fastest mean finish measured
(407 s vs 470–487 s for earlier configs — straighter driving is also faster
driving).

`minpos.enabled` and `select_nearest` are both required — setting either alone
changes nothing.

**How inter-robot avoidance was tuned**, in two eras, because the honest
history is the argument for the current values. Era one: raising the
avoidance *strength* did almost nothing (weight `1e+4` and 4×-finer sampling
measured inside counting noise); what worked was stopping robots from
converging at all — the 4.0 m frontier keep-out plus a 2.0 m clearance
envelope took inter-robot contact from 15 episodes to 0 across 20 runs. Era
two exposed what that 2.0 m envelope cost: it demands 2 m of center
separation in ~2 m corridors, which is geometrically impossible when two
robots pass — so close passes destabilized and shoved trajectories into
walls instead (the base config runs `stat_weight 0`; walls exist only as
corridor constraints, and the giant peer term outvoted them). The current
setting is the measured middle: `planner_Cw 1.4` (passing robots need
1.4 + 0.6 = 2.0 m — exactly what the corridors provide), `dynamic_weight`
back at `1e+3`, and the **peer-yield governor**: in a close pass the robot
that loses the namespace tie-break (the same order that settles claim
contests) crawls at 0.3× while the winner — who never slows, so no deadlock
— drives through. A planner-level variant of the yield and a 4 m yield
radius were both built, measured worse, and reverted (details in the overlay
comments).

**How duplicate targeting was eliminated**, because the numbers are the
argument for the claim/status design: before intent sharing existed, the
robots spent a measured **82–107 seconds per run** with at least two of them
targeting the same spot (time-averaged pair-duplication mean 0.083–0.225 over
a 3-run baseline) — two robots near-equidistant from a frontier both compute
rank 0, and the pursuit lock keeps both driving for the whole leg. With
claims and status sharing on, the 10-run campaign measured
**0 duplicated-targeting seconds in 10/10 runs**. A first design that made
verdicts *permanent* (each 2 Hz rebroadcast re-applied, plus insert
suppression) over-corrected: in 2 of 9 runs a single early verdict near a
doorway walled off an unexplored wing and the team went home at 0.64
coverage. The shipped rule is therefore **a verdict applies once** — new
records near an old verdict are new information, and the local
cooldown/strike lifecycle keeps governing retries.

**What a healthy run looks like.** Robots fan out within the first minute rather
than convoying. A **staggered finish is normal**: each robot announces
`No frontiers left` and returns on its own schedule, so one parking while
another is still exploring is expected, not a stall. An occasional late
resume trip — a robot leaving home again to check a leftover sliver — is
also normal and self-terminating.

Measured over 10 runs with claims + status sharing
(`--scenario ground-exploration-multi`; the pre-claims 20-run campaign in
parentheses where comparable):

| | 3 robots | 1 robot |
|---|---|---|
| Completed, all robots home | **10/10** (20/20) | 18/18 |
| Duplicated-targeting time per run | **0 s in 10/10** (was 82–107 s) | n/a |
| Coverage of the exploration bounds | 0.994 – 0.999 (0.995 – 0.998) | 0.991 – 0.996 |
| Wall-clock to finish | 368 – 528 s, mean **407**¹ | 622 – 887 s |
| Robot-robot contacts | **0 in 10 runs** (0 in 20) | n/a |
| Coordination verified live | 10/10 | n/a |

¹ Duration is from the 5-run campaign on the final configuration (with the
smoothness stack above); earlier configs measured 383–673 s, mean 470–487.
The coordination rows are from the 10-run claims campaign; every subsequent
5-run campaign reproduced them (0 contacts, 0 duplication, all home).

Roughly half the time for slightly better coverage. Coordination is checked at
runtime, not assumed: the harness requires each planner to log MinPos, map
sharing, claim sharing and status sharing as enabled *and* observes peer-pose,
visited-map, claim and frontier-status traffic from all N robots before it
will pass a run. See Known Issues for the wall contact rate.

### Docker (Linux)

```bash
cd mighty/docker

# Autonomous frontier exploration (ACL_office world)
make run-ground-exploration

# Ground robot in a forest world (explores by default; see note above for goal nav)
make run-ground-robot

# Forest navigation with a fixed goal
make run-ground-robot GOAL_X=90 GOAL_Y=0 GOAL_Z=0 ENV=easy_forest
```

> If you hit GPU errors, append `GPU=false` (see the Docker section above).

### Docker (Mac)

Uses Xpra — after launching, open **http://localhost:8080** in your browser.

```bash
cd mighty/docker

make run-mac-ground-exploration      # autonomous exploration
make run-mac-ground-robot            # forest world
```

### Native (Linux)

```bash
cd ~/code/mighty_ws

# Autonomous frontier exploration (ACL_office world)
python3 src/mighty/scripts/run_sim.py --mode exploration-singleagent-ground \
  -s ~/code/mighty_ws/install/setup.bash

# Multi-robot exploration: 3 robots dividing the map via MinPos
python3 src/mighty/scripts/run_sim.py --mode exploration-multiagent-ground \
  -s ~/code/mighty_ws/install/setup.bash

# Ground robot in a forest world (explores by default; see note above for goal nav)
python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot \
  -s ~/code/mighty_ws/install/setup.bash

# Forest navigation with a fixed goal
python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot \
  --env easy_forest --goal 30 0 0 -s ~/code/mighty_ws/install/setup.bash
```

<details>
  <summary><b>Ground Robot Make Targets Reference</b></summary>

  | Target | Description | Options |
  |--------|-------------|---------|
  | `make run-ground-exploration` | Single ground robot, autonomous frontier exploration (ACL_office) | `ENV` (default: ACL_office) |
  | `make run-multiagent-ground-exploration` | N ground robots (default 3) exploring together via MinPos | `NUM_AGENTS` (default: 3), `ENV` (default: ACL_office) |
  | `make run-ground-robot` | Single ground robot crossing a forest world to a fixed goal | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` — all optional; unset means the same defaults as the native command (hard_forest, goal 305 0 0) |
  | `make run-mac-ground-exploration` | Autonomous exploration on Mac (Xpra, browser at localhost:8080) | `ENV` (default: ACL_office) |
  | `make run-mac-multiagent-ground-exploration` | Multi-robot exploration on Mac | `NUM_AGENTS` (default: 3), `ENV` |
  | `make run-mac-ground-robot` | Forest world on Mac (Xpra, browser at localhost:8080) | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` |

</details>

---

## Known Issues

Measured with `scripts/exploration_test.py` (see below). Listed so you know what
you're looking at rather than discovering it mid-demo.

**Rare planner crash under rapid goal churn.** One `double free or
corruption` abort of a planner process was observed in roughly forty
10-robot swap runs — during a `random_slot` run, whose goal turnover is ~3×
the antipodal pattern's. The affected robot stops; the other nine continue.
Suspected heap race between the goal callback and the replan loop; needs an
AddressSanitizer session to pin down. If a robot freezes mid-demo, check its
pane for `process has died`.

**The ground robot occasionally brushes a wall.** Contacts occurred in 8 of 15
exploration runs (two campaigns of 10 and 5). The robot always recovered and
completed the mission — 15/15 reached full coverage and returned home — so this
costs appearance, not success.

With three robots the wall-contact rate rises to 10 of 10 runs (317-1136 events
per run, i.e. per-robot rates comparable to the single-robot case at 3x the
exposure) — the same phenomenon, not a multi-robot regression.

It is not uniformly distributed. Taking each run's closest approach to an
obstacle as the proxy for where it happened, they cluster at two spots in
`ACL_office`:

| Closest-approach point | Runs |
|------------------------|------|
| ~(12.6, −11.4)         | 4    |
| ~(2.6, −11.6)          | 4    |
| ~(−4.2, −23.1)         | 1    |
| ~(1.7, −21.3)          | 1    |

Note that which corner dominates shifts between campaigns — the first campaign
after the cornering fix hit (12.6, −11.4) repeatedly and never touched (2.6,
−11.6), while the next campaign did the opposite. Treat the specific
coordinates as "the tight corners of this map", not as a fixed list.

The cause is a speed/accuracy trade rather than a clearance one. A
differential-drive robot's tightest arc while moving is `v_max / w_max` (see
`mpc_sim.yaml`); at those two corners the achievable arc plus path-tracking error
still doesn't quite fit. Widening the planner's clearance does **not** help — we
measured that directly: raising the safe-corridor inflation from 0.30 m to 0.45 m
left contacts unchanged while dropping coverage from 0.991 to 0.948 and causing
the robot to stop returning home. Lowering `v_max` further would trade run time
for accuracy everywhere to fix two corners, which didn't seem worth it for a demo.

**Multi-robot wall contacts scale with the robot count.** Three robots produce
wall contacts in essentially every run (the single-robot entry above explains
why); per robot the rate is comparable to the solo case, at 3x the exposure.

**Occasional low-speed robot-robot contact in pinch encounters.** Two robots
meeting head-on in a corridor or doorway is the one encounter class no
clearance radius resolves — the geometry requires one robot to not be there.
The earlier "0 contacts in 40 runs" era achieved that by an oversized 2.0 m
avoidance envelope that paid at the walls instead (unstable, wall-shoved
trajectories in close passes). With the current corridor-feasible tuning
plus the yield governor, pooled measurement over 10 runs: **8 contact-free**
(minimum pair separations 0.61–1.00 m), 2 runs with one low-speed pinch
contact each (one ~0.1 s brush, one few-second shove; the robots slide past,
recover, and both runs completed with ≥0.99 coverage — no instability, no
wall violations, no failed runs). The distribution is heavy-tailed: roughly
one pinch meeting every 2–5 runs, and pinch geometry decides between brush
and shove. Eliminating the tail outright needs a *hard* mechanism — e.g.
priority-based reservation of narrow segments over the existing claims
channel — which is deliberate future work, not a tuning knob.

**Residual multi-robot driving waviness on a shared host — a simulation
artifact, not a planner or controller defect.** Two real contributors were
found and fixed (`agent_dyn_horizon_sec` and `mpc_path_publish_rate_hz`,
documented in Multi-Robot Exploration above); together they took driven-path
turning at three robots from 1.07–1.29 rad/m down to **0.83**, with the yaw
command profile now identical to solo. The remaining gap to the one-robot
floor (0.54 rad/m — measured even when that one robot runs the full
multi-robot config) tracks host contention: with three robots Gazebo's
real-time factor drops to 0.91 (1.00 solo), the three MPCs each run near a
full core, and the stack runs on wall clock — so controllers compute
commands for real-time dynamics that the slowed, jittery physics doesn't
quite deliver. On hardware each robot has its own computer, so this doesn't
transfer. On a shared sim host: use a beefier machine, or `NUM_AGENTS=2` for
smoother-looking runs. Migrating the multi sim to `use_sim_time` would remove
the artifact in principle but re-times the entire stack and its validation —
a deliberate future project, not a config flip.

**Multi-agent UAVs pass closer than their own bounding box.** In the antipodal
swap scenario the closest approach between two agents was 0.451 m with 5 agents
and 0.329 m with 10, against `drone_bbox: [0.5, 0.5, 0.5]` — two 0.5 m boxes
overlap below 0.5 m centre-to-centre. `fake_sim` has no collision physics so
nothing visibly goes wrong, but the margin is real and would matter on hardware,
where that 0.5 m figure is documented as already including prop guards.

The obvious knob is not the cause: raising `dyn_base_inflation_m` from 0.2 to 0.4
moved the mean closest approach from 0.394 m to 0.393 m over two runs each at 5
and 10 agents (10 agents improved, 5 agents got worse — i.e. noise), so it was
reverted. Whatever bounds inter-agent spacing lives in the local trajectory
optimisation rather than in the dynamic-obstacle inflation, and chasing it was
out of scope here. Everything else about the scenario is healthy: all agents
launch, each completes 20–21 crossings in 8 minutes, altitude holds, nothing
crashes.

**`spawn_entity` can fail at startup.** Occasionally the robot is never spawned
into Gazebo ("Spawn service failed. Exiting."), so there is no sensor data — yet
the run still looks successful, because `fake_sim` integrates the commanded
trajectory whether or not anything is perceiving. If a Gazebo run reaches its
goal implausibly smoothly, check the launch log for that message.

**`Default goal z is lower than the ground level`** is logged at ERROR severity
on healthy `--mode interactive` runs. It is one of the few messages loud enough
to survive the default `--log-level error`, so it is usually the first thing a
new user sees; it does not appear to affect the flight.

## Automated Simulation Tests

`scripts/exploration_test.py` runs a scenario headless, watches it over the ROS
graph, and judges it against explicit gates. It exists because these are
stochastic, minutes-long, closed-loop simulations: a single successful run
proves very little, so the interesting question is how often a scenario succeeds.

```bash
source ~/code/mighty_ws/install/setup.bash    # required: monitors dynus_interfaces topics

# What can be tested
python3 src/mighty/scripts/exploration_test.py scenarios

# One run
python3 src/mighty/scripts/exploration_test.py run --scenario ground-exploration

# Ten runs, with a summary table
python3 src/mighty/scripts/exploration_test.py campaign \
    --scenario ground-exploration-multi --runs 10 --out /tmp/campaign
```

Each run writes `run<NN>.json` (metrics, per-gate verdicts, and position/goal
traces for offline reconstruction) plus complete per-pane logs; a campaign also
writes `summary.md` and `summary.json`.

**A run can never hang.** Three independent stoppers bound every run: a per-run
wall-clock cap enforced by the harness (not by the simulation behaving itself),
early abort the moment a stall is detected, and a campaign-wide deadline after
which remaining runs are reported as *not run* rather than assumed passing.

Stalls are classified rather than lumped together, because the fix differs:

| Class | Meaning |
|-------|---------|
| `STUCK_NO_PROGRESS` | Not translating at all while a goal is active — wedged or motors stalled |
| `STUCK_NO_COVERAGE` | Learned nothing new for four minutes and hasn't met the coverage gate |
| `STUCK_OSCILLATING` | Same, but the evidence shows repeated retargeting and lots of driving with no net progress |
| `TIMEOUT` | Hit the wall-clock cap with no terminal state |

Stall detection is deliberately tuned to be *specific* rather than sensitive: a
false alarm sends you chasing a phantom bug, while a slow catch costs at most one
wall-clock cap. In particular, greedy nearest-frontier exploration legitimately
makes long round trips, so motion alone is never treated as evidence of a stall.

---

## Acknowledgments

The L-BFGS solver implementation in this repository (`src/mighty/lbfgs_solver.cpp`, `include/mighty/lbfgs_solver.hpp`) is adapted from [ZJU-FAST-Lab/GCOPTER](https://github.com/ZJU-FAST-Lab/GCOPTER/blob/main/gcopter/include/gcopter/lbfgs.hpp). We thank the authors for making their implementation publicly available.

The ground robot's MPC path-tracking controller ([kotakondo/mpc](https://github.com/kotakondo/mpc), pulled in via `mighty.repos`) was originally developed by **[Lucas Jia (@lucas-yyy000)](https://github.com/lucas-yyy000)**. The ground-robot simulation depends on it.

The `ACL_office` Gazebo environment used in the ground-robot exploration simulations was developed by **[Paul Leonhard Kohler (@p-kohler)](https://github.com/p-kohler)**.

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
