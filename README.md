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

## Running the Simulations

Every scenario, and the one command that starts it. Native commands assume the
workspace is built and run from `~/code/mighty_ws`; Docker commands run from
`mighty/docker`. See [Installation](#installation) first.

| Scenario | Native | Docker |
|----------|--------|--------|
| **Ground robot, autonomous exploration** (Pioneer 3-AT explores `ACL_office` and returns to its start) | `python3 src/mighty/scripts/run_sim.py --mode exploration-singleagent-ground -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-exploration` |
| Ground robot in a forest world | `python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot -s ~/code/mighty_ws/install/setup.bash` | `make run-ground-robot` |
| Single UAV, goals clicked in RViz | `python3 src/mighty/scripts/run_sim.py --mode interactive -s ~/code/mighty_ws/install/setup.bash` | `make run-interactive` |
| Multi-agent UAVs swapping on a circle | `python3 src/mighty/scripts/run_sim.py --mode multiagent -s ~/code/mighty_ws/install/setup.bash` | `make run-multiagent` |
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
  --goal, -g          Goal position X Y Z for gazebo mode (default: 105 0 3)
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

  > The Docker entrypoint defaults `GOAL_X` to 305, so `make run-gazebo` and
  > `run_sim.py --mode gazebo` do **not** start with the same goal.
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

## Ground Robot Simulation (Single Agent)

Besides the aerial vehicle, MIGHTY drives a single wheeled ground robot (Pioneer
3-AT) in Gazebo. Two single-agent ground-robot scenarios are available:

| Scenario | What it does | Goal |
|----------|--------------|------|
| **Autonomous exploration** | Robot explores an unknown `ACL_office` map on its own using frontier detection. | None — self-driven |
| **Forest navigation** | Robot drives through a forest world (e.g. `easy_forest`). | Frontier exploration by default; a fixed goal after a one-line config change (see note) |

> **Note on goal navigation.** `config/mighty_ground_robot.yaml` ships with
> `exploration.enabled: true`, so the forest scenario **autonomously explores**
> out of the box and the `--goal` / `GOAL_*` values are ignored. To make the
> robot drive to a fixed goal instead, set `exploration.enabled: false` in
> `config/mighty_ground_robot.yaml`, then provide a goal. In Docker, use `DEV=1`
> to edit the config live without rebuilding the image.

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

### Docker (Linux)

```bash
cd mighty/docker

# Autonomous frontier exploration (ACL_office world)
make run-ground-exploration

# Ground robot in a forest world (explores by default; see note above for goal nav)
make run-ground-robot

# Forest navigation with a fixed goal (requires exploration.enabled:false)
make run-ground-robot GOAL_X=30 GOAL_Y=0 GOAL_Z=0 ENV=easy_forest DEV=1
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

# Ground robot in a forest world (explores by default; see note above for goal nav)
python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot \
  -s ~/code/mighty_ws/install/setup.bash

# Forest navigation with a fixed goal (requires exploration.enabled:false)
python3 src/mighty/scripts/run_sim.py --mode gazebo --ground-robot \
  --env easy_forest --goal 30 0 0 -s ~/code/mighty_ws/install/setup.bash
```

<details>
  <summary><b>Ground Robot Make Targets Reference</b></summary>

  | Target | Description | Options |
  |--------|-------------|---------|
  | `make run-ground-exploration` | Single ground robot, autonomous frontier exploration (ACL_office) | `ENV` (default: ACL_office) |
  | `make run-ground-robot` | Single ground robot in a forest world (explores unless `exploration.enabled:false`) | `GOAL_X`, `GOAL_Y`, `GOAL_Z` (default: 30, 0, 0), `ENV` (default: easy_forest) |
  | `make run-mac-ground-exploration` | Autonomous exploration on Mac (Xpra, browser at localhost:8080) | `ENV` (default: ACL_office) |
  | `make run-mac-ground-robot` | Forest world on Mac (Xpra, browser at localhost:8080) | `GOAL_X`, `GOAL_Y`, `GOAL_Z`, `ENV` |

</details>

---

## Known Issues

Measured with `scripts/exploration_test.py` (see below). Listed so you know what
you're looking at rather than discovering it mid-demo.

**The ground robot occasionally brushes a wall.** Contacts occurred in 8 of 15
exploration runs (two campaigns of 10 and 5). The robot always recovered and
completed the mission — 15/15 reached full coverage and returned home — so this
costs appearance, not success.

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
    --scenario ground-exploration --runs 10 --out /tmp/campaign
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

The ground robot's MPC path-tracking controller ([kotakondo/mpc](https://github.com/kotakondo/mpc), pulled in via `mighty.repos`) was originally developed by **[Lucas Jia (@lucas-yyy000)](https://github.com/lucas-yyy000)**, with later contributions from Sera Ham and Elon Raya. The ground-robot simulation depends on it.

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
