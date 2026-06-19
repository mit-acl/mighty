#!/bin/bash
set -e

# MIGHTY Complete Setup Script
# Run this once to install all dependencies, clone repos, and build everything
# Combines system setup, dependency installation, and workspace building
#
# Usage: ./setup.sh [-j N] [--jetson]
#   -j N      Number of parallel jobs for building (default: all CPUs)
#   --jetson  Skip packages unavailable on Jetson (ARM64), e.g. Gazebo

# Parse arguments
NUM_JOBS=$(nproc)
JETSON=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            NUM_JOBS="$2"
            shift 2
            ;;
        --jetson)
            JETSON=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: ./setup.sh [-j N] [--jetson]"
            exit 1
            ;;
    esac
done

# Auto-detect Jetson if not explicitly set
if ! $JETSON && [ -f /etc/nv_tegra_release ]; then
    echo "Detected NVIDIA Jetson platform, enabling --jetson mode automatically."
    JETSON=true
fi

echo "============================================="
echo "MIGHTY Complete Setup Script"
echo "============================================="
echo "Using $NUM_JOBS parallel jobs for building"
if $JETSON; then
    echo "Jetson mode: skipping Gazebo and other x86-only packages"
fi
echo ""

# Prompt for sudo password once at the beginning and keep it cached
echo "This script requires sudo access for:"
echo "  - Installing system packages"
echo "  - Installing Livox-SDK2"
echo ""
sudo -v
# Keep sudo alive in background (updates timestamp every 60 seconds)
while true; do sudo -n true; sleep 60; kill -0 "$$" || exit; done 2>/dev/null &

# Determine workspace locations
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CODE_DIR="${CODE_DIR:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
MIGHTY_WS="$CODE_DIR/mighty_ws"
DECOMP_WS="$CODE_DIR/decomp_ws"
LIVOX_WS="$CODE_DIR/livox_ws"
LIVOX_SDK_DIR="$CODE_DIR/Livox-SDK2"

echo "Installation directories:"
echo "  CODE_DIR: $CODE_DIR"
echo "  MIGHTY_WS: $MIGHTY_WS"
echo "  DECOMP_WS: $DECOMP_WS"
echo "  LIVOX_WS: $LIVOX_WS"
echo "  LIVOX_SDK_DIR: $LIVOX_SDK_DIR"
echo ""

# ============================================
# STEP 1: Fix Repository Issues
# ============================================
echo "============================================="
echo "STEP 1: Fixing Repository Issues"
echo "============================================="

# Only add ROS 2 repository if none exists — never overwrite existing setup
if ls /etc/apt/sources.list.d/ros2* > /dev/null 2>&1; then
    echo "ROS 2 repository already configured, keeping existing setup."
else
    echo "No ROS 2 repository found, adding..."
    sudo curl -sSL https://raw.githubusercontent.com/ros/rosdistro/master/ros.key -o /usr/share/keyrings/ros-archive-keyring.gpg
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/ros-archive-keyring.gpg] http://packages.ros.org/ros2/ubuntu $(. /etc/os-release && echo $UBUNTU_CODENAME) main" | sudo tee /etc/apt/sources.list.d/ros2.list > /dev/null
fi

# ============================================
# STEP 2: System Updates and Basic Software
# ============================================
echo ""
echo "============================================="
echo "STEP 2: System Updates and Basic Software"
echo "============================================="
sudo apt update
sudo apt upgrade -y
sudo apt-get install -q -y --no-install-recommends \
    git tmux vim wget tmuxp make openssh-server net-tools g++ xterm python3-pip

pip install pymavlink
sudo apt install -y libomp-dev libpcl-dev libeigen3-dev

# ============================================
# STEP 3: ROS 2 Humble Installation
# ============================================
echo ""
echo "============================================="
echo "STEP 3: Installing ROS 2 Humble"
echo "============================================="

# Check if ROS 2 is already installed
if [ ! -d "/opt/ros/humble" ]; then
    echo "Installing ROS 2 Humble..."
    sudo apt update && sudo apt install -y locales
    sudo locale-gen en_US en_US.UTF-8
    sudo update-locale LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
    export LANG=en_US.UTF-8
    sudo apt install -y software-properties-common
    echo -e "\n" | sudo add-apt-repository universe
    sudo apt update && sudo apt install -y curl

    # ROS 2 repo already set up in Step 1 — just update
    sudo apt update
    sudo apt install -y ros-humble-desktop
    sudo apt install -y ros-dev-tools
else
    echo "ROS 2 Humble already installed, skipping..."
fi

export ROS_DISTRO=humble
source /opt/ros/humble/setup.bash

# ============================================
# STEP 4: ROS Dependencies
# ============================================
echo ""
echo "============================================="
echo "STEP 4: Installing ROS Dependencies"
echo "============================================="
ROS_DEPS=(
    ros-${ROS_DISTRO}-pcl-conversions
    ros-${ROS_DISTRO}-example-interfaces
    ros-${ROS_DISTRO}-pcl-ros
    ros-${ROS_DISTRO}-rviz2
    ros-${ROS_DISTRO}-rqt-gui
    ros-${ROS_DISTRO}-rqt-gui-py
    ros-${ROS_DISTRO}-tf2-tools
    ros-${ROS_DISTRO}-tf-transformations
    ros-${ROS_DISTRO}-turtlesim
    ros-${ROS_DISTRO}-rqt*
    ros-${ROS_DISTRO}-rviz-common
    nlohmann-json3-dev
    libpcl-dev
    build-essential
)

if ! $JETSON; then
    ROS_DEPS+=( ros-${ROS_DISTRO}-gazebo-* )
else
    echo "Jetson mode: skipping Gazebo packages (not available on ARM64)"
fi

sudo apt-get install -y "${ROS_DEPS[@]}"

# ============================================
# STEP 5: Create Workspaces and Import Repositories
# ============================================
echo ""
echo "============================================="
echo "STEP 5: Creating Workspaces and Importing Repos"
echo "============================================="

# Create all workspaces
mkdir -p "$MIGHTY_WS/src"
mkdir -p "$DECOMP_WS/src"
mkdir -p "$LIVOX_WS/src"

# Check if MIGHTY is already cloned
if [ ! -d "$MIGHTY_WS/src/mighty" ]; then
    echo "Cloning MIGHTY..."
    cd "$MIGHTY_WS/src"
    git clone --depth 1 --branch v0.0.6 https://github.com/mit-acl/mighty.git
else
    echo "MIGHTY already exists, updating..."
    cd "$MIGHTY_WS/src/mighty"
    git fetch
fi

# Import all dependencies using mighty.repos
echo "Importing dependencies from mighty.repos..."
cd "$MIGHTY_WS"
vcs import src < src/mighty/mighty.repos || true

# Move special dependencies to their workspaces
if [ -d "$MIGHTY_WS/src/DecompROS2" ]; then
    if [ -d "$DECOMP_WS/src/DecompROS2" ]; then
        echo "DecompROS2 already exists in decomp_ws, skipping move..."
        rm -rf "$MIGHTY_WS/src/DecompROS2"
    else
        echo "Moving DecompROS2 to decomp_ws..."
        mv "$MIGHTY_WS/src/DecompROS2" "$DECOMP_WS/src/"
    fi
fi

if [ -d "$MIGHTY_WS/src/livox_ros_driver2" ]; then
    if [ -d "$LIVOX_WS/src/livox_ros_driver2" ]; then
        echo "livox_ros_driver2 already exists in livox_ws, skipping move..."
        rm -rf "$MIGHTY_WS/src/livox_ros_driver2"
    else
        echo "Moving livox_ros_driver2 to livox_ws..."
        mv "$MIGHTY_WS/src/livox_ros_driver2" "$LIVOX_WS/src/"
    fi
fi

if [ -d "$MIGHTY_WS/src/Livox-SDK2" ]; then
    if [ -d "$LIVOX_SDK_DIR" ]; then
        echo "Livox-SDK2 already exists in parent directory, skipping move..."
        rm -rf "$MIGHTY_WS/src/Livox-SDK2"
    else
        echo "Moving Livox-SDK2 to parent directory..."
        mv "$MIGHTY_WS/src/Livox-SDK2" "$LIVOX_SDK_DIR"
    fi
fi

# ============================================
# STEP 5: Build DecompROS2
# ============================================
echo ""
echo "============================================="
echo "STEP 6: Building DecompROS2"
echo "============================================="
cd "$DECOMP_WS"
source /opt/ros/humble/setup.bash

echo "Building decomp_util first (required)..."
colcon build --packages-select decomp_util --parallel-workers "$NUM_JOBS" --cmake-args -DCMAKE_BUILD_TYPE=Release

echo "Building remaining decomp packages..."
source "$DECOMP_WS/install/setup.bash"
colcon build --parallel-workers "$NUM_JOBS" --cmake-args -DCMAKE_BUILD_TYPE=Release

# ============================================
# STEP 6: Build Livox-SDK2
# ============================================
echo ""
echo "============================================="
echo "STEP 7: Building Livox-SDK2"
echo "============================================="
cd "$LIVOX_SDK_DIR"
mkdir -p build
cd build
cmake .. && make -j"$NUM_JOBS"
echo "Installing Livox-SDK2..."
sudo make install

# ============================================
# STEP 7: Build livox_ros_driver2
# ============================================
echo ""
echo "============================================="
echo "STEP 8: Building livox_ros_driver2"
echo "============================================="
cd "$LIVOX_WS/src/livox_ros_driver2"
source /opt/ros/humble/setup.bash
./build.sh humble

# ============================================
# STEP 8: Build MIGHTY Workspace
# ============================================
echo ""
echo "============================================="
echo "STEP 9: Building MIGHTY Workspace"
echo "============================================="
cd "$MIGHTY_WS"
source /opt/ros/humble/setup.bash
source "$DECOMP_WS/install/setup.bash"
export CMAKE_PREFIX_PATH="$LIVOX_WS/install/livox_ros_driver2:$DECOMP_WS/install/decomp_util:$CMAKE_PREFIX_PATH"
export LD_LIBRARY_PATH="$LIVOX_WS/install/livox_ros_driver2/lib:$LD_LIBRARY_PATH"

CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release)
COLCON_BUILD_ARGS=(--parallel-workers "$NUM_JOBS")

if $JETSON; then
    COLCON_BUILD_ARGS+=(
        --packages-ignore
            gazebo_dev
            gazebo_ros
            gazebo_plugins
            gazebo_ros_pkgs
            realsense_gazebo_plugin
            ros2_livox_simulation
        --executor sequential
    )
    CMAKE_ARGS+=(-DBUILD_SIMULATION=OFF)
    # Limit compiler parallelism to avoid OOM on memory-constrained Jetson boards
    export CMAKE_BUILD_PARALLEL_LEVEL=1
    export MAKEFLAGS="-j1"
    echo "Jetson mode: skipping Gazebo simulation packages, building sequentially to save RAM"
fi

colcon build "${COLCON_BUILD_ARGS[@]}" --cmake-args "${CMAKE_ARGS[@]}"

# ============================================
# STEP 9: Setup Bash Configuration
# ============================================
echo ""
echo "============================================="
echo "STEP 10: Setting Up Bash Configuration"
echo "============================================="

# Add to bashrc if not already there
if ! grep -q "MIGHTY Setup" ~/.bashrc; then
    echo "" >> ~/.bashrc
    echo "# ============================================" >> ~/.bashrc
    echo "# MIGHTY Setup" >> ~/.bashrc
    echo "# ============================================" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "# Source ROS 2 Humble" >> ~/.bashrc
    echo "source /opt/ros/humble/setup.bash" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "# ROS 2 RTPS network" >> ~/.bashrc
    echo "export ROS_DOMAIN_ID=20" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "# Livox library path" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=$LIVOX_WS/install/livox_ros_driver2/lib:\$LD_LIBRARY_PATH:/usr/local/lib:/usr/local/include" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=/opt/ros/humble/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=$DECOMP_WS/install/decomp_ros_msgs/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=/opt/ros/humble/lib/$(dpkg --print-architecture | sed 's/amd64/x86_64-linux-gnu/' | sed 's/arm64/aarch64-linux-gnu/'):\$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "# Source MIGHTY workspaces" >> ~/.bashrc
    echo "source $DECOMP_WS/install/setup.bash" >> ~/.bashrc
    echo "source $MIGHTY_WS/install/setup.bash" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "Setup added to ~/.bashrc"
else
    echo "MIGHTY setup already exists in ~/.bashrc, skipping..."
fi

# ============================================
# STEP 11: Summary
# ============================================
echo ""
echo "============================================="
echo "MIGHTY Setup Complete!"
echo "============================================="
echo ""
echo "Workspaces:"
echo "  - MIGHTY:  $MIGHTY_WS"
echo "  - Decomp:  $DECOMP_WS"
echo "  - Livox:   $LIVOX_WS"
echo ""
echo "To get started:"
echo "  source ~/.bashrc"
echo "  cd $MIGHTY_WS"
echo ""
echo "  # Single-agent interactive simulation (click goals in RViz2 with \"2D Goal Pose\")"
echo "  python3 src/mighty/scripts/run_sim.py --mode interactive --setup-bash ~/code/mighty_ws/install/setup.bash"
echo ""
echo "See README.md for more simulation modes and options."
echo "============================================="
