#!/bin/bash
set -e

# MIGHTY Complete Setup Script
# Run this once to install all dependencies, clone repos, and build everything
# Combines system setup, dependency installation, and workspace building

echo "============================================="
echo "MIGHTY Complete Setup Script"
echo "============================================="
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
# STEP 1: System Updates and Basic Software
# ============================================
echo "============================================="
echo "STEP 1: System Updates and Basic Software"
echo "============================================="
sudo rm -rf /var/lib/apt/lists/*
sudo apt update
sudo apt upgrade -y
sudo apt-get update
sudo apt-get upgrade -y
sudo apt-get install -q -y --no-install-recommends \
    git tmux vim wget tmuxp make openssh-server net-tools g++ xterm python3-pip

pip install pymavlink
sudo apt install -y libomp-dev libpcl-dev libeigen3-dev

# ============================================
# STEP 2: ROS 2 Humble Installation
# ============================================
echo ""
echo "============================================="
echo "STEP 2: Installing ROS 2 Humble"
echo "============================================="

# Check if ROS 2 is already installed
if [ ! -d "/opt/ros/humble" ]; then
    echo "Installing ROS 2 Humble..."
    sudo apt update && sudo apt install locales
    sudo locale-gen en_US en_US.UTF-8
    sudo update-locale LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
    export LANG=en_US.UTF-8
    sudo apt install software-properties-common
    echo -e "\n" | sudo add-apt-repository universe
    sudo apt update && sudo apt install curl -y
    sudo curl -sSL https://raw.githubusercontent.com/ros/rosdistro/master/ros.key -o /usr/share/keyrings/ros-archive-keyring.gpg
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/ros-archive-keyring.gpg] http://packages.ros.org/ros2/ubuntu $(. /etc/os-release && echo $UBUNTU_CODENAME) main" | sudo tee /etc/apt/sources.list.d/ros2.list > /dev/null
    sudo apt install -y ros-humble-desktop
    sudo apt install -y ros-dev-tools
else
    echo "ROS 2 Humble already installed, skipping..."
fi

export ROS_DISTRO=humble
source /opt/ros/humble/setup.bash

# ============================================
# STEP 3: ROS Dependencies
# ============================================
echo ""
echo "============================================="
echo "STEP 3: Installing ROS Dependencies"
echo "============================================="
sudo apt-get install -y \
    ros-${ROS_DISTRO}-gazebo-* \
    ros-${ROS_DISTRO}-pcl-conversions \
    ros-${ROS_DISTRO}-example-interfaces \
    ros-${ROS_DISTRO}-pcl-ros \
    ros-${ROS_DISTRO}-rviz2 \
    ros-${ROS_DISTRO}-rqt-gui \
    ros-${ROS_DISTRO}-rqt-gui-py \
    ros-${ROS_DISTRO}-tf2-tools \
    ros-${ROS_DISTRO}-tf-transformations \
    ros-${ROS_DISTRO}-turtlesim \
    ros-${ROS_DISTRO}-rqt* \
    ros-${ROS_DISTRO}-rviz-common \
    nlohmann-json3-dev \
    libpcl-dev \
    build-essential

# ============================================
# STEP 4: Create Workspaces and Import Repositories
# ============================================
echo ""
echo "============================================="
echo "STEP 4: Creating Workspaces and Importing Repos"
echo "============================================="

# Create all workspaces
mkdir -p "$MIGHTY_WS/src"
mkdir -p "$DECOMP_WS/src"
mkdir -p "$LIVOX_WS/src"

# Check if MIGHTY is already cloned
if [ ! -d "$MIGHTY_WS/src/mighty" ]; then
    echo "Cloning MIGHTY..."
    cd "$MIGHTY_WS/src"
    git clone https://github.com/mit-acl/mighty.git
    cd mighty
    git switch multiagent_sim
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
    echo "Moving DecompROS2 to decomp_ws..."
    mv "$MIGHTY_WS/src/DecompROS2" "$DECOMP_WS/src/"
fi

if [ -d "$MIGHTY_WS/src/livox_ros_driver2" ]; then
    echo "Moving livox_ros_driver2 to livox_ws..."
    mv "$MIGHTY_WS/src/livox_ros_driver2" "$LIVOX_WS/src/"
fi

if [ -d "$MIGHTY_WS/src/Livox-SDK2" ]; then
    echo "Moving Livox-SDK2 to parent directory..."
    mv "$MIGHTY_WS/src/Livox-SDK2" "$LIVOX_SDK_DIR"
fi

# ============================================
# STEP 5: Build DecompROS2
# ============================================
echo ""
echo "============================================="
echo "STEP 5: Building DecompROS2"
echo "============================================="
cd "$DECOMP_WS"
source /opt/ros/humble/setup.bash

echo "Building decomp_util first (required)..."
colcon build --packages-select decomp_util --cmake-args -DCMAKE_BUILD_TYPE=Release

echo "Building remaining decomp packages..."
source "$DECOMP_WS/install/setup.bash"
colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

# ============================================
# STEP 6: Build Livox-SDK2
# ============================================
echo ""
echo "============================================="
echo "STEP 6: Building Livox-SDK2"
echo "============================================="
cd "$LIVOX_SDK_DIR"
mkdir -p build
cd build
cmake .. && make -j$(nproc)
echo "Installing Livox-SDK2..."
sudo make install

# ============================================
# STEP 7: Build livox_ros_driver2
# ============================================
echo ""
echo "============================================="
echo "STEP 7: Building livox_ros_driver2"
echo "============================================="
cd "$LIVOX_WS/src/livox_ros_driver2"
source /opt/ros/humble/setup.bash
./build.sh humble

# ============================================
# STEP 8: Build MIGHTY Workspace
# ============================================
echo ""
echo "============================================="
echo "STEP 8: Building MIGHTY Workspace"
echo "============================================="
cd "$MIGHTY_WS"
source /opt/ros/humble/setup.bash
source "$DECOMP_WS/install/setup.bash"
export CMAKE_PREFIX_PATH="$LIVOX_WS/install/livox_ros_driver2:$DECOMP_WS/install/decomp_util:$CMAKE_PREFIX_PATH"
export LD_LIBRARY_PATH="$LIVOX_WS/install/livox_ros_driver2/lib:$LD_LIBRARY_PATH"

colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

# ============================================
# STEP 9: Setup Bash Configuration
# ============================================
echo ""
echo "============================================="
echo "STEP 9: Setting Up Bash Configuration"
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
    echo "export ROS_DOMAIN_ID=10" >> ~/.bashrc
    echo "" >> ~/.bashrc
    echo "# Livox library path" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=$LIVOX_WS/install/livox_ros_driver2/lib:\$LD_LIBRARY_PATH:/usr/local/lib:/usr/local/include" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=/opt/ros/humble/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=$DECOMP_WS/install/decomp_ros_msgs/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
    echo "export LD_LIBRARY_PATH=/opt/ros/humble/lib/x86_64-linux-gnu:\$LD_LIBRARY_PATH" >> ~/.bashrc
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
# STEP 10: Summary
# ============================================
echo ""
echo "============================================="
echo "✅ MIGHTY Setup Complete!"
echo "============================================="
echo ""
echo "All dependencies installed and built:"
echo "  ✅ ROS 2 Humble"
echo "  ✅ System dependencies"
echo "  ✅ DecompROS2 (decomp_ws)"
echo "  ✅ Livox-SDK2"
echo "  ✅ livox_ros_driver2 (livox_ws)"
echo "  ✅ MIGHTY and all ROS packages (mighty_ws)"
echo ""
echo "Workspaces:"
echo "  - MIGHTY:  $MIGHTY_WS"
echo "  - Decomp:  $DECOMP_WS"
echo "  - Livox:   $LIVOX_WS"
echo ""
echo "To use MIGHTY immediately:"
echo "  source ~/.bashrc"
echo ""
echo "Or manually source:"
echo "  source $DECOMP_WS/install/setup.bash"
echo "  source $MIGHTY_WS/install/setup.bash"
echo ""
echo "All future terminal sessions will have MIGHTY ready to use!"
echo "============================================="
