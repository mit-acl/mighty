FROM ros:humble

LABEL Maintainer="Kota Kondo <kkondo@mit.edu>" \
      Description="Ubuntu 22.04, ROS2 humble for MIGHTY"

ENV ROS_DISTRO=humble
ENV DEBIAN_FRONTEND=noninteractive

# System packages
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -q -y --no-install-recommends \
        git tmux vim wget tmuxp make gdb && \
    rm -rf /var/lib/apt/lists/*

# ROS 2 dependencies
RUN apt-get update && apt-get install -y \
        ros-${ROS_DISTRO}-gazebo-* \
        ros-${ROS_DISTRO}-pcl-conversions \
        ros-${ROS_DISTRO}-example-interfaces \
        ros-${ROS_DISTRO}-pcl-ros \
        ros-${ROS_DISTRO}-rviz2 \
        ros-${ROS_DISTRO}-rviz-common \
        ros-${ROS_DISTRO}-rqt-gui \
        ros-${ROS_DISTRO}-rqt-gui-py \
        ros-${ROS_DISTRO}-tf2-tools \
        ros-${ROS_DISTRO}-tf-transformations \
        ros-${ROS_DISTRO}-desktop \
        ros-${ROS_DISTRO}-turtlesim \
        ros-dev-tools \
        python3-vcstool \
        libpcl-dev build-essential nlohmann-json3-dev \
        libeigen3-dev \
        libtbb-dev && \
    rm -rf /var/lib/apt/lists/*

# Workspace layout (paths must match mighty.sh / Makefile DEV_MOUNTS)
ENV CODE_DIR=/home/kkondo/code
ENV MIGHTY_WS=${CODE_DIR}/mighty_ws
ENV DECOMP_WS=${CODE_DIR}/decomp_ws
ENV LIVOX_WS=${CODE_DIR}/livox_ws
ENV LIVOX_SDK_DIR=${CODE_DIR}/Livox-SDK2

# Clone MIGHTY at a pinned release tag (override with --build-arg MIGHTY_VERSION=<tag>).
# All other dep versions live in mighty.repos — single source of truth.
RUN mkdir -p ${MIGHTY_WS}/src
WORKDIR ${MIGHTY_WS}/src
ARG DUMMY=unknown
ARG MIGHTY_VERSION=v0.0.5
RUN DUMMY=${DUMMY} git clone --depth 1 --branch ${MIGHTY_VERSION} https://github.com/mit-acl/mighty.git

# Import every dependency (pinned in mighty.repos)
WORKDIR ${MIGHTY_WS}
RUN vcs import src < src/mighty/mighty.repos

# Move special-build deps out of mighty_ws into their own trees
RUN mkdir -p ${DECOMP_WS}/src ${LIVOX_WS}/src && \
    mv ${MIGHTY_WS}/src/DecompROS2       ${DECOMP_WS}/src/ && \
    mv ${MIGHTY_WS}/src/livox_ros_driver2 ${LIVOX_WS}/src/ && \
    mv ${MIGHTY_WS}/src/Livox-SDK2        ${LIVOX_SDK_DIR}

SHELL ["/bin/bash", "-c"]

# Build parallelism — for amd64 emulation on Mac, limit with: --build-arg JOBS=4
ARG JOBS=""

# Build DecompROS2 (decomp_util first, then the rest)
WORKDIR ${DECOMP_WS}
RUN source /opt/ros/humble/setup.bash && \
    colcon build --packages-select decomp_util \
                 ${JOBS:+--parallel-workers $JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release
RUN source /opt/ros/humble/setup.bash && \
    source ${DECOMP_WS}/install/setup.bash && \
    colcon build ${JOBS:+--parallel-workers $JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release

# Build Livox-SDK2 (cmake install to /usr/local)
WORKDIR ${LIVOX_SDK_DIR}
RUN mkdir -p build && cd build && \
    cmake .. && make -j${JOBS:-$(nproc)} && make install

# Build livox_ros_driver2 (uses its own build.sh for package file setup)
WORKDIR ${LIVOX_WS}/src/livox_ros_driver2
RUN source /opt/ros/humble/setup.bash && ./build.sh humble

# Build the MIGHTY workspace
WORKDIR ${MIGHTY_WS}
RUN source /opt/ros/humble/setup.bash && \
    source ${DECOMP_WS}/install/setup.bash && \
    export CMAKE_PREFIX_PATH=${LIVOX_WS}/install/livox_ros_driver2:${DECOMP_WS}/install/decomp_util:${CMAKE_PREFIX_PATH} && \
    export MAKEFLAGS="-j${JOBS:-$(nproc)}" && \
    colcon build ${JOBS:+--parallel-workers $JOBS} \
                 --cmake-args -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_EIGEN3=ON

# Gazebo audio fix
RUN echo 'pcm.!default { type plug; slave.pcm "null" }' >> /usr/share/alsa/alsa.conf

# Livox library path
ENV LD_LIBRARY_PATH=${LIVOX_WS}/install/livox_ros_driver2/lib:/usr/local/lib:/usr/local/include

# Shell aliases
RUN echo 'alias tks="tmux kill-server"'                              >> ~/.bashrc && \
    echo 'alias sb="source ~/.bashrc"'                               >> ~/.bashrc && \
    echo 'alias gs="git status"'                                     >> ~/.bashrc && \
    echo 'alias gp="git push"'                                       >> ~/.bashrc && \
    echo "alias roscd='cd ${MIGHTY_WS}'"                             >> ~/.bashrc && \
    echo 'alias cb="roscd && colcon build && sb"'                    >> ~/.bashrc && \
    echo 'alias ss="roscd && source install/setup.bash"'             >> ~/.bashrc && \
    echo 'alias cbd="clear && roscd && colcon build && ss"'          >> ~/.bashrc && \
    echo 'alias cbsl="roscd && colcon build --symlink-install && sb"' >> ~/.bashrc && \
    echo 'alias cbps="roscd && colcon build --packages-select"'      >> ~/.bashrc && \
    echo 'alias tf_visualize="ros2 run rqt_tf_tree rqt_tf_tree"'     >> ~/.bashrc

# macOS visualization support
# Xpra: forwards individual windows efficiently (recommended for Mac)
# Foxglove: web-based topic viewer
RUN apt-get update && apt-get install -y --no-install-recommends \
        xpra xvfb libgl1-mesa-glx libgl1-mesa-dri \
        ros-${ROS_DISTRO}-foxglove-bridge \
        psmisc && \
    rm -rf /var/lib/apt/lists/*

# Entry point
COPY ./mighty.sh /
RUN chmod +x /mighty.sh
ENTRYPOINT ["/mighty.sh"]
