FROM ros:humble

LABEL Maintainer="Kota Kondo <kkondo@mit.edu>" \
      Description="Ubuntu 22.04, ROS2 humble for MIGHTY"

# Set environment variables
##########################################
ENV ROS_DISTRO=humble
ENV DEBIAN_FRONTEND=noninteractive

# Usuals
##########################################
RUN rm -rf /var/lib/apt/lists/*
RUN apt update
RUN apt upgrade -y
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -q -y --no-install-recommends git tmux vim wget tmuxp make gdb

# Making directories and setting up the environment
##########################################
WORKDIR /home/kkondo
RUN mkdir code

# Install dependencies
##########################################
RUN apt-get install -y  ros-${ROS_DISTRO}-gazebo-* \
                        ros-${ROS_DISTRO}-pcl-conversions \
                        ros-${ROS_DISTRO}-example-interfaces \
                        ros-${ROS_DISTRO}-pcl-ros \
                        ros-${ROS_DISTRO}-rviz2 \
                        ros-${ROS_DISTRO}-rqt-gui \
                        ros-${ROS_DISTRO}-rqt-gui-py \
                        ros-${ROS_DISTRO}-tf2-tools \
                        ros-${ROS_DISTRO}-tf-transformations

RUN apt install -y  ros-${ROS_DISTRO}-desktop \
                        ros-dev-tools \
                        ros-${ROS_DISTRO}-turtlesim \
                        ros-${ROS_DISTRO}-rqt* \
                        ros-${ROS_DISTRO}-rviz2 \
                        ros-${ROS_DISTRO}-rviz-common \ 
                        libpcl-dev \
                        build-essential \
                        nlohmann-json3-dev

# Install MIGHTY and dependencies
##########################################
# Version control: commits from mighty.repos
ARG MIGHTY_VERSION=6237c6524df309e50e34b5501cbef70b2e1c6d80
ARG ACL_MAPPING_VERSION=c3b27af1549ef1ab549acbe1143bf810a900da95
ARG DYNUS_INTERFACES_VERSION=ec381d5efeaaa5d25c73a3453d4dc78d59a6594c
ARG GAZEBO_ROS_PKGS_VERSION=ed1f7b886971d7bf898e6e0f59c622df929c0521
ARG LIVOX_LASER_SIM_VERSION=246040d802ec2d7374f38ffd07e5e7fc49e297fb
ARG REALSENSE_GAZEBO_VERSION=8c357e69da2160f54328db70915bd9b28eb975bf
ARG UAV_SIMULATOR_VERSION=b196208d3cbf0269b8e80220f7cc1cc992b6c7f8
ARG DECOMP_ROS2_VERSION=a2bc42ed46bf934323994706f17af8a3883860ea
ARG LIVOX_ROS_DRIVER2_VERSION=274aad204b5dfce15c648e88dfd375fa5d7b3026
ARG LIVOX_SDK2_VERSION=6a940156dd7151c3ab6a52442d86bc83613bd11b

RUN mkdir -p /home/kkondo/code/mighty_ws/src
WORKDIR /home/kkondo/code/mighty_ws/src
# dummy arg to invalidate cache and build from here
ARG DUMMY=unknown
RUN DUMMY=${DUMMY} git clone https://github.com/mit-acl/mighty.git
WORKDIR /home/kkondo/code/mighty_ws/src/mighty
RUN git checkout ${MIGHTY_VERSION}
WORKDIR /home/kkondo/code/mighty_ws/src
RUN git clone https://github.com/kotakondo/dynus_interfaces.git && \
    cd dynus_interfaces && git checkout ${DYNUS_INTERFACES_VERSION}
RUN git clone https://github.com/kotakondo/realsense_gazebo_plugin.git && \
    cd realsense_gazebo_plugin && git checkout ${REALSENSE_GAZEBO_VERSION}
RUN git clone https://github.com/kotakondo/livox_laser_simulation_ros2.git && \
    cd livox_laser_simulation_ros2 && git checkout ${LIVOX_LASER_SIM_VERSION}
RUN git clone https://github.com/kotakondo/gazebo_ros_pkgs.git && \
    cd gazebo_ros_pkgs && git checkout ${GAZEBO_ROS_PKGS_VERSION}
RUN git clone https://gitlab.com/mit-acl/lab/acl-mapping.git && \
    cd acl-mapping && git checkout ${ACL_MAPPING_VERSION}
RUN git clone https://github.com/kotakondo/uav_simulator.git && \
    cd uav_simulator && git checkout ${UAV_SIMULATOR_VERSION}
RUN mkdir -p /home/kkondo/code/decomp_ws/src
WORKDIR /home/kkondo/code/decomp_ws/src
RUN git clone https://github.com/kotakondo/DecompROS2.git && \
    cd DecompROS2 && git checkout ${DECOMP_ROS2_VERSION}
RUN mkdir -p /home/kkondo/code/livox_ws/src
WORKDIR /home/kkondo/code/livox_ws/src
RUN git clone https://github.com/kotakondo/livox_ros_driver2.git && \
    cd livox_ros_driver2 && git checkout ${LIVOX_ROS_DRIVER2_VERSION}
WORKDIR /home/kkondo/code
RUN git clone https://github.com/Livox-SDK/Livox-SDK2.git && \
    cd Livox-SDK2 && git checkout ${LIVOX_SDK2_VERSION}

# Build the workspace
##########################################
### for decomp
WORKDIR /home/kkondo/code/decomp_ws
RUN . /opt/ros/humble/setup.sh && colcon build --packages-select decomp_util
RUN . /home/kkondo/code/decomp_ws/install/setup.sh && colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

### for Livox-SDK2
WORKDIR /home/kkondo/code/Livox-SDK2
RUN mkdir build
WORKDIR /home/kkondo/code/Livox-SDK2/build
RUN cmake .. && make -j && sudo make install

### for livox_ros_driver2
WORKDIR /home/kkondo/code/livox_ws/src/livox_ros_driver2
RUN . /opt/ros/humble/setup.sh && ./build.sh humble

### for other packages
WORKDIR /home/kkondo/code/mighty_ws
RUN . /opt/ros/humble/setup.sh && \
    . /home/kkondo/code/decomp_ws/install/setup.sh && \
    export CMAKE_PREFIX_PATH=/home/kkondo/code/livox_ws/install/livox_ros_driver2:/home/kkondo/code/decomp_ws/install/decomp_util && colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

# Gazebo
##########################################
### handle ALSA-related error
RUN echo    'pcm.!default { \
                  type plug \
                  slave.pcm "null" \
            }' >> /usr/share/alsa/alsa.conf

# For livox lidar
ENV LD_LIBRARY_PATH=/home/kkondo/code/livox_ws/install/livox_ros_driver2/lib:${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/include

# Alias
RUN echo 'alias tks="tmux kill-server"' >> ~/.bashrc
RUN echo 'alias sb="source ~/.bashrc"' >> ~/.bashrc
RUN echo 'alias eb="code ~/.bashrc"' >> ~/.bashrc
RUN echo 'alias gs="git status"' >> ~/.bashrc
RUN echo 'alias gp="git push"' >> ~/.bashrc
RUN echo 'alias roscd="cd ~/code/mighty_ws"' >> ~/.bashrc
RUN echo 'alias cb="roscd && colcon build && sb"' >> ~/.bashrc
RUN echo 'alias ss="roscd && source install/setup.bash"' >> ~/.bashrc
RUN echo 'alias cbd="clear && roscd && colcon build && ss"' >> ~/.bashrc
RUN echo 'alias cbm="clear && roscd && colcon build --packages-select ros2_mapper && ss"' >> ~/.bashrc
RUN echo 'alias cbsl="roscd && colcon build --symlink-install && sb"' >> ~/.bashrc
RUN echo 'alias cbps="roscd && colcon build --packages-select"' >> ~/.bashrc
RUN echo 'alias tf_visualize="ros2 run rqt_tf_tree rqt_tf_tree"' >> ~/.bashrc
RUN echo 'alias tks="tmux kill-server"' >> ~/.bashrc

# Set up Entry Point
# ##########################################
COPY ./mighty.sh /
RUN chmod +x /mighty.sh
ENTRYPOINT ["/mighty.sh"]
