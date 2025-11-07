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
RUN mkdir -p /home/kkondo/code/mighty_ws/src
WORKDIR /home/kkondo/code/mighty_ws/src
# dummy arg to invalidate cache and build from here
ARG DUMMY=unknown
RUN DUMMY=${DUMMY} git clone https://github.com/mit-acl/mighty.git
WORKDIR /home/kkondo/code/mighty_ws/src/mighty
RUN git switch dev
WORKDIR /home/kkondo/code/mighty_ws/src
RUN git clone https://github.com/kotakondo/dynus_interfaces.git
RUN git clone https://github.com/kotakondo/realsense_gazebo_plugin.git
RUN git clone https://github.com/kotakondo/livox_laser_simulation_ros2.git
RUN git clone https://github.com/kotakondo/gazebo_ros_pkgs.git
RUN git clone https://gitlab.com/mit-acl/lab/acl-mapping.git
WORKDIR /home/kkondo/code/mighty_ws/src/acl-mapping
RUN git switch ros2
RUN mkdir -p /home/kkondo/code/decomp_ws/src
WORKDIR /home/kkondo/code/decomp_ws/src
RUN git clone https://github.com/kotakondo/DecompROS2.git
RUN mkdir -p /home/kkondo/code/livox_ws/src
WORKDIR /home/kkondo/code/livox_ws/src
RUN git clone https://github.com/kotakondo/livox_ros_driver2.git
WORKDIR /home/kkondo/code
RUN git clone https://github.com/Livox-SDK/Livox-SDK2.git

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
    export CMAKE_PREFIX_PATH=/home/kkondo/code/livox_ws/install/livox_ros_driver2:/home/kkondo/code/decomp_ws/install/decomp_util && \
    colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release

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
