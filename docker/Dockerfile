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
        libpcl-dev build-essential nlohmann-json3-dev && \
    rm -rf /var/lib/apt/lists/*

# Workspace layout (paths must match mighty.sh / Makefile DEV_MOUNTS)
ENV CODE_DIR=/home/kkondo/code
ENV MIGHTY_WS=${CODE_DIR}/mighty_ws
ENV DECOMP_WS=${CODE_DIR}/decomp_ws
ENV LIVOX_WS=${CODE_DIR}/livox_ws
ENV LIVOX_SDK_DIR=${CODE_DIR}/Livox-SDK2

# Clone MIGHTY at a pinned release tag (override with --build-arg MIGHTY_VERSION=<tag>).
# All other dep versions live in mighty.repos — single source of truth.
#
# RELEASING — this tag and mighty.repos are read from the CLONE, not from your
# working tree, so anything committed locally but not pushed and tagged is
# invisible here. That is not hypothetical: v0.0.6 predates the `mpc` pin, so
# building against it produced an image with no ground-robot controller at all,
# and it also predates the acl-mapping bump that carries
# global_mapper_ros/cfg/sim_ground_robot.yaml. Cutting a release means:
#   1. bump MIGHTY_VERSION below IN THE COMMIT YOU ARE ABOUT TO TAG — if you tag
#      first and bump after, a checkout of the new tag still builds the old one
#   2. push the branch, then verify the clone path before committing to a tag:
#          make build MIGHTY_VERSION=<branch>
#      (git clone --branch takes a branch name, so this exercises exactly this
#      RUN without spending a tag you would then have to supersede)
#   3. merge, tag, push the tag, and re-run plain `make build`
# `make build-local` does NOT exercise any of this — it copies the working tree
# and re-imports mighty.repos from it, so it always sees your local state.
RUN mkdir -p ${MIGHTY_WS}/src
WORKDIR ${MIGHTY_WS}/src
ARG DUMMY=unknown
ARG MIGHTY_VERSION=v0.0.7
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
                 --cmake-args -DCMAKE_BUILD_TYPE=Release

# Gazebo audio fix
RUN echo 'pcm.!default { type plug; slave.pcm "null" }' >> /usr/share/alsa/alsa.conf

# Livox library path
ENV LD_LIBRARY_PATH=${LIVOX_WS}/install/livox_ros_driver2/lib:/usr/local/lib:/usr/local/include

# Register the Livox libraries with the system linker as well.
#
# LD_LIBRARY_PATH above is not sufficient: scripts/run_sim.py starts every tmux
# pane with
#     unset ... LD_LIBRARY_PATH
#     . "$SETUP_BASH"
# and mighty_ws's setup.bash does NOT chain to livox_ws — this image builds the
# workspace via CMAKE_PREFIX_PATH instead of sourcing it, so colcon never
# recorded the parent prefix (a native workspace does chain, which is why this
# only ever broke in Docker). The result: gzserver runs without
# livox_ws/install/livox_ros_driver2/lib on its search path,
# libros2_livox.so cannot resolve liblivox_ros_driver2__rosidl_typesupport_cpp.so,
# and Gazebo silently declines to load the Livox SENSOR plugin. The model still
# spawns and the RealSense MODEL plugin still works, so the only visible symptom
# is a point cloud that never publishes — hence an empty map, no frontiers, and
# a ground robot that sits at its spawn point looking healthy.
#
# ldconfig survives the unset, so this fixes it regardless of environment.
RUN echo "${LIVOX_WS}/install/livox_ros_driver2/lib" > /etc/ld.so.conf.d/livox.conf && \
    echo "/usr/local/lib"                            >> /etc/ld.so.conf.d/livox.conf && \
    ldconfig

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

# Runtime Python dependencies of the `mpc` ground-robot controller.
#
# mpc/mpc_node.py imports casadi and do_mpc at module scope, but mpc's
# package.xml declares only ROS packages, so rosdep never installs these and the
# image BUILDS FINE without them — mpc is pure Python, nothing is compiled. The
# failure is runtime-only and total:
#
#     [mpc_node-5] ModuleNotFoundError: No module named 'casadi'
#     [ERROR] [mpc_node-5]: process has died [exit code 1]
#
# with the consequence that the ground robot never receives a velocity command,
# never moves, and exploration never starts — while every other node comes up
# clean, so the terminal looks healthy. Pinned to the versions the native
# 10-run validation campaign actually ran against.
#
# python3-pip is installed here too: none of the apt lists above pull it in, so
# pip3 does not exist in this image (the first attempt at this layer failed with
# exit 127). setup.sh installs it on the native path, which is why only Docker
# was affected.
# numpy is pinned below 2 deliberately. do-mpc pulls the newest numpy (2.2.6) if
# left alone, and `import do_mpc` then dies importing matplotlib:
#
#     ImportError: numpy.core.multiarray failed to import
#
# because this image's matplotlib comes from apt and is compiled against the
# numpy 1.x C ABI. The symptom is identical to having no casadi at all — a dead
# mpc_node and a motionless robot — so it is worth keeping this pin explicit.
# 1.26.4 is what the native validation campaign ran against.
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --no-cache-dir numpy==1.26.4 casadi==3.6.7 do-mpc==5.1.1

# Entry point
# NOTE: the build context is the REPOSITORY ROOT (see docker/Makefile), not
# docker/, so that the local-source override below can see the working tree.
COPY ./docker/mighty.sh /
RUN chmod +x /mighty.sh

# ---------------------------------------------------------------------------
# Optional local-source override (contributors / CI)
# ---------------------------------------------------------------------------
# By default the image is built from the pinned MIGHTY_VERSION tag above, which
# is the path a new user takes. With --build-arg MIGHTY_SOURCE=local the working
# tree replaces the cloned package and ONLY the mighty package is rebuilt:
#
#     make build-local        # see docker/Makefile
#
# This lives at the very end on purpose. Every dependency layer (DecompROS2,
# Livox, and the tag's own mighty build) stays cached, so iterating on local C++
# costs one package rebuild instead of the full multi-hour image. Note that
# mighty.repos still comes from the cloned tag, so a local change that needs a
# *new dependency pin* requires a normal `make build` at the updated tag.
ARG MIGHTY_SOURCE=git
# Cache-bust for the local path. Without it this layer can be served from cache
# even when mighty.repos has changed, silently producing an image built against
# stale dependency pins — observed shipping mpc at an older commit than the one
# pinned in the working tree. docker/Makefile passes a timestamp.
ARG MIGHTY_LOCAL_REV=none
COPY . /tmp/mighty_local
WORKDIR ${MIGHTY_WS}
RUN if [ "${MIGHTY_SOURCE}" = "local" ]; then \
      echo "MIGHTY_SOURCE=local (rev ${MIGHTY_LOCAL_REV}): rebuilding mighty from the working tree" && \
      rm -rf ${MIGHTY_WS}/src/mighty && \
      cp -a /tmp/mighty_local ${MIGHTY_WS}/src/mighty && \
      : "cp -a preserves source mtimes, which are OLDER than the build artifacts" \
        "baked in by the tag build above. Without discarding those artifacts, make" \
        "sees everything as up to date and silently ships the TAG binary with your" \
        "local sources and configs — a build that exits 0 and changes nothing." && \
      rm -rf ${MIGHTY_WS}/build/mighty ${MIGHTY_WS}/install/mighty && \
      : "Import WITHOUT --skip-existing. The local tree may both add dependencies" \
        "the tag lacks (mpc) and bump the pin of one the tag already cloned. This" \
        "path originally used --skip-existing to stay cheap, which silently kept" \
        "acl-mapping at the tag's older commit: global_mapper_ros then shipped" \
        "without cfg/sim_ground_robot.yaml, the mapper came up with no parameter" \
        "file, published no occupancy grid, and the ground robot never moved." \
        "Updating existing checkouts to the pinned versions is the whole point of" \
        "building from the local tree, so pay the extra fetch." && \
      vcs import ${MIGHTY_WS}/src < ${MIGHTY_WS}/src/mighty/mighty.repos && \
      : "vcs re-clones the three deps that were MOVED to their own workspaces" \
        "above (they are absent from mighty_ws/src, so --skip-existing does not" \
        "skip them). Drop the copies again, or colcon tries to build them here" \
        "and fails on the half-imported trees." && \
      rm -rf ${MIGHTY_WS}/src/DecompROS2 ${MIGHTY_WS}/src/livox_ros_driver2 \
             ${MIGHTY_WS}/src/Livox-SDK2 && \
      source /opt/ros/humble/setup.bash && \
      source ${DECOMP_WS}/install/setup.bash && \
      export CMAKE_PREFIX_PATH=${LIVOX_WS}/install/livox_ros_driver2:${DECOMP_WS}/install/decomp_util:${CMAKE_PREFIX_PATH} && \
      export MAKEFLAGS="-j${JOBS:-$(nproc)}" && \
      colcon build ${JOBS:+--parallel-workers $JOBS} \
                   --cmake-args -DCMAKE_BUILD_TYPE=Release ; \
    else \
      echo "MIGHTY_SOURCE=git: keeping the ${MIGHTY_VERSION} build" ; \
    fi && \
    rm -rf /tmp/mighty_local

ENTRYPOINT ["/mighty.sh"]
