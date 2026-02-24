# Platform detection: auto-configure for macOS vs Linux
# macOS uses Foxglove Studio for visualization (native app, no X11 needed)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  PLATFORM_FLAG = --platform linux/amd64
  GPU ?= false
  JOBS ?= 4
  XHOST_CMD = /opt/X11/bin/xhost +localhost
  DISPLAY_ENV = -e DISPLAY=host.docker.internal:0 -e LIBGL_ALWAYS_SOFTWARE=1 -e GALLIUM_DRIVER=llvmpipe
  NETWORK_FLAG =
else
  PLATFORM_FLAG =
  GPU ?= true
  JOBS ?=
  XHOST_CMD = xhost +
  DISPLAY_ENV = -e DISPLAY=${DISPLAY}
  NETWORK_FLAG = --network=host
endif

# GPU support: set GPU=false to run without NVIDIA GPU
# Usage: make run GPU=false
ifeq ($(GPU),true)
GPU_FLAGS = --gpus all --env="NVIDIA_DRIVER_CAPABILITIES=all"
else
GPU_FLAGS =
endif

# Mount local mighty source into container (picks up script/config/launch changes without rebuild)
# Source mount: for run_sim.py and other scripts called directly from src/
# Config mount: overlays installed config so ROS nodes pick up local mighty.yaml changes
MIGHTY_SRC_MOUNT = -v $(CURDIR)/..:/home/kkondo/code/mighty_ws/src/mighty \
		   -v $(CURDIR)/../config:/home/kkondo/code/mighty_ws/install/mighty/share/mighty/config

# Build args for parallelism (empty = unlimited, set JOBS=4 to limit)
ifneq ($(JOBS),)
BUILD_ARGS = --build-arg JOBS=$(JOBS)
else
BUILD_ARGS =
endif

.PHONY: build
build:
	docker build \
	  $(PLATFORM_FLAG) \
	  $(BUILD_ARGS) \
	  -t mighty .

.PHONY: build-no-cache
build-no-cache:
	docker build \
	$(PLATFORM_FLAG) \
	$(BUILD_ARGS) \
	--build-arg DUMMY=`date +%s` -t mighty .

# Incremental rebuild of just the mighty package inside Docker (much faster than full build)
# Mounts local source, runs colcon build --packages-select mighty, commits result back to image
.PHONY: build-mighty
build-mighty:
	docker run $(PLATFORM_FLAG) --name mighty_rebuild \
				$(MIGHTY_SRC_MOUNT) \
				--entrypoint /bin/bash \
				mighty -c \
				'. /opt/ros/humble/setup.sh && \
				 . /home/kkondo/code/decomp_ws/install/setup.sh && \
				 export CMAKE_PREFIX_PATH=/home/kkondo/code/livox_ws/install/livox_ros_driver2:/home/kkondo/code/decomp_ws/install/decomp_util && \
				 cd /home/kkondo/code/mighty_ws && \
				 colcon build --packages-select mighty --cmake-args -DCMAKE_BUILD_TYPE=Release'
	docker commit --change='ENTRYPOINT ["/mighty.sh"]' mighty_rebuild mighty
	docker rm mighty_rebuild

# Default run (multiagent with 10 agents)
run:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=multiagent \
				-e NUM_AGENTS=10 \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Run multiagent simulation with custom number of agents
# Usage: make run-multiagent NUM_AGENTS=5
run-multiagent:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=multiagent \
				-e NUM_AGENTS=$(or $(NUM_AGENTS),10) \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Run gazebo simulation with custom goal
# Usage: make run-gazebo GOAL_X=100 GOAL_Y=50 GOAL_Z=3
run-gazebo:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=gazebo \
				-e GOAL_X=$(or $(GOAL_X),305.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),3.0) \
				-e ENV=$(or $(ENV),hard_forest) \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Run with Xpra (macOS visualization via browser at http://localhost:8080)
# Usage: make run-mac [NUM_AGENTS=5]
run-mac:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e USE_XPRA=true \
				-e MODE=multiagent \
				-e NUM_AGENTS=$(or $(NUM_AGENTS),10) \
				-p 8080:8080 \
				-v $(CURDIR)/mighty.sh:/mighty.sh \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Run single-agent interactive mode with Xpra (manual 2D Goal Pose via RViz2)
# Usage: make run-mac-interactive
run-mac-interactive:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e USE_XPRA=true \
				-e MODE=multiagent \
				-e NUM_AGENTS=1 \
				-e NO_AUTO_GOAL=true \
				-p 8080:8080 \
				-v $(CURDIR)/mighty.sh:/mighty.sh \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Run Gazebo with Xpra (macOS visualization via browser at http://localhost:8080)
# Usage: make run-mac-gazebo [GOAL_X=100 GOAL_Y=50 GOAL_Z=3]
run-mac-gazebo:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				-e USE_XPRA=true \
				-e MODE=gazebo \
				-e GOAL_X=$(or $(GOAL_X),305.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),3.0) \
				-e ENV=$(or $(ENV),hard_forest) \
				-p 8080:8080 \
				-v $(CURDIR)/mighty.sh:/mighty.sh \
				$(MIGHTY_SRC_MOUNT) \
				--rm mighty

# Interactive shell (for debugging)
shell:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(MIGHTY_SRC_MOUNT) \
				--rm \
				--entrypoint /bin/bash \
				mighty
