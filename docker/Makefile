# Platform detection: auto-configure for macOS vs Linux
# macOS uses Foxglove Studio for visualization (native app, no X11 needed)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  PLATFORM_FLAG = --platform linux/amd64
  GPU ?= false
  JOBS ?= 2
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
				--rm mighty

# Interactive shell (for debugging)
shell:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				--rm \
				--entrypoint /bin/bash \
				mighty
