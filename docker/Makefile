# Platform detection: auto-configure for macOS vs Linux
# macOS uses Xpra for browser-based visualization (no X11 needed)
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

# Dev mode: bind-mount the local mighty source over the in-image one for fast
# parameter / launch / rviz iteration without rebuilding the image. C++ source
# changes still need a rebuild (the binaries are baked into the image).
# Usage: make run-interactive DEV=1
ifeq ($(DEV),1)
MIGHTY_SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
DEV_MOUNTS = -v $(MIGHTY_SRC):/home/kkondo/code/mighty_ws/src/mighty \
             -v $(MIGHTY_SRC)/config:/home/kkondo/code/mighty_ws/install/mighty/share/mighty/config \
             -v $(MIGHTY_SRC)/launch:/home/kkondo/code/mighty_ws/install/mighty/share/mighty/launch \
             -v $(MIGHTY_SRC)/rviz:/home/kkondo/code/mighty_ws/install/mighty/share/mighty/rviz
else
DEV_MOUNTS =
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

# Run single-agent interactive mode (manual 2D Goal Pose via RViz2)
# Usage: make run-interactive [DEV=1]
run-interactive:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e MODE=interactive \
				--rm mighty

# Run multiagent simulation with custom number of agents
# Usage: make run-multiagent NUM_AGENTS=5 [DEV=1]
run-multiagent:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e MODE=multiagent \
				-e NUM_AGENTS=$(or $(NUM_AGENTS),10) \
				--rm mighty

# Run gazebo simulation with custom goal
# Usage: make run-gazebo GOAL_X=100 GOAL_Y=50 GOAL_Z=3 [DEV=1]
run-gazebo:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e MODE=gazebo \
				-e GOAL_X=$(or $(GOAL_X),305.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),3.0) \
				-e ENV=$(or $(ENV),hard_forest) \
				--rm mighty

# Run with Xpra (macOS visualization via browser at http://localhost:8080)
# Usage: make run-mac [NUM_AGENTS=5] [DEV=1]
run-mac:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e USE_XPRA=true \
				-e MODE=multiagent \
				-e NUM_AGENTS=$(or $(NUM_AGENTS),10) \
				-p 8080:8080 \
				--rm mighty

# Run single-agent interactive mode with Xpra (manual 2D Goal Pose via RViz2)
# Usage: make run-mac-interactive [DEV=1]
run-mac-interactive:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e USE_XPRA=true \
				-e MODE=interactive \
				-p 8080:8080 \
				--rm mighty

# Run Gazebo with Xpra (macOS visualization via browser at http://localhost:8080)
# Usage: make run-mac-gazebo [GOAL_X=100 GOAL_Y=50 GOAL_Z=3] [DEV=1]
run-mac-gazebo:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
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


# ---------------------------------------------------------------------------
# Hardware (rover) image: mighty-hw — see README.hw.md. The targets above are
# the SIM image and are untouched.
# ---------------------------------------------------------------------------
HW_COMPOSE        = docker compose -f compose.hw.yaml
HW_IMAGE          = mighty-hw:local
HW_REGISTRY_IMAGE ?= ghcr.io/mit-acl/mighty-hw
HW_TAG            ?= $(shell date +%Y%m%d)

.PHONY: hw-build hw-push hw-up hw-down hw-attach hw-status
hw-build:
	$(HW_COMPOSE) build

hw-push: hw-build
	docker tag $(HW_IMAGE) $(HW_REGISTRY_IMAGE):$(HW_TAG)
	docker push $(HW_REGISTRY_IMAGE):$(HW_TAG)

hw-up:
	$(HW_COMPOSE) up -d

hw-down:
	$(HW_COMPOSE) down

hw-attach:
	docker exec -it hw-mighty tmux attach -t hw_mighty

hw-status:
	./mighty_hw.sh status
