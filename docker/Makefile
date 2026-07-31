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

# Override the git ref the image is built from. Accepts anything
# `git clone --branch` takes, so a BRANCH works as well as a tag:
#     make build MIGHTY_VERSION=dev-sim
# That is the way to verify the clone-from-GitHub path before cutting a tag —
# `make build-local` cannot, since it copies the working tree instead of
# cloning. Empty (the default) uses the pinned tag in the Dockerfile.
ifneq ($(MIGHTY_VERSION),)
VERSION_ARG = --build-arg MIGHTY_VERSION=$(MIGHTY_VERSION)
else
VERSION_ARG =
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

# The build context is the REPOSITORY ROOT (..), not docker/, so that
# `build-local` can copy the working tree into the image. .dockerignore at the
# repo root keeps the context small (excludes .git, build/, install/, log/).
.PHONY: build
build:
	docker build \
	  $(PLATFORM_FLAG) \
	  $(BUILD_ARGS) \
	  $(VERSION_ARG) \
	  -f Dockerfile -t mighty ..

.PHONY: build-no-cache
build-no-cache:
	docker build \
	$(PLATFORM_FLAG) \
	$(BUILD_ARGS) \
	$(VERSION_ARG) \
	-f Dockerfile \
	--build-arg DUMMY=`date +%s` -t mighty ..

# Build the image from the LOCAL working tree instead of the pinned release tag.
# Requires a prior `make build` (or any cached image layers): every dependency
# layer is reused and only the mighty package is recompiled, so this takes
# minutes rather than hours. Use it to test uncommitted C++ changes in Docker —
# DEV=1 cannot, because it only bind-mounts config/launch/rviz over an image
# whose binaries were compiled from the tag.
# Usage: make build-local [JOBS=8]
.PHONY: build-local
build-local:
	docker build \
	  $(PLATFORM_FLAG) \
	  $(BUILD_ARGS) \
	  -f Dockerfile \
	  --build-arg MIGHTY_SOURCE=local \
	  --build-arg MIGHTY_LOCAL_REV=`date +%s` -t mighty ..

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

# Run single-agent ground robot (Pioneer 3-AT) autonomous frontier exploration
# in the ACL_office world. The robot self-drives — no goal needed.
# Usage: make run-ground-exploration [ENV=ACL_office] [DEV=1]
run-ground-exploration:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e MODE=exploration-singleagent-ground \
				-e ENV=$(or $(ENV),ACL_office) \
				-e NO_RVIZ=$(or $(NO_RVIZ),false) \
				--rm mighty

# Run single-agent ground robot (Pioneer 3-AT) crossing a forest world to a
# fixed goal — the ground-robot counterpart of run-gazebo. Exploration is turned
# off at launch, so GOAL_X/GOAL_Y/GOAL_Z are honoured even though
# config/mighty_ground_robot.yaml enables exploration. For frontier exploration
# use run-ground-exploration instead (ENV selects the world there too).
# Usage: make run-ground-robot [GOAL_X=30 GOAL_Y=0 GOAL_Z=0] [ENV=easy_forest] [DEV=1]
run-ground-robot:
	$(XHOST_CMD) && docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(DISPLAY_ENV) \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e MODE=gazebo \
				-e GROUND_ROBOT=true \
				-e GOAL_X=$(or $(GOAL_X),30.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),0.0) \
				-e ENV=$(or $(ENV),easy_forest) \
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

# Run ground robot autonomous exploration on Mac (Xpra, browser at localhost:8080)
# Usage: make run-mac-ground-exploration [ENV=ACL_office] [DEV=1]
run-mac-ground-exploration:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e USE_XPRA=true \
				-e MODE=exploration-singleagent-ground \
				-e ENV=$(or $(ENV),ACL_office) \
				-p 8080:8080 \
				--rm mighty

# Run ground robot in a forest world on Mac (Xpra, browser at localhost:8080)
# See run-ground-robot notes re: exploration.enabled and goal navigation.
# Usage: make run-mac-ground-robot [GOAL_X=30 GOAL_Y=0 GOAL_Z=0] [ENV=easy_forest] [DEV=1]
run-mac-ground-robot:
	docker run $(PLATFORM_FLAG) -it \
				$(GPU_FLAGS) \
				--privileged \
				$(NETWORK_FLAG) \
				--env="QT_X11_NO_MITSHM=1" \
				$(DEV_MOUNTS) \
				-e USE_XPRA=true \
				-e MODE=gazebo \
				-e GROUND_ROBOT=true \
				-e GOAL_X=$(or $(GOAL_X),30.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),0.0) \
				-e ENV=$(or $(ENV),easy_forest) \
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
