# GPU support: set GPU=false to run without NVIDIA GPU
# Usage: make run GPU=false
GPU ?= true
ifeq ($(GPU),true)
GPU_FLAGS = --gpus all --env="NVIDIA_DRIVER_CAPABILITIES=all"
else
GPU_FLAGS =
endif

.PHONY: build
build:
	docker build \
	  -t mighty .

.PHONY: build-no-cache
build-no-cache:
	docker build \
	--build-arg DUMMY=`date +%s` -t mighty .

# Default run (multiagent with 10 agents)
run:
	xhost + && docker run -it \
				$(GPU_FLAGS) \
				--privileged \
				-e DISPLAY=${DISPLAY} \
				--network=host \
				--env="DISPLAY" \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=multiagent \
				-e NUM_AGENTS=10 \
				--rm mighty

# Run multiagent simulation with custom number of agents
# Usage: make run-multiagent NUM_AGENTS=5
run-multiagent:
	xhost + && docker run -it \
				$(GPU_FLAGS) \
				--privileged \
				-e DISPLAY=${DISPLAY} \
				--network=host \
				--env="DISPLAY" \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=multiagent \
				-e NUM_AGENTS=$(or $(NUM_AGENTS),10) \
				--rm mighty

# Run gazebo simulation with custom goal
# Usage: make run-gazebo GOAL_X=100 GOAL_Y=50 GOAL_Z=3
run-gazebo:
	xhost + && docker run -it \
				$(GPU_FLAGS) \
				--privileged \
				-e DISPLAY=${DISPLAY} \
				--network=host \
				--env="DISPLAY" \
				--env="QT_X11_NO_MITSHM=1" \
				-e MODE=gazebo \
				-e GOAL_X=$(or $(GOAL_X),305.0) \
				-e GOAL_Y=$(or $(GOAL_Y),0.0) \
				-e GOAL_Z=$(or $(GOAL_Z),3.0) \
				-e ENV=$(or $(ENV),hard_forest) \
				--rm mighty

# Interactive shell (for debugging)
shell:
	xhost + && docker run -it \
				$(GPU_FLAGS) \
				--privileged \
				-e DISPLAY=${DISPLAY} \
				--network=host \
				--env="DISPLAY" \
				--env="QT_X11_NO_MITSHM=1" \
				--rm \
				--entrypoint /bin/bash \
				mighty
