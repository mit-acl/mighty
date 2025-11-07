# at the top of your Makefile
.PHONY: build
build:
	docker build \
	  -t mighty .

.PHONY: build-no-cache
build-no-cache:
	docker build \
	--build-arg DUMMY=`date +%s` -t mighty .

run:
	xhost + && docker run  -it \
				--gpus all \
				--env="NVIDIA_DRIVER_CAPABILITIES=all" \
				--privileged \
				-e DISPLAY=${DISPLAY} \
				--network=host \
				--env="DISPLAY" \
				--env="QT_X11_NO_MITSHM=1" \
				--rm mighty