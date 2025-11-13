# MIGHTY Interactive Demo

<a target="_blank" href="https://youtu.be/Pvb-VPUdLvg"><img src="./imgs/mighty_click_example.gif" width="700" height="386" alt="Complex Benchmarks"></a>

This is a revised version for benchmarking based on Kumar Lab's repo (https://github.com/yuwei-wu/GCOPTER.git), which implements the GCOPTER algorithm (See the "About GCOPTER" section below for more details). 

## Paper

## Video

The full video is available [https://youtu.be/Pvb-VPUdLvg](https://youtu.be/Pvb-VPUdLvg).

## System Demo

If you are interested in the real-world system demo, please check out `main` branch: [https://github.com/mit-acl/mighty/tree/main]
```
bash
git checkout main
```
and follow the instructions in the README file there.

## Setup Instructions 

### Prerequisites

Before running the code, make sure you have ROS 2 (Humble) installed. You may also need to install some additional dependencies:

```bash
sudo apt update
sudo apt install ros-humble-pcl-conversions ros-humble-pcl-ros
sudo apt install libpcl-dev libompl-dev
```

### Build

```bash
mkdir -p ws/src
cd ws/src
git clone https://github.com/kotakondo/DecompROS2.git
git clone https://github.com/kotakondo/GCOPTER.git
cd ../
colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

### Run Simulation

```bash
ros2 launch mighty mighty.launch.py
```

Now RViz should open automatically. You can click the "2D Nav Goal" button or type "g", then the cursor turns into a green arrow. Click somewhere in the RViz window to set the start position and then click the "2D Nav Goal" button (or type "g") again and click somewhere else to set the goal position. The planner will start running automatically after you set both the start and goal positions.

## About GCOPTER

__Author__: [Zhepei Wang](https://zhepeiwang.github.io) and [Fei Gao](https://scholar.google.com/citations?hl=en&user=4RObDv0AAAAJ) from [ZJU FAST Lab](http://zju-fast.com).

__Paper__: [Geometrically Constrained Trajectory Optimization for Multicopters](https://arxiv.org/abs/2103.00190), Zhepei Wang, Xin Zhou, Chao Xu, and Fei Gao, <em>[IEEE Transactions on Robotics](https://doi.org/10.1109/TRO.2022.3160022)</em> (__T-RO__), Regular Paper.
```
@article{WANG2022GCOPTER,
    title={Geometrically Constrained Trajectory Optimization for Multicopters}, 
    author={Wang, Zhepei and Zhou, Xin and Xu, Chao and Gao, Fei}, 
    journal={IEEE Transactions on Robotics}, 
    year={2022}, 
    volume={38}, 
    number={5}, 
    pages={3259-3278}, 
    doi={10.1109/TRO.2022.3160022}
}
```

## Submodules
- [SDLP: Seidel's Algorithm](https://github.com/ZJU-FAST-Lab/SDLP) on Linear-Complexity Linear Programming for Computational Geometry.
- [VertexEnumeration3D](https://github.com/ZJU-FAST-Lab/VertexEnumeration3D): Highly Efficient Vertex Enumeration for 3D Convex Polytopes (Outperforms [cddlib](https://github.com/cddlib/cddlib) in 3D).
- [LBFGS-Lite](https://github.com/ZJU-FAST-Lab/LBFGS-Lite): An Easy-to-Use Header-Only L-BFGS Solver.
