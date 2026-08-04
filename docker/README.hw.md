# hw-mighty container (hardware stack)

The MIGHTY hardware autonomy stack (planner + MPC + global mapper + DLIO)
packaged as **one idling container** driven by **one script**: `mighty_hw.sh`
builds a HOST tmux session `hw_mighty` in which every pane `docker exec`s its
node into the container. Per-pane Ctrl-C / Up / Enter restarts a single node;
Ctrl-b d detaches and the stack keeps running.

Sources are **baked from the host checkouts** at their host paths
(`/home/swarm/code/...`): the build context is `/home/swarm/code` because the
image bakes four sibling trees — `mighty_ws/src`, `decomp_ws/src`,
`livox_ws/src`, `Livox-SDK2`. See `Dockerfile.hw.dockerignore` for exactly what
ships. The session builder lives on the host, so pane/flag changes need **no
image rebuild** — only source changes do (`make hw-build` or
`./mighty_hw.sh rebuild`, ~13 s warm thanks to the manifest cache pin).

## Usage

```bash
./mighty_hw.sh start [--odom-type dlio|dlio_in_mocap|mocap] [--two-d-only|--no-two-d]
./mighty_hw.sh attach      # or: tmux attach -t hw_mighty
./mighty_hw.sh status
./mighty_hw.sh stop
./mighty_hw.sh rebuild     # stop + image rebuild + start (start flags forwarded)
```

## Prerequisites (every mode)

- **`drive.service`** — hosts the zenoh router on `:7447`; every pane spin-waits
  on it, and `start` fails fast if it is absent.
- **`sensors.service`** — livox MID-360 + D455 + the `base_link->lidar` static
  TF. This stack launches **no sensors, ever** (the old `--no-sensors` flag is
  gone because it is the only behavior).
- **`dlio.service` must be disabled** — this stack runs its own DLIO panes (the
  same nodes; both running = double-publish). Rollback: re-enable
  `dlio.service`, use the native `mighty` alias (`run_hw_red_rover.py`).

## Odom types

| `--odom-type` | Panes beyond MIGHTY / RViz 2D goal / Global Mapper | Notes |
|---|---|---|
| `dlio` (default) | DLIO, seed pose, `tf map->odom` | Replicates what `dlio.service` ran: DLIO is told `initial_pose_topic:=world` and a 2 Hz constant `PoseStamped` spoof on `/<ns>/world` (z = `DLIO_SEED_Z`, default 0.4) anchors it. Works with no external infrastructure. |
| `dlio_in_mocap` | DLIO, `tf map->odom` | A **real** mocap publishes `/<ns>/world`. No spoof pane — DLIO anchors to the **first** pose it receives, so a spoof would race the mocap and win. |
| `mocap` | `tf mocap->base_link`, `tf world->map`, `tf map->odom` | No DLIO. Mighty runs `use_onboard_localization:=false` with `twist_topic:=mocap/twist`; the mapper consumes **raw** `livox/lidar` with `pose_topic:=world`. |

In the dlio modes the mapper consumes DLIO's deskewed cloud
(`dlio/odom_node/pointcloud/deskewed`) and pose (`dlio/odom_node/pose`).
DLIO runs a ~3 s stationary IMU calibration whenever its pane (re)starts — keep
the rover still. `--two-d-only` (default ON) pins DLIO's published z to the
seed; `--no-two-d` disables that.

## Env & networking

- Identity/transport come from **`/etc/rover/rover.env`** via compose
  `env_file` (`ROBOT_NAME`, `VEHTYPE`/`VEHNUM`, `ROS_DOMAIN_ID`,
  `RMW_IMPLEMENTATION`, optional `DLIO_SEED_Z`) — same as the
  drive/sensors/dlio siblings, so this directory is rover-agnostic.
- The container runs the mighty **nodes only** over `network_mode: host`,
  attaching to the host zenoh router at `localhost:7447` (owned by
  `drive.service`). No zenoh config is mounted; `ZENOH_ROUTER_CONFIG_URI` from
  rover.env is ignored by plain sessions.

## Startup ordering (/tf_static race)

The static TF panes latch one transient-local sample; over rmw_zenoh a
subscriber that forms before the publisher exists never receives it. Consumers
therefore gate themselves with `ros2 run mighty wait_for_tf.py <pairs>` before
launching (mighty waits for `map->odom`; the mapper additionally for
`base_link->lidar`, plus `world->map` in mocap mode). `wait_for_tf.py` warns
and continues on timeout, so the stack never deadlocks.
