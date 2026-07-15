# hw-mighty container (hardware stack) — RR04

The MIGHTY hardware autonomy stack (planner + MPC + global mapper + DLIO + Livox
driver) packaged as **one container** running the same attachable 7-pane tmux
`hw_mighty` session inside it. It runs as closely as possible to the native
`mighty` alias on RR04: same launcher
(`run_hw_red_rover.py --odom-type dlio_in_mocap --two-d-only`), same 7 panes,
per-pane Ctrl-C / Up / Enter restarts a single node.

Sources are **baked from the host checkouts** at their host paths
(`/home/swarm/code/...`), so `run_hw_red_rover.py` (which hardcodes those paths)
runs unmodified. The build context is `/home/swarm/code` because the image bakes
four trees that are siblings — `mighty_ws/src`, `decomp_ws/src`, `livox_ws/src`,
`Livox-SDK2` — not all under `mighty_ws`. See `Dockerfile.hw.dockerignore` for
exactly what ships.

## Env & networking (RR04 specifics)

- **No `/etc/rover/rover.env` on RR04** — identity/transport are set directly in
  `compose.hw.yaml`'s `environment:` block: `ROVER_NAME=RR04`, `VEHTYPE=RR`,
  `VEHNUM=04`, `ROS_DOMAIN_ID=22`, `RMW_IMPLEMENTATION=rmw_zenoh_cpp` (the native
  values, from `~/config/roman_rover.sh` + `RR04.sh`). The image itself is
  rover-agnostic; identity lives in the compose file.
- **Zenoh:** the container runs the mighty **nodes only** and attaches to the
  **host** Zenoh router over `network_mode: host`. **Start the host router first**
  (native workflow): run the `zenoh_route` alias
  (`ros2 run rmw_zenoh_cpp rmw_zenohd`) before starting the container. No config
  is mounted — the nodes use the default router at `localhost:7447`, shared via
  host networking.

## Operate (manual — the native `mighty` alias is untouched and is the rollback)

```bash
cd ~/code/mighty_ws/src/mighty/docker

./mighty_hw.sh start                 # compose up -d + attach (detach: Ctrl-b d)
./mighty_hw.sh start --odom-type mocap
./mighty_hw.sh start --no-two-d      # disable DLIO 2D-only z pinning
./mighty_hw.sh attach                # re-attach to a running stack
./mighty_hw.sh status                # container + pane status
./mighty_hw.sh stop                  # tear the whole stack down
docker logs hw-mighty                # PID-1 wrapper output
```

Default run is `--odom-type dlio_in_mocap --two-d-only` (matches the `mighty`
alias). Pane intervention is identical to native: Ctrl-C, ↑, Enter in any pane
restarts that node. Detaching never stops the stack; the container exits when the
`hw_mighty` session dies (or on `stop`).

## Build

```bash
make hw-build        # docker compose -f compose.hw.yaml build → mighty-hw:local
```

Layer caching: apt/pip and the dependency workspaces (decomp, Livox SDK + driver)
only rebuild when their sources change; a mighty/mpc edit re-runs just the final
colcon layer. The image is built with `-DBUILD_SIMULATION=OFF` (no Gazebo).

Note: RR04's `mighty_ws/src` still contains a **duplicate `DecompROS2`** (decomp is
built separately from `decomp_ws/src`); it is excluded in
`Dockerfile.hw.dockerignore` to avoid a double-build/conflict.

## Rollback to native

The native flow is untouched — after `./mighty_hw.sh stop`, just use the `mighty`
alias. Don't run both at once: two stacks would double-publish on the same topics.

## Notes

- `restart: "no"` — a dead node stays dead until hand-restarted in its pane
  (native semantics). The container itself only dies with the session.
- No `tty: true` in compose — a detached container pty reports width 0 and newer
  libtmux fails with `new-session: width too small`; the launch command pins
  `COLUMNS`/`LINES` instead.
- **No bind-mounts and no dev-override compose file** — the image is the single
  source of truth (control/repeatability). To pick up code edits, rebuild
  (`make hw-build`).
