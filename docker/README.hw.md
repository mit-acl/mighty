# hw-mighty container (hardware stack)

The MIGHTY hardware autonomy stack (planner + MPC + global mapper + DLIO +
Livox driver) packaged as **one container** running the same attachable 7-pane
tmux `hw_mighty` session inside it. Sibling of the `rover-drive` container
(rfc0001 `rover_bringup/docker/`) and follows its conventions: pinned,
reproducible image; identical observability — attach, watch each node live,
Ctrl-C / Up / Enter a pane to restart just that node.

Sources are **baked from the host checkouts** at their host paths
(`/home/swarm/code/...`), so `run_hw_red_rover.py` runs unmodified and the
image captures exactly what the rover runs natively (build context is
`/home/swarm/code`; see `Dockerfile.hw.dockerignore` for what ships). Only
`/etc/rover` is mounted at runtime. Rover identity (RR06/RR08/...) comes purely
from `/etc/rover/rover.env` — the same image runs on any rover.

The zenoh router is **not** ours: rover-drive owns `rmw_zenohd` on `:7447`;
mighty's nodes attach to it over `network_mode: host`. The Livox MID360
(UDP, host `192.168.1.100`) also rides on host networking — no devices needed.

## Operate (manual — the native `mighty` alias is untouched and is the rollback)

```bash
cd ~/code/mighty_ws/src/mighty/docker

./mighty_hw.sh start                 # compose up -d + attach (detach: Ctrl-b d)
./mighty_hw.sh start --odom-type dlio_in_mocap
./mighty_hw.sh start --no-two-d      # disable DLIO 2D-only z pinning
./mighty_hw.sh attach                # re-attach to a running stack
./mighty_hw.sh status                # container + pane status
./mighty_hw.sh stop                  # tear the whole stack down
docker logs hw-mighty                # PID-1 wrapper output
```

Pane intervention is identical to native: Ctrl-C, ↑, Enter in any pane
restarts that node. Detaching never stops the stack; the container exits when
the `hw_mighty` session dies (or on `stop`).

## Build / update

```bash
make hw-build        # bake current host checkouts into mighty-hw:local
make hw-push         # retag + push to ghcr.io/mit-acl/mighty-hw (needs ghcr login)
```

Layer caching: apt/pip and the dependency workspaces (decomp, Livox SDK +
driver) only rebuild when their sources change; a mighty/mpc edit re-runs just
the final colcon layer.

## Dev loop (edit on the rover without rebuilding the image)

```bash
./mighty_hw.sh start --dev
```

Mounts host `src/mighty` + `src/mpc` read-only, colcon-builds them inside the
container (into named volumes; incremental after the first run), then
launches. Config/launch/YAML edits land immediately this way. Edits to other
packages (acl-mapping, dlio, livox) need `make hw-build` instead.

## Rollback to native

The native flow is untouched — just use the `mighty` alias (or
`scripts/run_mighty_dlio.sh`) after `./mighty_hw.sh stop`. Don't run both at
once: two stacks would double-publish on the same topics.

## Notes

- `restart: "no"` — a dead node stays dead until hand-restarted in its pane
  (native semantics). The container itself only dies with the session.
- No `tty: true` in compose — a detached container pty reports width 0 and
  newer libtmux fails with `new-session: width too small`; the launch command
  pins `COLUMNS`/`LINES` instead (gotcha inherited from rover-drive).
- Image builds with `-DBUILD_SIMULATION=OFF` (no Gazebo); sim work stays on
  the existing sim image (`Dockerfile` + `Makefile` targets in this dir).
