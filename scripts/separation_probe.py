#!/usr/bin/env python3
"""Min-separation probe for the multi-robot ground swap sims.

The verifier behind the swap demos' >=1 m pairwise-separation guarantee.
Subscribes to every /NXxx/state, computes the pairwise min distance at 20 Hz,
tracks below-threshold episodes (with 0.5 m / 0.3 m severity tiers), and
counts goal legs + true arrivals so a "safe" run that secretly gridlocked is
visible (separation via parking is not success). Prints STATUS lines every
30 s, VIOL lines per below-threshold episode, a final single-line SUMMARY
json (also written to <outdir>/sep_<tag>.json), and a 4 Hz min-sep timeline
CSV for locating events.

Run it next to a live sim (same host / same ROS_DOMAIN_ID; inside the sim's
container for docker runs — separate containers silently drop FastDDS SHM
traffic):
  python3 separation_probe.py --agents 10 --threshold 1.0 --duration 900 \
      --tag myrun --outdir /tmp
"""
import argparse
import json
import math
import time

import rclpy
from rclpy.node import Node

from dynus_interfaces.msg import State
from geometry_msgs.msg import PoseStamped


class SepProbe(Node):
    def __init__(self, args):
        super().__init__('sep_probe')
        self.args = args
        self.names = [f'{args.prefix}{i:02d}' for i in range(1, args.agents + 1)]
        self.pos = {}            # ns -> (x, y, t_mono)
        self.goal = {}           # ns -> (x, y)
        self.reached_current = {ns: False for ns in self.names}
        self.legs = {ns: 0 for ns in self.names}
        self.arrivals = {ns: 0 for ns in self.names}
        self.run_min = float('inf')
        self.run_min_info = None          # (t, pair, d)
        self.episodes = []                # closed: dict(t0,t1,minv,pair)
        self.cur_ep = None
        self.below_time = 0.0
        self.last_tick = None
        self.t0 = time.monotonic()
        self.csv = open(f'{args.outdir}/sep_{args.tag}.csv', 'w', buffering=1)
        self.csv.write('t,min_dist,pair\n')
        self.last_csv = 0.0

        for ns in self.names:
            self.create_subscription(
                State, f'/{ns}/state',
                lambda m, n=ns: self.state_cb(n, m), 10)
            self.create_subscription(
                PoseStamped, f'/{ns}/term_goal',
                lambda m, n=ns: self.goal_cb(n, m), 10)
        self.create_timer(0.05, self.tick)
        self.create_timer(30.0, self.status)

    def state_cb(self, ns, msg):
        self.pos[ns] = (msg.pos.x, msg.pos.y, time.monotonic())
        g = self.goal.get(ns)
        if g is not None and not self.reached_current[ns]:
            if math.hypot(msg.pos.x - g[0], msg.pos.y - g[1]) < 1.05:
                self.reached_current[ns] = True
                self.arrivals[ns] += 1

    def goal_cb(self, ns, msg):
        newg = (msg.pose.position.x, msg.pose.position.y)
        old = self.goal.get(ns)
        if old is None or math.hypot(newg[0] - old[0], newg[1] - old[1]) >= 1.0:
            if old is not None:
                self.legs[ns] += 1
            self.goal[ns] = newg
            self.reached_current[ns] = False

    def tick(self):
        now = time.monotonic()
        dt = 0.0 if self.last_tick is None else now - self.last_tick
        self.last_tick = now
        t = now - self.t0

        fresh = [(ns, x, y) for ns, (x, y, ts) in self.pos.items()
                 if now - ts < 1.0]
        if len(fresh) >= 2:
            mind, pair = float('inf'), ''
            for i in range(len(fresh)):
                for j in range(i + 1, len(fresh)):
                    d = math.hypot(fresh[i][1] - fresh[j][1],
                                   fresh[i][2] - fresh[j][2])
                    if d < mind:
                        mind, pair = d, f'{fresh[i][0]}-{fresh[j][0]}'
            if mind < self.run_min:
                self.run_min = mind
                self.run_min_info = (round(t, 1), pair, round(mind, 3))
            thr = self.args.threshold
            if mind < thr:
                self.below_time += dt
                if self.cur_ep is None:
                    self.cur_ep = {'t0': round(t, 1), 'minv': mind, 'pair': pair}
                elif mind < self.cur_ep['minv']:
                    self.cur_ep.update(minv=mind, pair=pair)
            elif self.cur_ep is not None and mind > thr + 0.05:
                self.close_ep(t)
            if t - self.last_csv >= 0.25:
                self.csv.write(f'{t:.2f},{mind:.3f},{pair}\n')
                self.last_csv = t

        if t >= self.args.duration:
            self.finish(t)

    def close_ep(self, t):
        ep = self.cur_ep
        ep['t1'] = round(t, 1)
        ep['minv'] = round(ep['minv'], 3)
        self.episodes.append(ep)
        print(f"VIOL t={ep['t0']}-{ep['t1']}s pair={ep['pair']} "
              f"min={ep['minv']}", flush=True)
        self.cur_ep = None

    def counts(self):
        eps = self.episodes + ([self.cur_ep] if self.cur_ep else [])
        return (len(eps),
                sum(1 for e in eps if e['minv'] < 0.5),
                sum(1 for e in eps if e['minv'] < 0.3))

    def status(self):
        t = time.monotonic() - self.t0
        n1, n05, n03 = self.counts()
        print(f'STATUS t={t:.0f}s fresh={len(self.pos)} '
              f'run_min={self.run_min:.3f} eps={n1}(<0.5:{n05} <0.3:{n03}) '
              f'below_s={self.below_time:.1f} '
              f'arrivals={sum(self.arrivals.values())} '
              f'legs={sum(self.legs.values())}', flush=True)

    def finish(self, t):
        if self.cur_ep is not None:
            self.close_ep(t)
        n1, n05, n03 = self.counts()
        summary = {
            'tag': self.args.tag,
            'duration_s': round(t, 1),
            'agents_seen': len(self.pos),
            'run_min': round(self.run_min, 3),
            'run_min_at': self.run_min_info,
            'episodes_below_thr': n1,
            'episodes_below_0p5': n05,
            'episodes_below_0p3': n03,
            'time_below_thr_s': round(self.below_time, 1),
            'arrivals_total': sum(self.arrivals.values()),
            'legs_total': sum(self.legs.values()),
            'arrivals_per_robot': self.arrivals,
            'worst_episodes': sorted(self.episodes,
                                     key=lambda e: e['minv'])[:8],
        }
        line = 'SUMMARY ' + json.dumps(summary)
        print(line, flush=True)
        with open(f'{self.args.outdir}/sep_{self.args.tag}.json', 'w') as f:
            json.dump(summary, f, indent=1)
        self.csv.close()
        raise SystemExit(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--agents', type=int, default=10)
    ap.add_argument('--prefix', default='NX')
    ap.add_argument('--threshold', type=float, default=1.0)
    ap.add_argument('--duration', type=float, default=360.0)
    ap.add_argument('--tag', default='run')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()
    rclpy.init()
    node = SepProbe(args)
    try:
        rclpy.spin(node)
    except SystemExit:
        pass


if __name__ == '__main__':
    main()
