#!/usr/bin/env python3
"""
Dynamic obstacle collision + closest-distance analysis with caching.

Example:
  python3 src/dynus/benchmarking/dynamic_obstacle_collision_checker.py 
        --num_obstacles 50 
        --bag_folder /media/kkondo/kota_elements/mighty/dynamic_obstacle/bags/mighty 
        --drone_radius 0.1 
        --sample_interval 0.01 
        --hist_bins 60

Re-plot later from cache only (no bag I/O):
python3 src/dynus/benchmarking/dynamic_obstacle_collision_checker.py 
        --num_obstacles 50 
        --bag_folder /media/kkondo/kota_elements/mighty/dynamic_obstacle/bags/mighty 
        --drone_radius 0.1 
        --sample_interval 0.01 
        --read-cache-only
        --hist_bins 60
"""

import os
import argparse
import math
from pathlib import Path
from typing import Dict, Tuple, Optional, List

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# rosbag2_py imports
from rosbag2_py import SequentialReader, StorageOptions, ConverterOptions
from rclpy.serialization import deserialize_message
from tf2_msgs.msg import TFMessage

FrameSeries = Dict[str, Dict[str, np.ndarray]]  # {"frame": {"times": np.ndarray[N], "pos": np.ndarray[N,3]}}

def _euclidean(a: np.ndarray, b: np.ndarray) -> float:
    d = a - b
    return float(np.sqrt(np.dot(d, d)))


def _read_tf_series(bag_path: str, frames_of_interest: List[str]) -> Tuple[FrameSeries, Optional[float], Optional[float]]:
    """Read /tf and return per-frame time series restricted to frames_of_interest."""
    frames_set = set(frames_of_interest)
    series: FrameSeries = {}
    for f in frames_of_interest:
        series[f] = {"times": [], "pos": []}

    storage_options = StorageOptions(uri=bag_path, storage_id='sqlite3')
    converter_options = ConverterOptions(input_serialization_format='cdr', output_serialization_format='cdr')
    reader = SequentialReader()
    reader.open(storage_options, converter_options)

    t0 = None
    t1 = None

    while reader.has_next():
        topic, data, timestamp = reader.read_next()
        if topic != "/tf":
            continue

        t_sec = timestamp / 1e9
        if t0 is None or t_sec < t0:
            t0 = t_sec
        if t1 is None or t_sec > t1:
            t1 = t_sec

        try:
            msg = deserialize_message(data, TFMessage)
        except Exception as e:
            print(f"Error deserializing message: {e}")
            continue

        for tr in msg.transforms:
            frame = tr.child_frame_id.strip()
            if frame not in frames_set:
                continue
            pos = (tr.transform.translation.x, tr.transform.translation.y, tr.transform.translation.z)
            s = series[frame]
            s["times"].append(t_sec)
            s["pos"].append(pos)

    # Convert to numpy & sort just in case
    for f in list(series.keys()):
        times = np.asarray(series[f]["times"], dtype=np.float64)
        pos = np.asarray(series[f]["pos"], dtype=np.float64)
        if times.size == 0:
            # keep placeholders; will be handled by caller
            series[f]["times"] = times
            series[f]["pos"] = pos
            continue
        order = np.argsort(times)
        series[f]["times"] = times[order]
        series[f]["pos"] = pos[order]

    return series, t0, t1


def _get_latest(series_times: np.ndarray, series_pos: np.ndarray, t_sample: float) -> Optional[np.ndarray]:
    """Latest position at or before t_sample."""
    if series_times.size == 0:
        return None
    idx = np.searchsorted(series_times, t_sample, side="right") - 1
    if idx >= 0:
        return series_pos[idx]
    return None


def analyze_bag_for_distances(
    bag_path: str,
    num_obstacles: int,
    sample_interval: float,
    drone_frame: str = "NX01/base_link",
) -> Tuple[np.ndarray, np.ndarray, bool]:
    """
    Returns:
      rel_times: shape [K], times (sec) relative to bag start **only for valid samples**
      min_dists: shape [K], closest obstacle distance per valid sample time
      collided:  (unused here; kept for signature compatibility)
    """
    # ----- gating config -----
    START_POS = np.array([0.0, 0.0, 3.0], dtype=np.float64)
    START_RADIUS = 1.0  # start collecting once > 1.0 m from start
    GOAL_POS  = np.array([105.0, 0.0, 3.0], dtype=np.float64)
    GOAL_RADIUS  = 1.0  # stop when within 1.0 m of goal

    frames = [drone_frame] + [f"obstacle_{i}" for i in range(num_obstacles)]
    series, t0, t1 = _read_tf_series(bag_path, frames)

    if t0 is None or t1 is None:
        print("No TF data found in the bag.")
        return np.array([]), np.array([]), False

    total_steps = int((t1 - t0) / sample_interval) + 1
    t_samples = t0 + np.arange(total_steps, dtype=np.float64) * sample_interval

    # Only append when we actually have both the drone and at least one obstacle at that time.
    rel_times_valid: List[float] = []
    min_dists_valid: List[float] = []

    drone_times = series[drone_frame]["times"]
    drone_pos   = series[drone_frame]["pos"]

    # Pre-bind obstacle arrays for speed
    obs_times = [series[f"obstacle_{i}"]["times"] for i in range(num_obstacles)]
    obs_pos   = [series[f"obstacle_{i}"]["pos"]   for i in range(num_obstacles)]

    collecting = False  # flips to True once drone > START_RADIUS from START_POS

    for t in tqdm(t_samples, desc=f"Sampling {Path(bag_path).name}"):
        dp = _get_latest(drone_times, drone_pos, t)
        if dp is None:
            continue

        # ----- gating: wait until drone has moved away from the start -----
        if not collecting:
            if _euclidean(dp, START_POS) <= START_RADIUS:
                continue
            collecting = True

        # ----- gating: stop once near the goal -----
        if _euclidean(dp, GOAL_POS) <= GOAL_RADIUS:
            break

        # ----- compute closest obstacle distance only while collecting -----
        best = np.inf
        for i in range(num_obstacles):
            if obs_times[i].size == 0:
                continue
            op = _get_latest(obs_times[i], obs_pos[i], t)
            if op is None:
                continue
            d = _euclidean(dp, op)
            if d < best:
                best = d

        if np.isfinite(best):
            rel_times_valid.append(float(t - t0))
            min_dists_valid.append(float(best))

    return np.asarray(rel_times_valid), np.asarray(min_dists_valid), False

def save_cache(cache_dir: Path, bag_key: str, rel_times: np.ndarray, min_dists: np.ndarray,
               meta: Dict[str, float]) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{bag_key}_closest_dists.npz"
    np.savez_compressed(out, times=rel_times, min_dists=min_dists, **meta)
    return out


def load_cache(cache_dir: Path, bag_key: str) -> Optional[Tuple[np.ndarray, np.ndarray, Dict[str, float]]]:
    path = cache_dir / f"{bag_key}_closest_dists.npz"
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=False)
    rel_times = data["times"]
    min_dists = data["min_dists"]
    meta = {k: float(data[k]) for k in data.files if k not in ("times", "min_dists")}
    return rel_times, min_dists, meta


def plot_hist(min_dists_all: np.ndarray, bins: int, out_path: Path, title: str):
    valid = min_dists_all[np.isfinite(min_dists_all)]
    if valid.size == 0:
        print("No valid distances to plot.")
        return
    weights = np.ones(valid.shape[0], dtype=np.float64) * (100.0 / valid.shape[0])
    plt.figure(figsize=(7, 4.5))
    plt.hist(valid, bins=bins, weights=weights)  # percentages
    plt.xlabel("Closest obstacle distance [m]")
    plt.ylabel("Percentage of samples [%]")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Saved histogram to: {out_path}")



def main():
    parser = argparse.ArgumentParser(description="Collision + closest-distance analyzer with caching & histograms.")
    parser.add_argument('--num_obstacles', type=int, required=True,
                        help="Number of obstacles (expected indices: 0..num_obstacles-1)")
    parser.add_argument('--bag_folder', type=str, required=True,
                        help="Folder containing rosbag2 directories (e.g., num_0, num_1, ...)")
    parser.add_argument('--drone_radius', type=float, default=0.1, help="Collision radius [m].")
    parser.add_argument('--sample_interval', type=float, default=0.01, help="Sampling interval [s].")
    parser.add_argument('--hist_bins', type=int, default=50, help="Number of histogram bins.")
    parser.add_argument('--per_bag_hist', action='store_true', help="Also write one histogram per bag.")
    parser.add_argument('--read-cache-only', action='store_true',
                        help="Do not read bags; only read cached NPZ and plot/aggregate.")
    parser.add_argument('--reuse-cache', dest='reuse_cache', action='store_true', default=True,
                        help="Reuse cache when present (default: True).")
    parser.add_argument('--no-reuse-cache', dest='reuse_cache', action='store_false',
                        help="Force recomputation even if cache exists.")
    args = parser.parse_args()

    bag_folder = Path(args.bag_folder)
    cache_dir = bag_folder / "closest_dist_cache"
    out_dir = bag_folder / "plots_closest_dist"
    out_dir.mkdir(parents=True, exist_ok=True)

    # List bag directories like num_0, num_1, ...
    bag_dirs = sorted([d for d in os.listdir(bag_folder) if d.startswith("num_")
                       and (bag_folder / d).is_dir()])

    all_min_dists = []
    results = {}

    print(f"List of bags to process: {bag_dirs}")

    for bag_dir in bag_dirs:
        bag_key = bag_dir  # e.g., "num_3"
        # Bags are typically at <bag_folder>/num_X/num_X
        bag_path = (bag_folder / bag_dir / bag_dir).as_posix()

        rel_times = None
        min_dists = None

        if args.read_cache_only:
            cached = load_cache(cache_dir, bag_key)
            if cached is None:
                print(f"[WARN] No cache for {bag_key}, skipping (read-cache-only).")
                results[bag_dir] = "no-cache"
                continue
            rel_times, min_dists, meta = cached
        else:
            if args.reuse_cache:
                cached = load_cache(cache_dir, bag_key)
                if cached is not None:
                    rel_times, min_dists, meta = cached
                    # If sampling mismatch, recompute
                    if "sample_interval" in meta and abs(meta["sample_interval"] - args.sample_interval) > 1e-12:
                        print(f"[INFO] Cache sample_interval differs for {bag_key}; recomputing.")
                        rel_times = None
                        min_dists = None

            if rel_times is None or min_dists is None:
                print(f"Processing bag: {bag_path}")
                rel_times, min_dists, _ = analyze_bag_for_distances(
                    bag_path=bag_path,
                    num_obstacles=args.num_obstacles,
                    sample_interval=args.sample_interval,
                    drone_frame="NX01/base_link",
                )
                if rel_times.size == 0:
                    results[bag_dir] = "no tf"
                    continue
                meta = {
                    "sample_interval": float(args.sample_interval),
                    "num_obstacles": float(args.num_obstacles),
                }
                save_path = save_cache(cache_dir, bag_key, rel_times, min_dists, meta)
                print(f"Cached distances -> {save_path}")

        # Collision based on Euclidean distance, with empty-guard
        if min_dists.size > 0:
            dmin = float(min_dists.min())
            collided = (dmin < args.drone_radius)
            results[bag_dir] = "collision" if collided else "no collision"

            if collided:
                print(f"\033[91m{bag_key}: COLLISION (min distance: {dmin:.3f} m)\033[0m")
            else:
                print(f"\033[92m{bag_key}: no collision (min distance: {dmin:.3f} m)\033[0m")
        else:
            # No valid (drone, obstacle) pairs ever found at the same sampled time
            collided = False
            results[bag_dir] = "no data"
            print(f"\033[93m{bag_key}: no data (no valid pairs at sampled times; min distance: N/A)\033[0m")

        # Per-bag histogram (optional)
        if args.per_bag_hist:
            plot_hist(min_dists, args.hist_bins, out_dir / f"{bag_key}_closest_dist_hist.png",
                      title=f"Closest distances: {bag_key}")

        # Accumulate for global histogram
        all_min_dists.append(min_dists)

    # Write results file with stats
    output_file = bag_folder / "dynamic_collision_check.txt"
    with open(output_file, "w") as f:
        for bag_dir in bag_dirs:
            res = results.get(bag_dir, "skipped")
            line = f"{bag_dir}: {res}"
            # Add quick stats if cached/processed
            cached = load_cache(cache_dir, bag_dir)
            if cached is not None:
                _, md, _ = cached
                valid = md[np.isfinite(md)]
                if valid.size > 0:
                    stats = dict(
                        min=float(np.min(valid)),
                        p05=float(np.percentile(valid, 5)),
                        p50=float(np.percentile(valid, 50)),
                        p95=float(np.percentile(valid, 95)),
                    )
                    line += f" | min={stats['min']:.3f} m, p05={stats['p05']:.3f}, p50={stats['p50']:.3f}, p95={stats['p95']:.3f}"
            f.write(line + "\n")

        # Global summary
        if all_min_dists:
            stacked = np.concatenate([a[np.isfinite(a)] for a in all_min_dists if a.size > 0], axis=0)
            if stacked.size > 0:
                gmin = float(np.min(stacked))
                p05 = float(np.percentile(stacked, 5))
                p50 = float(np.percentile(stacked, 50))
                p95 = float(np.percentile(stacked, 95))
                f.write("\nSUMMARY (all bags)\n")
                f.write(f"min={gmin:.3f} m, p05={p05:.3f}, p50={p50:.3f}, p95={p95:.3f}\n")

    print(f"Results written to {output_file}")

    # Global histogram
    if all_min_dists:
        stacked = np.concatenate([a[np.isfinite(a)] for a in all_min_dists if a.size > 0], axis=0)
        if stacked.size > 0:
            plot_hist(stacked, args.hist_bins, out_dir / "closest_dist_hist_all.png",
                      title="Closest obstacle distance (all bags)")


if __name__ == '__main__':
    main()
