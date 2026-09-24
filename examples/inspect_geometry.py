"""Run one viewpoint through the existing sampled FOV implementation.

This is an execution example, not an independent correctness certificate.
The track-column convention and natural splines match precompute_fov_cache.py.
"""
from pathlib import Path
import argparse
import json
import sys

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fov_extended import compute_fov_at_s0


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Cannot encode {type(value).__name__}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--track", type=Path, default=ROOT / "Monza.csv")
    parser.add_argument("--s0", type=float, default=930.0)
    parser.add_argument("--look", type=float, default=200.0)
    args = parser.parse_args()
    data = np.loadtxt(args.track, delimiter=",")
    if data.ndim != 2 or data.shape[1] < 4:
        parser.error("Track columns must be x, y, right width, left width.")
    station = np.r_[0., np.cumsum(np.linalg.norm(np.diff(data[:, :2], axis=0), axis=1))]
    if np.any(np.diff(station) <= 0):
        parser.error("Consecutive centerline samples must have distinct positions.")
    if not 0 <= args.s0 < station[-1] or args.look <= 0:
        parser.error("Require 0 <= s0 < track length and look > 0.")
    x = CubicSpline(station, data[:, 0], bc_type="natural")
    y = CubicSpline(station, data[:, 1], bc_type="natural")
    result = compute_fov_at_s0(x, y, data[:, 3], data[:, 2], args.s0, look=args.look)
    summary = {
        "track": args.track.name,
        "s0_m": args.s0,
        "look_ahead_m": args.look,
        "sampled_front": result.get("final_fpv"),
        "interpretation": "Output of the current sampled implementation; see docs/implementation_notes.md.",
    }
    print(json.dumps(summary, indent=2, default=json_value))


if __name__ == "__main__":
    main()
