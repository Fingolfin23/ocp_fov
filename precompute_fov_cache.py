#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Precompute FOV/Unknown-Zone cache for Monza (or any track CSV in the same format),
so that OCP and visualization can load a single file quickly.

Inputs
------
--csv: path to Monza.csv. Expected columns (case tolerant):
    '# x_m' or 'x_m', 'y_m', 'w_tr_right_m', 'w_tr_left_m'
If width columns are missing, use --w to set a constant half-width.

Outputs
-------
- fov_cache.json : full per-s0 dictionaries for later visualization/OCP
- fov_cache.npz  : compact arrays for fast OCP lookup

Dependencies
------------
NumPy, Pandas, SciPy
This script uses functions from fov_extended.py in the same folder.

cd Retake
python precompute_fov_cache.py --csv Monza.csv \
  --look 200 --s0_step 5 --s_dense_step 0.5 \
  --theta_min 4 --peak_min 0 --min_dist 5 --min_arc 10 \
  --eps_dir 1e-3 --fov_border_mode segment

"""
from __future__ import annotations
import os, json, argparse
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from fov_extended import compute_fov_at_s0

# ---- IO utils ----
def compute_arclength(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ds = np.hypot(np.diff(x), np.diff(y))
    s = np.zeros_like(x, dtype=float)
    s[1:] = np.cumsum(ds)
    return s

def load_track(csv_path: str, default_half_width: float = 5.0):
    """Load Monza-like CSV and return (s, x_spline, y_spline, w_left, w_right)."""
    df = pd.read_csv(csv_path)
    def pick(*cands):
        for c in cands:
            if c in df.columns: return c
        return None
    xcol = pick('# x_m', 'x_m', 'x')
    ycol = pick('y_m', 'y')
    if xcol is None or ycol is None:
        raise ValueError(f"CSV missing x/y columns. Found: {list(df.columns)}")
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    s = compute_arclength(x, y)

    wr_col = pick('w_tr_right_m', 'wr_m', 'w_r')
    wl_col = pick('w_tr_left_m',  'wl_m', 'w_l')
    if wr_col is not None and wl_col is not None:
        w_right = df[wr_col].to_numpy(float)
        w_left  = df[wl_col].to_numpy(float)
    else:
        w_right = np.full_like(s, default_half_width, dtype=float)
        w_left  = np.full_like(s, default_half_width, dtype=float)

    x_spline = CubicSpline(s, x, bc_type='natural')
    y_spline = CubicSpline(s, y, bc_type='natural')
    return s, x_spline, y_spline, w_left, w_right

# ---- Main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, help='Path to Monza.csv')
    ap.add_argument('--out_json', default='fov_cache.json')
    ap.add_argument('--out_npz',  default='fov_cache.npz')
    ap.add_argument('--look', type=float, default=200.0, help='forward look distance [m] per s0')
    ap.add_argument('--s0_step', type=float, default=5.0, help='sampling step along arclength [m]')
    ap.add_argument('--w', type=float, default=5.0, help='half-width if width columns missing')

    # FOV params（离散对齐阈值）
    ap.add_argument('--s_dense_step', type=float, default=0.5)
    ap.add_argument('--theta_min', type=float, default=4.0)
    ap.add_argument('--peak_min', type=float, default=0.0)
    ap.add_argument('--min_dist', type=float, default=5.0)
    ap.add_argument('--min_arc', type=float, default=10.0)
    ap.add_argument('--kappa_min', type=float, default=1e-5)
    ap.add_argument('--eps_dir', type=float, default=1e-3)
    ap.add_argument('--fov_border_mode', choices=['segment','front'], default='segment')

    args = ap.parse_args()

    s, xsp, ysp, w_left, w_right = load_track(args.csv, default_half_width=args.w)
    total_len = float(s[-1])
    s0_list = np.arange(0.0, total_len, args.s0_step)

    entries = []             # 完整字典（给可视化/OCP都能直接用）
    s0_arr, s_fpv_arr, side_arr = [], [], []  # 紧凑数组（给OCP快速查）

    for s0 in s0_list:
        res = compute_fov_at_s0(
            xsp, ysp, w_left, w_right, s0,
            look=args.look,
            kappa_spline=None,
            s_dense_step=args.s_dense_step,
            theta_min_deg=args.theta_min,
            peak_min_deg=args.peak_min,
            min_dist=args.min_dist,
            min_arc=args.min_arc,
            kappa_min=args.kappa_min,
            epsilon_dir=args.eps_dir,
            fov_border_mode=args.fov_border_mode,
        )
        e = {
            's0': float(res['s0']),
            'tangent_points': [
                {'s': float(tp['s']), 'side': int(tp['side']),
                 'point': (float(tp['point'][0]), float(tp['point'][1]))}
                for tp in res.get('tangent_points', [])
            ],
            'unknown_zones': [
                {'s_start': float(uz['s_start']), 's_end': float(uz['s_end']),
                 'side': int(uz['side']),
                 'polygon': [(float(x), float(y)) for (x,y) in uz.get('polygon', [])]}
                for uz in res.get('unknown_zones', [])
            ],
            'final_fpv': None if res.get('final_fpv') is None else {
                's_itp': float(res['final_fpv']['s_itp']),
                's_fpv': float(res['final_fpv']['s_fpv']),
                'side':  int(res['final_fpv']['side'])
            },
            'fov_border': {
                'mode': res['fov_border'].get('mode','segment') if res.get('fov_border') else 'segment',
                'points': [ (float(x), float(y)) for (x,y) in res.get('fov_border',{}).get('points',[]) ],
                's_on_boundary': [ float(si) for si in res.get('fov_border',{}).get('s_on_boundary',[]) ],
                'side': [ int(si) for si in res.get('fov_border',{}).get('side',[]) ],
                'angles': [ float(a) for a in res.get('fov_border',{}).get('angles',[]) ],
            }
        }
        entries.append(e)

        s0_arr.append(float(res['s0']))
        fpv = res.get('final_fpv')
        if fpv is None:
            s_fpv_arr.append(np.nan)
            side_arr.append(-1)
        else:
            s_fpv_arr.append(float(fpv['s_fpv']))
            side_arr.append(int(fpv['side']))

    # JSON（全量）
    payload = {
        'source_csv': os.path.abspath(args.csv),
        'total_length': total_len,
        'params': {
            'look': args.look, 's0_step': args.s0_step,
            's_dense_step': args.s_dense_step,
            'theta_min': args.theta_min, 'peak_min': args.peak_min,
            'min_dist': args.min_dist, 'min_arc': args.min_arc,
            'kappa_min': args.kappa_min, 'eps_dir': args.eps_dir,
            'fov_border_mode': args.fov_border_mode,
        },
        'entries': entries,
    }
    with open(args.out_json, 'w', encoding='utf-8') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"[Saved] {args.out_json} with {len(entries)} frames")

    # NPZ（紧凑）
    np.savez_compressed(
        args.out_npz,
        s0=np.array(s0_arr, dtype=float),
        s_fpv=np.array(s_fpv_arr, dtype=float),
        side=np.array(side_arr, dtype=int),
        total_length=total_len,
    )
    print(f"[Saved] {args.out_npz}")

if __name__ == '__main__':
    main()