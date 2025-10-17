#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick visualization of FOV cache entries.
Usage:
  # 指定 s0=930m，显示线段模式，打标并保存
python plot_fov.py --csv Monza.csv --cache fov_cache.json \
  --s0 930 --mode segment --markers --fill \
  --show-uz-normals --show-tangent-rays --panx -40 --pany 15
"""
from __future__ import annotations
import json, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from fov_extended import build_boundaries, compute_visibility_border_discrete

def compute_arclength(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ds = np.hypot(np.diff(x), np.diff(y))
    s = np.zeros_like(x, dtype=float)
    s[1:] = np.cumsum(ds)
    return s

def load_track(csv_path: str):
    df = pd.read_csv(csv_path)
    xcol = '# x_m' if '# x_m' in df.columns else 'x_m'
    ycol = 'y_m'
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)
    s = compute_arclength(x,y)
    wR = df['w_tr_right_m'].to_numpy(float) if 'w_tr_right_m' in df.columns else np.full_like(s, 5.0)
    wL = df['w_tr_left_m' ].to_numpy(float) if 'w_tr_left_m'  in df.columns else np.full_like(s, 5.0)
    xsp = CubicSpline(s,x); ysp = CubicSpline(s,y)
    return s, xsp, ysp, wL, wR

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--cache', default='fov_cache.json')
    ap.add_argument('--k', type=int, default=None, help='index of frame to visualize')
    ap.add_argument('--s0', type=float, default=None, help='arclength to visualize (overrides --k)')
    ap.add_argument('--mode', choices=['segment','front'], default='segment')
    ap.add_argument('--save', default=None, help='path to save the figure (optional)')
    ap.add_argument('--no-legend', action='store_true', help='hide legend')
    ap.add_argument('--markers', action='store_true', help='show debug markers: ITP/FPV/Unknown-front')
    ap.add_argument('--fill', action='store_true', help='fill FOV area (green) using front border (computed on the fly if needed)')
    ap.add_argument('--alpha', type=float, default=0.22, help='fill transparency for FOV/Unknown')
    ap.add_argument('--show-uz-normals', action='store_true', help='draw normal half-planes of Unknown zones')
    ap.add_argument('--show-tangent-rays', action='store_true', help='draw tangent rays for ITP scanning (car→ITP, ITP→frontiers)')
    ap.add_argument('--back', type=float, default=40.0, help='meters to include *behind* s0 for drawing (visualization only)')
    ap.add_argument('--panx', type=float, default=0.0, help='pan the view window horizontally (meters, + moves right)')
    ap.add_argument('--pany', type=float, default=0.0, help='pan the view window vertically (meters, + moves up)')
    args = ap.parse_args()

    with open(args.cache,'r',encoding='utf-8') as f:
        cache = json.load(f)

    s, xsp, ysp, wL, wR = load_track(args.csv)
    entries = cache['entries']
    # choose frame by --s0 if provided, else by --k (default to 0)
    if args.s0 is not None:
        # pick the entry with nearest s0
        s0_array = np.array([en['s0'] for en in entries], dtype=float)
        k = int(np.argmin(np.abs(s0_array - args.s0)))
    else:
        k = 0 if args.k is None else int(max(0, min(args.k, len(entries)-1)))
    e = entries[k]
    s0 = float(e['s0'])

    # local window (extend backward for drawing if needed)
    look = cache['params'].get('look', 200.0)
    s_start = max(0.0, s0 - float(args.back))
    s_dense = np.arange(s_start, s0 + look, 0.5)
    try:
        wL_loc = np.interp(s_dense, xsp.x, wL)
        wR_loc = np.interp(s_dense, xsp.x, wR)
    except Exception:
        wL_loc = np.full_like(s_dense, wL[0])
        wR_loc = np.full_like(s_dense, wR[0])
    (xL,yL),(xR,yR) = build_boundaries(xsp,ysp,wL_loc,wR_loc,s_dense)

    def s_to_xy_on_side(s_val, side):
        j = int(np.clip(np.searchsorted(s_dense, s_val, 'left'), 0, len(s_dense)-1))
        return (xR[j], yR[j]) if side == 0 else (xL[j], yL[j])

    plt.figure(figsize=(10, 8))
    plt.plot(xL,yL, color='k', linewidth=2)
    plt.plot(xR,yR, color='k', linewidth=2)

    car = np.array([xsp(s0), ysp(s0)])
    plt.scatter([car[0]],[car[1]], marker='s', zorder=5)

    # heading arrow from spline derivative
    dt = 1e-3
    dx = float(xsp.derivative(1)(s0))
    dy = float(ysp.derivative(1)(s0))
    norm = (dx**2 + dy**2)**0.5 or 1.0
    ux, uy = dx/norm, dy/norm
    plt.arrow(car[0], car[1], 6*ux, 6*uy, length_includes_head=True, head_width=0.6, alpha=0.8)

    # Unknown zones (use boundary segment from s_start to s_end on the same side, closed to car)
    for uz in e.get('unknown_zones', []):
        if uz.get('s_start') is None or uz.get('s_end') is None:
            continue
        side = int(uz['side'])  # 0=right,1=left
        s_a, s_b = float(uz['s_start']), float(uz['s_end'])
        # ensure s_a <= s_b
        if s_a > s_b:
            s_a, s_b = s_b, s_a
        j0 = int(np.clip(np.searchsorted(s_dense, s_a, 'left'), 0, len(s_dense)-1))
        j1 = int(np.clip(np.searchsorted(s_dense, s_b, 'left'), 0, len(s_dense)-1))
        bx = (xR if side==0 else xL)[j0:j1+1]
        by = (yR if side==0 else yL)[j0:j1+1]
        if len(bx) >= 2:
            poly_x = [car[0]] + list(bx) + [car[0]]
            poly_y = [car[1]] + list(by) + [car[1]]
            plt.fill(poly_x, poly_y, color=(0.8,0,0,0.6), linewidth=0)  # dark red, fixed alpha
            plt.plot(bx, by, color='red', alpha=0.9, linewidth=1.5)

    # if args.show_uz_normals:
    #     for uz in e.get('unknown_zones', []):
    #         if uz.get('s_start') is not None and uz.get('s_end') is not None:
    #             # get endpoints XY
    #             p0 = s_to_xy_on_side(uz['s_start'], uz['side'])
    #             p1 = s_to_xy_on_side(uz['s_end'], uz['side'])
    #             # midpoint and edge vector
    #             mx, my = (0.5*(p0[0]+p1[0]), 0.5*(p0[1]+p1[1]))
    #             ex, ey = (p1[0]-p0[0], p1[1]-p0[1])
    #             nx, ny = -ey, ex
    #             norm_len = (nx**2+ny**2)**0.5 or 1.0
    #             nx/=norm_len; ny/=norm_len
    #             plt.arrow(mx,my, 8*nx, 8*ny, head_width=1.0, color='red', alpha=0.5)

    # FOV border
    fb = e.get('fov_border', {})
    if fb and fb.get('mode','segment') == 'segment' and len(fb.get('points',[]))==2:
        (ix,iy),(fx,fy) = fb['points']
        plt.plot([ix,fx],[iy,fy], linewidth=2, label='ITP→FPV')
    elif fb and fb.get('points'):
        px,py = zip(*fb['points'])
        plt.plot(px,py, linewidth=2, label='FOV front')

    # Optional: fill FOV visible area in green as a corridor (boundary segment + visible border)
    if args.fill:
        fpv = e.get('final_fpv')
        border_pts = []
        if fb and fb.get('points'):
            border_pts = fb['points']  # either segment [ITP,FPV] or front polyline
        else:
            # fallback: compute a discrete front
            tmp = compute_visibility_border_discrete(
                spline_x=xsp, spline_y=ysp, w_left=wL_loc, w_right=wR_loc,
                s0=s0, look=look, unknown_zones=e.get('unknown_zones', []),
                fpv_side=(fpv or {}).get('side'), s_fpv=(fpv or {}).get('s_fpv'),
                s_step=0.5, n_angles=721,
            )
            border_pts = tmp.get('points', [])
        if fpv is not None and border_pts:
            # near-side boundary segment from s0 to s_itp
            itp_side = int(fb['side'][0]) if (fb and fb.get('side')) else int(fpv['side'])
            j0 = int(np.clip(np.searchsorted(s_dense, s0, 'left'), 0, len(s_dense)-1))
            jI = int(np.clip(np.searchsorted(s_dense, fpv['s_itp'], 'left'), 0, len(s_dense)-1))

            # 1. ITP侧边界段 (s0→s_itp)
            seg_x = (xR if itp_side==0 else xL)[j0:jI+1]
            seg_y = (yR if itp_side==0 else yL)[j0:jI+1]

            # 2. 可见前沿
            front_x = [p[0] for p in border_pts]
            front_y = [p[1] for p in border_pts]

            # 3. FPV侧边界段 (s_fpv→s0, 另一侧边界)
            jF = int(np.clip(np.searchsorted(s_dense, fpv['s_fpv'], 'left'), 0, len(s_dense)-1))
            seg2_x = (xR if fpv['side']==0 else xL)[jF:j0-1:-1]
            seg2_y = (yR if fpv['side']==0 else yL)[jF:j0-1:-1]

            # Add opposite side boundary segment from car to FPV
            j0_other = int(np.clip(np.searchsorted(s_dense, s0, 'left'), 0, len(s_dense)-1))
            jF_other = int(np.clip(np.searchsorted(s_dense, fpv['s_fpv'], 'left'), 0, len(s_dense)-1))
            seg_other_x = (xL if itp_side==0 else xR)[j0_other:jF_other+1]
            seg_other_y = (yL if itp_side==0 else yR)[j0_other:jF_other+1]

            poly_x = [car[0]] + list(seg_x) + list(front_x) + list(seg_other_x[::-1]) + [car[0]]
            poly_y = [car[1]] + list(seg_y) + list(front_y) + list(seg_other_y[::-1]) + [car[1]]
            if len(poly_x) >= 3:
                plt.fill(poly_x, poly_y, color='green', alpha=args.alpha, linewidth=0, label='FOV area')

    # === Debug markers: ITP / FPV / Unknown-front (optional) ===
    if args.markers:
        # final FPV / ITP markers
        fpv = e.get('final_fpv')
        if fpv is not None:
            # try reading sides from fov_border first (more direct)
            itp_side = None
            if fb and fb.get('side') and len(fb['side']) >= 1:
                itp_side = int(fb['side'][0])
            if itp_side is None:
                itp_side = int(fpv['side'])  # fallback if not stored
            itp_xy = s_to_xy_on_side(fpv['s_itp'], itp_side)
            fpv_xy = s_to_xy_on_side(fpv['s_fpv'], fpv['side'])
            plt.scatter([itp_xy[0]],[itp_xy[1]], c='r', s=60, zorder=6, label='ITP')
            plt.scatter([fpv_xy[0]],[fpv_xy[1]], c='g', s=60, zorder=6, label='FPV')

            if args.show_tangent_rays:
                plt.plot([car[0], itp_xy[0]],[car[1], itp_xy[1]], 'b--', alpha=0.6)
                plt.plot([itp_xy[0], fpv_xy[0]],[itp_xy[1], fpv_xy[1]], 'g--', alpha=0.6)
                # and for unknowns: from ITP to front point
                for uz in e.get('unknown_zones', []):
                    uz_front_xy = s_to_xy_on_side(uz['s_end'], uz['side'])
                    plt.plot([itp_xy[0], uz_front_xy[0]],[itp_xy[1], uz_front_xy[1]], 'r--', alpha=0.6)

        # Unknown-front points
        for uz in e.get('unknown_zones', []):
            uz_front_xy = s_to_xy_on_side(uz['s_end'], uz['side'])
            plt.scatter([uz_front_xy[0]],[uz_front_xy[1]], c='orange', s=50, zorder=6, label='Unknown front')

    plt.axis('equal')
    plt.title(f"s0={s0:.1f}  (frame {k}/{len(entries)-1})  mode={fb.get('mode','segment')}")
    plt.tight_layout()
    if not args.no_legend:
        plt.legend(loc='upper left', frameon=False, ncol=2)
    if args.save:
        plt.savefig(args.save, dpi=200)
        print(f"[Saved] {args.save}")

    # zoom around the car position but with a longer forward view
    zoom_x_back = 120   # keep 40m behind
    zoom_x_front = 120 # extend forward view
    zoom_y = 80        # wider lateral view
    plt.xlim(car[0] - zoom_x_back + args.panx, car[0] + zoom_x_front + args.panx)
    plt.ylim(car[1] - zoom_y + args.pany,  car[1] + zoom_y + args.pany)

    plt.show()

if __name__ == '__main__':
    main()