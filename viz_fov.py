# viz_fov.py
# Visualize full-track FOV from a fov_cache.json over a track CSV.
# Usage examples:
#   Static (no GUI):
#     python viz_fov.py fov_cache.json --stride 5 --fill --save fov_monza.png --no_show
#   GIF (no GUI):
#     python viz_fov.py fov_cache.json --stride 3 --fill --gif fov_monza.gif --fps 20 --no_show
#   MP4 (requires ffmpeg):
#     python viz_fov.py fov_cache.json --stride 5 --fill --mp4 fov_monza.mp4 --fps 30 --no_show
#   Segment only:
#     python viz_fov.py fov_cache.json --range 0:500 --stride 2 --gif fov_seg.gif --fps 24 --no_show
'''
python viz_fov.py fov_cache.json \
  --track_fill --bounds_outline \
  --stride 2 --fps 20 --no_show \
  --follow_cam --cam_w 220 --cam_h 220 --cam_margin 12 \
  --sector --sector_r 120 --sector_half_deg 35 \
  --gif fov_follow.gif --gif_dpi 140 --out_px_w 1080 --out_px_h 1080
  '''

import argparse
import json
import sys
import numpy as np
import matplotlib
# Switch to non-interactive backend when exporting or asked to avoid GUI
if any(flag in sys.argv for flag in ['--gif', '--mp4', '--save', '--no_show', '--no-show']):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Tuple, List, Optional, Iterable, Deque
from matplotlib import animation
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D

# ----------------------------- Track utilities -----------------------------

def load_track(csv_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load track CSV: columns [x, y, right_width, left_width] (no header)."""
    arr = np.loadtxt(csv_path, delimiter=',')
    x, y, wr, wl = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    return x, y, wr, wl


def cumulative_s(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    dx = np.diff(x)
    dy = np.diff(y)
    s = np.zeros_like(x)
    s[1:] = np.cumsum(np.hypot(dx, dy))
    return s


def interp_xy(s_query: float, s_ref: np.ndarray, x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Linear interpolation of (x,y) at s_query along the polyline."""
    xq = np.interp(s_query, s_ref, x)
    yq = np.interp(s_query, s_ref, y)
    return float(xq), float(yq)


def heading_theta(x: np.ndarray, y: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Compute heading angle along the centerline."""
    dx_ds = np.gradient(x, s, edge_order=2)
    dy_ds = np.gradient(y, s, edge_order=2)
    return np.arctan2(dy_ds, dx_ds)


def build_boundaries(x: np.ndarray, y: np.ndarray, wl: np.ndarray, wr: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute left/right road boundaries from centerline and widths."""
    s = cumulative_s(x, y)
    th = heading_theta(x, y, s)
    # left = normal +90deg, right = -90deg
    xL = x - wl * np.sin(th)
    yL = y + wl * np.cos(th)
    xR = x + wr * np.sin(th)
    yR = y - wr * np.cos(th)
    return xL, yL, xR, yR

# -------------------------- Unknown-zone parsing ---------------------------

def _coerce_point(pt) -> Optional[Tuple[float, float]]:
    """Try to coerce a single point to (x, y) float tuple from various shapes."""
    if pt is None:
        return None
    # [x, y]
    if isinstance(pt, (list, tuple)) and len(pt) == 2 and all(isinstance(v, (int, float)) for v in pt):
        return float(pt[0]), float(pt[1])
    # {'x': ..., 'y': ...}
    if isinstance(pt, dict):
        if 'x' in pt and 'y' in pt:
            return float(pt['x']), float(pt['y'])
        if 'point' in pt and isinstance(pt['point'], (list, tuple)) and len(pt['point']) == 2:
            return _coerce_point(pt['point'])
    return None


def parse_unknown_zone_polys(zone_entry) -> List[np.ndarray]:
    """
    Parse an unknown-zone entry into a list of polygon arrays (N,2).
    Accept many shapes:
      - {'points': [[x,y], ...]}
      - {'vertices': [[x,y], ...]}
      - [{'x':..,'y':..}, ...]
      - [[x,y], ...]
      - {'type':'Polygon','coordinates': [[[x,y],...]] } (GeoJSON-like)
      - {'polygons': [ ...same as above... ]}  # multiple polys
    Returns a possibly empty list (no exception on malformed entries).
    """
    polys: List[np.ndarray] = []
    if zone_entry is None:
        return polys

    # Wrapper for multiple polygons
    if isinstance(zone_entry, dict) and 'polygons' in zone_entry:
        items = zone_entry.get('polygons', [])
        for it in items:
            polys.extend(parse_unknown_zone_polys(it))
        return polys

    # GeoJSON-like
    if isinstance(zone_entry, dict) and zone_entry.get('type', '').lower() == 'polygon' and 'coordinates' in zone_entry:
        coords = zone_entry['coordinates']
        if isinstance(coords, list) and len(coords) > 0:
            ring = coords[0]  # exterior ring
            pts = []
            for p in ring:
                cp = _coerce_point(p)
                if cp is not None:
                    pts.append(cp)
            if len(pts) >= 3:
                polys.append(np.asarray(pts, dtype=float))
        return polys

    # Common dict with "points" or "vertices"
    if isinstance(zone_entry, dict):
        pts_raw = zone_entry.get('points', zone_entry.get('vertices'))
        if isinstance(pts_raw, list) and len(pts_raw) >= 3:
            pts = []
            for p in pts_raw:
                cp = _coerce_point(p)
                if cp is not None:
                    pts.append(cp)
            if len(pts) >= 3:
                polys.append(np.asarray(pts, dtype=float))
            return polys
        # Single polygon nested under 'polygon' key
        if 'polygon' in zone_entry:
            polys.extend(parse_unknown_zone_polys(zone_entry['polygon']))
            return polys

    # List: either list of points, or list of polygons
    if isinstance(zone_entry, list):
        if len(zone_entry) == 0:
            return polys
        # Is it a list of polygons? (nested lists)
        is_poly_list = all(isinstance(el, (list, tuple, dict)) for el in zone_entry) and \
                       any(isinstance(el, list) and len(el) > 0 and isinstance(el[0], (list, tuple, dict)) for el in zone_entry)
        if is_poly_list:
            for sub in zone_entry:
                polys.extend(parse_unknown_zone_polys(sub))
            return polys
        # Otherwise: treat as a single polygon
        pts = []
        for p in zone_entry:
            cp = _coerce_point(p)
            if cp is not None:
                pts.append(cp)
        if len(pts) >= 3:
            polys.append(np.asarray(pts, dtype=float))
        return polys

    # Fallback
    return polys

# ----------------------------- Drawing helpers -----------------------------

def draw_fov(ax, p0: Tuple[float, float], tp: List[Tuple[float, float]], fill: bool,
             ray_lw: float = 1.2, fill_alpha: float = 0.12, ray_color: str = 'C0'):
    """Draw FOV from car point p0 to two tangent points tp (len==2)."""
    if len(tp) != 2:
        return
    x0, y0 = p0
    ax.plot([x0, tp[0][0]], [y0, tp[0][1]],
            linewidth=ray_lw, color=ray_color, solid_capstyle='round', zorder=4)
    ax.plot([x0, tp[1][0]], [y0, tp[1][1]],
            linewidth=ray_lw, color=ray_color, solid_capstyle='round', zorder=4)
    if fill:
        ax.fill([x0, tp[0][0], tp[1][0]],
                [y0, tp[0][1], tp[1][1]],
                alpha=fill_alpha, facecolor=ray_color, linewidth=0, zorder=1)

# --------------------------------- Main ------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Visualize full-track FOV overlay.")
    ap.add_argument("fov_json", type=str, help="Path to fov_cache.json")
    ap.add_argument("--track", type=str, default=None, help="Override track CSV path (x,y,wr,wl)")
    ap.add_argument("--stride", type=int, default=10, help="Plot every N-th FOV to reduce clutter")
    ap.add_argument("--fill", action="store_true", help="Fill FOV triangle in static overlay")
    ap.add_argument("--range", type=str, default=None, help="Limit s-range, e.g. '0:500'")
    ap.add_argument("--save", type=str, default=None, help="Path to save static figure (png/pdf)")
    ap.add_argument("--gif", type=str, default=None, help="Path to save animated GIF (PillowWriter)")
    ap.add_argument("--mp4", type=str, default=None, help="Path to save MP4 (FFMpegWriter)")
    ap.add_argument("--fps", type=int, default=20, help="Frames per second for animation")
    ap.add_argument("--car_radius", type=float, default=1.0, help="Marker radius for the animated car")
    ap.add_argument("--no_show", action="store_true", help="Do not open an interactive window (useful when saving)")
    ap.add_argument("--anim_no_uz", action="store_true", help="Disable unknown-zone polygons in animation (faster)")
    ap.add_argument("--anim_no_fill", action="store_true", help="Disable FOV triangle fill in animation (faster)")
    ap.add_argument("--max_frames", type=int, default=None, help="Hard cap for number of animation frames (subsample evenly)")
    ap.add_argument("--track_fill", action="store_true", help="Fill the road ribbon between left/right bounds")
    ap.add_argument("--track_color", type=str, default="#B0B7C3", help="Road fill color")
    ap.add_argument("--track_alpha", type=float, default=0.35, help="Road fill alpha")
    ap.add_argument("--bounds_lw", type=float, default=2.2, help="Line width for left/right bounds")
    ap.add_argument("--center_lw", type=float, default=1.0, help="Line width for centerline")
    ap.add_argument("--track_expand", type=float, default=0.0, help="Visual inflation of road half-widths in meters (for display only)")

    # Figure size / DPI / output pixel control
    ap.add_argument("--fig_w", type=float, default=12.0, help="Figure width in inches")
    ap.add_argument("--fig_h", type=float, default=8.0, help="Figure height in inches")
    ap.add_argument("--dpi", type=int, default=120, help="Base figure DPI for display/static export")
    ap.add_argument("--gif_dpi", type=int, default=None, help="DPI for GIF/MP4 export (overrides --dpi for animation)")
    ap.add_argument("--out_px_w", type=int, default=None, help="Target output width in pixels (overrides fig size for export)")
    ap.add_argument("--out_px_h", type=int, default=None, help="Target output height in pixels (overrides fig size for export)")

    # Visual styling
    ap.add_argument("--ray_lw", type=float, default=2.5, help="Line width for FOV rays")
    ap.add_argument("--border_lw", type=float, default=2.0, help="Line width for ITP–FPV border in static overlay")
    ap.add_argument("--car_ms", type=float, default=0.0, help="Marker size for the car (if 0, derive from --car_radius)")
    ap.add_argument("--bounds_outline", action="store_true", help="Draw a white outline under the road bounds for contrast")
    ap.add_argument("--bounds_outline_w", type=float, default=2.0, help="Width added to bounds lines for white outline")
    ap.add_argument("--fov_fill_alpha", type=float, default=0.15, help="Alpha for FOV triangle fill")
    ap.add_argument("--show_rays", action="store_true", help="Show ITP-to-FPV FOV rays in animation (off by default)")
    ap.add_argument("--uz_alpha", type=float, default=0.8, help="Alpha for unknown-zone polygons")
    ap.add_argument("--uz_color", type=str, default="#E76F51", help="Color for unknown-zone polygons")

    # Corner notation (legend) for animation
    ap.add_argument("--no_notation", action="store_true", help="Disable corner notation legend in animation")
    ap.add_argument("--notation_loc", type=str, default="upper right",
                    help="Legend corner: 'upper right','upper left','lower right','lower left'")

    # Animated trail
    ap.add_argument("--trail_len", type=int, default=0, help="Length of the motion trail (number of past points). 0 disables.")
    ap.add_argument("--trail_lw", type=float, default=2.0, help="Line width of the trail")
    ap.add_argument("--trail_color", type=str, default="#444444", help="Color of the motion trail")

    # Follow camera window
    ap.add_argument("--follow_cam", action="store_true", help="Enable moving camera window centered at the car")
    ap.add_argument("--cam_w", type=float, default=200.0, help="Follow-cam window width (m)")
    ap.add_argument("--cam_h", type=float, default=200.0, help="Follow-cam window height (m)")
    ap.add_argument("--cam_margin", type=float, default=10.0, help="Extra margin added to follow-cam window (m)")

    # Sector-shaped FOV (local, view-oriented)
    ap.add_argument("--sector", action="store_true", help="Draw a small sector FOV centered at the car heading")
    ap.add_argument("--sector_r", type=float, default=90.0, help="Sector radius (m)")
    ap.add_argument("--sector_half_deg", type=float, default=35.0, help="Half opening angle of the sector (deg)")
    ap.add_argument("--sector_color", type=str, default="#2D9CDB", help="Sector fill color")
    ap.add_argument("--sector_alpha", type=float, default=0.22, help="Sector fill alpha")
    args = ap.parse_args()

    # Load JSON
    with open(args.fov_json, "r", encoding="utf-8") as f:
        cache = json.load(f)

    # Choose DPI for animation vs static
    export_anim = bool('--gif' in sys.argv or '--mp4' in sys.argv or args.gif or args.mp4)
    dpi_anim = args.gif_dpi if (args.gif_dpi is not None) else args.dpi

    # If exact pixel size requested, compute inches from dpi used for export/display
    dpi_for_size = dpi_anim if export_anim else args.dpi
    fig_w_in = args.fig_w
    fig_h_in = args.fig_h
    if args.out_px_w is not None and args.out_px_h is not None:
        fig_w_in = max(1e-3, args.out_px_w / float(dpi_for_size))
        fig_h_in = max(1e-3, args.out_px_h / float(dpi_for_size))

    # Track CSV path
    track_csv = args.track or cache.get("source_csv")
    if track_csv is None:
        raise ValueError("Track CSV path not found. Provide --track or ensure 'source_csv' in JSON.")

    # Load track & derived geometry
    x, y, wr, wl = load_track(track_csv)
    s_track = cumulative_s(x, y)
    xL, yL, xR, yR = build_boundaries(x, y, wl, wr)
    # Heading along centerline (unwrap to ensure smooth interpolation)
    theta_track = heading_theta(x, y, s_track)
    theta_track = np.unwrap(theta_track)

    # Visual-only widening of the road (does not affect geometry/FOV)
    xL_vis, yL_vis, xR_vis, yR_vis = xL, yL, xR, yR
    if abs(args.track_expand) > 1e-9:
        xL_vis, yL_vis, xR_vis, yR_vis = build_boundaries(x, y, wl + args.track_expand, wr + args.track_expand)

    # Build a single closed polygon for the road ribbon (left forward, right backward)
    road_poly = np.concatenate([
        np.column_stack([xL_vis, yL_vis]),
        np.column_stack([xR_vis, yR_vis])[::-1, :]
    ], axis=0)

    # s-range filter
    s_min, s_max = float(s_track[0]), float(s_track[-1])
    if args.range:
        parts = args.range.split(":")
        if len(parts) == 2:
            if parts[0] != "":
                s_min = float(parts[0])
            if parts[1] != "":
                s_max = float(parts[1])

    entries = cache.get("entries", cache)  # support top-level list as well

    # Filter by s-range then stride for overlay density
    filtered = [e for e in entries if s_min <= float(e.get("s0", 0.0)) <= s_max]
    if args.stride > 1:
        filtered = filtered[::args.stride]

    if len(filtered) == 0:
        print("[WARN] No entries after range/stride filtering. A bare track will be shown/saved.")

    # ---------------- Pre-parse all frames once (avoid heavy per-frame parsing) ----------------
    frames_data: List[dict] = []
    for e in filtered:
        s0 = float(e.get("s0", 0.0))
        p0 = interp_xy(s0, s_track, x, y)
        # heading at s0
        theta0 = float(np.interp(s0, s_track, theta_track))
        # tangent points
        tps: List[Tuple[float, float]] = []
        for tp in e.get("tangent_points", []):
            pt = tp.get("point", None) if isinstance(tp, dict) else None
            if pt and isinstance(pt, (list, tuple)) and len(pt) == 2:
                tps.append((float(pt[0]), float(pt[1])))
            elif isinstance(tp, dict) and ('x' in tp and 'y' in tp):
                tps.append((float(tp['x']), float(tp['y'])))
        # unknown zones
        uz_polys: List[np.ndarray] = []
        uzs = e.get("unknown_zones", [])
        if isinstance(uzs, (list, dict)):
            uz_list = uzs if isinstance(uzs, list) else [uzs]
            for uz in uz_list:
                uz_polys.extend(parse_unknown_zone_polys(uz))

        # fov border (unchanged)
        fov_border_pts = None
        fb = e.get("fov_border", {})
        pts = fb.get("points", []) if isinstance(fb, dict) else []
        if isinstance(pts, list) and len(pts) >= 2:
            try:
                fov_border_pts = np.asarray(pts, dtype=float)
            except Exception:
                fov_border_pts = None

        # min distance to any tangent point
        min_tp_dist = None
        if len(tps) > 0:
            try:
                dists = [np.hypot(tp[0]-p0[0], tp[1]-p0[1]) for tp in tps]
                if len(dists) > 0:
                    min_tp_dist = float(min(dists))
            except Exception:
                min_tp_dist = None

        frames_data.append({
            's0': s0,
            'p0': p0,
            'theta0': theta0,
            'tps': tps,
            'uz_polys': uz_polys,
            'fov_border_pts': fov_border_pts,
            'min_tp_dist': min_tp_dist
        })

    print(f"[INFO] Frames prepared: {len(frames_data)} (range=({s_min:.1f},{s_max:.1f}), stride={args.stride})")

    # ------------------------------ Build figure ------------------------------
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=args.dpi)
    ax = plt.gca()
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Draw road: optional filled ribbon + boundary outlines + centerline
    if args.track_fill:
        road_patch = Polygon(road_poly, closed=True, facecolor=args.track_color,
                             edgecolor='none', alpha=args.track_alpha, zorder=0)
        ax.add_patch(road_patch)

    # Optional white outline under bounds for contrast
    if args.bounds_outline:
        ax.plot(xL_vis, yL_vis, linewidth=args.bounds_lw + args.bounds_outline_w,
                color='white', alpha=1.0, zorder=1, solid_capstyle='round')
        ax.plot(xR_vis, yR_vis, linewidth=args.bounds_lw + args.bounds_outline_w,
                color='white', alpha=1.0, zorder=1, solid_capstyle='round')

    # Boundary outlines on top of fill
    ax.plot(xL_vis, yL_vis, linewidth=args.bounds_lw, color='k', alpha=0.95,
            label='left/right bounds', zorder=2, solid_capstyle='round')
    ax.plot(xR_vis, yR_vis, linewidth=args.bounds_lw, color='k', alpha=0.95,
            zorder=2, solid_capstyle='round')

    # Centerline dashed
    ax.plot(x, y, linewidth=args.center_lw, linestyle='--', color='k', alpha=0.6, label='centerline', zorder=3)

    # Static sparse overlay of FOV (blue rays/borders disabled; keep unknown-zone shading only)
    for fd in frames_data:
        # Only render unknown-zone polygons for context; skip static FOV rays/borders
        for poly in fd['uz_polys']:
            ax.fill(poly[:, 0], poly[:, 1], alpha=args.uz_alpha, facecolor=args.uz_color,
                    linewidth=0, zorder=1)

    ax.autoscale_view()
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('Track with FOV Overlay')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(True)
    # ---------------- Legend/Notation handling ----------------

    want_anim = bool(args.gif or args.mp4)
    # Corner notation (legend)
    # If exporting animation and not disabled, show a compact legend in a corner
    if want_anim and not args.no_notation:
        # Remove any previous legend
        leg_prev = ax.get_legend()
        if leg_prev is not None:
            try:
                leg_prev.remove()
            except Exception:
                pass
        handles = []
        # Road bounds and centerline (proxies)
        handles.append(Line2D([0], [0], color='k', lw=args.bounds_lw, label='Bounds'))
        handles.append(Line2D([0], [0], color='k', lw=args.center_lw, ls='--', label='Centerline'))
        # Dynamic FOV rays (optional)
        if args.show_rays:
            handles.append(Line2D([0], [0], color='C0', lw=args.ray_lw, label='FOV rays'))
        # No FOV triangle fill in this mode; sector remains as the FOV depiction
        # Heading-oriented sector (optional)
        if args.sector:
            handles.append(Patch(facecolor=args.sector_color, alpha=args.sector_alpha, edgecolor='none', label='Heading sector'))
        # Motion trail (optional)
        if args.trail_len and args.trail_len > 0:
            handles.append(Line2D([0], [0], color=args.trail_color, lw=args.trail_lw, label='Trail'))
        # Car marker
        car_ms_eff = args.car_ms if (args.car_ms and args.car_ms > 0) else (args.car_radius * 2.5)
        handles.append(Line2D([0], [0], marker='o', color='C3', linestyle='None', markersize=car_ms_eff, label='Car'))

        ax.legend(handles=handles,
                  loc=args.notation_loc,
                  frameon=True, framealpha=0.92,
                  fancybox=True, borderpad=0.6,
                  fontsize=9)
    # For static export / interactive view (no animation), keep a simple legend if not present
    if (not want_anim) and (ax.get_legend() is None):
        ax.legend(loc='best')

    # ------------------------------ Animation ---------------------------------
    if want_anim:
        # Even subsample if too many frames requested
        frames_seq: Iterable[dict] = frames_data
        if args.max_frames is not None and len(frames_data) > args.max_frames:
            idxs = np.linspace(0, len(frames_data) - 1, args.max_frames, dtype=int)
            frames_seq = [frames_data[i] for i in idxs]
            print(f"[INFO] Animation frames subsampled to {len(frames_seq)} (max_frames)")

        if len(list(frames_seq)) == 0:
            print("[WARN] No frames available for animation after filtering/subsampling.")
        else:
            # Prepare dynamic artists
            car_size = args.car_ms if (args.car_ms and args.car_ms > 0) else (args.car_radius * 2.5)
            car_pt, = ax.plot([], [], marker='o', markersize=car_size, linestyle='None',
                              color='C3', zorder=5)
            ray1, = ax.plot([], [], linewidth=args.ray_lw, color='C0', zorder=5, solid_capstyle='round')
            ray2, = ax.plot([], [], linewidth=args.ray_lw, color='C0', zorder=5, solid_capstyle='round')
            tri_patch = Polygon([[0, 0], [0, 0], [0, 0]], closed=True,
                                facecolor='C0', edgecolor='none',
                                alpha=args.fov_fill_alpha, visible=False, zorder=2)
            ax.add_patch(tri_patch)
            # Optional heading-oriented sector FOV
            sector_patch = Polygon([[0, 0], [0, 0], [0, 0]], closed=True,
                                   facecolor=args.sector_color, edgecolor='none',
                                   alpha=args.sector_alpha, visible=False, zorder=3)
            ax.add_patch(sector_patch)
            # Optional motion trail
            trail_line, = ax.plot([], [], color=args.trail_color, linewidth=args.trail_lw,
                                  alpha=0.9, zorder=4, solid_capstyle='round')
            from collections import deque
            trail_buf: Deque[Tuple[float, float]] = deque(maxlen=max(0, args.trail_len))
            uz_artists: List[Polygon] = []

            # Hide rays/triangle by default unless --show_rays is set
            if not args.show_rays:
                ray1.set_visible(False)
                ray2.set_visible(False)
                tri_patch.set_visible(False)

            def init_anim():
                car_pt.set_data([], [])
                ray1.set_data([], [])
                ray2.set_data([], [])
                tri_patch.set_visible(False)
                sector_patch.set_visible(False)
                trail_line.set_data([], [])
                for a in uz_artists:
                    a.remove()
                uz_artists.clear()
                return car_pt, ray1, ray2, tri_patch, sector_patch, trail_line

            def update_anim(fd: dict):
                # Clear unknown zone artists from previous frame
                for a in uz_artists:
                    a.remove()
                uz_artists.clear()
                tri_patch.set_visible(False)

                (x0, y0) = fd['p0']
                car_pt.set_data([x0], [y0])
                # Update motion trail
                if trail_buf.maxlen and trail_buf.maxlen > 0:
                    trail_buf.append((x0, y0))
                    xs = [p[0] for p in trail_buf]
                    ys = [p[1] for p in trail_buf]
                    trail_line.set_data(xs, ys)
                else:
                    trail_line.set_data([], [])

                tps = fd['tps']
                # FOV rays/triangle (optional)
                if args.show_rays:
                    if len(tps) >= 1:
                        ray1.set_data([x0, tps[0][0]], [y0, tps[0][1]])
                    else:
                        ray1.set_data([], [])
                    if len(tps) >= 2:
                        ray2.set_data([x0, tps[1][0]], [y0, tps[1][1]])
                    else:
                        ray2.set_data([], [])
                    if (len(tps) == 2) and (not args.anim_no_fill):
                        tri_patch.set_xy(np.array([[x0, y0], tps[0], tps[1]], dtype=float))
                        tri_patch.set_visible(True)
                    else:
                        tri_patch.set_visible(False)
                    ray1.set_visible(True)
                    ray2.set_visible(True)
                else:
                    ray1.set_data([], [])
                    ray2.set_data([], [])
                    tri_patch.set_visible(False)
                    ray1.set_visible(False)
                    ray2.set_visible(False)

                if not args.anim_no_uz:
                    for poly in fd['uz_polys']:
                        p = Polygon(poly, closed=True, facecolor=args.uz_color, edgecolor='none',
                                    alpha=args.uz_alpha, zorder=1)
                        ax.add_patch(p)
                        uz_artists.append(p)

                # Sector FOV (heading-oriented)
                sector_patch.set_visible(False)
                if args.sector:
                    ang0 = fd.get('theta0', None)
                    if ang0 is not None:
                        # radius: limited by nearest tangent point when available
                        r = float(args.sector_r)
                        if fd.get('min_tp_dist') is not None:
                            try:
                                r = min(r, float(fd['min_tp_dist']))
                            except Exception:
                                pass
                        half = np.deg2rad(float(args.sector_half_deg))
                        # build arc (21 points)
                        angs = np.linspace(ang0 - half, ang0 + half, 21)
                        xs = x0 + r * np.cos(angs)
                        ys = y0 + r * np.sin(angs)
                        sector_xy = np.vstack([[x0, y0], np.column_stack([xs, ys]), [x0, y0]])
                        sector_patch.set_xy(sector_xy)
                        sector_patch.set_visible(True)

                # Follow camera window
                if args.follow_cam:
                    w = float(args.cam_w) * 0.5 + float(args.cam_margin)
                    h = float(args.cam_h) * 0.5 + float(args.cam_margin)
                    ax.set_xlim(x0 - w, x0 + w)
                    ax.set_ylim(y0 - h, y0 + h)

                return car_pt, ray1, ray2, tri_patch, sector_patch, trail_line

            anim = animation.FuncAnimation(
                fig, update_anim, init_func=init_anim,
                frames=frames_seq, interval=1000.0 / float(args.fps), blit=False, repeat=False
            )

            if args.gif:
                try:
                    save_dpi = dpi_anim
                    anim.save(args.gif, writer=animation.PillowWriter(fps=args.fps), dpi=save_dpi)
                    print(f"[OK] Saved GIF to {args.gif}")
                except Exception as ex:
                    print(f"[WARN] Failed to save GIF: {ex}")
            if args.mp4:
                try:
                    save_dpi = dpi_anim
                    anim.save(args.mp4, writer=animation.FFMpegWriter(fps=args.fps, bitrate=2400), dpi=save_dpi)
                    print(f"[OK] Saved MP4 to {args.mp4}")
                except Exception as ex:
                    print(f"[WARN] Failed to save MP4 (need ffmpeg?): {ex}")

    # ------------------------------ Save/Show ---------------------------------
    if args.save:
        try:
            plt.tight_layout()
            fig.savefig(args.save, dpi=args.dpi, bbox_inches='tight')
            print(f"[OK] Saved static figure to {args.save}")
        except Exception as ex:
            print(f"[WARN] Failed to save static figure: {ex}")

    # Show only if no export and not suppressed
    if (not want_anim) and (not args.no_show):
        plt.tight_layout()
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()