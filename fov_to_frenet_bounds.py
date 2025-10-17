# -*- coding: utf-8 -*-
"""FOV/Unknown → Frenet bounds & utilities

This module turns the precomputed FOV dictionary (fov_cache.json) into arrays that the
Frenet OCP can directly use:

A1) FOV fronts (ITP→FPV) → lateral bounds n_min^fov(s_k), n_max^fov(s_k)
    via half-space projection along the normal line p(s,n)=r(s)+n*n_r(s).

A2) Unknown zones (polygons + s-interval + side) → optional clipping on the same
    normal line (hard) and/or a mask for soft-penalty usage.

It also exposes small helpers used by visualization/solvers.
"""
from __future__ import annotations
from typing import List, Tuple, Dict, Any, Optional
import json
import math
import numpy as np

# --------------------------- basic geometry ---------------------------

def heading_normal(theta: float) -> np.ndarray:
    """Left normal of a heading angle theta."""
    return np.array([-math.sin(theta), math.cos(theta)], dtype=float)


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v.copy()


def _perp(v: np.ndarray) -> np.ndarray:
    return np.array([-v[1], v[0]], dtype=float)


def _solve_2x2(A: np.ndarray, b: np.ndarray) -> Optional[np.ndarray]:
    """Solve 2x2; return None if singular."""
    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    if abs(det) < 1e-12:
        return None
    inv = np.array([[ A[1,1], -A[0,1]],
                    [-A[1,0],  A[0,0]]], dtype=float) / det
    return inv @ b


def _line_segment_intersection_param(p0: np.ndarray, d: np.ndarray,
                                     a: np.ndarray, b: np.ndarray) -> Optional[float]:
    """Intersection of ray p(t)=p0+t d with segment a→b.
    Returns parameter t along the ray if it intersects within segment (u∈[0,1]); else None.
    """
    A = np.column_stack((d, -(b - a)))  # [d, -(b-a)]
    sol = _solve_2x2(A, a - p0)
    if sol is None:
        return None
    t, u = float(sol[0]), float(sol[1])
    if u < -1e-12 or u > 1+1e-12:
        return None
    return t


def _point_in_polygon(pt: np.ndarray, poly: List[List[float]]) -> bool:
    """Ray casting; poly is list of [x,y]."""
    x, y = float(pt[0]), float(pt[1])
    inside = False
    n = len(poly)
    if n < 3:
        return False
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i+1)%n]
        # Check if edge crosses the horizontal ray
        cond = ((y1 > y) != (y2 > y))
        if cond:
            xin = (x2 - x1) * (y - y1) / (y2 - y1 + 1e-16) + x1
            if xin > x:
                inside = not inside
    return inside


# --------------------------- FOV → half-space bound ---------------------------

def _halfspace_bound_at_s(pA: np.ndarray, pB: np.ndarray,
                          car_xy_at_s0: np.ndarray,
                          r_s: np.ndarray,
                          n_hat_s: np.ndarray) -> Tuple[float, bool]:
    """Return (n_bound, is_upper) induced by front segment A→B at normal line of (r_s, n_hat_s).

    Half-space orientation is set so that the car position is on the feasible side.
    If the normal is near-parallel to the front, returns (±inf, True) meaning no useful bound.

    is_upper=True means inequality n <= n_bound; False means n >= n_bound.
    """
    t = _unit(pB - pA)
    nN = _perp(t)
    if float(nN @ (car_xy_at_s0 - pA)) > 0:
        nN = -nN
    num = float(nN @ (pA - r_s))
    den = float(nN @ n_hat_s)
    if abs(den) < 1e-12:
        return (math.inf if num >= 0 else -math.inf), True
    n_bound = num / den
    is_upper = den > 0
    return n_bound, is_upper


# --------------------------- cache I/O ---------------------------

def load_cache(path: str) -> Dict[str, Any]:
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


# --------------------------- main constructor ---------------------------

def build_fov_bounds_arrays(*,
    s_grid: np.ndarray,
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    theta: np.ndarray,
    wr: np.ndarray,
    wl: np.ndarray,
    cache_json_path: str,
    blend_mode: str = 'nearest',
    use_unknown_to_clip: bool = True,
    s_tolerance: float = 2.5,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Project FOV fronts and (optionally) Unknown polygons onto the lateral normal at each s_k.

    Returns:
        nMinFov, nMaxFov (np.ndarray, length=len(s_grid)) and aux info dict.
    """
    cache = load_cache(cache_json_path)
    entries: List[Dict[str, Any]] = cache.get('entries', [])
    if not entries:
        # fall back to pure road bounds
        return -wr.copy(), wl.copy(), { 'mode':'road-only', 'picked_idx': np.full_like(s_grid, -1, int) }

    # sort by s0
    s0_list = np.array([float(e.get('s0', 0.0)) for e in entries], dtype=float)
    sort_idx = np.argsort(s0_list)
    s0_list = s0_list[sort_idx]
    entries = [entries[i] for i in sort_idx]

    # prepare outputs (init with road bounds)
    nMin = -wr.astype(float).copy()
    nMax =  wl.astype(float).copy()

    picked_idx = np.empty(len(s_grid), dtype=int)

    # helper: get car XY at a given s0 by nearest on the provided grid
    def car_xy_at(s0: float) -> np.ndarray:
        k0 = int(np.argmin(np.abs(s_grid - s0)))
        return np.array([x_ref[k0], y_ref[k0]], dtype=float)

    for k, s_k in enumerate(s_grid):
        # pick entry/entries
        if blend_mode == 'pairwise_intersect':
            # left-right neighbors
            ri = int(np.searchsorted(s0_list, s_k, side='right'))
            li = max(0, ri-1)
            cand_idx = np.unique([li, min(ri, len(entries)-1)])
        else:  # 'nearest'
            cand_idx = [int(np.argmin(np.abs(s0_list - s_k)))]
        picked_idx[k] = cand_idx[0]

        r_s = np.array([x_ref[k], y_ref[k]], dtype=float)
        n_hat = heading_normal(float(theta[k]))

        # 1) apply FOV fronts from candidate entries
        local_min = -wr[k]
        local_max =  wl[k]
        for j in cand_idx:
            e = entries[j]
            pts = e.get('fov_border', {}).get('points', [])
            if not pts or len(pts) < 2:
                continue
            pA = np.array(pts[0], dtype=float)
            pB = np.array(pts[-1], dtype=float)
            car_xy = car_xy_at(float(e.get('s0', s_k)))
            n_bound, is_upper = _halfspace_bound_at_s(pA, pB, car_xy, r_s, n_hat)
            if math.isfinite(n_bound):
                if is_upper:
                    local_max = min(local_max, n_bound)
                else:
                    local_min = max(local_min, n_bound)

        # 2) optionally clip with Unknown polygons on the corresponding side
        if use_unknown_to_clip:
            for j in cand_idx:
                e = entries[j]
                uz_list = e.get('unknown_zones', []) or []
                for uz in uz_list:
                    s_start = float(uz.get('s_start', -1e9))
                    s_end   = float(uz.get('s_end',   1e9))
                    if s_k < s_start - s_tolerance or s_k > s_end + s_tolerance:
                        continue
                    side = int(uz.get('side', 1))  # 1=left, 0=right
                    poly = uz.get('polygon', [])
                    if not poly or len(poly) < 3:
                        continue
                    # intersect normal line p(t)=r_s + t n_hat with polygon edges; pick nearest on that side
                    best_t: Optional[float] = None
                    for i in range(len(poly)):
                        a = np.array(poly[i], dtype=float)
                        b = np.array(poly[(i+1)%len(poly)], dtype=float)
                        t = _line_segment_intersection_param(r_s, n_hat, a, b)
                        if t is None:
                            continue
                        # keep only left(+) or right(-) side
                        if side == 1 and t <= 1e-9:
                            continue
                        if side == 0 and t >= -1e-9:
                            continue
                        if best_t is None:
                            best_t = t
                        else:
                            # choose the one closest to r_s in that direction
                            if side == 1:  # left: positive small t
                                if 0 < t < best_t:
                                    best_t = t
                            else:          # right: negative, choose the one closer to 0 (less negative)
                                if t < 0 and (t > best_t):
                                    best_t = t
                    if best_t is not None and math.isfinite(best_t):
                        if side == 1:
                            local_max = min(local_max, best_t)
                        else:
                            local_min = max(local_min, best_t)

        nMin[k] = local_min
        nMax[k] = local_max

        # safety: keep order
        if nMin[k] > nMax[k]:
            mid = 0.5*(nMin[k]+nMax[k])
            nMin[k] = nMax[k] = mid

    aux = {
        'mode': blend_mode,
        'picked_idx': picked_idx,
        's0_list': s0_list,
    }
    return nMin, nMax, aux


# --- Compatibility wrapper: historical external name ---
def build_planning_bounds_arrays(*,
    s_grid: np.ndarray,
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    theta: np.ndarray,
    wr: np.ndarray,
    wl: np.ndarray,
    cache_json_path: str,
    blend_mode: str = 'pairwise_intersect',
    include_unknown: bool = True,
    s_tolerance: float = 2.5,
):
    """Alias to build_fov_bounds_arrays with argument mapping.
    include_unknown → use_unknown_to_clip
    """
    return build_fov_bounds_arrays(
        s_grid=s_grid,
        x_ref=x_ref,
        y_ref=y_ref,
        theta=theta,
        wr=wr,
        wl=wl,
        cache_json_path=cache_json_path,
        blend_mode=blend_mode,
        use_unknown_to_clip=include_unknown,
        s_tolerance=s_tolerance,
    )


# --------------------------- unknown mask for a given trajectory ---------------------------

def unknown_mask_for_traj(*,
    s_grid: np.ndarray,
    n_traj: np.ndarray,
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    theta: np.ndarray,
    cache_json_path: str,
    s_tolerance: float = 2.5,
) -> np.ndarray:
    """Return boolean mask M[k]=True if p(s_k,n_k) lies inside any Unknown polygon whose
    s-interval covers s_k. Useful for soft-penalty visualization/evaluation.
    """
    cache = load_cache(cache_json_path)
    entries: List[Dict[str, Any]] = cache.get('entries', [])
    if not entries:
        return np.zeros_like(s_grid, dtype=bool)

    # Flatten unknowns with their s-interval for quick scan
    unk_list: List[Tuple[float,float,int,List[List[float]]]] = []
    for e in entries:
        uzs = e.get('unknown_zones', []) or []
        for uz in uzs:
            poly = uz.get('polygon', [])
            if not poly or len(poly) < 3:
                continue
            unk_list.append((float(uz.get('s_start', -1e9)),
                             float(uz.get('s_end', 1e9)),
                             int(uz.get('side', 1)),
                             poly))

    M = np.zeros(len(s_grid), dtype=bool)
    for k, s_k in enumerate(s_grid):
        r_s = np.array([x_ref[k], y_ref[k]], dtype=float)
        n_hat = heading_normal(float(theta[k]))
        p = r_s + float(n_traj[k]) * n_hat
        for s0, s1, _side, poly in unk_list:
            if s_k < s0 - s_tolerance or s_k > s1 + s_tolerance:
                continue
            if _point_in_polygon(p, poly):
                M[k] = True
                break
    return M

# --------------------------- direct OPTI integration helper ---------------------------

def apply_planning_band_to_opti(*,
    opti,
    n_var,
    s_grid: np.ndarray,
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    theta: np.ndarray,
    wr: np.ndarray,
    wl: np.ndarray,
    cache_json_path: str,
    mode: str = 'hard',
    lam_fov: float = 5e2,
    blend_mode: str = 'pairwise_intersect',
    include_unknown: bool = True,
):
    """Enforce the convex planning band on the Opti model.

    Args:
        opti: casadi.Opti() instance
        n_var: the decision variable array for lateral deviation n at all grid nodes (shape (N,) or (1,N))
        mode: 'hard' adds subject_to bounds; 'soft' adds quadratic penalty slack variables
        lam_fov: weight for the soft-wall penalty
    Returns:
        A dict with arrays (nMinPlan, nMaxPlan) and references to any created slack vars (soft mode).
    """
    try:
        import casadi as ca
    except Exception:
        raise RuntimeError("apply_planning_band_to_opti requires casadi to be installed.")

    n_var = ca.reshape(n_var, (-1,))  # flatten to (N,)
    nMin, nMax, aux = build_planning_bounds_arrays(
        s_grid=s_grid, x_ref=x_ref, y_ref=y_ref, theta=theta,
        wr=wr, wl=wl,
        cache_json_path=cache_json_path,
        blend_mode=blend_mode,
        include_unknown=include_unknown,
    )

    if mode == 'hard':
        for k in range(len(s_grid)):
            opti.subject_to(n_var[k] >= nMin[k])
            opti.subject_to(n_var[k] <= nMax[k])
        return { 'nMinPlan': nMin, 'nMaxPlan': nMax, 'aux': aux }

    elif mode == 'soft':
        # add slack variables and quadratic penalties
        s_plus  = opti.variable(len(s_grid))
        s_minus = opti.variable(len(s_grid))
        opti.subject_to(s_plus  >= 0)
        opti.subject_to(s_minus >= 0)
        for k in range(len(s_grid)):
            opti.subject_to(n_var[k] - nMax[k] <= s_plus[k])
            opti.subject_to(nMin[k] - n_var[k] <= s_minus[k])
        penalty = lam_fov * (ca.sumsqr(s_plus) + ca.sumsqr(s_minus))
        return { 'nMinPlan': nMin, 'nMaxPlan': nMax, 'aux': aux, 'soft_penalty': penalty,
                 'slack_plus': s_plus, 'slack_minus': s_minus }

    else:
        raise ValueError("mode must be 'hard' or 'soft'")
