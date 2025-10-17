"""
Extended field-of-view (FOV) computation utilities for high-speed driving.

This module builds upon the simple FOV estimation presented in the
`compute_fov_simple` function by providing a structured pipeline to
detect Interested Tangent Points (ITP) on the left and right track
boundaries, determine whether the current path segment is a bend,
search for the farthest point of view (FPV) by aligning with the
vehicle's heading, and construct unknown zones bounded by the
identified tangent and farthest points.

The functions defined here are intentionally modular so they can
be reused in a larger optimisation framework.  Each function
contains detailed documentation explaining the expected inputs and
algorithmic behaviour.  The implementation follows the pseudo-code
outlined in Section 8.2.4 and Algorithm 8.1 of the referenced
report, but avoids external dependencies beyond NumPy and SciPy.

python plot_fov_cache.py --csv Monza.csv --cache fov_cache.json --k 120 --mode segment
# 如果你用 'front' 生成了折线前沿：
# python plot_fov_cache.py --csv Monza.csv --cache fov_cache.json --k 120 --mode front
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline
from typing import Callable, List, Optional, Tuple, Dict

def build_boundaries(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    w_left: np.ndarray,
    w_right: np.ndarray,
    s_dense: np.ndarray,
) -> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:
    """Compute the left and right boundary coordinates of a track."""
    dx = spline_x.derivative(1)(s_dense)
    dy = spline_y.derivative(1)(s_dense)
    tang_norm = np.hypot(dx, dy)
    tang_norm[tang_norm == 0.0] = 1.0
    nx = -dy / tang_norm
    ny = dx / tang_norm
    xc = spline_x(s_dense)
    yc = spline_y(s_dense)
    x_left = xc - w_left * nx
    y_left = yc - w_left * ny
    x_right = xc + w_right * nx
    y_right = yc + w_right * ny
    return (x_left, y_left), (x_right, y_right)

def _angle_monotonicity_itp(
    xb: np.ndarray,
    yb: np.ndarray,
    s_segment: np.ndarray,
    p_car: np.ndarray,
    dir_v: np.ndarray,
    want_max: bool,
    ang_min: float,
    peak_min: float,
    min_dist: float,
    min_arc: float,
    s0: float
) -> float:
    """Internal helper to detect the interested tangent point (ITP)."""
    Vx = xb - p_car[0]
    Vy = yb - p_car[1]
    dist = np.hypot(Vx, Vy)
    dot = dir_v[0] * Vx + dir_v[1] * Vy
    mask = (dist > min_dist) & (dot > 0)
    if not np.any(mask):
        return np.inf
    Vx = Vx[mask]; Vy = Vy[mask]; dot = dot[mask]
    s_seg = s_segment[mask]
    theta = np.arctan2(dir_v[0] * Vy - dir_v[1] * Vx, dot)
    if np.ptp(theta) < ang_min:
        return np.inf
    dtheta = np.diff(theta)
    sign = np.sign(dtheta)
    flips = np.where(sign[:-1] * sign[1:] < 0)[0] + 1
    cand_idx = flips if flips.size else ([np.argmax(theta)] if want_max else [np.argmin(theta)])
    good = []
    for idx in cand_idx:
        left_gap = abs(theta[idx] - theta[max(idx-1,0)])
        right_gap = abs(theta[idx] - theta[min(idx+1,len(theta)-1)])
        if max(left_gap, right_gap) >= peak_min:
            good.append(idx)
    if not good:
        return np.inf
    idx_sel = min(good, key=lambda i: s_seg[i])
    s_tan = s_seg[idx_sel]
    if s_tan - s0 < min_arc:
        return np.inf
    return s_tan

def find_ITP_by_monotonicity(
    boundary: Tuple[np.ndarray, np.ndarray],
    s_window: np.ndarray,
    p_car: np.ndarray,
    dir_v: np.ndarray,
    want_max: bool,
    theta_min_deg: float,
    peak_min_deg: float,
    min_dist: float,
    min_arc: float,
    s0: float
) -> float:
    ang_min = np.deg2rad(theta_min_deg)
    peak_min = np.deg2rad(peak_min_deg)
    return _angle_monotonicity_itp(
        xb=boundary[0], yb=boundary[1],
        s_segment=s_window,
        p_car=p_car, dir_v=dir_v,
        want_max=want_max,
        ang_min=ang_min, peak_min=peak_min,
        min_dist=min_dist, min_arc=min_arc, s0=s0
    )           

def find_TPs_by_monotonicity(
    boundary: Tuple[np.ndarray, np.ndarray],
    s_window: np.ndarray,
    p_car: np.ndarray,
    dir_v: np.ndarray,
    want_max: bool,
    theta_min_deg: float,
    peak_min_deg: float,
    min_dist: float,
    min_arc: float,
    s0: float
) -> List[float]:
    """
    Find all tangent points (TPs) on the boundary by monotonicity of angle theta(s).
    Returns a sorted list of s values of candidate tangent points.
    """
    xb, yb = boundary
    Vx = xb - p_car[0]
    Vy = yb - p_car[1]
    dist = np.hypot(Vx, Vy)
    dot = dir_v[0] * Vx + dir_v[1] * Vy
    mask = (dist > min_dist) & (dot > 0)
    if not np.any(mask):
        return []
    Vx = Vx[mask]
    Vy = Vy[mask]
    dot = dot[mask]
    s_seg = s_window[mask]
    theta = np.arctan2(dir_v[0] * Vy - dir_v[1] * Vx, dot)
    if np.ptp(theta) < np.deg2rad(theta_min_deg):
        return []
    dtheta = np.diff(theta)
    sign = np.sign(dtheta)
    flips = np.where(sign[:-1] * sign[1:] < 0)[0] + 1
    # Determine candidate indices of local extrema according to want_max
    if flips.size == 0:
        # No sign flips: whole segment monotonic, take global extremum
        idxs = [np.argmax(theta)] if want_max else [np.argmin(theta)]
    else:
        # For each flip, check if it is local max or min
        idxs = []
        for idx in flips:
            left = theta[idx-1] if idx-1 >= 0 else theta[idx]
            center = theta[idx]
            right = theta[idx+1] if idx+1 < len(theta) else theta[idx]
            if want_max:
                if center > left and center > right:
                    idxs.append(idx)
            else:
                if center < left and center < right:
                    idxs.append(idx)
    # Filter by angular prominence and min_arc
    peak_min = np.deg2rad(peak_min_deg)
    good_idxs = []
    for idx in idxs:
        left_gap = abs(theta[idx] - theta[max(idx-1,0)])
        right_gap = abs(theta[idx] - theta[min(idx+1,len(theta)-1)])
        if max(left_gap, right_gap) >= peak_min:
            s_tan = s_seg[idx]
            if s_tan - s0 >= min_arc:
                good_idxs.append(s_tan)
    return sorted(good_idxs)

def is_bend_at_itp(
    s_itp: float,
    kappa_spline: Optional[CubicSpline],
    kappa_min: float
) -> bool:
    """Check whether the track has significant curvature at the ITP."""
    if kappa_spline is None or s_itp == np.inf:
        return True
    kappa_val = abs(kappa_spline(s_itp))
    return kappa_val >= kappa_min


def find_fpv_by_alignment(
    tangent_boundary: Tuple[np.ndarray, np.ndarray],
    opposite_boundary: Tuple[np.ndarray, np.ndarray],
    s_window: np.ndarray,
    p_car: np.ndarray,
    dir_v: np.ndarray,
    s_itp: float,
    epsilon_dir: float
) -> Tuple[float, bool]:
    """Search for the farthest point of view (FPV) by alignment."""
    if s_itp == np.inf:
        return np.inf, False
    itp_idx = np.searchsorted(s_window, s_itp, side="left")
    itp_idx = max(0, min(itp_idx, len(s_window) - 1))
    x_tan_itp = tangent_boundary[0][itp_idx]
    y_tan_itp = tangent_boundary[1][itp_idx]
    tan_dir = np.array([x_tan_itp - p_car[0], y_tan_itp - p_car[1]])
    norm = np.hypot(tan_dir[0], tan_dir[1])
    if norm == 0:
        return np.inf, False
    tan_dir /= norm
    x_opp, y_opp = opposite_boundary
    best_s = np.inf
    found = False
    for idx in range(itp_idx, len(s_window)):
        vx = x_opp[idx] - p_car[0]
        vy = y_opp[idx] - p_car[1]
        dot = tan_dir[0] * vx + tan_dir[1] * vy
        if dot <= 0:
            continue
        cross = abs(tan_dir[0] * vy - tan_dir[1] * vx)
        if cross / dot < epsilon_dir:
            best_s = s_window[idx]
            found = True
        elif found:
            break
    return best_s, found


def build_unknown_zone(
    p_car: np.ndarray,
    boundary: Tuple[np.ndarray, np.ndarray],
    s_start: float,
    s_end: float,
    s_window: np.ndarray
) -> List[Tuple[float, float]]:
    """Construct a polygon approximating the unknown zone."""
    if np.isinf(s_start) or np.isinf(s_end):
        return []
    idx_start = np.searchsorted(s_window, s_start, side="left")
    idx_start = max(0, min(idx_start, len(s_window) - 1))
    idx_end = np.searchsorted(s_window, s_end, side="left")
    idx_end = max(0, min(idx_end, len(s_window) - 1))
    x_b, y_b = boundary
    start_point = (x_b[idx_start], y_b[idx_start])
    end_point = (x_b[idx_end], y_b[idx_end])
    return [tuple(p_car), start_point, end_point, tuple(p_car)]

def _first_aligned_hit_on_boundary(
    x_b: np.ndarray,
    y_b: np.ndarray,
    s_window: np.ndarray,
    p_car: np.ndarray,
    ray_dir: np.ndarray,
    start_idx: int,
    epsilon_dir: float
) -> Optional[float]:
    """
    在 boundary 上从 start_idx 往前找第一个与 ray_dir 方向对齐的点；
    返回其弧长 s。若找不到，返回 None。
    """
    for idx in range(start_idx, len(s_window)):
        vx = x_b[idx] - p_car[0]
        vy = y_b[idx] - p_car[1]
        dot = ray_dir[0] * vx + ray_dir[1] * vy
        if dot <= 0:
            continue  # 车后方/切线反向
        cross = abs(ray_dir[0] * vy - ray_dir[1] * vx)
        if dot != 0 and cross / dot < epsilon_dir:
            return s_window[idx]
    return None


# --- New helper: alignment from ITP point, not car ---
def _first_aligned_from_itp(
    x_b: np.ndarray,
    y_b: np.ndarray,
    s_window: np.ndarray,
    itp_pt: np.ndarray,
    ray_dir: np.ndarray,
    start_idx: int,
    epsilon_dir: float
) -> Optional[float]:
    """
    从 start_idx+1 开始，遍历 boundary 离散点，对于每个点计算向量 v = [x_b[idx]-itp_pt[0], y_b[idx]-itp_pt[1]]。
    计算 dot = ray_dir @ v, cross = |ray_dir_x*v_y - ray_dir_y*v_x|。
    若 dot>0 且 cross/dot < epsilon_dir，则返回 s_window[idx]；否则继续。若找不到，返回 None。
    """
    for idx in range(start_idx + 1, len(s_window)):
        vx = x_b[idx] - itp_pt[0]
        vy = y_b[idx] - itp_pt[1]
        dot = ray_dir[0] * vx + ray_dir[1] * vy
        if dot <= 0:
            continue
        cross = abs(ray_dir[0] * vy - ray_dir[1] * vx)
        if dot != 0 and cross / dot < epsilon_dir:
            return s_window[idx]
    return None


def compute_fov_full(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    kappa_spline: Optional[CubicSpline],
    w_left: np.ndarray,
    w_right: np.ndarray,
    s0: float,
    delta_s: float = 200.0,
    s_dense_step: float = 0.5,
    theta_min_deg: float = 4.0,
    peak_min_deg: float = 0.0,
    min_dist: float = 5.0,
    min_arc: float = 10.0,
    kappa_min: float = 1e-5,
    epsilon_dir: float = 1e-3,
    fov_border_mode: str = 'segment'
) -> Dict[str, object]:
    """
    单次 FOV 计算：
    - 找两侧 ITP（角度单调极值）
    - 每个 ITP 仅在切点处做一次射线判别：先同侧=Unknown Zone，先对侧=潜在 FPV
    - 在所有潜在 FPV 中取弧长最小者为最终 FPV
    返回：
      {
        's0', 's_itp', 's_fpv', 'side', 'unknown_zones', 'fov_polygon', 's_fov'
      }
    参数 fov_border_mode: str, 可选，指定 fov_border 的类型：
      - 'segment'：仅用 ITP–FPV 的直线段作为可见域前沿
      - 'front'：通过角度扫描获得离散可见前沿
    """
    # 采样窗口
    s_dense = np.arange(s0, s0 + delta_s + s_dense_step, s_dense_step)
    # 插宽度
    try:
        s_global = spline_x.x
        w_left_dense = np.interp(s_dense, s_global, w_left)
        w_right_dense = np.interp(s_dense, s_global, w_right)
    except AttributeError:
        w_left_dense = np.full_like(s_dense, w_left if np.isscalar(w_left) else w_left[0])
        w_right_dense = np.full_like(s_dense, w_right if np.isscalar(w_right) else w_right[0])

    # 边界采样
    (x_left, y_left), (x_right, y_right) = build_boundaries(
        spline_x, spline_y, w_left_dense, w_right_dense, s_dense
    )

    # 车辆位置与朝向
    p_car = np.array([spline_x(s0), spline_y(s0)])
    t_vec = np.array([spline_x.derivative(1)(s0), spline_y.derivative(1)(s0)])
    t_norm = np.hypot(t_vec[0], t_vec[1]) or 1.0
    dir_v = t_vec / t_norm
    s_window = s_dense

    # 两侧 ITP - 修改为找所有切点
    tps_left = find_TPs_by_monotonicity(
        (x_left, y_left), s_window, p_car, dir_v,
        want_max=True, theta_min_deg=theta_min_deg, peak_min_deg=peak_min_deg,
        min_dist=min_dist, min_arc=min_arc, s0=s0
    )
    tps_right = find_TPs_by_monotonicity(
        (x_right, y_right), s_window, p_car, dir_v,
        want_max=False, theta_min_deg=theta_min_deg, peak_min_deg=peak_min_deg,
        min_dist=min_dist, min_arc=min_arc, s0=s0
    )

    # 结果骨架
    result: Dict[str, object] = {
        's0': s0,
        's_itp': None,
        's_fpv': None,
        'side': None,              # 0=右边界, 1=左边界
        'unknown_zones': [],       # 多个 unknown zone
        'fov_polygon': [],         # 车→ITP→最终FPV 三角示意
        'fov_border': {},          # 离散版可视域前沿（后续填充）
        's_fov': s0 + delta_s,
        'tangent_points': []       # 所有检测到的切点列表
    }

    pot_fpvs: List[Tuple[float, int, float, int, int]] = []  # (s_itp, fpv_side, s_fpv, itp_side, itp_idx)

    # 处理左侧所有切点
    for s_itp in tps_left:
        if not is_bend_at_itp(s_itp, kappa_spline, kappa_min):
            continue
        itp_idx = int(np.searchsorted(s_window, s_itp, side='left'))
        itp_idx = max(0, min(itp_idx, len(s_window)-1))
        itp_point = np.array([x_left[itp_idx], y_left[itp_idx]])
        ray_dir = itp_point - p_car
        n = np.hypot(ray_dir[0], ray_dir[1]) or 1.0
        ray_dir /= n

        result['tangent_points'].append({'s': s_itp, 'side': 1, 'point': (x_left[itp_idx], y_left[itp_idx])})

        # Alignment is now checked between the car→ITP reference ray and the vectors ITP→boundary-point (discrete)
        hit_same = _first_aligned_from_itp(x_left, y_left, s_window, itp_point, ray_dir, itp_idx, epsilon_dir)
        hit_opp  = _first_aligned_from_itp(x_right, y_right, s_window, itp_point, ray_dir, itp_idx, epsilon_dir)

        if hit_same is not None and (hit_opp is None or hit_same <= hit_opp):
            poly = build_unknown_zone(p_car, (x_left, y_left), s_itp, hit_same, s_window)
            result['unknown_zones'].append({
                's_start': s_itp,
                's_end': hit_same,
                'side': 1,
                'polygon': poly
            })
        elif hit_opp is not None:
            pot_fpvs.append((s_itp, 0, hit_opp, 1, itp_idx))

    # 处理右侧所有切点
    for s_itp in tps_right:
        if not is_bend_at_itp(s_itp, kappa_spline, kappa_min):
            continue
        itp_idx = int(np.searchsorted(s_window, s_itp, side='left'))
        itp_idx = max(0, min(itp_idx, len(s_window)-1))
        itp_point = np.array([x_right[itp_idx], y_right[itp_idx]])
        ray_dir = itp_point - p_car
        n = np.hypot(ray_dir[0], ray_dir[1]) or 1.0
        ray_dir /= n

        result['tangent_points'].append({'s': s_itp, 'side': 0, 'point': (x_right[itp_idx], y_right[itp_idx])})

        # Alignment is now checked between the car→ITP reference ray and the vectors ITP→boundary-point (discrete)
        hit_same = _first_aligned_from_itp(x_right, y_right, s_window, itp_point, ray_dir, itp_idx, epsilon_dir)
        hit_opp  = _first_aligned_from_itp(x_left, y_left, s_window, itp_point, ray_dir, itp_idx, epsilon_dir)

        if hit_same is not None and (hit_opp is None or hit_same <= hit_opp):
            poly = build_unknown_zone(p_car, (x_right, y_right), s_itp, hit_same, s_window)
            result['unknown_zones'].append({
                's_start': s_itp,
                's_end': hit_same,
                'side': 0,
                'polygon': poly
            })
        elif hit_opp is not None:
            pot_fpvs.append((s_itp, 1, hit_opp, 0, itp_idx))

    # 选择最终 FPV（潜在里取 s_fpv 最小）
    if pot_fpvs:
        pot_fpvs.sort(key=lambda t: t[2])
        s_itp_sel, fpv_side, s_fpv_sel, itp_side, itp_idx_sel = pot_fpvs[0]
        result['s_itp'] = s_itp_sel
        result['s_fpv'] = s_fpv_sel
        result['side'] = fpv_side
        result['s_fov'] = s_fpv_sel

        # 构造 ITP / FPV 坐标
        fpv_idx = int(np.searchsorted(s_window, s_fpv_sel, side='left'))
        fpv_idx = max(0, min(fpv_idx, len(s_window)-1))
        if fpv_side == 0:  # FPV 在右
            fpv_pt = (x_right[fpv_idx], y_right[fpv_idx])
        else:              # FPV 在左
            fpv_pt = (x_left[fpv_idx], y_left[fpv_idx])

        if itp_side == 1:
            itp_pt = (x_left[itp_idx_sel], y_left[itp_idx_sel])
        else:
            itp_pt = (x_right[itp_idx_sel], y_right[itp_idx_sel])

        # FOV 三角示意（车→ITP→FPV→车）
        result['fov_polygon'] = [tuple(p_car), itp_pt, fpv_pt, tuple(p_car)]

    # 填充 fov_border 按所选模式
    if result.get('s_fpv') is not None and result.get('s_itp') is not None and result.get('side') is not None:
        if fov_border_mode == 'segment':
            # 仅使用 ITP–FPV 的直线段作为 fov_border（符合你的定义：可见域前沿=这条射线）
            result['fov_border'] = {
                'mode': 'segment',
                'points': [itp_pt, fpv_pt],
                's_on_boundary': [float(result['s_itp']), float(result['s_fpv'])],
                'side': [int(itp_side), int(result['side'])],
                'angles': []
            }
        elif fov_border_mode == 'front':
            # 扫描角度得到非三角可见前沿（不依赖 ITP）
            result['fov_border'] = compute_visibility_border_discrete(
                spline_x=spline_x,
                spline_y=spline_y,
                w_left=w_left,
                w_right=w_right,
                s0=s0,
                look=delta_s,
                unknown_zones=result.get('unknown_zones', []),
                fpv_side=result.get('side'),
                s_fpv=result.get('s_fpv'),
                s_step=s_dense_step,
                n_angles=721,
            )
    else:
        result['fov_border'] = {}

    return result


def compute_fov_range(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    kappa_spline: Optional[CubicSpline],
    w_left: np.ndarray,
    w_right: np.ndarray,
    s0: float,
    search_range: float = 200.0,
    s_dense_step: float = 0.5,
    theta_min_deg: float = 4.0,
    peak_min_deg: float = 0.0,
    min_dist: float = 5.0,
    min_arc: float = 10.0,
    kappa_min: float = 1e-5,
    epsilon_dir: float = 1e-3
) -> Dict[str, object]:
    """
    在 [s0, s0+search_range] 内重复执行单次 FOV 计算，收集所有 Unknown Zone。
    返回：
      {
        's0': s0,
        'unknown_zones': [
            {'s_start','s_end','side','polygon'}, ...
        ]
      }
    """
    unknown_zones: List[Dict[str, object]] = []
    s_cur = s0
    s_end = s0 + search_range

    while s_cur < s_end:
        remaining = s_end - s_cur
        res = compute_fov_full(
            spline_x, spline_y, kappa_spline,
            w_left, w_right,
            s0=s_cur, delta_s=remaining,
            s_dense_step=s_dense_step,
            theta_min_deg=theta_min_deg,
            peak_min_deg=peak_min_deg,
            min_dist=min_dist, min_arc=min_arc,
            kappa_min=kappa_min, epsilon_dir=epsilon_dir,
            fov_border_mode='segment'
        )

        # 若没有 ITP，退出
        if res.get('s_itp') is None and (not res.get('unknown_zones')):
            break

        # 记录未知区
        for uz in res.get('unknown_zones', []):
            unknown_zones.append(uz)

        # 推进 s 指针：优先用 s_fpv，否则用未知区末端，否则退出
        s_next = res.get('s_fpv')
        if s_next is None or not np.isfinite(s_next):
            if res.get('unknown_zones'):
                s_next = max(uz['s_end'] for uz in res['unknown_zones'])
            else:
                break
        s_cur = s_next + max(min_arc, 1e-6)

    return {'s0': s0, 'unknown_zones': unknown_zones}


# ---------------- Discrete visibility border (alignment-threshold style) ----------------

def compute_visibility_border_discrete(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    w_left: np.ndarray,
    w_right: np.ndarray,
    s0: float,
    look: float,
    unknown_zones: List[Dict[str, object]],
    fpv_side: Optional[int],
    s_fpv: Optional[float],
    s_step: float = 0.5,
    n_angles: int = 721,
) -> Dict[str, object]:
    """离散版可视域前沿（非三角）。
    做法：
      1) 在 [s0, s0+look] 上离散出左右边界；
      2) 只看前向半平面，按角度等间隔采样；
      3) 对每个角度 a，在两条边界的离散点里找与射线方向最接近的点且在前向（dot>0），并且不落在 Unknown Zone 的同侧弧段里；
      4) 两侧候选取距离车更近的作为该角度的“首碰点”；
      5) 若存在最终 FPV，则把角度范围裁切在 [0, 车→FPV 的方向] 内；
      6) 返回 points/s_on_boundary/side/angles 四个数组。
    说明：完全采用离散采样与方向对齐近似，不做几何线段相交。
    """
    # 边界离散
    s_dense = np.arange(s0, s0 + look + s_step, s_step)
    try:
        s_global = spline_x.x
        wL = np.interp(s_dense, s_global, w_left)
        wR = np.interp(s_dense, s_global, w_right)
    except Exception:
        wL = np.full_like(s_dense, w_left if np.isscalar(w_left) else w_left[0])
        wR = np.full_like(s_dense, w_right if np.isscalar(w_right) else w_right[0])

    (xL, yL), (xR, yR) = build_boundaries(spline_x, spline_y, wL, wR, s_dense)

    car = np.array([spline_x(s0), spline_y(s0)])
    t = np.array([spline_x.derivative(1)(s0), spline_y.derivative(1)(s0)])
    t /= (np.hypot(*t) or 1.0)

    # 预计算角度/距离（前向）
    def rel(X, Y):
        Vx = X - car[0]; Vy = Y - car[1]
        dot = Vx * t[0] + Vy * t[1]
        fwd = dot > 0
        idx = np.nonzero(fwd)[0]
        if idx.size == 0:
            return np.empty((0,2)), np.empty(0, dtype=int), np.empty(0), np.empty(0)
        V = np.stack([Vx[idx], Vy[idx]], axis=1)
        dotv = V @ t
        cross = t[0]*V[:,1] - t[1]*V[:,0]
        ang = np.arctan2(cross, dotv)
        dist = np.hypot(V[:,0], V[:,1])
        return V, idx, ang, dist

    VL, idxL, angL, distL = rel(xL, yL)
    VR, idxR, angR, distR = rel(xR, yR)

    # Unknown Zone 区间（剔除同侧遮挡）
    uzL = [(uz['s_start'], uz['s_end']) for uz in unknown_zones if uz.get('side') == 1]
    uzR = [(uz['s_start'], uz['s_end']) for uz in unknown_zones if uz.get('side') == 0]

    # 角度范围：[-pi/2, +pi/2] 并可裁到 FPV 方向
    a_min, a_max = -np.pi/2, np.pi/2
    if (fpv_side is not None) and (s_fpv is not None):
        j = int(np.searchsorted(s_dense, s_fpv, side='left'))
        j = max(0, min(j, len(s_dense)-1))
        fpv_pt = (xR[j], yR[j]) if fpv_side == 0 else (xL[j], yL[j])
        v = np.array([fpv_pt[0] - car[0], fpv_pt[1] - car[1]])
        v /= (np.hypot(*v) or 1.0)
        a_fpv = np.arctan2(t[0]*v[1] - t[1]*v[0], t @ v)
        # 只保留 [min(0,a_fpv), max(0,a_fpv)]，通常你希望从车头方向到 FPV 的扇区
        lo, hi = (0.0, a_fpv) if a_fpv >= 0 else (a_fpv, 0.0)
        a_min, a_max = max(a_min, lo), min(a_max, hi)

    angles = np.linspace(a_min, a_max, n_angles)

    border_pts, border_s, border_side, border_ang = [], [], [], []

    for a in angles:
        # 目标射线方向：将车体坐标系绕车头旋转 a
        ray = np.array([np.cos(a) * t[0] - np.sin(a) * t[1],
                        np.sin(a) * t[0] + np.cos(a) * t[1]])
        # 在每侧找“角度最近且不在 Unknown 段”的离散点；取距离最近者
        def pick_side(V, ang, idxs, X, Y, uz_intervals, side_tag):
            if ang.size == 0:
                return None
            k = int(np.argmin(np.abs(ang - a)))
            i = idxs[k]
            # 前向已经由 rel() 保证，这里仅检查 Unknown Zone
            s_i = s_dense[i]
            if any(s0 <= s_i <= s1 for (s0, s1) in uz_intervals):
                return None
            d = np.hypot(X[i]-car[0], Y[i]-car[1])
            return (d, (X[i], Y[i]), s_i, side_tag)

        cands = []
        cL = pick_side(VL, angL, idxL, xL, yL, uzL, 1)
        if cL is not None: cands.append(cL)
        cR = pick_side(VR, angR, idxR, xR, yR, uzR, 0)
        if cR is not None: cands.append(cR)
        if not cands:
            continue
        dmin, pt, s_i, sid = min(cands, key=lambda t: t[0])
        # 去抖：相邻过近不重复
        if border_pts and np.hypot(pt[0]-border_pts[-1][0], pt[1]-border_pts[-1][1]) < 0.3:
            continue
        border_pts.append(pt)
        border_s.append(float(s_i))
        border_side.append(int(sid))
        border_ang.append(float(a))

    return {
        'points': [tuple(p) for p in border_pts],
        's_on_boundary': border_s,
        'side': border_side,
        'angles': border_ang,
    }

# ---------------- High-level helpers for downstream usage ----------------

def compute_fov_at_s0(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    w_left: np.ndarray,
    w_right: np.ndarray,
    s0: float,
    look: float = 200.0,
    **kwargs,
) -> Dict[str, object]:
    """便捷接口：一次性返回该 s0 的所有 FOV 相关信息。
    返回字典包含：
      - 'tangent_points'：全部切点
      - 'unknown_zones'：全部未知区
      - 'final_fpv'：{s_itp, s_fpv, side, itp_point, fpv_point} 或 None
      - 'fov_border'：可见域前沿。fov_border_mode='segment' 时为 ITP–FPV 直线段，'front' 时为角度扫描离散前沿。
    参数 fov_border_mode: str, 可选，'segment' 用 ITP–FPV 线段，'front' 用角度扫描前沿
    """
    fov_border_mode = kwargs.get('fov_border_mode', 'segment')
    res = compute_fov_full(
        spline_x=spline_x,
        spline_y=spline_y,
        kappa_spline=kwargs.get('kappa_spline', None),
        w_left=w_left,
        w_right=w_right,
        s0=s0,
        delta_s=look,
        s_dense_step=kwargs.get('s_dense_step', 0.5),
        theta_min_deg=kwargs.get('theta_min_deg', 4.0),
        peak_min_deg=kwargs.get('peak_min_deg', 0.0),
        min_dist=kwargs.get('min_dist', 5.0),
        min_arc=kwargs.get('min_arc', 10.0),
        kappa_min=kwargs.get('kappa_min', 1e-5),
        epsilon_dir=kwargs.get('epsilon_dir', 1e-3),
        fov_border_mode=fov_border_mode,
    )
    final_fpv = None
    if res.get('s_fpv') is not None and res.get('side') is not None:
        final_fpv = {
            's_itp': res.get('s_itp'),
            's_fpv': res.get('s_fpv'),
            'side': res.get('side'),
        }
    return {
        's0': s0,
        'tangent_points': res.get('tangent_points', []),
        'unknown_zones': res.get('unknown_zones', []),
        'final_fpv': final_fpv,
        'fov_border': res.get('fov_border', {}),
    }


def compute_fov_over_range(
    spline_x: CubicSpline,
    spline_y: CubicSpline,
    w_left: np.ndarray,
    w_right: np.ndarray,
    s0: float,
    search_range: float = 200.0,
    s0_step: float = 5.0,
    **kwargs,
) -> Dict[str, object]:
    """批量计算 [s0, s0+search_range) 区间内，以 s0_step 采样的 FOV 结果。
    返回：
      {
        's_list': [...],
        'final_fpv_s': [...],     # NaN 表示无 FPV
        'final_fpv_side': [...],  # -1 表示无 FPV
        'entries': [ {compute_fov_at_s0 的返回}, ... ]
      }
    参数 fov_border_mode: str, 可选。'segment'（默认）用 ITP–FPV 线段，'front' 用角度扫描离散前沿。
    """
    s_list = np.arange(s0, s0 + search_range, s0_step)
    final_fpv_s = []
    final_fpv_side = []
    entries = []
    for s in s_list:
        e = compute_fov_at_s0(
            spline_x, spline_y, w_left, w_right, s, look=kwargs.get('look', search_range), **kwargs
        )
        entries.append(e)
        fpv = e.get('final_fpv')
        if fpv is None:
            final_fpv_s.append(np.nan)
            final_fpv_side.append(-1)
        else:
            final_fpv_s.append(fpv['s_fpv'])
            final_fpv_side.append(fpv['side'])
    return {
        's_list': s_list.tolist(),
        'final_fpv_s': final_fpv_s,
        'final_fpv_side': final_fpv_side,
        'entries': entries,
    }