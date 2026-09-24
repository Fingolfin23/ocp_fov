"""Periodic Monza replay geometry, with the original sampled FOV algorithm.

The algorithm markers come directly from ``fov_extended.compute_fov_full``.
The blue visible-area reference is computed independently by first ray/segment
intersection with a sampled forward road window.  It is a rendering aid, not a
replacement for the algorithm's station comparisons or a certificate of their
completeness.  A missing algorithm front does not imply an entirely visible
window.  Core boundary names are preserved internally; because their historical
left/right convention differs from geometric left/right, use neutral A/B names
in a presentation.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import CubicSpline


def _core_module(track_path: Path, repo_path: str | Path | None = None):
    candidates = [Path(repo_path)] if repo_path else []
    candidates += [track_path.parent, Path(__file__).resolve().parents[1]]
    for root in candidates:
        path = root / "fov_extended.py"
        if path.is_file():
            spec = importlib.util.spec_from_file_location("monza_original_fov", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError("Cannot locate fov_extended.py; pass repo_path.")


def _area(points: np.ndarray) -> float:
    if len(points) < 3:
        return 0.0
    # Shift before applying shoelace to reduce cancellation at world coordinates.
    q = points - points[0]
    return 0.5 * abs(float(np.sum(q[:, 0] * np.roll(q[:, 1], -1)
                                  - q[:, 1] * np.roll(q[:, 0], -1))))


def _visible_polygon(p: np.ndarray, heading: float,
                     window: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    """First exits from P in the forward half-plane of a closed polyline.

    Boundary-vertex directions and their two one-sided limits capture narrow
    occlusion wedges.  Uniform directions improve visual smoothness.  Zero-range
    contact at the artificial starting cap is ignored; P is on that cap.
    """
    rel = window - p
    angles = np.arctan2(rel[:, 1], rel[:, 0]) - heading
    angles = np.arctan2(np.sin(angles), np.cos(angles))
    eps = 1e-7
    limit = np.pi / 2 - eps
    angles = angles[(angles >= -limit) & (angles <= limit)]
    angles = np.unique(np.round(np.concatenate([
        np.linspace(-limit, limit, 181), angles - eps, angles, angles + eps,
    ]), 12))
    angles = angles[(angles >= -limit) & (angles <= limit)]
    directions = np.column_stack((np.cos(angles + heading), np.sin(angles + heading)))
    a = rel
    e = np.roll(window, -1, axis=0) - window
    numerator = a[:, 0] * e[:, 1] - a[:, 1] * e[:, 0]
    # Chunk the query to keep temporary arrays small on typical laptops.
    ranges = []
    for start in range(0, len(directions), 256):
        d = directions[start:start + 256]
        denominator = d[:, None, 0] * e[None, :, 1] - d[:, None, 1] * e[None, :, 0]
        valid_denominator = abs(denominator) > 1e-12
        denominator = np.where(valid_denominator, denominator, np.nan)
        t = numerator[None, :] / denominator
        u = (a[None, :, 0] * d[:, None, 1] - a[None, :, 1] * d[:, None, 0]) / denominator
        valid = valid_denominator & (t > 1e-7) & (u >= -1e-9) & (u <= 1 + 1e-9)
        ranges.append(np.min(np.where(valid, t, np.inf), axis=1))
    distance = np.concatenate(ranges)
    finite = np.isfinite(distance)
    points = p + directions[finite] * distance[finite, None]
    return np.vstack((p, points, p)), {
        "rays": len(directions), "rays_without_exit": int(np.sum(~finite)),
    }


def _los_crosscheck(frames: list[dict[str, Any]], samples: int = 60) -> dict[str, Any]:
    """Independent segment-in-polygon checks at deterministic interior samples.

    Shapely checks whether the full P--Q segment stays inside the sampled road
    window.  The result is compared with membership in the ray-cast polygon.
    Points within 2 cm of either boundary are excluded from this floating-point
    illustration check; this is not a proof for the continuous road.
    """
    try:
        from shapely.geometry import LineString, Point, Polygon
    except ImportError:
        return {"status": "skipped", "reason": "optional shapely dependency unavailable"}
    rng = np.random.default_rng(20260924)
    tested = mismatches = invalid = 0
    selected = np.unique(np.linspace(0, len(frames) - 1, min(12, len(frames)), dtype=int))
    for idx in selected:
        frame = frames[int(idx)]
        window = Polygon(frame["window"])
        visible = Polygon(frame["visible"])
        if not window.is_valid or not visible.is_valid:
            invalid += 1
            continue
        near, far = frame["forward_left"], frame["forward_right"]
        heading = np.array([math.cos(frame["heading"]), math.sin(frame["heading"])])
        for _ in range(samples):
            station_idx = int(rng.integers(1, len(near) - 1))
            lam = float(rng.uniform(0.05, 0.95))
            q = lam * near[station_idx] + (1 - lam) * far[station_idx]
            point = Point(q)
            if not window.contains(point):
                continue
            if window.boundary.distance(point) < 0.02 or visible.boundary.distance(point) < 0.02:
                continue
            actual = (float((q - frame["p"]) @ heading) > 0
                      and window.buffer(1e-7).covers(LineString([frame["p"], q])))
            drawn = visible.covers(point)
            tested += 1
            mismatches += bool(actual != drawn)
    return {"status": "checked", "frames": len(selected), "tested_points": tested,
            "mismatches": mismatches, "invalid_sampled_polygons": invalid,
            "boundary_exclusion_m": 0.02, "random_seed": 20260924}


def _front_reference_check(frames: list[dict[str, Any]], tolerance: float = 0.25) -> dict[str, Any]:
    """Annotate sampled output against the independent drawing reference.

    This diagnostic checks the visible prefix, road-contained TP--FPV segment,
    and FPV membership in the ray-cast polygon.  It does not alter the core
    result, prove maximality, or establish correctness for the continuous road.
    None means no algorithm front or an unavailable optional diagnostic.
    """
    for frame in frames:
        frame["front_reference_agrees"] = None
    audit: dict[str, Any] = {
        "reference_front_check_tolerance_m": tolerance,
        "reference_front_check_scope": "Sampled-polygon diagnostic only; not a proof for the continuous road or of FPV maximality.",
    }
    try:
        from shapely.geometry import LineString, Point, Polygon
    except ImportError:
        return {**audit, "reference_front_check_status": "skipped",
                "reference_front_check_reason": "optional shapely dependency unavailable",
                "reference_front_checked_frames": None,
                "reference_front_mismatch_frames": None,
                "reference_front_mismatch_ids": None}
    mismatch_ids = []
    checked = 0
    for index, frame in enumerate(frames):
        if frame["fpv_delta"] is None:
            continue
        checked += 1
        tp, fpv = frame["border"]
        road = Polygon(frame["window"]).buffer(tolerance)
        visible = Polygon(frame["visible"]).buffer(tolerance)
        agrees = bool(road.covers(LineString([frame["p"], tp]))
                      and road.covers(LineString([tp, fpv]))
                      and visible.covers(Point(fpv)))
        frame["front_reference_agrees"] = agrees
        if not agrees:
            mismatch_ids.append(index)
    return {**audit, "reference_front_check_status": "checked",
            "reference_front_checked_frames": checked,
            "reference_front_mismatch_frames": len(mismatch_ids),
            "reference_front_mismatch_ids": mismatch_ids}


def prepare_replay(track_path: str | Path, frames: int = 720, look: float = 200.0,
                   step: float = 0.5, repo_path: str | Path | None = None) -> dict[str, Any]:
    """Compute an evenly spaced full lap without changing the FOV algorithm.

    The final ~5 m source gap is closed by appending the first sample.  Position,
    tangent, and widths repeat periodically, including windows across the finish
    line.  The inherited core boundary convention is used for every overlay.

    Return NumPy coordinate arrays and one frame per station in [0, L).
    ``window`` and ``visible`` use a separate <= 1 m boundary polyline.  ``pockets``
    are boundary-arc/chord candidates reconstructed from core unknown-zone
    station intervals, never the old nearly collinear P--TP--hit triangles.
    """
    if frames < 2 or look <= 0 or step <= 0:
        raise ValueError("frames >= 2, look > 0 and step > 0 are required")
    track_path = Path(track_path).resolve()
    core = _core_module(track_path, repo_path)
    data = np.loadtxt(track_path, delimiter=",")
    if data.ndim != 2 or data.shape[1] < 4 or not np.all(np.isfinite(data)):
        raise ValueError("Track must have finite x, y, right-width, left-width columns")
    delta = np.linalg.norm(np.diff(data[:, :2], axis=0), axis=1)
    if np.any(delta <= 0) or np.any(data[:, 2:4] <= 0):
        raise ValueError("Consecutive track points and positive widths are required")
    original_s = np.r_[0.0, np.cumsum(delta)]
    closing_gap = float(np.linalg.norm(data[-1, :2] - data[0, :2]))
    if closing_gap > 1e-9:
        closed = np.vstack((data, data[0]))
        station = np.r_[original_s, original_s[-1] + closing_gap]
    else:
        closed = data.copy()
        closed[-1] = closed[0]
        station = original_s
    length = float(station[-1])
    if look + 240 >= length:
        raise ValueError("Replay window is intended to cover much less than one lap")
    # Extending the knots is essential: the core interpolates widths with
    # np.interp(s, spline.x, widths), which does not itself wrap periodically.
    knots = np.r_[station[:-1] - length, station, station[1:] + length]
    values = np.concatenate((closed[:-1], closed, closed[1:]))
    x = CubicSpline(knots, values[:, 0], bc_type="periodic")
    y = CubicSpline(knots, values[:, 1], bc_type="periodic")
    w_left, w_right = values[:, 3], values[:, 2]

    def geometry(s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        left, right = core.build_boundaries(
            x, y, np.interp(s, knots, w_left), np.interp(s, knots, w_right), s)
        return np.column_stack((x(s), y(s))), np.column_stack(left), np.column_stack(right)

    whole_s = np.linspace(0, length, int(np.ceil(length / 2.0)) + 1)
    whole_center, whole_left, whole_right = geometry(whole_s)
    result_frames: list[dict[str, Any]] = []
    ray_misses = ray_count = 0
    render_step = min(1.0, max(step, 0.25))
    # At default step=0.5 the drawing and core boundary polyline coincide.
    render_offsets = np.linspace(0.0, look, int(np.ceil(look / render_step)) + 1)
    for s0 in np.linspace(0.0, length, frames, endpoint=False):
        result = core.compute_fov_full(
            spline_x=x, spline_y=y, kappa_spline=None,
            w_left=w_left, w_right=w_right, s0=float(s0), delta_s=look,
            s_dense_step=step, theta_min_deg=4.0,
            fov_border_mode="segment")
        forward_s = s0 + render_offsets
        forward_center, forward_left, forward_right = geometry(forward_s)
        context_s = np.arange(s0 - 90.0, s0 + 240.0 + 1.0, 1.0)
        context_center, context_left, context_right = geometry(context_s)
        p = np.array([x(s0), y(s0)], dtype=float)
        heading = float(np.arctan2(y.derivative()(s0), x.derivative()(s0)))
        window = np.concatenate((forward_left, forward_right[::-1]))
        visible, ray_audit = _visible_polygon(p, heading, window)
        ray_count += ray_audit["rays"]
        ray_misses += ray_audit["rays_without_exit"]
        pockets = []
        for zone in result["unknown_zones"]:
            zs = np.linspace(zone["s_start"], zone["s_end"],
                             max(3, int(np.ceil((zone["s_end"] - zone["s_start"]) / step)) + 1))
            _, zl, zr = geometry(zs)
            arc = zl if zone["side"] == 1 else zr
            poly = np.vstack((arc, arc[0]))
            if _area(poly) > 1e-4:
                pockets.append(poly)
        border = np.asarray(result.get("fov_border", {}).get("points", []), dtype=float).reshape(-1, 2)
        tangent_meta = result.get("tangent_points", [])
        tangents = np.asarray([tp["point"] for tp in tangent_meta], dtype=float).reshape(-1, 2)
        result_frames.append({
            "s0": float(s0), "p": p, "heading": heading,
            "forward_s": forward_s, "forward_center": forward_center,
            "forward_left": forward_left, "forward_right": forward_right,
            "context_center": context_center, "context_left": context_left,
            "context_right": context_right, "tangent_points": tangents,
            "tangent_meta": tangent_meta, "border": border,
            "fpv_delta": None if result["s_fpv"] is None else float(result["s_fpv"] - s0),
            "s_itp": result["s_itp"], "s_fpv": result["s_fpv"],
            "window": window, "visible": visible, "pockets": pockets,
            "unknown_zone_count": len(result["unknown_zones"]),
        })
    finite = all(np.all(np.isfinite(frame[key]))
                 for frame in result_frames
                 for key in ["p", "forward_left", "forward_right", "context_center",
                             "context_left", "context_right", "border", "tangent_points",
                             "window", "visible"])
    fpv_delta = [f["fpv_delta"] for f in result_frames if f["fpv_delta"] is not None]
    front_reference_audit = _front_reference_check(result_frames)
    audit = {
        "source_samples": int(len(data)), "source_open_length_m": float(original_s[-1]),
        "source_closure_gap_m": closing_gap, "closed_lap_length_m": length,
        "seam_position_error_m": float(np.hypot(x(length) - x(0), y(length) - y(0))),
        "seam_tangent_error": float(np.hypot(x.derivative()(length) - x.derivative()(0),
                                              y.derivative()(length) - y.derivative()(0))),
        "frame_count": int(frames), "station_step_m": length / frames,
        "core_look_m": look, "core_sample_step_m": step,
        "core_angle_range_threshold_deg": 4.0,
        "core_parameters": {"peak_min_deg": 0.0, "min_dist": 5.0, "min_arc": 10.0,
                            "kappa_min": 1e-5, "epsilon_dir": 1e-3,
                            "fov_border_mode": "segment", "kappa_spline": None},
        "reference_boundary_sample_step_m": look / (len(render_offsets) - 1),
        "front_frames": len(fpv_delta), "no_front_frames": frames - len(fpv_delta),
        "min_detected_front_distance_m": min(fpv_delta) if fpv_delta else None,
        "max_detected_front_distance_m": max(fpv_delta) if fpv_delta else None,
        "frames_with_unknown_zone_candidates": sum(bool(f["unknown_zone_count"]) for f in result_frames),
        "reference_ray_count": ray_count, "reference_rays_without_exit": ray_misses,
        "finite_coordinates": bool(finite), "line_of_sight_crosscheck": _los_crosscheck(result_frames),
        "scope": "Independent sampled-polygon forward visibility is a rendering reference; core FOV results are unchanged. Missing front is not full visibility.",
        "boundary_convention": "Original fov_extended.build_boundaries convention; use neutral A/B labels.",
        **front_reference_audit,
    }
    return {"length": length, "whole_s": whole_s, "whole_center": whole_center,
            "whole_left": whole_left, "whole_right": whole_right,
            "frames": result_frames, "audit": audit}
