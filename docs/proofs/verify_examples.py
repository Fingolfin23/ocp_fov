"""Independent analytic checks for the manuscript, not a test of the whole source pipeline.

Both roads have an injective Frenet map. The straight variable-width example is
an illustration of the geometry theorem, not a claim about a curvature-gated run.
Visibility oracles use polynomial minimization and exact circle distances;
predicted regions use the proposed boundary-arc/crosscut constructions.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon
from shapely.geometry import Point, Polygon

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260922)
N_TEST = 10000


def left_width(x):
    return 0.2*x + 0.01*(x-2)**2*(6-x)


def straight_visible_exact(x, y):
    # Minimize w(t) - (y/x)t on [0,x]. w is cubic; its derivative is quadratic.
    roots = np.roots([-0.03, 0.2, -0.08-y/x])
    candidates = [0., x]
    candidates.extend(float(z.real) for z in roots if abs(z.imag) < 1e-12 and 0 < z.real < x)
    margin = min(left_width(t) - (y/x)*t for t in candidates)
    return margin >= -1e-12 and y >= -1., margin


def circle_boundary(s, radius, R=10.):
    a = np.asarray(s)/R
    return np.column_stack([radius*np.sin(a), R-radius*np.cos(a)])


def draw_polygon(ax, geom, color, alpha=1.):
    if geom.is_empty:
        return
    if geom.geom_type == "Polygon":
        ax.add_patch(PlotPolygon(np.asarray(geom.exterior.coords), closed=True,
                                 facecolor=color, edgecolor="none", alpha=alpha))
    elif hasattr(geom, "geoms"):
        for p in geom.geoms:
            draw_polygon(ax, p, color, alpha)


def main():
    xs = np.linspace(0., 8., 4001)
    road1 = Polygon(np.vstack([np.column_stack([xs, left_width(xs)]), [[8., -1.], [0., -1.]]]))
    arc = np.linspace(2., 6., 4001)
    pocket = Polygon(np.column_stack([arc, left_width(arc)]))
    T1, U1 = np.array([2., .4]), np.array([6., 1.2])
    errors1, skipped1 = [], 0
    for _ in range(N_TEST):
        x = RNG.uniform(1e-4, 8.)
        y = RNG.uniform(-1., left_width(x))
        exact, margin = straight_visible_exact(x, y)
        if abs(margin) < 5e-6:
            skipped1 += 1
            continue
        predicted = not pocket.contains(Point(x, y))
        if predicted != exact:
            errors1.append([x, y, exact, predicted])

    R, w, H = 10., 1., 15.
    ri, ro = R-w, R+w
    alpha = np.arccos(ri/R)
    beta = alpha + np.arccos(ri/ro)
    st, su = R*alpha, R*beta
    T2 = circle_boundary([st], ri)[0]
    U2 = circle_boundary([su], ro)[0]
    road2 = Polygon(np.vstack([circle_boundary(np.linspace(0., H, 6001), ri),
                               circle_boundary(np.linspace(H, 0., 6001), ro)]))
    near = Polygon(np.vstack([circle_boundary(np.linspace(0., st, 4001), ri),
                              circle_boundary(np.linspace(su, 0., 4001), ro)]))
    errors2, skipped2 = [], 0
    O = np.array([0., R])
    for _ in range(N_TEST):
        s = RNG.uniform(1e-5, H)
        radius = RNG.uniform(ri, ro)
        q = circle_boundary([s], radius)[0]
        lam = np.clip(np.dot(O, q)/np.dot(q, q), 0., 1.)
        clearance = np.linalg.norm(lam*q-O) - ri
        if abs(clearance) < 5e-6:
            skipped2 += 1
            continue
        exact = bool(clearance >= 0.)
        predicted = near.covers(Point(*q))
        if predicted != exact:
            errors2.append([s, radius, exact, predicted])

    d = T2/np.linalg.norm(T2)
    self_hit = []
    for step in [.1, .01, .001, .0001]:
        vec = circle_boundary([st+step], ri)[0]-T2
        ratio = abs(d[0]*vec[1]-d[1]*vec[0])/np.dot(d, vec)
        self_hit.append({"station_step": step, "alignment_ratio": float(ratio),
                         "accepted_by_fixed_epsilon_1e-3": bool(ratio < 1e-3),
                         "is_exact_second_intersection": False})

    results = {
        "purpose": "Analytic sanity checks; not a proof by testing and not the complete report pipeline.",
        "seed": 20260922,
        "same_side_pocket": {
            "centerline": "c(s)=(s,0), 0<=s<=8",
            "left_width": "0.2*s+0.01*(s-2)^2*(6-s)",
            "right_width": 1., "T": T1.tolist(), "U": U1.tolist(),
            "oracle": "Exact cubic minimum along each sight segment",
            "samples": N_TEST, "skipped_near_boundary": skipped1,
            "mismatches": len(errors1), "first_mismatches": errors1[:5]},
        "opposite_side_front": {
            "center_radius": R, "half_width": w, "window": H,
            "s_T": float(st), "s_U": float(su), "T": T2.tolist(), "U": U2.tolist(),
            "oracle": "Exact minimum distance of sight segment to inner-circle center",
            "samples": N_TEST, "skipped_near_boundary": skipped2,
            "mismatches": len(errors2), "first_mismatches": errors2[:5]},
        "fixed_tolerance_near_tangent": self_hit}
    (OUT/"verification_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), constrained_layout=True)
    green, red = "#b7d9c5", "#edaaa3"
    for ax, road, unknown, T, U, title in [
            (axes[0], road1, pocket, T1, U1, "Same-side return: a local blind pocket"),
            (axes[1], road2, road2.difference(near), T2, U2, "Opposite-side return: a global front")]:
        draw_polygon(ax, road, green)
        draw_polygon(ax, unknown, red)
        coords = np.asarray(road.exterior.coords)
        ax.plot(coords[:, 0], coords[:, 1], color="#475569", lw=1.1)
        ax.plot([0, T[0]], [0, T[1]], "--", color="#62748a", lw=1.4)
        ax.plot([T[0], U[0]], [T[1], U[1]], color="#1c4568", lw=2.4)
        ax.scatter([0], [0], marker="*", s=100, color="#172f4d", zorder=5)
        ax.scatter([T[0], U[0]], [T[1], U[1]], s=45, color=["#ae640d", "#245d46"], zorder=5)
        for point, label, offset in [(np.zeros(2), "P", (-14, -15)), (T, "T", (5, -15)), (U, "U", (5, 7))]:
            ax.annotate(label, point, xytext=offset, textcoords="offset points", fontsize=12, fontweight="bold")
        ax.set_title(title, fontsize=12, pad=12)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.grid(alpha=.12)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylim(-1.7, 2.2)
    axes[1].set_ylim(-1.8, 10.)
    fig.suptitle("Analytic examples for the two-stage visibility construction", fontsize=15)
    fig.savefig(OUT/"proof_examples.png", dpi=190, facecolor="white")
    plt.close(fig)
    assert not errors1 and not errors2, results
    print(json.dumps({"pocket_checked": N_TEST-skipped1, "front_checked": N_TEST-skipped2,
                      "mismatches": len(errors1)+len(errors2), "s_T": st, "s_FPV": su}, indent=2))


if __name__ == "__main__":
    main()
