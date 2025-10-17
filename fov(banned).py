import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import casadi as ca
import pandas as pd





# 1. set file path
file_path = "../Laptimeoptimization/Monza.csv"

# 2. read csv file
df = pd.read_csv(file_path)

# 3. look head
print(df.head())

# 4. rename the first row
##df.rename(columns={'# x_m': 'x_m'}, inplace=True)

# 5. read the columns
x = df['# x_m'].values
y = df['y_m'].values
w_left = df['w_tr_left_m'].values
w_right = df['w_tr_right_m'].values








# Step 1: Compute arclengths
s = np.zeros(len(x))
s[1:] = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
s_max=s[-1] #the length of the track

# Step 2: Create cubic spline interpolators
x_spline = CubicSpline(s, x)
y_spline = CubicSpline(s, y)

# Step 3: Get the first and second derivative of the parameteriation
dx_ds = x_spline.derivative(1)
d2x_ds2 = x_spline.derivative(2)
dy_ds = y_spline.derivative(1)
d2y_ds2 = y_spline.derivative(2)

## Generate interpolated values for plotting
s_query = np.linspace(0, s_max, 500)
x_interp = x_spline(s_query)
y_interp = y_spline(s_query)

## Plot the spline
plt.figure(figsize=(8, 6))
plt.plot(x_interp, y_interp, label='Interpolated Track', color='blue')
##plt.plot(x, y, 'ro', label='Original Points')
plt.axis('equal')
plt.title('Race Track Centerline via Cubic Spline Interpolation')
plt.xlabel('x (meters)')
plt.ylabel('y (meters)')
plt.legend()
plt.grid(True)
plt.show()

# output parameterisation and curvature as a casadi instance
x_ref = ca.interpolant('x_ref', 'linear', [s_query], x_interp)
y_ref = ca.interpolant('y_ref', 'linear', [s_query], y_interp)
kappa_ref_ = dx_ds(s_query)*d2y_ds2(s_query)-dy_ds(s_query)*d2x_ds2(s_query)
kappa_ref = ca.interpolant('kappa', 'linear', [s_query], kappa_ref_)

















def compute_fov_simple(
        spline_x, spline_y,
        w_left, w_right,
        s_query,
        delta_s       = 200,   # search window length [m] – we look for tangent points
                           # in the range  s0 … s0 + delta_s.  Increase (→300-400 m)
                           # if you have very large-radius corners.

    dense_step    = 0.5,   # resampling step for track boundaries [m].
                           # Finer → smoother θ-curve but heavier computation.

    ang_min_deg   = 4,     # minimum total bearing span (θ.ptp) [deg] that qualifies
                           # as a “curved” segment.  If θ.ptp < ang_min_deg we treat
                           # the segment as straight and simply set sfov = s + delta_s.

    peak_min_deg  = 0,     # peak-gap threshold [deg] used to filter out
                           # ultra-flat “pseudo peaks”.  0 = no filtering;
                           # in production 0.5–2 deg is typical.

    min_dist      = 5,     # minimum distance from car to candidate point [m].
                           # Suppresses tangent points that are too close
                           # (numerical noise or inside track width).

    min_arc       = 10,    # minimum arc-length (s_tan − s0) [m] before a point
                           # can be accepted as tangent; prevents 1–2 m “false hits”.

    curv_min      = 1e-5   # curvature threshold – below this |κ| the centreline
                           # is considered straight, skipping tangent search.
):


    s_dense = np.arange(0, s_max + dense_step, dense_step)
    (xl, yl), (xr, yr) = build_boundaries(spline_x, spline_y,
                                          w_left, w_right, s_dense)

    sfov = np.empty_like(s_query)
    ang_min  = np.deg2rad(ang_min_deg)
    peak_min = np.deg2rad(peak_min_deg)

    # ---------------- 主循环 ----------------
    for k, s0 in enumerate(s_query):

        # 车辆位置 + 朝向
        p_car = np.array([spline_x(s0), spline_y(s0)])
        t     = np.array([spline_x.derivative(1)(s0),
                          spline_y.derivative(1)(s0)])
        dir_v = t / np.linalg.norm(t)

        # 直线段直接给最大 FOV
        kappa_now = (spline_x.derivative(1)(s0)*spline_y.derivative(2)(s0)
                     - spline_y.derivative(1)(s0)*spline_x.derivative(2)(s0))
        if abs(kappa_now) < curv_min:
            sfov[k] = s0 + delta_s
            continue

        idx0 = np.searchsorted(s_dense, s0)
        idx1 = np.searchsorted(s_dense, s0 + delta_s, side="right")
        s_win = s_dense[idx0:idx1]

        # ----------- 局部极值检测工具 -----------
        # ───── 把原来的 pick_extreme_local 整块替换为下面版本 ─────
        def pick_extreme_local(xb, yb, want_max):
            Vx = xb[idx0:idx1] - p_car[0]
            Vy = yb[idx0:idx1] - p_car[1]
            dist = np.hypot(Vx, Vy)
            dot  = dir_v[0]*Vx + dir_v[1]*Vy

            mask = (dist > min_dist) & (dot > 0)
            if not np.any(mask):
                return np.inf

            Vx, Vy, dot, s_seg = Vx[mask], Vy[mask], dot[mask], s_win[mask]
            theta = np.arctan2(dir_v[0]*Vy - dir_v[1]*Vx, dot)

            # ── 调试：始终打印一次 θ.ptp 及翻转数 ──
            if abs(s0 - dbg_s0) < dense_step/2:
                flips = np.where(np.sign(np.diff(theta))[:-1] *
                                np.sign(np.diff(theta))[1:] < 0)[0] + 1
                print(f"[dbg] s0={s0:.1f} θ.ptp={np.ptp(theta)*180/np.pi:5.2f}° "
                    f"flips={flips.size}")
            # ─────────────────────────────────────────

            # 1) 先找局部极值（符号翻转）
            dtheta = np.diff(theta)
            sign   = np.sign(dtheta)
            flips  = np.where(sign[:-1] * sign[1:] < 0)[0] + 1

            if flips.size:                       # A. 有翻转 → 用局部极值
                cand_idx = flips
            else:                                # B. 单调 → 用全局极值
                cand_idx = [np.argmax(theta)] if want_max else [np.argmin(theta)]

            # 2) 过滤尖峰
            peak_min = peak_min_deg * np.pi/180.0
            good = []
            for idx in cand_idx:
                left_gap  = abs(theta[idx] - theta[max(idx-1, 0)])
                right_gap = abs(theta[idx] - theta[min(idx+1, len(theta)-1)])
                if max(left_gap, right_gap) >= peak_min:
                    good.append(idx)
            if not good:
                return np.inf

            idx_sel = max(good, key=lambda i: theta[i]) if want_max else min(good, key=lambda i: theta[i])
            s_tan   = s_seg[idx_sel]

            # 3) 距离阈值
            if s_tan - s0 < min_arc:
                return np.inf

            # －－ 可选再打印一次最终筛选结果 －－
            if abs(s0 - dbg_s0) < dense_step/2:
                print(f"     arc={s_tan-s0:5.1f} m  gap={max(left_gap,right_gap)*180/np.pi:4.2f}°")

            return s_tan
        # ─────────────────────────────────────────────────────────────

        # 左右切点
        s_L = pick_extreme_local(xl, yl, want_max=True)   # 左边找最大 θ
        s_R = pick_extreme_local(xr, yr, want_max=False)  # 右边找最小 θ
        sfov[k] = min(s_L, s_R, s0 + delta_s)

    return sfov


# 1. 只接收一个返回值
sfov_vals = compute_fov_simple(
    x_spline, y_spline, w_left, w_right,
    s_query, delta_s=DELTA_S
)

# 2. CasADi 插值器
sfov_ca = ca.interpolant('sfov', 'linear', [s_query], sfov_vals)