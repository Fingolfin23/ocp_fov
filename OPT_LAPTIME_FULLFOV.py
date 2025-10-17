import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import pandas as pd
import json
from fov_to_frenet_bounds import build_planning_bounds_arrays

# ---------------- FOV helpers: FPV mapping & safe speed profile ----------------

def fpv_s_on_grid(cache_path: str, s_grid: np.ndarray, default_look: float = 200.0) -> np.ndarray:
    """Map each s_k to farthest-visible arclength s_fpv[k] using the FOV cache.
    If FPV is missing/behind, fall back to sk + default_look (meters).
    """
    try:
        with open(cache_path, 'r', encoding='utf-8') as f:
            cache = json.load(f)
    except Exception:
        return s_grid + float(default_look)
    entries = cache.get('entries', [])
    if not entries:
        return s_grid + float(default_look)
    s0  = np.array([float(e.get('s0', 0.0)) for e in entries], dtype=float)
    sfv = np.array([ (e.get('final_fpv', {}) or {}).get('s', np.nan) for e in entries ], dtype=float)
    s_fpv = np.empty_like(s_grid, dtype=float)
    for k, sk in enumerate(s_grid):
        j = int(np.argmin(np.abs(s0 - sk)))
        val = sfv[j]
        if np.isfinite(val) and (val >= sk):
            s_fpv[k] = val
        else:
            s_fpv[k] = sk + float(default_look)
    return s_fpv


def build_vsafe_array(mode: str,
                      s_grid: np.ndarray,
                      kappa: np.ndarray,
                      a_y_lim: float = 8.0,
                      v_const: float = 0.0,
                      s_fpv: np.ndarray | None = None) -> np.ndarray:
    """Construct per-node safe speed v_safe[k].
    mode:
      - 'constant': use v_const everywhere.
      - 'curvature': v_safe = sqrt(a_y_lim / max(|kappa_at_fpv|, eps)). If s_fpv is None, use kappa at s_k.
    """
    eps = 1e-6
    if mode == 'constant':
        return np.full_like(s_grid, float(v_const), dtype=float)
    if mode == 'curvature':
        if s_fpv is None:
            k_use = np.abs(kappa)
        else:
            idx = np.clip(np.searchsorted(s_grid, s_fpv), 0, len(s_grid)-1)
            k_use = np.abs(kappa[idx])
        v = np.sqrt(np.maximum(a_y_lim, 0.0) / np.maximum(k_use, eps))
        return v
    return np.full_like(s_grid, float(v_const), dtype=float)

# 车辆模型
def vehicle_dynamics(x, u, kappa):
    v    = x[0]
    beta = x[1]
    omega= x[2]
    n    = x[3]
    xi   = x[4]
    delta= u[0]
    F_dr = u[1]
    # 参数
    m = 1200.0
    Jzz = 1260.0
    lf, lr = 1.5, 1.4
    wf, wr_ = 1.6, 1.5
    g = 9.81
    Fz = m * g / 4
    delta_r = 0.0
    epsilon = 1e-4
    mu_roll = 0.015
    rho = 1.2041
    Cd = 0.3
    Af = 1


# 轮胎模型
    def slip_angle(delta, beta, l, w, omega, is_front):
        dy = v * ca.sin(beta) + l * omega if is_front else v * ca.sin(beta) - l * omega
        dx_left = v * ca.cos(beta) - w / 2 * omega
        dx_right = v * ca.cos(beta) + w / 2 * omega
        return delta - ca.atan2(dy, dx_left), delta - ca.atan2(dy, dx_right)
    alpha_fl, alpha_fr = slip_angle(delta, beta, lf, wf, omega, True)
    alpha_rl, alpha_rr = slip_angle(delta_r, beta, lr, wr_, omega, False)

    def Fy(alpha):
        #return ca.sin(1.9 * ca.atan(10 * ca.sin(alpha))) + 0.97 * ca.sin(ca.atan(10 * ca.sin(alpha)))
        return 6000 * alpha  # 线性化模型

    FY_fl = Fz * Fy(alpha_fl)
    FY_fr = Fz * Fy(alpha_fr)
    FY_rl = Fz * Fy(alpha_rl)
    FY_rr = Fz * Fy(alpha_rr)



    SF = (1 - n * kappa) / (v * ca.cos(xi + beta) + epsilon)
    F_roll = mu_roll * m * g
    F_aero = 0.5 * rho * Cd * Af * v**2

    dot_n = SF * v * ca.sin(xi + beta)
    dot_xi = SF * omega - kappa
    dot_beta = SF * (-omega + 1 / (m * v + epsilon) *
                    ((FY_fl + FY_fr) * ca.cos(delta - beta) +
                     (FY_rl + FY_rr) * ca.cos(delta_r - beta) -
                     F_dr * ca.sin(beta)))
    dot_omega = SF / Jzz * (
        FY_fl * (lf * ca.cos(delta) - wf / 2 * ca.sin(delta)) +
        FY_fr * (lf * ca.cos(delta) + wf / 2 * ca.sin(delta)) +
        FY_rl * (-lr * ca.cos(delta_r) - wr_ / 2 * ca.sin(delta_r)) +
        FY_rr * (-lr * ca.cos(delta_r) + wr_ / 2 * ca.sin(delta_r))
    )
    dot_v = SF / m * (
        (FY_fl + FY_fr) * ca.sin(beta - delta) +
        (FY_rl + FY_rr) * ca.sin(beta - delta_r) +
        F_dr * ca.cos(beta) - F_roll - F_aero
    )

    dx = ca.vertcat(dot_v, dot_beta, dot_omega, dot_n, dot_xi)
    return dx, SF, FY_fl, FY_fr, FY_rl, FY_rr

# CSV采样并等间距重采样
def generate_reference_path(csv_file, ds=0.5, max_points=None):
    data = np.loadtxt(csv_file, delimiter=',')
    # if max_points is not None:
    #     data = data[:max_points]
    #  在这里设置采样路段
    data = data[140:240]
    x_ref = data[:, 0]
    y_ref = data[:, 1]
    wr = data[:, 2]
    wl = data[:, 3]

    dx = np.diff(x_ref)
    dy = np.diff(y_ref)
    ds_raw = np.hypot(dx, dy)
    s_raw = np.insert(np.cumsum(ds_raw), 0, 0)
    s_max = s_raw[-1]

    s_uniform = np.arange(0, s_max, ds)
    if s_uniform[-1] < s_max:
        s_uniform = np.append(s_uniform, s_max)

    fx = interp1d(s_raw, x_ref, kind='cubic')
    fy = interp1d(s_raw, y_ref, kind='cubic')
    fwr = interp1d(s_raw, wr, kind='linear')
    fwl = interp1d(s_raw, wl, kind='linear')

    x_ref_new = fx(s_uniform)
    y_ref_new = fy(s_uniform)
    wr_new = fwr(s_uniform)
    wl_new = fwl(s_uniform)

    dx_ds = np.gradient(x_ref_new, s_uniform)
    dy_ds = np.gradient(y_ref_new, s_uniform)
    theta = np.arctan2(dy_ds, dx_ds)
    dtheta_ds = np.gradient(theta, s_uniform)
    kappa = dtheta_ds

    return s_uniform, x_ref_new, y_ref_new, theta, kappa, wr_new, wl_new

# 优化主程序
def solve_min_time(csv_file='Nuerburgring.csv', ds=2,
                   fov_cache_path='fov_cache.json',
                   blend_mode='pairwise_intersect',
                   include_unknown=True,
                   a_brake=10.0,
                   vsafe_mode='curvature',
                   vsafe_const=0.0,
                   ay_limit=8.0,
                   default_look=200.0):
    # === 生成等间距参考路径 ===
    s, x_ref, y_ref, theta, kappa, wr, wl = generate_reference_path(csv_file, ds=ds, max_points=1160)

    # === FOV→Frenet planning band (road ∩ FOV front ∩ Unknown) ===
    try:
        nMinPlan, nMaxPlan, _aux = build_planning_bounds_arrays(
            s_grid=s,
            x_ref=x_ref,
            y_ref=y_ref,
            theta=theta,
            wr=wr,
            wl=wl,
            cache_json_path=fov_cache_path,
            blend_mode=blend_mode,
            include_unknown=include_unknown,
        )
    except Exception:
        # fallback to pure road bounds if FOV cache not available
        nMinPlan = -wr.copy()
        nMaxPlan =  wl.copy()

    # === FPV map & per-node safe speed ===
    s_fpv = fpv_s_on_grid(fov_cache_path, s, default_look=default_look)
    v_safe_arr = build_vsafe_array(vsafe_mode, s, kappa,
                                   a_y_lim=ay_limit,
                                   v_const=vsafe_const,
                                   s_fpv=s_fpv)

    N = len(s) - 1

    # ---- 优化参数 ----
    lambda_ddelta = 10 /3 # 3.4
    lambda_dF = 0.0001 /3 #  600120
    lambda_dn = 1.0 /3 # 21
    lambda_dxi = 5 /3 #  0.4
    lambda_dv = 0.1/3 # 0.8
    lambda_jerk = 1.0 /9 # 6.38

    m = 1500.0
    g = 9.81
    F_z = m * g / 4
    mu = 0.9
    F_x_max = 7100.0 / 2
    F_y_max = mu * F_z
    a_x_max = mu * g
    a_y_max = mu * g

    nx, nu = 5, 2
    x0_val = [20.0, 0.0, 0.0, 0.0, 0.0]

    w, w0, lbw, ubw, g_list, lbg, ubg = [], [], [], [], [], [], []
    J = 0
    Xk = ca.MX.sym("X0", nx)
    w += [Xk]
    w0 += x0_val
    lbw += x0_val
    ubw += x0_val
    delta_prev = ca.MX(0)
    delta_prev_prev = ca.MX(0)
    F_prev = ca.MX(0)
    Xk_prev = Xk
    J_time = 0
    for k in range(N):
        Uk = ca.MX.sym(f"U_{k}", nu)
        w += [Uk]
        w0 += [0.0, 500.0]
        lbw += [-np.pi/6, -21000.0]
        ubw += [np.pi/6, 7100.0]

        kappa_k = float(kappa[k])
        Xk_next = ca.MX.sym(f"X_{k + 1}", nx)
        w += [Xk_next]
        w0 += [20.0, 0.0, 0.0, 0.0, 0.0]

        # ==== 绿色规划带（道路∩FOV∩Unknown）横向约束 ====
        idx_n = min(k+1, len(s)-1)  # bounds apply to Xk_next
        n_min_k = float(nMinPlan[idx_n])
        n_max_k = float(nMaxPlan[idx_n])
        lbw += [0, -ca.inf, -2, n_min_k, -ca.inf]
        ubw += [42.5, ca.inf,  2, n_max_k,  ca.inf]

        f_k, SF_k, FY_fl, FY_fr, FY_rl, FY_rr = vehicle_dynamics(Xk, Uk, kappa_k)
        g_list += [Xk + ds * f_k - Xk_next]
        lbg += [0.0] * nx
        ubg += [0.0] * nx

        # === 可见距离内可刹到安全速度： 2*a_brake*(s_fpv - s_k) + v_safe^2 - v_k^2 >= 0 ===
        g_list += [ 2.0*a_brake * ( s_fpv[k] - s[k] ) + (v_safe_arr[k]**2) - (Xk[0]**2) ]
        lbg    += [ 0.0 ]
        ubg    += [ ca.inf ]

        J += ds * SF_k
        J_time += ds * SF_k
        # ====================================
        # Regularization Term
        # ====================================
        # if k > 0:
        #     J += lambda_ddelta * (Uk[0] - delta_prev) ** 2
        #     # J += lambda_dF * (Uk[1] - F_prev) ** 2
        #     J += lambda_dn * (Xk[3] - Xk_prev[3]) ** 2
        #     J += lambda_dxi * (Xk[4] - Xk_prev[4]) ** 2
        #     J += lambda_dv * (Xk[0] - Xk_prev[0]) ** 2
        # if k > 1:
        #     J += lambda_jerk * (Uk[0] - 2 * delta_prev + delta_prev_prev) ** 2
        delta_prev_prev = delta_prev
        delta_prev = Uk[0]
        F_prev = Uk[1]
        Xk_prev = Xk

        # ----------- 轮胎椭圆限制 -----------
        F_dr = Uk[1]
        F_x_rl = F_dr / 2
        F_x_rr = F_dr / 2
        ellipse_rl = (F_x_rl / F_x_max) ** 2 + (FY_rl / F_y_max) ** 2
        ellipse_rr = (F_x_rr / F_x_max) ** 2 + (FY_rr / F_y_max) ** 2
        g_list += [ellipse_rl, ellipse_rr]
        lbg += [-ca.inf, -ca.inf]
        ubg += [1.0, 1.0]

        Xk = Xk_next

    nlp = {"f": J, "x": ca.vertcat(*w), "g": ca.vertcat(*g_list)}
    # 在这里设置求解器迭代次数
    solver = ca.nlpsol("solver", "ipopt", nlp, {"ipopt.print_level": 5, "print_time": False,"ipopt.max_iter": 3000})
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    # 在solve_min_time末尾（求解后）
    J_time_func = ca.Function("J_time", [ca.vertcat(*w)], [J_time])
    total_time = float(J_time_func(sol['x']))

    w_opt = sol['x'].full().flatten()
    X_opt, U_opt = [], []
    for k in range(N):
        base = k * (nx + nu)
        X_opt.append(w_opt[base: base + nx])
        U_opt.append(w_opt[base + nx: base + nx + nu])
    X_opt.append(w_opt[-nx:])
    X_opt = np.array(X_opt)
    U_opt = np.array(U_opt)

    n_opt = X_opt[:, 3]
    x_opt = x_ref - n_opt * np.sin(theta)
    y_opt = y_ref + n_opt * np.cos(theta)

    # === 两个子图 ===
    plt.figure(figsize=(12, 5))

    # 1. velocity profile
    plt.subplot(1, 2, 1)
    plt.plot(s, X_opt[:, 0],color='#2F5597', label='velocity (m/s)')
    plt.xlabel('s (m)')
    plt.ylabel('velocity (m/s)')
    plt.grid(True)
    plt.title("Velocity Profile")
    plt.legend()
    # 2. trajectory
    plt.subplot(1, 2, 2)
    # Left edge: r + wL * n_hat (n_hat = [-sin theta, cos theta])
    x_left  = x_ref - wl * np.sin(theta)
    y_left  = y_ref + wl * np.cos(theta)
    # Right edge: r - wR * n_hat
    x_right = x_ref + wr * np.sin(theta)
    y_right = y_ref - wr * np.cos(theta)
    plt.plot(x_ref, y_ref, 'k--', label='Reference Path')  # 黑色虚线
    plt.plot(x_left, y_left, color='black', linestyle='-')  # 8B0000
    plt.plot(x_right, y_right, color='black', linestyle='-')  # 2F5597
    plt.plot(x_opt, y_opt, color='red', label='Optimized Path')  # 红色实线
    plt.axis('equal')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title("Trajectory")
    plt.text(0.05, 0.95, f'Total Time: {total_time:.2f} s',
             transform=plt.gca().transAxes,
             fontsize=12, color='blue',
             verticalalignment='top')
    plt.grid(True)
    plt.legend()
    # s_list = s[:len(X_opt)]
    # df_data = {
    #     's': s_list,
    #     'v': X_opt[:, 0],
    #     'beta': X_opt[:, 1],
    #     'omega': X_opt[:, 2],
    #     'n': X_opt[:, 3],
    #     'xi': X_opt[:, 4],
    #     'delta': np.concatenate((U_opt[:, 0],[0])),
    #     'F_dr': np.concatenate((U_opt[:, 1],[0])),
    # }
    # df = pd.DataFrame(df_data)
    # df.to_csv("optimal_trajectory.csv", index=False)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    solve_min_time('Monza.csv', ds=3,
                   fov_cache_path='fov_cache.json',
                   blend_mode='pairwise_intersect', include_unknown=True,
                   a_brake=9.0, vsafe_mode='curvature', ay_limit=8.0, default_look=200.0)
