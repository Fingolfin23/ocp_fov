# compare_tire_models.py
# ------------------------------------------------------------
# Compare linear vs. nonlinear (Pacejka) tire models in the same
# time-optimal OCP. Adds switches to expose differences:
#   • --mu_cap {both,linear,nonlinear,none}
#       Control whether the per-wheel μ-ellipse constraint is enforced
#       for each model. If you set --mu_cap nonlinear, the linear model
#       can exceed μFz and will typically achieve unrealistically short
#       times; the nonlinear model remains capped by physics.
#   • --front_brake_frac f  (default 0.6)
#       During braking (F_dr<0) send fraction f of the braking force
#       to the front axle to induce combined slip on the front tires.
#   • --match_slope
#       Match linear small-angle stiffness to Pacejka via Cα ≈ Fz·B·C·D
#       so near-zero behavior is comparable, highlighting saturation.
#   • --Calpha_scale k
#       Optional multiplier on the matched Cα to stress-test sensitivity.
# Plots include velocity/trajectory, μ-utilization, Fy–α samples,
# and normalized friction-ellipse usage.
# ------------------------------------------------------------

# Example to reveal differences (linear uncapped, combined slip on fronts):
   

'''
python compare_tire_models.py --csv Nuerburgring.csv --idx_from 110 --idx_to 400 \
   --ds 1.0 --mu_cap nonlinear --front_brake_frac 0.7 --match_slope
python compare_tire_models.py --csv Monza.csv --idx_from 770 --idx_to 900 \
   --ds 1.0 --mu_cap nonlinear --front_brake_frac 0.7 --match_slope

python compare_tire_models.py --csv Nuerburgring.csv \
  --idx_from 110 --idx_to 250 --ds 1.0 \
  --ipopt_iter 3000 \
  --mu_cap nonlinear --front_brake_frac 0.7 --match_slope


python compare_tire_models.py --csv Nuerburgring.csv \
  --idx_from 110 --idx_to 250 --ds 1.0 \
  --ipopt_iter 2000 --mu_cap nonlinear --front_brake_frac 0.6 --match_slope
'''

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import argparse

# -------------------------------
# Utilities: reference path
# -------------------------------
def resample_reference(csv_file, ds=0.5, idx_from=None, idx_to=None):
    """Load CSV: x,y,wr,wl; slice by index (NOT meters!), resample to uniform ds."""
    data = np.loadtxt(csv_file, delimiter=',')
    if idx_from is not None or idx_to is not None:
        data = data[idx_from:idx_to]

    x_ref = data[:, 0]
    y_ref = data[:, 1]
    wr    = data[:, 2]
    wl    = data[:, 3]

    dx = np.diff(x_ref)
    dy = np.diff(y_ref)
    ds_raw = np.hypot(dx, dy)
    s_raw = np.insert(np.cumsum(ds_raw), 0, 0.0)
    s_max = s_raw[-1]

    s_u = np.arange(0.0, s_max, ds)
    if s_u[-1] < s_max:
        s_u = np.append(s_u, s_max)

    fx = interp1d(s_raw, x_ref, kind='cubic')
    fy = interp1d(s_raw, y_ref, kind='cubic')
    fwr = interp1d(s_raw, wr, kind='linear')
    fwl = interp1d(s_raw, wl, kind='linear')

    x_u  = fx(s_u)
    y_u  = fy(s_u)
    wr_u = fwr(s_u)
    wl_u = fwl(s_u)

    dx_ds = np.gradient(x_u, s_u)
    dy_ds = np.gradient(y_u, s_u)
    theta = np.arctan2(dy_ds, dx_ds)
    kappa = np.gradient(theta, s_u)

    return s_u, x_u, y_u, theta, kappa, wr_u, wl_u


# -------------------------------
# Tire models
# -------------------------------
def tire_linear(alpha, C_alpha):
    # per-tire linear; alpha [rad], C_alpha [N/rad]
    return C_alpha * alpha

def tire_pacejka(alpha, Fz, B=10.0, C=1.9, D=1.0, E=0.97, mu_scale=1.0):
    # Magic Formula lateral; returns Fy [N]
    # D effectively scales with mu; we include mu_scale if you want to < 1 on wet etc.
    return mu_scale * Fz * D * ca.sin(C * ca.atan(B*alpha - E*(B*alpha - ca.atan(B*alpha))))

def per_wheel_slip_angles(v, beta, omega, delta_f, lf, lr, wf, wr):
    # return alpha_fl, alpha_fr, alpha_rl, alpha_rr
    # Front axle (y positive to left), left uses -wf/2 in denominator, right +wf/2
    # dy = v*sin(beta) +/- l*omega; dx = v*cos(beta) +/- (w/2)*omega
    dy_f = v*ca.sin(beta) + lf*omega
    dx_fl = v*ca.cos(beta) - (wf/2.0)*omega
    dx_fr = v*ca.cos(beta) + (wf/2.0)*omega
    a_fl = delta_f - ca.atan2(dy_f, dx_fl)
    a_fr = delta_f - ca.atan2(dy_f, dx_fr)

    dy_r = v*ca.sin(beta) - lr*omega
    dx_rl = v*ca.cos(beta) - (wr/2.0)*omega
    dx_rr = v*ca.cos(beta) + (wr/2.0)*omega
    a_rl = 0.0 - ca.atan2(dy_r, dx_rl)  # rear steer delta_r = 0
    a_rr = 0.0 - ca.atan2(dy_r, dx_rr)
    return a_fl, a_fr, a_rl, a_rr

# -------------------------------
# Dynamics (Frenet with tire model switch)
# -------------------------------
def dynamics_frenet(x, u, kappa, params, tire_kind):
    """
    x = [v, beta, omega, n, xi]
    u = [delta_f, F_dr]
    tire_kind in {'linear','pacejka'}
    params must include 'front_brake_frac' (0..1)
    returns: f(x), SF, Fy tuple, Fx tuple, (alphas)
    """
    v, beta, omega, n, xi = x[0], x[1], x[2], x[3], x[4]
    delta_f, F_dr = u[0], u[1]

    # vehicle params
    m    = params['m']
    Jzz  = params['Jzz']
    lf   = params['lf']; lr = params['lr']
    wf   = params['wf']; wr = params['wr']
    mu_r = params['mu_roll']
    rho  = params['rho']; Cd = params['Cd']; Af = params['Af']
    mu   = params['mu']
    eps  = 1e-6

    # static load split (flat)
    Fz_per_wheel = (m*params['g'])/4.0

    # slip angles
    a_fl, a_fr, a_rl, a_rr = per_wheel_slip_angles(v, beta, omega, delta_f, lf, lr, wf, wr)

    # lateral forces by tire model
    if tire_kind == 'linear':
        Cf_tire_f = params['Calpha_f']/2.0  # per tire
        Cf_tire_r = params['Calpha_r']/2.0
        Fy_fl = tire_linear(a_fl, Cf_tire_f)
        Fy_fr = tire_linear(a_fr, Cf_tire_f)
        Fy_rl = tire_linear(a_rl, Cf_tire_r)
        Fy_rr = tire_linear(a_rr, Cf_tire_r)
    else:
        B,C,D,E = params['B'], params['C'], params['D'], params['E']
        Fy_fl = tire_pacejka(a_fl, Fz_per_wheel, B,C,D,E, mu_scale=1.0)
        Fy_fr = tire_pacejka(a_fr, Fz_per_wheel, B,C,D,E, mu_scale=1.0)
        Fy_rl = tire_pacejka(a_rl, Fz_per_wheel, B,C,D,E, mu_scale=1.0)
        Fy_rr = tire_pacejka(a_rr, Fz_per_wheel, B,C,D,E, mu_scale=1.0)

    # Longitudinal force distribution with combined slip on fronts during braking
    fb = params.get('front_brake_frac', 0.6)  # fraction of |brake| to front axle
    F_brake = -ca.fmin(F_dr, 0.0)  # ≥0 when braking
    F_drive =  ca.fmax(F_dr, 0.0)  # ≥0 when driving

    # Split equally per axle, apply signs (negative Fx for braking)
    Fx_fl = -0.5*fb*F_brake
    Fx_fr = -0.5*fb*F_brake
    Fx_rl =  0.5*F_drive - 0.5*(1.0-fb)*F_brake
    Fx_rr =  0.5*F_drive - 0.5*(1.0-fb)*F_brake

    # Frenet scale factor
    SF = (1.0 - n*kappa) / (v*ca.cos(xi + beta) + eps)

    # drags
    F_roll = mu_r * m * params['g']
    F_aero = 0.5 * rho * Cd * Af * v**2

    # state derivatives
    dot_n    = SF * v * ca.sin(xi + beta)
    dot_xi   = SF * omega - kappa
    dot_beta = SF * ( -omega + 1.0/(m*v + eps) *
                     ((Fy_fl + Fy_fr)*ca.cos(delta_f - beta) +
                      (Fy_rl + Fy_rr)*ca.cos(0.0      - beta) - F_dr*ca.sin(beta)) )
    dot_omega = SF / Jzz *(
        Fy_fl*(lf*ca.cos(delta_f) - (wf/2.0)*ca.sin(delta_f)) +
        Fy_fr*(lf*ca.cos(delta_f) + (wf/2.0)*ca.sin(delta_f)) +
        Fy_rl*(-lr*ca.cos(0.0)    - (wr/2.0)*ca.sin(0.0))     +
        Fy_rr*(-lr*ca.cos(0.0)    + (wr/2.0)*ca.sin(0.0))
    )
    dot_v = SF / m * (
        (Fy_fl + Fy_fr)*ca.sin(beta - delta_f) +
        (Fy_rl + Fy_rr)*ca.sin(beta - 0.0) + F_dr*ca.cos(beta) - F_roll - F_aero
    )

    f = ca.vertcat(dot_v, dot_beta, dot_omega, dot_n, dot_xi)
    Fy_tuple = (Fy_fl, Fy_fr, Fy_rl, Fy_rr)
    Fx_tuple = (Fx_fl, Fx_fr, Fx_rl, Fx_rr)
    return f, SF, Fy_tuple, Fx_tuple, (a_fl,a_fr,a_rl,a_rr)

# -------------------------------
# Minimal "apex-biased" lateral initial guess (keeps code small)
# -------------------------------
def build_apex_n_guess(kappa, wl, wr, gain=0.25, margin=0.20):
    """
    Create a tiny lateral offset toward the inside of each bend to avoid the
    optimizer sticking to the centerline when initialized with n=0.
    Positive curvature -> left turn -> bias to + (toward wl);
    Negative curvature -> right turn -> bias to - (toward wr).
    Returns an array n0 with same length as kappa.
    """
    n0 = np.zeros_like(kappa, dtype=float)
    left_mask  = (np.asarray(kappa) > 0.0)
    right_mask = (np.asarray(kappa) < 0.0)
    if np.any(left_mask):
        # stay within left width minus a small margin
        n0[left_mask] = np.minimum(gain * wl[left_mask], wl[left_mask] - margin)
    if np.any(right_mask):
        # stay within right width minus a small margin (negative toward right)
        n0[right_mask] = -np.minimum(gain * wr[right_mask], wr[right_mask] - margin)
    # light smoothing to avoid abrupt jumps
    if n0.size >= 4:
        kern = np.array([1.0, 2.0, 2.0, 1.0], dtype=float)
        kern /= kern.sum()
        n0 = np.convolve(n0, kern, mode='same')
    return n0

# -------------------------------
# Build & solve OCP once
# -------------------------------
def solve_ocp_once(s, theta, kappa, wr, wl, tire_kind, ds,
                   v0=20.0, params=None, ipopt_max_iter=400,
                   mu_cap_mode='both', match_slope=False):
    N = len(s)-1
    nx, nu = 5, 2
    n0_guess = build_apex_n_guess(kappa, wl, wr, gain=0.25, margin=0.20)

    # decision stacks
    w, lbw, ubw, w0 = [], [], [], []
    g, lbg, ubg = [], [], []

    # cost
    J = 0
    # parameters fallback
    p = dict(
        m=1200.0, Jzz=1260.0, lf=1.5, lr=1.4, wf=1.6, wr=1.5,
        g=9.81, mu_roll=0.015, rho=1.2041, Cd=0.30, Af=1.0,
        mu=0.95,
        Calpha_f=60000.0, Calpha_r=60000.0, # axle [N/rad], per tire/2
        B=10.0, C=1.9, D=1.0, E=0.97,       # Pacejka
        v_max=42.5, delta_max=np.deg2rad(20.0), F_dr_min=-21000.0, F_dr_max=7100.0,
        front_brake_frac=0.6,
    )
    if params is not None:
        p.update(params)

    # Optional: match linear small-angle stiffness to Pacejka near α≈0
    if tire_kind == 'linear' and match_slope:
        Fz_w = (p['m']*p['g'])/4.0
        C_linear_per_tire = Fz_w * p['B'] * p['C'] * p['D']  # from small-angle of Pacejka
        scale = float(p.get('Calpha_scale', 1.0))
        C_linear_per_tire *= scale
        p['Calpha_f'] = 2.0 * C_linear_per_tire
        p['Calpha_r'] = 2.0 * C_linear_per_tire

    # initial state
    Xk = ca.MX.sym("X0", nx)
    w += [Xk]
    x0 = [v0, 0.0, 0.0, float(n0_guess[0]), 0.0]
    lbw += x0; ubw += x0; w0 += x0

    # per-wheel normal load (flat)
    Fz_w = (p['m']*p['g'])/4.0
    mu = p['mu']

    # small regularization
    lam_delta = 1.0
    lam_F     = 1e-4
    lam_state = 1e-2

    for k in range(N):
        # input
        Uk = ca.MX.sym(f"U_{k}", nu)
        w += [Uk]
        lbw += [-p['delta_max'], p['F_dr_min']]
        ubw += [ p['delta_max'], p['F_dr_max']]
        w0  += [0.0,  1000.0]

        # next state
        Xk_next = ca.MX.sym(f"X_{k+1}", nx)
        w += [Xk_next]
        lbw += [0.0, -ca.inf, -ca.inf, -wr[k], -ca.inf]
        ubw += [p['v_max'], ca.inf, ca.inf,  wl[k], ca.inf]
        w0  += [v0, 0.0, 0.0, float(n0_guess[k+1]), 0.0]

        # dynamics
        f_k, SF_k, Fy_k, Fx_k, alphas = dynamics_frenet(Xk, Uk, float(kappa[k]), p, tire_kind)
        g += [Xk + ds*f_k - Xk_next]
        lbg += [0.0]*nx
        ubg += [0.0]*nx

        # time objective
        J += ds * SF_k

        # friction ellipse (all four wheels), only apply depending on mu_cap_mode and tire_kind
        Fx_fl, Fx_fr, Fx_rl, Fx_rr = Fx_k
        Fy_fl, Fy_fr, Fy_rl, Fy_rr = Fy_k
        apply_mu = (mu_cap_mode == 'both') or (mu_cap_mode == tire_kind)
        if apply_mu:
            for Fx_i, Fy_i in [(Fx_fl,Fy_fl), (Fx_fr,Fy_fr), (Fx_rl,Fy_rl), (Fx_rr,Fy_rr)]:
                mu_used_sq = (Fx_i/(mu*Fz_w))**2 + (Fy_i/(mu*Fz_w))**2
                g += [mu_used_sq]
                lbg += [-ca.inf]
                ubg += [1.0]

        # mild smoothness
        if k>0:
            J += lam_delta*(Uk[0] - U_prev[0])**2 + lam_F*(Uk[1]-U_prev[1])**2
            J += lam_state*ca.sumsqr(Xk - X_prev)

        U_prev = Uk; X_prev = Xk
        Xk = Xk_next

    nlp = dict(f=J, x=ca.vertcat(*w), g=ca.vertcat(*g))
    solver = ca.nlpsol("solver","ipopt", nlp,
        {"ipopt.max_iter": ipopt_max_iter, "ipopt.print_level":0, "print_time":0})

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    xopt = sol['x'].full().ravel()

    # unpack
    # pattern per stage: [Uk, Xk+1]; with initial X0
    blocks = []
    ptr = 0
    X_list = [xopt[ptr:ptr+nx]]; ptr += nx
    for k in range(N):
        U_k = xopt[ptr:ptr+nu]; ptr += nu
        Xn  = xopt[ptr:ptr+nx]; ptr += nx
        blocks.append((X_list[-1], U_k))
        X_list.append(Xn)

    X_arr = np.array(X_list)               # (N+1, nx)
    U_arr = np.array([b[1] for b in blocks])  # (N, nu)

    # recompute metrics along solution
    mu_used = {'fl':[], 'fr':[], 'rl':[], 'rr':[]}
    alpha_log = {'fl':[], 'fr':[], 'rl':[], 'rr':[]}
    Fy_log = {'fl':[], 'fr':[], 'rl':[], 'rr':[]}
    SF_sum = 0.0
    mu_series = {'fl':[], 'fr':[], 'rl':[], 'rr':[]}

    for k in range(N):
        Xk = X_arr[k,:]
        Uk = U_arr[k,:]
        f_k, SF_k, Fy_k, Fx_k, alphas = dynamics_frenet(ca.vertcat(*Xk), ca.vertcat(*Uk),
                                                        float(kappa[k]), p, tire_kind)
        SF_sum += (ds * float(SF_k))
        Fx_fl, Fx_fr, Fx_rl, Fx_rr = [float(v) for v in Fx_k]
        Fy_fl, Fy_fr, Fy_rl, Fy_rr = [float(v) for v in Fy_k]
        a_fl, a_fr, a_rl, a_rr = [float(v) for v in alphas]

        for key,(Fx,Fy) in zip(['fl','fr','rl','rr'],
                               [(Fx_fl,Fy_fl),(Fx_fr,Fy_fr),(Fx_rl,Fy_rl),(Fx_rr,Fy_rr)]):
            mu_val = np.sqrt((Fx/(mu*Fz_w))**2 + (Fy/(mu*Fz_w))**2)
            mu_used[key].append(mu_val)
            mu_series[key].append(mu_val)
        for key,val in zip(['fl','fr','rl','rr'], [a_fl,a_fr,a_rl,a_rr]):
            alpha_log[key].append(val)
        for key,val in zip(['fl','fr','rl','rr'], [Fy_fl,Fy_fr,Fy_rl,Fy_rr]):
            Fy_log[key].append(val)

    result = dict(
        X=X_arr, U=U_arr, total_time=SF_sum,
        mu_used=mu_used, alpha=alpha_log, Fy=Fy_log,
        params=p,
        mu_series=mu_series,
    )
    return result

# -------------------------------
# Plot helpers
# -------------------------------
def frenet_to_xy(x_ref, y_ref, theta, n_arr):
    x = x_ref - n_arr*np.sin(theta)
    y = y_ref + n_arr*np.cos(theta)
    return x,y

def plot_comparison(s, x_ref, y_ref, theta, wr, wl, res_lin, res_nl, title_prefix=""):
    # 1) Velocity + Trajectory
    x_l = x_ref - wr*np.sin(theta); y_l = y_ref + wr*np.cos(theta)
    x_r = x_ref + wl*np.sin(theta); y_r = y_ref - wl*np.cos(theta)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={'width_ratios': [1.0, 1.25]})
    # velocity
    ax[0].plot(s, res_lin['X'][:,0], label='Linear', lw=2)
    ax[0].plot(s, res_nl['X'][:,0], label='Nonlinear (Pacejka)', lw=2)
    ax[0].set_xlabel('s (m)'); ax[0].set_ylabel('velocity (m/s)')
    ax[0].set_title('Velocity Profile'); ax[0].grid(True); ax[0].legend()

    # trajectory
    x_lin, y_lin = frenet_to_xy(x_ref, y_ref, theta, res_lin['X'][:,3])
    x_nl,  y_nl  = frenet_to_xy(x_ref, y_ref, theta, res_nl['X'][:,3])
    ax[1].plot(x_ref, y_ref, 'k--', label='Reference')
    ax[1].plot(x_l, y_l, 'k-'); ax[1].plot(x_r, y_r, 'k-')
    ax[1].plot(x_lin, y_lin, color='#E24A33', label='Optimized (Linear)')
    ax[1].plot(x_nl,  y_nl,  color='#348ABD', label='Optimized (Nonlinear)')
    # annotate total Frenet arc length on the trajectory subplot
    path_len = float(s[-1] - s[0]) if len(s) > 0 else 0.0
    ax[1].text(0.02, 0.98,
               f"Path length ≈ {path_len:.1f} m",
               transform=ax[1].transAxes, va='top', ha='left',
               fontsize=10, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'))
    ax[1].set_aspect('equal', 'box')
    ax[1].set_title('Trajectory'); ax[1].grid(True); ax[1].legend()
    tt = f"{title_prefix}  T_linear={res_lin['total_time']:.2f}s,  T_nonlinear={res_nl['total_time']:.2f}s"
    fig.suptitle(tt)
    plt.tight_layout()

    # 2) mu-utilization (front-left as example)
    plt.figure(figsize=(12,4))
    plt.plot(s[:-1], res_lin['mu_used']['fl'], label='Linear')
    plt.plot(s[:-1], res_nl['mu_used']['fl'], label='Nonlinear')
    plt.axhline(1.0, color='r', ls='--', lw=1, label='μ limit')
    plt.xlabel('s (m)'); plt.ylabel('μ_used (FL)')
    plt.title('Front-Left wheel μ-utilization')
    plt.grid(True); plt.legend(); plt.tight_layout()

    # 2b) μ exceedance indicator (FL)
    plt.figure(figsize=(12,2.8))
    lin_mu = np.array(res_lin['mu_used']['fl'])
    nl_mu  = np.array(res_nl['mu_used']['fl'])
    plt.plot(s[:-1], np.maximum(0.0, lin_mu-1.0), label='Linear: max(0, μ-1)')
    plt.plot(s[:-1], np.maximum(0.0,  nl_mu-1.0), label='Nonlinear: max(0, μ-1)')
    plt.xlabel('s (m)'); plt.ylabel('exceedance')
    plt.title('Where the linear model would over-use friction if uncapped')
    plt.grid(True); plt.legend(); plt.tight_layout()

    # 3) Fy vs alpha (FL)
    plt.figure(figsize=(12,4))
    plt.plot(np.rad2deg(res_lin['alpha']['fl']), res_lin['Fy']['fl'], '.', label='Linear samples', alpha=0.7)
    plt.plot(np.rad2deg(res_nl['alpha']['fl']),  res_nl['Fy']['fl'],  '.', label='Nonlinear samples', alpha=0.7)
    plt.axhline( res_lin['params']['mu']*res_lin['params']['m']*res_lin['params']['g']/4.0, color='k', ls='--', lw=1, label='μFz')
    plt.axhline(-res_lin['params']['mu']*res_lin['params']['m']*res_lin['params']['g']/4.0, color='k', ls='--', lw=1)
    plt.xlabel('slip angle α (deg)'); plt.ylabel('Fy (N)')
    plt.title('Fy–α samples (Front-Left). Nonlinear saturates; Linear grows ~unbounded.')
    plt.grid(True); plt.legend(); plt.tight_layout()

    # 4) Normalized friction ellipse scatter (FL)
    def ellipse_points(res):
        Fz_w = (res['params']['m']*res['params']['g'])/4.0
        mu = res['params']['mu']
        Fx = np.zeros_like(res['Fy']['fl'])  # FL gets no drive/brake in this setup
        Fy = np.array(res['Fy']['fl'])
        return Fx/(mu*Fz_w), Fy/(mu*Fz_w)

    xl, yl = ellipse_points(res_lin)
    xn, yn = ellipse_points(res_nl)
    ang = np.linspace(0, 2*np.pi, 200)
    plt.figure(figsize=(4.5,4.5))
    plt.plot(np.cos(ang), np.sin(ang), 'k--', lw=1, label='unit ellipse')
    plt.scatter(xl, yl, s=10, alpha=0.5, label='Linear')
    plt.scatter(xn, yn, s=10, alpha=0.5, label='Nonlinear')
    plt.axis('equal'); plt.xlabel('Fx/(μFz)'); plt.ylabel('Fy/(μFz)')
    plt.title('Normalized friction usage (FL)')
    plt.grid(True); plt.legend(); plt.tight_layout()

# -------------------------------
# Main
# -------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', type=str, default='Monza.csv')
    ap.add_argument('--idx_from', type=int, default=770, help='row index start (inclusive)')
    ap.add_argument('--idx_to',   type=int, default=900, help='row index end (exclusive)')
    ap.add_argument('--ds', type=float, default=2.0)
    ap.add_argument('--ipopt_iter', type=int, default=400)
    ap.add_argument('--mu_cap', type=str, default='both', choices=['both','linear','nonlinear','none'])
    ap.add_argument('--front_brake_frac', type=float, default=0.6)
    ap.add_argument('--match_slope', action='store_true')
    ap.add_argument('--Calpha_scale', type=float, default=1.0)
    args = ap.parse_args()

    s, x, y, th, kap, wr, wl = resample_reference(args.csv, ds=args.ds,
                                                  idx_from=args.idx_from, idx_to=args.idx_to)
    # Solve both variants
    base_params = dict(mu=0.95, front_brake_frac=args.front_brake_frac,
                       B=10.0, C=1.9, D=1.0, E=0.97, Calpha_scale=args.Calpha_scale)

    res_lin = solve_ocp_once(s, th, kap, wr, wl, 'linear',  ds=args.ds,
                             params=dict(base_params, Calpha_f=60000.0, Calpha_r=60000.0),
                             ipopt_max_iter=args.ipopt_iter,
                             mu_cap_mode=('linear' if args.mu_cap=='linear' else
                                          ('both' if args.mu_cap=='both' else
                                           ('none' if args.mu_cap=='none' else 'none'))),
                             match_slope=args.match_slope)

    res_nl  = solve_ocp_once(s, th, kap, wr, wl, 'pacejka', ds=args.ds,
                             params=base_params,
                             ipopt_max_iter=args.ipopt_iter,
                             mu_cap_mode=('nonlinear' if args.mu_cap=='nonlinear' else
                                          ('both' if args.mu_cap=='both' else
                                           ('none' if args.mu_cap=='none' else 'none'))),
                             match_slope=False)

    # Visualize
    title = f"{args.csv}  idx[{args.idx_from}:{args.idx_to}]  ds={args.ds}"
    plot_comparison(s, x, y, th, wr, wl, res_lin, res_nl, title_prefix=title)
    plt.show()

if __name__ == "__main__":
    main()