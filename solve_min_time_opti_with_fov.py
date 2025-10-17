# -*- coding: utf-8 -*-
"""
Opti-based global OCP with FOV planning band and FPV→safe-speed braking constraint.

This version uses casadi.Opti() (high-level modeling) instead of the low-level
nlpsol with manual lbw/ubw stacks. It keeps the same vehicle model and reference
path utilities from solve_min_time_with_fov.py, but enforces the convex planning
band (road ∩ FOV ∩ Unknown) via apply_planning_band_to_opti().

Braking constraint (flexible): at each s_k,
    2 a_brake ( s_fpv[k] - s_k ) + v_safe[k]^2 - v(k)^2  >= 0
where v_safe can be a constant or curvature-based safe speed.
"""
import argparse
import numpy as np
import casadi as ca

from fov_to_frenet_bounds import apply_planning_band_to_opti
# reuse utilities and model from the existing script (only functions are imported)
from solve_min_time_with_fov import (
    generate_reference_path,
    vehicle_dynamics,
    visualize_track,
    fpv_s_on_grid,
    build_vsafe_array,
)

MONZA_REF_PATH_FILENAME = 'Monza.csv'
MONZA_OBSTACLES_FILENAME = 'MonzaObstacles.csv'
NUERBURGRING_REF_PATH_FILENAME = 'Nuerburgring.csv'
NUERBURGRING_OBSTACLES_FILENAME = 'NuerburgringObstacles.csv'


def solve_min_time_opti(
    ref_path_filename='Nuerburgring.csv',
    obstacles_filename='NuerburgringObstacles.csv',
    ds=3.0,
    fov_cache_path='fov_cache.json',
    fov_mode='hard',                   # 'hard' or 'soft'
    lam_fov=5e2,
    blend_mode='pairwise_intersect',
    include_unknown=True,
    a_brake=10.0,
    vsafe_mode='constant',             # 'constant' or 'curvature'
    vsafe_const=0.0,
    ay_limit=8.0,
    default_look=200.0,
    ipopt_max_iter=2000,
):
    # === Reference path & geometry ===
    s, theta, kappa, new_ref_path, init_ref_path, obstacles = generate_reference_path(
        ref_path_filename=ref_path_filename,
        obstacles_filename=obstacles_filename,
        ds=ds,
        max_points=1160,
    )
    N = len(s) - 1

    # === FPV map & safe speed field ===
    s_fpv = fpv_s_on_grid(fov_cache_path, s, default_look=default_look)
    v_safe_arr = build_vsafe_array(vsafe_mode, s, kappa, a_y_lim=ay_limit,
                                   v_const=vsafe_const, s_fpv=s_fpv)

    # === Opti model ===
    opti = ca.Opti()

    nx, nu = 5, 2   # [v, beta, omega, n, xi], [delta, F_dr]
    X = opti.variable(nx, N+1)
    U = opti.variable(nu, N)

    # initial state
    x0_val = np.array([10.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
    opti.subject_to(X[:, 0] == x0_val)

    # control bounds
    delta_max = np.pi/6
    F_min, F_max = -21000.0, 7100.0
    opti.subject_to(opti.bounded(-delta_max, U[0, :],  delta_max))
    opti.subject_to(opti.bounded(F_min,      U[1, :],  F_max))

    # state box bounds (except n handled by FOV band)
    v_max = 42.5
    opti.subject_to(opti.bounded(0.0,  X[0, :], v_max))     # v
    opti.subject_to(opti.bounded(-ca.inf, X[1, :], ca.inf)) # beta
    opti.subject_to(opti.bounded(-2.0, X[2, :], 2.0))       # omega
    opti.subject_to(opti.bounded(-ca.inf, X[4, :], ca.inf)) # xi

    # planning band (road ∩ FOV ∩ Unknown) as hard/soft constraints on n = X[3,:]
    J = 0
    ret = apply_planning_band_to_opti(
        opti=opti,
        n_var=X[3, :],
        s_grid=s,
        x_ref=new_ref_path.x_coordinates,
        y_ref=new_ref_path.y_coordinates,
        theta=theta,
        wr=new_ref_path.dist_to_right_bound,
        wl=new_ref_path.dist_to_left_bound,
        cache_json_path=fov_cache_path,
        mode=fov_mode,
        lam_fov=lam_fov,
        blend_mode=blend_mode,
        include_unknown=include_unknown,
    )
    if fov_mode == 'soft' and 'soft_penalty' in ret:
        J = J + ret['soft_penalty']

    # weights
    lambda_ddelta = 10.0
    lambda_dF     = 0.01
    lambda_dn     = 10.0
    lambda_dxi    = 5.0
    lambda_dv     = 1.0
    lambda_jerk   = 3.0

    # multiple-shooting dynamics and braking-to-safe-speed constraint
    for k in range(N):
        Xk   = X[:, k]
        Xk1  = X[:, k+1]
        Uk   = U[:, k]
        kap  = float(kappa[k])

        f_k, SF_k, FY_fl, FY_fr, FY_rl, FY_rr = vehicle_dynamics(Xk, Uk, kap)
        # explicit Euler
        opti.subject_to(Xk1 == Xk + ds * f_k)

        # time objective — keep same form as original (sum ds * SF)
        J = J + ds * SF_k

        # smoothness penalties
        if k >= 1:
            J = J + lambda_ddelta * (U[0, k] - U[0, k-1])**2
            J = J + lambda_dF     * (U[1, k] - U[1, k-1])**2
            J = J + lambda_dn     * (X[3, k] - X[3, k-1])**2
            J = J + lambda_dxi    * (X[4, k] - X[4, k-1])**2
            J = J + lambda_dv     * (X[0, k] - X[0, k-1])**2
        if k >= 2:
            J = J + lambda_jerk   * (U[0, k] - 2*U[0, k-1] + U[0, k-2])**2

        # braking-to-safe-speed within visibility
        # 2*a_brake*(s_fpv - s_k) + v_safe^2 - v_k^2 >= 0
        lhs = 2.0*a_brake*( s_fpv[k] - s[k] ) + (v_safe_arr[k]**2) - (X[0, k]**2)
        opti.subject_to(lhs >= 0)

        # tire ellipse (rear-left & rear-right share drive force)
        F_dr = U[1, k]
        m = 1500.0; g = 9.81; mu = 1.0
        F_z = m * g / 4.0
        F_x_max = 14000.0/2.0
        F_y_max = mu * F_z
        F_x_rl = F_dr / 2.0
        F_x_rr = F_dr / 2.0
        ellipse_rl = (F_x_rl / F_x_max)**2 + (FY_rl / F_y_max)**2
        ellipse_rr = (F_x_rr / F_x_max)**2 + (FY_rr / F_y_max)**2
        opti.subject_to(ellipse_rl <= 1.0)
        opti.subject_to(ellipse_rr <= 1.0)

    # solver
    p_opts = { }
    s_opts = { 'ipopt.print_level': 5, 'print_time': False, 'ipopt.max_iter': int(ipopt_max_iter) }
    opti.minimize(J)
    opti.solver('ipopt', p_opts, s_opts)

    sol = opti.solve()

    X_opt = sol.value(X).T  # shape (N+1, nx)
    U_opt = sol.value(U).T

    n_opt = X_opt[:, 3]
    x_opt = new_ref_path.x_coordinates - n_opt * np.sin(theta)
    y_opt = new_ref_path.y_coordinates + n_opt * np.cos(theta)

    visualize_track(X_opt, new_ref_path, x_opt, y_opt, theta, s, obstacles, init_ref_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--track', choices=['monza','nuerburgring'], default='nuerburgring')
    parser.add_argument('--ds', type=float, default=3.0)
    parser.add_argument('--fov_cache', type=str, default='fov_cache.json')
    parser.add_argument('--fov_mode', choices=['hard','soft'], default='hard')
    parser.add_argument('--lam_fov', type=float, default=5e2)
    parser.add_argument('--blend_mode', choices=['nearest','pairwise_intersect'], default='pairwise_intersect')
    parser.add_argument('--include_unknown', action='store_true', default=True)
    parser.add_argument('--a_brake', type=float, default=10.0)
    parser.add_argument('--vsafe_mode', choices=['constant','curvature'], default='constant')
    parser.add_argument('--vsafe_const', type=float, default=0.0)
    parser.add_argument('--ay_limit', type=float, default=8.0)
    parser.add_argument('--default_look', type=float, default=200.0)
    parser.add_argument('--ipopt_max_iter', type=int, default=2000)
    args = parser.parse_args()

    if args.track == 'monza':
        ref_path_filename = MONZA_REF_PATH_FILENAME
        obstacles_filename = MONZA_OBSTACLES_FILENAME
    else:
        ref_path_filename = NUERBURGRING_REF_PATH_FILENAME
        obstacles_filename = NUERBURGRING_OBSTACLES_FILENAME

    solve_min_time_opti(
        ref_path_filename=ref_path_filename,
        obstacles_filename=obstacles_filename,
        ds=args.ds,
        fov_cache_path=args.fov_cache,
        fov_mode=args.fov_mode,
        lam_fov=args.lam_fov,
        blend_mode=args.blend_mode,
        include_unknown=args.include_unknown,
        a_brake=args.a_brake,
        vsafe_mode=args.vsafe_mode,
        vsafe_const=args.vsafe_const,
        ay_limit=args.ay_limit,
        default_look=args.default_look,
        ipopt_max_iter=args.ipopt_max_iter,
    )


if __name__ == '__main__':
    main()