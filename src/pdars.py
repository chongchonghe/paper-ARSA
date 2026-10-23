"""Von Neumann analysis of the IMEX PD-ARS radiation update with a local
reduced speed of light, derived without any assumed form for phi.

Units: one grid cell is the unit of length (dx = 1), so the continuum
system  d_t W + A d_x W = tau_c * Stilde W  (tau_c = cell optical depth)
is the paper's linearized system (sec 2.3) with ell = dx, i.e. calL = tau_c.

Parameters of a frozen-coefficient cell:
  nu    = phi * C * dt      local CFL number at the local reduced speed
  nuc   = C * dt            CFL number at the true speed of light (= nu/phi)
  tau_c = cell optical depth
  P     = radiation-to-gas pressure parameter calP
  chi   = Eddington factor (1/3 isotropic ... 1 free streaming); the radiation
          signal speed is sigma = sqrt(chi) (fraction of phi*C).

The radiation transport is the HLL/Rusanov flux with wave speeds +-sigma*phi*C
(identical to characteristic upwinding for the symmetric radiation block).
Matter transport is done by the separate hydro step and is omitted here
(optionally included with `hydro=True`).

Scheme (Honghu.jl imex_substep!, Quokka, He et al. 2024a):
  stage 2:  U1 = (I - Q)^-1 (I + Z) U^n
  stage 3:  U2 = (I - Q/2)^-1 [ U^n/2 + (I + Z) U1 / 2 ]
with Z the transport symbol times dt and Q = dt * tau_c * Stilde.
"""
import numpy as np

I5 = np.eye(5)


def matrices(M, phi, C, P, gamma, chi=1.0 / 3.0):
    """Advection matrix A and source matrix per unit optical depth Stilde.

    For chi = 1/3 this is the paper's (A, S/tau) for any M. For chi != 1/3 it
    is valid at M = 0 only (the v terms use T^4 + chi E)."""
    g1 = gamma - 1.0
    if abs(chi - 1.0 / 3.0) > 1e-14 and M != 0.0:
        raise ValueError("chi != 1/3 only implemented for M = 0")
    A = np.array([
        [M, 1, 0, 0, 0],
        [1, M, 1, 0, 0],
        [0, g1, M, 0, 0],
        [0, 0, 0, 0, phi * C],
        [0, 0, 0, phi * C * chi, 0],
    ], dtype=float)
    if M == 0.0:
        S = np.array([
            [0, 0, 0, 0, 0],
            [0, -(1 + chi) * P / C, 0, 0, P],
            [0, 0, -4.0 * C * P * g1, C * P * g1, 0],
            [0, 0, 4.0 * phi * C, -phi * C, 0],
            [0, (1 + chi) * phi, 0, 0, -phi * C],
        ], dtype=float)
    else:
        D = 3.0 * C**2 - M**2
        S = np.array([
            [0, 0, 0, 0, 0],
            [0, -4.0 * C * P / D, -4.0 * P * M / C, -P * M / (3.0 * C), P],
            [0, 0, -4.0 * (C**2 - M**2) * P * g1 / C,
             (3.0 * C**2 + M**2) * P * g1 / (3.0 * C), -2.0 * P * M * g1],
            [0, 4.0 * phi * C * M / D, 4.0 * phi * C, -phi * C, phi * M],
            [0, 4.0 * phi * C**2 / D, 4.0 * phi * M, phi * M / 3.0, -phi * C],
        ], dtype=float)
    return A, S


def transport_Z(A, dt, theta, phi, C, chi, hydro=False, epsE=1.0):
    """dt * (spatial operator symbol) for the HLL radiation flux (speeds
    +-sqrt(chi) phi C) and, optionally, upwind matter transport."""
    Z = np.zeros((5, 5), dtype=complex)
    sig = np.sqrt(chi) * phi * C
    Ar = A[3:, 3:]
    Z[3:, 3:] = -dt * (1j * np.sin(theta) * Ar + (1 - np.cos(theta)) * sig * np.diag([epsE, 1.0]))
    if hydro:
        Am = A[:3, :3]
        lam, R = np.linalg.eig(Am)
        absA = (R @ np.diag(np.abs(lam)) @ np.linalg.inv(R)).real
        Z[:3, :3] = -dt * (1j * np.sin(theta) * Am + (1 - np.cos(theta)) * absA)
    return Z


def amplification(theta, tau_c, nu, nuc, P, C=1e4, M=0.0, gamma=5 / 3,
                  chi=1 / 3, scheme="pdars", hydro=False, ap=False):
    """Amplification matrix of one radiation substep.

    scheme: 'pdars'   -- the actual IMEX PD-ARS update (default)
            'euler'   -- forward Euler transport + backward Euler source
            'rk2mid'  -- SSP-RK2 with a full-dt source after each stage
                         (the model scheme of the earlier draft)
            'rk2end'  -- SSP-RK2 transport, source once at the end
    """
    phi = nu / nuc
    dt = nuc / C
    A, S = matrices(M, phi, C, P, gamma, chi)
    Q = dt * tau_c * S
    Z = transport_Z(A, dt, theta, phi, C, chi, hydro, min(1.0, 1.0 / tau_c) if ap else 1.0)
    E = I5 + Z
    if scheme == "pdars":
        G1 = np.linalg.solve(I5 - Q, E)
        return np.linalg.solve(I5 - 0.5 * Q, 0.5 * I5 + 0.5 * E @ G1)
    if scheme == "euler":
        return np.linalg.solve(I5 - Q, E)
    if scheme == "rk2mid":
        G1 = np.linalg.solve(I5 - Q, E)
        return np.linalg.solve(I5 - Q, 0.5 * I5 + 0.5 * E @ G1)
    if scheme == "rk2end":
        return np.linalg.solve(I5 - Q, I5 + Z + 0.5 * Z @ Z)
    if scheme == "lw":
        # Lax-Wendroff radiation transport followed by a backward-Euler source
        Ar = A[3:, 3:] * dt
        L = np.eye(5, dtype=complex)
        L[3:, 3:] = np.eye(2) - 1j * np.sin(theta) * Ar - (1 - np.cos(theta)) * Ar @ Ar
        return np.linalg.solve(I5 - Q, L)
    raise ValueError(scheme)


def growth(theta, tau_c, nu, nuc, P, **kw):
    return np.max(np.abs(np.linalg.eigvals(amplification(theta, tau_c, nu, nuc, P, **kw))))


THETA = np.pi * np.linspace(0.05, 1.0, 20)


def gmax_theta(tau_c, nu, nuc, P, thetas=THETA, **kw):
    return max(growth(t, tau_c, nu, nuc, P, **kw) for t in thetas)


def nu_max(tau_c, nuc, P, tol=1e-10, nu_hi=None, thetas=THETA, rel=2e-3, **kw):
    """Largest local CFL number nu = phi*C*dt (phi <= 1, i.e. nu <= nuc) for
    which the cell is stable, by a monotone scan then bisection. Returns nuc
    if phi = 1 is stable."""
    stable = lambda nu: gmax_theta(tau_c, nu, nuc, P, thetas, **kw) <= 1 + tol
    hi = nuc if nu_hi is None else min(nu_hi, nuc)
    if stable(hi):
        return hi
    # geometric scan from tiny nu upward to find the first unstable value
    grid = np.geomspace(1e-3, hi, 60)
    lo = 0.0
    for nu in grid:
        if stable(nu):
            lo = nu
        else:
            hi = nu
            break
    while hi - lo > rel * hi:
        mid = 0.5 * (lo + hi)
        if stable(mid):
            lo = mid
        else:
            hi = mid
    return lo
