"""Check that the 2D checkerboard sets eq:bound_md: gray, M = 0, radiation step only,
variables (rho, vx, vy, T, E, Fx, Fy). Random cases; for each, no mode grows at 0.999 phi_max
and the fastest mode at 1.001 phi_max is (pi, pi)."""
# 2D, M=0, gray, radiation step only: variables (rho, vx, vy, T, E, Fx, Fy)
import numpy as np
rng = np.random.default_rng(3)
def build(phi, C, P, g1, L, dt, chi, dx, dy):
    S = np.zeros((7, 7))
    for v, F in [(1, 5), (2, 6)]:
        S[v, v] = -(1+chi)*P/C; S[v, F] = P
        S[F, v] = (1+chi)*phi; S[F, F] = -phi*C
    S[3, 3] = -4*C*P*g1; S[3, 4] = C*P*g1; S[4, 3] = 4*phi*C; S[4, 4] = -phi*C
    S *= dt*L
    Ax = np.zeros((7, 7)); Ax[4, 5] = 1; Ax[5, 4] = chi
    Ay = np.zeros((7, 7)); Ay[4, 6] = 1; Ay[6, 4] = chi
    Pr = np.diag([0, 0, 0, 0, 1, 1, 1.0])
    s = np.sqrt(chi)
    nux, nuy = phi*C*dt/dx, phi*C*dt/dy
    return S, Ax, Ay, Pr, s, nux, nuy
def rho_max(args, th):
    S, Ax, Ay, Pr, s, nux, nuy = args
    I = np.eye(7)
    A = np.linalg.inv(I - S/2); Bi = np.linalg.inv(I - S)
    tx, ty = th[:, 0][:, None, None], th[:, 1][:, None, None]
    Z = -(nux*(1j*np.sin(tx)*Ax + s*(1-np.cos(tx))*Pr) + nuy*(1j*np.sin(ty)*Ay + s*(1-np.cos(ty))*Pr))
    IZ = I + Z
    G = A @ (I/2 + IZ @ Bi @ IZ / 2)
    r = np.abs(np.linalg.eigvals(G)).max(axis=1)
    return r
t = np.linspace(0, np.pi, 33); TH = np.array([(a, b) for a in t for b in t])[1:]
worst_below, n_pi, n_act = 0.0, 0, 0
for _ in range(150):
    chi = rng.choice([1/3, 1.0]); s = np.sqrt(chi); g1 = rng.choice([2/3, 0.4])
    C = 1e4; P = 10**rng.uniform(-6, 0); L = 10**rng.uniform(-2, 3)
    dx, dy = 1.0, 10**rng.uniform(-0.5, 0.5)
    nuc = 10**rng.uniform(1.5, 3)                 # c dt / dx
    dt = nuc*dx/C
    K = 4*g1*P*C*dt*L
    bE = C*dt*(s/dx + s/dy - L/(2+K)); bF = C*dt*(s/dx + s/dy - L/2)
    b = max(bE, bF)
    if b <= 1: continue                            # phi = 1 stable; boundary inactive
    phimax = 1/b; n_act += 1
    r_lo = rho_max(build(0.999*phimax, C, P, g1, L, dt, chi, dx, dy), TH)
    worst_below = max(worst_below, r_lo.max() - 1)
    r_hi = rho_max(build(1.001*phimax, C, P, g1, L, dt, chi, dx, dy), TH)
    n_pi += np.allclose(TH[r_hi.argmax()], [np.pi, np.pi])
print(f"active {n_act}, max(rho-1) below boundary {worst_below:.2e}, argmax at (pi,pi) above: {n_pi}/{n_act}")
