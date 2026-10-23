"""Energy block at theta = pi with NG groups sharing one gas temperature.
Variables (T, E_1..E_NG), background E_g = 4 w_g T (linearised).  Units as in
pdars.py: q_g = nu_g tau_g, gas heating dT = g1*P*C*dt*sum tau_g (E_g - 4 w_g T).
Question: is the per-cell phi = min_g phi_max(tau_g, Lambda) stable?"""
import numpy as np
C = 1e4; g1 = 2/3; sig = 1/np.sqrt(3)

def G(phis, taus, ws, nuc, P):
    n = len(taus); dt = nuc / C
    S = np.zeros((n+1, n+1))
    for g,(ph,t,w) in enumerate(zip(phis, taus, ws)):
        S[0,0] += -4*g1*P*C*t*w; S[0,g+1] += g1*P*C*t
        S[g+1,0] = 4*ph*C*t*w;  S[g+1,g+1] = -ph*C*t
    Q = dt*S; I = np.eye(n+1)
    Z = np.diag([0] + [-2*sig*ph*nuc for ph in phis])   # nu_g = phi_g*nuc
    E = I + Z
    G1 = np.linalg.solve(I - Q, E)
    return np.linalg.solve(I - Q/2, I/2 + E @ G1 / 2)

def phimax(tau, Lam, nuc):
    den = sig - tau/(2 + Lam*tau)
    return 1.0 if den <= 0 else min(1.0, 1/(nuc*den))

rng = np.random.default_rng(1)
worst = 0; nbad = 0
for trial in range(4000):
    NG = rng.integers(2, 6)
    taus = 10**rng.uniform(-3, 4, NG)
    ws = rng.dirichlet(np.ones(NG))
    nuc = 10**rng.uniform(1, 3.5); P = 10**rng.uniform(-6, 0)
    Lam = 4*g1*P*nuc
    for mode in ('cell', 'group'):
        if mode == 'cell':
            ph = min(phimax(t, Lam, nuc) for t in taus); phis = [ph]*NG
        else:
            phis = [phimax(t, Lam, nuc) for t in taus]
        rho = np.max(np.abs(np.linalg.eigvals(G(phis, taus, ws, nuc, P))))
        if mode == 'group': worst = max(worst, rho)
        if rho > 1 + 1e-9:
            nbad += 1
            if nbad < 8: print(mode, 'unstable', rho, np.round(np.log10(taus),2), np.round(ws,2), round(nuc), P)
print('max rho (per-group phi):', worst, ' unstable cases:', nbad)
