"""Independent check of the checkerboard theorem with the paper's own pdars.py matrices.
H is built from the formula in the agent's report (paper variables, density dropped)."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from pdars import amplification
rng = np.random.default_rng(7)
TH = np.pi * np.linspace(1e-3, 1, 400)

def a_star(q, K):
    return 2 * (2 + q + K) / (2 + K)

def setup(chi, gamma, nuc, P, C, tau):
    sig = np.sqrt(chi); g1 = gamma - 1
    # nu_max from a = 2 nu sigma <= min(a*_E, a*_F); both depend on nu through q = nu tau
    # energy: K = 4 g1 P nuc tau (no phi);  momentum: K_F = (1+chi) P nuc tau / C^2
    K, KF = 4 * g1 * P * nuc * tau, (1 + chi) * P * nuc * tau / C**2
    # 2 nu sig <= 2(2 + nu tau + K)/(2+K)  <=>  nu (sig - tau/(2+K)) <= 1
    def numax(Kc):
        d = sig - tau / (2 + Kc); return np.inf if d <= 0 else 1 / d
    return min(numax(K), numax(KF), nuc), K, KF

def H_of(nu, nuc, chi, gamma, P, C, tau):
    phi = nu / nuc; sig2 = chi; q = nu * tau
    w = np.array([0, (1 + chi) * phi / (sig2 * P), 4 * phi / ((gamma - 1) * P), 1.0, 1 / sig2])  # (rho,v,T,E,F)
    H = 2 * np.diag(w)[1:, 1:]
    uE = np.array([0, 4, -1, 0]); uF = np.array([(1 + chi) / C, 0, 0, -1])   # (v,T,E,F)
    H += q * np.outer(uE, uE) + (q / sig2) * np.outer(uF, uF)
    return H

worst_cert, worst_rho, n_in, fails_above = -np.inf, -np.inf, 0, 0
for trial in range(600):
    chi = rng.choice([1/3, 0.5, 1.0]); gamma = rng.choice([5/3, 1.4, 1.1])
    nuc = 10**rng.uniform(1, 3.5); P = 10**rng.uniform(-6, 1); C = 10**rng.uniform(2, 5)
    tau = 10**rng.uniform(-3, 4)
    nm, K, KF = setup(chi, gamma, nuc, P, C, tau)
    if not np.isfinite(nm) or nm >= nuc:
        continue
    for fac in (rng.uniform(0.3, 1.0), 1.0):
        nu = fac * nm
        H = H_of(nu, nuc, chi, gamma, P, C, tau)
        Hn = np.linalg.norm(H, 2)
        for th in TH:
            G = amplification(th, tau, nu, nuc, P, C=C, gamma=gamma, chi=chi)[1:, 1:]
            Mx = H - G.conj().T @ H @ G
            worst_cert = max(worst_cert, -np.linalg.eigvalsh((Mx + Mx.conj().T) / 2).min() / Hn)
            worst_rho = max(worst_rho, np.abs(np.linalg.eigvals(G)).max() - 1)
        n_in += 1
    # just above: rho(G(pi)) must exceed 1
    G = amplification(np.pi, tau, 1.001 * nm, nuc, P, C=C, gamma=gamma, chi=chi)
    fails_above += np.abs(np.linalg.eigvals(G)).max() > 1
print(f"cases at or below the bound: {n_in}")
print(f"worst normalized violation of H - G*HG >= 0 over theta: {worst_cert:.2e}")
print(f"worst rho(G(theta)) - 1 over theta: {worst_rho:.2e}")
print(f"cases unstable at theta = pi just above the bound (1.001 nu_max): {fails_above}")
