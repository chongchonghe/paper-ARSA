"""Two checks of the stability bound: (1) the closed form eq:numax against the numerical boundary of
pdars.py at Mach 0 and 3 in a worst case; (2) in 2D and 3D, the stage certificate D* B^-1 D <= B holds at the
bound (3000 random cases, fixed seed) and fails 1% above it at the all-pi checkerboard."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'src'))
from pdars import nu_max
rng = np.random.default_rng(3)

# --- 1. Mach 3 worst case reported by agent C, against pdars.py's own boundary search
def formula(tau, nuc, P, gamma, chi=1/3):
    K = 4 * (gamma - 1) * P * nuc * tau
    d = np.sqrt(chi) - tau / (2 + K)
    return min(nuc, np.inf if d <= 0 else 1 / d)
TH = np.pi * np.linspace(0.02, 1, 400)
for M in (0.0, 3.0):
    tau, gamma, P, nuc = 1.13, 1.19, 1e-6, 79.0
    n = nu_max(tau, nuc, P, M=M, gamma=gamma, thetas=TH, rel=1e-6)
    f = formula(tau, nuc, P, gamma)
    print(f"M={M}: numerical nu_max = {n:.6f}, closed form = {f:.6f}, relative (formula - true)/true = {(f - n) / n:+.2e}")

# --- 2. multi-D stage certificate D* B^-1 D <= B at mu_tot r_max = 1, own construction
def stage_min_eig(D, mu, q, K, qF, KF, th):
    # variables: matter (T, v_1..v_D), radiation (E, F_1..F_D); symmetrized exchange, rank one per channel
    n = 2 * (D + 1)
    Q = np.zeros((n, n))
    def chan(m, r, qc, Kc):
        Q[m, m] += Kc; Q[r, r] += qc; Q[m, r] -= np.sqrt(Kc * qc); Q[r, m] -= np.sqrt(Kc * qc)
    chan(0, D + 1, q, K)
    for d in range(D):
        chan(1 + d, D + 2 + d, qF, KF)
    B = np.eye(n) + Q
    Zr = np.zeros((D + 1, D + 1), dtype=complex)
    for d in range(D):
        X = np.zeros((D + 1, D + 1)); X[0, 1 + d] = X[1 + d, 0] = 1
        Zr += mu[d] * ((1 - np.cos(th[d])) * np.eye(D + 1) + 1j * np.sin(th[d]) * X)
    Dm = np.eye(n, dtype=complex); Dm[D + 1:, D + 1:] -= Zr
    Mx = B - Dm.conj().T @ np.linalg.solve(B, Dm)
    return np.linalg.eigvalsh((Mx + Mx.conj().T) / 2).min()

worst_at, viol_above = np.inf, 0
for trial in range(3000):
    D = rng.choice([2, 3])
    q, K = 10**rng.uniform(-3, 4), 10**rng.uniform(-4, 4)
    qF, KF = q * 10**rng.uniform(-2, 2), 10**rng.uniform(-8, 1)
    rmax = max((2 + K) / (2 + q + K), (2 + KF) / (2 + qF + KF))
    w = rng.dirichlet(np.ones(D)); mu = w / rmax                    # mu_tot r_max = 1 exactly
    th = rng.uniform(0, np.pi, D) if trial % 3 else np.pi * np.ones(D) - rng.uniform(0, 0.05, D)
    worst_at = min(worst_at, stage_min_eig(D, mu, q, K, qF, KF, th) / (1 + q + K + qF + KF))
    # above the bound, at the all-pi checkerboard it must fail
    viol_above += stage_min_eig(D, 1.01 * mu, q, K, qF, KF, np.pi * np.ones(D)) < 0
print(f"multi-D, at the bound: smallest normalized eigenvalue of B - D*B^-1D over 3000 draws = {worst_at:.2e}")
print(f"multi-D, 1% above the bound at the all-pi checkerboard: certificate fails in {viol_above}/3000")
