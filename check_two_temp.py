"""Two-temperature exchange under PD-ARS + HLL: which parameters set the bound?

Model (gray transport, non-gray mean opacities as in the paper): the emitter
has a temperature T_d = T_d(T, E), a function of the gas temperature T and of
the radiation energy E.  The net emission is N = rho c (kP a_r T_d^4 - kE E);
the radiation gains phi*N, the gas (heat capacity rho c_v, temperature T)
loses N.  Notation used in the printout:

  psi   = dT_d/dT at fixed E          (the parameter of the hypothesis)
  eta   = dT_d/dE at fixed T
  alpha = 1 - (kP/kE) d(a_r T_d^4)/dE|_T = 1 - 4 a_r T_d^3 (kP/kE) eta
  psiN  = dT_d/dT at fixed N  (= psi/alpha; the 'psi' of branch tmp-two-T)
  h_d   = 4 a_r T_d^3/(rho c_v)

Claimed closed form (D dimensions, static background):
  nu [ sum_d sigma_d - alpha tauE/(2 + K) ] <= 1,  K = h_d psi nuc tauP,
  nu [ sum_d sigma_d - tauF/(2 + K_F) ] <= 1      (momentum block unchanged).

Units of the numerical part are those of pdars.py: dx = 1, c = C, dt = nuc/C,
phi = nu/nuc, energy unit a_r T0^4 with T0 the background gas temperature
(so That = 1 in the background), h0 = 4 a_r T0^3/(rho c_v) = 4 (gamma-1) P.
The exchange rows are built from finite-difference derivatives of the
NONLINEAR net emission n(That, Ehat) = r That_d^4 - Ehat (r = kP/kE), so the
numerical boundary does not use the closed-form algebra.

Run:  ~/softwares/python-envs/base/bin/python check_two_temp.py
"""
import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdars import amplification as amp_pdars  # noqa: E402

C = 1e4
GAMMA = 5 / 3
G1 = GAMMA - 1


# ----------------------------------------------------------------------------
# 1. Symbolic checks
# ----------------------------------------------------------------------------
def symbolic():
    print("=" * 78)
    print("1. SYMBOLIC CHECKS (sympy)")
    print("=" * 78)
    T, E, rho, c, kP, kE, ar, cv, phi, dt, dx = sp.symbols(
        "T E rho c kappa_P kappa_E a_r c_v phi Delta_t Delta_x", positive=True)
    Td = sp.Function("T_d")(T, E)
    N = rho * c * (kP * ar * Td**4 - kE * E)
    NT, NE = sp.diff(N, T), -sp.diff(N, E)
    S = dt * sp.Matrix([[-NT / (rho * cv), NE / (rho * cv)],
                        [phi * NT, -phi * NE]])
    print("(1a) exchange over one step on (dT, dE): det S_Delta =", sp.simplify(S.det()), " (rank one)")
    psi, eta = sp.diff(Td, T), sp.diff(Td, E)
    tauE, tauP = rho * kE * dx, rho * kP * dx
    nu, nuc = phi * c * dt / dx, c * dt / dx
    hd = 4 * ar * Td**3 / (rho * cv)
    alpha = 1 - 4 * ar * Td**3 * (kP / kE) * eta
    K = hd * psi * nuc * tauP
    q = nu * tauE * alpha
    print("     S11 + K           =", sp.simplify(S[0, 0] + K), "  with K = h_d psi nuc tauP")
    print("     S22 + q           =", sp.simplify(S[1, 1] + q), "  with q = nu alpha tauE")
    print("     S12*S21 - S11*S22 =", sp.simplify(S[0, 1] * S[1, 0] - S[0, 0] * S[1, 1]))

    # (1b) energy block of PD-ARS for a general rank-one S with diagonal (-K, -q)
    a, Ks, qs, x = sp.symbols("a K q x", positive=True)
    S2 = sp.Matrix([[-Ks, x], [Ks * qs / x, -qs]])
    I2 = sp.eye(2)
    D = sp.diag(1, 1 - a)
    G = (I2 - S2 / 2).inv() * (I2 / 2 + D * (I2 - S2).inv() * D / 2)
    trG, detG = sp.simplify(G.trace()), sp.simplify(G.det())
    print("(1b) PD-ARS energy block with S = [[-K, x], [K q/x, -q]], D = diag(1, 1-a):")
    print("     x appears in tr G:", x in trG.free_symbols, "; in det G:", x in detG.free_symbols,
          "  -> eigenvalues depend on (a, K, q) only")
    astar = 2 * (2 + qs + Ks) / (2 + Ks)
    detB, detH = 1 + qs + Ks, 1 + (qs + Ks) / 2
    lhs = sp.simplify((I2 - G).det())
    rhs = Ks * (2 + Ks) * a * (astar - a) / (4 * detB * detH)
    print("     det(I - G) - K(2+K) a (a* - a)/(4 detB detH) =", sp.simplify(lhs - rhs),
          "  with a* = 2(2+q+K)/(2+K)")
    lam = Ks / qs  # the paper's lambda, so that lambda q = K
    Npar = (2 + Ks)**2 * a**2 - 2 * (lam**2 * qs**2 + lam * qs**2 + 4 * lam * qs + 4) * a + 4 * (2 + qs + Ks)**2
    print("     det(I + G) - N(a)/(4 detB detH) (paper eq:detGplusI) =",
          sp.simplify((I2 + G).det() - Npar / (4 * detB * detH)))

    # (1c) dust balance rho c (kP a_r T_d^4 - kE E) = Gamma(T, T_d)
    Tds = sp.symbols("T_d", positive=True)
    Gam = sp.Function("Gamma")(T, Tds)
    F = rho * c * (kP * ar * Tds**4 - kE * E) - Gam
    psi_b = sp.simplify(-sp.diff(F, T) / sp.diff(F, Tds))
    eta_b = sp.simplify(-sp.diff(F, E) / sp.diff(F, Tds))
    GT, Gd = sp.diff(Gam, T), -sp.diff(Gam, Tds)
    R = 4 * ar * Tds**3 * rho * kP * c
    alpha_old, psi_old = Gd / (Gd + R), GT / Gd
    alpha_b = 1 - 4 * ar * Tds**3 * (kP / kE) * eta_b
    print("(1c) dust balance (branch tmp-two-T): alpha_old = Gd/(Gd+R), psi_old = GT/Gd, R = 4 a_r T_d^3 rho kP c")
    print("     dT_d/dT|_E - alpha_old*psi_old  =", sp.simplify(psi_b - alpha_old * psi_old))
    print("     alpha(general) - alpha_old       =", sp.simplify(alpha_b - alpha_old))
    print("     dT_d/dE|_T - rho c kE/(Gd + R)   =", sp.simplify(eta_b - rho * c * kE / (Gd + R)))
    print("     -> tmp-two-T substitutions = general ones: alpha tauE and h_d psi_old alpha = h_d psi")

    # (1d) dT_d/dT at fixed N equals psi/alpha for any T_d(T, E)
    psiN = psi + eta * NT / NE
    print("(1d) dT_d/dT|_N - psi/alpha =", sp.simplify(psiN - psi / alpha), "  (any T_d(T, E))")

    print()


# ----------------------------------------------------------------------------
# 2. Emitter-temperature models T_d(That, Ehat) and the exchange Jacobian
# ----------------------------------------------------------------------------
class PowerLaw:
    """T_d = t0 * That^p * (Ehat/E0)^s : psi = p t0, eta = s t0/E0."""

    def __init__(self, t0, p, s, E0):
        self.t0, self.p, self.s, self.E0 = t0, p, s, E0

    def Td(self, T, E):
        return self.t0 * T**self.p * (E / self.E0)**self.s


class Balance:
    """Dust with no heat capacity: r T_d^4 - E = g T^m (T - T_d) + H, with
    H a constant extra dust heating (does not enter the linearization)."""

    def __init__(self, r, g, m, Td0, E0):
        self.r, self.g, self.m = r, g, m
        self.E0 = E0
        self.H = r * Td0**4 - E0 - g * (1 - Td0)  # background balance at That = 1
        self.Td0 = Td0

    def Td(self, T, E):
        f = lambda x: self.r * x**4 - E - self.g * T**self.m * (T - x) - self.H
        return brentq(f, 1e-6, 1e3, xtol=1e-15, rtol=1e-15)

    def old_params(self):
        """alpha and psi_old of branch tmp-two-T, from Gamma derivatives."""
        Td, g, m = self.Td0, self.g, self.m
        Gd = g
        GT = g * (m * (1 - Td) + 1)
        R = 4 * self.r * Td**3
        return Gd / (Gd + R), GT / Gd


def derivs(model, r, E0, hrel=1e-6):
    """Finite-difference derivatives at the background (That = 1, Ehat = E0).
    Returns T_d, psi, eta and the derivatives nT, nE of n = r T_d^4 - E."""
    Td = model.Td(1.0, E0)
    dT, dE = hrel, hrel * E0
    n = lambda T, E: r * model.Td(T, E)**4 - E
    psi = (model.Td(1 + dT, E0) - model.Td(1 - dT, E0)) / (2 * dT)
    eta = (model.Td(1.0, E0 + dE) - model.Td(1.0, E0 - dE)) / (2 * dE)
    nT = (n(1 + dT, E0) - n(1 - dT, E0)) / (2 * dT)
    nE = (n(1.0, E0 + dE) - n(1.0, E0 - dE)) / (2 * dE)
    return dict(Td=Td, psi=psi, eta=eta, nT=nT, nE=nE, alpha=1 - 4 * r * Td**3 * eta)


# ----------------------------------------------------------------------------
# 3. Amplification matrix of PD-ARS + HLL in D dimensions, M = 0
#    variables [T, v_1..v_D, E, F_1..F_D]
# ----------------------------------------------------------------------------
def build(D, nu, nuc, tauE, tauF, h0, nT, nE, Ebg, chi):
    n = 2 * D + 2
    iT, iE = 0, D + 1
    phi = nu / nuc
    dt = nuc / C
    P = h0 / (4 * G1)
    Q = np.zeros((n, n))
    # energy channel: gas loses N, radiation gains phi*N (per unit tauE, rho kE c)
    Q[iT, iT] = -(h0 / 4) * nuc * tauE * nT
    Q[iT, iE] = -(h0 / 4) * nuc * tauE * nE
    Q[iE, iT] = nu * tauE * nT
    Q[iE, iE] = nu * tauE * nE
    # momentum channels (pdars.matrices at M = 0, drag scaled by the background E)
    for d in range(D):
        iv, iF = 1 + d, D + 2 + d
        Q[iv, iv] = -dt * tauF * (1 + chi) * P * Ebg / C
        Q[iv, iF] = dt * tauF * P
        Q[iF, iv] = dt * tauF * (1 + chi) * phi * Ebg
        Q[iF, iF] = -nu * tauF
    return Q


def amplification(Q, D, nu, chi, thetas):
    n = Q.shape[0]
    iE = D + 1
    sig = np.sqrt(chi)
    Z = np.zeros((n, n), dtype=complex)
    for d, th in enumerate(thetas):
        iF = D + 2 + d
        X = np.zeros((n, n))
        X[iE, iF], X[iF, iE] = 1.0, chi
        Ir = np.zeros((n, n))
        Ir[iE:, iE:] = np.eye(D + 1)
        Z += -nu * (1j * np.sin(th) * X + sig * (1 - np.cos(th)) * Ir)
    I = np.eye(n)
    Em = I + Z
    G1s = np.linalg.solve(I - Q, Em)
    return np.linalg.solve(I - 0.5 * Q, 0.5 * I + 0.5 * Em @ G1s)


def rho_G(case, nu, thetas):
    Q = build(case["D"], nu, case["nuc"], case["tauE"], case["tauF"], case["h0"],
              case["nT"], case["nE"], case["Ebg"], case["chi"])
    return np.max(np.abs(np.linalg.eigvals(amplification(Q, case["D"], nu, case["chi"], thetas))))


def nu_checkerboard(case, rel=1e-10):
    """Largest nu <= nuc with rho(G(pi,...,pi)) <= 1, by bisection."""
    th = [np.pi] * case["D"]
    ok = lambda nu: rho_G(case, nu, th) <= 1 + 1e-12
    if ok(case["nuc"]):
        return case["nuc"]
    lo, hi = 1e-6, case["nuc"]
    assert ok(lo)
    while hi - lo > rel * hi:
        mid = np.sqrt(lo * hi) if hi / lo > 4 else 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


THETA_1D = np.concatenate([np.pi * np.geomspace(1e-4, 4e-3, 20), np.pi * np.linspace(5e-3, 1, 400)])
_g = np.pi * np.linspace(0.02, 1, 30)
THETA_2D = [(a, b) for a in _g for b in _g]


def theta_tests(case, nupi):
    """(i) below nu_pi no theta is unstable; (ii) just above, theta = pi grows fastest."""
    D = case["D"]
    grid = [(t,) for t in THETA_1D] if D == 1 else THETA_2D
    below = max(rho_G(case, (1 - 1e-6) * nupi, th) for th in grid)
    ok_below = below <= 1 + 1e-10
    ok_above = True
    if nupi < case["nuc"]:
        vals = [rho_G(case, (1 + 1e-3) * nupi, th) for th in grid]
        k = int(np.argmax(vals))
        ok_above = vals[k] > 1 and np.allclose(grid[k], np.pi)
    return ok_below, ok_above, below


def closed_form(case, alpha, psi, Td):
    """The two-temperature bound (energy and momentum conditions)."""
    sig = case["D"] * np.sqrt(case["chi"])
    hd = case["h0"] * Td**3
    K = hd * psi * case["nuc"] * case["r"] * case["tauE"]
    P = case["h0"] / (4 * G1)
    KF = (1 + case["chi"]) * P * case["Ebg"] * case["nuc"] * case["tauF"] / C**2
    out = case["nuc"]
    for den in (sig - alpha * case["tauE"] / (2 + K), sig - case["tauF"] / (2 + KF)):
        if den > 0:
            out = min(out, 1 / den)
    return out


def make_case(model, r, E0, D=1, nuc=1000.0, tauE=1.0, tauF=None, h0=1e-3, chi=1 / 3):
    dv = derivs(model, r, E0)
    case = dict(D=D, nuc=nuc, tauE=tauE, tauF=tauE if tauF is None else tauF, h0=h0, chi=chi,
                r=r, Ebg=E0, nT=dv["nT"], nE=dv["nE"])
    return case, dv


# ----------------------------------------------------------------------------
# 4. Numerical checks
# ----------------------------------------------------------------------------
def check_against_pdars():
    print("(2a) builder vs pdars.amplification, one temperature (T_d = T, r = 1, E0 = 1):")
    worst = 0.0
    rng = np.random.default_rng(0)
    for _ in range(50):
        tau, nuc, P = 10**rng.uniform(-2, 4), 10**rng.uniform(1, 3.5), 10**rng.uniform(-6, 1)
        chi = rng.choice([1 / 3, 1.0])
        nu = nuc * 10**rng.uniform(-3, 0)
        th = rng.uniform(0, np.pi)
        case, _ = make_case(PowerLaw(1.0, 1.0, 0.0, 1.0), 1.0, 1.0, nuc=nuc, tauE=tau, h0=4 * G1 * P, chi=chi)
        mine = np.linalg.eigvals(amplification(
            build(1, nu, nuc, tau, tau, 4 * G1 * P, case["nT"], case["nE"], 1.0, chi), 1, nu, chi, [th]))
        ref = np.linalg.eigvals(amp_pdars(th, tau, nu, nuc, P, C=C, gamma=GAMMA, chi=chi))
        # every eigenvalue of the 4x4 system is one of the 5x5 (which adds the density, g = 1)
        rel = max(np.min(np.abs(ref - g)) / max(1.0, abs(g)) for g in mine)
        worst = max(worst, rel)
    print(f"     largest relative eigenvalue difference over 50 random cases: {worst:.1e}")


def same_psi_pairs():
    print("(2b) SAME psi = dT_d/dT|_E, DIFFERENT eta = dT_d/dE|_T (power-law T_d, identical background)")
    print("     T_d/T = 0.5, r = kP/kE = 1, E0 = r T_d^4, nuc = 1000, chi = 1/3, tauF = tauE, psi = 0.3")
    print("     nu_num: largest nu with rho(G) <= 1 at theta = pi; theta-tests: (i) stable below for all theta,"
          " (ii) theta = pi grows fastest just above")
    print("     nu_2T: closed form with alpha;  nu_psi: psi-only form (alpha = 1);  phi = nu/nuc")
    t0, r = 0.5, 1.0
    E0 = r * t0**4
    regimes = [("K << 2", 1.0, 1e-5), ("K ~ 10", 5.0, 0.05), ("K >> 2", 1e4, 0.1)]
    psi = 0.3
    print(f"     {'regime':8s} {'tauE':>7s} {'h0':>7s} {'dlnTd/dlnE':>10s} {'alpha':>6s} {'K':>9s} "
          f"{'nu_num':>9s} {'nu_2T':>9s} {'nu_psi':>9s} {'thetaOK':>7s}")
    for name, tauE, h0 in regimes:
        for alpha_t in (1.0, 0.5, 0.1):
            s = (1 - alpha_t) * E0 / (4 * r * t0**4)
            model = PowerLaw(t0, psi / t0, s, E0)
            case, dv = make_case(model, r, E0, tauE=tauE, h0=h0)
            nupi = nu_checkerboard(case)
            okb, oka, _ = theta_tests(case, nupi)
            f2 = closed_form(case, dv["alpha"], dv["psi"], dv["Td"])
            fpsi = closed_form(case, 1.0, dv["psi"], dv["Td"])
            K = h0 * dv["Td"]**3 * dv["psi"] * case["nuc"] * r * tauE
            print(f"     {name:8s} {tauE:7.0e} {h0:7.0e} {dv['eta'] * E0 / dv['Td']:10.3f} {dv['alpha']:6.3f} {K:9.3g} "
                  f"{nupi:9.5g} {f2:9.5g} {fpsi:9.5g} {str(okb and oka):>7s}")


def limits():
    print("(2c) LIMITS: which combination matters (same setup as 2b)")
    t0, r = 0.5, 1.0
    E0 = r * t0**4

    def run(tauE, h0, psi, alpha_t):
        s = (1 - alpha_t) * E0 / (4 * r * t0**4)
        case, dv = make_case(PowerLaw(t0, psi / t0, s, E0), r, E0, tauE=tauE, h0=h0)
        return nu_checkerboard(case), closed_form(case, dv["alpha"], dv["psi"], dv["Td"])

    print("     K << 2 (tauE = 1, h0 = 1e-5, K <= 2.5e-3): vary psi at fixed alpha = 0.3 -> nu_num")
    for psi in (0.05, 0.3, 1.0, 2.0):
        n, f = run(1.0, 1e-5, psi, 0.3)
        print(f"        psi = {psi:4.2f}: nu_num = {n:.6g}  (closed form {f:.6g})")
    print("     K >> 2 (tauE = 1e4, h0 = 0.1): fixed psiN = psi/alpha = 0.6, vary the split -> nu_num")
    for psi, al in ((0.6, 1.0), (0.3, 0.5), (0.06, 0.1), (0.006, 0.01)):
        n, f = run(1e4, 0.1, psi, al)
        print(f"        psi = {psi:5.3f}, alpha = {al:4.2f}: nu_num = {n:.6g}  (closed form {f:.6g})")
    print("     thin (tauE = 1e-3, h0 = 0.1): any psi, alpha")
    for psi, al in ((1.0, 1.0), (0.3, 0.1)):
        n, f = run(1e-3, 0.1, psi, al)
        print(f"        psi = {psi:4.2f}, alpha = {al:4.2f}: nu_num = {n:.6g}  (closed form {f:.6g}; CFL 1/sigma = {np.sqrt(3):.6g})")


def balance_models():
    print("(2d) DUST BALANCE r T_d^4 - E = g T^m (T - T_d) + H (zero dust heat capacity)")
    r, Td0 = 1.0, 0.5
    R = 4 * r * Td0**3
    psi_t = 0.15
    print(f"     T_d/T = {Td0}, r = {r}; both models tuned to psi = dT_d/dT|_E = {psi_t}")
    print(f"     {'model':22s} {'psi':>7s} {'alpha':>7s} {'alpha_old':>9s} {'psi_old':>7s} "
          f"{'tauE':>6s} {'h0':>5s} {'nu_num':>9s} {'nu_2T':>9s} {'nu_old':>9s} {'nu_psi':>9s} {'nu_psiN':>9s} {'nu_1par':>9s}")
    for m in (0.0, 0.5):
        psi_old = 1 + m * (1 - Td0)
        al = psi_t / psi_old
        g = al * R / (1 - al)
        E0 = r * Td0**4 - g * (1 - Td0)  # H = 0: the gas heats the dust, which radiates
        model = Balance(r, g, m, Td0, E0)
        a_old, p_old = model.old_params()
        for tauE, h0 in ((5.0, 1e-3), (1e4, 0.02)):
            case, dv = make_case(model, r, E0, tauE=tauE, h0=h0)
            nupi = nu_checkerboard(case)
            okb, oka, _ = theta_tests(case, nupi)
            f2 = closed_form(case, dv["alpha"], dv["psi"], dv["Td"])
            # tmp-two-T: tauE -> alpha tauE, K = h nuc alpha tauP with h = h_d psi_old
            fold = closed_form(case, a_old, a_old * p_old, dv["Td"])
            fpsi = closed_form(case, 1.0, dv["psi"], dv["Td"])        # psi only, alpha = 1
            fpsiN = closed_form(case, 1.0, p_old, dv["Td"])           # psiN = GT/Gd only, alpha = 1
            f1p = closed_form(case, dv["psi"], dv["psi"], dv["Td"])   # psi only, alpha = psi
            lab = f"m={m} (Gamma ~ T^{m}(T-Td))"
            print(f"     {lab:22s} {dv['psi']:7.4f} {dv['alpha']:7.4f} {a_old:9.4f} {p_old:7.3f} "
                  f"{tauE:6.0e} {h0:5.0e} {nupi:9.5g} {f2:9.5g} {fold:9.5g} {fpsi:9.5g} {fpsiN:9.5g} {f1p:9.5g}"
                  + ("" if okb and oka else "  THETA-TEST FAILED"))
    print("     nu_old: tmp-two-T substitutions with (alpha_old, psi_old) from Gamma; nu_psi: alpha = 1, K = h_d psi nuc tauP;")
    print("     nu_psiN: alpha = 1, K = h_d psiN nuc tauP (psiN = dT_d/dT|_N = GT/Gd as the only parameter);"
          " nu_1par: tauE -> psi tauE, K = h_d psi nuc tauP")

    print("     equilibrium background (T_d = T, E = r T^4, H = 0), m = 1/2: is alpha = psi?")
    for g in (1e-3, 0.1, 10.0):
        model = Balance(r, g, 0.5, 1.0, r * 1.0)
        dv = derivs(model, r, r * 1.0)
        print(f"        g = {g:6.0e}: psi = {dv['psi']:.8f}, alpha = {dv['alpha']:.8f}")
    print("     Keto-Field form (m = 1/2) away from equilibrium: psiN = GT/Gd = 3/2 - T_d/(2T) =",
          f"{1.5 - Td0 / 2:.3f} at T_d/T = {Td0}")


def random_scan(n1=400, n2=40, seed=7):
    print(f"(2e) RANDOM SCAN: {n1} cases in 1D, {n2} in 2D, power-law and balance T_d models")
    print("     ranges: tauE 1e-3..1e5, tauF/tauE 0.3..3, kP/kE 0.1..10, h0 1e-6..10, nuc 10..3e3,"
          " T_d/T 0.1..2, psi 0.01..2 (power law), alpha 0.01..1 (power law),")
    print("     balance: Gamma = g T^m (T - T_d), m in {0, 1/2, 1}, g/(4 r T_d^3) 1e-3..1e2, E/(r T_d^4) 0.1..10;"
          " chi in {1/3, 1} (1D), {1/3, 1/2} (2D)")
    rng = np.random.default_rng(seed)
    worst = {1: 0.0, 2: 0.0}
    worst_psi = 0.0
    nfail = 0
    nunsafe_psi = 0
    nactive = {1: 0, 2: 0}
    ranges = dict(psi=[np.inf, 0], alpha=[np.inf, 0], K=[np.inf, 0])
    for k in range(n1 + n2):
        D = 1 if k < n1 else 2
        r = 10**rng.uniform(-1, 1)
        if rng.uniform() < 0.5:
            t0 = 10**rng.uniform(-1, 0.3)
            E0 = r * t0**4 * 10**rng.uniform(-1, 1)
            psi = 10**rng.uniform(-2, 0.3)
            al = rng.uniform(0.01, 1.0)
            s = (1 - al) * E0 / (4 * r * t0**4)
            model = PowerLaw(t0, psi / t0, s, E0)
        else:
            Td0 = 10**rng.uniform(-1, 0.3)
            m = rng.choice([0.0, 0.5, 1.0])
            g = 10**rng.uniform(-3, 2) * 4 * r * Td0**3
            E0 = r * Td0**4 * 10**rng.uniform(-1, 1)
            model = Balance(r, g, m, Td0, E0)
        nuc = 10**rng.uniform(1, 3.5)
        tauE = 10**rng.uniform(-3, 5)
        tauF = tauE * 10**rng.uniform(-0.5, 0.5)
        h0 = 10**rng.uniform(-6, 1)
        chi = rng.choice([1 / 3, 1.0]) if D == 1 else rng.choice([1 / 3, 1 / 2])
        case, dv = make_case(model, r, E0, D=D, nuc=nuc, tauE=tauE, tauF=tauF, h0=h0, chi=chi)
        if dv["psi"] < 0 or dv["alpha"] <= 0:
            continue
        nupi = nu_checkerboard(case)
        okb, oka, _ = theta_tests(case, nupi)
        nfail += not (okb and oka)
        f2 = closed_form(case, dv["alpha"], dv["psi"], dv["Td"])
        fpsi = closed_form(case, 1.0, dv["psi"], dv["Td"])
        err = abs(nupi - f2) / f2
        worst[D] = max(worst[D], err)
        if nupi < nuc:
            nactive[D] += 1
        K = h0 * dv["Td"]**3 * dv["psi"] * nuc * r * tauE
        for key, v in (("psi", dv["psi"]), ("alpha", dv["alpha"]), ("K", K)):
            ranges[key] = [min(ranges[key][0], v), max(ranges[key][1], v)]
        if fpsi > nupi * (1 + 1e-6):
            nunsafe_psi += 1
            worst_psi = max(worst_psi, fpsi / nupi)
    print("     realized: " + ", ".join(f"{k} {v[0]:.2g}..{v[1]:.2g}" for k, v in ranges.items()))
    print(f"     worst |nu_num - nu_2T|/nu_2T: 1D {worst[1]:.1e}, 2D {worst[2]:.1e};  theta-test failures: {nfail};"
          f"  active bound (nu < nuc): 1D {nactive[1]}, 2D {nactive[2]}")
    print(f"     psi-only form (alpha = 1) above the true boundary in {nunsafe_psi} cases;"
          f" largest overshoot nu_psi/nu_num = {worst_psi:.3g}")


def outside_dissipative():
    print("(2f) OUTSIDE THE DISSIPATIVE RANGE (alpha <= 0 or psi < 0)")
    print("     rho(G) at nu = 1e-3 and theta = 0 tests the exchange alone (no transport)")
    t0, r = 0.5, 1.0
    E0 = r * t0**4
    for al, psi, tauE, h0 in ((0.0, 0.3, 10.0, 1e-3), (-0.2, 0.3, 10.0, 1e-3), (-0.2, 0.3, 1e3, 1.0),
                              (0.5, -0.05, 10.0, 1e-2), (0.5, -0.3, 10.0, 0.1)):
        s = (1 - al) * E0 / (4 * r * t0**4)
        case, dv = make_case(PowerLaw(t0, psi / t0, s, E0), r, E0, tauE=tauE, h0=h0)
        K = h0 * dv["Td"]**3 * dv["psi"] * case["nuc"] * r * tauE
        rho0 = rho_G(case, 1e-3, [0.0])
        nupi = nu_checkerboard(case) if rho_G(case, 1e-3, [np.pi]) <= 1 + 1e-12 else float("nan")
        f2 = closed_form(case, dv["alpha"], dv["psi"], dv["Td"])
        alp = 0.0 if abs(dv["alpha"]) < 1e-9 else dv["alpha"]
        print(f"     alpha = {alp:+.2f}, psi = {dv['psi']:+.2f}, tauE = {tauE:g}, h0 = {h0:g}, K = {K:+.3g}:"
              f" rho(G(0)) = {rho0:.6f}; nu_num = {nupi:.5g}, formula = {f2:.5g}")


if __name__ == "__main__":
    symbolic()
    print("=" * 78)
    print("2. NUMERICAL CHECKS: PD-ARS + HLL amplification matrix, M = 0, C = 1e4, gamma = 5/3")
    print("=" * 78)
    check_against_pdars()
    same_psi_pairs()
    limits()
    balance_models()
    random_scan()
    outside_dissipative()
