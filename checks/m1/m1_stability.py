"""Check the full-M1 linearized HLL/PD-ARS derivation.

This is a NumPy reproduction of the first-order flux, not a Quokka run.
Run: python3 checks/m1/m1_stability.py (writes data/m1/checks.json; runs from any directory)
"""

import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


def closure(u):
    p = u[1:] / u[0]  # F/(c E); c=1 in these variables
    f = np.linalg.norm(p)
    if f == 0:
        return np.eye(3) / 3
    n = p / f
    chi = (5 - 2 * np.sqrt(4 - 3 * f * f)) / 3
    return (1 - chi) / 2 * np.eye(3) + (3 * chi - 1) / 2 * np.outer(n, n)


def flux(u, d):
    return np.r_[u[d + 1], u[0] * closure(u)[:, d]]


def hll(left, right, d):
    sl = -max(0.1, np.sqrt(closure(left)[d, d]))
    sr = max(0.1, np.sqrt(closure(right)[d, d]))
    return (sr * flux(left, d) - sl * flux(right, d)
            + sr * sl * (right - left)) / (sr - sl)


def matrices(beta, n):
    n = np.asarray(n, dtype=float)
    n /= np.linalg.norm(n)
    f = 4 * beta / (3 + beta * beta)
    u = np.r_[1., f * n]
    t = closure(u)
    chi = (1 + 3 * beta * beta) / (3 + beta * beta)
    cp = 4 * beta / (3 - beta * beta)
    nn = np.outer(n, n)
    dt_df = -cp / 2 * np.eye(3) + 3 * cp / 2 * nn
    deriv = np.zeros((3, 3, 3))
    if f > 0:
        for j in range(3):
            tangent = np.eye(3)[j] - n[j] * n
            deriv[:, :, j] = (n[j] * dt_df + (3 * chi - 1) / (2 * f)
                              * (np.outer(tangent, n) + np.outer(n, tangent)))
    a = np.zeros((3, 4, 4))
    for d in range(3):
        a[d, 0, d + 1] = 1
        a[d, 1:, 0] = (t - f * dt_df)[:, d]
        a[d, 1:, 1:] = deriv[:, d, :]
    den = (1 - beta * beta) * (3 - beta * beta)
    h = np.zeros((4, 4))
    h[0, 0] = (1 + 5 * beta * beta) / den
    h[0, 1:] = h[1:, 0] = -beta * (5 + beta * beta) / den * n
    h[1:, 1:] = np.eye(3) + (3 * (1 + beta * beta) / den - 1) * nn
    speeds = np.maximum(0.1, np.sqrt(np.diag(t)))
    return u, a, h, speeds


def symbolic_checks():
    b, mu, lam, x = sp.symbols('b mu lam x', real=True)
    den = (1 - b*b) * (3 - b*b)
    cp = 4*b/(3-b*b)
    f = 4*b/(3+b*b)
    chi = (1+3*b*b)/(3+b*b)
    h = sp.Matrix([[(1+5*b*b)/den, -b*(5+b*b)/den, 0],
                   [-b*(5+b*b)/den, 3*(1+b*b)/den, 0], [0, 0, 1]])
    apar = sp.Matrix([[0, 1, 0], [chi-f*cp, cp, 0], [0, 0, b]])
    aperp = sp.Matrix([[0, 0, 1], [0, 0, b], [(1-chi+f*cp)/2, -cp/2, 0]])
    for a in [apar, aperp]:
        assert (h*a - a.T*h).applyfunc(sp.simplify) == sp.zeros(3)
    assert sp.simplify(h[:2, :2].det() - 1/(3-b*b)) == 0
    # A unit direction with parallel component mu; the fourth eigenvalue is b*mu.
    a = mu*apar + sp.sqrt(1-mu*mu)*aperp
    expected = ((lam-b*mu) * ((3-b*b)*lam**2 - 4*b*mu*lam
                            + 2*b*b*mu*mu - (1-b*b))) / (3-b*b)
    assert sp.simplify(a.charpoly(lam).as_expr().subs({'lam': lam}) - expected) == 0
    # The Cauchy-Schwarz speed bound has a nonnegative residual.
    residual = 3*(3-b*b)**2 - (3+b*b)*(4-b*b-2*b*b*mu*mu)
    assert sp.expand(residual - ((1-b*b)*(15-4*b*b)
                                + 2*b*b*mu*mu*(3+b*b))) == 0
    y2 = 3*x*(2-x)
    mod2 = (1-x+x*x/2-y2/2)**2 + y2*(1-x)**2
    assert sp.factor(mod2-1) == x*(x-2)*(x-1)**2
    return {'common_symmetrizer': True, 'directional_characteristic_polynomial': True,
            'speed_bound_residual': True, 'rk2_ellipse_identity': True}


def numerical_checks():
    rng = np.random.default_rng(20261023)
    stats = dict(states=0, fourier_modes=0, maximum_symmetry_residual=0.,
                 maximum_flux_jacobian_error=0., maximum_hll_jacobian_error=0.,
                 maximum_H_operator_norm=0., maximum_speed_ratio=0.)
    for beta in [0., .1, .5, .8, .9, .99, .999]:
        for n in [np.array([1., 0, 0]), np.ones(3), *rng.normal(size=(8, 3))]:
            u, a, h, speeds = matrices(beta, n)
            chol = np.linalg.cholesky(h).T  # H = chol^T chol
            ci = np.linalg.inv(chol)
            ah = np.array([chol @ ad @ ci for ad in a])
            stats['states'] += 1
            stats['maximum_symmetry_residual'] = max(
                stats['maximum_symmetry_residual'], float(np.max(abs(ah-ah.transpose(0, 2, 1)))))
            stats['maximum_speed_ratio'] = max(stats['maximum_speed_ratio'],
                max(float(np.max(abs(np.linalg.eigvalsh(ad)))) / sd for ad, sd in zip(ah, speeds)))
            # Central differences are independent of the analytic derivatives.
            # Stay away from the realizability boundary for the differencing check.
            if beta <= .9:
                eps = 1e-6
                for d in range(3):
                    fd = np.column_stack([(flux(u+eps*e, d)-flux(u-eps*e, d))/(2*eps)
                                          for e in np.eye(4)])
                    stats['maximum_flux_jacobian_error'] = max(
                        stats['maximum_flux_jacobian_error'], float(np.max(abs(fd-a[d]))))
                    for side in [0, 1]:
                        def perturbed(v):
                            return hll(v, u, d) if side == 0 else hll(u, v, d)
                        fd = np.column_stack([(perturbed(u+eps*e)-perturbed(u-eps*e))/(2*eps)
                                              for e in np.eye(4)])
                        predicted = (a[d] + (-1)**side * speeds[d]*np.eye(4))/2
                        stats['maximum_hll_jacobian_error'] = max(
                            stats['maximum_hll_jacobian_error'], float(np.max(abs(fd-predicted))))
            for dim in [1, 2, 3]:
                nu = 1 / np.sum(speeds[:dim])
                phases = np.vstack([np.zeros(dim), np.full(dim, np.pi),
                                    rng.uniform(-np.pi, np.pi, (100, dim))])
                for theta in phases:
                    xx = nu * np.dot(speeds[:dim], 1-np.cos(theta))
                    yy = nu * np.einsum('d,dij->ij', np.sin(theta), ah[:dim])
                    w = xx*np.eye(4) + 1j*yy
                    g = np.eye(4)-w+w@w/2
                    op = float(np.linalg.svd(g, compute_uv=False)[0])
                    stats['maximum_H_operator_norm'] = max(stats['maximum_H_operator_norm'], op)
                    stats['fourier_modes'] += 1
    assert stats['maximum_symmetry_residual'] < 1e-7
    assert stats['maximum_flux_jacobian_error'] < 1e-7
    assert stats['maximum_hll_jacobian_error'] < 1e-7
    assert stats['maximum_H_operator_norm'] < 1+1e-7
    return stats


def absorption_counterexample():
    mp.mp.dps = 70
    b, mu, q = mp.mpf(9)/10, mp.mpf(1)/6, mp.mpf(10)
    sigma = mp.sqrt((1-b*b+4*b*b*mu*mu)/(3+b*b))
    lam = (2*b*mu+mp.sqrt((1-b*b)*(3-b*b-2*b*b*mu*mu)))/(3-b*b)
    theta = 3*mp.pi/4
    # Strictly below the checkerboard bound, avoiding a marginal-only example.
    nu = mp.mpf('0.999')*(2+q)/(2*sigma)
    w = nu*(sigma*(1-mp.cos(theta))+1j*lam*mp.sin(theta))
    g = (1+(1-w)**2/(1+q))/(2+q)
    checker = (1+(1-2*nu*sigma)**2/(1+q))/(2+q)
    assert abs(g) > 1 and abs(checker) < 1
    u, a, h, speeds = matrices(float(b), [float(mu), float(mp.sqrt(1-mu*mu)), 0])
    wmat = float(nu)*(speeds[0]*(1-np.cos(float(theta)))*np.eye(4)
                     + 1j*np.sin(float(theta))*a[0])
    gmat = (np.eye(4)+(np.eye(4)-wmat)@(np.eye(4)-wmat)/(1+float(q)))/(2+float(q))
    rho = float(max(abs(np.linalg.eigvals(gmat))))
    assert abs(rho-float(abs(g))) < 1e-12
    # Finite material response, with equal grey means and v_bg=0.
    # Variables: (g, v_x/c, v_y/c, v_z/c, delta E, delta F/c).
    # g is the perturbation in the LTE emission energy; opacities are frozen.
    # Constant external sources cancel the nonzero background exchange only.
    ke = km = 1e-4
    exchange = np.zeros((8, 8))  # Q = -dt * source Jacobian
    energy_force = np.r_[-1., -u[1:], 1., np.zeros(3)]
    momentum_force = np.zeros((3, 8))
    momentum_force[:, 1:4] = -(np.eye(3)+closure(u))
    momentum_force[:, 5:8] = np.eye(3)
    exchange[0] = -ke*energy_force
    exchange[1:4] = -km*momentum_force
    exchange[4] = float(q)*energy_force
    exchange[5:8] = float(q)*momentum_force
    def coupled_growth(wr):
        d = np.eye(8, dtype=complex)
        d[4:, 4:] -= wr
        gfull = np.linalg.solve(2*np.eye(8)+exchange,
                              np.eye(8)+d@np.linalg.solve(np.eye(8)+exchange, d))
        return float(max(abs(np.linalg.eigvals(gfull))))
    coupled_rho = coupled_growth(wmat)
    coupled_checker = coupled_growth(2*float(nu)*speeds[0]*np.eye(4))
    assert coupled_rho > 1 and coupled_checker < 1
    # Real periodic eight-cell evolution, using centered flux plus HLL diffusion.
    vals, vecs = np.linalg.eig(a[0])
    v = vecs[:, np.argmax(vals)].real
    phase = np.arange(8)*float(theta)
    initial = np.cos(phase)[:, None]*v[None, :]
    state = initial.copy()
    def explicit(data):
        right, left = np.roll(data, -1, axis=0), np.roll(data, 1, axis=0)
        return data - float(nu)/2*((right-left)@a[0].T
                                  - speeds[0]*(right-2*data+left))
    for _ in range(1000):
        stage = explicit(state)/(1+float(q))
        state = (state + explicit(stage))/(2+float(q))
    measured = np.linalg.norm(state)/np.linalg.norm(initial)
    predicted = float(abs(g)**1000)
    assert abs(measured/predicted-1) < 1e-9
    return {**{key: mp.nstr(value, 35) for key, value in dict(
                beta=b, reduced_flux=4*b/(3+b*b), normal_cosine=mu,
                sigma=sigma, lambda_plus=lam, q=q, nu=nu, tau_cell=q/nu,
                theta=theta, amplification_modulus=abs(g), checkerboard=checker,
                gain_after_1000_steps=abs(g)**1000).items()},
            'matrix_spectral_radius': rho, 'real_grid_gain_after_1000_steps': measured,
            'finite_material_response': {'K_energy': ke, 'K_momentum_prefactor': km,
                                         'spectral_radius': coupled_rho,
                                         'checkerboard_spectral_radius': coupled_checker},
            'scope': 'Grey absorption into a fixed reservoir, with constant forcing maintaining the anisotropic background.'}


if __name__ == '__main__':
    results = {'symbolic': symbolic_checks(), 'numerical': numerical_checks(),
               'absorption_counterexample': absorption_counterexample()}
    dest = Path(__file__).resolve().parents[2]/'data/m1'
    dest.mkdir(parents=True, exist_ok=True)
    (dest/'checks.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))
