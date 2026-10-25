"""All-mode M1 bounds with exchange; run from the workspace root.

Analytic certificates, numerical spot checks, and a continuous-system
obstruction to unconditional all-time nonamplification are kept separate.
"""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
import sympy as sp

from m1_stability import matrices, closure


ALPHA = np.sqrt(2) / 3


def ellipse_maximum(q, L, R):
    """Exact candidate set in real arithmetic; floating-point polynomial roots.

    Returns max |g|^2 on y^2=R*x*(2*L-x), including both endpoints.
    This numerical routine is not an interval-arithmetic certificate.
    """
    c = Polynomial([0., 1.])
    x = L * (1-c)
    v = R * L*L * (1-c*c)
    numerator = (1+q+(1-x)**2-v)**2 + 4*(1-x)**2*v
    candidates = [-1., 1.]
    for root in numerator.deriv().roots():
        if abs(root.imag) < 1e-8 and -1 < root.real < 1:
            candidates.append(float(root.real))
    vals = [float(numerator(t)/((1+q)*(2+q))**2) for t in candidates]
    j = int(np.argmax(vals))
    return vals[j], candidates[j]


def estimate_ellipse_limit(q, R):
    """Boundary estimate, with the returned point checked by ellipse_maximum."""
    low, high = 1., 1+q/2
    if ellipse_maximum(q, high, R)[0] <= 1+1e-13:
        return high
    for _ in range(60):
        middle = (low+high)/2
        if ellipse_maximum(q, middle, R)[0] <= 1:
            low = middle
        else:
            high = middle
    return low


def bernstein_coefficients(poly, variable, lo, hi):
    t = sp.Symbol('t', real=True)
    p = sp.Poly(sp.expand(poly.subs(variable, lo+(hi-lo)*t)), t)
    return [sp.simplify(sum(p.nth(k)*sp.binomial(j, k)/sp.binomial(4, k)
                            for k in range(j+1))) for j in range(5)]


def symbolic_certificate():
    x, v, q = sp.symbols('x v q', real=True)
    alpha = sp.sqrt(2)/3
    # |B(z)|^2, B(z)=1+2*alpha*(z^2-z), is convex in v=(Im z)^2.
    B2 = (1+2*alpha*(x*x-x-v))**2 + 4*alpha**2*v*(2*x-1)**2
    assert sp.diff(B2, v, 2) == 8*alpha**2
    upper = sp.expand(9-B2.subs(v, 3*x*(2-x)))
    coeffs = [bernstein_coefficients(upper, x, lo, hi) for lo, hi in [(0, 1), (1, 2)]]
    assert all(c.is_positive for row in coeffs for c in row)
    # At v=0, B is real and in [1-alpha/2, 1+4*alpha], a subset of (0,3).
    assert (3-(1+4*alpha)).is_positive
    # |z|^2 <= x^2+3*x*(2-x) <= 9/2, so alpha^2 |z|^2 <= 1.
    assert sp.expand(sp.Rational(9, 2)-(x*x+3*x*(2-x))-2*(x-sp.Rational(3, 2))**2) == 0
    # Numerator decomposition used in the all-q proof.
    z = sp.Symbol('z')
    Rz = 1-z+z*z/2
    assert sp.expand(1+q+(1-(1+alpha*q)*z)**2
                     -(2*Rz+q*(1+2*alpha*(z*z-z))+alpha**2*q*q*z*z)) == 0
    t = sp.Symbol('t', real=True)
    correction = 1/(2*(1-t/2)) + 3/(2*(1-t/2)*(1-t))
    assert sp.simplify(correction-(4-t)/((2-t)*(1-t))) == 0
    return {'bernstein_coefficients': [[str(c) for c in row] for row in coeffs],
            'scalar_bound_symbolic_certificate': True,
            'full_cfl_finite_time_constant': True}


def source_matrix(u, q, ke, km):
    C = np.eye(4)
    C[0, 1:] = u[1:]/u[0]
    C[1:, 1:] += closure(u)
    D = np.diag([ke, km, km, km])
    return np.block([[D@C, -D], [-q*C, q*np.eye(4)]])


def pdars(Q, W):
    eye = np.eye(len(Q))
    D = eye-W
    return np.linalg.solve(2*eye+Q, eye+D@np.linalg.solve(eye+Q, D))


def scalar_checks():
    rng = np.random.default_rng(2401)
    worst = 0.
    count = 0
    rows = []
    for q in [0., .01, .1, .5, 1., 2., 10., 100., 10000.]:
        L = 1+ALPHA*q
        max2, phase_cos = ellipse_maximum(q, L, 3.)
        worst = max(worst, np.sqrt(max2))
        assert max2 < 1+1e-12
        rows.append({'q': q, 'closed_bound_L': L, 'checkerboard_L': 1+q/2,
                     'universal_ellipse_limit_estimate': estimate_ellipse_limit(q, 3.),
                     'maximum_modulus_at_closed_bound': float(np.sqrt(max2)),
                     'maximizing_cosine': phase_cos})
        for _ in range(15):
            beta = rng.uniform(0, .995)
            u, a, h, speeds = matrices(beta, rng.normal(size=3))
            chol = np.linalg.cholesky(h).T
            ci = np.linalg.inv(chol)
            ah = np.array([chol@ad@ci for ad in a])
            for dim in [1, 2, 3]:
                nu = L/sum(speeds[:dim])
                for theta in rng.uniform(-np.pi, np.pi, (20, dim)):
                    x = nu*np.dot(speeds[:dim], 1-np.cos(theta))
                    y = nu*np.einsum('d,dij->ij', np.sin(theta), ah[:dim])
                    G = pdars(q*np.eye(4), x*np.eye(4)+1j*y)
                    op = float(np.linalg.svd(G, compute_uv=False)[0])
                    assert op < 1+1e-9
                    worst = max(worst, op)
                    count += 1
    u, a, h, speeds = matrices(.9, [1/6, np.sqrt(35)/6, 0])
    R = float(max(abs(np.linalg.eigvals(a[0])))**2/speeds[0]**2)
    actual_max2, _ = ellipse_maximum(10., .999*6, R)
    assert actual_max2 > 1
    state_bound = estimate_ellipse_limit(10., R)
    assert ellipse_maximum(10., state_bound, R)[0] < 1+1e-12
    return {'matrix_mode_checks': count, 'largest_checked_H_norm': worst, 'table': rows,
            'previous_counterexample': {'R': R, 'new_closed_bound_L': 1+10*ALPHA,
                                        'state_specific_limit_L_estimate': state_bound,
                                        'old_example_L': 5.994}}


def continuous_obstruction():
    # A two-variable parallel radiation block coupled to gas thermal/momentum response.
    b = sp.Rational(4, 5)
    f = 4*b/(3+b*b)
    chi = (1+3*b*b)/(3+b*b)
    C = sp.Matrix([[1, f], [0, 1+chi]])
    D = sp.diag(1, sp.Rational(1, 10000))
    A = sp.Matrix([[0, 1], [(1-3*b*b)/(3-b*b), 4*b/(3-b*b)]])
    P = C*D*(sp.eye(2)+C*D).inv()
    Aeq = P*A
    disc = sp.factor(sp.trace(Aeq)**2-4*Aeq.det())
    assert disc < 0  # Exact rational certificate of complex equilibrium characteristics.
    Q = np.array((D*C).row_join(-D).col_join((-C).row_join(sp.eye(2))), dtype=float)
    At = np.zeros((4, 4))
    At[2:, 2:] = np.array(A, dtype=float)
    k = .01
    continuous = -Q-1j*k*At
    growth = float(np.max(np.linalg.eigvals(continuous).real))
    assert growth > 0
    # Fixed dx, dt -> 0: the discrete growth tends to the semidiscrete growth.
    dx = .1
    theta = k*dx
    sigma = np.sqrt(float(chi))
    Lspace = (1j*np.sin(theta)*At)/dx
    Lspace[2:, 2:] += sigma*(1-np.cos(theta))/dx*np.eye(2)
    semidiscrete_growth = float(np.max(np.linalg.eigvals(-Q-Lspace).real))
    table = []
    for dt in [.01, .001, .0001]:
        g = pdars(dt*Q, dt*Lspace)
        rho = float(max(abs(np.linalg.eigvals(g))))
        table.append({'dt': dt, 'rho': rho, 'log_rho_over_dt': float(np.log(rho)/dt)})
        assert rho > 1
    assert abs(table[-1]['log_rho_over_dt']/semidiscrete_growth-1) < 1e-3
    return {'beta': str(b), 'f': str(f), 'q_rate': 1, 'K_energy_rate': 1,
            'K_momentum_rate': .0001, 'equilibrium_characteristic_discriminant': str(disc),
            'equilibrium_characteristics': [str(z) for z in np.linalg.eigvals(np.array(Aeq, dtype=float))],
            'k': k, 'continuous_growth_rate': growth,
            'semidiscrete_growth_rate': semidiscrete_growth, 'discrete_convergence': table}


def coupled_finite_time_checks():
    rng = np.random.default_rng(2402)
    worst_ratio = 0.
    count = 0
    for _ in range(40):
        u, a, h, speeds = matrices(rng.uniform(0, .98), rng.normal(size=3))
        Qrate = source_matrix(u, 1., 10**rng.uniform(-3, 2), 10**rng.uniform(-5, 0))
        H = np.zeros((8, 8))
        H[:4, :4] = np.eye(4)
        H[4:, 4:] = h
        chol = np.linalg.cholesky(H).T
        ci = np.linalg.inv(chol)
        Sh = -chol@Qrate@ci
        omega = max(0., float(max(np.linalg.eigvalsh((Sh+Sh.T)/2))))
        dt = .2/max(1., omega)
        R = max(1., max(float(max(abs(np.linalg.eigvals(ad))))**2/sd**2
                        for ad, sd in zip(a, speeds)))
        nu = .99/(R*sum(speeds))
        bound = 1/(1-dt*omega)
        for theta in np.vstack([np.zeros(3), np.full(3, np.pi), rng.uniform(-np.pi, np.pi, (30, 3))]):
            wr = nu*(np.dot(speeds, 1-np.cos(theta))*np.eye(4)
                     + 1j*np.einsum('d,dij->ij', np.sin(theta), a))
            W = np.zeros((8, 8), dtype=complex)
            W[4:, 4:] = wr
            gh = chol@pdars(dt*Qrate, W)@ci
            op = float(np.linalg.svd(gh, compute_uv=False)[0])
            assert op <= bound+1e-10
            worst_ratio = max(worst_ratio, op/bound)
            count += 1
    return {'matrix_mode_checks': count, 'maximum_norm_over_proved_bound': worst_ratio}


def full_transport_cfl_finite_time_checks():
    """Fixed source rates, L=.98, and mesh/step refinement; H is fixed per state."""
    rng = np.random.default_rng(2403)
    count = 0
    powers = 0
    worst_ratio = 0.
    identity_error = 0.
    largest_L = 0.
    refinement = []
    resolvent_checks = 0
    largest_resolvent_ratio = 0.
    unrestricted_step_checks = 0
    enlarged_checks = 0
    enlarged_identity_error = 0.
    for beta in [0., .8, .95]:
        for ke, km in [(1., 1e-4), (.1, .1), (10., .001)]:
            u, a, h, speeds = matrices(beta, [1., .3, .2])
            Qrate = source_matrix(u, 1., ke, km)
            H = np.zeros((8, 8))
            H[:4, :4] = np.eye(4)
            H[4:, 4:] = h
            chol = np.linalg.cholesky(H).T
            ci = np.linalg.inv(chol)
            Qh = chol@Qrate@ci
            Sh = -Qh
            M = float(np.linalg.svd(Sh, compute_uv=False)[0])
            omega = max(0., float(max(np.linalg.eigvalsh((Sh+Sh.T)/2))))
            Csource = np.eye(4)
            Csource[0, 1:] = u[1:]/u[0]
            Csource[1:, 1:] += closure(u)
            Dsource = np.diag([ke, km, km, km])
            transform = np.block([[np.eye(4), Dsource], [-Csource, np.eye(4)]])
            B = np.eye(4)+Csource@Dsource
            transformed_source = transform@Qrate@np.linalg.inv(transform)
            expected = np.zeros((8, 8))
            expected[4:, 4:] = B
            assert np.linalg.norm(transformed_source-expected) < 1e-10
            f = float(np.linalg.norm(u[1:]/u[0]))
            chi = (1+3*beta*beta)/(3+beta*beta)
            KB = 1+f*km/(np.sqrt(1+ke)+np.sqrt(1+km*(1+chi)))**2
            K = float(np.linalg.cond(transform@ci))*KB
            Cstar = K*M*(1+3*K)/2
            for source_step in np.logspace(-6, 6, 13):
                resolvent = np.linalg.inv(np.eye(8)+source_step*Qh)
                rop = float(np.linalg.svd(resolvent, compute_uv=False)[0])
                assert rop <= K+1e-8
                largest_resolvent_ratio = max(largest_resolvent_ratio, rop/K)
                resolvent_checks += 1
            # Large implicit steps test the version with no omega*dt restriction.
            for source_step in [1., 10., 100.]:
                nu_full = .98/sum(speeds)
                theta = np.array([.4, 1.2, 2.8])
                wr = nu_full*(np.dot(speeds, 1-np.cos(theta))*np.eye(4)
                              + 1j*np.einsum('d,dij->ij', np.sin(theta), a))
                W = np.zeros((8, 8), dtype=complex)
                W[4:, 4:] = wr
                G = pdars(source_step*Qh, chol@W@ci)
                assert np.linalg.svd(G, compute_uv=False)[0] <= 1+Cstar*source_step+1e-8
                unrestricted_step_checks += 1
            # The enlarged gray envelope L=1+q/2, with fixed q_rate=1.
            ell_rate = .5
            Aconstant = 3*np.sqrt(3/2)*K*K*ell_rate
            Bconstant = 9/4*K*K*ell_rate**2
            assert (Cstar+Aconstant)**2/2 >= Bconstant
            for dt in [.01, .1, 1., 10.]:
                theta = np.array([.7, 1.8, 2.5])
                normalized_nu = 1/sum(speeds)
                wr = normalized_nu*(np.dot(speeds, 1-np.cos(theta))*np.eye(4)
                                     + 1j*np.einsum('d,dij->ij', np.sin(theta), a))
                Wbar = np.zeros((8, 8), dtype=complex)
                Wbar[4:, 4:] = wr
                Wbar = chol@Wbar@ci
                Dbar = np.eye(8)-Wbar
                Gbar = pdars(dt*Qh, Wbar)
                G = pdars(dt*Qh, (1+ell_rate*dt)*Wbar)
                Rhalf = np.linalg.inv(np.eye(8)+dt*Qh/2)
                Rone = np.linalg.inv(np.eye(8)+dt*Qh)
                exact_difference = Rhalf@(-ell_rate*dt*(Wbar@Rone@Dbar+Dbar@Rone@Wbar)
                                           +(ell_rate*dt)**2*Wbar@Rone@Wbar)/2
                residual = float(np.linalg.norm(G-Gbar-exact_difference))
                enlarged_identity_error = max(enlarged_identity_error, residual)
                assert residual < 1e-10
                assert np.linalg.svd(G-Gbar, compute_uv=False)[0] <= Aconstant*dt+Bconstant*dt**2+1e-8
                op = float(np.linalg.svd(G, compute_uv=False)[0])
                assert np.log(op) <= (Cstar+Aconstant)*dt+1e-10
                enlarged_checks += 1
            delta = .5
            C = M*(4-delta)/((2-delta)*(1-delta))
            base_dt = .02/max(1., omega)
            T = 20*base_dt
            for dim in [1, 2, 3]:
                nu = .98/sum(speeds[:dim])
                largest_L = max(largest_L, nu*sum(speeds[:dim]))
                phases = np.vstack([np.zeros(dim), np.full(dim, np.pi),
                                    rng.uniform(-np.pi, np.pi, (12, dim))])
                for dt in [base_dt, base_dt/2, base_dt/4]:
                    # Grid spacing is refined with dt at fixed radiation CFL.
                    dx = dt/nu
                    assert omega*dt <= delta
                    for theta in phases:
                        wr = nu*(np.dot(speeds[:dim], 1-np.cos(theta))*np.eye(4)
                                 + 1j*np.einsum('d,dij->ij', np.sin(theta), a[:dim]))
                        W = np.zeros((8, 8), dtype=complex)
                        W[4:, 4:] = wr
                        Wh = chol@W@ci
                        eye = np.eye(8)
                        D = eye-Wh
                        G0 = (eye+D@D)/2
                        assert np.linalg.svd(G0, compute_uv=False)[0] <= 1+1e-10
                        assert np.linalg.svd(D, compute_uv=False)[0]**2 <= 3+1e-10
                        G = pdars(dt*Qh, Wh)
                        Rhalf = np.linalg.inv(eye+dt*Qh/2)
                        Rone = np.linalg.inv(eye+dt*Qh)
                        correction = (Rhalf-eye)@G0 + Rhalf@D@(Rone-eye)@D/2
                        residual = float(np.linalg.norm(G-G0-correction))
                        identity_error = max(identity_error, residual)
                        assert residual < 1e-11
                        assert np.linalg.svd(correction, compute_uv=False)[0] <= C*dt+1e-10
                        op = float(np.linalg.svd(G, compute_uv=False)[0])
                        assert op <= 1+C*dt+1e-10
                        worst_ratio = max(worst_ratio, op/(1+C*dt))
                        nsteps = int(round(T/dt))
                        power_norm = float(np.linalg.svd(np.linalg.matrix_power(G, nsteps), compute_uv=False)[0])
                        assert np.log(power_norm) <= C*T+1e-10
                        count += 1
                        powers += 1
                    if beta == .8 and ke == 1. and dim == 1:
                        refinement.append({'dt': dt, 'dx': dx, 'source_rate_norm': M,
                                           'omega': omega, 'T': T, 'C_delta': C,
                                           'nsteps': nsteps})
    return {'matrix_mode_checks': count, 'matrix_power_checks': powers,
            'largest_L': largest_L, 'maximum_norm_over_proved_bound': worst_ratio,
            'maximum_resolvent_decomposition_residual': identity_error,
            'gray_resolvent_checks': resolvent_checks,
            'maximum_gray_resolvent_norm_over_bound': largest_resolvent_ratio,
            'unrestricted_source_step_checks': unrestricted_step_checks,
            'enlarged_cfl_checks': enlarged_checks,
            'maximum_enlarged_transport_identity_residual': enlarged_identity_error,
            'example_fixed_source_refinement_family': refinement}


if __name__ == '__main__':
    result = {'symbolic': symbolic_certificate(), 'scalar_absorption': scalar_checks(),
              'continuous_obstruction': continuous_obstruction(),
              'coupled_finite_time': coupled_finite_time_checks(),
              'full_transport_cfl_finite_time': full_transport_cfl_finite_time_checks()}
    out = Path(__file__).resolve().parents[2]/'data/m1/exchange_bounds.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
