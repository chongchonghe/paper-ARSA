"""Exact certificates and independent matrix checks for the thick-cell extension.

These are checks of the frozen first-order radiation substep, not Quokka runs.
Run: python3 checks/m1/m1_thick_cells.py
"""

import json
from pathlib import Path

import numpy as np
import sympy as s

from m1_stability import matrices
from m1_exchange_bounds import pdars


def hurwitz_minors(coefficients):
    n = len(coefficients) - 1
    h = s.Matrix(n, n, lambda i, j: coefficients[2*j-i+1]
                 if 0 <= 2*j-i+1 <= n else 0)
    return [s.factor(h[:j, :j].det()) for j in range(1, n+1)]


def exact_counterexample():
    I = s.eye(3)
    mu = s.diag(s.Rational(1, 5000), s.Rational(1, 10000), s.Rational(1, 10000))
    C = s.Matrix([[1, s.Rational(2, 25), 0],
                  [0, s.Rational(48, 25), 0], [0, 0, s.Rational(26, 25)]])
    iA = s.Matrix([[0, 0, 1], [0, 0, 11], [-s.Rational(6, 7), s.Rational(1, 14), 0]])
    Q = (mu*C).row_join(-mu).col_join((-C).row_join(I))
    W = s.diag(s.zeros(3), I/5+iA)
    nu, tau = 20000, 1000
    Id = s.eye(6)
    D = Id-nu*W
    G = (2*Id+nu*tau*Q).inv()*(Id+D*(Id+nu*tau*Q).inv()*D)
    trace4 = s.trace(G**4)
    assert trace4 > 12  # |tr G^4| <= 6 rho(G)^4: exact instability certificate.

    # An independent, unscaled full 8-variable M1 matrix gives the same growth.
    beta = np.sqrt(11/13)
    u, a, _, speeds = matrices(beta, [0, 1, 0])
    from m1_exchange_bounds import source_matrix
    Q8 = source_matrix(u, 1., 1/5000, 1/10000)
    W8 = np.zeros((8, 8), dtype=complex)
    W8[4:, 4:] = nu*(speeds[0]*np.eye(4)+1j*a[0])
    rho = float(max(abs(np.linalg.eigvals(np.array(G, float)))))
    rho8 = float(max(abs(np.linalg.eigvals(pdars(nu*tau*Q8, W8)))))
    assert abs(rho-rho8) < 1e-6

    # The checkerboard criterion holds strictly, in exact arithmetic.
    L = s.Rational(nu, 5)
    xi = s.Rational(1, 5000)
    q = nu*tau
    checker_limit = 1+q/(2+xi*q)
    assert L < checker_limit
    Dc = Id-s.diag(s.zeros(3), 2*L*I)
    Gc = (2*Id+q*Q).inv()*(Id+Dc*(Id+q*Q).inv()*Dc)
    # Each checkerboard source channel is similar to a symmetric scalar channel;
    # the exact inequality above proves contraction in its channel norm.

    # Strict same-mode semidiscrete damping: all six Hurwitz determinants positive.
    semicoeff = (-tau*Q-W).charpoly().all_coeffs()
    semiminors = hurwitz_minors(semicoeff)
    assert all(v > 0 for v in semiminors)

    # Continuous stability for EVERY real nonzero k, not just the sampled mode.
    k = s.symbols('k', real=True)
    coefficients = (-tau*Q-s.diag(s.zeros(3), k*iA)).charpoly().all_coeffs()
    assert coefficients[-1] == 0
    assert coefficients[-2] != 0  # simple neutral eigenvalue when k != 0
    continuous_minors = hurwitz_minors(coefficients[:-1])
    for v in continuous_minors:
        poly = s.Poly(s.expand(v), k)
        assert all(c > 0 and power[0] % 2 == 0 for power, c in poly.terms())
    # k=0 is source-only: three zero eigenvalues are semisimple and others negative.
    F = C*mu*(I+C*mu).inv()
    transform = I.row_join(mu).col_join((-C).row_join(I))
    P = transform.inv()*s.diag(I, s.zeros(3))*transform
    assert P*P == P and Q*P == s.zeros(6) and P*Q == s.zeros(6)
    Ginf = (I+(I-nu*F*(I/5+iA))**2)/2
    assert s.trace(Ginf**4) > 11  # > dimension=3; persists at infinite optical depth.

    # Derive transformed flux/closure independently from the exact Levermore formulas.
    b = s.sqrt(143)/13
    f = 4*b/(3+b*b)
    chi = (1+3*b*b)/(3+b*b)
    cp = 4*b/(3-b*b)
    aperp = s.Matrix([[0, 0, 1], [0, 0, b], [(1-chi+f*cp)/2, -cp/2, 0]])
    S = s.diag(1, s.sqrt(143), s.I)
    Coriginal = s.Matrix([[1, f, 0], [0, 1+chi, 0], [0, 0, 1+(1-chi)/2]])
    assert (s.I*S*aperp*S.inv()-iA).applyfunc(s.simplify) == s.zeros(3)
    assert (S*Coriginal*S.inv()-C).applyfunc(s.simplify) == s.zeros(3)
    return dict(beta=str(b), f=str(s.simplify(f)), chi=str(s.simplify(chi)),
                nu=nu, tau_cell=tau, L=str(L), q=q, checker_limit=str(checker_limit),
                spectral_radius=rho, full_eight_variable_radius=rho8,
                checker_radius=float(max(abs(np.linalg.eigvals(np.array(Gc, float))))),
                trace_fourth_power=str(trace4), trace_fourth_power_decimal=float(trace4),
                semidiscrete_characteristic_coefficients=[str(v) for v in semicoeff],
                semidiscrete_hurwitz_minors=[str(v) for v in semiminors],
                continuous_characteristic_coefficients=[str(v) for v in coefficients],
                continuous_hurwitz_minors=[str(v) for v in continuous_minors],
                stiff_limit_trace_fourth_power=str(s.trace(Ginf**4)),
                stiff_limit_radius=float(max(abs(np.linalg.eigvals(np.array(Ginf, float))))),
                scope='Positive frozen responses; not a nonrelativistic ideal gas at energy balance.')


def compatible_positive_example():
    b, f, chi, cp, aa = map(s.Rational, ['1/3', '3/7', '3/7', '6/13', '3/13'])
    ue, um = s.Rational(1, 100), s.Rational(1, 10**8)
    c, ct = 1+chi, 1+(1-chi)/2
    h = f*um/(ue-c*um)
    Hr = s.diag(s.Matrix([[1, h], [h, (1+h*cp)/aa]]), 1, 1)
    C = s.diag(1, c, ct, ct)
    C[0, 1] = f
    M = s.diag(ue, um, um, um)
    Hg = C.T*Hr*M.inv()
    A = s.diag(s.Matrix([[0, 1], [aa, cp]]), b, b)
    Q = (M*C).row_join(-M).col_join((-C).row_join(s.eye(4)))
    H = s.diag(Hg, Hr)
    assert Hr*A == A.T*Hr and Hg == Hg.T and H*Q == Q.T*H
    assert all(H[:j, :j].det() > 0 for j in range(1, 9))
    assert H*Q == (C.T*Hr*C).row_join(-C.T*Hr).col_join((-Hr*C).row_join(Hr))
    # Independently verify analytic flux and actual full source code reproduction.
    u, an, _, speeds = matrices(float(b), [1, 0, 0])
    from m1_exchange_bounds import source_matrix
    assert np.max(abs(an[0]-np.array(A, float))) < 1e-14
    assert np.max(abs(source_matrix(u, 1., float(ue), float(um))-np.array(Q, float))) < 1e-14
    R = (133+56*s.sqrt(3))/169
    xi = ue
    Hn, Qn, An = np.array(H, float), np.array(Q, float), np.array(A, float)
    chol = np.linalg.cholesky(Hn).T
    ci = np.linalg.inv(chol)
    Qh = chol@Qn@ci
    rng = np.random.default_rng(20261024)
    count, worst, power_ratio = 0, 0., 0.
    for q in np.logspace(-6, 10, 17):
        L = .999*(1+q/(2+float(xi)*q))/float(R)
        B = np.eye(8)+q*Qh
        energy = np.eye(8)+B
        chB, chE = np.linalg.cholesky(B).T, np.linalg.cholesky(energy).T
        for theta in np.r_[0., np.pi, rng.uniform(0, np.pi, 24)]:
            W = np.zeros((8, 8), dtype=complex)
            W[4:, 4:] = L*((1-np.cos(theta))*np.eye(4)+1j*np.sin(theta)*An/speeds[0])
            Wh = chol@W@ci
            D = np.eye(8)-Wh
            stage = np.linalg.solve(B, D)
            G = pdars(q*Qh, Wh)
            stage_op = np.linalg.svd(chB@stage@np.linalg.inv(chB), compute_uv=False)[0]
            full_op = np.linalg.svd(chE@G@np.linalg.inv(chE), compute_uv=False)[0]
            assert stage_op <= 1+1e-5 and full_op <= 1+1e-5
            worst = max(worst, stage_op, full_op)
            # The opacity-independent first-step and power estimate.
            bound = (1+np.linalg.norm(D, 2)**2)/2
            for n in [1, 2, 10]:
                op = np.linalg.norm(np.linalg.matrix_power(G, n), 2)
                assert op <= bound*(1+1e-5)
                power_ratio = max(power_ratio, op/bound)
            count += 1
    tau = 1000.
    sigma = speeds[0]
    # Positive quadratic root of R L <= 1+q/(2+xi q), q=nu*tau.
    ac = float(R*xi)*sigma*tau
    bc = 2*float(R)*sigma-(1+float(xi))*tau
    nu = (-bc+np.sqrt(bc*bc+8*ac))/(2*ac)
    return dict(Hr=[[str(v) for v in row] for row in Hr.tolist()],
                mu_energy=str(ue), mu_momentum=str(um), xi_max=str(xi),
                R=str(R), R_decimal=float(R), tau_cell=tau, nu_bound=nu,
                L_bound=nu*sigma, asymptotic_L_bound=float((1+1/xi)/R),
                isothermal_c_over_sound_speed=float(s.sqrt(s.Rational(3, 8)*ue/um)),
                matrix_mode_checks=count, largest_contraction_norm=float(worst),
                largest_power_norm_over_proved_bound=float(power_ratio))


def diffusion_perturbation_checks():
    """Check the exact perturbation identity, including nonequilibrium sources.

    The theorem's compactness/Lipschitz proof supplies uniform constants;
    these tests independently exercise its matrix identities and norm bounds.
    """
    from m1_stability import closure
    rng = np.random.default_rng(87122)
    Id = np.eye(8)
    mu = np.diag([.01, 1e-8, 1e-8, 1e-8])
    count, worst_residual, worst_ratio = 0, 0., 0.
    largest_L0, largest_rho = 0., 0.
    for eta in [.2, 1., 3.]:
        C0 = np.diag([1., eta+1/3, eta+1/3, eta+1/3])
        Hr0 = np.diag([1., 3., 3., 3.])
        H0 = np.block([[C0@Hr0@np.linalg.inv(mu), np.zeros((4, 4))],
                       [np.zeros((4, 4)), Hr0]])
        chol = np.linalg.cholesky(H0).T
        ci = np.linalg.inv(chol)
        Q0 = chol@np.block([[mu@C0, -mu], [-C0, np.eye(4)]])@ci
        _, a0, _, _ = matrices(0., [1, 0, 0])
        xi = max(.01, (eta+1/3)*1e-8)
        for alpha in [10., 100., 1000.]:
            f = 1/alpha  # alpha*f = B_f = 1 at all stiffnesses.
            beta = 3*f/(2+np.sqrt(4-3*f*f))
            for direction in [np.array([1., 0, 0]), *rng.normal(size=(3, 3))]:
                u, a, _, speeds = matrices(beta, direction)
                C = np.eye(4)
                C[0, 1:] = u[1:]
                C[1:, 1:] = eta*np.eye(3)+closure(u)
                Qf = chol@np.block([[mu@C, -mu], [-C, np.eye(4)]])@ci
                for dim in [1, 2, 3]:
                    for tau in [10., 1000.]:
                        dx = tau/alpha
                        sigma_sum = dim/np.sqrt(3)
                        ac, bc = sigma_sum*xi*tau, 2*sigma_sum-(1+xi)*tau
                        nu = .9*(-bc+np.sqrt(bc*bc+8*ac))/(2*ac)
                        dt = nu*dx
                        q, L0 = alpha*dt, nu*sigma_sum
                        assert L0 <= 1+q/(2+xi*q)
                        assert alpha*f == 1. and alpha*dx >= 10.
                        largest_L0 = max(largest_L0, L0)
                        Rf = np.linalg.inv(Id+q*Qf)
                        R0 = np.linalg.inv(Id+q*Q0)
                        Af = np.linalg.inv(2*Id+q*Qf)
                        A0 = np.linalg.inv(2*Id+q*Q0)
                        K = max(1., np.linalg.norm(Rf, 2), 2*np.linalg.norm(Af, 2))
                        for theta in [np.zeros(dim), np.full(dim, np.pi),
                                      *rng.uniform(-np.pi, np.pi, (4, dim))]:
                            def symbol(aa, ss):
                                W = np.zeros((8, 8), dtype=complex)
                                W[4:, 4:] = nu*(np.dot(ss[:dim], 1-np.cos(theta))*np.eye(4)
                                                + 1j*np.einsum('d,dij->ij', np.sin(theta), aa[:dim]))
                                return chol@W@ci
                            Wf = symbol(a, speeds)
                            W0 = symbol(a0, np.ones(3)/np.sqrt(3))
                            Df, D0 = Id-Wf, Id-W0
                            Gf, G0 = pdars(q*Qf, Wf), pdars(q*Q0, W0)
                            deltaD = Df-D0
                            expansion = ((Af-A0)@(Id+D0@R0@D0)
                                         + Af@(deltaD@Rf@Df+D0@(Rf-R0)@Df+D0@R0@deltaD))
                            residual = np.linalg.norm((Gf-G0)-expansion, 2)
                            assert residual < 1e-8
                            worst_residual = max(worst_residual, residual)
                            d0, df = np.linalg.norm(D0, 2), np.linalg.norm(Df, 2)
                            bound = (q*np.linalg.norm(Qf-Q0, 2)
                                     *(K/4*(1+d0*d0)+K*K/2*d0*df)
                                     + np.linalg.norm(deltaD, 2)*K/2*(K*df+d0))
                            actual = np.linalg.norm(Gf-G0, 2)
                            assert actual <= bound+1e-9
                            worst_ratio = max(worst_ratio, actual/bound)
                            largest_rho = max(largest_rho, max(abs(np.linalg.eigvals(Gf))))
                            count += 1
    return dict(matrix_mode_checks=count, maximum_identity_residual=worst_residual,
                maximum_perturbation_norm_over_bound=worst_ratio,
                largest_reference_L=largest_L0, largest_observed_spectral_radius=largest_rho,
                scope='Matrix checks of the analytic perturbation proof; not a proof by sampling.')


if __name__ == '__main__':
    results = {'counterexample': exact_counterexample(),
               'positive_example': compatible_positive_example(),
               'diffusion_perturbation': diffusion_perturbation_checks()}
    dest = Path(__file__).resolve().parents[2]/'data/m1/thick_cells.json'
    dest.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps({k: {n: v for n, v in d.items() if not isinstance(v, list)
                         and len(str(v)) < 300} for k, d in results.items()}, indent=2))
