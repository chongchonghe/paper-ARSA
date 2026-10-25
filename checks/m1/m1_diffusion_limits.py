"""Checks for the focused M1 diffusion note's finite-time/power-bound distinction.

The symbolic checks certify the projected long-wave derivative and the negative
characteristic discriminant. The printed finite-phase eigenvalue is a separate
high-precision numerical illustration, not an exact eigenvalue certificate.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as s


def check():
    b = s.Rational(4, 5)
    f = 4*b/(3+b*b)
    chi = (1+3*b*b)/(3+b*b)
    cp = 4*b/(3-b*b)
    C = s.Matrix([[1, f], [0, 1+chi]])
    M = s.diag(1, s.Rational(1, 10000))
    A = s.Matrix([[0, 1], [chi-f*cp, cp]])
    I2, I4 = s.eye(2), s.eye(4)
    Qs = (M*C).row_join(-M).col_join((-C).row_join(I2))
    F = C*M*(I2+C*M).inv()
    Aeq = F*A
    discriminant = s.factor(s.trace(Aeq)**2-4*Aeq.det())
    assert discriminant == -s.Rational(25310094734, 180228450156561)

    # Source-conserved variables m+Mr, and departure r-Cm from equilibrium.
    transform = I2.row_join(M).col_join((-C).row_join(I2))
    P = transform.inv()*s.diag(I2, s.zeros(2))*transform
    embedding = C.inv().col_join(I2)  # equilibrium state parametrized by r
    extract = P[2:, :]
    assert P*P == P and Qs*P == s.zeros(4) and P*Qs == s.zeros(4)
    assert extract*embedding == I2 and embedding*extract == P

    # At every q>0: the eigenvalue 1 is semisimple and its phase derivative
    # is -i Aeq. No infinite-opacity limit is needed for this obstruction.
    q = s.symbols('q', positive=True)
    R = (I4+q*Qs).inv()
    H = (2*I4+q*Qs).inv()
    Gzero = H*(I4+R)
    dW = s.diag(s.zeros(2), s.I*A)  # nu=1; HLL diffusion has zero derivative
    dG = -H*(dW*R+R*dW)
    zero = lambda mat: mat.applyfunc(lambda x: s.simplify(s.factor(x))) == s.zeros(*mat.shape)
    assert zero(Gzero*P-P) and zero(P*Gzero-P)
    assert zero(extract*dG*embedding+s.I*Aeq)
    relaxing = (transform*Gzero*transform.inv())[2:, 2:]
    # det(relaxing-I) has a positive numerator and denominator for q>0.
    det_relax = s.factor((relaxing-I2).det())
    num, den = s.fraction(det_relax)
    assert all(v > 0 for v in s.Poly(num, q).all_coeffs() if v != 0)
    assert all(v > 0 for v in s.Poly(den, q).all_coeffs() if v != 0)
    assert det_relax != 0

    # A reproducible finite-q, finite-phase illustration of the analytic result.
    q0, nu, tau = 1000, 1, 1000
    t = s.Rational(1, 2000)  # theta=2 atan(t), making sin/cos rational
    sn, cs = 2*t/(1+t*t), (1-t*t)/(1+t*t)
    sigma = s.sqrt(chi)
    W = s.diag(s.zeros(2), nu*(sigma*(1-cs)*I2+s.I*sn*A))
    D = I4-W
    G = (2*I4+q0*Qs).inv()*(I4+D*(I4+q0*Qs).inv()*D)
    radii = []
    for digits in (60, 90):
        with mp.workdps(digits):
            vals = [[mp.mpc(str(s.re(x).evalf(digits+5)),
                             str(s.im(x).evalf(digits+5))) for x in G.row(i)]
                    for i in range(4)]
            eigenvalues = mp.eig(mp.matrix(vals), left=False, right=False)
            rho = max(abs(z) for z in eigenvalues)
            assert rho > 1+mp.mpf('0.0000058')
            radii.append(mp.nstr(rho, 50))
    assert radii[0] == radii[1]
    return dict(
        exact_discriminant=str(discriminant),
        source_equilibrium_projection=True,
        semisimple_unit_eigenvalue_for_positive_q=True,
        projected_phase_derivative='-i Aeq for all q>0 at nu=1',
        relaxing_block_det_minus_identity=str(det_relax),
        q=q0, nu=nu, tau_cell=tau, theta='2 atan(1/2000)',
        spectral_radius_60_and_90_digits=radii,
        L0=str(1/s.sqrt(3)), isotropic_CFL_limit=str(1+s.Rational(q0, 2+q0)),
        uniformity_obstruction='Keep q, nu, f, theta fixed; alpha=tau/dx and dt=dx/chat.',
        scope='Frozen radiation substep, not full hydrodynamics or a nonlinear Quokka run.')


if __name__ == '__main__':
    result = check()
    output = Path(__file__).resolve().parents[2]/'data/m1/diffusion_limits.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
