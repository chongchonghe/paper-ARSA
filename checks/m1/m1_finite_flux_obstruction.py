"""Finite-f M1: continuous-time obstruction to an opacity-uniform timestep bound.

Exact algebra certifies a complex pair of equilibrium propagation speeds.
Small-phase perturbation then gives a growing semidiscrete mode. Numerical
values at a fixed phase and shrinking timesteps supplement that proof.
This is the frozen radiation substep, not an evolving nonlinear background.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as s


def symbolic():
    b, e, m = s.symbols('b e m', positive=True)
    f = 4*b/(3+b*b)
    chi = (1+3*b*b)/(3+b*b)
    cp = 4*b/(3-b*b)
    a = s.factor(chi-f*cp)
    d = 1+chi
    C = s.Matrix([[1, f], [0, d]])
    M = s.diag(e, m)
    A = s.Matrix([[0, 1], [a, cp]])
    F = C*M*(s.eye(2)+C*M).inv()
    u, v = e/(1+e), d*m/(1+d*m)
    t = cp+f*a/(d*(1+e))
    disc = s.factor(s.trace(F*A)**2-4*(F*A).det())
    assert s.factor(disc-v*(v*t*t+4*u*a)) == 0
    assert s.factor(a-(1-3*b*b)/(3-b*b)) == 0
    derivative = s.factor(s.diff(disc, m).subs(m, 0))
    assert s.factor(derivative-4*u*d*a) == 0
    assert s.simplify(f.subs(b, 1/s.sqrt(3))-2*s.sqrt(3)/5) == 0

    values = {b:s.Rational(4,5), e:s.Integer(1), m:s.Rational(1,10000)}
    C, M, A, F = [x.subs(values) for x in [C, M, A, F]]
    Q = (M*C).row_join(-M).col_join((-C).row_join(s.eye(2)))
    J = C.inv().col_join(s.eye(2))
    T = s.eye(2).row_join(M).col_join((-C).row_join(s.eye(2)))
    Ti = T.inv()
    assert T.det() > 0
    Qtransformed = (T*Q*Ti).applyfunc(s.factor)
    assert Qtransformed[:2, :] == s.zeros(2,4)
    assert Qtransformed[2:, :2] == s.zeros(2)
    assert Qtransformed[2:, 2:] == s.eye(2)+C*M
    # Source kernel is semisimple; relaxing block has positive eigenvalues.
    assert Qtransformed[2, 2] > 0 and Qtransformed[3, 3] > 0
    P = Ti*s.diag(s.eye(2), s.zeros(2))*T
    E = P[2:, :]
    assert E*J == s.eye(2) and J*E == P
    Afull = s.diag(s.zeros(2), A)
    assert E*Afull*J == F*A
    assert s.factor(disc.subs(values)) == -s.Rational(25310094734,180228450156561)
    return dict(general_discriminant='v*(v*t^2+4*u*a)',
                a='(1-3*b^2)/(3-b^2)',
                u='mu_E/(1+mu_E)', v='(1+chi)*mu_M/(1+(1+chi)*mu_M)',
                t="chi_prime+f*a/((1+chi)*(1+mu_E))",
                derivative_at_muM_zero='4*u*(1+chi)*a',
                flux_threshold='2*sqrt(3)/5',
                particular_discriminant=str(s.factor(disc.subs(values))),
                source_kernel_semisimple=True, projected_generator_identity=True)


def numerical(dps):
    with mp.workdps(dps):
        b=mp.mpf(4)/5; f=4*b/(3+b*b); chi=(1+3*b*b)/(3+b*b)
        cp=4*b/(3-b*b)
        C=mp.matrix([[1,f],[0,1+chi]])
        M=mp.diag([1,mp.mpf(1)/10000])
        A=mp.matrix([[0,1],[chi-f*cp,cp]])
        Q=mp.matrix(4)
        for i in range(2):
            for j in range(2):
                Q[i,j]=(M*C)[i,j];Q[i,j+2]=-M[i,j]
                Q[i+2,j]=-C[i,j];Q[i+2,j+2]=int(i==j)
        tau=mp.mpf(1000); angle=mp.mpf(1)/2000
        sn=2*angle/(1+angle*angle);cs=(1-angle*angle)/(1+angle*angle)
        W=mp.matrix(4)
        wr=mp.sqrt(chi)*(1-cs)*mp.eye(2)+1j*sn*A
        for i in range(2):
            for j in range(2):W[i+2,j+2]=wr[i,j]
        B=-tau*Q-W
        spectrum=mp.eig(B,left=False,right=False)
        growth=max(mp.re(z) for z in spectrum)
        assert growth > mp.mpf('0.000005')
        rows=[]
        for exponent in [0,1,2,3,4,6,8]:
            step=mp.mpf(10)**(-exponent)
            D=mp.eye(4)-step*W
            G=(2*mp.eye(4)+step*tau*Q)**-1*(mp.eye(4)+D*(mp.eye(4)+step*tau*Q)**-1*D)
            vals=mp.eig(G,left=False,right=False)
            rho=max(abs(z) for z in vals)
            rate=mp.log(rho)/step
            assert rate > 0
            rows.append(dict(nu=mp.nstr(step,15),q=mp.nstr(step*tau,15),
                             spectral_radius=mp.nstr(rho,35),
                             log_rho_over_nu=mp.nstr(rate,35)))
        # Simultaneously shrink cell width h and nu=h, so dt=h^2, chat=1.
        # Binary matrix powers execute the discrete update to physical time T=1.
        refinement=[]
        for exponent in [4,5,6,7]:
            h=mp.mpf(10)**(-exponent);step=h;n=10**(2*exponent)
            D=mp.eye(4)-step*W
            G=(2*mp.eye(4)+step*tau*Q)**-1*(mp.eye(4)+D*(mp.eye(4)+step*tau*Q)**-1*D)
            vals, vectors=mp.eig(G,left=False,right=True)
            k=max(range(4),key=lambda k:abs(vals[k]));z=vectors[:,k]
            z=z/mp.norm(z)
            actual=mp.norm((G**n)*z)
            predicted=abs(vals[k])**n
            err=abs(actual/predicted-1)
            assert err < mp.mpf('1e-35')
            refinement.append(dict(dx=mp.nstr(h,15),alpha=mp.nstr(tau/h,15),
                dt=mp.nstr(h*h,15),steps=n,physical_time=1,
                log_amplitude_gain=mp.nstr(mp.log(actual),30),
                matrix_power_relative_error=mp.nstr(err,5)))
        return dict(digits=dps, generator_growth=mp.nstr(growth,35),
                    timestep_sequence=rows,refinement=refinement)


if __name__=='__main__':
    result=dict(symbolic=symbolic(),numerical=[numerical(70),numerical(100)],
        scope='Frozen static full-exchange M1, equal opacities, eta=1, f=80/91; no uniformity in opacity without extra restrictions.',
        limitation='Constant background powers; not a nonlinear streaming transient or full hydrodynamics.')
    assert result['numerical'][0]['generator_growth']==result['numerical'][1]['generator_growth']
    path=Path(__file__).resolve().parents[2]/'data/m1/finite_flux_obstruction.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
