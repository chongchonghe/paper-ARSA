import sympy as sp
a, q, lam = sp.symbols('a q lambda', positive=True)   # a = 2 nu sigma, q = nu tau, lam = 4 gamma1 P/phi
I = sp.eye(2)
Q = q*sp.Matrix([[-lam, lam/4],[4, -1]])
Z = sp.diag(0, -a)
E = I + Z
G1 = (I - Q).inv()*E
G = (I - Q/2).inv()*(sp.Rational(1,2)*I + sp.Rational(1,2)*E*G1)
G = sp.simplify(G)
for s in (1, -1):
    d = sp.factor(sp.simplify((G - s*I).det()))
    print('det(G -', s, 'I) =', d)
