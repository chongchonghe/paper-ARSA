const NOMAIN = true
include("checkerboard.jl")
for rec in (PCM(), PLM(), PPM()), D in (1, 2)
    n = 16; dx = 1/n; phi0 = 1e-3; cfl = 0.4; tau = 0.3
    dt = cfl / (phi0 * D / dx)
    den = dt * (D / sqrt(3) / dx - (tau/dx) / 2)
    for fac in (0.8, 0.9, 0.95, 1.0, 1.05)
        g, lhs = growth(; D, n, phi=fac/den, phi0, tau, Pcal=1e-6, cfl, rec)
        @printf("%s D=%d fac=%.2f g=%.5f\n", typeof(rec), D, fac, g)
    end
end
