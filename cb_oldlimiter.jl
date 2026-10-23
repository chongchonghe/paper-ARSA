const NOMAIN = true
include("checkerboard.jl")
println("Current StabilityBound (tau = rho kappa dx_min) on a uniform medium:")
for D in (1, 2, 3), tau in (0.6, 1.0)
    gg = growth(StabilityBound(); D, n=16, tau, Pcal=1e-6, nst=30)
    @printf("  D=%d tau_cell=%.1f  g/substep=%.5f\n", D, tau, gg)
end
