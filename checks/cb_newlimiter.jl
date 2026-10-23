const NOMAIN = true
include("checkerboard.jl")
println("New StabilityBound (beta = 1 and 0.9) on a uniform medium, E checkerboard:")
for beta in (1.0, 0.9), D in (1, 2, 3), (tau, Pcal) in ((0.6, 1e-6), (1.0, 1e-6), (0.3, 1e-2), (3.0, 1e-1))
    gg = growth(StabilityBound(; beta); D, n=16, tau, Pcal, nst=30)
    @printf("  beta=%.1f D=%d tau=%.1f Pcal=%.0e  g/substep=%.5f\n", beta, D, tau, Pcal, gg)
end
