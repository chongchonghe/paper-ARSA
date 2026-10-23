const NOMAIN = true
include("checkerboard.jl")
println("Flux checkerboard, hydro on, beta_order=1, calC=1e4:")
for lim in (ConstantChat(), FixedPhi(1.0)), Pcal in (1e-2, 1.0, 100.0), tau in (0.3, 10.0, 1e3)
    gg = growth(lim; D=1, n=16, tau, Pcal, hydro=true, nst=40, c=1e4, beta_order=1, pert=:F)
    @printf("  %-14s Pcal=%6.2f tau=%7.1f  g/step=%.5f\n", typeof(lim), Pcal, tau, gg)
end
