const NOMAIN = true
include("checkerboard.jl")
println("Standard RSLA (ConstantChat) at C_r = 0.4, calC = 1e4, chat = 10 a, large calP:")
for hydro in (false, true), Pcal in (1.0, 100.0), tau in (0.3, 10.0), D in (1, 2)
    gg = growth(ConstantChat(); D, n=16, tau, Pcal, hydro, nst=40, c=1e4, beta_order=1)
    @printf("  hydro=%-5s D=%d Pcal=%6.1f tau=%5.1f  g/step=%.5f\n", hydro, D, Pcal, tau, gg)
end
