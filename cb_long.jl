const NOMAIN = true
include("checkerboard.jl")
using Honghu, Printf
function series(lim; D, n=16, phi0=1e-3, tau, Pcal, nst=300)
    cv, gamma = 1.0, 5/3; a_rad = Pcal*(gamma-1)*cv; dx = 1/n
    rad = RadiationPhysics(; c=1.0, chat=phi0, a_rad, beta_order=0, opacity=ConstantOpacity(tau/dx), chat_limiter=lim, tolerance=1e-14)
    sim = Simulation(; n=ntuple(d -> d <= D ? n : 1, 3), eos=IdealGas(; gamma, cv), radiation=rad, hydro=false,
                     boundaries=Boundaries(Periodic()), cfl_radiation=0.4, radiation_reconstruction=PCM())
    cell(x) = floor(Int, x/dx)
    set_initial_conditions!(sim) do x, y, z
        s = iseven(cell(x) + (D>1)*cell(y) + (D>2)*cell(z)) ? 1.0 : -1.0
        (rho=1.0, T=1.0, Erad=a_rad*(1 + 1e-4*s))
    end
    g = sim.grid; I = (g.lo[1]:g.hi[1], g.lo[2]:g.hi[2], g.lo[3]:g.hi[3])
    sgn = [(-1.0)^((i-g.lo[1]) + (D>1)*(j-g.lo[2]) + (D>2)*(k-g.lo[3])) for i in I[1], j in I[2], k in I[3]]
    dt = Honghu.radiation_timestep(sim.radiation, g, sim.cfl_radiation)
    E = Float64[]; Tt = Float64[]; tot = Float64[]
    for s in 1:nst
        step!(sim, dt)
        push!(E, sum(sim.R[1, I...] .* sgn)/length(sgn)/a_rad)
        push!(Tt, sum((gas_temperature(sim)) .* sgn)/length(sgn))
        m = sim.R[1, I...] ./ a_rad; push!(tot, sqrt(sum(abs2, m .- sum(m)/length(m))/length(m)))
    end
    E, Tt, tot
end
E, Tt, tot = series(StabilityBound(); D=2, tau=1.0, Pcal=1e-6)
for s in (1, 10, 30, 60, 100, 200, 300)
    @printf("step %3d  E-checker %+.3e  T-checker %+.3e  rms(E) %.3e\n", s, E[s], Tt[s], tot[s])
end
println("--- fixed phi at the limiter's value, 2D")
sim0, _ = (nothing, nothing)
# the limiter's phi for this case: beta / (c dt D)
n=16; dx=1/n; dt=0.4/(1e-3*2*n); chi=1.0/dx
phi = 0.9 / (dt*(2*n/sqrt(3) - chi/2))
println("phi = ", phi)
for (lab, lim, D) in (("fixed 2D", FixedPhi(phi), 2), ("limiter 1D", StabilityBound(), 1), ("RSLA 2D", ConstantChat(), 2))
    E, Tt, tot = series(lim; D, tau=1.0, Pcal=1e-6)
    @printf("%-11s steps 30/100/300: E %+.3e %+.3e %+.3e   T %+.3e\n", lab, E[30], E[100], E[300], Tt[300])
end
