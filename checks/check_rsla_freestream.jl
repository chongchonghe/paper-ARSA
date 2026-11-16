# Checks Proposition 1 (paper Sec. 3.1.6) in the nonlinear code: the standard RSLA is stable at nu0 = 1/sqrt(D)
# for every radiation state, here free-streaming M1 radiation in directions not aligned with the grid.
#
# A uniform, optically thin (kappa ~ 0), periodic box of N^D cubic cells holds radiation E = 1 with reduced flux
# f in direction n, plus white noise of relative amplitude 1e-8 in E and F (it contains every Fourier mode,
# the checkerboard included). Gas and hydrodynamics are frozen; only the M1 transport (HLL, speeds
# max(0.1, sqrt(T_dd))) acts, with piecewise-constant reconstruction (PCM, as the analysis assumes) and with
# PPM (as in the tests). Each run takes NSTEP radiation steps at
# C_r = s sqrt(D), i.e. nu0 = C_r/D = s/sqrt(D): s = 1 is Proposition 1's step; s = 1.05 is a control above
# it, which should grow where sum_d sigma_d = sqrt(D) (radiation along a full diagonal) and stay stable where
# the sum is smaller. Prints, per run, sum_d sigma_d of the background, the predicted checkerboard factor
# nu0 sum_d sigma_d, and the growth of the noise, max|delta E| at the end over that at the start (and the
# largest ratio over the run). Writes no files.
#
#   julia --project=. checks/check_rsla_freestream.jl
using Distributed, Printf
addprocs(8; exeflags=`--project=$(Base.active_project()) -t 1`)

@everywhere begin
using Honghu, Random, Statistics
const NSTEP, AMP = 3000, 1e-8

function run_case(D, N, dir, f, s, rec)
    n = dir ./ sqrt(sum(abs2, dir))
    rad = RadiationPhysics(; c=1.0, chat=1.0, a_rad=1.0, opacity=ConstantOpacity(1e-12))
    dims = ntuple(d -> d <= D ? N : 1, 3)
    sim = Simulation(; n=dims, xlo=(0.0, 0.0, 0.0), xhi=Float64.(dims) ./ N,
                     eos=IdealGas(gamma=5 / 3, cv=1e10), radiation=rad, hydro=false,
                     boundaries=Boundaries(Periodic()), cfl_radiation=s * sqrt(D),
                     radiation_reconstruction=rec == :pcm ? PCM() : PPM())
    set_initial_conditions!(sim) do x, y, z
        (rho=1.0, p=1.0, Erad=1.0, Fx=f * n[1], Fy=f * n[2], Fz=f * n[3])
    end
    g = sim.grid; R = sim.R; b = Honghu.gbase(1)
    I = (g.lo[1]:g.hi[1], g.lo[2]:g.hi[2], g.lo[3]:g.hi[3])
    rng = MersenneTwister(1)
    for v in (Honghu.IERAD, Honghu.IFRX, Honghu.IFRY, Honghu.IFRZ)
        view(R, b + v, I...) .+= AMP .* randn(rng, length.(I)...)
    end
    amp() = (E = view(R, b + Honghu.IERAD, I...); maximum(abs.(E .- mean(E))))
    a0 = amp(); amax = a0
    cb = s -> (amax = max(amax, amp()); true)
    run!(sim; tmax=1e30, maxsteps=NSTEP, verbose=false, callback=cb)
    t = Honghu.eddington_tensor(LevermoreClosure(), f * n[1], f * n[2], f * n[3])
    sig = sum(max(0.1, sqrt(t[d, d])) for d in 1:D)
    (; D, dir, f, s, rec, sig, factor=s / sqrt(D) * sig, growth=amp() / a0, gmax=amax / a0,
       steps=sim.diag.steps, finite=all(isfinite, R))
end
end

dirs2 = ((1, 0, 0), (cosd(22.5), sind(22.5), 0), (1, 1, 0))
dirs3 = ((1, 0, 0), (1, 1, 0), (1, 2, 3), (1, 1, 1))
jobs = [(D, N, dir, f, s, rec) for (D, N, dirs) in ((2, 32, dirs2), (3, 16, dirs3)) for dir in dirs
        for f in (0.9, 0.999) for s in (1.0, 1.05) for rec in (:pcm, :ppm)]
res = pmap(j -> run_case(j...), jobs)

@printf("%d radiation steps per run; noise amplitude %g; nu0 = s/sqrt(D)\n", NSTEP, AMP)
@printf("%2s %-22s %6s %5s %4s %10s %16s %12s %12s %6s\n", "D", "direction", "f", "s", "rec", "sum sigma", "nu0 sum sigma", "growth end", "growth max", "steps")
for r in res
    dir = join((@sprintf("%.3g", x) for x in r.dir), ",")
    @printf("%2d %-22s %6.3f %5.2f %4s %10.4f %16.4f %12.3g %12.3g %6d%s\n", r.D, "(" * dir * ")", r.f, r.s, string(r.rec), r.sig,
            r.factor, r.growth, r.gmax, r.steps, r.finite ? "" : "  NON-FINITE")
end
bad = [r for r in res if r.rec == :pcm && r.s == 1.0 && !(r.growth < 1 && r.finite)]
println(isempty(bad) ? "\nPCM, nu0 = 1/sqrt(D): every run decays by the end (after a transient), as Proposition 1 states." :
        "\nPCM, nu0 = 1/sqrt(D): some runs grow: see the rows with rec = pcm, s = 1.00 and growth end > 1.")
