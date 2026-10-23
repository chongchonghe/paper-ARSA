# The advecting pulse in the dynamic diffusion limit (He, Wibking & Krumholz 2024a; Honghu.jl's
# radhydro_pulse_dyn) with a colder background T_lo, so that the gas holds the heat capacity: a
# Gaussian pulse, T from T_lo to 2 T_lo over a width of 24 cm, in pressure balance, carried by
# v0 = 3e7 cm/s through a periodic box of 1024 cm (128 cells), kappa = 500 cm^2/g, Levermore
# closure, v/c terms on, to t = 4.8e-5 s (the pulse moves 1440 cm, 1.4 box lengths). RSLA and ARSA at
# phi0 = v0/c (chat = v0), 2 v0/c and 0.01, against chat = c at the same resolution. Then the same
# pulse at rest (static diffusion), T_lo = 1e6 K, phi0 = 0.01, to t = 4.8e-4 s.
#
# Every run is a 3D slab: n x NT x NT cubic cells, periodic in y and z, the same in every column (a set of
# parallel slabs, as a 3D code sees a planar problem). The RSLA and the ARSA take the 3D RSLA step for
# isotropic radiation at beta = 0.9, C_r = 3 nu0 = 0.9 sqrt3; the references chat = c take nu0 = 0.4 per
# direction, as the 1D references did. phi_* and the threshold of eq:tau_one are for isotropic radiation
# in 3D at beta = 1: nu0 = 1/sqrt3, sum_d sigma_d = sqrt3, so phi_* = h/(1 + h).
#
# Prints, for each T_lo: the derived variables (h, tau_c, phi_*, the eq:tau_one threshold) from the
# initial state, h/phi0 at the peak, and for each run the L1 error of T relative to the pulse, the lag of the
# pulse centroid behind the reference, phi and the sub-steps; for the static run, the diffusion length
# and the error relative to the change of the solution. Writes data/diffusion_pulse_profiles.csv;
# run_all.sh saves the printout in data/diffusion_pulse.txt.
#
#   ./run_all.sh diffusion_pulse      (or: julia --project=. -t 1 src/sim_diffusion_pulse.jl)
using Honghu, Printf
const CGS = Honghu.CGS
const MU, GAMMA, KAPPA, RHO0, W, V0, TMAX, L = 2.33 * CGS.m_u, 5 / 3, 500.0, 1.2, 24.0, 3e7, 4.8e-5, 1024.0
const CV = CGS.k_B / (MU * (GAMMA - 1))
const NT = 2                                 # cells in y and in z
const CR = 0.9 * sqrt(3)                     # C_r = 3 nu0, nu0 = 0.9/sqrt3: the 3D RSLA step at beta = 0.9
const CR_REF = 1.2                           # references: nu0 = 0.4 per direction
const SIGSUM, NU1 = sqrt(3), 1 / sqrt(3)     # isotropic radiation in 3D; largest stable nu0 (beta = 1)

Tprof(x, T_lo) = T_lo * (1 + exp(-x^2 / (2 * W^2)))
rhoprof(x, T_lo) = (T = Tprof(x, T_lo);
    RHO0 * T_lo / T + (CGS.a_rad * MU / (3 * CGS.k_B)) * (T_lo^4 / T - T^3))   # gas + radiation pressure balance

function pulse(; T_lo, phi0, lim=ConstantChat(), n=128, cfl=CR, v0=V0)
    rad = RadiationPhysics(; c=CGS.c_light, chat=phi0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           beta_order=1, chat_limiter=lim, opacity=ConstantOpacity(KAPPA))
    sim = Simulation(; n=(n, NT, NT), xlo=(-L / 2, 0.0, 0.0), xhi=(L / 2, NT * L / n, NT * L / n),
                     eos=IdealGas(; gamma=GAMMA, cv=CV), radiation=rad,
                     boundaries=Boundaries(Periodic()), cfl=0.8, cfl_radiation=cfl)
    set_initial_conditions!(sim) do x, y, z
        T = Tprof(x, T_lo)
        (rho=rhoprof(x, T_lo), vx=v0, T=T, Erad=CGS.a_rad * T^4, Fx=0.0)
    end
    sim
end

# T and phi of the first column at t = TMAX (all columns are the same); phi from the last full step (the
# last step may be cut short to hit TMAX)
function score(sim; tmax=TMAX)
    g = sim.grid
    cur = zeros(g.n[1]); prev = copy(cur); dts = [0.0, 0.0]
    cb = function (s)
        prev .= cur; cur .= view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])
        dts .= (dts[2], s.dt)
        true
    end
    run!(sim; tmax, verbose=false, callback=cb)
    (T=gas_temperature(sim)[:, 1, 1], phi=(dts[2] >= dts[1] ? cur : prev) ./ sim.radiation.c,
     steps=sim.diag.steps, sub=sim.diag.rad_substeps)
end

n = 128
dx = L / n
xs = [-L / 2 + (i - 0.5) * dx for i in 1:n]
lp(T, ref, T_lo) = sum(abs, T .- ref) / sum(abs, ref .- T_lo)          # relative to the pulse
xc = mod.(xs .- V0 * TMAX .+ L / 2, L) .- L / 2                       # x - v0 t in the periodic box
centroid(T, T_lo) = sum(xc .* (T .- T_lo)) / sum(T .- T_lo)           # in the comoving frame
profiles = []
refs = Dict{Symbol,Any}()                                             # for the clipped runs
for T_lo in (1e5, 1e6)
    T = Tprof.(xs, T_lo); rho = rhoprof.(xs, T_lo)
    h = 4 * CGS.a_rad .* T .^ 3 ./ (rho * CV)
    tau = KAPPA .* rho .* dx
    cs = sqrt(GAMMA * CGS.k_B * T_lo / MU)
    @printf("\n== T_lo = %.0e K: h in [%.3g, %.3g], tau_c in [%.3g, %.3g], Mach %.2f, tau_pulse v/c = %.3g ==\n",
            T_lo, extrema(h)..., extrema(tau)..., V0 / cs, KAPPA * RHO0 * W * V0 / CGS.c_light)
    t0 = time()
    ref = score(pulse(; T_lo, phi0=1.0, cfl=CR_REF))
    @printf("  reference chat = c: %d steps, %d sub-steps (%.0f s); centroid at x - v0 t = %.2f cm\n",
            ref.steps, ref.sub, time() - t0, centroid(ref.T, T_lo))
    T_lo == 1e6 && (refs[:moving] = ref)
    for (lab, phi0) in (("v0/c", V0 / CGS.c_light), ("2 v0/c", 2V0 / CGS.c_light), ("", 0.01))
        ps = NU1 * SIGSUM .* h ./ (1 .+ h)
        ok = ps .< phi0
        thr = [ok[i] ? 2 * (SIGSUM - phi0 / NU1) * phi0 / ((1 + h[i]) * (phi0 - ps[i])) : Inf for i in 1:n]
        hp = maximum(h)
        @printf("  phi0 = %s %.4g: largest phi_* = %.3g, cells with phi_* < phi0 %d/%d, above the eq:tau_one threshold %d/%d (finite thresholds in [%.3g, %.3g]); h/phi0 at the peak %.3g\n",
                lab, phi0, maximum(ps), count(ok), n, count(tau .>= thr), n,
                minimum(thr[isfinite.(thr)]; init=Inf), maximum(thr[isfinite.(thr)]; init=-Inf), hp / phi0)
        a = score(pulse(; T_lo, phi0))
        b = score(pulse(; T_lo, phi0, lim=StabilityBound()))
        for (name, r) in (("RSLA", a), ("ARSA", b))
            @printf("    %s: L1/pulse %8.3g%%, centroid lag %7.2f cm, phi in [%.3g, %.3g], cells at phi = 1: %3d/%d, %d steps, %d sub-steps\n",
                    name, 100lp(r.T, ref.T, T_lo), centroid(ref.T, T_lo) - centroid(r.T, T_lo),
                    extrema(r.phi)..., count(r.phi .>= 0.999), n, r.steps, r.sub)
        end
        for i in 1:n
            push!(profiles, (T_lo, V0, TMAX, phi0, xs[i], T[i], ref.T[i], a.T[i], b.T[i], b.phi[i]))
        end
    end
end
# The same pulse at rest: the gas temperature diffuses with D = c h/(3 kappa rho), and under the RSLA
# the radiation's heat capacity grows by 1/phi0, which slows the diffusion by 1/(1 + h/phi0).
let T_lo = 1e6, phi0 = 0.01, tmax = 4.8e-4
    T = Tprof.(xs, T_lo); rho = rhoprof.(xs, T_lo)
    h = 4 * CGS.a_rad .* T .^ 3 ./ (rho * CV)
    i = argmax(T)
    D = CGS.c_light * h[i] / (3 * KAPPA * rho[i])
    @printf("\n== static: T_lo = %.0e K, v = 0, phi0 = %g, t = %.3g s: at the peak D = %.3g cm^2/s, sqrt(D t) = %.3g cm (width %g cm), h/phi0 = %.3g ==\n",
            T_lo, phi0, tmax, D, sqrt(D * tmax), W, h[i] / phi0)
    t0 = time()
    ref = score(pulse(; T_lo, phi0=1.0, cfl=CR_REF, v0=0.0); tmax)
    refs[:static] = ref
    @printf("  reference chat = c: %d steps, %d sub-steps (%.0f s); change of the solution %.3g%% of the pulse; peak T %.4g K (initial %.4g K)\n",
            ref.steps, ref.sub, time() - t0, 100lp(T, ref.T, T_lo), maximum(ref.T), maximum(T))
    a = score(pulse(; T_lo, phi0, v0=0.0); tmax)
    b = score(pulse(; T_lo, phi0, lim=StabilityBound(), v0=0.0); tmax)
    for (name, r) in (("RSLA", a), ("ARSA", b))
        @printf("    %s: L1/pulse %8.3g%% (%.3g of the change), peak T %.4g K, phi in [%.3g, %.3g], cells at phi = 1: %3d/%d, %d steps, %d sub-steps\n",
                name, 100lp(r.T, ref.T, T_lo), sum(abs, r.T .- ref.T) / sum(abs, ref.T .- T),
                maximum(r.T), extrema(r.phi)..., count(r.phi .>= 0.999), n, r.steps, r.sub)
    end
    for i in 1:n
        push!(profiles, (T_lo, 0.0, tmax, phi0, xs[i], T[i], ref.T[i], a.T[i], b.T[i], b.phi[i]))
    end
end

open(joinpath(@__DIR__, "..", "data", "diffusion_pulse_profiles.csv"), "w") do io
    println(io, "T_lo,v0,t,phi0,x,T_initial,T_ref,T_rsla,T_arsa,phi_arsa")
    foreach(r -> println(io, join(r, ",")), profiles)
end
println("wrote data/diffusion_pulse_profiles.csv")

# The ARSA with |F| clipped to cE before the exchange update, the full realizability repair that codes
# commonly apply there (T_lo = 1e6 K, phi0 = 0.01). This redefines a Honghu.jl method, so it runs last,
# through invokelatest; the references above are unclipped.
@eval Honghu function floor_rad_energy!(R::AbstractArray{T,4}, U::AbstractArray{T,4},
                                         chat::AbstractArray{T,4}, phys::RadiationPhysics{NG,M,T},
                                         g::CartesianGrid{T}, floors::Floors{T}) where {NG,M,T}
    amend_rad_state!(R, phys, g, floors)
end
println("\n== ARSA with |F| clipped to cE before the exchange update, T_lo = 1e6 K, phi0 = 0.01 ==")
for (lab, v0, tmax, ref) in (("static", 0.0, 4.8e-4, refs[:static]), ("moving", V0, TMAX, refs[:moving]))
    r = Base.invokelatest(() -> score(pulse(; T_lo=1e6, phi0=0.01, lim=StabilityBound(), v0); tmax))
    @printf("  %s: L1/pulse %8.3g%%, phi in [%.3g, %.3g], %d sub-steps\n", lab, 100lp(r.T, ref.T, 1e6),
            extrema(r.phi)..., r.sub)
end
