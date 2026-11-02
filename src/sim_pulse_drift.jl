# Under the RSLA a temperature structure carried by the flow drifts behind it (paper Sec. 2.2), even where the gas holds
# most of the heat capacity, h << phi0, the green region of Fig. 1. This test follows the advected pulse of
# sim_diffusion_pulse.jl (kappa rho w v/c = 14.4, dynamic diffusion) for NCROSS crossings of the periodic box at
# phi0 = 0.01 (chat = 10 v0), with T_lo chosen so that h/phi0 at the peak of the pulse is each of HPHI (rho, kappa,
# w and v fixed). The RSLA, the ARSA and the chat = c reference advance together, one box crossing at a time; after
# each crossing the L1 errors of the gas and radiation temperatures relative to the pulse, and the lag of the pulse
# centroid behind the reference, are recorded. The cases run in parallel, one task each. Writes
# data/pulse_drift.csv; run_all.sh saves the printout in data/pulse_drift.txt. About 15 minutes with 6 threads.
#
#   ./run_all.sh pulse_drift      (or: julia --project=. -t 6 src/sim_pulse_drift.jl)
const SRC = @__DIR__
src = read(joinpath(SRC, "sim_diffusion_pulse.jl"), String)
include_string(@__MODULE__, src[1:first(findfirst("\nn = 128\n", src))])    # constants, Tprof, rhoprof, pulse

const HPHI = (0.001, 0.003, 0.01, 0.03, 0.1, 0.3)    # h/phi0 at the peak of the pulse; 0.3 is in the amber region
const PHI0, NCROSS, NX = 0.01, 8, 128
const TCROSS = L / V0                                 # one crossing of the box
const XS = [-L / 2 + (i - 0.5) * L / NX for i in 1:NX]

hpeak(T_lo) = 4 * CGS.a_rad * Tprof(0.0, T_lo)^3 / (rhoprof(0.0, T_lo) * CV)
function find_Tlo(target)                             # bisection in log T_lo; h at the peak grows with T_lo
    lo, hi = 1e4, 1e7
    for _ in 1:100
        mid = sqrt(lo * hi)
        hpeak(mid) < target ? (lo = mid) : (hi = mid)
    end
    sqrt(lo * hi)
end

function drift_case(x)
    T_lo = find_Tlo(x * PHI0)
    lp(a, b) = sum(abs, a .- b) / sum(abs, b .- T_lo)                   # L1 relative to the pulse
    centroid(T, t) = (xc = mod.(XS .- V0 * t .+ L / 2, L) .- L / 2;     # x - v0 t in the periodic box
                      sum(xc .* (T .- T_lo)) / sum(T .- T_lo))
    sims = (ref=pulse(; T_lo, phi0=1.0, cfl=CR_REF), rsla=pulse(; T_lo, phi0=PHI0),
            arsa=pulse(; T_lo, phi0=PHI0, lim=StabilityBound()))
    rows = []
    for k in 1:NCROSS
        t = k * TCROSS
        st = Dict(name => (run!(s; tmax=t, verbose=false);
                           (T=gas_temperature(s)[:, 1, 1], Tr=radiation_temperature(s)[:, 1, 1]))
                  for (name, s) in pairs(sims))
        r, a, b = st[:ref], st[:rsla], st[:arsa]
        push!(rows, (x, T_lo, k, t, lp(a.T, r.T), lp(b.T, r.T), lp(a.Tr, r.Tr), lp(b.Tr, r.Tr),
                     centroid(r.T, t) - centroid(a.T, t), centroid(r.T, t) - centroid(b.T, t)))
    end
    rows
end

tasks = [Threads.@spawn drift_case(x) for x in HPHI]
results = fetch.(tasks)

@printf("phi0 = %g, kappa rho w v/c = %.3g, %d crossings of %g cm (%.3g s each); errors relative to the pulse\n",
        PHI0, KAPPA * RHO0 * W * V0 / CGS.c_light, NCROSS, L, TCROSS)
for (x, rows) in zip(HPHI, results)
    T_lo = rows[1][2]
    T0 = Tprof.(XS, T_lo); h = 4 * CGS.a_rad .* T0 .^ 3 ./ (rhoprof.(XS, T_lo) * CV)
    lag1 = rows[1][9] / L                              # RSLA lag per distance travelled, first crossing
    @printf("\n== h/phi0 at the peak %g: T_lo = %.4g K, h in [%.3g, %.3g], Mach %.2f; RSLA lag per distance after one crossing %.3g = %.2f x h_peak (1 - phi0)/(gamma phi0) (the drift speed of Sec. 2.2) ==\n",
            x, T_lo, extrema(h)..., V0 / sqrt(GAMMA * CGS.k_B * T_lo / MU), lag1, lag1 / (x * (1 - PHI0) / GAMMA))
    for r in rows
        @printf("  crossing %d: L1 T RSLA %7.3g%% ARSA %6.3g%% | L1 Tr RSLA %7.3g%% ARSA %6.3g%% | lag RSLA %7.2f cm ARSA %6.2f cm\n",
                r[3], 100r[5], 100r[6], 100r[7], 100r[8], r[9], r[10])
    end
end

open(joinpath(SRC, "..", "data", "pulse_drift.csv"), "w") do io
    println(io, "h_phi0_peak,T_lo,crossing,t,L1_T_rsla,L1_T_arsa,L1_Tr_rsla,L1_Tr_arsa,lag_rsla,lag_arsa")
    for rows in results, r in rows
        println(io, join(r, ","))
    end
end
println("\nwrote data/pulse_drift.csv")
