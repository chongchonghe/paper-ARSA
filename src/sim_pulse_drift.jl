# Under the RSLA a temperature structure carried by the flow falls behind it where phi0 is not >> D, even where
# phi0 >> P (paper Sec. 5.3, Fig. 7). This test follows the advected pulse of sim_diffusion_pulse.jl
# (kappa rho w v/c = D = 14.4, dynamic diffusion) for NCROSS crossings of the periodic box at phi0 = 0.01
# (chat = 10 v0), with T_lo chosen so that phi0/P at the peak of the pulse is each of PHI0_P (rho, kappa, w and v
# fixed); P = a_r T^4/(rho k_B T/mu) is the pressure ratio of the paper, P = h/(4 (gamma - 1)). For each phi0/P,
# the reference (chat = c), the RSLA and the ARSA run as separate tasks, twelve in all. Where h at the peak
# exceeds phi0 (phi0/P = 1), phi of the ARSA saturates below 1 at the peak; there the ARSA also runs at
# phi0 = PHI0_HI, above h, as Fig. 1 prescribes (the run arsa_hi; json keys *_arsa_hi). Each task records the
# position of the pulse NSAMPLE times per crossing: the centroid of T - T_lo in the periodic box (a circular mean,
# then the mean of x over the box centred on it), unwrapped into the real position x, which starts at 0. Writes
# data/pulse_drift.json: for each phi0/P, T_lo, P and h at the peak, the Mach number, the smallest phi of the ARSA,
# and at every integer crossing k = 0..NCROSS the real centroid x [cm] of each run, the lags (x_ref - x)/w of the
# RSLA and the ARSA in pulse widths w = 24 cm, and the profiles of T and T_rad [K] on the cells XS of the box;
# run_all.sh saves the printout in data/pulse_drift.txt. About an hour with 12 threads.
#
#   ./run_all.sh pulse_drift      (or: julia --project=. -t 12 src/sim_pulse_drift.jl)
const SRC = @__DIR__
src = read(joinpath(SRC, "sim_diffusion_pulse.jl"), String)
include_string(@__MODULE__, src[1:first(findfirst("\nn = 128\n", src))])    # constants, Tprof, rhoprof, pulse

const PHI0_P = (1, 10, 100, 1000)                     # phi0/P at the peak of the pulse
const PHI0, NCROSS, NSAMPLE, NX = 0.01, 32, 8, 128
const PHI0_HI = 0.03                                  # the ARSA where h at the peak > PHI0
const TCROSS = L / V0                                 # one crossing of the box
const XS = [-L / 2 + (i - 0.5) * L / NX for i in 1:NX]

Ppeak(T_lo) = CGS.a_rad * Tprof(0.0, T_lo)^3 * MU / (rhoprof(0.0, T_lo) * CGS.k_B)
function find_Tlo(target)                             # bisection in log T_lo; P at the peak grows with T_lo
    lo, hi = 1e4, 1e7
    for _ in 1:100
        mid = sqrt(lo * hi)
        Ppeak(mid) < target ? (lo = mid) : (hi = mid)
    end
    sqrt(lo * hi)
end

function centroid(T, T_lo)                            # in the periodic box
    w = T .- T_lo
    c = atan(sum(w .* sin.(2pi .* XS ./ L)), sum(w .* cos.(2pi .* XS ./ L))) * L / (2pi)
    d = mod.(XS .- c .+ L / 2, L) .- L / 2
    c + sum(d .* w) / sum(w)
end

function track(T_lo, phi0, cfl, lim)                  # one run: real centroid and profiles at every crossing
    s = pulse(; T_lo, phi0, cfl, lim)
    g = s.grid
    prof(s) = (T=gas_temperature(s)[:, 1, 1], Tr=radiation_temperature(s)[:, 1, 1])
    p = prof(s)
    xb = centroid(p.T, T_lo); x = xb                  # centroid in the box, and unwrapped
    pos, profs, phimin = [x], [p], 1.0
    for k in 1:NCROSS, j in 1:NSAMPLE
        run!(s; tmax=(k - 1 + j / NSAMPLE) * TCROSS, verbose=false)
        p = prof(s)
        c = centroid(p.T, T_lo)
        x += mod(c - xb + L / 2, L) - L / 2; xb = c
        phimin = min(phimin, minimum(view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])) / s.radiation.c)
        j == NSAMPLE && (push!(pos, x); push!(profs, p))
    end
    (; pos, profs, phimin, steps=s.diag.steps, sub=s.diag.rad_substeps)
end

const RUNS = ((:ref, 1.0, CR_REF, ConstantChat()), (:rsla, PHI0, CR, ConstantChat()),
              (:arsa, PHI0, CR, StabilityBound()))
T_los = Dict(r => find_Tlo(PHI0 / r) for r in PHI0_P)
hpeak(T_lo) = 4 * (GAMMA - 1) * Ppeak(T_lo)
tasks = Dict((r, name) => Threads.@spawn(track(T_los[r], phi0, cfl, lim))
             for r in PHI0_P for (name, phi0, cfl, lim) in RUNS)
for r in PHI0_P
    hpeak(T_los[r]) > PHI0 && (tasks[(r, :arsa_hi)] = Threads.@spawn(track(T_los[r], PHI0_HI, CR, StabilityBound())))
end
res = Dict(key => fetch(t) for (key, t) in tasks)

@printf("phi0 = %g, D = kappa rho w v/c = %.3g, %d crossings of %g cm (%.3g s each, %.3g w); w = %g cm\n",
        PHI0, KAPPA * RHO0 * W * V0 / CGS.c_light, NCROSS, L, TCROSS, L / W, W)
rows = []
for r in PHI0_P
    T_lo = T_los[r]
    ref, rsla, arsa = res[(r, :ref)], res[(r, :rsla)], res[(r, :arsa)]
    T0 = Tprof.(XS, T_lo); h = 4 * CGS.a_rad .* T0 .^ 3 ./ (rhoprof.(XS, T_lo) * CV)
    mach = V0 / sqrt(GAMMA * CGS.k_B * T_lo / MU)
    lag(a) = (ref.pos .- a.pos) ./ W
    @printf("\n== phi0/P at the peak %g: T_lo = %.4g K, P at the peak %.3g, h in [%.3g, %.3g], Mach %.2f; ARSA phi >= %.3g; steps ref %d, RSLA %d, ARSA %d ==\n",
            r, T_lo, Ppeak(T_lo), extrema(h)..., mach, arsa.phimin, ref.steps, rsla.steps, arsa.steps)
    hi = get(res, (r, :arsa_hi), nothing)
    hi === nothing || @printf("  ARSA at phi0 = %g: phi >= %.3g, %d steps\n", PHI0_HI, hi.phimin, hi.steps)
    for k in (1, 2, 4, 8, 16, 32)
        k > NCROSS && continue
        @printf("  crossing %2d (%6.0f w travelled): centroid lag RSLA %9.3f w, ARSA %7.3f w%s\n",
                k, (ref.pos[k+1] - ref.pos[1]) / W, lag(rsla)[k+1], lag(arsa)[k+1],
                hi === nothing ? "" : @sprintf(", ARSA at phi0 = %g %7.3f w", PHI0_HI, lag(hi)[k+1]))
    end
    row = (phi0_P=r, T_lo=T_lo, P_peak=Ppeak(T_lo), h_peak=maximum(h), mach=mach, arsa_phi_min=arsa.phimin,
           crossing=collect(0:NCROSS), x_ref=ref.pos, x_rsla=rsla.pos, x_arsa=arsa.pos,
           lag_rsla_w=lag(rsla), lag_arsa_w=lag(arsa),
           T_ref=[p.T for p in ref.profs], T_rsla=[p.T for p in rsla.profs], T_arsa=[p.T for p in arsa.profs],
           Tr_ref=[p.Tr for p in ref.profs], Tr_rsla=[p.Tr for p in rsla.profs], Tr_arsa=[p.Tr for p in arsa.profs])
    hi === nothing || (row = (; row..., phi0_arsa_hi=PHI0_HI, arsa_hi_phi_min=hi.phimin, x_arsa_hi=hi.pos,
                                lag_arsa_hi_w=lag(hi), T_arsa_hi=[p.T for p in hi.profs],
                                Tr_arsa_hi=[p.Tr for p in hi.profs]))
    push!(rows, row)
end

json(io, x::Real) = print(io, isfinite(x) ? @sprintf("%.7g", x) : "null")
json(io, x::Symbol) = print(io, '"', x, '"')
json(io, x::AbstractVector) = (print(io, '['); for (i, v) in enumerate(x); i > 1 && print(io, ','); json(io, v); end; print(io, ']'))
json(io, x::NamedTuple) = (print(io, '{'); for (i, (k, v)) in enumerate(pairs(x)); i > 1 && print(io, ",\n"); json(io, k); print(io, ':'); json(io, v); end; print(io, '}'))
open(joinpath(SRC, "..", "data", "pulse_drift.json"), "w") do io
    json(io, (phi0=PHI0, D=KAPPA * RHO0 * W * V0 / CGS.c_light, w=W, L=L, v0=V0, t_cross=TCROSS, x=XS, runs=rows))
end
println("\nwrote data/pulse_drift.json")
