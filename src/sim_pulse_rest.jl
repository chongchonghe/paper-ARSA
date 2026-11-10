# A weak pulse at rest (D = 0, the left edge of Fig. 1), in the upper part of the amber region: the pulse of
# sim_pulse_drift.jl with amplitude AMP = 0.1 instead of 1, T = T_lo [1 + AMP exp(-x^2/(2 w^2))], and phi0/P = 3 at
# the peak, so that h = 0.61-0.89 phi0 over the whole pulse and the RSLA slows the diffusion of all of it, by
# (1 + h)/(1 + h/phi0), about a factor of 2, with no drift. (Until 2026-11-08 the full pulse, AMP = 1, whose wings
# have h = 0.055 phi0.) The pulse diffuses on t_diff ~ 3 kappa rho w^2/(h c), so the run is long: to NW t_w,
# t_w = 3 kappa rho0 w^2/c, when the exact peak has fallen to about a third. kappa = 30 cm^2/g keeps the cells thick enough for the ARSA to
# reach phi = 1 (tau_c = 36-72 against tau_thresh <= 28) at the least cost: the radiation substep does not depend on
# kappa, so the number of substeps scales with it. Without advection the box only needs to hold the diffused pulse:
# L = 768 cm on NX = 384 cells (dx = 2 cm, as in sim_pulse_drift.jl).
# The RSLA and the ARSA at phi0 = 0.01, 3D slab NX x 2 x 2, C_r = 0.9 sqrt3, against the exact solution (the enthalpy
# diffusion of sim_pulse_drift.jl in column mass, on NM cells, without the shift by v t).
#
# Prints T_lo, P, h, tau_c and tau_thresh (eq:tau_one, beta = 1 and nu0 = 1/sqrt3 as Table 3, and with beta = 0.9 as
# the limiter), and for each run the smallest phi, the steps, and at every output the peak of (T - T_lo)/T_lo and the
# L1 error of T - T_lo relative to the exact pulse. Writes data/pulse_rest.json: phi0, kappa, w, L, t_w, the output
# times t [s], the cell centres x [cm], and T [K] of the exact solution, the RSLA and the ARSA at every output (the
# first is the initial pulse), and amp; run_all.sh saves the printout in data/pulse_rest.txt. Two single-threaded workers.
#
#   ./run_all.sh pulse_rest      (or: julia --project=. src/sim_pulse_rest.jl)
using Distributed, Printf
addprocs(2; exeflags=`--project=$(Base.active_project()) -t 1`)

@everywhere begin
using Honghu
const CGS = Honghu.CGS
const MU, GAMMA, RHO0, W = 2.33 * CGS.m_u, 5 / 3, 1.2, 24.0     # as in sim_diffusion_pulse.jl
const CV = CGS.k_B / (MU * (GAMMA - 1))
const PHI0, R, KAP, L, NX, NT = 0.01, 3, 30.0, 768.0, 384, 2   # phi0, phi0/P at the peak, kappa, box, cells
const AMP = 0.1                                                 # amplitude of the pulse
const CR = 0.9 * sqrt(3)
const TW = 3 * KAP * RHO0 * W^2 / CGS.c_light                    # t_w
const NW = 1000                                                 # run to NW t_w
const TOUTS = [0.0, 250TW, 500TW, 1000TW]
const XS = [-L / 2 + (i - 0.5) * L / NX for i in 1:NX]

Tprof(x, T_lo) = T_lo * (1 + AMP * exp(-x^2 / (2 * W^2)))
rhoprof(x, T_lo) = (T = Tprof(x, T_lo);
    RHO0 * T_lo / T + (CGS.a_rad * MU / (3 * CGS.k_B)) * (T_lo^4 / T - T^3))   # gas + radiation pressure balance
Ppeak(T_lo) = CGS.a_rad * Tprof(0.0, T_lo)^3 * MU / (rhoprof(0.0, T_lo) * CGS.k_B)

function pulse(T_lo, lim)
    rad = RadiationPhysics(; c=CGS.c_light, chat=PHI0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           beta_order=1, chat_limiter=lim, opacity=ConstantOpacity(KAP))
    sim = Simulation(; n=(NX, NT, NT), xlo=(-L / 2, 0.0, 0.0), xhi=(L / 2, NT * L / NX, NT * L / NX),
                     eos=IdealGas(; gamma=GAMMA, cv=CV), radiation=rad,
                     boundaries=Boundaries(Periodic()), cfl=0.8, cfl_radiation=CR)
    set_initial_conditions!(sim) do x, y, z
        T = Tprof(x, T_lo)
        (rho=rhoprof(x, T_lo), T=T, Erad=CGS.a_rad * T^4, Fx=0.0)
    end
    sim
end

function run_case(T_lo, lim)
    s = pulse(T_lo, lim)
    g = s.grid
    phimin = Ref(1.0)
    cb = s -> (phimin[] = min(phimin[], minimum(view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])) / s.radiation.c); true)
    Ts = [gas_temperature(s)[:, 1, 1]]
    e0 = conserved_totals(s); dE = [0.0]              # total energy in the periodic box
    for t in TOUTS[2:end]
        run!(s; tmax=t, verbose=false, callback=cb)
        push!(Ts, gas_temperature(s)[:, 1, 1]); push!(dE, conserved_totals(s).energy_total / e0.energy_total - 1)
    end
    (; Ts, phimin=phimin[], steps=s.diag.steps, sub=s.diag.rad_substeps, dE, erad0=e0.energy_rad / e0.energy_total)
end
end

function find_Tlo(target)                             # as in sim_pulse_drift.jl
    lo, hi = 1e4, 1e7
    for _ in 1:100
        mid = sqrt(lo * hi)
        Ppeak(mid) < target ? (lo = mid) : (hi = mid)
    end
    sqrt(lo * hi)
end

const NM = 6144                                       # cells in m of the exact solution (dm as in sim_pulse_drift.jl)
interp(xq, xs, ys) = (i = clamp(searchsortedlast(xs, xq), 1, length(xs) - 1);
                      ys[i] + (ys[i+1] - ys[i]) * clamp((xq - xs[i]) / (xs[i+1] - xs[i]), 0, 1))
rhoT(T, T_lo) = RHO0 * T_lo / T + (CGS.a_rad * MU / (3 * CGS.k_B)) * (T_lo^4 / T - T^3)
hT(T, T_lo) = GAMMA * CV * T + 4 * CGS.a_rad * T^4 / (3 * rhoT(T, T_lo))  # enthalpy per mass
function exact(T_lo)                                  # as in sim_pulse_drift.jl, at rest: T on XS at every output
    xf = collect(range(-L / 2, L / 2, length=64NM + 1))
    mf = [0.0; cumsum((rhoprof.(xf[1:end-1], T_lo) .+ rhoprof.(xf[2:end], T_lo)) ./ 2 .* (L / 64NM))]
    dm = mf[end] / NM
    T = [Tprof(interp((j - 0.5) * dm, mf, xf), T_lo) for j in 1:NM]
    h, K, t = hT.(T, T_lo), CGS.c_light / (3 * KAP), 0.0
    dhdT(T) = (hT(T * (1 + 1e-6), T_lo) - hT(T, T_lo)) / (1e-6 * T)
    profs = []
    for tout in TOUTS
        while t < tout
            E = CGS.a_rad .* T .^ 4
            dt = min(0.4 * dm^2 / (K * maximum(4 .* CGS.a_rad .* T .^ 3 ./ dhdT.(T))), tout - t)
            h .+= dt * K .* (circshift(E, 1) .- 2 .* E .+ circshift(E, -1)) ./ dm^2
            for _ in 1:3; T .-= (hT.(T, T_lo) .- h) ./ dhdT.(T); end
            t += dt
        end
        xe = [0.0; cumsum(dm ./ rhoT.(T, T_lo))]; xe .-= xe[NM ÷ 2 + 1]
        xc = (xe[1:end-1] .+ xe[2:end]) ./ 2
        push!(profs, [interp(x, xc, T) for x in XS])
    end
    profs
end

T_lo = find_Tlo(PHI0 / R)
T0 = Tprof.(XS, T_lo); rho = rhoprof.(XS, T_lo)
h = 4 * CGS.a_rad .* T0 .^ 3 ./ (rho * CV)
tauc = KAP .* rho .* (L / NX)
thr1 = 2 * (sqrt(3) - PHI0 * sqrt(3)) * PHI0 ./ ((1 .+ h) .* (PHI0 .- h ./ (1 .+ h)))               # Table 3
ps9 = 0.9 .* h ./ (1 .+ 0.9 .* h)
thr9 = 2 * (sqrt(3) - 0.9 * PHI0 / (0.9 / sqrt(3))) * PHI0 ./ ((1 .+ 0.9 .* h) .* (PHI0 .- ps9))    # the limiter
@printf("amplitude %g; phi0 = %g, phi0/P at the peak %g: T_lo = %.4g K, P at the peak %.3g, h in [%.3g, %.3g] (h/phi0 up to %.3f), kappa = %g cm^2/g\n",
        AMP, PHI0, R, T_lo, Ppeak(T_lo), extrema(h)..., maximum(h) / PHI0, KAP)
@printf("tau_c in [%.3g, %.3g]; tau_thresh in [%.3g, %.3g] (beta = 1, nu0 = 1/sqrt3, Table 3) and [%.3g, %.3g] (beta = 0.9, the limiter)\n",
        extrema(tauc)..., extrema(thr1)..., extrema(thr9)...)
@printf("box %g cm, %d cells; t_w = 3 kappa rho0 w^2/c = %.4g s; outputs at %s t_w\n", L, NX, TW, string(round.(Int, TOUTS ./ TW)))
t1 = time()
fut = [remotecall(run_case, w, T_lo, lim) for (w, lim) in zip(workers(), (ConstantChat(), StabilityBound()))]
ex = exact(T_lo)
rsla, arsa = fetch.(fut)
@printf("(%.0f s)\n", time() - t1)
for (name, r) in (("RSLA", rsla), ("ARSA", arsa))
    @printf("%s: phi >= %.3g, %d steps, %d sub-steps\n", name, r.phimin, r.steps, r.sub)
end
@printf("total energy E_gas + E_rad in the periodic box (exactly conserved by the RHD equations; E_rad/E_tot = %.3g at the start): relative change at %d t_w RSLA %+.3e, ARSA %+.3e\n",
        rsla.erad0, NW, rsla.dE[end], arsa.dE[end])
pk(T) = (maximum(T) - T_lo) / (AMP * T_lo)
l1(T, e) = sum(abs, T .- e) / sum(abs, e .- T_lo)
for k in eachindex(TOUTS)
    @printf("  t = %5d t_w: peak (T - T_lo)/(AMP T_lo) exact %.4f, RSLA %.4f, ARSA %.4f; L1 of T - T_lo RSLA %.4g, ARSA %.4g\n",
            round(Int, TOUTS[k] / TW), pk(ex[k]), pk(rsla.Ts[k]), pk(arsa.Ts[k]), l1(rsla.Ts[k], ex[k]), l1(arsa.Ts[k], ex[k]))
end
@printf("Table 3: h max %.3g, tau_c in [%.3g, %.3g], tau_thresh (beta = 1) in [%.3g, %.3g]; at %d t_w L1 RSLA %.4g, ARSA %.4g, ratio %.3g\n",
        maximum(h), extrema(tauc)..., extrema(thr1)..., NW, l1(rsla.Ts[end], ex[end]), l1(arsa.Ts[end], ex[end]),
        l1(rsla.Ts[end], ex[end]) / l1(arsa.Ts[end], ex[end]))

json(io, x::Real) = print(io, isfinite(x) ? @sprintf("%.7g", x) : "null")
json(io, x::Symbol) = print(io, '"', x, '"')
json(io, x::AbstractVector) = (print(io, '['); for (i, v) in enumerate(x); i > 1 && print(io, ','); json(io, v); end; print(io, ']'))
json(io, x::NamedTuple) = (print(io, '{'); for (i, (k, v)) in enumerate(pairs(x)); i > 1 && print(io, ",\n"); json(io, k); print(io, ':'); json(io, v); end; print(io, '}'))
open(joinpath(@__DIR__, "..", "data", "pulse_rest.json"), "w") do io
    json(io, (phi0=PHI0, phi0_P=R, amp=AMP, kappa=KAP, w=W, L=L, T_lo=T_lo, t_w=TW, t=TOUTS, x=XS,
              T_exact=ex, T_rsla=rsla.Ts, T_arsa=arsa.Ts, dEtot_rsla=rsla.dE, dEtot_arsa=arsa.dE))
end
println("wrote data/pulse_rest.json")
