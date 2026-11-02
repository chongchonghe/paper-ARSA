# The pulse of sim_pulse_drift.jl in the corner of small P and large D of Fig. 1 (paper Sec. 5.3), run until it has
# diffused away: since a structure drifts by about 3 D/phi of its size over its diffusive life, about 3 gamma D/h
# widths traveled, whatever h, the drift test needs a short life. Here D = kappa rho w v/c = 0.1 = 10 phi0, with
# v = 1e7 cm/s and kappa = 0.1 c/(rho w v) = 10.4 cm^2/g (cells of optical depth 12-25 on 512 cells), and
# phi0/P = 30 at the peak of the pulse (h = 8.9e-4 = 0.09 phi0 there), so that the pulse diffuses over its width
# in a few hundred widths traveled. The RSLA and the ARSA at phi0 = 0.01 run to NOUT DW = 2048 widths traveled
# (20 crossings of the box), with output every DW = 64 widths, against the exact solution (the pulse at rest,
# diffusing, shifted by v t; as in sim_pulse_drift.jl). Each run appends every output, as it reaches it, to
# logs/pulse_drift_long_<run>.jsonl. Writes data/pulse_drift_long.json, in the layout of a row of
# data/pulse_drift.json plus kappa and v0 at the top level; run_all.sh saves the printout in
# data/pulse_drift_long.txt. About 40 min on two single-threaded workers.
#
#   ./run_all.sh pulse_drift_long      (or: julia --project=. src/sim_pulse_drift_long.jl)
using Distributed, Printf
addprocs(2; exeflags=`--project=$(Base.active_project()) -t 1`)
const SRC = @__DIR__
@everywhere workers() const SRC = $SRC
@everywhere begin
    src = read(joinpath(SRC, "sim_diffusion_pulse.jl"), String)
    include_string(Main, src[1:first(findfirst("\nn = 128\n", src))])    # constants, Tprof, rhoprof
end

@everywhere begin
const PHI0, R, VEL = 0.01, 30, 1e7                    # phi0, phi0/P at the peak, flow speed [cm/s]
const KAP = 0.1 * CGS.c_light / (RHO0 * W * VEL)      # D = kappa rho w v/c = 0.1
const DW, NOUT, NSAMPLE, NX, WIN = 64, 32, 8, 512, 8  # output every DW widths, NOUT of them; centroid window [w]
const TOUT = DW * W / VEL
const XS = [-L / 2 + (i - 0.5) * L / NX for i in 1:NX]
const LOGS = joinpath(SRC, "..", "logs")

Ppeak(T_lo) = CGS.a_rad * Tprof(0.0, T_lo)^3 * MU / (rhoprof(0.0, T_lo) * CGS.k_B)
function find_Tlo(target)                             # as in sim_pulse_drift.jl
    lo, hi = 1e4, 1e7
    for _ in 1:100
        mid = sqrt(lo * hi)
        Ppeak(mid) < target ? (lo = mid) : (hi = mid)
    end
    sqrt(lo * hi)
end

function pulse_long(T_lo, lim)                        # pulse() of sim_diffusion_pulse.jl with kappa = KAP, v = VEL
    rad = RadiationPhysics(; c=CGS.c_light, chat=PHI0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           beta_order=1, chat_limiter=lim, opacity=ConstantOpacity(KAP))
    sim = Simulation(; n=(NX, NT, NT), xlo=(-L / 2, 0.0, 0.0), xhi=(L / 2, NT * L / NX, NT * L / NX),
                     eos=IdealGas(; gamma=GAMMA, cv=CV), radiation=rad, boundaries=Boundaries(Periodic()),
                     cfl=0.8, cfl_radiation=CR)
    set_initial_conditions!(sim) do x, y, z
        T = Tprof(x, T_lo)
        (rho=rhoprof(x, T_lo), vx=VEL, T=T, Erad=CGS.a_rad * T^4, Fx=0.0)
    end
    sim
end

function centroid(T, T_lo)                            # as in sim_pulse_drift.jl
    w = T .- T_lo
    c = atan(sum(w .* sin.(2pi .* XS ./ L)), sum(w .* cos.(2pi .* XS ./ L))) * L / (2pi)
    d = mod.(XS .- c .+ L / 2, L) .- L / 2
    near = abs.(d) .< WIN * W
    c + sum(d[near] .* w[near]) / sum(w[near])
end

vec_json(v) = "[" * join((@sprintf("%.7g", x) for x in v), ",") * "]"
function track(name, T_lo, lim)                       # one run; appends every output to logs/pulse_drift_long_<name>.jsonl
    s = pulse_long(T_lo, lim)
    g = s.grid
    prof(s) = (T=gas_temperature(s)[:, 1, 1], Tr=radiation_temperature(s)[:, 1, 1])
    p = prof(s)
    xb = centroid(p.T, T_lo); x = xb
    pos, profs, phimin = [x], [p], 1.0
    logf = joinpath(LOGS, "pulse_drift_long_$name.jsonl")
    rm(logf; force=true)
    for k in 1:NOUT, j in 1:NSAMPLE
        run!(s; tmax=(k - 1 + j / NSAMPLE) * TOUT, verbose=false)
        p = prof(s)
        c = centroid(p.T, T_lo)
        x += mod(c - xb + L / 2, L) - L / 2; xb = c
        phimin = min(phimin, minimum(view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])) / s.radiation.c)
        j == NSAMPLE || continue
        push!(pos, x); push!(profs, p)
        open(logf, "a") do io                         # one writer per file
            println(io, "{\"widths\":", k * DW, ",\"time\":\"", Libc.strftime("%H:%M:%S", time()), "\",\"x\":",
                    @sprintf("%.7g", x), ",\"phi_min\":", @sprintf("%.4g", phimin), ",\"T\":", vec_json(p.T),
                    ",\"Tr\":", vec_json(p.Tr), "}")
        end
    end
    (; pos, profs, phimin, steps=s.diag.steps, sub=s.diag.rad_substeps)
end
end

# the exact solution, as exact() of sim_pulse_drift.jl with kappa = KAP and v = VEL
const NM = 4096
interp(xq, xs, ys) = (i = clamp(searchsortedlast(xs, xq), 1, length(xs) - 1);
                      ys[i] + (ys[i+1] - ys[i]) * clamp((xq - xs[i]) / (xs[i+1] - xs[i]), 0, 1))
rhoT(T, T_lo) = RHO0 * T_lo / T + (CGS.a_rad * MU / (3 * CGS.k_B)) * (T_lo^4 / T - T^3)
hT(T, T_lo) = GAMMA * CV * T + 4 * CGS.a_rad * T^4 / (3 * rhoT(T, T_lo))
function exact(T_lo)
    xf = collect(range(-L / 2, L / 2, length=64NM + 1))
    mf = [0.0; cumsum((rhoprof.(xf[1:end-1], T_lo) .+ rhoprof.(xf[2:end], T_lo)) ./ 2 .* (L / 64NM))]
    dm = mf[end] / NM
    T = [Tprof(interp((j - 0.5) * dm, mf, xf), T_lo) for j in 1:NM]
    h, K, t = hT.(T, T_lo), CGS.c_light / (3 * KAP), 0.0
    dhdT(T) = (hT(T * (1 + 1e-6), T_lo) - hT(T, T_lo)) / (1e-6 * T)
    profs = []
    for k in 0:NOUT
        while t < k * TOUT
            E = CGS.a_rad .* T .^ 4
            dt = min(0.4 * dm^2 / (K * maximum(4 .* CGS.a_rad .* T .^ 3 ./ dhdT.(T))), k * TOUT - t)
            h .+= dt * K .* (circshift(E, 1) .- 2 .* E .+ circshift(E, -1)) ./ dm^2
            for _ in 1:3; T .-= (hT.(T, T_lo) .- h) ./ dhdT.(T); end
            t += dt
        end
        xe = [0.0; cumsum(dm ./ rhoT.(T, T_lo))]; xe .-= xe[NM ÷ 2 + 1]
        xc, s = (xe[1:end-1] .+ xe[2:end]) ./ 2, VEL * k * TOUT
        push!(profs, [interp(mod(x - s + L / 2, L) - L / 2, xc, T) for x in XS])
    end
    profs
end

T_lo = find_Tlo(PHI0 / R)
mkpath(LOGS)
jobs = [(:rsla, ConstantChat()), (:arsa, StabilityBound())]
out = pmap(j -> track(j[1], T_lo, j[2]), jobs)
rsla, arsa = out
ex = exact(T_lo)

T0 = Tprof.(XS, T_lo); h = 4 * CGS.a_rad .* T0 .^ 3 ./ (rhoprof.(XS, T_lo) * CV)
mach = VEL / sqrt(GAMMA * CGS.k_B * T_lo / MU)
xex = [VEL * k * TOUT for k in 0:NOUT]
lag(a) = (xex .- a.pos) ./ W
pk(T) = (maximum(T) - T_lo) / T_lo
@printf("phi0 = %g, D = kappa rho w v/c = %.3g (kappa = %.4g cm^2/g, v = %.3g cm/s), phi0/P at the peak %d: T_lo = %.4g K, P at the peak %.3g, h in [%.3g, %.3g], Mach %.2f, tau_c in [%.3g, %.3g]\n",
        PHI0, KAP * RHO0 * W * VEL / CGS.c_light, KAP, VEL, R, T_lo, Ppeak(T_lo), extrema(h)..., mach,
        extrema(KAP .* rhoprof.(XS, T_lo) .* (L / NX))...)
@printf("%d widths traveled, output every %d w (%.3g s); %d cells; ARSA phi >= %.3g; steps RSLA %d, ARSA %d\n",
        NOUT * DW, DW, TOUT, NX, arsa.phimin, rsla.steps, arsa.steps)
for k in 1:NOUT
    @printf("  %4d w traveled: centroid lag RSLA %8.3f w, ARSA %7.3f w; peak (T - T_lo)/T_lo exact %.4f, RSLA %.4f, ARSA %.4f\n",
            k * DW, lag(rsla)[k+1], lag(arsa)[k+1], pk(ex[k+1]), pk(rsla.profs[k+1].T), pk(arsa.profs[k+1].T))
end

row = (phi0_P=R, T_lo=T_lo, P_peak=Ppeak(T_lo), h_peak=maximum(h), mach=mach, arsa_phi_min=arsa.phimin,
       widths=collect(0:DW:NOUT * DW), x_exact=xex, x_rsla=rsla.pos, x_arsa=arsa.pos,
       lag_rsla_w=lag(rsla), lag_arsa_w=lag(arsa),
       T_exact=ex, T_rsla=[p.T for p in rsla.profs], T_arsa=[p.T for p in arsa.profs],
       Tr_rsla=[p.Tr for p in rsla.profs], Tr_arsa=[p.Tr for p in arsa.profs])
json(io, x::Real) = print(io, isfinite(x) ? @sprintf("%.7g", x) : "null")
json(io, x::Symbol) = print(io, '"', x, '"')
json(io, x::AbstractVector) = (print(io, '['); for (i, v) in enumerate(x); i > 1 && print(io, ','); json(io, v); end; print(io, ']'))
json(io, x::NamedTuple) = (print(io, '{'); for (i, (k, v)) in enumerate(pairs(x)); i > 1 && print(io, ",\n"); json(io, k); print(io, ':'); json(io, v); end; print(io, '}'))
open(joinpath(SRC, "..", "data", "pulse_drift_long.json"), "w") do io
    json(io, (phi0=PHI0, D=KAP * RHO0 * W * VEL / CGS.c_light, kappa=KAP, v0=VEL, w=W, L=L, t_out=TOUT, x=XS,
              runs=[row]))
end
println("\nwrote data/pulse_drift_long.json")
