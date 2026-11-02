# Under the RSLA a temperature structure carried by the flow falls behind it where phi0 is not >> D, even where
# phi0 >> P (paper Sec. 5.3, Fig. 7). This test follows the advected pulse of sim_diffusion_pulse.jl
# (kappa rho w v/c = D = 14.4, dynamic diffusion), on NX = 512 cells, for NOUT DW = 256 pulse widths traveled
# (6 crossings of the box) at phi0 = 0.01 (chat = 10 v0), with T_lo chosen so that phi0/P at the peak of the pulse is
# each of PHI0_P (rho, kappa, w and v fixed); P = a_r T^4/(rho k_B T/mu) is the pressure ratio of the paper,
# P = h/(4 (gamma - 1)), so h < phi0 at the peak in every row. The RSLA and the ARSA of each phi0/P run as separate
# jobs; the reference is the exact solution (phi = 1, no RSLA). By Galilean invariance it is the pulse at rest,
# shifted by v t. At rest the pulse stays in pressure balance and kappa rho dx >> 1, so in the column mass m
# (dm = rho dx) the enthalpy per mass h = c_p T + (4/3) a T^4/rho obeys dh/dt = (c/(3 kappa)) d^2(a T^4)/dm^2,
# with rho(T) from the pressure balance; exact() solves it explicitly on NM cells in m (converged to 2e-4 in T)
# and maps it back to x. Each run records the position of the pulse NSAMPLE times per DW = 32 widths: the centroid of
# T - T_lo in the periodic box (a circular mean, then the mean of x over the box centred on it), unwrapped into
# the real position x, which starts at 0; that of the exact solution is v t. Writes
# data/pulse_drift.json: for each phi0/P, T_lo, P and h at the peak, the Mach number, the smallest phi of the
# ARSA, and every DW widths traveled (0, 32, ..., 256) the real centroid x [cm] of each run, the lags (v t - x)/w of
# the RSLA and the ARSA in pulse widths w = 24 cm, and the profiles of T [K] of the exact solution and of T and
# T_rad [K] of the runs on the cells XS of the box; run_all.sh saves the printout in data/pulse_drift.txt. The runs
# use single-threaded worker processes, since Honghu.jl spreads its grid loops over the threads of a process;
# about 3 min with 8 workers. logs/pulse_drift.progress shows each output as it is reached.
#
#   ./run_all.sh pulse_drift      (or: julia --project=. src/sim_pulse_drift.jl)
using Distributed, Printf
addprocs(8; exeflags=`--project=$(Base.active_project()) -t 1`)
const SRC = @__DIR__
@everywhere workers() const SRC = $SRC
@everywhere begin
    src = read(joinpath(SRC, "sim_diffusion_pulse.jl"), String)
    include_string(Main, src[1:first(findfirst("\nn = 128\n", src))])    # constants, Tprof, rhoprof, pulse
end

@everywhere begin
const PHI0_P = (3, 30, 300, 3000)                     # phi0/P at the peak of the pulse
const PHI0, DW, NOUT, NSAMPLE, NX = 0.01, 32, 8, 4, 512   # outputs every DW widths, NOUT of them
const TOUT = DW * W / V0                              # the time to travel DW widths
const XS = [-L / 2 + (i - 0.5) * L / NX for i in 1:NX]
const PROGRESS = joinpath(SRC, "..", "logs", "pulse_drift.progress")

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

function track(label, T_lo, phi0, cfl, lim)           # one run: real centroid and profiles at every output
    s = pulse(; T_lo, phi0, cfl, lim, n=NX)
    g = s.grid
    prof(s) = (T=gas_temperature(s)[:, 1, 1], Tr=radiation_temperature(s)[:, 1, 1])
    p = prof(s)
    xb = centroid(p.T, T_lo); x = xb                  # centroid in the box, and unwrapped
    pos, profs, phimin = [x], [p], 1.0
    for k in 1:NOUT, j in 1:NSAMPLE
        run!(s; tmax=(k - 1 + j / NSAMPLE) * TOUT, verbose=false)
        p = prof(s)
        c = centroid(p.T, T_lo)
        x += mod(c - xb + L / 2, L) - L / 2; xb = c
        phimin = min(phimin, minimum(view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])) / s.radiation.c)
        j == NSAMPLE || continue
        push!(pos, x); push!(profs, p)
        f = Base.Filesystem.open(PROGRESS, Base.JL_O_WRONLY | Base.JL_O_CREAT | Base.JL_O_APPEND, 0o644)
        write(f, string(Libc.strftime("%H:%M:%S", time()), "  ", label, ": ", k * DW, " of ", NOUT * DW, " widths\n"))
        close(f)                                      # one O_APPEND write: the lines of the runs do not mix
    end
    (; pos, profs, phimin, steps=s.diag.steps, sub=s.diag.rad_substeps)
end
end

const NM = 8192                                       # cells in m of the exact solution
interp(xq, xs, ys) = (i = clamp(searchsortedlast(xs, xq), 1, length(xs) - 1);
                      ys[i] + (ys[i+1] - ys[i]) * clamp((xq - xs[i]) / (xs[i+1] - xs[i]), 0, 1))
rhoT(T, T_lo) = RHO0 * T_lo / T + (CGS.a_rad * MU / (3 * CGS.k_B)) * (T_lo^4 / T - T^3)    # pressure balance
hT(T, T_lo) = GAMMA * CV * T + 4 * CGS.a_rad * T^4 / (3 * rhoT(T, T_lo))                     # enthalpy per mass
function exact(T_lo)                                  # T of the exact solution on XS at every output
    xf = collect(range(-L / 2, L / 2, length=64NM + 1))   # m(x) of the initial pulse, m = 0 at x = -L/2
    mf = [0.0; cumsum((rhoprof.(xf[1:end-1], T_lo) .+ rhoprof.(xf[2:end], T_lo)) ./ 2 .* (L / 64NM))]
    dm = mf[end] / NM
    T = [Tprof(interp((j - 0.5) * dm, mf, xf), T_lo) for j in 1:NM]
    h, K, t = hT.(T, T_lo), CGS.c_light / (3 * KAPPA), 0.0
    dhdT(T) = (hT(T * (1 + 1e-6), T_lo) - hT(T, T_lo)) / (1e-6 * T)
    profs = []
    for k in 0:NOUT
        while t < k * TOUT
            E = CGS.a_rad .* T .^ 4
            dt = min(0.4 * dm^2 / (K * maximum(4 .* CGS.a_rad .* T .^ 3 ./ dhdT.(T))), k * TOUT - t)
            h .+= dt * K .* (circshift(E, 1) .- 2 .* E .+ circshift(E, -1)) ./ dm^2
            for _ in 1:3; T .-= (hT.(T, T_lo) .- h) ./ dhdT.(T); end    # T from h (Newton)
            t += dt
        end
        xe = [0.0; cumsum(dm ./ rhoT.(T, T_lo))]; xe .-= xe[NM ÷ 2 + 1]     # x of the cell edges, centre at 0
        xc, s = (xe[1:end-1] .+ xe[2:end]) ./ 2, V0 * k * TOUT                # shifted by v t in the box
        push!(profs, [interp(mod(x - s + L / 2, L) - L / 2, xc, T) for x in XS])
    end
    profs
end

const RUNS = ((:rsla, ConstantChat()), (:arsa, StabilityBound()))
T_los = Dict(r => find_Tlo(PHI0 / r) for r in PHI0_P)
jobs = [(r, name, lim) for r in PHI0_P for (name, lim) in RUNS]
mkpath(dirname(PROGRESS)); rm(PROGRESS; force=true)
out = pmap(j -> track("phi0/P = $(j[1]) $(j[2])", T_los[j[1]], PHI0, CR, j[3]), jobs)
res = Dict((j[1], j[2]) => o for (j, o) in zip(jobs, out))

@printf("phi0 = %g, D = kappa rho w v/c = %.3g, %d widths traveled (w = %g cm), output every %d w (%.3g s); box %g cm (%.3g w), %d cells\n",
        PHI0, KAPPA * RHO0 * W * V0 / CGS.c_light, NOUT * DW, W, DW, TOUT, L, L / W, NX)
rows = []
for r in PHI0_P
    T_lo = T_los[r]
    rsla, arsa, ex = res[(r, :rsla)], res[(r, :arsa)], exact(T_lo)
    T0 = Tprof.(XS, T_lo); h = 4 * CGS.a_rad .* T0 .^ 3 ./ (rhoprof.(XS, T_lo) * CV)
    mach = V0 / sqrt(GAMMA * CGS.k_B * T_lo / MU)
    xex = [V0 * k * TOUT for k in 0:NOUT]                  # the centroid of the exact solution, v t
    lag(a) = (xex .- a.pos) ./ W
    @printf("\n== phi0/P at the peak %g: T_lo = %.4g K, P at the peak %.3g, h in [%.3g, %.3g], Mach %.2f; ARSA phi >= %.3g; steps RSLA %d, ARSA %d; peak (T - T_lo)/T_lo of the exact solution %.4f at the start, %.4f at the end ==\n",
            r, T_lo, Ppeak(T_lo), extrema(h)..., mach, arsa.phimin, rsla.steps, arsa.steps,
            (maximum(ex[1]) - T_lo) / T_lo, (maximum(ex[end]) - T_lo) / T_lo)
    for k in 1:NOUT
        @printf("  %3d w traveled: centroid lag RSLA %8.3f w, ARSA %7.3f w\n", k * DW, lag(rsla)[k+1], lag(arsa)[k+1])
    end
    push!(rows, (phi0_P=r, T_lo=T_lo, P_peak=Ppeak(T_lo), h_peak=maximum(h), mach=mach, arsa_phi_min=arsa.phimin,
                 widths=collect(0:DW:NOUT * DW), x_exact=xex, x_rsla=rsla.pos, x_arsa=arsa.pos,
                 lag_rsla_w=lag(rsla), lag_arsa_w=lag(arsa),
                 T_exact=ex, T_rsla=[p.T for p in rsla.profs], T_arsa=[p.T for p in arsa.profs],
                 Tr_rsla=[p.Tr for p in rsla.profs], Tr_arsa=[p.Tr for p in arsa.profs]))
end

json(io, x::Real) = print(io, isfinite(x) ? @sprintf("%.7g", x) : "null")
json(io, x::Symbol) = print(io, '"', x, '"')
json(io, x::AbstractVector) = (print(io, '['); for (i, v) in enumerate(x); i > 1 && print(io, ','); json(io, v); end; print(io, ']'))
json(io, x::NamedTuple) = (print(io, '{'); for (i, (k, v)) in enumerate(pairs(x)); i > 1 && print(io, ",\n"); json(io, k); print(io, ':'); json(io, v); end; print(io, '}'))
open(joinpath(SRC, "..", "data", "pulse_drift.json"), "w") do io
    json(io, (phi0=PHI0, D=KAPPA * RHO0 * W * V0 / CGS.c_light, w=W, L=L, v0=V0, t_out=TOUT, x=XS, runs=rows))
end
println("\nwrote data/pulse_drift.json")
