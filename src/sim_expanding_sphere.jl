# The interior of the opaque expanding sphere of Roth & Kasen (2015, ApJS 217, 9, Sec. 6.3), with absorption:
# a uniform medium in homologous expansion, v = x/(t + t0), t0 = R0/VOUT, which is an exact solution of the Euler
# equations (rho ~ (t + t0)^-3, uniform pressure, no acceleration). As in the paper, rho = 4.75e-7 g/cm^3,
# T = 1e4 K for gas and radiation, mu = 0.5, gamma = 5/3, kappa = 0.4 cm^2/g, and the sphere grows from R0 = 1e13 cm
# to 11 R0; unlike the paper, the opacity is absorption (kappa_P = kappa_E = kappa_F), so the gas and the radiation
# share one temperature on the adiabat of gas plus radiation, and the expansion is slower, VOUT = 10 km/s (the paper:
# 1e4 km/s), so that chat = 30 VOUT at phi0 = 1e-3 and the flow stays subsonic within the box (the RSLA error depends
# on h/phi0 only). The pure-scattering case of the paper is not run: with kappa_E = 0 the limiter sets phi = phi0.
#
# The box [-R0, R0]^3 of N^3 cubic cells is the interior of the sphere. Its ghost cells copy the nearest interior
# cell and take the homologous velocity at their own position, with F shifted by (4/3) E dv (HomologousBC), so that
# each run sees only its own solution. Levermore closure, v/c terms on. The RSLA and the ARSA take the 3D RSLA step
# at beta = 0.9 (C_r = 0.9 sqrt3). The reference is exact: s = (3/2) ln T - ln rho + (4/3) a T^3 mu/(rho k_B),
# the entropy of gas plus radiation in units of k_B/mu, is constant.
#
# Prints the initial h, tau_c, D = rho kappa R0 VOUT/c and phi_* = h/(1 + h), and for each run at r_out = 11 R0 the
# box mean of T_gas, its error against the exact solution, its range in the box, phi and the steps. Writes
# data/expanding_sphere.csv: for each phi0 and method (rsla, arsa) at NOUT radii r_out/R0 from 1 to 11 the time
# [s], the mean, minimum and maximum of T_gas [K] in the box, the mean T_rad [K], the range of phi, and the exact
# T [K]. run_all.sh saves the printout in data/expanding_sphere.txt. The six runs use one worker process each,
# with THREADS threads (environment variable, default 1): Honghu.jl spreads its grid loops over the threads of a
# process. On a machine with 64 cores, THREADS=10.
#
#   ./run_all.sh expanding_sphere      (or: julia --project=. src/sim_expanding_sphere.jl)
#   THREADS=10 ./run_all.sh expanding_sphere
using Distributed, Printf
const PHI0S = (1e-3, 3e-3, 1e-2)
const THREADS = parse(Int, get(ENV, "THREADS", "1"))
addprocs(2 * length(PHI0S); exeflags=`--project=$(Base.active_project()) -t $THREADS`)

@everywhere begin
using Honghu
const CGS = Honghu.CGS
const MU, GAMMA, R0, RHO0, KAP, T0, VOUT = 0.5 * CGS.m_u, 5 / 3, 1e13, 4.75e-7, 0.4, 1e4, 1e6
const N, NOUT = 32, 13
const CV = CGS.k_B / (MU * (GAMMA - 1))
const TZERO = R0 / VOUT
const ROUT = exp.(range(0, log(11), length=NOUT))     # r_out/R0 at the outputs

struct HomologousBC <: Honghu.BoundaryCondition end
function Honghu.fill_face!(A::AbstractArray{T,4}, ::HomologousBC, g::Honghu.CartesianGrid, dir::Int, side::Int,
                           t) where {T}
    rad = size(A, 1) == Honghu.NRAD
    others = dir == 1 ? (2, 3) : dir == 2 ? (1, 3) : (1, 2)
    src = side == 1 ? g.lo[dir] : g.hi[dir]
    for m in 1:g.ng[dir], kk in 1:g.dims[others[2]], jj in 1:g.dims[others[1]]
        idx = side == 1 ? g.lo[dir] - m : g.hi[dir] + m
        sub = ntuple(d -> d == dir ? idx : (d == others[1] ? jj : kk), 3)
        ins = ntuple(d -> d == dir ? src : (d == others[1] ? jj : kk), 3)
        xg, xi = Honghu.cellcenter(g, sub...), Honghu.cellcenter(g, ins...)
        for n in axes(A, 1)
            A[n, sub...] = A[n, ins...]
        end
        if rad                                        # F = (4/3) v E + the comoving flux of the interior cell
            for d in 1:3
                A[d+1, sub...] += 4 / 3 * A[1, ins...] * (xg[d] - xi[d]) / (t + TZERO)
            end
        else                                          # same rho and internal energy, homologous velocity
            rho = A[1, ins...]
            ke = sum(A[d+1, ins...]^2 for d in 1:3) / (2rho)
            for d in 1:3
                A[d+1, sub...] = rho * xg[d] / (t + TZERO)
            end
            A[5, sub...] += sum(A[d+1, sub...]^2 for d in 1:3) / (2rho) - ke
        end
    end
    return nothing
end

function sphere(; phi0, lim)
    rad = RadiationPhysics(; c=CGS.c_light, chat=phi0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           beta_order=1, chat_limiter=lim, opacity=ConstantOpacity(KAP))
    sim = Simulation(; n=(N, N, N), xlo=(-R0, -R0, -R0), xhi=(R0, R0, R0), eos=IdealGas(; gamma=GAMMA, cv=CV),
                     radiation=rad, boundaries=Boundaries(HomologousBC()), cfl=0.8, cfl_radiation=0.9 * sqrt(3))
    E0 = CGS.a_rad * T0^4
    set_initial_conditions!(sim) do x, y, z
        (rho=RHO0, vx=x / TZERO, vy=y / TZERO, vz=z / TZERO, T=T0, Erad=E0,
         Fx=4 / 3 * x / TZERO * E0, Fy=4 / 3 * y / TZERO * E0, Fz=4 / 3 * z / TZERO * E0)
    end
    sim
end

function run_case(phi0, method)
    sim = sphere(; phi0, lim=method == "arsa" ? StabilityBound() : ConstantChat())
    g = sim.grid
    inner = A -> view(A, g.lo[1]:g.hi[1], g.lo[2]:g.hi[2], g.lo[3]:g.hi[3])
    prange = [Inf, -Inf]
    cb = function (s)
        lo, hi = extrema(inner(view(s.rw.chat, 1, :, :, :))) ./ CGS.c_light
        prange .= (min(prange[1], lo), max(prange[2], hi))
        true
    end
    rows = []
    for r in ROUT
        t = (r - 1) * TZERO
        t > 0 && run!(sim; tmax=t, verbose=false, callback=cb)
        Tg, Tr = gas_temperature(sim), radiation_temperature(sim)
        push!(rows, (; r, t, Tg=sum(Tg) / length(Tg), Tgmin=minimum(Tg), Tgmax=maximum(Tg), Tr=sum(Tr) / length(Tr),
                     phimin=t > 0 ? prange[1] : phi0, phimax=t > 0 ? prange[2] : phi0))
    end
    (; phi0, method, rows, steps=sim.diag.steps, sub=sim.diag.rad_substeps)
end
end

entropy(T, rho) = 1.5 * log(T) - log(rho) + 4 / 3 * CGS.a_rad * T^3 * MU / (rho * CGS.k_B)
function exact(r)                                     # T on the adiabat at rho = RHO0 r^-3, by bisection in log T
    rho, s = RHO0 / r^3, entropy(T0, RHO0)
    lo, hi = log(1.0), log(T0)
    for _ in 1:200
        m = (lo + hi) / 2
        entropy(exp(m), rho) > s ? (hi = m) : (lo = m)
    end
    exp((lo + hi) / 2)
end

h0 = 4 * CGS.a_rad * T0^3 / (RHO0 * CV)
@printf("initial state: h = %.3g, phi_* = h/(1 + h) = %.3g, tau_c = %.3g, D = rho kappa R0 VOUT/c = %.3g, Mach at the box edge %.2f\n",
        h0, h0 / (1 + h0), RHO0 * KAP * 2R0 / N, RHO0 * KAP * R0 * VOUT / CGS.c_light,
        VOUT / sqrt(GAMMA * CGS.k_B * T0 / MU))
@printf("exact T at r_out = 11 R0: %.4g K (gas alone: %.4g K)\n", exact(11.0), T0 / 121)
@printf("%d worker processes, %d threads each\n", nworkers(), THREADS)
t1 = time()
res = pmap(a -> run_case(a...), [(p, m) for p in PHI0S for m in ("rsla", "arsa")])
for r in res
    f = r.rows[end]
    @printf("phi0 = %-6g %s: T_gas %.4g K, error %+.4f, range in the box [%.4g, %.4g] K, T_rad %.4g K, phi in [%.3g, %.3g], %d steps, %d sub-steps\n",
            r.phi0, uppercase(r.method), f.Tg, f.Tg / exact(f.r) - 1, f.Tgmin, f.Tgmax, f.Tr, f.phimin, f.phimax, r.steps, r.sub)
end
@printf("(%.0f s)\n", time() - t1)

out = joinpath(@__DIR__, "..", "data", "expanding_sphere.csv")
open(out, "w") do io
    println(io, "phi0,method,r_out_R0,t_s,Tg_mean_K,Tg_min_K,Tg_max_K,Tr_mean_K,phi_min,phi_max,T_exact_K")
    for r in res, w in r.rows
        @printf(io, "%g,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", r.phi0, r.method, w.r, w.t, w.Tg, w.Tgmin,
                w.Tgmax, w.Tr, w.phimin, w.phimax, exact(w.r))
    end
end
println("wrote data/expanding_sphere.csv")
