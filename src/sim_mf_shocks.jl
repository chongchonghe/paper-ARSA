# Radiative shocks of Melon Fuksman et al. (2021, Sec. 4.1; Ensman 1994) on a coarse grid:
# gas at 10 K with radiation in LTE, rho = 7.78e-10 g/cm^3, mu = 1, gamma = 7/5, inflow onto a
# reflecting wall in a domain of 7e10 cm, M1 closure. Subcritical: u = 6 km/s to
# t = 3.8e4 s; supercritical: u = 20 km/s to t = 7.5e3 s. Opacity: their absorption coefficient
# kappa rho = 3.1e-10 /cm, fixed (PowerLawOpacity with m = -1), or `kfac` times their mass opacity
# (ConstantOpacity, as in the paper's coarse-shock test).
#
# Every run is a 3D slab: 64 x NT x NT cubic cells, periodic in y and z, the same in every column. The
# RSLA and the ARSA take the 3D RSLA step at beta = 0.9, C_r = 3 nu0 = 0.9 sqrt3; the references chat = c
# take nu = 0.4 (C_r = 3 nu = 1.2). phi_* and the threshold of eq:tau_one are for 3D at beta = 1,
# nu0 = 1/sqrt3, with sum_d sigma_d = sqrt3 (isotropic radiation) or 1 (streaming along x).
#
# For each case: the dimensionless numbers (P, h, phi_*, tau_c, the threshold of eq:tau_one)
# from the chat = c reference, the L1 error of T_gas of the RSLA and the ARSA at phi0, and the
# errors of test 3 at other phi0, the phi0 at which the ARSA matches the error of the RSLA of test 3
# (bisection), and the ARSA with E floored before the exchange update. Writes data/mf_shocks.csv (summary),
# data/mf_shocks_profiles.csv (profiles at the final time), data/mf_shocks_phi0.csv (test 3 at phi0 =
# 0.03, 0.01, 1e-3) and data/mf_shocks_match.csv (the equal-accuracy phi0).
#
#   ./run_all.sh mf_shocks      (or: julia --project=. -t 1 src/sim_mf_shocks.jl)
using Honghu, Printf
const CGS = Honghu.CGS
const RHO0, T0, MU, GAMMA, ALPHA, L = 7.78e-10, 10.0, 1.0, 1.4, 3.1e-10, 7e10
const NT = 2                                   # cells in y and in z
const CR, CR_REF = 0.9 * sqrt(3), 1.2          # C_r = 3 nu0: RSLA/ARSA (nu0 = 0.9/sqrt3), references (nu0 = 0.4)
const NU1 = 1 / sqrt(3)                        # the largest stable nu0 in 3D (beta = 1)
const CV = CGS.k_B / (MU * CGS.m_u * (GAMMA - 1))

function shock(; u, n, phi0, lim=ConstantChat(), opacity, cfl=CR)
    rad = RadiationPhysics(; c=CGS.c_light, chat=phi0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           opacity, chat_limiter=lim)
    sim = Simulation(; n=(n, NT, NT), xlo=(0.0, 0.0, 0.0), xhi=(L, NT * L / n, NT * L / n),
                     eos=IdealGas(; gamma=GAMMA, cv=CV), radiation=rad,
                     boundaries=Boundaries((Reflecting(), Outflow()), (Periodic(), Periodic()),
                                           (Periodic(), Periodic())),
                     floors=Floors(erad=CGS.a_rad * 1e-4), cfl=0.4, cfl_radiation=cfl)
    set_initial_conditions!(sim) do x, y, z
        (rho=RHO0, vx=-u, T=T0, Erad=CGS.a_rad * T0^4, Fx=0.0)
    end
    sim
end

density(sim) = (g = sim.grid;
    [Honghu.getstate(sim.U, Val(Honghu.NHYDRO), i + g.ng[1], 1 + g.ng[2], 1 + g.ng[3])[Honghu.IRHO]
     for i in 1:g.n[1]])
hratio(T, rho) = 4 * CGS.a_rad * T^3 / (rho * CV)
l1(T, ref) = sum(abs, T .- ref) / sum(abs, ref)

function score(sim, tmax; track=false)
    hmax = Ref(0.0); at = Ref((0.0, 0.0, 0.0))
    g = sim.grid
    cur = zeros(g.n[1]); prev = copy(cur); dts = [0.0, 0.0]
    cb = function (s)
        prev .= cur; cur .= view(s.rw.chat, 1, g.lo[1]:g.hi[1], g.lo[2], g.lo[3])     # first column
        dts .= (dts[2], s.dt)
        if track
            T = gas_temperature(s)[:, 1, 1]; rho = density(s); h = hratio.(T, rho); i = argmax(h)
            h[i] > hmax[] && (hmax[] = h[i]; at[] = (s.t, T[i], rho[i]))
        end
        true
    end
    run!(sim; tmax, verbose=false, callback=cb)
    phi = (dts[2] >= dts[1] ? cur : prev) ./ sim.radiation.c
    (T=gas_temperature(sim)[:, 1, 1], Tr=radiation_temperature(sim)[:, 1, 1], rho=density(sim), phi=phi, steps=sim.diag.steps,
     sub=sim.diag.rad_substeps, hmax=hmax[], at=at[])
end

# Predictions of Sec. 3.3 from a reference state (h and tau_c of each cell) for a signal-speed sum sigsum
# at the 3D step nu0 = 1/sqrt3 (beta = 1):
# the largest phi_* (from h_max), the cells with phi_* < phi0, the range of the eq:tau_one
# threshold over those cells, and the cells whose tau_c is above it (the cells predicted at phi = 1).
function predict(h, tau, hmax, phi0, sigsum)
    ps = NU1 * sigsum .* h ./ (1 .+ h)                            # phi_* of each cell
    ok = ps .< phi0
    thr = [ok[i] ? 2 * (sigsum - phi0 / NU1) * phi0 / ((1 + h[i]) * (phi0 - ps[i])) : Inf for i in eachindex(h)]
    fin = thr[isfinite.(thr)]
    (phistar=NU1 * sigsum * hmax / (1 + hmax), nstar=count(ok), thr_min=isempty(fin) ? NaN : minimum(fin),
     thr_max=isempty(fin) ? NaN : maximum(fin), nthr=count(tau .>= thr))
end

# test, label, u, t_end, opacity, phi0 (nothing: the phi0 at which the ARSA matches the error of the RSLA
# of test 3, found by the bisection in test 3, which therefore runs before test 2)
const OP_MF = PowerLawOpacity(; kappa0=ALPHA / RHO0, rho0=RHO0, m=-1.0)   # their kappa rho, fixed
const OP_10 = ConstantOpacity(10 * ALPHA / RHO0)                          # 10 x their mass opacity
const CASES = (
    (1, "subcritical, kappa rho fixed", 6e5, 3.8e4, OP_MF, 0.001),
    (3, "supercritical, 10 x kappa", 20e5, 7.5e3, OP_10, 0.03),
    (2, "supercritical, 10 x kappa", 20e5, 7.5e3, OP_10, nothing),
)
n = 64
dx = L / n
xs = [(i - 0.5) * dx for i in 1:n]
if !isdefined(Main, :NOMAIN)                                      # NOMAIN: include the definitions only
    summary, profiles, scan, match = [], [], [], []
    phimatch = Ref(NaN)                                              # set by the bisection of test 3
    for (test, name, u, tmax, op, phi0) in CASES
        phi0 === nothing && (phi0 = phimatch[])
        t0 = time()
        ref = score(shock(; u, n, phi0=1.0, opacity=op, cfl=CR_REF), tmax; track=true)
        tau = [Honghu.planck_opacity(op, r, 10.0) * r * dx for r in ref.rho]
        h = hratio.(ref.T, ref.rho)
        hm = max(ref.hmax, maximum(h))
        @printf("\n== %s: u = %g km/s, t = %.3g s, %d x %d x %d cells, phi0 = %g (reference %.0f s) ==\n",
                name, u / 1e5, tmax, n, NT, NT, phi0, time() - t0)
        @printf("  final: T max %.4g K, rho/rho0 max %.3g; tau_c in [%.3g, %.3g]\n",
                maximum(ref.T), maximum(ref.rho) / RHO0, minimum(tau), maximum(tau))
        @printf("  h max %.3g (over the run, at t = %.3g s, T = %.4g K, rho/rho0 = %.3g); P = h/(4(gamma-1)) = %.3g\n",
                hm, ref.at[1], ref.at[2], ref.at[3] / RHO0, hm / (4 * (GAMMA - 1)))
        for (lab, sigsum) in (("isotropic", sqrt(3)), ("streaming", 1.0))
            q = predict(h, tau, hm, phi0, sigsum)
            @printf("  %-10s largest phi_* = %.3g (phi0 = %g); cells with phi_* < phi0: %d/%d; cells above the eq:tau_one threshold: %d/%d (threshold in [%.3g, %.3g])\n",
                    lab, q.phistar, phi0, q.nstar, n, q.nthr, n, q.thr_min, q.thr_max)
        end
        a = score(shock(; u, n, phi0, opacity=op), tmax)
        b = score(shock(; u, n, phi0, lim=StabilityBound(), opacity=op), tmax)
        ea, eb = l1(a.T, ref.T), l1(b.T, ref.T)
        ef = l1(score(shock(; u, n, phi0, lim=StabilityBound(flux_form=true), opacity=op), tmax).T, ref.T)
        @printf("  ARSA with phi inside the flux divergence (flux form): %.3g%%\n", 100ef)
        @printf("  RSLA %.3g%%   ARSA %.3g%%   ARSA phi: min %.3g, max %.3g, cells at phi = 1: %d/%d;  T_rad L1: RSLA %.3g%%, ARSA %.3g%%\n",
                100ea, 100eb, minimum(b.phi), maximum(b.phi), count(b.phi .>= 0.999), n,
                100 * l1(a.Tr, ref.Tr), 100 * l1(b.Tr, ref.Tr))
        push!(summary, (test, u, tmax, phi0, hm, hm / (4 * (GAMMA - 1)), minimum(tau), maximum(tau), ea, eb,
                        l1(a.Tr, ref.Tr), l1(b.Tr, ref.Tr), ef, count(b.phi .>= 0.999)))
        for i in 1:n
            push!(profiles, (test, xs[i], ref.T[i], ref.Tr[i], a.T[i], a.Tr[i], b.T[i], b.Tr[i], b.phi[i]))
        end
        if test == 3                                                  # other phi0, errors only
            for p in (0.03, 0.01, 1e-3)
                ra = score(shock(; u, n, phi0=p, opacity=op), tmax)
                rb = score(shock(; u, n, phi0=p, lim=StabilityBound(), opacity=op), tmax)
                push!(scan, (p, l1(ra.T, ref.T), l1(rb.T, ref.T), count(rb.phi .>= 0.999)))
                @printf("  test 3 at phi0 = %g: RSLA %.3g%%, ARSA %.3g%%, ARSA cells at phi = 1: %d/%d\n",
                        p, 100 * scan[end][2], 100 * scan[end][3], scan[end][4], n)
            end
            # the phi0 at which the ARSA has the error of the RSLA at this phi0: bisection in log(phi0) between
            # 1e-4 and phi0, 12 steps; a failed run (non-finite error) counts as worse
            worse(p) = !(l1(score(shock(; u, n, phi0=p, lim=StabilityBound(), opacity=op), tmax).T, ref.T) <= ea)
            lo, hi = 1e-4, phi0
            if worse(lo)
                for _ in 1:12
                    mid = sqrt(lo * hi)
                    worse(mid) ? (lo = mid) : (hi = mid)
                end
            else
                hi = lo
            end
            phimatch[] = hi
            m = score(shock(; u, n, phi0=hi, lim=StabilityBound(), opacity=op), tmax)
            push!(match, (phi0, ea, a.sub, hi, l1(m.T, ref.T), m.sub, a.sub / m.sub))
            @printf("  the ARSA matches the RSLA at phi0 = %g (%.3g%%) at phi0 = %.3g (%.3g%%): %d against %d sub-steps, speedup %.2f\n",
                    phi0, 100ea, hi, 100 * l1(m.T, ref.T), m.sub, a.sub, a.sub / m.sub)
        end
    end

    function writecsv(path, header, rows)
        open(joinpath(@__DIR__, "..", path), "w") do io
            println(io, join(header, ","))
            for r in rows
                println(io, join(r, ","))
            end
        end
        println("wrote ", path)
    end
    writecsv("data/mf_shocks.csv", ("test", "u", "t_end", "phi0", "h_max", "P_max", "tau_c_min", "tau_c_max",
                                    "L1_rsla", "L1_arsa", "L1_Trad_rsla", "L1_Trad_arsa", "L1_arsa_fluxform",
                                    "cells_phi1"), summary)
    writecsv("data/mf_shocks_profiles.csv", ("test", "x", "T_ref", "Tr_ref", "T_rsla", "Tr_rsla", "T_arsa", "Tr_arsa",
                                             "phi_arsa"), profiles)
    writecsv("data/mf_shocks_phi0.csv", ("phi0", "L1_rsla", "L1_arsa", "cells_phi1"), scan)
    writecsv("data/mf_shocks_match.csv", ("phi0_rsla", "L1_rsla", "sub_rsla", "phi0_arsa", "L1_arsa", "sub_arsa",
                                          "speedup_substeps"), match)

    # The ARSA with E raised to the floor before the exchange update, as an earlier version did, on the
    # supercritical shock with 10 x kappa. This redefines a Honghu.jl method, so it runs last, through
    # invokelatest; the reference is unfloored.
    @eval Honghu function floor_rad_energy!(R::AbstractArray{T,4}, U::AbstractArray{T,4},
                                             chat::AbstractArray{T,4}, phys::RadiationPhysics{NG,M,T},
                                             g::CartesianGrid{T}, floors::Floors{T}) where {NG,M,T}
        ef = floors.erad / NG
        n = 0
        for k in g.lo[3]:g.hi[3], j in g.lo[2]:g.hi[2], i in g.lo[1]:g.hi[1], gr in 1:NG
            ib = gbase(gr) + IERAD
            R[ib, i, j, k] > ef || (R[ib, i, j, k] = ef; n += 1)
        end
        return n
    end
    ref10 = score(shock(; u=20e5, n, phi0=1.0, opacity=OP_10, cfl=CR_REF), 7.5e3)
    for p in (0.03, phimatch[])
        f = Base.invokelatest(() -> score(shock(; u=20e5, n, phi0=p, lim=StabilityBound(), opacity=OP_10), 7.5e3))
        @printf("\n  supercritical, 10 x kappa, phi0 = %g: ARSA with E floored before the exchange update %.3g%%\n",
                p, 100 * l1(f.T, ref10.T))
    end
end
