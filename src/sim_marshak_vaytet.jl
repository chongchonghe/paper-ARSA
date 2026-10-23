# The multigroup Marshak wave (Vaytet et al. 2011; the set-up of Honghu.jl's marshak_vaytet) as a 3D slab:
# 64 x NT x NT cubic cells, periodic in y and z, the same in every column. 4 groups. The RSLA and the ARSA
# take the 3D RSLA step at beta = 0.9, C_r = 3 nu0 = 0.9 sqrt3; the reference chat = c takes nu0 = 0.4 per
# direction (C_r = 1.2). phi_* and the threshold of eq:tau_one are for isotropic radiation in 3D at beta = 1:
# nu0 = 1/sqrt3, sum_d sigma_d = sqrt3, so phi_* = h/(1 + h).
#
#   1. phi0 = 0.003, 0.01, 0.03: the L1 errors of T_gas and the sub-steps, and the phi0 at which the
#      RSLA matches the error of the ARSA, by bisection in log(phi0) (14 steps, as in irsla_paper.jl,
#      section equal_accuracy); the speedup is the ratio of sub-steps. -> scores.csv
#   2. phi0 = 0.01: profiles of T_gas, of the radiation temperature of each group and of phi of each
#      group, for the reference, the RSLA and the ARSA. -> profiles.csv
#   3. For each group: the edges, the opacity at the geometric centre of the group (as the problem
#      evaluates it) and the cell optical depth. -> groups.csv
#   4. Which groups carry the energy, in the reference: the share of each group in the flux at t_end
#      at four positions, and in the energy given to the gas over the run, c rho kappa_g (E_g -
#      f_g a_r T^4) summed over the steps, behind and ahead of the front (T_gas = 650 K). -> budget.csv
#   5. Which groups cause the error of the RSLA: at the step of the reference, chat = 0.01 c in the thin
#      group only (rho kappa dx < 1, group 4), in the thick groups only (1-3), and in all groups; the
#      L1 error of T_gas of each. -> attribution.csv
#
# The radiation temperature of a group is its brightness temperature, the T at which a blackbody puts
# the energy E_g in the group, a_r T^4 f_g(T) = E_g, with f_g the Planck fraction of the group
# (planck_fractions); every group has T_g = T_gas in equilibrium. Writes data/marshak_vaytet_*.csv; run_all.sh
# saves the printout in data/marshak_vaytet.txt.
#
#   ./run_all.sh marshak_vaytet      (or: julia --project=. -t 1 src/sim_marshak_vaytet.jl)
using Honghu, Honghu.Problems, Printf
const CGS = Honghu.CGS
const NG, N, RHO, LX, KAPPA0, NU_PIVOT = 4, 64, 1e-3, 20.0, 2000.0, 4e13
const NT = 2                                         # cells in y and in z
const CR, CR_REF = 0.9 * sqrt(3), 1.2                # C_r = 3 nu0: RSLA/ARSA (nu0 = 0.9/sqrt3), reference (nu0 = 0.4)
const SIGSUM, NU1 = sqrt(3), 1 / sqrt(3)             # isotropic radiation in 3D; largest stable nu0 (beta = 1)
const PHI0S, PHI0_PROFILE = (0.003, 0.01, 0.03), 0.01
const RHOCV = 1e-3                                   # rho c_v of marshak_vaytet, erg cm^-3 K^-1

# The set-up of Honghu.jl's marshak_vaytet (problems/radiation.jl), on an N x nt x nt grid of cubic cells,
# periodic in y and z: a wall in radiative equilibrium at 1000 K on the left, outflow on the right, a static
# slab at 300 K with rho c_v = 1e-3 erg cm^-3 K^-1 and kappa = 2000 (nu/4e13 Hz)^-2 cm^2/g.
function vay(; nt=NT, chat_over_c=1.0, chat_limiter=ConstantChat(), cfl_radiation=CR_REF)
    rho0, T_initial, T_L = RHO, 300.0, 1000.0
    c, a_rad, k_B = CGS.c_light, CGS.a_rad, CGS.k_B
    groups = RadiationGroups(6.0e10, 6.0e14, NG; energy_unit=CGS.h_planck, spacing=:log)
    opacity = FrequencyOpacity((nu, rho, T) -> KAPPA0 * (NU_PIVOT / nu)^2)
    rad = RadiationPhysics(; c, chat=c * chat_over_c, a_rad, k_B, groups, opacity,
                           closure=LevermoreClosure(), coupling=QuokkaCoupling(), chat_limiter, beta_order=0)
    fracL = planck_fractions(groups, T_L, k_B)
    wall = Fixed(ntuple(m -> (gr = (m - 1) ÷ Honghu.NRAD;
                              (m - 1) % Honghu.NRAD == 0 ? a_rad * T_L^4 * fracL[gr+1] : 0.0),
                        Honghu.NRAD * NG))
    per = (Periodic(), Periodic())
    dx = LX / N
    sim = Simulation(; n=(N, nt, nt), xhi=(LX, nt * dx, nt * dx),
                     eos=IdealGas(; gamma=5 / 3, cv=RHOCV / rho0), radiation=rad, hydro=false,
                     boundaries=Boundaries((Outflow(), Outflow()), per, per),
                     radiation_boundaries=Boundaries((wall, Outflow()), per, per),
                     floors=Floors(erad=a_rad * T_initial^4 * 1e-20), cfl_radiation)
    set_initial_conditions!(sim) do x, y, z
        (rho=rho0, T=T_initial, Erad=a_rad * T_initial^4, Fx=0.0)
    end
    return sim, 1.36e-7
end

# A fixed phi per group: `thin` where rho kappa_E dx < 1 (group 4), `thick` elsewhere (groups 1-3)
struct GroupPhi <: Honghu.AbstractChatLimiter
    thin::Float64
    thick::Float64
end
Honghu.chat_factor(l::GroupPhi, phi0, w::Honghu.LocalCell) = w.chiE * w.dxmin < 1 ? l.thin : l.thick
Honghu._convert_limiter(l::GroupPhi, ::Type{T}) where {T} = l
Honghu.per_group(::GroupPhi) = true

kappa_g(groups, gr) = KAPPA0 * (NU_PIVOT / sqrt(groups.edges[gr] * groups.edges[gr+1]))^2   # as marshak_vaytet

# T_gas, E, F and phi of every group at t_end; phi from the last full step (the last step may be cut
# short to hit t_end). With `heat = true`, also the energy each group gives the gas per unit volume
# over the run, c rho kappa_g (E_g - f_g a_r T^4) dt summed over the steps.
function marshak(; heat=false, kw...)
    sim, tmax = vay(; kw...)
    g = sim.grid
    ii, j, k = g.lo[1]:g.hi[1], g.lo[2], g.lo[3]
    groups = sim.radiation.groups
    cur = zeros(NG, N); prev = copy(cur); dts = [0.0, 0.0]; Q = zeros(NG, N)
    cb = function (s)
        prev .= cur; cur .= view(s.rw.chat, 1:NG, ii, j, k)
        dts .= (dts[2], s.dt)
        if heat
            T = gas_temperature(s)[:, 1, 1]
            for (n, i) in enumerate(ii), gr in 1:NG
                f = planck_fractions(groups, T[n], CGS.k_B)[gr]
                E = s.R[Honghu.gbase(gr)+Honghu.IERAD, i, j, k]
                Q[gr, n] += s.dt * CGS.c_light * RHO * kappa_g(groups, gr) * (E - f * CGS.a_rad * T[n]^4)
            end
        end
        true
    end
    run!(sim; tmax, verbose=false, callback=cb)
    E = [sim.R[Honghu.gbase(gr)+Honghu.IERAD, i, j, k] for gr in 1:NG, i in ii]
    F = [sim.R[Honghu.gbase(gr)+Honghu.IFRX, i, j, k] for gr in 1:NG, i in ii]
    (T=gas_temperature(sim)[:, 1, 1], E=E, F=F, Q=Q, phi=(dts[2] >= dts[1] ? cur : prev) ./ sim.radiation.c,
     sub=sim.diag.rad_substeps, groups=groups)
end
rsla(phi0) = marshak(; chat_over_c=phi0, cfl_radiation=CR)
arsa(phi0) = marshak(; chat_over_c=phi0, chat_limiter=StabilityBound(), cfl_radiation=CR)
l1(T, ref) = sum(abs, T .- ref) / sum(abs, ref)

# brightness temperature of group gr: bisection in log T of a_r T^4 f_gr(T) = E (increasing in T)
function tgroup(E, gr, groups)
    f(T) = CGS.a_rad * T^4 * planck_fractions(groups, T, CGS.k_B)[gr] - E
    lo, hi = 1.0, 1e5
    for _ in 1:100
        mid = sqrt(lo * hi)
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    sqrt(lo * hi)
end

ref = marshak(; heat=true)
dx = LX / N

# 1. errors against phi0, and the RSLA phi0 of equal accuracy
@printf("== multigroup Marshak wave, %d x %d x %d cells, %d groups, C_r = %.3g; reference chat = c: %d sub-steps ==\n",
        N, NT, NT, NG, CR, ref.sub)
scores = []
for phi0 in PHI0S
    a, b = rsla(phi0), arsa(phi0)
    ea, eb = l1(a.T, ref.T), l1(b.T, ref.T)
    err(p) = l1(rsla(p).T, ref.T)
    lo, hi = phi0, min(1.0, 30phi0)
    if err(lo) >= eb >= err(hi)
        for _ in 1:14
            mid = sqrt(lo * hi)
            err(mid) > eb ? (lo = mid) : (hi = mid)
        end
        m = rsla(hi); pm, em, sm = hi, l1(m.T, ref.T), m.sub
    else
        pm, em, sm = NaN, NaN, 0
    end
    push!(scores, (phi0, ea, eb, a.sub, b.sub, pm, em, sm, sm / b.sub))
    @printf("  phi0 = %-5g RSLA %6.3g%%  ARSA %6.3g%%  (%d, %d sub-steps); the RSLA matches the ARSA at phi0 = %.4g (%.3g%%) in %d sub-steps: speedup %.2f\n",
            phi0, 100ea, 100eb, a.sub, b.sub, pm, 100em, sm, sm / b.sub)
end

# 2. profiles at PHI0_PROFILE, and 3. the groups
a, b = rsla(PHI0_PROFILE), arsa(PHI0_PROFILE)
groups = ref.groups
nu = groups.edges
info = []
# heat-capacity ratio h = 4 a_r T^3/(rho c_v) of the reference at t_end (the gas only heats, so its largest
# value is the largest over the run), and phi_* = s h/(1 + h) for isotropic radiation, s = nu0 sum sigma = 1
hs = 4 * CGS.a_rad .* ref.T .^ 3 ./ RHOCV
@printf("\n== reference at t_end: T_gas in [%.4g, %.4g] K, h = 4 a_r T^3/(rho c_v) in [%.3g, %.3g], largest phi_* (nu0 sum sigma = %.3g) = %.3g ==\n",
        extrema(ref.T)..., extrema(hs)..., NU1 * SIGSUM, NU1 * SIGSUM * maximum(hs) / (1 + maximum(hs)))
# the threshold of eq:tau_one in each cell, for isotropic radiation in 3D; infinite where phi_* >= phi0
for phi0 in PHI0S
    ps = NU1 * SIGSUM .* hs ./ (1 .+ hs)
    thr = [ps[i] < phi0 ? 2 * (SIGSUM - phi0 / NU1) * phi0 / ((1 + hs[i]) * (phi0 - ps[i])) : Inf for i in eachindex(hs)]
    fin = thr[isfinite.(thr)]
    @printf("  phi0 = %-5g eq:tau_one threshold: finite in %d/%d cells, in [%.3g, %.3g]\n", phi0, length(fin), N,
            minimum(fin; init=Inf), maximum(fin; init=-Inf))
end
@printf("\n== profiles at phi0 = %g ==\n", PHI0_PROFILE)
for gr in 1:NG
    kappa = kappa_g(groups, gr)
    push!(info, (gr, nu[gr], nu[gr+1], kappa, kappa * RHO * dx))
    @printf("  group %d: nu = %.1e-%.1e Hz, kappa = %.3g cm^2/g, tau_c = %.3g;  ARSA phi in [%.3g, %.3g], cells at phi = 1: %d/%d\n",
            gr, nu[gr], nu[gr+1], kappa, kappa * RHO * dx, extrema(b.phi[gr, :])..., count(b.phi[gr, :] .>= 0.999), N)
end

# 4. energy budget of the reference
xs = [(i - 0.5) * dx for i in 1:N]
ifr = findlast(ref.T .> 650.0)
share(v) = 100 .* v ./ sum(v)
budget = []
@printf("\n== which groups carry the energy, in the reference; front (T_gas = 650 K) at x = %.2f cm ==\n", xs[ifr])
for (lab, i) in (("flux, first cell", 1), ("flux, x = 2 cm", 7), ("flux, front", ifr), ("flux, x = 10 cm", 32))
    push!(budget, (lab, xs[i], share(ref.F[:, i])...))
end
for (lab, r) in (("heating, behind the front", 1:ifr), ("heating, ahead of the front", ifr+1:N),
                 ("heating, whole slab", 1:N))
    push!(budget, (lab, NaN, share(vec(sum(ref.Q[:, r]; dims=2)))...))
end
for b in budget
    @printf("  %-28s %s  share of groups 1-4: %s\n", b[1], isnan(b[2]) ? "         " : @sprintf("x = %5.2f", b[2]),
            join((@sprintf("%5.1f%%", v) for v in b[3:end]), " "))
end

# 5. attribution of the error of the RSLA, at the step of the reference
attribution = []
@printf("\n== which groups cause the error of the RSLA: chat = %g c at the step of the reference ==\n", PHI0_PROFILE)
for (lab, thin, thick) in (("group 4 only", PHI0_PROFILE, 1.0), ("groups 1-3 only", 1.0, PHI0_PROFILE),
                           ("all groups", PHI0_PROFILE, PHI0_PROFILE))
    r = marshak(; chat_limiter=GroupPhi(thin, thick))
    push!(attribution, (lab, thin, thick, l1(r.T, ref.T), r.sub))
    @printf("  reduced in %-16s L1 of T_gas %6.3g%%  (%d sub-steps)\n", lab, 100l1(r.T, ref.T), r.sub)
end

writerows(name, header, rows) = open(joinpath(@__DIR__, "..", "data", "marshak_vaytet_" * name), "w") do io
    println(io, join(header, ","))
    foreach(r -> println(io, join(r, ",")), rows)
end
writerows("scores.csv", ("phi0", "L1_rsla", "L1_arsa", "sub_rsla", "sub_arsa", "phi0_rsla_matched",
                         "L1_rsla_matched", "sub_rsla_matched", "speedup"), scores)
writerows("profiles.csv", ("run", "x", "T_gas", ("T_g$gr" for gr in 1:NG)..., ("phi_g$gr" for gr in 1:NG)...),
          [(name, xs[i], r.T[i], (tgroup(r.E[gr, i], gr, groups) for gr in 1:NG)..., r.phi[:, i]...)
           for (name, r) in (("ref", ref), ("rsla", a), ("arsa", b)) for i in 1:N])
writerows("groups.csv", ("group", "nu_lo", "nu_hi", "kappa", "tau_c"), info)
writerows("budget.csv", ("quantity", "x", ("share_g$gr" for gr in 1:NG)...), budget)
writerows("attribution.csv", ("reduced_in", "phi_thin", "phi_thick", "L1", "substeps"), attribution)
println("\nwrote data/marshak_vaytet_{scores,profiles,groups,budget,attribution}.csv")
