# The derived variables of the ARSA and its error as a function of phi0, for the supercritical radiative
# shock of sim_mf_shocks.jl (3D slab of 64 x 2 x 2 cells, 10 x the mass opacity of Melon Fuksman et al.). Uses the
# set-up and the predictions (`predict`) of sim_mf_shocks.jl. Prints a table, which run_all.sh saves in
# data/arsa_phi0.txt: a header with the numbers that do not depend on phi0, then one row per phi0,
# whitespace-separated. About 2 min.
#
#   ./run_all.sh arsa_phi0      (or: julia --project=. -t 1 src/sim_arsa_phi0.jl)
const NOMAIN = true
include("sim_mf_shocks.jl")

const SETUPS = ((2, "supercritical, kappa rho = 3.1e-9 rho/rho0 /cm (10 x mass opacity)", 20e5, 7.5e3, OP_10),)
# the sweep of phi0, in half decades; a phi0 is skipped where chat = phi0 c is below the fastest matter signal,
# the inflow speed plus the sound speed of the hottest gas of the reference
const PHI0S = (0.1, 0.03, 0.01, 0.003, 0.001)
const SIGS = (sqrt(3), 1.0)                    # sum_d sigma_d in 3D: isotropic, streaming along x

head = ["# arsa_phi0.jl: derived variables of the ARSA and errors against chat = c, 64 x $NT x $NT cells, C_r = 0.9 sqrt3 (beta = 0.9)",
        "# phi_*, its threshold and the cell counts come from the final chat = c state, at the 3D step nu0 = 1/sqrt3",
        "# (beta = 1); the suffix _a is isotropic radiation (sum sigma = sqrt3), _b streaming along x (sum sigma = 1). Columns:",
        "#   shock, phi0",
        "#   nstar     cells with phi_* < phi0",
        "#   thr_min, thr_max   smallest and largest eq:tau_one threshold of those cells (cells with tau_c above it should reach phi = 1)",
        "#   nthr      cells with phi_* < phi0 and tau_c above the threshold (predicted at phi = 1)",
        "#   phi_min, phi_max, nphi1   the ARSA phi at the end of the run: smallest, largest, cells at phi = 1",
        "#   substeps  radiation substeps of the ARSA (the RSLA takes the same)",
        "#   rsla_T%, arsa_T%, rsla_Tr%, arsa_Tr%   relative L1 errors of T_gas and T_rad, in percent",
        "#"]
cols = ("shock", "phi0", "nstar_a", "thr_min_a", "thr_max_a", "nthr_a", "nstar_b", "thr_min_b", "nthr_b",
        "phi_min", "phi_max", "nphi1", "substeps", "rsla_T%", "arsa_T%", "rsla_Tr%", "arsa_Tr%")
rows = String[]
for (k, name, u, tmax, op) in SETUPS
    ref = score(shock(; u, n, phi0=1.0, opacity=op, cfl=CR_REF), tmax; track=true)
    tau = [Honghu.planck_opacity(op, r, 10.0) * r * dx for r in ref.rho]
    h = hratio.(ref.T, ref.rho)
    hm = max(ref.hmax, maximum(h))
    push!(head, @sprintf("# shock %d: %s, u = %g km/s, t = %.3g s; h_max = %.3g (P = %.3g), tau_c in [%.3g, %.3g], largest phi_* = %.3g (isotropic), %.3g (streaming)",
                         k, name, u / 1e5, tmax, hm, hm / (4 * (GAMMA - 1)), minimum(tau), maximum(tau),
                         predict(h, tau, hm, 1.0, SIGS[1]).phistar, predict(h, tau, hm, 1.0, SIGS[2]).phistar))
    vmax = u + sqrt(GAMMA * CGS.k_B * maximum(ref.T) / (MU * CGS.m_u))
    push!(head, @sprintf("#   fastest matter signal u + c_s(T max) = %.3g km/s; chat = phi0 c is %.3g km/s at the smallest phi0",
                         vmax / 1e5, minimum(PHI0S) * CGS.c_light / 1e5))
    for phi0 in PHI0S
        if phi0 * CGS.c_light < vmax
            push!(head, @sprintf("#   shock %d: phi0 = %g skipped (chat below the matter speed)", k, phi0))
            continue
        end
        a = score(shock(; u, n, phi0, opacity=op), tmax)
        b = score(shock(; u, n, phi0, lim=StabilityBound(), opacity=op), tmax)
        qa, qb = (predict(h, tau, hm, phi0, sig) for sig in SIGS)
        push!(rows, @sprintf("%5d %8.3g %7d %9.3g %9.3g %6d %7d %9.3g %6d %8.3g %8.3g %5d %8d %8.4g %8.4g %8.4g %8.4g",
                             k, phi0, qa.nstar, qa.thr_min, qa.thr_max, qa.nthr, qb.nstar, qb.thr_min, qb.nthr,
                             minimum(b.phi), maximum(b.phi), count(b.phi .>= 0.999), b.sub,
                             100l1(a.T, ref.T), 100l1(b.T, ref.T), 100l1(a.Tr, ref.Tr), 100l1(b.Tr, ref.Tr)))
    end
end
widths = (5, 8, 7, 9, 9, 6, 7, 9, 6, 8, 8, 5, 8, 8, 8, 8, 8)
table = vcat(head, join((lpad(c, w) for (c, w) in zip(cols, widths)), " "), rows)
println(join(table, "\n"))
