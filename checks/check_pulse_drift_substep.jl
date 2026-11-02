# Checks that the deviation of the ARSA from the exact solution in the top row of Fig. 6 (phi0/P = 3, paper
# Sec. 5.3) comes from the length of its radiation substeps, not from the hydrodynamic step. At phi = 1 the
# radiation substeps run at a local CFL number chat dt/dx far above 1. The ARSA of sim_pulse_drift.jl (512 cells,
# phi0 = 0.01) runs to 32 widths traveled three times: as there; with radiation substeps ten times shorter and the
# same hydrodynamic step (cfl_radiation / 10, max_radiation_substeps 100); and with the hydrodynamic step about
# seven times shorter and the same substep (max_radiation_substeps 1). Prints for each the local CFL number of a
# radiation substep at phi = 1, c dt_sub/dx, the steps, and the L1 error of T - T_lo against the exact solution
# of data/pulse_drift.json (written by sim_pulse_drift.jl) at 32 widths, relative to the pulse. Reads that file;
# writes none. run_all.sh saves the printout in data/pulse_drift_substep.txt. About 3 min.
#
#   ./run_all.sh pulse_drift_substep      (or: julia --project=. checks/check_pulse_drift_substep.jl)
using Printf
const SRC = joinpath(@__DIR__, "..", "src")
src = read(joinpath(SRC, "sim_diffusion_pulse.jl"), String)
include_string(@__MODULE__, src[1:first(findfirst("\nn = 128\n", src))])    # constants, Tprof, rhoprof

const NX, PHI0, DW = 512, 0.01, 32
Ppeak(T_lo) = CGS.a_rad * Tprof(0.0, T_lo)^3 * MU / (rhoprof(0.0, T_lo) * CGS.k_B)
function find_Tlo(target)                             # as in sim_pulse_drift.jl
    lo, hi = 1e4, 1e7
    for _ in 1:100
        mid = sqrt(lo * hi)
        Ppeak(mid) < target ? (lo = mid) : (hi = mid)
    end
    sqrt(lo * hi)
end

function arsa(T_lo; cfl_rad, maxsub)                  # pulse() of sim_diffusion_pulse.jl with the substep cap exposed
    rad = RadiationPhysics(; c=CGS.c_light, chat=PHI0 * CGS.c_light, a_rad=CGS.a_rad, k_B=CGS.k_B,
                           beta_order=1, chat_limiter=StabilityBound(), opacity=ConstantOpacity(KAPPA))
    sim = Simulation(; n=(NX, NT, NT), xlo=(-L / 2, 0.0, 0.0), xhi=(L / 2, NT * L / NX, NT * L / NX),
                     eos=IdealGas(; gamma=GAMMA, cv=CV), radiation=rad, boundaries=Boundaries(Periodic()),
                     cfl=0.8, cfl_radiation=cfl_rad, max_radiation_substeps=maxsub)
    set_initial_conditions!(sim) do x, y, z
        T = Tprof(x, T_lo)
        (rho=rhoprof(x, T_lo), vx=V0, T=T, Erad=CGS.a_rad * T^4, Fx=0.0)
    end
    sim
end

# the exact solution of the first row (phi0/P = 3) at 32 widths: the json is valid Julia array syntax
json = read(joinpath(@__DIR__, "..", "data", "pulse_drift.json"), String)
row = json[first(findfirst("\"phi0_P\":3,", json)):end]
Tex = eval(Meta.parse(match(r"\"T_exact\":(\[\[.*?\]\])", row)[1]))[2]
widths = eval(Meta.parse(match(r"\"widths\":(\[.*?\])", row)[1]))
@assert widths[2] == DW && length(Tex) == NX

T_lo = find_Tlo(PHI0 / 3)
dx = L / NX
@printf("phi0/P = 3 (T_lo = %.4g K), %d cells, ARSA at phi0 = %g, to %d widths traveled; L1 of T - T_lo against the exact solution, relative to the pulse\n",
        T_lo, NX, PHI0, DW)
for (name, cfl_rad, maxsub) in (("as in Fig. 6", CR, 10), ("radiation substeps 10x shorter", CR / 10, 100),
                                ("hydro step shorter, same substeps", CR, 1))
    s = arsa(T_lo; cfl_rad, maxsub)
    run!(s; tmax=DW * W / V0, verbose=false)
    T = gas_temperature(s)[:, 1, 1]
    l1 = sum(abs.(T .- Tex)) / sum(Tex .- T_lo)
    phimin = minimum(view(s.rw.chat, 1, s.grid.lo[1]:s.grid.hi[1], s.grid.lo[2], s.grid.lo[3])) / s.radiation.c
    @printf("  %-34s c dt_sub/dx = %5.1f, %5d hydro steps, %6d radiation substeps; phi >= %.3g; L1 = %.2f%%\n",
            name, CGS.c_light * s.t / s.diag.rad_substeps / dx, s.diag.steps, s.diag.rad_substeps, phimin, 100l1)
end
