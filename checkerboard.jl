# Checkerboard stability experiment in Honghu.jl: a uniform, static medium with a
# grid-scale (checkerboard) perturbation of E_rad, run at a fixed local phi.
# The measured growth factor per substep is compared with the von Neumann
# prediction for the IMEX PD-ARS scheme,
#
#     stable  <=>  phi c dt [ sum_d sigma/dx_d - rho kappa/(2 + K) ] <= 1,
#     K = 4 a_r T^3 c dt rho kappa / (rho c_v),   sigma = 1/sqrt(3) here.
#
#   julia --project=. -t 4 checkerboard.jl
using Honghu, Printf

struct FixedPhi <: Honghu.AbstractChatLimiter
    phi::Float64
end
Honghu.chat_factor(l::FixedPhi, phi0, cell::Honghu.LocalCell) = l.phi
Honghu._convert_limiter(l::FixedPhi, ::Type{T}) where {T} = l

const SIG = 1 / sqrt(3)

"The CFL-limited RSLA step and the `phi` at which the bound's left side equals `fac`."
function phi_at(fac; D, n, phi0, cfl, tau, Pcal, gamma=5 / 3)
    dx = 1 / n
    dt = cfl / (phi0 * D / dx)
    kappa = tau / dx
    K = 4 * Pcal * (gamma - 1) * dt * kappa        # 4 a T^3 c dt rho kappa/(rho cv)
    den = dt * (D * SIG / dx - kappa / (2 + K))
    return den <= 0 ? 1.0 : min(1.0, fac / den)
end

"""
Growth factor per substep of the checkerboard projection of E_rad, measured
over the last `nmeas` of `nst` substeps, for a limiter `lim` (a `FixedPhi` or
a real limiter).
"""
function growth(lim; D=1, n=16, phi0=1e-3, tau=0.5, Pcal=1e-6, cfl=0.4, gamma=5 / 3,
                nst=40, nmeas=20, amp=1e-3, rec=PCM(), hydro=false, cfl_hydro=0.4, c=1.0, beta_order=0, pert=:E)
    cv = 1.0
    a_rad = Pcal * (gamma - 1) * cv                # calP = a T^4/(rho a_iso^2), T = rho = 1
    dx = 1 / n
    kappa = tau / dx
    rad = RadiationPhysics(; c, chat=phi0 * c, a_rad, beta_order,
                           opacity=ConstantOpacity(kappa), chat_limiter=lim)
    dims = ntuple(d -> d <= D ? n : 1, 3)
    sim = Simulation(; n=dims, eos=IdealGas(; gamma, cv), radiation=rad, hydro, cfl=cfl_hydro,
                     boundaries=Boundaries(Periodic()), cfl_radiation=cfl,
                     radiation_reconstruction=rec)
    cell(x) = floor(Int, x / dx)
    set_initial_conditions!(sim) do x, y, z
        s = iseven(cell(x) + (D > 1 ? cell(y) : 0) + (D > 2 ? cell(z) : 0)) ? 1.0 : -1.0
        pert === :E ? (rho=1.0, T=1.0, Erad=a_rad * (1 + amp * s)) :
            (rho=1.0, T=1.0, Erad=a_rad, Fx=amp * s * c * a_rad)
    end
    g = sim.grid
    I = (g.lo[1]:g.hi[1], g.lo[2]:g.hi[2], g.lo[3]:g.hi[3])
    sgn = [(-1.0)^((i - g.lo[1]) + (D > 1) * (j - g.lo[2]) + (D > 2) * (k - g.lo[3]))
           for i in I[1], j in I[2], k in I[3]]
    iv = pert === :E ? 1 : 2
    proj() = abs(sum(sim.R[iv, I...] .* sgn)) / length(sgn) / a_rad
    dt = Honghu.radiation_timestep(sim.radiation, g, sim.cfl_radiation)  # hydro CFL is looser when chat >> a
    a1 = 0.0
    for s in 1:nst
        step!(sim, dt)
        s == nst - nmeas && (a1 = proj())
        isfinite(proj()) || return Inf
    end
    return (proj() / a1)^(1 / nmeas)
end

if !isdefined(Main, :NOMAIN)
    println("Table 1. Fixed phi at a fraction of the predicted boundary.")
    println("  D   tau   Pcal    phi/phi_b  phi       g/substep")
    for D in (1, 2, 3), (tau, Pcal) in ((0.3, 1e-6), (0.3, 1e-2), (3.0, 1e-2), (1.0, 1e-1))
        for fac in (0.9, 1.1)
            n = 16
            phi = phi_at(fac; D, n, phi0=1e-3, cfl=0.4, tau, Pcal)
            gg = growth(FixedPhi(phi); D, n, tau, Pcal, nst=D == 3 ? 30 : 40)
            @printf("  %d  %4.1f  %6.0e   %4.2f     %.5f   %.5f\n", D, tau, Pcal, fac, phi, gg)
        end
    end
end
