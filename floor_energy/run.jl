# Does raising E to the floor before the exchange update (the earlier algorithm) help or hurt the ARSA,
# beyond the supercritical shock? For each problem, the ARSA with Honghu.jl's conservative rule (floor only
# where rho e + sum_g E_g c/chat_g <= 0) and with the floor everywhere E is below it. Counts the
# (cell, group, call) events with E at or below the floor after transport, and the energy the floor adds
# to the conserved budget, sum (E_floor - E) c/chat dV, as a fraction of a reference energy: the thermal
# energy of the pulse, the heat the Marshak slab gains, or the kinetic energy that flows into the shock over
# the run. Errors against the stored chat = c references of fig_diffusion_pulse/, fig_marshak_vaytet/ and
# data/mf_shocks_profiles.csv (3D slabs, as there). Writes results.csv; run.sh saves the printout in
# metadata.txt.
#
#   ./run.sh      (or: julia --project=.. -t 1 run.jl)
using Honghu, Printf
const V2 = joinpath(@__DIR__, "..")
const RESULTS = []
const MODE = Ref(:conservative)                 # :conservative or :floor
const NLOW, EADD = Ref(0), Ref(0.0)             # events with E <= floor; added energy (budget units x volume)

@eval Honghu function floor_rad_energy!(R::AbstractArray{T,4}, U::AbstractArray{T,4},
                                         chat::AbstractArray{T,4}, phys::RadiationPhysics{NG,M,T},
                                         g::CartesianGrid{T}, floors::Floors{T}) where {NG,M,T}
    ef = floors.erad / NG
    c = phys.c
    dV = prod(g.dx)
    nbad = 0
    for k in g.lo[3]:g.hi[3], j in g.lo[2]:g.hi[2], i in g.lo[1]:g.hi[1]
        rho = U[IRHO, i, j, k]
        mx, my, mz = U[IMX, i, j, k], U[IMY, i, j, k], U[IMZ, i, j, k]
        budget = U[IEGAS, i, j, k] - (mx * mx + my * my + mz * mz) / (2 * rho)
        low = false
        for gr in 1:NG
            er = R[gbase(gr)+IERAD, i, j, k]
            budget += er * c / chat[gr, i, j, k]
            if !(er > ef)
                low = true
                Main.NLOW[] += 1
            end
        end
        if Main.MODE[] === :floor || (low && !(budget > 0))
            for gr in 1:NG
                b = gbase(gr)
                er = R[b+IERAD, i, j, k]
                if !(er > ef)
                    Main.EADD[] += (ef - er) * c / chat[gr, i, j, k] * dV
                    R[b+IERAD, i, j, k] = ef
                end
            end
            low && (nbad += 1)
        end
    end
    return nbad
end

# each problem's set-up in its own module (their constants clash)
module P   # the pulse
    using Honghu, Printf
    src = read(joinpath(Main.V2, "fig_diffusion_pulse", "run.jl"), String)
    include_string(@__MODULE__, src[1:first(findfirst("\nn = 128\n", src))])
end
module Mw  # the Marshak wave
    using Honghu, Printf
    src = read(joinpath(Main.V2, "fig_marshak_vaytet", "run.jl"), String)
    include_string(@__MODULE__, src[1:first(findfirst("ref = marshak(; heat=true)", src))-1])
end
const NOMAIN = true
module S   # the shocks
    using Honghu, Printf
    include(joinpath(Main.V2, "mf_shocks.jl"))
end

readcsv(f) = (lines = readlines(f); h = split(lines[1], ","); [Dict(zip(h, split(l, ","))) for l in lines[2:end]])
num(r, k) = parse(Float64, r[k])

function both(label, run, err, scale)
    out = []
    for mode in (:conservative, :floor)
        MODE[] = mode; NLOW[] = 0; EADD[] = 0.0
        r = run()
        push!(out, (err(r), NLOW[], EADD[] / scale, r.sub))
    end
    push!(RESULTS, (label, out[1][1], out[2][1], out[1][2], out[1][4], out[2][3]))
    @printf("%-44s conservative %8.3g%%   floor %8.3g%%   events with E <= floor: %d (%d sub-steps); floor adds %.3g of the reference energy\n",
            label, 100out[1][1], 100out[2][1], out[1][2], out[1][4], out[2][3])
end

# --- the pulse: T_lo = 1e6 K, references from profiles.csv ---------------------------------------------
pr = readcsv(joinpath(V2, "fig_diffusion_pulse", "profiles.csv"))
nt2 = P.NT^2
for (lab, v0, phi0, tmax) in (("pulse, advected, phi0 = 0.01", P.V0, 0.01, P.TMAX),
                              ("pulse, advected, phi0 = v0/c", P.V0, P.V0 / P.CGS.c_light, P.TMAX),
                              ("pulse, advected, phi0 = 2 v0/c", P.V0, 2P.V0 / P.CGS.c_light, P.TMAX),
                              ("pulse, at rest, phi0 = 0.01", 0.0, 0.01, 4.8e-4))
    rows = [r for r in pr if num(r, "T_lo") == 1e6 && num(r, "v0") == v0 && isapprox(num(r, "phi0"), phi0; rtol=1e-6)]
    ref = [num(r, "T_ref") for r in rows]
    dx = P.L / 128
    rho = P.rhoprof.([num(r, "x") for r in rows], 1e6)
    epulse = sum(rho .* P.CV .* ([num(r, "T_initial") for r in rows] .- 1e6)) * dx * (P.NT * dx)^2   # pulse thermal energy, whole slab
    both(lab, () -> P.score(P.pulse(; T_lo=1e6, phi0, lim=StabilityBound(), v0); tmax),
         r -> sum(abs, r.T .- ref) / sum(abs, ref .- 1e6), epulse)
end

# --- the Marshak wave: reference T_gas from profiles.csv -------------------------------------------------
mp = [r for r in readcsv(joinpath(V2, "fig_marshak_vaytet", "profiles.csv")) if r["run"] == "ref"]
mref = [num(r, "T_gas") for r in mp]
dxm = Mw.LX / Mw.N
emar = sum(Mw.RHOCV .* (mref .- 300.0)) * dxm * (Mw.NT * dxm)^2                 # heat gained by the slab
for phi0 in (0.03, 0.01, 0.003)
    both(@sprintf("Marshak, phi0 = %g", phi0), () -> Mw.arsa(phi0), r -> sum(abs, r.T .- mref) / sum(abs, mref), emar)
end

# --- the supercritical shock: references from data/mf_shocks_profiles.csv --------------------------------
sp = readcsv(joinpath(V2, "data", "mf_shocks_profiles.csv"))
sm = Dict(r["test"] => num(r, "phi0") for r in readcsv(joinpath(V2, "data", "mf_shocks.csv")))
for (t, lab) in (("3", "shock, supercritical"), ("2", "shock, supercritical"))
    ref = [num(r, "T_ref") for r in sp if r["test"] == t]
    phi0 = sm[t]
    dxs = S.L / S.n
    eshock = 0.5 * S.RHO0 * (20e5)^2 * 20e5 * 7.5e3 * (S.NT * dxs)^2                # kinetic energy brought in over the run
    both(@sprintf("%s, phi0 = %.3g", lab, phi0),
         () -> S.score(S.shock(; u=20e5, n=S.n, phi0, lim=StabilityBound(), opacity=S.OP_10), 7.5e3),
         r -> sum(abs, r.T .- ref) / sum(abs, ref), eshock)
end

open(joinpath(@__DIR__, "results.csv"), "w") do io
    println(io, "case,L1_conservative,L1_floor,events_below_floor,substeps,floor_energy_fraction")
    foreach(r -> println(io, join((replace(r[1], "," => ";"), r[2:end]...), ",")), RESULTS)
end
println("wrote results.csv")
