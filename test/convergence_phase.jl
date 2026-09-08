"""
Phase/modulus convergence of Q (both normalisations), per sector at 3 frequencies, against
(a) rmax, (b) quad_rmin, (c) ODE/quadrature tolerances, (d) cfree.
Usage: julia --project=. test/convergence_phase.jl [sector ...]
Writes data/convphase_{sector}_{tag}.dat, columns: kind value w |Qout| argQout |Qin| argQin.
"""
nothing
include(joinpath(@__DIR__, "..", "jl", "QFactor.jl"))
using Printf

const QNM_RE = Dict(2 => 0.37367, 3 => 0.59944, 4 => 0.80918, 5 => 1.01233)
const CHANNELS = Dict(
    "eee" => (2, 2, 2, -2, 2), "eeo" => (2, 2, 3, -2, 2), "eoe" => (2, 2, 2, -2, 3),
    "eoo" => (2, 2, 2, -2, 2), "ooe" => (2, 2, 2, -2, 2), "ooo" => (2, 2, 3, -2, 2))
DATA_DIR = joinpath(@__DIR__, "..", "data"); mkpath(DATA_DIR)
sectors = isempty(ARGS) ? collect(SECTORS) : ARGS

fmt(Q) = @sprintf("%.10e %.8f", abs(Q), angle(Q))
function row(f, kind, val, w, R)
    println(f, kind, " ", val, " ", w, " ", fmt(R.Qout), " ", fmt(R.Qin)); flush(f)
    @printf("  %-6s %-10.4g w=%.4f |Qout|=%.8e arg=%.6f |Qin|=%.8e arg=%.6f\n", kind, val, w, abs(R.Qout), angle(R.Qout), abs(R.Qin), angle(R.Qin))
end
safe(f; kw...) = try f(; kw...) catch e; e isa QuadratureNotConverged ? (println("  SKIP ", kw, " ", e); nothing) : rethrow() end

for sector in sectors
    l1, m1, l2, m2, l = CHANNELS[sector]
    ws = (0.1, 0.5 * QNM_RE[l] / 2 + 0.05, 0.75 * QNM_RE[l])   # low, intermediate, top of range (2ω = 1.5 ω_QNM)
    tag = "$(sector)_$(l1)$(m1)_$(l2)$(m2)_$(l)"
    println("\n=== $sector ($l1,$m1)x($l2,$m2)->$l ===")
    open(joinpath(DATA_DIR, "convphase_$(tag).dat"), "w") do f
        println(f, "# kind value w |Qout| argQout |Qin| argQin")
        println(f, "# sector=$sector l1=$l1 m1=$m1 l2=$l2 m2=$m2 l=$l  rstar=r+2log(r/2-1)  default: rmax=10^3/(2w) quad_rmin=5 tol=1e-10 cfree=0")
        for w in ws
            Qf(; kw...) = qfactor_full(sector, (l1, m1, w), (l2, m2, w), l; kw...)
            rdef = 10^3 / (2w); period = 2π / (2w)
            for k in (0.25, 0.5, 1, 2, 4, 8, 16)
                rmax = round(k * rdef / period) * period
                R = safe(Qf; rmax=rmax); R === nothing || row(f, "rmax", rmax, w, R)
            end
            for qr in (2, 3, 4, 5)
                R = safe(Qf; rmax=rdef, quad_rmin=qr); R === nothing || row(f, "rmin", 10.0^(-qr), w, R)
            end
            for tol in (1e-6, 1e-7, 1e-9, 1e-11)
                R = safe(Qf; rmax=rdef, atol=tol, rtol=tol, quad_rtol=tol); R === nothing || row(f, "tol", tol, w, R)
            end
            for (i, cf) in enumerate(((;), (c1m1=0.5+0.2im,), (c1m1=-1.0, c2m1=0.3im)))
                R = safe(Qf; rmax=rdef, cfree=cf); R === nothing || row(f, "cfree", i - 1, w, R)
            end
        end
    end
end
