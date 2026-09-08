"""
Convergence test: does |Q_h| depend on the inner/outer integration radius?
Main entry: run as a script (julia --project=. test/convergence_test.jl).

For each sector, sweeps (a) inner cutoff r_min at fixed r_max and (b) outer cutoff r_max at
fixed r_min, near a channel's near-resonance frequency pair, against the package's default
settings as reference (quad_rmin=3, rmax at 2x the λ=4 default) — not the sweep's most extreme
point, since r_min-2M below ~1e-5 hits near-horizon cancellation noise and very large rmax hits
an oscillatory-quadrature convergence limit. Outer-sweep r_max points are snapped to an integer
multiple of the output mode's oscillation period 2π/ω_out so all points sample the same phase.

Writes, per sector:
    data/convergence_{sector}_{l1}{m1}_{l2}{m2}_{l}_inner.dat  # (rmin-2M)/M  rel_err
    data/convergence_{sector}_{l1}{m1}_{l2}{m2}_{l}_outer.dat  # rmax/M       rel_err
"""

include(joinpath(@__DIR__, "..", "jl", "QFactor.jl"))

# QNM fundamental frequencies (real part, M=1), same table as plots/compute_q_allsectors.py
const QNM_RE = Dict(2 => 0.37367, 3 => 0.59944, 4 => 0.80918, 5 => 1.01233)

# one representative non-vanishing channel per sector (satisfies the l+l1+l2 selection rule,
# doesn't vanish by l1=l2 exchange symmetry) — picked as in plots/compute_q_allsectors.py
const CHANNELS = Dict(
    "eee" => (2, 2, 2, -2, 2),
    "eeo" => (2, 2, 3, -2, 2),
    "eoe" => (2, 2, 2, -2, 3),
    "eoo" => (2, 2, 2, -2, 2),
    "ooe" => (2, 2, 2, -2, 2),
    "ooo" => (2, 2, 3, -2, 2),
)

DATA_DIR = joinpath(@__DIR__, "..", "data")
mkpath(DATA_DIR)

function logspaced(lo, hi, n)
    return exp10.(range(log10(lo), log10(hi), length=n))
end

for sector in SECTORS
    l1, m1, l2, m2, l = CHANNELS[sector]
    ω_qnm = QNM_RE[l]
    ω1 = ω_qnm / 2
    ω2 = ω_qnm / 2
    ω_out = ω1 + ω2
    tag = "$(sector)_$(l1)$(m1)_$(l2)$(m2)_$(l)"

    println("\n=== $sector  channel=($l1,$m1)×($l2,$m2)→($l)  ω1=$ω1 ω2=$ω2 ω_out=$ω_out (QNM=$ω_qnm) ===")

    λ_default = 4.0
    rmax_default = 10^(λ_default - 1) / (2 * ω_out)
    quad_rmin_default = 3   # r_min = 2.001, package default

    function compute_Q(; quad_rmin, rmax)
        qfactor(sector, (l1, m1, ω1), (l2, m2, ω2), l; rmax=rmax, quad_rmin=quad_rmin)
    end

    # catches QuadratureNotConverged so one bad point doesn't block the whole sweep
    function safe_Q(; quad_rmin, rmax)
        try
            return compute_Q(quad_rmin=quad_rmin, rmax=rmax)
        catch e
            e isa QuadratureNotConverged || rethrow()
            println("  [SKIPPED] quad_rmin=$quad_rmin rmax=$rmax : $(sprint(showerror, e))")
            flush(stdout)
            return nothing
        end
    end

    # (a) inner-cutoff sweep, rmax fixed at the default; 8 points log-spaced over (rmin-2M) in
    # [1e-5, 1e-1] (validated free of near-horizon cancellation noise); reference quad_rmin=3
    # (rmin-2M=1e-3) sits inside that range.
    inner_plotted = logspaced(1e-5, 1e-1, 8)

    println("-- inner-cutoff sweep (rmax=$(round(rmax_default,digits=1)) fixed) --"); flush(stdout)
    Q_ref_in = safe_Q(quad_rmin=quad_rmin_default, rmax=rmax_default)
    println("  [reference] rmin-2M=$(10.0^(-quad_rmin_default))  |Q|=$(Q_ref_in === nothing ? "N/A" : abs(Q_ref_in))"); flush(stdout)

    open(joinpath(DATA_DIR, "convergence_$(tag)_inner.dat"), "w") do f
        println(f, "# (rmin-2M)/M  rel_err")
        println(f, "# sector=$sector channel=($l1,$m1)x($l2,$m2)->($l) w1=$ω1 w2=$ω2 rmax=$rmax_default reference_rmin-2M=$(10.0^(-quad_rmin_default))")
        if Q_ref_in !== nothing
            for off in inner_plotted
                Q = safe_Q(quad_rmin=-log10(off), rmax=rmax_default)
                Q === nothing && continue
                rel_err = abs(abs(Q) - abs(Q_ref_in)) / abs(Q_ref_in)
                println(f, off, "  ", rel_err)
                println("  rmin-2M=$(off)  |Q|=$(abs(Q))  rel_err=$rel_err"); flush(stdout)
            end
        end
    end

    # (b) outer-cutoff sweep, quad_rmin fixed; 8 points log-spaced over rmax in
    # [rmax_default/8, rmax_default*1.5], reference at rmax_default*2 (validated to converge
    # cleanly here; rmax_default*8 does not converge for eee/eeo within quad_maxevals). All
    # points snapped to an integer multiple of the oscillation period to sample the same phase.
    period = 2π / ω_out
    snap(x) = round(x / period) * period

    outer_plotted = snap.(logspaced(rmax_default / 8, rmax_default * 1.5, 8))
    outer_ref = snap(rmax_default * 2)

    println("-- outer-cutoff sweep (quad_rmin=$quad_rmin_default fixed, period=$(round(period,digits=2))) --"); flush(stdout)
    Q_ref_out = safe_Q(quad_rmin=quad_rmin_default, rmax=outer_ref)
    println("  [reference] rmax=$(round(outer_ref,digits=1))  |Q|=$(Q_ref_out === nothing ? "N/A" : abs(Q_ref_out))"); flush(stdout)

    open(joinpath(DATA_DIR, "convergence_$(tag)_outer.dat"), "w") do f
        println(f, "# rmax/M  rel_err")
        println(f, "# sector=$sector channel=($l1,$m1)x($l2,$m2)->($l) w1=$ω1 w2=$ω2 quad_rmin=$quad_rmin_default reference_rmax=$outer_ref")
        if Q_ref_out !== nothing
            for rmx in outer_plotted
                Q = safe_Q(quad_rmin=quad_rmin_default, rmax=rmx)
                Q === nothing && continue
                rel_err = abs(abs(Q) - abs(Q_ref_out)) / abs(Q_ref_out)
                println(f, rmx, "  ", rel_err)
                println("  rmax=$(round(rmx,digits=1))  |Q|=$(abs(Q))  rel_err=$rel_err"); flush(stdout)
            end
        end
    end
end

println("\nWrote data/convergence_{sector}_{...}_inner.dat and _outer.dat for all sectors.")
