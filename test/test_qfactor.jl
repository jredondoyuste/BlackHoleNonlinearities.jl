include(joinpath(@__DIR__, "..", "jl", "QFactor.jl"))

# cfree independence (test a) is limited by a horizon boundary term of the c-dependent source:
# at r0-2=1e-3 the spread is ~1% (ooo), ~2% (eeo), independent of ODE tolerance/rmax (see
# test/convergence_phase.jl). Proper fix: add the boundary term analytically in the regulator.
const TOL = 2.5e-2
const ω1 = 0.4
const ω2 = 0.3
const CASES = Dict("eee" => (2, 2, 2, -2, 2), "eeo" => (2, 2, 2, -2, 3), "eoe" => (2, 2, 2, -2, 3), "ooo" => (2, 2, 2, -2, 3),
                   "eoo" => (2, 2, 2, -2, 2), "ooe" => (2, 2, 2, -2, 2))
const DOM1 = (rmin=5.0, rmax=1250.0)
const DOM2 = (rmin=6.0, rmax=12500.0)

allpass = true
println("sector  test               rel_diff      time(s)  status")
for sector in SECTORS
    l1, m1, l2, m2, l = CASES[sector]
    mode1, mode2 = (l1, m1, ω1), (l2, m2, ω2)
    t0 = @elapsed Qa1 = qfactor(sector, mode1, mode2, l; DOM2...)
    t1 = @elapsed Qa2 = qfactor(sector, mode1, mode2, l; DOM2..., cfree=(c1m1=0.5+0.2im, c2m1=-0.3))
    rel_a = abs(Qa1 - Qa2) / abs(Qa1)
    global allpass &= rel_a < TOL
    println("$sector    (a) cfree indep    $(round(rel_a, sigdigits=3))      $(round(t0+t1, digits=1))  $(rel_a < TOL ? "PASS" : "FAIL")"); flush(stdout)
    t2 = @elapsed Qb1 = qfactor(sector, mode1, mode2, l; DOM1...)
    rel_b = abs(abs(Qb1) - abs(Qa1)) / abs(Qa1)
    global allpass &= rel_b < TOL
    println("$sector    (b) |Q| domain     $(round(rel_b, sigdigits=3))      $(round(t2, digits=1))  $(rel_b < TOL ? "PASS" : "FAIL")   |Q|=$(round(abs(Qa1), sigdigits=6))"); flush(stdout)
end
println(allpass ? "FINAL: PASS" : "FINAL: FAIL")
allpass || exit(1)
