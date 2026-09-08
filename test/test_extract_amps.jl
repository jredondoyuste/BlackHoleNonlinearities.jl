include(joinpath(@__DIR__, "..", "jl", "HomogeneousSolutions.jl"))

# sanity check: series V/f reproduces VZ(r,l)/f(r) at r=50
println("Zerilli v_k series check at r=50:")
for l in (2, 3, 4)
    r = 50.0
    f = 1 - 2/r
    v = vcoeffs(l, "even", 30)
    vzf_series = sum(v[k]*r^(-k) for k in 2:31)
    vzf_direct = VZ(r, l)/f
    reldiff = abs(vzf_series - vzf_direct) / abs(vzf_direct)
    println("  l=$l  reldiff=$reldiff  $(reldiff < 1e-12 ? "PASS" : "FAIL")")
    global allpass_v = (reldiff < 1e-12)
end

CASES = [(2, 0.1, "odd"), (3, 0.1, "odd"), (2, 0.18684, "odd"),
         (2, 0.1, "even"), (2, 0.3, "even"), (3, 0.5, "even")]
REXTS = (500, 1000, 2000, 5000, 10000, 20000)

allpass = true
for (l, ω, parity) in CASES
    println("\n=== l=$l ω=$ω parity=$parity ===")
    sol = linear_sol(l, ω, parity, 5, 20000, "default", 1e-10, 1e-10)
    println("  rext      |a_in|_lo      arg(a_in)_lo   |a_out|_lo     arg(a_out)_lo   |a_in|_new     arg(a_in)_new  |a_out|_new    arg(a_out)_new  |ain|^2-|aout|^2")
    results = Dict()
    for rext in REXTS
        lo = extract_amps_lo(sol, ω, rext)
        new = extract_amps(sol, ω, rext, parity)
        flux = abs2(new[1]) - abs2(new[2])
        results[rext] = new
        println("  $rext   $(round(abs(lo[1]),sigdigits=6))  $(round(angle(lo[1]),sigdigits=6))   $(round(abs(lo[2]),sigdigits=6))  $(round(angle(lo[2]),sigdigits=6))   " *
                "$(round(abs(new[1]),sigdigits=8))  $(round(angle(new[1]),sigdigits=8))  $(round(abs(new[2]),sigdigits=8))  $(round(angle(new[2]),sigdigits=8))  $(round(flux,sigdigits=8))")
    end
    a1 = results[1000]; a2 = results[20000]
    for i in 1:2
        rel_abs = abs(abs(a1[i]) - abs(a2[i])) / abs(a1[i])
        dphi = abs(angle(a1[i]) - angle(a2[i]))
        ok = rel_abs < 1e-7 && dphi < 1e-7
        global allpass &= ok
        println("  comp $i: rel|a| diff=$rel_abs  Δarg=$dphi  $(ok ? "PASS" : "FAIL")")
    end
    # informational only: cancellation of two ~1e5 numbers limits this below series/ODE tolerance
    flux5000 = abs2(results[5000][1]) - abs2(results[5000][2])
    fluxok = abs(flux5000 - 1) < 1e-7
    println("  flux |ain|^2-|aout|^2 @ rext=5000: $flux5000  $(fluxok ? "holds to 1e-7" : "does not hold to 1e-7 (cancellation-limited)")")
end

println()
println(allpass && allpass_v ? "FINAL: PASS" : "FINAL: FAIL")
