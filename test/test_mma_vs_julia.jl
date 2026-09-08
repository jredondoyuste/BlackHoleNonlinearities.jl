"""
Compares the Julia-translated source terms (jl/generated/SourceXYZ.jl) against the Mathematica
reference (tools/mma_reference.wls -> test/reference/mma_{source,coeffs}_XYZ.dat).
Usage: julia --project=. test/test_mma_vs_julia.jl (nonzero exit on any mismatch).
"""

using DelimitedFiles
using Printf

const HERE = @__DIR__
const REFDIR = joinpath(HERE, "reference")

include(joinpath(HERE, "..", "jl", "Source.jl"))

# must match tools/mma_reference.wls exactly (same rationals, same order c1m2,c1m1,c2m2,c2m1,c3m2,c4m2)
const C1M2 = 3/10 + (1/5)im
const C1M1 = -9/20 + (13/20)im
const C2M2 = 11/20 - (7/20)im
const C2M1 = -1/4 - (1/2)im
const C3M2 = 2/5 + (3/20)im
const C4M2 = -3/5 + (1/10)im

const RTOL = 1e-9
const ATOL = 1e-12

const SECTORS = ["EEE", "EEO", "EOE", "EOO", "OOE", "OOO"]

function read_dat(path)
    raw = readdlm(path, ' '; comments=true, comment_char='#')
    header_line = ""
    open(path) do io
        for line in eachline(io)
            if startswith(line, "#")
                header_line = line
                break
            end
        end
    end
    cols = split(strip(lstrip(header_line, '#')))
    colidx = Dict{String,Int}(c => i for (i, c) in enumerate(cols))
    return raw, colidx
end

# "0", "1", "m1", "m2" -> 0, 1, -1, -2
function parse_signed_idx(tok::AbstractString)
    if startswith(tok, "m")
        return -parse(Int, tok[2:end])
    else
        return parse(Int, tok)
    end
end

relerr(a, b) = abs(a - b) <= ATOL ? 0.0 : abs(a - b) / max(abs(b), eps())

function check(name, got, expected, maxerr_ref, failures)
    e = relerr(got, expected)
    maxerr_ref[] = max(maxerr_ref[], e)
    ok = abs(got - expected) <= ATOL || e <= RTOL
    if !ok
        push!(failures, @sprintf("%-40s got=%-40s expected=%-40s relerr=%.3e",
                                  name, string(got), string(expected), e))
    end
    return ok
end

function test_coeffs_file(sector, failures)
    path = joinpath(REFDIR, "mma_coeffs_$(sector).dat")
    isfile(path) || return 0.0, 0
    raw, colidx = read_dat(path)
    nrows = size(raw, 1)
    maxerr_ref = Ref(0.0)
    ncomparisons = 0

    # base names, without re_/im_ prefix, in column order
    basenames = String[]
    for c in sort(collect(keys(colidx)); by = k -> colidx[k])
        if startswith(c, "re_")
            push!(basenames, c[4:end])
        end
    end

    sector_lc = lowercase(sector)

    for row in 1:nrows
        l1 = Int(raw[row, colidx["l1"]])
        m1 = Int(raw[row, colidx["m1"]])
        l2 = Int(raw[row, colidx["l2"]])
        m2 = Int(raw[row, colidx["m2"]])
        l  = Int(raw[row, colidx["l"]])
        w1 = raw[row, colidx["re_w1"]] + im * raw[row, colidx["im_w1"]]
        w2 = raw[row, colidx["re_w2"]] + im * raw[row, colidx["im_w2"]]

        cf = coefficients(sector_lc, l1, l2, l, m1, m2, w1, w2;
                           c1m2=C1M2, c1m1=C1M1, c2m2=C2M2, c2m1=C2M1,
                           c3m2=C3M2, c4m2=C4M2)

        for base in basenames
            expected = raw[row, colidx["re_"*base]] + im * raw[row, colidx["im_"*base]]
            parts = split(base, "_")
            got = nothing
            label = "$(sector) case=($l1,$m1,$l2,$m2,$l) $base"
            if parts[1] == "c" && length(parts) == 3
                i = parse(Int, parts[2]); j = parse_signed_idx(parts[3])
                fname = Symbol("c_$(i)_$(j < 0 ? "m$(-j)" : string(j))")
                got = getfield(cf, fname)
            elseif parts[1] == "clog" && length(parts) == 2
                i = parse(Int, parts[2])
                fname = Symbol("clog_$(i)")
                got = getfield(cf, fname)
            elseif parts[1] == "zq" && length(parts) == 4
                i = parse(Int, parts[2]); j = parse(Int, parts[3]); k = parse_signed_idx(parts[4])
                zqd = cf.zq(1, -1)
                got = zqd[(i, j, k)]
            else
                continue
            end
            ncomparisons += 1
            check(label, got, expected, maxerr_ref, failures)
        end
    end
    return maxerr_ref[], ncomparisons
end

function test_source_file(sector, failures)
    path = joinpath(REFDIR, "mma_source_$(sector).dat")
    isfile(path) || return 0.0, 0, true
    raw, colidx = read_dat(path)
    nrows = size(raw, 1)
    maxerr_ref = Ref(0.0)
    ncomparisons = 0
    sector_lc = lowercase(sector)

    # grouped by (l1,m1,l2,m2,l,r) to check theta-independence of the reference itself
    theta_groups = Dict{NTuple{6,Float64},Vector{Tuple{Float64,ComplexF64}}}()

    for row in 1:nrows
        l1 = Int(raw[row, colidx["l1"]])
        m1 = Int(raw[row, colidx["m1"]])
        l2 = Int(raw[row, colidx["l2"]])
        m2 = Int(raw[row, colidx["m2"]])
        l  = Int(raw[row, colidx["l"]])
        w1 = raw[row, colidx["re_w1"]] + im * raw[row, colidx["im_w1"]]
        w2 = raw[row, colidx["re_w2"]] + im * raw[row, colidx["im_w2"]]
        theta = raw[row, colidx["theta"]]
        r = raw[row, colidx["r"]]
        phi1 = raw[row, colidx["re_phi1"]] + im * raw[row, colidx["im_phi1"]]
        psi1 = raw[row, colidx["re_psi1"]] + im * raw[row, colidx["im_psi1"]]
        phi2 = raw[row, colidx["re_phi2"]] + im * raw[row, colidx["im_phi2"]]
        psi2 = raw[row, colidx["re_psi2"]] + im * raw[row, colidx["im_psi2"]]
        expected = raw[row, colidx["re_S"]] + im * raw[row, colidx["im_S"]]

        c1m2 = raw[row, colidx["re_c1m2"]] + im * raw[row, colidx["im_c1m2"]]
        c1m1 = raw[row, colidx["re_c1m1"]] + im * raw[row, colidx["im_c1m1"]]
        c2m2 = raw[row, colidx["re_c2m2"]] + im * raw[row, colidx["im_c2m2"]]
        c2m1 = raw[row, colidx["re_c2m1"]] + im * raw[row, colidx["im_c2m1"]]
        c3m2 = raw[row, colidx["re_c3m2"]] + im * raw[row, colidx["im_c3m2"]]
        c4m2 = raw[row, colidx["re_c4m2"]] + im * raw[row, colidx["im_c4m2"]]

        cf = coefficients(sector_lc, l1, l2, l, m1, m2, w1, w2;
                           c1m2=c1m2, c1m1=c1m1, c2m2=c2m2, c2m1=c2m1,
                           c3m2=c3m2, c4m2=c4m2)
        got = source(sector_lc, cf, phi1, psi1, phi2, psi2, r)

        ncomparisons += 1
        label = "$(sector) case=($l1,$m1,$l2,$m2,$l) r=$r theta=$theta"
        check(label, got, expected, maxerr_ref, failures)

        key = (Float64(l1), Float64(m1), Float64(l2), Float64(m2), Float64(l), Float64(r))
        push!(get!(theta_groups, key, Vector{Tuple{Float64,ComplexF64}}()), (theta, expected))
    end

    theta_ok = true
    for (key, vals) in theta_groups
        length(vals) < 2 && continue
        base = vals[1][2]
        for (th, v) in vals[2:end]
            if abs(v - base) > ATOL && abs(v - base) / max(abs(base), eps()) > RTOL
                theta_ok = false
                push!(failures, @sprintf("%s theta-independence FAILED: %s vs %s",
                                          sector, string(base), string(v)))
            end
        end
    end

    return maxerr_ref[], ncomparisons, theta_ok
end

function main()
    overall_ok = true
    println(rpad("sector", 8), rpad("coeffs max relerr", 22), rpad("source max relerr", 22),
            rpad("theta-consistent", 18), "n_coeff_cmp  n_src_cmp")
    println("-"^100)

    for sector in SECTORS
        failures = String[]
        coeffs_maxerr, n_coeff = test_coeffs_file(sector, failures)
        src_maxerr, n_src, theta_ok = test_source_file(sector, failures)

        sector_ok = isempty(failures) && theta_ok
        overall_ok &= sector_ok

        println(rpad(sector, 8), rpad(@sprintf("%.3e", coeffs_maxerr), 22),
                rpad(@sprintf("%.3e", src_maxerr), 22),
                rpad(theta_ok ? "yes" : "NO", 18),
                rpad(string(n_coeff), 13), n_src)

        if !isempty(failures)
            println("  Failures for $sector:")
            for f in first(failures, 20)
                println("    ", f)
            end
            if length(failures) > 20
                println("    ... and ", length(failures) - 20, " more")
            end
        end
    end

    println("-"^100)
    if overall_ok
        println("PASS: Julia source/coefficients match the Mathematica reference (rtol=$RTOL, atol=$ATOL)")
    else
        println("FAIL: mismatches found -- see above")
    end
    return overall_ok
end

ok = main()
if !ok
    exit(1)
end
# no exit(0) here: explicit exit() after using WignerSymbols crashes on shutdown in this env
