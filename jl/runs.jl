"""
General runner script for nonlinear QNM computations (all parity sectors).
Usage: julia --project=<path> runs.jl <mode> <parity> [args...]

Modes:
  source   parity l1 l2 l m1 m2 omega rmin rmax npts [source]
      → evaluate source term S(r) along a radial grid

  qfactor  parity l1 l2 l m1 m2 omega [lambda1 lambda2] [source]
      → compute Q_h at two integration domains, print both + difference

  qscan    parity l1 l2 l m1 m2 omega_min omega_max npts [lambda1 lambda2] [source]
      → 1D frequency scan: |Q_h| at two lambda values per ω; err = ||Q1| - |Q2||

  qscan2d  parity l1 l2 l m1 m2 w1_min w1_max n1 w2_min w2_max n2 [lambda1 lambda2] [source]
      → 2D scan over independent (ω1, ω2): Q_h at two lambda values, columns include ΔQ
      sol1(l1,ω1) is cached across the inner ω2 loop (n1 × n2 − n1 fewer ODE solves)

Arguments:
  parity          "ooo" or "ooe" (error if not implemented)
  lambda1/lambda2 integration-domain scale: rmin = λ+1, rmax = 10^(λ-1)/(2ω1)
                  defaults: lambda1=4.0  lambda2=5.0
  source          "bruno" (default) or "adrien" (OOO only)

Output columns:
  source   → # r  Re(S)  Im(S)  Re(phi1)  Im(phi1)  Re(psi1)  Im(psi1)
  qfactor  → # Re(Q1)  Im(Q1)  |Q1|  Re(Q2)  Im(Q2)  |Q2|  Re(ΔQ)  Im(ΔQ)  |ΔQ|
  qscan    → # omega  |Q1|  |Q2|  err    (err = ||Q1|-|Q2||)
  qscan2d  → # omega1  omega2  |Q1|  |Q2|  err    (err = ||Q1|-|Q2||)
"""

include("QFactor.jl")  # pulls in Source.jl, SourceOO*.jl, HomogeneousSolutions.jl
using ProgressMeter

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
    lambda_domain(λ, ω) → (rmin, rmax)

Convert a domain-scale parameter λ and frequency ω to integration bounds.
"""
function lambda_domain(λ::Float64, ω::Float64)
    rmin = λ + 1.0
    rmax = 10^(λ - 1) / (2 * ω)
    return rmin, rmax
end

"""
    source_fn(parity, src) → function

Return the source function for the given parity and source variant.
"""
function source_fn(parity::String, src::String)
    if parity == "ooo"
        return src == "adrien" ? S_OOO_adrien : S_OOO
    elseif parity == "ooe"
        return S_OOE
    else
        error("Source for $parity not yet implemented")
    end
end

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

mode   = ARGS[1]
parity = ARGS[2]

# ---------------------------------------------------------------------------
# source mode
# ---------------------------------------------------------------------------

if mode == "source"
    l1   = parse(Int,     ARGS[3])
    l2   = parse(Int,     ARGS[4])
    l    = parse(Int,     ARGS[5])
    m1   = parse(Int,     ARGS[6])
    m2   = parse(Int,     ARGS[7])
    ω    = parse(Float64, ARGS[8])
    rmin = parse(Float64, ARGS[9])
    rmax = parse(Float64, ARGS[10])
    npts = parse(Int,     ARGS[11])
    src  = length(ARGS) >= 12 ? ARGS[12] : "bruno"

    S_fn = source_fn(parity, src)
    sol1 = linear_sol(l1, ω, "odd", 5, Int(rmax), "default", 1e-7, 1e-7)
    sol2 = (l1 == l2) ? sol1 : linear_sol(l2, ω, "odd", 5, Int(rmax), "default", 1e-7, 1e-7)

    println("# r  Re(S)  Im(S)  Re(phi1)  Im(phi1)  Re(psi1)  Im(psi1)")
    for r in range(rmin, rmax, length=npts)
        S    = S_fn(sol1, sol2, ω, ω, l1, l2, l, m1, m2, r)
        u, v = sol1(r)
        println(r, " ", real(S), " ", imag(S), " ",
                real(u), " ", imag(u), " ", real(v), " ", imag(v))
    end

# ---------------------------------------------------------------------------
# qfactor mode
# ---------------------------------------------------------------------------

elseif mode == "qfactor"
    l1  = parse(Int,     ARGS[3])
    l2  = parse(Int,     ARGS[4])
    l   = parse(Int,     ARGS[5])
    m1  = parse(Int,     ARGS[6])
    m2  = parse(Int,     ARGS[7])
    ω   = parse(Float64, ARGS[8])
    λ1  = length(ARGS) >= 9  ? parse(Float64, ARGS[9])  : 4.0
    λ2  = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 5.0
    src = length(ARGS) >= 11 ? ARGS[11] : "bruno"

    rmin1, rmax1 = lambda_domain(λ1, ω)
    rmin2, rmax2 = lambda_domain(λ2, ω)

    Q1 = qfactor(parity, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin1, rmax=rmax1, source=src)
    Q2 = qfactor(parity, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin2, rmax=rmax2, source=src)
    ΔQ = Q2 - Q1

    println("# Re(Q1)  Im(Q1)  |Q1|  Re(Q2)  Im(Q2)  |Q2|  Re(ΔQ)  Im(ΔQ)  |ΔQ|")
    println(real(Q1), " ", imag(Q1), " ", abs(Q1), " ",
            real(Q2), " ", imag(Q2), " ", abs(Q2), " ",
            real(ΔQ), " ", imag(ΔQ), " ", abs(ΔQ))

# ---------------------------------------------------------------------------
# qscan mode
# ---------------------------------------------------------------------------

elseif mode == "qscan"
    l1    = parse(Int,     ARGS[3])
    l2    = parse(Int,     ARGS[4])
    l     = parse(Int,     ARGS[5])
    m1    = parse(Int,     ARGS[6])
    m2    = parse(Int,     ARGS[7])
    ω_min = parse(Float64, ARGS[8])
    ω_max = parse(Float64, ARGS[9])
    npts  = parse(Int,     ARGS[10])
    λ1    = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 4.0
    λ2    = length(ARGS) >= 12 ? parse(Float64, ARGS[12]) : 5.0
    src   = length(ARGS) >= 13 ? ARGS[13] : "bruno"

    println("# omega  |Q1|  |Q2|  err")
    for ω in range(ω_min, ω_max, length=npts)
        local rmin1, rmax1, rmin2, rmax2, Q1, Q2, err
        rmin1, rmax1 = lambda_domain(λ1, ω)
        rmin2, rmax2 = lambda_domain(λ2, ω)
        Q1 = qfactor(parity, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin1, rmax=rmax1, source=src)
        Q2 = qfactor(parity, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin2, rmax=rmax2, source=src)
        err = abs(abs(Q1) - abs(Q2))
        println(ω, " ", abs(Q1), " ", abs(Q2), " ", err)
    end

# ---------------------------------------------------------------------------
# qscan2d mode  — independent ω1, ω2 grid
# ---------------------------------------------------------------------------

elseif mode == "qscan2d"
    l1     = parse(Int,     ARGS[3])
    l2     = parse(Int,     ARGS[4])
    l      = parse(Int,     ARGS[5])
    m1     = parse(Int,     ARGS[6])
    m2     = parse(Int,     ARGS[7])
    ω1_min = parse(Float64, ARGS[8])
    ω1_max = parse(Float64, ARGS[9])
    n1     = parse(Int,     ARGS[10])
    ω2_min = parse(Float64, ARGS[11])
    ω2_max = parse(Float64, ARGS[12])
    n2     = parse(Int,     ARGS[13])
    λ1     = length(ARGS) >= 14 ? parse(Float64, ARGS[14]) : 4.0
    λ2     = length(ARGS) >= 15 ? parse(Float64, ARGS[15]) : 5.0
    src    = length(ARGS) >= 16 ? ARGS[16] : "bruno"

    println("# omega1  omega2  |Q1|  |Q2|  err")
    prog = Progress(n1 * n2; desc="qscan2d: ", output=stderr, showspeed=true)
    for ω1 in range(ω1_min, ω1_max, length=n1)
        # ── Compute and cache sol1 at both lambda domains ──────────────────
        rmin1_λ1, rmax1_λ1 = lambda_domain(λ1, ω1)
        sol1_λ1  = linear_sol(l1, ω1, "odd", rmin1_λ1, rmax1_λ1, "default", 1e-7, 1e-7)
        _, aout1_λ1 = extract_amps(sol1_λ1, ω1, rmax1_λ1)

        rmin1_λ2, rmax1_λ2 = lambda_domain(λ2, ω1)
        sol1_λ2  = linear_sol(l1, ω1, "odd", rmin1_λ2, rmax1_λ2, "default", 1e-7, 1e-7)
        _, aout1_λ2 = extract_amps(sol1_λ2, ω1, rmax1_λ2)

        # ── Inner loop over ω2 ────────────────────────────────────────────
        for ω2 in range(ω2_min, ω2_max, length=n2)
            local Q1, Q2, err
            Q1 = qfactor_row(parity, sol1_λ1, aout1_λ1,
                             (l1, m1, ω1), (l2, m2, ω2), l;
                             rmin=rmin1_λ1, rmax=rmax1_λ1, source=src)
            Q2 = qfactor_row(parity, sol1_λ2, aout1_λ2,
                             (l1, m1, ω1), (l2, m2, ω2), l;
                             rmin=rmin1_λ2, rmax=rmax1_λ2, source=src)
            err = abs(abs(Q1) - abs(Q2))
            println(ω1, " ", ω2, " ", abs(Q1), " ", abs(Q2), " ", err)
            next!(prog)
        end
    end
    finish!(prog)

else
    error("Unknown mode: $mode. Use source, qfactor, qscan, or qscan2d.")
end
