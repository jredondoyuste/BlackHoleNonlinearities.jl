"""
Runner script for OOO sector computations.
Usage: julia --project=<path> run_ooo.jl <mode> [args...]

Modes:
  source   l1 l2 l m1 m2 omega rmin rmax npts [source]    → evaluate source along r
  qfactor  l1 l2 l m1 m2 omega [source]                   → compute single Q_h
  qscan    l1 l2 l m1 m2 omega_min omega_max npts [source] → frequency scan of Q_h

`source` is optional: "bruno" (default) or "adrien".
Output: whitespace-separated columns written to stdout. Q values are Q_h (h-amplitude normalized).
"""

include("QFactor.jl")  # also pulls in Source.jl, SourceOOO.jl, HomogeneousSolutions.jl

mode = ARGS[1]

if mode == "source"
    l1   = parse(Int, ARGS[2])
    l2   = parse(Int, ARGS[3])
    l    = parse(Int, ARGS[4])
    m1   = parse(Int, ARGS[5])
    m2   = parse(Int, ARGS[6])
    ω    = parse(Float64, ARGS[7])
    rmin = parse(Float64, ARGS[8])
    rmax = parse(Float64, ARGS[9])
    npts = parse(Int, ARGS[10])
    src  = length(ARGS) >= 11 ? ARGS[11] : "bruno"

    sol1 = linear_sol(l1, ω, "odd", 5, Int(rmax), "default", 1e-7, 1e-7)
    sol2 = (l1 == l2) ? sol1 : linear_sol(l2, ω, "odd", 5, Int(rmax), "default", 1e-7, 1e-7)

    S_fn = src == "adrien" ? S_OOO_adrien : S_OOO

    rs = range(rmin, rmax, length=npts)
    println("# r  Re(S)  Im(S)  Re(phi1)  Im(phi1)  Re(psi1)  Im(psi1)")
    for r in rs
        S = S_fn(sol1, sol2, ω, ω, l1, l2, l, m1, m2, r)
        u, v = sol1(r)
        println(r, " ", real(S), " ", imag(S), " ", real(u), " ", imag(u), " ", real(v), " ", imag(v))
    end

elseif mode == "qfactor"
    l1 = parse(Int, ARGS[2])
    l2 = parse(Int, ARGS[3])
    l  = parse(Int, ARGS[4])
    m1 = parse(Int, ARGS[5])
    m2 = parse(Int, ARGS[6])
    ω  = parse(Float64, ARGS[7])
    src = length(ARGS) >= 8 ? ARGS[8] : "bruno"

    local Qval = qfactor("ooo", (l1, m1, ω), (l2, m2, ω), l; source=src)
    println("# Re(Q_h)  Im(Q_h)  |Q_h|")
    println(real(Qval), " ", imag(Qval), " ", abs(Qval))

elseif mode == "qscan"
    l1      = parse(Int, ARGS[2])
    l2      = parse(Int, ARGS[3])
    l       = parse(Int, ARGS[4])
    m1      = parse(Int, ARGS[5])
    m2      = parse(Int, ARGS[6])
    ω_min   = parse(Float64, ARGS[7])
    ω_max   = parse(Float64, ARGS[8])
    npts    = parse(Int, ARGS[9])
    src     = length(ARGS) >= 10 ? ARGS[10] : "bruno"

    ωs = range(ω_min, ω_max, length=npts)
    println("# omega  Re(Q_h)  Im(Q_h)  |Q_h|")
    for ω in ωs
        if ω < 0.5
            Λ = 4.0
        else
            Λ = 5.0
        end
        local Qval = qfactor("ooo", (l1, m1, ω), (l2, m2, ω), l;
                             rmin=Λ+1, rmax=10^(Λ-1)/(2*ω), source=src)
        println(ω, " ", real(Qval), " ", imag(Qval), " ", abs(Qval))
    end

else
    error("Unknown mode: $mode. Use source, qfactor, or qscan.")
end
