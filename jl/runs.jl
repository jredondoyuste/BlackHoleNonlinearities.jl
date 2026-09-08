"""
CLI driver for nonlinear QNM computations (all sectors). Main entry: run as a script.
Usage: julia --project=<path> runs.jl <mode> <sector> [args...]

  source   sector l1 l2 l m1 m2 omega rmin rmax npts             → S(r) along a radial grid
  qfactor  sector l1 l2 l m1 m2 omega [lambda1 lambda2] [c1m1]   → Q_h at two domains + diff
  qscan    sector l1 l2 l m1 m2 omega_min omega_max npts [lambda1 lambda2]   → 1D ω scan
  qscan2d  sector l1 l2 l m1 m2 w1_min w1_max n1 w2_min w2_max n2 [lambda1 lambda2]  → 2D (ω1,ω2) scan

sector: one of SECTORS; XYZ = parity of mode1, mode2, output ("e"=even/Zerilli, "o"=odd/RW).
lambda1/lambda2: rmin=λ+1, rmax=10^(λ-1)/(2ω1); defaults 4.0/5.0.
c1m1: qfactor only, free regularization constant passed as cfree=(c1m1=c1m1,); default 0.0.

Output columns (1/2 = λ1/λ2 domain; Qout by outgoing amplitudes, r*-origin independent;
Qin by ingoing amplitudes, r*-dependent, r* = r+2log(r/2-1), M=1):
  source   → # r  Re(S)  Im(S)  Re(phi1)  Im(phi1)  Re(psi1)  Im(psi1)
  qfactor  → # ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2
  qscan    → # w  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2
  qscan2d  → # w1 w2 ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2
"""

include("QFactor.jl")
using ProgressMeter

function lambda_domain(λ::Float64, ω::Float64)
    rmin = λ + 1.0
    rmax = 10^(λ - 1) / (2 * ω)
    return rmin, rmax
end

mode   = ARGS[1]
sector = ARGS[2]
sector in SECTORS || error("Unknown sector: $sector. Use one of $(SECTORS).")

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

    sol1 = linear_sol(l1, ω, _parity(sector[1]), 5, Int(rmax), "default", 1e-10, 1e-10)
    sol2 = (l1 == l2 && sector[1] == sector[2]) ? sol1 :
           linear_sol(l2, ω, _parity(sector[2]), 5, Int(rmax), "default", 1e-10, 1e-10)
    S_fn = make_source(sector, sol1, sol2, l1, l2, l, m1, m2, ω, ω)

    println("# r  Re(S)  Im(S)  Re(phi1)  Im(phi1)  Re(psi1)  Im(psi1)")
    for r in range(rmin, rmax, length=npts)
        S    = S_fn(r)
        u, v = sol1(r)
        println(r, " ", real(S), " ", imag(S), " ",
                real(u), " ", imag(u), " ", real(v), " ", imag(v))
    end

elseif mode == "qfactor"
    l1   = parse(Int,     ARGS[3])
    l2   = parse(Int,     ARGS[4])
    l    = parse(Int,     ARGS[5])
    m1   = parse(Int,     ARGS[6])
    m2   = parse(Int,     ARGS[7])
    ω    = parse(Float64, ARGS[8])
    λ1   = length(ARGS) >= 9  ? parse(Float64, ARGS[9])  : 4.0
    λ2   = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 5.0
    c1m1 = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 0.0

    rmin1, rmax1 = lambda_domain(λ1, ω)
    rmin2, rmax2 = lambda_domain(λ2, ω)

    R1 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin1, rmax=rmax1, cfree=(c1m1=c1m1,))
    R2 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin2, rmax=rmax2, cfree=(c1m1=c1m1,))

    println("# ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2")
    println(real(R1.Qout), " ", imag(R1.Qout), " ", real(R2.Qout), " ", imag(R2.Qout), " ",
            real(R1.Qin), " ", imag(R1.Qin), " ", real(R2.Qin), " ", imag(R2.Qin))

    for (name, Q1, Q2) in (("Qout", R1.Qout, R2.Qout), ("Qin", R1.Qin, R2.Qin))
        println("$name: |1|=$(abs(Q1)) arg1/pi=$(angle(Q1)/pi)  |2|=$(abs(Q2)) arg2/pi=$(angle(Q2)/pi)  " *
                "Δ|.|=$(abs(Q2)-abs(Q1)) Δarg/pi=$(angle(Q2)/pi-angle(Q1)/pi)")
    end

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

    println("# w  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2")
    for ω in range(ω_min, ω_max, length=npts)
        local rmin1, rmax1, rmin2, rmax2, R1, R2
        rmin1, rmax1 = lambda_domain(λ1, ω)
        rmin2, rmax2 = lambda_domain(λ2, ω)
        R1 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin1, rmax=rmax1)
        R2 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin2, rmax=rmax2)
        println(ω, " ", real(R1.Qout), " ", imag(R1.Qout), " ", real(R2.Qout), " ", imag(R2.Qout), " ",
                real(R1.Qin), " ", imag(R1.Qin), " ", real(R2.Qin), " ", imag(R2.Qin))
    end

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

    p1 = _parity(sector[1])

    println("# w1 w2 ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2")
    prog = Progress(n1 * n2; desc="qscan2d: ", output=stderr, showspeed=true)
    for ω1 in range(ω1_min, ω1_max, length=n1)
        rmin1_λ1, rmax1_λ1 = lambda_domain(λ1, ω1)
        sol1_λ1  = linear_sol(l1, ω1, p1, rmin1_λ1, rmax1_λ1, "default", 1e-10, 1e-10)
        amps1_λ1 = extract_amps(sol1_λ1, ω1, rmax1_λ1, p1)

        rmin1_λ2, rmax1_λ2 = lambda_domain(λ2, ω1)
        sol1_λ2  = linear_sol(l1, ω1, p1, rmin1_λ2, rmax1_λ2, "default", 1e-10, 1e-10)
        amps1_λ2 = extract_amps(sol1_λ2, ω1, rmax1_λ2, p1)

        for ω2 in range(ω2_min, ω2_max, length=n2)
            local R1, R2
            R1 = qfactor_row_full(sector, sol1_λ1, amps1_λ1,
                             (l1, m1, ω1), (l2, m2, ω2), l;
                             rmin=rmin1_λ1, rmax=rmax1_λ1)
            R2 = qfactor_row_full(sector, sol1_λ2, amps1_λ2,
                             (l1, m1, ω1), (l2, m2, ω2), l;
                             rmin=rmin1_λ2, rmax=rmax1_λ2)
            println(ω1, " ", ω2, " ", real(R1.Qout), " ", imag(R1.Qout), " ", real(R2.Qout), " ", imag(R2.Qout), " ",
                    real(R1.Qin), " ", imag(R1.Qin), " ", real(R2.Qin), " ", imag(R2.Qin))
            next!(prog)
        end
    end
    finish!(prog)

else
    error("Unknown mode: $mode. Use source, qfactor, qscan, or qscan2d.")
end
