"""
Batch Q-factor scanner, one Julia process per sector (see slurm/run_sector.slurm),
multithreaded over the mode combinations in <modes_file>.

  julia --project=. --threads=N jl/scan_modes.jl <sector> <modes_file> <outdir> \
        <same_freq> <w1_min> <qnm_frac> <n1> [w2_min n2] [lambda1 lambda2]

  modes_file  lines of "l1 m1 l2 m2 l m" (# comments, blank lines ok)
  same_freq   true → ω1=ω2, scan [w1_min, qnm_frac*Re(ω_QNM(l,n=0))], n1 pts
              false → independent (ω1,ω2) grid, both capped the same way, n1×n2 pts
  w2_min/n2   only needed if same_freq=false
  lambda1/2   integration-domain scale, default 4.0/5.0 (two domains → built-in
              convergence check, Q1 vs Q2 columns in the output)

Output: outdir/q_<sector>_<l1><m1>_<l2><m2>_<l><m>.dat, columns
  # w [w2]  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2
"""

include("QFactor.jl")
using ProgressMeter, Base.Threads

# Schwarzschild gravitational QNM, n=0, Re(Mω) (Berti-Cardoso-Starinets tables, M=1).
# Extend if you need l > 7.
const QNM_RE_OMEGA = Dict(2=>0.37367, 3=>0.59943, 4=>0.80918, 5=>1.01230, 6=>1.21200, 7=>1.40970)

function lambda_domain(λ::Float64, ω::Float64)
    rmin = λ + 1.0
    rmax = 10^(λ - 1) / (2 * ω)
    return rmin, rmax
end

function read_modes(path)
    modes = NTuple{6,Int}[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        f = split(s)
        push!(modes, ntuple(i -> parse(Int, f[i]), 6))
    end
    return modes
end

function scan_1d(sector, l1, l2, l, m1, m2, ω_min, ω_max, npts, λ1, λ2, outpath)
    open(outpath, "w") do io
        println(io, "# w  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2")
        for ω in range(ω_min, ω_max, length=npts)
            rmin1, rmax1 = lambda_domain(λ1, ω)
            rmin2, rmax2 = lambda_domain(λ2, ω)
            R1 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin1, rmax=rmax1)
            R2 = qfactor_full(sector, (l1, m1, ω), (l2, m2, ω), l; rmin=rmin2, rmax=rmax2)
            println(io, ω, " ", real(R1.Qout), " ", imag(R1.Qout), " ", real(R2.Qout), " ", imag(R2.Qout), " ",
                    real(R1.Qin), " ", imag(R1.Qin), " ", real(R2.Qin), " ", imag(R2.Qin))
            flush(io)
        end
    end
end

function scan_2d(sector, l1, l2, l, m1, m2, ω1_min, ω1_max, n1, ω2_min, ω2_max, n2, λ1, λ2, outpath)
    p1 = _parity(sector[1])
    open(outpath, "w") do io
        println(io, "# w1 w2  ReQout1 ImQout1 ReQout2 ImQout2 ReQin1 ImQin1 ReQin2 ImQin2")
        for ω1 in range(ω1_min, ω1_max, length=n1)
            rmin1_λ1, rmax1_λ1 = lambda_domain(λ1, ω1)
            sol1_λ1  = linear_sol(l1, ω1, p1, rmin1_λ1, rmax1_λ1, "default", 1e-10, 1e-10)
            amps1_λ1 = extract_amps(sol1_λ1, ω1, rmax1_λ1, p1)

            rmin1_λ2, rmax1_λ2 = lambda_domain(λ2, ω1)
            sol1_λ2  = linear_sol(l1, ω1, p1, rmin1_λ2, rmax1_λ2, "default", 1e-10, 1e-10)
            amps1_λ2 = extract_amps(sol1_λ2, ω1, rmax1_λ2, p1)

            for ω2 in range(ω2_min, ω2_max, length=n2)
                R1 = qfactor_row_full(sector, sol1_λ1, amps1_λ1, (l1, m1, ω1), (l2, m2, ω2), l;
                                       rmin=rmin1_λ1, rmax=rmax1_λ1)
                R2 = qfactor_row_full(sector, sol1_λ2, amps1_λ2, (l1, m1, ω1), (l2, m2, ω2), l;
                                       rmin=rmin1_λ2, rmax=rmax1_λ2)
                println(io, ω1, " ", ω2, " ", real(R1.Qout), " ", imag(R1.Qout), " ", real(R2.Qout), " ", imag(R2.Qout), " ",
                        real(R1.Qin), " ", imag(R1.Qin), " ", real(R2.Qin), " ", imag(R2.Qin))
            end
            flush(io)
        end
    end
end

sector     = ARGS[1]
sector in SECTORS || error("Unknown sector: $sector. Use one of $(SECTORS).")
modes_file = ARGS[2]
outdir     = ARGS[3]
same_freq  = lowercase(ARGS[4]) == "true"

w1_min   = parse(Float64, ARGS[5])
qnm_frac = parse(Float64, ARGS[6])
n1       = parse(Int,     ARGS[7])

if same_freq
    λ1 = length(ARGS) >= 8  ? parse(Float64, ARGS[8])  : 4.0
    λ2 = length(ARGS) >= 9  ? parse(Float64, ARGS[9])  : 5.0
else
    w2_min = parse(Float64, ARGS[8])
    n2     = parse(Int,     ARGS[9])
    λ1     = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 4.0
    λ2     = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 5.0
end

mkpath(outdir)
modes = read_modes(modes_file)
isempty(modes) && error("No mode combinations found in $modes_file")

println("sector=$sector  same_freq=$same_freq  n_modes=$(length(modes))  nthreads=$(nthreads())")
prog = Progress(length(modes); desc="$sector modes: ", output=stderr, showspeed=true)
@threads for i in eachindex(modes)
    l1, m1, l2, m2, l, m = modes[i]
    tag = "$(l1)$(m1)_$(l2)$(m2)_$(l)$(m)"
    outpath = joinpath(outdir, "q_$(sector)_$(tag).dat")
    try
        haskey(QNM_RE_OMEGA, l) || error("no QNM Re(ω) tabulated for l=$l — add it to QNM_RE_OMEGA")
        w_max = qnm_frac * QNM_RE_OMEGA[l]
        if same_freq
            scan_1d(sector, l1, l2, l, m1, m2, w1_min, w_max, n1, λ1, λ2, outpath)
        else
            scan_2d(sector, l1, l2, l, m1, m2, w1_min, w_max, n1, w2_min, w_max, n2, λ1, λ2, outpath)
        end
    catch e
        println(stderr, "FAILED sector=$sector mode=$tag: $e")
    end
    next!(prog)
end
finish!(prog)
