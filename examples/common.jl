using Plots

const BHN_EXAMPLE_COMMON = true
const EXAMPLE_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = joinpath(@__DIR__, "output")
const QNM_RE = Dict(2 => 0.37367, 3 => 0.59944, 4 => 0.80918,
                    5 => 1.01230, 6 => 1.21201, 7 => 1.40974, 8 => 1.60619)

isdefined(Main, :qfactor_full) || include(joinpath(EXAMPLE_ROOT, "jl", "QFactor.jl"))
mkpath(OUTPUT_DIR)
default(fontfamily="Computer Modern", linewidth=1.8, framestyle=:box,
        grid=false, legendfontsize=7, guidefontsize=10, tickfontsize=8)

function quick_q(sector, mode1, mode2, l; λ=3.0, kw...)
    rmax = 10.0^(λ - 1) / (2min(mode1[3], mode2[3]))
    return qfactor_full(sector, mode1, mode2, l;
                        rmin=4.0, rmax=rmax, atol=1e-8, rtol=1e-8,
                        quad_rtol=1e-7, kw...).Qin
end

function write_plot(p, name)
    path = joinpath(OUTPUT_DIR, name)
    savefig(p, path)
    println("wrote $path")
    return path
end
