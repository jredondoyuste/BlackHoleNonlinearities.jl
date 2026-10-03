isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function self_coupling()
    p = plot(xlabel="M omegaout / Re(M omegaQNM)", ylabel="|Qin|",
             yscale=:log10, legend=:topleft)
    x = range(0.45, 1.08, length=8)
    for l in 2:4
        lout = 2l
        q = [quick_q("eee", (l, l, z*QNM_RE[lout]/2),
                     (l, l, z*QNM_RE[lout]/2), lout) for z in x]
        plot!(p, x, abs.(q), marker=:circle, label="l = $l")
    end
    vline!(p, [1.0], linestyle=:dot, color=:gray, label=false)
    return write_plot(p, "self_coupling.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && self_coupling()
