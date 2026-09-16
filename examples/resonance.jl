isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function resonance()
    channels = [
        ("eee", 2, 2, 2, -2, 2),
        ("eeo", 2, 2, 4, -2, 3),
        ("eee", 2, 2, 2, 2, 4),
    ]
    p1 = plot(yscale=:log10, ylabel="|Qin|", legend=:topleft)
    p2 = plot(xlabel="M omegaout", ylabel="arg(Qin) / pi", legend=false)
    for (sector, l1, m1, l2, m2, l) in channels
        ωout = range(0.45QNM_RE[l], 1.35QNM_RE[l], length=10)
        q = [quick_q(sector, (l1, m1, w/2), (l2, m2, w/2), l) for w in ωout]
        label = "lout = $l"
        plot!(p1, ωout, abs.(q), marker=:circle, label=label)
        plot!(p2, ωout, angle.(q) ./ π, marker=:circle)
        vline!(p1, [QNM_RE[l]], linestyle=:dot, color=:gray, label=false)
        vline!(p2, [QNM_RE[l]], linestyle=:dot, color=:gray, label=false)
    end
    return write_plot(plot(p1, p2, layout=(2, 1), size=(500, 620)), "resonance.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && resonance()
