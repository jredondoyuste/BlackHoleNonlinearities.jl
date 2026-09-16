isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function low_frequency()
    channels = [
        ("eee", 2, 2, 2, -2, 2), ("eeo", 2, 2, 3, -2, 2),
        ("eoe", 2, 2, 2, -2, 3), ("eoo", 2, 2, 2, -2, 2),
        ("ooe", 2, 2, 2, -2, 2), ("ooo", 2, 2, 3, -2, 2),
    ]
    ωout = exp10.(range(log10(0.05), log10(0.30), length=5))
    p1 = plot(xscale=:log10, yscale=:log10, ylabel="|Qin|", legend=:topleft)
    p2 = plot(xscale=:log10, xlabel="M omegaout", ylabel="arg(Qin) / pi", legend=false)
    for (sector, l1, m1, l2, m2, l) in channels
        q = [quick_q(sector, (l1, m1, w/2), (l2, m2, w/2), l) for w in ωout]
        plot!(p1, ωout, abs.(q), marker=:circle, label=uppercase(sector))
        plot!(p2, ωout, angle.(q) ./ π, marker=:circle)
    end
    ref = abs(quick_q("eee", (2, 2, ωout[1]/2), (2, -2, ωout[1]/2), 2)) .* (ωout ./ ωout[1]).^2
    plot!(p1, ωout, ref, linestyle=:dash, color=:black, label="omegaout^2")
    return write_plot(plot(p1, p2, layout=(2, 1), size=(500, 620)), "low_frequency.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && low_frequency()
