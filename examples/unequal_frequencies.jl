isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function unequal_frequencies()
    ω = range(0.10, 0.42, length=7)
    # exchange-summed kernel, smooth across the diagonal (no identical-parent 1/2)
    q = [quick_q("eee", (2, 2, w1), (2, 2, w2), 4; symmetry_factor=false) for w2 in ω, w1 in ω]
    p1 = heatmap(ω, ω, log10.(abs.(q)), xlabel="M omega1", ylabel="M omega2",
                 colorbar_title="log10 |Qin|", aspect_ratio=:equal)
    p2 = heatmap(ω, ω, angle.(q) ./ π, xlabel="M omega1", ylabel="M omega2",
                 colorbar_title="arg(Qin) / pi", aspect_ratio=:equal)
    return write_plot(plot(p1, p2, layout=(1, 2), size=(900, 380)), "unequal_frequencies.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && unequal_frequencies()
