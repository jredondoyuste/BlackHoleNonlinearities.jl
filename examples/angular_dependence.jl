isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function angular_dependence()
    ms = -2:2
    panels = Any[]
    for lout in (2, 3, 4)
        sector = iseven(lout) ? "eee" : "eeo"
        values = zeros(length(ms), length(ms))
        for (i, m1) in enumerate(ms), (j, m2) in enumerate(ms)
            abs(m1 + m2) > lout && continue
            q = quick_q(sector, (2, m1, 0.15), (2, m2, 0.15), lout)
            values[j, i] = abs(q)
        end
        push!(panels, heatmap(ms, ms, values, xlabel="m1", ylabel="m2",
                              title="lout = $lout", colorbar=false, aspect_ratio=:equal,
                              bottom_margin=7Plots.mm))
    end
    return write_plot(plot(panels..., layout=(1, 3), size=(900, 340)), "angular_dependence.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && angular_dependence()
