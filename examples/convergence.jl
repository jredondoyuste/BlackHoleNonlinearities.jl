isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function convergence()
    cases = [
        ("eee", 2, 2, 2, -2, 2), ("eeo", 2, 2, 3, -2, 2),
        ("eoe", 2, 2, 2, -2, 3), ("eoo", 2, 2, 2, -2, 2),
        ("ooe", 2, 2, 2, -2, 2), ("ooo", 2, 2, 3, -2, 2),
    ]
    p1 = plot(xscale=:log10, yscale=:log10, xlabel="(rmin - 2M) / M",
              ylabel="relative change", legend=false)
    p2 = plot(xscale=:log10, yscale=:log10, xlabel="rmax / M",
              ylabel="relative change", legend=:topright)
    offsets = exp10.(range(-4, -1, length=4))
    radii = exp10.(range(log10(180), log10(1200), length=4))
    for (sector, l1, m1, l2, m2, l) in cases
        w = QNM_RE[l] / 2
        modes = ((l1, m1, w), (l2, m2, w))
        reference = qfactor_full(sector, modes[1], modes[2], l;
                                 rmin=5.0, rmax=1800.0, atol=1e-8, rtol=1e-8,
                                 quad_rtol=1e-7).Qin
        qinner = [qfactor_full(sector, modes[1], modes[2], l;
                  rmin=5.0, rmax=1800.0, quad_rmin=-log10(x), atol=1e-8,
                  rtol=1e-8, quad_rtol=1e-7).Qin for x in offsets]
        qouter = [qfactor_full(sector, modes[1], modes[2], l;
                  rmin=5.0, rmax=r, atol=1e-8, rtol=1e-8,
                  quad_rtol=1e-7).Qin for r in radii]
        err_inner = max.(abs.(qinner .- reference) ./ abs(reference), 1e-12)
        err_outer = max.(abs.(qouter .- reference) ./ abs(reference), 1e-12)
        plot!(p1, offsets, err_inner, marker=:circle)
        plot!(p2, radii, err_outer,
              marker=:circle, label=uppercase(sector))
    end
    hline!(p1, [1e-2], linestyle=:dot, color=:black)
    hline!(p2, [1e-2], linestyle=:dot, color=:black, label=false)
    return write_plot(plot(p1, p2, layout=(2, 1), size=(500, 620)), "convergence.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && convergence()
