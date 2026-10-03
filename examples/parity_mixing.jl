isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function parity_mixing()
    mode = (2, 2, 0.19)
    a = quick_q("eee", mode, mode, 4)
    b = quick_q("eoo", mode, mode, 4)
    c = quick_q("ooe", mode, mode, 4)
    x = range(-2, 2, length=61)
    y = range(-2, 2, length=61)
    # incident +m strain ∝ C+ + i C-, so the pole sits at p = +i
    z = [abs(a + 2im*(xr + im*yi)*b - (xr + im*yi)^2*c) /
         abs(1 + im*(xr + im*yi))^2 for yi in y, xr in x]
    p = heatmap(x, y, log10.(z), xlabel="Re(p)", ylabel="Im(p)",
                colorbar_title="log10 |Qin(p)|", aspect_ratio=:equal)
    scatter!(p, [0.0], [0.0], color=:white, markerstrokecolor=:black, label=false)
    return write_plot(p, "parity_mixing.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && parity_mixing()
