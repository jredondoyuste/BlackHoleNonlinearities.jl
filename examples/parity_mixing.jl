isdefined(Main, :BHN_EXAMPLE_COMMON) || include("common.jl")

function parity_mixing()
    mode = (2, 2, 0.19)
    qeee = quick_q("eee", mode, mode, 4) / 2
    qeoo = quick_q("eoo", mode, mode, 4) / 2
    qooe = quick_q("ooe", mode, mode, 4) / 2
    x = range(-1.2, 1.2, length=45)
    y = range(-1.2, 1.2, length=45)
    z = [abs(qeee - 2im*(xr + im*yi)*qeoo + (xr + im*yi)^2*qooe) /
         abs(1 - im*(xr + im*yi))^2 for yi in y, xr in x]
    p = heatmap(x, y, log10.(z), xlabel="Re(p)", ylabel="Im(p)",
                colorbar_title="log10 |Qin(p)|", aspect_ratio=:equal)
    scatter!(p, [0.0], [0.0], color=:white, markerstrokecolor=:black, label=false)
    return write_plot(p, "parity_mixing.pdf")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && parity_mixing()
