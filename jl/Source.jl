include("generated/SourceEEE.jl")
include("generated/SourceEEO.jl")
include("generated/SourceEOE.jl")
include("generated/SourceEOO.jl")
include("generated/SourceOOE.jl")
include("generated/SourceOOO.jl")

const SECTORS = ("eee", "eeo", "eoe", "eoo", "ooe", "ooo")

function coefficients(sector::String, l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "eee" && return coefficients_eee(l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "eeo" && return coefficients_eeo(l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "eoe" && return coefficients_eoe(l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "eoo" && return coefficients_eoo(l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "ooe" && return coefficients_ooe(l1, l2, l, m1, m2, ω1, ω2; kw...)
    sector == "ooo" && return coefficients_ooo(l1, l2, l, m1, m2, ω1, ω2; kw...)
    error("Source for $sector not implemented")
end

function source(sector::String, cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "eee" && return source_eee(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "eeo" && return source_eeo(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "eoe" && return source_eoe(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "eoo" && return source_eoo(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "ooe" && return source_ooe(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    sector == "ooo" && return source_ooo(cf, ϕ1, ψ1, ϕ2, ψ2, r)
    error("Source for $sector not implemented")
end

function make_source(sector, sol1, sol2, l1, l2, l, m1, m2, ω1, ω2; kw...)
    cf = coefficients(sector, l1, l2, l, m1, m2, ω1, ω2; kw...)
    # The generated sources use KK with the bare 3j symbol; projecting onto
    # Y^{l m} adds (-1)^m, m = m1 + m2 (Eq. (39) of the paper).
    sgn = iseven(m1 + m2) ? 1.0 : -1.0
    # float(r): integer r overflows Int64 in the r^n terms of the generated sources
    r -> (r = float(r); u1 = sol1(r); u2 = sol2(r); sgn * source(sector, cf, u1[1], u1[2], u2[1], u2[2], r))
end
