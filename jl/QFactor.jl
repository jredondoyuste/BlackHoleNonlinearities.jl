import QuadGK
include("HomogeneousSolutions.jl")
include("Source.jl")

atoh(l::Int) = sqrt(l * (l+1) * (l+2) * (l-1)) / 2.0
atoh_factor(l1::Int, l2::Int, l::Int) = atoh(l) / (atoh(l1) * atoh(l2))

# Phase of an odd master amplitude in the strain amplitude, Eq. (19) of the paper:
# -i for outgoing waves and +i for incoming ones (the TT gauge vector flips sign).
strain_phase(parity::Char, dir::Symbol) =
    parity == 'e' ? one(ComplexF64) : (dir == :out ? -im : im) * one(ComplexF64)

const QUAD_RMIN_DEFAULT = Dict("eee"=>3, "eeo"=>3, "eoe"=>5, "eoo"=>3, "ooe"=>3, "ooo"=>3)
const QUAD_PANEL_MAXEVALS = 150
const QUAD_MAX_PANELS = 500_000

_parity(c::Char) = c == 'o' ? "odd" : c == 'e' ? "even" : error("invalid parity letter '$c'")

struct QuadratureNotConverged <: Exception
    n::Int
end
Base.showerror(io::IO, e::QuadratureNotConverged) =
    print(io, "QuadratureNotConverged: did not converge within $(e.n) integrand evaluations")

function _wavelength_panels(r0::Float64, rmax::Float64, ω_out::Float64)
    λ_star = ω_out > 0 ? π / ω_out : Inf
    npanel_est = isfinite(λ_star) ? ceil(Int, (rmax - r0) / (λ_star * (1 - 2/rmax))) : 1
    edges = Float64[r0]
    sizehint!(edges, npanel_est + 1)
    r = r0
    while r < rmax
        dr = isfinite(λ_star) ? λ_star * (1 - 2 / r) : rmax - r
        r = min(r + dr, rmax)
        push!(edges, r)
    end
    length(edges) - 1 > QUAD_MAX_PANELS && throw(QuadratureNotConverged(length(edges) - 1))
    return edges
end

function _window_average(edges::Vector{Float64}, cum::Vector{ComplexF64},
                         k_hi::Int, W::Float64)
    B = edges[k_hi]
    A = B - W
    k_lo = k_hi
    while k_lo > 1 && edges[k_lo - 1] >= A
        k_lo -= 1
    end
    k_lo >= k_hi && return cum[k_hi], B
    acc = zero(ComplexF64)
    for k in k_lo:(k_hi - 1)
        acc += (cum[k] + cum[k+1]) * 0.5 * (edges[k+1] - edges[k])
    end
    span = edges[k_hi] - edges[k_lo]
    return acc / span, (edges[k_hi] + edges[k_lo]) * 0.5
end

function _qfactor(sector::String,
                  sol1, amps1,
                  mode1::Tuple{Int,Int,Float64},
                  mode2::Tuple{Int,Int,Float64},
                  l::Int;
                  rmin=5, rmax=10^5,
                  solver_scd="verne",
                  oo=8, atol=1e-10, rtol=1e-10,
                  quad_rmin=nothing, quad_rtol=1e-10,
                  quad_maxevals=20_000_000,
                  symmetry_factor=true,
                  cfree::NamedTuple=(;))
    sector in SECTORS || error("sector \"$sector\" not implemented (valid: $(SECTORS))")
    quad_rmin = something(quad_rmin, QUAD_RMIN_DEFAULT[sector])

    ain1, aout1 = amps1
    l1, m1, ω1 = mode1
    l2, m2, ω2 = mode2

    if l1 == l2 && ω1 == ω2 && sector[1] == sector[2]
        sol2 = sol1
        ain2, aout2 = ain1, aout1
    else
        sol2 = linear_sol(l2, ω2, _parity(sector[2]), rmin, rmax, solver_scd, atol, rtol)
        ain2, aout2 = extract_amps(sol2, ω2, rmax, _parity(sector[2]))
    end

    ω_out = ω1 + ω2
    sol_hom = linear_sol(l, ω_out, _parity(sector[3]), rmin, rmax, solver_scd, atol, rtol)
    ain_hom, aout_hom = extract_amps(sol_hom, ω_out, rmax, _parity(sector[3]))

    S_fn = make_source(sector, sol1, sol2, l1, l2, l, m1, m2, ω1, ω2; cfree...)

    n_evals = Ref(0)
    function integrand(x)
        n_evals[] += 1
        n_evals[] > quad_maxevals && throw(QuadratureNotConverged(n_evals[]))
        return sol_hom(x)[1] * S_fn(x) / (1 - 2/x)
    end

    r0 = 2 + 10.0^(-quad_rmin)
    edges = _wavelength_panels(r0, Float64(rmax), ω_out)
    npanel = length(edges) - 1

    cum = Vector{ComplexF64}(undef, npanel + 1)
    cum[1] = zero(ComplexF64)
    for i in 1:npanel
        v, _ = QuadGK.quadgk(integrand, edges[i], edges[i+1];
                             rtol=quad_rtol, order=oo, maxevals=QUAD_PANEL_MAXEVALS)
        cum[i+1] = cum[i] + v
    end
    total = cum[end]

    # Remove the leading 1/R tail using wavelength averages at rmax/2 and rmax.
    if ω_out > 0
        W = π / ω_out
        k_hi = npanel + 1
        k_lo = searchsortedlast(edges, rmax / 2)
        if k_lo >= 2 && edges[k_lo] > edges[1] + W && edges[k_hi] - edges[k_lo] > W
            Ihi, Rhi = _window_average(edges, cum, k_hi, W)
            Ilo, Rlo = _window_average(edges, cum, k_lo, W)
            if Rhi > Rlo
                total = Ihi + (Ihi - Ilo) * Rlo / (Rhi - Rlo)
            end
        end
    end

    # W = 2iω' a_in for the convention used by the Mathematica source.
    # The source is the exchange-summed bilinear kernel; a single parent mode
    # (equal l, m, ω) carries the usual 1/2.
    sym = symmetry_factor && (l1, m1, ω1) == (l2, m2, ω2) ? 0.5 : 1.0
    A2 = sym * atoh_factor(l1, l2, l) * strain_phase(sector[3], :out) *
         (-im) / (2 * ω_out * ain_hom) * total
    Qout = A2 / (strain_phase(sector[1], :out) * aout1 * strain_phase(sector[2], :out) * aout2)
    Qin  = A2 / (strain_phase(sector[1], :in) * ain1 * strain_phase(sector[2], :in) * ain2)
    return (Qout=Qout, Qin=Qin, A2=A2, ain1=ain1, aout1=aout1,
            ain2=ain2, aout2=aout2, ain_hom=ain_hom, aout_hom=aout_hom)
end

function qfactor_full(sector::String,
                 mode1::Tuple{Int,Int,Float64},
                 mode2::Tuple{Int,Int,Float64},
                 l::Int;
                 rmin=5, rmax=10^5,
                 solver_lin="default", solver_scd="verne",
                 oo=8, atol=1e-10, rtol=1e-10,
                 quad_rmin=nothing, quad_rtol=1e-10,
                 quad_maxevals=20_000_000,
                 symmetry_factor=true,
                 cfree::NamedTuple=(;))
    sector in SECTORS || error("sector \"$sector\" not implemented (valid: $(SECTORS))")

    l1, m1, ω1 = mode1
    sol1 = linear_sol(l1, ω1, _parity(sector[1]), rmin, rmax, solver_lin, atol, rtol)
    amps1 = extract_amps(sol1, ω1, rmax, _parity(sector[1]))

    return _qfactor(sector, sol1, amps1, mode1, mode2, l;
                    rmin=rmin, rmax=rmax, solver_scd=solver_scd,
                    oo=oo, atol=atol, rtol=rtol,
                    quad_rmin=quad_rmin, quad_rtol=quad_rtol,
                    quad_maxevals=quad_maxevals,
                    symmetry_factor=symmetry_factor,
                    cfree=cfree)
end

function qfactor(sector::String,
                 mode1::Tuple{Int,Int,Float64},
                 mode2::Tuple{Int,Int,Float64},
                 l::Int;
                 kwargs...)
    return qfactor_full(sector, mode1, mode2, l; kwargs...).Qout
end
