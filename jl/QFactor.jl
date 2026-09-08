"""
Q-factor computation for black hole nonlinear perturbations.
Main entry: `qfactor(sector, mode1, mode2, l; kwargs...)` / `qfactor_full(...)` → NamedTuple.

sector = XYZ: X,Y = parity of mode1/mode2, Z = parity of output ("e"=Zerilli/even, "o"=RW/odd).
mode1/mode2 = (l,m,ω)::Tuple{Int,Int,Float64}. `qfactor` returns Qout = A2/(aout1*aout2) (old
normalization). `qfactor_full` also returns Qin = A2/(ain1*ain2) (depends on the r* origin,
r* = r+2log(r/2-1), M=1), and A2, ain1, aout1, ain2, aout2, ain_hom, aout_hom. Returns Q_h
(h-amplitude, via atoh_factor).

Keyword arguments:
    rmin=5, rmax=10^5 — ODE near-horizon seed r0 = 2+10^(-rmin); keep rmin large for a good series IC.
    solver_lin="default", solver_scd="verne"
    oo=8, atol=1e-10, rtol=1e-10 — ODE solver tolerances (unrelated to quadrature).
    quad_rmin=nothing → QUAD_RMIN_DEFAULT[sector] (eoe 5, else 3) — quadrature inner cutoff
        r=2+10^(-quad_rmin); eoe needs 5 (cfree shifts Q 15% at 3, 0.2% at 5), eeo/eee must stay
        at 3 (near-horizon cancellation is noise below that).
    quad_rtol=1e-10, quad_atol=0.0 — per-panel tolerances (panel = one ODE step, see
        `_quad_panels`); quad_rtol matches the ~1e-7 noise floor inherited from the ODE solves.
        quad_atol=0 means atol_i = QUAD_ATOL_FRAC·(panel contribution + mean contribution),
        bounding total error by 2·QUAD_ATOL_FRAC·Σ|∫panel|; pass quad_atol>0 to set that bound directly.
    quad_nwave=1.0 — no panel wider than 1/quad_nwave of the local 2ω_out beat wavelength.
    quad_extrap=true — correct the rmax truncation: average over one beat period (kills the
        oscillatory tail), then Richardson-extrapolate rmax/2↔rmax (kills the 1/R algebraic tail).
        Free (reuses the panel loop's partial sums).
    quad_maxevals=20_000_000 — hard cap on Green's-function integrand evaluations; raises
        QuadratureNotConverged past this.
    cfree=(;) — free regularization constants, forwarded to make_source.
"""

import QuadGK
include("HomogeneousSolutions.jl")
include("Source.jl")

# ψ-amplitude → h-amplitude: h² = ψ²/√(l(l+1)), atoh(l) = √(l(l+1)(l+2)(l-1))/2, Q_h = atoh_factor·Q_ψ.
atoh(l::Int) = sqrt(l * (l+1) * (l+2) * (l-1)) / 2.0
atoh_factor(l1::Int, l2::Int, l::Int) = atoh(l) / (atoh(l1) * atoh(l2))

# QUAD_ATOL_FRAC: knee of the cost/accuracy curve — 1e-6 settles a panel in one Gauss-Kronrod
# rule (~17 evals at oo=8) and matches a rtol=1e-10 reference to ~5e-8 rel on |Q|; 1e-8 costs
# ~300x more evaluations for only the 8th digit (below ~1e-7 it's resolving ODE-solver noise).
const QUAD_ATOL_FRAC = 1e-6

# QUAD_PANEL_MAXEVALS=85: a panel (one ODE step) is resolved to near machine precision by the
# first Gauss-Kronrod rule (17 evals at oo=8); this caps the cost of noise-floor panels that
# can never converge further — without it, eee drove whole sweeps past 20M evals at rmax=25000.
const QUAD_PANEL_MAXEVALS = 85

# Per-sector inner cutoff r0 = 2+10^(-quad_rmin): eoe needs quad_rmin=5 (cfree shifts Q up to
# 15% at quad_rmin=3, a horizon boundary term, down to 0.2% at 5); eeo/eee must stay at 3 (the
# near-horizon source there is a cancellation, pure noise below 1e-3, Q off by 6%-100%). Other
# sectors are insensitive (<1e-4) between 1e-2 and 1e-5.
const QUAD_RMIN_DEFAULT = Dict("eee"=>3, "eeo"=>3, "eoe"=>5, "eoo"=>3, "ooe"=>3, "ooo"=>3)

_parity(c::Char) = c == 'o' ? "odd" : c == 'e' ? "even" : error("invalid parity letter '$c'")

# Thrown when the Green's-function quadrature would exceed quad_maxevals evaluations or
# _quad_panels' max_panels — a hard wall-clock safety net, not a warning.
struct QuadratureNotConverged <: Exception
    n::Int
end
Base.showerror(io::IO, e::QuadratureNotConverged) =
    print(io, "QuadratureNotConverged: did not converge within $(e.n) integrand evaluations")

# Panel edges = union of sol1/sol2/sol_hom's own ODE step grids: the integrand is their
# piecewise-poly dense output, smooth within a step but kinked at every step boundary, so a
# single global Gauss-Kronrod call over all steps used to exhaust quad_maxevals on nearly every
# channel. Panels still wider than n_per_wave⁻¹ of the local 2ω_out beat wavelength (π f(r)/ω_out
# in r, f=1-2/r, from r*=r+2log(r/2-1)) are further subdivided as a safety net.
function _quad_panels(r0::Float64, rmax::Float64, sols, ω_out::Float64;
                      n_per_wave::Float64=1.0, max_panels::Int=400_000)
    pts = Float64[r0, rmax]
    for s in sols, t in s.t
        r0 < t < rmax && push!(pts, t)
    end
    sort!(pts)
    unique!(pts)

    λ_star = ω_out > 0 ? π / ω_out : Inf     # wavelength of the 2ω_out beat in r_*
    edges = Float64[pts[1]]
    sizehint!(edges, length(pts))
    for i in 2:length(pts)
        a, b = edges[end], pts[i]
        b > a || continue
        wmax = λ_star * (1 - 2 / a) / n_per_wave    # a > 2 always
        nsub = (isfinite(wmax) && wmax > 0) ? clamp(ceil(Int, (b - a) / wmax), 1, 4096) : 1
        if nsub == 1
            push!(edges, b)
        else
            for j in 1:nsub
                push!(edges, a + (b - a) * (j / nsub))
            end
        end
    end
    length(edges) - 1 > max_panels && throw(QuadratureNotConverged(length(edges) - 1))
    return edges
end


# Trapezoid-average the partial-sum curve I(edges[k]) over [edges[k_hi]-W, edges[k_hi]] (snapped
# to panel edges, mismatch < one panel width). Averaging over one beat period kills the
# oscillatory truncation error, leaving the algebraic part for _qfactor's Richardson step.
function _window_average(edges::Vector{Float64}, cum::Vector{ComplexF64},
                         k_hi::Int, W::Float64)
    B = edges[k_hi]
    A = B - W
    k_lo = k_hi
    while k_lo > 1 && edges[k_lo - 1] >= A
        k_lo -= 1
    end
    k_lo >= k_hi && return cum[k_hi], B      # window narrower than one panel
    acc = zero(ComplexF64)
    for k in k_lo:(k_hi - 1)
        acc += (cum[k] + cum[k+1]) * 0.5 * (edges[k+1] - edges[k])
    end
    span = edges[k_hi] - edges[k_lo]
    return acc / span, (edges[k_hi] + edges[k_lo]) * 0.5
end

# Shared impl behind qfactor/qfactor_row. sol1/amps1=(ain1,aout1) must already be solved
# (parity sector[1]); solves sol2 (parity sector[2], reusing sol1 if modes coincide) and
# sol_hom (parity sector[3], ω_out=ω1+ω2), then evaluates the Green's-function integral.
function _qfactor(sector::String,
                  sol1, amps1,
                  mode1::Tuple{Int,Int,Float64},
                  mode2::Tuple{Int,Int,Float64},
                  l::Int;
                  rmin=5, rmax=10^5,
                  solver_scd="verne",
                  oo=8, atol=1e-10, rtol=1e-10,
                  quad_rmin=nothing, quad_atol=0.0, quad_rtol=1e-10, quad_nwave=1.0, quad_extrap=true,
                  quad_maxevals=20_000_000,
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

    # quad_rmin/quad_atol/quad_rtol are decoupled from the ODE near-horizon seed (rmin/atol/rtol);
    # quad_maxevals is the wall-clock safety net (see QuadratureNotConverged).
    n_evals = Ref(0)
    function integrand(x)
        n_evals[] += 1
        n_evals[] > quad_maxevals && throw(QuadratureNotConverged(n_evals[]))
        return sol_hom(x)[1] * S_fn(x) / (1 - 2/x)
    end

    r0 = 2 + 10.0^(-quad_rmin)
    edges = _quad_panels(r0, Float64(rmax), (sol1, sol2, sol_hom), ω_out;
                         n_per_wave=quad_nwave)
    npanel = length(edges) - 1

    # Panel contributions cancel heavily (running partial sum >> final answer), so per-panel
    # atol_i = QUAD_ATOL_FRAC·(w_i + mean(w)) — bounding total error by 2·QUAD_ATOL_FRAC·L1-scale —
    # instead of a flat floor, which tightens as npanel grows (doubling rmax blew eee past 20M
    # evals at rmax=25000). w_i from one midpoint sample per panel (~6% overhead); the local term
    # is a relative tolerance (rmax-independent), the mean term keeps near-zero panels from
    # bisecting forever.
    wmid = Vector{Float64}(undef, npanel)
    scale = 0.0
    for i in 1:npanel
        a, b = edges[i], edges[i+1]
        wmid[i] = abs(integrand((a + b) / 2)) * (b - a)
        scale += wmid[i]
    end
    mean_w = scale / npanel
    atol_frac = quad_atol > 0 ? quad_atol / (2 * scale) : QUAD_ATOL_FRAC

    # QUAD_PANEL_MAXEVALS caps a panel that can't meet quad_rtol at O(1) cost instead of stalling
    # (below the ~1e-7 noise floor inherited from the ODE solves, G-K is resolving solver noise).
    # Partial sums are kept at every edge for the truncation control below, at no extra cost.
    cum = Vector{ComplexF64}(undef, npanel + 1)
    cum[1] = zero(ComplexF64)
    total = zero(ComplexF64)
    for i in 1:npanel
        v, _ = QuadGK.quadgk(integrand, edges[i], edges[i+1];
                             atol=atol_frac * (wmid[i] + mean_w), rtol=quad_rtol,
                             order=oo, maxevals=QUAD_PANEL_MAXEVALS)
        total += v
        cum[i+1] = total
    end

    # rmax truncation has two errors: (1) smooth ~1/r² algebraic tail (ingoing sol_hom x outgoing
    # source beat) — ooo test |Q| drifts 0.0085177→0.0082664 (rmax 600→50000, 1.5% over
    # 1250-12500), pure truncation not quadrature error (identical to 13 digits at quad_atol 0 vs
    # 1e-12). (2) oscillatory tail at the 2ω_out beat (decays only in envelope), dominates the
    # residual cfree-dependence. Richardson alone cancels (1) but amplifies (2) by up to 3x, so
    # average over one beat period first (kills (2)), then Richardson the smoothed curve (kills
    # (1)) — both reuse the panel loop's partial sums, no extra evaluations.
    if quad_extrap
        W = π / ω_out                      # one period of the 2ω_out beat
        k_half = searchsortedlast(edges, rmax / 2)
        if k_half > 1 && edges[k_half] - W > edges[1]
            Ihi, Rhi = _window_average(edges, cum, npanel + 1, W)
            Ilo, Rlo = _window_average(edges, cum, k_half, W)
            if Rhi > Rlo
                total = Ihi + (Ihi - Ilo) * Rlo / (Rhi - Rlo)
            end
        end
    end

    # Second-order outgoing amplitude for unit-horizon-normalised linear solutions.
    A2 = atoh_factor(l1, l2, l) * im / (2 * ω_out * ain_hom) * total
    Qout = A2 / (aout1 * aout2)
    Qin  = A2 / (ain1 * ain2)
    return (Qout=Qout, Qin=Qin, A2=A2, ain1=ain1, aout1=aout1,
            ain2=ain2, aout2=aout2, ain_hom=ain_hom, aout_hom=aout_hom)
end

# Full Q-factor NamedTuple for the given sector; solves sol1 then delegates to _qfactor.
function qfactor_full(sector::String,
                 mode1::Tuple{Int,Int,Float64},
                 mode2::Tuple{Int,Int,Float64},
                 l::Int;
                 rmin=5, rmax=10^5,
                 solver_lin="default", solver_scd="verne",
                 oo=8, atol=1e-10, rtol=1e-10,
                 quad_rmin=nothing, quad_atol=0.0, quad_rtol=1e-10, quad_nwave=1.0, quad_extrap=true,
                  quad_maxevals=20_000_000,
                 cfree::NamedTuple=(;))
    sector in SECTORS || error("sector \"$sector\" not implemented (valid: $(SECTORS))")

    l1, m1, ω1 = mode1
    sol1 = linear_sol(l1, ω1, _parity(sector[1]), rmin, rmax, solver_lin, atol, rtol)
    amps1 = extract_amps(sol1, ω1, rmax, _parity(sector[1]))

    return _qfactor(sector, sol1, amps1, mode1, mode2, l;
                    rmin=rmin, rmax=rmax, solver_scd=solver_scd,
                    oo=oo, atol=atol, rtol=rtol,
                    quad_rmin=quad_rmin, quad_atol=quad_atol, quad_rtol=quad_rtol,
                    quad_nwave=quad_nwave, quad_extrap=quad_extrap,
                    quad_maxevals=quad_maxevals,
                    cfree=cfree)
end

# Qout = A2/(aout1*aout2) for the given sector (backward-compatible with old Q_h); see
# qfactor_full for the complete NamedTuple, including Qin normalised by the ingoing amplitudes.
function qfactor(sector::String,
                 mode1::Tuple{Int,Int,Float64},
                 mode2::Tuple{Int,Int,Float64},
                 l::Int;
                 kwargs...)
    return qfactor_full(sector, mode1, mode2, l; kwargs...).Qout
end

# Like qfactor_full but skips resolving sol1 — pass a pre-solved sol1 and amps1=(ain1,aout1);
# useful when scanning ω2 at fixed (l1,m1,ω1): compute sol1 once, call this in the inner loop.
function qfactor_row_full(sector::String,
                     sol1, amps1,
                     mode1::Tuple{Int,Int,Float64},
                     mode2::Tuple{Int,Int,Float64},
                     l::Int;
                     rmin=5, rmax=10^5,
                     solver_lin="default", solver_scd="verne",
                     oo=8, atol=1e-10, rtol=1e-10,
                     quad_rmin=nothing, quad_atol=0.0, quad_rtol=1e-10, quad_nwave=1.0, quad_extrap=true,
                  quad_maxevals=20_000_000,
                     cfree::NamedTuple=(;))
    sector in SECTORS || error("sector \"$sector\" not implemented (valid: $(SECTORS))")

    return _qfactor(sector, sol1, amps1, mode1, mode2, l;
                    rmin=rmin, rmax=rmax, solver_scd=solver_scd,
                    oo=oo, atol=atol, rtol=rtol,
                    quad_rmin=quad_rmin, quad_atol=quad_atol, quad_rtol=quad_rtol,
                    quad_nwave=quad_nwave, quad_extrap=quad_extrap,
                    quad_maxevals=quad_maxevals,
                    cfree=cfree)
end

# Like qfactor but with a pre-solved sol1/amps1 — see qfactor_row_full.
function qfactor_row(sector::String,
                     sol1, amps1,
                     mode1::Tuple{Int,Int,Float64},
                     mode2::Tuple{Int,Int,Float64},
                     l::Int;
                     kwargs...)
    return qfactor_row_full(sector, sol1, amps1, mode1, mode2, l; kwargs...).Qout
end
