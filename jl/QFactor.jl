"""
Q-factor computation for black hole nonlinear perturbations.

    qfactor(parity, mode1, mode2, l; kwargs...) → Complex{Float64}

    parity  ::  String            — "ooo", "ooe", etc. (error if not implemented)
    mode1   ::  (l1, m1, ω1)     — Tuple{Int,Int,Float64}
    mode2   ::  (l2, m2, ω2)     — Tuple{Int,Int,Float64}
    l       ::  Int               — output angular multipole

Keyword arguments:
    rmin=5, rmax=10^5
    solver_lin="default", solver_scd="verne"
    oo=8, atol=1e-7, rtol=1e-7
    source="bruno"    # "bruno" (canonical) or "adrien" (legacy); OOO only

Implemented parities:
    "ooo" — odd×odd→odd   (Regge-Wheeler output; Bruno's source, or "adrien" legacy)
    "ooe" — odd×odd→even  (Zerilli output;       Bruno's source)

Returns Q_h (h-amplitude normalized via atoh conversion).
"""

import Integrals
include("HomogeneousSolutions.jl")
include("Source.jl")

# ---------------------------------------------------------------------------
# Normalization: ψ-amplitude → h-amplitude
# h² = ψ²/√(l(l+1))  →  atoh(l) = √(l(l+1)(l+2)(l-1)) / 2
# Q_h = atoh_factor(l1,l2,l) * Q_ψ
# ---------------------------------------------------------------------------

atoh(l::Int) = sqrt(l * (l+1) * (l+2) * (l-1)) / 2.0
atoh_factor(l1::Int, l2::Int, l::Int) = atoh(l) / (atoh(l1) * atoh(l2))

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    qfactor(parity, mode1, mode2, l; rmin, rmax, solver_lin, solver_scd, oo, atol, rtol, source)

Compute Q_h for the given parity sector. Errors if the sector is not yet implemented.
"""
function qfactor(parity::String,
                 mode1::Tuple{Int,Int,Float64},
                 mode2::Tuple{Int,Int,Float64},
                 l::Int;
                 rmin=5, rmax=10^5,
                 solver_lin="default", solver_scd="verne",
                 oo=8, atol=1e-7, rtol=1e-7,
                 source::String="bruno")
    if parity == "ooo"
        return _qfactor_ooo(mode1, mode2, l;
                            rmin=rmin, rmax=rmax,
                            solver_lin=solver_lin, solver_scd=solver_scd,
                            oo=oo, atol=atol, rtol=rtol, source=source)
    elseif parity == "ooe"
        return _qfactor_ooe(mode1, mode2, l;
                            rmin=rmin, rmax=rmax,
                            solver_lin=solver_lin, solver_scd=solver_scd,
                            oo=oo, atol=atol, rtol=rtol)
    else
        error("Source for $parity not yet implemented")
    end
end

# ---------------------------------------------------------------------------
# OOO (private)
# ---------------------------------------------------------------------------

function _qfactor_ooo(mode1::Tuple{Int,Int,Float64},
                      mode2::Tuple{Int,Int,Float64},
                      l::Int;
                      rmin=5, rmax=10^5,
                      solver_lin="default", solver_scd="verne",
                      oo=8, atol=1e-7, rtol=1e-7,
                      source::String="bruno")
    l1, m1, ω1 = mode1
    l2, m2, ω2 = mode2

    # Linear odd-parity solutions
    sol1 = linear_sol(l1, ω1, "odd", rmin, rmax, solver_lin, atol, rtol)
    _, aout1 = extract_amps(sol1, ω1, rmax)

    if l1 == l2 && ω1 == ω2
        sol2  = sol1
        aout2 = aout1
    else
        sol2 = linear_sol(l2, ω2, "odd", rmin, rmax, solver_lin, atol, rtol)
        _, aout2 = extract_amps(sol2, ω2, rmax)
    end

    # Homogeneous odd solution at output frequency
    ω_out = ω1 + ω2
    sol_hom = linear_sol(l, ω_out, "odd", rmin, rmax, solver_scd, atol, rtol)
    ain_hom, _ = extract_amps(sol_hom, ω_out, rmax)

    # Source selection
    S_fn = source == "adrien" ? S_OOO_adrien : S_OOO

    # Green's function integral
    integrand(x, p) = sol_hom(x)[1] * S_fn(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, x) / (1 - 2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol_int = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=atol, reltol=rtol)

    # Q_ψ → Q_h
    Q_psi = im / (2 * ω_out * ain_hom) * sol_int.u / (aout1 * aout2)
    return atoh_factor(l1, l2, l) * Q_psi
end

# ---------------------------------------------------------------------------
# OOE (private)  — odd×odd→even, Zerilli output
# ---------------------------------------------------------------------------

function _qfactor_ooe(mode1::Tuple{Int,Int,Float64},
                      mode2::Tuple{Int,Int,Float64},
                      l::Int;
                      rmin=5, rmax=10^5,
                      solver_lin="default", solver_scd="verne",
                      oo=8, atol=1e-7, rtol=1e-7)
    l1, m1, ω1 = mode1
    l2, m2, ω2 = mode2

    # Linear odd-parity (Regge-Wheeler) input solutions
    sol1 = linear_sol(l1, ω1, "odd", rmin, rmax, solver_lin, atol, rtol)
    _, aout1 = extract_amps(sol1, ω1, rmax)

    if l1 == l2 && ω1 == ω2
        sol2  = sol1
        aout2 = aout1
    else
        sol2 = linear_sol(l2, ω2, "odd", rmin, rmax, solver_lin, atol, rtol)
        _, aout2 = extract_amps(sol2, ω2, rmax)
    end

    # Homogeneous even-parity (Zerilli) solution at output frequency
    ω_out = ω1 + ω2
    sol_hom = linear_sol(l, ω_out, "even", rmin, rmax, solver_scd, atol, rtol)
    ain_hom, _ = extract_amps(sol_hom, ω_out, rmax)

    # Green's function integral with Bruno's OOE source
    integrand(x, p) = sol_hom(x)[1] * S_OOE(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, x) / (1 - 2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol_int = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=atol, reltol=rtol)

    # Q_ψ → Q_h (same atoh normalization as OOO)
    Q_psi = im / (2 * ω_out * ain_hom) * sol_int.u / (aout1 * aout2)
    return atoh_factor(l1, l2, l) * Q_psi
end

# ---------------------------------------------------------------------------
# Row-cached variants for 2D scans
# These accept a pre-computed (sol1, aout1) pair so the caller can reuse
# the l1,ω1 solution across an entire inner loop over ω2.
# ---------------------------------------------------------------------------

"""
    qfactor_row(parity, sol1, aout1, mode1, mode2, l; kwargs...) → Complex{Float64}

Like `qfactor` but accepts a pre-solved `sol1` and its outgoing amplitude
`aout1`, skipping that ODE integration. Useful when scanning over ω2 at
fixed (l1, m1, ω1): compute sol1 once, then call this in the inner loop.
"""
function qfactor_row(parity::String,
                     sol1, aout1::Complex,
                     mode1::Tuple{Int,Int,Float64},
                     mode2::Tuple{Int,Int,Float64},
                     l::Int;
                     rmin=5, rmax=10^5,
                     solver_lin="default", solver_scd="verne",
                     oo=8, atol=1e-7, rtol=1e-7,
                     source::String="bruno")
    if parity == "ooo"
        S_fn = source == "adrien" ? S_OOO_adrien : S_OOO
        out_parity = "odd"
    elseif parity == "ooe"
        S_fn = S_OOE
        out_parity = "even"
    else
        error("Source for $parity not yet implemented")
    end

    l1, m1, ω1 = mode1
    l2, m2, ω2 = mode2

    sol2 = linear_sol(l2, ω2, "odd", rmin, rmax, solver_lin, atol, rtol)
    _, aout2 = extract_amps(sol2, ω2, rmax)

    ω_out   = ω1 + ω2
    sol_hom = linear_sol(l, ω_out, out_parity, rmin, rmax, solver_scd, atol, rtol)
    ain_hom, _ = extract_amps(sol_hom, ω_out, rmax)

    integrand(x, p) = sol_hom(x)[1] * S_fn(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, x) / (1 - 2/x)
    domain  = (2 + 10.0^(-rmin), rmax)
    prob    = Integrals.IntegralProblem(integrand, domain)
    sol_int = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=atol, reltol=rtol)

    Q_psi = im / (2 * ω_out * ain_hom) * sol_int.u / (aout1 * aout2)
    return atoh_factor(l1, l2, l) * Q_psi
end
