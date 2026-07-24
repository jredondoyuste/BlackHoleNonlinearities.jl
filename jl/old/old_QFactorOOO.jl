# Q-factor for OOO sector (odd×odd→odd)
import Integrals
include("HomogeneousSolutions.jl")
include("SourceOOO.jl")

"""
    qfactor_ooo(ω, l1, l2, l, m1, m2; rmin, rmax, solver_lin, solver_scd, oo, atol, rtol, source)

Compute Q-factor for the odd×odd→odd sector at frequency ω.
The 2nd-order output harmonic has angular numbers (l, m1+m2) and frequency 2ω.
The GW strain is h² = ψ²/√(l(l+1)).
`source` selects the regulated source term: "adrien" (default) or "bruno".
"""
function qfactor_ooo(ω1::Float64, ω2::Float64, l1::Int, l2::Int, l::Int, m1::Int, m2::Int;
                     rmin=5, rmax=10^5, solver_lin="default", solver_scd="verne",
                     oo=8, atol=1e-7, rtol=1e-7, source::String="adrien")
    # Linear odd-parity solutions at frequency ω
    sol1 = linear_sol(l1, ω1, "odd", rmin, rmax, solver_lin, atol, rtol)
    _, aout1 = extract_amps(sol1, ω1, rmax)
    if l1 == l2
        sol2 = sol1
        aout2 = aout1
    else
        sol2 = linear_sol(l2, ω2, "odd", rmin, rmax, solver_lin, atol, rtol)
        _, aout2 = extract_amps(sol2, ω2, rmax)
    end

    # Homogeneous ODD solution at 2ω (Regge-Wheeler, same parity as input)
    sol_hom = linear_sol(l, ω1+ω2, "odd", rmin, rmax, solver_scd, atol, rtol)
    ain_hom, _ = extract_amps(sol_hom, ω1+ω2, rmax)

    # Green's function integral
    S_fn = source == "bruno" ? SOOO_bruno : SOOO
    integrand(x, p) = sol_hom(x)[1] * S_fn(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, x) / (1 - 2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=atol, reltol=rtol)

    prefactor = im / (2 * (ω1+ω2) * ain_hom)
    return prefactor * sol.u / (aout1 * aout2)
end
