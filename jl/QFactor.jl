"""
Module for computing Q-factor of black hole nonlinear mode.
"""
# import DifferentialEquations as DE
import Integrals
include("HomogeneousSolutions.jl")
include("Source.jl")

"""
    qfactor(ω, rmin, rmax, solver_lin, solver_scd, oo, atol, rtol)

Compute Q-factor for frequency ω.
Solves linear problem, second-order problem, and integrates source term.
"""
function qfactor(ω::Float64, rmin=5, rmax=10^5, solver_lin="default", solver_scd="verne", oo=8,atol=1e-7,rtol=1e-7)
    sol1 = linear_sol(2, ω, "odd", rmin, rmax, solver_lin,atol,rtol)
    _, aout1 = extract_amps(sol1, ω, rmax)
    sol2_hom = linear_sol(4, 2*ω, "even", rmin, rmax, solver_scd,atol,rtol)
    ain2, _ = extract_amps(sol2_hom, 2*ω, rmax)
    integrand(x, p) = sol2_hom(x)[1] * SOOE(sol1, ω, x) / (1-2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=1e-6, reltol=1e-6)
    prefactor = im / (4*ω*ain2)
    return prefactor*sol.u/(aout1^2)
end

function qfactor3(ω::Float64, rmin=5, rmax=10^5, solver_lin="default", solver_scd="verne", oo=8,atol=1e-7,rtol=1e-7)
    sol1 = linear_sol(3, ω, "odd", rmin, rmax, solver_lin,atol,rtol)
    _, aout1 = extract_amps(sol1, ω, rmax)
    sol2_hom = linear_sol(6, 2*ω, "even", rmin, rmax, solver_scd,atol,rtol)
    ain2, _ = extract_amps(sol2_hom, 2*ω, rmax)
    integrand(x, p) = sol2_hom(x)[1] * SOOE3(sol1, ω, x) / (1-2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=1e-6, reltol=1e-6)
    prefactor = im / (4*ω*ain2)
    return prefactor*sol.u/(aout1^2)
end

function qfactor4(ω::Float64, rmin=5, rmax=10^5, solver_lin="default", solver_scd="verne", oo=8,atol=1e-7,rtol=1e-7)
    sol1 = linear_sol(4, ω, "odd", rmin, rmax, solver_lin,atol,rtol)
    _, aout1 = extract_amps(sol1, ω, rmax)
    sol2_hom = linear_sol(8, 2*ω, "even", rmin, rmax, solver_scd,atol,rtol)
    ain2, _ = extract_amps(sol2_hom, 2*ω, rmax)
    integrand(x, p) = sol2_hom(x)[1] * SOOE4(sol1, ω, x) / (1-2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=1e-6, reltol=1e-6)
    prefactor = im / (4*ω*ain2)
    return prefactor*sol.u/(aout1^2)
end