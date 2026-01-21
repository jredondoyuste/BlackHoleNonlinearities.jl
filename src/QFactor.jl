import DifferentialEquations as DE
import Integrals
include("HomogeneousSolutions.jl")
include("Source.jl")

function qfactor(ω::Float64, rmin, rmax, solver_lin="default", solver_scd="verne", oo=7,atol=1e-12,rtol=1e-12)
    sol1 = linear_sol(2, ω, "odd", rmin, rmax, solver_lin,atol,rtol)
    ain1, aout1 = extract_amps(sol1, ω, rmax)
    sol2_hom = linear_sol(4, 2*ω, "even", rmin, rmax, solver_scd,atol,rtol)
    ain2, aout2 = extract_amps(sol2_hom, 2*ω, rmax)
    integrand(x, p) = sol2_hom(x)[1] * SOOE(sol1, ω, x) / (1-2/x)
    integrand_2(x) = sol2_hom(x)[1] * SOOE(sol1, ω, x) / (1-2/x)
    domain = (2 + 10.0^(-rmin), rmax)
    prob = Integrals.IntegralProblem(integrand, domain)
    sol = Integrals.solve(prob, Integrals.QuadGKJL(order=oo), abstol=1e-6, reltol=1e-6)
    prefactor = im / (4*ω*ain2)
    uncorrected=abs(prefactor*sol.u/(aout1^2))
    sol_good=sol.u+asymptotics(ω,rmax,ain1,aout1,ain2,aout2)
    corrected=abs(prefactor*sol_good/(aout1^2))
    return [uncorrected, abs(corrected-uncorrected)]
end