"""
Module for nonlinear source terms in black hole perturbation equations.
"""

"""
    SOOE(sol, w, r)

Second-order source from odd-odd coupling for even parity.
Computed from linear solution sol at frequency w and radius r.
"""
function SOOE(sol, w, r)
    ϕ, ψ = sol(r)
    f1 = @. (1728 + 5184*r - 15984*r^2 + r^5*(25254 - 2032*w^2) + r^4*(42174 - 696*w^2) - 12*r^8*w^2*(594 + w^2) - 4*r^3*(10269 + 20*w^2) - 6*r^6*(18801 + 40*w^2) + r^7*(47790 + 2187*w^2 - 4*w^4) - 756*r^2*(1 + 3*r)^2*(-16 + 26*r - 13*r^2 + 2*r^3)*log(r))/(1890*sqrt(pi)*r^8*(1 + 3*r)^2)
    f2 = @. ((-2 + r)^2*(-288 - 1404*r - 468*r^2 - 2106*r^5 + r^3*(3546 + 4*w^2) + 3*r^4*(-81 + 8*w^2) + 126*r^2*(-13 - 69*r - 63*r^2 + 18*r^3)*log(r)))/(945*sqrt(pi)*r^7*(1 + 3*r)^2)
    f3 = @. -1/1890*((-2 + r)^2*(-168 - 852*r - 108*r^2 + 6*r^4*(855 + 2*w^2) + r^3*(2817 + 4*w^2)))/(sqrt(pi)*r^6*(1 + 3*r)^2)
    return @. f1*ϕ^2 + f2*ϕ*ψ + f3*ψ^2
end

"""
    NONREGSOURCE(sol, w, r)

Nonregularized source term for second-order equation.
Computed from linear solution sol at frequency w and radius r.
"""
function NONREGSOURCE(sol, w, r)
    ϕ, ψ = sol(r)
    f1 = @. (-56 - 176*r - 730*r^2 - 282*r^5*w^2 + 33*r^6*w^2 + 4*r^7*w^4 + 12*r^8*w^4 + r^4*(1275 - 86*w^2) + r^3*(-2134 + 56*w^2))/(70*sqrt(pi)*r^5*(1 + 3*r)^2)
    f2 = @. ((-2 + r)^2*(-14 - 53*r - 54*r^2 + 8*r^3*w^2 + 12*r^4*w^2))/(35*sqrt(pi)*r^4*(1 + 3*r)^2)
    f3 = @. ((-2 + r)^2*(-14 - 64*r - 87*r^2 + 4*r^3*w^2 + 12*r^4*w^2))/(70*sqrt(pi)*r^3*(1 + 3*r)^2)
    return @. f1*ϕ^2 + f2*ϕ*ψ + f3*ψ^2
end
