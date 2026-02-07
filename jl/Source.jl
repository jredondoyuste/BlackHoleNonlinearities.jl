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
    f1 = @. (26244*r^4 - 80*w^2 - 240*r*w^2 - 216*r^3*(262 + r^2*w^2) - 9*r^2*(322 + r^2*w^2) - 4*(-927 - 29*r^2*w^2 + r^4*w^4) - 6*r*(-3393 - 2*r^2*w^2 + 2*r^4*w^4))/(1890*sqrt(pi)*r^5*(1 + 3*r)^2)
    f2 = @. -1/945*((-2 + r)^2*(639 + 3186*r + 4374*r^2 - 4*w^2 - 24*r*w^2))/(sqrt(pi)*r^4*(1 + 3*r)^2)
    f3 = @. -1/1890*((-2 + r)^2*(117 + 162*r + 4*w^2 + 12*r*w^2))/(sqrt(pi)*r^3*(1 + 3*r)^2)
    return @. f1*ϕ^2 + f2*ϕ*ψ + f3*ψ^2
end

function SOOE3(sol, w, r)
    ϕ, ψ = sol(r)
    f1 = @. (2769600 + 128000000*r^4 - 1620*w^2 - 12204*r*w^2 + 27*r^2*w^2 - 81*r^4*w^4 - 400*r^2*(-129391 + 9*r^2*w^2) - 16000*r^3*(18235 + 23*r^2*w^2) + 20*r*(1907120 + 162*r^2*w^2 - 27*r^4*w^4))/(12320*sqrt(13)*sqrt(pi)*r^5*(3 + 20*r)^2)
    f2 = @. -1/6160*((-2 + r)^2*(495600 + 5360000*r + 16000000*r^2 - 81*w^2 - 1080*r*w^2))/(sqrt(13)*sqrt(pi)*r^4*(3 + 20*r)^2)
    f3 = @. -1/12320*((-2 + r)^2*(99600 + 320000*r + 81*w^2 + 540*r*w^2))/(sqrt(13)*sqrt(pi)*r^3*(3 + 20*r)^2)
    return @. f1*ϕ^2 + f2*ϕ*ψ + f3*ψ^2
end

function SOOE4(sol, w, r)
    ϕ, ψ = sol(r)
    f1 = @. -1/25025*(-2401000000*r^4 + 720*w^2 + 9888*r*w^2 + 3675*r^2*(-314942 + 3*r^2*w^2) + 171500*r^3*(31884 + 19*r^2*w^2) + 12*(-1133125 + 483*r^2*w^2 + 3*r^4*w^4) + 70*r*(-4897165 - 78*r^2*w^2 + 6*r^4*w^4))/(sqrt(17)*sqrt(pi)*r^5*(3 + 35*r)^2)
    f2 = @. (-2*(-2 + r)^2*(2487975 + 46476500*r + 240100000*r^2 - 36*w^2 - 840*r*w^2))/(25025*sqrt(17)*sqrt(pi)*r^4*(3 + 35*r)^2)
    f3 = @. -1/25025*((-2 + r)^2*(525525 + 3001250*r + 36*w^2 + 420*r*w^2))/(sqrt(17)*sqrt(pi)*r^3*(3 + 35*r)^2)
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
