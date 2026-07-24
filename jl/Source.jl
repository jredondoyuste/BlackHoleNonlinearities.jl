"""
Source terms for black hole perturbation equations at 2nd order.

Public API:
  S_OOO(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r)        # OOO, Bruno's (canonical)
  S_OOO_adrien(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r) # OOO, Adrien's (legacy)
  S_OOE(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r)        # OOE, Bruno's (canonical)

Other sectors (OEO, OEE, EOO, EOE, EEE) are not yet implemented.
"""

# Internal implementations — do NOT edit math expressions inside these files
include("SourceOOO.jl")
include("SourceOOE.jl")

"""
    S_OOO(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r)

Odd×odd→odd source term (canonical, Bruno's derivation).
"""
S_OOO(args...) = SOOO_bruno(args...)

"""
    S_OOO_adrien(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r)

Odd×odd→odd source term (Adrien's legacy derivation).
"""
S_OOO_adrien(args...) = SOOO(args...)

"""
    S_OOE(sol1, sol2, ω1, ω2, l1, l2, l, m1, m2, r)

Odd×odd→even source term (canonical, Bruno's derivation).
Output drives the Zerilli equation.
"""
S_OOE(args...) = SOOE_bruno(args...)
