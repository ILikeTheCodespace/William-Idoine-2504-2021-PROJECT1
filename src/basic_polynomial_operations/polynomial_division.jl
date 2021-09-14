#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::PolynomialModP, den::PolynomialModP)
    @assert num.prime == den.prime "Both polynomials must share the same prime divisor"
    f, g = mod(num.terms,num.prime), mod(den.terms,den.prime)
    degree(f) < degree(num.terms) && return nothing 
    iszero(g) && throw(DivideError())
    q = Polynomial()
    prev_degree = degree(f)
    while degree(f) ≥ degree(g) 
        h = Polynomial( (leading(f) ÷ leading(g))(num.prime) )  #syzergy 
        f = mod((f - h*g), num.prime)
        q = mod((q + h), num.prime)  
        prev_degree == degree(f) && break
        prev_degree = degree(f)
    end
    @assert iszero( mod((num.terms  - (q*g + f)),num.prime))
    return q, f
end

"""
The quotient from polynomial division. Returns a function of an integer.
"""

÷(num::PolynomialModP, den::PolynomialModP)  = first(divide(num,den))

"""
The remainder from polynomial division. Returns a function of an integer.
"""

rem(num::PolynomialModP, den::PolynomialModP)  = last(divide(num,den))