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
function divide(num::Polynomial, den::Polynomial)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Polynomial()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = Polynomial( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

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
÷(num::Polynomial, den::Polynomial)  = (p::Int) -> first(divide(num,den)(p))
÷(num::PolynomialModP, den::PolynomialModP)  = first(divide(num,den))

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::Polynomial, den::Polynomial)  = (p::Int) -> last(divide(num,den)(p))