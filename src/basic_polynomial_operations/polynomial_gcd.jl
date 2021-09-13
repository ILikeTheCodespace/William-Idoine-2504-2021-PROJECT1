#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a::Polynomial, b::Polynomial, prime::Int)
    old_r, r = mod(a,prime), mod(b,prime)
    old_s, s = one(Polynomial), zero(Polynomial)
    old_t, t = zero(Polynomial), one(Polynomial)

    while !iszero(mod(r,prime))
        q = divide(old_r, r)(prime) |> first
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::Polynomial, b::Polynomial, prime::Int) = extended_euclid_alg(a,b,prime) |> first

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a::PolynomialModP, b::PolynomialModP)
    @assert a.prime == b.prime
    old_r, r = mod(a.terms,a.prime), mod(b.terms,b.prime)
    old_s, s = one(Polynomial), zero(Polynomial)
    old_t, t = zero(Polynomial), one(Polynomial)

    while !iszero(mod(r,a.prime))
        q = divide(PolynomialModP(old_r, a.prime), PolynomialModP(r, a.prime)) |> first
        old_r, r = r, mod(old_r - q*r, a.prime)
        old_s, s = s, mod(old_s - q*s, a.prime)
        old_t, t = t, mod(old_t - q*t, a.prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a.terms + t*b.terms - g, a.prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialModP, b::PolynomialModP) = extended_euclid_alg(a,b) |> first