#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. Irreducibles are fixed points on the function factor.
"""

function factor(f::PolynomialModP)::Vector{Tuple{Polynomial,Int}}
    #Cantor Zassenhaus factorization
    f.prime == 3 && return [(mod(f.terms, f.prime),1)]
    f_modp = PolynomialModP(mod(f.terms, f.prime), f.prime)
    degree(f_modp.terms) ≤ 1 && return [(f_modp.terms,1)]

    # make f primitive
    ff = PolynomialModP(prim_part(f_modp.terms)(f.prime), f.prime)      
    # @show "after prim:", ff

    ff_der = PolynomialModP(derivative(prim_part(f_modp.terms)(f.prime)), f.prime)      

     # make f square-free
    squares_poly = gcd(f, ff_der)
    ff = PolynomialModP(ff ÷ PolynomialModP(squares_poly, f.prime), f.prime)
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(ff.terms).coeff
    ff = ff ÷ old_coeff       
    # @show "after monic:", ff

    dds = dd_factor(ff)

    ret_val = Tuple{Polynomial,Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k, f.prime)
        sp = map((p)->(p ÷ leading(p).coeff)(f.prime),sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f_modp,PolynomialModP(mp,f.prime))))
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f_modp.terms).coeff* one(Polynomial), 1) )
    # if f.prime == 3 && degree(mod(expand_factorization(ret_val),f.prime)) == 0 && mod(f.terms-expand_factorization(ret_val),f.prime) != 0
    #     ret_val = [(ret_val[1][1]*(Term(1,3)), ret_val[1][2])]
    # end
    return ret_val
end

"""
Expand a factorization.
"""
function expand_factorization(factorization::Vector{Tuple{Polynomial,Int}})::Polynomial 
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f
"""
function multiplicity(f::PolynomialModP, g::PolynomialModP)::Int
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(PolynomialModP((f ÷ g),g.prime), g)
end

"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) is equal to `f` (mod `prime`).
"""
function dd_factor(f::PolynomialModP)::Array{Polynomial}
    x = x_poly()
    w = deepcopy(x)
    g = Array{Polynomial}(undef,degree(f.terms)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f.terms)
        w = rem(PolynomialModP(^(w,f.prime, f.prime), f.prime), f)
        g[k] = gcd(PolynomialModP(w - x, f.prime), f)
        f = PolynomialModP((f ÷ PolynomialModP(g[k], f.prime)), f.prime)
    end

    #edge case for final factor
    f.terms != one(Polynomial) && push!(g,f.terms)
    
    return g
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""

function dd_split(f::Polynomial, d::Int, prime::Int)::Vector{Polynomial}
    f = PolynomialModP(mod(f,prime), prime)
    degree(f.terms) == d && return [f.terms]
    degree(f.terms) == 0 && return []
    w = PolynomialModP(mod(rand(Polynomial, degree = d, monic = true), prime), prime)
    n_power = (prime^d-1) ÷ 2
    g = PolynomialModP(gcd(PolynomialModP(^(w,n_power) - one(Polynomial), prime), f), prime)
    ḡ = (f ÷ g) # g\bar + [TAB]
    return vcat(dd_split(g.terms, d, prime), dd_split(ḡ, d, prime) )
end