#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################


"""
Test factorization of polynomials.
"""
function factor_test_poly(;N::Int = 25, seed::Int = 0, primes::Vector{Int} = [5,7,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(Polynomial)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly - PASSED")
end

"""
Test factorization of polynomials modulo p.
"""
function factor_test_poly_mod_p(;N::Int = 25, seed::Int = 0, primes::Vector{Int} = [5,7,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(Polynomial)
            pmod = PolynomialModP(p, prime)
            factorization = factor(pmod)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(pmod.terms-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_mod_p - PASSED")
end

