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
function factor_test_poly(;N::Int = 25, seed::Int = 0, primes::Vector{Int} = [3,5,7,17,19])
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
function factor_test_poly_mod_p(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = primes(3,25))
    Random.seed!(seed)
    x = x_poly()
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            # Limiting to a polynomial of degree 4 since the ddfactor function scales poorly with an increase in the degree of the input polynomial
            p = rand(Polynomial, degree = rand(1:6))
            pmod = PolynomialModP(p, prime)
            factorization = factor(pmod)
            pr = mod(expand_factorization(factorization),prime)
            # if length(pr.terms) == 1 && pr.terms[1] == Term(1,0)
            #     pr = pmod.terms
            # end
            @assert mod(pmod.terms-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_mod_p - PASSED")
end