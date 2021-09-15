#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations
#                                                                               
#############################################################################
#############################################################################


"""
Test product of polynomials.
"""
function prod_test_poly(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end

    for _ in 1:N
        p_base = Polynomial(Term(1,0))
        for _ in 1:N_prods
            p = rand(Polynomial)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly - PASSED")
end

"""
Benchmark multiplication 
"""
function prod_benchmark(;N::Int = 500, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    p1 = 0
    p2 = 0
    for i in 1:N
        x = x_poly()
        p1 += rand(1:100)*x^(rand(1:100))
        p2 += rand(1:100)*x^(rand(1:100))
        if i%100 == 0
            print("Time taken to multiply two polynomials with $i terms:")
            @time prod = p1*p2
            @assert leading(prod) == leading(p1)*leading(p2)
        end
    end
    println("prod_benchmark COMPLETED")
end

function crt_benchmark(;N::Int = 500, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    p1 = 0
    p2 = 0
    for i in 1:N
        x = x_poly()
        p1 += rand(1:100)*x^(rand(1:100))
        p2 += rand(1:100)*x^(rand(1:100))
        if i%100 == 0
            print("Time taken to multiply two polynomials with $i terms:")
            @time prod = p1*p2
            @assert leading(prod) == leading(p1)*leading(p2)
        end
    end
    println("prod_benchmark COMPLETED")
end


"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly - PASSED")
end


"""
Test division of polynomials modulo p.
"""
## REMOVE
function division_test_poly(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p_prod = p1*p2
        q, r = Polynomial(), Polynomial()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_poly - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly - PASSED")
end

"""
Test division of polynomials modulo p using PolynomialModP data type.
"""
function division_test_poly_mod_p(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p1mod = PolynomialModP(p1, prime)
        p2mod = PolynomialModP(p2, prime)
        p_prod = PolynomialModP(p1mod.terms*p2mod.terms, prime)
        q, r = Polynomial(), Polynomial()
        try
            q, r = divide(p_prod, p2mod)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2mod.terms, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero(mod(q*p2mod.terms+r - p_prod.terms, prime))
    end
    println("division_test_poly_mod_p - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly_mod_p(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p1mod = PolynomialModP(p1, prime)
        p2mod = PolynomialModP(p2, prime)
        g, s, t = extended_euclid_alg(p1mod, p2mod)
        @assert mod(s*p1mod.terms + t*p2mod.terms - g, prime) == 0
    end
    println("ext_euclid_test_poly_mod_p - PASSED")
end
