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
function prod_test_poly(;N::Int = 10, N_prods::Int = 20, seed::Int = 0)
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
Test product of polynomials modulo P.
"""
function prod_test_poly_mod_p(;N::Int = 10, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    sample_primes_array = [3,5,7,11,13,17,19]
    for _ in 1:N
        rand_prime = rand(sample_primes_array)
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p1mod = PolynomialModP(p1, rand_prime)
        p2mod = PolynomialModP(p2, rand_prime)
        prod = p1mod*p2mod
        @assert mod(leading(prod), rand_prime) == leading(mod(p1mod.terms*p2mod.terms, rand_prime))
    end

    for _ in 1:N
        rand_prime = rand(sample_primes_array)
        p_base = PolynomialModP(Polynomial(Term(1,0)), rand_prime)
        for _ in 1:N_prods
            p = PolynomialModP(rand(Polynomial), rand_prime)
            prod = p_base*p
            prod = PolynomialModP(prod, rand_prime)
            @assert leading(prod.terms) == leading(mod(p_base.terms*p.terms, rand_prime))
            p_base = prod
        end
    end
    println("prod_test_poly_mod_p - PASSED")
end

"""
Test power of polynomial/integer mod P
"""
function power_test_poly(;N::Int = 10, seed::Int = 0)
    sample_primes_array = [3,5,7,11,13,17,19]
    Random.seed!(seed)
    x = x_poly()
    print("Proceeding with power tests ~ [")
    for i in 1:N
        exponent = rand(10:20)
        rand_prime = rand(sample_primes_array)
        p1 = rand(Polynomial)
        out = 1
        for i in 1:exponent
            out = mod(out*p1, rand_prime)
        end
        if ^(p1, exponent, rand_prime) - out != 0
            println("--- TEST FAILED ---")
            break
        end 
        print("-")
    end
    print("]\n")
    println("power_test_poly - COMPLETED\n")
end

"""
Test power of polynomial mod P
"""
function power_test_poly_modp(;N::Int = 10, seed::Int = 0)
    sample_primes_array = [3,5,7,11,13,17,19]
    Random.seed!(seed)
    x = x_poly()
    print("Proceeding with power tests using PolynomialModP data type ~ [")
    for i in 1:N
        exponent = rand(10:20)
        rand_prime = rand(sample_primes_array)
        p1 = rand(Polynomial)
        p1mod = PolynomialModP(p1, rand_prime)
        out = 1
        for i in 1:exponent
            out = mod(out*p1mod.terms, p1mod.prime)
        end
        if ^(p1mod, exponent) - out != 0
            println("--- TEST FAILED ---")
            break
        end 
        print("-")
    end
    print("]\n")
    println("power_test_poly_modp - COMPLETED\n")
end

"""
Benchmark multiplication 
"""
function prod_benchmark(;N::Int = 50, seed::Int = 0)
    println("Begin prod_benchmark: Demonstrating performance of polynomial multiplication algorithm")
    Random.seed!(seed)
    p1 = 0
    p2 = 0
    x = x_poly()
    for i in 1:N
        p1 += rand(1:N)*x^(i)
        p2 += rand(1:N)*x^(i)
        if i%10 == 0
            print("Time taken to multiply two polynomials with $i terms:")
            @time prod = p1*p2
            @assert leading(prod) == leading(p1)*leading(p2)
        end
    end
    println("prod_benchmark COMPLETED\n")
end

function prod_modP_benchmark(;N::Int = 50, seed::Int = 0)
    println("Begin prod_modP_benchmark: Demonstrating performance of polynomial multiplication algorithm in modulo P")
    Random.seed!(seed)
    sample_primes_array = [3,5,7,11,13,17,19]
    p1 = 0
    p2 = 0
    x = x_poly()
    for i in 1:N
        rand_prime = rand(sample_primes_array)
        p1 += rand(1:N)*x^(i)
        p2 += rand(1:N)*x^(i)
        p1mod = PolynomialModP(p1, rand_prime)
        p2mod = PolynomialModP(p2, rand_prime)
        if i%10 == 0
            print("Time taken to multiply two polynomials with $i terms:")
            @time prod = p1mod*p2mod
            @assert leading(prod) == leading(mod((mod(p1mod.terms, rand_prime) * mod(p2mod.terms, rand_prime)), rand_prime))
        end
    end
    println("prod_modP_benchmark COMPLETED\n")
end

function crt_benchmark(;N::Int = 50, seed::Int = 0)
    println("Begin crt_benchmark: Demonstrating performance of polynomial multiplication using the Chinese Remainder Theorem")
    Random.seed!(seed)
    sample_primes_array = [3,5,7,11,13,17,19]
    p1 = 0
    p2 = 0
    x = x_poly()
    for i in 1:N
        x = x_poly()
        rand_prime = rand(sample_primes_array)
        p1 += rand(1:N)*x^(i)
        p2 += rand(1:N)*x^(i)
        if i%10 == 0
            print("Time taken to multiply two polynomials with $i terms:")
            @time prod = mult_poly_with_crt(p1, p2)
            @assert prod == p1 * p2 "Error: CRT does not net p1*p2"
        end
    end
    println("crt_benchmark COMPLETED\n")
end


"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly(;N::Int = 10^3,  seed::Int = 0)
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
Test derivative of polynomials modulo P (as well as product).
"""
function prod_derivative_test_poly_mod_p(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    sample_primes_array = [3,5,7,11,13,17,19]
    for _ in 1:N
        rand_prime = rand(sample_primes_array)
        p1 = rand(Polynomial)
        p2 = rand(Polynomial)
        p1mod = PolynomialModP(p1, rand_prime)
        p2mod = PolynomialModP(p2, rand_prime)
        p1d = derivative(p1)
        p2d = derivative(p2)
        p1mod_d = PolynomialModP(p1d, rand_prime)
        p2mod_d = PolynomialModP(p2d, rand_prime)
        @assert mod((p1d*p2) + (p1*p2d), rand_prime) == mod(derivative(p1mod.terms*p2mod.terms), rand_prime)
    end
    println("prod_derivative_test_poly_mod_p - PASSED")
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
function division_test_poly_mod_p(;prime::Int = 101, N::Int = 10^3, seed::Int = 0)
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
