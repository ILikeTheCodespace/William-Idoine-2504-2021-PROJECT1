#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials.
"""

function *(p1::Polynomial, p2::Polynomial)::Polynomial
    p_out = Polynomial()
    for t in p1
        p_out = p_out + (t * p2)
    end
    return p_out
end

"""
Multiply two polynomials modulo P.
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::Polynomial
    @assert p1.prime == p2.prime "Primes must be the same"
    p_out = Polynomial()
    for t in p1.terms
        p_out = p_out + (t * p2.terms)
    end
    return mod(p_out, p1.prime)
end


"""
Power of a polynomial.
"""
# function ^(p::Polynomial, n::Int)
#     n < 0 && error("No negative power")
#     out = one(p)
#     for _ in 1:n
#         out *= p
#     end
#     return out
# end

function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    binary_arr = digits(n, base=2)
    length(binary_arr) == 1 && return p*binary_arr[1]
    for i in 1:length(binary_arr)
        binary_arr[i] == 1 ? out = out*p : out *= 1
        p = p*p
    end
    return out
end

# function ^(p, n::Int, prime::Int)
#     n < 0 && error("No negative power")
#     out = one(p)
#     binary_arr = digits(n, base=2)
#     length(binary_arr) == 1 && return mod(p*binary_arr[1], prime)
#     for i in 1:length(binary_arr)
#         binary_arr[i] == 1 ? out = mod(out*p, prime) : out *= 1
#         p = mod(p*p, prime)
#     end
#     return out
# end

# I decided to keep the ^ as a general function and not restrict it to only polynomials since this also works with raising integers and other types to a power n mod p.
function ^(p, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    binary_arr = digits(n, base=2)
    length(binary_arr) == 1 && return p*binary_arr[1]
    p = mod(p, prime)
    for i in 1:length(binary_arr)
        binary_arr[i] == 1 ? out = mod(out*p, prime) : out *= 1
        p = mod(p*p, prime)
    end
    return mod(out, prime)
end

^(p::PolynomialModP, n::Int) = ^(p.terms, n, p.prime)

function CRT(p1::PolynomialModP, p2::Vararg{PolynomialModP, N} where N)   
    x = x_poly()
    function inner_CRT_smod(u::Vector{Int64}, m::Vector{Int64}) 
        @assert length(u) == length(m) "Number of elements in input vectors not equal"
        v = Vector{Int}(undef,length(m))
        v[1] = u[1]
        for i in 1:length(m)-1
            token = u[i+1]
            for j in 1:i
                token -= v[j]*mult_vec_el_up_to(m, j-1)
            end
            v[i+1] = mod(token*inverse_mod(mult_vec_el_up_to(m, i) , m[i+1]), m[i+1])
        end

        sol = 0
        for i in 1:length(u)
            sol += v[i]*mult_vec_el_up_to(m, i-1)
        end
        for i in 1:length(m)
            @assert sol % m[i] == u[i] "Did not find solution, terminating now." 
        end
        return smod(sol, prod(m))
    end

    m = [3,5,7]
    prods = Vector{Polynomial}(undef,length(m))

    @time for i in 1:length(m)
        prod = mod(p1.terms, m[i])
        for arg in p2
            prod *= mod(arg.terms, m[i])
        end
        prods[i] = mod(prod, m[i])
    end

    poly_data = vec([vec(zeros(Int, 1, maximum(degree, prods)+1)) for _ in 1:1, _ in 1:length(m)])

    for i in 1:length(prods)
        for j in prods[i]
            poly_data[i][j.degree + 1] = j.coeff
        end
    end

    output_poly = 0
    for i in 1:maximum(degree, prods)+1
        CRT_input = vec(zeros(Int, 1, length(m)))
        for j in 1:length(m)
            CRT_input[j] = poly_data[j][i]
        end
        output_poly += inner_CRT_smod(CRT_input, m)x^(i-1)
    end

    # for i in 1:length(poly_data)
    #     output_poly += inner_CRT_smod(poly_data[i], m)x^(i-1)
    # end
    return output_poly
end

function mult_vec_el_up_to(arr::Vector, el::Int)
    @assert length(arr) > 0 "Array must contain at least 1 element." 
    el == 0 ? 1 : arr[el]*mult_vec_el_up_to(arr, el-1) 
end

function i_ext_euclid_alg(a,b)
    a == 0 && return b, 0, 1
    g, s, t = i_ext_euclid_alg(b % a, a)
    return g, t - (b ÷ a)*s, s
end

smod(a::Int,m::Int)::Int = mod(a,m) > m ÷ 2 ?  mod(a,m) - m : mod(a,m) 

inverse_mod(a,m) = mod(i_ext_euclid_alg(a,m)[2],m);
