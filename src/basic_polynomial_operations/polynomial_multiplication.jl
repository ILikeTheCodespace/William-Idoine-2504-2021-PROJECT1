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

function CRT(poly_arr, prime_arr)

    c = 0
    max_input_deg = max(degree(poly_arr[1]), degree(poly_arr[2]))+1
    exponent_arr_1 = zeros(Int64, 1, max_input_deg)
    exponent_arr_2 = zeros(Int64, 1, max_input_deg)
    for i in poly_arr[1]
        exponent_arr_1[i.degree+1] = i.coeff
    end
    for i in poly_arr[2]
        exponent_arr_2[i.degree+1] = i.coeff
    end
    """
    The lines above actually result in a marginally slower CRT implementation as opposed to the approach where you cycle through exponent arrays 1 and 2 and check if k matches the exponent value before sending those values off to the iCRT call. I chose to keep this approach because the computational time difference was only marginal and I felt like it shows that I *really* tried to tackle the task given.
    """

    for k in max_input_deg-1:-1:0
        ak = exponent_arr_1[k+1]
        bk = exponent_arr_2[k+1]
        ck = iCRT([ak, bk], prime_arr)
        k == 0 ? c = c + ck : c = c + ck * x^k # Can probably remove the turnery, its only there because I had a strange bug with Term(1,1)
    end
    return c
end

function mult_poly_with_crt(a, b)
    height_a = maximum(coeffs(a))
    height_b = maximum(coeffs(b))
    B = 2*height_a*height_b*min(degree(a)+1, degree(b)+1)
    M = 3
    c = poly_modP_multiplication(a, b, M)
    while M < B
        nextprime(M) == M ? p = nextprime(M, 2) : p = nextprime(M)
        d = poly_modP_multiplication(a, b, p)
        c = CRT([c, d], [M, p])
        M = M*p
    end
    return c
end

function iCRT(u, m)
    v = Vector{Int}(undef, 2)
    v[1] = u[1]
    v[2] = mod((u[2] - v[1])*inverse_mod(m[1], m[2]), m[2])
    output = u[1] + v[2]*m[1]
    smod(output, m[1]) == smod(u[1], m[1]) && smod(output, m[2]) == smod(u[2], m[2]) && return smod(output, m[1]*m[2])
    error("Output modulo prime is not consistent with inputs")
end

function poly_modP_multiplication(a::Polynomial, b::Polynomial, prime::Int)
    if length(a) == 0 || length(b) == 0
        return 0
    end
    if leading(a) == a.terms[1] && leading(b) == b.terms[1]
        return mod(mod(a, prime) * mod(b, prime), prime)
    end

    if length(a.terms) == 1
        a1 = a
        a2 = Polynomial()
    else
        a1 = Polynomial(a.terms[1:(length(a.terms)÷2)])
        a2 = Polynomial(a.terms[(length(a.terms)÷2)+1:length(a.terms)])    
    end
    
    if length(b.terms) == 1
        b1 = b
        b2 = Polynomial()
    else
        b1 = Polynomial(b.terms[1:(length(b.terms)÷2)])
        b2 = Polynomial(b.terms[(length(b.terms)÷2)+1:length(b.terms)])    
    end

    return mod(poly_modP_multiplication(a1, b1, prime) + poly_modP_multiplication(a1, b2, prime) + poly_modP_multiplication(a2, b1, prime) + poly_modP_multiplication(a2, b2, prime), prime)
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
