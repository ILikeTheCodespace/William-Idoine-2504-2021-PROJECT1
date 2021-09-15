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

# """
# Multiply two polynomials.
# """
# function *(p1::PolynomialModP, p2::PolynomialModP)::Polynomial
#     @assert p1.prime == p2.prime
#     p_out = Polynomial()
#     for t in p1.terms
#         p_out = p_out + (t * p2.terms)
#     end
#     return mod(p_out, p1.prime)
# end


"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

# function CRT(p1::PolynomialModP, p2::Vararg{PolynomialModP, N} where N)   
#     u = [p1.terms]
#     m = [p1.prime]
#     for i in 1:length(p2)
#         push!(u, p2[i].terms)
#         push!(m, p2[i].prime)
#     end
#     v = Vector{Polynomial}(undef,length(u))

#     v[1] = u[1]
#     for i in 1:length(p2)
#         token = u[i+1]
#         for j in 1:i
#             token -= v[j]*mult_vec_el_up_to(m, j-1)
#         end
#         v[i+1] = mod(token*inverse_mod(mult_vec_el_up_to(m, i) , m[i+1]), m[i+1])
#         # v[i+1] = token*inverse_mod(mult_vec_el_up_to(m, i) , m[i+1]) % m[i+1] 
#     end

#     sol = 0
#     for i in 1:length(u)
#         sol += v[i]*mult_vec_el_up_to(m, i-1)
#     end
#     if mod(sol,m[1]) == u[1]
#         if mod(sol, m[2]) == u[2]
#             println("HAHAYES!")
#         end
#     end
#     return sol
# end

function CRT(p1::PolynomialModP, p2::Vararg{PolynomialModP, N} where N)   
    prods = Vector{Polynomial}(undef,3)
    m = [3,5,7]
    for i in 1:length(m)
        prod = mod(p1.terms, m[i])
        for arg in p2
            prod *= mod(arg.terms, m[i])
        end
        prods[i] = mod(prod, m[i])
    end

    # SO I NOW HAV PRODUCTS IN THE ARRAY PRODS IN THIS FORM, TOMORROW I NEED TO GET THE COEFFS OUT OF THEM AND DO THE FUNNY CRT ON THEM. https://gyazo.com/f28d7371de813b0700907f5a4f88bfe0 

    v = Vector{Int}(undef,3)
    v[1] = u[1]
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

inverse_mod(a,m) = mod(i_ext_euclid_alg(a,m)[2],m);
