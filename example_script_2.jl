include("poly_factorization_project.jl")

"""
Takes some input polynomial (of type Polynomial) and converts it to a string whilst also simplifying it to make it more human readable.
"""
function pretty_printer(input::Polynomial)::String
    manipulated_string = replace(string(input), "+ -" => "- ")
    manipulated_string = replace(manipulated_string, "x^1" => "x")
    return replace(manipulated_string, "⋅x^0" => "")
end

x = x_poly()

p1 = x^2 + x^5

function testbed(p1, exponent, prime)
    out = 1
    for i in 1:exponent
        out = mod(out*p1, prime)
    end
    return mod(out, prime)
end

@time lmao = testbed(p1, 50, 13)
@time ^(p1, 50, 13)
println(^(p1, 50, 13) - lmao)

p2 = x^5 + x^3

# Class example
p1 = 3x + -4
p2 = 6x + 5

p1mod = PolynomialModP(p1, 5)
p2mod = PolynomialModP(p2, 5)

@show p1 * p2
@show CRT(p1mod, p2mod)
# p3 = x^8 + x^5
# p1+p2
# p1 * p2 + p3




# println("To showcase some of the functions within this repo, we will use the two polynomials p1 = ", pretty_printer(p1), " and p2 = ", pretty_printer(p2), "\n")

# println("The scripts in this repo allow you to:\n")
# println("Add polynomials together:")
# println("p1 + p2 = ", pretty_printer(p1+p2), "\n")

# println("Multiply polynomials and raise them by a power:")
# println("p1 * p2 = ", pretty_printer(p1*p2))
# println("p1^3 = ", pretty_printer(p1^3), "\n")

# println("Divide one polynomial by another:")
# println("p1/p2 = ", pretty_printer(divide(p2, p1)(4)[1]), "\n")

# println("Take the derivative of a polynomial:")
# println("p1' * p2 = ", pretty_printer(derivative(p1*p2)))
# println("p1' * p2 + p1 * p2' = ", derivative(p1)*p2 + p1*derivative(p2), "\n")

# println("Find the GCD of two polynomials modulo prime:")
# println("GCD(p1,p2)(mod 3) = ", pretty_printer(gcd(p1,p2,3)), "\n")

# println("Factorize polynomials:")
# prime = 5
# p = mod(p1,prime)
# println("Will factor p1(mod $prime): ", p)
# factorization = factor(p,prime)
# println("Here is the factorization: ", factorization)

# pr = mod(expand_factorization(factorization),prime)
# println("Reconstructing: ", pr)

# println("Aside from the functions showcased, this repo does much much more!")