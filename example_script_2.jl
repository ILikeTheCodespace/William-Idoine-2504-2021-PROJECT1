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

p1 = 11x^2 + 52x^5
p2 = 103x^5 + 7x^3
prime = 5

println("To showcase some of the functions within this repo, we will use the two polynomials p1 = ", pretty_printer(p1), " and p2 = ", pretty_printer(p2), "\n")

println("The scripts in this repo allow you to:\n")
println("Add polynomials together:")
println("p1 + p2 = ", pretty_printer(p1+p2), "\n")

println("Multiply polynomials and raise them by a power:")
println("p1 * p2 = ", pretty_printer(p1*p2))
println("p1^3 = ", pretty_printer(p1^prime), "\n")

println("Divide one polynomial by another:")
println("p1/p2 = ", pretty_printer((p2÷p1)(prime)), "\n")

println("Find the remainder of this polynomial division")
println(rem(p1,p2)(prime), "\n")

println("Take the derivative of a polynomial:")
println("p1' * p2 = ", pretty_printer(derivative(p1*p2)))
println("p1' * p2 + p1 * p2' = ", derivative(p1)*p2 + p1*derivative(p2), "\n")

println("Find the GCD of two polynomials modulo prime:")
println("GCD(p1,p2)(mod 3) = ", pretty_printer(gcd(p1,p2,prime)), "\n")

println("Factorize polynomials:")
p = mod(p1,prime)
println("Will factor p1(mod $prime): ", p)
factorization = factor(p,prime)
println("Here is the factorization: ", factorization)

pr = mod(expand_factorization(factorization),prime)
println("Reconstructing: ", pretty_printer(pr), "\n")

println("Make a polynomial square free:")
println("The square free version of the polynomial, ", pretty_printer(p1), " is, ", pretty_printer(square_free(p1, prime)), "\n")

println("Construct a cyclotonic polynomial given a prime integer input:")
println("The following is a cyclotonic polynomial, ", pretty_printer(cyclotonic_polynomial(prime)), "\n")

println("Construct a linear monic polynomial given an integer input:")
println("The following is a linear monic polynomial generated with input $prime, ", pretty_printer(linear_monic_polynomial(prime)), "\n")

println("Make a polynomial primitive:")
println("The primitive part of the polynomial ", pretty_printer(p1), " is ", pretty_printer(prim_part(p1)(prime)), "\n")

println("Take the modulus of a polynomial:")
println(pretty_printer(p1), " modulo $prime is congruent to ", pretty_printer(mod(p1, prime)), "\n")

println("Find the degree of a polynomial:")
println("The degree of the polynomial ", pretty_printer(p1), " is equal to ", degree(p1), "\n")

println("Find the leading term of any given polynomial:")
println("The leading term of the polynomial ", pretty_printer(p1), " is equal to ", pretty_printer(Polynomial(leading(p1))), "\n")

println("Aside from the functions showcased, this software package has much more functionality! \nIt can also perform arithmetic on individual polynomial terms, manipulate polynomials \nprogrammatically with pop! and push! and much, much more!")