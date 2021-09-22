#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

include("../poly_factorization_project.jl")

println("Starting tests")
x = x_poly()
####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
prod_test_poly()
prod_test_poly_mod_p()
prod_derivative_test_poly()
prod_derivative_test_poly_mod_p()
ext_euclid_test_poly_mod_p()
division_test_poly_mod_p()
power_test_poly()
power_test_poly_modp()

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
@time factor_test_poly_mod_p()

###################################
# Benchmark tests for polynomials #
###################################
include("polynomials_test.jl")
prod_benchmark()
prod_modP_benchmark()
crt_benchmark()