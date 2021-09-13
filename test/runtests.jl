#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

include("../poly_factorization_project.jl")

####
# Execute unit tests for integers
###
include("integers_test.jl")
# test_euclid_ints()
# test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
# TAKEN OUT A TEST FOR TIME PURPOSES
# prod_test_poly()
# prod_derivative_test_poly()
# ext_euclid_test_poly()
# division_test_poly()
ext_euclid_test_poly_mod_p()
division_test_poly_mod_p()

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
# factor_test_poly()
factor_test_poly_mod_p()