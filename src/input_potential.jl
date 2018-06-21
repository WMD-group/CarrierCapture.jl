# Reading in potential data and finding polynomial fit

# To make Julia call its own mini-python library (PyCall), set:
# ENV["PYTHON"] = "" then run Pkg.build("PyCall") in Julia REPL

using DataFrames
using Polynomials

######################### Import data #########################
# import data from two files
# Q (first column, amu^{1/2}Å) and E (second column, eV)

# import first potential curve and convert to matrix
potential_data_1 = readtable("Potential_1.csv")
potential_matrix_1 = convert(Matrix, potential_data_1)

# import second potential curve and convert to matrix
potential_data_2 = readtable("Potential_2.csv")
potential_matrix_2 = convert(Matrix, potential_data_2)

######################### Setting up variables #########################
# Put in configuration file later?
poly_order = 4 # order of polynomial for fitting potential

######################### Defining data ##########################
Q1 = potential_matrix_1[:,1]
E1 = potential_matrix_1[:,2]

Q2 = potential_matrix_2[:,1]
E2 = potential_matrix_2[:,2]

######################### Polynomial fit #########################
# polynomial fitting for first potential
V1 = polyfit(Q1, E1, poly_order)
# polynomial coefficients
c1 = Polynomials.coeffs(V1)

# polynomial fitting for second potential
V2 = polyfit(Q2, E2, poly_order)
# polynomial coefficients
c2 = Polynomials.coeffs(V2)

######################### Function for polynomial fit #########################

# This is now in Phonon.jl

# Inputs:
# x, e.g. x = linspace(-10,10,100)
# coefficients, polynomial coefficients
# poly_order, order of polynomial (integer)
#
# function polyfunc(x, coefficients, poly_order)
#     y_terms = zeros(length(x),poly_order + 1)
#     for i = 1:poly_order + 1
#         y_terms[:,i] = coefficients[i].*x.^(i-1)
#     end
#     return sum(y_terms,2)
# end
