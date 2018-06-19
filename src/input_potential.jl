# Reading in potential data and finding polynomial fit

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
Q1_data = potential_matrix_1[:,1]
E1_data = potential_matrix_1[:,2]

Q2_data = potential_matrix_2[:,1]
E2_data = potential_matrix_2[:,2]

######################### Polynomial fit #########################
# polynomial fitting for first potential
potential_energy_surface_1 = polyfit(Q1_data, E1_data, poly_order)
# polynomial coefficients
c_i_1 = Polynomials.coeffs(potential_energy_surface_1)

# polynomial fitting for second potential
potential_energy_surface_2 = polyfit(Q2_data, E2_data, poly_order)
# polynomial coefficients
c_i_2 = Polynomials.coeffs(potential_energy_surface_2)

######################### Function for polynomial fit #########################

# Inputs:
# x, e.g. x = linspace(-10,10,100)
# coefficients, polynomial coefficients
# poly_order, order of polynomial (integer)

function polyfunc(x, coefficients, poly_order)
    y_terms = zeros(length(x),poly_order + 1)
    for i = 1:poly_order + 1
        y_terms[:,i] = coefficients[i].*x.^(i-1)
    end
    return sum(y_terms,2)
end
