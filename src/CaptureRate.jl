__precompile__()

push!(LOAD_PATH,"../src/")
module CaptureRate
using Potential
using Polynomials
# using Plots

export conf_coord, cc_from_dict, calc_overlap!, calc_capt_coeff!

ħ = 6.582119514E-16 # eV⋅s
kB = 8.6173303E-5 # eV⋅K⁻¹

occ_cut_off = 1E-5


mutable struct conf_coord
    # Configuration coordinate
    name::String
    # potentials
    V1::potential; V2::potential
    # e-ph coupling matrix element; degeneracy
    W::Float64; g::Int
    # vibrational wave function overlap integral
    # (initial state) phonon eigenvalue; phonon overlap; Gaussian function energy 
    ϵ_list::Array{Float64,1}; overlap_list::Array{Float64,1}; δ_list::Array{Float64,1}
    temp::Array{Float64,1}; capt_coeff::Array{Float64,1}
end
conf_coord(pot_i::potential, pot_f::potential) = conf_coord("", pot_i, pot_f, Inf, 1, [], [], [], [], [])


function cc_from_dict(pot_i, pot_f, cfg::Dict)::conf_coord
    name = cfg["initial"]*" => "*cfg["final"]
    cc = conf_coord(name, pot_i, pot_f, cfg["W"], cfg["g"], [], [], [], [], [])
    return cc
end

"""
function calc_harm_wave_func(ħω1, ħω2, ΔQ, ΔE; Qi=-10, Qf=10, NQ=100, nev=20, nev2=Nothing)
    if nev2 == Nothing
        nev2 = nev
    end

    # potentials
    x = range(Qi, stop=Qf, length=NQ)

    # define potential
    # Ground state
    E1 = harmonic(x, ħω1)
    V1 = potential(x, E1)
    # Excited state
    E2 = harmonic(x .- ΔQ, ħω2) .+ ΔE
    V2 = potential(x, E2)

    # solve Schrödinger equation
    # Ground state
    ϵ1, χ1 = solve1D_ev_amu(x->harmonic(x, ħω1), NQ=NQ, Qi=Qi, Qf=Qf, nev=nev)
    # Excited state
    ϵ2, χ2 = solve1D_ev_amu(x->harmonic(x .- ΔQ, ħω2), NQ=NQ, Qi=Qi, Qf=Qf, nev=nev2)
    ϵ2 = ϵ2 .+ ΔE

    # Assign
    cc = CC()
    cc.V1 = V1; cc.V2 = V2
    cc.ϵ1 = ϵ1; cc.χ1 = χ1
    cc.ϵ2 = ϵ2; cc.χ2 = χ2
    cc
end 


function calc_poly_wave_func(potential_matrix_1, potential_matrix_2, poly_order, Qi=-10, Qf=10, NQ=100, nev=10, nev2=Nothing)
    if nev2 == Nothing
        nev2 = nev
    end
    ######################### Defining data ##########################
    #     Q, Configuration coordinate
    #     E, Energy

    # data for first potential
    Q1 = potential_matrix_1[:,1]
    Energies_1 = potential_matrix_1[:,2]

    # data for second potential
    Q2 = potential_matrix_2[:,1]
    Energies_2 = potential_matrix_2[:,2]

    ######################### Polynomial fit #########################
    # polynomial fitting for first potential
    poly1 = polyfit(Q1, Energies_1, poly_order)
    # polynomial coefficients
    c1 = Polynomials.coeffs(poly1)

    # polynomial fitting for second potential
    poly2 = polyfit(Q2, Energies_2, poly_order)
    # polynomial coefficients
    c2 = Polynomials.coeffs(poly2)

    # x coordinate vector
    x = LinRange(Qi, Qf, NQ)

    # define potential using coefficients (calls polyfunc from Potential)
    E1 = polyfunc(x, c1, poly_order)
    V1 = potential(x, E1)

    E2 = polyfunc(x, c2, poly_order)
    V2 = potential(x, E2)

    # calculate ΔQ and ΔE here
    Q1min = Q1[argmin(Energies_1)]
    Q2min = Q2[argmin(Energies_2)]
    ΔQ = Q2min - Q1min
    ΔE = minimum(Energies_2) - minimum(Energies_1)

    # solve Schrödinger equation
    # Ground state
    ϵ1, χ1 = solve1D_ev_amu(x->polyfunc(x, c1, poly_order); NQ=NQ, Qi=Qi, Qf=Qf, nev=nev)
    # Excited state
    ϵ2, χ2 = solve1D_ev_amu(x->polyfunc(x, c2, poly_order); NQ=NQ, Qi=Qi, Qf=Qf, nev=nev2)

    # Assign
    cc = CC()
    cc.V1 = V1; cc.V2 = V2
    cc.ϵ1 = ϵ1; cc.χ1 = χ1
    cc.ϵ2 = ϵ2; cc.χ2 = χ2
    cc
end


function plot_potentials(cc::CC; plt=Nothing, scale_factor=2e-2)
    if plt == Nothing
       plt = plot()
    end
    # Initial state
    plot_potential(cc.V1.Q, cc.V1.E, cc.ϵ1, cc.χ1, plt=plt, color="#bd0026", scale_factor=scale_factor)

    # Final state
    plot_potential(cc.V2.Q, cc.V2.E, cc.ϵ2, cc.χ2, plt=plt, color="#4575b4", scale_factor=scale_factor)
end
"""

function calc_overlap!(cc::conf_coord; cut_off=0.25, σ=0.025)
    ΔL = (maximum(cc.V1.Q) - minimum(cc.V1.Q))/length(cc.V1.Q)
    cc.ϵ_list = []
    cc.overlap_list = []
    cc.δ_list = []
    for i in UnitRange(1, length(cc.V1.ϵ))
        for j in UnitRange(1, length(cc.V2.ϵ))
            Δϵ = abs(cc.V1.ϵ[i] - cc.V2.ϵ[j])
            if  Δϵ < cut_off
                integrand = (cc.V1.χ[i] .* cc.V1.Q .* cc.V2.χ[j]) 
                overlap = sum(integrand)*ΔL

                append!(cc.ϵ_list, cc.V1.ϵ[i])
                append!(cc.overlap_list, overlap)

                append!(cc.δ_list, exp(-(Δϵ/σ)^2/2)/(σ*sqrt(2*π)))
            end
        end
    end
end

function calc_capt_coeff!(cc::conf_coord, V, temp)
    # TODO:       convergence over σ
    capt_coeff = zeros(length(temp))
    Z = 0
    β = 1 ./ (kB .* temp)
    for ϵ in cc.V1.ϵ
        Z = Z .+ exp.(-β*ϵ)
    end
    for var in zip(cc.ϵ_list, cc.overlap_list, cc.δ_list)
        ϵ, overlap, δ = var
        occ = exp.(-β*ϵ) ./ Z
        capt_coeff += occ * overlap .* overlap * δ
    end
    capt_coeff = V*2*π/ħ*cc.g*cc.W^2 * capt_coeff

    occ_high = exp(-β[length(Z)]*cc.V1.ϵ[length(cc.V1.ϵ)] ./ Z[length(Z)])
    
    @assert occ_high < occ_cut_off "occ(ϵ_max, T_max): $occ_high should be less than $occ_cut_off"
    
    cc.capt_coeff = capt_coeff
    cc.temp = temp
end


end # module
