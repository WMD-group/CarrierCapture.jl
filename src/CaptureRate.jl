module CaptureRate

push!(LOAD_PATH, ".")

using Potential
using Potential: polyfunc, harmonic, solve1D_ev_amu
using Polynomials
using Plots

# defining constants
ħ = 6.582119514E-16 # eV⋅s
kB = 8.6173303E-5 # eV⋅K⁻¹

mutable struct potential
    Q # Configuration coordinate
    E # Energy
end

mutable struct CC
    # Configuration coordinate
    # potentials
    V1; V2
    # eigenvalue and eigenvectors.
    ϵ1; ϵ2
    χ1; χ2
    # vibrational wave function overlap integral
    # (initial state) phonon eigenvalue; phonon overlap; Gaussian function energy
    ϵ_list; overlap_list; δ_list
end
CC() = CC([], [], [], [], [], [], [], [], [])

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
    plot!(plt, cc.V1.Q, cc.V1.E, lw=4, color="#bd0026", label="")
    for i = 1:length(cc.ϵ1)
        plot!(cc.V1.Q, cc.χ1[i]*scale_factor .+ cc.ϵ1[i], fillrange=[cc.χ1[i]*0 .+ cc.ϵ1[i], cc.χ1[i]*scale_factor .+ cc.ϵ1[i]], c="#bd0026", alpha=0.5, label="")
        # plot!(cc.V1.Q, cc.χ1[i]*scale_factor .+ cc.ϵ1[i], color="#d73027", label="")
    end

    # Final state
    plot!(cc.V2.Q, cc.V2.E, lw=4, color="#4575b4", label="")
    for i = 1:length(cc.ϵ2)
        plot!(cc.V2.Q, cc.χ2[i]*scale_factor .+ cc.ϵ2[i], fillrange=[cc.χ2[i]*0 .+ cc.ϵ2[i], cc.χ2[i]*scale_factor .+ cc.ϵ2[i]], c="#4575b4", alpha=0.5, label="")
        # plot!(cc.V2.Q, cc.χ2[i]*scale_factor .+ cc.ϵ2[i], color="#4575b4", label="")
    end
end

function calc_overlap!(cc::CC; plt=Nothing, cut_off=0.25, σ=0.025, lplot=false)
    if lplot == true
        plt = plot_potentials(cc; plt=plt)
    end

    ΔL = (maximum(cc.V1.Q) - minimum(cc.V1.Q))/length(cc.V1.Q)
    cc.ϵ_list = []
    cc.overlap_list = []
    cc.δ_list = []
    for i in UnitRange(1, length(cc.ϵ1))
        for j in UnitRange(1, length(cc.ϵ2))
            Δϵ = abs(cc.ϵ1[i] - cc.ϵ2[j])
            if  Δϵ < cut_off
                integrand = (cc.χ1[i] .* cc.V1.Q .* cc.χ2[j])
                overlap = sum(integrand)*ΔL

                append!(cc.ϵ_list, cc.ϵ1[i])
                append!(cc.overlap_list, overlap)

                append!(cc.δ_list, exp(-(Δϵ/σ)^2/2)/(σ*sqrt(2*π)))

                # plot
                # alpha = (cut_off-Δϵ)/cut_off
                # plot!(cc.V1.Q, cc.ϵ1[i] .+ integrand*1E-1, color="#31a354", lw=3, alpha=alpha, label="")
            end
        end
    end
    if lplot
        return(plt)
    end
end

function calc_capt_coeff(W, V, T_range, cc::CC)
    # TODO:       convergence over σ
    capt_coeff = zeros(length(T_range))
    Z = 0
    β = 1 ./ (kB .* T_range)
    for ϵ in cc.ϵ1
        Z = Z .+ exp.(-β*ϵ)
    end
    for summand in zip(cc.ϵ_list, cc.overlap_list, cc.δ_list)
        ϵ, overlap, δ = summand
        occ = exp.(-β*ϵ) ./ Z
        capt_coeff += occ * overlap .* overlap * δ
    end
    capt_coeff = V*2*π/ħ*g*W^2 * capt_coeff

    occ_high = exp(-β[length(Z)]*cc.ϵ1[length(cc.ϵ1)] ./ Z[length(Z)])

    println("occupation(ϵ_max, T_max): ", occ_high)

    return capt_coeff
end


end # module
