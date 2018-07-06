module AnharmCaptureRate

push!(LOAD_PATH, ".")

using Phonon: polyfunc, solve1D_ev_amu
# using Plots
using Polynomials
using PyPlot
######################### Defining constants #########################
ħ = 6.582119514E-16 # eV⋅s
kB = 8.6173303E-5 # eV⋅K⁻¹

######################### Defining variables #########################

V = 1.1E-21 # Å³ volume
g = 1       # degeneracy
W = 0.204868962802   # ev/(amu^(1/2)*Å)
# poly_order = 4 # order of polynomial for potenital fittings

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

function calc_anharm_wave_func(potential_matrix_1, potential_matrix_2, poly_order, Qi=-10, Qf=10, NQ=100, nev=10)
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
    x = linspace(Qi, Qf, NQ)

    # define potential using coefficients (calls polyfunc from Phonon)
    E1 = polyfunc(x, c1, poly_order)
    V1 = potential(x, E1)

    E2 = polyfunc(x, c2, poly_order)
    V2 = potential(x, E2)

    # calculate ΔQ and ΔE here
    Q1min = Q1[indmin(Energies_1)]
    Q2min = Q2[indmin(Energies_2)]
    ΔQ = Q2min - Q1min
    ΔE = minimum(Energies_2) - minimum(Energies_1)

    # solve Schrödinger equation
    # function solve1D_ev_amu(pot_ev_amu; NQ=1000, Qi=-10, Qf=10, nev=30, maxiter=10000)
    # Ground state

    ϵ1, χ1 = solve1D_ev_amu(x->polyfunc(x, c1, poly_order); NQ=NQ, Qi=Qi, Qf=Qf, nev=nev)
    # Excited state
    ϵ2, χ2 = solve1D_ev_amu(x->polyfunc(x, c2, poly_order); NQ=NQ, Qi=Qi, Qf=Qf, nev=nev)

    # Assign
    cc = CC()
    cc.V1 = V1; cc.V2 = V2
    cc.ϵ1 = ϵ1; cc.χ1 = χ1
    cc.ϵ2 = ϵ2; cc.χ2 = χ2
    cc
end

function plot_potentials(cc::CC)
    # close("all")
    fig = figure()
    ax = gca()
    scaling = 1E-2
    # plot potentials
    PyPlot.plot(cc.V1.Q, cc.V1.E, lw=1, color="red")
    PyPlot.plot(cc.V2.Q, cc.V2.E, lw=1, color="blue")

    # plot wavefunctions
    for i = 1:10
        PyPlot.plot(cc.V1.Q, cc.χ1[i]*scaling+cc.ϵ1[i], color="red", linewidth=0.5, alpha = 0.3)
        fill_between(cc.V1.Q, cc.ϵ1[i], cc.χ1[i]*scaling+cc.ϵ1[i], color="red", alpha=0.2)
        PyPlot.plot(cc.V2.Q, cc.χ2[i]*scaling+cc.ϵ2[i], color="blue", linewidth=0.5, alpha = 0.3)
        fill_between(cc.V2.Q, cc.ϵ2[i], cc.χ2[i]*scaling+cc.ϵ2[i], color="blue", alpha=0.2)
    end

    # # plot overlap
    # for i = 1:10
    #     PyPlot.plot(cc.V1.Q, cc.ϵ1[i]+integrand*1E-1, color="green", linewidth=0.5, alpha = 0.3)
    #     fill_between(cc.V1.Q, cc.ϵ1[i], cc.ϵ1[i]+integrand*1E-1, color="green", alpha=0.3)
    # end
    xlabel("Q (amu)"); ylabel("Energy (eV)");
    PyPlot.grid("off")
    ax[:grid](color="gray", linestyle=":", linewidth=0.5)
    tight_layout()
end

function calc_overlap!(cc::CC; cut_off=0.25, σ=0.025)
    p = plot_potentials(cc)
    ΔL = (maximum(cc.V1.Q) - minimum(cc.V1.Q))/length(cc.V1.Q)
    cc.ϵ_list = []
    cc.overlap_list = []
    cc.δ_list = []
    for i in range(1, length(cc.ϵ1))
        for j in range(1, length(cc.ϵ2))
            Δϵ = abs(cc.ϵ1[i] - cc.ϵ2[j])
            if  Δϵ < cut_off
                integrand = (cc.χ1[i] .* cc.V1.Q .* cc.χ2[j])
                overlap = sum(integrand)*ΔL

                append!(cc.ϵ_list, cc.ϵ1[i])
                append!(cc.overlap_list, overlap)

                append!(cc.δ_list, exp(-(Δϵ/σ)^2/2)/(σ*sqrt(2*π)))

                # plot
                alpha = (cut_off-Δϵ)/cut_off
                # plot!(cc.V1.Q, cc.ϵ1[i]+integrand*1E-1, color="#31a354", lw=2, alpha=0.5)
                PyPlot.plot(cc.V1.Q, cc.ϵ1[i]+integrand*1E-1, color="orange", linewidth=0.5, alpha=0.1)
                fill_between(cc.V1.Q, cc.ϵ1[i], cc.ϵ1[i]+integrand*1E-1, color="orange", linewidth=0.5, alpha=0.1)
                # #31a354 is green
            end
        end
    end
    return(p)
end

function calc_capt_coeff(W, V, T_range, cc::CC)
    capt_coeff = zeros(length(T_range))
    # TODO: CHECK convergence of the partition function Z(number of eigenvalue)
    # TODO:       convergence over σ
    Z = 0
    β = 1 ./ (kB .* T_range)
    for ϵ in cc.ϵ1
        Z+=exp.(-β*ϵ)
    end
    # println('Z', Z)
    for summand in zip(cc.ϵ_list, cc.overlap_list, cc.δ_list)
        ϵ, overlap, δ = summand
        occ = exp.(-β*ϵ) ./ Z
        capt_coeff += occ * overlap .* overlap * δ
    end
    capt_coeff = V*2*π/ħ*g*W^2 * capt_coeff
    return(capt_coeff)
end

end # module
