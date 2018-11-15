module Potential
using Brooglie # Atomic unit
using Plots

amu = 931.4940954E6   # eV / c^2
ħc = 0.19732697E-6    # eV m

function solve1D_ev_amu(pot_ev_amu; NQ=100, Qi=-10, Qf=10, nev=30, maxiter=nev*NQ)
    factor = (1/amu) * (ħc*1E10)^2

    ϵ1, χ1 = solve1D(x->pot_ev_amu(x*factor^0.5);
                  N=NQ, a=Qi/factor^0.5, b=Qf/factor^0.5, m=1, nev=nev, maxiter=maxiter)
    return ϵ1, χ1/factor^0.25
end

function plot_potential(Q, E, ϵ1, χ1; plt=Nothing, color="#bd0026", label="", scale_factor=2e-2)
    if plt == Nothing
       plt = plot()
    end
    plot!(plt, Q, E, lw=4, color=color, label=label)
    for i = 1:length(ϵ1)
        plot!(plt, Q, χ1[i]*scale_factor .+ ϵ1[i],
              fillrange=[χ1[i]*0 .+ ϵ1[i], χ1[i]*scale_factor .+ ϵ1[i]],
<<<<<<< HEAD
              color=color, alpha=0.5, label="")    
=======
              c="#bd0026", alpha=0.5, label="")
>>>>>>> master
    end
end


# Set of potentials
function sqwell(x, width, depth; x0=0.)
    x = x .- x0
    pot_well = (x .< -width/2.) .| (x .> width/2.)
    return pot_well .* depth
end

function harmonic(x, ħω)
    # ev(amu)
    a = amu / 2. * (ħω/ħc/1E10)^2
    return a*x.*x
end

function double(x, ħω1, ħω2)
    a1 = amu / 2. * (ħω1/ħc/1E10)^2
    a2 = amu / 2. * (ħω2/ħc/1E10)^2
    return - a1*x.*x + a2*x.*x.*x.*x
end

function polyfunc(x, coeffs, poly_order)
    y = 0 .* x
    for i = 1:poly_order + 1
        y += coeffs[i].*x.^(i-1)
    end
    return y
end

end # module
