"""
Submodule providing helper functions for generating plots.
"""
module Plotter
using ..CarrierCapture: Potential, conf_coord

export plot_pots, plot_pot!, plot_ccs, plot_cc!

using Plots, LaTeXStrings


"""
    plot_pot!(pot; lplt_wf = false, plt = nothing, color = Nothing, label = "", scale_factor = 2e-2)

Plots `Potential` with its' wave functions if `lplt_wf` is `true`.
"""
function plot_pot!(pot::Potential; lplt_wf = false, plt = nothing, color = Nothing, label = "", scale_factor = 2e-2)
    if plt == nothing
        plt = plot()
    end
    label = if label == "" pot.name else label end
    color = if color == nothing "black" else color end

    # plot wave functions
    if lplt_wf
        ϵ = pot.ϵ; χ = pot.χ
        for i = 1:length(ϵ)
            plot!(plt, pot.Q, χ[i, :] * scale_factor .+ ϵ[i], lw = 0,
                  fillrange = [χ[i, :] * 0 .+ ϵ[i], χ[i, :] * scale_factor .+ ϵ[i]],
                  color = color, alpha = 0.5, label = "")
        end
    end

    # plot function
    plot!(plt, pot.Q, pot.E, lw = 4, color = color, label = label)

    # plot data
    if size(pot.QE_data)[1] > 1 
        scatter!(plt, pot.QE_data.Q, pot.QE_data.E, color = color, label = "")
    end
    xaxis!(L"\ Q (amu^{1\/2} Å) \ (^{}$$"); yaxis!(L"\ Energy  (eV) \ (^{}$$")

    return plt
end


"""
    plot_cc!(cc; plt = nothing, color = nothing, label = "")

Plots capture coefficient `cc.capt_coeff` as a function of `1000/T`.
"""
function plot_cc!(cc::conf_coord; plt = nothing, color = nothing, label = "")
    if plt == nothing
        plt = plot()
    end
    label = if label == "" cc.name else label end
    # color = if color==Nothing "black" else color end
    plot!(plt, 1000 ./ cc.temperature, cc.capt_coeff, lw = 4, label = label)
    xaxis!(L"\ 1000\/T (K^{-1}) \ (^{}$$"); yaxis!(L"C (cm^{3}\/s) \ (^{}$$", :log10)
    return plt
end



end
