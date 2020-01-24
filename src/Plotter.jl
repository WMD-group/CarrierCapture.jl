"""
submodule providing helper functions for generating plots.
"""
module Plotter
using ..CarrierCapture: potential, conf_coord

export plot_pots, plot_pot!, plot_ccs, plot_cc!

using Plots, LaTeXStrings

function plot_pots(pots::Dict{String,potential}, plot_cfg; output_fig = "potential.pdf")
    plt = plot(legend = :bottomleft)

    for (name, pot) in pots
        plot_pot!(pot, lplt_wf = true, plt = plt)
    end

    Emin = get(plot_cfg, "Emin", minimum(hcat([pot.E for (name, pot) in pots]...)))
    Emax = get(plot_cfg, "Emax", maximum(hcat([pot.E for (name, pot) in pots]...)))
    ylims!(Emin, Emax)

    Qmin = get(plot_cfg, "Qmin", minimum(hcat([pot.Q for (name, pot) in pots]...)))
    Qmax = get(plot_cfg, "Qmax", maximum(hcat([pot.Q for (name, pot) in pots]...)))
    xlims!(Qmin, Qmax)

    savefig(plt, output_fig)
end

function plot_pot!(pot::potential; lplt_wf = false, plt = Nothing, color = Nothing, label = "", scale_factor = 2e-2)
    if plt == Nothing
        plt = plot()
    end
    label = if label == "" pot.name else label end
    color = if color == Nothing "black" else color end

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


function plot_ccs(ccs::Array{conf_coord,1}, plot_cfg::Dict; output_fig = "captcoeff.pdf")
    plt = plot(legend = :bottomleft)
    for cc in ccs
        plot_cc!(cc; plt = plt)
    end
    
    Cmin = get(plot_cfg, "Cmin", minimum(hcat([cc.capt_coeff for cc in ccs]...)) / 2.)
    Cmax = get(plot_cfg, "Cmax", maximum(hcat([cc.capt_coeff for cc in ccs]...)) * 2.)
    ylims!(Cmin, Cmax)

    invTmin = get(plot_cfg, "invTmin", 1000 / maximum(hcat([cc.temperature for cc in ccs]...)))
    invTmax = get(plot_cfg, "invTmax", 1000 / minimum(hcat([cc.temperature for cc in ccs]...)))
    xlims!(invTmin, invTmax)

    savefig(plt, output_fig)
end

function plot_cc!(cc; plt = Nothing, color = Nothing, label = "")
    if plt == Nothing
        plt = plot()
    end
    label = if label == "" cc.name else label end
    # color = if color==Nothing "black" else color end
    plot!(plt, 1000 ./ cc.temperature, cc.capt_coeff, lw = 4, label = label)
    xaxis!(L"\ 1000\/T (K^{-1}) \ (^{}$$"); yaxis!(L"C (cm^{3}\/s) \ (^{}$$", :log10)
    return plt
end

end