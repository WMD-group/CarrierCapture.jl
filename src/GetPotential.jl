#!/usr/bin/env julia

push!(LOAD_PATH,"../src/")
module GetPotential
using Potential
using Plots, LaTeXStrings
using ArgParse, YAML
using CSV, DataFrames # consider CSVFiles
# Native serialization in Julia works fine.
# using JLD2, FileIO
using Serialization


function plot_pot!(pot::potential; lplt_wf=false, plt=Nothing, color=Nothing, label="", scale_factor=2e-2)
    if plt == Nothing
        plt = plot()
    end
    label = if label=="" pot.name else label end
    color = if color==Nothing pot.color else color end
    
    # plot wave functions
    if lplt_wf
        ϵ = pot.ϵ; χ = pot.χ
        for i = 1:length(ϵ)
            plot!(plt, pot.Q, χ[i]*scale_factor .+ ϵ[i], lw=0,
                  fillrange=[χ[i]*0 .+ ϵ[i], χ[i]*scale_factor .+ ϵ[i]],
                  color=color, alpha=0.5, label="")    
        end
    end
    
    # plot function
    plot!(plt, pot.Q, pot.E, lw=4, color=color, label=label)
    # plot data
    scatter!(plt, pot.QE_data.Q, pot.QE_data.E, color=color, label="")

    xaxis!(L"\ Q (amu^{1\/2} Å) \ (^{}$$"); yaxis!(L"\ Energy  (eV) \ (^{}$$")
    return plt
end

s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "input file in a yaml format (default:pot_input.yaml)"
        default = "pot_input.yaml"
        arg_type = String
end

input = YAML.load(open(parse_args(ARGS, s)["input"]))
# input = YAML.load(open("pot_input.yaml"))

Qi, Qf, NQ = input["Qi"], input["Qf"], input["NQ"]
Q = range(Qi, stop=Qf, length=NQ)

pots = Dict{String, potential}()
plt = plot()
for potential in input["potentials"]
	pot_cfg = potential["potential"]
	# TODO: make CSV.read compatable with a format "-.30073981E+03"
	#       consider <CSVFiles> package
	QE_data = names!(CSV.read(pot_cfg["data"]; allowmissing=:none), [:Q, :E])
	pot = pot_from_dict(QE_data, pot_cfg)
	pots[pot.name] = pot
	fit_pot!(pot, Q)
	solve_pot!(pot)
	plot_pot!(pot, lplt_wf=true, plt=plt)
end

open("potential.jld", "w") do file
    serialize(file, pots)
end

# plot
plot_cfg = get(input, "plot", Nothing)
if plot_cfg == Nothing
    ylims!(0., 3)
    xlims!(Qi, Qf)
else
    ylims!(plot_cfg["Emin"], plot_cfg["Emax"])
    xlims!(plot_cfg["Qmin"], plot_cfg["Qmax"])
end
savefig(plt, "pot.pdf")

# write capture coefficient cvs
for (name, pot) in pots
    CSV.write("eigval_$(name).csv", DataFrame([pot.ϵ .- pot.E0, pot.ϵ], [:ϵ_E0, :ϵ]))
end

end # module
