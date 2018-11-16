push!(LOAD_PATH,"../src/")
module GetPotential
using Potential
using Plots
using ArgParse
using YAML
using CSV
using JLD2, FileIO


function plot_pot!(pot::potential; lplt_wf=false, plt=Nothing, 
			       color=Nothing, label="", scale_factor=2e-2)
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
    # plot
    scatter!(plt, pot.QE_data[:,1], pot.QE_data[:,2], color=color, label="")

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

pots = []
plt = plot()
for potential in input["potentials"]
	pot_cfg = potential["potential"]
	QE_data = convert(Array{Float64, 2}, CSV.read(pot_cfg["data"]))
	pot = pot_from_dict(QE_data, pot_cfg)
	append!(pots, [pot])
	fit_pot!(pot, Q)
	solve_pot!(pot)
	plot_pot!(pot, lplt_wf=true, plt=plt)
end

# Dump data for Rate
jldopen("potential.jld2", "w") do file
    # addrequire(file, Potential)
    for pot in pots
    	file[pot.name] = pot
    end
end

plot_cfg = get(input, "plot", Nothing)
if plot_cfg == Nothing
    ylims!(0., 3)
    xlims!(Qi, Qf)
else
    ylims!(plot_cfg["Emin"], plot_cfg["Emax"])
    xlims!(plot_cfg["Qmin"], plot_cfg["Qmax"])
end
savefig(plt, "pot.pdf")

end # module
