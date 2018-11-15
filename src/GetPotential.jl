push!(LOAD_PATH,"../src/")
module GetPotential
using Potential
using Plots
using ArgParse
using YAML
using CSV


function plot_pot!(pot::potential; lplt_wf=false, plt=Nothing, 
			       color="#bd0026", label="", scale_factor=2e-2)
    if plt == Nothing
        plt = plot()
    end
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

# Global
Qi, Qf, NQ = input["Qi"], input["Qf"], input["NQ"]
Q = range(Qi, stop=Qf, length=NQ)

# potential_i
pot_i_cfg = input["potential_i"]
# import data 
QE_data = convert(Array{Float64, 2}, CSV.read(pot_i_cfg["data"]))
pot_i = pot_from_dict(QE_data, pot_i_cfg)
fit_pot!(pot_i, Q)
solve_pot!(pot_i)

# potential_f
pot_f_cfg = input["potential_f"]
# import data 
QE_data = convert(Array{Float64, 2}, CSV.read(pot_f_cfg["data"]))
pot_f = pot_from_dict(QE_data, pot_f_cfg)
fit_pot!(pot_f, Q)
solve_pot!(pot_f)


# plot
plt = plot_pot!(pot_i, lplt_wf=true)
plot_pot!(pot_f, lplt_wf=true, color="#005dff", plt=plt)

if get(input, "Egap", Nothing) != Nothing
	# potential_i + Egap
	pot_g_cfg = input["potential_i"]
	# import data 
	QE_data = convert(Array{Float64, 2}, CSV.read(pot_g_cfg["data"]))
	QE_data[:, 2] .+= input["Egap"]
	pot_g = pot_from_dict(QE_data, pot_g_cfg)
	fit_pot!(pot_g, Q)
	solve_pot!(pot_g)
	plot_pot!(pot_g, lplt_wf=true, plt=plt)
end

ylims!(0., 3)
xlims!(Qi, Qf)

# savefig(plt, "pot.png")
savefig(plt, "pot.pdf")
end # module
