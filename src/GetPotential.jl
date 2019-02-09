#!/usr/bin/env julia

push!(LOAD_PATH,"../src/")
module GetPotential
using Potential
using ArgParse, YAML
using CSV, DataFrames # consider CSVFiles
# Native serialization in Julia works fine.
# using JLD2, FileIO
using Serialization
using Printf

# read arguments and input file
s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "Input file in a yaml format (default:input.yaml)"
        default = "input.yaml"
        arg_type = String
    "--dryrun"
    	help = "Turn off the Srödinger equation solver"
    	action = :store_true
    "--plot", "-p"
    	help = "Plot potentials"
    	action = :store_true
end
args = parse_args(ARGS, s)

input = YAML.load(open(args["input"]))

# Set up global variables
Qi, Qf, NQ = input["Qi"], input["Qf"], input["NQ"]
Q = range(Qi, stop=Qf, length=NQ)

# solve potentials
pots = Dict{String, potential}()
for potential in input["potentials"]
	pot_cfg = potential["potential"]
	# TODO: make CSV.read compatable with a format "-.30073981E+03"
	#       consider <CSVFiles> package
	QE_data = names!(CSV.read(pot_cfg["data"]; allowmissing=:none), [:Q, :E])
	pot = pot_from_dict(QE_data, pot_cfg)
	pots[pot.name] = pot
	fit_pot!(pot, Q)
	if args["dryrun"] == false
		solve_pot!(pot)
	end
end

# print zero-phonon frequencies
if args["dryrun"] == false
	println("===========ħω0===========")
	for (name, pot) in pots
		ϵ0 = pot.ϵ[2] - pot.ϵ[1]
		ħω0 = @sprintf("%.1f", ϵ0*1000)
		println("$(name): $(ħω0) meV")
	end
	println("=========================")
end

# save potential structs
open("potential.jld", "w") do file
    serialize(file, pots)
end

# plot
if args["plot"] == true
	using Plotter
	plot_cfg = get(input, "plot", Nothing)
	plot_pots(pots, plot_cfg)
end

# write potential cvs
for (name, pot) in pots
    CSV.write("pot_fit_$(name).csv", DataFrame([pot.Q, pot.E], [:Q, :E]))
end

# write capture coefficient cvs
for (name, pot) in pots
    CSV.write("eigval_$(name).csv", DataFrame([pot.ϵ .- pot.E0, pot.ϵ], [:ϵ_E0, :ϵ]))
end

end # module
