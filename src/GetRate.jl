push!(LOAD_PATH,"../src/")
module GetRate
using Potential
using CaptureRate
using Plots
using ArgParse
using YAML
using CSV
using JLD2, FileIO


s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "input file in a yaml format (default: pot_input.yaml)"
        default = "pot_input.yaml"
        arg_type = String
    "--potential", "-p"
        help = "potential file in JLD2 format (default: potential.jld2"
        default = "potential.jld2"
        arg_type = String
end

input = YAML.load(open(parse_args(ARGS, s)["input"]))
pot_path = open(parse_args(ARGS, s)["potential"])

Qi, Qf, NQ = input["Qi"], input["Qf"], input["NQ"]
Q = range(Qi, stop=Qf, length=NQ)

f = jldopen(pot_path, "r")
println(f)


end # module
