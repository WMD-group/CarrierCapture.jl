push!(LOAD_PATH,"../src/")
module GetRate
using CaptureRate
using Plots
using ArgParse
using YAML
# using JLD2, FileIO
using Serialization
using Potential
using LsqFit
using LaTeXStrings

function plot_cc!(cc::conf_coord; plt=Nothing, color=Nothing, label="")
    if plt == Nothing
        plt = plot()
    end
    # label = if label=="" pot.name else label end
    # color = if color==Nothing pot.color else color end
    # plot!(plt, 1000 ./ cc.temp, cc.capt_coeff, lw=4, color=color, label=label)
    plot!(plt, 1000 ./ cc.temp, cc.capt_coeff, lw=4)
    xlims!(0, 5)
    ylims!(1E-10, 1E-5)
    xaxis!(L"\ 1000\/T (K^{-1}) \ (^{}$$")
    yaxis!(L"C_{n}  (cm^{3}\/s) \ (^{}$$", :log10)
end

s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "input file in a yaml format"
        default = "pot_input.yaml"
        arg_type = String
    "--potential", "-p"
        help = "potential file serialized"
        default = "potential.jld"
        arg_type = String
end
input = YAML.load(open(parse_args(ARGS, s)["input"]))
pot_path = parse_args(ARGS, s)["potential"]

# Global variables
# TODO: ? assert Q == Q in pots
# Qi, Qf, NQ = input["Qi"], input["Qf"], input["NQ"]
# Q = range(Qi, stop=Qf, length=NQ)

# capture variables
input_capt = input["captures"]
V = input_capt["Volume"]
Tmin, Tmax, NT = input_capt["Tmin"], input_capt["Tmax"], input_capt["NT"]
temp = range(Tmin, stop=Tmax, length=NT)

# JLD2 has some problem with type
# pots = jldopen(pot_path, "r")
println("==============================")
println("Read potential from $pot_path")
file = open("potential.jld", "r")
pots = deserialize(file)
close(file)
println(keys(pots))
println("==============================")


ccs = []
for cc in input_capt["ccs"]
    cc_cfg = cc["cc"]
    pot_i = pots[cc_cfg["initial"]]
    pot_f = pots[cc_cfg["final"]]
    # cc = conf_coord(pot_i, pot_f)
    cc = cc_from_dict(pot_i, pot_f, cc_cfg)
    calc_overlap!(cc)
    calc_capt_coeff!(cc, V, temp)
    append!(ccs, [cc])
end
println(length(ccs))

plt = plot()
for cc in ccs
    plot_cc!(cc; plt=plt)
end
savefig(plt, "capt.pdf")


end # module
