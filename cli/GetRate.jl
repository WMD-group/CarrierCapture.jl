#!/usr/bin/env julia

module GetRate

using CarrierCapture

# read arguments and input file
s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "input file in a yaml format"
        default = "input.yaml"
        arg_type = String
    "--potential", "-p"
        help = "potential file serialized"
        default = "potential.jld"
        arg_type = String
    "--verbose", "-v"
        help = "write verbose capture coefficient"
        action = :store_true
end

input = YAML.load(open(parse_args(ARGS, s)["input"]))
pot_path = parse_args(ARGS, s)["potential"]
is_verbose = parse_args(ARGS, s)["verbose"]

println(raw"
      _____                _            _____            _
     / ____|              (_)          / ____|          | |
    | |     __ _ _ __ _ __ _  ___ _ __| |     __ _ _ __ | |_ _   _ _ __ ___
    | |    / _` | '__| '__| |/ _ \ '__| |    / _` | '_ \| __| | | | '__/ _ \
    | |___| (_| | |  | |  | |  __/ |  | |___| (_| | |_) | |_| |_| | | |  __/
     \_____\__,_|_|  |_|  |_|\___|_|   \_____\__,_| .__/ \__|\__,_|_|  \___|
                                                  | |
                                                  |_| v 0.1
      ____      _   ____       _
     / ___| ___| |_|  _ \ __ _| |_ ___
    | |  _ / _ \ __| |_) / _` | __/ _ \
    | |_| |  __/ |_|  _ < (_| | ||  __/
     \____|\___|\__|_| \_\__,_|\__\___|
")

# capture variables
input_capt = input["captures"]
V = input_capt["Volume"]
Tmin, Tmax, NT = input_capt["Tmin"], input_capt["Tmax"], input_capt["NT"]
temperature = range(Tmin, stop=Tmax, length=NT)

# JLD2 has some problem with type
# pots = jldopen(pot_path, "r")
println("==============================")
println("Read potential from <$(pot_path)>")
file = open(pot_path, "r"); pots = deserialize(file); close(file)
println(keys(pots))
println("==============================")

# TODO: raise no eigenvalue error!
#       after dryrun of GetPotential
# calculate capture coefficient
ccs = Array{conf_coord,1}()
for cc in input_capt["ccs"]
    cc_cfg = cc["cc"]
    pot_i = pots[cc_cfg["initial"]]
    pot_f = pots[cc_cfg["final"]]
    cc = cc_from_dict(pot_i, pot_f, cc_cfg)
    calc_overlap!(cc)
    calc_capt_coeff!(cc, V, temperature)
    append!(ccs, [cc])
end

# plot
plot_cfg = get(input, "plot", Dict())
plot_ccs(ccs, plot_cfg)

# write capture coefficient into cvs files
for (i, cc) in enumerate(ccs)
    CSV.write("capt_coeff_$(i).csv", DataFrame([cc.temperature, cc.capt_coeff], [:T, :C]))
end

# write verbose ouput
if is_verbose
    for (i, cc) in enumerate(ccs)
        h5open("partial_capt_coeff_$(i).hdf5","w") do file
            write(file, "partial_capt_coeff", collect(cc.partial_capt_coeff))
        end
    end
end

end # module
