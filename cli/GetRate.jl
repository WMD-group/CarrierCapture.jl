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
    "--wave", "-w"
        help = "wave function file serialized"
        default = "wave.jld"
        arg_type = String
    "--plot", "-p"
        help = "Plot potentials"
        action = :store_true
    "--verbose", "-v"
        help = "write verbose capture coefficient"
        action = :store_true
end

args = parse_args(ARGS, s)
input = YAML.load(open(args["input"]))
wave_path = args["wave"]
is_verbose = args["verbose"]

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

println("==============================")
println("Read wave functions from <$(wave_path)>")
file = open(wave_path, "r"); pots = deserialize(file); close(file)
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
if args["plot"] == true
    println("====Ploting capture rate====\n")
    plot_cfg = get(input, "plot", Dict())
    plot_ccs(ccs, plot_cfg)
end

# write capture coefficient into cvs files
for (i, cc) in enumerate(ccs)
    CSV.write("capt_coeff_$(i).csv", DataFrame([cc.temperature, cc.capt_coeff], [:T, :C]))
end

# write verbose ouput
if is_verbose
    for (i, cc) in enumerate(ccs)
        h5open("overlap_$(i).hdf5", "w") do file
            write(file, "overlap", collect(cc.overlap_matrix))
        end
        h5open("partial_capt_coeff_$(i).hdf5","w") do file
            write(file, "partial_capt_coeff", collect(cc.partial_capt_coeff))
        end
    end
end

end # module
