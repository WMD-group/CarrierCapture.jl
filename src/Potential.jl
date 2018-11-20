__precompile__()

push!(LOAD_PATH,"../src/")

module Potential
using Brooglie # Atomic unit
# using Plots
using DataFrames
using LsqFit
# using JLD2

amu = 931.4940954E6   # eV / c^2
ħc = 0.19732697E-6    # eV m

export potential, pot_from_dict, fit_pot!, solve_pot!
export solve1D_ev_amu
export sqwell, harmonic, double, polyfunc

mutable struct potential
    name::String
    color::String
    QE_data::DataFrame
    E0::Float64
    func_type::String
    func::Function
    params::Dict{String, Number}
    p0::Array{Float64,1}
    Q::Array{Float64,1}; E::Array{Float64,1}
    nev::Int
    # ϵ includes E0
    ϵ::Array{Float64,1}; χ::Array{Array{Float64,1},1}
    # TODO: JLD2
    #       Don't blame S. Kim.
    #       Blame JLD2
    potential() = new("", "black", DataFrame([0. 0.]), Inf, 
                      "polyfunc", x->0, Dict(), [0.],  
                      [], [], 
                      0, [], [[]])
end


function pot_from_dict(QE_data::DataFrame, cfg::Dict)::potential
    pot = potential()
    pot.name = cfg["name"]
    pot.color = cfg["color"]
    pot.nev = cfg["nev"]
    pot.func_type = cfg["function"]["type"]
    pot.p0 =  parse.(Float64, split(cfg["function"]["p0"]))
    pot.params = convert(Dict{String, Number}, cfg["function"]["params"])
    pot.E0 = get(cfg, "E0", Inf)
    
    pot.QE_data = QE_data 
    if pot.E0 < Inf
        pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
    end

    return pot
end


function fit_pot!(pot::potential, Q)
    # find func
    func = @eval $(Symbol(pot.func_type))
    params = pot.params
    fit = curve_fit((x,p) -> func.(x, Ref(p); param=params), pot.QE_data.Q, pot.QE_data.E, pot.p0)
    pot.Q = Q
    pot.E = func.(Q, Ref(fit.param); param=params)
    pot.func = x -> func.(x, Ref(fit.param); param=params)
end


function solve_pot!(pot::potential)
    # solve
    pot.ϵ, pot.χ = solve1D_ev_amu(pot.func; 
        NQ=length(pot.Q), Qi=minimum(pot.Q), Qf=maximum(pot.Q), nev=pot.nev)
end


function solve1D_ev_amu(pot_ev_amu; NQ=100, Qi=-10, Qf=10, nev=30, maxiter=nev*NQ)
    factor = (1/amu) * (ħc*1E10)^2

    ϵ1, χ1 = solve1D(x->pot_ev_amu(x*factor^0.5);
                  N=NQ, a=Qi/factor^0.5, b=Qf/factor^0.5, m=1, nev=nev, maxiter=maxiter)
    return ϵ1, χ1/factor^0.25
end


# function plot_potential(Q, E, ϵ1, χ1; plt=Nothing, color="#bd0026", label="", scale_factor=2e-2)
#     if plt == Nothing
#        plt = plot()
#     end
#     plot!(plt, Q, E, lw=4, color=color, label=label)
#     for i = 1:length(ϵ1)
#         plot!(plt, Q, χ1[i]*scale_factor .+ ϵ1[i],
#               fillrange=[χ1[i]*0 .+ ϵ1[i], χ1[i]*scale_factor .+ ϵ1[i]],
#               color=color, alpha=0.5, label="")    
#     end
# end


# Set of potentials
function sqwell(x, width, depth; x0=0.)
    x = x .- x0
    pot_well = (x .< -width/2.) .| (x .> width/2.)
    return pot_well .* depth
end

function harmonic(x, ħω)
    # ev(amu)
    a = amu / 2. * (ħω/ħc/1E10)^2
    return a*x.*x
end

function double(x, ħω1, ħω2)
    a1 = amu / 2. * (ħω1/ħc/1E10)^2
    a2 = amu / 2. * (ħω2/ħc/1E10)^2
    return - a1*x.*x + a2*x.*x.*x.*x
end

function polyfunc(x, coeffs; param)
    poly_order = Int(param["poly_order"])
    y = 0 .* x
    for i = 1:poly_order + 1
        y += coeffs[i].*x.^(i-1)
    end
    return y
end

end # module
