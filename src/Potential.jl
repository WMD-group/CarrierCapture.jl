__precompile__()

push!(LOAD_PATH,"../src/")

module Potential
using Brooglie # Atomic unit
# using Plots
using DataFrames
using LsqFit
# using JLD2
using Polynomials

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
                      "func_type", x->0, Dict(), [0.],  
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
    pot.E0 = get(cfg, "E0", Inf)
    pot.QE_data = QE_data 

    pot.params = get(cfg["function"], "params",  Dict("E0"=>pot.E0))
    pot.params["E0"] = pot.E0
    pot.params["Q0"] = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    if pot.E0 < Inf
        pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
    end

    return pot
end


function fit_pot!(pot::potential, Q)
    """
    """
    pot.params["E0"] = pot.E0
    pot.params["Q0"] = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    E_CUT = 3
    e_cut_ind = pot.QE_data.E .< E_CUT

    func = @eval $(Symbol(pot.func_type))
    params = pot.params
    fit = curve_fit((x,p) -> func.(x, Ref(p); param=params), pot.QE_data.Q[e_cut_ind], pot.QE_data.E[e_cut_ind], pot.p0)
    pot.Q = Q
    pot.E = func.(Q, Ref(fit.param); param=params)
    pot.func = x -> func.(x, Ref(fit.param); param=params)

    println("===========Fit===========")
    println("Function: $(pot.func_type)")
    println("Best fit: $(fit.param)")
    println("=========================")
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


# Set of potentials
function sqwell(x, width, depth; x0=0.)
    x = x .- x0
    pot_well = (x .< -width/2.) .| (x .> width/2.)
    return pot_well .* depth
end

function harmonic(x, ħω; param)
    # ev(amu)
    a = amu / 2. * (ħω/ħc/1E10)^2
    return a*x.*x
end

function double(x, ħω1, ħω2; param)
    a1 = amu / 2. * (ħω1/ħc/1E10)^2
    a2 = amu / 2. * (ħω2/ħc/1E10)^2
    return - a1*x.*x + a2*x.*x.*x.*x
end

function polyfunc(x, coeffs; param)
    E0 = param["E0"]
    Q0 = param["Q0"]
    poly_order = Int(param["poly_order"])
    y = 0 .* x
    y += E0
    for i = 2:poly_order + 1
        y += coeffs[i].* (x - Q0) .^(i-1)
    end
    return y
end

function morse(x, coeffs; param)
    E0 = param["E0"]
    Q0 = param["Q0"]
    return coeffs[1].*((1-exp.(-coeffs[2].*(x-Q0))).^2)+E0
end

function morse_quatic(x, coeffs; param)
    E0 = param["E0"]
    Q0 = param["Q0"]
    m = Int(param["poly_order"])
    
    A = coeffs[1]^2
    a = coeffs[2]
    r0 = coeffs[3]
    r1 = coeffs[4]
    B = coeffs[5]^2

    morse = A*(1-exp.(-a.*(x-Q0-r0))).^2 .- A*(1-exp.(a*r0)).^2
    poly = B*(x-Q0-r1).^m - B*(-r1)^m
    derv = 2*A*a*exp(a*r0)*(1-exp(a*r0)) + m*B*(-r1)^(m-1)
    linear = - derv.*x + derv*Q0
    return morse + poly + linear + E0
end

function morse_poly(x, coeffs; param)
    E0 = param["E0"]
    Q0 = param["Q0"]
    m = Int(param["poly_order"])
    A = coeffs[1]
    a = coeffs[2]

    r0 = coeffs[3]
    r1 = coeffs[4]

    morse = x -> A*(exp.(-2*a.*(x-Q0-r0)) - 2*exp.(-a.*(x-Q0-r0)))
    
    X = Poly([-Q0-r1, 1])
    poly = coeffs[3+m]^2*X^m
    slope = A*(-2a)*(exp(2*a*r0) - exp(a*r0)) + m*coeffs[3+m]*(-r1)^(m-1)
    for i = 2:m-1
        poly += coeffs[3+i]*X^i
        slope += i*coeffs[3+i]*(-r1)^(i-1)
    end
    linear = -slope.*x .+ slope*Q0
    return morse(x) + poly(x) + linear .+ E0 .- morse(Q0) .- poly(Q0)
end

# function morse_harmonic(x, coeffs; param)
#     E0 = param["E0"]
#     Q0 = param["Q0"]
#     if coeffs[2]*(x - Q0) > 0
#         return coeffs[1].*((1-exp.(-coeffs[2].*(x-Q0))).^2)+E0
#     else
#         return coeffs[5].*(x-Q0)^2+E0
#     end
# end



end # module
