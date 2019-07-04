# defining constants
amu = 931.4940954E6   # eV / c^2
ħc = 0.19732697E-6    # eV m

export potential, pot_from_dict, fit_pot!, solve_pot!, find_crossing
export solve1D_ev_amu
export sqwell, harmonic, double, polyfunc
export get_bspline, get_spline

# set up potential structure
mutable struct potential
    name::String
    color::String
    QE_data::DataFrame
    E0::Float64; Q0::Float64
    func_type::String
    func
    params::Dict{String, Any}
    p0::Array{Float64,1}
    Q::Array{Float64,1}; E::Array{Float64,1}
    nev::Int
    # ϵ includes E0
    ϵ::Array{Float64,1}; χ::Array{Float64,2}
    # TODO: JLD2 doesn't work
    #       Don't blame S. Kim.
    #       Blame JLD2
    potential() = new("", "black", DataFrame(), Inf, 0,
                      "func_type", x->0, Dict(), [0],
                      [], [],
                      0, [], Array{Float64}(undef, 0, 2))
end

# read potential
function pot_from_dict(QE_data::DataFrame, cfg::Dict)::potential
    pot = potential()
    pot.name = cfg["name"]
    pot.color = cfg["color"]
    pot.nev = cfg["nev"]
    pot.func_type = cfg["function"]["type"]
    pot.p0 =  parse.(Float64, split(get(cfg["function"], "p0", "0")))
    pot.E0 = get(cfg, "E0", Inf)
    pot.QE_data = QE_data

    pot.params = get(cfg["function"], "params",  Dict("E0" => pot.E0))

    pot.Q0 = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    if pot.E0 < Inf
        pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
    end

    return pot
end

# fitting potential
function fit_pot!(pot::potential, Q)
    """
    """
    # pot.params["E0"] = pot.E0
    # pot.params["Q0"] = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    E_CUT = 2 # defines upper energy limit on data used for fitting
    e_cut_ind = pot.QE_data.E .< E_CUT + pot.E0

    params = pot.params
    pot.Q = Q

    println("Potential fitting: $(pot.name)")

    if pot.func_type == "bspline"
        println("=========bspline=========\n")
        bspline = get_bspline(pot.QE_data.Q, pot.QE_data.E)
        pot.E = bspline(Q)
        pot.func = bspline
    elseif pot.func_type == "spline"
        println("=========spline==========\n")
        weight = get(params, "weight", nothing)
        weight = parse.(Float64, split(weight))
        smoothness = get(params, "smoothness", 0)
        order = get(params, "order", 2)
        spline = get_spline(pot.QE_data.Q, pot.QE_data.E; weight = weight, smoothness = smoothness, order = order)
        pot.E = spline(Q)
        pot.func = spline
    elseif pot.func_type == "harmonic"
        println("========harmonic=========\n")
        # func = x -> harmonic(x, params["hw"]; param = params)
        pot.E = harmonic.(Q, params["hw"]; E₀ = pot.E0, Q₀ = pot.Q0)
        pot.func = x -> harmonic(x, params["hw"]; E₀ = pot.E0, Q₀ = pot.Q0)
    else
        if pot.func_type == "polyfunc"
            func = polyfunc
        elseif pot.func_type == "morse_poly"
            func = (x, p) -> morse_poly(x, Ref(p); E₀ = pot.E0, Q₀ = pot.Q0, poly_order = params["poly_order"])
        else
            func = (x, p) -> polyfunc(x, Ref(p); E₀ = pot.E0, Q₀ = pot.Q0, poly_order = params["poly_order"])
        end
        fit = curve_fit((x, p) -> func.(x, p), pot.QE_data.Q[e_cut_ind], pot.QE_data.E[e_cut_ind], pot.p0)
        pot.E = func.(Q, Ref(fit.param); param = params)
        pot.func = x -> func(x, fit.param; param = params)

        println("===========Fit===========")
        println("Function: $(pot.func_type)")
        println("Best fit: $(fit.param)")
        println("=========================\n")
    end
end


function solve_pot!(pot::potential)
    pot.ϵ, pot.χ = solve1D_ev_amu(pot.func;
        NQ=length(pot.Q), Qi=minimum(pot.Q), Qf=maximum(pot.Q), nev=pot.nev)
end


function solve1D_ev_amu(func; NQ=100, Qi=-10, Qf=10, nev=30, maxiter=nev*NQ)
    factor = (1/amu) * (ħc*1E10)^2

    ϵ1, χ1 = Brooglie.solve(x -> func.(x*factor^0.5);
                  N=NQ, a=Qi/factor^0.5, b=Qf/factor^0.5, m=1, nev=nev, maxiter=maxiter)
    return ϵ1, vcat(χ1'...)/factor^0.25
end

function find_crossing(pot_1::potential, pot_2::potential)
    # find root of pot_1 - pot_2 = 0
    diff_func = x-> pot_1.func(x) - pot_2.func(x)
    Q = pot_1.Q
    roots = find_zero(diff_func, Q[length(Q)÷2])
    return roots, pot_1.func.(roots)
end


############################################################
# Set of potentials
function sqwell(x, width, depth; x0=0.)
    x = x .- x0
    pot_well = (x .< -width/2.) .| (x .> width/2.)
    return pot_well .* depth
end

function harmonic(x, ħω; E₀, Q₀)
    a = amu/2 * (ħω/ħc/1E10)^2
    return a*(x - Q₀)^2 + E₀
end

function double(x, ħω1, ħω2; param)
    a1 = amu / 2. * (ħω1/ħc/1E10)^2
    a2 = amu / 2. * (ħω2/ħc/1E10)^2
    return - a1*x.*x + a2*x.*x.*x.*x
end

# Set of potentials for fitting
function polyfunc(x, coeffs; E₀, Q₀, poly_order)
    # E0 = param["E0"]
    # Q₀ = param["Q0"]
    # poly_order = Int(param["poly_order"])
    y = 0 .* x
    y += E₀
    for i = 2:poly_order + 1
        y += coeffs[i].* (x - Q₀) .^(i-1)
    end
    return y
end

function morse(x, coeffs; E₀, Q₀)
    return coeffs[1].*((1-exp.(-coeffs[2].*(x-Q₀))).^2)+E₀
end

function morse_poly(x, coeffs; E₀, Q₀, poly_order)
    # TODO: CHECK # of coeffs and param["poly_order"], and give a WARNING
    # E₀ = param["E0"]
    # Q₀ = param["Q0"]
    # po = param["poly_order"]

    orders = if (typeof(po) == Int64) Array([po]) else parse.(Int, split(po)) end

    A = abs(coeffs[1])
    a = coeffs[2]

    r₀ = coeffs[3]

    morse = x -> A*(exp(-2*a*(x-Q₀-r₀)) - 2*exp(-a*(x-Q₀-r₀)))

    coeff = zeros(maximum(orders)+1)
    for (index, poly_order) in enumerate(orders)
        coeff[poly_order+1] = coeffs[3+index]
    end

    coeff[end] = abs(coeff[end])

    poly = Poly(coeff)
    tot_derv = polyder(poly, 1) -  A*(-2a)*(exp(2a*r₀) - exp(a*r₀))

    r₁ = Polynomials.roots(tot_derv)
    r₁ = real.(r₁[imag.(r₁) .== 0])
    r₁ = r₁[argmin(abs.(r₁))]

    return morse(x) + poly(x-Q₀-r₁) + E₀ - morse(Q₀) - poly(-r₁)
end

function get_spline(Qs, Es; weight = nothing, smoothness = 0, order = 2)
    if Qs[1] > Qs[end]
        Qs = reverse(Qs)
        Es = reverse(Es)
    end

    if weight == nothing
        spl = Spline1D(Qs, Es, k = order, bc = "extrapolate", s = smoothness)
    else
        spl = Spline1D(Qs, Es, w = weight, k = order, bc = "extrapolate", s = smoothness)
    end
    return spl
end

function get_bspline(Qs, Es)
    # BSplines assume your data is uniformly spaced on the grid
    # Qs, Es have to be eqaully-spaced (Range)
    if Qs[1] > Qs[end]
        Qs = reverse(Qs)
        Es = reverse(Es)
    end
    Qs_range = range(Qs[1], stop = Qs[end], length = length(Qs))

    itp = interpolate(Es, BSpline(Quadratic(Line(OnCell()))))
    etpf = extrapolate(itp, Line())
    setpf = scale(etpf, Qs_range)

    return setpf
end
