# defining constants
amu = 931.4940954E6   # eV / c^2
ħc = 0.19732697E-6    # eV m
boltz = 8.617333262E-5 # eV/K

export Potential, pot_from_dict, filter_sample_points!, fit_pot!, solve_pot!, find_crossing, potential_from_file, cleave_pot
export Plotter
# export solve1D_ev_amu
# export sqwell, harmonic, double_well, polyfunc, morse
# export get_bspline, get_spline

"""
Stores a potential in one-dimensional space Q, with discreet points (E0, Q0) and fitting function func.

## Fields

- `name` -- the name of Potential.
- `QE_data`   -- the (n X 2) DataFrame of data points (Q vs Energy).
- `E0`, `Q0`  -- the minimum point of the Potential [`Q0`, `E0`].
- `func_type`     -- the type of fitting function ("bspline", "spline", "harmonic", "polyfunc", "morse_poly", "morse").
- `params`    -- the list of hyper parameters for the fitting function.
- `Q`, `E`  -- `Q` and `E`=`func(Q, p_opt; params)`. They differ from QE_data
in that they are data points of the fitted function.
- `nev`  -- the number of eigenvalues to be evaluated.
- `ϵ` -- the list of eigenvalues
- `T` -- temperature (only used for filtering sample points)

## Constructor
    Potential()

"""
mutable struct Potential
    name::String
    QE_data::DataFrame
    E0::Float64; Q0::Float64
    func_type::String
    func # can be either Function or Interpolations.ScaledInterpolation
    params::Dict{String, Any}
    Q::Array{Float64,1}; E::Array{Float64,1}
    nev::Int64
    # ϵ includes E0
    ϵ::Array{Float64,1}; χ::Array{Float64,2}
    # options for thermally accessible fit
    T::Float64
    Potential() = new("", DataFrame([0 0], [:Q, :E]), 0, 0,
                      "harmonic_fittable", x->0, Dict(),
                      [], [],
                      0, [], Array{Float64}(undef, 0, 2), 293.15)
end

# extend copy to work with potentials, a cleaner implementation would be welcome
function Base.copy(pot::Potential)
    copy_pot = Potential()
    copy_pot.name = pot.name
    copy_pot.QE_data = pot.QE_data
    copy_pot.E0 = pot.E0
    copy_pot.Q0 = pot.Q0
    copy_pot.func_type = pot.func_type
    copy_pot.func = pot.func
    copy_pot.params = pot.params
    copy_pot.Q = pot.Q
    copy_pot.E = pot.E
    copy_pot.nev = pot.nev
    copy_pot.ϵ = pot.ϵ
    copy_pot.χ = pot.χ
    copy_pot.T = pot.T
    return copy_pot
end

"""
    potential_from_file(filename::String, resolution::Int64 = 3001)

Parse a two column file with data for reaction coordinate and energy.
Lines beginning with a hash are ignored. The optional argument `resolution`
determines how many points are included in the interpolation between the
read-in data points, pot.Q.

Additionally, any pot.QE_data.Q not included in pot.Q are added.

"""
function potential_from_file(filename::String, resolution::Int64 = 3001)
    Q_dat = Float64[]
    E_dat = Float64[]
    # read each line
    for line in eachline(filename)
        # skip comments
        if line[begin] != '#'
            line_split=split(line)
            append!(Q_dat,parse(Float64,line_split[1]))
            append!(E_dat,parse(Float64,line_split[2]))
        end
    end
    pot = Potential()
    pot.QE_data = DataFrame(Q = Q_dat[:], E = E_dat[:])
    # range of points
    interpol_points = range(Q_dat[1], stop=Q_dat[end], length=resolution)
    # add sample points from QE_data, sort and remove duplicates
    pot.Q = unique(sort(vcat(interpol_points, Q_dat[:])))

    return pot
end


"""
    filter_sample_points!(pot::Potential)

Remove all fit points after a thermally unsurmountable point in increasing order.
This function is only robust if

"""
function filter_sample_points!(pot::Potential, thermal_energy::Float64)
    thermal_energy = pot.T*boltz

    # find minima
    min_ind_array = findall(pot.QE_data.E .== minimum(pot.QE_data.E))
    len_min_ind = length(min_ind_array)
    # check if there are multiple maxima
    if len_min_ind > 1
        throw(ArgumentError("The fitted potential has several minima."))
    end

    min_ind = len_min_ind[1]

    # select the interconnected data points below the energy threshold and
    # containing the minimum energy point we call these points an island
    below_thresh = pot.QE_data.E .<= thermal_energy + pot.E0

    island = []

    for (ind, bool) in enumerate(below_thresh)
        # if the data point is not part of an island
        if bool == 0
            # if there is a non-0 island in memory, add it and clear memory
            if length(island) > 0
                if min_ind in island
                    break
                else
                island = []
                end
            end
        # if the point is part of an island
        else
            # add point to island
            append!(island, ind)
        end
    end

    pot.QE_data = DataFrame(Q = pot.QE_data.Q[island], E = pot.QE_data.E[island])
    println(island)

end

"""
    function filter_sample_points!(pot::Potential, thermal_energy::Float64)

Same as above but use the temperature associated with the potential to filter points

"""
function filter_sample_points!(pot::Potential)
    filter_sample_points!(pot, pot.T*boltz)
end

"""
    fit_pot!(pot::Potential)

Fit a function `pot.func_type` to `pot.QE_data` on the domain `pot.Q`.

## parameters

- `pot`: `Potential`
    - `pot.func_type`: the fitting function; `{"spline" (preferred), "bspline", "harmonic", "polyfunc", "morse_poly", "morse"}`.
    - `pot.params`: the hyperparameters.

### Hyperparameters `params`

- spline (preferred)

    Spline interpolation. See more detail in [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl).
    - `order`: spline order (between 1 and 5; default 2).
    - `smoothness`: the amount of smoothness is determined by the condition that `sum((w[i]*(y[i]-spline(x[i])))**2) <= s`
    - `weight`: the weight applied to each `QE_data` point (length m 1-d array).


- bspline

    Basic spline interpolation. See more detail in [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).


- harmonic

    Harmonic function whose minimum is at [`pot.Q0`, `pot.E0`].
    - `hw`: the energy quanta of the harmonic oscillator.


- polyfunc

    Polynomial function;
        `y = E₀ + Σ coeffs[i].* (x .- Q₀) .^(i-1)`.
    - `poly_order`: the maximum order of polynomials.
    - `p0`: the initial parameters for the fitting function.

## Example
- Spline fit
```julia
nev = 60
name = "DX-+h"

Q1 = [30.66721918 29.860133 29.05306268 28.24612992 27.43911044 26.63197583 25.82507425 25.01797262 24.21096115 23.40387798 22.59690043 21.78991117 20.98281579 20.17570172 19.36884219 18.56169249 17.75463179 16.94772679 16.14061031 15.33347439 14.52663309 13.71945696 12.91240658 12.10544441 11.29841968 10.4913338 9.684370388 8.877289725 8.070184138]
E1 = [8.0902 7.5970 7.0749 6.5242 5.9461 5.3451 4.7300 4.1147 3.5182 2.9637 2.4769 2.0819 1.7972 1.6315 1.5800 1.6237 1.7301 1.8586 1.9687 2.0283 2.0190 2.0673 1.9910 2.0528 1.9730 2.0857 2.4550 3.1653 4.3448]

pot = Potential(); pot.name = name
pot.nev = nev
pot.Q0 = Q1[findmin(E1)[2]]; pot.E0 = 1.69834
pot.QE_data = DataFrame(Q = Q1[:], E = E1[:])
pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
pot.Q = Q

pot.func_type = "spline"
# spline fitting parameters
params = Dict()
params["order"] = 4
params["smoothness"] = 0.001
params["weight"] = [1 1 1 1 1 0.5 0.4 0.4 0.5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
fit_pot!(pot)
```

- Harmonic
```julia
nev = 40
name = "SnIII+h"

pot = Potential(); pot.name = name
pot.Q0 = 1.90291674728; pot.E0 = 0.585005
pot.nev = nev
pot.func_type = "harmonic"
params = Dict()
params["hw"] = 0.0281812646475
fit_pot!(pot)
```

"""
function fit_pot!(pot::Potential)
    # pot.params["E0"] = pot.E0
    # pot.params["Q0"] = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    E_CUT = 2 # defines upper energy limit on data used for fitting
    e_cut_ind = pot.QE_data.E .< E_CUT + pot.E0 # boolean

    println("Potential fitting: $(pot.name)")

    if pot.func_type == "bspline"
        println("=========bspline=========\n")

        bspline = get_bspline(pot.QE_data.Q, pot.QE_data.E)
        pot.E = bspline(pot.Q)
        pot.func = bspline

    elseif pot.func_type == "spline"
        println("=========spline==========\n")

        weight = get(pot.params, "weight", nothing)
        if isa(weight, String)
            weight = parse.(Float64, split(weight))
        end
        if weight ≠ nothing
            weight = vec(weight)
        end
        smoothness = get(pot.params, "smoothness", 0)
        order = get(pot.params, "order", 2)

        spline = get_spline(pot.QE_data.Q, pot.QE_data.E; weight = weight, smoothness = smoothness, order = order)
        pot.E = spline(pot.Q)
        pot.func = spline

    elseif pot.func_type == "harmonic"
        println("========harmonic=========\n")
        # func = x -> harmonic(x, params["hw"]; param = params)
        pot.E = harmonic.(pot.Q, pot.params["hw"]; E₀ = pot.E0, Q₀ = pot.Q0)
        pot.func = x -> harmonic(x, pot.params["hw"]; E₀ = pot.E0, Q₀ = pot.Q0)

    # potentials that need fitting
    else
        # self-consistent fitting
        if pot.func_type == "sc_harmonic"
        println("========self-consistent harmonic fit========\n")
        func = (x, p) -> harmonic_fittable(x, p; E₀ = pot.E0, Q₀ = pot.Q0)
        fit = sc_fit(func, pot.QE_data.Q[e_cut_ind], pot.QE_data.E[e_cut_ind], Float64.(pot.params["p0"]), pot.T)

        # one-shot fitting
        else
            if pot.func_type == "polyfunc"
                println("========polynomial========\n")
                func = (x, p) -> polyfunc(x, p; E₀ = pot.E0, Q₀ = pot.Q0, poly_order = pot.params["poly_order"])
            elseif pot.func_type == "morse_poly"
                println("========morse polynomial========\n")
                func = (x, p) -> morse_poly(x, p; E₀ = pot.E0, Q₀ = pot.Q0, poly_order = pot.params["poly_order"])
            elseif pot.func_type == "morse"
                func = (x, p) -> morse(x, p; E₀ = pot.E0, Q₀ = pot.Q0)
            elseif pot.func_type == "harmonic_fittable"
                println("========harmonic fit========\n")
                func = (x, p) -> harmonic_fittable(x, p; E₀ = pot.E0, Q₀ = pot.Q0)
            end
            fit = curve_fit(func, pot.QE_data.Q[e_cut_ind], pot.QE_data.E[e_cut_ind], Float64.(pot.params["p0"])) # curve_fit : model, x data, y data, p0
        end

        pot.E = func.(pot.Q, Ref(fit.param)) # when calling function, need to Ref() so that it can deal with array mismatch
        pot.func = x -> func(x, fit.param)

        println("===========Fit===========")
        println("Function: $(pot.func_type)")
        println("Best fit: $(fit.param)")
        println("=========================\n")
    end
end


"""
    solve_pot!(pot::Potential)

Solve 1D Shrödinger equation for `Potential`. `pot.ϵ` and `pot.χ` store the eigenvalues and eigenvectors, respectively.
"""
function solve_pot!(pot::Potential)
    # pot.ϵ, pot.χ = solve1D_ev_amu(pot.func, pot.Q; nev=pot.nev)
    # Q = pot.Q
    # NQ=length(Q); Qi=minimum(Q); Qf=maximum(Q)
    pot.ϵ, pot.χ = solve1D_ev_amu(pot.func, pot.Q; nev=pot.nev)
end


"""
Solve 1D Shrödinger equation. The Brooglie wrapper in the unit of `eV` and `amu`.
"""
function solve1D_ev_amu(func, Q::AbstractArray{<:Real,1}; nev=30, maxiter=nev*length(Q))
    NQ=length(Q); Qi=minimum(Q); Qf=maximum(Q)
    return solve1D_ev_amu(func, NQ, Qi, Qf; nev=nev, maxiter=maxiter)
end

function solve1D_ev_amu(func, NQ::Int, Qi::Real, Qf::Real; nev=30, maxiter=nev*length(Q))
    factor = (1/amu) * (ħc*1E10)^2

    ϵ1, χ1 = Brooglie.solve(x -> func.(x*factor^0.5);
                  N=NQ, a=Qi/factor^0.5, b=Qf/factor^0.5, m=1, nev=nev, maxiter=maxiter)
    return ϵ1, vcat(χ1'...)/factor^0.25
end


"""
    find_crossing(pot_1::Potential, pot_2::Potential)

Find the crossing point between two potential energy surfaces `pot1` and `pot2`.
`Qx, Ex = find_crossing(pot1, pot2)`.
"""
function find_crossing(pot_1::Potential, pot_2::Potential)
    # find root of pot_1 - pot_2 = 0
    diff_func = x-> pot_1.func(x) - pot_2.func(x)
    Q = pot_1.Q
    rts = find_zero(diff_func, Q[length(Q)÷2]) # start search at the midpoint of the Q vector
    return rts, pot_1.func.(rts)
end


# read potential
"""
Depreciated.
Construct `Potential` from `QE_data` and configure dictionary `cfg`.
"""
function pot_from_dict(QE_data::DataFrame, cfg::Dict)::Potential
    pot = Potential()
    pot.name = cfg["name"]
    pot.nev = cfg["nev"]
    pot.func_type = cfg["function"]["type"]
    # pot.p0 =  parse.(Float64, split(get(cfg["function"], "p0", "1 1")))
    pot.E0 = get(cfg, "E0", Inf)
    pot.QE_data = QE_data

    pot.params = get(cfg["function"], "params",  Dict("E0" => pot.E0))
    pot.params["p0"] =  parse.(Float64, split(get(cfg["function"], "p0", "1 1")))

    pot.Q0 = pot.QE_data.Q[findmin(pot.QE_data.E)[2]]

    if pot.E0 < Inf
        pot.QE_data.E .+= - minimum(pot.QE_data.E) + pot.E0
    end

    return pot
end

############################################################
# Set of potentials
function sqwell(x, width, depth; x0=0.)
    x = x .- x0
    pot_well = (x .< -width/2.) .| (x .> width/2.)
    return pot_well .* depth
end

function harmonic(x, ħω::Real; E₀, Q₀)
    a = amu/2 * (ħω/ħc/1E10)^2
    return a*(x - Q₀)^2 + E₀
end

function harmonic_fittable(x, coeff; E₀, Q₀)
    # an alternate syntax for harmonic() in order to make it fittable
    y = (0 .* x) .+ E₀ .+ coeff[1] * (x .- Q₀) .^2
    return y
end

function double_well(x, ħω1::Real, ħω2::Real; param)
    a1 = amu / 2. * (ħω1/ħc/1E10)^2
    a2 = amu / 2. * (ħω2/ħc/1E10)^2
    return - a1*x.*x + a2*x.*x.*x.*x
end

# Set of potentials for fitting
function polyfunc(x, coeffs; E₀::Real, Q₀::Real, poly_order::Int)
    # E0 = param["E0"]
    # Q₀ = param["Q0"]
    # poly_order = Int(param["poly_order"])
    y = 0 .* x # array of size x
    y = y .+ E₀

    for i = 2:poly_order + 1
        y = y .+ coeffs[i].* (x .- Q₀) .^(i-1)
    end

    return y
end

function morse(x, coeffs; E₀, Q₀)
    y = 0 .* x # array of size x
    y = y .+ E₀

    y = y .+ coeffs[1] .* ((1 .- exp.(-coeffs[2].*(x.-Q₀))).^2)

    return y
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
    # Qs, Es have to be equally-spaced (Range)
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

"""
    function cleave_pot(pot::Potential)

Return two potentials with the `QE_data` of the original potential split by the apex.
If there is more than one maximum, throw an error

"""
function cleave_pot(pot::Potential)
    # find maxima
    max_ind_array = findall(pot.E .== maximum(pot.E))
    len_max_ind = length(max_ind_array)
    # check if there are multiple maxima
    if len_max_ind > 1
        throw(ArgumentError("The fitted potential has several maxima."))
    end

    # cast unique maximum to int
    max_ind = max_ind_array[1]

    # find the corresponding Q value
    Q_limit = pot.Q[max_ind]

    # define two sub-potentials
    pot_a = copy(pot)
    pot_a.E = []
    pot_a.func = x->0

    pot_b = copy(pot_a)

    # assign corresponding data points
    Q_data = pot.QE_data.Q
    E_data = pot.QE_data.E

    pot_a.QE_data = DataFrame(Q = Q_data[Q_data.<=Q_limit], E = E_data[Q_data.<=Q_limit])
    pot_b.QE_data = DataFrame(Q = Q_data[Q_data.>=Q_limit], E = E_data[Q_data.>=Q_limit])

    return pot_a, pot_b
end

