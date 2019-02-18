const H2eV = 27.21138602 # Hartree to electronvolt

function integrate(φ::Array{T}, L) where {T<:Number}
    @assert φ|>size|>unique|>length == 1 "The array is has different widths."
    D = length(size(φ))
    N = size(φ,1)
    basehyperarea = (L/N)^D
    return sum(φ) * basehyperarea
end

normalizewf(φ, L) = φ / √integrate(φ.^2, L)

function solve1D(V; N=500, a=-1, b=1, m=1, nev=N÷20, maxiter=1000)
    #    2 + V'₁   -1        0     ⋯
    #     -1     2 + V'₂     -1    ⋯
    #      0       -1     2 + V'₃  ⋯
    #      ⋮      ⋮        ⋮     ⋱
    nev = trunc(Int, nev)
    ε = (b-a)/N
    Θ(x) = V(x) * 2*m*ε^2
    xx = range(a, stop=b, length=N)
    d  = 2 .+ Θ.(xx)
    dd = -ones(N-1)
    M  = spdiagm(-1 => dd, 0 => d, 1 => dd)
    try
        λ, v = eigs(M, nev=nev, which=:SR, maxiter=maxiter)
        E = λ / (2*m*ε^2)
        return E, normalizewf.([v[:,i] for i in 1:nev], b-a)
    catch e
        @warn("Eigenvalue computation threw an error. Try the following:")
        @warn("    - Use the official binaries if there ", prefix="")
        @warn("      is an 'unexpected error'", prefix="")
        @warn("    - Increase the number of iterations or decrease N", prefix="")
        @warn("      if you get unespecified ARPACK error 1.", prefix="")
        throw(e)
    end
end

export solve1D