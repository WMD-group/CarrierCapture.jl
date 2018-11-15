"""
Solve the time independent Schödinger equation for a particle in an
arbitrary potential. See the documentation for `Solve1D`, `Solve2D`
and `Solve3D` and the PDF documentation in the `/doc` folder.
"""
module Brooglie
export solve1D, solve2D, solve3D

using Base.Iterators
using SparseArrays, Arpack

const H2eV = 27.21138602 # Hartree to electronvolt
"""
    integrate(φ, L)

Multidimensional integration of the D-dimensional array `φ` in an
hypercube of side L.

Is assumed that the array has the same number of elements in each
dimension.
"""
function integrate(φ::Array{T}, L) where {T<:Number}
    @assert φ|>size|>unique|>length == 1 "The array is has different widths."
    D = length(size(φ))
    N = size(φ,1)
    basehyperarea = (L/N)^D
    return sum(φ) * basehyperarea
end
"""
    normalizewf(φ,L)

    Normalize the wavefunction vector `φ` in a
    hypercube of side `L` so its modulus squared has unit area.

    The dimension of the hypercube is that of the wavefunction.
"""
normalizewf(φ, L) = φ / √integrate(φ.^2, L)
"""
    solve1D(V; N=500, a=-1, b=1,
               m=1, nev=N÷20, maxiter=1000)

Solve the 1D potential `V`(x) in a grid x ∈ [`a`,`b`], discretized in
`N` steps.

The particle is asumed to have mass `m` (by default the electron
mass).

The function will return the `nev` first energy levels (in Hartree[^1])
and its normalized eigenfunctions.

[^1] A Hartree is equivalent to 27.21 eV. The global variable `H2eV`
can be accessed under Brooglie.H2eV for convenience.

"""
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
"""
    solve2D(V; N=250, a=-1, b=1,
               m=1, nev=N÷20, maxiter=1000)

Identical to `solve1D`, but for a 2D grid `x`,`y` ∈ [`a`,`b`].
"""
function solve2D(V; N=250, a=-1, b=1, m=1, nev=N÷20, maxiter=1000)
    nev = trunc(Int, nev)
    ε = (b-a)/N
    Θ(x,y) = V(x,y) * 2*m*ε^2
    ll = range(a, stop=b, length=N)
    rr = vec([(x,y) for x in ll, y in ll])
    d  = 4 .+ [Θ(r...) for r in rr]
    dd = take(cycle([-ones(N-1); 0]), N^2-1) |> collect
    ddd = -ones(N*(N-1))
    M  = spdiagm(-N => ddd, -1 => dd, 0 => d, 1 => dd, -N => ddd)
    λ, v = eigs(M, nev=nev, which=:SR, maxiter=maxiter)
    E = λ / (2*m*ε^2)
    return E, [normalizewf(reshape(v[:,i], (N,N)), b-a) for i in 1:nev]
end
"""
    solve3D(V; N=100, a=-1, b=1,
               m=1, nev=N÷20, maxiter=1000)

Like `solve1D`, but for a 3D grid `x`,`y`,`z` ∈ [`a`,`b`].
"""
function solve3D(V; N=100, a=-1, b=1, m=1, nev=N÷20, maxiter=1000)
    nev = trunc(Int, nev)
    ε = (b-a)/N
    Θ(x,y,z) = V(x,y,z) * 2*m*ε^2
    ll = range(a, stop=b, length=N)
    rr = vec([(x,y,z) for x in ll, y in ll, z in ll])
    d  = 6 .+ [Θ(r...) for r in rr]
    dd = take(cycle([-ones(N-1); 0]), N^3-1) |> collect
    ddd = take(cycle([-ones(N*(N-1)); zeros(N)]), N^3-N) |> collect
    dddd = -ones(N^2*(N-1))
    M  = spdiagm(-N^2 => dddd, -N => ddd, -1 => dd, 0 => d, 1 => dd, N => ddd, N^2 => dddd)
    λ, v = eigs(M, nev=nev, which=:SR, maxiter=maxiter)
    E = λ / (2*m*ε^2)
    return E, [reshape(v[:,i], (N,N,N)) for i in 1:nev]
end

end#module
