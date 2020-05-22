"""

A copy of [https://github.com/RedPointyJackson/Brooglie](https://github.com/RedPointyJackson/Brooglie)

"""
module Brooglie

export buildH, solve

using LinearAlgebra
using SparseArrays
using Base.Iterators
using Arpack


const H2eV = 27.21138602 # Hartree to electronvolt

"""
    combinations(els, N)

Return an array of all the `N` element combinations with repetition of `els`,
ordered lexicographically.
"""
function combinations(els, N)
    @assert N ≥ 0
    if N == 0
        return typeof(first(els))[]
    elseif N == 1
        return vec([[el] for el in els])
    else
        return vec([vcat(el, c) for c in combinations(els, N-1), el in els])
    end
end

"""
    numberofarguments(f)

Number of arguments of function `f`. Fails if `f` has more than one method.
"""
function numberofarguments(f)
    m = methods(f).ms
    if length(m) > 1
        error("Function $f has more than one method.")
    elseif length(m) < 1
        error("Function $f has no methods!.") # Happens for constants
    end
    return m[1].nargs - 1
end

"""
    buildskeleton(D, N)

Helper for building the Hamiltonian in `D` dimensions with `N` base
elements. Returns the Hamiltonian but with a diagonal of zeros.
"""
function buildskeleton(D, N)

    function ith_diag(i)
        motif = [-ones(N^(i-1)*(N-1)); zeros(N^(i-1))]
        return take(cycle(motif), N^D - N^(i-1)) |> collect
    end

    diags = ith_diag.(1:D)
    S = spdiagm([N^(i-1) => diags[i] for i in 1:D]...)

    return S + S'

end

"""
    buildH(V; N=20, a=-1, b=1, m=1)

Hamiltonian of a particle of mass `m` in a box spanning from `a` to `b` in all
D dimensions with a basis of `N` elements (number of partitions of the space in
each coordinate). The potential `V`(x, y, z, ...) is a function of D arguments.
"""
function buildH(V; N=20, a=-1, b=1, m=1)
    D = numberofarguments(V)
    xx = combinations(range(a, stop=b, length=N), D)
    h = (b-a)/N
    θ(ζ) = V(ζ...) * 2*m*h^2
    return buildskeleton(D, N) + spdiagm(0 => map(ζ->2D+θ(ζ), xx))
end

"""
    integrate(φ, L) = ∫φ dA = ∫∫⋯∫ φ(x,y,z,...) dx⋅dy⋅dz...

Multidimensional integration of the D-dimensional array `φ` in an
hypercube of side L.

Is assumed that the array has the same number of elements in each
dimension.
"""
function integrate(φ, L)
    @assert length(unique(size(φ))) == 1 "The array is has different widths."
    D = length(size(φ))
    N = size(φ,1)
    dA = (L/N)^D
    return sum(φ) * dA
end

"""
    normalizewf(φ, L) = φ / ∫|φ|² dA

    Normalize the wavefunction vector `φ` in a
    hypercube of side `L` so its modulus squared has unit area.

    The dimension of the hypercube is that of the wavefunction.
"""
normalizewf(φ, L) = φ / √integrate(φ.^2, L)

"""
    solve(V; N=500, a=-1, b=1, m=1, nev=N÷20, maxiter=1000)

Solve the potential `V`(x,y,z,...) in a grid xᵢ ∈ [`a`,`b`], discretized in `N`
steps.

The particle is assumed to have mass `m` (by default 1, a electron mass).

The function will return the `nev` first energy levels (in Hartree[^1])
and its normalized eigenfunctions.

[^1] A Hartree is equivalent to 27.21… eV. The global variable `H2eV`,
equal to that value, can be accessed under Brooglie.H2eV for
convenience.

"""
function solve(V; N=500, a=-1, b=1, m=1, nev=N÷20, maxiter=1000)
    D = numberofarguments(V)
    H = buildH(V; N=N, a=a, b=b, m=m)
    λλ = zeros(nev)
    v = zeros(N^D, nev)
    try
        λλ, v = eigs(H, nev=nev, which=:SR, maxiter=maxiter)
    catch ex
        if typeof(ex) == ARPACKException && ex.info == 1
            error("Arpack runned out of iters (maxiter = $maxiter)")
        end
    end
    h = (b-a)/N
    E = λλ / (2*m*h^2)
    vv = [reshape(v[:,i], Tuple(repeat([N], D))) for i in 1:nev]
    φφ = normalizewf.([vv[i] for i in 1:nev], b-a)
    return E, φφ
end

end#module