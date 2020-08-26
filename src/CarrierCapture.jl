"""
Main module for `CarrierCapture.jl` -- A set of codes to compute carrier capture and recombination rates in semiconducting compounds.

Two structs are exported from this module for public use:

- [`Potential`](@ref): Potential.
- [`conf_coord`](@ref): Configuration coordinate.
- [`Plotter`](@ref): Auxiliary submodule for plotting.
- [`Brooglie`](@ref): 1D Shrödinger equation solver.

"""
module CarrierCapture
    using DocStringExtensions
    # using Reexport
    using LsqFit, Polynomials, Dierckx, Interpolations, Roots

    # using Plots, LaTeXStrings
    using DataFrames
    # @reexport using Plots, LaTeXStrings
    # @reexport using DataFrames

    include("Brooglie.jl")

    include("Potential.jl")
    include("CaptureRate.jl")
    include("Plotter.jl")
    include("paramScan.jl")
end # module
