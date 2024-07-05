using Documenter
using CarrierCapture

makedocs(
    sitename = "CarrierCapture",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    checkdocs = :none,  # Disable strict checking for missing docstrings
    modules = [CarrierCapture],
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Library" => Any[
            "Public" => "lib/public.md",
            "Brooglie" => "lib/brooglie.md",
            "Plotter" => "lib/plotter.md"
        ],
    ]
)

deploydocs(
    target = "build",
    repo   = "github.com/WMD-group/CarrierCapture.jl.git",
)
