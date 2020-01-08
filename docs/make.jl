using Documenter
using CarrierCapture

makedocs(
    sitename = "CarrierCapture",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [CarrierCapture],
    pages = [
        "Home" => "index.md",
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/public.md"
        ],
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
