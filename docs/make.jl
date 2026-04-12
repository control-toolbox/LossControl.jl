# to run the documentation generation:
# julia --project=. docs/make.jl
pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, @__DIR__)

using Documenter

# ===============================================
# --- Copy Project Assets for reproducibility ---
# ===============================================
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force=true,
)
cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force=true,
)

# Repository URL (used for links in docs)
repo_url = "github.com/control-toolbox/LossControl.jl"

# ==============================
# --- Generate Documentation ---
# ==============================
# If draft is true, the Julia code in markdown is not executed.
# To disable draft mode in a specific markdown file, add:
#=
```@meta
Draft = false
```
=#
makedocs(;
    draft=false,
    warnonly=:cross_references,
    sitename="Loss control",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=[
            "statement.md", "numerical.md", "zermelo1.md", "zermelo2.md", "ho.md"
        ],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Introduction" => "index.md",
        "Mathematical background" => [
            "Optimal control and loss control" => "statement.md",
            "Numerical approach" => "numerical.md",
        ],
        "Examples" => [
            "Zermelo navigation: Example 1" => "zermelo1.md",
            "Zermelo navigation: Example 2" => "zermelo2.md",
            "Harmonic oscillator problem" => "ho.md",
        ],
    ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main", push_preview=true)
# push_preview: use https://control-toolbox.org/LossControl.jl/previews/PRXXX where XXX is the pull request number
