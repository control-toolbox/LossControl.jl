using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Control loss",
    format = Documenter.HTML(
        prettyurls = false, 
        size_threshold_ignore = ["zermelo1.md", "zermelo2.md", "ho.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
        inventory_version = "0.1.0"  # Add this line
    ),
    pages = [
        "Introduction"              => "index.md",
        "Zermelo example 1"         => "zermelo1.md",
        "Zermelo example 2"         => "zermelo2.md",
        "Harmonic oscillator"       => "ho.md",

    ]
)

deploydocs(
    repo = "https://github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)
