using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Control loss",
    format = Documenter.HTML(
        prettyurls = false, 
        size_threshold_ignore = ["statement.md", "numerical.md","zermelo.md", "ho.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages = [
        "Introduction"                       => "index.md",
        "Statement of the problem"           => "statement.md",
        "Numerical approach"                 => "numerical.md",
        "Zermelo navigation problem"         => "zermelo.md",
        "Harmonic oscillator problem"        => "ho.md",

    ]
)

deploydocs(
    repo = "github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)
