using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Control loss",
    format = Documenter.HTML(
        prettyurls = false, 
        size_threshold_ignore = ["zermelo.md","statement.md", "ho.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
        inventory_version = "0.1.0"  # Add this line
    ),
    pages = [
        "Introduction"              => "index.md",
        "PMP with loss control regions and numerical methods" => "statement.md",
        "Zermelo with loss control regions" => "zermelo.md",
        "Harmonic oscillator with loss control regions" => "ho.md",

    ]
)

deploydocs(
    repo = "https://github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)


