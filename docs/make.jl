using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Control loss",
    format = Documenter.HTML(
        prettyurls = false, 
        #size_threshold_ignore = [""],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages = [
        "Introduction"       => "index.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)