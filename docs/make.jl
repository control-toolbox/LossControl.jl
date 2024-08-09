using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Control loss",
    format = Documenter.HTML(
        prettyurls = false, 
<<<<<<< HEAD
        #size_threshold_ignore = [""],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages = [
        "Introduction"       => "index.md",
=======
        size_threshold_ignore = ["zermelo1.md", "zermelo2.md"],
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

>>>>>>> adjoint
    ]
)

deploydocs(
<<<<<<< HEAD
    repo = "github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)
=======
    repo = "https://github.com/control-toolbox/control-loss.git",
    devbranch = "main"
)


>>>>>>> adjoint
