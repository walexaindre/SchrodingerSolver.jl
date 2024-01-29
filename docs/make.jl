using SchrodingerSolver
using Documenter

DocMeta.setdocmeta!(SchrodingerSolver, :DocTestSetup, :(using SchrodingerSolver); recursive=true)

makedocs(;
    modules=[SchrodingerSolver],
    authors="walexaindre <walexaindre@hotmail.com
> and contributors",
    sitename="SchrodingerSolver.jl",
    format=Documenter.HTML(;
        canonical="https://walexaindre.github.io/SchrodingerSolver.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/walexaindre/SchrodingerSolver.jl",
    devbranch="master",
)
