using MRSI
using Documenter

DocMeta.setdocmeta!(MRSI, :DocTestSetup, :(using MRSI); recursive=true)

makedocs(;
    modules=[MRSI],
    authors="Korbinian Eckstein korbinian90@gmail.com",
    repo="https://github.com/korbinian90/MRSI.jl/blob/{commit}{path}#{line}",
    sitename="MRSI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://korbinian90.github.io/MRSI.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/korbinian90/MRSI.jl",
    devbranch="main",
)
