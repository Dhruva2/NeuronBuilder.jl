using NBuilder
using Documenter

DocMeta.setdocmeta!(NBuilder, :DocTestSetup, :(using NBuilder); recursive=true)

makedocs(;
    modules=[NBuilder],
    authors="Dhruva V. Raman, University of Sussex",
    repo="https://github.com/Dhruva2/NBuilder.jl/blob/{commit}{path}#{line}",
    sitename="NBuilder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dhruva2.github.io/NBuilder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Dhruva2/NBuilder.jl",
    devbranch="main",
)
