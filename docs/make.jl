using Documenter, PAPA

makedocs(;
    modules=[PAPA],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://qiplab.bbn.com/matthewware/PAPA.jl/blob/{commit}{path}#L{line}",
    sitename="PAPA.jl",
    authors="Luke C. Govia",
    assets=String[],
)
