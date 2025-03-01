using DihedralGroups
using Documenter

DocMeta.setdocmeta!(DihedralGroups, :DocTestSetup, :(using DihedralGroups); recursive=true)

makedocs(;
    modules=[DihedralGroups],
    authors="LauBMo <laurea987@gmail.com> and contributors",
    sitename="DihedralGroups.jl",
    format=Documenter.HTML(;
        canonical="https://LauraBMo.github.io/DihedralGroups.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/DihedralGroups.jl",
    devbranch="main",
)
