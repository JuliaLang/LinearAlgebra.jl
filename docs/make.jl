docs_dir = @__DIR__
root_dir = dirname(docs_dir)
test_dir = joinpath(root_dir, "test")
include(joinpath(test_dir, "util.jl"))
    
# Make sure to run this before we do `import LinearAlgebra`
delete_all_methods()

# Load our version of LinearAlgebra
using LinearAlgebra

using Documenter: DocMeta, makedocs, deploydocs

DocMeta.setdocmeta!(LinearAlgebra, :DocTestSetup, :(using LinearAlgebra); recursive=true)

makedocs(
    modules = [LinearAlgebra],
    sitename = "LinearAlgebra",
    pages = Any[
        "LinearAlgebra" => "index.md",
        ];
    # strict = true,
    strict = Symbol[:doctest],
    )

deploydocs(repo = "github.com/JuliaLang/LinearAlgebra.jl.git")
