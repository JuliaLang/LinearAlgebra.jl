# This file is a part of Julia. License is MIT: https://julialang.org/license

include("util.jl")

# Make sure to run this before we do `import LinearAlgebra`
delete_all_methods()

# Load our version of LinearAlgebra
import LinearAlgebra

for file in readlines(joinpath(@__DIR__, "testgroups"))
    include(file * ".jl")
end
