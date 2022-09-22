# This file is a part of Julia. License is MIT: https://julialang.org/license

include("util.jl")

# Make sure to run this before we do `import LinearAlgebra`
delete_methods_from("LinearAlgebra")

# Load our version of LinearAlgebra
import LinearAlgebra
using Test: @testset

@testset "LinearAlgebra" begin # TODO: delete this line
for file in readlines(joinpath(@__DIR__, "testgroups"))
    @testset "$(file)" begin
        include(file * ".jl")
    end
end
end # TODO: delete this line
