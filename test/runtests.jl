# This file is a part of Julia. License is MIT: https://julialang.org/license

for file in readlines(joinpath(@__DIR__, "testgroups"))
    @info "Testing $file"
    include(file * ".jl")
end
