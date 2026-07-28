# This file is a part of Julia. License is MIT: https://julialang.org/license

include("prune_old_LA.jl")

using Test, LinearAlgebra

testfiles = String[]
for file in readlines(joinpath(@__DIR__, "testgroups"))
    push!(testfiles, file * ".jl")
end

# ParallelTestRunner comes from the Pkg.test target; Julia base CI runs this
# file without it and falls back to the serial path.
if Base.find_package("ParallelTestRunner") !== nothing
    using ParallelTestRunner
    # Auto CPU thread count detection in ParallelTestRunner is bad
    push!(ARGS, "--jobs=$(Sys.CPU_THREADS)")
    testsuite = Dict{String,Expr}(splitext(f)[1] => :(include($(joinpath(@__DIR__, f))))
                                  for f in testfiles)
    runtests(LinearAlgebra, ARGS; testsuite,
         init_worker_code = :(include($(joinpath(@__DIR__, "prune_old_LA.jl")))))
else
    foreach(include, testfiles)
end

@testset "Docstrings" begin
    @test isempty(Docs.undocumented_names(LinearAlgebra))
end
