# This file is a part of Julia. License is MIT: https://julialang.org/license
using Test, LinearAlgebra

for file in readlines(joinpath(@__DIR__, "testgroups"))
    @info "Testing $file"
    include(file * ".jl")
end

@testset "Docstrings" begin
    @test isempty(Docs.undocumented_names(LinearAlgebra))
end

@testset "versioninfo" begin
    vinfo = sprint(LinearAlgebra.versioninfo)
    @test occursin("Threading:", vinfo)
    @test occursin(r"Threads.threadpoolsize\(\) = [0-9]+", vinfo)
    @test occursin(r"Threads.maxthreadid\(\) = [0-9]+", vinfo)
    @test occursin(r"LinearAlgebra.BLAS.get_num_threads\(\) = [0-9]+", vinfo)
    @test occursin("Relevant environment variables:", vinfo)
    vars = strip(split(vinfo, "Relevant environment variables:")[end])
    @test any(occursin(vars), [r"JULIA_NUM_THREADS = [0-9]+", r"MKL_DYNAMIC = [0-9]+",
                r"MKL_NUM_THREADS = [0-9]+",
                r"OPENBLAS_NUM_THREADS = [0-9]+",
                r"GOTO_NUM_THREADS = [0-9]+",
                r"OMP_NUM_THREADS = [0-9]+", r"\[none\]"])

    withenv("MKL_NUM_THREADS" => 1) do
        vinfo = sprint(LinearAlgebra.versioninfo)
        vars = strip(split(vinfo, "Relevant environment variables:")[end])
        @test occursin("MKL_NUM_THREADS = 1", vars)
    end
end
