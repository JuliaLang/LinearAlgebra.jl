# This file is a part of Julia. License is MIT: https://julialang.org/license

include("prune_old_LA.jl")

using Test, LinearAlgebra

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

testfiles = String[]
for file in readlines(joinpath(@__DIR__, "testgroups"))
    push!(testfiles, file * ".jl")
end

# Optionally run only a subset of the test files, e.g. LINALG_TEST_SHARD="1/2".
# `testgroups` is sorted by decreasing runtime, and files are dealt out to the
# shards in a snake order (1, 2, 2, 1, 1, 2, ...) to balance the total runtime.
let shard = get(ENV, "LINALG_TEST_SHARD", "")
    if !isempty(shard)
        i, n = parse.(Int, split(shard, '/'))
        1 <= i <= n || error("invalid LINALG_TEST_SHARD=$shard")
        snake(p) = (r = fld(p - 1, n); k = mod(p - 1, n) + 1; iseven(r) ? k : n + 1 - k)
        keepat!(testfiles, [p for p in eachindex(testfiles) if snake(p) == i])
    end
end

# ParallelTestRunner comes from the Pkg.test target; Julia base CI runs this
# file without it and falls back to the serial path.
if Base.find_package("ParallelTestRunner") !== nothing
    using ParallelTestRunner
    LinearAlgebra.versioninfo()
    push!(ARGS, "--verbose")
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
