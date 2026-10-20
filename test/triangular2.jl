# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTriangularReal

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Random, LinearAlgebra

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers: COMPREHENSIVE

Random.seed!(123)

include("testtriag.jl") # test_approx_eq_modphase

# Float32 takes the same BLAS paths as Float64
test_triangular((Float64, BigFloat, Int, (@static COMPREHENSIVE ? (Float32,) : ())...))

end # module TestTriangularReal
