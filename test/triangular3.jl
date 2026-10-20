# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTriangularComplex

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Random, LinearAlgebra

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers: COMPREHENSIVE

Random.seed!(123)

include("testtriag.jl") # test_approx_eq_modphase

# ComplexF32 is also the element type used for the mixed-eltype tests, and ComplexF64
# takes the same BLAS paths
test_triangular((ComplexF32, Complex{BigFloat}, (@static COMPREHENSIVE ? (ComplexF64,) : ())...))

end # module TestTriangularComplex
