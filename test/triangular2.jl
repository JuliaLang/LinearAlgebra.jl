# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTriangularReal

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Random

Random.seed!(123)

include("testtriag.jl") # test_approx_eq_modphase

test_triangular((Float32, Float64, BigFloat, Int))

end # module TestTriangularReal
