# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTriangularComplex

using Random

Random.seed!(123)

include("testtriag.jl") # test_approx_eq_modphase

test_triangular((ComplexF32, ComplexF64, Complex{BigFloat}))

end # module TestTriangularComplex
