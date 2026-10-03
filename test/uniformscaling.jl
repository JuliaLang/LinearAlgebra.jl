# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestUniformscaling

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.Quaternions
using Main.LinearAlgebraTestHelpers.OffsetArrays

Random.seed!(1234543)

@testset "basic functions" begin
    @test I === I' # transpose
    @test ndims(I) == 2
    @test one(UniformScaling{Float32}) == UniformScaling(one(Float32))
    @test zero(UniformScaling{Float32}) == UniformScaling(zero(Float32))
    @test eltype(one(UniformScaling{Float32})) == Float32
    @test zero(UniformScaling(rand(ComplexF64))) == zero(UniformScaling{ComplexF64})
    @test one(UniformScaling(rand(ComplexF64))) == one(UniformScaling{ComplexF64})
    @test eltype(one(UniformScaling(rand(ComplexF64)))) == ComplexF64
    @test -one(UniformScaling(2)) == UniformScaling(-1)
    @test opnorm(UniformScaling(1+im)) ≈ sqrt(2)
    @test convert(UniformScaling{Float64}, 2I) === 2.0I
    @test float(2I) === 2.0*I
end

@testset "getindex" begin
    @test I[1,1] == I[CartesianIndex(1,1)] == 1
    @test I[1,2] == I[CartesianIndex(1,2)] == 0

    J = I(15)
    for (a, b) in [
        # indexing that returns a Vector
        (1:10, 1),
        (4, 1:10),
        (11, 1:10),
        # indexing that returns a Matrix
        (1:2, 1:2),
        (1:2:3, 1:2:3),
        (1:2:8, 2:2:9),
        (1:2:8, 9:-4:1),
        (9:-4:1, 1:2:8),
        (2:3, 1:2),
        (2:-1:1, 1:2),
        (1:2:9, 5:2:13),
        (1, [1,2,5]),
        (1, [1,10,5,2]),
        (10, [10]),
        ([1], 1),
        ([15,1,5,2], 6),
        ([2], [2]),
        ([2,9,8,2,1], [2,8,4,3,1]),
        ([8,3,5,3], 2:9),
    ]
        @test I[a,b] == J[a,b]
        ndims(a) == 1 && @test I[OffsetArray(a,-10),b] == J[OffsetArray(a,-10),b]
        ndims(b) == 1 && @test I[a,OffsetArray(b,-9)] == J[a,OffsetArray(b,-9)]
        ndims(a) == ndims(b) == 1 && @test I[OffsetArray(a,-7),OffsetArray(b,-8)] == J[OffsetArray(a,-7),OffsetArray(b,-8)]
    end
end

@testset "sqrt, exp, log, and trigonometric functions" begin
    # convert to a dense matrix with random size
    M(J) = (N = rand(1:10); Matrix(J, N, N))

    # on complex plane
    J = UniformScaling(randn(ComplexF64))
    for f in ( exp,   log, cis,
               sqrt,  abs,
               sin,   cos,   tan,
               asin,  acos,  atan,
               csc,   sec,   cot,
               acsc,  asec,  acot,
               sinh,  cosh,  tanh,
               asinh, acosh, atanh,
               csch,  sech,  coth,
               acsch, asech, acoth )
        @test f(J) ≈ f(M(J))
    end

    for f in (sincos, sincosd)
        @test all(splat(≈), zip(f(J), f(M(J))))
    end

    # on real axis
    for (λ, fs) in (
        # functions defined for x ∈ ℝ
        (()->randn(),           (exp,   cbrt,
                                 sin,   cos,   tan,
                                 csc,   sec,   cot,
                                 atan,  acot,
                                 sinh,  cosh,  tanh,
                                 csch,  sech,  coth,
                                 asinh, acsch)),
        # functions defined for x ≥ 0
        (()->abs(randn()),      (log,   sqrt)),
        # functions defined for -1 ≤ x ≤ 1
        (()->2rand()-1,         (asin,  acos,  atanh)),
        # functions defined for x ≤ -1 or x ≥ 1
        (()->1/(2rand()-1),     (acsc,  asec,  acoth)),
        # functions defined for 0 ≤ x ≤ 1
        (()->rand(),            (asech,)),
        # functions defined for x ≥ 1
        (()->1/rand(),          (acosh,))
    )
        for f in fs
            J = UniformScaling(λ())
            @test f(J) ≈ f(M(J))
        end
    end
end

@testset "conjugation of UniformScaling" begin
    @test conj(UniformScaling(1))::UniformScaling{Int} == UniformScaling(1)
    @test conj(UniformScaling(1.0))::UniformScaling{Float64} == UniformScaling(1.0)
    @test conj(UniformScaling(1+1im))::UniformScaling{Complex{Int}} == UniformScaling(1-1im)
    @test conj(UniformScaling(1.0+1.0im))::UniformScaling{ComplexF64} == UniformScaling(1.0-1.0im)
end

@testset "isdiag, istriu, istril, issymmetric, ishermitian, isposdef, isapprox" begin
    @test isdiag(I)
    @test istriu(I)
    @test istril(I)
    @test issymmetric(I)
    @test issymmetric(UniformScaling(complex(1.0,1.0)))
    @test ishermitian(I)
    @test !ishermitian(UniformScaling(complex(1.0,1.0)))
    @test isposdef(UniformScaling(rand()))
    @test !isposdef(UniformScaling(-rand()))
    @test !isposdef(UniformScaling(randn(ComplexF64)))
    @test !isposdef(UniformScaling(NaN))
    @test isposdef(I)
    @test !isposdef(-I)
    @test isposdef(UniformScaling(complex(1.0, 0.0)))
    @test !isposdef(UniformScaling(complex(1.0, 1.0)))
    @test UniformScaling(4.00000000000001) ≈ UniformScaling(4.0)
    @test UniformScaling(4.32) ≈ UniformScaling(4.3) rtol=0.1 atol=0.01
    @test UniformScaling(4.32) ≈ 4.3 * [1 0; 0 1] rtol=0.1 atol=0.01
    @test UniformScaling(4.32) ≈ 4.3 * [1 0; 0 1] rtol=0.1 atol=0.01 norm=norm
    @test 4.3 * [1 0; 0 1] ≈ UniformScaling(4.32) rtol=0.1 atol=0.01
    @test [4.3201 0.002;0.001 4.32009] ≈ UniformScaling(4.32) rtol=0.1 atol=0.
    @test UniformScaling(4.32) ≉ fill(4.3,2,2) rtol=0.1 atol=0.01
    @test UniformScaling(4.32) ≈ 4.32 * [1 0; 0 1]
end

@testset "arithmetic with Number" begin
    α = rand()
    @test α + I == α + 1
    @test I + α == α + 1
    @test α - I == α - 1
    @test I - α == 1 - α
    @test α .* UniformScaling(1.0) == UniformScaling(1.0) .* α
    @test UniformScaling(α)./α == UniformScaling(1.0)
    @test α.\UniformScaling(α) == UniformScaling(1.0)
    @test α * UniformScaling(1.0) == UniformScaling(1.0) * α
    @test UniformScaling(α)/α == UniformScaling(1.0)
    @test 2I//3 == (2//3)*I
    @test (2I)^α == (2I).^α == (2^α)I

    β = rand()
    @test (α*I)^2    == UniformScaling(α^2)
    @test (α*I)^(-2) == UniformScaling(α^(-2))
    @test (α*I)^(.5) == UniformScaling(α^(.5))
    @test (α*I)^β    == UniformScaling(α^β)

    @test (α * I) .^ 2 == UniformScaling(α^2)
    @test (α * I) .^ β == UniformScaling(α^β)
end

@testset "unary" begin
    @test +I === +1*I
    @test -I === -1*I
end

@testset "tr, det and logdet" begin
    for T in (Int, Float64, ComplexF64, Bool)
        @test tr(UniformScaling(zero(T))) === zero(T)
    end
    @test_throws ArgumentError tr(UniformScaling(1))
    @test det(I) === true
    @test det(1.0I) === 1.0
    @test det(0I) === 0
    @test det(0.0I) === 0.0
    @test logdet(I) == 0
    @test_throws ArgumentError det(2I)
end

@test copy(UniformScaling(one(Float64))) == UniformScaling(one(Float64))
@test sprint(show,MIME"text/plain"(),UniformScaling(one(ComplexF64))) == "$(LinearAlgebra.UniformScaling){ComplexF64}\n(1.0 + 0.0im)*I"
@test sprint(show,MIME"text/plain"(),UniformScaling(one(Float32))) == "$(LinearAlgebra.UniformScaling){Float32}\n1.0*I"
@test sprint(show,UniformScaling(one(ComplexF64))) == "$(LinearAlgebra.UniformScaling){ComplexF64}(1.0 + 0.0im)"
@test sprint(show,UniformScaling(one(Float32))) == "$(LinearAlgebra.UniformScaling){Float32}(1.0f0)"

let
    λ = complex(randn(),randn())
    J = UniformScaling(λ)
    @testset "transpose, conj, inv, pinv, cond" begin
        @test ndims(J) == 2
        @test transpose(J) == J
        @test J * [1 0; 0 1] == conj(*(adjoint(J), [1 0; 0 1])) # ctranpose (and A(c)_mul_B)
        @test I + I === UniformScaling(2) # +
        @test inv(I) == I
        @test inv(J) == UniformScaling(inv(λ))
        @test pinv(J) == UniformScaling(inv(λ))
        @test @inferred(pinv(0.0I)) == 0.0I
        @test @inferred(pinv(0I)) == 0.0I
        @test @inferred(pinv(false*I)) == 0.0I
        @test @inferred(pinv(0im*I)) == 0im*I
        @test cond(I) == 1
        @test cond(J) == (λ ≠ zero(λ) ? one(real(λ)) : oftype(real(λ), Inf))
    end

    @testset "real, imag, reim" begin
        @test real(J) == UniformScaling(real(λ))
        @test imag(J) == UniformScaling(imag(λ))
        @test reim(J) == (UniformScaling(real(λ)), UniformScaling(imag(λ)))
    end

    @testset "copyto!" begin
        A = Matrix{Int}(undef, (3,3))
        @test copyto!(A, I) == one(A)
        B = Matrix{ComplexF64}(undef, (1,2))
        @test copyto!(B, J) == [λ zero(λ)]
    end

    @testset "copy!" begin
        A = Matrix{Int}(undef, (3,3))
        @test copy!(A, I) == one(A)
        B = Matrix{ComplexF64}(undef, (1,2))
        @test copy!(B, J) == [λ zero(λ)]
    end

    @testset "binary ops with vectors" begin
        v = complex.(randn(3), randn(3))
        # As shown in #20423@GitHub, vector acts like x1 matrix when participating in linear algebra
        @test v  * J ≈ v  * λ
        @test v' * J ≈ v' * λ
        @test J * v  ≈ λ * v
        @test J * v' ≈ λ * v'
        @test v  / J ≈ v  / λ
        @test v' / J ≈ v' / λ
        @test J \ v  ≈ λ \ v
        @test J \ v' ≈ λ \ v'
    end

    @testset "binary ops with matrices" begin
        B = bitrand(2, 2)
        @test B + I == B + Matrix(I, size(B))
        @test I + B == B + Matrix(I, size(B))
        AA = randn(2, 2)
        for A in (AA, view(AA, 1:2, 1:2))
            I22 = Matrix(I, size(A))
            @test @inferred(A + I) == A + I22
            @test @inferred(I + A) == A + I22
            @test @inferred(I - I) === UniformScaling(0)
            @test @inferred(B - I) == B - I22
            @test @inferred(I - B) == I22 - B
            @test @inferred(A - I) == A - I22
            @test @inferred(I - A) == I22 - A
            @test @inferred(I*J) === UniformScaling(λ)
            @test @inferred(B*J) == B*λ
            @test @inferred(J*B) == B*λ
            @test @inferred(I*A) !== A # Don't alias
            @test @inferred(A*I) !== A # Don't alias

            @test @inferred(A*J) == A*λ
            @test @inferred(J*A) == A*λ
            @test @inferred(J*fill(1, 3)) == fill(λ, 3)
            @test @inferred(λ*J) === UniformScaling(λ*J.λ)
            @test @inferred(J*λ) === UniformScaling(λ*J.λ)
            @test @inferred(J/I) === J
            @test @inferred(I/A) == inv(A)
            @test @inferred(A/I) == A
            @test @inferred(I/λ) === UniformScaling(1/λ)
            @test @inferred(I\J) === J

            if isa(A, Array)
                T = LowerTriangular(randn(3,3))
            else
                T = LowerTriangular(view(randn(3,3), 1:3, 1:3))
            end
            @test @inferred(T + J) == Array(T) + J
            @test @inferred(J + T) == J + Array(T)
            @test @inferred(T - J) == Array(T) - J
            @test @inferred(J - T) == J - Array(T)
            @test @inferred(T\I) == inv(T)

            if isa(A, Array)
                T = LinearAlgebra.UnitLowerTriangular(randn(3,3))
            else
                T = LinearAlgebra.UnitLowerTriangular(view(randn(3,3), 1:3, 1:3))
            end
            @test @inferred(T + J) == Array(T) + J
            @test @inferred(J + T) == J + Array(T)
            @test @inferred(T - J) == Array(T) - J
            @test @inferred(J - T) == J - Array(T)
            @test @inferred(T\I) == inv(T)

            if isa(A, Array)
                T = UpperTriangular(randn(3,3))
            else
                T = UpperTriangular(view(randn(3,3), 1:3, 1:3))
            end
            @test @inferred(T + J) == Array(T) + J
            @test @inferred(J + T) == J + Array(T)
            @test @inferred(T - J) == Array(T) - J
            @test @inferred(J - T) == J - Array(T)
            @test @inferred(T\I) == inv(T)

            if isa(A, Array)
                T = LinearAlgebra.UnitUpperTriangular(randn(3,3))
            else
                T = LinearAlgebra.UnitUpperTriangular(view(randn(3,3), 1:3, 1:3))
            end
            @test @inferred(T + J) == Array(T) + J
            @test @inferred(J + T) == J + Array(T)
            @test @inferred(T - J) == Array(T) - J
            @test @inferred(J - T) == J - Array(T)
            @test @inferred(T\I) == inv(T)

            for elty in (Float64, ComplexF64)
                if isa(A, Array)
                    T = Hermitian(randn(elty, 3,3))
                else
                    T = Hermitian(view(randn(elty, 3,3), 1:3, 1:3))
                end
                @test @inferred(T + J) == Array(T) + J
                @test @inferred(J + T) == J + Array(T)
                @test @inferred(T - J) == Array(T) - J
                @test @inferred(J - T) == J - Array(T)
            end

            @testset for f in (transpose, adjoint)
                if isa(A, Array)
                    T = f(randn(ComplexF64,3,3))
                else
                    T = f(view(randn(ComplexF64,3,3), 1:3, 1:3))
                end
                TA = Array(T)
                @test @inferred(T + J) == TA + J
                @test @inferred(J + T) == J + TA
                @test @inferred(T - J) == TA - J
                @test @inferred(J - T) == J - TA
            end

            @test @inferred(I\A) == A
            @test @inferred(A\I) == inv(A)
            @test @inferred(λ\I) === UniformScaling(1/λ)
        end
    end
end

@testset "hcat and vcat" begin
    @test_throws ArgumentError hcat(I)
    @test_throws ArgumentError [I I]
    @test_throws ArgumentError vcat(I)
    @test_throws ArgumentError [I; I]
    @test_throws ArgumentError [I I; I]

    A = rand(3,4)
    B = rand(3,3)
    C = rand(0,3)
    D = rand(2,0)
    E = rand(1,3)
    F = rand(3,1)
    α = rand()
    @test (hcat(A, 2I))::Matrix == hcat(A, Matrix(2I, 3, 3))
    @test (hcat(E, α))::Matrix == hcat(E, [α])
    @test (hcat(E, α, 2I))::Matrix == hcat(E, [α], fill(2, 1, 1))
    @test (vcat(A, 2I))::Matrix == vcat(A, Matrix(2I, 4, 4))
    @test (vcat(F, α))::Matrix == vcat(F, [α])
    @test (vcat(F, α, 2I))::Matrix == vcat(F, [α], fill(2, 1, 1))
    @test (hcat(C, 2I))::Matrix == C
    @test_throws DimensionMismatch hcat(C, α)
    @test (vcat(D, 2I))::Matrix == D
    @test_throws DimensionMismatch vcat(D, α)
    @test (hcat(I, 3I, A, 2I))::Matrix == hcat(Matrix(I, 3, 3), Matrix(3I, 3, 3), A, Matrix(2I, 3, 3))
    @test (vcat(I, 3I, A, 2I))::Matrix == vcat(Matrix(I, 4, 4), Matrix(3I, 4, 4), A, Matrix(2I, 4, 4))
    @test (hvcat((2,1,2), B, 2I, I, 3I, 4I))::Matrix ==
        hvcat((2,1,2), B, Matrix(2I, 3, 3), Matrix(I, 6, 6), Matrix(3I, 3, 3), Matrix(4I, 3, 3))
    @test hvcat((3,1), C, C, I, 3I)::Matrix == hvcat((2,1), C, C, Matrix(3I, 6,6))
    @test hvcat((2,2,2), I, 2I, 3I, 4I, C, C)::Matrix ==
        hvcat((2,2,2), Matrix(I, 3, 3), Matrix(2I, 3,3 ), Matrix(3I, 3,3), Matrix(4I, 3,3), C, C)
    @test hvcat((2,2,4), C, C, I, 2I, 3I, 4I, 5I, D)::Matrix ==
        hvcat((2,2,4), C, C, Matrix(I, 3, 3), Matrix(2I,3,3),
            Matrix(3I, 2, 2), Matrix(4I, 2, 2), Matrix(5I,2,2), D)
    @test (hvcat((2,3,2), B, 2I, C, C, I, 3I, 4I))::Matrix ==
        hvcat((2,2,2), B, Matrix(2I, 3, 3), C, C, Matrix(3I, 3, 3), Matrix(4I, 3, 3))
    @test hvcat((3,2,1), C, C, I, B ,3I, 2I)::Matrix ==
        hvcat((2,2,1), C, C, B, Matrix(3I,3,3), Matrix(2I,6,6))
    @test (hvcat((1,2), A, E, α))::Matrix == hvcat((1,2), A, E, [α]) == hvcat((1,2), A, E, α*I)
    @test (hvcat((2,2), α, E, F, 3I))::Matrix == hvcat((2,2), [α], E, F, Matrix(3I, 3, 3))
    @test (hvcat((2,2), 3I, F, E, α))::Matrix == hvcat((2,2), Matrix(3I, 3, 3), F, E, [α])
end

@testset "Matrix/Array construction from UniformScaling" begin
    I2_33 = [2 0 0; 0 2 0; 0 0 2]
    I2_34 = [2 0 0 0; 0 2 0 0; 0 0 2 0]
    I2_43 = [2 0 0; 0 2 0; 0 0 2; 0 0 0]
    for ArrType in (Matrix, Array)
        @test ArrType(2I, 3, 3)::Matrix{Int} == I2_33
        @test ArrType(2I, 3, 4)::Matrix{Int} == I2_34
        @test ArrType(2I, 4, 3)::Matrix{Int} == I2_43
        @test ArrType(2.0I, 3, 3)::Matrix{Float64} == I2_33
        @test ArrType{Real}(2I, 3, 3)::Matrix{Real} == I2_33
        @test ArrType{Float64}(2I, 3, 3)::Matrix{Float64} == I2_33
    end
end

@testset "Diagonal construction from UniformScaling" begin
    @test Diagonal(2I, 3)::Diagonal{Int} == Matrix(2I, 3, 3)
    @test Diagonal(2.0I, 3)::Diagonal{Float64} == Matrix(2I, 3, 3)
    @test Diagonal{Real}(2I, 3)::Diagonal{Real} == Matrix(2I, 3, 3)
    @test Diagonal{Float64}(2I, 3)::Diagonal{Float64} == Matrix(2I, 3, 3)
end

@testset "equality comparison of matrices with UniformScaling" begin
    # AbstractMatrix methods
    diagI = Diagonal(fill(1, 3))
    rdiagI = view(diagI, 1:2, 1:3)
    bidiag = Bidiagonal(fill(2, 3), fill(2, 2), :U)
    @test diagI  ==  I == diagI  # test isone(I) path / equality
    @test 2diagI !=  I != 2diagI # test isone(I) path / inequality
    @test 0diagI == 0I == 0diagI # test iszero(I) path / equality
    @test 2diagI != 0I != 2diagI # test iszero(I) path / inequality
    @test 2diagI == 2I == 2diagI # test generic path / equality
    @test 0diagI != 2I != 0diagI # test generic path / inequality on diag
    @test bidiag != 2I != bidiag # test generic path / inequality off diag
    @test rdiagI !=  I != rdiagI # test square matrix check
    # StridedMatrix specialization
    denseI = [1 0 0; 0 1 0; 0 0 1]
    rdenseI = [1 0 0 0; 0 1 0 0; 0 0 1 0]
    alltwos = fill(2, (3, 3))
    @test denseI  ==  I == denseI  # test isone(I) path / equality
    @test 2denseI !=  I != 2denseI # test isone(I) path / inequality
    @test 0denseI == 0I == 0denseI # test iszero(I) path / equality
    @test 2denseI != 0I != 2denseI # test iszero(I) path / inequality
    @test 2denseI == 2I == 2denseI # test generic path / equality
    @test 0denseI != 2I != 0denseI # test generic path / inequality on diag
    @test alltwos != 2I != alltwos # test generic path / inequality off diag
    @test rdenseI !=  I != rdenseI # test square matrix check

    # isequal
    @test !isequal(I, I(3))
    @test !isequal(I(1), I)
    @test !isequal([1], I)
    @test isequal(I, 1I)
    @test !isequal(2I, 3I)
end

@testset "operations involving I should preserve eltype" begin
    @test isa(Int8(1) + I, Int8)
    @test isa(Float16(1) + I, Float16)
    @test eltype(Int8(1)I) == Int8
    @test eltype(Float16(1)I) == Float16
    @test eltype(fill(Int8(1), 2, 2)I) == Int8
    @test eltype(fill(Float16(1), 2, 2)I) == Float16
    @test eltype(fill(Int8(1), 2, 2) + I) == Int8
    @test eltype(fill(Float16(1), 2, 2) + I) == Float16
end

@testset "test that UniformScaling is applied correctly for matrices of matrices" begin
    LL = Bidiagonal(fill(0*I, 3), fill(1*I, 2), :L)
    @test (I - LL')\[[0], [0], [1]] == (I - LL)'\[[0], [0], [1]] == fill([1], 3)
end

@testset "broadcasting (#23197)" begin
    @testset "with dense arrays" begin
        A = [1 2; 3 4]
        @test @inferred(I .+ A) == @inferred(A .+ I) == A + I == [2 2; 3 5]
        @test @inferred(I .- A) == I - A == [0 -2; -3 -3]
        @test A .- I == A - I
        @test A .* I == Diagonal([1, 4])
        @test I .* A == Diagonal([1, 4])
        @test A .+ 2I == A + 2I == [3 2; 3 6]
        @test 2I .+ A == [3 2; 3 6]
        @test A .+ 2 .* I == [3 2; 3 6]
        @test I .- ones(4, 4) == Matrix(I, 4, 4) - ones(4, 4) # the motivating example
        @test A .+ I .+ 1 == A .+ (I .+ 1) == [3 3; 4 6]
        @test (I .+ A) .+ I == [3 2; 3 6]
        @test (A .+ I) .* 2 == [4 4; 6 10]
        @test (I .* 2) .+ (I .+ 1) .+ A == [5 3; 4 8]
        @test I .+ A .+ A' == [3 5; 5 9]
        @test I .+ [1 1; 1 1] == [2 1; 1 2]
        @test max.(I, A) == [1 2; 3 4]
        @test (x -> x^2).(A .+ I) == (A + I).^2
        @test (I .+ A) isa Matrix{Int}
        @test (1.0I .+ A) isa Matrix{Float64}
        @test (I .+ [1.0 2.0; 3.0 4.0]) isa Matrix{Float64}
        @test (I .+ complex(A)) isa Matrix{Complex{Int}}
        @test @inferred(I .& trues(2, 2)) == [true false; false true]
        @test @inferred(I .& trues(2, 2)) isa BitMatrix
        # 1×1 matrices
        @test I .+ ones(1, 1) == fill(2.0, 1, 1)
    end

    @testset "along trailing dimensions" begin
        A = ones(2, 2, 3)
        B = I .+ A
        @test size(B) == (2, 2, 3)
        for k in axes(A, 3)
            @test B[:, :, k] == ones(2, 2) + I
        end
        @test I .* A == cat(fill(Matrix(I, 2, 2), 3)...; dims=3)
    end

    @testset "with only UniformScalings and scalars" begin
        @test @inferred(I .+ I) === 2I
        @test @inferred(.-I) === -I
        @test I .- I === 0I
        @test I .* I === I
        @test I .* 2I === 2I
        @test I .+ 2I .- I === 2I
        @test 2 .* I === I .* 2 === 2I
        @test 2 .* I .+ I === 3I
        @test 2 .* I === Broadcast.broadcast(*, 2, I)
        @test I .* 2.0 === 2.0I
        @test I ./ 2 === 0.5I
        @test 2 .\ I === 0.5I
        @test I .^ 2 === I
        @test (2I) .^ 2 === 4I
        @test (2I) .^ 2.0 === 4.0I
        @test sqrt.(4I) === 2.0I
        @test abs.(-I) === UniformScaling(1)
        @test float.(2I) === 2.0I
        @test (x -> 2x).(I) === 2I
        @test (2I) .* Ref(2) === 4I
        @test (2I) .* fill(2) === 4I
        @test Base.literal_pow.(^, 2I, Val(3)) === 8I
        @test identity.(I) === I
        @test 3.0I .* (I .+ I) === 6.0I
        # functions that don't map the off-diagonal zeros to zero
        @test_throws ArgumentError I .+ 1
        @test_throws ArgumentError 1 .+ I
        @test_throws ArgumentError exp.(I)
        @test_throws ArgumentError I .- 2
        @test_throws ArgumentError (I .+ 1) .* 2
        @test (I .+ 1) .* I === 2I # the off-diagonal elements are (0 + 1) * 0 == 0
        # the same expressions are fine once a shape is provided
        @test (I .+ 1) .+ zeros(2, 2) == [2 1; 1 2]
        @test exp.(I) .* ones(2, 2) == [ℯ 1; 1 ℯ]
    end

    @testset "in-place broadcasting" begin
        A = [1 2; 3 4]
        B = copy(A)
        @test (B .+= I) === B
        @test B == A + I
        B = copy(A)
        @test (B .-= 2I) === B
        @test B == A - 2I
        B = copy(A)
        @test (B .*= I) === B
        @test B == Diagonal(A)
        B = zeros(Int, 3, 3)
        @test (B .= I) === B
        @test B == Matrix(I, 3, 3)
        B .= 2I .+ 1
        @test B == [3 1 1; 1 3 1; 1 1 3]
        B .= I .+ B .+ 1
        @test B == [5 2 2; 2 5 2; 2 2 5]
        B = zeros(2, 2, 2)
        B .= I
        @test B == cat(fill(Matrix(I, 2, 2), 2)...; dims=3)
        B = zeros(3, 3)
        @test Broadcast.broadcast!(+, B, I, 1) === B
        @test B == [2 1 1; 1 2 1; 1 1 2]
        B = zeros(3, 3)
        @test @inferred(Broadcast.broadcast!(-, B, I, ones(3, 3))) == I - ones(3, 3)
        # non-square destinations
        @test_throws DimensionMismatch zeros(2, 3) .= I
        @test_throws DimensionMismatch zeros(2, 3) .+= I
        @test_throws DimensionMismatch zeros(2, 3) .= I .+ 1
        @test_throws DimensionMismatch zeros(3) .= I
        @test_throws DimensionMismatch zeros(3) .= 2 .* I
        @test_throws DimensionMismatch zeros(3) .= I .+ I
        @test_throws DimensionMismatch fill(1) .= I
        @test_throws DimensionMismatch fill(1) .= I .+ I
    end

    @testset "shape mismatches" begin
        @test_throws DimensionMismatch I .+ ones(2, 3)
        @test_throws DimensionMismatch ones(2, 3) .+ I
        @test_throws DimensionMismatch I .+ ones(3, 1)
        @test_throws DimensionMismatch I .+ ones(1, 3)
        @test_throws DimensionMismatch I .+ ones(3)
        @test_throws DimensionMismatch ones(3) .+ I
        @test_throws DimensionMismatch I .+ ones(2, 3, 2)
        @test_throws DimensionMismatch I .+ (1, 2)
        @test_throws DimensionMismatch (I .+ 1) .+ ones(2, 3)
        @test_throws DimensionMismatch (I .* 2) .+ ones(3)
        @test_throws DimensionMismatch ones(2, 2) .+ I .+ ones(2, 3)
        @test_throws DimensionMismatch ones(2, 2) .+ I .+ ones(3, 3)
        @test_throws DimensionMismatch Broadcast.broadcast(+, I, ones(3))
    end

    @testset "with offset arrays" begin
        A = OffsetArray(ones(2, 2), 0:1, 0:1)
        B = A .+ I
        @test axes(B) == axes(A)
        @test B == A + I
        @test B[0, 0] == B[1, 1] == 2
        @test B[0, 1] == B[1, 0] == 1
        # same size, but different axes: the diagonal is where the indices are equal
        A = OffsetArray(ones(2, 2), 0:1, 1:2)
        B = A .+ I
        @test axes(B) == axes(A)
        @test B == A + I
        @test B[1, 1] == 2
        @test B[0, 1] == B[0, 2] == B[1, 2] == 1
        @test_throws DimensionMismatch OffsetArray(ones(2, 3), 0:1, 0:2) .+ I
    end

    @testset "with structured matrices" begin
        D = Diagonal([1, 2])
        @test @inferred(D .+ I) == D + I
        @test @inferred(D .+ I) isa Diagonal
        @test @inferred(I .+ D) isa Diagonal
        @test @inferred(D .* I) == Diagonal([1, 2])
        @test (D .* I) isa Diagonal
        @test D .+ I .+ 1 isa Matrix
        @test D .+ I .+ 1 == D + I .+ 1
        @test D .+ 2I .- I == D + I
        @test D ./ I isa Matrix # division by the zeros
        @test isequal(D ./ I, [1.0 NaN; NaN 2.0])
        D2 = copy(D)
        @test (D2 .+= I) === D2
        @test D2 == D + I
        D2 .= D .* I .+ 2I
        @test D2 == Diagonal([3, 4])
        @test_throws ArgumentError D2 .= I .+ 1
        # mixed dense and structured
        @test D .+ I .+ ones(2, 2) == [3 1; 1 4]
        @test (D .+ I .+ ones(2, 2)) isa Matrix
        B = Bidiagonal([1, 2, 3], [4, 5], :U)
        @test @inferred(B .+ I) == B + I
        @test (B .+ I) isa Bidiagonal
        @test @inferred(I .- B) == I - B
        @test (I .- B) isa Bidiagonal
        T = Tridiagonal([1, 2], [3, 4, 5], [6, 7])
        @test @inferred(T .+ I) == T + I
        @test (T .+ I) isa Tridiagonal
        S = SymTridiagonal([1, 2, 3], [4, 5])
        @test @inferred(S .+ I) == S + I
        @test (S .+ I) isa SymTridiagonal
        S2 = copy(S)
        @test (S2 .= I) === S2
        @test S2 == I(3)
        S2 .= S .+ 2I
        @test S2 == S + 2I
        S2 .+= I
        @test S2 == S + 3I
        A = [1 2; 3 4]
        for TT in (UpperTriangular, LowerTriangular, UpperHessenberg, Symmetric, Hermitian)
            M = TT(A)
            @test @inferred(M .+ I) == M + I
            @test (M .+ I) isa TT
            @test M .+ I .+ 1 == M + I .+ 1
            # adding a scalar preserves symmetry, but not the triangular structures
            @test M .+ I .+ 1 isa (TT <: Union{Symmetric,Hermitian} ? TT : Matrix)
        end
        @test (UnitUpperTriangular(A) .+ I) == UnitUpperTriangular(A) + I
        @test (UnitUpperTriangular(A) .+ I) isa UpperTriangular
        @test (UnitLowerTriangular(A) .+ I) == UnitLowerTriangular(A) + I
        @test (UnitLowerTriangular(A) .+ I) isa LowerTriangular
        @test Hermitian(A) .+ im*I == Hermitian(A) + im*I
        @test (Hermitian(A) .+ im*I) isa Matrix
        @test (Hermitian(A) .+ (1+0im)*I) isa Matrix
        @test (Hermitian(complex(A)) .+ I) isa Hermitian
        @test (Symmetric(A) .+ im*I) isa Symmetric
    end

    @testset "lazy Broadcasted objects" begin
        A = [1 2; 3 4]
        bc = Broadcast.broadcasted(+, I, A)
        @test axes(bc) == axes(A)
        @test size(bc) == size(A)
        @test bc[1, 1] == 2 && bc[1, 2] == 2 && bc[2, 1] == 3 && bc[2, 2] == 5
        @test collect(bc) == A + I
        bci = Broadcast.instantiate(bc)
        @test axes(bci) == axes(A)
        @test copy(bci) == A + I
        @test Broadcast.materialize(bc) == A + I
        bc = Broadcast.broadcasted(+, I, 2I)
        @test axes(bc) == ()
        @test copy(bc) === 3I
        @test Broadcast.materialize(bc) === 3I
        @test Broadcast.materialize(I) === I
        @test Broadcast.broadcastable(I) === I
    end

    @testset "with other number types" begin
        A = fill(1//2, 2, 2)
        @test @inferred(A .+ I) == A + I
        @test (A .+ I) isa Matrix{Rational{Int}}
        @test (A .+ (1//3)I) == A + (1//3)I
        J = Quaternion(1.0, 2.0, 3.0, 4.0) * I
        A = fill(Quaternion(1.0, 1.0, 1.0, 1.0), 2, 2)
        @test A .+ J == A + J
        @test J .+ J === Quaternion(2.0, 4.0, 6.0, 8.0) * I
        @test_throws ArgumentError J .+ Quaternion(1.0, 2.0, 3.0, 4.0)
    end
end

@testset "in-place mul! and div! methods" begin
    J = randn()*I
    A = randn(4, 3)
    C = similar(A)
    target_mul = J * A
    target_div = A / J
    @test mul!(C, J, A) == target_mul
    @test mul!(C, A, J) == target_mul
    @test lmul!(J, copyto!(C, A)) == target_mul
    @test rmul!(copyto!(C, A), J) == target_mul
    @test ldiv!(J, copyto!(C, A)) == target_div
    @test ldiv!(C, J, A) == target_div
    @test rdiv!(copyto!(C, A), J) == target_div

    A = randn(4, 3)
    C = randn!(similar(A))
    alpha = randn()
    beta = randn()
    target = J * A * alpha + C * beta
    @test mul!(copy(C), J, A, alpha, beta) ≈ target
    @test mul!(copy(C), A, J, alpha, beta) ≈ target

    a = randn()
    C = randn(3, 3)
    target_5mul = a*alpha*J + beta*C
    @test mul!(copy(C), a, J, alpha, beta) ≈ target_5mul
    @test mul!(copy(C), J, a, alpha, beta) ≈ target_5mul
    target_5mul = beta*C # alpha = 0
    @test mul!(copy(C), a, J, 0, beta) ≈ target_5mul
    target_5mul = a*alpha*Matrix(J, 3, 3) # beta = 0
    @test mul!(copy(C), a, J, alpha, 0) ≈ target_5mul

end

@testset "Construct Diagonal from UniformScaling" begin
    @test size(I(3)) === (3,3)
    @test I(3) isa Diagonal
    @test I(3) == [1 0 0; 0 1 0; 0 0 1]
end

@testset "dot" begin
    A = randn(3, 3)
    λ = randn()
    J = UniformScaling(λ)
    @test dot(A, J) ≈ dot(J, A)
    @test dot(A, J) ≈ tr(A' * J)

    A = rand(ComplexF64, 3, 3)
    λ = randn() + im * randn()
    J = UniformScaling(λ)
    @test dot(A, J) ≈ conj(dot(J, A))
    @test dot(A, J) ≈ tr(A' * J)
end

@testset "generalized dot" begin
    x = rand(-10:10, 3)
    y = rand(-10:10, 3)
    λ = rand(-10:10)
    J = UniformScaling(λ)
    @test dot(x, J, y) == λ*dot(x, y)
    λ = Quaternion(0.44567, 0.755871, 0.882548, 0.423612)
    x, y = Quaternion(rand(4)...), Quaternion(rand(4)...)
    @test dot([x], λ*I, [y]) ≈ dot(x, λ, y) ≈ dot(x, λ*y)
end

@testset "Factorization solutions" begin
    J = complex(randn(),randn()) * I
    qrp = A -> qr(A, ColumnNorm())

    # thin matrices
    X = randn(3,2)
    Z = pinv(X)
    for fac in (qr,qrp,svd)
        F = fac(X)
        @test @inferred(F \ I) ≈ Z
        @test @inferred(F \ J) ≈ Z * J
    end

    # square matrices
    X = randn(3,3)
    X = X'X + rand()I # make positive definite for cholesky
    Z = pinv(X)
    for fac in (bunchkaufman,cholesky,lu,qr,qrp,svd)
        F = fac(X)
        @test @inferred(F \ I) ≈ Z
        @test @inferred(F \ J) ≈ Z * J
    end

    # fat matrices - only rank-revealing variants
    X = randn(2,3)
    Z = pinv(X)
    for fac in (qrp,svd)
        F = fac(X)
        @test @inferred(F \ I) ≈ Z
        @test @inferred(F \ J) ≈ Z * J
    end
end

@testset "offset arrays" begin
    A = OffsetArray(zeros(4,4), -1:2, 0:3)
    @test sum(I + A) ≈ 3.0
    @test sum(A + I) ≈ 3.0
    @test sum(I - A) ≈ 3.0
    @test sum(A - I) ≈ -3.0
end

@testset "type promotion when dividing UniformScaling by matrix" begin
    A = randn(5,5)
    cA = complex(A)
    J = (5+2im)*I
    @test J/A ≈ J/cA
    @test A\J ≈ cA\J
end

@testset "block matrix equality" begin
    A = Diagonal(fill(I(2), 4))
    @test isone(A)
    @test A != I
    @test 0*A != 0*I
end

end # module TestUniformscaling
