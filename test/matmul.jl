# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestMatmul

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Base: rtoldefault
using Test, LinearAlgebra, Random
using LinearAlgebra: mul!, Symmetric, Hermitian

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.SizedArrays

## Test Julia fallbacks to BLAS routines

mul_wrappers = [
    m -> m,
    m -> Symmetric(m, :U),
    m -> Symmetric(m, :L),
    m -> Hermitian(m, :U),
    m -> Hermitian(m, :L),
    m -> adjoint(m),
    m -> transpose(m)]

@testset "wrap" begin
    f(A) = LinearAlgebra.wrap(A, 'N')
    A = ones(1,1)
    @test @inferred(f(A)) === A
    g(A) = LinearAlgebra.wrap(A, 'T')
    @test @inferred(g(A)) === transpose(A)
    # https://github.com/JuliaLang/julia/issues/52202
    @test Base.infer_return_type((Vector{Float64},)) do v
        LinearAlgebra.wrap(v, 'N')
    end == Vector{Float64}
    h(A) = LinearAlgebra.wrap(LinearAlgebra._unwrap(A), LinearAlgebra.wrapper_char(A))
    @test @inferred(h(transpose(A))) === transpose(A)
    @test @inferred(h(adjoint(A))) === transpose(A)

    M = rand(2,2)
    for S in (Symmetric(M), Hermitian(M))
        @test @inferred((A -> LinearAlgebra.wrap(parent(A), LinearAlgebra.wrapper_char(A)))(S)) === Symmetric(M)
    end
    M = rand(ComplexF64,2,2)
    for S in (Symmetric(M), Hermitian(M))
        @test @inferred((A -> LinearAlgebra.wrap(parent(A), LinearAlgebra.wrapper_char(A)))(S)) === S
    end

    @testset "WrapperChar" begin
        @test LinearAlgebra.WrapperChar('c') == 'c'
        @test LinearAlgebra.WrapperChar('C') == 'C'
        @testset "constant propagation in uppercase/lowercase" begin
            v = @inferred (() -> Val(uppercase(LinearAlgebra.WrapperChar('C'))))()
            @test v isa Val{'C'}
            v = @inferred (() -> Val(uppercase(LinearAlgebra.WrapperChar('s'))))()
            @test v isa Val{'S'}
            v = @inferred (() -> Val(lowercase(LinearAlgebra.WrapperChar('C'))))()
            @test v isa Val{'c'}
            v = @inferred (() -> Val(lowercase(LinearAlgebra.WrapperChar('s'))))()
            @test v isa Val{'s'}
        end
    end
end

@testset "matrices with zero dimensions" begin
    for (dimsA, dimsB, dimsC) in (
        ((0, 5), (5, 3), (0, 3)),
        ((3, 5), (5, 0), (3, 0)),
        ((3, 0), (0, 4), (3, 4)),
        ((0, 5), (5, 0), (0, 0)),
        ((0, 0), (0, 4), (0, 4)),
        ((3, 0), (0, 0), (3, 0)),
        ((0, 0), (0, 0), (0, 0)))
        @test Matrix{Float64}(undef, dimsA) * Matrix{Float64}(undef, dimsB) == zeros(dimsC)
    end
    @test Matrix{Float64}(undef, 5, 0) |> t -> t't == zeros(0, 0)
    @test Matrix{Float64}(undef, 5, 0) |> t -> t * t' == zeros(5, 5)
    @test Matrix{ComplexF64}(undef, 5, 0) |> t -> t't == zeros(0, 0)
    @test Matrix{ComplexF64}(undef, 5, 0) |> t -> t * t' == zeros(5, 5)
end
@testset "2x2 matmul" begin
    AA = [1 2; 3 4]
    BB = [5 6; 7 8]
    AAi = AA + (0.5 * im) .* BB
    BBi = BB + (2.5 * im) .* AA[[2, 1], [2, 1]]
    for A in (copy(AA), view(AA, 1:2, 1:2)), B in (copy(BB), view(BB, 1:2, 1:2))
        @test A * B == [19 22; 43 50]
        @test *(transpose(A), B) == [26 30; 38 44]
        @test *(A, transpose(B)) == [17 23; 39 53]
        @test *(transpose(A), transpose(B)) == [23 31; 34 46]
    end
    for Ai in (copy(AAi), view(AAi, 1:2, 1:2)), Bi in (copy(BBi), view(BBi, 1:2, 1:2))
        @test Ai * Bi == [-21+53.5im -4.25+51.5im; -12+95.5im 13.75+85.5im]
        @test *(adjoint(Ai), Bi) == [68.5-12im 57.5-28im; 88-3im 76.5-25im]
        @test *(Ai, adjoint(Bi)) == [64.5+5.5im 43+31.5im; 104-18.5im 80.5+31.5im]
        @test *(adjoint(Ai), adjoint(Bi)) == [-28.25-66im 9.75-58im; -26-89im 21-73im]
        @test_throws DimensionMismatch [1 2; 0 0; 0 0] * [1 2]
    end
    for wrapper_a in mul_wrappers, wrapper_b in mul_wrappers
        @test wrapper_a(AA) * wrapper_b(BB) == Array(wrapper_a(AA)) * Array(wrapper_b(BB))
    end
    @test_throws DimensionMismatch mul!(Matrix{Float64}(undef, 3, 3), AA, BB)
end
@testset "3x3 matmul" begin
    AA = [1 2 3; 4 5 6; 7 8 9] .- 5
    BB = [1 0 5; 6 -10 3; 2 -4 -1]
    AAi = AA + (0.5 * im) .* BB
    BBi = BB + (2.5 * im) .* AA[[2, 1, 3], [2, 3, 1]]
    for A in (copy(AA), view(AA, 1:3, 1:3)), B in (copy(BB), view(BB, 1:3, 1:3))
        @test A * B == [-26 38 -27; 1 -4 -6; 28 -46 15]
        @test *(adjoint(A), B) == [-6 2 -25; 3 -12 -18; 12 -26 -11]
        @test *(A, adjoint(B)) == [-14 0 6; 4 -3 -3; 22 -6 -12]
        @test *(adjoint(A), adjoint(B)) == [6 -8 -6; 12 -9 -9; 18 -10 -12]
    end
    for Ai in (copy(AAi), view(AAi, 1:3, 1:3)), Bi in (copy(BBi), view(BBi, 1:3, 1:3))
        @test Ai * Bi == [-44.75+13im 11.75-25im -38.25+30im; -47.75-16.5im -51.5+51.5im -56+6im; 16.75-4.5im -53.5+52im -15.5im]
        @test *(adjoint(Ai), Bi) == [-21+2im -1.75+49im -51.25+19.5im; 25.5+56.5im -7-35.5im 22+35.5im; -3+12im -32.25+43im -34.75-2.5im]
        @test *(Ai, adjoint(Bi)) == [-20.25+15.5im -28.75-54.5im 22.25+68.5im; -12.25+13im -15.5+75im -23+27im; 18.25+im 1.5+94.5im -27-54.5im]
        @test *(adjoint(Ai), adjoint(Bi)) == [1+2im 20.75+9im -44.75+42im; 19.5+17.5im -54-36.5im 51-14.5im; 13+7.5im 11.25+31.5im -43.25-14.5im]
        @test_throws DimensionMismatch [1 2 3; 0 0 0; 0 0 0] * [1 2 3]
    end
    for wrapper_a in mul_wrappers, wrapper_b in mul_wrappers
        @test wrapper_a(AA) * wrapper_b(BB) == Array(wrapper_a(AA)) * Array(wrapper_b(BB))
    end
    @test_throws DimensionMismatch mul!(Matrix{Float64}(undef, 4, 4), AA, BB)
end

# Generic AbstractArrays
module MyArray15367
using Test, Random

struct MyArray{T,N} <: AbstractArray{T,N}
    data::Array{T,N}
end

Base.size(A::MyArray) = size(A.data)
Base.getindex(A::MyArray, indices...) = A.data[indices...]

A = MyArray(rand(4, 5))
b = rand(5)
@test A * b ≈ A.data * b
end

@testset "Generic integer matrix multiplication" begin
    AA = [1 2 3; 4 5 6] .- 3
    BB = [2 -2; 3 -5; -4 7]
    for A in (copy(AA), view(AA, 1:2, 1:3)), B in (copy(BB), view(BB, 1:3, 1:2))
        @test A * B == [-7 9; -4 9]
        @test *(transpose(A), transpose(B)) == [-6 -11 15; -6 -13 18; -6 -15 21]
    end
    AA = fill(1, 2, 100)
    BB = fill(1, 100, 3)
    for A in (copy(AA), view(AA, 1:2, 1:100)), B in (copy(BB), view(BB, 1:100, 1:3))
        @test A * B == [100 100 100; 100 100 100]
    end
    AA = rand(1:20, 5, 5) .- 10
    BB = rand(1:20, 5, 5) .- 10
    CC = Matrix{Int}(undef, size(AA, 1), size(BB, 2))
    for A in (copy(AA), view(AA, 1:5, 1:5)), B in (copy(BB), view(BB, 1:5, 1:5)), C in (copy(CC), view(CC, 1:5, 1:5))
        @test *(transpose(A), B) == A' * B
        @test *(A, transpose(B)) == A * B'
        # Preallocated
        @test mul!(C, A, B) == A * B
        @test mul!(C, transpose(A), B) == A' * B
        @test mul!(C, A, transpose(B)) == A * B'
        @test mul!(C, transpose(A), transpose(B)) == A' * B'
        @test mul!(C, adjoint(A), transpose(B)) == A' * transpose(B)

        # Inplace multiply-add
        α = rand(-10:10)
        β = rand(-10:10)
        rand!(C, -10:10)
        βC = β * C
        _C0 = copy(C)
        C0() = (C .= _C0; C)  # reset C but don't change the container type
        @test mul!(C0(), A, B, α, β) == α * A * B .+ βC
        @test mul!(C0(), transpose(A), B, α, β) == α * A' * B .+ βC
        @test mul!(C0(), A, transpose(B), α, β) == α * A * B' .+ βC
        @test mul!(C0(), transpose(A), transpose(B), α, β) == α * A' * B' .+ βC
        @test mul!(C0(), adjoint(A), transpose(B), α, β) == α * A' * transpose(B) .+ βC

        #test DimensionMismatch for generic_matmatmul
        @test_throws DimensionMismatch mul!(C, adjoint(A), transpose(fill(1, 4, 4)))
        @test_throws DimensionMismatch mul!(C, adjoint(fill(1, 4, 4)), transpose(B))
    end
    vv = [1, 2]
    CC = Matrix{Int}(undef, 2, 2)
    for v in (copy(vv), view(vv, 1:2)), C in (copy(CC), view(CC, 1:2, 1:2))
        @test @inferred(mul!(C, v, adjoint(v))) == [1 2; 2 4]

        C .= [1 0; 0 1]
        @test @inferred(mul!(C, v, adjoint(v), 2, 3)) == [5 4; 4 11]
    end
end

@testset "generic_matvecmul" begin
    AA = rand(5, 5)
    BB = rand(5)
    for A in (copy(AA), view(AA, 1:5, 1:5)), B in (copy(BB), view(BB, 1:5))
        @test_throws DimensionMismatch LinearAlgebra.generic_matvecmul!(zeros(6), 'N', A, B)
        @test_throws DimensionMismatch LinearAlgebra.generic_matvecmul!(B, 'N', A, zeros(6))
    end
    vv = [1, 2, 3]
    CC = Matrix{Int}(undef, 3, 3)
    for v in (copy(vv), view(vv, 1:3)), C in (copy(CC), view(CC, 1:3, 1:3))
        @test mul!(C, v, transpose(v)) == v * v'
        C .= C0 = rand(-10:10, size(C))
        @test mul!(C, v, transpose(v), 2, 3) == 2v * v' .+ 3C0
    end
    vvf = map(Float64, vv)
    CC = Matrix{Float64}(undef, 3, 3)
    for vf in (copy(vvf), view(vvf, 1:3)), C in (copy(CC), view(CC, 1:3, 1:3))
        @test mul!(C, vf, transpose(vf)) == vf * vf'
        C .= C0 = rand(eltype(C), size(C))
        @test mul!(C, vf, transpose(vf), 2, 3) ≈ 2vf * vf' .+ 3C0
    end

    @testset "zero stride" begin
        for AAv in (view(AA, StepRangeLen(2,0,size(AA,1)), :),
                    view(AA, StepRangeLen.(2,0,size(AA))...),
                    view(complex.(AA, AA), StepRangeLen.(2,0,size(AA))...),)
            for BB2 in (BB, complex.(BB, BB))
                C = AAv * BB2
                @test allequal(C)
                @test C ≈ Array(AAv) * BB2
            end
        end
    end
end

@testset "generic_matvecmul for vectors of vectors" begin
    @testset "matrix of scalars" begin
        u = [[1, 2], [3, 4]]
        A = [1 2; 3 4]
        v = [[0, 0], [0, 0]]
        Au = [[7, 10], [15, 22]]
        @test A * u == Au
        mul!(v, A, u)
        @test v == Au
        mul!(v, A, u, 2, -1)
        @test v == Au
    end

    @testset "matrix of matrices" begin
        u = [[1, 2], [3, 4]]
        A = Matrix{Matrix{Int}}(undef, 2, 2)
        A[1, 1] = [1 2; 3 4]
        A[1, 2] = [5 6; 7 8]
        A[2, 1] = [9 10; 11 12]
        A[2, 2] = [13 14; 15 16]
        v = [[0, 0], [0, 0]]
        Au = [[44, 64], [124, 144]]
        @test A * u == Au
        mul!(v, A, u)
        @test v == Au
        mul!(v, A, u, 2, -1)
        @test v == Au
    end
end

@testset "generic_matvecmul for vectors of matrices" begin
    x = [1 2 3; 4 5 6]
    A = reshape([x,2x,3x,4x],2,2)
    b = [x, 2x]
    for f in (adjoint, transpose)
        c = f(A) * b
        for i in eachindex(c)
            @test c[i] == sum(f(A)[i, j] * b[j] for j in eachindex(b))
        end
    end
end

@testset "generic_matmatmul for matrices of vectors" begin
    B = Matrix{Vector{Int}}(undef, 2, 2)
    B[1, 1] = [1, 2]
    B[2, 1] = [3, 4]
    B[1, 2] = [5, 6]
    B[2, 2] = [7, 8]
    A = [1 2; 3 4]
    C = Matrix{Vector{Int}}(undef, 2, 2)
    AB = Matrix{Vector{Int}}(undef, 2, 2)
    AB[1, 1] = [7, 10]
    AB[2, 1] = [15, 22]
    AB[1, 2] = [19, 22]
    AB[2, 2] = [43, 50]
    @test A * B == AB
    mul!(C, A, B)
    @test C == AB
    mul!(C, A, B, 2, -1)
    @test C == AB
    LinearAlgebra.generic_matmatmul!(C, 'N', 'N', A, B, LinearAlgebra.MulAddMul(2, -1))
    @test C == AB
end

@testset "fallbacks & such for BlasFloats" begin
    AA = rand(Float64, 6, 6)
    BB = rand(Float64, 6, 6)
    CC = zeros(Float64, 6, 6)
    for A in (copy(AA), view(AA, 1:6, 1:6)), B in (copy(BB), view(BB, 1:6, 1:6)), C in (copy(CC), view(CC, 1:6, 1:6))
        @test mul!(C, transpose(A), transpose(B)) == transpose(A) * transpose(B)
        @test mul!(C, A, adjoint(B)) == A * transpose(B)
        @test mul!(C, adjoint(A), B) == transpose(A) * B

        # Inplace multiply-add
        α = rand(Float64)
        β = rand(Float64)
        rand!(C)
        βC = β * C
        _C0 = copy(C)
        C0() = (C .= _C0; C)  # reset C but don't change the container type
        @test mul!(C0(), transpose(A), transpose(B), α, β) ≈ α * transpose(A) * transpose(B) .+ βC
        @test mul!(C0(), A, adjoint(B), α, β) ≈ α * A * transpose(B) .+ βC
        @test mul!(C0(), adjoint(A), B, α, β) ≈ α * transpose(A) * B .+ βC
    end
end

@testset "allocations in BLAS-mul" begin
    for n in (2, 3, 6)
        A = rand(Float64, n, n)
        B = rand(Float64, n, n)
        C = zeros(Float64, n, n)
        # gemm
        for t in (identity, adjoint, transpose)
            At = t(A)
            Bt = t(B)
            mul!(C, At, B)
            @test 0 == @allocations mul!(C, At, B)
            mul!(C, A, Bt)
            @test 0 == @allocations mul!(C, A, Bt)
            mul!(C, At, Bt)
            @test 0 == @allocations mul!(C, At, Bt)
        end
        # syrk/herk
        mul!(C, transpose(A), A)
        mul!(C, adjoint(A), A)
        mul!(C, A, transpose(A))
        mul!(C, A, adjoint(A))
        @test 0 == @allocations mul!(C, transpose(A), A)
        @test 0 == @allocations mul!(C, adjoint(A), A)
        @test 0 == @allocations mul!(C, A, transpose(A))
        @test 0 == @allocations mul!(C, A, adjoint(A))
        # complex times real
        Cc = complex(C)
        Ac = complex(A)
        for t in (identity, adjoint, transpose)
            Bt = t(B)
            mul!(Cc, Ac, Bt)
            @test 0 == @allocations mul!(Cc, Ac, Bt)
        end
    end
end

@testset "mixed Blas-non-Blas matmul" begin
    AA = rand(-10:10, 6, 6)
    BB = ones(Float64, 6, 6)
    CC = zeros(Float64, 6, 6)
    for A in (copy(AA), view(AA, 1:6, 1:6)), B in (copy(BB), view(BB, 1:6, 1:6)), C in (copy(CC), view(CC, 1:6, 1:6))
        @test mul!(C, A, B) == A * B
        @test mul!(C, transpose(A), transpose(B)) == transpose(A) * transpose(B)
        @test mul!(C, A, adjoint(B)) == A * transpose(B)
        @test mul!(C, adjoint(A), B) == transpose(A) * B
    end
end

@testset "allocations in mixed Blas-non-Blas matmul" begin
    for n in (2, 3, 6)
        A = rand(-10:10, n, n)
        B = ones(Float64, n, n)
        C = zeros(Float64, n, n)
        mul!(C, A, B)
        mul!(C, A, transpose(B))
        mul!(C, adjoint(A), B)
        @test 0 == @allocations mul!(C, A, B)
        @test 0 == @allocations mul!(C, A, transpose(B))
        @test 0 == @allocations mul!(C, adjoint(A), B)
    end
end

@testset "matrix algebra with subarrays of floats (stride != 1)" begin
    A = reshape(map(Float64, 1:20), 5, 4)
    Aref = A[1:2:end, 1:2:end]
    Asub = view(A, 1:2:5, 1:2:4)
    b = [1.2, -2.5]
    @test (Aref * b) == (Asub * b)
    @test *(transpose(Asub), Asub) == *(transpose(Aref), Aref)
    @test *(Asub, transpose(Asub)) == *(Aref, transpose(Aref))
    Ai = A .+ im
    Aref = Ai[1:2:end, 1:2:end]
    Asub = view(Ai, 1:2:5, 1:2:4)
    @test *(adjoint(Asub), Asub) == *(adjoint(Aref), Aref)
    @test *(Asub, adjoint(Asub)) == *(Aref, adjoint(Aref))
end

@testset "matrix x matrix with negative stride" begin
    M = reshape(map(Float64, 1:77), 7, 11)
    N = reshape(map(Float64, 1:63), 9, 7)
    U = view(M, 7:-1:1, 11:-2:1)
    V = view(N, 7:-1:2, 7:-1:1)
    @test U * V ≈ Matrix(U) * Matrix(V)
end

@testset "dot product of subarrays of vectors (floats, negative stride, issue #37767)" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        a = Vector{T}(3:2:7)
        b = Vector{T}(1:10)
        v = view(b, 7:-2:3)
        @test dot(a, Vector(v)) ≈ 67.0
        @test dot(a, v) ≈ 67.0
        @test dot(v, a) ≈ 67.0
        @test dot(Vector(v), Vector(v)) ≈ 83.0
        @test dot(v, v) ≈ 83.0
    end
end

@testset "dot product of stride-vector like input" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        a = randn(T, 10)
        b = view(a, 1:10)
        c = reshape(b, 5, 2)
        d = view(c, :, 1:2)
        r = sum(abs2, a)
        for x in (a,b,c,d), y in (a,b,c,d)
            @test dot(x, y) ≈ r
        end
    end
end

@testset "Complex matrix x real MatOrVec etc (issue #29224)" for T in (Float32, Float64)
    A0 = randn(complex(T), 10, 10)
    B0 = randn(T, 10, 10)
    @testset "Combination Mat{$(complex(T))} Mat{$T}" for Bax1 in (1:5, 2:2:10), Bax2 in (1:5, 2:2:10)
        B = view(A0, Bax1, Bax2)
        tB = transpose(B)
        Bd, tBd = copy(B), copy(tB)
        for Aax1 in (1:5, 2:2:10, (:)), Aax2 in (1:5, 2:2:10)
            A = view(A0, Aax1, Aax2)
            AB_correct = copy(A) * Bd
            AtB_correct = copy(A) * tBd
            @test A*Bd ≈ AB_correct # view times matrix
            @test A*B ≈ AB_correct # view times view
            @test A*tBd ≈ AtB_correct # view times transposed matrix
            @test A*tB ≈ AtB_correct # view times transposed view
        end
    end
    x = randn(T, 10)
    y0 = similar(A0, 20)
    @testset "Combination Mat{$(complex(T))} Vec{$T}" for Aax1 in (1:5, 2:2:10, (:)), Aax2 in (1:5, 2:2:10)
        A = view(A0, Aax1, Aax2)
        Ad = copy(A)
        for indx in (1:5, 1:2:10, 6:-1:2)
            vx = view(x, indx)
            dx = x[indx]
            Ax_correct = Ad*dx
            @test A*vx ≈ A*dx ≈ Ad*vx ≈ Ax_correct # view/matrix times view/vector
            for indy in (1:2:2size(A,1), size(A,1):-1:1)
                y = view(y0, indy)
                @test mul!(y, A, vx) ≈ mul!(y, A, dx) ≈ mul!(y, Ad, vx) ≈
                    mul!(y, Ad, dx) ≈ Ax_correct   # test for uncontiguous dest
            end
        end
    end
end

@testset "real matrix x complex vec" begin
    _matmulres(M, v) = [mapreduce(*, +, row, v) for row in eachrow(M)]
    testmatmul(M, v) = @test M * v ≈ _matmulres(M, v)

    @testset for T in (Float32, Float64), n = (4, 5)
        M1 = reshape(Vector{T}(1:n^2), n, n)
        M2 = reinterpret(reshape, T, [Tuple(T(i + j) for j in 1:n) for i in 1:n])
        v = convert(Vector{Complex{T}}, (1:n) .+ im .* (4 .+ (1:n)))

        for M in (M1, M2)
            M_view_cont = @view M[:, :]
            v_view_cont = @view v[:]
            for _M in (M, M_view_cont), _v in (v, v_view_cont)
                testmatmul(_M, _v)
            end

            # construct a view with strides(M, 1) == 1 and strides(M, 2) != 1
            ax_noncont = 1:2:n
            n1 = length(ax_noncont)
            M_view_noncont = @view M[1:n1, ax_noncont]
            v_view_noncont = @view v[ax_noncont]
            testmatmul(M_view_noncont, v_view_noncont)

            @testset for op in (transpose, adjoint)
                for _M in (M, M_view_cont), _v in (v, v_view_cont)
                    _M2 = op(_M)
                    testmatmul(_M2, _v)
                end
                _M2 = op(M_view_noncont)
                testmatmul(_M2, v_view_noncont)
            end
        end
    end
end

@testset "matrix x vector with negative lda or 0 stride" for T in (Float32, Float64)
    for TA in (T, complex(T)), TB in (T, complex(T))
        A = view(randn(TA, 10, 10), 1:10, 10:-1:1) # negative lda
        v = view([randn(TB)], 1 .+ 0(1:10)) # 0 stride
        Ad, vd = copy(A), copy(v)
        @test Ad * vd ≈ A * vd ≈ Ad * v ≈ A * v
    end
end

@testset "issue #15286" begin
    A = reshape(map(Float64, 1:20), 5, 4)
    C = zeros(8, 8)
    sC = view(C, 1:2:8, 1:2:8)
    B = reshape(map(Float64, -9:10), 5, 4)
    @test mul!(sC, transpose(A), A) == A' * A
    @test mul!(sC, transpose(A), B) == A' * B

    Aim = A .- im
    C = zeros(ComplexF64, 8, 8)
    sC = view(C, 1:2:8, 1:2:8)
    B = reshape(map(Float64, -9:10), 5, 4) .+ im
    @test mul!(sC, adjoint(Aim), Aim) == Aim' * Aim
    @test mul!(sC, adjoint(Aim), B) == Aim' * B
end

@testset "syrk & herk" begin
    AA = reshape(1:1503, 501, 3) .- 750.0
    res = Float64[135228751 9979252 -115270247; 9979252 10481254 10983256; -115270247 10983256 137236759]
    for A in (copy(AA), view(AA, 1:501, 1:3))
        @test *(transpose(A), A) == res
        @test *(adjoint(A), transpose(copy(A'))) == res
    end
    cutoff = 501
    A = reshape(1:6*cutoff, 2 * cutoff, 3) .- (6 * cutoff) / 2
    Asub = view(A, 1:2:2*cutoff, 1:3)
    Aref = A[1:2:2*cutoff, 1:3]
    @test *(transpose(Asub), Asub) == *(transpose(Aref), Aref)
    Ai = A .- im
    Asub = view(Ai, 1:2:2*cutoff, 1:3)
    Aref = Ai[1:2:2*cutoff, 1:3]
    @test *(adjoint(Asub), Asub) == *(adjoint(Aref), Aref)

    A5x5, A6x5 = Matrix{Float64}.(undef, ((5, 5), (6, 5)))
    @test_throws DimensionMismatch LinearAlgebra.syrk_wrapper!(A5x5, 'N', A6x5)
    @test_throws DimensionMismatch LinearAlgebra.herk_wrapper!(complex(A5x5), 'N', complex(A6x5))
end

@testset "5-arg syrk! & herk!" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64), A in (randn(T, 5), randn(T, 5, 5))
        B = A' * A
        C = B isa Number ? [B;;] : Matrix(Hermitian(B))
        @test mul!(copy(C), A', A, true, 2) ≈ 3C
        D = Matrix(Hermitian(A * A'))
        @test mul!(copy(D), A, A', true, 3) ≈ 4D
        if T <: Complex
            @test mul!(2C, A', A, im, 2) ≈ (4 + im) * C
            @test mul!(2D, A, A', im, 3) ≈ (6 + im) * D
        end
    end
end

@testset "matmul for types w/o sizeof (issue #1282)" begin
    AA = fill(complex(1, 1), 10, 10)
    for A in (copy(AA), view(AA, 1:10, 1:10))
        A2 = A^2
        @test A2[1, 1] == 20im
    end
end

@testset "mul! (scaling)" begin
    A5x5, b5, C5x6 = Array{Float64}.(undef, ((5, 5), 5, (5, 6)))
    for A in (A5x5, view(A5x5, :, :)), b in (b5, view(b5, :)), C in (C5x6, view(C5x6, :, :))
        @test_throws DimensionMismatch mul!(A, Diagonal(b), C)
    end
end

@testset "muladd" begin
    A23 = reshape(1:6, 2, 3) .+ 0
    B34 = reshape(1:12, 3, 4) .+ im
    u2 = [10, 20]
    v3 = [3, 5, 7] .+ im
    w4 = [11, 13, 17, 19im]

    @testset "matrix-matrix" begin
        @test muladd(A23, B34, 0) == A23 * B34
        @test muladd(A23, B34, 100) == A23 * B34 .+ 100
        @test muladd(A23, B34, u2) == A23 * B34 .+ u2
        @test muladd(A23, B34, w4') == A23 * B34 .+ w4'
        @test_throws DimensionMismatch muladd(B34, A23, 1)
        @test muladd(ones(1, 3), ones(3, 4), ones(1, 4)) == fill(4.0, 1, 4)
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3, 4), ones(9, 4))

        # broadcasting fallback method allows trailing dims
        @test muladd(A23, B34, ones(2, 4, 1)) == A23 * B34 + ones(2, 4, 1)
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3, 4), ones(9, 4, 1))
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3, 4), ones(1, 4, 9))
        # and catches z::Array{T,0}
        @test muladd(A23, B34, fill(0)) == A23 * B34
    end
    @testset "matrix-vector" begin
        @test muladd(A23, v3, 0) == A23 * v3
        @test muladd(A23, v3, 100) == A23 * v3 .+ 100
        @test muladd(A23, v3, u2) == A23 * v3 .+ u2
        @test muladd(A23, v3, im) isa Vector{Complex{Int}}
        @test muladd(ones(1, 3), ones(3), ones(1)) == [4]
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3), ones(7))

        # fallback
        @test muladd(A23, v3, ones(2, 1, 1)) == A23 * v3 + ones(2, 1, 1)
        @test_throws DimensionMismatch muladd(A23, v3, ones(2, 2))
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3), ones(7, 1))
        @test_throws DimensionMismatch muladd(ones(1, 3), ones(3), ones(1, 7))
        @test muladd(A23, v3, fill(0)) == A23 * v3
    end
    @testset "adjoint-matrix" begin
        @test muladd(v3', B34, 0) isa Adjoint
        @test muladd(v3', B34, 2im) == v3' * B34 .+ 2im
        @test muladd(v3', B34, w4') == v3' * B34 .+ w4'

        # via fallback
        @test muladd(v3', B34, ones(1, 4)) == (B34' * v3 + ones(4, 1))'
        @test_throws DimensionMismatch muladd(v3', B34, ones(7, 4))
        @test_throws DimensionMismatch muladd(v3', B34, ones(1, 4, 7))
        @test muladd(v3', B34, fill(0)) == v3' * B34 # does not make an Adjoint
    end
    @testset "vector-adjoint" begin
        @test muladd(u2, v3', 0) isa Matrix
        @test muladd(u2, v3', 99) == u2 * v3' .+ 99
        @test muladd(u2, v3', A23) == u2 * v3' .+ A23

        @test muladd(u2, v3', ones(2, 3, 1)) == u2 * v3' + ones(2, 3, 1)
        @test_throws DimensionMismatch muladd(u2, v3', ones(2, 3, 4))
        @test_throws DimensionMismatch muladd([1], v3', ones(7, 3))
        @test muladd(u2, v3', fill(0)) == u2 * v3'
    end
    @testset "dot" begin # all use muladd(::Any, ::Any, ::Any)
        @test muladd(u2', u2, 0) isa Number
        @test muladd(v3', v3, im) == dot(v3, v3) + im
        @test muladd(u2', u2, [1]) == [dot(u2, u2) + 1]
        @test_throws DimensionMismatch muladd(u2', u2, [1, 1]) == [dot(u2, u2) + 1]
        @test muladd(u2', u2, fill(0)) == dot(u2, u2)
    end
    @testset "arrays of arrays" begin
        vofm = [rand(1:9, 2, 2) for _ in 1:3]
        Mofm = [rand(1:9, 2, 2) for _ in 1:3, _ in 1:3]

        @test muladd(vofm', vofm, vofm[1]) == vofm' * vofm .+ vofm[1] # inner
        @test muladd(vofm, vofm', Mofm) == vofm * vofm' .+ Mofm       # outer
        @test muladd(vofm', Mofm, vofm') == vofm' * Mofm .+ vofm'     # bra-mat
        @test muladd(Mofm, Mofm, vofm) == Mofm * Mofm .+ vofm         # mat-mat
        @test muladd(Mofm, vofm, vofm) == Mofm * vofm .+ vofm         # mat-vec
    end
end

@testset "muladd & structured matrices" begin
    A33 = reshape(1:9, 3, 3) .+ im
    v3 = [3, 5, 7im]

    # no special treatment
    @test muladd(Symmetric(A33), Symmetric(A33), 1) == Symmetric(A33) * Symmetric(A33) .+ 1
    @test muladd(Hermitian(A33), Hermitian(A33), v3) == Hermitian(A33) * Hermitian(A33) .+ v3
    @test muladd(adjoint(A33), transpose(A33), A33) == A33' * transpose(A33) .+ A33

    u1 = muladd(UpperTriangular(A33), UpperTriangular(A33), Diagonal(v3))
    @test u1 isa UpperTriangular
    @test u1 == UpperTriangular(A33) * UpperTriangular(A33) + Diagonal(v3)

    # diagonal
    @test muladd(Diagonal(v3), Diagonal(A33), Diagonal(v3)).diag == ([1, 5, 9] .+ im .+ 1) .* v3

    # uniformscaling
    @test muladd(Diagonal(v3), I, I).diag == v3 .+ 1
    @test muladd(2 * I, 3 * I, I).λ == 7
    @test muladd(A33, A33', I) == A33 * A33' + I

    # https://github.com/JuliaLang/julia/issues/38426
    @test @evalpoly(A33, 1.0 * I, 1.0 * I) == I + A33
    @test @evalpoly(A33, 1.0 * I, 1.0 * I, 1.0 * I) == I + A33 + A33^2
end

# issue #6450
@test dot(Any[1.0, 2.0], Any[3.5, 4.5]) === 12.5

@testset "dot" for elty in (Float32, Float64, ComplexF32, ComplexF64)
    x = convert(Vector{elty}, [1.0, 2.0, 3.0])
    y = convert(Vector{elty}, [3.5, 4.5, 5.5])
    @test_throws DimensionMismatch dot(x, 1:2, y, 1:3)
    @test_throws BoundsError dot(x, 1:4, y, 1:4)
    @test_throws BoundsError dot(x, 1:3, y, 2:4)
    @test dot(x, 1:2, y, 1:2) == convert(elty, 12.5)
    @test transpose(x) * y == convert(elty, 29.0)
    X = convert(Matrix{elty}, [1.0 2.0; 3.0 4.0])
    Y = convert(Matrix{elty}, [1.5 2.5; 3.5 4.5])
    @test dot(X, Y) == convert(elty, 35.0)
    Z = Matrix{elty}[reshape(1:4, 2, 2), fill(1, 2, 2)]
    @test dot(Z, Z) == convert(elty, 34.0)
end

dot1(x, y) = invoke(dot, Tuple{Any,Any}, x, y)
dot2(x, y) = invoke(dot, Tuple{AbstractArray,AbstractArray}, x, y)
@testset "generic dot" begin
    AA = [1+2im 3+4im; 5+6im 7+8im]
    BB = [2+7im 4+1im; 3+8im 6+5im]
    for A in (copy(AA), view(AA, 1:2, 1:2)), B in (copy(BB), view(BB, 1:2, 1:2))
        @test dot(A, B) == dot(vec(A), vec(B)) == dot1(A, B) == dot2(A, B) == dot(float.(A), float.(B))
        @test dot(Int[], Int[]) == 0 == dot1(Int[], Int[]) == dot2(Int[], Int[])
        @test_throws MethodError dot(Any[], Any[])
        @test_throws MethodError dot1(Any[], Any[])
        @test_throws MethodError dot2(Any[], Any[])
        for n1 = 0:2, n2 = 0:2, d in (dot, dot1, dot2)
            if n1 != n2
                @test_throws DimensionMismatch d(1:n1, 1:n2)
            else
                @test d(1:n1, 1:n2) ≈ norm(1:n1)^2
            end
        end
    end
end

@testset "Issue 11978" begin
    A = Matrix{Matrix{Float64}}(undef, 2, 2)
    A[1, 1] = Matrix(1.0I, 3, 3)
    A[2, 2] = Matrix(1.0I, 2, 2)
    A[1, 2] = Matrix(1.0I, 3, 2)
    A[2, 1] = Matrix(1.0I, 2, 3)
    b = Vector{Vector{Float64}}(undef, 2)
    b[1] = fill(1.0, 3)
    b[2] = fill(1.0, 2)
    @test A * b == Vector{Float64}[[2, 2, 1], [2, 2]]
end

@test_throws ArgumentError LinearAlgebra.copytri!(Matrix{Float64}(undef, 10, 10), 'Z')

@testset "Issue 30055" begin
    B = [1+im 2+im 3+im; 4+im 5+im 6+im; 7+im 9+im im]
    A = UpperTriangular(B)
    @test copy(transpose(A)) == transpose(A)
    @test copy(A') == A'
    A = LowerTriangular(B)
    @test copy(transpose(A)) == transpose(A)
    @test copy(A') == A'
    B = Matrix{Matrix{Complex{Int}}}(undef, 2, 2)
    B[1, 1] = [1+im 2+im; 3+im 4+im]
    B[2, 1] = [1+2im 1+3im; 1+3im 1+4im]
    B[1, 2] = [7+im 8+2im; 9+3im 4im]
    B[2, 2] = [9+im 8+im; 7+im 6+im]
    A = UpperTriangular(B)
    @test copy(transpose(A)) == transpose(A)
    @test copy(A') == A'
    A = LowerTriangular(B)
    @test copy(transpose(A)) == transpose(A)
    @test copy(A') == A'
end

@testset "gemv! and gemm_wrapper for $elty" for elty in [Float32, Float64, ComplexF64, ComplexF32]
    A10x10, x10, x11 = Array{elty}.(undef, ((10, 10), 10, 11))
    @test_throws DimensionMismatch LinearAlgebra.gemv!(x10, 'N', A10x10, x11)
    @test_throws DimensionMismatch LinearAlgebra.gemv!(x11, 'N', A10x10, x10)
    @test LinearAlgebra.gemv!(elty[], 'N', Matrix{elty}(undef, 0, 0), elty[]) == elty[]
    @test LinearAlgebra.gemv!(x10, 'N', Matrix{elty}(undef, 10, 0), elty[]) == zeros(elty, 10)

    I0x0 = Matrix{elty}(I, 0, 0)
    I10x10 = Matrix{elty}(I, 10, 10)
    I10x11 = Matrix{elty}(I, 10, 11)
    @test LinearAlgebra.gemm_wrapper('N', 'N', I10x10, I10x10) == I10x10
    @test_throws DimensionMismatch LinearAlgebra.gemm_wrapper!(I10x10, 'N', 'N', I10x11, I10x10)
    @test_throws DimensionMismatch LinearAlgebra.gemm_wrapper!(I10x10, 'N', 'N', I0x0, I0x0)

    A = rand(elty, 3, 3)
    @test LinearAlgebra.matmul3x3('T', 'N', A, Matrix{elty}(I, 3, 3)) == transpose(A)
end

@testset "#13593, #13488" begin
    aa = rand(3, 3)
    bb = rand(3, 3)
    for a in (copy(aa), view(aa, 1:3, 1:3)), b in (copy(bb), view(bb, 1:3, 1:3))
        @test_throws ArgumentError mul!(a, a, b)
        @test_throws ArgumentError mul!(a, b, a)
        @test_throws ArgumentError mul!(a, a, a)
    end
end

@testset "#35163" begin
    # typemax(Int32) * Int32(1) + Int32(1) * Int32(1) should wrap around
    # not promote to Int64, convert to Int32 and throw inexacterror
    val = mul!(Int32[1], fill(typemax(Int32), 1, 1), Int32[1], Int32(1), Int32(1))
    @test val[1] == typemin(Int32)
end

# Number types that lack conversion to the destination type
struct RootInt
    i::Int
end
import Base: *, adjoint, transpose
import LinearAlgebra: Adjoint, Transpose
(*)(x::RootInt, y::RootInt) = x.i * y.i
(*)(x::RootInt, y::Integer) = x.i * y
adjoint(x::RootInt) = x
transpose(x::RootInt) = x
Base.zero(::RootInt) = RootInt(0)

@test Base.promote_op(*, RootInt, RootInt) === Int

@testset "#14293" begin
    a = [RootInt(3)]
    C = [0;;]
    mul!(C, a, transpose(a))
    @test C[1] == 9
    C = [1;;]
    mul!(C, a, transpose(a), 2, 3)
    @test C[1] == 21
    a = [RootInt(2), RootInt(10)]
    @test a * adjoint(a) == [4 20; 20 100]
    A = [RootInt(3) RootInt(5)]
    @test A * a == [56]
end

function test_mul(C, A, B, S)
    mul!(C, A, B)
    @test Array(A) * Array(B) ≈ C
    @test A * B ≈ C

    # This is similar to how `isapprox` choose `rtol` (when `atol=0`)
    # but consider all number types involved:
    rtol = max(rtoldefault.(real.(eltype.((C, A, B))))...)

    rand!(C, S)
    T = promote_type(eltype.((A, B))...)
    α = T <: AbstractFloat ? rand(T) : rand(T(-10):T(10))
    β = T <: AbstractFloat ? rand(T) : rand(T(-10):T(10))
    βArrayC = β * Array(C)
    βC = β * C
    mul!(C, A, B, α, β)
    @test α * Array(A) * Array(B) .+ βArrayC ≈ C rtol = rtol
    @test α * A * B .+ βC ≈ C rtol = rtol
end

@testset "mul! vs * for special types" begin
    eltypes = [Float32, Float64, Int64(-100):Int64(100)]
    for k in [3, 4, 10]
        T = rand(eltypes)
        bi1 = Bidiagonal(rand(T, k), rand(T, k - 1), rand([:U, :L]))
        bi2 = Bidiagonal(rand(T, k), rand(T, k - 1), rand([:U, :L]))
        tri1 = Tridiagonal(rand(T, k - 1), rand(T, k), rand(T, k - 1))
        tri2 = Tridiagonal(rand(T, k - 1), rand(T, k), rand(T, k - 1))
        stri1 = SymTridiagonal(rand(T, k), rand(T, k - 1))
        stri2 = SymTridiagonal(rand(T, k), rand(T, k - 1))
        C = rand(T, k, k)
        specialmatrices = (bi1, bi2, tri1, tri2, stri1, stri2)
        for A in specialmatrices
            B = specialmatrices[rand(1:length(specialmatrices))]
            test_mul(C, A, B, T)
        end
        for S in specialmatrices
            l = rand(1:6)
            B = randn(k, l)
            C = randn(k, l)
            test_mul(C, S, B, T)
            A = randn(l, k)
            C = randn(l, k)
            test_mul(C, A, S, T)
        end
    end
    for T in eltypes
        A = Bidiagonal(rand(T, 2), rand(T, 1), rand([:U, :L]))
        B = Bidiagonal(rand(T, 2), rand(T, 1), rand([:U, :L]))
        C = randn(2, 2)
        test_mul(C, A, B, T)
        B = randn(2, 9)
        C = randn(2, 9)
        test_mul(C, A, B, T)
    end
    let
        tri44 = Tridiagonal(randn(3), randn(4), randn(3))
        tri33 = Tridiagonal(randn(2), randn(3), randn(2))
        full43 = randn(4, 3)
        full24 = randn(2, 4)
        full33 = randn(3, 3)
        full44 = randn(4, 4)
        @test_throws DimensionMismatch mul!(full43, tri44, tri33)
        @test_throws DimensionMismatch mul!(full44, tri44, tri33)
        @test_throws DimensionMismatch mul!(full44, tri44, full43)
        @test_throws DimensionMismatch mul!(full43, tri33, full43)
        @test_throws DimensionMismatch mul!(full43, full43, tri44)
    end
end

# #18218
module TestPR18218
using Test
import Base.*, Base.+, Base.zero
struct TypeA
    x::Int
end
Base.convert(::Type{TypeA}, x::Int) = TypeA(x)
struct TypeB
    x::Int
end
struct TypeC
    x::Int
end
Base.convert(::Type{TypeC}, x::Int) = TypeC(x)
zero(c::TypeC) = TypeC(0)
zero(::Type{TypeC}) = TypeC(0)
(*)(x::Int, a::TypeA) = TypeB(x * a.x)
(*)(a::TypeA, x::Int) = TypeB(a.x * x)
(+)(a::Union{TypeB,TypeC}, b::Union{TypeB,TypeC}) = TypeC(a.x + b.x)
A = TypeA[1 2; 3 4]
b = [1, 2]
d = A * b
@test typeof(d) == Vector{TypeC}
@test d == TypeC[5, 11]
end

@testset "VecOrMat of Vectors" begin
    X = rand(ComplexF64, 3, 3)
    Xv1 = [X[:, j] for i in 1:1, j in 1:3]
    Xv2 = [transpose(X[i, :]) for i in 1:3]
    Xv3 = [transpose(X[i, :]) for i in 1:3, j in 1:1]

    XX = X * X
    XtX = transpose(X) * X
    XcX = X' * X
    XXt = X * transpose(X)
    XtXt = transpose(XX)
    XcXt = X' * transpose(X)
    XXc = X * X'
    XtXc = transpose(X) * X'
    XcXc = X' * X'

    @test (Xv1*Xv2)[1] ≈ XX
    @test (Xv1*Xv3)[1] ≈ XX
    @test transpose(Xv1) * Xv1 ≈ XtX
    @test transpose(Xv2) * Xv2 ≈ XtX
    @test (transpose(Xv3)*Xv3)[1] ≈ XtX
    @test Xv1' * Xv1 ≈ XcX
    @test Xv2' * Xv2 ≈ XcX
    @test (Xv3'*Xv3)[1] ≈ XcX
    @test (Xv1*transpose(Xv1))[1] ≈ XXt
    @test Xv2 * transpose(Xv2) ≈ XXt
    @test Xv3 * transpose(Xv3) ≈ XXt
    @test transpose(Xv1) * transpose(Xv2) ≈ XtXt
    @test transpose(Xv1) * transpose(Xv3) ≈ XtXt
    @test Xv1' * transpose(Xv2) ≈ XcXt
    @test Xv1' * transpose(Xv3) ≈ XcXt
    @test (Xv1*Xv1')[1] ≈ XXc
    @test Xv2 * Xv2' ≈ XXc
    @test Xv3 * Xv3' ≈ XXc
    @test transpose(Xv1) * Xv2' ≈ XtXc
    @test transpose(Xv1) * Xv3' ≈ XtXc
    @test Xv1' * Xv2' ≈ XcXc
    @test Xv1' * Xv3' ≈ XcXc
end

@testset "copyto! for matrices of matrices" begin
    A = [randn(ComplexF64, 2,3) for _ in 1:2, _ in 1:3]
    for (tfun, tM) in ((identity, 'N'), (transpose, 'T'), (adjoint, 'C'))
        At = copy(tfun(A))
        B = zero.(At)
        copyto!(B, axes(B, 1), axes(B, 2), tM, A, axes(A, tM == 'N' ? 1 : 2), axes(A, tM == 'N' ? 2 : 1))
        @test B == At
    end
end

@testset "method ambiguity" begin
    # Ambiguity test is run inside a clean process.
    # https://github.com/JuliaLang/julia/issues/28804
    script = joinpath(@__DIR__, "ambiguous_exec.jl")
    cmd = `$(Base.julia_cmd()) --startup-file=no $script`
    @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
end

struct A32092
    x::Float64
end
Base.:+(x::Float64, a::A32092) = x + a.x
Base.:*(x::Float64, a::A32092) = x * a.x
@testset "Issue #32092" begin
    @test ones(2, 2) * [A32092(1.0), A32092(2.0)] == fill(3.0, (2,))
end

@testset "strong zero" begin
    @testset for α in Any[false, 0.0, 0], n in 1:4
        C = ones(n, n)
        A = fill!(zeros(n, n), NaN)
        B = ones(n, n)
        @test mul!(copy(C), A, B, α, 1.0) == C
    end
end

@testset "CartesianIndex handling in _modify!" begin
    C = rand(10, 10)
    A = rand(10, 10)
    @test mul!(view(C, 1:10, 1:10), A, 0.5) == A * 0.5
end

@testset "Issue #33214: tiled generic mul!" begin
    n = 100
    A = rand(n, n)
    B = rand(n, n)
    C = zeros(n, n)
    mul!(C, A, B, -1 + 0im, 0)
    D = -A * B
    @test D ≈ C

    # Just in case dispatching on the surface API `mul!` is changed in the future,
    # let's test the function where the tiled multiplication is defined.
    fill!(C, 0)
    LinearAlgebra.generic_matmatmul!(C, 'N', 'N', A, B, LinearAlgebra.MulAddMul(-1, 0))
    @test D ≈ C
end

@testset "size zero types in matrix mult (see issue 39362)" begin
    A = [missing missing; missing missing]
    v = [missing, missing]
    @test (A * v == v) === missing
    M = fill(1.0, 2, 2)
    a = fill(missing, 2, 1)
    @test (a' * M * a == fill(missing, 1, 1)) === missing
end


@testset "multiplication of empty matrices without calling zero" begin
    r, c = rand(0:9, 2)
    A = collect(Number, rand(r, c))
    B = rand(c, 0)
    C = A * B
    @test size(C) == (r, 0)
    @test_throws MethodError zero(eltype(C))
end

@testset "Issue #33873: genmatmul! with empty operands" begin
    @test Matrix{Any}(undef, 0, 2) * Matrix{Any}(undef, 2, 3) == Matrix{Any}(undef, 0, 3)
    @test_throws MethodError Matrix{Any}(undef, 2, 0) * Matrix{Any}(undef, 0, 3)
    @test Matrix{Int}(undef, 2, 0) * Matrix{Int}(undef, 0, 3) == zeros(Int, 2, 3)
end

@testset "3-arg *, order by type" begin
    x = [1, 2im]
    y = [im, 20, 30 + 40im]
    z = [-1, 200 + im, -3]
    A = [1 2 3im; 4 5 6+im]
    B = [-10 -20; -30 -40]
    a = 3 + im * round(Int, 10^6 * (pi - 3))
    b = 123

    @test x' * A * y == (x' * A) * y == x' * (A * y)
    @test y' * A' * x == (y' * A') * x == y' * (A' * x)
    @test y' * transpose(A) * x == (y' * transpose(A)) * x == y' * (transpose(A) * x)

    @test B * A * y == (B * A) * y == B * (A * y)

    @test a * A * y == (a * A) * y == a * (A * y)
    @test A * y * a == (A * y) * a == A * (y * a)

    @test a * B * A == (a * B) * A == a * (B * A)
    @test B * A * a == (B * A) * a == B * (A * a)

    @test a * y' * z == (a * y') * z == a * (y' * z)
    @test y' * z * a == (y' * z) * a == y' * (z * a)

    @test a * y * z' == (a * y) * z' == a * (y * z')
    @test y * z' * a == (y * z') * a == y * (z' * a)

    @test a * x' * A == (a * x') * A == a * (x' * A)
    @test x' * A * a == (x' * A) * a == x' * (A * a)
    @test a * x' * A isa Adjoint{<:Any,<:Vector}

    @test a * transpose(x) * A == (a * transpose(x)) * A == a * (transpose(x) * A)
    @test transpose(x) * A * a == (transpose(x) * A) * a == transpose(x) * (A * a)
    @test a * transpose(x) * A isa Transpose{<:Any,<:Vector}

    @test x' * B * A == (x' * B) * A == x' * (B * A)
    @test x' * B * A isa Adjoint{<:Any,<:Vector}

    @test y * x' * A == (y * x') * A == y * (x' * A)
    y31 = reshape(y, 3, 1)
    @test y31 * x' * A == (y31 * x') * A == y31 * (x' * A)

    vm = [rand(1:9, 2, 2) for _ in 1:3]
    Mm = [rand(1:9, 2, 2) for _ in 1:3, _ in 1:3]

    @test vm' * Mm * vm == (vm' * Mm) * vm == vm' * (Mm * vm)
    @test Mm * Mm' * vm == (Mm * Mm') * vm == Mm * (Mm' * vm)
    @test vm' * Mm * Mm == (vm' * Mm) * Mm == vm' * (Mm * Mm)
    @test Mm * Mm' * Mm == (Mm * Mm') * Mm == Mm * (Mm' * Mm)
end

@testset "3-arg *, order by size" begin
    M44 = randn(4, 4)
    M24 = randn(2, 4)
    M42 = randn(4, 2)
    @test M44 * M44 * M44 ≈ (M44 * M44) * M44 ≈ M44 * (M44 * M44)
    @test M42 * M24 * M44 ≈ (M42 * M24) * M44 ≈ M42 * (M24 * M44)
    @test M44 * M42 * M24 ≈ (M44 * M42) * M24 ≈ M44 * (M42 * M24)
end

@testset "4-arg *, by type" begin
    y = [im, 20, 30 + 40im]
    z = [-1, 200 + im, -3]
    a = 3 + im * round(Int, 10^6 * (pi - 3))
    b = 123
    M = rand(vcat(1:9, im .* [1, 2, 3]), 3, 3)
    N = rand(vcat(1:9, im .* [1, 2, 3]), 3, 3)

    @test a * b * M * y == (a * b) * (M * y)
    @test a * b * M * N == (a * b) * (M * N)
    @test a * M * N * y == (a * M) * (N * y)
    @test a * y' * M * z == (a * y') * (M * z)
    @test a * y' * M * N == (a * y') * (M * N)

    @test M * y * a * b == (M * y) * (a * b)
    @test M * N * a * b == (M * N) * (a * b)
    @test M * N * y * a == (a * M) * (N * y)
    @test y' * M * z * a == (a * y') * (M * z)
    @test y' * M * N * a == (a * y') * (M * N)

    @test M * N * conj(M) * y == (M * N) * (conj(M) * y)
    @test y' * M * N * conj(M) == (y' * M) * (N * conj(M))
    @test y' * M * N * z == (y' * M) * (N * z)
end

@testset "4-arg *, by size" begin
    for shift in 1:5
        s1, s2, s3, s4, s5 = circshift(3:7, shift)
        a = randn(s1, s2)
        b = randn(s2, s3)
        c = randn(s3, s4)
        d = randn(s4, s5)

        # _quad_matmul
        @test *(a, b, c, d) ≈ (a * b) * (c * d)

        # _tri_matmul(A,B,B,δ)
        @test *(11.1, b, c, d) ≈ (11.1 * b) * (c * d)
        @test *(a, b, c, 99.9) ≈ (a * b) * (c * 99.9)
    end
end

#46865
@testset "mul!() with non-const alpha, beta" begin
    f!(C,A,B,alphas,betas) = mul!(C, A, B, alphas[1], betas[1])
    alphas = [1.0]
    betas = [0.5]
    for d in [2,3,4]  # test native small-matrix cases as well as BLAS
        A = rand(d,d)
        B = copy(A)
        C = copy(A)
        f!(C, A, B, alphas, betas)
        @test (@allocated f!(C, A, B, alphas, betas)) == 0
    end
end

@testset "vector-matrix multiplication" begin
    a = [1,2]
    A = reshape([1,2], 2, 1)
    B = [1 2]
    @test a * B ≈ A * B
    B = reshape([1,2], 2, 1)
    @test a * B' ≈ A * B'
    @test a * transpose(B) ≈ A * transpose(B)
end

@testset "issue #56085" begin
    struct Thing
        data::Float64
    end

    Base.zero(::Type{Thing}) = Thing(0.)
    Base.zero(::Thing)       = Thing(0.)
    Base.one(::Type{Thing})  = Thing(1.)
    Base.one(::Thing)        = Thing(1.)
    Base.:+(t1::Thing, t::Thing...) = +(getfield.((t1, t...), :data)...)
    Base.:*(t1::Thing, t::Thing...) = *(getfield.((t1, t...), :data)...)

    M = Float64[1 2; 3 4]
    A = Thing.(M)

    @test A * A ≈ M * M
end

@testset "issue #1147: error messages in matmul" begin
    for T in (Int, Float64, ComplexF64)
        for f in (identity, Symmetric)
            @test_throws "incompatible dimensions for matrix multiplication" f(zeros(T,0,0)) * zeros(T,1,5)
            @test_throws "incompatible dimensions for matrix multiplication" f(zeros(T,0,0)) * zeros(T,1)
            @test_throws "incompatible dimensions for matrix multiplication" zeros(T,0) * f(zeros(T,2,2))
            @test_throws "incompatible dimensions for matrix multiplication" mul!(zeros(T,0,0), zeros(T,5), zeros(T,5))
            @test_throws "incompatible dimensions for matrix multiplication" mul!(zeros(T,0,0), f(zeros(T,1,1)), zeros(T,0,0))
            @test_throws "incompatible destination size" mul!(zeros(T,0,2), f(zeros(T,1,1)), zeros(T,1,2))
            @test_throws "incompatible destination size" mul!(zeros(T,1,0), f(zeros(T,1,1)), zeros(T,1,2))
            @test_throws "incompatible destination size" mul!(zeros(T,0,0), f(zeros(T,1,1)), zeros(T,1))
            @test_throws "incompatible destination size" mul!(zeros(T,0), f(zeros(T,1,1)), zeros(T,1))
        end

        @test_throws "expected size: (2, 2)" LinearAlgebra.matmul2x2!(zeros(T,2,2), 'N', 'N', zeros(T,2,3), zeros(T,3,2))
        @test_throws "expected size: (2, 2)" LinearAlgebra.matmul2x2!(zeros(T,2,3), 'N', 'N', zeros(T,2,2), zeros(T,2,3))
    end
end

@testset "zero-length generic matvec" begin
    m = SizedArrays.SizedArray{(2,2)}(ones(2,2))
    A = fill(m, 2, 0)
    v = fill(m, size(A,2))
    w = similar(v, size(A,1))
    mul!(w, A, v)
    @test all(iszero, w)
    A = fill(m, 0, 2)
    mul!(w, A', v)
    @test all(iszero, w)
end

@testset "zero-size matmul" begin
    A = zeros(0,2)
    S = Symmetric(zeros(0,0))
    @test S * A == A
    @test A' * S == A'
    S = Symmetric(zeros(2,2))
    @test S * A' == A'
    @test A * S == A
end

@testset "BlasFlag.NONE => generic_matmatmul!" begin
    A = ones(2,2)
    S = Symmetric(ones(2,2))
    @test mul!(similar(A), S, A, big(1), big(0)) ≈ S * A
    C1 = mul!(one(A), S, A, big(2), big(1))
    C2 = mul!(one(A), S, A, 2, 1)
    @test C1 ≈ C2
end

end # module TestMatmul
