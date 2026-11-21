# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestFactorization

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)
using Main.LinearAlgebraTestHelpers: COMPREHENSIVE

@testset "equality for factorizations - $f" for f in Any[
    bunchkaufman,
    cholesky,
    x -> cholesky(x, RowMaximum()),
    eigen,
    hessenberg,
    lq,
    lu,
    qr,
    x -> qr(x, ColumnNorm()),
    svd,
    schur,
]
    A = randn(3, 3)
    A = A * A' # ensure A is pos. def. and symmetric
    F, G = f(A), f(A)

    @test F == G
    @test isequal(F, G)
    @test hash(F) == hash(G)

    f === hessenberg && continue

    # change all arrays in F to have eltype Float32
    F = typeof(F).name.wrapper(Base.mapany(1:nfields(F)) do i
        x = getfield(F, i)
        return x isa AbstractArray{Float64} ? Float32.(x) : x
    end...)
    # round all arrays in G to the nearest Float64 representable as Float32
    G = typeof(G).name.wrapper(Base.mapany(1:nfields(G)) do i
        x = getfield(G, i)
        return x isa AbstractArray{Float64} ? Float64.(Float32.(x)) : x
    end...)

    @test F == G broken=!(f === eigen || f === qr || f == bunchkaufman || f == cholesky || F isa CholeskyPivoted)
    @test isequal(F, G) broken=!(f === eigen || f === qr || f == bunchkaufman || f == cholesky || F isa CholeskyPivoted)
    @test hash(F) == hash(G)
end

@testset "size for factorizations - $f" for f in Any[
    bunchkaufman,
    cholesky,
    x -> cholesky(x, RowMaximum()),
    hessenberg,
    lq,
    lu,
    qr,
    x -> qr(x, ColumnNorm()),
    svd,
]
    A = randn(3, 3)
    A = A * A' # ensure A is pos. def. and symmetric
    F = f(A)
    @test size(F) == size(A)
    @test size(F') == size(A')
end

@testset "size for transpose factorizations - $f" for f in Any[
    bunchkaufman,
    cholesky,
    x -> cholesky(x, RowMaximum()),
    hessenberg,
    lq,
    lu,
    svd,
]
    A = randn(3, 3)
    A = A * A' # ensure A is pos. def. and symmetric
    F = f(A)
    @test size(F) == size(A)
    @test size(transpose(F)) == size(transpose(A))
end

@testset "equality of QRCompactWY" begin
    A = rand(100, 100)
    F, G = qr(A), qr(A)

    @test F == G
    @test isequal(F, G)
    @test hash(F) == hash(G)

    G.T[28, 100] = 42

    @test F != G
    @test !isequal(F, G)
    @test hash(F) != hash(G)
end

# A dense array type that is not an `Array`, like FixedSizeArrays.jl
struct WrappedDenseArray{T,N} <: DenseArray{T,N}
    data::Array{T,N}
end
Base.size(A::WrappedDenseArray) = size(A.data)
Base.IndexStyle(::Type{<:WrappedDenseArray}) = IndexLinear()
Base.getindex(A::WrappedDenseArray, i::Int) = A.data[i]
Base.setindex!(A::WrappedDenseArray, v, i::Int) = setindex!(A.data, v, i)
Base.similar(::WrappedDenseArray, ::Type{T}, dims::Dims) where {T} = WrappedDenseArray(Array{T}(undef, dims))
# strided array interface, required by LAPACK
Base.unsafe_convert(::Type{Ptr{T}}, A::WrappedDenseArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, A.data)
Base.elsize(::Type{WrappedDenseArray{T,N}}) where {T,N} = Base.elsize(Array{T,N})

@testset "\\ with a factorization keeps the array type of the rhs (#1742)" begin
    for T in (Float64, ComplexF64)
        A = T[4 1 0.5; 1 3 0.2; 0.5 0.2 2]
        Awide, Atall = A[1:2, :], [A; 1 2 3]
        Fs = Any[lu(A), cholesky(Hermitian(A)), bunchkaufman(Hermitian(A)),
                 lq(A), lq(Awide), qr(A), qr(Atall)]
        for Ã in (A, Awide, Atall)
            push!(Fs, qr(Ã, ColumnNorm()), svd(Ã))
        end
        for F in Fs, cols in ((), (2,))
            B = rand(T, size(F, 1), cols...)
            @test (F \ WrappedDenseArray(B))::WrappedDenseArray{T} ≈ F \ B
        end
        @test (WrappedDenseArray(A) \ WrappedDenseArray(T[1, 2, 3]))::WrappedDenseArray{T} ≈ A \ T[1, 2, 3]
        # right hand sides that are not strided still give an `Array`
        @test (lu(A) \ (1:3))::Vector{T} ≈ lu(A) \ [1, 2, 3]
    end
end

@testset "deprecated conversion of factorizations to arrays" begin
    A = [4.0 1.0; 1.0 3.0]
    for F in (lu(A), cholesky(A), qr(A), svd(A), eigen(A), schur(A), hessenberg(A), lq(A))
        @test_deprecated convert(Matrix, F)
        @test_deprecated convert(AbstractMatrix, F)
    end
end

@testset "ldiv!/rdiv! into an array that shares memory with the factorization" begin
    n = 4
    @testset for T in (Float64, BigFloat, (@static COMPREHENSIVE ? (ComplexF64,) : ())...)
        spd() = (X = rand(T, n, n); X * X' + n * I)
        # factorized in place, so that the factorization is stored in M
        factorizations = Any[lu!, M -> cholesky!(Hermitian(M)), M -> cholesky!(Hermitian(M), RowMaximum()),
            qr!, M -> qr!(M, ColumnNorm()), lq!, M -> bunchkaufman!(Hermitian(M))]
        T <: LinearAlgebra.BlasFloat && push!(factorizations, hessenberg!)
        for fact in factorizations
            M = spd(); F = fact(M); M0 = copy(M)
            @test_throws ArgumentError ldiv!(F, M)
            @test_throws ArgumentError ldiv!(F, view(M, :, 1))
            @test_throws ArgumentError ldiv!(F', M)
            @test_throws ArgumentError ldiv!(M, F, rand(T, n, n))
            if F isa Union{LU,Cholesky,CholeskyPivoted,Hessenberg}
                @test_throws ArgumentError rdiv!(M, F)
            end
            @test M == M0 # the check comes before anything is written
        end
        # factorizations of banded matrices, stored in their bands
        for F in (ldlt!(SymTridiagonal(rand(T, n) .+ 10, rand(T, n-1))),
                  lu!(Tridiagonal(rand(T, n-1), rand(T, n) .+ 10, rand(T, n-1))))
            d = F isa LDLt ? F.data.dv : F.factors.d
            d0 = copy(d)
            @test_throws ArgumentError ldiv!(F, d)
            @test_throws ArgumentError ldiv!(d, F, rand(T, n))
            @test d == d0
        end
        # these solve on a copy, but overwriting the factorization would make it unusable
        if T <: LinearAlgebra.BlasFloat # no generic svd
            F = svd!(spd()); U0 = copy(F.U)
            @test_throws ArgumentError ldiv!(F, F.U)
            @test F.U == U0
        end
        M = spd(); F = bunchkaufman!(Hermitian(M)); M0 = copy(M)
        @test_throws ArgumentError rdiv!(M, F)
        @test M == M0
    end
end

end
