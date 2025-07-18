# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestDiagonal

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: BlasFloat, BlasComplex

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.OffsetArrays
using Main.LinearAlgebraTestHelpers.InfiniteArrays
using Main.LinearAlgebraTestHelpers.FillArrays
using Main.LinearAlgebraTestHelpers.SizedArrays
using Main.LinearAlgebraTestHelpers.ImmutableArrays

const n=12 # Size of matrix problem to test
Random.seed!(1)

# wrapper to avoid dispatching to diagonal methods
struct NotDiagonal{T,A<:AbstractMatrix{T}} <: AbstractMatrix{T}
    a :: A
end
Base.size(N::NotDiagonal) = size(N.a)
Base.getindex(N::NotDiagonal, i::Int, j::Int) = N.a[i, j]
LinearAlgebra.isdiag(N::NotDiagonal) = false # this contradicts `getindex`
LinearAlgebra.ishermitian(N::NotDiagonal) = ishermitian(N.a)
LinearAlgebra.istriu(N::NotDiagonal) = istriu(N.a)
LinearAlgebra.istril(N::NotDiagonal) = istril(N.a)

@testset for relty in (Float32, Float64, BigFloat), elty in (relty, Complex{relty})
    dd=convert(Vector{elty}, randn(n))
    vv=convert(Vector{elty}, randn(n))
    UU=convert(Matrix{elty}, randn(n,n))
    if elty <: Complex
        dd+=im*convert(Vector{elty}, randn(n))
        vv+=im*convert(Vector{elty}, randn(n))
        UU+=im*convert(Matrix{elty}, randn(n,n))
    end
    D = Diagonal(dd)
    M = Matrix(D)
    # we can't directly compare with a Matrix, since the dense methods often dispatch
    # to Diagonal ones. We therefore compare with other structured matrix types
    # which have their own implementations.
    # We wrap the complex matrices in NotDiagonal to avoid falling back to Diagonal methods
    DM = elty <: Real ? Hermitian(M) : NotDiagonal(UpperTriangular(M))

    @testset "constructor" begin
        for x in (dd, GenericArray(dd))
            @test Diagonal(x)::Diagonal{elty,typeof(x)} == M
            @test Diagonal(x).diag === x
            @test Diagonal{elty}(x)::Diagonal{elty,typeof(x)} == M
            @test Diagonal{elty}(x).diag === x
            @test Diagonal{elty}(D) === D
        end
        @test eltype(Diagonal{elty}([1,2,3,4])) == elty
        @test isa(Diagonal{elty,Vector{elty}}(GenericArray([1,2,3,4])), Diagonal{elty,Vector{elty}})
        @test isa(Diagonal{elty}(rand(Int,n,n)), Diagonal{elty,Vector{elty}})
        DI = Diagonal([1,2,3,4])
        @test Diagonal(DI) === DI
        @test isa(Diagonal{elty}(DI), Diagonal{elty})

        # diagonal matrices may be converted to Diagonal
        local A = [1 0; 0 2]
        local DA = convert(Diagonal{Float32,Vector{Float32}}, A)
        @test DA isa Diagonal{Float32,Vector{Float32}}
        @test DA == A

        # issue #26178
        @test_throws MethodError convert(Diagonal, [1,2,3,4])
        @test_throws DimensionMismatch convert(Diagonal, [1 2 3 4])
        @test_throws InexactError convert(Diagonal, ones(2,2))

        # Test reversing
        # Test reversing along rows
        @test reverse(D, dims=1) == reverse(Matrix(D), dims=1)

        # Test reversing along columns
        @test reverse(D, dims=2) == reverse(Matrix(D), dims=2)

        # Test reversing the entire matrix
        @test reverse(D)::Diagonal == reverse(Matrix(D)) == reverse!(copy(D))
    end

    @testset "Basic properties" begin
        @test_throws BoundsError size(D,0)
        @test size(D,1) == size(D,2) == length(dd)
        @test size(D,3) == 1
        @test typeof(convert(Diagonal{ComplexF32},D)) <: Diagonal{ComplexF32}
        @test typeof(convert(AbstractMatrix{ComplexF32},D)) <: Diagonal{ComplexF32}

        @test convert(Array, real(D)) == real(M)
        @test convert(Array, abs.(D)) == abs.(M)
        @test convert(Array, imag(D)) == imag(M)

        @test parent(D) == dd
        @test D[1,1] == dd[1]
        @test D[1,2] == 0

        @test issymmetric(D)
        @test isdiag(D)
        @test isdiag(Diagonal([[1 0; 0 1], [1 0; 0 1]]))
        @test !isdiag(Diagonal([[1 0; 0 1], [1 0; 1 1]]))
        @test istriu(D)
        @test istriu(D, -1)
        @test !istriu(D, 1)
        @test istriu(Diagonal(zero(diag(D))), 1)
        @test istril(D)
        @test !istril(D, -1)
        @test istril(D, 1)
        @test istril(Diagonal(zero(diag(D))), -1)
        @test Base.isstored(D,1,1)
        @test !Base.isstored(D,1,2)
        @test_throws BoundsError Base.isstored(D, n + 1, 1)
        if elty <: Real
            @test ishermitian(D)
        end
    end

    @testset "diag" begin
        @test isempty(@inferred diag(D,  n+1))
        @test isempty(@inferred diag(D, -n-1))
        @test (@inferred diag(D))::typeof(dd) == dd
        @test (@inferred diag(D, 0))::typeof(dd) == dd
        @test (@inferred diag(D, 1))::typeof(dd) == zeros(elty, n-1)
        DG = Diagonal(GenericArray(dd))
        @test (@inferred diag(DG))::typeof(GenericArray(dd)) == GenericArray(dd)
        @test (@inferred diag(DG, 1))::typeof(GenericArray(dd)) == GenericArray(zeros(elty, n-1))
    end


    @testset "Simple unary functions" begin
        for op in (-,)
            @test op(D)==op(DM)
        end

        for func in (det, tr)
            @test func(D) ≈ func(DM) atol=n^2*eps(relty)*(1+(elty<:Complex))
        end

        if eltype(D) <: Real
            @test minimum(D) ≈ minimum(DM)
            @test maximum(D) ≈ maximum(DM)
        end

        if relty <: BlasFloat
            for func in (exp, cis, sinh, cosh, tanh, sech, csch, coth)
                @test func(D) ≈ func(DM) atol=n^3*eps(relty)
            end
            @test log(Diagonal(abs.(D.diag))) ≈ log(abs.(DM)) atol=n^3*eps(relty)
        end
        if elty <: BlasComplex
            for func in (logdet, sqrt, sin, cos, tan, sec, csc, cot,
                         asin, acos, atan, asec, acsc, acot,
                         asinh, acosh, atanh, asech, acsch, acoth)
                @test func(D) ≈ func(DM) atol=n^2*eps(relty)*2
            end
        end
    end

    @testset "Two-dimensional Euler formula for Diagonal" begin
        @test cis(Diagonal([π, π])) ≈ -I
    end

    @testset "Linear solve" begin
        for (v, U) in ((vv, UU), (view(vv, 1:n), view(UU, 1:n, 1:2)))
            @test D*v ≈ DM*v atol=n*eps(relty)*(1+(elty<:Complex))
            @test D*U ≈ DM*U atol=n^2*eps(relty)*(1+(elty<:Complex))

            @test transpose(U)*D ≈ transpose(U)*M
            @test U'*D ≈ U'*M

            if relty != BigFloat
                atol_two = 2n^2 * eps(relty) * (1 + (elty <: Complex))
                atol_three = 2n^3 * eps(relty) * (1 + (elty <: Complex))
                @test D\v ≈ DM\v atol=atol_two
                @test D\U ≈ DM\U atol=atol_three
                @test ldiv!(D, copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(transpose(D), copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(adjoint(conj(D)), copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(D, copy(U)) ≈ DM\U atol=atol_three
                @test ldiv!(transpose(D), copy(U)) ≈ DM\U atol=atol_three
                @test ldiv!(adjoint(conj(D)), copy(U)) ≈ DM\U atol=atol_three
                # this method tests AbstractMatrix/AbstractVec for second arg
                Usym_bad = Symmetric(ones(elty, n+1, n+1))
                @test_throws DimensionMismatch ldiv!(D, copy(Usym_bad))

                @test ldiv!(zero(v), D, copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(zero(v), transpose(D), copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(zero(v), adjoint(conj(D)), copy(v)) ≈ DM\v atol=atol_two
                @test ldiv!(zero(U), D, copy(U)) ≈ DM\U atol=atol_three
                @test ldiv!(zero(U), transpose(D), copy(U)) ≈ DM\U atol=atol_three
                @test ldiv!(zero(U), adjoint(conj(D)), copy(U)) ≈ DM\U atol=atol_three

                Uc = copy(U')
                target = rmul!(Uc, Diagonal(inv.(D.diag)))
                @test rdiv!(Uc, D) ≈ target atol=atol_three
                @test_throws DimensionMismatch rdiv!(Matrix{elty}(I, n-1, n-1), D)
                @test_throws SingularException rdiv!(Uc, Diagonal(fill!(similar(D.diag), 0)))
                @test rdiv!(Uc, transpose(D)) ≈ target atol=atol_three
                @test rdiv!(Uc, adjoint(conj(D))) ≈ target atol=atol_three
                @test ldiv!(D, Matrix{eltype(D)}(I, size(D))) ≈ D \ Matrix{eltype(D)}(I, size(D)) atol=atol_three
                @test_throws DimensionMismatch ldiv!(D, fill(elty(1), n + 1))
                @test_throws SingularException ldiv!(Diagonal(zeros(relty, n)), copy(v))
                b = rand(elty, n, n)
                @test ldiv!(D, copy(b)) ≈ M\b
                @test_throws SingularException ldiv!(Diagonal(zeros(elty, n)), copy(b))
                b = view(rand(elty, n), Vector(1:n))
                b2 = copy(b)
                c = ldiv!(D, b)
                d = M\b2
                @test c ≈ d
                @test_throws SingularException ldiv!(Diagonal(zeros(elty, n)), b)
                b = rand(elty, n+1, n+1)
                @test_throws DimensionMismatch ldiv!(D, copy(b))
                b = view(rand(elty, n+1), Vector(1:n+1))
                @test_throws DimensionMismatch ldiv!(D, b)
            end
        end
    end
    d = convert(Vector{elty}, randn(n))
    D2 = Diagonal(d)
    DM2= Matrix(Diagonal(d))
    @testset "Binary operations" begin
        for op in (+, -, *)
            @test Array(op(D, D2)) ≈ op(DM, DM2)
        end
        @testset "with plain numbers" begin
            a = rand()
            @test Array(a*D) ≈ a*DM
            @test Array(D*a) ≈ DM*a
            @test Array(D/a) ≈ DM/a
            if elty <: Real
                @test convert(Array, abs.(D)^a) ≈ abs.(DM)^a
            else
                @test convert(Array, D^a) ≈ DM^a rtol=max(eps(relty), 1e-15) # TODO: improve precision
            end
            @test Diagonal(1:100)^2 == Diagonal((1:100).^2)
            p = 3
            @test Diagonal(1:100)^p == Diagonal((1:100).^p)
            @test Diagonal(1:100)^(-1) == Diagonal(inv.(1:100))
            @test Diagonal(1:100)^2.0 == Diagonal((1:100).^2.0)
            @test Diagonal(1:100)^(2.0+0im) == Diagonal((1:100).^(2.0+0im))
        end

        if relty <: BlasFloat
            for b in (rand(elty,n,n), rand(elty,n))
                @test lmul!(copy(D), copy(b)) ≈ M*b
                @test lmul!(transpose(copy(D)), copy(b)) ≈ transpose(M)*b
                @test lmul!(adjoint(copy(D)), copy(b)) ≈ M'*b
            end
        end

        #a few missing mults
        bd = Bidiagonal(D2)
        @test D*transpose(D2) ≈ M*transpose(DM2)
        @test D2*transpose(D) ≈ DM2*transpose(M)
        @test D2*D' ≈ DM2*M'

        #division of two Diagonals
        @test D/D2 ≈ Diagonal(D.diag./D2.diag)
        @test D\D2 ≈ Diagonal(D2.diag./D.diag)

        # QR \ Diagonal
        A = rand(elty, n, n)
        qrA = qr(A)
        @test qrA \ D ≈ A \ D

        # HermOrSym
        A     = rand(elty, n, n)
        Asym  = Symmetric(A + transpose(A), :U)
        Aherm = Hermitian(A + adjoint(A), :U)
        Msym = Array(Asym)
        Mherm = Array(Aherm)
        for op in (+, -)
            @test op(Asym, D) isa Symmetric
            @test convert(Array, op(Asym, D)) ≈ Array(Symmetric(op(Msym, M)))
            @test op(D, Asym) isa Symmetric
            @test convert(Array, op(D, Asym)) ≈ Array(Symmetric(op(M, Msym)))
            if !(elty <: Real)
                Dr = real(D)
                Mr = Array(Dr)
                @test op(Aherm, Dr) isa Hermitian
                @test convert(Array, op(Aherm, Dr)) ≈ Array(Hermitian(op(Mherm, Mr)))
                @test op(Dr, Aherm) isa Hermitian
                @test convert(Array, op(Dr, Aherm)) ≈ Array(Hermitian(op(Mr, Mherm)))
            end
        end
        Msym = Array(Asym)
        @test convert(Array, D*transpose(Asym)) ≈ M * convert(Array, transpose(Msym))
        @test convert(Array, D*adjoint(Asym)) ≈ M * convert(Array, adjoint(Asym))
        @test convert(Array, D*transpose(Aherm)) ≈ M * convert(Array, transpose(Aherm))
        @test convert(Array, D*adjoint(Aherm)) ≈ M * convert(Array, adjoint(Aherm))
        @test convert(Array, Asym*transpose(D)) ≈ Msym * convert(Array, transpose(D))
        @test convert(Array, transpose(D)*Asym) ≈ convert(Array, transpose(D)) * Msym
        @test convert(Array, adjoint(Aherm)*adjoint(D)) ≈ convert(Array, adjoint(Aherm)) * convert(Array, adjoint(D))
        @test convert(Array, adjoint(D)*adjoint(Aherm)) ≈ convert(Array, adjoint(D)) * convert(Array, adjoint(Aherm))

        # Performance specialisations for A*_mul_B!
        vvv = similar(vv)
        @test (r = M * vv   ; mul!(vvv, D, vv)  ≈ r ≈ vvv)
        @test (r = M' * vv  ; mul!(vvv, adjoint(D), vv) ≈ r ≈ vvv)
        @test (r = transpose(M) * vv ; mul!(vvv, transpose(D), vv) ≈ r ≈ vvv)

        UUU = similar(UU)
        for transformA in (identity, adjoint, transpose)
            for transformD in (identity, adjoint, transpose)
                @test mul!(UUU, transformA(UU), transformD(D)) ≈  transformA(UU) * Matrix(transformD(D))
                @test mul!(UUU, transformD(D), transformA(UU)) ≈  Matrix(transformD(D)) * transformA(UU)
            end
        end

        alpha = elty(randn())  # randn(elty) does not work with BigFloat
        beta = elty(randn())
        @testset begin
            vvv = similar(vv)
            vvv .= randn(size(vvv))  # randn!(vvv) does not work with BigFloat
            r = alpha * M * vv + beta * vvv
            @test mul!(vvv, D, vv, alpha, beta) === vvv
            @test r ≈ vvv
        end
        @testset begin
            vvv = similar(vv)
            vvv .= randn(size(vvv))  # randn!(vvv) does not work with BigFloat
            r = alpha * M' * vv + beta * vvv
            @test mul!(vvv, adjoint(D), vv, alpha, beta) === vvv
            @test r ≈ vvv
        end
        @testset begin
            vvv = similar(vv)
            vvv .= randn(size(vvv))  # randn!(vvv) does not work with BigFloat
            r = alpha * transpose(M) * vv + beta * vvv
            @test mul!(vvv, transpose(D), vv, alpha, beta) === vvv
            @test r ≈ vvv
        end

        @testset begin
            UUU = similar(UU)
            UUU .= randn(size(UUU))  # randn!(UUU) does not work with BigFloat
            r = alpha * M * UU + beta * UUU
            @test mul!(UUU, D, UU, alpha, beta) === UUU
            @test r ≈ UUU
        end
        @testset begin
            UUU = similar(UU)
            UUU .= randn(size(UUU))  # randn!(UUU) does not work with BigFloat
            r = alpha * M' * UU + beta * UUU
            @test mul!(UUU, adjoint(D), UU, alpha, beta) === UUU
            @test r ≈ UUU
        end
        @testset begin
            UUU = similar(UU)
            UUU .= randn(size(UUU))  # randn!(UUU) does not work with BigFloat
            r = alpha * transpose(M) * UU + beta * UUU
            @test mul!(UUU, transpose(D), UU, alpha, beta) === UUU
            @test r ≈ UUU
        end

        # make sure that mul!(A, {Adj|Trans}(B)) works with B as a Diagonal
        VV = Array(D)
        r  = VV * M
        @test rmul!(VV, D) ≈ r ≈ M*M
        if transpose(D) !== D
            r  = VV * transpose(M)
            @test rmul!(VV, transpose(D)) ≈ r
        end
        if adjoint(D) !== D
            r  = VV * M'
            @test rmul!(VV, adjoint(D)) ≈ r
        end

        # kron
        D3 = Diagonal(convert(Vector{elty}, rand(n÷2)))
        DM3= Matrix(D3)
        @test Matrix(kron(D, D3)) ≈ kron(DM, DM3)
        M4 = rand(elty, size(D3,1) + 1, size(D3,2) + 2) # choose a different size from D3
        @test kron(D3, M4) ≈ kron(DM3, M4)
        @test kron(M4, D3) ≈ kron(M4, DM3)
        X = [ones(1,1) for i in 1:2, j in 1:2]
        @test kron(I(2), X)[1,3] == zeros(1,1)
        X = [ones(2,2) for i in 1:2, j in 1:2]
        @test kron(I(2), X)[1,3] == zeros(2,2)
    end
    @testset "iszero, isone, triu, tril" begin
        Dzero = Diagonal(zeros(elty, 10))
        Done = Diagonal(ones(elty, 10))
        Dmix = Diagonal(zeros(elty, 10))
        Dmix[end,end] = one(elty)
        @test iszero(Dzero)
        @test !isone(Dzero)
        @test !iszero(Done)
        @test isone(Done)
        @test !iszero(Dmix)
        @test !isone(Dmix)
        @test istriu(D)
        @test istril(D)
        @test iszero(triu(D,1))
        @test triu(D,0)  == D
        @test triu(D,-1) == D
        @test tril(D,1)  == D
        @test iszero(tril(D,-1))
        @test tril(D,0)  == D
        @test_throws ArgumentError tril(D, -n - 2)
        @test_throws ArgumentError tril(D, n)
        @test_throws ArgumentError triu(D, -n)
        @test_throws ArgumentError triu(D, n + 2)
    end

    # factorize
    @test factorize(D) == D

    @testset "Eigensystem" begin
        eigD = eigen(D)
        @test Diagonal(eigD.values) == D
        @test eigD.vectors == Matrix(I, size(D))
        eigsortD = eigen(D, sortby=LinearAlgebra.eigsortby)
        @test eigsortD.values !== D.diag
        @test eigsortD.values == sort(D.diag, by=LinearAlgebra.eigsortby)
        @test Matrix(eigsortD) == D
    end

    @testset "ldiv" begin
        v = rand(n + 1)
        @test_throws DimensionMismatch D\v
        v = rand(n)
        @test D\v ≈ DM\v
        V = rand(n + 1, n)
        @test_throws DimensionMismatch D\V
        V = rand(n, n)
        @test D\V ≈ DM\V
    end

    @testset "conj and transpose" begin
        @test transpose(D) == D
        if elty <: Real
            @test transpose(D) === D
            @test adjoint(D) === D
        elseif elty <: BlasComplex
            @test Array(conj(D)) ≈ conj(DM)
            @test adjoint(D) == conj(D)
            local D2 = copy(D)
            local D2adj = adjoint(D2)
            D2adj[1,1] = rand(eltype(D2adj))
            @test D2[1,1] == adjoint(D2adj[1,1])
            @test D2adj' === D2
        end
        # Translates to Ac/t_mul_B, which is specialized after issue 21286
        @test(D' * vv == conj(D) * vv)
        @test(transpose(D) * vv == D * vv)
    end

    # logdet and logabsdet
    if relty <: Real
        lD = Diagonal(convert(Vector{relty}, rand(n)))
        lM = Matrix(lD)
        @test logdet(lD) ≈ logdet(lM)
        d1, s1 = @inferred logabsdet(lD)
        d2, s2 = logabsdet(lM)
        @test d1 ≈ d2
        @test s1 == s2
        @test logdet(Diagonal(relty[-1,-2])) ≈ log(2)
        @test_throws DomainError logdet(Diagonal(relty[-1,-2,-3]))
    end

    @testset "similar" begin
        @test isa(similar(D), Diagonal{elty})
        @test isa(similar(D, Int), Diagonal{Int})
        @test isa(similar(D, (3,2)), Matrix{elty})
        @test isa(similar(D, Int, (3,2)), Matrix{Int})
    end

    # Issue number 10036
    # make sure issymmetric/ishermitian work for
    # non-real diagonal matrices
    @testset "issymmetric/hermitian for complex Diagonal" begin
        @test issymmetric(D2)
        @test ishermitian(D2)
        if elty <: Complex
            dc = d .+ elty(1im)
            D3 = Diagonal(dc)
            @test issymmetric(D3)
            @test !ishermitian(D3)
        end
    end

    @testset "svd (#11120/#11247/#1149)" begin
        D[1] = 0
        U, s, V = svd(D)
        @test (U*Diagonal(s))*V' ≈ D
        @test svdvals(D) == s
        @test svd(D).V == V
    end
end

@testset "axes" begin
    v = OffsetArray(1:3)
    D = Diagonal(v)
    @test axes(D) isa NTuple{2,typeof(axes(v,1))}
end

@testset "rdiv! (#40887)" begin
    @test rdiv!(Matrix(Diagonal([2.0, 3.0])), Diagonal(2:3)) == Diagonal([1.0, 1.0])
    @test rdiv!(fill(3.0, 3, 3), 3.0I(3)) == ones(3,3)
end

@testset "kron (issue #40595)" begin
    # custom array type to test that kron on Diagonal matrices preserves types of the parents if possible
    struct KronTestArray{T, N, AT} <: AbstractArray{T, N}
        data::AT
    end
    KronTestArray(data::AbstractArray) = KronTestArray{eltype(data), ndims(data), typeof(data)}(data)
    Base.size(A::KronTestArray) = size(A.data)
    LinearAlgebra.kron(A::KronTestArray, B::KronTestArray) = KronTestArray(kron(A.data, B.data))
    Base.getindex(K::KronTestArray{<:Any,N}, i::Vararg{Int,N}) where {N} = K.data[i...]

    A = KronTestArray([1, 2, 3]);
    @test kron(A, A) isa KronTestArray
    Ad = Diagonal(A);
    @test kron(Ad, Ad).diag isa KronTestArray
    @test kron(Ad, Ad).diag == kron([1, 2, 3], [1, 2, 3])
end

# Define a vector type that does not support `deleteat!`, to ensure that `kron` handles this
struct SimpleVector{T} <: AbstractVector{T}
    vec::Vector{T}
end
SimpleVector(x::SimpleVector) = SimpleVector(Vector(x.vec))
SimpleVector{T}(::UndefInitializer, n::Integer) where {T} = SimpleVector(Vector{T}(undef, n))
Base.:(==)(x::SimpleVector, y::SimpleVector) = x == y
Base.axes(x::SimpleVector) = axes(x.vec)
Base.convert(::Type{Vector{T}}, x::SimpleVector) where {T} = convert(Vector{T}, x.vec)
Base.convert(::Type{Vector}, x::SimpleVector{T}) where {T} = convert(Vector{T}, x)
Base.convert(::Type{Array{T}}, x::SimpleVector) where {T} = convert(Vector{T}, x)
Base.convert(::Type{Array}, x::SimpleVector) = convert(Vector, x)
Base.copyto!(x::SimpleVector, y::SimpleVector) = (copyto!(x.vec, y.vec); x)
Base.eltype(::Type{SimpleVector{T}}) where {T} = T
Base.getindex(x::SimpleVector, ind...) = getindex(x.vec, ind...)
Base.kron(x::SimpleVector, y::SimpleVector) = SimpleVector(kron(x.vec, y.vec))
Base.promote_rule(::Type{<:AbstractVector{T}}, ::Type{SimpleVector{U}}) where {T,U} = Vector{promote_type(T, U)}
Base.promote_rule(::Type{SimpleVector{T}}, ::Type{SimpleVector{U}}) where {T,U} = SimpleVector{promote_type(T, U)}
Base.setindex!(x::SimpleVector, val, ind...) = (setindex!(x.vec, val, ind...), x)
Base.similar(x::SimpleVector, ::Type{T}) where {T} = SimpleVector(similar(x.vec, T))
Base.similar(x::SimpleVector, ::Type{T}, dims::Dims{1}) where {T} = SimpleVector(similar(x.vec, T, dims))
Base.size(x::SimpleVector) = size(x.vec)

@testset "kron (issue #46456)" for repr in Any[identity, SimpleVector]
    A = Diagonal(repr(randn(10)))
    M = Array(A)
    BL = Bidiagonal(repr(randn(10)), repr(randn(9)), :L)
    BU = Bidiagonal(repr(randn(10)), repr(randn(9)), :U)
    C = SymTridiagonal(repr(randn(10)), repr(randn(9)))
    Cl = SymTridiagonal(repr(randn(10)), repr(randn(10)))
    D = Tridiagonal(repr(randn(9)), repr(randn(10)), repr(randn(9)))
    @test kron(A, BL)::Bidiagonal == kron(M, Array(BL))
    @test kron(A, BU)::Bidiagonal == kron(M, Array(BU))
    @test kron(A, C)::SymTridiagonal == kron(M, Array(C))
    @test kron(A, Cl)::SymTridiagonal == kron(M, Array(Cl))
    @test kron(A, D)::Tridiagonal == kron(M, Array(D))
end

@testset "svdvals and eigvals (#11120/#11247)" begin
    D = Diagonal(Matrix{Float64}[randn(3,3), randn(2,2)])
    @test sort([svdvals(D)...;], rev = true) ≈ svdvals([D.diag[1] zeros(3,2); zeros(2,3) D.diag[2]])
    @test sort([eigvals(D)...;], by=LinearAlgebra.eigsortby) ≈ eigvals([D.diag[1] zeros(3,2); zeros(2,3) D.diag[2]])
end

@testset "eigvals should return a copy of the diagonal" begin
    D = Diagonal([1, 2, 3])
    lam = eigvals(D)
    D[3,3] = 4 # should not affect lam
    @test lam == [1, 2, 3]
end

@testset "eigmin (#27847)" begin
    for _ in 1:100
        d = randn(rand(1:10))
        D = Diagonal(d)
        @test eigmin(D) == minimum(d)
    end
end

@testset "isposdef" begin
    @test isposdef(Diagonal(1.0 .+ rand(n)))
    @test !isposdef(Diagonal(-1.0 * rand(n)))
    @test isposdef(Diagonal(complex(1.0, 0.0) .+ rand(n)))
    @test !isposdef(Diagonal(complex(1.0, 1.0) .+ rand(n)))
    @test isposdef(Diagonal([[1 0; 0 1], [1 0; 0 1]]))
    @test !isposdef(Diagonal([[1 0; 0 1], [1 0; 1 1]]))
end

@testset "getindex" begin
    d = randn(n)
    D = Diagonal(d)
    # getindex bounds checking
    @test_throws BoundsError D[0, 0]
    @test_throws BoundsError D[-1, -2]
    @test_throws BoundsError D[n, n + 1]
    @test_throws BoundsError D[n + 1, n]
    @test_throws BoundsError D[n + 1, n + 1]
    # getindex on and off the diagonal
    for i in 1:n, j in 1:n
        @test D[i, j] == (i == j ? d[i] : 0)
    end
end

@testset "setindex!" begin
    d = randn(n)
    D = Diagonal(d)
    # setindex! bounds checking
    @test_throws BoundsError D[0, 0] = 0
    @test_throws BoundsError D[-1 , -2] = 0
    @test_throws BoundsError D[n, n + 1] = 0
    @test_throws BoundsError D[n + 1, n] = 0
    @test_throws BoundsError D[n + 1, n + 1] = 0
    for i in 1:n, j in 1:n
        if i == j
            # setindex on! the diagonal
            @test ((D[i, j] = i) == i; D[i, j] == i)
        else
            # setindex! off the diagonal
            @test ((D[i, j] = 0) == 0; iszero(D[i, j]))
            @test_throws ArgumentError D[i, j] = 1
        end
    end
    # setindex should return the destination
    @test setindex!(D, 1, 1, 1) === D
end

@testset "Test reverse" begin
    D = Diagonal(randn(5))
    M = Matrix(D)
    @test reverse(D, dims=1) == reverse(M, dims=1)
    @test reverse(D, dims=2) == reverse(M, dims=2)
    @test reverse(D)::Diagonal == reverse(M)
end

@testset "inverse" begin
    for d in Any[randn(n), Int[], [1, 2, 3], [1im, 2im, 3im], [1//1, 2//1, 3//1], [1+1im//1, 2//1, 3im//1]]
        D = Diagonal(d)
        @test inv(D) ≈ inv(Array(D))
    end
    @test_throws SingularException inv(Diagonal(zeros(n)))
    @test_throws SingularException inv(Diagonal([0, 1, 2]))
    @test_throws SingularException inv(Diagonal([0im, 1im, 2im]))
end

@testset "pseudoinverse" begin
    for d in Any[randn(n), zeros(n), Int[], [0, 2, 0.003], [0im, 1+2im, 0.003im], [0//1, 2//1, 3//100], [0//1, 1//1+2im, 3im//100]]
        D = Diagonal(d)
        M = Array(D)
        @test pinv(D) ≈ pinv(M)
        @test pinv(D, 1.0e-2) ≈ pinv(M, 1.0e-2)
    end
end

# allow construct from range
@test all(Diagonal(range(1, stop=3, length=3)) .== Diagonal([1.0,2.0,3.0]))

# Issue 12803
for t in (Float32, Float64, Int, ComplexF64, Rational{Int})
    @test Diagonal(Matrix{t}[fill(t(1), 2, 2), fill(t(1), 3, 3)])[2,1] == zeros(t, 3, 2)
end

# Issue 15401
@test Matrix(1.0I, 5, 5) \ Diagonal(fill(1.,5)) == Matrix(I, 5, 5)

@testset "Triangular and Diagonal" begin
    function _randomarray(type, ::Val{N} = Val(2)) where {N}
        sz = ntuple(_->5, N)
        if type == Int
            return rand(1:9, sz...)
        else
            return randn(type, sz...)
        end
    end
    types = (Float64, Int, ComplexF64)
    for ta in types
        D = Diagonal(_randomarray(ta, Val(1)))
        M = Matrix(D)
        for tb in types
            B = _randomarray(tb, Val(2))
            Tmats = (LowerTriangular(B), UnitLowerTriangular(B), UpperTriangular(B), UnitUpperTriangular(B))
            restypes = (LowerTriangular, LowerTriangular, UpperTriangular, UpperTriangular)
            for (T, rtype) in zip(Tmats, restypes)
                adjtype = (rtype == LowerTriangular) ? UpperTriangular : LowerTriangular

                # Triangular * Diagonal
                R = T * D
                TA = Array(T)
                @test R ≈ TA * M
                @test isa(R, rtype)

                # Diagonal * Triangular
                R = D * T
                @test R ≈ M * TA
                @test isa(R, rtype)

                # Adjoint of Triangular * Diagonal
                R = T' * D
                @test R ≈ TA' * M
                @test isa(R, adjtype)

                # Diagonal * Adjoint of Triangular
                R = D * T'
                @test R ≈ M * TA'
                @test isa(R, adjtype)

                # Transpose of Triangular * Diagonal
                R = transpose(T) * D
                @test R ≈ transpose(TA) * M
                @test isa(R, adjtype)

                # Diagonal * Transpose of Triangular
                R = D * transpose(T)
                @test R ≈ M * transpose(TA)
                @test isa(R, adjtype)
            end
        end
    end
end

let D1 = Diagonal(rand(5)), D2 = Diagonal(rand(5))
    @test LinearAlgebra.rmul!(copy(D1),D2) == D1*D2
    @test LinearAlgebra.lmul!(D1,copy(D2)) == D1*D2
    @test LinearAlgebra.rmul!(copy(D1),transpose(D2)) == D1*transpose(D2)
    @test LinearAlgebra.lmul!(transpose(D1),copy(D2)) == transpose(D1)*D2
    @test LinearAlgebra.rmul!(copy(D1),adjoint(D2)) == D1*adjoint(D2)
    @test LinearAlgebra.lmul!(adjoint(D1),copy(D2)) == adjoint(D1)*D2
end

@testset "multiplication of a Diagonal with a Matrix" begin
    A = collect(reshape(1:8, 4, 2));
    B = BigFloat.(A);
    DL = Diagonal(collect(axes(A, 1)));
    DR = Diagonal(Float16.(collect(axes(A, 2))));

    @test DL * A == collect(DL) * A
    @test A * DR == A * collect(DR)
    @test DL * B == collect(DL) * B
    @test B * DR == B * collect(DR)

    A = reshape([ones(2,2), ones(2,2)*2, ones(2,2)*3, ones(2,2)*4], 2, 2)
    Ac = collect(A)
    D = Diagonal([collect(reshape(1:4, 2, 2)), collect(reshape(5:8, 2, 2))])
    Dc = collect(D)
    @test A * D == Ac * Dc
    @test D * A == Dc * Ac
    @test D * D == Dc * Dc

    AS = similar(A)
    mul!(AS, A, D, true, false)
    @test AS == A * D

    D2 = similar(D)
    mul!(D2, D, D)
    @test D2 == D * D

    copyto!(D2, D)
    lmul!(D, D2)
    @test D2 == D * D
    copyto!(D2, D)
    rmul!(D2, D)
    @test D2 == D * D
end

@testset "multiplication of 2 Diagonal and a Matrix (#46400)" begin
    A = randn(10, 10)
    D = Diagonal(randn(10))
    D2 = Diagonal(randn(10))
    @test D * A * D2 ≈ D * (A * D2)
    @test D * A * D2 ≈ (D * A) * D2
    @test_throws DimensionMismatch Diagonal(ones(9)) * A * D2
    @test_throws DimensionMismatch D * A * Diagonal(ones(9))
end

@testset "multiplication of QR Q-factor and Diagonal (#16615 spot test)" begin
    D = Diagonal(randn(5))
    Q = qr(randn(5, 5)).Q
    @test D * Q' == Array(D) * Q'
    Q = qr(randn(5, 5), ColumnNorm()).Q
    @test_throws ArgumentError lmul!(Q, D)
end

@testset "block diagonal matrices" begin
    D = Diagonal([[1 2; 3 4], [1 2; 3 4]])
    Dherm = Diagonal([[1 1+im; 1-im 1], [1 1+im; 1-im 1]])
    Dsym = Diagonal([[1 1+im; 1+im 1], [1 1+im; 1+im 1]])
    @test adjoint(D) == Diagonal([[1 3; 2 4], [1 3; 2 4]])
    @test transpose(D) == Diagonal([[1 3; 2 4], [1 3; 2 4]])
    @test adjoint(Dherm) == Dherm
    @test transpose(Dherm) == Diagonal([[1 1-im; 1+im 1], [1 1-im; 1+im 1]])
    @test adjoint(Dsym) == Diagonal([[1 1-im; 1-im 1], [1 1-im; 1-im 1]])
    @test transpose(Dsym) == Dsym
    @test diag(D, 0) == diag(D) == [[1 2; 3 4], [1 2; 3 4]]
    @test diag(D, 1) == diag(D, -1) == [zeros(Int,2,2)]
    @test diag(D, 2) == diag(D, -2) == []

    v = [[1, 2], [3, 4]]
    @test Dherm' * v == Dherm * v
    @test transpose(D) * v == [[7, 10], [15, 22]]

    @test issymmetric(D) == false
    @test issymmetric(Dherm) == false
    @test issymmetric(Dsym) == true

    @test ishermitian(D) == false
    @test ishermitian(Dherm) == true
    @test ishermitian(Dsym) == false

    @test exp(D) == Diagonal([exp([1 2; 3 4]), exp([1 2; 3 4])])
    @test cis(D) == Diagonal([cis([1 2; 3 4]), cis([1 2; 3 4])])
    @test log(D) == Diagonal([log([1 2; 3 4]), log([1 2; 3 4])])
    @test sqrt(D) == Diagonal([sqrt([1 2; 3 4]), sqrt([1 2; 3 4])])

    @test tr(D) == 10
    @test det(D) == 4

    M = [1 2; 3 4]
    for n in 0:1
        D = Diagonal(fill(M, n))
        @test D == Matrix{eltype(D)}(D)
    end

    S = SizedArray{(2,3)}(reshape([1:6;],2,3))
    D = Diagonal(fill(S,3))
    @test D * fill(S,2,3)' == fill(S * S', 3, 2)
    @test fill(S,3,2)' * D == fill(S' * S, 2, 3)

    @testset "indexing with non-standard-axes" begin
        s = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
        D = Diagonal(fill(s,3))
        @test @inferred(D[1,2]) isa typeof(s)
        @test all(iszero, D[1,2])
    end

    @testset "mul!" begin
        D1 = Diagonal(fill(ones(2,3), 2))
        D2 = Diagonal(fill(ones(3,2), 2))
        C = similar(D1, size(D1))
        mul!(C, D1, D2)
        @test all(x -> size(x) == (2,2), C)
        @test C == D1 * D2
        D = similar(D1)
        mul!(D, D1, D2)
        @test all(x -> size(x) == (2,2), D)
        @test D == D1 * D2
    end
end

@testset "Eigensystem for block diagonal (issue #30681)" begin
    I2 = Matrix(I, 2,2)
    D = Diagonal([2.0*I2, 3.0*I2])
    eigD = eigen(D)
    evals = [ 2.0, 2.0, 3.0, 3.0 ]
    evecs = [ [[ 1.0, 0.0 ]]  [[ 0.0, 1.0 ]]  [[ 0.0, 0.0 ]]  [[ 0.0, 0.0 ]];
              [[ 0.0, 0.0 ]]  [[ 0.0, 0.0 ]]  [[ 1.0, 0.0 ]]  [[ 0.0, 1.0 ]] ]
    @test eigD.values == evals
    @test eigD.vectors == evecs
    @test D * eigD.vectors ≈ eigD.vectors * Diagonal(eigD.values)

    I3 = Matrix(I, 3,3)
    D = Diagonal([[0.0 -1.0; 1.0 0.0], 2.0*I3])
    eigD = eigen(D)
    evals = [ -1.0im, 1.0im, 2.0, 2.0, 2.0 ]
    evecs = [ [[ 1/sqrt(2)+0im, 1/sqrt(2)*im ]]  [[ 1/sqrt(2)+0im, -1/sqrt(2)*im ]]  [[ 0.0, 0.0 ]]       [[ 0.0, 0.0 ]]      [[ 0.0, 0.0]];
              [[ 0.0, 0.0, 0.0 ]]                [[ 0.0, 0.0, 0.0 ]]                 [[ 1.0, 0.0, 0.0 ]]  [[ 0.0, 1.0, 0.0 ]] [[ 0.0, 0.0, 1.0]] ]
    @test eigD.values == evals
    @test eigD.vectors ≈ evecs
    @test D * eigD.vectors ≈ eigD.vectors * Diagonal(eigD.values)

    # test concrete types
    D = Diagonal([I2 for _ in 1:4])
    @test eigen(D) isa Eigen{Vector{Float64}, Float64, Matrix{Vector{Float64}}, Vector{Float64}}
end

@testset "linear solve for block diagonal matrices" begin
    D = Diagonal([rand(2,2) for _ in 1:5])
    b = [rand(2,2) for _ in 1:5]
    B = [rand(2,2) for _ in 1:5, _ in 1:5]
    @test ldiv!(D, copy(b)) ≈ Diagonal(inv.(D.diag)) * b
    @test ldiv!(D, copy(B)) ≈ Diagonal(inv.(D.diag)) * B
    @test rdiv!(copy(B), D) ≈ B * Diagonal(inv.(D.diag))
end

@testset "multiplication/division with Symmetric/Hermitian" begin
    for T in (Float64, ComplexF64)
        D = Diagonal(randn(T, n))
        A = randn(T, n, n); A = A'A
        S = Symmetric(A)
        H = Hermitian(A)
        for (transform1, transform2) in ((identity,  identity),
                (identity,  adjoint  ), (adjoint,   identity ), (adjoint,   adjoint  ),
                (identity,  transpose), (transpose, identity ), (transpose, transpose) )
            @test *(transform1(D), transform2(S)) ≈ *(transform1(Matrix(D)), transform2(Matrix(S)))
            @test *(transform1(D), transform2(H)) ≈ *(transform1(Matrix(D)), transform2(Matrix(H)))
            @test *(transform1(S), transform2(D)) ≈ *(transform1(Matrix(S)), transform2(Matrix(D)))
            @test *(transform1(S), transform2(H)) ≈ *(transform1(Matrix(S)), transform2(Matrix(H)))
            @test (transform1(H)/D) * D ≈ transform1(H)
            @test (transform1(S)/D) * D ≈ transform1(S)
            @test D * (D\transform2(H)) ≈ transform2(H)
            @test D * (D\transform2(S)) ≈ transform2(S)
        end
    end
end

@testset "multiplication of transposes of Diagonal (#22428)" begin
    for T in (Float64, ComplexF64)
        D = Diagonal(randn(T, 5, 5))
        B = Diagonal(randn(T, 5, 5))
        DD = Diagonal([randn(T, 2, 2), rand(T, 2, 2)])
        BB = Diagonal([randn(T, 2, 2), rand(T, 2, 2)])
        fullDD = copyto!(Matrix{Matrix{T}}(undef, 2, 2), DD)
        fullBB = copyto!(Matrix{Matrix{T}}(undef, 2, 2), BB)
        for (transform1, transform2) in ((identity,  identity),
                (identity,  adjoint  ), (adjoint,   identity ), (adjoint,   adjoint  ),
                (identity,  transpose), (transpose, identity ), (transpose, transpose))
            @test *(transform1(D), transform2(B))::typeof(D) ≈ *(transform1(Matrix(D)), transform2(Matrix(B))) atol=2 * eps()
            @test *(transform1(DD), transform2(BB))::typeof(DD) == *(transform1(fullDD), transform2(fullBB))
        end
        M = randn(T, 5, 5)
        MM = [randn(T, 2, 2) for _ in 1:2, _ in 1:2]
        for transform in (identity, adjoint, transpose)
            @test lmul!(transform(D), copy(M)) ≈ *(transform(Matrix(D)), M)
            @test rmul!(copy(M), transform(D)) ≈ *(M, transform(Matrix(D)))
            @test lmul!(transform(DD), copy(MM)) ≈ *(transform(fullDD), MM)
            @test rmul!(copy(MM), transform(DD)) ≈ *(MM, transform(fullDD))
        end
    end
end

@testset "Diagonal of adjoint/transpose vectors (#23649)" begin
    @test Diagonal(adjoint([1, 2, 3])) == Diagonal([1 2 3])
    @test Diagonal(transpose([1, 2, 3])) == Diagonal([1 2 3])
end

@testset "Multiplication with adjoint and transpose vectors (#26863)" begin
    x = collect(1:2)
    xt = transpose(x)
    A = reshape([[1 2; 3 4], zeros(Int,2,2), zeros(Int, 2, 2), [5 6; 7 8]], 2, 2)
    D = Diagonal(A)
    @test x'*D == x'*A == collect(x')*D == collect(x')*A
    @test xt*D == xt*A == collect(xt)*D == collect(xt)*A
    outadjxD = similar(x'*D); outtrxD = similar(xt*D);
    mul!(outadjxD, x', D)
    @test outadjxD == x'*D
    mul!(outtrxD, xt, D)
    @test outtrxD == xt*D

    D1 = Diagonal([[1 2; 3 4]])
    @test D1 * x' == D1 * collect(x') == collect(D1) * collect(x')
    @test D1 * xt == D1 * collect(xt) == collect(D1) * collect(xt)
    outD1adjx = similar(D1 * x'); outD1trx = similar(D1 * xt);
    mul!(outadjxD, D1, x')
    @test outadjxD == D1*x'
    mul!(outtrxD, D1, xt)
    @test outtrxD == D1*xt

    y = [x, x]
    yt = transpose(y)
    @test y'*D*y == (y'*D)*y == (y'*A)*y
    @test yt*D*y == (yt*D)*y == (yt*A)*y
    outadjyD = similar(y'*D); outtryD = similar(yt*D);
    outadjyD2 = similar(collect(y'*D)); outtryD2 = similar(collect(yt*D));
    mul!(outadjyD, y', D)
    mul!(outadjyD2, y', D)
    @test outadjyD == outadjyD2 == y'*D
    mul!(outtryD, yt, D)
    mul!(outtryD2, yt, D)
    @test outtryD == outtryD2 == yt*D
end

@testset "Multiplication of single element Diagonal (#36746, #40726)" begin
    @test_throws DimensionMismatch Diagonal(randn(1)) * randn(5)
    @test_throws DimensionMismatch Diagonal(randn(1)) * Diagonal(randn(3, 3))
    A = [1 0; 0 2]
    v = [3, 4]
    @test Diagonal(A) * v == A * v
    @test Diagonal(A) * Diagonal(A) == A * A
    @test_throws DimensionMismatch [1 0;0 1] * Diagonal([2 3])   # Issue #40726
    @test_throws DimensionMismatch lmul!(Diagonal([1]), [1,2,3]) # nearby
end

@testset "Multiplication of a Diagonal with an OffsetArray" begin
    # Offset indices should throw
    D = Diagonal(1:4)
    A = OffsetArray(rand(4,4), 2, 2)
    @test_throws ArgumentError D * A
    @test_throws ArgumentError A * D
    @test_throws ArgumentError mul!(similar(A, size(A)), A, D)
    @test_throws ArgumentError mul!(similar(A, size(A)), D, A)
end

@testset "Triangular division by Diagonal #27989" begin
    K = 5
    for elty in (Float32, Float64, ComplexF32, ComplexF64)
        U = UpperTriangular(randn(elty, K, K))
        L = LowerTriangular(randn(elty, K, K))
        D = Diagonal(randn(elty, K))
        @test (U / D)::UpperTriangular{elty} == UpperTriangular(Matrix(U) / Matrix(D))
        @test (L / D)::LowerTriangular{elty} == LowerTriangular(Matrix(L) / Matrix(D))
        @test (D \ U)::UpperTriangular{elty} == UpperTriangular(Matrix(D) \ Matrix(U))
        @test (D \ L)::LowerTriangular{elty} == LowerTriangular(Matrix(D) \ Matrix(L))
    end
end

@testset "(Sym)Tridiagonal division by Diagonal" begin
    for K in (5, 1), elty in (Float64, ComplexF32), overlength in (1, 0)
        S = SymTridiagonal(randn(elty, K), randn(elty, K-overlength))
        T = Tridiagonal(randn(elty, K-1), randn(elty, K), randn(elty, K-1))
        D = Diagonal(randn(elty, K))
        D0 = Diagonal(zeros(elty, K))
        @test (D \ S)::Tridiagonal{elty} == Tridiagonal(Matrix(D) \ Matrix(S))
        @test (D \ T)::Tridiagonal{elty} == Tridiagonal(Matrix(D) \ Matrix(T))
        @test (S / D)::Tridiagonal{elty} == Tridiagonal(Matrix(S) / Matrix(D))
        @test (T / D)::Tridiagonal{elty} == Tridiagonal(Matrix(T) / Matrix(D))
        @test_throws SingularException D0 \ S
        @test_throws SingularException D0 \ T
        @test_throws SingularException S / D0
        @test_throws SingularException T / D0
    end
    # 0-length case
    S = SymTridiagonal(Float64[], Float64[])
    T = Tridiagonal(Float64[], Float64[], Float64[])
    D = Diagonal(Float64[])
    @test (D \ S)::Tridiagonal{Float64} == T
    @test (D \ T)::Tridiagonal{Float64} == T
    @test (S / D)::Tridiagonal{Float64} == T
    @test (T / D)::Tridiagonal{Float64} == T
    # matrix eltype case
    K = 5
    for elty in (Float64, ComplexF32), overlength in (1, 0)
        S = SymTridiagonal([rand(elty, 2, 2) for _ in 1:K], [rand(elty, 2, 2) for _ in 1:K-overlength])
        T = Tridiagonal([rand(elty, 2, 2) for _ in 1:K-1], [rand(elty, 2, 2) for _ in 1:K], [rand(elty, 2, 2) for _ in 1:K-1])
        D = Diagonal(randn(elty, K))
        SM = fill(zeros(elty, 2, 2), K, K)
        TM = copy(SM)
        SM[1,1] = S[1,1]; TM[1,1] = T[1,1]
        for j in 2:K
            SM[j,j-1] = S[j,j-1]; SM[j,j] = S[j,j]; SM[j-1,j] = S[j-1,j]
            TM[j,j-1] = T[j,j-1]; TM[j,j] = T[j,j]; TM[j-1,j] = T[j-1,j]
        end
        for (M, Mm) in ((S, SM), (T, TM))
            DS = D \ M
            @test DS isa Tridiagonal
            DM = D \ Mm
            for i in -1:1; @test diag(DS, i) ≈ diag(DM, i) end
            DS = M / D
            @test DS isa Tridiagonal
            DM = Mm / D
            for i in -1:1; @test diag(DS, i) ≈ diag(DM, i) end
        end
    end
    # eltype promotion case
    S = SymTridiagonal(rand(-20:20, K), rand(-20:20, K-1))
    T = Tridiagonal(rand(-20:20, K-1), rand(-20:20, K), rand(-20:20, K-1))
    D = Diagonal(rand(1:20, K))
    @test (D \ S)::Tridiagonal{Float64} == Tridiagonal(Matrix(D) \ Matrix(S))
    @test (D \ T)::Tridiagonal{Float64} == Tridiagonal(Matrix(D) \ Matrix(T))
    @test (S / D)::Tridiagonal{Float64} == Tridiagonal(Matrix(S) / Matrix(D))
    @test (T / D)::Tridiagonal{Float64} == Tridiagonal(Matrix(T) / Matrix(D))
end

@testset "eigenvalue sorting" begin
    D = Diagonal([0.4, 0.2, -1.3])
    @test eigvals(D) == eigen(D).values == [0.4, 0.2, -1.3] # not sorted by default
    @test eigvals(Matrix(D)) == eigen(Matrix(D)).values == [-1.3, 0.2, 0.4] # sorted even if diagonal special case is detected
    E = eigen(D, sortby=abs) # sortby keyword supported for eigen(::Diagonal)
    @test E.values == [0.2, 0.4, -1.3]
    @test E.vectors == [0 1 0; 1 0 0; 0 0 1]
end

@testset "sum, mapreduce" begin
    D = Diagonal([1,2,3])
    Ddense = Matrix(D)
    @test sum(D) == 6
    @test_throws ArgumentError sum(D, dims=0)
    @test sum(D, dims=1) == sum(Ddense, dims=1)
    @test sum(D, dims=2) == sum(Ddense, dims=2)
    @test sum(D, dims=3) == sum(Ddense, dims=3)
    @test typeof(sum(D, dims=1)) == typeof(sum(Ddense, dims=1))
    @test mapreduce(one, min, D, dims=1) == mapreduce(one, min, Ddense, dims=1)
    @test mapreduce(one, min, D, dims=2) == mapreduce(one, min, Ddense, dims=2)
    @test mapreduce(one, min, D, dims=3) == mapreduce(one, min, Ddense, dims=3)
    @test typeof(mapreduce(one, min, D, dims=1)) == typeof(mapreduce(one, min, Ddense, dims=1))
    @test mapreduce(zero, max, D, dims=1) == mapreduce(zero, max, Ddense, dims=1)
    @test mapreduce(zero, max, D, dims=2) == mapreduce(zero, max, Ddense, dims=2)
    @test mapreduce(zero, max, D, dims=3) == mapreduce(zero, max, Ddense, dims=3)
    @test typeof(mapreduce(zero, max, D, dims=1)) == typeof(mapreduce(zero, max, Ddense, dims=1))

    D = Diagonal(Int[])
    Ddense = Matrix(D)
    @test sum(D) == 0
    @test_throws ArgumentError sum(D, dims=0)
    @test sum(D, dims=1) == sum(Ddense, dims=1)
    @test sum(D, dims=2) == sum(Ddense, dims=2)
    @test sum(D, dims=3) == sum(Ddense, dims=3)
    @test typeof(sum(D, dims=1)) == typeof(sum(Ddense, dims=1))

    D = Diagonal(Int[2])
    Ddense = Matrix(D)
    @test sum(D) == 2
    @test_throws ArgumentError sum(D, dims=0)
    @test sum(D, dims=1) == sum(Ddense, dims=1)
    @test sum(D, dims=2) == sum(Ddense, dims=2)
    @test sum(D, dims=3) == sum(Ddense, dims=3)
    @test typeof(sum(D, dims=1)) == typeof(sum(Ddense, dims=1))
end

@testset "logabsdet for generic eltype" begin
    d = Any[1, -2.0, -3.0]
    D = Diagonal(d)
    d1, s1 = logabsdet(D)
    @test d1 ≈ sum(log ∘ abs, d)
    @test s1 == prod(sign, d)
end

@testset "Empty (#35424) & size checks (#47060)" begin
    @test zeros(0)'*Diagonal(zeros(0))*zeros(0) === 0.0
    @test transpose(zeros(0))*Diagonal(zeros(Complex{Int}, 0))*zeros(0) === 0.0 + 0.0im
    @test dot(zeros(Int32, 0), Diagonal(zeros(Int, 0)), zeros(Int16, 0)) === 0
    @test_throws DimensionMismatch zeros(2)' * Diagonal(zeros(2)) * zeros(3)
    @test_throws DimensionMismatch zeros(3)' * Diagonal(zeros(2)) * zeros(2)
    @test_throws DimensionMismatch dot(zeros(2), Diagonal(zeros(2)), zeros(3))
    @test_throws DimensionMismatch dot(zeros(3), Diagonal(zeros(2)), zeros(2))
end

@testset "Diagonal(undef)" begin
    d = Diagonal{Float32}(undef, 2)
    @test length(d.diag) == 2
end

@testset "permutedims (#39447)" begin
    for D in (Diagonal(zeros(5)), Diagonal(zeros(5) .+ 1im), Diagonal([[1,2],[3,4]]))
        @test permutedims(D) === permutedims(D,(1,2)) === permutedims(D,(2,1)) === D
        @test_throws ArgumentError permutedims(D,(1,3))
    end
end

@testset "Inner product" begin
    A = Diagonal(rand(10) .+ im)
    B = Diagonal(rand(10) .+ im)
    @test dot(A, B) ≈ dot(Matrix(A), B)
    @test dot(A, B) ≈ dot(A, Matrix(B))
    @test dot(A, B) ≈ dot(Matrix(A), Matrix(B))
    @test dot(A, B) ≈ conj(dot(B, A))
end

@testset "eltype relaxation(#41015)" begin
    A = rand(3,3)
    for trans in (identity, adjoint, transpose)
        @test ldiv!(trans(I(3)), A) == A
        @test rdiv!(A, trans(I(3))) == A
    end
end

@testset "Conversion to AbstractArray" begin
    # tests corresponding to #34995
    d = ImmutableArray([1, 2, 3, 4])
    D = Diagonal(d)

    @test convert(AbstractArray{Float64}, D)::Diagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == D
    @test convert(AbstractMatrix{Float64}, D)::Diagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == D
end

@testset "divisions functionality" for elty in (Int, Float64, ComplexF64)
    B = Diagonal(rand(elty,5,5))
    x = rand(elty)
    @test \(x, B) == /(B, x)
end

@testset "promotion" begin
    for (v1, v2) in (([true], [1]), ([zeros(2,2)], [zeros(Int, 2,2)]))
        T = promote_type(eltype(v1), eltype(v2))
        V = promote_type(typeof(v1), typeof(v2))
        d1 = Diagonal(v1)
        d2 = Diagonal(v2)
        v = [d1, d2]
        @test (@inferred eltype(v)) == Diagonal{T, V}
    end
    # test for a type for which promote_type doesn't lead to a concrete eltype
    struct MyArrayWrapper{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
       a :: A
    end
    Base.size(M::MyArrayWrapper) = size(M.a)
    Base.axes(M::MyArrayWrapper) = axes(M.a)
    Base.length(M::MyArrayWrapper) = length(M.a)
    Base.getindex(M::MyArrayWrapper, i::Int...) = M.a[i...]
    Base.setindex!(M::MyArrayWrapper, v, i::Int...) = M.a[i...] = v
    d1 = Diagonal(MyArrayWrapper(1:3))
    d2 = Diagonal(MyArrayWrapper(1.0:3.0))
    c = [d1, d2]
    @test c[1] == d1
    @test c[2] == d2
end

@testset "zero and one" begin
    D1 = Diagonal(rand(3))
    @test D1 + zero(D1) == D1
    @test D1 * one(D1) == D1
    @test D1 * oneunit(D1) == D1
    @test oneunit(D1) isa typeof(D1)
    D2 = Diagonal([collect(reshape(1:4, 2, 2)), collect(reshape(5:8, 2, 2))])
    @test D2 + zero(D2) == D2
    @test D2 * one(D2) == D2
    @test D2 * oneunit(D2) == D2
    @test oneunit(D2) isa typeof(D2)
    D3 = Diagonal([D2, D2]);
    @test D3 + zero(D3) == D3
    @test D3 * one(D3) == D3
    @test D3 * oneunit(D3) == D3
    @test oneunit(D3) isa typeof(D3)
end

@testset "$Tri" for (Tri, UTri) in ((UpperTriangular, UnitUpperTriangular), (LowerTriangular, UnitLowerTriangular))
    A = randn(4, 4)
    TriA = Tri(A)
    UTriA = UTri(A)
    D = Diagonal(1.0:4.0)
    DM = Matrix(D)
    DMF = factorize(DM)
    outTri = similar(TriA)
    out = similar(A)
    # 2 args
    @testset for fun in (*, rmul!, rdiv!, /)
        @test fun(copy(TriA), D)::Tri == fun(Matrix(TriA), D)
        @test fun(copy(UTriA), D)::Tri == fun(Matrix(UTriA), D)
    end
    @testset for fun in (*, lmul!, ldiv!, \)
        @test fun(D, copy(TriA))::Tri == fun(D, Matrix(TriA))
        @test fun(D, copy(UTriA))::Tri == fun(D, Matrix(UTriA))
    end
    # 3 args
    @test outTri === ldiv!(outTri, D, TriA)::Tri == ldiv!(out, D, Matrix(TriA))
    @test outTri === ldiv!(outTri, D, UTriA)::Tri == ldiv!(out, D, Matrix(UTriA))
    @test outTri === mul!(outTri, D, TriA)::Tri == mul!(out, D, Matrix(TriA))
    @test outTri === mul!(outTri, D, UTriA)::Tri == mul!(out, D, Matrix(UTriA))
    @test outTri === mul!(outTri, TriA, D)::Tri == mul!(out, Matrix(TriA), D)
    @test outTri === mul!(outTri, UTriA, D)::Tri == mul!(out, Matrix(UTriA), D)
    # 5 args
    @test outTri === mul!(outTri, D, TriA, 2, 1)::Tri == mul!(out, D, Matrix(TriA), 2, 1)
    @test outTri === mul!(outTri, D, UTriA, 2, 1)::Tri == mul!(out, D, Matrix(UTriA), 2, 1)
    @test outTri === mul!(outTri, TriA, D, 2, 1)::Tri == mul!(out, Matrix(TriA), D, 2, 1)
    @test outTri === mul!(outTri, UTriA, D, 2, 1)::Tri == mul!(out, Matrix(UTriA), D, 2, 1)

    # we may write to a Unit triangular if the diagonal is preserved
    ID = Diagonal(ones(size(UTriA,2)))
    @test mul!(copy(UTriA), UTriA, ID) == UTriA
    @test mul!(copy(UTriA), ID, UTriA) == UTriA

    @testset "partly filled parents" begin
        M = Matrix{BigFloat}(undef, 2, 2)
        M[1,1] = M[2,2] = 3
        isupper = Tri == UpperTriangular
        M[1+!isupper, 1+isupper] = 3
        D = Diagonal(1:2)
        T = Tri(M)
        TA = Array(T)
        @test T * D == TA * D
        @test D * T == D * TA
        @test mul!(copy(T), T, D, 2, 3) == 2T * D + 3T
        @test mul!(copy(T), D, T, 2, 3) == 2D * T + 3T

        U = UTri(M)
        UA = Array(U)
        @test U * D == UA * D
        @test D * U == D * UA
        @test mul!(copy(T), U, D, 2, 3) == 2 * UA * D + 3TA
        @test mul!(copy(T), D, U, 2, 3) == 2 * D * UA + 3TA

        M2 = Matrix{BigFloat}(undef, 2, 2)
        M2[1+!isupper, 1+isupper] = 3
        U = UTri(M2)
        UA = Array(U)
        @test U * D == UA * D
        @test D * U == D * UA
        ID = Diagonal(ones(size(U,2)))
        @test mul!(copy(U), U, ID) == U
        @test mul!(copy(U), ID, U) == U
        @test mul!(copy(U), U, ID, 2, -1) == U
        @test mul!(copy(U), ID, U, 2, -1) == U
    end
end

@testset "rmul!/lmul! for adj/trans" begin
    for T in (Float64, ComplexF64)
        A = rand(T,5,4); B = similar(A)
        for f in (adjoint, transpose)
            D = Diagonal(rand(T, size(A,1)))
            B .= A
            rmul!(f(B), D)
            @test f(B) == f(A) * D
            D = Diagonal(rand(T, size(A,2)))
            B .= A
            lmul!(D, f(B))
            @test f(B) == D * f(A)
        end
    end
end

struct SMatrix1{T} <: AbstractArray{T,2}
    elt::T
end
Base.:(==)(A::SMatrix1, B::SMatrix1) = A.elt == B.elt
Base.zero(::Type{SMatrix1{T}}) where {T} = SMatrix1(zero(T))
Base.iszero(A::SMatrix1) = iszero(A.elt)
Base.getindex(A::SMatrix1, inds...) = A.elt
Base.size(::SMatrix1) = (1, 1)
@testset "map for Diagonal matrices (#46292)" begin
    A = Diagonal([1])
    @test A isa Diagonal{Int,Vector{Int}}
    @test 2*A isa Diagonal{Int,Vector{Int}}
    @test A.+1 isa Matrix{Int}
    # Numeric element types remain diagonal
    B = map(SMatrix1, A)
    @test B == fill(SMatrix1(1), 1, 1)
    @test B isa Diagonal{SMatrix1{Int},Vector{SMatrix1{Int}}}
    # Non-numeric element types become dense
    C = map(a -> SMatrix1(string(a)), A)
    @test C == fill(SMatrix1(string(1)), 1, 1)
    @test C isa Matrix{SMatrix1{String}}
end

@testset "show" begin
    @test repr(Diagonal([1,2])) == "Diagonal([1, 2])"  # 2-arg show
    @test contains(repr(MIME"text/plain"(), Diagonal([1,2])), "⋅  2")  # 3-arg show
end

@testset "copyto! with UniformScaling" begin
    @testset "Fill" begin
        for len in (4, InfiniteArrays.Infinity())
            d = FillArrays.Fill(1, len)
            D = Diagonal(d)
            @test copyto!(D, I) === D
        end
    end
    D = Diagonal(fill(2, 2))
    copyto!(D, I)
    @test all(isone, diag(D))
end

@testset "diagonal triple multiplication (#49005)" begin
    local n = 10
    @test *(Diagonal(ones(n)), Diagonal(1:n), Diagonal(ones(n))) isa Diagonal
    @test_throws DimensionMismatch (*(Diagonal(ones(n)), Diagonal(1:n), Diagonal(ones(n+1))))
    @test_throws DimensionMismatch (*(Diagonal(ones(n)), Diagonal(1:n+1), Diagonal(ones(n+1))))
    @test_throws DimensionMismatch (*(Diagonal(ones(n+1)), Diagonal(1:n), Diagonal(ones(n))))

    # currently falls back to two-term *
    @test *(Diagonal(ones(n)), Diagonal(1:n), Diagonal(ones(n)), Diagonal(1:n)) isa Diagonal
end

@testset "triple multiplication with a sandwiched BandedMatrix" begin
    D = Diagonal(StepRangeLen(NaN, 0, 4));
    B = Bidiagonal(1:4, 1:3, :U)
    C = D * B * D
    @test iszero(diag(C, 2))
    # test associativity
    C1 = (D * B) * D
    C2 = D * (B * D)
    @test diag(C,2) == diag(C1,2) == diag(C2,2)
end

@testset "diagind" begin
    D = Diagonal(1:4)
    M = Matrix(D)
    @testset for k in -4:4
        @test D[diagind(D,k)] == M[diagind(M,k)]
    end
end

@testset "avoid matmul ambiguities with ::MyMatrix * ::AbstractMatrix" begin
    A = [i+j for i in 1:2, j in 1:2]
    S = SizedArrays.SizedArray{(2,2)}(A)
    D = Diagonal([1:2;])
    @test S * D == A * D
    @test D * S == D * A
    C1, C2 = zeros(2,2), zeros(2,2)
    @test mul!(C1, S, D) == mul!(C2, A, D)
    @test mul!(C1, S, D, 1, 2) == mul!(C2, A, D, 1 ,2)
    @test mul!(C1, D, S) == mul!(C2, D, A)
    @test mul!(C1, D, S, 1, 2) == mul!(C2, D, A, 1 ,2)

    v = [i for i in 1:2]
    sv = SizedArrays.SizedArray{(2,)}(v)
    @test D * sv == D * v
    C1, C2 = zeros(2), zeros(2)
    @test mul!(C1, D, sv) == mul!(C2, D, v)
    @test mul!(C1, D, sv, 1, 2) == mul!(C2, D, v, 1 ,2)
end

@testset "copy" begin
    @test copy(Diagonal(1:5)) === Diagonal(1:5)
end

@testset "kron! for Diagonal" begin
    a = Diagonal([1, 2])
    b = Diagonal([3, 4])
    # Diagonal out
    c = Diagonal([0, 0, 0, 0])
    kron!(c, b, a)
    @test c == Diagonal([3, 6, 4, 8])
    @test c == kron!(fill(0, 4, 4), Matrix(b), Matrix(a)) # against dense kron!
    c = Diagonal(Vector{Float64}(undef, 4))
    kron!(c, a, b)
    @test c == Diagonal([3.0, 4.0, 6.0, 8.0])

    # AbstractArray out
    c = fill(0, 4, 4)
    kron!(c, b, a) 
    @test c == diagm([3, 6, 4, 8])
    @test c == kron!(fill(0, 4, 4), Matrix(b), Matrix(a)) # against dense kron!
    c = Matrix{Float64}(undef, 4, 4)
    kron!(c, a, b)
    @test c == diagm([3.0, 4.0, 6.0, 8.0])
    @test_throws DimensionMismatch kron!(Diagonal(zeros(5)), Diagonal(zeros(2)), Diagonal(zeros(2)))
end

@testset "uppertriangular/lowertriangular" begin
    D = Diagonal([1,2])
    @test LinearAlgebra.uppertriangular(D) === D
    @test LinearAlgebra.lowertriangular(D) === D
end

@testset "mul/div with an adjoint vector" begin
    A = [1.0;;]
    x = [1.0]
    yadj = Diagonal(A) \ x'
    @test typeof(yadj) == typeof(x')
    @test yadj == x'
    yadj = Diagonal(A) * x'
    @test typeof(yadj) == typeof(x')
    @test yadj == x'
end

@testset "Matrix conversion for non-numeric" begin
    D = Diagonal(fill(Diagonal([1,3]), 2))
    M = Matrix{eltype(D)}(D)
    @test M isa Matrix{eltype(D)}
    @test M == D
end

@testset "rmul!/lmul! with banded matrices" begin
    @testset "$(nameof(typeof(B)))" for B in (
                            Bidiagonal(rand(4), rand(3), :L),
                            Tridiagonal(rand(3), rand(4), rand(3))
                    )
        BA = Array(B)
        D = Diagonal(rand(size(B,1)))
        DA = Array(D)
        @test rmul!(copy(B), D) ≈ B * D ≈ BA * DA
        @test lmul!(D, copy(B)) ≈ D * B ≈ DA * BA
    end
end

@testset "rmul!/lmul! with numbers" begin
    D = Diagonal(rand(4))
    @test rmul!(copy(D), 0.2) ≈ rmul!(Array(D), 0.2)
    @test lmul!(0.2, copy(D)) ≈ lmul!(0.2, Array(D))
    @test_throws ArgumentError rmul!(D, NaN)
    @test_throws ArgumentError lmul!(NaN, D)
    D = Diagonal(rand(1))
    @test all(isnan, rmul!(copy(D), NaN))
    @test all(isnan, lmul!(NaN, copy(D)))
end

@testset "+/- with block Symmetric/Hermitian" begin
    for p in ([1 2; 3 4], [1 2+im; 2-im 4+2im])
        m = SizedArrays.SizedArray{(2,2)}(p)
        D = Diagonal(fill(m, 2))
        M = Matrix(D)
        for T in (Symmetric, Hermitian)
            S = T(fill(m, 2, 2))
            SA = Array(S)
            @test D + S == M + SA
            @test S + D == SA + M
        end
    end
end

@testset "bounds-check with CartesianIndex ranges" begin
    D = Diagonal(1:typemax(Int))
    @test checkbounds(Bool, D, diagind(D, IndexCartesian()))
end

@testset "zeros in kron with block matrices" begin
    D = Diagonal(1:4)
    M = Matrix(D)
    B = reshape([ones(2,2), ones(3,2), ones(2,3), ones(3,3)], 2, 2)
    @test kron(D, B) == kron(M, B)
    @test kron(B, D) == kron(B, M)
    D2 = Diagonal([ones(2,2), ones(3,3)])
    M2 = Array{eltype(D2)}(D2)
    @test kron(D, D2) == kron(D, M2)
    @test kron(D2, D) == kron(M2, D)
end

@testset "opnorms" begin
    D = Diagonal([1,-2,3,-4])

    @test opnorm(D, 1) == opnorm(Matrix(D), 1)
    @test opnorm(D, 2) ≈ opnorm(Matrix(D), 2)
    @test opnorm(D, Inf) == opnorm(Matrix(D), Inf)

    D = Diagonal([-11])
    @test opnorm(D, 1) == opnorm(Matrix(D), 1)
    @test opnorm(D, 2) ≈ opnorm(Matrix(D), 2)
    @test opnorm(D, Inf) == opnorm(Matrix(D), Inf)

    # block diagonal matrices
    D = Diagonal([[1 2; 3 4], [5 6; 7 8]])
    A = [1 2 0 0; 3 4 0 0; 0 0 5 6; 0 0 7 8] # full matrix of D
    @test opnorm(D, 1) == opnorm(A, 1)
    @test opnorm(D, 2) ≈ opnorm(A, 2)
    @test opnorm(D, Inf) == opnorm(A, Inf)
end

@testset "issymmetric with NaN" begin
    D = Diagonal(fill(NaN,3))
    A = Array(D)
    @test issymmetric(D) == issymmetric(A)
    @test ishermitian(D) == ishermitian(A)
end

@testset "isreal" begin
    D = Diagonal(ones(2))
    @test @inferred((D -> Val(isreal(D)))(D)) == Val(true)
    D = complex.(D)
    @test isreal(D)
    @test !isreal(im*D)
end

@testset "setindex! with BandIndex" begin
    D = Diagonal(zeros(2))
    D[LinearAlgebra.BandIndex(0,2)] = 1
    @test D[2,2] == 1
    @test_throws "cannot set off-diagonal entry $((1,2))" D[LinearAlgebra.BandIndex(1,1)] = 1
    @test_throws BoundsError D[LinearAlgebra.BandIndex(size(D,1),1)]
    @test_throws BoundsError D[LinearAlgebra.BandIndex(0,size(D,1)+1)]
end

@testset "lazy adjtrans" begin
    D = Diagonal(fill([1 2; 3 4], 3))
    m = [2 4; 6 8]
    for op in (transpose, adjoint)
        C = op(D)
        el = op(m)
        C[1,1] = el
        @test D[1,1] == m
        @test (@allocated op(D)) == 0
        @test (@allocated op(op(D))) == 0
    end
end

@testset "fillband!" begin
    D = Diagonal(zeros(4))
    LinearAlgebra.fillband!(D, 2, 0, 0)
    @test all(==(2), diagview(D,0))
    @test all(==(0), diagview(D,-1))
    @test_throws ArgumentError LinearAlgebra.fillband!(D, 3, -2, 2)

    LinearAlgebra.fillstored!(D, 1)
    LinearAlgebra.fillband!(D, 0, -3, 3)
    @test iszero(D)
    LinearAlgebra.fillstored!(D, 1)
    LinearAlgebra.fillband!(D, 0, -10, 10)
    @test iszero(D)

    LinearAlgebra.fillstored!(D, 1)
    D2 = copy(D)
    LinearAlgebra.fillband!(D, 0, -1, -3)
    @test D == D2
    LinearAlgebra.fillband!(D, 0, 10, 10)
    @test D == D2
end

end # module TestDiagonal
