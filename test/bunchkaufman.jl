# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestBunchKaufman

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: BlasComplex, BlasFloat, BlasReal, QRPivoted
using Base: getproperty

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers: COMPREHENSIVE
using Main.LinearAlgebraTestHelpers.Quaternions

n = 10

# Split n into 2 parts for tests needing two matrices
n1 = div(n, 2)
n2 = 2*n1

Random.seed!(12343212)

areal = randn(n,n)/2
aimg  = randn(n,n)/2
a2real = randn(n,n)/2
a2img  = randn(n,n)/2
breal = randn(n,2)/2
bimg  = randn(n,2)/2

areint = rand(1:7, n, n)
aimint = rand(1:7, n, n)
a2reint = rand(1:7, n, n)
a2imint = rand(1:7, n, n)
breint = rand(1:5, n, 2)
bimint = rand(1:5, n, 2)

# Standard mode: BLAS types, Int, generic real and complex floats, and exact rationals
@testset "$eltya argument A" for eltya in (Float64, ComplexF64, Int, BigFloat, Complex{BigFloat},
    Rational{BigInt}, (@static COMPREHENSIVE ? (Float32, ComplexF32, Float16, Complex{Float16},
    Complex{Int}, BigInt, Complex{BigInt}, Complex{Rational{BigInt}}) : ())...)
    # a = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex.(areal, aimg) : areal)
    # a2 = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex.(a2real, a2img) : a2real)
    a = convert(Matrix{eltya}, eltya <: Complex ? (real(eltya) <: AbstractFloat ?
        complex.(areal, aimg) : complex.(areint, aimint)) : (eltya <: AbstractFloat ?
        areal : areint))
    a2 = convert(Matrix{eltya}, eltya <: Complex ? (real(eltya) <: AbstractFloat ?
        complex.(a2real, a2img) : complex.(a2reint, a2imint)) : (eltya <: AbstractFloat ?
        a2real : a2reint))
    asym = transpose(a) + a                  # symmetric indefinite
    aher = a' + a                  # Hermitian indefinite
    apd  = a' * a                  # Positive-definite
    for (a, a2, aher, apd) in ((a, a2, aher, apd),
                               (view(a, 1:n, 1:n),
                                view(a2, 1:n, 1:n),
                                view(aher, 1:n, 1:n),
                                view(apd , 1:n, 1:n)))
        ε = εa = eps(abs(float(one(eltya))))

        # Inertia tests
        @testset "$uplo Bunch-Kaufman factor inertia" for uplo in (:L, :U)
            @testset "rook pivoting: $rook" for rook in (false, true)
                test_list = eltya <: Complex ? (Hermitian(aher, uplo), Hermitian(apd, uplo)) :
                    (Symmetric(transpose(a) + a, uplo), Hermitian(aher, uplo),
                Hermitian(apd, uplo))
                ελ = n*max(eps(Float64), εa) # zero-eigenvalue threshold
                ελ = typeof(Integer(one(real(eltya)))) <: Signed ? Rational{BigInt}(ελ) :
                    real(eltya(ελ))
                for M in test_list
                    bc = bunchkaufman(M, rook)
                    D = bc.D
                    λ = real(eltya <: Complex ? eigen(ComplexF64.(D)).values :
                        eigen(Float64.(D)).values)
                    σ₁ = norm(λ, Inf)
                    np = sum(λ .> ελ*σ₁)
                    nn = sum(λ .< -ελ*σ₁)
                    nz = n - np - nn
                    if real(eltya) <: AbstractFloat
                        @test inertia(bc) == (np, nn, nz)
                    else
                        @test inertia(bc; rtol=ελ) == (np, nn, nz)
                    end
                end
            end
        end

        # check that factorize gives a Bunch-Kaufman
        if eltya <: Union{Float32, Float64, ComplexF32, ComplexF64, Int}
            # Default behaviour only uses Bunch-Kaufman for these types, for now.
            @test isa(factorize(asym), LinearAlgebra.BunchKaufman)
            @test isa(factorize(aher), LinearAlgebra.BunchKaufman)
        end
        @testset "$uplo Bunch-Kaufman factor of indefinite matrix" for uplo in (:L, :U)
            bc1 = bunchkaufman(Hermitian(aher, uplo))
            @test LinearAlgebra.issuccess(bc1)
            @test logabsdet(bc1)[1] ≈ log(abs(det(bc1)))
            if eltya <: Real
                @test logabsdet(bc1)[2] == sign(det(bc1))
            else
                @test logabsdet(bc1)[2] ≈ sign(det(bc1))
            end
            @test inv(bc1)*aher ≈ Matrix(I, n, n)
            @testset for rook in (false, true)
                @test inv(bunchkaufman(Symmetric(transpose(a) + a, uplo), rook))*(transpose(a) + a) ≈ Matrix(I, n, n)
                if eltya <: BlasFloat
                    # test also bunchkaufman! without explicit type tag
                    # no bunchkaufman! method for Int ... yet
                    @test inv(bunchkaufman!(transpose(a) + a, rook))*(transpose(a) + a) ≈ Matrix(I, n, n)
                end
                @test size(bc1) == size(bc1.LD)
                @test size(bc1, 1) == size(bc1.LD, 1)
                @test size(bc1, 2) == size(bc1.LD, 2)
                if eltya <: BlasReal
                    @test_throws ArgumentError bunchkaufman(a)
                end
                # Test extraction of factors
                if eltya <: Real
                    @test getproperty(bc1, uplo)*bc1.D*getproperty(bc1, uplo)' ≈ aher[bc1.p, bc1.p]
                    @test getproperty(bc1, uplo)*bc1.D*getproperty(bc1, uplo)' ≈ bc1.P*aher*bc1.P'
                end

                bc1 = bunchkaufman(Symmetric(asym, uplo))
                @test getproperty(bc1, uplo)*bc1.D*transpose(getproperty(bc1, uplo)) ≈ asym[bc1.p, bc1.p]
                @test getproperty(bc1, uplo)*bc1.D*transpose(getproperty(bc1, uplo)) ≈ bc1.P*asym*transpose(bc1.P)
                @test_throws FieldError bc1.Z
                @test_throws ArgumentError uplo === :L ? bc1.U : bc1.L
            end
            # test Base.iterate
            ref_objs = (bc1.D, uplo === :L ? bc1.L : bc1.U, bc1.p)
            for (bki, bkobj) in enumerate(bc1)
                @test bkobj == ref_objs[bki]
            end
            if eltya <: BlasFloat
                @test convert(LinearAlgebra.BunchKaufman{eltya}, bc1) === bc1
                @test convert(LinearAlgebra.Factorization{eltya}, bc1) === bc1
                if eltya <: BlasReal
                    @test convert(LinearAlgebra.Factorization{Float16}, bc1) == convert(LinearAlgebra.BunchKaufman{Float16}, bc1)
                elseif eltya <: BlasComplex
                    @test convert(LinearAlgebra.Factorization{ComplexF16}, bc1) == convert(LinearAlgebra.BunchKaufman{ComplexF16}, bc1)
                end
            end
            @test Base.propertynames(bc1) == (:p, :P, :L, :U, :D)
        end

        @testset "$eltyb argument B" for eltyb in (Float64, ComplexF64, Int, BigFloat,
            (@static COMPREHENSIVE ? (Float32, ComplexF32, Float16, Complex{Float16}, Complex{BigFloat},
            Complex{Int}, BigInt, Complex{BigInt}, Rational{BigInt}, Complex{Rational{BigInt}}) : ())...)
            # b = eltyb == Int ? rand(1:5, n, 2) : convert(Matrix{eltyb}, eltyb <: Complex ? complex.(breal, bimg) : breal)
            b = convert(Matrix{eltyb}, eltyb <: Complex ? (real(eltyb) <: AbstractFloat ?
                complex.(breal, bimg) : complex.(breint, bimint)) : (eltyb <: AbstractFloat ?
                breal : breint))
            for b in (b, view(b, 1:n, 1:2))
                εb = eps(abs(float(one(eltyb))))
                ε = max(εa,εb)
                epsc = eltya <: Complex ? sqrt(2)*n : n # tolerance scale

                @testset "$uplo Bunch-Kaufman factor of indefinite matrix" for uplo in (:L, :U)
                    bc1 = bunchkaufman(Hermitian(aher, uplo))
                    # @test aher*(bc1\b) ≈ b atol=1000ε
                    cda = eltya <: Complex ? cond(ComplexF64.(aher)) : cond(Float64.(aher))
                    cda = real(eltya) <: AbstractFloat ? real(eltya(cda)) : cda
                    @test norm(aher*(bc1\b) - b) <= epsc*sqrt(eps(cda))*max(
                            norm(aher*(bc1\b)), norm(b))
                end

                @testset "$uplo Bunch-Kaufman factors of a pos-def matrix" for uplo in (:U, :L)
                    @testset "rook pivoting: $rook" for rook in (false, true)
                        bc2 = bunchkaufman(Hermitian(apd, uplo), rook)
                        @test LinearAlgebra.issuccess(bc2)
                        bks = split(sprint(show, "text/plain", bc2), "\n")
                        @test bks[1] == summary(bc2)
                        @test bks[2] == "D factor:"
                        @test bks[4+n] == "$uplo factor:"
                        @test bks[6+2n] == "permutation:"
                        @test logdet(bc2) ≈ log(det(bc2))
                        @test logabsdet(bc2)[1] ≈ log(abs(det(bc2)))
                        @test logabsdet(bc2)[2] == sign(det(bc2))
                        # @test inv(bc2)*apd ≈ Matrix(I, n, n) rtol=Base.rtoldefault(real(eltya))
                        # @test apd*(bc2\b) ≈ b rtol=eps(cond(apd))
                        @test norm(inv(bc2)*apd - Matrix(I, n, n)) <= epsc*Base.rtoldefault(
                            real(eltya))*max(norm(inv(bc2)*apd), norm(Matrix(I, n, n)))
                        cda = eltya <: Complex ? cond(ComplexF64.(apd)) : cond(Float64.(apd))
                        cda = real(eltya) <: AbstractFloat ? real(eltya(cda)) : cda
                        @test norm(apd*(bc2\b) - b) <= epsc*sqrt(eps(cda))*max(
                            norm(apd*(bc2\b)), norm(b))
                        @test ishermitian(bc2)
                        @test !issymmetric(bc2) || eltya <: Real
                    end
                end
            end
        end
    end
end

@testset "Singular matrices" begin
    R = Float64[1 0; 0 0]
    C = ComplexF64[1 0; 0 0]
    for A in (R, Symmetric(R), C, Hermitian(C))
        @test_throws SingularException bunchkaufman(A)
        @test_throws SingularException bunchkaufman!(copy(A))
        @test_throws SingularException bunchkaufman(A; check = true)
        @test_throws SingularException bunchkaufman!(copy(A); check = true)
        @test !issuccess(bunchkaufman(A; check = false))
        @test !issuccess(bunchkaufman!(copy(A); check = false))
    end
    F = bunchkaufman(R; check = false)
    @test sprint(show, "text/plain", F) == "Failed factorization of type $(typeof(F))"
end

@testset "test example due to @timholy in PR 15354" begin
    A = rand(6,5); A = complex(A'*A) # to avoid calling the real-lhs-complex-rhs method
    F = cholesky(A);
    v6 = rand(ComplexF64, 6)
    v5 = view(v6, 1:5)
    @test F\v5 == F\v6[1:5]
end

@testset "issue #32080" begin
    A = Symmetric([-5 -9 9; -9 4 1; 9 1 2])
    B = bunchkaufman(A, true)
    @test B.U * B.D * B.U' ≈ A[B.p, B.p]
end

@test_throws DomainError logdet(bunchkaufman([-1 -1; -1 1]))
@test logabsdet(bunchkaufman([8 4; 4 2]; check = false))[1] == -Inf

@testset "0x0 matrix" begin
    for ul in (:U, :L)
        B = bunchkaufman(Symmetric(ones(0, 0), ul))
        @test isa(B, BunchKaufman)
        @test B.D == Tridiagonal([], [], [])
        @test B.P == ones(0, 0)
        @test B.p == []
        if ul === :U
            @test B.U == UnitUpperTriangular(ones(0, 0))
            @test_throws ArgumentError B.L
        else
            @test B.L == UnitLowerTriangular(ones(0, 0))
            @test_throws ArgumentError B.U
        end
    end
end

@testset "adjoint of BunchKaufman" begin
    Ar = randn(5, 5)
    Ar = Ar + Ar'
    Actmp = complex.(randn(5, 5), randn(5, 5))
    Ac1 = Actmp + Actmp'
    Ac2 = Actmp + transpose(Actmp)
    b = ones(size(Ar, 1))

    F = bunchkaufman(Ar)
    @test F\b == F'\b

    F = bunchkaufman(Ac1)
    @test F\b == F'\b

    F = bunchkaufman(Ac2)
    @test_throws ArgumentError("adjoint not implemented for complex symmetric matrices") F'
end

@testset "BunchKaufman for AbstractMatrix" begin
    S = SymTridiagonal(fill(2.0, 4), ones(3))
    B = bunchkaufman(S)
    @test B.U * B.D * B.U' ≈ S
end

@testset "BunchKaufman array constructors #1461" begin
    a = randn(5,5)
    A = a'a
    for ul in (:U, :L)
        B = bunchkaufman(Symmetric(A, ul))
        @test A ≈ Array(B) ≈ Matrix(B) ≈ AbstractArray(B) ≈ AbstractMatrix(B)
    end
end

@testset "generic Bunch-Kaufman solve" begin
    for T in (BigFloat, Rational{BigInt}, Complex{Rational{BigInt}}), uplo in (:U, :L), rook in (false, true)
        A = T.([0 1 2 0; 1 0 3 1; 2 3 1 0; 0 1 0 2])
        F = bunchkaufman(Symmetric(A, uplo), rook)
        # ldiv! returns its argument, also for a vector
        b = T.([1, 2, 3, 4])
        x = copy(b)
        @test ldiv!(F, x) === x
        @test A * x == b || A * x ≈ b
        X = T.([1 2; 3 4; 5 6; 7 8])
        Y = copy(X)
        @test ldiv!(F, Y) === Y
        @test A * Y == X || A * Y ≈ X
        @test_throws DimensionMismatch F \ T.([1, 2, 3])
        @test_throws DimensionMismatch ldiv!(F, zeros(T, 5, 2))
    end
end

@testset "matrices with an abstract element type" begin
    # the element type is determined by the stored values
    for (A, T) in ((Union{Float32,Float64}[2.0 1.0; 1.0 1.5], Float64),
                   (AbstractFloat[2.0 1.0; 1.0 1.5], Float64),
                   (Real[2 1; 1 1.5], Float64),
                   (Real[2//1 1//2; 1//2 3//2], Rational{Int}),
                   (Hermitian(Number[2 im; -im 3]), Complex{Rational{BigInt}}))
        F = bunchkaufman(A)
        @test F isa BunchKaufman{T}
        @test Matrix(A) * (F \ [1, 2]) ≈ [1, 2]
    end
    # without stored values, the element type is promoted as by division
    @test bunchkaufman(Matrix{Real}(undef, 0, 0)) isa BunchKaufman{Float64}
end

@testset "generic Bunch-Kaufman for plain matrices" begin
    # a complex symmetric, non-Hermitian matrix is factorized as symmetric
    A = Complex{BigFloat}.([1 2im 0; 2im 3 1; 0 1 2im])
    F = bunchkaufman(A)
    @test issymmetric(F)
    @test F.U * F.D * transpose(F.U) ≈ A[F.p, F.p]
    @test_throws ArgumentError bunchkaufman(BigFloat.([1 2; 3 4]))
end

@testset "rook pivoting with a pivot below the safe minimum" begin
    x = Float16(1e-6) # a subnormal number
    for uplo in (:U, :L)
        A = Symmetric(Float16[x x 0; x 1 1; 0 1 3], uplo)
        F = bunchkaufman(A, true)
        @test F.ipiv == 1:3
        M = uplo === :U ? F.U * F.D * F.U' : F.L * F.D * F.L'
        @test M ≈ Matrix(A)[F.p, F.p]
    end
end

@testset "inertia with zero eigenvalues" begin
    # 2x2 block with two eigenvalues below the tolerance
    @test inertia(bunchkaufman(Symmetric([0 1e-3; 1e-3 0])); atol=1.0) == (0, 0, 2)
    # 2x2 block with one eigenvalue below the tolerance, whose sign is that of the trace
    F = bunchkaufman(Symmetric([0.5 1; 1 0.1]))
    @test F.ipiv == [-1, -1]
    @test inertia(F; atol=0.7) == (1, 0, 1)
    @test inertia(bunchkaufman(Symmetric(-[0.5 1; 1 0.1])); atol=0.7) == (0, 1, 1)
    # 1x1 block below the default tolerance
    @test inertia(bunchkaufman(Symmetric([1 0; 0 1e-20]))) == (1, 0, 1)
    # only the number of zero eigenvalues is defined for complex symmetric matrices
    @test inertia(bunchkaufman(Symmetric(ComplexF64[1 0; 0 1e-20]))) == (nothing, nothing, 1)
end

@testset "right division (#1032)" begin
    # the example from the issue
    A = rand(3, 3)
    A = A' * A
    @test A / bunchkaufman(A) ≈ I
    for T in (Float64, ComplexF64, BigFloat, Complex{BigFloat}, Rational{BigInt}, Complex{Rational{BigInt}}),
            W in (Symmetric, Hermitian), uplo in (:U, :L), rook in (false, true)
        T <: Real && W === Hermitian && continue
        M = T.([0 1 2 0; 1 0 3 1; 2 3 1 0; 0 1 0 2] .+ (T <: Complex ? im .* [0 1 0 2; 1 0 1 0; 0 1 0 1; 2 0 1 0] : 0))
        M = W === Hermitian ? M + M' : M + transpose(M)
        A = W(M, uplo)
        F = bunchkaufman(A, rook)
        B = T.([1 2 0 1; 3 4 1 0; 5 6 2 2])
        X = B / F
        @test X * Matrix(A) == B || X * Matrix(A) ≈ B
        @test X == B / Matrix(A) || X ≈ B / Matrix(A)
        Y = copy(B)
        @test rdiv!(Y, F) === Y
        @test Y == X
        b = T.([1, 2, 3, 4])
        # a row vector gives a 1×n matrix, as for other factorizations
        @test size(b' / F) == (1, 4)
        @test (b' / F) * Matrix(A) == b' || (b' / F) * Matrix(A) ≈ b'
        @test transpose(b) / F * Matrix(A) == transpose(b) || transpose(b) / F * Matrix(A) ≈ transpose(b)
        @test_throws DimensionMismatch T.(ones(2, 3)) / F
    end
    # complex right-hand side with a real factorization
    A = Symmetric([2.0 1 0; 1 0 3; 0 3 1])
    B = [1.0+2im 0 1; 3 4im 1]
    @test B / bunchkaufman(A) ≈ B / Matrix(A)
end

@testset "sign of logabsdet" begin
    M = [0 1+2im 3; 1-2im 0 1im; 3 -1im 0]
    for T in (ComplexF64, Complex{BigFloat}), uplo in (:U, :L), rook in (false, true)
        A = T.(M)
        # the sign has the element type, and for a Hermitian matrix it is real-valued, also with
        # 2-by-2 pivots
        F = bunchkaufman(Hermitian(A, uplo), rook)
        @test any(<(0), F.ipiv)
        l, s = @inferred logabsdet(F)
        @test s isa T && s == -1
        @test l ≈ log(12)
        S = Symmetric(A, uplo)
        l, s = @inferred logabsdet(bunchkaufman(S, rook))
        @test l ≈ log(abs(det(Matrix(S))))
        @test s ≈ sign(det(Matrix(S)))
        # a singular matrix
        F0 = bunchkaufman(Hermitian(zeros(T, 2, 2), uplo), rook; check=false)
        @test isequal(@inferred(logabsdet(F0)), (real(T)(-Inf), zero(T)))
    end
    # the logarithm of a rational determinant is a floating-point number, and so is the sign of a
    # complex rational number
    for (T, TS) in ((Rational{BigInt}, Rational{BigInt}), (Complex{Rational{BigInt}}, Complex{BigFloat}))
        F = bunchkaufman(Hermitian(T.([2 1; 1 -3])))
        l, s = @inferred logabsdet(F)
        @test l isa BigFloat && l ≈ log(7)
        @test s isa TS && s == -1
        F0 = bunchkaufman(Hermitian(zeros(T, 2, 2)); check=false)
        @test isequal(@inferred(logabsdet(F0)), (big(-Inf), zero(TS)))
    end
end

@testset "Hermitian matrices of quaternions" begin
    # the real 4n×4n matrix of the map x -> A*x has each eigenvalue of a Hermitian A four times
    realmat(q) = [q.s -q.v1 -q.v2 -q.v3; q.v1 q.s -q.v3 q.v2; q.v2 q.v3 q.s -q.v1; q.v3 -q.v2 q.v1 q.s]
    realmat(A::AbstractMatrix) = hvcat(size(A, 2), (realmat(A[i,j]) for i in axes(A, 1) for j in axes(A, 2))...)
    for n in (1, 2, 5, 8), zerodiag in (false, true), uplo in (:U, :L), rook in (false, true)
        n == 1 && zerodiag && continue
        B = randn(Quaternion{Float64}, n, n)
        M = B + B'
        # a zero diagonal enforces 2-by-2 pivots
        zerodiag && foreach(i -> M[i,i] = 0, 1:n)
        F = bunchkaufman(Hermitian(M, uplo), rook)
        @test F isa BunchKaufman{Quaternion{Float64}}
        @test ishermitian(F) && !issymmetric(F)
        zerodiag && @test any(<(0), F.ipiv)
        @test Matrix(F) ≈ M
        b = randn(Quaternion{Float64}, n)
        X = randn(Quaternion{Float64}, 2, n)
        @test M * (F \ b) ≈ b
        @test (X / F) * M ≈ X
        @test inv(F) * M ≈ I
        λ = eigvals(Symmetric(realmat(M)))
        np, nn = count(>(0), λ) ÷ 4, count(<(0), λ) ÷ 4
        @test inertia(F) == (np, nn, 0)
        l, s = @inferred logabsdet(F)
        @test l ≈ sum(log ∘ abs, λ) / 4
        @test s === (-1.0)^nn
    end
    M = let B = randn(Quaternion{Float64}, 4, 4); B + B'; end
    # a plain Hermitian matrix is factorized as Hermitian
    @test ishermitian(bunchkaufman(M))
    # integer quaternions are promoted as by division
    Mi = let B = Quaternion.(rand(-3:3, 4, 4), rand(-3:3, 4, 4), 0, 1); B + B'; end
    @test bunchkaufman(Hermitian(Mi); check=false) isa BunchKaufman{Quaternion{Float64}}
    # the symmetric factorization requires commuting entries
    msg = "the Bunch-Kaufman decomposition of a symmetric, non-Hermitian matrix requires " *
        "real or complex numbers, but the element type is $(Quaternion{Float64})"
    @test_throws ArgumentError(msg) bunchkaufman(Symmetric(M))
    S = let B = randn(Quaternion{Float64}, 4, 4); B + transpose(B); end
    @test_throws ArgumentError(msg) bunchkaufman(S)
    # the real type of the generic code must be floating-point or rational
    msg = "the generic Bunch-Kaufman factorization requires a floating-point or rational " *
        "real type, but the element type is $Int"
    @test_throws ArgumentError(msg) bunchkaufman!([2 1; 1 2])
end

end # module TestBunchKaufman
