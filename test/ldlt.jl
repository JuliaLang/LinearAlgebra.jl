# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestLDLT

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random

Random.seed!(123)

@testset "Factorization conversions of LDLT" begin
    S = SymTridiagonal(randn(5), randn(4))
    F = ldlt(S)
    @test Factorization{eltype(S)}(F) === F
    @test Array(Factorization{complex(eltype(S))}(F)) ≈ Array(ldlt(complex(S)))
    @test eltype(Factorization{complex(eltype(S))}) == complex(eltype(S))
end

@testset "eltype conversions of LDLT" begin
    S = SymTridiagonal(randn(5), randn(4))
    F = ldlt(S)
    Fc = LDLt{ComplexF32}(F.data)
    @test Fc isa LDLt{ComplexF32}
    @test Array(Fc) ≈ ComplexF32.(Array(S))
end

@testset "Accessing fields of LDLT" begin
    S = SymTridiagonal(randn(5), randn(4))
    F = ldlt(S)
    @test getproperty(F, :L) == transpose(getproperty(F, :Lt))
    @test getproperty(F, :d) == diag(getproperty(F, :D), 0)
end

@testset "REPL printing of LDLT" begin
    S = SymTridiagonal(randn(5), randn(4))
    F = ldlt(S)
    ldltstring = sprint((t, s) -> show(t, "text/plain", s), F)
    lstring = sprint((t, s) -> show(t, "text/plain", s), F.L)
    dstring = sprint((t, s) -> show(t, "text/plain", s), F.D)
    @test ldltstring == "$(summary(F))\nL factor:\n$lstring\nD factor:\n$dstring"
end

@testset "Pivoted LDLt (BunchPivot)" begin
    @testset "unsupported pivoting strategy" begin
        S = SymTridiagonal(randn(5), randn(4))
        @test_throws ArgumentError ldlt(S, RowMaximum())
        @test_throws ArgumentError ldlt!(copy(S), ColumnNorm())
    end

    @testset "factorize and \\ use BunchPivot" begin
        S = SymTridiagonal(zeros(4), [1.0, 0.0, 1.0])
        @test factorize(S) isa LDLtPivoted
        @test factorize(S).blocks == [1, 3]
        b = [1.0, 2.0, 3.0, 4.0]
        @test S \ b == [2.0, 1.0, 4.0, 3.0]
        @test factorize(S) \ b == S \ b
        @test_throws SingularException SymTridiagonal([1.0, 1.0], [1.0]) \ [1.0, 1.0]
        # inv uses factorize, which failed with a zero pivot before (#601)
        A = SymTridiagonal([1.0, 0.0, 0.0], [0.0, 1.0])
        @test inv(A) == inv(Matrix(A))
        @test A \ [1.0, 2.0, 3.0] == [1.0, 3.0, 2.0]
        # ldiv! and rdiv!, also with a shift that produces a zero diagonal
        x = [2.0, 1.0, 4.0, 3.0]
        @test ldiv!(S, copy(b)) == x
        @test rdiv!(copy(b'), S) == x'
        S1 = SymTridiagonal(ones(4), [1.0, 0.0, 1.0])
        @test_throws ZeroPivotException ldlt(S1; shift=-1)
        @test ldiv!(S1, copy(b); shift=-1) == x
        @test rdiv!(copy(b'), S1; shift=-1) == x'
        B = [b 2b]
        @test ldiv!(S1, copy(B); shift=-1.0) == [x 2x]
        @test rdiv!(copy(B'), S1; shift=-1.0) == [x 2x]'
        @test_throws SingularException ldiv!(S1, copy(b); shift=-2)
        # logabsdet, also for singular matrices and zero pivots
        @test logabsdet(S) == (0.0, 1.0)
        @test logabsdet(S1; shift=-1) == (0.0, 1.0)
        for Ssing in (SymTridiagonal([1.0, 1.0], [1.0]), SymTridiagonal([0.0, 1.0, 1.0], [0.0, 1.0]),
                      SymTridiagonal(zeros(5), ones(4)))
            @test logabsdet(Ssing) == (-Inf, 0.0)
            @test logabsdet(Ssing) == logabsdet(Matrix(Ssing))
        end
        @test logabsdet(S1; shift=-2) == (-Inf, 0.0)
        for elty in (Float64, ComplexF64)
            Sr = SymTridiagonal(randn(elty, 20) .* rand((0, 1), 20), randn(elty, 19))
            @test logabsdet(Sr)[1] ≈ logabsdet(Matrix(Sr))[1]
            @test logabsdet(Sr)[2] ≈ logabsdet(Matrix(Sr))[2]
            @test logabsdet(Sr; shift=0.5)[1] ≈ logabsdet(Matrix(Sr) + 0.5I)[1]
        end
        # shifted solves with the Hessenberg factorization of an indefinite Hermitian matrix
        A = Matrix(Hermitian(randn(ComplexF64, 6, 6)))
        F = hessenberg(Hermitian(A))
        for μ in (0.3, -0.5 + 0.2im, F.H.dv[1] + 0)
            bc = randn(ComplexF64, 6)
            @test (F + μ*I) \ bc ≈ (A + μ*I) \ bc
            @test bc' / (F + μ*I) ≈ bc' / (A + μ*I)
        end
    end

    @testset "NoPivot is the default and returns an LDLt" begin
        S = SymTridiagonal(randn(5), randn(4))
        @test ldlt(S, NoPivot()) == ldlt(S)
        @test ldlt(S, NoPivot()) isa LDLt
        @test ldlt(S, BunchPivot()) isa LDLtPivoted
    end

    @testset "indefinite matrices with zero diagonal, eltype $elty" for elty in (Float64, ComplexF64, BigFloat, Rational{BigInt})
        # zero pivot without pivoting, 2×2 pivot at the first position
        S = SymTridiagonal(elty[0, 1, 2], elty[1, 1])
        @test_throws ZeroPivotException ldlt(S)
        F = ldlt(S, BunchPivot())
        @test F isa LDLtPivoted{elty, SymTridiagonal{elty, Vector{elty}}, Vector{Int}}
        @test F.blocks == [1]
        @test size(F) == (3, 3)
        @test size(F, 1) == size(F, 2) == 3
        @test F.d == F.D.dv
        @test F.D == SymTridiagonal(elty[0, 1, 2], elty[1, 0])
        @test F.L == [1 0 0; 0 1 0; 1 0 1]
        @test F.Lt == transpose(F.L)
        @test F.L isa UnitLowerTriangular
        @test F.Lt isa UnitUpperTriangular
        @test F.L * F.D * F.Lt == S
        @test SymTridiagonal(F) == S
        @test Array(F) == Array(S)
        @test Matrix(F) == Array(S)
        @test AbstractMatrix(F) == S
        b = elty[1, 2, 3]
        x = F \ b
        @test x == elty[0, 1, 1]
        @test S * x == b
        @test det(F) ≈ -2
        @test logabsdet(F)[1] ≈ log(2)
        @test logabsdet(F)[2] == -1
        @test inv(F) ≈ inv(Array(S))
        # all-zero diagonal, 2×2 pivot at the last position
        S = SymTridiagonal(zeros(elty, 4), ones(elty, 3))
        @test_throws ZeroPivotException ldlt(S)
        F = ldlt(S, BunchPivot())
        @test F.blocks == [1, 3]
        @test F.L * F.D * F.Lt == S
        @test F.D == SymTridiagonal(zeros(elty, 4), elty[1, 0, 1])
        @test F.L == [1 0 0 0; 0 1 0 0; 1 0 1 0; 0 0 0 1]
        B = elty[1 2; 3 4; 5 6; 7 8]
        @test F \ B == Array(S) \ B
        @test B' / F == B' / Array(S)
        @test det(F) ≈ 1
    end

    @testset "random matrices, eltype $elty, size $n" for elty in (Float32, Float64, ComplexF64), n in (1, 2, 3, 4, 7, 20)
        for pattern in (ones(n), rand((0, 1), n), zeros(n))
            # zero out parts of the diagonal to produce 2×2 pivots
            dv = randn(elty, n) .* pattern
            ev = randn(elty, n - 1)
            S = SymTridiagonal(dv, ev)
            A = Array(S)
            # skip singular and ill-conditioned matrices
            κ = cond(A)
            (isfinite(κ) && κ < 1e4) || continue
            F = ldlt(S, BunchPivot())
            @test issorted(F.blocks) && all(k -> 1 <= k < n, F.blocks)
            @test all(diff(F.blocks) .>= 2)
            @test F.L * F.D * F.Lt ≈ A
            @test Array(F) ≈ A
            @test SymTridiagonal(F) ≈ S
            b = randn(elty, n)
            B = randn(elty, n, 3)
            @test F \ b ≈ A \ b
            @test F \ B ≈ A \ B
            @test B' / F ≈ B' / A
            @test transpose(B) / F ≈ transpose(B) / A
            @test F' \ b ≈ A' \ b
            @test transpose(F) \ b ≈ transpose(A) \ b
            @test inv(F) ≈ inv(A)
            @test det(F) ≈ det(A)
            @test logabsdet(F)[1] ≈ logabsdet(A)[1] atol=(1e-5 * (1 + abs(logabsdet(A)[1])))
            @test logabsdet(F)[2] ≈ logabsdet(A)[2]
            # in-place factorization
            Sc = copy(S)
            Fc = ldlt!(Sc, BunchPivot())
            @test Fc.data === Sc
            @test Fc == F
            # shift
            s = 2one(elty)
            Fs = ldlt(S, BunchPivot(); shift=s)
            @test Array(Fs) ≈ A + s*I
            @test Fs \ b ≈ (A + s*I) \ b
            # dimension mismatch
            @test_throws DimensionMismatch F \ randn(elty, n + 1)
            @test_throws DimensionMismatch ldiv!(F, randn(elty, n + 1, 2))
        end
    end

    @testset "comparison with bunchkaufman" begin
        # the standard example: two 2×2 blocks, which do not require any interchanges,
        # so the factors coincide with those of the (dense) Bunch-Kaufman factorization
        for ev in ([1., 0., 1.], [1., 1., 1.], [1., 2., 3.])
            S = SymTridiagonal(zeros(4), ev)
            @test_throws ZeroPivotException ldlt(S)
            F = ldlt(S, BunchPivot())
            B = bunchkaufman(Symmetric(Matrix(S), :L))
            @test F.blocks == [1, 3]
            @test B.p == 1:4
            @test F.D ≈ B.D
            @test F.L ≈ B.L
            @test F.L * F.D * F.Lt ≈ S
            b = [1., 2., 3., 4.]
            @test F \ b ≈ B \ b
            @test det(F) ≈ det(B)
            @test inv(F) ≈ inv(B)
        end
        # random indefinite matrices: the pivots may differ, but the solutions agree
        for n in (5, 20), elty in (Float64, ComplexF64)
            S = SymTridiagonal(randn(elty, n) .* rand((0, 1), n), randn(elty, n - 1))
            A = Matrix(S)
            κ = cond(A)
            (isfinite(κ) && κ < 1e4) || continue
            F = ldlt(S, BunchPivot())
            B = bunchkaufman(Symmetric(A, :L))
            b = randn(elty, n, 2)
            @test F \ b ≈ B \ b
            @test det(F) ≈ det(B)
            @test logabsdet(F)[1] ≈ logabsdet(B)[1]
            @test logabsdet(F)[2] ≈ logabsdet(B)[2]
        end
    end

    @testset "numerical stability" begin
        # a tiny leading pivot would blow up the unpivoted factorization
        ε = 1e-14
        S = SymTridiagonal([ε, 1.0, 1.0], [1.0, 1.0])
        A = Array(S)
        b = [1.0, 2.0, 3.0]
        x = A \ b
        F = ldlt(S, BunchPivot())
        @test F.blocks == [1]
        @test F \ b ≈ x rtol=1e-12
        @test maximum(abs, F.d) <= 3 * maximum(abs, A)
        # a 1×1 pivot is chosen if the diagonal entry is large enough
        S = SymTridiagonal([1.0, -2.0, 3.0], [1.0, 1.0])
        F = ldlt(S, BunchPivot())
        @test isempty(F.blocks)
        @test F.L == ldlt(S).L
        @test F.d == ldlt(S).d
    end

    @testset "singular matrices" begin
        S = SymTridiagonal([0.0, 1.0, 1.0], [0.0, 1.0])
        @test_throws SingularException(1) ldlt(S, BunchPivot())
        S = SymTridiagonal([1.0, 0.0, 1.0], [0.0, 0.0])
        @test_throws SingularException(2) ldlt(S, BunchPivot())
        # zero final pivot
        S = SymTridiagonal([1.0, 1.0], [1.0])
        @test_throws SingularException(2) ldlt(S, BunchPivot())
        S = SymTridiagonal([0.0], Float64[])
        @test_throws SingularException(1) ldlt(S, BunchPivot())
        S = SymTridiagonal([1//1, 2, 1], [1//1, 1])
        @test_throws SingularException(3) ldlt(S, BunchPivot())
    end

    @testset "empty matrix" begin
        S = SymTridiagonal(Float64[], Float64[])
        F = ldlt(S, BunchPivot())
        @test size(F) == (0, 0)
        @test isempty(F.blocks)
        @test F \ Float64[] == Float64[]
        @test Array(F) == zeros(0, 0)
    end

    @testset "conversions" begin
        S = SymTridiagonal([0., 1., 2., 0., 1.], [1., 1., 1., 1.])
        F = ldlt(S, BunchPivot())
        @test Factorization{Float64}(F) === F
        @test LDLtPivoted{Float64}(F) === F
        @test typeof(F)(F) === F
        Fc = Factorization{ComplexF64}(F)
        @test Fc isa LDLtPivoted{ComplexF64, SymTridiagonal{ComplexF64, Vector{ComplexF64}}, Vector{Int}}
        @test Array(Fc) ≈ Array(S)
        @test Fc.blocks == F.blocks
        @test eltype(Fc) == ComplexF64
        Fc32 = LDLtPivoted{ComplexF32}(F.data, F.blocks)
        @test Fc32 isa LDLtPivoted{ComplexF32}
        @test Array(Fc32) ≈ ComplexF32.(Array(S))
        @test convert(LDLtPivoted{Float32, SymTridiagonal{Float32, Vector{Float32}}, Vector{Int32}}, F) isa
            LDLtPivoted{Float32, SymTridiagonal{Float32, Vector{Float32}}, Vector{Int32}}
        @test Array(convert(LDLtPivoted{Float32}, F)) ≈ Float32.(Array(S))
    end

    @testset "adjoint" begin
        Sr = SymTridiagonal([0., 1., 2.], [1., 1.])
        F = ldlt(Sr, BunchPivot())
        @test F' === F
        @test transpose(F) === F
        Sc = SymTridiagonal(ComplexF64[0, 1, 2], ComplexF64[1 + im, 1 - im])
        F = ldlt(Sc, BunchPivot())
        @test F' isa LDLtPivoted{ComplexF64}
        @test Array(F') ≈ Array(Sc)'
        @test transpose(F) === F
        b = randn(ComplexF64, 3)
        @test F' \ b ≈ Array(Sc)' \ b
        @test transpose(F) \ b ≈ transpose(Array(Sc)) \ b
    end

    @testset "REPL printing" begin
        S = SymTridiagonal([0., 1., 2.], [1., 1.])
        F = ldlt(S, BunchPivot())
        ldltstring = sprint((t, s) -> show(t, "text/plain", s), F)
        lstring = sprint((t, s) -> show(t, "text/plain", s), F.L)
        dstring = sprint((t, s) -> show(t, "text/plain", s), F.D)
        @test ldltstring == "$(summary(F))\nL factor:\n$lstring\nD factor:\n$dstring"
    end
end

end # module TestLDLT
