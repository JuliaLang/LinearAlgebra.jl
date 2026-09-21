# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestLQ

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: BlasComplex, BlasFloat, BlasReal, rmul!, lmul!
using LinearAlgebra: LQPackedQ, QRPackedQ, _lmul_lq!, _rmul_lq!, lqfactUnblocked!

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)
using Main.LinearAlgebraTestHelpers.Quaternions

m = 10

Random.seed!(1234321)

asquare = randn(ComplexF64, m, m) / 2
awide = randn(ComplexF64, m, m+3) / 2
bcomplex = randn(ComplexF64, m, 2) / 2

# helper functions to unambiguously recover explicit forms of an LQPackedQ
squareQ(Q::LinearAlgebra.LQPackedQ) = (n = size(Q.factors, 2); lmul!(Q, Matrix{eltype(Q)}(I, n, n)))
rectangularQ(Q::LinearAlgebra.LQPackedQ) = convert(Array, Q)

@testset for eltya in (Float32, Float64, ComplexF32, ComplexF64), n in (m, size(awide, 2))
    adata = m == n ? asquare : awide
    a = convert(Matrix{eltya}, eltya <: Complex ? adata : real(adata))
    ε = εa = eps(abs(float(one(eltya))))
    n1 = n ÷ 2

    α = rand(eltya)
    aα = fill(α,1,1)
    @test lq(α).L*lq(α).Q ≈ lq(aα).L*lq(aα).Q
    @test abs(lq(α).Q[1,1]) ≈ one(eltya)

    @testset for eltyb in (Float32, Float64, ComplexF32, ComplexF64, Int)
        b = eltyb == Int ? rand(1:5, m, 2) : convert(Matrix{eltyb}, eltyb <: Complex ? bcomplex : real(bcomplex))
        εb = eps(abs(float(one(eltyb))))
        ε = max(εa,εb)

        tab = promote_type(eltya,eltyb)

        @testset for isview in (false,true)
            let a = isview ? view(a, 1:m - 1, 1:n - 1) : a, b = isview ? view(b, 1:m - 1) : b, m = m - isview, n = n - isview
                lqa = lq(a)
                x = lqa\b
                l, q = lqa.L, lqa.Q
                qra = qr(a, ColumnNorm())
                @testset "Basic ops" begin
                    @test size(lqa,1) == size(a,1)
                    @test size(lqa,3) == 1
                    @test size(lqa.Q,3) == 1
                    @test Base.propertynames(lqa) == (:L, :Q)
                    ref_obs = (l, q)
                    for (ii, lq_obj) in enumerate(lqa)
                        @test ref_obs[ii] == lq_obj
                    end
                    @test_throws FieldError lqa.Z
                    @test Array(copy(adjoint(lqa))) ≈ a'
                    @test q*squareQ(q)' ≈ Matrix(I, n, n)
                    @test l*q ≈ a
                    @test Array(lqa) ≈ a
                    @test Array(copy(lqa)) ≈ a
                    @test LinearAlgebra.Factorization{eltya}(lqa) === lqa
                    @test Matrix{eltya}(q) isa Matrix{eltya}
                    # test Array{T}(LQPackedQ{T})
                    @test Array{eltya}(q) ≈ Matrix(q)
                end
                @testset "Binary ops" begin
                    k = size(a, 2)
                    T = Tridiagonal(rand(eltya, k-1), rand(eltya, k), rand(eltya, k-1))
                    @test lq(T) * T ≈ T * T rtol=3000ε
                    @test lqa * T ≈ a * T rtol=3000ε
                    @test a*x ≈ b rtol=3000ε
                    @test x ≈ qra \ b rtol=3000ε
                    @test lqa*x ≈ a*x rtol=3000ε
                    @test (sq = size(q.factors, 2); *(Matrix{eltyb}(I, sq, sq), adjoint(q))*squareQ(q)) ≈ Matrix(I, n, n) rtol=5000ε
                    if eltya != Int
                        @test Matrix{eltyb}(I, n, n)*q ≈ Matrix(I, n, n) * convert(LinearAlgebra.AbstractQ{tab}, q)
                    end
                    @test q*x ≈ squareQ(q)*x rtol=100ε
                    @test q'*x ≈ squareQ(q)'*x rtol=100ε
                    @test a*q ≈ a*squareQ(q) rtol=100ε
                    @test a*q' ≈ a*squareQ(q)' rtol=100ε
                    @test q*a'≈ squareQ(q)*a' rtol=100ε
                    @test q'*a' ≈ squareQ(q)'*a' rtol=100ε
                    @test_throws DimensionMismatch q*x[1:n1 + 1]
                    @test_throws DimensionMismatch adjoint(q) * Matrix{eltya}(undef,m+2,m+2)
                    @test_throws DimensionMismatch Matrix{eltyb}(undef,m+2,m+2)*q
                    if isa(a, DenseArray) && isa(b, DenseArray)
                        # use this to test 2nd branch in mult code
                        pad_a = vcat(I, a)
                        pad_x = hcat(I, x)
                        @test pad_a*q ≈ pad_a*squareQ(q) rtol=100ε
                        @test q'*pad_x ≈ squareQ(q)'*pad_x rtol=100ε
                    end
                end
            end
        end

        @testset "Matmul with LQ factorizations" begin
            lqa = lq(a[:,1:n1])
            l,q = lqa.L, lqa.Q
            @test rectangularQ(q)*rectangularQ(q)' ≈ Matrix(I, n1, n1)
            @test squareQ(q)'*squareQ(q) ≈ Matrix(I, n1, n1)
            @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1),q)
            @test lmul!(adjoint(q), rectangularQ(q)) ≈ Matrix(I, n1, n1)
            @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1), adjoint(q))
            @test_throws BoundsError size(q,-1)
        end
    end
end

@testset "getindex on LQPackedQ (#23733)" begin
    local m, n
    function getqs(F::LinearAlgebra.LQ)
        implicitQ = F.Q
        sq = size(implicitQ.factors, 2)
        explicitQ = lmul!(implicitQ, Matrix{eltype(implicitQ)}(I, sq, sq))
        return implicitQ, explicitQ
    end

    m, n = 3, 3 # reduced Q 3-by-3, full Q 3-by-3
    implicitQ, explicitQ = getqs(lq(randn(m, n)))
    @test implicitQ[1, 1] == explicitQ[1, 1]
    @test implicitQ[m, 1] == explicitQ[m, 1]
    @test implicitQ[1, n] == explicitQ[1, n]
    @test implicitQ[m, n] == explicitQ[m, n]

    m, n = 3, 4 # reduced Q 3-by-4, full Q 4-by-4
    implicitQ, explicitQ = getqs(lq(randn(m, n)))
    @test implicitQ[1, 1] == explicitQ[1, 1]
    @test implicitQ[m, 1] == explicitQ[m, 1]
    @test implicitQ[1, n] == explicitQ[1, n]
    @test implicitQ[m, n] == explicitQ[m, n]
    @test implicitQ[m+1, 1] == explicitQ[m+1, 1]
    @test implicitQ[m+1, n] == explicitQ[m+1, n]

    m, n = 4, 3 # reduced Q 3-by-3, full Q 3-by-3
    implicitQ, explicitQ = getqs(lq(randn(m, n)))
    @test implicitQ[1, 1] == explicitQ[1, 1]
    @test implicitQ[n, 1] == explicitQ[n, 1]
    @test implicitQ[1, n] == explicitQ[1, n]
    @test implicitQ[n, n] == explicitQ[n, n]
end

@testset "size on LQPackedQ (#23780)" begin
    # size(Q::LQPackedQ) yields the shape of Q's full/square form
    for ((mA, nA), nQ) in (
        ((3, 3), 3), # A 3-by-3 => full/square Q 3-by-3
        ((3, 4), 4), # A 3-by-4 => full/square Q 4-by-4
        ((4, 3), 3) )# A 4-by-3 => full/square Q 3-by-3
        @test size(lq(randn(mA, nA)).Q) == (nQ, nQ)
    end
end

@testset "postmultiplication with / right-application of LQPackedQ (#23779)" begin
    function getqs(F::LinearAlgebra.LQ)
        implicitQ = F.Q
        explicitQ = lmul!(implicitQ, Matrix{eltype(implicitQ)}(I, size(implicitQ)...))
        return implicitQ, explicitQ
    end
    # for any shape m-by-n of LQ-factored matrix, where Q is an LQPackedQ
    # A_mul_B*(C, Q) (Ac_mul_B*(C, Q)) operations should work for
    # *-by-n (n-by-*) C, which we test below via n-by-n C
    for (mA, nA) in ((3, 3), (3, 4), (4, 3))
        implicitQ, explicitQ = getqs(lq(randn(mA, nA)))
        C = randn(nA, nA)
        @test *(C, implicitQ) ≈ *(C, explicitQ)
        @test *(C, adjoint(implicitQ)) ≈ *(C, adjoint(explicitQ))
        @test *(adjoint(C), implicitQ) ≈ *(adjoint(C), explicitQ)
        @test *(adjoint(C), adjoint(implicitQ)) ≈ *(adjoint(C), adjoint(explicitQ))
    end
    # where the LQ-factored matrix has at least as many rows m as columns n,
    # Q's full/square and reduced/rectangular forms have the same shape (n-by-n). hence we expect
    # _only_ *-by-n (n-by-*) C to work in A_mul_B*(C, Q) (Ac_mul_B*(C, Q)) ops.
    # and hence the n-by-n C tests above suffice.
    #
    # where the LQ-factored matrix has more columns n than rows m,
    # Q's full/square form is n-by-n whereas its reduced/rectangular form is m-by-n.
    # hence we need also test *-by-m C with
    # A*_mul_B(C, Q) ops, as below via m-by-m C.
    mA, nA = 3, 4
    implicitQ, explicitQ = getqs(lq(randn(mA, nA)))
    C = randn(mA, mA)
    zeroextCright = hcat(C, zeros(eltype(C), mA))
    zeroextCdown = vcat(C, zeros(eltype(C), (1, mA)))
    @test *(C, implicitQ) ≈ *(zeroextCright, explicitQ)
    @test *(adjoint(C), implicitQ) ≈ *(adjoint(zeroextCdown), explicitQ)
    @test_throws DimensionMismatch C * adjoint(implicitQ)
    @test_throws DimensionMismatch adjoint(C) * adjoint(implicitQ)
end

@testset "det(Q::LQPackedQ)" begin
    @testset for n in 1:3, m in 1:3
        @testset "real" begin
            _, Q = lq(randn(n, m))
            @test det(Q) ≈ det(Q*I)
            @test abs(det(Q)) ≈ 1
        end
        @testset "complex" begin
            _, Q = lq(randn(ComplexF64, n, m))
            @test det(Q) ≈ det(Q*I)
            @test abs(det(Q)) ≈ 1
        end
    end
end

@testset "REPL printing" begin
    bf = IOBuffer()
    show(bf, "text/plain", lq(Matrix(I, 4, 4)))
    seekstart(bf)
    @test String(take!(bf)) == """
$(LinearAlgebra.LQ){Float64, Matrix{Float64}, Vector{Float64}}
L factor:
4×4 Matrix{Float64}:
 1.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
Q factor: 4×4 $(LinearAlgebra.LQPackedQ){Float64, Matrix{Float64}, Vector{Float64}}"""
end

@testset "adjoint of LQ" begin
    n = 5

    for b in (ones(n), ones(n, 2), ones(Complex{Float64}, n, 2))
        for A in (
            randn(n, n),
            # Tall problems become least squares problems similarly to QR
            randn(n - 2, n),
            complex.(randn(n, n), randn(n, n)))

            F = lq(A)
            @test A'\b ≈ F'\b
        end
        @test_throws DimensionMismatch lq(randn(n, n + 2))'\b
    end

end

@testset "LQ factorization of Q" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        L1, Q1 = lq(randn(T, 5, 5))
        L2, Q2 = lq(Q1)
        @test Matrix(Q1) ≈ Matrix(Q2)
        @test L2 ≈ I
    end
end

# The reflectors of a `Float64` LQ factorization, promoted to `BigFloat` and with τ
# recomputed so that every `I - τᵢvᵢvᵢ'` is orthogonal to full `BigFloat` precision.
function bigfloat_lq_reflectors(m, n)
    factors = big.(lq(randn(m, n)).factors)
    nQ, k = size(factors, 2), min(size(factors)...)
    τ = map(1:k) do i
        v = zeros(BigFloat, nQ)
        v[i] = 1
        for l in i+1:nQ
            v[l] = conj(factors[i,l])
        end
        return 2 / (v'v)
    end
    return LQPackedQ(factors, τ)
end

# `Q = H_k' ⋯ H_1'` with `Hᵢ = I - vᵢτᵢvᵢ'`, built from explicit matrix products so that
# every operand order is fixed. Used to check the kernels for a non-commutative element type.
function quaternion_lq_refQ(factors, τ, n, k)
    QT = eltype(factors)
    Q = Matrix{QT}(I, n, n)
    for i in 1:k
        v = zeros(QT, n)
        v[i] = one(Float64)
        for l in i+1:n
            v[l] = conj(factors[i,l])
        end
        # Hᵢ' = I - vᵢ conj(τᵢ) vᵢ', and Q accumulates so that H_k' ends up leftmost
        Q = (Matrix{QT}(I, n, n) - (v .* conj(τ[i])) * v') * Q
    end
    return Q
end

@testset "generic lmul!/rmul! with LQPackedQ and its adjoint" begin
    @testset "matches LAPACK: $elty, ($m,$n)" for
            elty in (Float32, Float64, ComplexF32, ComplexF64),
            (m, n) in ((4, 6), (6, 4), (5, 5), (1, 1), (1, 7), (8, 1), (13, 6))
        A = elty <: Complex ? complex.(randn(m, n), randn(m, n)) : randn(m, n)
        Q = lq(convert(Matrix{elty}, A)).Q
        nQ = size(Q, 1)
        for p in (1, 3)
            B = elty <: Complex ? complex.(randn(nQ, p), randn(nQ, p)) : randn(nQ, p)
            B = convert(Matrix{elty}, B)
            @test _lmul_lq!(Q, copy(B), Val(false)) ≈ lmul!(Q, copy(B)) ≈ Q * B
            @test _lmul_lq!(Q, copy(B), Val(true)) ≈ lmul!(Q', copy(B)) ≈ Q' * B
            C = elty <: Complex ? complex.(randn(p, nQ), randn(p, nQ)) : randn(p, nQ)
            C = convert(Matrix{elty}, C)
            @test _rmul_lq!(copy(C), Q, Val(false)) ≈ rmul!(copy(C), Q) ≈ C * Q
            @test _rmul_lq!(copy(C), Q, Val(true)) ≈ rmul!(copy(C), Q') ≈ C * Q'
        end
        b = elty <: Complex ? complex.(randn(nQ), randn(nQ)) : randn(nQ)
        b = convert(Vector{elty}, b)
        @test _lmul_lq!(Q, copy(b), Val(false)) ≈ lmul!(Q, copy(b)) ≈ Q * b
        @test _lmul_lq!(Q, copy(b), Val(true)) ≈ lmul!(Q', copy(b)) ≈ Q' * b
    end

    @testset "dispatch for non-BLAS eltypes" begin
        for (m, n) in ((4, 6), (6, 4), (5, 5))
            Q = lq(randn(m, n)).Q
            Qsq = squareQ(Q)
            nQ = size(Q, 1)
            for B in (big.(randn(nQ, 3)), big.(randn(nQ)))
                @test lmul!(Q, copy(B)) ≈ Qsq * B rtol=1e-12
                @test lmul!(Q', copy(B)) ≈ Qsq' * B rtol=1e-12
            end
            C = big.(randn(3, nQ))
            @test rmul!(copy(C), Q) ≈ C * Qsq rtol=1e-12
            @test rmul!(copy(C), Q') ≈ C * Qsq' rtol=1e-12
        end
    end

    @testset "full BigFloat precision: ($m,$n)" for (m, n) in ((4, 6), (6, 4), (5, 5), (9, 3))
        Q = bigfloat_lq_reflectors(m, n)
        nQ = size(Q, 1)
        Id = Matrix{BigFloat}(I, nQ, nQ)
        @test lmul!(Q', lmul!(Q, copy(Id))) ≈ Id
        @test lmul!(Q, lmul!(Q', copy(Id))) ≈ Id
        @test rmul!(rmul!(copy(Id), Q), Q') ≈ Id
        @test rmul!(rmul!(copy(Id), Q'), Q) ≈ Id
        # left- and right-multiplication have to be consistent: A*Q == (Q'*A')'
        A = big.(randn(4, nQ))
        @test rmul!(copy(A), Q) ≈ collect(lmul!(Q', collect(A'))')
        @test rmul!(copy(A), Q') ≈ collect(lmul!(Q, collect(A'))')
    end

    # `qsize_check` lets `Q'*B` and `A*Q` also take the operand with `size(Q.factors, 1)`
    # rows resp. columns, which `mul!` zero-extends to the full `nQ`. That is only
    # reachable through `*`, never through `lmul!`/`rmul!`.
    @testset "flexible operand size via *: $elty, ($m,$n)" for
            elty in (Float64, ComplexF64), (m, n) in ((4, 6), (2, 7), (1, 5), (5, 6))
        A = elty <: Complex ? complex.(randn(m, n), randn(m, n)) : randn(m, n)
        Q = lq(convert(Matrix{elty}, A)).Q   # Q is n×n, factors are m×n with m < n
        Qsq = squareQ(Q)
        p = 3
        B = elty <: Complex ? complex.(randn(m, p), randn(m, p)) : randn(m, p)
        B = convert(Matrix{elty}, B)
        Bext = [B; zeros(elty, n - m, p)]
        @test Q' * B ≈ Qsq' * Bext
        @test size(Q' * B) == (n, p)
        C = elty <: Complex ? complex.(randn(p, m), randn(p, m)) : randn(p, m)
        C = convert(Matrix{elty}, C)
        Cext = [C zeros(elty, p, n - m)]
        @test C * Q ≈ Cext * Qsq
        @test size(C * Q) == (p, n)
        b = elty <: Complex ? complex.(randn(m), randn(m)) : randn(m)
        b = convert(Vector{elty}, b)
        @test Q' * b ≈ Qsq' * [b; zeros(elty, n - m)]
        # the other two directions admit only the full size
        @test_throws DimensionMismatch Q * B
        @test_throws DimensionMismatch C * Q'
    end

    # The scalar τᵢ sits between the vectors in `Hᵢ = I - vᵢτᵢvᵢ'`, so `lmul!` must apply it
    # as `vᵢ*(τᵢ*(vᵢ'B))` while `rmul!` must apply it as `((A*vᵢ)*τᵢ)*vᵢ'`. Commutative
    # element types cannot distinguish the two, hence the quaternion check.
    @testset "non-commutative element type: ($m,$n)" for (m, n) in ((4, 6), (6, 4), (5, 5), (3, 7))
        QT = Quaternion{Float64}
        k = min(m, n)
        factors = [randn(QT) for _ in CartesianIndices((m, n))]
        τ = [randn(QT) for _ in 1:k]
        Q = LQPackedQ(factors, τ)
        Qref = quaternion_lq_refQ(factors, τ, n, k)
        # Anchor the hand-built reference against the pre-existing `QRPackedQ` kernel: the
        # adjoint undoes the conjugated row storage, so `LQPackedQ(factors, τ)` acts as
        # `QRPackedQ(factors', τ)'`. That kernel predates this branch and is independent of
        # everything under test here, so a consistent mistake in the helper cannot make a
        # wrong kernel look right.
        let P = QRPackedQ(factors', τ), B = [randn(QT) for _ in CartesianIndices((n, 2))]
            @test lmul!(P', copy(B)) ≈ Qref * B
            @test lmul!(P, copy(B)) ≈ Qref' * B
            A = [randn(QT) for _ in CartesianIndices((2, n))]
            @test rmul!(copy(A), P') ≈ A * Qref
            @test rmul!(copy(A), P) ≈ A * Qref'
        end
        for p in (1, 3)
            B = [randn(QT) for _ in CartesianIndices((n, p))]
            @test _lmul_lq!(Q, copy(B), Val(false)) ≈ Qref * B
            @test _lmul_lq!(Q, copy(B), Val(true)) ≈ Qref' * B
            A = [randn(QT) for _ in CartesianIndices((p, n))]
            @test _rmul_lq!(copy(A), Q, Val(false)) ≈ A * Qref
            @test _rmul_lq!(copy(A), Q, Val(true)) ≈ A * Qref'
        end
    end

    @testset "dimension mismatch" begin
        Q = lq(randn(4, 6)).Q   # Q is 6×6
        @test_throws DimensionMismatch lmul!(Q, big.(randn(5, 2)))
        @test_throws DimensionMismatch lmul!(Q', big.(randn(5, 2)))
        @test_throws DimensionMismatch lmul!(Q, big.(randn(5)))
        @test_throws DimensionMismatch lmul!(Q', big.(randn(5)))
        @test_throws DimensionMismatch rmul!(big.(randn(2, 5)), Q)
        @test_throws DimensionMismatch rmul!(big.(randn(2, 5)), Q')
    end
end

@testset "generic (unblocked) LQ factorization" begin
    @testset "reproduces LAPACK gelqf!: $elty, ($m,$n)" for
            elty in (Float32, Float64, ComplexF32, ComplexF64),
            (m, n) in ((4, 6), (6, 4), (5, 5), (1, 1), (1, 7), (8, 1), (13, 6), (2, 9))
        A = elty <: Complex ? complex.(randn(m, n), randn(m, n)) : randn(m, n)
        A = convert(Matrix{elty}, A)
        Flap = lq!(copy(A))                  # LAPACK path
        Fgen = lqfactUnblocked!(copy(A))     # generic path
        @test Fgen.factors ≈ Flap.factors
        @test Fgen.τ ≈ Flap.τ
        @test Fgen.L ≈ Flap.L
        @test Matrix(Fgen.Q) ≈ Matrix(Flap.Q)
        @test Fgen.L * Fgen.Q ≈ A
    end

    @testset "lq for non-BLAS eltypes: ($m,$n)" for (m, n) in ((4, 6), (6, 4), (5, 5), (1, 7), (8, 1), (7, 3))
        for A in (big.(randn(m, n)), complex.(big.(randn(m, n)), big.(randn(m, n))))
            F = lq(A)
            @test F isa LQ{eltype(A)}
            @test F.L * F.Q ≈ A
            @test istril(F.L)
            nQ = size(F.Q, 1)
            Id = Matrix{eltype(A)}(I, nQ, nQ)
            # the reflectors are orthogonal to full BigFloat precision
            @test lmul!(F.Q', lmul!(F.Q, copy(Id))) ≈ Id
            @test rmul!(rmul!(copy(Id), F.Q), F.Q') ≈ Id
        end
    end

    @testset "lq of exact and unusual element types" begin
        @test lq(3).L * lq(3).Q ≈ fill(3.0, 1, 1)
        Ai = [1 2 3; 4 5 6]
        @test lq(Ai).L * lq(Ai).Q ≈ Ai
        Ar = Rational{BigInt}[1//2 1//3; 1//5 1//7]
        @test lq(Ar).L * lq(Ar).Q ≈ float.(Ar)
    end

    @testset "minimum-norm solve for an underdetermined BigFloat system" begin
        A = big.(randn(3, 6))
        b = big.(randn(3))
        x = lq(A) \ b
        @test A * x ≈ b
    end
end

end # module TestLQ
