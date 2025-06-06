# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestQR

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: BlasComplex, BlasFloat, BlasReal, QRPivoted, rmul!, lmul!

n = 10

# Split n into 2 parts for tests needing two matrices
n1 = div(n, 2)
n2 = 2*n1

Random.seed!(1234325)

areal = randn(n,n)/2
aimg  = randn(n,n)/2
a2real = randn(n,n)/2
a2img  = randn(n,n)/2
breal = randn(n,2)/2
bimg  = randn(n,2)/2

# helper functions to unambiguously recover explicit forms of an implicit QR Q
squareQ(Q::LinearAlgebra.AbstractQ) = Q*I
rectangularQ(Q::LinearAlgebra.AbstractQ) = Matrix(Q)

@testset for eltya in (Float32, Float64, ComplexF32, ComplexF64, BigFloat, Int)
    raw_a = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex.(areal, aimg) : areal)
    raw_a2 = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex.(a2real, a2img) : a2real)
    asym = raw_a' + raw_a                  # symmetric indefinite
    apd  = raw_a' * raw_a                 # symmetric positive-definite
    ε = εa = eps(abs(float(one(eltya))))

    @testset for eltyb in (Float32, Float64, ComplexF32, ComplexF64, Int)
        raw_b = eltyb == Int ? rand(1:5, n, 2) : convert(Matrix{eltyb}, eltyb <: Complex ? complex.(breal, bimg) : breal)
        εb = eps(abs(float(one(eltyb))))
        ε = max(εa, εb)
        tab = promote_type(eltya, eltyb)

        @testset "QR decomposition of a Number" begin
            α = rand(eltyb)
            aα = fill(α, 1, 1)
            @test qr(α).Q * qr(α).R ≈ qr(aα).Q * qr(aα).R
            @test abs(qr(α).Q[1,1]) ≈ one(eltyb)
        end

        for (a, b) in ((raw_a, raw_b),
               (view(raw_a, 1:n-1, 1:n-1), view(raw_b, 1:n-1, 1)))
            a_1 = size(a, 1)
            @testset "QR decomposition (without pivoting)" begin
                qra   = @inferred qr(a)
                q, r  = qra.Q, qra.R
                @test_throws FieldError qra.Z
                @test q'*squareQ(q) ≈ Matrix(I, a_1, a_1)
                @test q*squareQ(q)' ≈ Matrix(I, a_1, a_1)
                @test q'*Matrix(1.0I, a_1, a_1)' ≈ squareQ(q)'
                @test squareQ(q)'q ≈ Matrix(I, a_1, a_1)
                @test Matrix(1.0I, a_1, a_1)'q' ≈ squareQ(q)'
                @test q*r ≈ a
                @test a*(qra\b) ≈ b atol=3000ε
                @test Array(qra) ≈ a
                sq = size(q.factors, 2)
                @test *(Matrix{eltyb}(I, sq, sq), adjoint(q)) * squareQ(q) ≈ Matrix(I, sq, sq) atol=5000ε
                if eltya != Int
                    @test Matrix{eltyb}(I, a_1, a_1)*q ≈ squareQ(convert(LinearAlgebra.AbstractQ{tab}, q))
                    ac = copy(a)
                    @test qr!(a[:, 1:5])\b == qr!(view(ac, :, 1:5))\b
                end
                qrstring = sprint((t, s) -> show(t, "text/plain", s), qra)
                rstring  = sprint((t, s) -> show(t, "text/plain", s), r)
                qstring  = sprint((t, s) -> show(t, "text/plain", s), q)
                @test qrstring == "$(summary(qra))\nQ factor: $qstring\nR factor:\n$rstring"
                # iterate
                q, r = qra
                @test q*r ≈ a
                # property names
                @test Base.propertynames(qra)       == (:R, :Q)
            end
            @testset "Thin QR decomposition (without pivoting)" begin
                qra   = @inferred qr(a[:, 1:n1], NoPivot())
                q,r   = qra.Q, qra.R
                @test_throws FieldError qra.Z
                @test q'*squareQ(q) ≈ Matrix(I, a_1, a_1)
                @test q'*rectangularQ(q) ≈ Matrix(I, a_1, n1)
                @test q*r ≈ a[:, 1:n1]
                @test q*b[1:n1] ≈ rectangularQ(q)*b[1:n1] atol=100ε
                @test q*b ≈ squareQ(q)*b atol=100ε
                if eltya != Int
                    @test Array{eltya}(q) ≈ rectangularQ(q)
                end
                @test_throws DimensionMismatch q*b[1:n1 + 1]
                @test_throws DimensionMismatch b[1:n1 + 1]*q'
                sq = size(q.factors, 2)
                @test *(UpperTriangular(Matrix{eltyb}(I, sq, sq)), adjoint(q))*squareQ(q) ≈ Matrix(I, n1, a_1) atol=5000ε
                if eltya != Int
                    @test Matrix{eltyb}(I, a_1, a_1)*q ≈ squareQ(convert(LinearAlgebra.AbstractQ{tab},q))
                end
                # iterate
                q, r = qra
                @test q*r ≈ a[:, 1:n1]
                # property names
                @test Base.propertynames(qra)       == (:R, :Q)
            end
            @testset "(Automatic) Fat (pivoted) QR decomposition" begin
                @inferred qr(a, ColumnNorm())

                qrpa  = factorize(a[1:n1,:])
                q,r = qrpa.Q, qrpa.R
                @test_throws FieldError qrpa.Z
                p = qrpa.p
                @test q'*squareQ(q) ≈ Matrix(I, n1, n1)
                @test q*squareQ(q)' ≈ Matrix(I, n1, n1)
                sq = size(q, 2);
                @test (UpperTriangular(Matrix{eltya}(I, sq, sq))*q')*squareQ(q) ≈ Matrix(I, n1, n1)
                @test q*r ≈ (isa(qrpa,QRPivoted) ? a[1:n1,p] : a[1:n1,:])
                @test q*r[:,invperm(p)] ≈ a[1:n1,:]
                @test q*r*transpose(qrpa.P) ≈ a[1:n1,:]
                @test a[1:n1,:]*(qrpa\b[1:n1]) ≈ b[1:n1] atol=5000ε
                @test Array(qrpa) ≈ a[1:5,:]
                if eltya != Int
                    @test Array{eltya}(q) ≈ Matrix(q)
                end
                @test_throws DimensionMismatch q*b[1:n1+1]
                @test_throws DimensionMismatch b[1:n1+1]*q'
                if eltya != Int
                    @test Matrix{eltyb}(I, n1, n1)*q ≈ squareQ(convert(LinearAlgebra.AbstractQ{tab},q))
                end
                # iterate
                q, r, p = qrpa
                @test q*r[:,invperm(p)] ≈ a[1:n1,:]
                # property names
                @test Base.propertynames(qrpa)       == (:R, :Q, :p, :P)
            end
            @testset "(Automatic) Thin (pivoted) QR decomposition" begin
                qrpa  = factorize(a[:,1:n1])
                q,r = qrpa.Q, qrpa.R
                @test_throws FieldError qrpa.Z
                p = qrpa.p
                @test q'*squareQ(q) ≈ Matrix(I, a_1, a_1)
                @test q*squareQ(q)' ≈ Matrix(I, a_1, a_1)
                @test q*r ≈ a[:,p]
                @test q*r[:,invperm(p)] ≈ a[:,1:n1]
                @test Array(qrpa) ≈ a[:,1:5]
                if eltya != Int
                    @test Array{eltya}(q) ≈ Matrix(q)
                end
                @test_throws DimensionMismatch q*b[1:n1+1]
                @test_throws DimensionMismatch b[1:n1+1]*q'
                sq = size(q.factors, 2)
                @test *(UpperTriangular(Matrix{eltyb}(I, sq, sq)), adjoint(q))*squareQ(q) ≈ Matrix(I, n1, a_1) atol=5000ε
                if eltya != Int
                    @test Matrix{eltyb}(I, a_1, a_1)*q ≈ squareQ(convert(LinearAlgebra.AbstractQ{tab},q))
                end
                qrstring = sprint((t, s) -> show(t, "text/plain", s), qrpa)
                rstring  = sprint((t, s) -> show(t, "text/plain", s), r)
                qstring  = sprint((t, s) -> show(t, "text/plain", s), q)
                pstring  = sprint((t, s) -> show(t, "text/plain", s), p)
                @test qrstring == "$(summary(qrpa))\nQ factor: $qstring\nR factor:\n$rstring\npermutation:\n$pstring"
                # iterate
                q, r, p = qrpa
                @test q*r[:,invperm(p)] ≈ a[:,1:n1]
                # property names
                @test Base.propertynames(qrpa)       == (:R, :Q, :p, :P)
            end
        end
        if eltya != Int
            @testset "Matmul with QR factorizations" begin
                a = raw_a
                qrpa = factorize(a[:,1:n1])
                q, r = qrpa.Q, qrpa.R
                @test rmul!(copy(squareQ(q)'), q) ≈ Matrix(I, n, n)
                @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1),q)
                @test rmul!(squareQ(q), adjoint(q)) ≈ Matrix(I, n, n)
                @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1), adjoint(q))
                @test_throws ErrorException size(q,-1)
                @test_throws DimensionMismatch LinearAlgebra.lmul!(q,zeros(eltya,n1+1))
                @test_throws DimensionMismatch LinearAlgebra.lmul!(adjoint(q), zeros(eltya,n1+1))

                b = similar(a); rand!(b)
                c = similar(a)
                d = similar(a[:,1:n1])
                @test mul!(c, q, b) ≈ q*b
                @test mul!(d, q, r) ≈ q*r ≈ a[:,qrpa.p]
                @test mul!(c, q', b) ≈ q'*b
                @test mul!(d, q', a[:,qrpa.p])[1:n1,:] ≈ r
                @test all(x -> abs(x) < ε*norm(a), d[n1+1:end,:])
                @test mul!(c, b, q) ≈ b*q
                @test mul!(c, b, q') ≈ b*q'
                @test_throws DimensionMismatch mul!(Matrix{eltya}(I, n+1, n), q, b)

                qra = qr(a[:,1:n1], NoPivot())
                q, r = qra.Q, qra.R
                @test rmul!(copy(squareQ(q)'), q) ≈ Matrix(I, n, n)
                @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1),q)
                @test rmul!(squareQ(q), adjoint(q)) ≈ Matrix(I, n, n)
                @test_throws DimensionMismatch rmul!(Matrix{eltya}(I, n+1, n+1),adjoint(q))
                @test_throws ErrorException size(q,-1)
                @test_throws DimensionMismatch q * Matrix{Int8}(I, n+4, n+4)

                @test mul!(c, q, b) ≈ q*b
                @test mul!(d, q, r) ≈ a[:,1:n1]
                @test mul!(c, q', b) ≈ q'*b
                @test mul!(d, q', a[:,1:n1])[1:n1,:] ≈ r
                @test all(x -> abs(x) < ε*norm(a), d[n1+1:end,:])
                @test mul!(c, b, q) ≈ b*q
                @test mul!(c, b, q') ≈ b*q'
                @test_throws DimensionMismatch mul!(Matrix{eltya}(I, n+1, n), q, b)

                b = similar(a[:,1]); rand!(b)
                c = similar(a[:,1])
                d = similar(a[:,1])
                @test mul!(c, q, b) ≈ q*b
                @test mul!(c, q', b) ≈ q'*b
                @test_throws DimensionMismatch mul!(Vector{eltya}(undef, n+1), q, b)
            end
        end
    end
end

@testset "transpose errors" begin
    @test_throws ArgumentError transpose(qr(randn(ComplexF64,3,3)))
    @test_throws ArgumentError transpose(qr(randn(ComplexF64,3,3), NoPivot()))
    @test_throws ArgumentError transpose(qr(big.(randn(ComplexF64,3,3))))
end

@testset "Issue 7304" begin
    A = [-√.5 -√.5; -√.5 √.5]
    Q = rectangularQ(qr(A).Q)
    @test norm(A-Q) < eps()
end

@testset "qr on AbstractVector" begin
    vr = [3.0, 4.0]
    for Tr in (Float32, Float64)
        for T in (Tr, Complex{Tr})
            v = convert(Vector{T}, vr)
            nv, nm = qr(v)
            @test norm(nv*Matrix(I, (2,2)) - [-0.6 -0.8; -0.8 0.6], Inf) < eps(Tr)
            @test nm == fill(-5.0, 1, 1)
        end
    end
end

@testset "QR on Ints" begin
    # not sure what to do about this edge case now that we build decompositions
    # for qr(...), so for now just commenting this out
    # @test qr(Int[]) == (Int[],1)

    B = rand(7,2)
    @test (1:7)\B ≈ Vector(1:7)\B
end

@testset "Issue 16520" begin
    @test_throws DimensionMismatch rand(3,2)\(1:5)
end

@testset "Issue 22810" begin
    A = zeros(1, 2)
    B = zeros(1, 1)
    @test A \ B == zeros(2, 1)
    @test qr(A, ColumnNorm()) \ B == zeros(2, 1)
end

@testset "Issue 24107" begin
    A = rand(200,2)
    @test A \ range(0, stop=1, length=200) == A \ Vector(range(0, stop=1, length=200))
end

@testset "Issue 24589. Promotion of rational matrices" begin
    A = rand(1//1:5//5, 4,3)
    @test Matrix(first(qr(A))) == Matrix(first(qr(float(A))))
end

@testset "Issue Test Factorization fallbacks for rectangular problems" begin
    A  = randn(3,2)
    Ac = copy(A')
    b  = randn(3)
    b0 = copy(b)
    c  = randn(2)
    B  = randn(3,3)
    B0 = copy(B)
    C  = randn(2,3)
    @test A \b ≈ ldiv!(c, qr(A ), b)
    @test b == b0
    @test A \B ≈ ldiv!(C, qr(A ), B)
    @test B == B0
    c0 = copy(c)
    C0 = copy(C)
    @test Ac\c ≈ ldiv!(b, qr(Ac, ColumnNorm()), c)
    @test c0 == c
    @test Ac\C ≈ ldiv!(B, qr(Ac, ColumnNorm()), C)
    @test C0 == C
end

@testset "Issue reflector of zero-length vector" begin
    a = [2.0]
    x = view(a,1:0)
    τ = LinearAlgebra.reflector!(view(x,1:0))
    @test τ == 0.0

    b = reshape([3.0],1,1)
    @test isempty(LinearAlgebra.reflectorApply!(x, τ, view(b,1:0,:)))
    @test b[1] == 3.0
end

@testset "det(Q::Union{QRCompactWYQ, QRPackedQ})" begin
    # 40 is the number larger than the default block size 36 of QRCompactWY
    @testset for n in [1:3; 40], m in [1:3; 40], pivot in (NoPivot(), ColumnNorm())
        @testset "real" begin
            @testset for k in 0:min(n, m, 5)
                A = cat(Array(I(k)), randn(n - k, m - k); dims=(1, 2))
                Q, = qr(A, pivot)
                @test det(Q) ≈ det(Q*Matrix(I, size(Q, 1), size(Q, 1)))
                @test abs(det(Q)) ≈ 1
            end
        end
        @testset "complex" begin
            @testset for k in 0:min(n, m, 5)
                A = cat(Array(I(k)), randn(ComplexF64, n - k, m - k); dims=(1, 2))
                Q, = qr(A, pivot)
                @test det(Q) ≈ det(Q*Matrix(I, size(Q, 1), size(Q, 1)))
                @test abs(det(Q)) ≈ 1
            end
        end
    end
end

@testset "inv(::AbstractQ)" begin
    for T in (Float64, ComplexF64)
        Q = qr(randn(T,5,5)).Q
        @test inv(Q) === Q'
        @test inv(Q)' === inv(Q') === Q
    end
end

@testset "QR factorization of Q" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        Q1, R1 = qr(randn(T,5,5))
        Q2, R2 = qr(Q1)
        @test Matrix(Q1) ≈ Matrix(Q2)
        @test R2 ≈ I
    end
end

@testset "Generation of orthogonal matrices" begin
    for T in (Float32, Float64)
        n = 5
        Q, R = qr(randn(T,n,n))
        O = Q * Diagonal(sign.(diag(R)))
        @test O' * O ≈ I
    end
end

@testset "Multiplication of Q by special matrices" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        n = 5
        Q, R = qr(randn(T,n,n))
        Qmat = Matrix(Q)
        D = Diagonal(randn(T,n))
        @test Q * D ≈ Qmat * D
        @test D * Q ≈ D * Qmat
        J = 2*I
        @test Q * J ≈ Qmat * J
        @test J * Q ≈ J * Qmat
    end
end

@testset "copyto! for Q" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        n = 5
        Q, R = qr(randn(T,n,n))
        Qmat = Matrix(Q)
        dest1 = Matrix{T}(undef, size(Q))
        copyto!(dest1, Q)
        @test dest1 ≈ Qmat
        dest2 = PermutedDimsArray(Matrix{T}(undef, size(Q)), (1, 2))
        copyto!(dest2, Q)
        @test dest2 ≈ Qmat
        dest3 = PermutedDimsArray(Matrix{T}(undef, size(Q)), (2, 1))
        copyto!(dest3, Q)
        @test dest3 ≈ Qmat
    end
end

@testset "adjoint of QR" begin
    n = 5
    B = randn(5, 2)

    @testset "size(b)=$(size(b))" for b in (B[:, 1], B)
        @testset "size(A)=$(size(A))" for A in (
            randn(n, n),
            # Wide problems become minimum norm (in x) problems similarly to LQ
            randn(n + 2, n),
            complex.(randn(n, n), randn(n, n)))

            @testset "QRCompactWY" begin
                F = qr(A)
                x = F'\b
                @test x ≈ A'\b
                @test length(size(x)) == length(size(b))
            end

            @testset "QR" begin
                F = LinearAlgebra.qrfactUnblocked!(copy(A))
                x = F'\b
                @test x ≈ A'\b
                @test length(size(x)) == length(size(b))
            end

            @testset "QRPivoted" begin
                F = LinearAlgebra.qr(A, ColumnNorm())
                x = F'\b
                @test x ≈ A'\b
                @test length(size(x)) == length(size(b))
            end
        end
        @test_throws DimensionMismatch("overdetermined systems are not supported")    qr(randn(n - 2, n))'\b
        @test_throws DimensionMismatch("arguments must have the same number of rows") qr(randn(n, n + 1))'\b
        @test_throws DimensionMismatch("overdetermined systems are not supported")    LinearAlgebra.qrfactUnblocked!(randn(n - 2, n))'\b
        @test_throws DimensionMismatch("arguments must have the same number of rows") LinearAlgebra.qrfactUnblocked!(randn(n, n + 1))'\b
        @test_throws DimensionMismatch("overdetermined systems are not supported")    qr(randn(n - 2, n), ColumnNorm())'\b
        @test_throws DimensionMismatch("arguments must have the same number of rows") qr(randn(n, n + 1), ColumnNorm())'\b
    end
end

@testset "issue #38974" begin
    A = qr(ones(3, 1))
    B = I(3)
    C = B*A.Q'
    @test C ≈ A.Q * Matrix(I, 3, 3)
    @test A.Q' * B ≈ A.Q * Matrix(I, 3, 3)
end

@testset "convert between eltypes" begin
    a = rand(Float64, 10, 5)
    qra = qr(a)
    qrwy = LinearAlgebra.QRCompactWY{Float32}(qra.factors, qra.T)
    @test Array(qrwy) ≈ Array(qr(Float32.(a)))
    @test eltype(qrwy.factors) == eltype(qrwy.T) == Float32
    qra = qr(a, ColumnNorm())
    qrp = QRPivoted{Float32}(qra.factors, qra.τ, qra.jpvt)
    @test Array(qrp) ≈ Array(qr(Float32.(a), ColumnNorm()))
    @test eltype(qrp.factors) == eltype(qrp.τ) == Float32
    a = rand(Float16, 10, 5)
    qra = qr(a)
    qrnonblas = QR{ComplexF16}(qra.factors, qra.τ)
    @test Array(qrnonblas) ≈ Array(qr(ComplexF16.(a)))
    @test eltype(qrnonblas.factors) == eltype(qrnonblas.τ) == ComplexF16
end

# We use approximate equals to get MKL.jl tests to pass.
@testset "optimized getindex for an AbstractQ" begin
    for T in [Float64, ComplexF64]
        Q = qr(rand(T, 4, 4))
        Q2 = Q.Q
        M = Matrix(Q2)
        for j in axes(M, 2)
            @test Q2[:, j] ≈ M[:, j]
            for i in axes(M, 1)
                @test Q2[i, :] ≈ M[i, :]
                @test Q2[i, j] ≈ M[i, j]
            end
        end
        @test Q2[:] ≈ M[:]
        @test Q2[:, :] ≈ M[:, :]
        @test Q2[:, :, :] ≈ M[:, :, :]
    end
    # Check that getindex works if copy returns itself (#44729)
    struct MyIdentity{T} <: LinearAlgebra.AbstractQ{T} end
    Base.size(::MyIdentity, dim::Integer) = dim in (1,2) ? 2 : 1
    Base.size(::MyIdentity) = (2, 2)
    Base.copy(J::MyIdentity) = J
    LinearAlgebra.lmul!(::MyIdentity{T}, M::Array{T}) where {T} = M
    @test MyIdentity{Float64}()[1,:] == [1.0, 0.0]
end

@testset "issue #48911" begin
    # testcase in the original issue
    # test ldiv!(::QRPivoted, ::AbstractVector)
    A = Complex{BigFloat}[1+im 1-im]
    b = Complex{BigFloat}[3+im]
    x = A\b
    AF = Complex{Float64}[1+im 1-im]
    bf = Complex{Float64}[3+im]
    xf = AF\bf
    @test x ≈ xf

    # test ldiv!(::QRPivoted, ::AbstractVector)
    A = Complex{BigFloat}[1+im 2-2im 3+3im; 4-4im 5+5im 6-6im]
    b = Complex{BigFloat}[1+im; 0]
    x = A\b
    AF = Complex{Float64}[1+im 2-2im 3+3im; 4-4im 5+5im 6-6im]
    bf = Complex{Float64}[1+im; 0]
    xf = AF\bf
    @test x ≈ xf

    # test ldiv!(::QRPivoted, ::AbstractMatrix)
    C = Complex{BigFloat}[1+im 2-2im 3+3im; 4-4im 5+5im 6-6im]
    D = Complex{BigFloat}[1+im 1-im; 0 0]
    x = C\D
    CF = Complex{Float64}[1+im 2-2im 3+3im; 4-4im 5+5im 6-6im]
    DF = Complex{Float64}[1+im 1-im; 0 0]
    xf = CF\DF
    @test x ≈ xf
end

@testset "issue #53451" begin
    # in the issue it was noted that QR factorizations of zero-column matrices
    # were possible, but zero row-matrices errored, because LAPACK does not
    # accept these empty matrices. now, the `geqrt!` call should be forwarded only
    # if both matrix dimensions are positive.

    for dimA in (0, 1, 2, 4)
        for F in (Float32, Float64, ComplexF32, ComplexF64, BigFloat)
            # this should have worked before, Q is square, and R is 0 × 0:
            A_zero_cols = rand(F, dimA, 0)
            qr_zero_cols = qr(A_zero_cols)
            @test size(qr_zero_cols.Q) == (dimA, dimA)
            @test size(qr_zero_cols.R) == (0, 0)
            @test qr_zero_cols.Q == LinearAlgebra.I(dimA)

            # this should work now, Q is 0 × 0, and R has `dimA` columns:
            A_zero_rows = rand(F, 0, dimA)
            qr_zero_rows = qr(A_zero_rows)
            @test size(qr_zero_rows.Q) == (0, 0)
            @test size(qr_zero_rows.R) == (0, dimA)
        end
    end
end

@testset "issue #53214" begin
    # Test that the rank of a QRPivoted matrix is computed correctly
    @test rank(qr([1.0 0.0; 0.0 1.0], ColumnNorm())) == 2
    @test rank(qr([1.0 0.0; 0.0 0.9], ColumnNorm()), rtol=0.95) == 1
    @test rank(qr([1.0 0.0; 0.0 0.9], ColumnNorm()), atol=0.95) == 1
    @test rank(qr([1.0 0.0; 0.0 1.0], ColumnNorm()), rtol=1.01) == 0
    @test rank(qr([1.0 0.0; 0.0 1.0], ColumnNorm()), atol=1.01) == 0

    @test rank(qr([1.0 2.0; 2.0 4.0], ColumnNorm())) == 1
    @test rank(qr([1.0 2.0 3.0; 4.0 5.0 6.0 ; 7.0 8.0 9.0], ColumnNorm())) == 2
end

end # module TestQR
