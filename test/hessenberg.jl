# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestHessenberg

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.SizedArrays
using Main.LinearAlgebraTestHelpers.ImmutableArrays

# for tuple tests below
≅(x,y) = all(p -> p[1] ≈ p[2], zip(x,y))

let n = 10
    Random.seed!(1234321)

    Areal  = randn(n,n)/2
    Aimg   = randn(n,n)/2
    b_ = randn(n)
    B_ = randn(n,3)

    # UpperHessenberg methods not covered by the tests below
    @testset "UpperHessenberg" begin
        A = Areal
        H = UpperHessenberg(A)
        AH = triu(A,-1)
        for k in -2:2
            @test istril(H, k) == istril(AH, k)
            @test istriu(H, k) == istriu(AH, k)
            @test (k <= -1 ? istriu(H, k) : !istriu(H, k))
        end
        @test UpperHessenberg(H) === H
        @test parent(H) === A
        @test Matrix(H) == Array(H) == H == AH
        @test real(H) == real(AH)
        @test real(UpperHessenberg{ComplexF64}(A)) == H
        @test real(UpperHessenberg{ComplexF64}(H)) == H
        sim = similar(H, ComplexF64)
        @test sim isa UpperHessenberg{ComplexF64}
        @test size(sim) == size(H)
        for x in (2,2+3im)
            @test x*H == H*x == x*AH
            for op in (+,-)
                @test op(H,x*I) == op(AH,x*I) == op(op(x*I,H))
                @test op(H,x*I)*x == op(AH,x*I)*x == x*op(H,x*I)
            end
        end
        @test [H[i,j] for i=1:size(H,1), j=1:size(H,2)] == triu(A,-1)
        H1 = LinearAlgebra.fillstored!(copy(H), 1)
        @test H1 == triu(fill(1, n,n), -1)
        @test tril(H1.data,-2) == tril(H.data,-2)
        A2, H2 = copy(A), copy(H)
        A2[1:4,3]=H2[1:4,3]=1:4
        H2[5,3]=0
        @test H2 == triu(A2,-1)
        @test_throws ArgumentError H[5,3]=1
        Hc = UpperHessenberg(Areal + im .* Aimg)
        AHc = triu(Areal + im .* Aimg,-1)
        @test real(Hc) == real(AHc)
        @test imag(Hc) == imag(AHc)
        @test Array(copy(adjoint(Hc))) == adjoint(Array(Hc))
        @test Array(copy(transpose(Hc))) == transpose(Array(Hc))
        @test rmul!(copy(Hc), 2.0) == lmul!(2.0, copy(Hc))
        H = UpperHessenberg(Areal)
        @test Array(Hc + H) == Array(Hc) + Array(H)
        @test Array(Hc - H) == Array(Hc) - Array(H)
        @testset "ldiv and rdiv" begin
            for b in (b_, B_), H in (H, Hc, H', Hc', transpose(Hc))
                @test H * (H \ b) ≈ b
                @test (b' / H) * H ≈ (Matrix(b') / H) * H ≈ b'
                @test (transpose(b) / H) * H ≈ (Matrix(transpose(b)) / H) * H ≈ transpose(b)
            end
        end
        @testset "Preserve UpperHessenberg shape (issue #39388)" begin
            H = UpperHessenberg(Areal)
            A = rand(n,n)
            d = rand(n)
            dl = rand(n-1)
            du = rand(n-1)
            us = 1*I
            @testset "$op" for op = (+,-)
                for x = (us, Diagonal(d), Bidiagonal(d,dl,:U), Bidiagonal(d,dl,:L),
                            Tridiagonal(dl,d,du), SymTridiagonal(d,dl),
                            UpperTriangular(A), UnitUpperTriangular(A))
                    @test op(H,x) == op(Array(H),x)
                    @test op(x,H) == op(x,Array(H))
                    @test op(H,x) isa UpperHessenberg
                    @test op(x,H) isa UpperHessenberg
                end
            end
            H = UpperHessenberg(Areal)
            A = randn(n,n)
            d = randn(n)
            dl = randn(n-1)
            @testset "Multiplication/division" begin
                for x = (5, 5I, Diagonal(d), Bidiagonal(d,dl,:U),
                            UpperTriangular(A), UnitUpperTriangular(A))
                    @test (H*x)::UpperHessenberg ≈ Array(H)*x
                    @test (x*H)::UpperHessenberg ≈ x*Array(H)
                    @test H/x ≈ Array(H)/x
                    @test x\H ≈ x\Array(H)
                    @test H/x isa UpperHessenberg
                    @test x\H isa UpperHessenberg
                end
                x = Bidiagonal(d, dl, :L)
                @test H*x == Array(H)*x
                @test x*H == x*Array(H)
                @test H/x == Array(H)/x
                @test x\H == x\Array(H)
            end
        end
    end

    @testset for eltya in (Float32, Float64, ComplexF32, ComplexF64, Int), herm in (false, true)
        A_ = eltya == Int ?
                rand(1:7, n, n) :
                convert(Matrix{eltya}, eltya <: Complex ?
                    complex.(Areal, Aimg) :
                    Areal)
        A = herm ? Hermitian(A_ + A_') : A_

        H = hessenberg(A)
        @test Hessenberg(H) === H
        eltyh = eltype(H)
        @test size(H.Q, 1) == size(A, 1)
        @test size(H.Q, 2) == size(A, 2)
        @test size(H.Q) == size(A)
        @test size(H) == size(A)
        @test_throws FieldError H.Z
        @test convert(Array, H) ≈ A
        @test (H.Q * H.H) * H.Q' ≈ A ≈ (Matrix(H.Q) * Matrix(H.H)) * Matrix(H.Q)'
        @test (H.Q' * A) * H.Q ≈ H.H
        #getindex for HessenbergQ
        @test H.Q[1,1] ≈ Array(H.Q)[1,1]
        @test det(H.Q) ≈ det(Matrix(H.Q))
        @test logabsdet(H.Q)[1] ≈ logabsdet(Matrix(H.Q))[1] atol=2n*eps(float(real(eltya)))

        # REPL show
        hessstring = sprint((t, s) -> show(t, "text/plain", s), H)
        qstring = sprint((t, s) -> show(t, "text/plain", s), H.Q)
        hstring = sprint((t, s) -> show(t, "text/plain", s), H.H)
        @test hessstring == "$(summary(H))\nQ factor: $qstring\nH factor:\n$hstring"

        #iterate
        q,h = H
        @test q == H.Q
        @test h == H.H

        @test convert(Array, 2 * H) ≈ 2 * A ≈ convert(Array, H * 2)
        @test convert(Array, H + 2I) ≈ A + 2I ≈ convert(Array, 2I + H)
        @test convert(Array, H + (2+4im)I) ≈ A + (2+4im)I ≈ convert(Array, (2+4im)I + H)
        @test convert(Array, H - 2I) ≈ A - 2I ≈ -convert(Array, 2I - H)
        @test convert(Array, -H) == -convert(Array, H)
        @test convert(Array, 2*(H + (2+4im)I)) ≈ 2A + (4+8im)I

        b = convert(Vector{eltype(H)}, b_)
        B = convert(Matrix{eltype(H)}, B_)
        @test H \ b ≈ A \ b ≈ H \ complex(b)
        @test H \ B ≈ A \ B ≈ H \ complex(B)
        @test (H - I) \ B ≈ (A - I) \ B
        @test (H - (3+4im)I) \ B ≈ (A - (3+4im)I) \ B
        @test b' / H ≈ b' / A ≈ complex(b') / H
        @test transpose(b) / H ≈ transpose(b) / A ≈ transpose(complex(b)) / H
        @test B' / H ≈ B' / A ≈ complex(B') / H
        @test b' / H' ≈ complex(b)' / H'
        @test B' / (H - I) ≈ B' / (A - I)
        @test B' / (H - (3+4im)I) ≈ B' / (A - (3+4im)I)
        @test (H - (3+4im)I)' \ B ≈ (A - (3+4im)I)' \ B
        @test B' / (H - (3+4im)I)' ≈ B' / (A - (3+4im)I)'

        for shift in (0,1,3+4im)
            @test det(H + shift*I) ≈ det(A + shift*I)
            @test logabsdet(H + shift*I) ≅ logabsdet(A + shift*I)
        end

        HM = Matrix(h)
        @test dot(b, h, b) ≈ dot(h'b, b) ≈ dot(b, HM, b) ≈ dot(HM'b, b)
        c = b .+ 1
        @test dot(b, h, c) ≈ dot(h'b, c) ≈ dot(b, HM, c) ≈ dot(HM'b, c)
    end
end

@testset "Reverse operation on UpperHessenberg" begin
    A = UpperHessenberg(randn(5, 5))
    @test reverse(A, dims=1) == reverse(Matrix(A), dims=1)
    @test reverse(A, dims=2) == reverse(Matrix(A), dims=2)
    @test reverse(A) == reverse(Matrix(A))
end

@testset "hessenberg(::AbstractMatrix)" begin
    n = 10
    A = Tridiagonal(rand(n-1), rand(n), rand(n-1))
    H = hessenberg(A)
    @test convert(Array, H) ≈ A
end

# check logdet on a matrix that has a positive determinant
let A = [0.5 0.1 0.9 0.4; 0.9 0.7 0.5 0.4; 0.3 0.4 0.9 0.0; 0.4 0.0 0.0 0.5]
    @test logdet(hessenberg(A)) ≈ logdet(A) ≈ -3.5065578973199822
end

@testset "Base.propertynames" begin
    F =  hessenberg([4. 9. 7.; 4. 4. 1.; 4. 3. 2.])
    @test Base.propertynames(F) == (:Q, :H, :μ)
    @test Base.propertynames(F, true) == (:Q, :H, :μ, :τ, :factors, :uplo)
end

@testset "adjoint of Hessenberg" begin
    Ar = randn(5, 5)
    Ac = complex.(randn(5, 5), randn(5, 5))
    b = ones(size(Ar, 1))

    for A in (Ar, Ac)
        F = hessenberg(A)
        @test A'\b ≈ F'\b
    end
end

@testset "Conversion to AbstractArray" begin
    # tests corresponding to #34995
    A = ImmutableArray([1 2 3; 4 5 6; 7 8 9])
    H = UpperHessenberg(A)

    @test convert(AbstractArray{Float64}, H)::UpperHessenberg{Float64,ImmutableArray{Float64,2,Array{Float64,2}}} == H
    @test convert(AbstractMatrix{Float64}, H)::UpperHessenberg{Float64,ImmutableArray{Float64,2,Array{Float64,2}}} == H
end

@testset "custom axes" begin
    SZA = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
    S = UpperHessenberg(SZA)
    r = SizedArrays.SOneTo(2)
    @test axes(S) === (r,r)
end

@testset "copyto! with aliasing (#39460)" begin
    M = Matrix(reshape(1:36, 6, 6))
    A = UpperHessenberg(view(M, 1:5, 1:5))
    A2 = copy(A)
    B = UpperHessenberg(view(M, 2:6, 2:6))
    @test copyto!(B, A) == A2
end

@testset "getindex with Integers" begin
    M = reshape(1:9, 3, 3)
    S = UpperHessenberg(M)
    @test_throws "invalid index" S[3, true]
    @test S[1,2] == S[Int8(1),UInt16(2)] == S[big(1), Int16(2)]
end

@testset "complex Symmetric" begin
    D = diagm(0=>ComplexF64[1,2])
    S = Symmetric(D)
    H = hessenberg(S)
    @test H.H == D
end

@testset "istriu/istril forwards to parent" begin
    n = 10
    @testset "$(nameof(typeof(M)))" for M in [Tridiagonal(rand(n-1), rand(n), rand(n-1)),
                Tridiagonal(zeros(n-1), zeros(n), zeros(n-1)),
                Diagonal(randn(n)),
                Diagonal(zeros(n)),
                ]
        U = UpperHessenberg(M)
        A = Array(U)
        for k in -n:n
            @test istriu(U, k) == istriu(A, k)
            @test istril(U, k) == istril(A, k)
        end
    end
    z = zeros(n,n)
    P = Matrix{BigFloat}(undef, n, n)
    copytrito!(P, z, 'U')
    P[diagind(P,-1)] .= 0
    U = UpperHessenberg(P)
    A = Array(U)
    @testset for k in -n:n
        @test istriu(U, k) == istriu(A, k)
        @test istril(U, k) == istril(A, k)
    end
end

@testset "Hessenberg factorization of Q" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        Q1, H1 = hessenberg(randn(T,5,5))
        Q2, H2 = hessenberg(Q1)
        @test Q2'Q2 ≈ I
    end
end

@testset "multiplication with empty HessenbergQ" begin
    @test ones(2, 0)*hessenberg(zeros(0,0)).Q == zeros(2,0)
    @test_throws DimensionMismatch ones(2, 1)*hessenberg(zeros(0,0)).Q
    @test hessenberg(zeros(0,0)).Q * ones(0, 2) == zeros(0,2)
    @test_throws DimensionMismatch hessenberg(zeros(0,0)).Q * ones(1, 2)
end

@testset "fillband" begin
    U = UpperHessenberg(zeros(4,4))
    @test_throws ArgumentError LinearAlgebra.fillband!(U, 1, -2, 1)
    @test iszero(U)

    LinearAlgebra.fillband!(U, 10, -1, 2)
    @test all(==(10), diagview(U,-1))
    @test all(==(10), diagview(U,2))
    @test all(==(0), diagview(U,3))

    LinearAlgebra.fillband!(U, 0, -5, 5)
    @test iszero(U)

    U2 = copy(U)
    LinearAlgebra.fillband!(U, -10, 1, -2)
    @test U == U2
    LinearAlgebra.fillband!(U, -10, 10, 10)
    @test U == U2
end

end # module TestHessenberg
