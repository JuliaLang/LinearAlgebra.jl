# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTridiagonal

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.Quaternions
using Main.LinearAlgebraTestHelpers.InfiniteArrays
using Main.LinearAlgebraTestHelpers.FillArrays
using Main.LinearAlgebraTestHelpers.OffsetArrays
using Main.LinearAlgebraTestHelpers.SizedArrays
using Main.LinearAlgebraTestHelpers.ImmutableArrays

include("testutils.jl") # test_approx_eq_modphase

#Test equivalence of eigenvectors/singular vectors taking into account possible phase (sign) differences
function test_approx_eq_vecs(a::StridedVecOrMat{S}, b::StridedVecOrMat{T}, error=nothing) where {S<:Real,T<:Real}
    n = size(a, 1)
    @test n==size(b,1) && size(a,2)==size(b,2)
    error===nothing && (error=n^3*(eps(S)+eps(T)))
    for i=1:n
        ev1, ev2 = a[:,i], b[:,i]
        deviation = min(abs(norm(ev1-ev2)),abs(norm(ev1+ev2)))
        if !isnan(deviation)
            @test deviation ≈ 0.0 atol=error
        end
    end
end

@testset for elty in (Float32, Float64, ComplexF32, ComplexF64, Int)
    n = 12 #Size of matrix problem to test
    Random.seed!(123)
    if elty == Int
        Random.seed!(61516384)
        d = rand(1:100, n)
        dl = -rand(0:10, n-1)
        du = -rand(0:10, n-1)
        v = rand(1:100, n)
        B = rand(1:100, n, 2)
        a = rand(1:100, n-1)
        b = rand(1:100, n)
        c = rand(1:100, n-1)
    else
        d = convert(Vector{elty}, 1 .+ randn(n))
        dl = convert(Vector{elty}, randn(n - 1))
        du = convert(Vector{elty}, randn(n - 1))
        v = convert(Vector{elty}, randn(n))
        B = convert(Matrix{elty}, randn(n, 2))
        a = convert(Vector{elty}, randn(n - 1))
        b = convert(Vector{elty}, randn(n))
        c = convert(Vector{elty}, randn(n - 1))
        if elty <: Complex
            a += im*convert(Vector{elty}, randn(n - 1))
            b += im*convert(Vector{elty}, randn(n))
            c += im*convert(Vector{elty}, randn(n - 1))
        end
    end
    @test_throws DimensionMismatch SymTridiagonal(dl, fill(elty(1), n+1))
    @test_throws ArgumentError SymTridiagonal(rand(n, n))
    @test_throws ArgumentError Tridiagonal(dl, dl, dl)
    @test_throws ArgumentError convert(SymTridiagonal{elty}, Tridiagonal(dl, d, du))

    if elty != Int
        @testset "issue #1490" begin
            @test det(fill(elty(1),3,3)) ≈ zero(elty) atol=3*eps(real(one(elty)))
            @test det(SymTridiagonal(elty[],elty[])) == one(elty)
        end
    end

    @testset "constructor" begin
        for (x, y) in ((d, dl), (GenericArray(d), GenericArray(dl)))
            ST = (SymTridiagonal(x, y))::SymTridiagonal{elty, typeof(x)}
            @test ST == Matrix(ST)
            @test ST.dv === x
            @test ST.ev === y
            @test typeof(ST)(ST) === ST
            TT = (Tridiagonal(y, x, y))::Tridiagonal{elty, typeof(x)}
            @test TT == Matrix(TT)
            @test TT.dl === y
            @test TT.d  === x
            @test TT.du == y
            @test typeof(TT)(TT) === TT
        end
        ST = SymTridiagonal{elty}([1,2,3,4], [1,2,3])
        @test eltype(ST) == elty
        @test SymTridiagonal{elty, Vector{elty}}(ST) === ST
        @test SymTridiagonal{Int64, Vector{Int64}}(ST) isa SymTridiagonal{Int64, Vector{Int64}}
        TT = Tridiagonal{elty}([1,2,3], [1,2,3,4], [1,2,3])
        @test eltype(TT) == elty
        ST = SymTridiagonal{elty,Vector{elty}}(d, GenericArray(dl))
        @test isa(ST, SymTridiagonal{elty,Vector{elty}})
        TT = Tridiagonal{elty,Vector{elty}}(GenericArray(dl), d, GenericArray(dl))
        @test isa(TT, Tridiagonal{elty,Vector{elty}})
        @test_throws ArgumentError SymTridiagonal(d, GenericArray(dl))
        @test_throws ArgumentError SymTridiagonal(GenericArray(d), dl)
        @test_throws ArgumentError Tridiagonal(GenericArray(dl), d, GenericArray(dl))
        @test_throws ArgumentError Tridiagonal(dl, GenericArray(d), dl)
        @test_throws ArgumentError SymTridiagonal{elty}(d, GenericArray(dl))
        @test_throws ArgumentError Tridiagonal{elty}(GenericArray(dl), d,GenericArray(dl))
        STI = SymTridiagonal([1,2,3,4], [1,2,3])
        TTI = Tridiagonal([1,2,3], [1,2,3,4], [1,2,3])
        TTI2 = Tridiagonal([1,2,3], [1,2,3,4], [1,2,3], [1,2])
        @test SymTridiagonal(STI) === STI
        @test Tridiagonal(TTI)    === TTI
        @test Tridiagonal(TTI2)   === TTI2
        @test isa(SymTridiagonal{elty}(STI), SymTridiagonal{elty})
        @test isa(Tridiagonal{elty}(TTI), Tridiagonal{elty})
        TTI2y = Tridiagonal{elty}(TTI2)
        @test isa(TTI2y, Tridiagonal{elty})
        @test TTI2y.du2 == convert(Vector{elty}, [1,2])
    end
    @testset "interconversion of Tridiagonal and SymTridiagonal" begin
        @test Tridiagonal(dl, d, dl) == SymTridiagonal(d, dl)
        @test SymTridiagonal(d, dl) == Tridiagonal(dl, d, dl)
        @test Tridiagonal(dl, d, du) + Tridiagonal(du, d, dl) == SymTridiagonal(2d, dl+du)
        @test SymTridiagonal(d, dl) + Tridiagonal(dl, d, du) == Tridiagonal(dl + dl, d+d, dl+du)
        @test convert(SymTridiagonal,Tridiagonal(SymTridiagonal(d, dl))) == SymTridiagonal(d, dl)
        @test Array(convert(SymTridiagonal{ComplexF32},Tridiagonal(SymTridiagonal(d, dl)))) == convert(Matrix{ComplexF32}, SymTridiagonal(d, dl))
    end
    @testset "tril/triu" begin
        zerosd = fill!(similar(d), 0)
        zerosdl = fill!(similar(dl), 0)
        zerosdu = fill!(similar(du), 0)
        @test_throws ArgumentError tril!(SymTridiagonal(d, dl), -n - 2)
        @test_throws ArgumentError tril!(SymTridiagonal(d, dl), n)
        @test_throws ArgumentError tril!(Tridiagonal(dl, d, du), -n - 2)
        @test_throws ArgumentError tril!(Tridiagonal(dl, d, du), n)
        @test @inferred(tril(SymTridiagonal(d,dl)))    == Tridiagonal(dl,d,zerosdl)
        @test @inferred(tril(SymTridiagonal(d,dl),1))  == Tridiagonal(dl,d,dl)
        @test @inferred(tril(SymTridiagonal(d,dl),-1)) == Tridiagonal(dl,zerosd,zerosdl)
        @test @inferred(tril(SymTridiagonal(d,dl),-2)) == Tridiagonal(zerosdl,zerosd,zerosdl)
        @test @inferred(tril(Tridiagonal(dl,d,du)))    == Tridiagonal(dl,d,zerosdu)
        @test @inferred(tril(Tridiagonal(dl,d,du),1))  == Tridiagonal(dl,d,du)
        @test @inferred(tril(Tridiagonal(dl,d,du),-1)) == Tridiagonal(dl,zerosd,zerosdu)
        @test @inferred(tril(Tridiagonal(dl,d,du),-2)) == Tridiagonal(zerosdl,zerosd,zerosdu)
        @test @inferred(tril!(copy(SymTridiagonal(d,dl))))    == Tridiagonal(dl,d,zerosdl)
        @test @inferred(tril!(copy(SymTridiagonal(d,dl)),1))  == Tridiagonal(dl,d,dl)
        @test @inferred(tril!(copy(SymTridiagonal(d,dl)),-1)) == Tridiagonal(dl,zerosd,zerosdl)
        @test @inferred(tril!(copy(SymTridiagonal(d,dl)),-2)) == Tridiagonal(zerosdl,zerosd,zerosdl)
        @test @inferred(tril!(copy(Tridiagonal(dl,d,du))))    == Tridiagonal(dl,d,zerosdu)
        @test @inferred(tril!(copy(Tridiagonal(dl,d,du)),1))  == Tridiagonal(dl,d,du)
        @test @inferred(tril!(copy(Tridiagonal(dl,d,du)),-1)) == Tridiagonal(dl,zerosd,zerosdu)
        @test @inferred(tril!(copy(Tridiagonal(dl,d,du)),-2)) == Tridiagonal(zerosdl,zerosd,zerosdu)

        @test_throws ArgumentError triu!(SymTridiagonal(d, dl), -n)
        @test_throws ArgumentError triu!(SymTridiagonal(d, dl), n + 2)
        @test_throws ArgumentError triu!(Tridiagonal(dl, d, du), -n)
        @test_throws ArgumentError triu!(Tridiagonal(dl, d, du), n + 2)
        @test @inferred(triu(SymTridiagonal(d,dl)))    == Tridiagonal(zerosdl,d,dl)
        @test @inferred(triu(SymTridiagonal(d,dl),-1)) == Tridiagonal(dl,d,dl)
        @test @inferred(triu(SymTridiagonal(d,dl),1))  == Tridiagonal(zerosdl,zerosd,dl)
        @test @inferred(triu(SymTridiagonal(d,dl),2))  == Tridiagonal(zerosdl,zerosd,zerosdl)
        @test @inferred(triu(Tridiagonal(dl,d,du)))    == Tridiagonal(zerosdl,d,du)
        @test @inferred(triu(Tridiagonal(dl,d,du),-1)) == Tridiagonal(dl,d,du)
        @test @inferred(triu(Tridiagonal(dl,d,du),1))  == Tridiagonal(zerosdl,zerosd,du)
        @test @inferred(triu(Tridiagonal(dl,d,du),2))  == Tridiagonal(zerosdl,zerosd,zerosdu)
        @test @inferred(triu!(copy(SymTridiagonal(d,dl))))    == Tridiagonal(zerosdl,d,dl)
        @test @inferred(triu!(copy(SymTridiagonal(d,dl)),-1)) == Tridiagonal(dl,d,dl)
        @test @inferred(triu!(copy(SymTridiagonal(d,dl)),1))  == Tridiagonal(zerosdl,zerosd,dl)
        @test @inferred(triu!(copy(SymTridiagonal(d,dl)),2))  == Tridiagonal(zerosdl,zerosd,zerosdl)
        @test @inferred(triu!(copy(Tridiagonal(dl,d,du))))    == Tridiagonal(zerosdl,d,du)
        @test @inferred(triu!(copy(Tridiagonal(dl,d,du)),-1)) == Tridiagonal(dl,d,du)
        @test @inferred(triu!(copy(Tridiagonal(dl,d,du)),1))  == Tridiagonal(zerosdl,zerosd,du)
        @test @inferred(triu!(copy(Tridiagonal(dl,d,du)),2))  == Tridiagonal(zerosdl,zerosd,zerosdu)

        @test !istril(SymTridiagonal(d,dl))
        @test istril(SymTridiagonal(d,zerosdl))
        @test !istril(SymTridiagonal(d,dl),-2)
        @test !istriu(SymTridiagonal(d,dl))
        @test istriu(SymTridiagonal(d,zerosdl))
        @test !istriu(SymTridiagonal(d,dl),2)
        @test istriu(Tridiagonal(zerosdl,d,du))
        @test !istriu(Tridiagonal(dl,d,zerosdu))
        @test istriu(Tridiagonal(zerosdl,zerosd,du),1)
        @test !istriu(Tridiagonal(dl,d,zerosdu),2)
        @test istril(Tridiagonal(dl,d,zerosdu))
        @test !istril(Tridiagonal(zerosdl,d,du))
        @test istril(Tridiagonal(dl,zerosd,zerosdu),-1)
        @test !istril(Tridiagonal(dl,d,zerosdu),-2)

        @test isdiag(SymTridiagonal(d,zerosdl))
        @test !isdiag(SymTridiagonal(d,dl))
        @test isdiag(Tridiagonal(zerosdl,d,zerosdu))
        @test !isdiag(Tridiagonal(dl,d,zerosdu))
        @test !isdiag(Tridiagonal(zerosdl,d,du))
        @test !isdiag(Tridiagonal(dl,d,du))

        # Test methods that could fail due to dv and ev having the same length
        # see #41089

        badev = zero(d)
        badev[end] = 1
        S = SymTridiagonal(d, badev)

        @test istriu(S, -2)
        @test istriu(S, 0)
        @test !istriu(S, 2)

        @test isdiag(S)
    end

    @testset "iszero and isone" begin
        Tzero = Tridiagonal(zeros(elty, 9), zeros(elty, 10), zeros(elty, 9))
        Tone = Tridiagonal(zeros(elty, 9), ones(elty, 10), zeros(elty, 9))
        Tmix = Tridiagonal(zeros(elty, 9), zeros(elty, 10), zeros(elty, 9))
        Tmix[end, end] = one(elty)

        Szero = SymTridiagonal(zeros(elty, 10), zeros(elty, 9))
        Sone = SymTridiagonal(ones(elty, 10), zeros(elty, 9))
        Smix = SymTridiagonal(zeros(elty, 10), zeros(elty, 9))
        Smix[end, end] = one(elty)

        @test iszero(Tzero)
        @test !isone(Tzero)
        @test !iszero(Tone)
        @test isone(Tone)
        @test !iszero(Tmix)
        @test !isone(Tmix)

        @test iszero(Szero)
        @test !isone(Szero)
        @test !iszero(Sone)
        @test isone(Sone)
        @test !iszero(Smix)
        @test !isone(Smix)

        badev = zeros(elty, 3)
        badev[end] = 1

        @test isone(SymTridiagonal(ones(elty, 3), badev))
        @test iszero(SymTridiagonal(zeros(elty, 3), badev))
    end

    @testset for mat_type in (Tridiagonal, SymTridiagonal)
        A = mat_type == Tridiagonal ? mat_type(dl, d, du) : mat_type(d, dl)
        fA = map(elty <: Complex ? ComplexF64 : Float64, Array(A))
        @testset "similar, size, and copyto!" begin
            B = similar(A)
            @test size(B) == size(A)
            copyto!(B, A)
            @test B == A
            @test isa(similar(A), mat_type{elty})
            @test isa(similar(A, Int), mat_type{Int})
            @test isa(similar(A, (3, 2)), Matrix)
            @test isa(similar(A, Int, (3, 2)), Matrix{Int})
            @test size(A, 3) == 1
            @test size(A, 1) == n
            @test size(A) == (n, n)
            @test_throws BoundsError size(A, 0)
        end
        @testset "getindex" begin
            @test_throws BoundsError A[n + 1, 1]
            @test_throws BoundsError A[1, n + 1]
            @test A[1, n] == convert(elty, 0.0)
            @test A[1, 1] == d[1]
        end
        @testset "setindex!" begin
            @test_throws BoundsError A[n + 1, 1] = 0 # test bounds check
            @test_throws BoundsError A[1, n + 1] = 0 # test bounds check
            @test_throws ArgumentError A[1, 3]   = 1 # test assignment off the main/sub/super diagonal
            if mat_type == Tridiagonal
                @test (A[3, 3] = A[3, 3]; A == fA) # test assignment on the main diagonal
                @test (A[3, 2] = A[3, 2]; A == fA) # test assignment on the subdiagonal
                @test (A[2, 3] = A[2, 3]; A == fA) # test assignment on the superdiagonal
                @test ((A[1, 3] = 0) == 0; A == fA) # test zero assignment off the main/sub/super diagonal
            else # mat_type is SymTridiagonal
                @test ((A[3, 3] = A[3, 3]) == A[3, 3]; A == fA) # test assignment on the main diagonal
                @test_throws ArgumentError A[3, 2] = 1 # test assignment on the subdiagonal
                @test_throws ArgumentError A[2, 3] = 1 # test assignment on the superdiagonal
            end
            # setindex! should return the destination
            @test setindex!(A, A[2,2], 2, 2) === A
        end
        @testset "diag" begin
            @test (@inferred diag(A))::typeof(d) == d
            @test (@inferred diag(A, 0))::typeof(d) == d
            @test (@inferred diag(A, 1))::typeof(d) == (mat_type == Tridiagonal ? du : dl)
            @test (@inferred diag(A, -1))::typeof(d) == dl
            @test (@inferred diag(A, n-1))::typeof(d) == zeros(elty, 1)
            @test isempty(@inferred diag(A, -n - 1))
            @test isempty(@inferred diag(A, n + 1))
            GA = mat_type == Tridiagonal ? mat_type(GenericArray.((dl, d, du))...) : mat_type(GenericArray.((d, dl))...)
            @test (@inferred diag(GA))::typeof(GenericArray(d)) == GenericArray(d)
            @test (@inferred diag(GA, -1))::typeof(GenericArray(d)) == GenericArray(dl)
        end
        @testset "trace" begin
            if real(elty) <: Integer
                @test tr(A) == tr(fA)
            else
                @test tr(A) ≈ tr(fA) rtol=2eps(real(elty))
            end
        end
        @testset "Idempotent tests" begin
            for func in (conj, transpose, adjoint)
                @test func(func(A)) == A
                if func ∈ (transpose, adjoint)
                    @test func(func(A)) === A
                end
            end
        end
        @testset "permutedims(::[Sym]Tridiagonal)" begin
            @test permutedims(permutedims(A)) === A
            @test permutedims(A) == transpose.(transpose(A))
            @test permutedims(A, [1, 2]) === A
            @test permutedims(A, (2, 1)) == permutedims(A)
        end
        if elty != Int
            @testset "Simple unary functions" begin
                for func in (det, inv)
                    @test func(A) ≈ func(fA) atol=n^2*sqrt(eps(real(one(elty))))
                end
            end
        end
        ds = mat_type == Tridiagonal ? (dl, d, du) : (d, dl)
        for f in (real, imag)
            @test f(A)::mat_type == mat_type(map(f, ds)...)
        end
        if elty <: Real
            for f in (round, trunc, floor, ceil)
                fds = [f.(d) for d in ds]
                @test f.(A)::mat_type == mat_type(fds...)
                @test f.(Int, A)::mat_type == f.(Int, fA)
            end
        end
        fds = [abs.(d) for d in ds]
        @test abs.(A)::mat_type == mat_type(fds...)
        @testset "Multiplication with strided matrix/vector" begin
            @test (x = fill(1.,n); A*x ≈ Array(A)*x)
            @test (X = fill(1.,n,2); A*X ≈ Array(A)*X)
        end
        @testset "Binary operations" begin
            B = mat_type == Tridiagonal ? mat_type(a, b, c) : mat_type(b, a)
            fB = map(elty <: Complex ? ComplexF64 : Float64, Array(B))
            for op in (+, -, *)
                @test Array(op(A, B)) ≈ op(fA, fB)
            end
            α = rand(elty)
            @test Array(α*A) ≈ α*Array(A)
            @test Array(A*α) ≈ Array(A)*α
            @test Array(A/α) ≈ Array(A)/α

            @testset "Matmul with Triangular types" begin
                @test A*LinearAlgebra.UnitUpperTriangular(Matrix(1.0I, n, n)) ≈ fA
                @test A*LinearAlgebra.UnitLowerTriangular(Matrix(1.0I, n, n)) ≈ fA
                @test A*UpperTriangular(Matrix(1.0I, n, n)) ≈ fA
                @test A*LowerTriangular(Matrix(1.0I, n, n)) ≈ fA
            end
            @testset "mul! errors" begin
                Cnn, Cnm, Cmn = Matrix{elty}.(undef, ((n,n), (n,n+1), (n+1,n)))
                @test_throws DimensionMismatch LinearAlgebra.mul!(Cnn,A,Cnm)
                @test_throws DimensionMismatch LinearAlgebra.mul!(Cnn,A,Cmn)
                @test_throws DimensionMismatch LinearAlgebra.mul!(Cnn,B,Cmn)
                @test_throws DimensionMismatch LinearAlgebra.mul!(Cmn,B,Cnn)
                @test_throws DimensionMismatch LinearAlgebra.mul!(Cnm,B,Cnn)
            end
        end
        @testset "Negation" begin
            mA = -A
            @test mA isa mat_type
            @test -mA == A
        end
        if mat_type == SymTridiagonal
            @testset "Tridiagonal/SymTridiagonal mixing ops" begin
                B = convert(Tridiagonal{elty}, A)
                @test B == A
                @test B + A == A + B
                @test B - A == A - B
            end
            if elty <: LinearAlgebra.BlasReal
                @testset "Eigensystems" begin
                    zero, infinity = convert(elty, 0), convert(elty, Inf)
                    @testset "stebz! and stein!" begin
                        w, iblock, isplit = LAPACK.stebz!('V', 'B', -infinity, infinity, 0, 0, zero, b, a)
                        evecs = LAPACK.stein!(b, a, w)

                        (e, v) = eigen(SymTridiagonal(b, a))
                        @test e ≈ w
                        test_approx_eq_vecs(v, evecs)
                    end
                    @testset "stein! call using iblock and isplit" begin
                        w, iblock, isplit = LAPACK.stebz!('V', 'B', -infinity, infinity, 0, 0, zero, b, a)
                        evecs = LAPACK.stein!(b, a, w, iblock, isplit)
                        test_approx_eq_vecs(v, evecs)
                    end
                    @testset "stegr! call with index range" begin
                        F = eigen(SymTridiagonal(b, a),1:2)
                        fF = eigen(Symmetric(Array(SymTridiagonal(b, a))),1:2)
                        test_approx_eq_modphase(F.vectors, fF.vectors)
                        @test F.values ≈ fF.values
                    end
                    @testset "stegr! call with value range" begin
                        F = eigen(SymTridiagonal(b, a),0.0,1.0)
                        fF = eigen(Symmetric(Array(SymTridiagonal(b, a))),0.0,1.0)
                        test_approx_eq_modphase(F.vectors, fF.vectors)
                        @test F.values ≈ fF.values
                    end
                    @testset "eigenvalues/eigenvectors of symmetric tridiagonal" begin
                        if elty === Float32 || elty === Float64
                            DT, VT = @inferred eigen(A)
                            @inferred eigen(A, 2:4)
                            @inferred eigen(A, 1.0, 2.0)
                            D, Vecs = eigen(fA)
                            @test DT ≈ D
                            @test abs.(VT'Vecs) ≈ Matrix(elty(1)I, n, n)
                            test_approx_eq_modphase(eigvecs(A), eigvecs(fA))
                            #call to LAPACK.stein here
                            test_approx_eq_modphase(eigvecs(A,eigvals(A)),eigvecs(A))
                        elseif elty != Int
                            # check that undef is determined accurately even if type inference
                            # bails out due to the number of try/catch blocks in this code.
                            @test_throws UndefVarError fA
                        end
                    end
                end
            end
            if elty <: Real
                Ts = SymTridiagonal(d, dl)
                Fs = Array(Ts)
                Tldlt = factorize(Ts)
                @testset "symmetric tridiagonal" begin
                    @test_throws DimensionMismatch Tldlt\rand(elty,n+1)
                    @test size(Tldlt) == size(Ts)
                    if elty <: AbstractFloat
                        @test LinearAlgebra.LDLt{elty,SymTridiagonal{elty,Vector{elty}}}(Tldlt) === Tldlt
                        @test LinearAlgebra.LDLt{elty}(Tldlt) === Tldlt
                        @test typeof(convert(LinearAlgebra.LDLt{Float32,Matrix{Float32}},Tldlt)) ==
                            LinearAlgebra.LDLt{Float32,Matrix{Float32}}
                        @test typeof(convert(LinearAlgebra.LDLt{Float32},Tldlt)) ==
                            LinearAlgebra.LDLt{Float32,SymTridiagonal{Float32,Vector{Float32}}}
                    end
                    for vv in (copy(v), view(v, 1:n))
                        invFsv = Fs\vv
                        x = Ts\vv
                        @test x ≈ invFsv
                        @test Array(Tldlt) ≈ Fs
                    end

                    @testset "similar" begin
                        @test isa(similar(Ts), SymTridiagonal{elty})
                        @test isa(similar(Ts, Int), SymTridiagonal{Int})
                        @test isa(similar(Ts, (3, 2)), Matrix)
                        @test isa(similar(Ts, Int, (3, 2)), Matrix{Int})
                    end

                    @test first(logabsdet(Tldlt)) ≈ first(logabsdet(Fs))
                    @test last(logabsdet(Tldlt))  ≈ last(logabsdet(Fs))
                    # just test that the det method exists. The numerical value of the
                    # determinant is unreliable
                    det(Tldlt)
                end
            end
        else # mat_type is Tridiagonal
            @testset "tridiagonal linear algebra" begin
                for vv in (copy(v), view(copy(v), 1:n))
                    @test A*vv ≈ fA*vv
                    invFv = fA\vv
                    @test A\vv ≈ invFv
                    Tlu = factorize(A)
                    x = Tlu\vv
                    @test x ≈ invFv
                end
                elty != Int && @test A \ v ≈ ldiv!(copy(A), copy(v))
            end
            F = lu(A)
            L1, U1, p1 = F
            G = lu!(F, 2A)
            L2, U2, p2 = F
            @test L1 ≈ L2
            @test 2U1 ≈ U2
            @test p1 == p2
        end
        @testset "generalized dot" begin
            x = fill(convert(elty, 1), n)
            y = fill(convert(elty, 1), n)
            @test dot(x, A, y) ≈ dot(A'x, y) ≈ dot(x, A*y)
            @test dot([1], SymTridiagonal([1], Int[]), [1]) == 1
            @test dot([1], Tridiagonal(Int[], [1], Int[]), [1]) == 1
            @test dot(Int[], SymTridiagonal(Int[], Int[]), Int[]) === 0
            @test dot(Int[], Tridiagonal(Int[], Int[], Int[]), Int[]) === 0
        end
    end
end

@testset "SymTridiagonal/Tridiagonal block matrix" begin
    M = [1 2; 3 4]
    n = 5
    A = SymTridiagonal(fill(M, n), fill(M, n-1))
    @test @inferred A[1,1] == Symmetric(M)
    @test @inferred A[1,2] == M
    @test @inferred A[2,1] == transpose(M)
    @test @inferred diag(A, 1) == fill(M, n-1)
    @test @inferred diag(A, 0) == fill(Symmetric(M), n)
    @test @inferred diag(A, -1) == fill(transpose(M), n-1)
    @test_broken diag(A, -2) == fill(M, n-2)
    @test_broken diag(A, 2) == fill(M, n-2)
    @test isempty(@inferred diag(A, n+1))
    @test isempty(@inferred diag(A, -n-1))

    A[1,1] = Symmetric(2M)
    @test A[1,1] == Symmetric(2M)
    @test_throws ArgumentError A[1,1] = M

    @test tr(A) == sum(diag(A))
    @test issymmetric(tr(A))

    A = Tridiagonal(fill(M, n-1), fill(M, n), fill(M, n-1))
    @test @inferred A[1,1] == M
    @test @inferred A[1,2] == M
    @test @inferred A[2,1] == M
    @test @inferred diag(A, 1) == fill(M, n-1)
    @test @inferred diag(A, 0) == fill(M, n)
    @test @inferred diag(A, -1) == fill(M, n-1)
    @test_broken diag(A, -2) == fill(M, n-2)
    @test_broken diag(A, 2) == fill(M, n-2)
    @test isempty(@inferred diag(A, n+1))
    @test isempty(@inferred diag(A, -n-1))

    for n in 0:2
        dv, ev = fill(M, n), fill(M, max(n-1,0))
        A = SymTridiagonal(dv, ev)
        @test A == Matrix{eltype(A)}(A)

        A = Tridiagonal(ev, dv, ev)
        @test A == Matrix{eltype(A)}(A)
    end

    M = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
    S = SymTridiagonal(fill(M,4), fill(M,3))
    @test diag(S,2) == fill(zero(M), 2)
    @test diag(S,-2) == fill(zero(M), 2)
    @test isempty(diag(S,4))
    @test isempty(diag(S,-4))
end

@testset "Issue 12068" begin
    @test SymTridiagonal([1, 2], [0])^3 == [1 0; 0 8]
end

@testset "Issue #48505" begin
    @test SymTridiagonal([1,2,3],[4,5.0]) == [1.0 4.0 0.0; 4.0 2.0 5.0; 0.0 5.0 3.0]
    @test Tridiagonal([1, 2], [4, 5, 1], [6.0, 7]) == [4.0 6.0 0.0; 1.0 5.0 7.0; 0.0 2.0 1.0]
end

@testset "convert for SymTridiagonal" begin
    STF32 = SymTridiagonal{Float32}(fill(1f0, 5), fill(1f0, 4))
    @test convert(SymTridiagonal{Float64}, STF32)::SymTridiagonal{Float64} == STF32
    @test convert(AbstractMatrix{Float64}, STF32)::SymTridiagonal{Float64} == STF32
end

@testset "constructors from matrix" begin
    @test SymTridiagonal([1 2 3; 2 5 6; 0 6 9]) == [1 2 0; 2 5 6; 0 6 9]
    @test Tridiagonal([1 2 3; 4 5 6; 7 8 9]) == [1 2 0; 4 5 6; 0 8 9]
end

@testset "constructors with range and other abstract vectors" begin
    @test SymTridiagonal(1:3, 1:2) == [1 1 0; 1 2 2; 0 2 3]
    @test Tridiagonal(4:5, 1:3, 1:2) == [1 1 0; 4 2 2; 0 5 3]
end

@testset "Prevent off-diagonal aliasing in Tridiagonal" begin
    e = ones(4)
    f = e[1:end-1]
    T = Tridiagonal(f, 2e, f)
    T ./= 10
    @test all(==(0.1), f)
end

@testset "Issue #26994 (and the empty case)" begin
    T = SymTridiagonal([1.0],[3.0])
    x = ones(1)
    @test T*x == ones(1)
    @test SymTridiagonal(ones(0), ones(0)) * ones(0, 2) == ones(0, 2)
end

@testset "Issue 29630" begin
    function central_difference_discretization(N; dfunc = x -> 12x^2 - 2N^2,
                                               dufunc = x -> N^2 + 4N*x,
                                               dlfunc = x -> N^2 - 4N*x,
                                               bfunc = x -> 114ℯ^-x * (1 + 3x),
                                               b0 = 0, bf = 57/ℯ,
                                               x0 = 0, xf = 1)
        h = 1/N
        d, du, dl, b = map(dfunc, (x0+h):h:(xf-h)), map(dufunc, (x0+h):h:(xf-2h)),
                       map(dlfunc, (x0+2h):h:(xf-h)), map(bfunc, (x0+h):h:(xf-h))
        b[1] -= dlfunc(x0)*b0     # subtract the boundary term
        b[end] -= dufunc(xf)*bf   # subtract the boundary term
        Tridiagonal(dl, d, du), b
    end

    A90, b90 = central_difference_discretization(90)

    @test A90\b90 ≈ inv(A90)*b90
end

@testset "singular values of SymTridiag" begin
    @test svdvals(SymTridiagonal([-4,2,3], [0,0])) ≈ [4,3,2]
    @test svdvals(SymTridiagonal(collect(0.:10.), zeros(10))) ≈ reverse(0:10)
    @test svdvals(SymTridiagonal([1,2,1], [1,1])) ≈ [3,1,0]
    # test that dependent methods such as `cond` also work
    @test cond(SymTridiagonal([1,2,3], [0,0])) ≈ 3
end

@testset "sum, mapreduce" begin
    T = Tridiagonal([1,2], [1,2,3], [7,8])
    Tdense = Matrix(T)
    S = SymTridiagonal([1,2,3], [1,2])
    Sdense = Matrix(S)
    @test sum(T) == 24
    @test sum(S) == 12
    @test_throws ArgumentError sum(T, dims=0)
    @test sum(T, dims=1) == sum(Tdense, dims=1)
    @test sum(T, dims=2) == sum(Tdense, dims=2)
    @test sum(T, dims=3) == sum(Tdense, dims=3)
    @test typeof(sum(T, dims=1)) == typeof(sum(Tdense, dims=1))
    @test mapreduce(one, min, T, dims=1) == mapreduce(one, min, Tdense, dims=1)
    @test mapreduce(one, min, T, dims=2) == mapreduce(one, min, Tdense, dims=2)
    @test mapreduce(one, min, T, dims=3) == mapreduce(one, min, Tdense, dims=3)
    @test typeof(mapreduce(one, min, T, dims=1)) == typeof(mapreduce(one, min, Tdense, dims=1))
    @test mapreduce(zero, max, T, dims=1) == mapreduce(zero, max, Tdense, dims=1)
    @test mapreduce(zero, max, T, dims=2) == mapreduce(zero, max, Tdense, dims=2)
    @test mapreduce(zero, max, T, dims=3) == mapreduce(zero, max, Tdense, dims=3)
    @test typeof(mapreduce(zero, max, T, dims=1)) == typeof(mapreduce(zero, max, Tdense, dims=1))
    @test_throws ArgumentError sum(S, dims=0)
    @test sum(S, dims=1) == sum(Sdense, dims=1)
    @test sum(S, dims=2) == sum(Sdense, dims=2)
    @test sum(S, dims=3) == sum(Sdense, dims=3)
    @test typeof(sum(S, dims=1)) == typeof(sum(Sdense, dims=1))
    @test mapreduce(one, min, S, dims=1) == mapreduce(one, min, Sdense, dims=1)
    @test mapreduce(one, min, S, dims=2) == mapreduce(one, min, Sdense, dims=2)
    @test mapreduce(one, min, S, dims=3) == mapreduce(one, min, Sdense, dims=3)
    @test typeof(mapreduce(one, min, S, dims=1)) == typeof(mapreduce(one, min, Sdense, dims=1))
    @test mapreduce(zero, max, S, dims=1) == mapreduce(zero, max, Sdense, dims=1)
    @test mapreduce(zero, max, S, dims=2) == mapreduce(zero, max, Sdense, dims=2)
    @test mapreduce(zero, max, S, dims=3) == mapreduce(zero, max, Sdense, dims=3)
    @test typeof(mapreduce(zero, max, S, dims=1)) == typeof(mapreduce(zero, max, Sdense, dims=1))

    T = Tridiagonal(Int[], Int[], Int[])
    Tdense = Matrix(T)
    S = SymTridiagonal(Int[], Int[])
    Sdense = Matrix(S)
    @test sum(T) == 0
    @test sum(S) == 0
    @test_throws ArgumentError sum(T, dims=0)
    @test sum(T, dims=1) == sum(Tdense, dims=1)
    @test sum(T, dims=2) == sum(Tdense, dims=2)
    @test sum(T, dims=3) == sum(Tdense, dims=3)
    @test typeof(sum(T, dims=1)) == typeof(sum(Tdense, dims=1))
    @test_throws ArgumentError sum(S, dims=0)
    @test sum(S, dims=1) == sum(Sdense, dims=1)
    @test sum(S, dims=2) == sum(Sdense, dims=2)
    @test sum(S, dims=3) == sum(Sdense, dims=3)
    @test typeof(sum(S, dims=1)) == typeof(sum(Sdense, dims=1))

    T = Tridiagonal(Int[], Int[2], Int[])
    Tdense = Matrix(T)
    S = SymTridiagonal(Int[2], Int[])
    Sdense = Matrix(S)
    @test sum(T) == 2
    @test sum(S) == 2
    @test_throws ArgumentError sum(T, dims=0)
    @test sum(T, dims=1) == sum(Tdense, dims=1)
    @test sum(T, dims=2) == sum(Tdense, dims=2)
    @test sum(T, dims=3) == sum(Tdense, dims=3)
    @test typeof(sum(T, dims=1)) == typeof(sum(Tdense, dims=1))
    @test_throws ArgumentError sum(S, dims=0)
    @test sum(S, dims=1) == sum(Sdense, dims=1)
    @test sum(S, dims=2) == sum(Sdense, dims=2)
    @test sum(S, dims=3) == sum(Sdense, dims=3)
    @test typeof(sum(S, dims=1)) == typeof(sum(Sdense, dims=1))
end

@testset "Issue #28994 (sum of Tridigonal and UniformScaling)" begin
    dl = [1., 1.]
    d = [-2., -2., -2.]
    T = Tridiagonal(dl, d, dl)
    S = SymTridiagonal(T)

    @test diag(T + 2I) == zero(d)
    @test diag(S + 2I) == zero(d)
end

@testset "convert Tridiagonal to SymTridiagonal error" begin
    du = rand(Float64, 4)
    d  = rand(Float64, 5)
    dl = rand(Float64, 4)
    T = Tridiagonal(dl, d, du)
    @test_throws ArgumentError SymTridiagonal{Float32}(T)
end

# Issue #38765
@testset "Eigendecomposition with different lengths" begin
    # length(A.ev) can be either length(A.dv) or length(A.dv) - 1
    A = SymTridiagonal(fill(1.0, 3), fill(-1.0, 3))
    F = eigen(A)
    A2 = SymTridiagonal(fill(1.0, 3), fill(-1.0, 2))
    F2 = eigen(A2)
    test_approx_eq_modphase(F.vectors, F2.vectors)
    @test F.values ≈ F2.values ≈ eigvals(A) ≈ eigvals(A2)
    @test eigvecs(A) ≈ eigvecs(A2)
    @test eigvecs(A, eigvals(A)[1:1]) ≈ eigvecs(A2, eigvals(A2)[1:1])
end

@testset "non-commutative algebra (#39701)" begin
    for A in (SymTridiagonal(Quaternion.(randn(5), randn(5), randn(5), randn(5)), Quaternion.(randn(4), randn(4), randn(4), randn(4))),
              Tridiagonal(Quaternion.(randn(4), randn(4), randn(4), randn(4)), Quaternion.(randn(5), randn(5), randn(5), randn(5)), Quaternion.(randn(4), randn(4), randn(4), randn(4))))
        c = Quaternion(1,2,3,4)
        @test A * c ≈ Matrix(A) * c
        @test A / c ≈ Matrix(A) / c
        @test c * A ≈ c * Matrix(A)
        @test c \ A ≈ c \ Matrix(A)
    end
end

@testset "adjoint of LDLt" begin
    Sr = SymTridiagonal(randn(5), randn(4))
    Sc = SymTridiagonal(complex.(randn(5)) .+ 1im, complex.(randn(4), randn(4)))
    b = ones(size(Sr, 1))

    F = ldlt(Sr)
    @test F\b == F'\b

    F = ldlt(Sc)
    @test copy(Sc')\b == F'\b
end

@testset "symmetric and hermitian tridiagonals" begin
    A = [im 0; 0 -im]
    @test issymmetric(A)
    @test !ishermitian(A)

    # real
    A = SymTridiagonal(randn(5), randn(4))
    @test issymmetric(A)
    @test ishermitian(A)

    A = Tridiagonal(A.ev, A.dv, A.ev .+ 1)
    @test !issymmetric(A)
    @test !ishermitian(A)

    # complex
    # https://github.com/JuliaLang/julia/pull/41037#discussion_r645524081
    S = SymTridiagonal(randn(5) .+ 0im, randn(5) .+ 0im)
    S.ev[end] = im
    @test issymmetric(S)
    @test ishermitian(S)

    S = SymTridiagonal(randn(5) .+ 1im, randn(4) .+ 1im)
    @test issymmetric(S)
    @test !ishermitian(S)

    S = Tridiagonal(S.ev, S.dv, adjoint.(S.ev))
    @test !issymmetric(S)
    @test !ishermitian(S)

    S = Tridiagonal(S.dl, real.(S.d) .+ 0im, S.du)
    @test !issymmetric(S)
    @test ishermitian(S)
end

@testset "Conversion to AbstractArray" begin
    # tests corresponding to #34995
    v1 = ImmutableArray([1, 2])
    v2 = ImmutableArray([3, 4, 5])
    v3 = ImmutableArray([6, 7])
    T = Tridiagonal(v1, v2, v3)
    Tsym = SymTridiagonal(v2, v1)

    @test convert(AbstractArray{Float64}, T)::Tridiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == T
    @test convert(AbstractMatrix{Float64}, T)::Tridiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == T
    @test convert(AbstractArray{Float64}, Tsym)::SymTridiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Tsym
    @test convert(AbstractMatrix{Float64}, Tsym)::SymTridiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Tsym
end

@testset "dot(x,A,y) for A::Tridiagonal or SymTridiagonal" begin
    for elty in (Float32, Float64, ComplexF32, ComplexF64, Int)
        x = fill(convert(elty, 1), 0)
        T = Tridiagonal(x, x, x)
        Tsym = SymTridiagonal(x, x)
        @test dot(x, T, x) == 0.0
        @test dot(x, Tsym, x) == 0.0
    end
end

@testset "non-number eltype" begin
    @testset "sum for SymTridiagonal" begin
        dv = [SizedArray{(2,2)}(rand(1:2048,2,2)) for i in 1:10]
        ev = [SizedArray{(2,2)}(rand(1:2048,2,2)) for i in 1:10]
        S = SymTridiagonal(dv, ev)
        Sdense = Matrix(S)
        @test Sdense == collect(S)
        @test sum(S) == sum(Sdense)
        @test sum(S, dims = 1) == sum(Sdense, dims = 1)
        @test sum(S, dims = 2) == sum(Sdense, dims = 2)
    end
    @testset "issymmetric/ishermitian for Tridiagonal" begin
        @test !issymmetric(Tridiagonal([[1 2;3 4]], [[1 2;2 3], [1 2;2 3]], [[1 2;3 4]]))
        @test !issymmetric(Tridiagonal([[1 3;2 4]], [[1 2;3 4], [1 2;3 4]], [[1 2;3 4]]))
        @test issymmetric(Tridiagonal([[1 3;2 4]], [[1 2;2 3], [1 2;2 3]], [[1 2;3 4]]))

        @test ishermitian(Tridiagonal([[1 3;2 4].+im], [[1 2;2 3].+0im, [1 2;2 3].+0im], [[1 2;3 4].-im]))
        @test !ishermitian(Tridiagonal([[1 3;2 4].+im], [[1 2;2 3].+0im, [1 2;2 3].+0im], [[1 2;3 4].+im]))
        @test !ishermitian(Tridiagonal([[1 3;2 4].+im], [[1 2;2 3].+im, [1 2;2 3].+0im], [[1 2;3 4].-im]))
    end
    @testset "== between Tridiagonal and SymTridiagonal" begin
        dv = [SizedArray{(2,2)}([1 2;3 4]) for i in 1:4]
        ev = [SizedArray{(2,2)}([3 4;1 2]) for i in 1:4]
        S = SymTridiagonal(dv, ev)
        Sdense = Matrix(S)
        @test S == Tridiagonal(diag(Sdense, -1), diag(Sdense),  diag(Sdense, 1)) == S
        @test S !== Tridiagonal(diag(Sdense, 1), diag(Sdense),  diag(Sdense, 1)) !== S
    end
end

@testset "copyto! between SymTridiagonal and Tridiagonal" begin
    ev, dv = [1:4;], [1:5;]
    S = SymTridiagonal(dv, ev)
    T = Tridiagonal(zero(ev), zero(dv), zero(ev))
    @test copyto!(T, S) == S
    @test copyto!(zero(S), T) == T

    ev2 = [1:5;]
    S = SymTridiagonal(dv, ev2)
    T = Tridiagonal(zeros(length(ev2)-1), zero(dv), zeros(length(ev2)-1))
    @test copyto!(T, S) == S
    @test copyto!(zero(S), T) == T

    T2 = Tridiagonal(ones(length(ev)), zero(dv), zero(ev))
    @test_throws "cannot copy an asymmetric Tridiagonal matrix to a SymTridiagonal" copyto!(zero(S), T2)

    @testset "mismatched sizes" begin
        dv2 = [4; @view dv[2:end]]
        @test copyto!(S, SymTridiagonal([4], Int[])) == SymTridiagonal(dv2, ev)
        @test copyto!(T, SymTridiagonal([4], Int[])) == Tridiagonal(ev, dv2, ev)
        @test copyto!(S, Tridiagonal(Int[], [4], Int[])) == SymTridiagonal(dv2, ev)
        @test copyto!(T, Tridiagonal(Int[], [4], Int[])) == Tridiagonal(ev, dv2, ev)
        @test copyto!(S, SymTridiagonal(Int[], Int[])) == SymTridiagonal(dv, ev)
        @test copyto!(T, SymTridiagonal(Int[], Int[])) == Tridiagonal(ev, dv, ev)
        @test copyto!(S, Tridiagonal(Int[], Int[], Int[])) == SymTridiagonal(dv, ev)
        @test copyto!(T, Tridiagonal(Int[], Int[], Int[])) == Tridiagonal(ev, dv, ev)
    end
end

@testset "copyto! with UniformScaling" begin
    @testset "Tridiagonal" begin
        @testset "Fill" begin
            for len in (4, InfiniteArrays.Infinity())
                d = FillArrays.Fill(1, len)
                ud = FillArrays.Fill(0, len-1)
                T = Tridiagonal(ud, d, ud)
                @test copyto!(T, I) === T
            end
        end
        T = Tridiagonal(fill(3, 3), fill(2, 4), fill(3, 3))
        copyto!(T, I)
        @test all(isone, diag(T))
        @test all(iszero, diag(T, 1))
        @test all(iszero, diag(T, -1))
    end
    @testset "SymTridiagonal" begin
        @testset "Fill" begin
            for len in (4, InfiniteArrays.Infinity())
                d = FillArrays.Fill(1, len)
                ud = FillArrays.Fill(0, len-1)
                ST = SymTridiagonal(d, ud)
                @test copyto!(ST, I) === ST
            end
        end
        ST = SymTridiagonal(fill(2, 4), fill(3, 3))
        copyto!(ST, I)
        @test all(isone, diag(ST))
        @test all(iszero, diag(ST, 1))
        @test all(iszero, diag(ST, -1))
    end
end

@testset "custom axes" begin
    dv, uv = OffsetArray(1:4), OffsetArray(1:3)
    B = Tridiagonal(uv, dv, uv)
    ax = axes(dv, 1)
    @test axes(B) === (ax, ax)
    B = SymTridiagonal(dv, uv)
    @test axes(B) === (ax, ax)
end

@testset "Reverse operation on Tridiagonal" begin
    for n in 5:6
        d = randn(n)
        dl = randn(n - 1)
        du = randn(n - 1)
        T = Tridiagonal(dl, d, du)
        @test reverse(T, dims=1) == reverse(Matrix(T), dims=1)
        @test reverse(T, dims=2) == reverse(Matrix(T), dims=2)
        @test reverse(T)::Tridiagonal == reverse(Matrix(T)) == reverse!(copy(T))
    end
end

@testset "Reverse operation on SymTridiagonal" begin
    n = 5
    d = randn(n)
    dl = randn(n - 1)
    ST = SymTridiagonal(d, dl)
    @test reverse(ST, dims=1) == reverse(Matrix(ST), dims=1)
    @test reverse(ST, dims=2) == reverse(Matrix(ST), dims=2)
    @test reverse(ST)::SymTridiagonal == reverse(Matrix(ST))
end

@testset "getindex with Integers" begin
    dv, ev = 1:4, 1:3
    for S in (Tridiagonal(ev, dv, ev), SymTridiagonal(dv, ev))
        @test_throws "invalid index" S[3, true]
        @test S[1,2] == S[Int8(1),UInt16(2)] == S[big(1), Int16(2)]
    end
end

@testset "rmul!/lmul! with banded matrices" begin
    dl, d, du = rand(3), rand(4), rand(3)
    A = Tridiagonal(dl, d, du)
    D = Diagonal(d)
    @test rmul!(copy(A), D) ≈ A * D
    @test lmul!(D, copy(A)) ≈ D * A

    @testset "non-commutative" begin
        S32 = SizedArrays.SizedArray{(3,2)}(rand(3,2))
        S33 = SizedArrays.SizedArray{(3,3)}(rand(3,3))
        S22 = SizedArrays.SizedArray{(2,2)}(rand(2,2))
        T = Tridiagonal(fill(S32,3), fill(S32, 4), fill(S32, 3))
        D = Diagonal(fill(S22, size(T,2)))
        @test rmul!(copy(T), D) ≈ T * D
        D = Diagonal(fill(S33, size(T,1)))
        @test lmul!(D, copy(T)) ≈ D * T
    end
end

@testset "rmul!/lmul! with numbers" begin
    for T in (SymTridiagonal(rand(4), rand(3)), Tridiagonal(rand(3), rand(4), rand(3)))
        @test rmul!(copy(T), 0.2) ≈ rmul!(Array(T), 0.2)
        @test lmul!(0.2, copy(T)) ≈ lmul!(0.2, Array(T))
        @test_throws ArgumentError rmul!(T, NaN)
        @test_throws ArgumentError lmul!(NaN, T)
    end
    for T in (SymTridiagonal(rand(2), rand(1)), Tridiagonal(rand(1), rand(2), rand(1)))
        @test all(isnan, rmul!(copy(T), NaN))
        @test all(isnan, lmul!(NaN, copy(T)))
    end
end

@testset "mul with empty arrays" begin
    A = zeros(5,0)
    T = Tridiagonal(zeros(0), zeros(0), zeros(0))
    TL = Tridiagonal(zeros(4), zeros(5), zeros(4))
    @test size(A * T) == size(A)
    @test size(TL * A) == size(A)
    @test size(T * T) == size(T)
    C = similar(A)
    @test mul!(C, A, T) == A * T
    @test mul!(C, TL, A) == TL * A
    @test mul!(similar(T), T, T) == T * T
    @test mul!(similar(T, size(T)), T, T) == T * T

    v = zeros(size(T,2))
    @test size(T * v) == size(v)
    @test mul!(similar(v), T, v) == T * v

    D = Diagonal(zeros(size(T,2)))
    @test size(T * D) == size(D * T) == size(D)
    @test mul!(similar(D), T, D) == mul!(similar(D), D, T) == T * D
end

@testset "show" begin
    T = Tridiagonal(1:3, 1:4, 1:3)
    @test sprint(show, T) == "Tridiagonal(1:3, 1:4, 1:3)"
    S = SymTridiagonal(1:4, 1:3)
    @test sprint(show, S) == "SymTridiagonal(1:4, 1:3)"

    m = SizedArrays.SizedArray{(2,2)}(reshape([1:4;],2,2))
    T = Tridiagonal(fill(m,2), fill(m,3), fill(m,2))
    @test sprint(show, T) == "Tridiagonal($(repr(diag(T,-1))), $(repr(diag(T))), $(repr(diag(T,1))))"
    S = SymTridiagonal(fill(m,3), fill(m,2))
    @test sprint(show, S) == "SymTridiagonal($(repr(diag(S))), $(repr(diag(S,1))))"
end

@testset "mul for small matrices" begin
    @testset for n in 0:6
        for T in (
                Tridiagonal(rand(max(n-1,0)), rand(n), rand(max(n-1,0))),
                SymTridiagonal(rand(n), rand(max(n-1,0))),
                )
            M = Matrix(T)
            @test T * T ≈ M * M
            @test mul!(similar(T, size(T)), T, T) ≈ M * M
            @test mul!(ones(size(T)), T, T, 2, 4) ≈ M * M * 2 .+ 4

            for m in 0:6
                AR = rand(n,m)
                AL = rand(m,n)
                @test AL * T ≈ AL * M
                @test T * AR ≈ M * AR
                @test mul!(similar(AL), AL, T) ≈ AL * M
                @test mul!(similar(AR), T, AR) ≈ M * AR
                @test mul!(ones(size(AL)), AL, T, 2, 4) ≈ AL * M * 2 .+ 4
                @test mul!(ones(size(AR)), T, AR, 2, 4) ≈ M * AR * 2 .+ 4
            end

            v = rand(n)
            @test T * v ≈ M * v
            @test mul!(similar(v), T, v) ≈ M * v

            D = Diagonal(rand(n))
            @test T * D ≈ M * D
            @test D * T ≈ D * M
            @test mul!(Tridiagonal(similar(T)), D, T) ≈ D * M
            @test mul!(Tridiagonal(similar(T)), T, D) ≈ M * D
            @test mul!(similar(T, size(T)), D, T) ≈ D * M
            @test mul!(similar(T, size(T)), T, D) ≈ M * D
            @test mul!(ones(size(T)), D, T, 2, 4) ≈ D * M * 2 .+ 4
            @test mul!(ones(size(T)), T, D, 2, 4) ≈ M * D * 2 .+ 4

            for uplo in (:U, :L)
                B = Bidiagonal(rand(n), rand(max(0, n-1)), uplo)
                @test T * B ≈ M * B
                @test B * T ≈ B * M
                if n <= 2
                    @test mul!(Tridiagonal(similar(T)), B, T) ≈ B * M
                    @test mul!(Tridiagonal(similar(T)), T, B) ≈ M * B
                end
                @test mul!(similar(T, size(T)), B, T) ≈ B * M
                @test mul!(similar(T, size(T)), T, B) ≈ M * B
                @test mul!(ones(size(T)), B, T, 2, 4) ≈ B * M * 2 .+ 4
                @test mul!(ones(size(T)), T, B, 2, 4) ≈ M * B * 2 .+ 4
            end
        end
    end

    n = 4
    arr = SizedArrays.SizedArray{(2,2)}(reshape([1:4;],2,2))
    for T in (
            SymTridiagonal(fill(arr,n), fill(arr,n-1)),
            Tridiagonal(fill(arr,n-1), fill(arr,n), fill(arr,n-1)),
            )
        @test T * T ≈ Matrix(T) * Matrix(T)
        BL = Bidiagonal(fill(arr,n), fill(arr,n-1), :L)
        BU = Bidiagonal(fill(arr,n), fill(arr,n-1), :U)
        @test BL * T ≈ Matrix(BL) * Matrix(T)
        @test BU * T ≈ Matrix(BU) * Matrix(T)
        @test T * BL ≈ Matrix(T) * Matrix(BL)
        @test T * BU ≈ Matrix(T) * Matrix(BU)
        D = Diagonal(fill(arr,n))
        @test D * T ≈ Matrix(D) * Matrix(T)
        @test T * D ≈ Matrix(T) * Matrix(D)
    end
end

@testset "diagview" begin
    A = Tridiagonal(rand(3), rand(4), rand(3))
    for k in -5:5
        @test diagview(A,k) == diag(A,k)
    end
    v = diagview(A,1)
    v .= 0
    @test all(iszero, diag(A,1))
end

@testset "opnorms" begin
    for T in (Tridiagonal([1,2,3], [1,-2,3,-4], [1,2,3]),
                SymTridiagonal([1,-2,3,-4], [1,2,3]),
                Tridiagonal(Int[], [-5], Int[]),
                SymTridiagonal([-5], Int[]),
                Tridiagonal([1], [1,-2], [3]),
                SymTridiagonal([1,-2], [3])
             )
        @test opnorm(T, 1) == opnorm(Matrix(T), 1)
        @test_skip opnorm(T, 2) ≈ opnorm(Matrix(T), 2) # currently missing
        @test opnorm(T, Inf) == opnorm(Matrix(T), Inf)
    end
end

@testset "block-bidiagonal matrix indexing" begin
    dv = [ones(4,3), ones(2,2).*2, ones(2,3).*3, ones(4,4).*4]
    evu = [ones(4,2), ones(2,3).*2, ones(2,4).*3]
    evl = [ones(2,3), ones(2,2).*2, ones(4,3).*3]
    T = Tridiagonal(evl, dv, evu)
    # check that all the matrices along a column have the same number of columns,
    # and the matrices along a row have the same number of rows
    for j in axes(T, 2), i in 2:size(T, 1)
        @test size(T[i,j], 2) == size(T[1,j], 2)
        @test size(T[i,j], 1) == size(T[i,1], 1)
        if j < i-1 || j > i + 1
            @test iszero(T[i,j])
        end
    end

    @testset "non-standard axes" begin
        s = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
        T = Tridiagonal(fill(s,3), fill(s,4), fill(s,3))
        @test @inferred(T[3,1]) isa typeof(s)
        @test all(iszero, T[3,1])
    end

    # SymTridiagonal requires square diagonal blocks
    dv = [fill(i, i, i) for i in 1:3]
    ev = [ones(Int,1,2), ones(Int,2,3)]
    S = SymTridiagonal(dv, ev)
    @test S == Array{Matrix{Int}}(S)
end

@testset "convert to Tridiagonal/SymTridiagonal" begin
    @testset "Tridiagonal" begin
        for M in [diagm(0 => [1,2,3], 1=>[4,5]),
                diagm(0 => [1,2,3], 1=>[4,5], -1=>[6,7]),
                diagm(-1 => [1,2], 1=>[4,5])]
            B = convert(Tridiagonal, M)
            @test B == Tridiagonal(M)
            B = convert(Tridiagonal{Int8}, M)
            @test B == M
            @test B isa Tridiagonal{Int8}
            B = convert(Tridiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}, M)
            @test B == M
            @test B isa Tridiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}
        end
        @test_throws InexactError convert(Tridiagonal, fill(5, 4, 4))
    end
    @testset "SymTridiagonal" begin
        for M in [diagm(0 => [1,2,3], 1=>[4,5], -1=>[4,5]),
            diagm(0 => [1,2,3]),
            diagm(-1 => [1,2], 1=>[1,2])]
            B = convert(SymTridiagonal, M)
            @test B == SymTridiagonal(M)
            B = convert(SymTridiagonal{Int8}, M)
            @test B == M
            @test B isa SymTridiagonal{Int8}
            B = convert(SymTridiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}, M)
            @test B == M
            @test B isa SymTridiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}
        end
        @test_throws InexactError convert(SymTridiagonal, fill(5, 4, 4))
        @test_throws InexactError convert(SymTridiagonal, diagm(0=>fill(NaN,4)))
    end
end

@testset "isreal" begin
    for M in (SymTridiagonal(ones(2), ones(1)),
            Tridiagonal(ones(2), ones(3), ones(2)))
        @test @inferred((M -> Val(isreal(M)))(M)) == Val(true)
        M = complex.(M)
        @test isreal(M)
        @test !isreal(im*M)
    end
end

@testset "SymTridiagonal from Symmetric" begin
    S = Symmetric(reshape(1:9, 3, 3))
    @testset "helper functions" begin
        @test LinearAlgebra._issymmetric(S)
        @test !LinearAlgebra._issymmetric(Array(S))
    end
    ST = SymTridiagonal(S)
    @test ST == SymTridiagonal(diag(S), diag(S,1))
    S = Symmetric(Tridiagonal(1:3, 1:4, 1:3))
    @test convert(SymTridiagonal, S) == S
end

@testset "setindex! with BandIndex" begin
    T = Tridiagonal(zeros(3), zeros(4), zeros(3))
    T[LinearAlgebra.BandIndex(0,2)] = 1
    @test T[2,2] == 1
    T[LinearAlgebra.BandIndex(1,2)] = 2
    @test T[2,3] == 2
    T[LinearAlgebra.BandIndex(-1,2)] = 3
    @test T[3,2] == 3

    @test_throws "cannot set entry $((1,3)) off the tridiagonal band" T[LinearAlgebra.BandIndex(2,1)] = 1
    @test_throws "cannot set entry $((3,1)) off the tridiagonal band" T[LinearAlgebra.BandIndex(-2,1)] = 1
    @test_throws BoundsError T[LinearAlgebra.BandIndex(size(T,1),1)]
    @test_throws BoundsError T[LinearAlgebra.BandIndex(0,size(T,1)+1)]

    S = SymTridiagonal(zeros(4), zeros(3))
    S[LinearAlgebra.BandIndex(0,2)] = 1
    @test S[2,2] == 1

    @test_throws "cannot set off-diagonal entry $((1,3))" S[LinearAlgebra.BandIndex(2,1)] = 1
    @test_throws BoundsError S[LinearAlgebra.BandIndex(size(S,1),1)]
    @test_throws BoundsError S[LinearAlgebra.BandIndex(0,size(S,1)+1)]
end

@testset "fillband!" begin
    @testset "Tridiagonal" begin
        T = Tridiagonal(zeros(3), zeros(4), zeros(3))
        LinearAlgebra.fillband!(T, 2, 1, 1)
        @test all(==(2), diagview(T,1))
        @test all(==(0), diagview(T,0))
        @test all(==(0), diagview(T,-1))
        LinearAlgebra.fillband!(T, 3, 0, 0)
        @test all(==(3), diagview(T,0))
        @test all(==(2), diagview(T,1))
        @test all(==(0), diagview(T,-1))
        LinearAlgebra.fillband!(T, 4, -1, 1)
        @test all(==(4), diagview(T,-1))
        @test all(==(4), diagview(T,0))
        @test all(==(4), diagview(T,1))
        @test_throws ArgumentError LinearAlgebra.fillband!(T, 3, -2, 2)

        LinearAlgebra.fillstored!(T, 1)
        LinearAlgebra.fillband!(T, 0, -3, 3)
        @test iszero(T)
        LinearAlgebra.fillstored!(T, 1)
        LinearAlgebra.fillband!(T, 0, -10, 10)
        @test iszero(T)

        LinearAlgebra.fillstored!(T, 1)
        T2 = copy(T)
        LinearAlgebra.fillband!(T, 0, -1, -3)
        @test T == T2
        LinearAlgebra.fillband!(T, 0, 10, 10)
        @test T == T2
    end
    @testset "SymTridiagonal" begin
        S = SymTridiagonal(zeros(4), zeros(3))
        @test_throws ArgumentError LinearAlgebra.fillband!(S, 2, -1, -1)
        @test_throws ArgumentError LinearAlgebra.fillband!(S, 2, -2, 2)

        LinearAlgebra.fillband!(S, 1, -1, 1)
        @test all(==(1), diagview(S,-1))
        @test all(==(1), diagview(S,0))
        @test all(==(1), diagview(S,1))

        LinearAlgebra.fillstored!(S, 1)
        LinearAlgebra.fillband!(S, 0, -3, 3)
        @test iszero(S)
        LinearAlgebra.fillstored!(S, 1)
        LinearAlgebra.fillband!(S, 0, -10, 10)
        @test iszero(S)

        LinearAlgebra.fillstored!(S, 1)
        S2 = copy(S)
        LinearAlgebra.fillband!(S, 0, -1, -3)
        @test S == S2
        LinearAlgebra.fillband!(S, 0, 10, 10)
        @test S == S2
    end
end

@testset "issue #1491" begin
    # example where dstemr fails and we need to fall back to dstebz + dstein
    mat = SymTridiagonal([1.9081958235744878e-17, 2.4936649967166602e-17, 0.07337398524939442, -0.037808419232506454, -0.019817057861955256, 0.02058130717549595, 0.0954455291925313, -0.14462978163628976, -0.03603900926534216, 0.10574604437162191, 0.0032290029813105, -0.060091988785585734, 0.15705171399550666, -0.14949686924737385, -0.07781372420855004, 0.03553277340671564, -0.016676626898563265, 0.041446049523777104, 0.10027425296079334, -0.10113668506693083, -0.1489473358740452, 0.1254312751044567, 0.34898062637027333, -0.33073386964350204, 0.2590576535805837, -0.25350583802618526, 0.43822601856881893, -0.4307646088223485, 0.779796594155107, -0.8224956491946511, 0.9089301791662714, -0.8619090557448831, 0.2318453071607881, -0.23320912386320142, 1.3140530211788515, -1.3139256900796155, -0.08390903751033216, 0.08392205661902721, 0.22240033452483454, -0.22246103504490297, 0.029954811490418787, -0.027812830389636778, 0.5263757772376259, -0.5365156423835895, 0.39737424822653156, -0.37936748640222284, 0.9333689384817099, -0.9703849367869426, 1.4044213669711443, -1.3377550567558987, 0.13172978466921414, -0.13519350037871353, 0.4999843222780838, -0.5567430191155477, 1.581258524837276, -1.5341564150622489, 0.8586746358146186, -0.8739294540224967, -0.17581567809216825, 0.15338653603265703, 0.780052785514151, -0.7342698870670235, -0.3173535579981542, 1.3364975241689332, -2.007075402840063, 1.0579643498230522, 0.31139215377118973, -0.39386762936936015, 1.8984071979199968, -1.920544791889951, 0.49358795540233863, -0.4613352087464951, 2.455277289025621, -2.4875372089056578, 0.9570609479497718, -0.9570638796447994, -1.1779359208299722, 1.163074811976796, -1.9301179244054274, 1.9450298141233324, -1.1378646951262403, 0.6350150906475753, 0.08802515616268328, 0.44163280306776015, -0.32773792909920113, 0.32952024534152274, 1.4740728008471815, -1.5038699938845244, 0.6847048634935279, -0.6874625497241057, 1.8751566802154134, -1.8547286700155103, 1.2640098482403956, -1.2174383863708917], [4.0208830900619335, 2.571546145080502e-15, 3.6768867672571894, 2.1275700268385056, 2.6591525597707077, 1.5630051292577454, 3.169596676339297, 1.4476618215762946, 3.1868304175740785, 1.424991114261858, 3.4173704890401693, 1.081122203108007, 3.8388891688820155, 0.7188503200238101, 3.628959037595643, 1.2527365479664625, 3.655577511598652, 0.5795107784623137, 3.9154918603711404, 0.7305602983115349, 4.204984775263665, 0.4667273034707952, 4.086056779436833, 0.7605098745710153, 3.765995021766351, 0.684123305522116, 4.180477278972855, 0.22109649225044853, 4.1892495574511095, 0.5187290599147804, 4.267105890594963, 0.38956213309842774, 4.395771702058166, 0.012968327575103515, 4.370023869671557, 9.368149211179955e-5, 1.7478905204265782, 0.012088046413693718, 2.5400351007261417, 0.02355280591484052, 2.966532136162213, 0.1157714633509318, 3.4610631638791225, 0.21831668966320836, 3.6404645575284627, 0.2911241355881879, 4.444637765978322, 0.46127544435962, 3.9217404428871023, 0.4673958058890463, 3.6185829343499414, 0.6185854819553929, 4.063721954940776, 0.4953177806242663, 4.017151353166678, 0.4035551805802303, 4.318876887238373, 0.17260927157837835, 4.41629817034136, 0.251496875950196, 4.017922104317157, 0.47982184361329844, 3.1370471883575424, 2.4334652498110247, 1.9981994742982454, 0.8095271305424269, 4.40837268034087, 0.18738580231835025, 4.147189662985785, 0.02664936509207341, 2.521497611392047, 0.4752292595225756, 3.795983764313679, 0.004636199088478572, 2.8022242947838416, 0.0004261020805263317, 3.2721199558188827, 0.16661903081014354, 3.1196340292912215, 0.03715487568175004, 2.6328978561492336, 1.5366796979500048, 2.1884275437040994, 0.7509912208995596, 4.056762760327159, 0.3470101947911033, 4.267342895676449, 0.4303372808790385, 4.051346888184243, 0.14621115143865943, 4.2045273832895145, 0.2540956144073191, 4.064260328749705])
    λ = [-4.683759859765017, -4.683759859765017, -4.680622745270949, -4.680066202102797, -4.680066202102788, -4.642387544636312, -4.642387544636312, -4.62850079279404, -4.5645761799521845, -4.564576179952184, -4.563095222112828, -4.563095222112824, -4.4746776824461, -4.474677682446098, -4.372795904675404, -4.372795904675403, -4.32974284978971, -4.329742849789709, -4.190416629535501, -4.160130115632619, -4.160130115632617, -4.158746679334386, -4.020883090061934, -4.0208830900619335, -4.0208830900619335, -3.992984903781014, -3.992984903781012, -3.992659934264744, -3.699861925653249, -3.699861925653247, -3.699861925580689, -3.4509279319955604, -3.4509279319955577, -3.4509279319623327, -3.4082484759685445, -3.4082484759685436, -3.4082484739031704, -2.9611516414045025, -2.9611516414045007, -2.9611516414045003, -2.5495001812022124, -2.549500181202207, -2.5495001812022027, -1.7498663309027835, -1.7498663308995552, -1.749866330899554, -1.749866330899551, 1.7498663308995548, 1.7498663308995575, 1.749866330899561, 1.7498663309055233, 2.549500181202202, 2.5495001812022067, 2.5495001812022093, 2.961151641404497, 2.9611516414045007, 2.9611516414045043, 3.408248467695482, 3.4082484759685476, 3.408248475968552, 3.450927931921516, 3.45092793199556, 3.4509279319955635, 3.6998619253791105, 3.6998619256532437, 3.6998619256532463, 3.99283345851548, 3.9929849037810126, 3.9929849037810143, 4.0208830900619335, 4.0208830900619335, 4.020883090061935, 4.159323562240434, 4.160130115632611, 4.160130115632616, 4.245716444792351, 4.32974284978971, 4.329742849789711, 4.372795904675407, 4.37279590467541, 4.474677682446099, 4.474677682446103, 4.563095222112824, 4.563095222112828, 4.564576179952181, 4.564576179952186, 4.6363803854054515, 4.642387544636308, 4.642387544636311, 4.679778864396294, 4.6800662021027755, 4.680066202102799, 4.6837598597650185, 4.683759859765019]
    densemat = Symmetric(Matrix(mat))
    F = eigen(mat)
    Fdense = eigen(densemat)
    @test λ ≈ eigvals(mat) ≈ eigvals(densemat) ≈ F.values ≈ Fdense.values
    @test F.vectors ≈ Fdense.vectors
end

end # module TestTridiagonal
