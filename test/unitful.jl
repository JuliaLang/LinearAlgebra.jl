module TestUnitfulLinAlg

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random

Random.seed!(1234321)

const TESTDIR = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")
const TESTHELPERS = joinpath(TESTDIR, "testhelpers", "testhelpers.jl")
isdefined(Main, :LinearAlgebraTestHelpers) || Base.include(Main, TESTHELPERS)

using Main.LinearAlgebraTestHelpers.Furlongs

LinearAlgebra.sylvester(a::Furlong,b::Furlong,c::Furlong) = -c / (a + b)

getval(x) = x
getval(x::Furlong) = x.val

# specialized unitful multiplication/solves
@testset "unitful specialized mul/div" begin
    n = 10 #Size of test matrix
    function _bidiagdivmultest(T,
        x,
        typemul=T.uplo == 'U' ? UpperTriangular : Matrix,
        typediv=T.uplo == 'U' ? UpperTriangular : Matrix,
        typediv2=T.uplo == 'U' ? UpperTriangular : Matrix)
        TM = Matrix(T)
        @test map(getval, (T * x)::typemul) ≈ map(getval, TM * x)
        @test map(getval, (x * T)::typemul) ≈ map(getval, x * TM)
        @test map(getval, (x \ T)::typediv) ≈ map(getval, x \ TM)
        @test map(getval, (T / x)::typediv) ≈ map(getval, TM / x)
        if !isa(x, Number)
            @test map(getval, Array((T \ x)::typediv2)) ≈ map(getval, Array(TM \ x))
            @test map(getval, Array((x / T)::typediv2)) ≈ map(getval, Array(x / TM))
        end
        return nothing
    end
    @testset for relty in (Int, Float64, BigFloat), elty in (relty, Complex{relty})
        if relty <: AbstractFloat
            dv = convert(Vector{elty}, randn(n))
            ev = convert(Vector{elty}, randn(n - 1))
            if (elty <: Complex)
                dv += im * convert(Vector{elty}, randn(n))
                ev += im * convert(Vector{elty}, randn(n - 1))
            end
        elseif relty <: Integer
            dv = convert(Vector{elty}, rand(1:10, n))
            ev = convert(Vector{elty}, rand(1:10, n - 1))
            if (elty <: Complex)
                dv += im * convert(Vector{elty}, rand(1:10, n))
                ev += im * convert(Vector{elty}, rand(1:10, n - 1))
            end
        end
        @testset for uplo in (:U, :L)
            T = Bidiagonal(dv, ev, uplo)
            A = Matrix(T)
            for t in (T, Furlong.(T)), (A, dv, ev) in ((A, dv, ev), (Furlong.(A), Furlong.(dv), Furlong.(ev)))
                any(x -> x <: Furlong, (eltype(t), eltype(A))) || continue
                _bidiagdivmultest(t, 5, Bidiagonal, Bidiagonal)
                _bidiagdivmultest(t, 5I, Bidiagonal, Bidiagonal, t.uplo == 'U' ? UpperTriangular : LowerTriangular)
                _bidiagdivmultest(t, Diagonal(dv), Bidiagonal, Bidiagonal, t.uplo == 'U' ? UpperTriangular : LowerTriangular)
                _bidiagdivmultest(t, UpperTriangular(A))
                _bidiagdivmultest(t, UnitUpperTriangular(A))
                _bidiagdivmultest(t, LowerTriangular(A), t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix)
                _bidiagdivmultest(t, UnitLowerTriangular(A), t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix)
                _bidiagdivmultest(t, Bidiagonal(dv, ev, :U), Matrix, Matrix, Matrix)
                _bidiagdivmultest(t, Bidiagonal(dv, ev, :L), Matrix, Matrix, Matrix)
            end
        end
    end
end

# diagonal
@testset for relty in (Float32, Float64, BigFloat), elty in (relty, Complex{relty})
    n = 12
    dd = convert(Vector{elty}, randn(n))
    D = Diagonal(dd)
    @testset "svd/eigen with Diagonal{Furlong}" begin
        Du = Furlong.(D)
        @test Du isa Diagonal{<:Furlong{1}}
        F = svd(Du)
        U, s, V = F
        @test map(getval, Matrix(F)) ≈ map(getval, Du)
        @test svdvals(Du) == s
        @test U isa AbstractMatrix{<:Union{Real,Complex}}
        @test V isa AbstractMatrix{<:Union{Real,Complex}}
        @test s isa AbstractVector{<:Furlong{1}}
        E = eigen(Du; sortby=nothing)
        vals, vecs = E
        @test Matrix(E) == Du
        @test vals isa AbstractVector{<:Furlong{1}}
        @test vecs isa AbstractMatrix{<:Union{Real,Complex}}
        # default `sortby` (`eigsortby` needs a fallback for numbers that are neither `Real` nor `Complex`)
        Es = eigen(Du)
        p = sortperm(dd, by=LinearAlgebra.eigsortby)
        @test Es.values == Furlong.(dd[p])
        @test Es.values isa AbstractVector{<:Furlong{1}}
        @test Es.vectors == Matrix(I, n, n)[:, p]
        @test Du * Es.vectors == Es.vectors * Diagonal(Es.values)
        @test eigvals(Du) == Es.values
    end
end

@testset "pinv with dimensionful matrices" begin
    d = [3.0, 0.0, -2.0]
    A = Matrix(Diagonal(Furlong.(d)))  # a `Matrix`, so that the generic `pinv(::AbstractMatrix)` is exercised
    P = pinv(A)
    @test P isa Matrix{<:Furlong{-1}}
    @test P == Matrix(Diagonal([Furlong{-1}(inv(3.0)), Furlong{-1}(0.0), Furlong{-1}(inv(-2.0))]))
    @test map(getval, A * P * A) ≈ map(getval, A)
    # `rtol` is relative (dimensionless), `atol` defaults to a plain 0
    @test pinv(A; rtol = 0.5) == P
    @test pinv(A; rtol = 0.7) == Matrix(Diagonal([Furlong{-1}(inv(3.0)), Furlong{-1}(0.0), Furlong{-1}(0.0)]))
    @test pinv(Matrix{Furlong{1,Float64}}(undef, 0, 2)) isa Matrix{<:Furlong{-1}}
    @test pinv(Furlong(2.0)) == Furlong{-1}(0.5)
    # dimensionless matrices are unaffected
    @test pinv(Matrix(Diagonal(d))) == Matrix(Diagonal([inv(3.0), 0.0, inv(-2.0)]))
    @test pinv(Matrix(Diagonal([2, 0]))) == [0.5 0.0; 0.0 0.0]
    @test pinv(Matrix(Diagonal([2//1, 0//1]))) == [0.5 0.0; 0.0 0.0]
end

@testset "isapprox for dimensionful arrays" begin
    x = Furlong.([1.0, 2.0])
    # the default `atol = 0` must not be compared with the dimensionful `rtol*max(norm(x), norm(y))`
    @test isapprox(x, Furlong.([1.0, 2.0 + 1e-10]); rtol=1e-8)
    @test !isapprox(x, Furlong.([1.0, 2.0 + 1e-6]); rtol=1e-8)
    @test isapprox(x, x; rtol=0)
    @test !isapprox(x, Furlong.([1.0, 3.0]); rtol=0)
    @test isapprox(Furlong.([0.0 1.0; 2.0 3.0]), Furlong.([0.0 1.0; 2.0 3.0 + 1e-10]); rtol=1e-8)
end

@testset "solves with a dimensionful right-hand side" begin
    A0 = [4.0 1.0 0.0; 1.0 5.0 2.0; 0.0 2.0 6.0]
    b0 = [1.0, 2.0, 3.0]
    b = Furlong.(b0)
    x0 = A0 \ b0
    # a dimensionless matrix and a dimensionful right-hand side: the solution has the units of `b`
    for x in (A0 \ b, lu(A0) \ b, cholesky(A0) \ b)
        @test x isa Vector{<:Furlong{1}}
        @test map(getval, x) ≈ x0
    end
    X = A0 \ [b b]
    @test X isa Matrix{<:Furlong{1}}
    @test map(getval, X) ≈ [x0 x0]
    # the factorization is promoted to the numeric type of the solution (Float32 -> Float64)
    # while keeping its (here trivial) units
    F32 = lu(Float32.(A0))
    x32 = F32 \ b
    @test x32 isa Vector{Furlong{1,Float64}}
    @test map(getval, x32) == LU{Float64}(F32) \ b0
    @test F32 \ Float64.(b0) isa Vector{Float64}
    @test F32 \ Float32.(b0) isa Vector{Float32}
    # numeric promotion of the matrix is unchanged
    @test [4 1; 1 5] \ [1.0, 2.0] isa Vector{Float64}
    @test [4 1; 1 5] \ [1, 2] isa Vector{Float64}
    @test [4.0 1; 1 5] \ [1.0im, 2.0] isa Vector{ComplexF64}
    # a dimensionful matrix still fails inside `ldiv` (the factors and the solution have different units)
    @test_broken Furlong.(A0) \ b isa Vector{<:Furlong{0}}
end

# givens
@testset "testing dimensions with Furlongs #36430" begin
    @test_throws MethodError givens(Furlong(1.0), Furlong(2.0), 1, 2)
end

# hessenberg
@testset "dimensional Hessenberg" begin
    n = 10
    Random.seed!(1234321)

    Areal = randn(n, n) / 2
    @testset "Preserve UpperHessenberg shape (issue #39388)" begin
        H = UpperHessenberg(Furlong.(Areal))
        A = Furlong.(rand(n, n))
        d = Furlong.(rand(n))
        dl = Furlong.(rand(n - 1))
        du = Furlong.(rand(n - 1))
        us = Furlong(1) * I

        @testset "$op" for op = (+, -)
            for x = (us, Diagonal(d), Bidiagonal(d, dl, :U), Bidiagonal(d, dl, :L),
                Tridiagonal(dl, d, du), SymTridiagonal(d, dl),
                UpperTriangular(A), UnitUpperTriangular(A))
                @test op(H, x) == op(Array(H), x)
                @test op(x, H) == op(x, Array(H))
                @test op(H, x) isa UpperHessenberg
                @test op(x, H) isa UpperHessenberg
            end
        end
    end
    H = UpperHessenberg(Furlong.(Areal))
    A = randn(n, n)
    d = randn(n)
    dl = randn(n - 1)
    for A in (A, Furlong.(A))
        @testset "Multiplication/division Furlong" begin
            for x = (5, 5I, Diagonal(d), Bidiagonal(d, dl, :U),
                UpperTriangular(A), UnitUpperTriangular(A))
                @test map(getval, (H * x)::UpperHessenberg) ≈ map(getval, Array(H) * x)
                @test map(getval, (x * H)::UpperHessenberg) ≈ map(getval, x * Array(H))
                @test map(getval, (H / x)::UpperHessenberg) ≈ map(getval, Array(H) / x)
                @test map(getval, (x \ H)::UpperHessenberg) ≈ map(getval, x \ Array(H))
            end
            x = Bidiagonal(d, dl, :L)
            @test H * x == Array(H) * x
            @test x * H == x * Array(H)
            @test H / x == Array(H) / x
            @test x \ H == x \ Array(H)
        end
    end
end

# lu
@testset "lu factorization with dimension type" begin
    n = 4
    A = Matrix(Furlong(1.0) * I, n, n)
    F = lu(A).factors
    @test Diagonal(F) == Diagonal(A)
    # upper triangular part has a unit Furlong{1}
    @test all(x -> typeof(x) == Furlong{1,Float64}, F[i, j] for j = 1:n for i = 1:j)
    # lower triangular part is unitless Furlong{0}
    @test all(x -> typeof(x) == Furlong{0,Float64}, F[i, j] for j = 1:n for i = j+1:n)
end

# special
@testset "zero and one for unitful structured matrices" begin
    # eltype with dimensions
    D0 = Diagonal{Furlong{0,Int64}}([1, 2, 3, 4])
    Bu0 = Bidiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3], 'U')
    Bl0 = Bidiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3], 'L')
    T0 = Tridiagonal{Furlong{0,Int64}}([1, 2, 3], [1, 2, 3, 4], [1, 2, 3])
    S0 = SymTridiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3])
    F2 = Furlongs.Furlong{2}(1)
    D2 = Diagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2)
    Bu2 = Bidiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2, 'U')
    Bl2 = Bidiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2, 'L')
    T2 = Tridiagonal{Furlong{2,Int64}}([1, 2, 3] .* F2, [1, 2, 3, 4] .* F2, [1, 2, 3] .* F2)
    S2 = SymTridiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2)
    mats = Any[D0, Bu0, Bl0, T0, S0, D2, Bu2, Bl2, T2, S2]
    for A in mats
        @test iszero(zero(A))
        @test isone(one(A))
        @test zero(A) == zero(Matrix(A))
        @test one(A) == one(Matrix(A))
        @test eltype(one(A)) == typeof(one(eltype(A)))
    end
end

# triangular
@testset "triangular: dimensional correctness" begin
    A = UpperTriangular([Furlong(1) Furlong(4); Furlong(0) Furlong(1)])
    @test sqrt(A)::UpperTriangular == Furlong{1 // 2}.(UpperTriangular([1 2; 0 1]))
    @test inv(A)::UpperTriangular == Furlong{-1}.(UpperTriangular([1 -4; 0 1]))
    B = UnitUpperTriangular([Furlong(1) Furlong(4); Furlong(0) Furlong(1)])
    @test sqrt(B)::UnitUpperTriangular == Furlong{1 // 2}.(UpperTriangular([1 2; 0 1]))
    @test inv(B)::UnitUpperTriangular == Furlong{-1}.(UpperTriangular([1 -4; 0 1]))
    b = [Furlong(5), Furlong(8)]
    @test (A \ b)::Vector{<:Furlong{0}} == (B \ b)::Vector{<:Furlong{0}} == Furlong{0}.([-27, 8])
    C = LowerTriangular([Furlong(1) Furlong(0); Furlong(4) Furlong(1)])
    @test sqrt(C)::LowerTriangular == Furlong{1 // 2}.(LowerTriangular([1 0; 2 1]))
    @test inv(C)::LowerTriangular == Furlong{-1}.(LowerTriangular([1 0; -4 1]))
    D = UnitLowerTriangular([Furlong(1) Furlong(0); Furlong(4) Furlong(1)])
    @test sqrt(D)::UnitLowerTriangular == Furlong{1 // 2}.(UnitLowerTriangular([1 0; 2 1]))
    @test inv(D)::UnitLowerTriangular == Furlong{-1}.(UnitLowerTriangular([1 0; -4 1]))
    b = [Furlong(5), Furlong(8)]
    @test (C \ b)::Vector{<:Furlong{0}} == (D \ b)::Vector{<:Furlong{0}} == Furlong{0}.([5, -12])
end

@testset "unitful 3-arg *" begin
    for n in (2, 3, 5)
        λ = 5
        A = rand(-10:10, n, n)
        b = rand(-10:10, n)
        λu = Furlong{1}(λ)
        Au = Furlong{1}.(A)
        bu = Furlong{1}.(b)
        @test Furlong{3}.(A*A*λ) == Au*Au*λu
        @test Furlong{3}.(A'*A*λ) == Au'*Au*λu
        @test Furlong{3}.(A*A'*λ) == Au*Au'*λu
        @test Furlong{3}.(A'*A'*λ) == Au'*Au'*λu
        @test Furlong{3}.(A*b*λ) == Au*bu*λu
        @test Furlong{3}.(A'*b*λ) == Au'*bu*λu
    end
end

end # module TestUnitfulLinAlg
