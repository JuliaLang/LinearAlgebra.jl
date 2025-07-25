# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestBidiagonal

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: BlasReal, BlasFloat

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

n = 10 #Size of test matrix
Random.seed!(1)

@testset for relty in (Int, Float32, Float64, BigFloat), elty in (relty, Complex{relty})
    if relty <: AbstractFloat
        dv = convert(Vector{elty}, randn(n))
        ev = convert(Vector{elty}, randn(n-1))
        if (elty <: Complex)
            dv += im*convert(Vector{elty}, randn(n))
            ev += im*convert(Vector{elty}, randn(n-1))
        end
    elseif relty <: Integer
        dv = convert(Vector{elty}, rand(1:10, n))
        ev = convert(Vector{elty}, rand(1:10, n-1))
        if (elty <: Complex)
            dv += im*convert(Vector{elty}, rand(1:10, n))
            ev += im*convert(Vector{elty}, rand(1:10, n-1))
        end
    end
    dv0 = zeros(elty, 0)
    ev0 = zeros(elty, 0)

    @testset "Constructors" begin
        for (x, y) in ((dv0, ev0), (dv, ev), (GenericArray(dv), GenericArray(ev)))
            # from vectors
            ubd = Bidiagonal(x, y, :U)
            lbd = Bidiagonal(x, y, :L)
            @test ubd != lbd || x === dv0
            @test ubd.dv === x
            @test lbd.ev === y
            @test_throws ArgumentError Bidiagonal(x, y, :R)
            @test_throws ArgumentError Bidiagonal(x, y, 'R')
            x == dv0 || @test_throws DimensionMismatch Bidiagonal(x, x, :U)
            @test_throws MethodError Bidiagonal(x, y)
            # from matrix
            @test Bidiagonal(ubd, :U) == Bidiagonal(Matrix(ubd), :U) == ubd
            @test Bidiagonal(lbd, :L) == Bidiagonal(Matrix(lbd), :L) == lbd
            # from its own type
            @test typeof(ubd)(ubd) === ubd
            @test typeof(lbd)(lbd) === lbd
        end
        @test eltype(Bidiagonal{elty}([1,2,3,4], [1.0f0,2.0f0,3.0f0], :U)) == elty
        @test eltype(Bidiagonal([1,2,3,4], [1.0f0,2.0f0,3.0f0], :U)) == Float32 # promotion test
        @test isa(Bidiagonal{elty,Vector{elty}}(GenericArray(dv), ev, :U), Bidiagonal{elty,Vector{elty}})
        @test_throws MethodError Bidiagonal(dv, GenericArray(ev), :U)
        @test_throws MethodError Bidiagonal(GenericArray(dv), ev, :U)
        BI = Bidiagonal([1,2,3,4], [1,2,3], :U)
        @test Bidiagonal(BI) === BI
        @test isa(Bidiagonal{elty}(BI), Bidiagonal{elty})
    end

    @testset "getindex, setindex!, size, and similar" begin
        ubd = Bidiagonal(dv, ev, :U)
        lbd = Bidiagonal(dv, ev, :L)
        # bidiagonal getindex / upper & lower
        @test_throws BoundsError ubd[n + 1, 1]
        @test_throws BoundsError ubd[1, n + 1]
        @test ubd[2, 2] == dv[2]
        # bidiagonal getindex / upper
        @test ubd[2, 3] == ev[2]
        @test iszero(ubd[3, 2])
        # bidiagonal getindex / lower
        @test lbd[3, 2] == ev[2]
        @test iszero(lbd[2, 3])
        # bidiagonal setindex! / upper
        cubd = copy(ubd)
        @test_throws ArgumentError ubd[2, 1] = 1
        @test_throws ArgumentError ubd[3, 1] = 1
        @test (cubd[2, 1] = 0; cubd == ubd)
        @test ((cubd[1, 2] = 10) == 10; cubd[1, 2] == 10)
        # bidiagonal setindex! / lower
        clbd = copy(lbd)
        @test_throws ArgumentError lbd[1, 2] = 1
        @test_throws ArgumentError lbd[1, 3] = 1
        @test (clbd[1, 2] = 0; clbd == lbd)
        @test ((clbd[2, 1] = 10) == 10; clbd[2, 1] == 10)
        # bidiagonal setindex! / upper & lower
        @test_throws BoundsError ubd[n + 1, 1] = 1
        @test_throws BoundsError ubd[1, n + 1] = 1
        @test ((cubd[2, 2] = 10) == 10; cubd[2, 2] == 10)
        # bidiagonal size
        @test_throws BoundsError size(ubd, 0)
        @test size(ubd, 1) == size(ubd, 2) == n
        @test size(ubd, 3) == 1
        # bidiagonal similar
        @test isa(similar(ubd), Bidiagonal{elty})
        @test similar(ubd).uplo == ubd.uplo
        @test isa(similar(ubd, Int), Bidiagonal{Int})
        @test similar(ubd, Int).uplo == ubd.uplo
        @test isa(similar(ubd, (3, 2)), Matrix)
        @test isa(similar(ubd, Int, (3, 2)), Matrix{Int})

        # setindex! when off diagonal is zero bug
        Bu = Bidiagonal(rand(elty, 10), zeros(elty, 9), 'U')
        Bl = Bidiagonal(rand(elty, 10), zeros(elty, 9), 'L')
        @test_throws ArgumentError Bu[5, 4] = 1
        @test_throws ArgumentError Bl[4, 5] = 1

        # setindex should return the destination
        @test setindex!(ubd, 1, 1, 1) === ubd
    end

    @testset "isstored" begin
        ubd = Bidiagonal(dv, ev, :U)
        lbd = Bidiagonal(dv, ev, :L)
        # bidiagonal isstored / upper & lower
        @test_throws BoundsError Base.isstored(ubd, n + 1, 1)
        @test_throws BoundsError Base.isstored(ubd, 1, n + 1)
        @test Base.isstored(ubd, 2, 2)
        # bidiagonal isstored / upper
        @test Base.isstored(ubd, 2, 3)
        @test !Base.isstored(ubd, 3, 2)
        # bidiagonal isstored / lower
        @test Base.isstored(lbd, 3, 2)
        @test !Base.isstored(lbd, 2, 3)
    end

    @testset "show" begin
        BD = Bidiagonal(dv, ev, :U)
        @test sprint(show,BD) == "Bidiagonal($(repr(dv)), $(repr(ev)), :U)"
        BD = Bidiagonal(dv,ev,:L)
        @test sprint(show,BD) == "Bidiagonal($(repr(dv)), $(repr(ev)), :L)"
    end

    @testset for uplo in (:U, :L)
        T = Bidiagonal(dv, ev, uplo)

        @testset "Constructor and basic properties" begin
            @test size(T, 1) == size(T, 2) == n
            @test size(T) == (n, n)
            @test Array(T) == diagm(0 => dv, (uplo === :U ? 1 : -1) => ev)
            @test Bidiagonal(Array(T), uplo) == T
            @test big.(T) == T
            @test Array(abs.(T)) == abs.(diagm(0 => dv, (uplo === :U ? 1 : -1) => ev))
            @test Array(real(T)) == real(diagm(0 => dv, (uplo === :U ? 1 : -1) => ev))
            @test Array(imag(T)) == imag(diagm(0 => dv, (uplo === :U ? 1 : -1) => ev))
        end

        @testset for func in (conj, transpose, adjoint)
            @test func(func(T)) == T
            if func ∈ (transpose, adjoint)
                @test func(func(T)) === T
            end
        end

        @testset "permutedims(::Bidiagonal)" begin
            @test permutedims(permutedims(T)) === T
            @test permutedims(T) == transpose.(transpose(T))
            @test permutedims(T, [1, 2]) === T
            @test permutedims(T, (2, 1)) == permutedims(T)
        end

        @testset "triu and tril" begin
            zerosdv = zeros(elty, length(dv))
            zerosev = zeros(elty, length(ev))
            bidiagcopy(dv, ev, uplo) = Bidiagonal(copy(dv), copy(ev), uplo)

            @test istril(Bidiagonal(dv,ev,:L))
            @test istril(Bidiagonal(dv,ev,:L), 1)
            @test !istril(Bidiagonal(dv,ev,:L), -1)
            @test istril(Bidiagonal(zerosdv,ev,:L), -1)
            @test !istril(Bidiagonal(zerosdv,ev,:L), -2)
            @test istril(Bidiagonal(zerosdv,zerosev,:L), -2)
            @test !istril(Bidiagonal(dv,ev,:U))
            @test istril(Bidiagonal(dv,ev,:U), 1)
            @test !istril(Bidiagonal(dv,ev,:U), -1)
            @test !istril(Bidiagonal(zerosdv,ev,:U), -1)
            @test istril(Bidiagonal(zerosdv,zerosev,:U), -1)
            @test tril!(bidiagcopy(dv,ev,:U),-1) == Bidiagonal(zerosdv,zerosev,:U)
            @test tril!(bidiagcopy(dv,ev,:L),-1) == Bidiagonal(zerosdv,ev,:L)
            @test tril!(bidiagcopy(dv,ev,:U),-2) == Bidiagonal(zerosdv,zerosev,:U)
            @test tril!(bidiagcopy(dv,ev,:L),-2) == Bidiagonal(zerosdv,zerosev,:L)
            @test tril!(bidiagcopy(dv,ev,:U),1)  == Bidiagonal(dv,ev,:U)
            @test tril!(bidiagcopy(dv,ev,:L),1)  == Bidiagonal(dv,ev,:L)
            @test tril!(bidiagcopy(dv,ev,:U))    == Bidiagonal(dv,zerosev,:U)
            @test tril!(bidiagcopy(dv,ev,:L))    == Bidiagonal(dv,ev,:L)
            @test_throws ArgumentError tril!(bidiagcopy(dv, ev, :U), -n - 2)
            @test_throws ArgumentError tril!(bidiagcopy(dv, ev, :U), n)

            @test istriu(Bidiagonal(dv,ev,:U))
            @test istriu(Bidiagonal(dv,ev,:U), -1)
            @test !istriu(Bidiagonal(dv,ev,:U), 1)
            @test istriu(Bidiagonal(zerosdv,ev,:U), 1)
            @test !istriu(Bidiagonal(zerosdv,ev,:U), 2)
            @test istriu(Bidiagonal(zerosdv,zerosev,:U), 2)
            @test !istriu(Bidiagonal(dv,ev,:L))
            @test istriu(Bidiagonal(dv,ev,:L), -1)
            @test !istriu(Bidiagonal(dv,ev,:L), 1)
            @test !istriu(Bidiagonal(zerosdv,ev,:L), 1)
            @test istriu(Bidiagonal(zerosdv,zerosev,:L), 1)
            @test triu!(bidiagcopy(dv,ev,:L),1)  == Bidiagonal(zerosdv,zerosev,:L)
            @test triu!(bidiagcopy(dv,ev,:U),1)  == Bidiagonal(zerosdv,ev,:U)
            @test triu!(bidiagcopy(dv,ev,:U),2)  == Bidiagonal(zerosdv,zerosev,:U)
            @test triu!(bidiagcopy(dv,ev,:L),2)  == Bidiagonal(zerosdv,zerosev,:L)
            @test triu!(bidiagcopy(dv,ev,:U),-1) == Bidiagonal(dv,ev,:U)
            @test triu!(bidiagcopy(dv,ev,:L),-1) == Bidiagonal(dv,ev,:L)
            @test triu!(bidiagcopy(dv,ev,:L))    == Bidiagonal(dv,zerosev,:L)
            @test triu!(bidiagcopy(dv,ev,:U))    == Bidiagonal(dv,ev,:U)
            @test_throws ArgumentError triu!(bidiagcopy(dv, ev, :U), -n)
            @test_throws ArgumentError triu!(bidiagcopy(dv, ev, :U), n + 2)
            @test !isdiag(Bidiagonal(dv,ev,:U))
            @test !isdiag(Bidiagonal(dv,ev,:L))
            @test isdiag(Bidiagonal(dv,zerosev,:U))
            @test isdiag(Bidiagonal(dv,zerosev,:L))
        end

        @testset "iszero and isone" begin
            for uplo in (:U, :L)
                BDzero = Bidiagonal(zeros(elty, 10), zeros(elty, 9), uplo)
                BDone = Bidiagonal(ones(elty, 10), zeros(elty, 9), uplo)
                BDmix = Bidiagonal(zeros(elty, 10), zeros(elty, 9), uplo)
                BDmix[end,end] = one(elty)

                @test iszero(BDzero)
                @test !isone(BDzero)
                @test !iszero(BDone)
                @test isone(BDone)
                @test !iszero(BDmix)
                @test !isone(BDmix)
            end
        end

        @testset "trace" begin
            for uplo in (:U, :L)
                B = Bidiagonal(dv, ev, uplo)
                if relty <: Integer
                    @test tr(B) == tr(Matrix(B))
                else
                    @test tr(B) ≈ tr(Matrix(B)) rtol=2eps(relty)
                end
            end
        end

        Tfull = Array(T)
        @testset "Linear solves" begin
            if relty <: AbstractFloat
                c = convert(Matrix{elty}, randn(n,n))
                b = convert(Matrix{elty}, randn(n, 2))
                if (elty <: Complex)
                    b += im*convert(Matrix{elty}, randn(n, 2))
                end
            elseif relty <: Integer
                c = convert(Matrix{elty}, rand(1:10, n, n))
                b = convert(Matrix{elty}, rand(1:10, n, 2))
                if (elty <: Complex)
                    b += im*convert(Matrix{elty}, rand(1:10, n, 2))
                end
            end
            condT = cond(map(ComplexF64,Tfull))
            promty = typeof((zero(relty)*zero(relty) + zero(relty)*zero(relty))/one(relty))
            if relty != BigFloat
                x = transpose(T)\transpose(c)
                tx = transpose(Tfull) \ transpose(c)
                elty <: AbstractFloat && @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                @test_throws DimensionMismatch transpose(T)\transpose(b)
                x = T'\copy(transpose(c))
                tx = Tfull'\copy(transpose(c))
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                @test_throws DimensionMismatch T'\copy(transpose(b))
                x = T\transpose(c)
                tx = Tfull\transpose(c)
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                @test_throws DimensionMismatch T\transpose(b)
            end
            offsizemat = Matrix{elty}(undef, n+1, 2)
            @test_throws DimensionMismatch T \ offsizemat
            @test_throws DimensionMismatch transpose(T) \ offsizemat
            @test_throws DimensionMismatch T' \ offsizemat

            if elty <: BigFloat
                @test_throws SingularException ldiv!(Bidiagonal(zeros(elty, n), ones(elty, n-1), :U), rand(elty, n))
                @test_throws SingularException ldiv!(Bidiagonal(zeros(elty, n), ones(elty, n-1), :L), rand(elty, n))
            end
            let bb = b, cc = c
                for atype in ("Array", "SubArray")
                    if atype == "Array"
                        b = bb
                        c = cc
                    else
                        b = view(bb, 1:n)
                        c = view(cc, 1:n, 1:2)
                    end
                end
                x = T \ b
                tx = Tfull \ b
                @test_throws DimensionMismatch ldiv!(T, Vector{elty}(undef, n+1))
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = transpose(T) \ b
                tx = transpose(Tfull) \ b
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = copy(transpose(b)) / T
                tx = copy(transpose(b)) / Tfull
                @test_throws DimensionMismatch rdiv!(Matrix{elty}(undef, 1, n+1), T)
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = copy(transpose(b)) / transpose(T)
                tx = copy(transpose(b)) / transpose(Tfull)
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                @testset "Generic Mat-vec ops" begin
                    @test T*b ≈ Tfull*b
                    @test T'*b ≈ Tfull'*b
                    if relty != BigFloat # not supported by pivoted QR
                        @test T/b' ≈ Tfull/b'
                    end
                end
            end
            zdv = Vector{elty}(undef, 0)
            zev = Vector{elty}(undef, 0)
            zA  = Bidiagonal(zdv, zev, :U)
            zb  = Vector{elty}(undef, 0)
            @test ldiv!(zA, zb) === zb
            @testset "linear solves with abstract matrices" begin
                diag = b[:,1]
                D = Diagonal(diag)
                x = T \ D
                tx = Tfull \ D
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = D / T
                tx = D / Tfull
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = transpose(T) \ D
                tx = transpose(Tfull) \ D
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
                x = D / transpose(T)
                tx = D / transpose(Tfull)
                @test norm(x-tx,Inf) <= 4*condT*max(eps()*norm(tx,Inf), eps(promty)*norm(x,Inf))
            end
            @testset "Specialized multiplication/division" begin
                function _bidiagdivmultest(T,
                        x,
                        typemul=T.uplo == 'U' ? UpperTriangular : Matrix,
                        typediv=T.uplo == 'U' ? UpperTriangular : Matrix,
                        typediv2=T.uplo == 'U' ? UpperTriangular : Matrix)
                    TM = Matrix(T)
                    @test (T*x)::typemul ≈ TM*x
                    @test (x*T)::typemul ≈ x*TM
                    @test (x\T)::typediv ≈ x\TM
                    @test (T/x)::typediv ≈ TM/x
                    if !isa(x, Number)
                        U = T.uplo == 'U' ? UpperTriangular : LowerTriangular
                        @test Array((T\x)::typediv2) ≈ Array(U(TM)\x)
                        @test Array((x/T)::typediv2) ≈ Array(x/U(TM))
                    end
                    return nothing
                end
                A = Matrix(T)
                t = T
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

        if elty <: BlasReal
            @testset "$f" for f in (floor, trunc, round, ceil)
                @test (f.(Int, T))::Bidiagonal == Bidiagonal(f.(Int, T.dv), f.(Int, T.ev), T.uplo)
                @test (f.(T))::Bidiagonal == Bidiagonal(f.(T.dv), f.(T.ev), T.uplo)
            end
        end

        @testset "diag" begin
            @test (@inferred diag(T))::typeof(dv) == dv
            @test (@inferred diag(T, uplo === :U ? 1 : -1))::typeof(dv) == ev
            @test (@inferred diag(T,2))::typeof(dv) == zeros(elty, n-2)
            @test isempty(@inferred diag(T, -n - 1))
            @test isempty(@inferred diag(T,  n + 1))
            # test diag with another wrapped vector type
            gdv, gev = GenericArray(dv), GenericArray(ev)
            G = Bidiagonal(gdv, gev, uplo)
            @test (@inferred diag(G))::typeof(gdv) == gdv
            @test (@inferred diag(G, uplo === :U ? 1 : -1))::typeof(gdv) == gev
            @test (@inferred diag(G,2))::typeof(gdv) == GenericArray(zeros(elty, n-2))
        end

        @testset "Eigensystems" begin
            if relty <: AbstractFloat
                d1, v1 = eigen(T)
                d2, v2 = eigen(map(elty<:Complex ? ComplexF64 : Float64,Tfull), sortby=nothing)
                @test (uplo === :U ? d1 : reverse(d1)) ≈ d2
                if elty <: Real
                    test_approx_eq_modphase(v1, uplo === :U ? v2 : v2[:,n:-1:1])
                end
            end
        end

        @testset "Singular systems" begin
            if (elty <: BlasReal)
                @test AbstractArray(svd(T)) ≈ AbstractArray(svd!(copy(Tfull)))
                @test svdvals(Tfull) ≈ svdvals(T)
                u1, d1, v1 = svd(Tfull)
                u2, d2, v2 = svd(T)
                @test d1 ≈ d2
                if elty <: Real
                    test_approx_eq_modphase(u1, u2)
                    test_approx_eq_modphase(copy(v1), copy(v2))
                end
                @test 0 ≈ norm(u2*Diagonal(d2)*v2'-Tfull) atol=n*max(n^2*eps(relty),norm(u1*Diagonal(d1)*v1'-Tfull))
                @inferred svdvals(T)
                @inferred svd(T)
            end
        end

        @testset "Binary operations" begin
            @test -T == Bidiagonal(-T.dv,-T.ev,T.uplo)
            @test convert(elty,-1.0) * T == Bidiagonal(-T.dv,-T.ev,T.uplo)
            @test T / convert(elty,-1.0) == Bidiagonal(-T.dv,-T.ev,T.uplo)
            @test T * convert(elty,-1.0) == Bidiagonal(-T.dv,-T.ev,T.uplo)
            @testset for uplo2 in (:U, :L)
                dv = convert(Vector{elty}, relty <: AbstractFloat ? randn(n) : rand(1:10, n))
                ev = convert(Vector{elty}, relty <: AbstractFloat ? randn(n-1) : rand(1:10, n-1))
                T2 = Bidiagonal(dv, ev, uplo2)
                Tfull2 = Array(T2)
                for op in (+, -, *)
                    @test Array(op(T, T2)) ≈ op(Tfull, Tfull2)
                end
                A = kron(T.dv, T.dv')
                @test T * A ≈ lmul!(T, copy(A))
                @test A * T ≈ rmul!(copy(A), T)
            end
            # test pass-through of mul! for SymTridiagonal*Bidiagonal
            TriSym = SymTridiagonal(T.dv, T.ev)
            @test Array(TriSym*T) ≈ Array(TriSym)*Array(T)
            # test pass-through of mul! for AbstractTriangular*Bidiagonal
            Tri = UpperTriangular(diagm(1 => T.ev))
            Dia = Diagonal(T.dv)
            @test Array(Tri*T) ≈ Array(Tri)*Array(T) ≈ rmul!(copy(Tri), T)
            @test Array(T*Tri) ≈ Array(T)*Array(Tri) ≈ lmul!(T, copy(Tri))
            # test mul! itself for these types
            for AA in (Tri, Dia)
                for f in (identity, transpose, adjoint)
                    C = rand(elty, n, n)
                    D = copy(C) + 2.0 * Array(f(AA) * T)
                    mul!(C, f(AA), T, 2.0, 1.0) ≈ D
                end
            end
            # test mul! for BiTrySym * adjoint/transpose AbstractMat
            for f in (identity, transpose, adjoint)
                C = relty == Int ? rand(float(elty), n, n) : rand(elty, n, n)
                B = rand(elty, n, n)
                D = C + 2.0 * Array(T*f(B))
                @test mul!(C, T, f(B), 2.0, 1.0) ≈ D
                @test lmul!(T, copy(f(B))) ≈ T * f(B)
                @test rmul!(copy(f(B)), T) ≈ f(B) * T
            end

            # Issue #31870
            # Bi/Tri/Sym times Diagonal
            Diag = Diagonal(rand(elty, 10))
            BidiagU = Bidiagonal(rand(elty, 10), rand(elty, 9), 'U')
            BidiagL = Bidiagonal(rand(elty, 10), rand(elty, 9), 'L')
            Tridiag = Tridiagonal(rand(elty, 9), rand(elty, 10), rand(elty, 9))
            SymTri = SymTridiagonal(rand(elty, 10), rand(elty, 9))

            mats = Any[Diag, BidiagU, BidiagL, Tridiag, SymTri]
            for a in mats
                for b in mats
                    @test a*b ≈ Matrix(a)*Matrix(b)
                end
            end

            @test typeof(BidiagU*Diag) <: Bidiagonal
            @test typeof(BidiagL*Diag) <: Bidiagonal
            @test typeof(Tridiag*Diag) <: Tridiagonal
            @test typeof(SymTri*Diag)  <: Tridiagonal

            @test typeof(BidiagU*Diag) <: Bidiagonal
            @test typeof(Diag*BidiagL) <: Bidiagonal
            @test typeof(Diag*Tridiag) <: Tridiagonal
            @test typeof(Diag*SymTri)  <: Tridiagonal
        end

        @test inv(T)*Tfull ≈ Matrix(I, n, n)
        @test factorize(T) === T
    end
    BD = Bidiagonal(dv, ev, :U)
    @test Matrix{ComplexF64}(BD) == BD
end

@testset "Constructors with Char uplo" begin
    @test Bidiagonal(Int8[1,2], [1], 'U') == Bidiagonal(Int8[1,2], [1], :U)
end

# Issue 10742 and similar
let A = Bidiagonal([1,2,3], [0,0], :U)
    @test istril(A)
    @test isdiag(A)
end

# test construct from range
@test Bidiagonal(1:3, 1:2, :U) == [1 1 0; 0 2 2; 0 0 3]

@testset "promote_rule" begin
    A = Bidiagonal(fill(1f0,10),fill(1f0,9),:U)
    B = rand(Float64,10,10)
    C = Tridiagonal(rand(Float64,9),rand(Float64,10),rand(Float64,9))
    @test promote_rule(Matrix{Float64}, Bidiagonal{Float64}) == Matrix{Float64}
    @test promote(B,A) == (B, convert(Matrix{Float64}, A))
    @test promote(B,A) isa Tuple{Matrix{Float64}, Matrix{Float64}}
    @test promote(C,A) == (C,Tridiagonal(zeros(Float64,9),convert(Vector{Float64},A.dv),convert(Vector{Float64},A.ev)))
    @test promote(C,A) isa Tuple{Tridiagonal, Tridiagonal}
end

using LinearAlgebra: fillstored!, UnitLowerTriangular
@testset "fill! and fillstored!" begin
    let # fillstored!
        A = Tridiagonal(randn(2), randn(3), randn(2))
        @test fillstored!(A, 3) == Tridiagonal([3, 3], [3, 3, 3], [3, 3])
        B = Bidiagonal(randn(3), randn(2), :U)
        @test fillstored!(B, 2) == Bidiagonal([2,2,2], [2,2], :U)
        S = SymTridiagonal(randn(3), randn(2))
        @test fillstored!(S, 1) == SymTridiagonal([1,1,1], [1,1])
        Ult = UnitLowerTriangular(randn(3,3))
        @test fillstored!(Ult, 3) == UnitLowerTriangular([1 0 0; 3 1 0; 3 3 1])
    end
    let # fill!(exotic, 0)
        exotic_arrays = Any[Tridiagonal(randn(3), randn(4), randn(3)),
        Bidiagonal(randn(3), randn(2), rand([:U,:L])),
        SymTridiagonal(randn(3), randn(2)),
        Diagonal(randn(5)),
        # LowerTriangular(randn(3,3)), # AbstractTriangular fill! deprecated, see below
        # UpperTriangular(randn(3,3)) # AbstractTriangular fill! deprecated, see below
        ]
        for A in exotic_arrays
            @test iszero(fill!(A, 0))
        end

        # Diagonal fill! is no longer deprecated. See #29780
        # AbstractTriangular fill! was defined as fillstored!,
        # not matching the general behavior of fill!, and so it has been deprecated.
        # In a future dev cycle, this fill! methods should probably be reintroduced
        # with behavior matching that of fill! for other structured matrix types.
        # In the interim, equivalently test fillstored! below
        @test iszero(fillstored!(Diagonal(fill(1, 3)), 0))
        @test iszero(fillstored!(LowerTriangular(fill(1, 3, 3)), 0))
        @test iszero(fillstored!(UpperTriangular(fill(1, 3, 3)), 0))
    end
    let # fill!(small, x)
        val = randn()
        b = Bidiagonal(randn(1,1), :U)
        st = SymTridiagonal(randn(1,1))
        d = Diagonal(rand(1))
        for x in (b, st, d)
            @test Array(fill!(x, val)) == fill!(Array(x), val)
        end
        b = Bidiagonal(randn(2,2), :U)
        st = SymTridiagonal(randn(3), randn(2))
        t = Tridiagonal(randn(3,3))
        d = Diagonal(rand(3))
        for x in (b, t, st, d)
            @test_throws ArgumentError fill!(x, val)
            @test Array(fill!(x, 0)) == fill!(Array(x), 0)
        end
    end
end

@testset "pathological promotion (#24707)" begin
    @test promote_type(Matrix{Int}, Bidiagonal{Tuple{S}} where S<:Integer) <: Matrix
    @test promote_type(Matrix{Tuple{T}} where T<:Integer, Bidiagonal{Tuple{S}} where S<:Integer) <: Matrix
    @test promote_type(Matrix{Tuple{T}} where T<:Integer, Bidiagonal{Int}) <: Matrix
    @test promote_type(Tridiagonal{Int}, Bidiagonal{Tuple{S}} where S<:Integer) <: Tridiagonal
    @test promote_type(Tridiagonal{Tuple{T}} where T<:Integer, Bidiagonal{Tuple{S}} where S<:Integer) <: Tridiagonal
    @test promote_type(Tridiagonal{Tuple{T}} where T<:Integer, Bidiagonal{Int}) <: Tridiagonal
end

@testset "solve with matrix elements" begin
    A = triu(tril(randn(9, 9), 3), -3)
    b = randn(9)
    Alb = Bidiagonal(Any[tril(A[1:3,1:3]), tril(A[4:6,4:6]), tril(A[7:9,7:9])],
                     Any[triu(A[4:6,1:3]), triu(A[7:9,4:6])], 'L')
    Aub = Bidiagonal(Any[triu(A[1:3,1:3]), triu(A[4:6,4:6]), triu(A[7:9,7:9])],
                     Any[tril(A[1:3,4:6]), tril(A[4:6,7:9])], 'U')
    bb = Any[b[1:3], b[4:6], b[7:9]]
    @test vcat((Alb\bb)...) ≈ LowerTriangular(A)\b
    @test vcat((Aub\bb)...) ≈ UpperTriangular(A)\b
    Alb = Bidiagonal([tril(A[1:3,1:3]), tril(A[4:6,4:6]), tril(A[7:9,7:9])],
                     [triu(A[4:6,1:3]), triu(A[7:9,4:6])], 'L')
    Aub = Bidiagonal([triu(A[1:3,1:3]), triu(A[4:6,4:6]), triu(A[7:9,7:9])],
                     [tril(A[1:3,4:6]), tril(A[4:6,7:9])], 'U')
    d = [randn(3,3) for _ in 1:3]
    dl = [randn(3,3) for _ in 1:2]
    B = [randn(3,3) for _ in 1:3, _ in 1:3]
    for W in (UpperTriangular, LowerTriangular), t in (identity, adjoint, transpose)
        @test Matrix(t(Alb) \ W(B)) ≈ t(Alb) \ Matrix(W(B))
        @test Matrix(t(Aub) \ W(B)) ≈ t(Aub) \ Matrix(W(B))
        @test Matrix(W(B) / t(Alb)) ≈ Matrix(W(B)) / t(Alb)
        @test Matrix(W(B) / t(Aub)) ≈ Matrix(W(B)) / t(Aub)
    end
end

@testset "sum, mapreduce" begin
    Bu = Bidiagonal([1,2,3], [1,2], :U)
    Budense = Matrix(Bu)
    Bl = Bidiagonal([1,2,3], [1,2], :L)
    Bldense = Matrix(Bl)
    @test sum(Bu) == 9
    @test sum(Bl) == 9
    @test_throws ArgumentError sum(Bu, dims=0)
    @test sum(Bu, dims=1) == sum(Budense, dims=1)
    @test sum(Bu, dims=2) == sum(Budense, dims=2)
    @test sum(Bu, dims=3) == sum(Budense, dims=3)
    @test typeof(sum(Bu, dims=1)) == typeof(sum(Budense, dims=1))
    @test mapreduce(one, min, Bu, dims=1) == mapreduce(one, min, Budense, dims=1)
    @test mapreduce(one, min, Bu, dims=2) == mapreduce(one, min, Budense, dims=2)
    @test mapreduce(one, min, Bu, dims=3) == mapreduce(one, min, Budense, dims=3)
    @test typeof(mapreduce(one, min, Bu, dims=1)) == typeof(mapreduce(one, min, Budense, dims=1))
    @test mapreduce(zero, max, Bu, dims=1) == mapreduce(zero, max, Budense, dims=1)
    @test mapreduce(zero, max, Bu, dims=2) == mapreduce(zero, max, Budense, dims=2)
    @test mapreduce(zero, max, Bu, dims=3) == mapreduce(zero, max, Budense, dims=3)
    @test typeof(mapreduce(zero, max, Bu, dims=1)) == typeof(mapreduce(zero, max, Budense, dims=1))
    @test_throws ArgumentError sum(Bl, dims=0)
    @test sum(Bl, dims=1) == sum(Bldense, dims=1)
    @test sum(Bl, dims=2) == sum(Bldense, dims=2)
    @test sum(Bl, dims=3) == sum(Bldense, dims=3)
    @test typeof(sum(Bl, dims=1)) == typeof(sum(Bldense, dims=1))
    @test mapreduce(one, min, Bl, dims=1) == mapreduce(one, min, Bldense, dims=1)
    @test mapreduce(one, min, Bl, dims=2) == mapreduce(one, min, Bldense, dims=2)
    @test mapreduce(one, min, Bl, dims=3) == mapreduce(one, min, Bldense, dims=3)
    @test typeof(mapreduce(one, min, Bl, dims=1)) == typeof(mapreduce(one, min, Bldense, dims=1))
    @test mapreduce(zero, max, Bl, dims=1) == mapreduce(zero, max, Bldense, dims=1)
    @test mapreduce(zero, max, Bl, dims=2) == mapreduce(zero, max, Bldense, dims=2)
    @test mapreduce(zero, max, Bl, dims=3) == mapreduce(zero, max, Bldense, dims=3)
    @test typeof(mapreduce(zero, max, Bl, dims=1)) == typeof(mapreduce(zero, max, Bldense, dims=1))

    Bu = Bidiagonal([2], Int[], :U)
    Budense = Matrix(Bu)
    Bl = Bidiagonal([2], Int[], :L)
    Bldense = Matrix(Bl)
    @test sum(Bu) == 2
    @test sum(Bl) == 2
    @test_throws ArgumentError sum(Bu, dims=0)
    @test sum(Bu, dims=1) == sum(Budense, dims=1)
    @test sum(Bu, dims=2) == sum(Budense, dims=2)
    @test sum(Bu, dims=3) == sum(Budense, dims=3)
    @test typeof(sum(Bu, dims=1)) == typeof(sum(Budense, dims=1))
end

@testset "empty sub-diagonal" begin
    # `mul!` must use non-specialized method when sub-diagonal is empty
    A = [1 2 3 4]'
    @test A * Tridiagonal(ones(1, 1)) == A
end

@testset "generalized dot" begin
    for elty in (Float64, ComplexF64), n in (5, 1)
        dv = randn(elty, n)
        ev = randn(elty, n-1)
        x = randn(elty, n)
        y = randn(elty, n)
        for uplo in (:U, :L)
            B = Bidiagonal(dv, ev, uplo)
            @test dot(x, B, y) ≈ dot(B'x, y) ≈ dot(x, B*y) ≈ dot(x, Matrix(B), y)
        end
        dv = Vector{elty}(undef, 0)
        ev = Vector{elty}(undef, 0)
        x = Vector{elty}(undef, 0)
        y = Vector{elty}(undef, 0)
        for uplo in (:U, :L)
            B = Bidiagonal(dv, ev, uplo)
            @test dot(x, B, y) === zero(elty)
        end
    end
end

@testset "multiplication of bidiagonal and triangular matrix" begin
    n = 5
    for eltyB in (Int, ComplexF64)
        if eltyB == Int
            BU = Bidiagonal(rand(1:7, n), rand(1:7, n - 1), :U)
            BL = Bidiagonal(rand(1:7, n), rand(1:7, n - 1), :L)
        else
            BU = Bidiagonal(randn(eltyB, n), randn(eltyB, n - 1), :U)
            BL = Bidiagonal(randn(eltyB, n), randn(eltyB, n - 1), :L)
        end
        for eltyT in (Int, ComplexF64)
            for TriT in (LowerTriangular, UnitLowerTriangular, UpperTriangular, UnitUpperTriangular)
                if eltyT == Int
                    T = TriT(rand(1:7, n, n))
                else
                    T = TriT(randn(eltyT, n, n))
                end
                for B in (BU, BL)
                    MB = Matrix(B)
                    MT = Matrix(T)
                    for transB in (identity, adjoint, transpose), transT in (identity, adjoint, transpose)
                        @test transB(B) * transT(T) ≈ transB(MB) * transT(MT)
                        @test transT(T) * transB(B) ≈ transT(MT) * transB(MB)
                    end
                end
            end
        end
    end
end

struct MyNotANumberType
    n::Float64
end
Base.zero(n::MyNotANumberType)      = MyNotANumberType(zero(Float64))
Base.zero(T::Type{MyNotANumberType}) = MyNotANumberType(zero(Float64))
Base.copy(n::MyNotANumberType)      = MyNotANumberType(copy(n.n))
Base.transpose(n::MyNotANumberType) = n

@testset "transpose for a non-numeric eltype" begin
    @test !(MyNotANumberType(1.0) isa Number)
    a = [MyNotANumberType(1.0), MyNotANumberType(2.0), MyNotANumberType(3.0)]
    b = [MyNotANumberType(5.0), MyNotANumberType(6.0)]
    B = Bidiagonal(a, b, :U)
    tB = transpose(B)
    @test tB == Bidiagonal(a, b, :L)
    @test transpose(copy(tB)) == B
end

@testset "empty bidiagonal matrices" begin
    dv0 = zeros(0)
    ev0 = zeros(0)
    zm = zeros(0, 0)
    ubd = Bidiagonal(dv0, ev0, :U)
    lbd = Bidiagonal(dv0, ev0, :L)
    @test size(ubd) == (0, 0)
    @test_throws BoundsError getindex(ubd, 1, 1)
    @test_throws BoundsError setindex!(ubd, 0.0, 1, 1)
    @test similar(ubd) == ubd
    @test similar(lbd, Int) == zeros(Int, 0, 0)
    @test ubd == zm
    @test lbd == zm
    @test ubd == lbd
    @test ubd * ubd == ubd
    @test lbd + lbd == lbd
    @test lbd' == ubd
    @test ubd' == lbd
    @test triu(ubd, 1) == ubd
    @test triu(lbd, 1) == ubd
    @test tril(ubd, -1) == ubd
    @test tril(lbd, -1) == ubd
    @test_throws ArgumentError triu(ubd)
    @test_throws ArgumentError tril(ubd)
    @test sum(ubd) == 0.0
    @test reduce(+, ubd) == 0.0
    @test reduce(+, ubd, dims=1) == zeros(1, 0)
    @test reduce(+, ubd, dims=2) == zeros(0, 1)
    @test hcat(ubd, ubd) == zm
    @test vcat(ubd, lbd) == zm
    @test hcat(lbd, ones(0, 3)) == ones(0, 3)
    @test fill!(copy(ubd), 1.0) == ubd
    @test map(abs, ubd) == zm
    @test lbd .+ 1 == zm
    @test lbd + ubd isa Bidiagonal
    @test lbd .+ ubd isa Bidiagonal
    @test ubd * 5 == ubd
    @test ubd .* 3 == ubd
end

@testset "non-commutative algebra (#39701)" begin
    A = Bidiagonal(Quaternion.(randn(5), randn(5), randn(5), randn(5)), Quaternion.(randn(4), randn(4), randn(4), randn(4)), :U)
    c = Quaternion(1,2,3,4)
    @test A * c ≈ Matrix(A) * c
    @test A / c ≈ Matrix(A) / c
    @test c * A ≈ c * Matrix(A)
    @test c \ A ≈ c \ Matrix(A)
end

@testset "Conversion to AbstractArray" begin
    # tests corresponding to #34995
    dv = ImmutableArray([1, 2, 3, 4])
    ev = ImmutableArray([7, 8, 9])
    Bu = Bidiagonal(dv, ev, :U)
    Bl = Bidiagonal(dv, ev, :L)

    @test convert(AbstractArray{Float64}, Bu)::Bidiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Bu
    @test convert(AbstractMatrix{Float64}, Bu)::Bidiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Bu
    @test convert(AbstractArray{Float64}, Bl)::Bidiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Bl
    @test convert(AbstractMatrix{Float64}, Bl)::Bidiagonal{Float64,ImmutableArray{Float64,1,Array{Float64,1}}} == Bl
end

@testset "block-bidiagonal matrix" begin
    dv = [ones(4,3), ones(2,2).*2, ones(2,3).*3, ones(4,4).*4]
    evu = [ones(4,2), ones(2,3).*2, ones(2,4).*3]
    evl = [ones(2,3), ones(2,2).*2, ones(4,3).*3]
    BU = Bidiagonal(dv, evu, :U)
    BL = Bidiagonal(dv, evl, :L)
    # check that all the matrices along a column have the same number of columns,
    # and the matrices along a row have the same number of rows
    for j in axes(BU, 2), i in 2:size(BU, 1)
        @test size(BU[i,j], 2) == size(BU[1,j], 2)
        @test size(BU[i,j], 1) == size(BU[i,1], 1)
        if j < i || j > i + 1
            @test iszero(BU[i,j])
        end
    end
    for j in axes(BL, 2), i in 2:size(BL, 1)
        @test size(BL[i,j], 2) == size(BL[1,j], 2)
        @test size(BL[i,j], 1) == size(BL[i,1], 1)
        if j < i-1 || j > i
            @test iszero(BL[i,j])
        end
    end

    @test diag(BU, -1) == [zeros(size(dv[i+1], 1), size(dv[i],2)) for i in 1:length(dv)-1]
    @test diag(BL, 1) == [zeros(size(dv[i], 1), size(dv[i+1],2)) for i in 1:length(dv)-1]

    M = ones(2,2)
    for n in 0:1
        dv = fill(M, n)
        ev = fill(M, 0)
        B = Bidiagonal(dv, ev, :U)
        @test B == Matrix{eltype(B)}(B)
    end

    @testset "non-standard axes" begin
        s = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
        B = Bidiagonal(fill(s,4), fill(s,3), :U)
        @test @inferred(B[2,1]) isa typeof(s)
        @test all(iszero, B[2,1])
    end

    @testset "adjoint/transpose" begin
        m = rand(Int, 2, 2)
        for uplo in [:U, :L]
            B = Bidiagonal(fill(m,4), fill(m,3), uplo)
            A = Array{Matrix{Int}}(B)
            @testset for f in (adjoint, transpose)
                @test f(B) == f(A)
                @test f(f(B)) == B
            end
        end
    end
end

@testset "copyto!" begin
    ev, dv = [1:4;], [1:5;]
    B = Bidiagonal(dv, ev, :U)
    B2 = copyto!(zero(B), B)
    @test B2 == B
    for (ul1, ul2) in ((:U, :L), (:L, :U))
        B3 = Bidiagonal(dv, zero(ev), ul1)
        B2 = Bidiagonal(zero(dv), zero(ev), ul2)
        @test copyto!(B2, B3) == B3
    end

    @testset "mismatched sizes" begin
        dv2 = [4; @view dv[2:end]]
        @test copyto!(B, Bidiagonal([4], Int[], :U)) == Bidiagonal(dv2, ev, :U)
        @test copyto!(B, Bidiagonal([4], Int[], :L)) == Bidiagonal(dv2, ev, :U)
        @test copyto!(B, Bidiagonal(Int[], Int[], :U)) == Bidiagonal(dv, ev, :U)
        @test copyto!(B, Bidiagonal(Int[], Int[], :L)) == Bidiagonal(dv, ev, :U)
    end
end

@testset "copyto! with UniformScaling" begin
    @testset "Fill" begin
        for len in (4, InfiniteArrays.Infinity())
            d = FillArrays.Fill(1, len)
            ud = FillArrays.Fill(0, len-1)
            B = Bidiagonal(d, ud, :U)
            @test copyto!(B, I) === B
        end
    end
    B = Bidiagonal(fill(2, 4), fill(3, 3), :U)
    copyto!(B, I)
    @test all(isone, diag(B))
    @test all(iszero, diag(B, 1))
end

@testset "diagind" begin
    B = Bidiagonal(1:4, 1:3, :U)
    M = Matrix(B)
    @testset for k in -4:4
        @test B[diagind(B,k)] == M[diagind(M,k)]
    end
end

@testset "custom axes" begin
    dv, uv = OffsetArray(1:4), OffsetArray(1:3)
    B = Bidiagonal(dv, uv, :U)
    ax = axes(dv, 1)
    @test axes(B) === (ax, ax)
end

@testset "avoid matmul ambiguities with ::MyMatrix * ::AbstractMatrix" begin
    A = [i+j for i in 1:2, j in 1:2]
    S = SizedArrays.SizedArray{(2,2)}(A)
    B = Bidiagonal([1:2;], [1;], :U)
    @test S * B == A * B
    @test B * S == B * A
    C1, C2 = zeros(2,2), zeros(2,2)
    @test mul!(C1, S, B) == mul!(C2, A, B)
    @test mul!(C1, S, B, 1, 2) == mul!(C2, A, B, 1 ,2)
    @test mul!(C1, B, S) == mul!(C2, B, A)
    @test mul!(C1, B, S, 1, 2) == mul!(C2, B, A, 1 ,2)

    v = [i for i in 1:2]
    sv = SizedArrays.SizedArray{(2,)}(v)
    @test B * sv == B * v
    C1, C2 = zeros(2), zeros(2)
    @test mul!(C1, B, sv) == mul!(C2, B, v)
    @test mul!(C1, B, sv, 1, 2) == mul!(C2, B, v, 1 ,2)
end

@testset "Reverse operation on Bidiagonal" begin
    n = 5
    d = randn(n)
    e = randn(n - 1)
    for uplo in (:U, :L)
        B = Bidiagonal(d, e, uplo)
        @test reverse(B, dims=1) == reverse(Matrix(B), dims=1)
        @test reverse(B, dims=2) == reverse(Matrix(B), dims=2)
        @test reverse(B)::Bidiagonal == reverse(Matrix(B))
    end
end

@testset "Matrix conversion for non-numeric" begin
    B = Bidiagonal(fill(Diagonal([1,3]), 3), fill(Diagonal([1,3]), 2), :U)
    M = Matrix{eltype(B)}(B)
    @test M isa Matrix{eltype(B)}
    @test M == B
end

@testset "getindex with Integers" begin
    dv, ev = 1:4, 1:3
    B = Bidiagonal(dv, ev, :U)
    @test_throws "invalid index" B[3, true]
    @test B[1,2] == B[Int8(1),UInt16(2)] == B[big(1), Int16(2)]
end

@testset "rmul!/lmul! with banded matrices" begin
    dv, ev = rand(4), rand(3)
    for A in (Bidiagonal(dv, ev, :U), Bidiagonal(dv, ev, :L))
        @testset "$(nameof(typeof(B)))" for B in (
                                Bidiagonal(dv, ev, :U),
                                Bidiagonal(dv, ev, :L),
                                Diagonal(dv)
                        )
            @test_throws ArgumentError rmul!(B, A)
            @test_throws ArgumentError lmul!(A, B)
        end
    end
    @testset "non-commutative" begin
        S32 = SizedArrays.SizedArray{(3,2)}(rand(3,2))
        S33 = SizedArrays.SizedArray{(3,3)}(rand(3,3))
        S22 = SizedArrays.SizedArray{(2,2)}(rand(2,2))
        for uplo in (:L, :U)
            B = Bidiagonal(fill(S32, 4), fill(S32, 3), uplo)
            D = Diagonal(fill(S22, size(B,2)))
            @test rmul!(copy(B), D) ≈ B * D
            D = Diagonal(fill(S33, size(B,1)))
            @test lmul!(D, copy(B)) ≈ D * B
        end

        B = Bidiagonal(fill(S33, 4), fill(S33, 3), :U)
        D = Diagonal(fill(S32, 4))
        @test lmul!(B, Array(D)) ≈ B * D
        B = Bidiagonal(fill(S22, 4), fill(S22, 3), :U)
        @test rmul!(Array(D), B) ≈ D * B
    end
end

@testset "rmul!/lmul! with numbers" begin
    for T in (Bidiagonal(rand(4), rand(3), :U), Bidiagonal(rand(4), rand(3), :L))
        @test rmul!(copy(T), 0.2) ≈ rmul!(Array(T), 0.2)
        @test lmul!(0.2, copy(T)) ≈ lmul!(0.2, Array(T))
        @test_throws ArgumentError rmul!(T, NaN)
        @test_throws ArgumentError lmul!(NaN, T)
    end
    for T in (Bidiagonal(rand(1), rand(0), :U), Bidiagonal(rand(1), rand(0), :L))
        @test all(isnan, rmul!(copy(T), NaN))
        @test all(isnan, lmul!(NaN, copy(T)))
    end
end

@testset "mul with Diagonal" begin
    for n in 0:4
        dv, ev = rand(n), rand(max(n-1,0))
        d = rand(n)
        for uplo in (:U, :L)
            A = Bidiagonal(dv, ev, uplo)
            D = Diagonal(d)
            M = Matrix(A)
            S = similar(A, size(A))
            @test A * D ≈ mul!(S, A, D) ≈ M * D
            @test D * A ≈ mul!(S, D, A) ≈ D * M
            @test mul!(copy(S), D, A, 2, 2) ≈ D * M * 2 + S * 2
            @test mul!(copy(S), D, A, 0, 2) ≈ D * M * 0 + S * 2
            @test mul!(copy(S), A, D, 2, 2) ≈ M * D * 2 + S * 2
            @test mul!(copy(S), A, D, 0, 2) ≈ M * D * 0 + S * 2

            A2 = Bidiagonal(dv, zero(ev), uplo)
            M2 = Array(A2)
            S2 = Bidiagonal(copy(dv), copy(ev), uplo == (:U) ? (:L) : (:U))
            MS2 = Array(S2)
            @test mul!(copy(S2), D, A2) ≈ D * M2
            @test mul!(copy(S2), A2, D) ≈ M2 * D
            @test mul!(copy(S2), A2, D, 2, 2) ≈ M2 * D * 2 + MS2 * 2
            @test mul!(copy(S2), D, A2, 2, 2) ≈ D * M2 * 2 + MS2 * 2
        end
    end

    t1 = SizedArrays.SizedArray{(2,3)}([1 2 3; 3 4 5])
    t2 = SizedArrays.SizedArray{(3,2)}([1 2; 3 4; 5 6])
    dv, ev, d = fill(t1, 4), fill(2t1, 3), fill(t2, 4)
    for uplo in (:U, :L)
        A = Bidiagonal(dv, ev, uplo)
        D = Diagonal(d)
        @test A * D ≈ Array(A) * Array(D)
        @test D * A ≈ Array(D) * Array(A)
    end
end

@testset "conversion to Tridiagonal for immutable bands" begin
    n = 4
    dv = FillArrays.Fill(3, n)
    ev = FillArrays.Fill(2, n-1)
    z = FillArrays.Fill(0, n-1)
    dvf = FillArrays.Fill(Float64(3), n)
    evf = FillArrays.Fill(Float64(2), n-1)
    zf = FillArrays.Fill(Float64(0), n-1)
    B = Bidiagonal(dv, ev, :U)
    @test Tridiagonal{Int}(B) === Tridiagonal(B) === Tridiagonal(z, dv, ev)
    @test Tridiagonal{Float64}(B) === Tridiagonal(zf, dvf, evf)
    B = Bidiagonal(dv, ev, :L)
    @test Tridiagonal{Int}(B) === Tridiagonal(B) === Tridiagonal(ev, dv, z)
    @test Tridiagonal{Float64}(B) === Tridiagonal(evf, dvf, zf)
end

@testset "off-band indexing error" begin
    B = Bidiagonal(Vector{BigInt}(undef, 4), Vector{BigInt}(undef,3), :L)
    @test_throws "cannot set entry" B[1,2] = 4
end

@testset "mul with empty arrays" begin
    A = zeros(5,0)
    B = Bidiagonal(zeros(0), zeros(0), :U)
    BL = Bidiagonal(zeros(5), zeros(4), :U)
    @test size(A * B) == size(A)
    @test size(BL * A) == size(A)
    @test size(B * B) == size(B)
    C = similar(A)
    @test mul!(C, A, B) == A * B
    @test mul!(C, BL, A) == BL * A
    @test mul!(similar(B), B, B) == B * B
    @test mul!(similar(B, size(B)), B, B) == B * B

    v = zeros(size(B,2))
    @test size(B * v) == size(v)
    @test mul!(similar(v), B, v) == B * v

    D = Diagonal(zeros(size(B,2)))
    @test size(B * D) == size(D * B) == size(D)
    @test mul!(similar(D), B, D) == mul!(similar(D), D, B) == B * D
end

@testset "mul for small matrices" begin
    @testset for n in 0:6
        D = Diagonal(rand(n))
        v = rand(n)
        @testset for uplo in (:L, :U)
            B = Bidiagonal(rand(n), rand(max(n-1,0)), uplo)
            M = Matrix(B)

            @test B * v ≈ M * v
            @test mul!(similar(v), B, v) ≈ M * v
            @test mul!(ones(size(v)), B, v, 2, 3) ≈ M * v * 2 .+ 3
            @test mul!(ones(size(v)), B, v, 0, 3) ≈ M * v * 0 .+ 3

            @test B * B ≈ M * M
            @test mul!(similar(B, size(B)), B, B) ≈ M * M
            @test mul!(ones(size(B)), B, B, 2, 4) ≈ M * M * 2 .+ 4
            @test mul!(ones(size(B)), B, B, 0, 4) ≈ M * M * 0 .+ 4

            for m in 0:6
                AL = rand(m,n)
                AR = rand(n,m)
                @test AL * B ≈ AL * M
                @test B * AR ≈ M * AR
                @test mul!(similar(AL), AL, B) ≈ AL * M
                @test mul!(similar(AR), B, AR) ≈ M * AR
                @test mul!(ones(size(AL)), AL, B, 2, 4) ≈ AL * M * 2 .+ 4
                @test mul!(ones(size(AR)), B, AR, 2, 4) ≈ M * AR * 2 .+ 4
            end

            @test B * D ≈ M * D
            @test D * B ≈ D * M
            @test mul!(similar(B), B, D) ≈ M * D
            @test mul!(similar(B), B, D) ≈ M * D
            @test mul!(similar(B, size(B)), D, B) ≈ D * M
            @test mul!(similar(B, size(B)), B, D) ≈ M * D
            @test mul!(ones(size(B)), D, B, 2, 4) ≈ D * M * 2 .+ 4
            @test mul!(ones(size(B)), B, D, 2, 4) ≈ M * D * 2 .+ 4
        end
        BL = Bidiagonal(rand(n), rand(max(0, n-1)), :L)
        ML = Matrix(BL)
        BU = Bidiagonal(rand(n), rand(max(0, n-1)), :U)
        MU = Matrix(BU)
        T = Tridiagonal(zeros(max(0, n-1)), zeros(n), zeros(max(0, n-1)))
        @test mul!(T, BL, BU) ≈ ML * MU
        @test mul!(T, BU, BL) ≈ MU * ML
        T = Tridiagonal(ones(max(0, n-1)), ones(n), ones(max(0, n-1)))
        @test mul!(copy(T), BL, BU, 2, 3) ≈ ML * MU * 2 + T * 3
        @test mul!(copy(T), BU, BL, 2, 3) ≈ MU * ML * 2 + T * 3
    end

    n = 4
    arr = SizedArrays.SizedArray{(2,2)}(reshape([1:4;],2,2))
    for B in (
            Bidiagonal(fill(arr,n), fill(arr,n-1), :L),
            Bidiagonal(fill(arr,n), fill(arr,n-1), :U),
            )
        @test B * B ≈ Matrix(B) * Matrix(B)
        BL = Bidiagonal(fill(arr,n), fill(arr,n-1), :L)
        BU = Bidiagonal(fill(arr,n), fill(arr,n-1), :U)
        @test BL * B ≈ Matrix(BL) * Matrix(B)
        @test BU * B ≈ Matrix(BU) * Matrix(B)
        @test B * BL ≈ Matrix(B) * Matrix(BL)
        @test B * BU ≈ Matrix(B) * Matrix(BU)
        D = Diagonal(fill(arr,n))
        @test D * B ≈ Matrix(D) * Matrix(B)
        @test B * D ≈ Matrix(B) * Matrix(D)
    end
end

@testset "opnorms" begin
    B = Bidiagonal([1,-2,3,-4], [1,2,3], 'U')

    @test opnorm(B, 1) == opnorm(Matrix(B), 1)
    @test opnorm(B, 2) ≈ opnorm(Matrix(B), 2)
    @test opnorm(B, Inf) == opnorm(Matrix(B), Inf)

    B = Bidiagonal([1,-2,3,-4], [1,2,3], 'L')

    @test opnorm(B, 1) == opnorm(Matrix(B), 1)
    @test opnorm(B, 2) ≈ opnorm(Matrix(B), 2)
    @test opnorm(B, Inf) == opnorm(Matrix(B), Inf)

    B = Bidiagonal([2], Int[], 'L')

    @test opnorm(B, 1) == opnorm(Matrix(B), 1)
    @test opnorm(B, 2) ≈ opnorm(Matrix(B), 2)
    @test opnorm(B, Inf) == opnorm(Matrix(B), Inf)

    B = Bidiagonal([2], Int[], 'U')

    @test opnorm(B, 1) == opnorm(Matrix(B), 1)
    @test opnorm(B, 2) ≈ opnorm(Matrix(B), 2)
    @test opnorm(B, Inf) == opnorm(Matrix(B), Inf)
end

@testset "convert to Bidiagonal" begin
    M = diagm(0 => [1,2,3], 1=>[4,5])
    B = convert(Bidiagonal, M)
    @test B == Bidiagonal(M, :U)
    M = diagm(0 => [1,2,3], -1=>[4,5])
    B = convert(Bidiagonal, M)
    @test B == Bidiagonal(M, :L)
    B = convert(Bidiagonal{Int8}, M)
    @test B == M
    @test B isa Bidiagonal{Int8, Vector{Int8}}
    B = convert(Bidiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}, M)
    @test B == M
    @test B isa Bidiagonal{Int8, OffsetVector{Int8, Vector{Int8}}}
    M = diagm(-1 => [1,2], 1=>[4,5])
    @test_throws InexactError convert(Bidiagonal, M)
end

@testset "isreal" begin
    M = Bidiagonal(ones(2), ones(1), :U)
    @test @inferred((M -> Val(isreal(M)))(M)) == Val(true)
    M = complex.(M)
    @test isreal(M)
    @test !isreal(im*M)
end

@testset "ldiv! error message" begin
    C = zeros(2)
    B = Bidiagonal(1:0, 1:0, :U)
    msg = "size of result, (2,), does not match the size of b, (0, 1)"
    @test_throws msg ldiv!(C, B, zeros(0,1))
    msg = "the first dimension of the Bidiagonal matrix, 0, does not match the length of the right-hand-side, 2"
    @test_throws msg ldiv!(C, B, zeros(2))
    msg = "the first dimension of the Bidiagonal matrix, 0, does not match the first dimension of the right-hand-side, 2"
    @test_throws msg ldiv!(C, B, zeros(2,1))
end

@testset "l/rmul with 0-sized matrices" begin
    n = 0
    B = Bidiagonal(ones(n), ones(max(n-1,0)), :U)
    B2 = copy(B)
    D = Diagonal(ones(n))
    @test lmul!(D, B) == B2
    @test rmul!(B, D) == B2
end

@testset "setindex! with BandIndex" begin
    B = Bidiagonal(zeros(3), zeros(2), :U)
    B[LinearAlgebra.BandIndex(0,2)] = 1
    @test B[2,2] == 1
    B[LinearAlgebra.BandIndex(1,1)] = 2
    @test B[1,2] == 2
    @test_throws "cannot set entry $((1,3)) off the upper bidiagonal band" B[LinearAlgebra.BandIndex(2,1)] = 2

    B = Bidiagonal(zeros(3), zeros(2), :L)
    B[LinearAlgebra.BandIndex(-1,1)] = 2
    @test B[2,1] == 2
    @test_throws "cannot set entry $((3,1)) off the lower bidiagonal band" B[LinearAlgebra.BandIndex(-2,1)] = 2

    @test_throws BoundsError B[LinearAlgebra.BandIndex(size(B,1),1)]
    @test_throws BoundsError B[LinearAlgebra.BandIndex(0,size(B,1)+1)]
end

@testset "lazy adjtrans" begin
    B = Bidiagonal(fill([1 2; 3 4], 3), fill([5 6; 7 8], 2), :U)
    m = [2 4; 6 8]
    for op in (transpose, adjoint)
        C = op(B)
        el = op(m)
        C[1,1] = el
        @test B[1,1] == m
        C[2,1] = el
        @test B[1,2] == m
        @test (@allocated op(B)) == 0
        @test (@allocated op(op(B))) == 0
    end
end

@testset "fillband!" begin
    @testset "uplo = :U" begin
        B = Bidiagonal(zeros(4), zeros(3), :U)
        LinearAlgebra.fillband!(B, 2, 1, 1)
        @test all(==(2), diagview(B,1))
        LinearAlgebra.fillband!(B, 3, 0, 0)
        @test all(==(3), diagview(B,0))
        @test all(==(2), diagview(B,1))
        LinearAlgebra.fillband!(B, 4, 0, 1)
        @test all(==(4), diagview(B,0))
        @test all(==(4), diagview(B,1))
        @test_throws ArgumentError LinearAlgebra.fillband!(B, 3, -1, 0)

        LinearAlgebra.fillstored!(B, 1)
        LinearAlgebra.fillband!(B, 0, -3, 3)
        @test iszero(B)
        LinearAlgebra.fillband!(B, 0, -10, 10)
        @test iszero(B)
        LinearAlgebra.fillstored!(B, 1)
        B2 = copy(B)
        LinearAlgebra.fillband!(B, 0, -1, -3)
        @test B == B2
        LinearAlgebra.fillband!(B, 0, 10, 10)
        @test B == B2
    end

    @testset "uplo = :L" begin
        B = Bidiagonal(zeros(4), zeros(3), :L)
        LinearAlgebra.fillband!(B, 2, -1, -1)
        @test all(==(2), diagview(B,-1))
        LinearAlgebra.fillband!(B, 3, 0, 0)
        @test all(==(3), diagview(B,0))
        @test all(==(2), diagview(B,-1))
        LinearAlgebra.fillband!(B, 4, -1, 0)
        @test all(==(4), diagview(B,0))
        @test all(==(4), diagview(B,-1))
        @test_throws ArgumentError LinearAlgebra.fillband!(B, 3, 0, 1)

        LinearAlgebra.fillstored!(B, 1)
        LinearAlgebra.fillband!(B, 0, -3, 3)
        @test iszero(B)
        LinearAlgebra.fillband!(B, 0, -10, 10)
        @test iszero(B)
        LinearAlgebra.fillstored!(B, 1)
        B2 = copy(B)
        LinearAlgebra.fillband!(B, 0, -1, -3)
        @test B == B2
        LinearAlgebra.fillband!(B, 0, 10, 10)
        @test B == B2
    end
end

end # module TestBidiagonal
