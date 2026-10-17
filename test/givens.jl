# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestGivens

isdefined(Main, :pruned_old_LA) || @eval Main include("prune_old_LA.jl")

using Test, LinearAlgebra, Random
using LinearAlgebra: Givens, Rotation, givensAlgorithm

# Test givens rotations
@testset "Test Givens for $elty" for elty in (Float32, Float64, ComplexF32, ComplexF64)
    if elty <: Real
        raw_A = convert(Matrix{elty}, randn(10,10))
    else
        raw_A = convert(Matrix{elty}, complex.(randn(10,10),randn(10,10)))
    end
    @testset for A in (raw_A, view(raw_A, 1:10, 1:10))
        Ac = copy(A)
        R = Rotation(Givens{elty}[])
        T = Rotation(Givens{elty}[])
        for j = 1:8
            for i = j+2:10
                G, _ = givens(A, j+1, i, j)
                lmul!(G, A)
                rmul!(A, adjoint(G))
                lmul!(G, R)
                rmul!(T, G)

                @test lmul!(G, Matrix{elty}(I, 10, 10)) == [G[i,j] for i=1:10,j=1:10]

                @testset "transposes" begin
                    @test (@inferred G'*G)*Matrix(elty(1)I, 10, 10) ≈ Matrix(I, 10, 10)
                    @test (G*Matrix(elty(1)I, 10, 10))*G' ≈ Matrix(I, 10, 10)
                    @test (@inferred copy(R'))*(R*Matrix(elty(1)I, 10, 10)) ≈ Matrix(I, 10, 10)
                    @test_throws ErrorException transpose(G)
                    @test_throws ErrorException transpose(R)
                end
            end
        end
        @test (R')' === R
        # test products of Givens and Rotations
        for r in (R, T, *(R.rotations...), *(R.rotations[1], *(R.rotations[2:end]...)))
            @test r * A ≈ (A' * r')' ≈ lmul!(r, copy(A))
            @test A * r ≈ (r' * A')' ≈ rmul!(copy(A), r)
            @test r' * A ≈ lmul!(r', copy(A))
            @test A * r' ≈ rmul!(copy(A), r')
        end
        @test_throws ArgumentError givens(A, 3, 3, 2)
        @test_throws ArgumentError givens(one(elty),zero(elty),2,2)
        G, _ = givens(one(elty),zero(elty),11,12)
        @test_throws DimensionMismatch lmul!(G, A)
        @test_throws DimensionMismatch rmul!(A, adjoint(G))
        @test abs.(A) ≈ abs.(hessenberg(Ac).H)
        @test opnorm(R*Matrix{elty}(I, 10, 10)) ≈ one(elty)

        I10 = Matrix{elty}(I, 10, 10)
        G, _ = givens(one(elty),zero(elty),9,10)
        @test (G*I10)' * (G*I10) ≈ I10
        K, _ = givens(zero(elty),one(elty),9,10)
        @test (K*I10)' * (K*I10) ≈ I10
    end

    @testset "Givens * vectors" begin
        for x in (raw_A[:,1], view(raw_A, :, 1))
            G, r = @inferred  givens(x[2], x[4], 2, 4)
            @test (G*x)[2] ≈ r
            @test abs((G*x)[4]) < eps(real(elty))

            G, r = @inferred givens(x, 2, 4)
            @test (G*x)[2] ≈ r
            @test abs((G*x)[4]) < eps(real(elty))

            G, r = givens(x, 4, 2)
            @test (G*x)[4] ≈ r
            @test abs((G*x)[2]) < eps(real(elty))
        end
        d = rand(4)
        l = d[1]
        g2, l = givens(l, d[2], 1, 2)
        g3, l = givens(l, d[3], 1, 3)
        g4, l = givens(l, d[4], 1, 4)
        @test g2*(g3*d) ≈ g2*g3*d ≈ (g2*g3)*d
        @test g2*g3*g4 isa Rotation
    end
end

const TNumber = Union{Float64,ComplexF64}
struct MockUnitful{T<:TNumber} <: Number
    data::T
    MockUnitful(data::T) where T<:TNumber = new{T}(data)
end
import Base: *, /, one, oneunit
*(a::MockUnitful{T}, b::T) where T<:TNumber = MockUnitful(a.data * b)
*(a::T, b::MockUnitful{T}) where T<:TNumber = MockUnitful(a * b.data)
*(a::MockUnitful{T}, b::MockUnitful{T}) where T<:TNumber = MockUnitful(a.data * b.data)
/(a::MockUnitful{T}, b::MockUnitful{T}) where T<:TNumber = a.data / b.data
one(::Type{<:MockUnitful{T}}) where T = one(T)
oneunit(::Type{<:MockUnitful{T}}) where T = MockUnitful(one(T))

@testset "unitful givens rotation unitful $T " for T in (Float64, ComplexF64)
    g, r = givens(MockUnitful(T(3)), MockUnitful(T(4)), 1, 2)
    @test g.c ≈ 3/5
    @test g.s ≈ 4/5
    @test r.data ≈ 5.0
end

struct MockMeasurement{T<:AbstractFloat} <: AbstractFloat
    data::T
end
#these methods are only needed for preventing ambiguities
MockMeasurement{T}(x::MockMeasurement{T}) where {T<:AbstractFloat} = x
MockMeasurement{T}(z::Complex) where {T<:AbstractFloat} = MockMeasurement(T(z))
MockMeasurement{T}(r::Rational{P}) where {P,T<:AbstractFloat} = MockMeasurement(T(r))
MockMeasurement{T}(c::AbstractChar) where {T<:AbstractFloat} = MockMeasurement(T(c))
MockMeasurement{T}(x::Base.TwicePrecision) where {T<:AbstractFloat} = MockMeasurement(T(x))
#these methods are actually needed
import Base: promote_rule, floatmin, eps, ==, <, -, +, <=, sqrt
one(::Type{<:MockMeasurement{T}}) where {T<:Real} = one(T)
promote_rule(::Type{MockMeasurement{T}}, ::Type{<:Real}) where {T} = MockMeasurement{T}
for f in (:floatmin, :eps, :oneunit)
    @eval $f(::Type{MockMeasurement{T}}) where {T<:AbstractFloat} = $f(T)
end
for f in (:-, :sqrt)
    @eval $f(x::MockMeasurement) = MockMeasurement($f(x.data))
end
for f in (:*, :+, :/)
    @eval $f(x::MockMeasurement, y::MockMeasurement) = MockMeasurement($f(x.data,y.data))
end
for f in (:<, :(==), :(<=))
    @eval $f(x::MockMeasurement, y::MockMeasurement) = $f(x.data,y.data)
end

@testset "measurement givens rotation unitful $T " for T in (Float64, ComplexF64)
    g, r = givens(MockMeasurement(T(3)), MockMeasurement(T(4)), 1, 2)
    @test g.c ≈ 3/5
    @test g.s ≈ 4/5
    @test r.data ≈ 5.0
end

# 51554
# avoid infinite loop on Inf inputs
@testset "givensAlgorithm - Inf inputs" for T in (Float64, ComplexF64)
    cs, sn, r = givensAlgorithm(T(Inf), T(1.0))
    @test !isfinite(r)
    cs, sn, r = givensAlgorithm(T(1.0), T(Inf))
    @test !isfinite(r)
end

# exceptional cases, see Bindel et al., "On Computing Givens Rotations Reliably and Efficiently"
@testset "givensAlgorithm - non-finite inputs" begin
    for T in (Float32, Float64, BigFloat)
        for (f, g, cs, sn, r) in ((Inf, 2, 1, 0, Inf), (-Inf, 2, 1, 0, -Inf),
                                  (2, Inf, 0, 1, Inf), (-2, -Inf, 0, -1, Inf))
            res = givensAlgorithm(T(f), T(g))
            @test res isa NTuple{3,T}
            @test res == (cs, sn, r)
        end
        for (f, g) in ((Inf, Inf), (-Inf, Inf), (NaN, 1), (1, NaN), (NaN, Inf))
            @test all(isnan, givensAlgorithm(T(f), T(g)))
        end
    end
    for T in (Float32, Float64)
        CT = Complex{T}
        for (f, g, cs, sn, r) in ((CT(Inf, 1), CT(2, 3), 1, 0, CT(Inf, 1)),
                                  (CT(-Inf, Inf), CT(0, 1), 1, 0, CT(-Inf, Inf)),
                                  (CT(1, 2), CT(Inf, 0), 0, CT(1, 2)/sqrt(T(5)), CT(Inf, Inf)),
                                  (CT(0, -3), CT(1, -Inf), 0, 1, CT(0, -Inf)),
                                  (CT(-2, 0), CT(-Inf, -Inf), 0, CT(1, -1)/sqrt(T(2)), CT(-Inf, 0)),
                                  (CT(0), CT(1, -Inf), 0, CT(0, 1), CT(Inf, 0)))
            c, s, ρ = givensAlgorithm(f, g)
            @test c isa T && s isa CT && ρ isa CT
            @test c == cs && s ≈ sn && ρ == r
        end
        for (f, g) in ((CT(Inf), CT(0, Inf)), (CT(NaN), CT(1)), (CT(1), CT(1, NaN)))
            c, s, r = givensAlgorithm(f, g)
            @test isnan(c) && isnan(s) && isnan(r)
        end
    end
end

# The exceptional cases should be the limits of the finite cases: compare the result for
# infinite inputs to the result where each infinite component is replaced by a huge value.
@testset "givensAlgorithm - non-finite inputs are limits of large inputs" begin
    # replace ±Inf components by ±M
    finitize(x::Real, M) = isinf(x) ? copysign(M, x) : x
    finitize(z::Complex, M) = complex(finitize(real(z), M), finitize(imag(z), M))
    # does x (computed with huge values M) approximate the limit y (computed with Inf)?
    # The finite inputs are O(1), so the finite results converge like O(1/M).
    function islimit(x::Real, y::Real, M)
        isinf(y) ? sign(x) == sign(y) && abs(x) >= M/10 :
            isapprox(x, y; atol=sqrt(eps(typeof(x))) + 10/M)
    end
    islimit(x::Complex, y::Complex, M) =
        islimit(real(x), real(y), M) && islimit(imag(x), imag(y), M)
    islimit(x::Tuple, y::Tuple, M) = all(islimit(a, b, M) for (a, b) in zip(x, y))

    function check(f, g, T)
        res = givensAlgorithm(f, g)
        for M in (T(1e10), sqrt(floatmax(T)), floatmax(T)/4)
            if isinf(f) && isinf(g)
                # no limit: the rotation for huge inputs depends on their ratio
                @test all(isnan, res)
                c1, = givensAlgorithm(finitize(f, M), finitize(g, M))
                c2, = givensAlgorithm(finitize(f, M), finitize(g, M/2))
                @test isfinite(c1) && isfinite(c2) && !(c1 ≈ c2)
            else
                resM = givensAlgorithm(finitize(f, M), finitize(g, M))
                @test all(isfinite, resM)
                @test islimit(resM, res, M)
            end
        end
    end

    vals = (Inf, -Inf, 3, -2, 0)
    @testset for T in (Float32, Float64, BigFloat)
        for f in vals, g in vals
            (isinf(f) || isinf(g)) && check(T(f), T(g), T)
        end
    end
    @testset for T in (Float32, Float64)
        cvals = [complex(T(a), T(b)) for a in vals for b in vals]
        for f in cvals, g in cvals
            (isinf(f) || isinf(g)) && check(f, g, T)
        end
    end
end

# ordering of compositions
@testset "givens compositions ordering" begin
    R1, R2 = givens(1.,1.,1,2)[1], givens(1.,1.,2,3)[1]
    R = R1 * R2
    I3 = I(3)
    @test R1 * (R2 * I3) ≈ R * I3
    @test R2' * (R1' * I3) ≈ R' * I3
    @test (I3 * R1) * R2 ≈ I3 * R
    @test (I3 * R2') * R1' ≈ I3 * R'
    @test R' * (R * I3) ≈ R * (R' * I3) ≈ I
end

end # module TestGivens
