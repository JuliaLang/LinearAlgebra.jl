# This file is a part of Julia. License is MIT: https://julialang.org/license

import Base: copy, adjoint, getindex, show, transpose, one, zero, inv, float,
             hcat, vcat, hvcat, ^

"""
    UniformScaling{T<:Number}

Generically sized uniform scaling operator defined as a scalar times
the identity operator, `λ*I`. Although without an explicit `size`, it
acts similarly to a matrix in many cases and includes support for some
indexing. See also [`I`](@ref).

!!! compat "Julia 1.6"
     Indexing using ranges is available as of Julia 1.6.

# Examples
```jldoctest
julia> J = UniformScaling(2.)
UniformScaling{Float64}
2.0*I

julia> A = [1. 2.; 3. 4.]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> J*A
2×2 Matrix{Float64}:
 2.0  4.0
 6.0  8.0

julia> J[1:2, 1:2]
2×2 Matrix{Float64}:
 2.0  0.0
 0.0  2.0
```
"""
struct UniformScaling{T<:Number}
    λ::T
end

"""
    I

An object of type [`UniformScaling`](@ref), representing an identity matrix of any size.

# Examples
```jldoctest
julia> fill(1, (5,6)) * I == fill(1, (5,6))
true

julia> [1 2im 3; 1im 2 3] * I
2×3 Matrix{Complex{Int64}}:
 1+0im  0+2im  3+0im
 0+1im  2+0im  3+0im
```
"""
const I = UniformScaling(true)

"""
    (I::UniformScaling)(n::Integer)

Construct a `Diagonal` matrix from a `UniformScaling`.

!!! compat "Julia 1.2"
     This method is available as of Julia 1.2.

# Examples
```jldoctest
julia> I(3)
3×3 Diagonal{Bool, Vector{Bool}}:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1

julia> (0.7*I)(3)
3×3 Diagonal{Float64, Vector{Float64}}:
 0.7   ⋅    ⋅
  ⋅   0.7   ⋅
  ⋅    ⋅   0.7
```
"""
(I::UniformScaling)(n::Integer) = Diagonal(fill(I.λ, n))

eltype(::Type{UniformScaling{T}}) where {T} = T
ndims(J::UniformScaling) = 2
Base.has_offset_axes(::UniformScaling) = false
getindex(J::UniformScaling, ind::CartesianIndex{2}) = J[Tuple(ind)...]
getindex(J::UniformScaling, i::Integer,j::Integer) = ifelse(i==j,J.λ,zero(J.λ))

getindex(J::UniformScaling, n::Integer, m::AbstractVector{<:Integer}) = getindex(J, m, n)
function getindex(J::UniformScaling{T}, n::AbstractVector{<:Integer}, m::Integer) where T
    v = zeros(T, axes(n))
    @inbounds for (i,ii) in pairs(n)
        if ii == m
            v[i] = J.λ
        end
    end
    return v
end

function getindex(J::UniformScaling{T}, n::AbstractVector{<:Integer}, m::AbstractVector{<:Integer}) where T
    A = zeros(T, axes(n)..., axes(m)...)
    @inbounds for (j,jj) in pairs(m), (i,ii) in pairs(n)
        if ii == jj
            A[i,j] = J.λ
        end
    end
    return A
end

function show(io::IO, ::MIME"text/plain", J::UniformScaling)
    s = "$(J.λ)"
    if occursin(r"\w+\s*[\+\-]\s*\w+", s)
        s = "($s)"
    end
    print(io, typeof(J), "\n$s*I")
end
copy(J::UniformScaling) = UniformScaling(J.λ)

Base.convert(::Type{UniformScaling{T}}, J::UniformScaling) where {T} = UniformScaling(convert(T, J.λ))::UniformScaling{T}

conj(J::UniformScaling) = UniformScaling(conj(J.λ))
real(J::UniformScaling) = UniformScaling(real(J.λ))
imag(J::UniformScaling) = UniformScaling(imag(J.λ))

float(J::UniformScaling) = UniformScaling(float(J.λ))

transpose(J::UniformScaling) = J
adjoint(J::UniformScaling) = UniformScaling(conj(J.λ))

one(::Type{UniformScaling{T}}) where {T} = UniformScaling(one(T))
one(J::UniformScaling{T}) where {T} = one(UniformScaling{T})
oneunit(::Type{UniformScaling{T}}) where {T} = UniformScaling(oneunit(T))
oneunit(J::UniformScaling{T}) where {T} = oneunit(UniformScaling{T})
zero(::Type{UniformScaling{T}}) where {T} = UniformScaling(zero(T))
zero(J::UniformScaling{T}) where {T} = zero(UniformScaling{T})

isdiag(::UniformScaling) = true
istriu(::UniformScaling) = true
istril(::UniformScaling) = true
issymmetric(::UniformScaling) = true
ishermitian(J::UniformScaling) = isreal(J.λ)
isposdef(J::UniformScaling) = isposdef(J.λ)

(+)(J::UniformScaling, x::Number) = J.λ + x
(+)(x::Number, J::UniformScaling) = x + J.λ
(-)(J::UniformScaling, x::Number) = J.λ - x
(-)(x::Number, J::UniformScaling) = x - J.λ

(+)(J::UniformScaling)                      = UniformScaling(+J.λ)
(+)(J1::UniformScaling, J2::UniformScaling) = UniformScaling(J1.λ+J2.λ)
(+)(B::BitArray{2}, J::UniformScaling)      = Array(B) + J
(+)(J::UniformScaling, B::BitArray{2})      = J + Array(B)
(+)(J::UniformScaling, A::AbstractMatrix)   = A + J

(-)(J::UniformScaling)                      = UniformScaling(-J.λ)
(-)(J1::UniformScaling, J2::UniformScaling) = UniformScaling(J1.λ-J2.λ)
(-)(B::BitArray{2}, J::UniformScaling)      = Array(B) - J
(-)(J::UniformScaling, B::BitArray{2})      = J - Array(B)
(-)(A::AbstractMatrix, J::UniformScaling)   = A + (-J)

# matrix functions
for f in ( :exp,   :log, :cis,
           :expm1, :log1p,
           :sqrt,  :cbrt,
           :sin,   :cos,   :tan,
           :asin,  :acos,  :atan,
           :csc,   :sec,   :cot,
           :acsc,  :asec,  :acot,
           :sinh,  :cosh,  :tanh,
           :asinh, :acosh, :atanh,
           :csch,  :sech,  :coth,
           :acsch, :asech, :acoth )
    @eval Base.$f(J::UniformScaling) = UniformScaling($f(J.λ))
end
for f in (:sincos, :sincosd)
    @eval Base.$f(J::UniformScaling) = map(UniformScaling, $f(J.λ))
end

# Unit{Lower/Upper}Triangular matrices become {Lower/Upper}Triangular under
# addition with a UniformScaling
for (t1, t2) in ((:UnitUpperTriangular, :UpperTriangular),
                 (:UnitLowerTriangular, :LowerTriangular))
    @eval begin
        function (+)(UL::$t1, J::UniformScaling)
            ULnew = copymutable_oftype(UL.data, Base.promote_op(+, eltype(UL), typeof(J)))
            for i in axes(ULnew, 1)
                ULnew[i,i] = one(ULnew[i,i]) + J
            end
            return ($t2)(ULnew)
        end
    end
end

# Adding a complex UniformScaling to the diagonal of a Hermitian
# matrix breaks the hermiticity, if the UniformScaling is non-real.
# However, to preserve type stability, we do not special-case a
# UniformScaling{<:Complex} that happens to be real.
function (+)(A::Hermitian, J::UniformScaling{<:Complex})
    TS = Base.promote_op(+, eltype(A), typeof(J))
    B = copytri!(copymutable_oftype(parent(A), TS), A.uplo, true)
    for i in diagind(B, IndexStyle(B))
        B[i] = A[i] + J
    end
    return B
end

function (-)(J::UniformScaling{<:Complex}, A::Hermitian)
    TS = Base.promote_op(+, eltype(A), typeof(J))
    B = copytri!(copymutable_oftype(parent(A), TS), A.uplo, true)
    B .= .-B
    for i in diagind(B, IndexStyle(B))
        B[i] = J - A[i]
    end
    return B
end

function (+)(A::AdjOrTransAbsMat, J::UniformScaling)
    checksquare(A)
    op = wrapperop(A)
    op(op(A) + op(J))
end
function (-)(J::UniformScaling, A::AdjOrTransAbsMat)
    checksquare(A)
    op = wrapperop(A)
    op(op(J) - op(A))
end

function (+)(A::AbstractMatrix, J::UniformScaling)
    checksquare(A)
    B = copymutable_oftype(A, Base.promote_op(+, eltype(A), typeof(J)))
    for i in intersect(axes(A,1), axes(A,2))
        @inbounds B[i,i] += J
    end
    return B
end

function (-)(J::UniformScaling, A::AbstractMatrix)
    checksquare(A)
    B = convert(AbstractMatrix{Base.promote_op(+, eltype(A), typeof(J))}, -A)
    for i in intersect(axes(A,1), axes(A,2))
        @inbounds B[i,i] += J
    end
    return B
end

inv(J::UniformScaling) = UniformScaling(inv(J.λ))
opnorm(J::UniformScaling, p::Real=2) = opnorm(J.λ, p)

pinv(J::UniformScaling) = ifelse(iszero(J.λ),
                          UniformScaling(zero(inv(J.λ))),  # type stability
                          UniformScaling(inv(J.λ)))

function det(J::UniformScaling{T}) where T
    if isone(J.λ)
        one(T)
    elseif iszero(J.λ)
        zero(T)
    else
        throw(ArgumentError("Determinant of UniformScaling is only well-defined when λ = 0 or 1."))
    end
end

function tr(J::UniformScaling{T}) where T
    if iszero(J.λ)
        zero(T)
    else
        throw(ArgumentError("Trace of UniformScaling is only well-defined when λ = 0"))
    end
end

*(J1::UniformScaling, J2::UniformScaling) = UniformScaling(J1.λ*J2.λ)
*(B::BitArray{2}, J::UniformScaling) = *(Array(B), J::UniformScaling)
*(J::UniformScaling, B::BitArray{2}) = *(J::UniformScaling, Array(B))
*(A::AbstractMatrix, J::UniformScaling) = A*J.λ
*(v::AbstractVector, J::UniformScaling) = reshape(v, length(v), 1) * J
*(J::UniformScaling, A::AbstractVecOrMat) = J.λ*A
*(x::Number, J::UniformScaling) = UniformScaling(x*J.λ)
*(J::UniformScaling, x::Number) = UniformScaling(J.λ*x)

/(J1::UniformScaling, J2::UniformScaling) = J2.λ == 0 ? throw(SingularException(1)) : UniformScaling(J1.λ/J2.λ)
/(J::UniformScaling, A::AbstractMatrix) =
    (invA = inv(A); lmul!(J.λ, convert(AbstractMatrix{promote_type(eltype(J),eltype(invA))}, invA)))
/(A::AbstractMatrix, J::UniformScaling) = J.λ == 0 ? throw(SingularException(1)) : A/J.λ
/(v::AbstractVector, J::UniformScaling) = reshape(v, length(v), 1) / J

/(J::UniformScaling, x::Number) = UniformScaling(J.λ/x)
//(J::UniformScaling, x::Number) = UniformScaling(J.λ//x)

\(J1::UniformScaling, J2::UniformScaling) = J1.λ == 0 ? throw(SingularException(1)) : UniformScaling(J1.λ\J2.λ)
\(J::UniformScaling, A::AbstractVecOrMat) = J.λ == 0 ? throw(SingularException(1)) : J.λ\A
\(A::AbstractMatrix, J::UniformScaling) =
    (invA = inv(A); rmul!(convert(AbstractMatrix{promote_type(eltype(invA),eltype(J))}, invA), J.λ))
\(F::Factorization, J::UniformScaling) = F \ J(size(F,1))

\(x::Number, J::UniformScaling) = UniformScaling(x\J.λ)

@inline mul!(C::AbstractMatrix, A::AbstractMatrix, J::UniformScaling, alpha::Number, beta::Number) =
    mul!(C, A, J.λ, alpha, beta)
@inline mul!(C::AbstractVecOrMat, J::UniformScaling, B::AbstractVecOrMat, alpha::Number, beta::Number) =
    mul!(C, J.λ, B, alpha, beta)

function mul!(out::AbstractMatrix{T}, a::Number, B::UniformScaling, α::Number, β::Number) where {T}
    checksquare(out)
    if iszero(β)  # zero contribution of the out matrix
        fill!(out, zero(T))
    elseif !isone(β)
        rmul!(out, β)
    end
    s = convert(T, a*B.λ*α)
    if !iszero(s)
        @inbounds for i in diagind(out, IndexStyle(out))
            out[i] += s
        end
    end
    return out
end
@inline mul!(out::AbstractMatrix, A::UniformScaling, b::Number, α::Number, β::Number)=
    mul!(out, A.λ, UniformScaling(b), α, β)
rmul!(A::AbstractMatrix, J::UniformScaling) = rmul!(A, J.λ)
lmul!(J::UniformScaling, B::AbstractVecOrMat) = lmul!(J.λ, B)
rdiv!(A::AbstractMatrix, J::UniformScaling) = rdiv!(A, J.λ)
ldiv!(J::UniformScaling, B::AbstractVecOrMat) = ldiv!(J.λ, B)
ldiv!(Y::AbstractVecOrMat, J::UniformScaling, B::AbstractVecOrMat) = (Y .= J.λ .\ B)

Broadcast.broadcasted(::typeof(*), x::Number,J::UniformScaling) = UniformScaling(x*J.λ)
Broadcast.broadcasted(::typeof(*), J::UniformScaling,x::Number) = UniformScaling(J.λ*x)

Broadcast.broadcasted(::typeof(/), J::UniformScaling,x::Number) = UniformScaling(J.λ/x)

Broadcast.broadcasted(::typeof(\), x::Number,J::UniformScaling) = UniformScaling(x\J.λ)

(^)(J::UniformScaling, x::Number) = UniformScaling((J.λ)^x)
Base.literal_pow(::typeof(^), J::UniformScaling, x::Val) = UniformScaling(Base.literal_pow(^, J.λ, x))

Broadcast.broadcasted(::typeof(^), J::UniformScaling, x::Number) = UniformScaling(J.λ^x)
function Broadcast.broadcasted(::typeof(Base.literal_pow), ::typeof(^), J::UniformScaling, x::Val)
    UniformScaling(Base.literal_pow(^, J.λ, x))
end

==(J1::UniformScaling,J2::UniformScaling) = (J1.λ == J2.λ)

## equality comparison with UniformScaling
==(J::UniformScaling, A::AbstractMatrix) = A == J
function ==(A::AbstractMatrix, J::UniformScaling)
    require_one_based_indexing(A)
    size(A, 1) == size(A, 2) || return false
    isempty(A) && return true
    # Check that the elements of A are equal to those of J,
    # this ensures that if A == J, their elements are equal as well
    iszero(J.λ) && return first(A) == J.λ && iszero(A)
    isone(J.λ) && return first(A) == J.λ && isone(A)
    return _isequalto_uniformscaling(A, J)
end
function _isequalto_uniformscaling(A::AbstractMatrix, J::UniformScaling)
    return isdiag(A) && all(==(J.λ), diagview(A))
end
function _isequalto_uniformscaling(A::StridedMatrix, J::UniformScaling)
    for j in axes(A, 2), i in axes(A, 1)
        ifelse(i == j, A[i, j] == J.λ, iszero(A[i, j])) || return false
    end
    return true
end

isequal(A::AbstractMatrix, J::UniformScaling) = false
isequal(J::UniformScaling, A::AbstractMatrix) = false

function isapprox(J1::UniformScaling{T}, J2::UniformScaling{S};
            atol::Real=0, rtol::Real=Base.rtoldefault(T,S,atol), nans::Bool=false) where {T<:Number,S<:Number}
    isapprox(J1.λ, J2.λ, rtol=rtol, atol=atol, nans=nans)
end
function isapprox(J::UniformScaling, A::AbstractMatrix;
                  atol::Real = 0,
                  rtol::Real = Base.rtoldefault(promote_leaf_eltypes(A), eltype(J), atol),
                  nans::Bool = false, norm::Function = norm)
    n = checksquare(A)
    normJ = norm === opnorm             ? abs(J.λ) :
            norm === LinearAlgebra.norm ? abs(J.λ) * sqrt(n) :
                                          norm(Diagonal(fill(J.λ, n)))
    return norm(A - J) <= max(atol, rtol * max(norm(A), normJ))
end
isapprox(A::AbstractMatrix, J::UniformScaling; kwargs...) = isapprox(J, A; kwargs...)

"""
    copyto!(dest::AbstractMatrix, src::UniformScaling)

Copies a [`UniformScaling`](@ref) onto a matrix.

!!! compat "Julia 1.1"
    In Julia 1.0 this method only supported a square destination matrix. Julia 1.1. added
    support for a rectangular matrix.
"""
function copyto!(A::AbstractMatrix, J::UniformScaling)
    require_one_based_indexing(A)
    fill!(A, 0)
    λ = J.λ
    for i = 1:min(size(A,1),size(A,2))
        @inbounds A[i,i] = λ
    end
    return A
end

function copyto!(A::Diagonal, J::UniformScaling)
    A.diag .= J.λ
    return A
end
function copyto!(A::Union{Bidiagonal, SymTridiagonal}, J::UniformScaling)
    A.ev .= 0
    A.dv .= J.λ
    return A
end
function copyto!(A::Tridiagonal, J::UniformScaling)
    A.dl .= 0
    A.du .= 0
    A.d .= J.λ
    return A
end

"""
    copy!(dest::AbstractMatrix, src::UniformScaling)

Copies a [`UniformScaling`](@ref) onto a matrix.

!!! compat "Julia 1.12"
    This method is available as of Julia 1.12.
"""
Base.copy!(A::AbstractMatrix, J::UniformScaling) = copyto!(A, J)

function cond(J::UniformScaling{T}) where T
    onereal = inv(one(real(J.λ)))
    return J.λ ≠ zero(T) ? onereal : oftype(onereal, Inf)
end

## Matrix construction from UniformScaling
function Matrix{T}(s::UniformScaling, dims::Dims{2}) where {T}
    A = zeros(T, dims)
    v = T(s.λ)
    for i in diagind(dims...)
        @inbounds A[i] = v
    end
    return A
end
Matrix{T}(s::UniformScaling, m::Integer, n::Integer) where {T} = Matrix{T}(s, Dims((m, n)))
Matrix(s::UniformScaling, m::Integer, n::Integer) = Matrix(s, Dims((m, n)))
Matrix(s::UniformScaling, dims::Dims{2}) = Matrix{eltype(s)}(s, dims)
Array{T}(s::UniformScaling, dims::Dims{2}) where {T} = Matrix{T}(s, dims)
Array{T}(s::UniformScaling, m::Integer, n::Integer) where {T} = Matrix{T}(s, m, n)
Array(s::UniformScaling, m::Integer, n::Integer) = Matrix(s, m, n)
Array(s::UniformScaling, dims::Dims{2}) = Matrix(s, dims)

dot(A::AbstractMatrix, J::UniformScaling) = dot(tr(A), J.λ)
dot(J::UniformScaling, A::AbstractMatrix) = dot(J.λ, tr(A))

dot(x::AbstractVector, J::UniformScaling, y::AbstractVector) = dot(x, J.λ, y)
dot(x::AbstractVector, a::Number, y::AbstractVector) = sum(t -> dot(t[1], a, t[2]), zip(x, y))
dot(x::AbstractVector, a::Union{Real,Complex}, y::AbstractVector) = a*dot(x, y)

# muladd
Base.muladd(A::UniformScaling, B::UniformScaling, z::UniformScaling) =
    UniformScaling(A.λ * B.λ + z.λ)
