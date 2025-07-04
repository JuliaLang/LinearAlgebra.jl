# This file is a part of Julia. License is MIT: https://julialang.org/license

## Diagonal matrices

struct Diagonal{T,V<:AbstractVector{T}} <: AbstractMatrix{T}
    diag::V

    function Diagonal{T,V}(diag) where {T,V<:AbstractVector{T}}
        require_one_based_indexing(diag)
        new{T,V}(diag)
    end
end
Diagonal(v::AbstractVector{T}) where {T} = Diagonal{T,typeof(v)}(v)
Diagonal{T}(v::AbstractVector) where {T} = Diagonal(convert(AbstractVector{T}, v)::AbstractVector{T})

function Base.promote_rule(A::Type{<:Diagonal{<:Any,V}}, B::Type{<:Diagonal{<:Any,W}}) where {V,W}
    X = promote_type(V, W)
    T = eltype(X)
    isconcretetype(T) && return Diagonal{T,X}
    return typejoin(A, B)
end

"""
    Diagonal(V::AbstractVector)

Construct a lazy matrix with `V` as its diagonal.

See also [`UniformScaling`](@ref) for the lazy identity matrix `I`,
[`diagm`](@ref) to make a dense matrix, and [`diag`](@ref) to extract diagonal elements.

# Examples
```jldoctest
julia> d = Diagonal([1, 10, 100])
3×3 Diagonal{$Int, Vector{$Int}}:
 1   ⋅    ⋅
 ⋅  10    ⋅
 ⋅   ⋅  100

julia> diagm([7, 13])
2×2 Matrix{$Int}:
 7   0
 0  13

julia> ans + I
2×2 Matrix{Int64}:
 8   0
 0  14

julia> I(2)
2×2 Diagonal{Bool, Vector{Bool}}:
 1  ⋅
 ⋅  1
```

!!! note
    A one-column matrix is not treated like a vector, but instead calls the
    method `Diagonal(A::AbstractMatrix)` which extracts 1-element `diag(A)`:

```jldoctest
julia> A = transpose([7.0 13.0])
2×1 transpose(::Matrix{Float64}) with eltype Float64:
  7.0
 13.0

julia> Diagonal(A)
1×1 Diagonal{Float64, Vector{Float64}}:
 7.0
```
"""
Diagonal(V::AbstractVector)

"""
    Diagonal(A::AbstractMatrix)

Construct a matrix from the principal diagonal of `A`.
The input matrix `A` may be rectangular, but the output will
be square.

# Examples
```jldoctest
julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> D = Diagonal(A)
2×2 Diagonal{Int64, Vector{Int64}}:
 1  ⋅
 ⋅  4

julia> A = [1 2 3; 4 5 6]
2×3 Matrix{Int64}:
 1  2  3
 4  5  6

julia> Diagonal(A)
2×2 Diagonal{Int64, Vector{Int64}}:
 1  ⋅
 ⋅  5
```
"""
Diagonal(A::AbstractMatrix) = Diagonal(diag(A))
Diagonal{T}(A::AbstractMatrix) where T = Diagonal{T}(diag(A))
Diagonal{T,V}(A::AbstractMatrix) where {T,V<:AbstractVector{T}} = Diagonal{T,V}(diag(A))
function convert(::Type{T}, A::AbstractMatrix) where T<:Diagonal
    checksquare(A)
    isdiag(A) ? T(A) : throw(InexactError(:convert, T, A))
end

Diagonal(D::Diagonal) = D
Diagonal{T}(D::Diagonal{T}) where {T} = D
Diagonal{T}(D::Diagonal) where {T} = Diagonal{T}(D.diag)

AbstractMatrix{T}(D::Diagonal) where {T} = Diagonal{T}(D)
AbstractMatrix{T}(D::Diagonal{T}) where {T} = copy(D)
Matrix(D::Diagonal{T}) where {T} = Matrix{promote_type(T, typeof(zero(T)))}(D)
Matrix(D::Diagonal{Any}) = Matrix{Any}(D)
Array(D::Diagonal{T}) where {T} = Matrix(D)
function Matrix{T}(D::Diagonal) where {T}
    B = Matrix{T}(undef, size(D))
    if haszero(T) # optimized path for types with zero(T) defined
        size(B,1) > 1 && fill!(B, zero(T))
        copyto!(diagview(B), D.diag)
    else
        copyto!(B, D)
    end
    return B
end

"""
    Diagonal{T}(undef, n)

Construct an uninitialized `Diagonal{T}` of length `n`. See `undef`.
"""
Diagonal{T}(::UndefInitializer, n::Integer) where T = Diagonal(Vector{T}(undef, n))

similar(D::Diagonal, ::Type{T}) where {T} = Diagonal(similar(D.diag, T))
similar(D::Diagonal, ::Type{T}, dims::Union{Dims{1},Dims{2}}) where {T} = similar(D.diag, T, dims)

# copyto! for matching axes
_copyto_banded!(D1::Diagonal, D2::Diagonal) = (copyto!(D1.diag, D2.diag); D1)

size(D::Diagonal) = (n = length(D.diag); (n,n))

axes(D::Diagonal) = (ax = axes(D.diag, 1); (ax, ax))

@inline function Base.isassigned(D::Diagonal, i::Int, j::Int)
    @boundscheck checkbounds(Bool, D, i, j) || return false
    if i == j
        @inbounds r = isassigned(D.diag, i)
    else
        r = true
    end
    r
end

@inline function Base.isstored(D::Diagonal, i::Int, j::Int)
    @boundscheck checkbounds(D, i, j)
    if i == j
        @inbounds r = Base.isstored(D.diag, i)
    else
        r = false
    end
    r
end

function Base.minimum(D::Diagonal{T}) where T <: Number
    mindiag = minimum(D.diag)
    size(D, 1) > 1 && return (min(zero(T), mindiag))
    return mindiag
end

function Base.maximum(D::Diagonal{T}) where T <: Number
    maxdiag = Base.maximum(D.diag)
    size(D, 1) > 1 && return (max(zero(T), maxdiag))
    return maxdiag
end

@inline function getindex(D::Diagonal, i::Int, j::Int)
    @boundscheck checkbounds(D, i, j)
    if i == j
        @inbounds r = D.diag[i]
    else
        r = diagzero(D, i, j)
    end
    r
end
"""
    diagzero(A::AbstractMatrix, i, j)

Return the appropriate zero element `A[i, j]` corresponding to a banded matrix `A`.
"""
diagzero(A::AbstractMatrix, i, j) = zero(eltype(A))
@propagate_inbounds diagzero(A::AbstractMatrix{M}, i, j) where {M<:AbstractMatrix} =
    zeroslike(M, axes(A[BandIndex(0, i)], 1), axes(A[BandIndex(0, j)], 2))
@propagate_inbounds diagzero(A::AbstractMatrix, inds...) = diagzero(A, to_indices(A, inds)...)
# dispatching on the axes permits specializing on the axis types to return something other than an Array
zeroslike(M::Type, ax::Vararg{Union{AbstractUnitRange, Integer}}) = zeroslike(M, ax)
"""
    zeroslike(::Type{M}, ax::Tuple{AbstractUnitRange, Vararg{AbstractUnitRange}}) where {M<:AbstractMatrix}
    zeroslike(::Type{M}, sz::Tuple{Integer, Vararg{Integer}}) where {M<:AbstractMatrix}

Return an appropriate zero-ed array similar to `M`, with either the axes `ax` or the size `sz`.
This will be used as a structural zero element of a matrix-valued banded matrix.
By default, `zeroslike` falls back to using the size along each axis to construct the array.
"""
zeroslike(M::Type, ax::Tuple{AbstractUnitRange, Vararg{AbstractUnitRange}}) = zeroslike(M, map(length, ax))
zeroslike(M::Type, sz::Tuple{Integer, Vararg{Integer}}) = zeros(M, sz)
zeroslike(::Type{M}, sz::Tuple{Integer, Vararg{Integer}}) where {M<:AbstractMatrix} = zeros(eltype(M), sz)

@inline function getindex(D::Diagonal, b::BandIndex)
    @boundscheck checkbounds(D, b)
    if b.band == 0
        @inbounds r = D.diag[b.index]
    else
        r = diagzero(D, b)
    end
    r
end

@inline function setindex!(D::Diagonal, v, i::Int, j::Int)
    @boundscheck checkbounds(D, i, j)
    if i == j
        @inbounds D.diag[i] = v
    elseif !iszero(v)
        throw(ArgumentError(lazy"cannot set off-diagonal entry ($i, $j) to a nonzero value ($v)"))
    end
    return D
end

@inline function setindex!(D::Diagonal, v, b::BandIndex)
    @boundscheck checkbounds(D, b)
    if b.band == 0
        @inbounds D.diag[b.index] = v
    elseif !iszero(v)
        throw(ArgumentError(lazy"cannot set off-diagonal entry $(to_indices(D, (b,))) to a nonzero value ($v)"))
    end
    return D
end

## structured matrix methods ##
function Base.replace_in_print_matrix(A::Diagonal,i::Integer,j::Integer,s::AbstractString)
    i==j ? s : Base.replace_with_centered_mark(s)
end
function Base.show(io::IO, A::Diagonal)
    print(io, "Diagonal(")
    show(io, A.diag)
    print(io, ")")
end

parent(D::Diagonal) = D.diag

copy(D::Diagonal) = Diagonal(copy(D.diag))

Base._reverse(A::Diagonal, dims) = reverse!(Matrix(A); dims)
Base._reverse(A::Diagonal, ::Colon) = Diagonal(reverse(A.diag))
Base._reverse!(A::Diagonal, ::Colon) = (reverse!(A.diag); A)

ishermitian(D::Diagonal) = all(ishermitian, D.diag)
issymmetric(D::Diagonal) = all(issymmetric, D.diag)
isposdef(D::Diagonal) = all(isposdef, D.diag)

factorize(D::Diagonal) = D

real(D::Diagonal) = Diagonal(real(D.diag))
imag(D::Diagonal) = Diagonal(imag(D.diag))

isreal(D::Diagonal) = isreal(D.diag)

iszero(D::Diagonal) = all(iszero, D.diag)
isone(D::Diagonal) = all(isone, D.diag)
isdiag(D::Diagonal) = all(isdiag, D.diag)
isdiag(D::Diagonal{<:Number}) = true
Base.@constprop :aggressive istriu(D::Diagonal, k::Integer=0) = k <= 0 || iszero(D.diag) ? true : false
Base.@constprop :aggressive istril(D::Diagonal, k::Integer=0) = k >= 0 || iszero(D.diag) ? true : false
function triu!(D::Diagonal{T}, k::Integer=0) where T
    n = size(D,1)
    if !(-n + 1 <= k <= n + 1)
        throw(ArgumentError(string("the requested diagonal, $k, must be at least ",
            "$(-n + 1) and at most $(n + 1) in an $n-by-$n matrix")))
    elseif k > 0
        fill!(D.diag, zero(T))
    end
    return D
end

function tril!(D::Diagonal{T}, k::Integer=0) where T
    n = size(D,1)
    if !(-n - 1 <= k <= n - 1)
        throw(ArgumentError(LazyString(lazy"the requested diagonal, $k, must be at least ",
            lazy"$(-n - 1) and at most $(n - 1) in an $n-by-$n matrix")))
    elseif k < 0
        fill!(D.diag, zero(T))
    end
    return D
end

(==)(Da::Diagonal, Db::Diagonal) = Da.diag == Db.diag
(-)(A::Diagonal) = Diagonal(-A.diag)
(+)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag + Db.diag)
(-)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag - Db.diag)

(*)(x::Number, D::Diagonal) = Diagonal(x * D.diag)
(*)(D::Diagonal, x::Number) = Diagonal(D.diag * x)
function lmul!(x::Number, D::Diagonal)
    if size(D,1) > 1
        # ensure that zeros are preserved on scaling
        y = D[2,1] * x
        iszero(y) || throw(ArgumentError(LazyString("cannot set index (2, 1) off ",
            lazy"the tridiagonal band to a nonzero value ($y)")))
    end
    lmul!(x, D.diag)
    return D
end
function rmul!(D::Diagonal, x::Number)
    if size(D,1) > 1
        # ensure that zeros are preserved on scaling
        y = x * D[2,1]
        iszero(y) || throw(ArgumentError(LazyString("cannot set index (2, 1) off ",
            lazy"the tridiagonal band to a nonzero value ($y)")))
    end
    rmul!(D.diag, x)
    return D
end
(/)(D::Diagonal, x::Number) = Diagonal(D.diag / x)
(\)(x::Number, D::Diagonal) = Diagonal(x \ D.diag)
(^)(D::Diagonal, a::Number) = Diagonal(D.diag .^ a)
(^)(D::Diagonal, a::Real) = Diagonal(D.diag .^ a) # for disambiguation
(^)(D::Diagonal, a::Integer) = Diagonal(D.diag .^ a) # for disambiguation
Base.literal_pow(::typeof(^), D::Diagonal, valp::Val) =
    Diagonal(Base.literal_pow.(^, D.diag, valp)) # for speed
Base.literal_pow(::typeof(^), D::Diagonal, ::Val{-1}) = inv(D) # for disambiguation

function mul(Da::Diagonal, Db::Diagonal)
    matmul_size_check(size(Da), size(Db))
    return Diagonal(Da.diag .* Db.diag)
end

function (*)(D::Diagonal, V::AbstractVector)
    matmul_size_check(size(D), size(V))
    return D.diag .* V
end

function mul(A::AdjOrTransAbsMat, D::Diagonal)
    adj = _wrapperop(A)
    copy(adj(adj(D) * adj(A)))
end
function mul(A::AdjOrTransAbsMat{<:Number, <:StridedMatrix}, D::Diagonal{<:Number})
    @invoke mul(A::AbstractMatrix, D::AbstractMatrix)
end
function mul(D::Diagonal, A::AdjOrTransAbsMat)
    adj = _wrapperop(A)
    copy(adj(adj(A) * adj(D)))
end
function mul(D::Diagonal{<:Number}, A::AdjOrTransAbsMat{<:Number, <:StridedMatrix})
    @invoke mul(D::AbstractMatrix, A::AbstractMatrix)
end

function rmul!(A::AbstractMatrix, D::Diagonal)
    matmul_size_check(size(A), size(D))
    for I in CartesianIndices(A)
        row, col = Tuple(I)
        @inbounds A[row, col] *= D.diag[col]
    end
    return A
end
# A' = A' * D => A = D' * A
# This uses the fact that D' is a Diagonal
function rmul!(A::AdjOrTransAbsMat, D::Diagonal)
    f = _wrapperop(A)
    lmul!(f(D), f(A))
    A
end
# T .= T * D
function rmul!(T::Tridiagonal, D::Diagonal)
    matmul_size_check(size(T), size(D))
    (; dl, d, du) = T
    d[1] *= D.diag[1]
    for i in axes(dl,1)
        dl[i] *= D.diag[i]
        du[i] *= D.diag[i+1]
        d[i+1] *= D.diag[i+1]
    end
    return T
end
for T in [:UpperTriangular, :UnitUpperTriangular,
        :LowerTriangular, :UnitLowerTriangular]
    @eval rmul!(A::$T{<:Any, <:StridedMatrix}, D::Diagonal) = _rmul!(A, D)
    @eval lmul!(D::Diagonal, A::$T{<:Any, <:StridedMatrix}) = _lmul!(D, A)
end
function _rmul!(A::UpperOrLowerTriangular, D::Diagonal)
    P = parent(A)
    isunit = A isa UnitUpperOrUnitLowerTriangular
    isupper = A isa UpperOrUnitUpperTriangular
    for col in axes(A,2)
        rowstart = isupper ? firstindex(A,1) : col+isunit
        rowstop = isupper ? col-isunit : lastindex(A,1)
        for row in rowstart:rowstop
            P[row, col] *= D.diag[col]
        end
    end
    isunit && _setdiag!(P, identity, D.diag)
    TriWrapper = isupper ? UpperTriangular : LowerTriangular
    return TriWrapper(P)
end

function lmul!(D::Diagonal, B::AbstractVecOrMat)
    matmul_size_check(size(D), size(B))
    for I in CartesianIndices(B)
        row = I[1]
        @inbounds B[I] = D.diag[row] * B[I]
    end
    return B
end
# A' = D * A' => A = A * D'
# This uses the fact that D' is a Diagonal
function lmul!(D::Diagonal, A::AdjOrTransAbsMat)
    f = _wrapperop(A)
    rmul!(f(A), f(D))
    A
end

# in-place multiplication with a diagonal
# T .= D * T
function lmul!(D::Diagonal, T::Tridiagonal)
    matmul_size_check(size(D), size(T))
    (; dl, d, du) = T
    d[1] = D.diag[1] * d[1]
    for i in axes(dl,1)
        dl[i] = D.diag[i+1] * dl[i]
        du[i] = D.diag[i] * du[i]
        d[i+1] = D.diag[i+1] * d[i+1]
    end
    return T
end
function _lmul!(D::Diagonal, A::UpperOrLowerTriangular)
    P = parent(A)
    isunit = A isa UnitUpperOrUnitLowerTriangular
    isupper = A isa UpperOrUnitUpperTriangular
    for col in axes(A,2)
        rowstart = isupper ? firstindex(A,1) : col+isunit
        rowstop = isupper ? col-isunit : lastindex(A,1)
        for row in rowstart:rowstop
            P[row, col] = D.diag[row] * P[row, col]
        end
    end
    isunit && _setdiag!(P, identity, D.diag)
    TriWrapper = isupper ? UpperTriangular : LowerTriangular
    return TriWrapper(P)
end

@propagate_inbounds _modify_nonzeroalpha!(x, out, ind, alpha, beta) = @stable_muladdmul _modify!(MulAddMul(alpha,beta), x, out, ind)
@propagate_inbounds _modify_nonzeroalpha!(x, out, ind, ::Bool, beta) = @stable_muladdmul _modify!(MulAddMul(true,beta), x, out, ind)

@inline function __muldiag_nonzeroalpha!(out, D::Diagonal, B, alpha::Number, beta::Number)
    @inbounds for j in axes(B, 2)
        @simd for i in axes(B, 1)
            _modify_nonzeroalpha!(D.diag[i] * B[i,j], out, (i,j), alpha, beta)
        end
    end
    return out
end
_has_matching_zeros(out::UpperOrUnitUpperTriangular, A::UpperOrUnitUpperTriangular) = true
_has_matching_zeros(out::LowerOrUnitLowerTriangular, A::LowerOrUnitLowerTriangular) = true
_has_matching_zeros(out, A) = false
function _rowrange_tri_stored(B::UpperOrUnitUpperTriangular, col)
    isunit = B isa UnitUpperTriangular
    1:min(col-isunit, size(B,1))
end
function _rowrange_tri_stored(B::LowerOrUnitLowerTriangular, col)
    isunit = B isa UnitLowerTriangular
    col+isunit:size(B,1)
end
_rowrange_tri_zeros(B::UpperOrUnitUpperTriangular, col) = col+1:size(B,1)
_rowrange_tri_zeros(B::LowerOrUnitLowerTriangular, col) = 1:col-1
function __muldiag_nonzeroalpha!(out, D::Diagonal, B::UpperOrLowerTriangular, alpha::Number, beta::Number)
    isunit = B isa UnitUpperOrUnitLowerTriangular
    out_maybeparent, B_maybeparent = _has_matching_zeros(out, B) ? (parent(out), parent(B)) : (out, B)
    for j in axes(B, 2)
        # store the diagonal separately for unit triangular matrices
        if isunit
            @inbounds _modify_nonzeroalpha!(D.diag[j] * B[j,j], out, (j,j), alpha, beta)
        end
        # The indices of out corresponding to the stored indices of B
        rowrange = _rowrange_tri_stored(B, j)
        @inbounds @simd for i in rowrange
            _modify_nonzeroalpha!(D.diag[i] * B_maybeparent[i,j], out_maybeparent, (i,j), alpha, beta)
        end
        # Fill the indices of out corresponding to the zeros of B
        # we only fill these if out and B don't have matching zeros
        if !_has_matching_zeros(out, B)
            rowrange = _rowrange_tri_zeros(B, j)
            @inbounds @simd for i in rowrange
                _modify_nonzeroalpha!(D.diag[i] * B[i,j], out, (i,j), alpha, beta)
            end
        end
    end
    return out
end

@inline _djalpha_nonzero(dj, alpha) = @stable_muladdmul MulAddMul(alpha,false)(dj)
@inline _djalpha_nonzero(dj, ::Bool) = dj

@inline function __muldiag_nonzeroalpha_right!(out, A, D::Diagonal, alpha::Number, beta::Number)
    @inbounds for j in axes(A, 2)
        dja = _djalpha_nonzero(D.diag[j], alpha)
        @simd for i in axes(A, 1)
            @stable_muladdmul _modify!(MulAddMul(true,beta), A[i,j] * dja, out, (i,j))
        end
    end
    return out
end

function __muldiag_nonzeroalpha!(out, A, D::Diagonal, alpha::Number, beta::Number)
    __muldiag_nonzeroalpha_right!(out, A, D, alpha, beta)
end
function __muldiag_nonzeroalpha!(out, A::UpperOrLowerTriangular, D::Diagonal, alpha::Number, beta::Number)
    isunit = A isa UnitUpperOrUnitLowerTriangular
    # if both A and out have the same upper/lower triangular structure,
    # we may directly read and write from the parents
    out_maybeparent, A_maybeparent = _has_matching_zeros(out, A) ? (parent(out), parent(A)) : (out, A)
    for j in axes(A, 2)
        dja = @inbounds _djalpha_nonzero(D.diag[j], alpha)
        # store the diagonal separately for unit triangular matrices
        if isunit
            # since alpha is multiplied to the diagonal element of D,
            # we may skip alpha in the second multiplication by setting ais1 to true
            @inbounds @stable_muladdmul _modify!(MulAddMul(true,beta), A[j,j] * dja, out, (j,j))
        end
        # indices of out corresponding to the stored indices of A
        rowrange = _rowrange_tri_stored(A, j)
        @inbounds @simd for i in rowrange
            # since alpha is multiplied to the diagonal element of D,
            # we may skip alpha in the second multiplication by setting ais1 to true
            @stable_muladdmul _modify!(MulAddMul(true,beta), A_maybeparent[i,j] * dja, out_maybeparent, (i,j))
        end
        # Fill the indices of out corresponding to the zeros of A
        # we only fill these if out and A don't have matching zeros
        if !_has_matching_zeros(out, A)
            rowrange = _rowrange_tri_zeros(A, j)
            @inbounds @simd for i in rowrange
                @stable_muladdmul _modify!(MulAddMul(true,beta), A[i,j] * dja, out, (i,j))
            end
        end
    end
    return out
end

# ambiguity resolution
function __muldiag_nonzeroalpha!(out, D1::Diagonal, D2::Diagonal, alpha::Number, beta::Number)
    __muldiag_nonzeroalpha_right!(out, D1, D2, alpha, beta)
end

@inline function __muldiag_nonzeroalpha!(out::Diagonal, D1::Diagonal, D2::Diagonal, alpha::Number, beta::Number)
    d1 = D1.diag
    d2 = D2.diag
    outd = out.diag
    @inbounds @simd for i in eachindex(d1, d2, outd)
        _modify_nonzeroalpha!(d1[i] * d2[i], outd, i, alpha, beta)
    end
    return out
end

# muldiag handles the zero-alpha case, so that we need only
# specialize the non-trivial case
function _mul_diag!(out, A, B, alpha, beta)
    require_one_based_indexing(out, A, B)
    matmul_size_check(size(out), size(A), size(B))
    if iszero(alpha)
        _rmul_or_fill!(out, beta)
    else
        __muldiag_nonzeroalpha!(out, A, B, alpha, beta)
    end
    return out
end

_mul!(out::AbstractVector, D::Diagonal, V::AbstractVector, alpha::Number, beta::Number) =
    _mul_diag!(out, D, V, alpha, beta)
_mul!(out::AbstractMatrix, D::Diagonal, V::AbstractVector, alpha::Number, beta::Number) =
    _mul_diag!(out, D, V, alpha, beta)
for MT in (:AbstractMatrix, :AbstractTriangular)
    @eval begin
        _mul!(out::AbstractMatrix, D::Diagonal, B::$MT, alpha::Number, beta::Number) =
            _mul_diag!(out, D, B, alpha, beta)
        _mul!(out::AbstractMatrix, A::$MT, D::Diagonal, alpha::Number, beta::Number) =
            _mul_diag!(out, A, D, alpha, beta)
    end
end
_mul!(C::AbstractMatrix, Da::Diagonal, Db::Diagonal, alpha::Number, beta::Number) =
    _mul_diag!(C, Da, Db, alpha, beta)

function (*)(Da::Diagonal, A::AbstractMatrix, Db::Diagonal)
    matmul_size_check(size(Da), size(A))
    matmul_size_check(size(A), size(Db))
    return broadcast(*, Da.diag, A, permutedims(Db.diag))
end

function (*)(Da::Diagonal, Db::Diagonal, Dc::Diagonal)
    matmul_size_check(size(Da), size(Db))
    matmul_size_check(size(Db), size(Dc))
    return Diagonal(Da.diag .* Db.diag .* Dc.diag)
end

/(A::AbstractVecOrMat, D::Diagonal) = _rdiv!(matprod_dest(A, D, promote_op(/, eltype(A), eltype(D))), A, D)

rdiv!(A::AbstractVecOrMat, D::Diagonal) = @inline _rdiv!(A, A, D)
# avoid copy when possible via internal 3-arg backend
function _rdiv!(B::AbstractVecOrMat, A::AbstractVecOrMat, D::Diagonal)
    require_one_based_indexing(A)
    dd = D.diag
    m, n = size(A, 1), size(A, 2)
    if (k = length(dd)) != n
        throw(DimensionMismatch(lazy"left hand side has $n columns but D is $k by $k"))
    end
    @inbounds for j in 1:n
        ddj = dd[j]
        iszero(ddj) && throw(SingularException(j))
        for i in 1:m
            B[i, j] = A[i, j] / ddj
        end
    end
    B
end

function \(D::Diagonal, B::AbstractVector)
    j = findfirst(iszero, D.diag)
    isnothing(j) || throw(SingularException(j))
    return D.diag .\ B
end
\(D::Diagonal, B::AbstractMatrix) = ldiv!(matprod_dest(D, B, promote_op(\, eltype(D), eltype(B))), D, B)

ldiv!(D::Diagonal, B::AbstractVecOrMat) = @inline ldiv!(B, D, B)
function ldiv!(B::AbstractVecOrMat, D::Diagonal, A::AbstractVecOrMat)
    require_one_based_indexing(A, B)
    dd = D.diag
    d = length(dd)
    m, n = size(A, 1), size(A, 2)
    m′, n′ = size(B, 1), size(B, 2)
    m == d || throw(DimensionMismatch(lazy"right hand side has $m rows but D is $d by $d"))
    (m, n) == (m′, n′) || throw(DimensionMismatch(lazy"expect output to be $m by $n, but got $m′ by $n′"))
    j = findfirst(iszero, D.diag)
    isnothing(j) || throw(SingularException(j))
    @inbounds for j = 1:n, i = 1:m
        B[i, j] = dd[i] \ A[i, j]
    end
    B
end

function _rdiv!(Dc::Diagonal, Db::Diagonal, Da::Diagonal)
    n, k = length(Db.diag), length(Da.diag)
    n == k || throw(DimensionMismatch(lazy"left hand side has $n columns but D is $k by $k"))
    j = findfirst(iszero, Da.diag)
    isnothing(j) || throw(SingularException(j))
    Dc.diag .= Db.diag ./ Da.diag
    Dc
end
ldiv!(Dc::Diagonal, Da::Diagonal, Db::Diagonal) = Diagonal(ldiv!(Dc.diag, Da, Db.diag))

# optimizations for (Sym)Tridiagonal and Diagonal
@propagate_inbounds _getudiag(T::Tridiagonal, i) = T.du[i]
@propagate_inbounds _getudiag(S::SymTridiagonal, i) = S.ev[i]
@propagate_inbounds _getdiag(T::Tridiagonal, i) = T.d[i]
@propagate_inbounds _getdiag(S::SymTridiagonal, i) = symmetric(S.dv[i], :U)::symmetric_type(eltype(S.dv))
@propagate_inbounds _getldiag(T::Tridiagonal, i) = T.dl[i]
@propagate_inbounds _getldiag(S::SymTridiagonal, i) = transpose(S.ev[i])

function (\)(D::Diagonal, S::SymTridiagonal)
    T = promote_op(\, eltype(D), eltype(S))
    du = similar(S.ev, T, max(length(S.dv)-1, 0))
    d  = similar(S.dv, T, length(S.dv))
    dl = similar(S.ev, T, max(length(S.dv)-1, 0))
    ldiv!(Tridiagonal(dl, d, du), D, S)
end
(\)(D::Diagonal, T::Tridiagonal) = ldiv!(similar(T, promote_op(\, eltype(D), eltype(T))), D, T)
function ldiv!(T::Tridiagonal, D::Diagonal, S::Union{SymTridiagonal,Tridiagonal})
    m = size(S, 1)
    dd = D.diag
    if (k = length(dd)) != m
        throw(DimensionMismatch(lazy"diagonal matrix is $k by $k but right hand side has $m rows"))
    end
    if length(T.d) != m
        throw(DimensionMismatch(lazy"target matrix size $(size(T)) does not match input matrix size $(size(S))"))
    end
    m == 0 && return T
    j = findfirst(iszero, dd)
    isnothing(j) || throw(SingularException(j))
    ddj = dd[1]
    T.d[1] = ddj \ _getdiag(S, 1)
    @inbounds if m > 1
        T.du[1] = ddj \ _getudiag(S, 1)
        for j in 2:m-1
            ddj = dd[j]
            T.dl[j-1] = ddj \ _getldiag(S, j-1)
            T.d[j]  = ddj \ _getdiag(S, j)
            T.du[j] = ddj \ _getudiag(S, j)
        end
        ddj = dd[m]
        T.dl[m-1] = ddj \ _getldiag(S, m-1)
        T.d[m] = ddj \ _getdiag(S, m)
    end
    return T
end

function (/)(S::SymTridiagonal, D::Diagonal)
    T = promote_op(\, eltype(D), eltype(S))
    du = similar(S.ev, T, max(length(S.dv)-1, 0))
    d  = similar(S.dv, T, length(S.dv))
    dl = similar(S.ev, T, max(length(S.dv)-1, 0))
    _rdiv!(Tridiagonal(dl, d, du), S, D)
end
(/)(T::Tridiagonal, D::Diagonal) = _rdiv!(matprod_dest(T, D, promote_op(/, eltype(T), eltype(D))), T, D)
function _rdiv!(T::Tridiagonal, S::Union{SymTridiagonal,Tridiagonal}, D::Diagonal)
    n = size(S, 2)
    dd = D.diag
    if (k = length(dd)) != n
        throw(DimensionMismatch(lazy"left hand side has $n columns but D is $k by $k"))
    end
    if length(T.d) != n
        throw(DimensionMismatch(lazy"target matrix size $(size(T)) does not match input matrix size $(size(S))"))
    end
    n == 0 && return T
    j = findfirst(iszero, dd)
    isnothing(j) || throw(SingularException(j))
    ddj = dd[1]
    T.d[1] = _getdiag(S, 1) / ddj
    @inbounds if n > 1
        T.dl[1] = _getldiag(S, 1) / ddj
        for j in 2:n-1
            ddj = dd[j]
            T.dl[j] = _getldiag(S, j) / ddj
            T.d[j] = _getdiag(S, j) / ddj
            T.du[j-1] = _getudiag(S, j-1) / ddj
        end
        ddj = dd[n]
        T.d[n] = _getdiag(S, n) / ddj
        T.du[n-1] = _getudiag(S, n-1) / ddj
    end
    return T
end

# Optimizations for [l/r]mul!, l/rdiv!, *, / and \ between Triangular and Diagonal.
# These functions are generally more efficient if we calculate the whole data field.
# The following code implements them in a unified pattern to avoid missing.
@inline function _setdiag!(data, f, diag, diag′ = nothing)
    @inbounds for i in 1:length(diag)
        data[i,i] = isnothing(diag′) ? f(diag[i]) : f(diag[i],diag′[i])
    end
    data
end
for Tri in (:UpperTriangular, :LowerTriangular)
    UTri = Symbol(:Unit, Tri)
    # 2 args
    for (fun, f) in zip((:mul, :rmul!, :rdiv!, :/), (:identity, :identity, :inv, :inv))
        g = fun == :mul ? :* : fun
        @eval $fun(A::$Tri, D::Diagonal) = $Tri($g(A.data, D))
        @eval $fun(A::$UTri, D::Diagonal) = $Tri(_setdiag!($g(A.data, D), $f, D.diag))
    end
    @eval mul(A::$Tri{<:Any, <:StridedMaybeAdjOrTransMat}, D::Diagonal) =
            @invoke mul(A::AbstractMatrix, D::Diagonal)
    @eval mul(A::$UTri{<:Any, <:StridedMaybeAdjOrTransMat}, D::Diagonal) =
            @invoke mul(A::AbstractMatrix, D::Diagonal)
    for (fun, f) in zip((:mul, :lmul!, :ldiv!, :\), (:identity, :identity, :inv, :inv))
        g = fun == :mul ? :* : fun
        @eval $fun(D::Diagonal, A::$Tri) = $Tri($g(D, A.data))
        @eval $fun(D::Diagonal, A::$UTri) = $Tri(_setdiag!($g(D, A.data), $f, D.diag))
    end
    @eval mul(D::Diagonal, A::$Tri{<:Any, <:StridedMaybeAdjOrTransMat}) =
            @invoke mul(D::Diagonal, A::AbstractMatrix)
    @eval mul(D::Diagonal, A::$UTri{<:Any, <:StridedMaybeAdjOrTransMat}) =
            @invoke mul(D::Diagonal, A::AbstractMatrix)
    # 3-arg ldiv!
    @eval ldiv!(C::$Tri, D::Diagonal, A::$Tri) = $Tri(ldiv!(C.data, D, A.data))
    @eval ldiv!(C::$Tri, D::Diagonal, A::$UTri) = $Tri(_setdiag!(ldiv!(C.data, D, A.data), inv, D.diag))
end

@inline function kron!(C::AbstractMatrix, A::Diagonal, B::Diagonal)
    valA = A.diag; mA, nA = size(A)
    valB = B.diag; mB, nB = size(B)
    nC = checksquare(C)
    @boundscheck nC == nA*nB ||
        throw(DimensionMismatch(lazy"expect C to be a $(nA*nB)x$(nA*nB) matrix, got size $(nC)x$(nC)"))
    zerofilled = false
    if !(isempty(A) || isempty(B))
        z = A[1,1] * B[1,1]
        if haszero(typeof(z))
            # in this case, the zero is unique
            fill!(C, zero(z))
            zerofilled = true
        end
    end
    for i in eachindex(valA), j in eachindex(valB)
        idx = (i-1)*nB+j
        @inbounds C[idx, idx] = valA[i] * valB[j]
    end
    if !zerofilled
        for j in axes(A,2), i in axes(A,1)
            Δrow, Δcol = (i-1)*mB, (j-1)*nB
            for k in axes(B,2), l in axes(B,1)
                i == j && k == l && continue
                @inbounds C[Δrow + l, Δcol + k] = A[i,j] * B[l,k]
            end
        end
    end
    return C
end

kron(A::Diagonal, B::Diagonal) = Diagonal(kron(A.diag, B.diag))

function kron!(C::Diagonal, A::Diagonal, B::Diagonal)
    kron!(C.diag, A.diag, B.diag)
    return C
end

function kron(A::Diagonal, B::SymTridiagonal)
    kdv = kron(A.diag, B.dv)
    # We don't need to drop the last element
    kev = kron(A.diag, _pushzero(_evview(B)))
    SymTridiagonal(kdv, kev)
end
function kron(A::Diagonal, B::Tridiagonal)
    # `_droplast!` is only guaranteed to work with `Vector`
    kd = convert(Vector, kron(A.diag, B.d))
    kdl = _droplast!(convert(Vector, kron(A.diag, _pushzero(B.dl))))
    kdu = _droplast!(convert(Vector, kron(A.diag, _pushzero(B.du))))
    Tridiagonal(kdl, kd, kdu)
end

@inline function kron!(C::AbstractMatrix, A::Diagonal, B::AbstractMatrix)
    require_one_based_indexing(B)
    (mA, nA) = size(A)
    (mB, nB) = size(B)
    (mC, nC) = size(C)
    @boundscheck (mC, nC) == (mA * mB, nA * nB) ||
        throw(DimensionMismatch(lazy"expect C to be a $(mA * mB)x$(nA * nB) matrix, got size $(mC)x$(nC)"))
    zerofilled = false
    if !(isempty(A) || isempty(B))
        z = A[1,1] * B[1,1]
        if haszero(typeof(z))
            # in this case, the zero is unique
            fill!(C, zero(z))
            zerofilled = true
        end
    end
    m = 1
    for j in axes(A,2)
        A_jj = @inbounds A[j,j]
        for k in axes(B,2)
            for l in axes(B,1)
                @inbounds C[m] = A_jj * B[l,k]
                m += 1
            end
            m += (nA - 1) * mB
        end
        if !zerofilled
            # populate the zero elements
            for i in axes(A,1)
                i == j && continue
                A_ij = @inbounds A[i, j]
                Δrow, Δcol = (i-1)*mB, (j-1)*nB
                for k in axes(B,2), l in axes(B,1)
                    B_lk = @inbounds B[l, k]
                    @inbounds C[Δrow + l, Δcol + k] = A_ij * B_lk
                end
            end
        end
        m += mB
    end
    return C
end

@inline function kron!(C::AbstractMatrix, A::AbstractMatrix, B::Diagonal)
    require_one_based_indexing(A)
    (mA, nA) = size(A)
    (mB, nB) = size(B)
    (mC, nC) = size(C)
    @boundscheck (mC, nC) == (mA * mB, nA * nB) ||
        throw(DimensionMismatch(lazy"expect C to be a $(mA * mB)x$(nA * nB) matrix, got size $(mC)x$(nC)"))
    zerofilled = false
    if !(isempty(A) || isempty(B))
        z = A[1,1] * B[1,1]
        if haszero(typeof(z))
            # in this case, the zero is unique
            fill!(C, zero(z))
            zerofilled = true
        end
    end
    m = 1
    for j in axes(A,2)
        for l in axes(B,1)
            Bll = @inbounds B[l,l]
            for i in axes(A,1)
                @inbounds C[m] = A[i,j] * Bll
                m += nB
            end
            m += 1
        end
        if !zerofilled
            for i in axes(A,1)
                A_ij = @inbounds A[i, j]
                Δrow, Δcol = (i-1)*mB, (j-1)*nB
                for k in axes(B,2), l in axes(B,1)
                    l == k && continue
                    B_lk = @inbounds B[l, k]
                    @inbounds C[Δrow + l, Δcol + k] = A_ij * B_lk
                end
            end
        end
        m -= nB
    end
    return C
end

conj(D::Diagonal) = Diagonal(conj(D.diag))
transpose(D::Diagonal) = Diagonal(_vectranspose(D.diag))
adjoint(D::Diagonal) = Diagonal(_vecadjoint(D.diag))
permutedims(D::Diagonal) = D
permutedims(D::Diagonal, perm) = (Base.checkdims_perm(axes(D), axes(D), perm); D)

function diag(D::Diagonal, k::Integer=0)
    # every branch call similar(..., ::Int) to make sure the
    # same vector type is returned independent of k
    dinds = diagind(D, k, IndexStyle(D))
    v = similar(D.diag, length(dinds))
    if k == 0
        copyto!(v, D.diag)
    else
        for i in eachindex(v, dinds)
            @inbounds v[i] = D[dinds[i]]
        end
    end
    return v
end
tr(D::Diagonal) = sum(tr, D.diag)
det(D::Diagonal) = prod(det, D.diag)
function logdet(D::Diagonal{<:Complex}) # make sure branch cut is correct
    z = sum(log, D.diag)
    complex(real(z), rem2pi(imag(z), RoundNearest))
end

# Matrix functions
for f in (:exp, :cis, :log, :sqrt,
          :cos, :sin, :tan, :csc, :sec, :cot,
          :cosh, :sinh, :tanh, :csch, :sech, :coth,
          :acos, :asin, :atan, :acsc, :asec, :acot,
          :acosh, :asinh, :atanh, :acsch, :asech, :acoth)
    @eval $f(D::Diagonal) = Diagonal($f.(D.diag))
end

# Cube root of a real-valued diagonal matrix
cbrt(A::Diagonal{<:Real}) = Diagonal(cbrt.(A.diag))

function inv(D::Diagonal{T}) where T
    Di = similar(D.diag, typeof(inv(oneunit(T))))
    for i = 1:length(D.diag)
        if iszero(D.diag[i])
            throw(SingularException(i))
        end
        Di[i] = inv(D.diag[i])
    end
    Diagonal(Di)
end

function pinv(D::Diagonal{T}) where T
    Di = similar(D.diag, typeof(inv(oneunit(T))))
    for i = 1:length(D.diag)
        if !iszero(D.diag[i])
            invD = inv(D.diag[i])
            if isfinite(invD)
                Di[i] = invD
                continue
            end
        end
        # fallback
        Di[i] = zero(T)
    end
    Diagonal(Di)
end
function pinv(D::Diagonal{T}, tol::Real) where T
    Di = similar(D.diag, typeof(inv(oneunit(T))))
    if !isempty(D.diag)
        maxabsD = maximum(abs, D.diag)
        for i = 1:length(D.diag)
            if abs(D.diag[i]) > tol*maxabsD
                invD = inv(D.diag[i])
                if isfinite(invD)
                    Di[i] = invD
                    continue
                end
            end
            # fallback
            Di[i] = zero(T)
        end
    end
    Diagonal(Di)
end

_ortho_eltype(T) = Base.promote_op(/, T, T)
_ortho_eltype(T::Type{<:Number}) = typeof(one(T)/one(T))

# TODO Docstrings for eigvals, eigvecs, eigen all mention permute, scale, sortby as keyword args
# but not all of them below provide them. Do we need to fix that?
#Eigensystem
eigvals(D::Diagonal{<:Number}; permute::Bool=true, scale::Bool=true) = copy(D.diag)
eigvals(D::Diagonal; permute::Bool=true, scale::Bool=true) =
    reduce(vcat, eigvals(x) for x in D.diag) #For block matrices, etc.
function eigvecs(D::Diagonal{T}) where {T<:AbstractMatrix}
    diag_vecs = [ eigvecs(x) for x in D.diag ]
    matT = promote_type(map(typeof, diag_vecs)...)
    ncols_diag = [ size(x, 2) for x in D.diag ]
    nrows = size(D, 1)
    vecs = Matrix{Vector{eltype(matT)}}(undef, nrows, sum(ncols_diag))
    for j in axes(D, 2), i in axes(D, 1)
        jj = sum(view(ncols_diag,1:j-1))
        if i == j
            for k in 1:ncols_diag[j]
                vecs[i,jj+k] = diag_vecs[i][:,k]
            end
        else
            for k in 1:ncols_diag[j]
                vecs[i,jj+k] = zeros(eltype(T), ncols_diag[i])
            end
        end
    end
    return vecs
end
function eigen(D::Diagonal; permute::Bool=true, scale::Bool=true, sortby::Union{Function,Nothing}=nothing)
    if any(!isfinite, D.diag)
        throw(ArgumentError("matrix contains Infs or NaNs"))
    end
    Td = _ortho_eltype(eltype(D))
    λ = eigvals(D)
    if !isnothing(sortby)
        p = sortperm(λ; alg=QuickSort, by=sortby)
        λ = λ[p]
        evecs = zeros(Td, size(D))
        @inbounds for i in eachindex(p)
            evecs[p[i],i] = one(Td)
        end
    else
        evecs = Diagonal(ones(Td, length(λ)))
    end
    Eigen(λ, evecs)
end
function eigen(D::Diagonal{<:AbstractMatrix}; permute::Bool=true, scale::Bool=true, sortby::Union{Function,Nothing}=nothing)
    if any(any(!isfinite, x) for x in D.diag)
        throw(ArgumentError("matrix contains Infs or NaNs"))
    end
    λ = eigvals(D)
    evecs = eigvecs(D)
    if !isnothing(sortby)
        p = sortperm(λ; alg=QuickSort, by=sortby)
        λ = λ[p]
        evecs = evecs[:,p]
    end
    Eigen(λ, evecs)
end
function eigen(Da::Diagonal, Db::Diagonal; sortby::Union{Function,Nothing}=nothing)
    if any(!isfinite, Da.diag) || any(!isfinite, Db.diag)
        throw(ArgumentError("matrices contain Infs or NaNs"))
    end
    if any(iszero, Db.diag)
        throw(ArgumentError("right-hand side diagonal matrix is singular"))
    end
    return GeneralizedEigen(eigen(Db \ Da; sortby)...)
end
function eigen(A::AbstractMatrix, D::Diagonal; sortby::Union{Function,Nothing}=nothing)
    if any(iszero, D.diag)
        throw(ArgumentError("right-hand side diagonal matrix is singular"))
    end
    if size(A, 1) == size(A, 2) && isdiag(A)
        return eigen(Diagonal(A), D; sortby)
    elseif all(isposdef, D.diag)
        S = promote_type(eigtype(eltype(A)), eltype(D))
        return eigen(A, cholesky(Diagonal{S}(D)); sortby)
    else
        return eigen!(D \ A; sortby)
    end
end

#Singular system
svdvals(D::Diagonal{<:Number}) = sort!(abs.(D.diag), rev = true)
svdvals(D::Diagonal) = [svdvals(v) for v in D.diag]
function svd(D::Diagonal{T}) where {T<:Number}
    d = D.diag
    s = abs.(d)
    piv = sortperm(s, rev = true)
    S = s[piv]
    Td  = _ortho_eltype(T)
    U = zeros(Td, size(D))
    Vt = copy(U)
    for i in 1:length(d)
        j = piv[i]
        U[j,i] = iszero(d[j]) ? one(Td) : d[j] / S[i]
        Vt[i,j] = one(Td)
    end
    return SVD(U, S, Vt)
end

*(x::AdjointAbsVec,   D::Diagonal, y::AbstractVector) = _mapreduce_prod(*, x, D, y)
*(x::TransposeAbsVec, D::Diagonal, y::AbstractVector) = _mapreduce_prod(*, x, D, y)
/(u::AdjointAbsVec, D::Diagonal) = (D' \ u')'
/(u::TransposeAbsVec, D::Diagonal) = transpose(transpose(D) \ transpose(u))

_opnorm1(A::Diagonal) = maximum(norm(x) for x in A.diag)
_opnormInf(A::Diagonal) = maximum(norm(x) for x in A.diag)
_opnorm12Inf(A::Diagonal, p) = maximum(opnorm(x, p) for x in A.diag)

function opnorm(A::Diagonal, p::Real=2)
    if p == 1 || p == Inf || p == 2
        return _opnorm12Inf(A, p)
    else
        throw(ArgumentError("invalid p-norm p=$p. Valid: 1, 2, Inf"))
    end
end

dot(x::AbstractVector, D::Diagonal, y::AbstractVector) = _mapreduce_prod(dot, x, D, y)

dot(A::Diagonal, B::Diagonal) = dot(A.diag, B.diag)
function dot(D::Diagonal, B::AbstractMatrix)
    size(D) == size(B) || throw(DimensionMismatch(lazy"Matrix sizes $(size(D)) and $(size(B)) differ"))
    return dot(D.diag, diagview(B))
end

dot(A::AbstractMatrix, B::Diagonal) = conj(dot(B, A))

function _mapreduce_prod(f, x, D::Diagonal, y)
    if !(length(x) == length(D.diag) == length(y))
        throw(DimensionMismatch(lazy"x has length $(length(x)), D has size $(size(D)), and y has $(length(y))"))
    end
    if isempty(x) && isempty(D) && isempty(y)
        return zero(promote_op(f, eltype(x), eltype(D), eltype(y)))
    else
        return mapreduce(t -> f(t[1], t[2], t[3]), +, zip(x, D.diag, y))
    end
end

function cholesky!(A::Diagonal, ::NoPivot = NoPivot(); check::Bool = true)
    info = 0
    for (i, di) in enumerate(A.diag)
        if isreal(di) && real(di) > 0
            A.diag[i] = √di
        elseif check
            throw(PosDefException(i))
        else
            info = i
            break
        end
    end
    Cholesky(A, 'U', convert(BlasInt, info))
end
@deprecate cholesky!(A::Diagonal, ::Val{false}; check::Bool = true) cholesky!(A::Diagonal, NoPivot(); check) false
@deprecate cholesky(A::Diagonal, ::Val{false}; check::Bool = true) cholesky(A::Diagonal, NoPivot(); check) false

function cholesky!(A::Diagonal, ::RowMaximum; tol=0.0, check=true)
    if !ishermitian(A)
        C = CholeskyPivoted(A, 'U', Vector{BlasInt}(), convert(BlasInt, 1),
                            tol, convert(BlasInt, -1))
        check && checkpositivedefinite(convert(BlasInt, -1))
    else
        d = A.diag
        n = length(d)
        info = 0
        rank = n
        p = sortperm(d, rev = true, by = real)
        tol = tol < 0 ? n*eps(eltype(A))*real(d[p[1]]) : tol # LAPACK behavior
        permute!(d, p)
        @inbounds for i in eachindex(d)
            di = d[i]
            rootdi, j = _cholpivoted!(di, tol)
            if j == 0
                d[i] = rootdi
            else
                rank = i - 1
                info = 1
                break
            end
        end
        C = CholeskyPivoted(A, 'U', p, convert(BlasInt, rank), tol, convert(BlasInt, info))
        check && chkfullrank(C)
    end
    return C
end

inv(C::Cholesky{<:Any,<:Diagonal}) = Diagonal(map(inv∘abs2, C.factors.diag))

cholcopy(A::Diagonal) = copymutable_oftype(A, choltype(A))
cholcopy(A::RealHermSymComplexHerm{<:Any,<:Diagonal}) = Diagonal(copy_similar(diag(A), choltype(A)))

function getproperty(C::Cholesky{<:Any,<:Diagonal}, d::Symbol)
    Cfactors = getfield(C, :factors)
    if d in (:U, :L, :UL)
        return Cfactors
    else
        return getfield(C, d)
    end
end

Base._sum(A::Diagonal, ::Colon) = sum(A.diag)
function Base._sum(A::Diagonal, dims::Integer)
    res = Base.reducedim_initarray(A, dims, zero(eltype(A)))
    if dims <= 2
        for i = 1:length(A.diag)
            @inbounds res[i] = A.diag[i]
        end
    else
        for i = 1:length(A.diag)
            @inbounds res[i,i] = A.diag[i]
        end
    end
    res
end

function logabsdet(A::Diagonal)
     mapreduce(x -> (log(abs(x)), sign(x)), ((d1, s1), (d2, s2)) -> (d1 + d2, s1 * s2),
               A.diag)
end

function Base.muladd(A::Diagonal, B::Diagonal, z::Diagonal)
    Diagonal(A.diag .* B.diag .+ z.diag)
end

uppertriangular(D::Diagonal) = D
lowertriangular(D::Diagonal) = D

throw_fillband_error(l, u, x) = throw(ArgumentError(lazy"cannot set bands $l:$u to a nonzero value ($x)"))

function fillband!(D::Diagonal, x, l, u)
    if l > u
        return D
    end
    if (l < 0 || u > 0) && !iszero(x)
        throw_fillband_error(l, u, x)
    end
    if l <= 0 <= u
        fill!(D.diag, x)
    end
    return D
end
