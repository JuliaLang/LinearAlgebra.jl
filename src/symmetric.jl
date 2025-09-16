# This file is a part of Julia. License is MIT: https://julialang.org/license

# Symmetric and Hermitian matrices
struct Symmetric{T,S<:AbstractMatrix{<:T}} <: AbstractMatrix{T}
    data::S
    uplo::Char

    function Symmetric{T,S}(data, uplo::Char) where {T,S<:AbstractMatrix{<:T}}
        require_one_based_indexing(data)
        (uplo != 'U' && uplo != 'L') && throw_uplo()
        new{T,S}(data, uplo)
    end
end
"""
    Symmetric(A::AbstractMatrix, uplo::Symbol=:U)

Construct a `Symmetric` view of the upper (if `uplo = :U`) or lower (if `uplo = :L`)
triangle of the matrix `A`.

`Symmetric` views are mainly useful for real-symmetric matrices, for which
specialized algorithms (e.g. for eigenproblems) are enabled for `Symmetric` types.
More generally, see also [`Hermitian(A)`](@ref) for Hermitian matrices `A == A'`, which
is effectively equivalent to `Symmetric` for real matrices but is also useful for
complex matrices.  (Whereas complex `Symmetric` matrices are supported but have few
if any specialized algorithms.)

To compute the symmetric part of a real matrix, or more generally the Hermitian part `(A + A') / 2` of
a real or complex matrix `A`, use [`hermitianpart`](@ref).

The `uplo` symbol corresponding to the triangular half of `A` that is shared by the symmetric view may be
fetched by using the function [`LinearAlgebra.uplo`](@ref). The underlying matrix `A` may be fetched from the symmetric
view by using `parent`.

# Examples
```jldoctest
julia> A = [1 2 3; 4 5 6; 7 8 9]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> Supper = Symmetric(A)
3×3 Symmetric{Int64, Matrix{Int64}}:
 1  2  3
 2  5  6
 3  6  9

julia> Slower = Symmetric(A, :L)
3×3 Symmetric{Int64, Matrix{Int64}}:
 1  4  7
 4  5  8
 7  8  9

julia> hermitianpart(A)
3×3 Hermitian{Float64, Matrix{Float64}}:
 1.0  3.0  5.0
 3.0  5.0  7.0
 5.0  7.0  9.0
```

Note that `Supper` will not be equal to `Slower` unless `A` is itself symmetric (e.g. if
`A == transpose(A)`).
"""
function Symmetric(A::AbstractMatrix, uplo::Symbol=:U)
    checksquare(A)
    return symmetric_type(typeof(A))(A, char_uplo(uplo))
end

"""
    symmetric(A, uplo::Symbol=:U)

Construct a symmetric view of `A`. If `A` is a matrix, `uplo` controls whether the upper
(if `uplo = :U`) or lower (if `uplo = :L`) triangle of `A` is used to implicitly fill the
other one. If `A` is a `Number`, it is returned as is.

If a symmetric view of a matrix is to be constructed of which the elements are neither
matrices nor numbers, an appropriate method of `symmetric` has to be implemented. In that
case, `symmetric_type` has to be implemented, too.
"""
symmetric(A::AbstractMatrix, uplo::Symbol=:U) = Symmetric(A, uplo)
symmetric(A::Number, ::Symbol=:U) = A

"""
    symmetric_type(T::Type)

The type of the object returned by `symmetric(::T, ::Symbol)`. For matrices, this is an
appropriately typed `Symmetric`, for `Number`s, it is the original type. If `symmetric` is
implemented for a custom type, so should be `symmetric_type`, and vice versa.
"""
function symmetric_type(::Type{T}) where {S, T<:AbstractMatrix{S}}
    return Symmetric{Union{S, promote_op(transpose, S), symmetric_type(S)}, T}
end
function symmetric_type(::Type{T}) where {S<:Number, T<:AbstractMatrix{S}}
    return Symmetric{S, T}
end
function symmetric_type(::Type{T}) where {S<:AbstractMatrix, T<:AbstractMatrix{S}}
    return Symmetric{AbstractMatrix, T}
end
symmetric_type(::Type{T}) where {T<:Number} = T

struct Hermitian{T,S<:AbstractMatrix{<:T}} <: AbstractMatrix{T}
    data::S
    uplo::Char

    function Hermitian{T,S}(data, uplo::Char) where {T,S<:AbstractMatrix{<:T}}
        require_one_based_indexing(data)
        (uplo != 'U' && uplo != 'L') && throw_uplo()
        new{T,S}(data, uplo)
    end
end
"""
    Hermitian(A::AbstractMatrix, uplo::Symbol=:U)

Construct a `Hermitian` view of the upper (if `uplo = :U`) or lower (if `uplo = :L`)
triangle of the matrix `A`.

To compute the Hermitian part of `A`, use [`hermitianpart`](@ref).

The `uplo` symbol corresponding to the triangular half of `A` that is shared by the hermitian view may be
fetched by using the function [`LinearAlgebra.uplo`](@ref). The underlying matrix `A` may be fetched from the hermitian
view by using `parent`.

# Examples
```jldoctest
julia> A = [1 2+2im 3-3im; 4 5 6-6im; 7 8+8im 9]
3×3 Matrix{Complex{Int64}}:
 1+0im  2+2im  3-3im
 4+0im  5+0im  6-6im
 7+0im  8+8im  9+0im

julia> Hupper = Hermitian(A)
3×3 Hermitian{Complex{Int64}, Matrix{Complex{Int64}}}:
 1+0im  2+2im  3-3im
 2-2im  5+0im  6-6im
 3+3im  6+6im  9+0im

julia> Hlower = Hermitian(A, :L)
3×3 Hermitian{Complex{Int64}, Matrix{Complex{Int64}}}:
 1+0im  4+0im  7+0im
 4+0im  5+0im  8-8im
 7+0im  8+8im  9+0im

julia> hermitianpart(A)
3×3 Hermitian{ComplexF64, Matrix{ComplexF64}}:
 1.0+0.0im  3.0+1.0im  5.0-1.5im
 3.0-1.0im  5.0+0.0im  7.0-7.0im
 5.0+1.5im  7.0+7.0im  9.0+0.0im
```

Note that `Hupper` will not be equal to `Hlower` unless `A` is itself Hermitian (e.g. if `A == adjoint(A)`).

All non-real parts of the diagonal will be ignored.

```julia
Hermitian(fill(complex(1,1), 1, 1)) == fill(1, 1, 1)
```
"""
function Hermitian(A::AbstractMatrix, uplo::Symbol=:U)
    n = checksquare(A)
    return hermitian_type(typeof(A))(A, char_uplo(uplo))
end

"""
    hermitian(A, uplo::Symbol=:U)

Construct a hermitian view of `A`. If `A` is a matrix, `uplo` controls whether the upper
(if `uplo = :U`) or lower (if `uplo = :L`) triangle of `A` is used to implicitly fill the
other one. If `A` is a `Number`, its real part is returned converted back to the input
type.

If a hermitian view of a matrix is to be constructed of which the elements are neither
matrices nor numbers, an appropriate method of `hermitian` has to be implemented. In that
case, `hermitian_type` has to be implemented, too.
"""
hermitian(A::AbstractMatrix, uplo::Symbol=:U) = Hermitian(A, uplo)
hermitian(A::Number, ::Symbol=:U) = convert(typeof(A), real(A))

"""
    hermitian_type(T::Type)

The type of the object returned by `hermitian(::T, ::Symbol)`. For matrices, this is an
appropriately typed `Hermitian`, for `Number`s, it is the original type. If `hermitian` is
implemented for a custom type, so should be `hermitian_type`, and vice versa.
"""
function hermitian_type(::Type{T}) where {S, T<:AbstractMatrix{S}}
    return Hermitian{Union{S, promote_op(adjoint, S), hermitian_type(S)}, T}
end
function hermitian_type(::Type{T}) where {S<:Number, T<:AbstractMatrix{S}}
    return Hermitian{S, T}
end
function hermitian_type(::Type{T}) where {S<:AbstractMatrix, T<:AbstractMatrix{S}}
    return Hermitian{AbstractMatrix, T}
end
hermitian_type(::Type{T}) where {T<:Number} = T

_unwrap(A::Hermitian) = parent(A)
_unwrap(A::Symmetric) = parent(A)

for (S, H) in ((:Symmetric, :Hermitian), (:Hermitian, :Symmetric))
    @eval begin
        $S(A::$S) = A
        function $S(A::$S, uplo::Symbol)
            if A.uplo == char_uplo(uplo)
                return A
            else
                throw(ArgumentError("Cannot construct $($S); uplo doesn't match"))
            end
        end
        $S(A::$H) = $S(A, _sym_uplo(A.uplo))
        function $S(A::$H, uplo::Symbol)
            if A.uplo == char_uplo(uplo)
                if $H === Hermitian && !(eltype(A) <: Real) &&
                    any(!isreal, A.data[i] for i in diagind(A.data, IndexStyle(A.data)))

                    throw(ArgumentError("Cannot construct $($S)($($H))); diagonal contains complex values"))
                end
                return $S(A.data, _sym_uplo(A.uplo))
            else
                throw(ArgumentError("Cannot construct $($S); uplo doesn't match"))
            end
        end
    end
end

convert(::Type{T}, m::Union{Symmetric,Hermitian}) where {T<:Symmetric} = m isa T ? m : T(m)::T
convert(::Type{T}, m::Union{Symmetric,Hermitian}) where {T<:Hermitian} = m isa T ? m : T(m)::T

const HermOrSym{T,        S} = Union{Hermitian{T,S}, Symmetric{T,S}}
const RealHermSym{T<:Real,S} = Union{Hermitian{T,S}, Symmetric{T,S}}
const SymSymTri{T} = Union{Symmetric{T}, SymTridiagonal{T}}
const RealHermSymSymTri{T<:Real} = Union{RealHermSym{T}, SymTridiagonal{T}}
const RealHermSymComplexHerm{T<:Real,S} = Union{Hermitian{T,S}, Symmetric{T,S}, Hermitian{Complex{T},S}}
const RealHermSymComplexSym{T<:Real,S} = Union{Hermitian{T,S}, Symmetric{T,S}, Symmetric{Complex{T},S}}
const RealHermSymSymTriComplexHerm{T<:Real} = Union{RealHermSymComplexSym{T}, SymTridiagonal{T}}
const RealSymHermitian{S} = Union{Symmetric{<:Real,S}, Hermitian{<:Any,S}}
const SelfAdjoint = Union{SymTridiagonal{<:Real}, Symmetric{<:Real}, Hermitian}

wrappertype(::Union{Symmetric, SymTridiagonal}) = Symmetric
wrappertype(::Hermitian) = Hermitian

hswrapperop(::Symmetric) = symmetric
hswrapperop(::Hermitian) = hermitian

nonhermitianwrappertype(::SymSymTri{<:Real}) = Symmetric
nonhermitianwrappertype(::Hermitian{<:Real}) = Symmetric
nonhermitianwrappertype(::Hermitian) = identity

"""
    LinearAlgebra.uplo(S::Union{Symmetric, Hermitian})::Symbol

Return a `Symbol` corresponding to the stored triangular half (`:U` or `:L`) in the matrix `S`,
that is, the elements are common between `S` and `parent(S)` for that triangular half.

# Example
```jldoctest
julia> S = Symmetric([1 2; 3 4], :U)
2×2 Symmetric{Int64, Matrix{Int64}}:
 1  2
 2  4

julia> LinearAlgebra.uplo(S)
:U

julia> H = Hermitian([1 2; 3 4], :L)
2×2 Hermitian{Int64, Matrix{Int64}}:
 1  3
 3  4

julia> LinearAlgebra.uplo(H)
:L
```
"""
uplo(S::HermOrSym) = sym_uplo(S.uplo)

size(A::HermOrSym) = size(A.data)
axes(A::HermOrSym) = axes(A.data)
@inline function Base.isassigned(A::HermOrSym, i::Int, j::Int)
    @boundscheck checkbounds(Bool, A, i, j) || return false
    @inbounds if i == j || ((A.uplo == 'U') == (i < j))
        return isassigned(A.data, i, j)
    else
        return isassigned(A.data, j, i)
    end
end

@inline function getindex(A::Symmetric, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    @inbounds if i == j
        return symmetric(A.data[i, j], _sym_uplo(A.uplo))::symmetric_type(eltype(A.data))
    elseif (A.uplo == 'U') == (i < j)
        return A.data[i, j]
    else
        return transpose(A.data[j, i])
    end
end
@inline function getindex(A::Hermitian, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    @inbounds if i == j
        return hermitian(A.data[i, j], _sym_uplo(A.uplo))::hermitian_type(eltype(A.data))
    elseif (A.uplo == 'U') == (i < j)
        return A.data[i, j]
    else
        return adjoint(A.data[j, i])
    end
end

Base._reverse(A::Symmetric, dims::Integer) = reverse!(Matrix(A); dims)
Base._reverse(A::Symmetric, ::Colon) = Symmetric(reverse(A.data), A.uplo == 'U' ? :L : :U)

@propagate_inbounds function setindex!(A::Symmetric, v, i::Integer, j::Integer)
    i == j || throw(ArgumentError("Cannot set a non-diagonal index in a symmetric matrix"))
    issymmetric(v) || throw(ArgumentError("cannot set a diagonal element of a symmetric matrix to an asymmetric value"))
    setindex!(A.data, v, i, j)
    return A
end

Base._reverse(A::Hermitian, dims) = reverse!(Matrix(A); dims)
Base._reverse(A::Hermitian, ::Colon) = Hermitian(reverse(A.data), A.uplo == 'U' ? :L : :U)

@propagate_inbounds function setindex!(A::Hermitian, v, i::Integer, j::Integer)
    if i != j
        throw(ArgumentError("Cannot set a non-diagonal index in a Hermitian matrix"))
    elseif !ishermitian(v)
        throw(ArgumentError("cannot set a diagonal element of a hermitian matrix to a non-hermitian value"))
    else
        setindex!(A.data, v, i, j)
    end
    return A
end

Base.dataids(A::HermOrSym) = Base.dataids(parent(A))
Base.unaliascopy(A::Hermitian) = Hermitian(Base.unaliascopy(parent(A)), _sym_uplo(A.uplo))
Base.unaliascopy(A::Symmetric) = Symmetric(Base.unaliascopy(parent(A)), _sym_uplo(A.uplo))

_conjugation(::Union{Symmetric, Hermitian{<:Real}}) = transpose
_conjugation(::Hermitian) = adjoint

diag(A::Symmetric) = symmetric.(diag(parent(A)), _sym_uplo(A.uplo))
diag(A::Hermitian) = hermitian.(diag(parent(A)), _sym_uplo(A.uplo))

function applytri(f, A::HermOrSym)
    if A.uplo == 'U'
        f(uppertriangular(A.data))
    else
        f(lowertriangular(A.data))
    end
end

function applytri(f, A::HermOrSym, B::HermOrSym)
    if A.uplo == B.uplo == 'U'
        f(uppertriangular(A.data), uppertriangular(B.data))
    elseif A.uplo == B.uplo == 'L'
        f(lowertriangular(A.data), lowertriangular(B.data))
    elseif A.uplo == 'U'
        f(uppertriangular(A.data), uppertriangular(_conjugation(B)(B.data)))
    else # A.uplo == 'L'
        f(uppertriangular(_conjugation(A)(A.data)), uppertriangular(B.data))
    end
end
_parent_tri(U::UpperOrLowerTriangular) = parent(U)
_parent_tri(U) = U
parentof_applytri(f, args...) = applytri(_parent_tri∘f, args...)

isdiag(A::HermOrSym) = applytri(isdiag, A)

# For A<:Union{Symmetric,Hermitian}, similar(A[, neweltype]) should yield a matrix with the same
# symmetry type, uplo flag, and underlying storage type as A. The following methods cover these cases.
similar(A::Symmetric, ::Type{T}) where {T} = Symmetric(similar(parent(A), T), sym_uplo(A.uplo))
similar(A::Hermitian, ::Type{T}) where {T} = Hermitian(similar(parent(A), T), sym_uplo(A.uplo))
# On the other hand, similar(A, [neweltype,] shape...) should yield a matrix of the underlying
# storage type of A (not wrapped in a symmetry type). The following method covers these cases.
similar(A::Union{Symmetric,Hermitian}, ::Type{T}, dims::Dims{N}) where {T,N} = similar(parent(A), T, dims)

parent(A::HermOrSym) = A.data
Symmetric{T,S}(A::Symmetric{T,S}) where {T,S<:AbstractMatrix{T}} = A
Symmetric{T,S}(A::Symmetric) where {T,S<:AbstractMatrix{T}} = Symmetric{T,S}(convert(S,A.data),A.uplo)
AbstractMatrix{T}(A::Symmetric) where {T} = Symmetric(convert(AbstractMatrix{T}, A.data), _sym_uplo(A.uplo))
AbstractMatrix{T}(A::Symmetric{T}) where {T} = copy(A)
Hermitian{T,S}(A::Hermitian{T,S}) where {T,S<:AbstractMatrix{T}} = A
Hermitian{T,S}(A::Hermitian) where {T,S<:AbstractMatrix{T}} = Hermitian{T,S}(convert(S,A.data),A.uplo)
AbstractMatrix{T}(A::Hermitian) where {T} = Hermitian(convert(AbstractMatrix{T}, A.data), _sym_uplo(A.uplo))
AbstractMatrix{T}(A::Hermitian{T}) where {T} = copy(A)

copy(A::Symmetric) = (Symmetric(parentof_applytri(copy, A), _sym_uplo(A.uplo)))
copy(A::Hermitian) = (Hermitian(parentof_applytri(copy, A), _sym_uplo(A.uplo)))

function copyto!(dest::Symmetric, src::Symmetric)
    if axes(dest) != axes(src)
        @invoke copyto!(dest::AbstractMatrix, src::AbstractMatrix)
    elseif src.uplo == dest.uplo
        copytrito!(dest.data, src.data, src.uplo)
    else
        copytrito!(dest.data, transpose(Base.unalias(dest.data, src.data)), dest.uplo)
    end
    return dest
end

function copyto!(dest::Hermitian, src::Hermitian)
    if axes(dest) != axes(src)
        @invoke copyto!(dest::AbstractMatrix, src::AbstractMatrix)
    elseif src.uplo == dest.uplo
        copytrito!(dest.data, src.data, src.uplo)
    else
        copytrito!(dest.data, adjoint(Base.unalias(dest.data, src.data)), dest.uplo)
    end
    return dest
end

@propagate_inbounds function copyto!(dest::StridedMatrix, A::HermOrSym)
    if axes(dest) != axes(A)
        @invoke copyto!(dest::StridedMatrix, A::AbstractMatrix)
    else
        _copyto!(dest, Base.unalias(dest, A))
    end
    return dest
end
@propagate_inbounds function _copyto!(dest::StridedMatrix, A::HermOrSym)
    copytrito!(dest, parent(A), A.uplo)
    conjugate = A isa Hermitian
    copytri!(dest, A.uplo, conjugate)
    _symmetrize_diagonal!(dest, A)
    return dest
end
@inline function _symmetrize_diagonal!(B, A::Symmetric)
    for i = 1:size(A, 1)
        B[i,i] = symmetric(A[i,i], _sym_uplo(A.uplo))::symmetric_type(eltype(A.data))
    end
    return B
end
@inline function _symmetrize_diagonal!(B, A::Hermitian)
    for i = 1:size(A, 1)
        B[i,i] = hermitian(A[i,i], _sym_uplo(A.uplo))::hermitian_type(eltype(A.data))
    end
    return B
end

# fill[stored]!
fill!(A::HermOrSym, x) = fillstored!(A, x)
function fillstored!(A::HermOrSym{T}, x) where T
    xT = convert(T, x)
    if isa(A, Hermitian)
        ishermitian(xT) || throw(ArgumentError("cannot fill Hermitian matrix with a non-hermitian value"))
    elseif isa(A, Symmetric)
        issymmetric(xT) || throw(ArgumentError("cannot fill Symmetric matrix with an asymmetric value"))
    end
    applytri(A -> fillstored!(A, xT), A)
    return A
end

function fillband!(A::HermOrSym, x, l, u)
    if isa(A, Hermitian)
        ishermitian(x) || throw(ArgumentError("cannot fill Hermitian matrix with a non-hermitian value"))
    elseif isa(A, Symmetric)
        issymmetric(x) || throw(ArgumentError("cannot fill Symmetric matrix with an asymmetric value"))
    end
    l == -u || throw(ArgumentError(lazy"lower and upper bands must be equal in magnitude and opposite in sign, got l=$(l), u=$(u)"))
    lp = A.uplo == 'U' ? 0 : l
    up = A.uplo == 'U' ? u : 0
    applytri(A -> fillband!(A, x, lp, up), A)
    return A
end

Base.isreal(A::HermOrSym{<:Real}) = true
function Base.isreal(A::HermOrSym)
    n = size(A, 1)
    @inbounds if A.uplo == 'U'
        for j in 1:n
            for i in 1:(j - (A isa Hermitian))
                if !isreal(A.data[i,j])
                    return false
                end
            end
        end
    else
        for j in 1:n
            for i in (j + (A isa Hermitian)):n
                if !isreal(A.data[i,j])
                    return false
                end
            end
        end
    end
    return true
end

ishermitian(A::Hermitian) = true
ishermitian(A::Symmetric{<:Real}) = true
ishermitian(A::Symmetric{<:Complex}) = isreal(A)
issymmetric(A::Hermitian{<:Real}) = true
issymmetric(A::Hermitian{<:Complex}) = isreal(A)
issymmetric(A::Symmetric) = true

# check if the symmetry is known from the type
_issymmetric(::Union{SymSymTri, Hermitian{<:Real}}) = true
_issymmetric(::Any) = false

adjoint(A::Hermitian) = A
transpose(A::Symmetric) = A
adjoint(A::Symmetric{<:Real}) = A
transpose(A::Hermitian{<:Real}) = A

real(A::Symmetric{<:Real}) = A
real(A::Hermitian{<:Real}) = A
real(A::Symmetric) = Symmetric(parentof_applytri(real, A), _sym_uplo(A.uplo))
real(A::Hermitian) = Hermitian(parentof_applytri(real, A), _sym_uplo(A.uplo))
imag(A::Symmetric) = Symmetric(parentof_applytri(imag, A), _sym_uplo(A.uplo))

Base.copy(A::Adjoint{<:Any,<:Symmetric}) =
    Symmetric(copy(adjoint(A.parent.data)), ifelse(A.parent.uplo == 'U', :L, :U))
Base.copy(A::Transpose{<:Any,<:Hermitian}) =
    Hermitian(copy(transpose(A.parent.data)), ifelse(A.parent.uplo == 'U', :L, :U))

tr(A::Symmetric{<:Number}) = tr(A.data) # to avoid AbstractMatrix fallback (incl. allocations)
tr(A::Hermitian{<:Number}) = real(tr(A.data))

Base.conj(A::Symmetric) = Symmetric(parentof_applytri(conj, A), _sym_uplo(A.uplo))
Base.conj(A::Hermitian) = Hermitian(parentof_applytri(conj, A), _sym_uplo(A.uplo))
Base.conj!(A::HermOrSym) = typeof(A)(parentof_applytri(conj!, A), A.uplo)

# tril/triu
function tril(A::Hermitian, k::Integer=0)
    if A.uplo == 'U' && k <= 0
        return tril!(copy(A.data'),k)
    elseif A.uplo == 'U' && k > 0
        return tril!(copy(A.data'),-1) + tril!(triu(A.data),k)
    elseif A.uplo == 'L' && k <= 0
        return tril(A.data,k)
    else
        return tril(A.data,-1) + tril!(triu!(copy(A.data')),k)
    end
end

function tril(A::Symmetric, k::Integer=0)
    if A.uplo == 'U' && k <= 0
        return tril!(copy(transpose(A.data)),k)
    elseif A.uplo == 'U' && k > 0
        return tril!(copy(transpose(A.data)),-1) + tril!(triu(A.data),k)
    elseif A.uplo == 'L' && k <= 0
        return tril(A.data,k)
    else
        return tril(A.data,-1) + tril!(triu!(copy(transpose(A.data))),k)
    end
end

function triu(A::Hermitian, k::Integer=0)
    if A.uplo == 'U' && k >= 0
        return triu(A.data,k)
    elseif A.uplo == 'U' && k < 0
        return triu(A.data,1) + triu!(tril!(copy(A.data')),k)
    elseif A.uplo == 'L' && k >= 0
        return triu!(copy(A.data'),k)
    else
        return triu!(copy(A.data'),1) + triu!(tril(A.data),k)
    end
end

function triu(A::Symmetric, k::Integer=0)
    if A.uplo == 'U' && k >= 0
        return triu(A.data,k)
    elseif A.uplo == 'U' && k < 0
        return triu(A.data,1) + triu!(tril!(copy(transpose(A.data))),k)
    elseif A.uplo == 'L' && k >= 0
        return triu!(copy(transpose(A.data)),k)
    else
        return triu!(copy(transpose(A.data)),1) + triu!(tril(A.data),k)
    end
end

for (T, trans, real) in [(:Symmetric, :transpose, :identity), (:(Hermitian{<:Union{Real,Complex}}), :adjoint, :real)]
    @eval begin
        function dot(A::$T, B::$T)
            n = size(A, 2)
            if n != size(B, 2)
                throw(DimensionMismatch(lazy"A has dimensions $(size(A)) but B has dimensions $(size(B))"))
            end

            dotprod = $real(zero(dot(first(A), first(B))))
            @inbounds if A.uplo == 'U' && B.uplo == 'U'
                for j in 1:n
                    for i in 1:(j - 1)
                        dotprod += 2 * $real(dot(A.data[i, j], B.data[i, j]))
                    end
                    dotprod += $real(dot(A[j, j], B[j, j]))
                end
            elseif A.uplo == 'L' && B.uplo == 'L'
                for j in 1:n
                    dotprod += $real(dot(A[j, j], B[j, j]))
                    for i in (j + 1):n
                        dotprod += 2 * $real(dot(A.data[i, j], B.data[i, j]))
                    end
                end
            elseif A.uplo == 'U' && B.uplo == 'L'
                for j in 1:n
                    for i in 1:(j - 1)
                        dotprod += 2 * $real(dot(A.data[i, j], $trans(B.data[j, i])))
                    end
                    dotprod += $real(dot(A[j, j], B[j, j]))
                end
            else
                for j in 1:n
                    dotprod += $real(dot(A[j, j], B[j, j]))
                    for i in (j + 1):n
                        dotprod += 2 * $real(dot(A.data[i, j], $trans(B.data[j, i])))
                    end
                end
            end
            return dotprod
        end
    end
end

function kron(A::Hermitian{<:Union{Real,Complex},<:StridedMatrix}, B::Hermitian{<:Union{Real,Complex},<:StridedMatrix})
    resultuplo = A.uplo == 'U' || B.uplo == 'U' ? :U : :L
    C = Hermitian(Matrix{promote_op(*, eltype(A), eltype(B))}(undef, _kronsize(A, B)), resultuplo)
    return kron!(C, A, B)
end
function kron(A::Symmetric{<:Number,<:StridedMatrix}, B::Symmetric{<:Number,<:StridedMatrix})
    resultuplo = A.uplo == 'U' || B.uplo == 'U' ? :U : :L
    C = Symmetric(Matrix{promote_op(*, eltype(A), eltype(B))}(undef, _kronsize(A, B)), resultuplo)
    return kron!(C, A, B)
end

function kron!(C::Hermitian{<:Union{Real,Complex},<:StridedMatrix}, A::Hermitian{<:Union{Real,Complex},<:StridedMatrix}, B::Hermitian{<:Union{Real,Complex},<:StridedMatrix})
    size(C) == _kronsize(A, B) || throw(DimensionMismatch("kron!"))
    if ((A.uplo == 'U' || B.uplo == 'U') && C.uplo != 'U') || ((A.uplo == 'L' && B.uplo == 'L') && C.uplo != 'L')
        throw(ArgumentError("C.uplo must match A.uplo and B.uplo, got $(C.uplo) $(A.uplo) $(B.uplo)"))
    end
    _hermkron!(C.data, A.data, B.data, conj, real, A.uplo, B.uplo)
    return C
end
function kron!(C::Symmetric{<:Number,<:StridedMatrix}, A::Symmetric{<:Number,<:StridedMatrix}, B::Symmetric{<:Number,<:StridedMatrix})
    size(C) == _kronsize(A, B) || throw(DimensionMismatch("kron!"))
    if ((A.uplo == 'U' || B.uplo == 'U') && C.uplo != 'U') || ((A.uplo == 'L' && B.uplo == 'L') && C.uplo != 'L')
        throw(ArgumentError("C.uplo must match A.uplo and B.uplo, got $(C.uplo) $(A.uplo) $(B.uplo)"))
    end
    _hermkron!(C.data, A.data, B.data, identity, identity, A.uplo, B.uplo)
    return C
end

function _hermkron!(C, A, B, conj, real, Auplo, Buplo)
    n_A = size(A, 1)
    n_B = size(B, 1)
    @inbounds if Auplo == 'U' && Buplo == 'U'
        for j = 1:n_A
            jnB = (j - 1) * n_B
            for i = 1:(j-1)
                Aij = A[i, j]
                inB = (i - 1) * n_B
                for l = 1:n_B
                    for k = 1:(l-1)
                        C[inB+k, jnB+l] = Aij * B[k, l]
                        C[inB+l, jnB+k] = Aij * conj(B[k, l])
                    end
                    C[inB+l, jnB+l] = Aij * real(B[l, l])
                end
            end
            Ajj = real(A[j, j])
            for l = 1:n_B
                for k = 1:(l-1)
                    C[jnB+k, jnB+l] = Ajj * B[k, l]
                end
                C[jnB+l, jnB+l] = Ajj * real(B[l, l])
            end
        end
    elseif Auplo == 'U' && Buplo == 'L'
        for j = 1:n_A
            jnB = (j - 1) * n_B
            for i = 1:(j-1)
                Aij = A[i, j]
                inB = (i - 1) * n_B
                for l = 1:n_B
                    C[inB+l, jnB+l] = Aij * real(B[l, l])
                    for k = (l+1):n_B
                        C[inB+l, jnB+k] = Aij * conj(B[k, l])
                        C[inB+k, jnB+l] = Aij * B[k, l]
                    end
                end
            end
            Ajj = real(A[j, j])
            for l = 1:n_B
                C[jnB+l, jnB+l] = Ajj * real(B[l, l])
                for k = (l+1):n_B
                    C[jnB+l, jnB+k] = Ajj * conj(B[k, l])
                end
            end
        end
    elseif Auplo == 'L' && Buplo == 'U'
        for j = 1:n_A
            jnB = (j - 1) * n_B
            Ajj = real(A[j, j])
            for l = 1:n_B
                for k = 1:(l-1)
                    C[jnB+k, jnB+l] = Ajj * B[k, l]
                end
                C[jnB+l, jnB+l] = Ajj * real(B[l, l])
            end
            for i = (j+1):n_A
                conjAij = conj(A[i, j])
                inB = (i - 1) * n_B
                for l = 1:n_B
                    for k = 1:(l-1)
                        C[jnB+k, inB+l] = conjAij * B[k, l]
                        C[jnB+l, inB+k] = conjAij * conj(B[k, l])
                    end
                    C[jnB+l, inB+l] = conjAij * real(B[l, l])
                end
            end
        end
    else #if Auplo == 'L' && Buplo == 'L'
        for j = 1:n_A
            jnB = (j - 1) * n_B
            Ajj = real(A[j, j])
            for l = 1:n_B
                C[jnB+l, jnB+l] = Ajj * real(B[l, l])
                for k = (l+1):n_B
                    C[jnB+k, jnB+l] = Ajj * B[k, l]
                end
            end
            for i = (j+1):n_A
                Aij = A[i, j]
                inB = (i - 1) * n_B
                for l = 1:n_B
                    C[inB+l, jnB+l] = Aij * real(B[l, l])
                    for k = (l+1):n_B
                        C[inB+k, jnB+l] = Aij * B[k, l]
                        C[inB+l, jnB+k] = Aij * conj(B[k, l])
                    end
                end
            end
        end
    end
end

(-)(A::Symmetric) = Symmetric(parentof_applytri(-, A), _sym_uplo(A.uplo))
(-)(A::Hermitian) = Hermitian(parentof_applytri(-, A), _sym_uplo(A.uplo))

## Addition/subtraction
for f ∈ (:+, :-), Wrapper ∈ (:Hermitian, :Symmetric)
    @eval function $f(A::$Wrapper, B::$Wrapper)
        uplo = A.uplo == B.uplo ? _sym_uplo(A.uplo) : (:U)
        $Wrapper(parentof_applytri($f, A, B), uplo)
    end
end

for f in (:+, :-)
    @eval begin
        $f(A::Hermitian, B::Symmetric{<:Real}) = $f(A, Hermitian(parent(B), _sym_uplo(B.uplo)))
        $f(A::Symmetric{<:Real}, B::Hermitian) = $f(Hermitian(parent(A), _sym_uplo(A.uplo)), B)
        $f(A::SymTridiagonal, B::Symmetric) = $f(Symmetric(A, _sym_uplo(B.uplo)), B)
        $f(A::Symmetric, B::SymTridiagonal) = $f(A, Symmetric(B, _sym_uplo(A.uplo)))
        $f(A::SymTridiagonal{<:Real}, B::Hermitian) = $f(Hermitian(A, _sym_uplo(B.uplo)), B)
        $f(A::Hermitian, B::SymTridiagonal{<:Real}) = $f(A, Hermitian(B, _sym_uplo(A.uplo)))
    end
end

mul(A::HermOrSym, B::HermOrSym) = A * copyto!(similar(parent(B)), B)
# catch a few potential BLAS-cases
function mul(A::HermOrSym{<:BlasFloat,<:StridedMatrix}, B::AdjOrTrans{<:BlasFloat,<:StridedMatrix})
    matmul_size_check(size(A), size(B))
    T = promote_type(eltype(A), eltype(B))
    mul!(similar(B, T, (size(A, 1), size(B, 2))),
            convert(AbstractMatrix{T}, A),
            copy_oftype(B, T)) # make sure the AdjOrTrans wrapper is resolved
end
function mul(A::AdjOrTrans{<:BlasFloat,<:StridedMatrix}, B::HermOrSym{<:BlasFloat,<:StridedMatrix})
    matmul_size_check(size(A), size(B))
    T = promote_type(eltype(A), eltype(B))
    mul!(similar(B, T, (size(A, 1), size(B, 2))),
            copy_oftype(A, T), # make sure the AdjOrTrans wrapper is resolved
            convert(AbstractMatrix{T}, B))
end

function dot(x::AbstractVector, A::HermOrSym, y::AbstractVector)
    require_one_based_indexing(x, y)
    n = length(x)
    (n == length(y) == size(A, 1)) || throw(DimensionMismatch())
    data = A.data
    s = dot(first(x), first(A), first(y))
    r = zero(s+s)
    iszero(n) && return r
    if A.uplo == 'U'
        @inbounds for j = 1:length(y)
            r += dot(x[j], hswrapperop(A)(data[j,j], :U), y[j])
            @simd for i = 1:j-1
                Aij = data[i,j]
                r += dot(x[i], Aij, y[j]) + dot(x[j], _conjugation(A)(Aij), y[i])
            end
        end
    else # A.uplo == 'L'
        @inbounds for j = 1:length(y)
            r += dot(x[j], hswrapperop(A)(data[j,j], :L), y[j])
            @simd for i = j+1:length(y)
                Aij = data[i,j]
                r += dot(x[i], Aij, y[j]) + dot(x[j], _conjugation(A)(Aij), y[i])
            end
        end
    end
    return r
end

# Scaling with Number
*(A::Symmetric, x::Number) = Symmetric(parentof_applytri(y -> y * x, A), _sym_uplo(A.uplo))
*(x::Number, A::Symmetric) = Symmetric(parentof_applytri(y -> x * y, A), _sym_uplo(A.uplo))
*(A::Hermitian, x::Real) = Hermitian(parentof_applytri(y -> y * x, A), _sym_uplo(A.uplo))
*(x::Real, A::Hermitian) = Hermitian(parentof_applytri(y -> x * y, A), _sym_uplo(A.uplo))
/(A::Symmetric, x::Number) = Symmetric(parentof_applytri(y -> y/x, A), _sym_uplo(A.uplo))
/(A::Hermitian, x::Real) = Hermitian(parentof_applytri(y -> y/x, A), _sym_uplo(A.uplo))

factorize(A::HermOrSym) = _factorize(A)
function _factorize(A::HermOrSym{T}; check::Bool=true) where T
    TT = typeof(sqrt(oneunit(T)))
    if isdiag(A)
        return Diagonal(A)
    elseif TT <: BlasFloat
        return bunchkaufman(A; check=check)
    else # fallback
        return lu(A; check=check)
    end
end

logabsdet(A::RealHermSymComplexHerm) = ((l, s) = logabsdet(_factorize(A; check=false)); return real(l), s)
logabsdet(A::Symmetric{<:Real}) = logabsdet(_factorize(A; check=false))
logabsdet(A::Symmetric) = logabsdet(_factorize(A; check=false))
logdet(A::RealHermSymComplexHerm) = real(logdet(_factorize(A; check=false)))
logdet(A::Symmetric{<:Real}) = logdet(_factorize(A; check=false))
logdet(A::Symmetric) = logdet(_factorize(A; check=false))
det(A::RealHermSymComplexHerm) = real(det(_factorize(A; check=false)))
det(A::Symmetric{<:Real}) = det(_factorize(A; check=false))
det(A::Symmetric) = det(_factorize(A; check=false))

\(A::HermOrSym, B::AbstractVector) = \(factorize(A), B)
# Bunch-Kaufman solves can not utilize BLAS-3 for multiple right hand sides
# so using LU is faster for AbstractMatrix right hand side
\(A::HermOrSym, B::AbstractMatrix) = \(isdiag(A) ? Diagonal(A) : lu(A), B)

function _inv(A::HermOrSym)
    n = checksquare(A)
    B = inv!(lu(A))
    conjugate = isa(A, Hermitian)
    # symmetrize
    if A.uplo == 'U' # add to upper triangle
        @inbounds for i = 1:n, j = i:n
            B[i,j] = conjugate ? (B[i,j] + conj(B[j,i])) / 2 : (B[i,j] + B[j,i]) / 2
        end
    else # A.uplo == 'L', add to lower triangle
        @inbounds for i = 1:n, j = i:n
            B[j,i] = conjugate ? (B[j,i] + conj(B[i,j])) / 2 : (B[j,i] + B[i,j]) / 2
        end
    end
    B
end
# StridedMatrix restriction seems necessary due to inv! call in _inv above
inv(A::Hermitian{<:Any,<:StridedMatrix}) = Hermitian(_inv(A), _sym_uplo(A.uplo))
inv(A::Symmetric{<:Any,<:StridedMatrix}) = Symmetric(_inv(A), _sym_uplo(A.uplo))

function svd(A::RealHermSymComplexHerm; full::Bool=false)
    vals, vecs = eigen(A)
    I = sortperm(vals; by=abs, rev=true)
    permute!(vals, I)
    Base.permutecols!!(vecs, I)         # left-singular vectors
    V = copy(vecs)                      # right-singular vectors
    # shifting -1 from singular values to right-singular vectors
    @inbounds for i = 1:length(vals)
        if vals[i] < 0
            vals[i] = -vals[i]
            for j = 1:size(V,1); V[j,i] = -V[j,i]; end
        end
    end
    return SVD(vecs, vals, V')
end
function svd(A::RealHermSymComplexHerm{Float16}; full::Bool = false)
    T = eltype(A)
    F = svd(eigencopy_oftype(A, eigtype(T)); full)
    return SVD{T}(F)
end

function svdvals!(A::RealHermSymComplexHerm)
    vals = eigvals!(A)
    for i = 1:length(vals)
        vals[i] = abs(vals[i])
    end
    return sort!(vals, rev = true)
end

#computes U * Diagonal(abs2.(v)) * U'
function _psd_spectral_product(v, U)
    Uv = U * Diagonal(v)
    return Uv * Uv' # often faster than generic matmul by calling BLAS.herk
end

# Matrix functions
^(A::SymSymTri{<:Complex}, p::Integer) = sympow(A, p)
^(A::SelfAdjoint, p::Integer) = sympow(A, p)
function sympow(A, p::Integer)
    if p < 0
        retmat = Base.power_by_squaring(inv(A), -p)
    else
        retmat = Base.power_by_squaring(A, p)
    end
    return wrappertype(A)(retmat)
end
function ^(A::SelfAdjoint, p::Real)
    isinteger(p) && return integerpow(A, p)
    F = eigen(A)
    if all(λ -> λ ≥ 0, F.values)
        rootpower = map(λ -> λ^(p / 2), F.values)
        retmat = _psd_spectral_product(rootpower, F.vectors)
        return wrappertype(A)(retmat)
    else
        retmat = (F.vectors * Diagonal(complex.(F.values).^p)) * F.vectors'
        return nonhermitianwrappertype(A)(retmat)
    end
end
function ^(A::SymSymTri{<:Complex}, p::Real)
    isinteger(p) && return integerpow(A, p)
    return Symmetric(schurpow(A, p))
end

for func in (:cos, :sin, :tan, :cosh, :sinh, :tanh, :atan, :asinh, :cbrt)
    @eval begin
        function ($func)(A::SelfAdjoint)
            F = eigen(A)
            retmat = (F.vectors * Diagonal(($func).(F.values))) * F.vectors'
            return wrappertype(A)(retmat)
        end
    end
end

function exp(A::SelfAdjoint)
    F = eigen(A)
    rootexp = map(λ -> exp(λ / 2), F.values)
    retmat = _psd_spectral_product(rootexp, F.vectors)
    return wrappertype(A)(retmat)
end

function cis(A::SelfAdjoint)
    F = eigen(A)
    retmat = F.vectors .* cis.(F.values') * F.vectors'
    return nonhermitianwrappertype(A)(retmat)
end

for func in (:acos, :asin, :atanh)
    @eval begin
        function ($func)(A::SelfAdjoint)
            F = eigen(A)
            if all(λ -> -1 ≤ λ ≤ 1, F.values)
                retmat = (F.vectors * Diagonal(($func).(F.values))) * F.vectors'
                return wrappertype(A)(retmat)
            else
                retmat = (F.vectors * Diagonal(($func).(complex.(F.values)))) * F.vectors'
                return nonhermitianwrappertype(A)(retmat)
            end
        end
    end
end

function acosh(A::SelfAdjoint)
    F = eigen(A)
    if all(λ -> λ ≥ 1, F.values)
        retmat = (F.vectors * Diagonal(acosh.(F.values))) * F.vectors'
        return wrappertype(A)(retmat)
    else
        retmat = (F.vectors * Diagonal(acosh.(complex.(F.values)))) * F.vectors'
        return nonhermitianwrappertype(A)(retmat)
    end
end

function sincos(A::SelfAdjoint)
    n = checksquare(A)
    F = eigen(A)
    T = float(eltype(F.values))
    S, C = Diagonal(similar(A, T, (n,))), Diagonal(similar(A, T, (n,)))
    for i in eachindex(S.diag, C.diag, F.values)
        S.diag[i], C.diag[i] = sincos(F.values[i])
    end
    return wrappertype(A)((F.vectors * S) * F.vectors'), wrappertype(A)((F.vectors * C) * F.vectors')
end

function log(A::SelfAdjoint)
    F = eigen(A)
    if all(λ -> λ ≥ 0, F.values)
        retmat = (F.vectors * Diagonal(log.(F.values))) * F.vectors'
        return wrappertype(A)(retmat)
    else
        retmat = (F.vectors * Diagonal(log.(complex.(F.values)))) * F.vectors'
        return nonhermitianwrappertype(A)(retmat)
    end
end

# sqrt has rtol kwarg to handle matrices that are semidefinite up to roundoff errors
function sqrt(A::SelfAdjoint; rtol = eps(real(float(eltype(A)))) * size(A, 1))
    F = eigen(A)
    λ₀ = -maximum(abs, F.values) * rtol # treat λ ≥ λ₀ as "zero" eigenvalues up to roundoff
    if all(λ -> λ ≥ λ₀, F.values)
        rootroot = map(λ -> λ < 0 ? zero(λ) : fourthroot(λ), F.values)
        retmat = _psd_spectral_product(rootroot, F.vectors)
        return wrappertype(A)(retmat)
    else
        retmat = (F.vectors * Diagonal(sqrt.(complex.(F.values)))) * F.vectors'
        return nonhermitianwrappertype(A)(retmat)
    end
end

"""
    hermitianpart(A::AbstractMatrix, uplo::Symbol=:U) -> Hermitian
    hermitianpart(x::Number) -> Number

Return the Hermitian part of the square matrix `A`, defined as `(A + A') / 2`, as a
[`Hermitian`](@ref) matrix. For real matrices `A`, this is also known as the symmetric part
of `A`; it is also sometimes called the "operator real part". The optional argument `uplo` controls the corresponding argument of the
[`Hermitian`](@ref) view. For real matrices, the latter is equivalent to a
[`Symmetric`](@ref) view.

For scalar inputs `x`, the function returns the real part of `x`. Standard integer scalars
are promoted to `Float64` to align with the behavior of 1×1 Hermitian matrices, while other numeric types
(`Float32`, `BigInt`, `Rational`, `BigFloat`) are preserved.

See also [`hermitianpart!`](@ref) for the corresponding in-place operation.

!!! compat "Julia 1.10"
    This function requires Julia 1.10 or later.
"""

hermitianpart(A::AbstractMatrix, uplo::Symbol=:U) = Hermitian(_hermitianpart(A), uplo)
hermitianpart(x::Number) = float(real(x))

"""
    hermitianpart!(A::AbstractMatrix, uplo::Symbol=:U) -> Hermitian

Overwrite the square matrix `A` in-place with its Hermitian part `(A + A') / 2`, and return
[`Hermitian(A, uplo)`](@ref). For real matrices `A`, this is also known as the symmetric
part of `A`.

See also [`hermitianpart`](@ref) for the corresponding out-of-place operation.

!!! compat "Julia 1.10"
    This function requires Julia 1.10 or later.
"""
hermitianpart!(A::AbstractMatrix, uplo::Symbol=:U) = Hermitian(_hermitianpart!(A), uplo)

_hermitianpart(A::AbstractMatrix) = _hermitianpart!(copy_similar(A, Base.promote_op(/, eltype(A), Int)))
_hermitianpart(a::Number) = real(a)

function _hermitianpart!(A::AbstractMatrix)
    require_one_based_indexing(A)
    n = checksquare(A)
    @inbounds for j in 1:n
        A[j, j] = _hermitianpart(A[j, j])
        for i in 1:j-1
            A[i, j] = val = (A[i, j] + adjoint(A[j, i])) / 2
            A[j, i] = adjoint(val)
        end
    end
    return A
end

## structured matrix printing ##
function Base.replace_in_print_matrix(A::HermOrSym,i::Integer,j::Integer,s::AbstractString)
    ijminmax = minmax(i, j)
    inds = A.uplo == 'U' ? ijminmax : reverse(ijminmax)
    Base.replace_in_print_matrix(parent(A), inds..., s)
end