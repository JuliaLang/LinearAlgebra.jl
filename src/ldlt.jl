# This file is a part of Julia. License is MIT: https://julialang.org/license

"""
    LDLt <: Factorization

Matrix factorization type of the `LDLt` factorization of a real [`SymTridiagonal`](@ref)
matrix `S` such that `S = L*Diagonal(d)*L'`, where `L` is a [`UnitLowerTriangular`](@ref)
matrix and `d` is a vector. The main use of an `LDLt` factorization `F = ldlt(S)`
is to solve the linear system of equations `Sx = b` with `F\\b`. This is the
return type of [`ldlt`](@ref), the corresponding matrix factorization function.

The individual components of the factorization `F::LDLt` can be accessed via `getproperty`:

| Component | Description                                 |
|:---------:|:--------------------------------------------|
| `F.L`     | `L` (unit lower triangular) part of `LDLt`  |
| `F.D`     | `D` (diagonal) part of `LDLt`               |
| `F.Lt`    | `Lt` (unit upper triangular) part of `LDLt` |
| `F.d`     | diagonal values of `D` as a `Vector`        |

# Examples
```jldoctest
julia> S = SymTridiagonal([3., 4., 5.], [1., 2.])
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 3.0  1.0   ⋅
 1.0  4.0  2.0
  ⋅   2.0  5.0

julia> F = ldlt(S)
LDLt{Float64, SymTridiagonal{Float64, Vector{Float64}}}
L factor:
3×3 UnitLowerTriangular{Float64, SymTridiagonal{Float64, Vector{Float64}}}:
 1.0        ⋅         ⋅
 0.333333  1.0        ⋅
 0.0       0.545455  1.0
D factor:
3×3 Diagonal{Float64, Vector{Float64}}:
 3.0   ⋅        ⋅
  ⋅   3.66667   ⋅
  ⋅    ⋅       3.90909
```
"""
struct LDLt{T,S<:AbstractMatrix{T}} <: Factorization{T}
    data::S

    function LDLt{T,S}(data) where {T,S<:AbstractMatrix{T}}
        require_one_based_indexing(data)
        new{T,S}(data)
    end
end
LDLt(data::AbstractMatrix{T}) where {T} = LDLt{T,typeof(data)}(data)
LDLt{T}(data::AbstractMatrix) where {T} = LDLt(convert(AbstractMatrix{T}, data)::AbstractMatrix{T})

size(S::LDLt) = size(S.data)
size(S::LDLt, i::Integer) = size(S.data, i)

LDLt{T,S}(F::LDLt{T,S}) where {T,S<:AbstractMatrix{T}} = F
LDLt{T,S}(F::LDLt) where {T,S<:AbstractMatrix{T}} = LDLt{T,S}(convert(S, F.data)::S)
LDLt{T}(F::LDLt{T}) where {T} = F
LDLt{T}(F::LDLt) where {T} = LDLt(convert(AbstractMatrix{T}, F.data)::AbstractMatrix{T})

Factorization{T}(F::LDLt{T}) where {T} = F
Factorization{T}(F::LDLt) where {T} = LDLt{T}(F)

function getproperty(F::LDLt{<:Any, <:SymTridiagonal}, d::Symbol)
    Fdata = getfield(F, :data)
    if d === :d
        return Fdata.dv
    elseif d === :D
        return Diagonal(Fdata.dv)
    elseif d === :L
        return UnitLowerTriangular(Fdata)
    elseif d === :Lt
        return UnitUpperTriangular(Fdata)
    else
        return getfield(F, d)
    end
end

adjoint(F::LDLt{<:Real,<:SymTridiagonal}) = F
adjoint(F::LDLt) = LDLt(copy(adjoint(F.data)))

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::LDLt)
    summary(io, F); println(io)
    println(io, "L factor:")
    show(io, mime, F.L)
    println(io, "\nD factor:")
    show(io, mime, F.D)
end

"""
    LDLtPivoted <: Factorization

Matrix factorization type of the block `LDLt` factorization of a real [`SymTridiagonal`](@ref)
matrix `S` such that `S = L*D*L'`, where `L` is a unit lower triangular matrix and `D` is a
block diagonal matrix with `1×1` and `2×2` blocks (stored as a `SymTridiagonal`). The pivot
blocks are chosen according to the [`BunchPivot`](@ref) strategy, which does not require any
interchanges of rows or columns, so that the tridiagonal structure is preserved. In contrast
to [`LDLt`](@ref), this factorization exists and is numerically stable for any nonsingular
symmetric tridiagonal matrix, in particular for indefinite ones. For complex `S`, the
factorization is `S = L*D*transpose(L)`, since a `SymTridiagonal` matrix is complex symmetric,
not Hermitian. The main use of an `LDLtPivoted`
factorization `F = ldlt(S, BunchPivot())` is to solve the linear system of equations
`Sx = b` with `F\\b`. This is the return type of [`ldlt(_, BunchPivot())`](@ref ldlt),
the corresponding matrix factorization function.

The individual components of the factorization `F::LDLtPivoted` can be accessed via `getproperty`:

| Component  | Description                                                        |
|:----------:|:-------------------------------------------------------------------|
| `F.L`      | `L` (unit lower triangular) part of `LDLtPivoted`                  |
| `F.D`      | `D` (block diagonal) part of `LDLtPivoted` as a `SymTridiagonal`   |
| `F.Lt`     | `Lt` (unit upper triangular) part of `LDLtPivoted`                 |
| `F.d`      | diagonal values of `D` as a `Vector`                               |
| `F.blocks` | indices `k` such that `D[k:k+1,k:k+1]` is a `2×2` block of `D`     |

Internally, the factorization is stored in a `SymTridiagonal` `F.data` of the same size as `S`,
whose diagonal `F.data.dv` holds the diagonal of `D`. For a `1×1` pivot at index `k`, the
off-diagonal entry `F.data.ev[k]` holds the multiplier `L[k+1,k]`. For a `2×2` pivot at
index `k`, the entry `F.data.ev[k]` holds the off-diagonal entry `D[k+1,k]` of the pivot
block, and `F.data.ev[k+1]` (if present) holds `S[k+2,k+1]/det(D[k:k+1,k:k+1])`, from which
the multipliers `L[k+2,k] = -D[k+1,k]*F.data.ev[k+1]` and `L[k+2,k+1] = D[k,k]*F.data.ev[k+1]`
are obtained.

# Examples
```jldoctest
julia> S = SymTridiagonal([1., 2., 3.], [4., 1.])
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 1.0  4.0   ⋅
 4.0  2.0  1.0
  ⋅   1.0  3.0

julia> F = ldlt(S, BunchPivot())
LDLtPivoted{Float64, SymTridiagonal{Float64, Vector{Float64}}, Vector{Int64}}
L factor:
3×3 UnitLowerTriangular{Float64, Matrix{Float64}}:
 1.0         ⋅          ⋅
 0.0        1.0         ⋅
 0.285714  -0.0714286  1.0
D factor:
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 1.0  4.0   ⋅
 4.0  2.0  0.0
  ⋅   0.0  3.07143

julia> F.blocks
1-element Vector{Int64}:
 1

julia> F.L * F.D * F.Lt ≈ S
true
```
"""
struct LDLtPivoted{T,S<:AbstractMatrix{T},P<:AbstractVector{<:Integer}} <: Factorization{T}
    data::S
    blocks::P

    function LDLtPivoted{T,S,P}(data, blocks) where {T,S<:AbstractMatrix{T},P<:AbstractVector{<:Integer}}
        require_one_based_indexing(data, blocks)
        new{T,S,P}(data, blocks)
    end
end
LDLtPivoted(data::AbstractMatrix{T}, blocks::AbstractVector{<:Integer}) where {T} =
    LDLtPivoted{T,typeof(data),typeof(blocks)}(data, blocks)
LDLtPivoted{T}(data::AbstractMatrix, blocks::AbstractVector{<:Integer}) where {T} =
    LDLtPivoted(convert(AbstractMatrix{T}, data)::AbstractMatrix{T}, blocks)

size(F::LDLtPivoted) = size(getfield(F, :data))
size(F::LDLtPivoted, i::Integer) = size(getfield(F, :data), i)

LDLtPivoted{T,S,P}(F::LDLtPivoted{T,S,P}) where {T,S<:AbstractMatrix{T},P<:AbstractVector{<:Integer}} = F
LDLtPivoted{T,S,P}(F::LDLtPivoted) where {T,S<:AbstractMatrix{T},P<:AbstractVector{<:Integer}} =
    LDLtPivoted{T,S,P}(convert(S, getfield(F, :data))::S, convert(P, getfield(F, :blocks))::P)
LDLtPivoted{T}(F::LDLtPivoted{T}) where {T} = F
LDLtPivoted{T}(F::LDLtPivoted) where {T} =
    LDLtPivoted(convert(AbstractMatrix{T}, getfield(F, :data))::AbstractMatrix{T}, getfield(F, :blocks))

Factorization{T}(F::LDLtPivoted{T}) where {T} = F
Factorization{T}(F::LDLtPivoted) where {T} = LDLtPivoted{T}(F)

function getproperty(F::LDLtPivoted{T,<:SymTridiagonal}, d::Symbol) where {T}
    Fdata = getfield(F, :data)
    if d === :d
        return Fdata.dv
    elseif d === :D
        return _ldlt_blockdiagonal(F)
    elseif d === :L
        return UnitLowerTriangular(_ldlt_unitlower(F))
    elseif d === :Lt
        return UnitUpperTriangular(transpose(_ldlt_unitlower(F)))
    else
        return getfield(F, d)
    end
end

# the block diagonal factor, i.e., the stored data with the multipliers removed
function _ldlt_blockdiagonal(F::LDLtPivoted{<:Any,<:SymTridiagonal})
    Fdata = getfield(F, :data)
    e = zero(Fdata.ev)
    for k in getfield(F, :blocks)
        e[k] = Fdata.ev[k]
    end
    return SymTridiagonal(copy(Fdata.dv), e)
end

# the unit lower triangular factor, which has (at most) two nonzero subdiagonals
function _ldlt_unitlower(F::LDLtPivoted{T,<:SymTridiagonal}) where {T}
    Fdata = getfield(F, :data)
    d, e, blocks = Fdata.dv, Fdata.ev, getfield(F, :blocks)
    n = length(d)
    L = Matrix{T}(I, n, n)
    k, j = 1, 1
    @inbounds while k < n
        if j <= length(blocks) && blocks[j] == k
            if k + 1 < n
                L[k+2,k] = -e[k]*e[k+1]
                L[k+2,k+1] = d[k]*e[k+1]
            end
            k += 2
            j += 1
        else
            L[k+1,k] = e[k]
            k += 1
        end
    end
    return L
end

adjoint(F::LDLtPivoted{<:Real,<:SymTridiagonal}) = F
# the factorized matrix is symmetric, so its transpose is the matrix itself
transpose(F::LDLtPivoted{<:Any,<:SymTridiagonal}) = F
adjoint(F::LDLtPivoted) = LDLtPivoted(copy(adjoint(getfield(F, :data))), getfield(F, :blocks))

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::LDLtPivoted)
    summary(io, F); println(io)
    println(io, "L factor:")
    show(io, mime, F.L)
    println(io, "\nD factor:")
    show(io, mime, F.D)
end

# SymTridiagonal
"""
    ldlt!(S::SymTridiagonal, pivot::PivotingStrategy = NoPivot()) -> LDLt or LDLtPivoted

Same as [`ldlt`](@ref), but saves space by overwriting the input `S`, instead of creating a copy.

# Examples
```jldoctest
julia> S = SymTridiagonal([3., 4., 5.], [1., 2.])
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 3.0  1.0   ⋅
 1.0  4.0  2.0
  ⋅   2.0  5.0

julia> ldltS = ldlt!(S);

julia> ldltS === S
false

julia> S
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 3.0       0.333333   ⋅
 0.333333  3.66667   0.545455
  ⋅        0.545455  3.90909
```
"""
ldlt!(::SymTridiagonal{T,V}, pivot::PivotingStrategy) where {T,V} =
    throw(ArgumentError(lazy"pivoting strategy $pivot is not supported by ldlt; use NoPivot() or BunchPivot()"))

function ldlt!(S::SymTridiagonal{T,V}, ::NoPivot = NoPivot()) where {T,V}
    n = size(S,1)
    d = S.dv
    e = S.ev
    @inbounds for i in 1:n-1
        iszero(d[i]) && throw(ZeroPivotException(i))
        new_ei = e[i] / d[i]
        d[i+1] -= new_ei * e[i]
        e[i] = new_ei
    end
    return LDLt{T,SymTridiagonal{T,V}}(S)
end

# Block LDLt factorization with 1×1 and 2×2 diagonal pivots, chosen by the simplified
# pivoting strategy of Bunch and Marcia (Algorithm A1 in Numer. Linear Algebra Appl. 13
# (2006), 865-867), which requires neither interchanges nor the largest entry of the matrix.
# See the docstrings of `BunchPivot` and `LDLtPivoted` for the pivot choice and the storage
# format, respectively.
function ldlt!(S::SymTridiagonal{T,V}, ::BunchPivot) where {T,V}
    n = size(S,1)
    d = S.dv
    e = S.ev
    blocks = Int[]
    α = (sqrt(5) - 1)/2
    k = 1
    @inbounds while k < n
        a1, a2, b2 = d[k], d[k+1], e[k]
        b3 = k + 1 < n ? e[k+1] : zero(b2)
        Δ = a1 * a2 - b2 * b2
        if abs(a1 * a2) >= α * abs2(b2) || abs(Δ) <= α * abs(a1 * b3) ||
                abs(b2 * Δ) <= α * abs(a1 * a1 * b3)
            # 1×1 pivot: by the pivot choice, a zero pivot implies b2 == 0, so that the
            # matrix is singular
            iszero(a1) && throw(SingularException(k))
            lk = b2 / a1
            d[k+1] -= lk * b2
            e[k] = lk
            k += 1
        else
            # 2×2 pivot [a1 b2; b2 a2], nonsingular by construction; b2 is kept as the
            # off-diagonal entry of D
            push!(blocks, k)
            if k + 1 < n
                m = b3 / Δ
                d[k+2] -= a1 * b3 * m
                e[k+1] = m
            end
            k += 2
        end
    end
    # a final 1×1 pivot
    k == n && iszero(d[n]) && throw(SingularException(n))
    return LDLtPivoted{T,SymTridiagonal{T,V},typeof(blocks)}(S, blocks)
end

"""
    ldlt(S::SymTridiagonal, pivot::PivotingStrategy = NoPivot(); shift = 0) -> LDLt or LDLtPivoted

Compute an `LDLt` (i.e., ``LDL^T``) factorization of the real symmetric tridiagonal matrix `S` such that `S = L*Diagonal(d)*L'`
where `L` is a unit lower triangular matrix and `d` is a vector. The main use of an `LDLt`
factorization `F = ldlt(S)` is to solve the linear system of equations `Sx = b` with `F\\b`.

If `shift` is nonzero, the factorization of `S + shift*I` is computed instead.

The second argument `pivot` selects the pivoting strategy. The default `NoPivot()` computes
the factorization without pivoting, which fails with a [`ZeroPivotException`](@ref) if a zero
pivot is encountered and may be numerically unstable for indefinite matrices. With
`BunchPivot()`, a block factorization `S = L*D*L'` with a block diagonal `D` consisting
of `1×1` and `2×2` blocks is computed and returned as an [`LDLtPivoted`](@ref). The pivot
blocks are chosen according to the [`BunchPivot`](@ref) strategy, which is numerically stable
and does not require any interchanges, so that the tridiagonal structure is preserved. With
this strategy, a zero pivot can only occur for a singular matrix, in which case a
[`SingularException`](@ref) is thrown.

See also [`bunchkaufman`](@ref) for a similar, but pivoted, factorization of arbitrary symmetric or Hermitian matrices.

# Examples
```jldoctest
julia> S = SymTridiagonal([3., 4., 5.], [1., 2.])
3×3 SymTridiagonal{Float64, Vector{Float64}}:
 3.0  1.0   ⋅
 1.0  4.0  2.0
  ⋅   2.0  5.0

julia> ldltS = ldlt(S);

julia> b = [6., 7., 8.];

julia> ldltS \\ b
3-element Vector{Float64}:
 1.7906976744186047
 0.627906976744186
 1.3488372093023255

julia> S \\ b
3-element Vector{Float64}:
 1.7906976744186047
 0.627906976744186
 1.3488372093023255
```

For an indefinite matrix, the factorization without pivoting may not exist, whereas the
pivoted one does:

```jldoctest
julia> S = SymTridiagonal([0., 1., 2.], [1., 1.]);

julia> ldlt(S)
ERROR: ZeroPivotException: factorization encountered one or more zero pivots. Consider switching to a pivoted LU factorization.
[...]

julia> F = ldlt(S, BunchPivot());

julia> F \\ [2., 3., 5.]
3-element Vector{Float64}:
 -0.5
  2.0
  1.5
```
"""
function ldlt(M::SymTridiagonal{T}, pivot::PivotingStrategy = NoPivot(); shift::Number=false) where T
    S = typeof((zero(T)+shift)/one(T))
    Mₛ = SymTridiagonal{S}(copymutable_oftype(M.dv, S), copymutable_oftype(M.ev, S))
    if !iszero(shift)
        Mₛ.dv .+= shift
    end
    return ldlt!(Mₛ, pivot)
end

factorize(S::SymTridiagonal) = ldlt(S, BunchPivot())

function ldiv!(S::LDLt{<:Any,<:SymTridiagonal}, B::AbstractVecOrMat)
    require_one_based_indexing(B)
    n, nrhs = size(B, 1), size(B, 2)
    if size(S,1) != n
        throw(DimensionMismatch(lazy"Matrix has dimensions $(size(S)) but right hand side has first dimension $n"))
    end
    d = S.data.dv
    l = S.data.ev
    @inbounds begin
        for i = 2:n
            li1 = l[i-1]
            @simd for j = 1:nrhs
                B[i,j] -= li1*B[i-1,j]
            end
        end
        dn = d[n]
        @simd for j = 1:nrhs
            B[n,j] /= dn
        end
        for i = n-1:-1:1
            di = d[i]
            li = l[i]
            @simd for j = 1:nrhs
                B[i,j] /= di
                B[i,j] -= li*B[i+1,j]
            end
        end
    end
    return B
end

rdiv!(B::AbstractVecOrMat, S::LDLt{<:Any,<:SymTridiagonal}) =
    transpose(ldiv!(S, transpose(B)))

function ldiv!(F::LDLtPivoted{<:Any,<:SymTridiagonal}, B::AbstractVecOrMat)
    require_one_based_indexing(B)
    n, nrhs = size(B, 1), size(B, 2)
    if size(F,1) != n
        throw(DimensionMismatch(lazy"Matrix has dimensions $(size(F)) but right hand side has first dimension $n"))
    end
    Fdata = getfield(F, :data)
    d = Fdata.dv
    e = Fdata.ev
    blocks = getfield(F, :blocks)
    nb = length(blocks)
    # The 2×2 blocks split the rows into runs of 1×1 pivots, which are processed in tight
    # loops as in the unpivoted case, interrupted by the 2×2 blocks.
    @inbounds begin
        # forward substitution: solve L*y = b
        k = 1
        for j in 1:nb+1
            kb = j <= nb ? blocks[j] : n # start of the next 2×2 block
            # 1×1 blocks, L[i+1,i] = e[i]
            for i in k:kb-1
                li = e[i]
                @simd for c = 1:nrhs
                    B[i+1,c] -= li*B[i,c]
                end
            end
            j <= nb || break
            # 2×2 block, L[kb+2,kb] = -e[kb]*e[kb+1] and L[kb+2,kb+1] = d[kb]*e[kb+1]
            if kb + 1 < n
                dk, ek, m = d[kb], e[kb], e[kb+1]
                @simd for c = 1:nrhs
                    B[kb+2,c] -= m*(dk*B[kb+1,c] - ek*B[kb,c])
                end
            end
            k = kb + 2
        end
        # block diagonal solve D*z = y, combined with backward substitution L'*x = z
        k = n
        for j in nb:-1:0
            kb = j >= 1 ? blocks[j] + 1 : 0 # last row of the previous 2×2 block
            # 1×1 blocks
            if k == n && kb < n
                dn = d[n]
                @simd for c = 1:nrhs
                    B[n,c] /= dn
                end
                k -= 1
            end
            for i in k:-1:kb+1
                di, li = d[i], e[i]
                @simd for c = 1:nrhs
                    B[i,c] = B[i,c]/di - li*B[i+1,c]
                end
            end
            j >= 1 || break
            # 2×2 block [d1 e1; e1 d2] in rows kb-1 and kb, which is solved with the
            # explicit inverse (see the reference in the docstring of `BunchPivot`):
            # for a 2×2 pivot, the pivoting strategy guarantees |d1*d2| < α*e1^2, so
            # that the determinant Δ is computed without cancellation.
            d1, d2, e1 = d[kb-1], d[kb], e[kb-1]
            Δ⁻¹ = inv(d1*d2 - e1*e1)
            i11, i12, i22 = d2*Δ⁻¹, -e1*Δ⁻¹, d1*Δ⁻¹
            if kb < n
                m = e[kb]
                @simd for c = 1:nrhs
                    y1, y2, xk1 = B[kb-1,c], B[kb,c], B[kb+1,c]
                    B[kb-1,c] = i11*y1 + i12*y2 + e1*m*xk1
                    B[kb,c]   = i12*y1 + i22*y2 - d1*m*xk1
                end
            else
                @simd for c = 1:nrhs
                    y1, y2 = B[kb-1,c], B[kb,c]
                    B[kb-1,c] = i11*y1 + i12*y2
                    B[kb,c]   = i12*y1 + i22*y2
                end
            end
            k = kb - 2
        end
    end
    return B
end

rdiv!(B::AbstractVecOrMat, F::LDLtPivoted{<:Any,<:SymTridiagonal}) =
    transpose(ldiv!(F, transpose(B)))

function logabsdet(F::LDLt{<:Any,<:SymTridiagonal})
    it = (F.data[i,i] for i in 1:size(F, 1))
    return sum(log∘abs, it), prod(sign, it)
end

# determinants of the diagonal blocks of D
function _ldlt_blockdets(F::LDLtPivoted{T,<:SymTridiagonal}) where {T}
    Fdata = getfield(F, :data)
    d, e, blocks = Fdata.dv, Fdata.ev, getfield(F, :blocks)
    n = length(d)
    dets = Vector{T}(undef, n - length(blocks))
    k, j, i = 1, 1, 1
    @inbounds while k <= n
        if j <= length(blocks) && blocks[j] == k
            dets[i] = d[k]*d[k+1] - e[k]*e[k]
            k += 2
            j += 1
        else
            dets[i] = d[k]
            k += 1
        end
        i += 1
    end
    return dets
end

function logabsdet(F::LDLtPivoted{<:Any,<:SymTridiagonal})
    dets = _ldlt_blockdets(F)
    return sum(log∘abs, dets), prod(sign, dets)
end

# Conversion methods
function SymTridiagonal(F::LDLt{<:Any, <:SymTridiagonal})
    e = copy(F.data.ev)
    d = copy(F.data.dv)
    e .*= d[1:end-1]
    d[2:end] += e .* F.data.ev
    SymTridiagonal(d, e)
end
AbstractMatrix(F::LDLt) = SymTridiagonal(F)
AbstractArray(F::LDLt) = AbstractMatrix(F)
Matrix(F::LDLt) = Array(AbstractArray(F))
Array(F::LDLt) = Matrix(F)

function SymTridiagonal(F::LDLtPivoted{<:Any,<:SymTridiagonal})
    Fdata = getfield(F, :data)
    d, e, blocks = Fdata.dv, Fdata.ev, getfield(F, :blocks)
    n = length(d)
    dA = copy(d)
    eA = copy(e)
    k, j = 1, 1
    @inbounds while k < n
        if j <= length(blocks) && blocks[j] == k
            if k + 1 < n
                Δ = d[k]*d[k+1] - e[k]*e[k]
                eA[k+1] = e[k+1]*Δ
                dA[k+2] += d[k]*eA[k+1]*e[k+1]
            end
            k += 2
            j += 1
        else
            eA[k] = e[k]*d[k]
            dA[k+1] += e[k]*eA[k]
            k += 1
        end
    end
    return SymTridiagonal(dA, eA)
end
AbstractMatrix(F::LDLtPivoted) = SymTridiagonal(F)
AbstractArray(F::LDLtPivoted) = AbstractMatrix(F)
Matrix(F::LDLtPivoted) = Array(AbstractArray(F))
Array(F::LDLtPivoted) = Matrix(F)
