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
to [`LDLt`](@ref), this factorization exists and is numerically stable for any symmetric
tridiagonal matrix, in particular for indefinite ones. The main use of an `LDLtPivoted`
factorization `F = ldlt(S; alg=BunchPivot())` is to solve the linear system of equations
`Sx = b` with `F\\b`. This is the return type of [`ldlt(_; alg=BunchPivot())`](@ref ldlt),
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

julia> F = ldlt(S; alg=BunchPivot())
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
    ldlt!(S::SymTridiagonal; alg::PivotingStrategy = NoPivot()) -> LDLt or LDLtPivoted

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
ldlt!(S::SymTridiagonal; alg::PivotingStrategy = NoPivot()) = _ldlt!(S, alg)

_ldlt!(::SymTridiagonal{T,V}, alg::PivotingStrategy) where {T,V} =
    throw(ArgumentError(lazy"pivoting strategy $alg is not supported by ldlt; use NoPivot() or BunchPivot()"))

function _ldlt!(S::SymTridiagonal{T,V}, ::NoPivot) where {T,V}
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

# Block LDLt factorization with 1×1 and 2×2 diagonal pivots, chosen by Bunch's strategy
# (J. R. Bunch, SIAM J. Numer. Anal. 11 (1974), 521-528), which requires no interchanges
# for tridiagonal matrices. See the docstrings of `BunchPivot` and `LDLtPivoted` for the
# pivot choice and the storage format, respectively.
function _ldlt!(S::SymTridiagonal{T,V}, ::BunchPivot) where {T,V}
    n = size(S,1)
    d = S.dv
    e = S.ev
    blocks = Int[]
    # α minimizes the bound on the element growth, see the reference above
    α = (sqrt(5) - 1)/2
    σ = max(maximum(abs, d; init=abs(zero(T))), maximum(abs, e; init=abs(zero(T))))
    k = 1
    @inbounds while k < n
        dk, ek = d[k], e[k]
        if σ*abs(dk) >= α*abs2(ek)
            # 1×1 pivot: a zero pivot implies ek == 0, i.e., a structurally singular matrix
            iszero(dk) && throw(ZeroPivotException(k))
            lk = ek / dk
            d[k+1] -= lk * ek
            e[k] = lk
            k += 1
        else
            # 2×2 pivot [dk ek; ek d[k+1]], nonsingular by construction; ek is kept as the
            # off-diagonal entry of D
            push!(blocks, k)
            if k + 1 < n
                ek1 = e[k+1]
                Δ = dk * d[k+1] - ek * ek
                m = ek1 / Δ
                d[k+2] -= dk * ek1 * m
                e[k+1] = m
            end
            k += 2
        end
    end
    return LDLtPivoted{T,SymTridiagonal{T,V},typeof(blocks)}(S, blocks)
end

"""
    ldlt(S::SymTridiagonal; shift = 0, alg::PivotingStrategy = NoPivot()) -> LDLt or LDLtPivoted

Compute an `LDLt` (i.e., ``LDL^T``) factorization of the real symmetric tridiagonal matrix `S` such that `S = L*Diagonal(d)*L'`
where `L` is a unit lower triangular matrix and `d` is a vector. The main use of an `LDLt`
factorization `F = ldlt(S)` is to solve the linear system of equations `Sx = b` with `F\\b`.

If `shift` is nonzero, the factorization of `S + shift*I` is computed instead.

The keyword argument `alg` selects the pivoting strategy. The default `NoPivot()` computes
the factorization without pivoting, which fails with a [`ZeroPivotException`](@ref) if a zero
pivot is encountered and may be numerically unstable for indefinite matrices. With
`alg = BunchPivot()`, a block factorization `S = L*D*L'` with a block diagonal `D` consisting
of `1×1` and `2×2` blocks is computed and returned as an [`LDLtPivoted`](@ref). The pivot
blocks are chosen according to the [`BunchPivot`](@ref) strategy, which is numerically stable
and does not require any interchanges, so that the tridiagonal structure is preserved.

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

julia> F = ldlt(S; alg=BunchPivot());

julia> F \\ [2., 3., 5.]
3-element Vector{Float64}:
 -0.5
  2.0
  1.5
```
"""
function ldlt(M::SymTridiagonal{T}; shift::Number=false, alg::PivotingStrategy=NoPivot()) where T
    S = typeof((zero(T)+shift)/one(T))
    Mₛ = SymTridiagonal{S}(copymutable_oftype(M.dv, S), copymutable_oftype(M.ev, S))
    if !iszero(shift)
        Mₛ.dv .+= shift
    end
    return ldlt!(Mₛ; alg)
end

factorize(S::SymTridiagonal) = ldlt(S)

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
    @inbounds begin
        # forward substitution: solve L*y = b
        k, j = 1, 1
        while k < n
            if j <= nb && blocks[j] == k
                # 2×2 block, L[k+2,k] = -e[k]*e[k+1] and L[k+2,k+1] = d[k]*e[k+1]
                if k + 1 < n
                    dk, ek, m = d[k], e[k], e[k+1]
                    @simd for c = 1:nrhs
                        B[k+2,c] -= m*(dk*B[k+1,c] - ek*B[k,c])
                    end
                end
                k += 2
                j += 1
            else
                # 1×1 block, L[k+1,k] = e[k]
                lk = e[k]
                @simd for c = 1:nrhs
                    B[k+1,c] -= lk*B[k,c]
                end
                k += 1
            end
        end
        # block diagonal solve D*z = y, combined with backward substitution L'*x = z
        k, j = n, nb
        while k >= 1
            if j >= 1 && blocks[j] == k - 1
                # 2×2 block in rows k-1 and k
                d1, d2, e1 = d[k-1], d[k], e[k-1]
                Δ = d1*d2 - e1*e1
                if k < n
                    m = e[k]
                    @simd for c = 1:nrhs
                        y1, y2, xk1 = B[k-1,c], B[k,c], B[k+1,c]
                        B[k-1,c] = (d2*y1 - e1*y2)/Δ + e1*m*xk1
                        B[k,c]   = (d1*y2 - e1*y1)/Δ - d1*m*xk1
                    end
                else
                    @simd for c = 1:nrhs
                        y1, y2 = B[k-1,c], B[k,c]
                        B[k-1,c] = (d2*y1 - e1*y2)/Δ
                        B[k,c]   = (d1*y2 - e1*y1)/Δ
                    end
                end
                k -= 2
                j -= 1
            else
                # 1×1 block
                dk = d[k]
                if k < n
                    lk = e[k]
                    @simd for c = 1:nrhs
                        B[k,c] = B[k,c]/dk - lk*B[k+1,c]
                    end
                else
                    @simd for c = 1:nrhs
                        B[k,c] /= dk
                    end
                end
                k -= 1
            end
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
