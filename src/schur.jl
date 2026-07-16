# This file is a part of Julia. License is MIT: https://julialang.org/license

# Schur decomposition
"""
    Schur <: Factorization

Matrix factorization type of the Schur factorization of a matrix `A`. This is the
return type of [`schur(_)`](@ref), the corresponding matrix factorization function.

If `F::Schur` is the factorization object, the (quasi) triangular Schur factor can
be obtained via either `F.Schur` or `F.T` and the orthogonal/unitary Schur vectors
via `F.vectors` or `F.Z` such that `A = F.vectors * F.Schur * F.vectors'`. The
eigenvalues of `A` can be obtained with `F.values`.

Iterating the decomposition produces the components `F.T`, `F.Z`, and `F.values`.

# Examples
```jldoctest
julia> A = [5. 7.; -2. -4.]
2×2 Matrix{Float64}:
  5.0   7.0
 -2.0  -4.0

julia> F = schur(A)
Schur{Float64, Matrix{Float64}, Vector{Float64}}
T factor:
2×2 Matrix{Float64}:
 3.0   9.0
 0.0  -2.0
Z factor:
2×2 Matrix{Float64}:
  0.961524  0.274721
 -0.274721  0.961524
eigenvalues:
2-element Vector{Float64}:
  3.0
 -2.0

julia> F.vectors * F.Schur * F.vectors'
2×2 Matrix{Float64}:
  5.0   7.0
 -2.0  -4.0

julia> t, z, vals = F; # destructuring via iteration

julia> t == F.T && z == F.Z && vals == F.values
true
```
"""
struct Schur{Ty,S<:AbstractMatrix,C<:AbstractVector} <: Factorization{Ty}
    T::S
    Z::S
    values::C
    Schur{Ty,S,C}(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty},
                  values::AbstractVector) where {Ty,S,C} = new(T, Z, values)
end
Schur(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}, values::AbstractVector) where {Ty} =
    Schur{Ty, typeof(T), typeof(values)}(T, Z, values)
# backwards-compatible constructors (remove with Julia 2.0)
@deprecate(Schur{Ty,S}(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty},
                       values::AbstractVector) where {Ty,S},
           Schur{Ty,S,typeof(values)}(T, Z, values))

# iteration for destructuring into components
Base.iterate(S::Schur) = (S.T, Val(:Z))
Base.iterate(S::Schur, ::Val{:Z}) = (S.Z, Val(:values))
Base.iterate(S::Schur, ::Val{:values}) = (S.values, Val(:done))
Base.iterate(S::Schur, ::Val{:done}) = nothing

"""
    schur!(A) -> F::Schur

Same as [`schur`](@ref) but uses the input argument `A` as workspace.

# Examples
```jldoctest
julia> A = [5. 7.; -2. -4.]
2×2 Matrix{Float64}:
  5.0   7.0
 -2.0  -4.0

julia> F = schur!(A)
Schur{Float64, Matrix{Float64}, Vector{Float64}}
T factor:
2×2 Matrix{Float64}:
 3.0   9.0
 0.0  -2.0
Z factor:
2×2 Matrix{Float64}:
  0.961524  0.274721
 -0.274721  0.961524
eigenvalues:
2-element Vector{Float64}:
  3.0
 -2.0

julia> A
2×2 Matrix{Float64}:
 3.0   9.0
 0.0  -2.0
```
"""
schur!(A::StridedMatrix{<:BlasFloat}) = Schur(LinearAlgebra.LAPACK.gees!('V', A)...)

"""
    schur(A) -> F::Schur

Computes the Schur factorization of the matrix `A`. The (quasi) triangular Schur factor can
be obtained from the `Schur` object `F` with either `F.Schur` or `F.T` and the
orthogonal/unitary Schur vectors can be obtained with `F.vectors` or `F.Z` such that
`A = F.vectors * F.Schur * F.vectors'`. The eigenvalues of `A` can be obtained with `F.values`.

For real `A`, the Schur factorization is "quasitriangular", which means that it
is upper-triangular except with 2×2 diagonal blocks for any conjugate pair
of complex eigenvalues; this allows the factorization to be purely real even
when there are complex eigenvalues.  To obtain the (complex) purely upper-triangular
Schur factorization from a real quasitriangular factorization, you can use
`Schur{Complex}(schur(A))`.

Iterating the decomposition produces the components `F.T`, `F.Z`, and `F.values`.

# Examples
```jldoctest
julia> A = [5. 7.; -2. -4.]
2×2 Matrix{Float64}:
  5.0   7.0
 -2.0  -4.0

julia> F = schur(A)
Schur{Float64, Matrix{Float64}, Vector{Float64}}
T factor:
2×2 Matrix{Float64}:
 3.0   9.0
 0.0  -2.0
Z factor:
2×2 Matrix{Float64}:
  0.961524  0.274721
 -0.274721  0.961524
eigenvalues:
2-element Vector{Float64}:
  3.0
 -2.0

julia> F.vectors * F.Schur * F.vectors'
2×2 Matrix{Float64}:
  5.0   7.0
 -2.0  -4.0

julia> t, z, vals = F; # destructuring via iteration

julia> t == F.T && z == F.Z && vals == F.values
true
```
"""
schur(A::AbstractMatrix{T}) where {T} = schur!(eigencopy_oftype(A, eigtype(T)))
function schur(A::RealHermSymComplexHerm)
    F = eigen(A; sortby=nothing)
    return Schur(typeof(F.vectors)(Diagonal(F.values)), F.vectors, F.values)
end
function schur(A::Union{UnitUpperTriangular{T},UpperTriangular{T}}) where {T}
    t = eigtype(T)
    Z = copy_similar(A, t)
    return Schur(Z, Matrix{t}(I, size(A)), convert(Vector{t}, diag(A)))
end
function schur(A::Union{UnitLowerTriangular{T},LowerTriangular{T}}) where {T}
    t = eigtype(T)
    # double flip the matrix A
    Z = copy_similar(A, t)
    reverse!(reshape(Z, :))
    # construct "reverse" identity
    n = size(A, 1)
    J = zeros(t, n, n)
    for i in axes(J, 2)
       J[n+1-i, i] = oneunit(t)
    end
    return Schur(Z, J, reverse!(convert(Vector{t}, diag(A))))
end
function schur(A::Bidiagonal{T}) where {T}
    t = eigtype(T)
    if A.uplo == 'U'
        return Schur(Matrix{t}(A), Matrix{t}(I, size(A)), Vector{t}(A.dv))
    else # A.uplo == 'L'
        # construct "reverse" identity
        n = size(A, 1)
        J = zeros(t, n, n)
        for i in axes(J, 2)
            J[n+1-i, i] = oneunit(t)
        end
        dv = reverse!(Vector{t}(A.dv))
        ev = reverse!(Vector{t}(A.ev))
        return Schur(Matrix{t}(Bidiagonal(dv, ev, 'U')), J, dv)
    end
end

function getproperty(F::Schur, d::Symbol)
    if d === :Schur
        return getfield(F, :T)
    elseif d === :vectors
        return getfield(F, :Z)
    else
        getfield(F, d)
    end
end

Base.propertynames(F::Schur) =
    (:Schur, :vectors, fieldnames(typeof(F))...)

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::Schur)
    summary(io, F); println(io)
    println(io, "T factor:")
    show(io, mime, F.T)
    println(io, "\nZ factor:")
    show(io, mime, F.Z)
    println(io, "\neigenvalues:")
    show(io, mime, F.values)
end

# convert a (standard-form) quasi-triangular real Schur factorization into a
# triangular complex Schur factorization.
#
# Based on the "triangularize" function from GenericSchur.jl,
# released under the MIT "Expat" license by @RalphAS
function Schur{CT}(S::Schur{<:Real}) where {CT<:Complex}
    Tr = S.T
    T = CT.(Tr)
    Z = CT.(S.Z)
    n = size(T,1)
    for j=n:-1:2
        if !iszero(Tr[j,j-1])
            # We want a unitary similarity transform from
            # ┌   ┐      ┌     ┐
            # │a b│      │w₁  x│
            # │c a│ into │0  w₂│ where bc < 0 (a,b,c real)
            # └   ┘      └     ┘
            # If we write it as
            # ┌     ┐
            # │u  v'│
            # │-v u'│
            # └     ┘
            # and make the Ansatz that u is real (so v is imaginary),
            # we arrive at a Givens rotation:
            # θ = atan(sqrt(-Tr[j,j-1]/Tr[j-1,j]))
            # s,c = sin(θ), cos(θ)
            s = sqrt(abs(Tr[j,j-1]))
            c = sqrt(abs(Tr[j-1,j]))
            r = hypot(s,c)
            G = Givens(j-1,j,complex(c/r),im*(-s/r))
            lmul!(G,T)
            rmul!(T,G')
            rmul!(Z,G')
        end
    end
    return Schur(triu!(T),Z,diag(T))
end

Schur{Complex}(S::Schur{<:Complex}) = S
Schur{T}(S::Schur{T}) where {T} = S
Schur{T}(S::Schur) where {T} = Schur(T.(S.T), T.(S.Z), T <: Real && !(eltype(S.values) <: Real) ? complex(T).(S.values) : T.(S.values))

"""
    ordschur!(F::Schur, select::Union{Vector{Bool},BitVector}) -> F::Schur

Same as [`ordschur`](@ref) but overwrites the factorization `F`.
"""
function ordschur!(schur::Schur, select::Union{Vector{Bool},BitVector})
    _, _, vals = _ordschur!(schur.T, schur.Z, select)
    schur.values[:] = vals
    return schur
end

_ordschur(T::StridedMatrix{Ty}, Z::StridedMatrix{Ty}, select::Union{Vector{Bool},BitVector}) where {Ty<:BlasFloat} =
    _ordschur!(copy(T), copy(Z), select)

_ordschur!(T::StridedMatrix{Ty}, Z::StridedMatrix{Ty}, select::Union{Vector{Bool},BitVector}) where {Ty<:BlasFloat} =
    LinearAlgebra.LAPACK.trsen!(convert(Vector{BlasInt}, select), T, Z)[1:3]

"""
    ordschur(F::Schur, select::Union{Vector{Bool},BitVector}) -> F::Schur

Reorders the Schur factorization `F` of a matrix `A = Z*T*Z'` according to the logical array
`select` returning the reordered factorization `F` object. The selected eigenvalues appear
in the leading diagonal of `F.Schur` and the corresponding leading columns of
`F.vectors` form an orthogonal/unitary basis of the corresponding right invariant
subspace. In the real case, a complex conjugate pair of eigenvalues must be either both
included or both excluded via `select`.
"""
ordschur(schur::Schur, select::Union{Vector{Bool},BitVector}) =
    Schur(_ordschur(schur.T, schur.Z, select)...)

"""
    ordschur!(S::Schur, p::AbstractVector{<:Integer}) -> F::Schur

Reorders the Schur factorization `F` of a matrix `A = Z*T*Z'` according to the integer array
`p` returning the reordered factorization `F` object using the algorithm in [^BD93], see also
[^DK01]. The `i`-th diagonal entry of `F` is the `p[i]`-th diagonal element of `S`.  In the
real case, a complex conjugate pair of eigenvalues must be either both included in `p` successively.

[^BD93]   Zhaojun Bai and James W. Demmel, "On swapping diagonal blocks in real Schur form,"
        Linear Algebra and its Applications Volume 186 pp 75-95, 1993
        https://doi.org/10.1016/0024-3795(93)90286-W

[^DK01]   Daniel Kressner, "Block algorithms for reordering standard and generalized Schur forms,"
        ACM Transactions on Mathematical Software Volume 32 Issue 4 pp 521-532
        https://doi.org/10.1145/1186785.1186787
"""
@views @inbounds function ordschur!(S::Schur, p::AbstractVector{<:Integer})
    T = S.T
    Z = S.Z
    vals = S.values

    n = checksquare(T)
    size(Z,1) == n && size(Z,2) == n ||
        throw(ArgumentError("ordschur!: S.Z has incompatible size."))
    length(vals) == n ||
        throw(ArgumentError("ordschur!: S.values has incompatible length."))
    length(p) == n ||
        throw(ArgumentError("ordschur!: permutation p must have length equal to size(S.T,1)."))
    sort!(collect(p)) == collect(1:n) ||
        throw(ArgumentError("ordschur!: permutation p must contain each index 1:n exactly once."))

    # Identify 1x1 and 2x2 blocks
    sizes = ones(Int, n)
    i = 1
    while i < n
        if eltype(T) <: Real && !iszero(T[i+1,i])
            sizes[i] = 2
            sizes[i+1] = 0
            i += 2
        else
            sizes[i] = 1
            i += 1
        end
    end
    if i == n
        sizes[n] = 1
    end

    # A 2x2 real Schur block must move as a unit
    if eltype(T) <: Real
        pinv = similar(p)
        for i in 1:n
            pinv[p[i]] = i
        end
        i = 1
        while i < n
            if sizes[i] == 2
                abs(pinv[i] - pinv[i+1]) == 1 ||
                    throw(ArgumentError("ordschur!: indices $i and $(i+1) belong to the same 2x2 real Schur block and must remain adjacent in p."))
                i += 2
            else
                i += 1
            end
        end
    end

    # Initial block starts and block ids per entry
    nb = count(!iszero, sizes)
    blocks = Vector{Int}(undef, nb)
    block_id = zeros(Int, n)

    b = 0
    for i in 1:n
        if !iszero(sizes[i])
            b += 1
            blocks[b] = i
            block_id[i:i+sizes[i]-1] .= b
        end
    end

    # Desired block permutation induced by p
    pb = Vector{Int}(undef, nb)
    k = 0
    lastb = 0
    for idx in p
        b = block_id[idx]
        if b != lastb
            k += 1
            pb[k] = b
            lastb = b
        end
    end
    k == nb || throw(ArgumentError("ordschur!: invalid block permutation induced by p."))

    # current[pos] = original block id currently sitting at block position pos
    # pos_of[b] = current position of original block id b
    current = collect(1:nb)
    pos_of = collect(1:nb)

    # Workspace
    Δ = Matrix{eltype(T)}(I, 2, 2)
    M_K0  = zeros(eltype(T), 4, 4)
    M_K1  = zeros(eltype(T), 4, 4)
    M_rhs = zeros(eltype(T), 4)
    M_X   = zeros(eltype(T), 2, 2)
    M_Q   = zeros(eltype(T), 4, 4)

    # Realize desired block order by adjacent swaps
    for pos in 1:nb
        want = pb[pos]
        curpos = pos_of[want]

        while curpos > pos
            i1 = blocks[curpos-1]
            j1 = blocks[curpos]
            s1, s2 = sizes[i1], sizes[j1]
            i2, j2 = i1+s1-1, j1+s2-1

            # Swap blocks: Calls LAPACK or pure Julia version based on eltype(S.T)
            _swap_adj_schur_blocks!(T, Z, i1, i2, j1, j2, s1, s2, n, Δ, M_K0, M_K1, M_rhs, M_X, M_Q)

            # Update bookkeeping after swapping adjacent blocks of sizes s1 and s2
            sizes[i1:j2] .= 0
            sizes[i1] = s2
            sizes[i1+s2] = s1
            blocks[curpos-1] = i1
            blocks[curpos] = i1 + s2

            b1 = current[curpos-1]
            b2 = current[curpos]
            current[curpos-1] = b2
            current[curpos] = b1
            pos_of[b1] = curpos
            pos_of[b2] = curpos - 1

            curpos -= 1
        end
    end

    permute!(vals, p)
    return S
end

# Generic Julia routine
@views @inline function _swap_adj_schur_blocks!(T::AbstractMatrix, Z::AbstractMatrix,
    i1, i2, j1, j2, s1, s2, n, Δ, M_K0, M_K1, M_rhs, M_X, M_Q)
    m = s1 + s2
    rind = i1:j2
    K0  = M_K0[1:s1*s2, 1:s1*s2]
    K1  = M_K1[1:s1*s2, 1:s1*s2]
    rhs = M_rhs[1:s1*s2]
    X   = M_X[1:s1, 1:s2]
    Q   = M_Q[1:m, 1:m]

    # Solve A*X - X*B = -C
    kron!(K0, Δ[1:s2,1:s2], T[i1:i2,i1:i2])
    K0 .-= kron!(K1, transpose(T[j1:j2,j1:j2]), Δ[1:s1,1:s1])
    rhs .= -vec(T[i1:i2,j1:j2])
    ldiv!(lu!(K0), rhs)
    X .= reshape(rhs, s1, s2)

    # Build orthogonal/unitary similarity
    fill!(Q, zero(eltype(T)))
    Q[1:s1, 1:s2] .= X
    Q[s1+1:m, 1:s2] .= Δ[1:s2,1:s2]
    Q[1:s1, s2+1:m] .= Δ[1:s1,1:s1]
    F = qr!(Q)

    lmul!(adjoint(F.Q), T[rind,i1:n])
    rmul!(T[1:n,rind], F.Q)
    rmul!(Z[:,rind], F.Q)

    # Restore Schur structure
    for k in 1:m-2
        T[rind[k+2:end], rind[k]] .= zero(eltype(T))
    end
    T[i1+s2:j2, i1:i1+s2-1] .= zero(eltype(T))
    if !(eltype(T) <: Real)
        T[rind[2:end], rind[1:end-1]] .= zero(eltype(T))
    end
end

# LAPACK routine for BLAS types
@views @inline function _swap_adj_schur_blocks!(T::StridedMatrix{<:LinearAlgebra.BlasFloat}, Z::StridedMatrix{<:LinearAlgebra.BlasFloat},
    i1, _, j1, _, _, _, _, _, _, _, _, _, _)
    LAPACK.trexc!(j1, i1, T, Z)
end

"""
    GeneralizedSchur <: Factorization

Matrix factorization type of the generalized Schur factorization of two matrices
`A` and `B`. This is the return type of [`schur(_, _)`](@ref), the corresponding
matrix factorization function.

If `F::GeneralizedSchur` is the factorization object, the (quasi) triangular Schur
factors can be obtained via `F.S` and `F.T`, the left unitary/orthogonal Schur
vectors via `F.left` or `F.Q`, and the right unitary/orthogonal Schur vectors can
be obtained with `F.right` or `F.Z` such that `A=F.left*F.S*F.right'` and
`B=F.left*F.T*F.right'`. The generalized eigenvalues of `A` and `B` can be obtained
with `F.α./F.β`.

Iterating the decomposition produces the components `F.S`, `F.T`, `F.Q`, `F.Z`,
`F.α`, and `F.β`.
"""
struct GeneralizedSchur{Ty,M<:AbstractMatrix,A<:AbstractVector,B<:AbstractVector{Ty}} <: Factorization{Ty}
    S::M
    T::M
    α::A
    β::B
    Q::M
    Z::M
    function GeneralizedSchur{Ty,M,A,B}(S::AbstractMatrix{Ty}, T::AbstractMatrix{Ty},
                                        alpha::AbstractVector, beta::AbstractVector{Ty},
                                        Q::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}) where {Ty,M,A,B}
        new{Ty,M,A,B}(S, T, alpha, beta, Q, Z)
    end
end
function GeneralizedSchur(S::AbstractMatrix{Ty}, T::AbstractMatrix{Ty},
                          alpha::AbstractVector, beta::AbstractVector{Ty},
                          Q::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}) where Ty
    GeneralizedSchur{Ty, typeof(S), typeof(alpha), typeof(beta)}(S, T, alpha, beta, Q, Z)
end
# backwards-compatible constructors (remove with Julia 2.0)
@deprecate(GeneralizedSchur{Ty,M}(S::AbstractMatrix{Ty}, T::AbstractMatrix{Ty},
                                 alpha::AbstractVector, beta::AbstractVector{Ty},
                                 Q::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}) where {Ty,M},
           GeneralizedSchur{Ty,M,typeof(alpha),typeof(beta)}(S, T, alpha, beta, Q, Z))

# iteration for destructuring into components
Base.iterate(S::GeneralizedSchur) = (S.S, Val(:T))
Base.iterate(S::GeneralizedSchur, ::Val{:T}) = (S.T, Val(:Q))
Base.iterate(S::GeneralizedSchur, ::Val{:Q}) = (S.Q, Val(:Z))
Base.iterate(S::GeneralizedSchur, ::Val{:Z}) = (S.Z, Val(:α))
Base.iterate(S::GeneralizedSchur, ::Val{:α}) = (S.α, Val(:β))
Base.iterate(S::GeneralizedSchur, ::Val{:β}) = (S.β, Val(:done))
Base.iterate(S::GeneralizedSchur, ::Val{:done}) = nothing

"""
    schur!(A::StridedMatrix, B::StridedMatrix) -> F::GeneralizedSchur

Same as [`schur`](@ref) but uses the input matrices `A` and `B` as workspace.
"""
function schur!(A::StridedMatrix{T}, B::StridedMatrix{T}) where {T<:BlasFloat}
    if LAPACK.version() < v"3.6.0"
        GeneralizedSchur(LinearAlgebra.LAPACK.gges!('V', 'V', A, B)...)
    else
        GeneralizedSchur(LinearAlgebra.LAPACK.gges3!('V', 'V', A, B)...)
    end
end

"""
    schur(A, B) -> F::GeneralizedSchur

Computes the Generalized Schur (or QZ) factorization of the matrices `A` and `B`. The
(quasi) triangular Schur factors can be obtained from the `Schur` object `F` with `F.S`
and `F.T`, the left unitary/orthogonal Schur vectors can be obtained with `F.left` or
`F.Q` and the right unitary/orthogonal Schur vectors can be obtained with `F.right` or
`F.Z` such that `A=F.left*F.S*F.right'` and `B=F.left*F.T*F.right'`. The
generalized eigenvalues of `A` and `B` can be obtained with `F.α./F.β`.

Iterating the decomposition produces the components `F.S`, `F.T`, `F.Q`, `F.Z`,
`F.α`, and `F.β`.
"""
function schur(A::AbstractMatrix{TA}, B::AbstractMatrix{TB}) where {TA,TB}
    S = promote_type(eigtype(TA), TB)
    return schur!(copy_similar(A, S), copy_similar(B, S))
end

"""
    ordschur!(F::GeneralizedSchur, select::Union{Vector{Bool},BitVector}) -> F::GeneralizedSchur

Same as `ordschur` but overwrites the factorization `F`.
"""
function ordschur!(gschur::GeneralizedSchur, select::Union{Vector{Bool},BitVector})
    _, _, α, β, _, _ = _ordschur!(gschur.S, gschur.T, gschur.Q, gschur.Z, select)
    gschur.α[:] = α
    gschur.β[:] = β
    return gschur
end

_ordschur(S::StridedMatrix{Ty}, T::StridedMatrix{Ty}, Q::StridedMatrix{Ty},
    Z::StridedMatrix{Ty}, select::Union{Vector{Bool},BitVector}) where {Ty<:BlasFloat} =
        _ordschur!(copy(S), copy(T), copy(Q), copy(Z), select)

_ordschur!(S::StridedMatrix{Ty}, T::StridedMatrix{Ty}, Q::StridedMatrix{Ty},
    Z::StridedMatrix{Ty}, select::Union{Vector{Bool},BitVector}) where {Ty<:BlasFloat} =
        LinearAlgebra.LAPACK.tgsen!(convert(Vector{BlasInt}, select), S, T, Q, Z)

"""
    ordschur(F::GeneralizedSchur, select::Union{Vector{Bool},BitVector}) -> F::GeneralizedSchur

Reorders the Generalized Schur factorization `F` of a matrix pair `(A, B) = (Q*S*Z', Q*T*Z')`
according to the logical array `select` and returns a GeneralizedSchur object `F`. The
selected eigenvalues appear in the leading diagonal of both `F.S` and `F.T`, and the
left and right orthogonal/unitary Schur vectors are also reordered such that
`(A, B) = F.Q*(F.S, F.T)*F.Z'` still holds and the generalized eigenvalues of `A`
and `B` can still be obtained with `F.α./F.β`.
"""
ordschur(gschur::GeneralizedSchur, select::Union{Vector{Bool},BitVector}) =
    GeneralizedSchur(_ordschur(gschur.S, gschur.T, gschur.Q, gschur.Z, select)...)

function getproperty(F::GeneralizedSchur, d::Symbol)
    if d === :values
        return getfield(F, :α) ./ getfield(F, :β)
    elseif d === :alpha
        return getfield(F, :α)
    elseif d === :beta
        return getfield(F, :β)
    elseif d === :left
        return getfield(F, :Q)
    elseif d === :right
        return getfield(F, :Z)
    else
        getfield(F, d)
    end
end

Base.propertynames(F::GeneralizedSchur) =
    (:values, :left, :right, fieldnames(typeof(F))...)

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::GeneralizedSchur)
    summary(io, F); println(io)
    println(io, "S factor:")
    show(io, mime, F.S)
    println(io, "\nT factor:")
    show(io, mime, F.T)
    println(io, "\nQ factor:")
    show(io, mime, F.Q)
    println(io, "\nZ factor:")
    show(io, mime, F.Z)
    println(io, "\nα:")
    show(io, mime, F.α)
    println(io, "\nβ:")
    show(io, mime, F.β)
end

# Conversion
AbstractMatrix(F::Schur) = (F.Z * F.T) * F.Z'
AbstractArray(F::Schur) = AbstractMatrix(F)
Matrix(F::Schur) = Array(AbstractArray(F))
Array(F::Schur) = Matrix(F)

copy(F::Schur) = Schur(copy(F.T), copy(F.Z), copy(F.values))
copy(F::GeneralizedSchur) = GeneralizedSchur(copy(F.S), copy(F.T), copy(F.α), copy(F.β), copy(F.Q), copy(F.Z))
