"""
    nonzeroinds(v::AbstractVector)

Return an iterable collection of the indices of the nonzero entries of the vector `v`.
"""
nonzeroinds(v::AbstractVector) = eachindex(v)

"""
    nzrows(A, col::Integer)

Return an iterable collection of the row indices of the nonzero entries in column `col` of the matrix `A`.
The returned indices should be sorted.
"""
function nzrows(A, col)
    require_one_based_indexing(A)
    checkbounds(A, axes(A,1), col)
    _nzrows(A, col)
end

_nzrows(A, col) = axes(A, 1)

function nzrows(A::AbstractVector, col)
    require_one_based_indexing(A)
    col == 1 ? nonzeroinds(A) : throw(BoundsError(A, (":", col)))
end

"""
    nzcols(A, row::Integer)

Return an iterable collection of the column indices of the nonzero entries in row `row` of the matrix `A`.
The returned indices should be sorted.
"""
function nzcols(A, row)
    require_one_based_indexing(A)
    checkbounds(A, row, axes(A,2))
    _nzcols(A, row)
end

_nzcols(A, row) = axes(A, 2)

# special matrices

_nzrows(::Diagonal, col) = col:col

_nzcols(::Diagonal, row) = row:row

function _nzrows(B::Bidiagonal, col)
    if B.uplo == 'U'
        return max(1, col-1):col
    else
        return col:min(length(B.dv), col+1)
    end
end

function _nzcols(B::Bidiagonal, row)
    if B.uplo == 'U'
        return row:min(length(B.dv), row+1)
    else
        return max(1, row-1):row
    end
end

function _nzrows(T::Union{Tridiagonal, SymTridiagonal}, col)
    return max(1, col-1):min(length(T.d), col+1)
end

function _nzcols(T::Union{Tridiagonal, SymTridiagonal}, row)
    return max(1, row-1):min(length(T.d), row+1)
end

function _nzrows(U::UpperOrUnitUpperTriangular, col)
    nzrows_parent = nzrows(parent(U), col)
    ind_parent = findlast(<=(col), nzrows_parent)
    return @view nzrows_parent[begin:(isnothing(ind_parent) ? (begin-1) : ind_parent)]
end

function _nzcols(U::UpperOrUnitUpperTriangular, row)
    nzcols_parent = nzcols(parent(U), row)
    ind_parent = findfirst(>=(row), nzcols_parent)
    return @view nzcols_parent[(isnothing(ind_parent) ? (end+1) : ind_parent):end]
end

function _nzrows(L::LowerOrUnitLowerTriangular, col)
    nzrows_parent = nzrows(parent(L), col)
    inds_parent = findfirst(>=(col), nzrows_parent)
    return @view nzrows_parent[(isnothing(inds_parent) ? (end+1) : inds_parent):end]
end

function _nzcols(L::LowerOrUnitLowerTriangular, row)
    nzcols_parent = nzcols(parent(L), row)
    ind_parent = findlast(<=(row), nzcols_parent)
    return @view nzcols_parent[begin:(isnothing(ind_parent) ? (begin-1) : ind_parent)]
end

function _nzrows(H::UpperHessenberg, col)
    nzrows_parent = nzrows(parent(H), col)
    inds_parent = findlast(<=(col+1), nzrows_parent)
    return @view nzrows_parent[begin:(isnothing(inds_parent) ? (begin-1) : inds_parent)]
end

function _nzcols(H::UpperHessenberg, row)
    nzcols_parent = nzcols(parent(H), row)
    ind_parent = findfirst(>=(row-1), nzcols_parent)
    return @view nzcols_parent[(isnothing(ind_parent) ? (end+1) : ind_parent):end]
end
