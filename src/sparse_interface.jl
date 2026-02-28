"""
    nonzeroinds(v::AbstractVector)

Return an iterable collection of the indices of the nonzero entries of the vector `v`.
"""
nonzeroinds(v::AbstractVector) = eachindex(v)

"""
    nzrows(A, col::Integer)

Return an iterable collection of the row indices of the nonzero entries in column `col` of the matrix `A`. The indices must be sorted.
"""
function nzrows(A, col)
    checkbounds(A, axes(A,1), col)
    axes(A, 1)
end

nzrows(A::AbstractVector, col) = col == 1 ? nonzeroinds(A) : throw(BoundsError(A, (":", col)))

"""
    nzcols(A, row::Integer)

Return an iterable collection of the column indices of the nonzero entries in row `row` of the matrix `A`. The indices must be sorted.
"""
function nzcols(A, row)
    checkbounds(A, axes(A,2), row)
    axes(A, 2)
end

# special matrices

function nzrows(D::Diagonal, col)
    checkbounds(D, axes(D,1), col)
    return col:col
end

function nzcols(D::Diagonal, row)
    checkbounds(D, row, axes(D,2))
    return row:row
end

function nzrows(B::Bidiagonal, col)
    checkbounds(B, axes(B,1), col)
    if B.uplo == 'U'
        return max(1, col-1):col
    else
        return col:min(length(B.dv), col+1)
    end
end

function nzcols(B::Bidiagonal, row)
    checkbounds(B, row, axes(B,2))
    if B.uplo == 'U'
        return row:min(length(B.dv), row+1)
    else
        return max(1, row-1):row
    end
end

function nzrows(T::Union{Tridiagonal, SymTridiagonal}, col)
    checkbounds(T, axes(T,1), col)
    return max(1, col-1):min(length(T.d), col+1)
end

function nzcols(T::Union{Tridiagonal, SymTridiagonal}, row)
    checkbounds(T, row, axes(T,2))
    return max(1, row-1):min(length(T.d), row+1)
end

function nzrows(U::UpperOrUnitUpperTriangular, col)
    checkbounds(U, axes(U,1), col)
    nzrows_parent = nzrows(parent(U), col)
    ind_parent = findlast(<=(col), nzrows_parent)
    return @view nzrows_parent[begin:(isnothing(ind_parent) ? (begin-1) : ind_parent)]
end

function nzcols(U::UpperOrUnitUpperTriangular, row)
    checkbounds(U, row, axes(U,2))
    nzcols_parent = nzcols(parent(U), row)
    ind_parent = findfirst(>=(row), nzcols_parent)
    return @view nzcols_parent[(isnothing(ind_parent) ? (end+1) : ind_parent):end]
end

function nzrows(L::LowerOrUnitLowerTriangular, col)
    checkbounds(L, axes(L,1), col)
    nzrows_parent = nzrows(parent(L), col)
    inds_parent = findfirst(>=(col), nzrows_parent)
    return @view nzrows_parent[(isnothing(inds_parent) ? (end+1) : inds_parent):end]
end

function nzcols(L::LowerOrUnitLowerTriangular, row)
    checkbounds(L, row, axes(L,2))
    nzcols_parent = nzcols(parent(L), row)
    ind_parent = findlast(<=(row), nzcols_parent)
    return @view nzcols_parent[begin:(isnothing(ind_parent) ? (begin-1) : ind_parent)]
end

function nzrows(H::UpperHessenberg, col)
    checkbounds(H, axes(H,1), col)
    nzrows_parent = nzrows(parent(H), col)
    inds_parent = findlast(<=(col+1), nzrows_parent)
    return @view nzrows_parent[begin:(isnothing(inds_parent) ? (begin-1) : inds_parent)]
end

function nzcols(H::UpperHessenberg, row)
    checkbounds(H, row, axes(H,2))
    nzcols_parent = nzcols(parent(H), row)
    ind_parent = findfirst(>=(row-1), nzcols_parent)
    return @view nzcols_parent[(isnothing(ind_parent) ? (end+1) : ind_parent):end]
end