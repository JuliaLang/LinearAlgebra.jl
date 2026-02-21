"""
    nonzeroinds(v::AbstractVector)

Return an iterable collection of the indices of the nonzero entries of the vector `v`.
"""
nonzeroinds(v::AbstractVector) = eachindex(v)

"""
    nzrows(A, col::Integer)

Return an iterable collection of the row indices of the nonzero entries in column `col` of the matrix `A`.
"""
function nzrows(A, col)
    checkbounds(A, axes(A,1), col)
    axes(A, 1)
end

nzrows(A::AbstractVector, col) = col == 1 ? nonzeroinds(A) : throw(BoundsError(A, (":", col)))

"""
    nzcols(A, row::Integer)

Return an iterable collection of the column indices of the nonzero entries in row `row` of the matrix `A`.
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
