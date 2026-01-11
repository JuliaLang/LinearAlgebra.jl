"""
    nonzeroinds(v::AbstractVector)

Return an iterable collection of the indices of the nonzero entries of the vector `v`.
"""
nonzeroinds(v::AbstractVector) = eachindex(v)

"""
    nzrows(A, col::Integer)

Return an iterable collection of the row indices of the nonzero entries in column `col` of the matrix `A`.
"""
nzrows(A, col) = axes(A, 1)

"""
    nzcols(A, row::Integer)

Return an iterable collection of the column indices of the nonzero entries in row `row` of the matrix `A`.
"""
nzcols(A, row) = axes(A, 2)