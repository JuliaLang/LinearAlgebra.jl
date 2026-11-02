# This file is a part of Julia. License is MIT: https://julialang.org/license

# StridedArrays

# This test file defines array types and functions
# useful for testing the strided array interface.

module StridedArrays

export strided_ptr
export check_strided_get

function strided_ptr(f, a::AbstractArray{T}) where {T}
    a_cconv = Base.cconvert(Ptr{T}, a)
    GC.@preserve a_cconv begin
        f(Base.unsafe_convert(Ptr{T}, a_cconv))
    end
end

"""
    check_strided_get(a::AbstractArray{T,N})

Test that array `a` implements the strided array interface for reading.
Checks stride consistency and that `unsafe_load` matches `getindex`.
"""
function check_strided_get(a::AbstractArray{T,N})::Nothing where {T, N}
    if !isbitstype(eltype(a))
        error("a doesn't have isbits elements")
    end
    # Putting strided_ptr before the loop means that strided_ptr shouldn't error for empty arrays
    strided_ptr(a) do a_ptr
        for d in 1:N
            if stride(a, d) != strides(a)[d]
                error("stride(a, d) doesn't equal strides(a)[d] for dimension $(d)")
            end
        end
        for i in CartesianIndices(a)
            el_ptr = a_ptr
            for d in 1:N
                stride_in_bytes = stride(a, d) * Base.elsize(typeof(a))
                first_idx = first(axes(a, d))
                el_ptr += (i[d] - first_idx) * stride_in_bytes
            end
            if unsafe_load(el_ptr) !== a[i]
                error("getindex and unsafe_load mismatch at index $(i)")
            end
        end
    end
    nothing
end

end
