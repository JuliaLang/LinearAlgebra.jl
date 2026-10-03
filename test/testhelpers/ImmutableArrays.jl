# This file is a part of Julia. License is MIT: https://julialang.org/license

# ImmutableArrays (arrays that implement getindex but not setindex!)

# This test file defines an array wrapper that is immutable. It can be used to
# test the action of methods on immutable arrays.

module ImmutableArrays

using LinearAlgebra

export ImmutableArray

"An immutable wrapper type for arrays."
struct ImmutableArray{T,N,A<:AbstractArray} <: AbstractArray{T,N}
    data::A
end

ImmutableArray(data::AbstractArray{T,N}) where {T,N} = ImmutableArray{T,N,typeof(data)}(data)

# Minimal AbstractArray interface
Base.size(A::ImmutableArray) = size(A.data)
Base.size(A::ImmutableArray, d) = size(A.data, d)
Base.getindex(A::ImmutableArray, i...) = getindex(A.data, i...)

# The immutable array remains immutable after conversion to AbstractArray
AbstractArray{T}(A::ImmutableArray) where {T} = ImmutableArray(AbstractArray{T}(A.data))
AbstractArray{T,N}(A::ImmutableArray{S,N}) where {S,T,N} = ImmutableArray(AbstractArray{T,N}(A.data))

Base.copy(A::ImmutableArray) = ImmutableArray(copy(A.data))
Base.zero(A::ImmutableArray) = ImmutableArray(zero(A.data))

Base.:(-)(A::ImmutableArray) = ImmutableArray(-A.data)
Base.:(+)(A::ImmutableArray, B::ImmutableArray) = ImmutableArray(A.data + B.data)
Base.:(-)(A::ImmutableArray, B::ImmutableArray) = ImmutableArray(A.data - B.data)

Base.:(*)(A::ImmutableArray, x::Number) = ImmutableArray(A.data * x)
Base.:(*)(x::Number, A::ImmutableArray) = ImmutableArray(x * A.data)

Base.:(*)(A::ImmutableArray, B::ImmutableArray) = ImmutableArray(A.data * B.data)

function LinearAlgebra.eigen(S::SymTridiagonal{T, <:ImmutableArray{T,1}}) where {T}
    # Use the underlying data for the eigen computation
    S2 = SymTridiagonal(diag(S), diag(S,1))
    eigvals, eigvecs = eigen(S2)
    return Eigen(ImmutableArray(eigvals), ImmutableArray(eigvecs))
end

function LinearAlgebra.eigen(S::Symmetric{T, <:ImmutableArray{T,2}}) where {T<:Real}
    # Use the underlying data for the eigen computation
    S2 = Symmetric(parent(S).data)
    eigvals, eigvecs = eigen(S2)
    return Eigen(ImmutableArray(eigvals), ImmutableArray(eigvecs))
end

function LinearAlgebra.eigen(S::Hermitian{T, <:ImmutableArray{T,2}}) where {T<:Union{Real,Complex}}
    # Use the underlying data for the eigen computation
    S2 = Hermitian(parent(S).data)
    eigvals, eigvecs = eigen(S2)
    return Eigen(ImmutableArray(eigvals), ImmutableArray(eigvecs))
end

end
