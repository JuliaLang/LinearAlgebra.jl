module SimpleUnits
# a minimal unitful number type
struct SimpleUnit{pow,T<:Real} <: Number
    x::T
end
SimpleUnit{pow}(x::T) where {pow,T} = SimpleUnit{pow,T}(x)
SimpleUnit{pow,T}(a::SimpleUnit{pow,S}) where {pow,T<:Real,S<:Real} = SimpleUnit{pow}(convert(T, a.x))
Base.:+(a::SimpleUnit{pow}, b::SimpleUnit{pow}) where {pow} = SimpleUnit{pow}(a.x + b.x)
Base.:*(a::SimpleUnit{apow}, b::SimpleUnit{bpow}) where {apow,bpow} = SimpleUnit{apow+bpow}(a.x * b.x)
Base.inv(a::SimpleUnit{pow}) where {pow} = SimpleUnit{-pow}(inv(a.x))
Base.:/(a::SimpleUnit, b::SimpleUnit) = a * inv(b)
Base.zero(::Type{SimpleUnit{pow,T}}) where {pow,T} = SimpleUnit{pow}(zero(T))
Base.one(::Type{SimpleUnit{pow,T}}) where {pow,T} = SimpleUnit{0}(one(T))
Base.zero(a::SimpleUnit{pow}) where {pow} = SimpleUnit{pow}(zero(a.x))
Base.one(a::SimpleUnit) = SimpleUnit{0}(a.x)
Base.isapprox(a::SimpleUnit{pow}, b::SimpleUnit{pow}; kws...) where {pow} = isapprox(a.x, b.x; kws...)
end
