methods_to_delete =
[
:adjoint
:transpose
:inv
:literal_pow
:\
:/
:isapprox
:copyto!
:*
:muladd
:copyto!
:isone
:kron!
:kron
:^
:exp
:cis
:log
:sqrt
:cbrt
:inv
:cos
:sin
:sincos
:tan
:cosh
:sinh
:tanh
:acos
:asin
:atan
:acosh
:asinh
:atanh
:sec
:sech
:csc
:csch
:cot
:coth
:asec
:asech
:acsc
:acot
:acoth
:acsch
]

let
    LA = get(Base.loaded_modules, Base.PkgId(Base.UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e"), "LinearAlgebra"), nothing)
    if LA !== nothing
        @assert hasmethod(*, Tuple{Matrix{Float64}, Matrix{Float64}})
        for methss in methods_to_delete
            meths = getglobal(Base, methss)
            for meth in methods(meths)
                if meth.module === LA
                    Base.delete_method(meth)
                end
            end
        end
    end
    Base.unreference_module(Base.PkgId(Base.UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e"), "LinearAlgebra"))
end

@assert !hasmethod(*, Tuple{Matrix{Float64}, Matrix{Float64}})

pruned_old_LA = true
