let
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

    prune_old_LA = parse(Bool, get(ENV, "JULIA_PRUNE_OLD_LA", "false"))
    LinalgSysImg = Base.PkgId(Base.UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e"), "LinearAlgebra")
    LA = get(Base.loaded_modules, LinalgSysImg, nothing)
    if LA !== nothing && prune_old_LA
        @assert hasmethod(*, Tuple{Matrix{Float64}, Matrix{Float64}})
        ftypes_to_delete = Set{Any}()
        for methss in methods_to_delete
            meths = getglobal(Base, methss)
            push!(ftypes_to_delete, typeof(meths))
            for meth in methods(meths)
                if meth.module === LA
                    Base.delete_method(meth)
                end
            end
        end

        for meth in methods(Core.kwcall)
            meth.module === LA || continue
            sig = Base.unwrap_unionall(meth.sig)
            # signature: Tuple{typeof(Core.kwcall), NamedTuple, typeof(f), args...}
            length(sig.parameters) >= 3 || continue
            if sig.parameters[3] in ftypes_to_delete
                Base.delete_method(meth)
            end
        end
    end
end

# check in a separate block to ensure that the latest world age is used
let
    prune_old_LA = parse(Bool, get(ENV, "JULIA_PRUNE_OLD_LA", "false"))
    LinalgSysImg = Base.PkgId(Base.UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e"), "LinearAlgebra")
    LA = get(Base.loaded_modules, LinalgSysImg, nothing)
    if LA !== nothing && prune_old_LA
        @assert !hasmethod(*, Tuple{Matrix{Float64}, Matrix{Float64}})
        @assert isempty(methods(Core.kwcall,
            Tuple{NamedTuple, typeof(isapprox), Vector{Float64}, Vector{Float64}}))
    end
    prune_old_LA && Base.unreference_module(LinalgSysImg)
end

pruned_old_LA = true
