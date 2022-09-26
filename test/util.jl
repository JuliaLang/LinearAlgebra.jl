function delete_methods_from(
        from_module_name::AbstractString,
        modules_to_search = values(Base.loaded_modules),
    )
    num_deleted_methods = 0
    for m in modules_to_search
        m_names = names(m; all = true, imported = true)
        for name in m_names
            if !isdefined(m, name)
                continue
            end
            generic_function = getproperty(m, name)
            for method in methods(generic_function)
                if string(method.module) == from_module_name
                    @debug "$(num_deleted_methods). Deleting method: $(method)"
                    Base.delete_method(method)
                    num_deleted_methods += 1
                end
            end
        end
    end
    return num_deleted_methods
end

function delete_all_methods()
    delete_methods_from("LinearAlgebra")
    delete_methods_from("SparseArrays")
    return nothing
end
