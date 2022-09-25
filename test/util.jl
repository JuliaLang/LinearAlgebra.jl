function delete_methods_from(
        from_module_name::AbstractString,
        modules_to_search = values(Base.loaded_modules),
    )
    num_deleted_methods = 0
    for m in modules_to_search
        m_names = names(m; all = true, imported = true)
        for name in m_names
            generic_function = try
                getproperty(m, name)
            catch
                nothing
            end
            if generic_function === nothing
                continue
            end
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

# Make sure to run this before we do `import LinearAlgebra`
delete_methods_from("LinearAlgebra")
