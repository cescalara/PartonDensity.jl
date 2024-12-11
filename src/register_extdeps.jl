# This file is a part of JuliaPackageTemplate.jl, licensed under the MIT License (MIT).

struct _WeakDepNotLoadedErrorHint
    func::Function
    weakdeps::Vector{Symbol}
end

function (hint::_WeakDepNotLoadedErrorHint)(io, exc, _, _)
    thismodule = @__MODULE__
    func = hint.func
    ext_name = Symbol("$(thismodule)$(join(hint.weakdeps))Ext")
    if exc.f === func && (Base.get_extension(thismodule, ext_name) === nothing)
        deplist = join(hint.weakdeps, ", ")
        quoted_deplist = join(["`$dep`" for dep in hint.weakdeps], ", ")
        print(io, "\nNOTE: Using `$func` requires $quoted_deplist to be loaded to activate `$ext_name`, e.g. via `import $deplist`.")
    end
end
 
function _register_extension_deps(func_deps::Pair{<:Function,<:Union{Symbol,AbstractVector{<:Symbol}}}...) 
    if isdefined(Base.Experimental, :register_error_hint)
        for (func, deps) in func_deps
            dep_vec = deps isa Symbol ? [deps] : deps
            Base.Experimental.register_error_hint(_WeakDepNotLoadedErrorHint(func, dep_vec), MethodError)
        end
    end
end
