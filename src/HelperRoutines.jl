# HelperRoutines.jl

function extract_options(args)
    o = Dict{Symbol, Any}(args)
    if haskey(o, :zi)
        zi = o[:zi]
    else
        zi = 0.0
    end
    if haskey(o, :reltol)
        reltol = o[:reltol]
    else
        reltol = 1e-8
    end
    if haskey(o, :abstol)
        abstol = o[:abstol]
    else
        abstol = 1e-6
    end
    if haskey(o, :dtmax)
        dtmax = o[:dtmax]
    else
        dtmax = 1e-3
    end
    return (zi = zi, reltol = reltol, abstol = abstol, dtmax = dtmax)
end

