# CoulombModels.jl

"""
Applies Coulomb's model to compute the friction force at the defined state.
"""
function friction(mdl::CoulombModel, nrm::Number, vel::Number; p...)
    if vel == 0.0
        F = mdl.coefficient * nrm
    else
        F = mdl.coefficient * nrm * sign(vel)
    end
    return (f = F, dzdt = zero(typeof(nrm)))
end

function friction(
    mdl::CoulombModel,
    t::Array{T},
    nrm,
    vel;
    p...
) where T<: Number
    
    npts = length(t)
    F = zeros(T, npts)
    for i in (1:npts)
        rsp = friction(mdl, nrm(t[i]), vel(t[i]))
        F[i] = rsp.f
    end
    return (f = F, t = t, z = zeros(T, length(t)), dzdt = zeros(T, length(t)))
end

# ------------------------------------------------------------------------------
# Required routines to support friction model fitting.

function model_from_array(mdl::CoulombModel, x::Array{T}) where T <: Number
    CoulombModel(x[1])
end

function model_to_array(x::CoulombModel)
    [x.coefficient]
end
