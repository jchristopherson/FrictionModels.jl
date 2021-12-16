# CoulombModels.jl

"""
Applies Coulomb's model to compute the friction force at the defined state.
"""
function friction(mdl::CoulombModel, nrm::Number, vel::Number)
    if vel == 0.0
        F = mdl.coefficient * nrm
    else
        F = mdl.coefficient * nrm * sign(vel)
    end
    return F
end

# ------------------------------------------------------------------------------
# Required routines to support friction model fitting.

function model_from_array(mdl::CoulombModel, x::Array{T}) where T <: Number
    CoulombModel(x[1])
end

function model_to_array(x::CoulombModel)
    [x.coefficient]
end
