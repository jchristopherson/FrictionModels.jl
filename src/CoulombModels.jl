# CoulombModels.jl

"""
Applies Coulomb's model to compute the friction force at the defined state.
"""
function friction(mdl::CoulombModel, nrm::Number, vel::Number, p::Number...)
    if vel == 0.0
        F = mdl.coefficient * nrm
    else
        F = -mdl.coefficient * nrm * sign(vel)
    end
    return (force = F, params = ())
end
