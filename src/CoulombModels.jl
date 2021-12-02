# CoulombModels.jl

"""
Applies Coulomb's model to compute the friction force at the defined state.
"""
function friction(mdl::CoulombModel, nrm::Number, vel::Number)
    if vel == 0.0
        mdl.coefficient * nrm
    else
        -mdl.coefficient * nrm * sign(vel)
    end
end
