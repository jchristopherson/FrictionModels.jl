# HyperbolicModels.jl

"""
Applies the Hyperbolic model proposed by Rodriguez et al.
"""
function friction(mdl::HyperbolicModel, nrm::Number, vel::Number)
    F = mdl.friction_coefficient * nrm * mdl.normalization_coefficient * 
        tanh(mdl.dissipation_coefficient * vel)^(2 * mdl.hysteresis_coefficient - 1) / 
        (1.0 + atan(mdl.dissipation_coefficient * vel)^(2 * mdl.stribeck_velocity)) + 
        mdl.viscous_damping * vel
    return (f = F, params = ())
end
