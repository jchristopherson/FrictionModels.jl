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

function friction(
    mdl::HyperbolicModel,
    t::Array{T},
    nrm,
    vel;
    p...
) where T<: Number

    npts = length(t)
    F = zeros(T, npts)
    for i in (1:npts)
        F[i] = friction(mdl, nrm(t[i]), vel(t[i]))
    end
end

function model_from_array(mdl::HyperbolicModel, x::Array{T}) where T <: Number
    HyperbolicModel(x[1], x[2], x[3], x[4], x[5], x[6])
end

function model_to_array(x::HyperbolicModel)
    [
        x.friction_coefficient, 
        x.normalization_coefficient, 
        x.dissipation_coefficient, 
        x.hysteresis_coefficient, 
        x.stribeck_velocity, 
        x.viscous_damping
    ]
end
