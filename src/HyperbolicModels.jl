# HyperbolicModels.jl

"""
Applies the Hyperbolic model proposed by Rodriguez et al.
"""
function friction(mdl::HyperbolicModel, nrm::Number, vel::Number; p...)
    F = mdl.friction_coefficient * nrm * mdl.normalization_coefficient * 
        sign(mdl.dissipation_coefficient * vel) * 
        tanh(abs(mdl.dissipation_coefficient * vel))^(2 * mdl.hysteresis_coefficient - 1) / 
        (1.0 + atan(abs(mdl.dissipation_coefficient * vel))^(2 * mdl.stribeck_velocity)) + 
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
        rsp = friction(mdl, nrm(t[i]), vel(t[i]))
        F[i] = rsp.f
    end
    return (f = F, t = t)
end

# ------------------------------------------------------------------------------
# Required routines to support friction model fitting.

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
