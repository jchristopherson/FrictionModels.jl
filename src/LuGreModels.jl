# LuGreModels.jl

# The Lu-Gre state equation (bristle velocity)
function heuristic_state_equation(
    mdl::LuGreModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T}
) where T <: Number
    z0 = z[1]
    Fc = mdl.coulomb_coefficient * nrm
    Fs = mdl.static_coefficient * nrm
    g = Fc + (Fs - Fc) * exp(-(vel / mdl.stribeck_velocity)^2)
    return [vel - mdl.bristle_stiffness * abs(vel) * z0 / g]
end

# The Lu-Gre friction force equation.
function heuristic_force_equation(
    mdl::LuGreModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T},
    dzdt::Array{T}
) where T <: Number
    return mdl.bristle_stiffness * z[1] + mdl.bristle_damping * dzdt[1] + 
        mdl.viscous_damping * vel
end

# ------------------------------------------------------------------------------
# Required routines to support friction model fitting.

function model_from_array(mdl::LuGreModel, x::Array{T}) where T <: Number
    LuGreModel(x[1], x[2], x[3], x[4], x[5], x[6])
end

function model_to_array(x::LuGreModel)
    [
        x.static_coefficient,
        x.coulomb_coefficient,
        x.stribeck_velocity,
        x.bristle_stiffness,
        x.bristle_damping,
        x.viscous_damping
    ]
end
