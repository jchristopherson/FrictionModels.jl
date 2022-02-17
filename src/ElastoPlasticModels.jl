# ElastoPlasticModels.jl

# The Elasto-Plastic alpha function.
function ep_alpha(
    mdl::ElastoPlasticModel,
    z::Number
)
    zba = mdl.breakaway_displacement
    zmax = mdl.max_bristle_displacement
    if abs(z) < zba
        alpha = 0.0
    elseif  abs(z) >= zba && abs(z) <= zmax
        half = 0.5
        arg = pi * (z - ((zmax + zba) / 2.0)) / (zmax - zba)
        alpha = half * (1.0 + sin(arg))
    else
        alpha = 1.0
    end
    return alpha
end

# The Elasto-Plastic model Stribeck equation.
function ep_stribeck(
    mdl::ElastoPlasticModel,
    v::T,
    N::T
) where T <: Number
    fc = mdl.coulomb_coefficient * N
    fs = mdl.static_coefficient * N
    return fc + (fs - fc) * exp(-(v / mdl.stribeck_velocity)^2)
end

# The Elasto-Plastic state equation (bristle velocity)
function heuristic_state_equation(
    mdl::ElastoPlasticModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T}
) where T <: Number
    alpha = ep_alpha(mdl, z[1])
    fss = ep_stribeck(mdl, vel, nrm)
    zdot = vel * (1.0 - 
        alpha * mdl.bristle_stiffness * sign(vel) * z[1] / fss)
    return [zdot]
end

# The Elasto-Plastic friction force equation
function heuristic_force_equation(
    mdl::ElastoPlasticModel,
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

function model_from_array(mdl::ElastoPlasticModel, x::Array{T}) where T <: Number
    ElastoPlasticModel(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8])
end

function model_to_array(x::ElastoPlasticModel)
    [
        x.static_coefficient,
        x.coulomb_coefficient,
        x.stribeck_velocity,
        x.bristle_stiffness,
        x.bristle_damping,
        x.breakaway_displacement,
        x.max_bristle_displacement,
        x.viscous_damping
    ]
end

