# LuGreModels.jl

# Compute the bristle velocity.
function lugrevelocity(
    mdl::LuGreModel,
    nrm::Number,
    vel::Number,
    z::Number
)
    Fc = mdl.coulomb_coefficient * nrm
    Fs = mdl.static_coefficient * nrm
    g = Fc + (Fs - Fc) * exp(-(vel / mdl.stribeck_velocity)^2)
    dzdt = v - mdl.bristle_stiffness * abs(vel) * z / g
    return dzdt
end

# Compute the Lu-Gre friction force when the bristle deformation and bristle
# velocity are known.
function lugrefriction(
    mdl::LuGreModel,
    vel::Number,
    z::Number,
    dzdt::Number
)
    c = mdl.bristle_damping * exp(-(v / mdl.stribeck_velocity)^2)
    F = mdl.bristle_stiffness * z + c * dzdt + mdl.viscous_damping * vel
    return F
end

function friction(
    mdl::LuGreModel,
    nrm::Number,
    vel::Float64,
    z::Number = 0.0
)
    dzdt = lugrevelocity(mdl, nrm, vel, z)
    F = lugrefriction(mdl, vel, z, dzdt)
    return (force = F, bristle_velocity = dzdt)
end