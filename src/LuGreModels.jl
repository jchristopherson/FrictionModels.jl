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
    dzdt = vel - mdl.bristle_stiffness * abs(vel) * z / g
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
    c = mdl.bristle_damping * exp(-(vel / mdl.stribeck_velocity)^2)
    F = mdl.bristle_stiffness * z + c * dzdt + mdl.viscous_damping * vel
    return F
end

"""
Applies the Lu-Gre model to compute the friction force at the defined state.
The parameters argument 'p' is meant to accept the current bristle deformation;
however, if this parameter is not specified, a value of 0 will be utilized.
"""
function friction(
    mdl::LuGreModel,
    nrm::Number,
    vel::Float64,
    p::Number...
)
    if isempty(p)
        z = zero(typeof(nrm))
    else
        z = p[1]
    end
    dzdt = lugrevelocity(mdl, nrm, vel, z)
    F = lugrefriction(mdl, vel, z, dzdt)
    return (force = F, params = (dzdt))
end

"""
"""
function friction(
    mdl::LuGreModel,
    t::Array{T},
    vel,
    nrm,
    p::T...
) where T <: Number
    # Solve the differential equation describing bristle deformation
    function diffeq(u_, p_, t_)
        v = vel(t_)
        N = nrm(t_)
        dzdt = lugrevelocity(mdl, N, v, u_)
        return dzdt
    end
    if isempty(p)
        zi = zero(T)
    else
        zi = p[1]
    end
    tspan = (first(t), last(t))
    prob = ODEProblem(diffeq, zi, tspan)
    sol = solve(prob)   # z = sol.u
    v = vel.(sol.t)
    N = nrm.(sol.t)
    dzdt = zeros(T, length(v))
    F = zeros(T, length(v))
    for i in 1:length(v)
        dzdt[i] = lugrevelocity(mdl, N[i], v[i], sol.u[i])
        F[i] = lugrefriction(mdl, v[i], sol.u[i], dzdt[i])
    end
    return (force = F, t = sol.t, z = sol.u, dzdt = dzdt)
end
