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
however, if this parameter is not specified, a value of 0 will be utilized.  The
parameter is expected to be specified as 'z'.

The output of the function is a named tuple containing the following fields.
- f: The friction force.
- dzdt: The bristle deformation velocity.
"""
function friction(
    mdl::LuGreModel,
    nrm::Number,
    vel::Number;
    p...
)
    args = Dict{Symbol, Any}(p)
    if haskey(args, :z)
        z = args[:z]
    else
        z = zero(typeof(nrm))
    end

    dzdt = lugrevelocity(mdl, nrm, vel, z)
    F = lugrefriction(mdl, vel, z, dzdt)

    return (f = F, dzdt = dzdt)
end

"""
Applies the Lu-Gre model to compute the friction force at the defined state by
performing the integration of the differential equation describing bristle
deformation.  The user can specify that the solver choose the time step by only
supplying a start and end simulation time, or the user can specify the time 
points the solver returns by specifying more than 2 time points.  Notice, if
supplying time points at which the solver should return the results, the time
points do not need to be evenly spaced; however, they do need to be either
monotonically increasing or decreasing.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'p' input parameter.  The following
inputs are currently recognized.
- zi: The initial condition or initial bristle deformation.  The default is 0.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.

The output of the function is a named tuple containing the following fields.
- t: The solution time points.
- f: The friction force.
- z: The bristle deformation.
- dzdt: The bristle deformation velocity.
"""
function friction(
    mdl::LuGreModel,
    t::Array{T},
    nrm,
    vel;
    p...
) where T <: Number

    function diffeq(u_, p_, t_)
        v = vel(t_)
        N = nrm(t_)
        return lugrevelocity(mdl, N, v, u_)
    end
    
    args = extract_options(p)

    tspan = (first(t), last(t))
    prob = ODEProblem(diffeq, args.zi, tspan)

    if length(t) == 2
        sol = solve(
            prob, 
            alg_hints = [:stiff], 
            reltol = args.reltol, 
            abstol = args.abstol,
            dtmax = args.dtmax
        )
    elseif length(t) > 2
        sol = solve(
            prob, 
            alg_hints = [:stiff], 
            saveat = t,
            reltol = args.reltol, 
            abstol = args.abstol,
            dtmax = args.dtmax
        )
    else
        error("No time range for the solution has been given.  Ensure at least two values are provided in the time array input.")
    end

    dzdt = zeros(T, length(sol.t))
    F = zeros(T, length(sol.t))
    for i in 1:length(sol.t)
        v = vel(sol.t[i])
        N = nrm(sol.t[i])
        dzdt[i] = lugrevelocity(mdl, N, v, sol.u[i])
        F[i] = lugrefriction(mdl, v, sol.u[i], dzdt[i])
    end

    return (f = F, t = sol.t, z = sol.u, dzdt = dzdt)
end

function model_from_array(mdl::LuGreModel, x::Array{T}) where T <: Number
    LuGreModel(x[1], x[2], x[3], x[4], x[5], x[6])
end

function model_to_array(x::LuGreModel)
    [
        x.static_coefficient,
        x.coulomb_coefficient,
        x.static_coefficient,
        x.bristle_stiffness,
        x.bristle_damping,
        x.viscous_damping
    ]
end
