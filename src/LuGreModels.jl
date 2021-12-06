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
Applies the Lu-Gre model to compute the friction force at the defined state by
performing the integration of the differential equation describing bristle
deformation.  The user can specify that the solver choose the time step by only
supplying a start and end simulation time, or the user can specify the time 
points the solver returns by specifying more than 2 time points.  Notice, if
supplying time points at which the solver should return the results, the time
points do not need to be evenly spaced; however, they do need to be either
monotonically increasing or decreasing.
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
    if length(t) == 2
        sol = solve(prob, alg_hints = [:stiff])
    elseif length(t) > 2
        sol = solve(prob, alg_hints = [:stiff], saveat = t)
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
    return (force = F, t = sol.t, z = sol.u, dzdt = dzdt)
end

"""
Uses a Levenberg-Marquardt solver to compute the best fit of a Lu-Gre model to
a predefined data set.  As the supplied data is discretely sampled but the 
solver may require points between those supplied, the solver will utilize linear
interpolation to estimate values between supplied data points.

The routine returns the fitted model along with the LsqFitResults type returned
from the Levenberg-Marquardt solver.  The LsqFitResults allow exploration of 
error margins, confidence intervals, etc.  Both results are returned as a 
named tuple with the fitted model available as 'model' and the LsqFitResults
available as 'fit'.
"""
function fit_model(
    mdl::LuGreModel,
    timedata::Array{T},
    frictiondata::Array{T},
    normaldata::Array{T},
    velocitydata::Array{T}
) where T <: Number
    # Set up the interpolation for each data set
    normal_interp = LinearInterpolation(
        timedata, normaldata,
        extrapolation_bc = Line()
    )
    velocity_interp = LinearInterpolation(
        timedata, velocitydata,
        extrapolation_bc = Line()
    )

    # Define the model.
    # - t: The time data.
    # - y: An array of parameter values
    function model(t_, p_)  # t_ is an array of time points
        # Define the model
        m = LuGreModel(p_[1], p_[2], p_[3], p_[4], p_[5], p_[6])
        
        # Define functions to compute the velocity and normal force at the
        # appropriate time.  This routine will utilize the interpolation
        # objects previously defined.
        vel(ti) = velocity_interp(ti)
        nrm(ti) = normal_interp(ti)

        # Solve the model
        rsp = friction(m, t_, vel, nrm)
        return rsp.force
    end

    # Fit the model
    p0 = [
        mdl.static_coefficient,
        mdl.coulomb_coefficient,
        mdl.stribeck_velocity,
        mdl.bristle_stiffness,
        mdl.bristle_damping,
        mdl.viscous_damping
    ]
    fit = curve_fit(model, timedata, frictiondata, p0)
    p = coef(fit)
    return (model = LuGreModel(p[1], p[2], p[3], p[4], p[5], p[6]), fit = fit)
end
