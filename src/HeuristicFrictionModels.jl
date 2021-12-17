# HeuristicFrictionModels.jl

# Notes:
# - For the routines in this section it is expected that any heuristic model 
#   provided implement two routines.  These routines are as follows:
#   #1: dzdt = heuristic_state_equation(mdl, t, nrm, pos, vel, z)
#   #2: F = heuristic_force_equation(mdl, t, nrm, pos, vel, z, dzdt)
#
# The inputs to these functions are as follows:
#   - mdl: A type based upon HeuristicFrictionModel
#   - t: A numeric value defining the time point at which to evaluate the model.
#   - nrm: The compressive normal force.
#   - pos: The relative position of the two contacting bodies at time t.
#   - vel: The relative (sliding) velocity between the two contacting bodies.
#   - z: An N-element state vector at time t.
#   - dzdt: An N-element vector of the state variable time derivatives.

"""
Evaluates a heuristic-type friction model at the given state.

The output of the function is a named tuple containing the following fields.
- f: The friction force.
- dzdt: The time derivative(s) of the state variable(s).
"""
function friction(
    mdl::HeuristicFrictionModel,
    t::T,
    nrm::T,
    pos::T,
    vel::T,
    z::Array{T}
) where T <: Number

    # Evaluate the state equation of the heuristic system
    heuristic_state_equation!(mdl, t, nrm, pos, vel, z, dzdt)

    # Evaluate the friction force of the heuristic system
    F = heuristic_force_equation(mdl, t, nrm, pos, vel, z, dzdt)

    # Output
    return (f = F, dzdt = dzdt)
end

"""
Evaluates a heuristic-type friction model by integrating the differential 
equations describing the model state.  The user can specify that the solver 
choose the time step by only supplying a start and end simulation time, or the 
user can specify the time points the solver returns by specifying more than 2 
time points.  Notice, if supplying time points at which the solver should 
return the results, the time points do not need to be evenly spaced; however, 
they do need to be either monotonically increasing or decreasing.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'p' input parameter.  The following
inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.

The output of the function is a named tuple containing the following fields.
- t: The solution time points.
- f: The friction force.
- z: The state variable(s).
- dzdt: The time derivative(s) of the state variable(s).
"""
function friction(
    mdl::HeuristicFrictionModel,
    t::Array{T},
    nrm::Function,
    pos::Function,
    vel::Function,
    z0::Array{T};
    p...
) where T <: Number
    
    function diffeq!(du_, u_, p_, t_)
        x = pos(t_)
        v = vel(t_)
        N = nrm(t_)
        heuristic_state_equation!(mdl, t_, N, x, v, u_, du_)
    end

    args = extract_options(p)

    prob = ODEProblem(diffeq!, z0, [first(t), last(t)])

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

    z = Array(sol)
    npts = length(sol.t)
    nargs = length(z0)
    dzdt = zeros(T, nargs, npts)
    F = zeros(T, npts)
    for i in 1:npts
        x = pos(sol.t[i])
        v = vel(sol.t[i])
        N = nrm(sol.t[i])
        heuristic_state_equation!(mdl, sol.t[i], N, x, v, z[:,i], dzdt[:,i])
        F[i] = heuristic_force_equation(mdl, sol.t[i], N, x, v, z[:,i], dzdt[:,i])
    end

    return (t = sol.t, f = F, z = z, dzdt = dzdt)
end

"""
Evaluates a heuristic-type friction model, that is not dependent upon position, 
by integrating the differential equations describing the model state.  The user 
can specify that the solver choose the time step by only supplying a start and 
end simulation time, or the user can specify the time points the solver returns 
by specifying more than 2 time points.  Notice, if supplying time points at 
which the solver should return the results, the time points do not need to be 
evenly spaced; however, they do need to be either monotonically increasing or 
decreasing.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'p' input parameter.  The following
inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.

The output of the function is a named tuple containing the following fields.
- t: The solution time points.
- f: The friction force.
- z: The state variable(s).
- dzdt: The time derivative(s) of the state variable(s).
"""
function friction(
    mdl::HeuristicFrictionModel,
    t::Array{T},
    nrm::Function,
    vel::Function,
    z0::Array{T};
    p...
) where T <: Number

    pos(ti) = zero(T)
    return friction(mdl, t, nrm, pos, vel, z0; p...)
end

"""
Evaluates a heuristic friction model by integrating the state equations using 
the supplied data.  If the solver requires points lying between provided data 
points linear interpolation will be used.  The user can specify that the solver 
choose the time step by only supplying a start and end simulation time, or the 
user can specify the time points the solver returns by specifying more than 2 
time points.  Notice, if supplying time points at which the solver should return 
the results, the time points do not need to be evenly spaced; however, they do 
need to be either monotonically increasing or decreasing.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'p' input parameter.  The following
inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.

The output of the function is a named tuple containing the following fields.
- t: The solution time points.
- f: The friction force.
- z: The state variable(s).
- dzdt: The time derivative(s) of the state variable(s).
"""
function friction(
    mdl::HeuristicFrictionModel,
    t::Array{T},
    nrm::Array{T},
    pos::Array{T},
    vel::Array{T},
    z0::Array{T};
    p...
) where T <: Number

    normal_interp = LinearInterpolation(
        t, nrm,
        extrapolation_bc = Line()
    )
    position_interp = LinearInterpolation(
        t, pos,
        extrapolation_bc = Line()
    )
    velocity_interp = LinearInterpolation(
        t, vel,
        extrapolation_bc = Line()
    )

    N(ti) = normal_interp(ti)
    X(ti) = position_interp(ti)
    V(ti) = velocity_interp(ti)

    return friction(mdl, t, N, X, V, z0; p...)
end

"""
Evaluates a heuristic friction model that doesn't depend upon position data by 
integrating the state equations using the supplied data.  If the solver requires
points lying between provided data points linear interpolation will be used.  
The user can specify that the solver choose the time step by only supplying a 
start and end simulation time, or the user can specify the time points the 
solver returns by specifying more than 2 time points.  Notice, if supplying time
points at which the solver should return the results, the time points do not 
need to be evenly spaced; however, they do need to be either monotonically 
increasing or decreasing.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'p' input parameter.  The following
inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.

The output of the function is a named tuple containing the following fields.
- t: The solution time points.
- f: The friction force.
- z: The state variable(s).
- dzdt: The time derivative(s) of the state variable(s).
"""
function friction(
    mdl::HeuristicFrictionModel,
    t::Array{T},
    nrm::Array{T},
    vel::Array{T},
    z0::Array{T};
    p...
) where T <: Number

    normal_interp = LinearInterpolation(
        t, nrm,
        extrapolation_bc = Line()
    )
    velocity_interp = LinearInterpolation(
        t, vel,
        extrapolation_bc = Line()
    )

    N(ti) = normal_interp(ti)
    X(ti) = zero(T)
    V(ti) = velocity_interp(ti)

    return friction(mdl, t, N, X, V, z0; p...)
end
