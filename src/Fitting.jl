"""
Uses a Levenberg-Marquardt solver to compute the best fit of a friction model to
a predefined data set.  As the supplied data is discretely sampled but the 
solver may require points between those supplied, the solver will utilize linear
interpolation to estimate values between supplied data points.

The routine returns the fitted model along with the LsqFitResults type returned
from the Levenberg-Marquardt solver.  The LsqFitResults allow exploration of 
error margins, confidence intervals, etc.  Both results are returned as a 
named tuple with the fitted model available as 'model' and the LsqFitResults
available as 'fit'.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'pargs' input parameter.  The 
following inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.
- lower: An N-element array of lower bounds on each of the N coefficients.
- upper: An N-element array of upper bounds on each of the N coefficients.
"""
function fit_model(
    mdl::HeuristicFrictionModel,
    timedata::Array{T},
    frictiondata::Array{T},
    normaldata::Array{T},
    positiondata::Array{T},
    velocitydata::Array{T},
    z0::Array{T};
    pargs...
) where T <: Number

    # Set up the interpolation for each data set
    normal_interp = LinearInterpolation(
        timedata, normaldata,
        extrapolation_bc = Line()
    )
    position_interp = LinearInterpolation(
        timedata, positiondata,
        extrapolation_bc = Line()
    )
    velocity_interp = LinearInterpolation(
        timedata, velocitydata,
        extrapolation_bc = Line()
    )

    # Define the model.
    # - t_: The time data.
    # - p_: An array of parameter values
    function model(t_, p_) # t_ is an array of time points
        # Define the model
        m = model_from_array(mdl, p_)

        # Define functions to compute the velocity and normal force at the
        # appropriate time.  This routine will utilize the interpolation
        # objects previously defined.
        pos(ti) = position_interp(ti)
        vel(ti) = velocity_interp(ti)
        nrm(ti) = normal_interp(ti)

        # Solve the model
        rsp = friction(m, t_, nrm, pos, vel, z0; pargs)
        return rsp.f
    end

    # Convert the model to array form
    p0 = model_to_array(mdl)

    # Extract the model options
    opt = extract_options(pargs, length(p0))

    # Fit the model
    fit = curve_fit(
        model,
        timedata,
        frictiondata,
        p0,
        lower = opt.lower,
        upper = opt.upper
    )
    p = coef(fit)

    return (model = model_from_array(mdl, p), fit = fit)
end

"""
Uses a Levenberg-Marquardt solver to compute the best fit of a friction model to
a predefined data set.  As the supplied data is discretely sampled but the 
solver may require points between those supplied, the solver will utilize linear
interpolation to estimate values between supplied data points.  This instance
is used when the model does not depend upon position data.

The routine returns the fitted model along with the LsqFitResults type returned
from the Levenberg-Marquardt solver.  The LsqFitResults allow exploration of 
error margins, confidence intervals, etc.  Both results are returned as a 
named tuple with the fitted model available as 'model' and the LsqFitResults
available as 'fit'.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'pargs' input parameter.  The 
following inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.
- lower: An N-element array of lower bounds on each of the N coefficients.
- upper: An N-element array of upper bounds on each of the N coefficients.
"""
function fit_model(
    mdl::HeuristicFrictionModel,
    timedata::Array{T},
    frictiondata::Array{T},
    normaldata::Array{T},
    velocitydata::Array{T},
    z0::Array{T};
    pargs...
) where T <: Number

    pos = zeros(T, length(timedata))
    return fit_model(mdl, timedata, frictiondata, normaldata, pos, velocitydata, z0; pargs...)
end

"""
Uses a Levenberg-Marquardt solver to compute the best fit of a friction model to
a predefined data set.  As the supplied data is discretely sampled but the 
solver may require points between those supplied, the solver will utilize linear
interpolation to estimate values between supplied data points.

The routine returns the fitted model along with the LsqFitResults type returned
from the Levenberg-Marquardt solver.  The LsqFitResults allow exploration of 
error margins, confidence intervals, etc.  Both results are returned as a 
named tuple with the fitted model available as 'model' and the LsqFitResults
available as 'fit'.

Control over the ODE solver tolerances and maximum allowable time step are
exposed to the calling code via the varargs 'pargs' input parameter.  The 
following inputs are currently recognized.
- reltol: The relative ODE solver tolerance. The defaults is 1e-8.
- abstol: The absolute ODE solver tolerance.  The default is 1e-6.
- dtmax: The maximum allowable step size.  The default is 1e-3.
- lower: An N-element array of lower bounds on each of the N coefficients.
- upper: An N-element array of upper bounds on each of the N coefficients.
"""
function fit_model(
    mdl::FrictionModel,
    timedata::Array{T},
    frictiondata::Array{T},
    normaldata::Array{T},
    velocitydata::Array{T};
    pargs...
) where T <: Number

    normal_interp = LinearInterpolation(
        timedata, normaldata,
        extrapolation_bc = Line()
    )
    velocity_interp = LinearInterpolation(
        timedata, velocitydata,
        extrapolation_bc = Line()
    )

    function model(t_, p_)
        m = model_from_array(mdl, p_)
        return friction.(m, normal_interp.(t_), velocity_interp.(t_))
    end

    p0 = model_to_array(mdl)
    opt = extract_options(pargs, length(p0))

    fit = curve_fit(
        model,
        timedata,
        frictiondata,
        p0,
        lower = opt.lower,
        upper = opt.upper
    )
    p = coef(fit)
    return (model = model_from_array(mdl, p), fit = fit)
end
