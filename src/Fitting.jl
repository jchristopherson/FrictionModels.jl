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
- zi: The initial condition or initial bristle deformation.  The default is 0.
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
    function model(t_, p_) # t_ is an array of time points
        # Define the model
        m = model_from_array(mdl, p_)

        # Define functions to compute the velocity and normal force at the
        # appropriate time.  This routine will utilize the interpolation
        # objects previously defined.
        vel(ti) = velocity_interp(ti)
        nrm(ti) = normal_interp(ti)

        # Solve the model
        rsp = friction(m, t_, nrm, vel)
        return rsp.f
    end

    # Extract the model options
    opt = extract_options(pargs)

    # Fit the model
    p0 = model_to_array(mdl)
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