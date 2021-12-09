using Plots
using CSV
using Printf
using DataFrames
using LaTeXStrings
using Interpolations    # Required for the interpolation routines
using LsqFit            # Required for the statistics on the fitting process
using FrictionModels

# Read in a CSV file containing the data to fit.
#
# The file is structured as follows:
# Column 1: Time
# Column 2: Velocity
# Column 3: Normal force
# Column 4: Friction Force
println(pwd())
data = CSV.read("examples//test_data_1.csv", DataFrame)

# Construct an initial estimate for model parameters
init_guess = LuGreModel(
    0.25,
    0.15,
    0.01,
    1.0e6,
    6.0e3,
    0.0
)

# We can specify limits on each parameter.  If no limit is desired for A
# specific parameter, we can simply input -Inf or Inf for either a lower or
# upper constraint respectively.
lb = [0.2, 0.05, 0.0, 1e3, 0.0, 0.0]
ub = [1.0, 0.2, Inf, Inf, Inf, 1.0]

# Fit the model
results = fit_model(
    init_guess,     # Initial Guess
    data[:,1],      # Time
    data[:,4],      # Friction Force
    data[:,3],      # Normal Force
    data[:,2],      # Velocity,
    lower = lb,     # Lower Bounds (optional)
    upper = ub      # Upper Bounds (optional)
)

# Evaluate the model and plot against the data.  Use linear interpolation to 
# allow the solver to estimate normal force and velocity values at time points
# not aligned with the measured values.
normal_interp = LinearInterpolation(
    data[:,1], data[:,3],
    extrapolation_bc = Line()
)
velocity_interp = LinearInterpolation(
    data[:,1], data[:,2],
    extrapolation_bc = Line()
)

# The friction routine requires a function to describe the velocity and normal
# force; therefore, define functions that use the interpolation objects.
nrm(ti) = normal_interp(ti)
vel(ti) = velocity_interp(ti)

# Solve the model
tspan = [first(data[:,1]), last(data[:,1])]
rsp = friction(results.model, tspan, nrm, vel)

# Plot the data
plt = plot(
    vel.(rsp.t), rsp.f,
    xlabel = L"v(t)",
    ylabel = L"F(t)",
    label = "Fitted Model",
    lw = 2,
    legend = :bottomright
)
plot!(
    plt,
    data[:,2], data[:,4],
    st = :scatter,
    label = "Measured Data"
)
display(plt)

# Print out the coefficients of the fitted model
@printf("Fitted Model Coefficients:\n")

# Print out statistitics related to the fit
sigma = stderror(results.fit)
confidence = confidence_interval(results.fit)

# Print out the results
@printf(
    "static_coefficient: %f\n\tError: %f\n\tConfidence Interval: (%f, %f)\n",
    results.model.static_coefficient,
    sigma[1],
    confidence[1][1],
    confidence[1][2]
)
@printf(
    "coulomb_coefficient: %f\n\tError: %f\n\tConfidence Interval: (%f, %f)\n",
    results.model.coulomb_coefficient,
    sigma[2],
    confidence[2][1],
    confidence[2][2]
)
@printf(
    "stribeck_velocity: %f\n\tError: %f\n\tConfidence Interval: (%f, %f)\n",
    results.model.stribeck_velocity,
    sigma[3],
    confidence[3][1],
    confidence[3][2]
)
@printf(
    "bristle_stiffness: %e\n\tError: %e\n\tConfidence Interval: (%e, %e)\n",
    results.model.bristle_stiffness,
    sigma[4],
    confidence[4][1],
    confidence[4][2]
)
@printf(
    "bristle_damping: %e\n\tError: %e\n\tConfidence Interval: (%e, %e)\n",
    results.model.bristle_damping,
    sigma[5],
    confidence[5][1],
    confidence[5][2]
)
@printf(
    "viscous_damping: %f\n\tError: %f\n\tConfidence Interval: (%f, %f)\n",
    results.model.viscous_damping,
    sigma[6],
    confidence[6][1],
    confidence[6][2]
)
