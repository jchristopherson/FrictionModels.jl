# FrictionModels
FrictionModels is a library that provides various friction models suitable for use in mechanical simulations.  The library also provide a mechanism for fitting of the models to existing data.

[![Build Status](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jchristopherson/FrictionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jchristopherson/FrictionModels.jl)

## Basic Usage:
The following example illustrates the use of the Lu-Gre model by computing and plotting the force-velocity relationship.  This example utilizes the Lu-Gre model; however, the same approach may be utilized for any of the other models by simply declaring the appropriate model type.
```julia
using Plots
using LaTeXStrings
using FrictionModels

# Develop a function describing the relative velocity between the two bodies.
velocity(t) = 1.5 * sin(10.0 * pi * t)

# Develop a function describing the normal force between the two bodies.
normal(t) = 100.0

# Define the friction model
mdl = LuGreModel(
    0.25,
    0.15,
    0.01,
    1.0e6,
    6.1e3,
    0.0
)

# Compute the friction model over a segment of time
tspan = [0.0, 1.0]
rsp = friction(mdl, tspan, normal, velocity)

# Plot the results
plot(
    velocity.(rsp.t), 
    rsp.force ./ (mdl.static_coefficient .* normal.(rsp.t)),
    xlabel = L"v(t)",
    ylabel = L"\frac{F}{\mu_s N}",
    label = false
)
```
![](images/lu_gre_example_plot_1.png?raw=true)

## Model Fitting
The following example illustrates fitting data stored in a CSV file to a Lu-Gre model.
```julia
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
```
```text
Fitted Model Coefficients:
static_coefficient: 0.200000
        Error: 0.012948
        Confidence Interval: (0.174464, 0.225536)
coulomb_coefficient: 0.146012
        Error: 0.004643
        Confidence Interval: (0.136855, 0.155170)
stribeck_velocity: 0.152533
        Error: 0.002468
        Confidence Interval: (0.147666, 0.157401)
bristle_stiffness: 4.277057e+04
        Error: 9.193762e+03
        Confidence Interval: (2.463859e+04, 6.090254e+04)
bristle_damping: 3.109026e+02
        Error: 1.354454e+01
        Confidence Interval: (2.841900e+02, 3.376152e+02)
viscous_damping: 0.350666
        Error: 0.411170
        Confidence Interval: (-0.460245, 1.161576)
```
![](images/lu_gre_fit_example_plot_1.png?raw=true)

For reference, the CSV file looked like this:
```csv
Time,Velocity,Normal,Friction
0,0,100.0,0.417029813
0.001,0.047116139,100.0,12.19570558
0.002,0.094185779,100.0,23.82451708
0.003,0.14116247,100.0,15.38126635
...
0.2,6.29379E-15,100.0,0.424359046
```

## References
- Quinn, D.D., "A New Regularization of Coulomb Friction." Journal of Vibration and Acoustics, 126, 2004, 391-397.
- McMillan, A.J., "A Non-Linear Friction Model for Self-Excited Vibrations." Journal of Sound and Vibration, 205(3), 1997, 323-335.
- Lopez, I.; Busturia, J.M.; Nijmeijer, H., "Energy Dissipation of a Friction Damper." Journal of Sound and Vibration, 278, 2004, 539-561.
- Dupont, P.; Armstrong, B.; Hayward, V., "Elasto-Plastic Friction Model: Contact Compliance and Stiction."Proceedinge 2000 American Control Conference, 2000.
- Blaha, P. (2001). "Coulomb Friction Identification Using Harmonic Balance Method of Two-Relay System." PhD Thesis, BRNO University of Technology.
- Rodriguez, E.D.; Garcia, B.S.; Cortes, F.R., (2020). "New Friction Model for Manipulator Robots based on Hyperbolic Functions." 19th National Mechantronics Congress, Queretaro, Mexico, 2020.

