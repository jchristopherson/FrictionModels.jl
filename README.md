# FrictionModels
FrictionModels is a library that provides various friction models suitable for use in mechanical simulations.  The library also provide a mechanism for fitting of the models to existing data.

[![Build Status](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jchristopherson/FrictionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jchristopherson/FrictionModels.jl)

## Implemented Models:
The following models are currently implemented.
- Coulomb
- LuGre
- Hyperbolic [6]

## Basic Usage:
The following example illustrates the use of the Lu-Gre model by computing and plotting the force-velocity relationship.  This example utilizes the Lu-Gre model, however, the same approach may be utilized for any of the other models by simply declaring the appropriate model type.
```julia
using Plots
using LaTeXStrings
using FrictionModels

# Develop a function describing the relative velocity between the two bodies.
velocity(t) = 1.5 * sin(10.0 * pi * t)

# Develop a function describing the normal force between the two bodies.
normal(t) = 100.0

# A function describing the relative position between two bodies.  Notice, this
# is not needed by the Lu-Gre model, and as such, any function can be defined
# that only takes a single argument.
position(t) = 0.0

# Define the friction model
mdl = LuGreModel(
    0.25,
    0.15,
    0.01,
    1.0e6,
    7.5e2,
    0.0
)

# Compute the friction model over a segment of time
tspan = [0.0, 1.0]
rsp = friction(mdl, tspan, normal, position, velocity, [0.0])

# Plot the results
plot(
    velocity.(rsp.t), 
    rsp.f ./ (mdl.static_coefficient .* normal.(rsp.t)),
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
    data[:,2],      # Velocity
    [0.0],          # Initial Condition Vector
    lower = lb,     # Lower Bounds (optional)
    upper = ub      # Upper Bounds (optional)
)

# Solve the model
tspan = [first(data[:,1]), last(data[:,1])]
rsp = friction(results.model, data[:,1], data[:,3], data[:,2], [0.0])

# Plot the data
plt = plot(
    data[:,2], rsp.f,
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
static_coefficient: 0.200003
        Error: 0.000182
        Confidence Interval: (0.199645, 0.200361)
coulomb_coefficient: 0.132016
        Error: 0.000127
        Confidence Interval: (0.131766, 0.132266)
stribeck_velocity: 0.232961
        Error: 0.000175
        Confidence Interval: (0.232615, 0.233307)
bristle_stiffness: 9.528510e+05
        Error: 1.747849e+02
        Confidence Interval: (9.525062e+05, 9.531957e+05)
bristle_damping: 9.182547e+01
        Error: 7.017297e+01
        Confidence Interval: (-4.656995e+01, 2.302209e+02)
viscous_damping: 1.000000
        Error: 0.163265
        Confidence Interval: (0.678009, 1.321991)
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
1. Quinn, D.D., "A New Regularization of Coulomb Friction." Journal of Vibration and Acoustics, 126, 2004, 391-397.
2. McMillan, A.J., "A Non-Linear Friction Model for Self-Excited Vibrations." Journal of Sound and Vibration, 205(3), 1997, 323-335.
3. Lopez, I., Busturia, J.M., Nijmeijer, H., "Energy Dissipation of a Friction Damper." Journal of Sound and Vibration, 278, 2004, 539-561.
4. Dupont, P., Armstrong, B., Hayward, V., "Elasto-Plastic Friction Model: Contact Compliance and Stiction."Proceedinge 2000 American Control Conference, 2000.
5. Blaha, P. (2001). "Coulomb Friction Identification Using Harmonic Balance Method of Two-Relay System." PhD Thesis, BRNO University of Technology.
6. Rodriguez, E.D., Garcia, B.S., Cortes, F.R., "New Friction Model for Manipulator Robots based on Hyperbolic Functions." 19th National Mechantronics Congress, Queretaro, Mexico, 2020.
7. Al-Bender, F., Lampaert, V., Swevers, J., "Modeling of Dry Sliding Friction Dynamics: From Heuristic Models to Physically Motivated Models and Back." Americal Institute of Physics, 14(2): 445-460, 2004.
8. Al-Bender, F., Lampaert, V., Swevers, J., "A Novel Generic Model at Asperity Level for Dry Friction Force Dynamics." Tribology Letters, 16: 81-93, 2004.
9. Lampaert, V., Swevers, J., Al-Bender, F., "Modification of the Leuven Integrated Friction Model Structure." IEEE Transactions on Automatic Control, 47(4): 683-687, 2002.
10. Al-Bender, F., Lampaert, V., Swevers, J., "The Generalized Maxwell-Slip Model: A Novel Model for Friction Simulation and Compensation." IEEE Transactions on Automatic Control, 50(11):1883-1887, 2005.