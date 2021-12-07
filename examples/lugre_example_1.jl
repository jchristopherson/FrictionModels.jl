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
    rsp.f ./ (mdl.static_coefficient .* normal.(rsp.t)),
    xlabel = L"v(t)",
    ylabel = L"\frac{F}{\mu_s N}",
    label = false
)
