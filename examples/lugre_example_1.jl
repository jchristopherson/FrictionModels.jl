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
