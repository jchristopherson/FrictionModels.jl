# FrictionModels

[![Build Status](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jchristopherson/FrictionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jchristopherson/FrictionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jchristopherson/FrictionModels.jl)


This is a work in progress.  Stay tuned as development is active.  Contributions are welcome; however, I'm still piecing together what the API should look like.

# Example 1:
The following example illustrates the use of the Lu-Gre model by computing and plotting the force-velocity relationship.
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
rsp = friction(mdl, tspan, velocity, normal)

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

# References
- Quinn, D.D., "A New Regularization of Coulomb Friction." Journal of Vibration and Acoustics, 126, 2004, 391-397.
- McMillan, A.J., "A Non-Linear Friction Model for Self-Excited Vibrations." Journal of Sound and Vibration, 205(3), 1997, 323-335.
- Lopez, I.; Busturia, J.M.; Nijmeijer, H., "Energy Dissipation of a Friction Damper." Journal of Sound and Vibration, 278, 2004, 539-561.
- Dupont, P.; Armstrong, B.; Hayward, V., "Elasto-Plastic Friction Model: Contact Compliance and Stiction."Proceedinge 2000 American Control Conference, 2000.
- Blaha, P. (2001). "Coulomb Friction Identification Using Harmonic Balance Method of Two-Relay System." PhD Thesis, BRNO University of Technology.
