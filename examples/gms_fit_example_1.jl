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
e1 = MaxwellElement(
    1e5,
    1e2,
    0.25
)
e2 = MaxwellElement(
    1e5,
    1e2,
    0.25
)
e3 = MaxwellElement(
    1e5,
    1e2,
    0.25
)
e4 = MaxwellElement(
    1e6,
    1e2,
    0.25
)
init_guess = GeneralizedMaxwellSlipModel(
    [e1, e2, e3, e4],
    0.25,
    0.15,
    0.5,
    0.1,
    0.0
)

# Fit the model
init_condition = [0.0, 0.0, 0.0, 0.0]
results = fit_model(
    init_guess,     # Initial Guess
    data[:,1],      # Time
    data[:,4],      # Friction Force
    data[:,3],      # Normal Force
    data[:,2],      # Velocity
    init_condition  # Initial Condition Vector
)

# Solve the model
rsp = friction(results.model, data[:,1], data[:,3], data[:,2], init_condition)

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
