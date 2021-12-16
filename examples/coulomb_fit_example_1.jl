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
init_guess = CoulombModel(
    0.25
)

# Fit the model
results = fit_model(
    init_guess,     # Initial Guess
    data[:,1],      # Time
    data[:,4],      # Friction Force
    data[:,3],      # Normal Force
    data[:,2]       # Velocity
)

# Solve the model
F = friction.(results.model, data[:,3], data[:,2])

# Plot the data
plt = plot(
    data[:,2], F,
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
    "coefficient: %f\n\tError: %f\n\tConfidence Interval:: (%f, %f)\n",
    results.model.coefficient,
    sigma[1],
    confidence[1][1],
    confidence[1][2]
)
