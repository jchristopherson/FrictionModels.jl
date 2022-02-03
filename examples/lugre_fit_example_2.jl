using Plots
using CSV
using DataFrames
using LaTeXStrings
using FrictionModels

# Read in the data file containing the data to fit.
#
# The file is structured as follows.
# Column 1: Time [s]
# Column 2: Normal Force [uN]
# Column 3: Position [mm]
# Column 4: Friction Force [uN]
data = CSV.read("examples//test_data_2.txt", DataFrame)

# Extract data
t_data = data[:,1]
N_data = data[:,2]
x_data = data[:,3]
F_data = data[:,4]

# Tare position data
tare = x_data[1]
x_data = x_data .- tare

# Change the sign on the friction force data as the experimental data was 
# obtained using a different convention than used by the model.
F_data = -F_data

# Compute velocity via finite differences
n = length(t_data)
v_data = Vector{Float64}(undef, n)
v_data[1] = (x_data[2] - x_data[1]) / (t_data[2] - t_data[1])
v_data[n] = (x_data[n] - x_data[n-1]) / (t_data[n] - t_data[n-1])
for i in 2:n-1
    v_data[i] = (x_data[i+1] - x_data[i-1]) / (t_data[i+1] - t_data[i-1])
end

# Construct an initial estimate for model parameters
init_guess = LuGreModel(
    0.65,
    0.15,
    0.01,
    1e6,
    1e7,
    0.0
)

# Define limits on each parameters
lb = [0.1, 0.1, 0.0, 0.0, 0.0, 0.0]
ub = [1.0, 1.0, Inf, Inf, Inf, 0.01]

# Fit the model
results = fit_model(
    init_guess,
    t_data,
    F_data,
    N_data,
    v_data,
    [0.0],
    lower = lb,
    upper = ub
)

# Solve the model
rsp = friction(results.model, t_data, N_data, v_data, [0.0])

# Plot the force-displacement results
plt = plot(
    1e3 .* x_data, rsp.f ./ N_data,
    xlabel = L"x(t)",
    ylabel = L"\mu",
    label = "Fitted Model",
    lw = 2,
    ylims = (-1.0, 1.0)
)
plot!(
    plt,
    1e3 .* x_data, F_data ./ N_data,
    ls = :dot,
    label = "Measured Data"
)
display(plt)

# Plot the time history results
nt = 10000
plt_time = plot(
    t_data[1:nt], 1e-6 .* rsp.f[1:nt],
    xlabel = L"t [s]",
    ylabel = L"F [N]",
    label = "Fitted Model",
    lw = 2,
    legend = :bottomleft
)
plot!(
    plt_time,
    t_data[1:nt], 1e-6 .* F_data[1:nt],
    ls = :dot,
    lw = 2,
    label = "Measured Data"
)
display(plt_time)
