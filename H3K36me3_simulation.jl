using Pkg
using Catalyst
using JumpProcesses
using Plots
using DifferentialEquations
using Optim
using BlackBoxOptim
using Distributions
using Pipe
using StatsBase
using Statistics
using DifferentialEquations.EnsembleAnalysis

t = default_t()

histone_nr = 60
gene1_start = 1
gene1_end = 60

# unmethylated histones
@species um(t)[1:histone_nr]

# methylated histones
@species m(t)[1:histone_nr]

# RNA POLII location(s)
@species rna_pol(t)[1:histone_nr]

# SETD2 location(s)
@species setd2(t)[1:histone_nr]



@parameters rna_pol_bind, setd2_bind
@parameters rna_pol_diss, setd2_diss 
@parameters met_rate, rna_pol_sliding, demet_rate

rxs = []

# Add binding of rna pol 2
append!(rxs, [(@reaction rna_pol_bind * (1-$(rna_pol[gene1_start])), 0 --> $(rna_pol[gene1_start]))])
# Add binding of setd2
append!(rxs, [(@reaction setd2_bind  * $(rna_pol[i]) * (1-$(setd2[i])), 0 --> $(setd2[i])) for i in gene1_start:gene1_end])


# Add sliding of rna pol 2 (no setd2 bound)
append!(rxs, [@reaction rna_pol_sliding * (1-$(setd2[i])) * (1-$(rna_pol[i+1])), $(rna_pol[i]) --> $(rna_pol[i+1]) for i in gene1_start:gene1_end-1])
# Add sliding of rna pol 2 (with setd2 bound)
append!(rxs, [@reaction rna_pol_sliding * $(setd2[i]) * (1-$(rna_pol[i+1])), $(rna_pol[i]) + $(setd2[i]) --> $(rna_pol[i+1]) + $(setd2[i+1]) for i in gene1_start:gene1_end-1])


# Add setd2 methylation
append!(rxs, [(@reaction met_rate * $(um[i]) * $(setd2[i]), $(um[i]) --> $(m[i])) for i in gene1_start:gene1_end])
# Add kdm4a demethylation
append!(rxs, [(@reaction demet_rate * $(m[i]) * (1-$(rna_pol[i])), $(m[i]) --> $(um[i])) for i in gene1_start:gene1_end])


# Add dissociation of rna pol 2 (with setd2 bound)
append!(rxs, [(@reaction rna_pol_diss * (1-$(setd2[gene1_end])), $(rna_pol[gene1_end]) --> 0)])
# Add dissociation of rna pol 2 (no setd2 bound)
append!(rxs, [(@reaction rna_pol_diss * $(setd2[gene1_end]), $(rna_pol[gene1_end]) + $(setd2[gene1_end]) --> 0)])


# Add dissociation of setd2
append!(rxs, [(@reaction setd2_diss, $(setd2[i]) --> 0) for i in gene1_start:gene1_end])

# Model
@named auto_methylation = ReactionSystem(rxs, t)
auto_methylation = complete(auto_methylation)

# Paramters
param = [
:rna_pol_bind => 0.05,
:setd2_bind => 0.033333,
:rna_pol_diss => 0.9,
:setd2_diss => 0.005,
:met_rate => 0.4,
:rna_pol_sliding => 0.6,
:demet_rate => 0.04
]

# Run a very short simulation to determine the indeces of each histone
# Initial state
u0 = vcat(
    [um[i] => 1 for i in 1:histone_nr],
    [m[i] => 2 for i in 1:histone_nr],
    [rna_pol[i] => 3 for i in 1:histone_nr],
    [setd2[i] => 4 for i in 1:histone_nr])

# How long the simulation should run for
tspan = (0.0, 0.1)

# Run the simulation
jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput)
sol = solve(jprob, JumpProcesses.SSAStepper())
start_pos = sol.u[1]

# Find the histone indeces
um_index = findall(x -> x == 1, start_pos)
m_index = findall(x -> x == 2, start_pos)
rna_pol_index = findall(x -> x == 3, start_pos)
setd2_index = findall(x -> x == 4, start_pos)


#################################### Simulate single trajectory ####################################################


# Set inital state of all histones to unmethylated
u0 = vcat(
    [um[i] => 1 for i in 1:histone_nr],
    [m[i] => 0 for i in 1:histone_nr],
    [rna_pol[i] => 0 for i in 1:histone_nr],
    [setd2[i] => 0 for i in 1:histone_nr])

# Run one longer simulation
tspan = (0.0, 1000)

jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput, save_positions=(false, false))
sol = solve(jprob, JumpProcesses.SSAStepper(), saveat=tspan[begin]:1:tspan[end])

# using JLD2
# save_object("123big_run_solo.jld2", sol)
# sol = load_object("big_run_solo.jld2")

# Extract state changes over time
changes_over_time = transpose(hcat(sol.u...))

um_index_over_time = changes_over_time[:, um_index]
m_index_over_time = changes_over_time[:, m_index]
rna_pol_index_over_time = changes_over_time[:, rna_pol_index]
setd2_index_over_time = changes_over_time[:, setd2_index]

# Plot the change of a species over time
plot([mean(m_index_over_time, dims=2) |> vec |> x -> x[1:end], 
        mean(um_index_over_time, dims=2) |> vec |> x -> x[1:end]])

t_max = 1:500

# Calculate the autocorrelation using the methylation change over time for one position
autocorr_start = StatsBase.autocor(getindex.(sol.u[400:end], m_index[1]), 0:100)
autocorr_end = StatsBase.autocor(getindex.(sol.u[400:end], m_index[end]), 0:100)
autocorr_near_end = StatsBase.autocor(getindex.(sol.u[400:end], m_index[end-5]), 0:100)

# Plot the autocorrelation
plot(autocorr_start[1:end])
plot(autocorr_end[1:end])
plot(autocorr_near_end[1:end])

# Plot average methylation over time
plot(sol.t, mean.(getindex.(sol.u, Ref(m_index))))

# Plot a heatmap of the start of the simulation
heatmap(
    m_index_over_time[t_max,:] .|
    (rna_pol_index_over_time[t_max,:] .<< 1) .|
    (setd2_index_over_time[t_max,:] .<< 2),
    color = cgrad([:lightskyblue,  # Non-methylated    |   Nothing bound
            :purple,         # Methylated        |   Nothing bound
            :forestgreen,           # Non-methylated    |   RNA POLII bound
            :forestgreen,           # Methylated        |   RNA POLII bound
            :black,             # Non-methyalted    |   SETD2 bound (NOT possible)
            :black,             # Methylated        |   SETD2 bound (NOT possible)
            :darkorange3,            # Non-methylated    |   RNA POLII + SETD2 bound
            :darkorange3             # Methylated        |   RNA POLII + SETD2 bound
            ]),
    xlabel = "Position",
    ylabel = "Time",
    title = "Simulation over time"
)


# Calculate the crosscorrelation between methylation and RNA POL II at a specific location
crosscorr_start = StatsBase.crosscor(
    getindex.(sol.u[400:end], m_index[1]), 
    getindex.(sol.u[400:end], rna_pol_index[1]),
    -100:100)
crosscorr_start_reverse = StatsBase.crosscor(
    getindex.(sol.u[400:end], rna_pol_index[1]), 
    getindex.(sol.u[400:end], m_index[1]),
    -100:100)
crosscorr_end = StatsBase.crosscor(
    getindex.(sol.u[400:end], m_index[end]), 
    getindex.(sol.u[400:end], rna_pol_index[end]),
    -100:100)
crosscorr_near_end = StatsBase.crosscor(
    getindex.(sol.u[400:end], m_index[end-5]), 
    getindex.(sol.u[400:end], rna_pol_index[end-5]),
    -100:100)

# Plot the crosscorrelation between methylation and RNA POL II at a specific location
plot(-100:100, crosscorr_start, xlabel="lags", ylabel="crosscorrelation")
plot(-100:100, crosscorr_start_reverse, xlabel="lags", ylabel="crosscorrelation")
plot(-100:100, crosscorr_end, xlabel="lags", ylabel="crosscorrelation")
plot(-100:100, crosscorr_near_end, xlabel="lags", ylabel="crosscorrelation")


# bitwise comparison of binary matrices
# t1 = [0 1 0 1 ; 0 1 0 1]
# t2 = [0 0 1 1 ; 0 0 1 1]
# t3 = [0 0 0 0 ; 1 1 1 1]

# t1 .| (t2 .<< 1) .| (t3 .<< 2)


#################################### Simulate multiple trajectories ####################################################

# initialize parameters
param = [
:rna_pol_bind => 0.05,
:setd2_bind => 0.033333,
:rna_pol_diss => 0.9,
:setd2_diss => 0.005,
:met_rate => 0.4,
:rna_pol_sliding => 0.6,
:demet_rate => 0.04
]

# Set simulation run time
tspan = (0.0, 1500)

# Initialize problem
jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput, save_positions=(false, false))

jprob = remake(jprob, p=param)
ens_prob = EnsembleProblem(jprob)
# Set number of trajectories and run the simulations
res = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=100, saveat=tspan[begin]:50:tspan[end]) # Save every 50 steps

# all results
res
# result of one simulation
res[1]
# time steps (same for all simulations)
res[1].t
# state space for 1 time step for 1 simulation
res[1].u[1]
res[1].u[31]

# Plot average of species over time
# State space over time
avg_um_per_sim = [];
avg_m_per_sim = [];
avg_rna_pol_per_sim = [];
avg_setd2_per_sim = [];

std_um_per_sim = [];
std_m_per_sim = [];
std_rna_pol_per_sim = [];
std_setd2_per_sim = [];

for sim in res
    state_change_over_time = hcat(sim.u...) |> transpose;

    # For every simulation get the state over time. Filter by species
    um_index_over_time = state_change_over_time[:, um_index];
    m_index_over_time = state_change_over_time[:, m_index];
    rna_pol_index_over_time = state_change_over_time[:, rna_pol_index];
    setd2_index_over_time = state_change_over_time[:, setd2_index];

    # Save the average of a species at the different time points (for one simulation)
    push!(avg_um_per_sim, mean(um_index_over_time, dims=2) |> vec);
    push!(avg_m_per_sim, mean(m_index_over_time, dims=2) |> vec);
    push!(avg_rna_pol_per_sim, mean(rna_pol_index_over_time, dims=2) |> vec);
    push!(avg_setd2_per_sim, mean(setd2_index_over_time, dims=2) |> vec);

    # Save the standard deviation of a species at the different time points (for one simulation)
    push!(std_um_per_sim, std(um_index_over_time, dims=2) |> vec);
    push!(std_m_per_sim, std(m_index_over_time, dims=2) |> vec);
    push!(std_rna_pol_per_sim, std(rna_pol_index_over_time, dims=2) |> vec);
    push!(std_setd2_per_sim, std(setd2_index_over_time, dims=2) |> vec);
end


# Get the average of a species at the different time points (for all simulations)
# Combine the simulations (one simulation per column)
avg_um_matrix = hcat(avg_um_per_sim...);
avg_m_matrix = hcat(avg_m_per_sim...);
avg_rna_pol_matrix = hcat(avg_rna_pol_per_sim...);
avg_setd2_matrix = hcat(avg_setd2_per_sim...);
# Calculate the average specie number per time point
avg_um = mean(avg_um_matrix, dims=2) |> vec;
avg_m = mean(avg_m_matrix, dims=2) |> vec;
avg_rna_pol = mean(avg_rna_pol_matrix, dims=2) |> vec;
avg_setd2 = mean(avg_setd2_matrix, dims=2) |> vec;

plot(res[1].t,
    [avg_um, avg_m, avg_rna_pol, avg_setd2],
    label = ["Unmethylated" "Methylated" "RNA Pol II" "SETD2"])


# Plot species over time for all individual trajectories
using DifferentialEquations.EnsembleAnalysis

m_series, v_series = timeseries_steps_meanvar(res)
colors = fill(:blue, length(m_series.u[1]))
colors[m_index] .= :green
colors[rna_pol_index] .= :yellow
colors[setd2_index] .= :red

p = plot()
for i in 1:length(m_series.u[1])
    plot!(p, m_series.t, [m_series.u[j][i] for j in 1:length(m_series.u)], color=colors[i], label="")
end

plot!(p, [NaN], [NaN], color=:blue, label="unmethylated")
plot!(p, [NaN], [NaN], color=:green, label="methylated")
plot!(p, [NaN], [NaN], color=:yellow, label="rna pol II")
plot!(p, [NaN], [NaN], color=:red, label="setd2")
display(p)




# Plot the average methylation across the gene of all trajectories (last step)
# Extract the state of the last time step for all trajectories 
last_step = [sim.u[end] for sim in res]
last_step = hcat(last_step...) |> transpose

# Extract the species for all trajectories at the last step
um_last_step = last_step[:, um_index]
m_last_step = last_step[:, m_index]
rna_pol_last_step = last_step[:, rna_pol_index]
setd2_last_step = last_step[:, setd2_index]

# Calculate the average methylation of each histone
avg_m_pos = mean(m_last_step, dims=1)
std_m_pos = mean(m_last_step, dims=1)

# Plot the average methylation
plot(avg_m_pos |> vec, label="Mean methylated", 
    ribbon=std_m_pos)


#################################### Quickly modify and perform multiple runs ####################################################


function quick_run(p, t)
    # Create a problem and run the simulations
    jinput = JumpInputs(auto_methylation, u0, t, p)
    jprob = JumpProcesses.JumpProblem(jinput, save_positions=(false, false))

    jprob = remake(jprob, p=p)
    ens_prob = EnsembleProblem(jprob)
    res = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=2500, saveat=t[begin]:25:t[end])
  
    # Plot species over time
    m_series, _ = timeseries_steps_meanvar(res)
    colors = fill(:blue, length(m_series.u[1]))
    colors[m_index] .= :green
    colors[rna_pol_index] .= :yellow
    colors[setd2_index] .= :red
    
    p = plot()
    for i in 1:length(m_series.u[1])
        plot!(p, m_series.t, [m_series.u[j][i] for j in 1:length(m_series.u)], color=colors[i], label="")
    end
    
    plot!(p, [NaN], [NaN], color=:blue, label="unmethylated")
    plot!(p, [NaN], [NaN], color=:green, label="methylated")
    plot!(p, [NaN], [NaN], color=:yellow, label="rna pol II")
    plot!(p, [NaN], [NaN], color=:red, label="setd2")
    title!("Species over time")
    xlabel!("Time step")
    ylabel!("Ratio histones occupied by species")
    display(p)
    
    # Extract the methylation state at the last time step and plot the average methylation over the gene body
    last_step = [sim.u[end] for sim in res]
    last_step = hcat(last_step...) |> transpose

    m_last_step = last_step[:, m_index]

    avg_m_pos = mean(m_last_step, dims=1)
    std_m_pos = std(m_last_step, dims=1)

    plot(avg_m_pos |> vec, label=false, xlabel="Histone position", ylabel="Methylation intensity", title="Methylation intensity steady state") |> display

    # plot(avg_m_pos |> vec, label="Mean methylated", 
    #     ribbon=std_m_pos)
    
    return res
end

param = [
:rna_pol_bind => 0.05,
:setd2_bind => 0.033333,
:rna_pol_diss => 0.9,
:setd2_diss => 0.005,
:met_rate => 0.4,
:rna_pol_sliding => 0.6,
:demet_rate => 0.04
];

tspan = (0.0, 400);

# Run function that runs the simulation and plots the results
tres2 = quick_run(param, tspan)
