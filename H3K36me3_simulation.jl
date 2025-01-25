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

t = default_t()

histone_nr = 20
gene1_start = 5
gene1_end = 15

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


@named auto_methylation = ReactionSystem(rxs, t)
auto_methylation = complete(auto_methylation)


param = [
:rna_pol_bind => 0.05,
:setd2_bind => 0.40,
:rna_pol_diss => 0.40,
:setd2_diss => 0.005,
:met_rate => 0.80,
:rna_pol_sliding => 0.20,
:demet_rate => 0.015
]

u0 = vcat(
    [um[i] => 1 for i in 1:histone_nr],
    [m[i] => 2 for i in 1:histone_nr],
    [rna_pol[i] => 3 for i in 1:histone_nr],
    [setd2[i] => 4 for i in 1:histone_nr])

tspan = (0.0, 0.1)
jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput)
sol = solve(jprob, JumpProcesses.SSAStepper())
start_pos = sol.u[1]
um_index = findall(x -> x == 1, start_pos)
m_index = findall(x -> x == 2, start_pos)
rna_pol_index = findall(x -> x == 3, start_pos)
setd2_index = findall(x -> x == 4, start_pos)

u0 = vcat(
    [um[i] => 1 for i in 1:histone_nr],
    [m[i] => 0 for i in 1:histone_nr],
    [rna_pol[i] => 0 for i in 1:histone_nr],
    [setd2[i] => 0 for i in 1:histone_nr])

tspan = (0.0, 2000)

jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput)
sol = solve(jprob, JumpProcesses.SSAStepper())

changes_over_time = transpose(hcat(sol.u...))

um_index_over_time = changes_over_time[:, um_index]
m_index_over_time = changes_over_time[:, m_index]
rna_pol_index_over_time = changes_over_time[:, rna_pol_index]
setd2_index_over_time = changes_over_time[:, setd2_index]

print(start_pos)

t_max = 1:300

#plot(sol.t, m_index_over_time)


# heatmap(
#     m_index_over_time[t_max,:] .|
#     (rna_pol_index_over_time[t_max,:] .<< 1) .|
#     (setd2_index_over_time[t_max,:] .<< 2),
#     color = cgrad([:lightblue,  # Non-methylated    |   Nothing bound
#             :lightgreen,        # Methylated        |   Nothing bound
#             :orange,            # Non-methylated    |   RNA POLII bound
#             :olive,             # Methylated        |   RNA POLII bound
#             :black,             # Non-methyalted    |   SETD2 bound (NOT possible)
#             :black,             # Methylated        |   SETD2 bound (NOT possible)
#             :red,               # Non-methylated    |   RNA POLII + SETD2 bound
#             :darkgoldenrod      # Methylated        |   RNA POLII + SETD2 bound
#             ])
# )

heatmap(
    m_index_over_time[t_max,:] .|
    (rna_pol_index_over_time[t_max,:] .<< 1) .|
    (setd2_index_over_time[t_max,:] .<< 2),
    color = cgrad([:lightblue,  # Non-methylated    |   Nothing bound
            :olivedrab,         # Methylated        |   Nothing bound
            :sienna2,           # Non-methylated    |   RNA POLII bound
            :sienna2,           # Methylated        |   RNA POLII bound
            :black,             # Non-methyalted    |   SETD2 bound (NOT possible)
            :black,             # Methylated        |   SETD2 bound (NOT possible)
            :brown3,            # Non-methylated    |   RNA POLII + SETD2 bound
            :brown3             # Methylated        |   RNA POLII + SETD2 bound
            ])
)


m_index_over_time


# bitwise comparison of binary matrices
# t1 = [0 1 0 1 ; 0 1 0 1]
# t2 = [0 0 1 1 ; 0 0 1 1]
# t3 = [0 0 0 0 ; 1 1 1 1]

# t1 .| (t2 .<< 1) .| (t3 .<< 2)
