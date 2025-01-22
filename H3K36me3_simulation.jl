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

# unmethylated
@species um(t)[1:histone_nr]

# unmethylated + RNA POLII bound
@species um_rna(t)[1:histone_nr]

# unmethylated + RNA POL II bound + SETD2 bound
@species um_rna_setd2(t)[1:histone_nr]

# unmethylated + kdm4a bound
@species um_kdm4a(t)[1:histone_nr]

# methylated
@species m(t)[1:histone_nr]

# methylated + RNA POLII bound
@species m_rna(t)[1:histone_nr]

# methylated + RNA POL II bound + SETD2 bound
@species m_rna_setd2(t)[1:histone_nr]

# methylated + KDM4A
@species m_kdm4a(t)[1:histone_nr]


@parameters rna_poll_bind, setd2_bind, kdm4a_bind
@parameters rna_poll_diss, setd2_diss, rna_poll_setd2_diss, kdm4a_diss 
@parameters met_rate, rna_pol_sliding, unmet_rate

rxs = []


######## Later add restriction to be binary only ########
# It's probably fine like it is
# ifelse($m_rna[gene1_start] == 1 | $m_rna_setd2[gene1_start] == 1, 0, 1)
########## end ###############


# Add binding of rna pol 2
append!(rxs, [(@reaction rna_poll_bind, $(m[gene1_start]) --> $(m_rna[gene1_start]))])
append!(rxs, [(@reaction rna_poll_bind, $(um[gene1_start]) --> $(um_rna[gene1_start]))])
# Add binding of setd2
append!(rxs, [(@reaction setd2_bind, $(m_rna[i]) --> $(m_rna_setd2[i])) for i in gene1_start:gene1_end])
append!(rxs, [(@reaction setd2_bind, $(um_rna[i]) --> $(um_rna_setd2[i])) for i in gene1_start:gene1_end])
# Add binding of kdm4a
append!(rxs, [(@reaction kdm4a_bind, $(m[i]) --> $(m_kdm4a[i])) for i in gene1_start:gene1_end])

# Add sliding of rna pol 2
append!(rxs, [@reaction rna_pol_sliding, $(m_rna[i]) --> $(m_rna[i+1]) for i in gene1_start:gene1_end-1])
append!(rxs, [@reaction rna_pol_sliding, $(um_rna[i]) --> $(um_rna[i+1]) for i in gene1_start:gene1_end-1])
# Add sliding of rna pol 2 + setd2
append!(rxs, [@reaction rna_pol_sliding, $(m_rna_setd2[i]) --> $(m_rna_setd2[i+1]) for i in gene1_start:gene1_end-1])
append!(rxs, [@reaction rna_pol_sliding, $(um_rna_setd2[i]) --> $(um_rna_setd2[i+1]) for i in gene1_start:gene1_end-1])

# Add setd2 methylation
append!(rxs, [(@reaction met_rate, $(um_rna_setd2[i]) --> $(m_rna_setd2[i])) for i in gene1_start:gene1_end])
# Add kdm4a demethylation
append!(rxs, [(@reaction unmet_rate, $(m_kdm4a[i]) --> $(um_kdm4a[i])) for i in gene1_start:gene1_end])

# Add dissociation of rna pol 2
append!(rxs, [(@reaction rna_poll_diss, $(m_rna[gene1_end]) --> $(m[gene1_end]))])
append!(rxs, [(@reaction rna_poll_diss, $(um_rna[gene1_end]) --> $(um[gene1_end]))])
# Add dissociation of setd2
append!(rxs, [(@reaction setd2_diss, $(m_rna_setd2[i]) --> $(m_rna[i])) for i in gene1_start:gene1_end])
append!(rxs, [(@reaction setd2_diss, $(um_rna_setd2[i]) --> $(um_rna[i])) for i in gene1_start:gene1_end])
# Add dissociation of rna pol 2 + setd2
append!(rxs, [(@reaction rna_poll_setd2_diss, $(m_rna_setd2[gene1_end]) --> $(m[gene1_end]))])
append!(rxs, [(@reaction rna_poll_setd2_diss, $(um_rna_setd2[gene1_end]) --> $(um[gene1_end]))])
# Add dissociation of kdm4a
append!(rxs, [(@reaction kdm4a_diss, $(um_kdm4a[i]) --> $(um[i])) for i in gene1_start:gene1_end])

# I feel like rna pol 2 is able to slide on a nucleosome while kdm4a is bound. And kdm4a is able to bind while rna pol 2 is bound
# And two rna pol 2 can be on top of each other.


@named auto_methylation = ReactionSystem(rxs, t)
auto_methylation = complete(auto_methylation)


param = [
:rna_poll_bind => 0.20,
:setd2_bind => 0.65,
:kdm4a_bind => 0.10,
:rna_poll_diss => 0.65,
:setd2_diss => 0.05,
:rna_poll_setd2_diss => 0.65, 
:kdm4a_diss => 0.7, 
:met_rate => 0.9,
:rna_pol_sliding => 0.35,
:unmet_rate => 0.4
]

# u0 = vcat(
#     [um[i] => 1 for i in 1:histone_nr],    
#     [um_rna[i] => 2 for i in 1:histone_nr],
#     [um_rna_setd2[i] => 3 for i in 1:histone_nr],
#     [um_kdm4a[i] => 4 for i in 1:histone_nr],
#     [m[i] => 5 for i in 1:histone_nr],
#     [m_rna[i] => 6 for i in 1:histone_nr],
#     [m_rna_setd2[i] => 7 for i in 1:histone_nr],
#     [m_kdm4a[i] => 8 for i in 1:histone_nr])
# tspan = (0.0, 0.1)
# jinput = JumpInputs(auto_methylation, u0, tspan, param)
# jprob = JumpProcesses.JumpProblem(jinput)
# sol = solve(jprob, JumpProcesses.SSAStepper())
# start_pos = sol.u[1]
# um_index = findall(x -> x == 1, start_pos)
# um_rna_index = findall(x -> x == 2, start_pos)
# um_rna_setd2_index = findall(x -> x == 3, start_pos)
# um_kdm4a_index = findall(x -> x == 4, start_pos)
# m_index = findall(x -> x == 5, start_pos)
# m_rna_index = findall(x -> x == 6, start_pos)
# m_rna_setd2_index = findall(x -> x == 7, start_pos)
# m_kdm4a_index = findall(x -> x == 8, start_pos)

u0 = vcat(
    [um[i] => 1 for i in 1:histone_nr],    
    [um_rna[i] => 0 for i in 1:histone_nr],
    [um_rna_setd2[i] => 0 for i in 1:histone_nr],
    [um_kdm4a[i] => 0 for i in 1:histone_nr],
    [m[i] => 0 for i in 1:histone_nr],
    [m_rna[i] => 0 for i in 1:histone_nr],
    [m_rna_setd2[i] => 0 for i in 1:histone_nr],
    [m_kdm4a[i] => 0 for i in 1:histone_nr])

tspan = (0.0, 1000)

jinput = JumpInputs(auto_methylation, u0, tspan, param)
jprob = JumpProcesses.JumpProblem(jinput)
sol = solve(jprob, JumpProcesses.SSAStepper())

changes_over_time = transpose(hcat(sol.u...))

um_over_time = changes_over_time[:, um_index]
um_rna_over_time = changes_over_time[:, um_rna_index]
um_rna_setd2_over_time = changes_over_time[:, um_rna_setd2_index]
um_kdm4a_over_time = changes_over_time[:, um_kdm4a_index]
m_over_time = changes_over_time[:, m_index]
m_rna_over_time = changes_over_time[:, m_rna_index]
m_rna_setd2_over_time = changes_over_time[:, m_rna_setd2_index]
m_kdm4a_over_time = changes_over_time[:, m_kdm4a_index]

print(start_pos)

