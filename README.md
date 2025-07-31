# histone_project

This repository contains data and scripts used for the Major Internship of Lars Maas named "A Framework for Mechanistic Stochastic Modelling of H3K36me3 kinetics".

The folder "data_h3k36me3" contains files containing H3K36me3 upregulated genes and regions. The file with H3K36me3 densities across the whole genome is too big for Github, but can be shared upon request.

H3K36me3_simulation.jl is the Julia script used to create the model, run the model and process simulation results.
get_chip_seq_region_h3k26me3.jl is the Julia script used to process experimental data and plot the average methylation pattern on a gene body.
