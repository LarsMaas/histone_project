using Pkg
using CSV
using DataFrames
using Clustering

# This file contains enriched regions for different cell types, determined using an HMM
regions = CSV.File("data_h3k36me3/h3k36me3_regions.csv", header=3, skipto=4) |> DataFrame
# This file contains densities (how enriched H3K36me3 is) for ES cells, with a resolution of 25 bp
densities_esc = CSV.File("data_h3k36me3/GSM307620_ES.H3K36me3.densities.txt", header=["chromosome", "position", "density"], delim="\t") |> DataFrame
first(densities_esc, 5)

# Get the regions for the H3K36me3 modifications in ES cells
k36_regions = regions[:,5:6]
k36_regions = dropmissing(k36_regions, 2)

# Split the column containing the start and end of regions into seperate columns
k36_regions_esc = k36_regions[:,1] |> skipmissing |> collect
k36_regions_esc = split.(k36_regions_esc, ":")
k36_regions_esc = map(x -> [x[Not(2)]; parse.(Int, split(x[2], "-"))], k36_regions_esc)

# Split the densities by chromosome so we only have to do that once, this improves filtering efficiency a lot in the next steps
split_densities_esc = Dict(chrom => filter(:chromosome => chr -> chr == chrom, densities_esc) for chrom in unique(densities_esc.chromosome))

# For every region, return the densities
region_densities = [filter(:position => pos -> pos >= region[2] && pos <= region[3], split_densities_esc[region[1]]).density for region in k36_regions_esc]

# filter(:position => pos -> pos > 197069800 && pos < 197069900, split_densities_esc["chr1"]).density

# Number of regions
length(region_densities)
# Smallest region
length.(region_densities) |> minimum
# Largest region
length.(region_densities) |> maximum
# Histogram of the length of the regions
histogram(length.(region_densities))


####################################################
# GET SAME LENGTH AND AVERAGE REGION
# Make sure every region has an equal number of data points.

# create a list with x equally spaced numbers from 0 to 1
divide_by = 0:1/41:1 |> collect

# For every region get x equally spaced indices 
equal_len_region_densities = []
for region in region_densities
    region_indices = divide_by * length(region)
    region_indices[1] = 1

    # Get the densities for the x indices of the region. Get an average if the number has a decimal.
    equal_len_region = []
    for i in region_indices
        if isinteger(i)
            push!(equal_len_region, region[Int(i)])
        else
            push!(equal_len_region, region[floor(Int64, i)] + region[ceil(Int64, i)] / 2)
        end
    end
    push!(equal_len_region_densities, equal_len_region)
end

# Plot the average region density (y-axis) by relative position (x-axis)
plot(mean(equal_len_region_densities, dims=1), ylim=(0,6))

equal_len_region_densities

####################################################
# GET -150 AND +150 FROM THE MIDDLE
# Look at the density 200 positions up and down from the center of each region
region_width = 300

# for every region, if region > 300 do region[middle-150 : middle+150]
filter_reg_dens = [region[floor(Int64, length(region)/2-region_width/2):floor(Int64, length(region)/2+region_width/2)] for region in region_densities if length(region) > region_width+2]

length(filter_reg_dens)
avg_filter_reg_dens = mean(filter_reg_dens, dims=1)
plot(avg_filter_reg_dens, ylim=(0,5))

####################################################
# Look for k means clustering of different lenghts online (What they do is that they take a peak, take a region around that (for example - and + 200) and then cluster)
kmeans(avg_filter_reg_dens, 3)

# Look at different genes
# set interpelation to max and lowish
# Perform k means clustering
# Plot means
# 

####################################################
# Get regions per gene

gene_inf = CSV.File("data_h3k36me3/h3k36me3_genes.csv", header=3, skipto=4) |> DataFrame

gene_inf = filter(:("K36 DHMS") => dhms_region -> dhms_region == "ES", gene_inf)
first(gene_inf, 5)

gene_pos = Matrix(gene_inf[:,["txStart", "txEnd", "chrom"]])
gene_pos = [gene_pos[i, :] for i in 1:size(gene_pos,1)]
gene_pos[1]
####
# Extract genes + or - 100 bp (4 datapoints, since density is by 25bp)
gene_region_densities = [filter(:position => pos -> pos >= gene[1]-100 && pos <= gene[2]+100, split_densities_esc[gene[3]]).density for gene in gene_pos]
first(gene_region_densities, 10)

gene_region_densities
# Number of enriched genes
length(gene_region_densities)
# Smallest gene
length.(gene_region_densities) |> minimum
# Largest gene
length.(gene_region_densities) |> maximum
# Gene length distribution
length.(gene_region_densities) |> histogram
histogram(length.(gene_region_densities)   )#   , xlim=(0,3000), bins=0:50:3000)

plot(gene_region_densities[1])
plot(gene_region_densities[2])
plot(gene_region_densities[3])
plot(gene_region_densities[4])
plot(gene_region_densities[5])

####

plot(gene_region_densities[6], ylim=(0,4))

divide_by = 0:1/minimum(length.(gene_region_densities)):1 |> collect
divide_by = 0:1/40:1 |> collect

equal_len_gene_region_densities = []
for region in [reg for reg in gene_region_densities if length(reg) > 16]
    region_indices = divide_by * length(region)
    region_indices[1] = 1

    equal_len_region = []
    for i in region_indices
        if floor(i) == 0.0
            i = 1
        end
        if isinteger(i)
            push!(equal_len_region, region[Int(i)])
        else
            push!(equal_len_region, region[floor(Int64, i)] + region[ceil(Int64, i)] / 2)
        end
    end
    push!(equal_len_gene_region_densities, equal_len_region)
end

length.(equal_len_gene_region_densities) |> maximum

avg_density_genes = mean(equal_len_gene_region_densities, dims=1)
plot(avg_density_genes   )#   , ylim=(0,6))
plot([mean(avg_density_genes[1][i-5:i]) for i in 6:length(avg_density_genes[1])], label=false, xlabel="Relative histone position", ylabel="Methylation intensity", title="Methylation intensity ChIP-Seq")

avg_density_genes
length(equal_len_gene_region_densities)


####################################################

plot(region_densities[4])

region_densities[4]

# Moving average
mean1 = [mean(region_densities[4][i-10:i]) for i in 11:length(region_densities[4])]
plot(mean1)


###############################
# K means clustering using genes

Matrix(equal_len_gene_region_densities)

kmeans(equal_len_gene_region_densities, 3)