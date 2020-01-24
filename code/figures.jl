import Pkg; Pkg.activate(".")

using Plots
using StatsPlots
import CSV
using Distributions
using Statistics
using StatsBase

# get the data and filter for predation only
d = CSV.read(joinpath("data", "network_data.dat"))
d = d[d.predation .> 0 , :]

# Figure 1 -- Links - species
S = 2:1:1000
m(s) = s-1
M(s) = s*s
plot(S, m.(S), fill=(M.(S), :grey, 0.2), c=:transparent, lab="")
@df d scatter!(:nodes, :links, c=:black,
 ms = 2, alpha = 0.6,dpi = 300, lab = "",
 frame = :box)
xaxis!(:log, "Species richness", (2, 1000))
yaxis!(:log, "Number of links", (1, 1000000))
savefig(joinpath("figures", "fig_01_link_species"))


# Figure 2

bb_posterior = CSV.read(joinpath("data", "posterior_distributions", "bb_posterior.csv"))
@df bb_posterior density(:p)
plot!(Beta(3,7))
xaxis!((0.06, 0.13))

## compare density of posterior to pencil trick
## add value of mean to text

# Figure 3
function draw_posterior(row)
    α = row[:p]*row[:theta]
    β = (1.0-row[:p])*row[:theta]
    return (n) -> BetaBinomial(n, α, β)
end

distributions = draw_posterior.(eachrow(bb_posterior))
number_of_trials = S.*S.-(S.-1)

L_predict = zeros(Int64, (length(number_of_trials), length(d)))
for (i, t) in enumerate(number_of_trials), (j, d) in enumerate(distributions)
    L_predict[j,i] = rand(d(t))
end

L(s,p) = (s*s-(s-1))*p+(s-1)

L_predict = L.(S', random_p)
quantile_functions = vec(mapslices(ecdf, L_predict, dims=2))
L_possible = 1:100:1_000_000

Q99 = zeros(Float64, length(quantile_functions))

for (i, q) in enumerate(quantile_functions)
    temp_q = q.(L_possible)
    diff_q = temp_q .- 0.50
    Q99[i] = L_possible[findfirst(x -> x>0.0, diff_q)]
end

scatter(S, Q99)
yaxis!(:log)
xaxis!(:log)



# Figure 4

# 4A - connectance - species

# Minimum number of species
ms = (S .- 1) ./ (S .^2)

# Median p and phi (for beta distribution)
p_median = median(bb_posterior[:p])
phi_median = exp(median(bb_posterior[:theta]))

# Beta distribution
beta_dist = Beta.(p_median .* phi_median, (1 .- p_median) .* phi_median)

# Quantiles to plot
beta015 = quantile(beta_dist, 0.015) .* (1 .- ms) .+ ms
beta985 = quantile(beta_dist, 0.985) .* (1 .- ms) .+ ms
beta055 = quantile(beta_dist, 0.055) .* (1 .- ms) .+ ms
beta945 = quantile(beta_dist, 0.945) .* (1 .- ms) .+ ms
beta165 = quantile(beta_dist, 0.165) .* (1 .- ms) .+ ms
beta835 = quantile(beta_dist, 0.835) .* (1 .- ms) .+ ms
beta500 = quantile(beta_dist, 0.5) .* (1 .- ms) .+ ms


# Empirical connectance
links = d[:links]
species = d[:nodes]
co_emp = links ./ (species .^2)

# Connectance vs species
plot(S, range(beta015, stop=beta985, length=300), color=:lightgreen, fill=:lightgreen, label="") # 97% PI
plot!(S, range(beta055, stop=beta945, length=300), color=:green, fill=:green, label="") # 89% PI
plot!(S, range(beta165, stop=beta835, length=300), color=:darkgreen, fill=:darkgreen, label="") # 67% PI
plot!(S, beta500, linecolor=:black, linewidth=4, label="") # Median connectance
scatter!(species, co_emp, label="") # Empirical connectance
plot!(S, ms, linecolor=:black, label="") # Minimum connectance
xaxis!(:log, "Species richness")
yaxis!("Connectance")
savefig(joinpath("figures", "fig_04a_connectance_species"))


## 4B - Extent to which the relationship gets closer to a power law (k)

# Values of k function of s and posterior p
k_predict = zeros(Float64, (length(S), size(bb_posterior)[1]))

for (i,s) in enumerate(S), (j,p) in enumerate(bb_posterior[:p])
    k_predict[i, j] = ((1 - p) * s + (p - 1)) / (p * s^2)
end

# Replace the values of k by their row-wise quantiles
k_quantiles = zeros(Float64, size(k_predict))

n = 1000 # Round quantiles in n classes
for s in 1:length(S)
    qfinder = ecdf(k_predict[s,:]) # Create ECDF function

    k_quantiles[s,:] = qfinder(k_predict[s,:]) # Replace scores by quantiles

    k_quantiles[s,:] = round.(k_quantiles[s,:] * n) / n # Round quantiles
end

# Function for the mean values of k at the quantile q for each value of s
function get_quantile(q)
    k_quant = zeros(Float64, (length(S), 1))
    for s in 1:length(S)
        k_quant[s] = mean(k_predict[s, k_quantiles[s,:] .== q])
    end
    return(k_quant)
end

# Empirical k (to be removed)
links = d[:links]
species = d[:nodes]
p_emp = (links .- (species .- 1 )) ./ (species .^2 .- (species .- 1))
k_emp = ((1 .- p_emp) .* species .+ (p_emp .- 1)) ./ (p_emp .* species .^2)

# k - species plot with quantiles:
# 67% percentile interval: quantiles 0.165 and 0.835
# 89% percentile interval: quantiles 0.055 and 0.945
# 97% percentile interval: quantiles 0.015 and 0.985

plot(S, range(get_quantile(0.015), stop=get_quantile(0.985), length=300), color=:lightgreen, fill=:lightgreen, label="") # 97% PI
plot!(S, range(get_quantile(0.055), stop=get_quantile(0.945), length=300), color=:green, fill=:green, label="") # 89% PI
plot!(S, range(get_quantile(0.165), stop=get_quantile(0.835), length=300), color=:darkgreen, fill=:darkgreen, label="") # 67% PI
plot!(S, mean(k_predict, dims = 2), linecolor = :black, lab = "Mean")
xaxis!(:log, "Species richness")
yaxis!("k")
savefig(joinpath("figures", "fig_04b_k_species"))







###### Tim's problem ######
Sx = 11367
Lx = 7062647
p = median(bb_posterior[:p])
z = ((Lx - (Sx - 1)) - p * (Sx ^2 - (Sx - 1))) / (sqrt(p * (1 - p) * (Sx^2 - (Sx -1))))
Lx - (Sx - 1)) / (Sx ^2 - (Sx - 1)
minimum(bb_posterior[:p])
###### End of Tim's problem ######
