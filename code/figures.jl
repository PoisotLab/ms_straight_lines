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


# Figure 2 - Links estimate
min_species = 3
max_species = 750

# Counterfactuals of lssl model
lssl_posterior = CSV.read(joinpath("data", "posterior_distributions", "lssl.csv"))
lssl_cf_links = lssl_posterior[r"counterfactual_links"]
lssl_cf_links = lssl_cf_links[:, min_species:max_species]

# Counterfactuals of const connectance model
const_posterior = CSV.read(joinpath("data", "posterior_distributions", "const_posterior.csv"))
const_cf_links = const_posterior[r"counterfactual_links"]
const_cf_links = const_cf_links[:, min_species:max_species]

# Counterfactuals of power law model
powerlaw_posterior = CSV.read(joinpath("data", "posterior_distributions", "powerlaw_posterior.csv"))
powerlaw_cf_links = powerlaw_posterior[r"counterfactual_links"]
powerlaw_cf_links = powerlaw_cf_links[:, min_species:max_species]

# Counterfactuals of beta binomial model
betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))
betab_cf_links = betab_posterior[r"counterfactual_links"]
betab_cf_links = betab_cf_links[:, min_species:max_species]


log_zeros(quant) = quant > 0 ? log(quant) : 0

cf_species = min_species:max_species

# Function to plot the quantiles of the counterfactuals links of each model
# we use log_zeros function to plot the y-axis in log (to account for log(0))
function plot_links_quantile(model)
    quant_015 = log_zeros.(quantile.(eachcol(model), 0.015))
    quant_055 = log_zeros.(quantile.(eachcol(model), 0.055))
    quant_165 = log_zeros.(quantile.(eachcol(model), 0.165))
    quant_835 = log_zeros.(quantile.(eachcol(model), 0.835))
    quant_945 = log_zeros.(quantile.(eachcol(model), 0.945))
    quant_985 = log_zeros.(quantile.(eachcol(model), 0.985))

    quant_500 = log_zeros.(quantile.(eachcol(model), 0.5))

    plot(cf_species, quant_985, fill = quant_015, color=:lightgreen, label="") # 97% PI
    plot!(cf_species, quant_945, fill = quant_055, color=:green, label="") # 89% PI
    plot!(cf_species, quant_835, fill = quant_165, color=:darkgreen, label="") # 69% PI
    plot!(cf_species, quant_500, linecolor=:black, linewidth=4, label="") # Median link number
    scatter!(d[:nodes], log.(d[:links]), label="") # Empirical links
    xaxis!(:log, "Species richness")
    yaxis!("Number of links (log)")
end

plot_links_quantile(lssl_cf_links)
plot_links_quantile(const_cf_links)
plot_links_quantile(powerlaw_cf_links)
plot_links_quantile(betab_cf_links)

savefig(joinpath("figures", "fig_02a_link_species_lssl"))







# First trial
bb_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))

density(bb_posterior[:mu])
plot!(Beta(3,7))
xaxis!((0.06, 0.13))

bb_cf_links = bb_posterior[:,r"counterfactual"] # to change lssl_posterior[r"counterfactual_links"]

bb_cf_dropone = bb_cf_links[:,2:750]

bb_05 = quantile.(eachcol(bb_cf_dropone), 0.05)
bb_95 = quantile.(eachcol(bb_cf_dropone), 0.95)

plot(2:750, bb_95, fill = bb_05, color=:lightgreen, label="") # 97% PI
scatter!(d[:nodes], d[:links], label="", color = :orange) # Empirical links
xaxis!(:log, "Species richness")
yaxis!(:log, "Number of links (log)")


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

# Beta distribution -- median parameters
beta_dist = Beta(p_median * phi_median, (1 - p_median) * phi_median)

## scale the "expected" distribution according to the mimum value:
bquant = LocationScale.(ms, 1 .- ms, beta_dist)

# Quantiles to plot
beta015 = quantile.(bquant, 0.015)
beta985 = quantile.(bquant, 0.985)
beta945 = quantile.(bquant, 0.945)

beta055 = quantile.(bquant, 0.055)
beta165 = quantile.(bquant, 0.165)
beta835 = quantile.(bquant, 0.835)

beta500 = quantile.(bquant, 0.5)


# Empirical connectance
links = d[:links]
species = d[:nodes]
co_emp = links ./ (species .^2)
S = 2:1:1000
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


# Beta distribution
beta_dist = Beta(p_median .* phi_median, (1 .- p_median) .* phi_median)

BetaBinomial(100, p_median * phi_median, (1 - p_median) * phi_median)

max = 24
S =2:1:max# 2:1:1000
BBin_S = BetaBinomial.((S.^2 .- S .+ 1), p_median * phi_median, (1 - p_median) * phi_median)

# need to shift this by S-1 to be the right value

# when the quantiles are calculated

plot(S,quantile.(BBin_S, 1-0.1) .+ S .- 1, fill = quantile.(BBin_S, 0.1) .+ S .- 1,color=:lightgreen, label="")
plot!(S,quantile.(BBin_S, 1-0.3) .+ S .- 1, fill = quantile.(BBin_S, 0.3) .+ S .- 1,color=:green, label="" )
plot!(S,quantile.(BBin_S, 1-0.4) .+ S .- 1, fill = quantile.(BBin_S, 0.4) .+ S .- 1,color=:darkgreen, label="")
ds = d[d.nodes .< max,:   ]
scatter!(ds[:nodes], ds[:links], lab="")
plot!(S,S.-1, lab="")
plot!(S,S.^2, lab ="")
xaxis!(:log, "Species richness", (2, 35))
yaxis!(:log, "Number of links", (1, 200))
