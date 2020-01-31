import Pkg; Pkg.activate(".")

using Plots
using StatsPlots
import CSV
using Distributions
using Statistics
using StatsBase
using Random



# get the data and filter for predation only
d = CSV.read(joinpath("data", "network_data.dat"))
d = d[d.predation .> 0 , :]

# posterior samples for beta binomial model
betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))

# posterior samples from previous models
lssl_posterior = CSV.read(joinpath("data", "posterior_distributions", "lssl.csv"))
const_posterior = CSV.read(joinpath("data", "posterior_distributions", "const_posterior.csv"))
powerlaw_posterior = CSV.read(joinpath("data", "posterior_distributions", "powerlaw_posterior.csv"))

# map estimates
mu_map = median(betab_posterior[:mu])
phi_map = exp(median(betab_posterior[:phi]))

α = mu_map*phi_map
β = (1.0-mu_map)*phi_map

# number of species
S = 3:750
ms = S.-1 # min links
Ms = S.^2 # max links

# counterfactuals
betab_cf_links = betab_posterior[r"counterfactual_links"]
betab_cf_links = betab_cf_links[:, S]

lssl_cf_links = lssl_posterior[r"counterfactual_links"]
lssl_cf_links = lssl_cf_links[:, S]

const_cf_links = const_posterior[r"counterfactual_links"]
const_cf_links = const_cf_links[:, S]

powerlaw_cf_links = powerlaw_posterior[r"counterfactual_links"]
powerlaw_cf_links = powerlaw_cf_links[:, S]

# MLE of p
include("02_fit_p.jl")



# Figure 1 -- Beta fit with posterior samples

Random.seed!(1234)
index = rand(1:size(betab_posterior,2), 20) # 20 posterior samples
mu_rdm = betab_posterior[index, :mu]
phi_rdm = exp.(betab_posterior[index, :phi])

betab_random = Beta.(mu_rdm .* phi_rdm, (1 .- mu_rdm) .* phi_rdm)

density(d.pex, c=:lightgrey, fill=(:lightgrey, 0), frame=:semi, dpi=300, size=(400,400), lab="Empirical data")
density!(rand(p, 100_000), c=:black, ls=:dash, linewidth=2, lab="MLE fit")
plot!(betab_random[1], c=:darkgreen, linewidth=1, alpha=0.4, lab="Posterior samples")
for i in 1:length(index)
    plot!(betab_random[i], c=:darkgreen, linewidth=1, alpha=0.2, lab="")
end
xaxis!((0, 0.5), "p")
yaxis!((0, 9), "Density")
savefig(joinpath("figures", "beta_fit"))


# Figure 2 -- Links estimate from counterfactuals of the 4 models

log_zeros(quant) = quant > 0 ? log10(quant) : 0 # to deal with negative quantiles

# Function to plot the quantiles of the counterfactuals links of each model
# we use log_zeros function to plot the y-axis in log (to account for log(0))
function plot_links_quantile(model; title = "", xlabel = "", ylabel = "")
    quant_015 = log_zeros.(quantile.(eachcol(model), 0.015))
    quant_110 = log_zeros.(quantile.(eachcol(model), 0.11))
    quant_890 = log_zeros.(quantile.(eachcol(model), 0.89))
    quant_985 = log_zeros.(quantile.(eachcol(model), 0.985))

    quant_500 = log_zeros.(quantile.(eachcol(model), 0.5))

    plot(cf_species, quant_985, fill=quant_015, color=:"#03a1fc", label="",
        title=title, xlabel=xlabel, ylabel=ylabel) # 97% PI
    plot!(cf_species, quant_890, fill=quant_110, color=:lightblue, label="") # 89% PI
    plot!(cf_species, quant_500, linecolor=:black, linewidth=2, label="") # Median link number
    plot!(cf_species, log10.(ms), linecolor=:black, label="") # Minimum number of links
    plot!(cf_species, log10.(Ms), linecolor=:black, label="") # Maximum number of links
    scatter!(d[:nodes], log10.(d[:links]), color=:grey, markersize=3, label="") # Empirical links
    xaxis!(:log, xlabel=xlabel)
    yaxis!(ylabel=ylabel)
end

plot_lssl = plot_links_quantile(lssl_cf_links, title="A - lssl",
    ylabel="Number of links")
plot_const = plot_links_quantile(const_cf_links, title="B - constant connect")
plot_powerlaw = plot_links_quantile(powerlaw_cf_links, title="C - power law",
    xlabel="Species richness", ylabel="Number of links")
plot_betab = plot_links_quantile(betab_cf_links, title="D - beta binomial",
    xlabel="Species richness")

plot(plot_lssl, plot_const, plot_powerlaw, plot_betab, layout=(2,2))
savefig(joinpath("figures", "models_estimate_links"))



# Figure 3 -- Links estimate from MAP and normal approximation

# 3A BetaBinomial predictions from map values

beta_map = Beta(α, β)

beta_89 = quantile.(beta, 0.89) .* (Ms .- ms) .+ S .- 1
beta_11 = quantile.(beta, 0.11) .* (Ms .- ms) .+ S .- 1
beta_50 = quantile.(beta, 0.5) .* (Ms .- ms) .+ S .- 1

links_betab = plot(S, beta_89, fill=beta_11,label="", colour=:grey)
plot!(S, beta_50, color=:black, label="", linewidth=2)
scatter!(d[:nodes], d[:links], label="", color=:orange) # Empirical links
xaxis!(:log, "Species richness", label="")
yaxis!(:log, "Number of links")

# 3B Normal approximation of BetaBinomial
means = (Ms .- ms) .* mu_map .+ S .- 1
vars = (Ms .- ms) .* mu_map .* (1 .- mu_map) .* (1 .+ S .* (S .- 1) .* (1 / (1 + phi_map)))

approxs = Normal.(means, sqrt.(vars))

approx_89 = quantile.(approxs, 0.89)
approx_11 = quantile.(approxs, 0.11)

plot(S, approx_89, fill = approx_11,label = "", colour = :grey)
scatter!(d[:nodes], d[:links], label="", color = :orange) # Empirical links
xaxis!(:log, "Species richness")
yaxis!(:log, "Number of links (log)")
savefig(joinpath("figures", "ink_species_betab_normal"))





# Tim's question
Sx = 11367
Lx = 7062647

means = (Sx^2 - Sx + 1) * betab_mu_map + Sx - 1
sds = sqrt((Sx^2 - Sx + 1) * betab_mu_map * (1 - betab_mu_map) * (1 + Sx* (Sx - 1) / (1 + betab_phi_map)))

(Lx - means) / sds








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
S = 3:750
# Minimum number of species
ms = (S .- 1) ./ (S .^2)

# Median p and phi (for beta distribution)
p_median = median(betab_posterior[:mu])
phi_median = exp(median(betab_posterior[:phi]))

# Beta distribution -- median parameters
beta_dist = Beta(p_median * phi_median, (1 - p_median) * phi_median)

## scale the "expected" distribution according to the mimum value:
bquant = LocationScale.(ms, 1 .- ms, beta_dist)

# Quantiles to plot
beta015 = quantile.(bquant, 0.015)
beta985 = quantile.(bquant, 0.985)
beta11 = quantile.(bquant, 0.110)
beta89 = quantile.(bquant, 0.890)

beta500 = quantile.(bquant, 0.5)


# Empirical connectance
links = d[:links]
species = d[:nodes]
co_emp = links ./ (species .^2)

# Connectance vs species
plotA = plot(S, range(beta015, stop=beta985, length=1000), color=:"#03a1fc",  fill=:"#03a1fc", label="") # 97% PI
plot!(S, range(beta11, stop=beta89, length=1000), color=:lightblue, fill=:lightblue, label="") # 89% PI
plot!(S, beta500, linecolor=:black, linewidth=2, label="") # Median connectance
scatter!(species, co_emp, color=:grey, label="") # Empirical connectance
plot!(S, ms, label="", linecolor=:black) # Minimum connectance
xaxis!(:log, "Species richness")
yaxis!("Connectance")
savefig(joinpath("figures", "fig_04a_connectance_species"))


# 4B L/ S distribution of average degree

ms = (S .- 1) ./ (S)

## scale the "expected" distribution according to the mimum value:
bquant_LS = LocationScale.(ms, S .- ms, beta_dist)
beta015_LS = quantile.(bquant_LS, 0.015)
beta985_LS = quantile.(bquant_LS, 0.985)
beta11_LS = quantile.(bquant_LS, 0.11)
beta89_LS = quantile.(bquant_LS, 0.89)
beta50_LS = quantile.(bquant_LS, 0.50)


plotB = plot(S, beta985_LS, fill = beta015_LS, color="#03a1fc", lab = "")
plot!(S, beta11_LS, fill = beta89_LS, colour = :lightblue, lab = "")
plot!(S, beta50_LS, linecolor=:black, linewidth=2, lab = "")
xaxis!(:log, "Species richness")
yaxis!(:log, "Average degree, L/S")
## add real points
scatter!(d.nodes, d.links ./ d.nodes, colour = :grey, lab = "")

plot!(S, (S .- 1)./S, linecolor=:black, lab="")
plot!(S, S, linecolor=:black, lab="")
savefig(joinpath("figures", "L_S_distribution"))


#histogram(d.links ./ d.nodes, lab = "")

minimum(d.links ./ d.nodes)


## 4C - Extent to which the relationship gets closer to a power law (k)

# Values of k function of s and posterior p

mu = median(betab_posterior[:mu])
phi =median(exp.(betab_posterior[:phi]))
beta_p = Beta(mu * phi, (1 - mu) * phi)

k_predict = zeros(Float64, (length(S), 10000))

for (i,s) in enumerate(S), j in 1:size(k_predict, 2)
    p = rand(beta_p)
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


# k - species plot with quantiles:
# 67% percentile interval: quantiles 0.165 and 0.835
# 89% percentile interval: quantiles 0.055 and 0.945
# 97% percentile interval: quantiles 0.015 and 0.985

plotC = plot(S, get_quantile(0.015), fill=get_quantile(0.985), color=:"#03a1fc", label="") # 97% PI
plot!(S, get_quantile(0.11), fill=get_quantile(0.89), color=:lightblue, label="") # 89% PI
plot!(S, get_quantile(0.5), linecolor=:black, linewidth=2, lab = "")
xaxis!(:log, "Species richness")
yaxis!("k")
savefig(joinpath("figures", "fig_04b_k_species"))

median_k = get_quantile(0.5)
findall(x -> x < 0.1, median_k)


plot(plotA, plotB, plotC, layout = (1,3))
savefig(joinpath("figures", "fig_04_linkdens_connect_k"))





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
