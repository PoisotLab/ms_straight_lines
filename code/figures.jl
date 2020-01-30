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


## Figure 1 - Beta fit

d.exr = d.links .- (d.nodes .- 1)
d.exp = d.nodes.^2 .- (d.nodes .- 1)
d.pex = d.exr ./ d.exp

p = fit(Beta, d.pex)

betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))

using Random
Random.seed!(1234)
index = rand(1:size(betab_posterior,2), 20)
mu_random = betab_posterior[index, :mu]
phi_random = exp.(betab_posterior[index, :phi])


betab_random = Beta.(mu_random .* phi_random, (1 .- mu_random) .* phi_random)

density(d.pex, c=:lightgrey, ls=:dash, fill=(:lightgrey, 0), frame=:semi, dpi=300, size=(400,400), lab="Empirical data")
density!(rand(p, 100_000), c=:black, ls=:dash, linewidth=2, lab="Fit")
plot!(betab_random[1], c=:darkgreen, linewidth=0.5, alpha=0.4, lab="Posterior")
for i in 1:length(index)
    plot!(betab_random[i], c=:darkgreen, linewidth=0.5, alpha=0.2, lab="")
end
xaxis!((0, 0.5), "p")
yaxis!((0, 9), "Density")
savefig(joinpath("figures", "penciltrick"))





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

log_zeros(quant) = quant > 0 ? log10(quant) : 0

cf_species = min_species:max_species
min_links_log = log10.(cf_species .- 1)
max_links_log = log10.(cf_species .^2)

# Function to plot the quantiles of the counterfactuals links of each model
# we use log_zeros function to plot the y-axis in log (to account for log(0))
function plot_links_quantile(model; title = "", xlabel = "", ylabel = "")
    quant_015 = log_zeros.(quantile.(eachcol(model), 0.015))
    quant_055 = log_zeros.(quantile.(eachcol(model), 0.055))
    quant_165 = log_zeros.(quantile.(eachcol(model), 0.165))
    quant_835 = log_zeros.(quantile.(eachcol(model), 0.835))
    quant_945 = log_zeros.(quantile.(eachcol(model), 0.945))
    quant_985 = log_zeros.(quantile.(eachcol(model), 0.985))

    quant_500 = log_zeros.(quantile.(eachcol(model), 0.5))

    plot(cf_species, quant_985, fill = quant_015, color="#C0C0C0", label="",
        title = title, xlabel = xlabel, ylabel = ylabel) # 97% PI
    plot!(cf_species, quant_945, fill = quant_055, color="#A0A0A0", label="") # 89% PI
    plot!(cf_species, quant_835, fill = quant_165, color="#808080", label="") # 69% PI
    plot!(cf_species, quant_500, linecolor=:black, linewidth=4, label="") # Median link number
    plot!(cf_species, min_links_log, linecolor=:black, label="") # Minimum number of links
    plot!(cf_species, max_links_log, linecolor=:black, label="") # Maximum number of links
    scatter!(d[:nodes], log10.(d[:links]), color=:yellow, alpha=0.6, markersize = 3, label="") # Empirical links
    xaxis!(:log, xlabel = xlabel)
    yaxis!(ylabel = ylabel)
end

plot_lssl = plot_links_quantile(lssl_cf_links, title = "A",
    ylabel = "Number of links")
plot_const = plot_links_quantile(const_cf_links, title = "B")
plot_powerlaw = plot_links_quantile(powerlaw_cf_links, title = "C",
    xlabel = "Species richness", ylabel = "Number of links")
plot_betab = plot_links_quantile(betab_cf_links, title = "D",
    xlabel = "Species richness")

plot(plot_lssl, plot_const, plot_powerlaw, plot_betab, layout = (2,2))
savefig(joinpath("figures", "fig_02_link_species_4models"))




# First trial
bb_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))

density(bb_posterior[:mu])
plot!(Beta(3,7))
xaxis!((0.06, 0.13))

bb_cf_dropone = bb_cf_links[:,2:750]

bb_05 = quantile.(eachcol(bb_cf_dropone), 0.05)
bb_95 = quantile.(eachcol(bb_cf_dropone), 0.95)

plot(2:750, bb_95, fill = bb_05, color=:lightgreen, label="") # 97% PI
scatter!(d[:nodes], d[:links], label="", color = :orange) # Empirical links
xaxis!(:log, "Species richness")
yaxis!(:log, "Number of links (log)")


## compare density of posterior to pencil trick
## add value of mean to text


# Figure 3 - BetaBinomial predictions

betab_mu_map = median(betab_posterior[:mu])
betab_phi_map = exp(median(betab_posterior[:phi]))

α = betab_mu_map*betab_phi_map
β = (1.0-betab_mu_map)*betab_phi_map

S = 3:750
n = S.^2 .- (S .- 1)

betab_draw = zeros(Int64, (nb_draw, length(S)))

# counterfactuals
# (i+2) because S starts at 3 species
for i in 1:nb_draw, (j,n) in enumerate(n)
    betab_draw[i,j] = rand(BetaBinomial(n, α, β)) + (j+2) - 1
end


plot_links_quantile(betab_draw)
savefig(joinpath("figures", "fig_03_link_species_betab"))



# Normal approximation of BetaBinomial

means = n .* betab_mu_map .+ S .- 1
vars = n .* betab_mu_map .* (1 .- betab_mu_map) .* (1 .+ S .* (S .- 1) .* (1 / (1 + betab_phi_map)))

approxs = Normal.(means, sqrt.(vars))

approx_89 = quantile.(approxs, 0.89)
approx_11 = quantile.(approxs, 0.11)

plot(S, approx_89, fill = approx_11,label = "", colour = :grey)
scatter!(d[:nodes], d[:links], label="", color = :orange) # Empirical links
xaxis!(:log, "Species richness")
yaxis!(:log, "Number of links (log)")
savefig(joinpath("figures", "fig_03b_link_species_betab_normal"))





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
S = 1:200
# Minimum number of species
ms = (S .- 1) ./ (S .^2)

# Median p and phi (for beta distribution)
p_median = median(bb_posterior[:mu])
phi_median = exp(median(bb_posterior[:phi]))

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


# L / S distribution of average degree

S = 1.5:.1:800

ms = (S .- 1) ./ (S)
## scale the "expected" distribution according to the mimum value:
bquant_LS = LocationScale.(ms, S .- ms, beta_dist)
beta015_LS = quantile.(bquant_LS, 0.015)
beta985_LS = quantile.(bquant_LS, 0.985)
beta11_LS = quantile.(bquant_LS, 0.11)
beta89_LS = quantile.(bquant_LS, 0.89)

plot(S, beta985_LS, fill = beta015_LS, lab = "")
plot!(S, beta11_LS, fill = beta89_LS, colour = :lightblue, lab = "")
xaxis!(:log, "Species richness")
yaxis!(:log, "Average degree, L/S")
## add real points
scatter!(d.nodes, d.links ./ d.nodes, colour = :grey, lab = "")


plot!(S, S, lab="")
plot!(S, (S .- 1)./S, lab="")

savefig(joinpath("figures", "L_S_distribution"))


#histogram(d.links ./ d.nodes, lab = "")

minimum(d.links ./ d.nodes)

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
