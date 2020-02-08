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

# posterior samples for the flexible links model
betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))

# posterior samples from previous models
lssl_posterior = CSV.read(joinpath("data", "posterior_distributions", "lssl.csv"))
const_posterior = CSV.read(joinpath("data", "posterior_distributions", "const_posterior.csv"))
powerlaw_posterior = CSV.read(joinpath("data", "posterior_distributions", "powerlaw_posterior.csv"))

# color palette for models
pal = (
    lssl=RGB(230/255,159/255,0/255),
    cc=RGB(86/255,190/255,233/255),
    pl=RGB(0/255,158/255,115/255),
    fl=RGB(213/255,94/255,0/255)
    )

sym = (
    lssl = :dashdot
)

# number of species
S = 3:750
mms = S.-1 # min links
Ms = S.^2 # max links
ms = mms ./ Ms  # min connectance
msl = (S .- 1) ./ (S) # min link density
total_flex = S .^ 2 .-S .+1

# map estimates
mu_map = median(betab_posterior[:mu])
phi_map = exp(median(betab_posterior[:phi]))

α = mu_map*phi_map
β = (1.0-mu_map)*phi_map

beta_map = Beta(α, β)
betabin_map = BetaBinomial.(total_flex, α, β)


# counterfactuals
betab_cf_links = betab_posterior[r"counterfactual_links"]
betab_cf_links = betab_cf_links[:, S]

lssl_cf_links = lssl_posterior[r"counterfactual_links"]
lssl_cf_links = lssl_cf_links[:, S]

const_cf_links = const_posterior[r"counterfactual_links"]
const_cf_links = const_cf_links[:, S]

powerlaw_cf_links = powerlaw_posterior[r"counterfactual_links"]
powerlaw_cf_links = powerlaw_cf_links[:, S]

# Figure 1 -- Beta fit with posterior samples

# generate posterior draws of the Beta distribution
Random.seed!(1234)
index = rand(1:size(betab_posterior,2), 20) # 20 posterior samples
mu_rdm = betab_posterior[index, :mu]
phi_rdm = exp.(betab_posterior[index, :phi])

betab_random = Beta.(mu_rdm .* phi_rdm, (1 .- mu_rdm) .* phi_rdm)

# MLE fit
pex = (d.links .- (d.nodes .- 1)) ./  (d.nodes.^2 .- (d.nodes .- 1))

p = fit(Beta, pex)

# calculate for text
phi_MLE = p.α + p.β
mu_MLE = p.α / phi_MLE

density(pex, c=:lightgrey, fill=(:lightgrey, 0, 0.5), dpi=120, size=(800,500), lab="Empirical data",
    foreground_color_legend=nothing, background_color_legend=:white, framestyle=:box)
plot!(betab_random[1], c=pal.fl, linewidth=1, alpha=0.3, lab="Posterior samples")
for i in 1:length(index)
    plot!(betab_random[i], c=pal.fl, linewidth=1, alpha=0.3, lab="")
end
density!(rand(p, 100_000), c=:black, ls=:dash, linewidth=2, lab="MLE fit")
xaxis!((0, 0.5), "p")
yaxis!((0, 9.5), "Density")
savefig(joinpath("figures", "beta_fit"))


# Figure 2 -- Links estimate from counterfactuals of the 4 models

neg_to_zeros(quant) = quant > 0 ? quant : 0.01 # to deal with negative quantiles

# Function to plot the quantiles of the counterfactuals links of each model
# we use log_zeros function to plot the y-axis in log (to account for log(0))
function plot_links_quantile(model; title="", xlabel="", ylabel="", linecolor="")
    quant_015 = neg_to_zeros.(quantile.(eachcol(model), 0.015))
    quant_110 = neg_to_zeros.(quantile.(eachcol(model), 0.11))
    quant_890 = neg_to_zeros.(quantile.(eachcol(model), 0.89))
    quant_985 = neg_to_zeros.(quantile.(eachcol(model), 0.985))
    quant_500 = neg_to_zeros.(quantile.(eachcol(model), 0.5))

    plot(S, quant_985, fill=quant_015, color=:grey, alpha=0.15, label="",
        title=title, title_location=:left, titlefontsize=11,
        xlabel=xlabel, ylabel=ylabel, framestyle=:box) # 97% PI
    plot!(S, quant_890, fill=quant_110, color=:grey, alpha=0.15, label="") # 89% PI
    scatter!(d[:nodes], d[:links], c=:grey, msw=0, markersize=5, label="") # Empirical links
    plot!(S, quant_500, linecolor=linecolor, linewidth=3, label="") # Median link number
    plot!(S, mms, linecolor=:black, lw=1, label="") # Minimum number of links
    plot!(S, Ms, linecolor=:black, lw=2, label="") # Maximum number of links
    xaxis!(:log, xlabel=xlabel, xlims=(minimum(S), maximum(S)))
    yaxis!(:log, ylims = (1,100000), ylabel=ylabel)
end

plot_lssl = plot_links_quantile(lssl_cf_links, title="A. LSSL",
    ylabel="Number of links", linecolor=pal.lssl)
plot_const = plot_links_quantile(const_cf_links, title="B. Constant connectance",
    linecolor=pal.cc)
plot_powerlaw = plot_links_quantile(powerlaw_cf_links, title="C. Power law",
    xlabel="Species richness", ylabel="Number of links", linecolor=pal.pl)
plot_betab = plot_links_quantile(betab_cf_links, title="D. Flexible links",
    xlabel="Species richness", linecolor=pal.fl)

plot(plot_lssl, plot_const, plot_powerlaw, plot_betab, layout=(2,2), size=(700,700), dpi=120)
savefig(joinpath("figures", "models_links"))


# Figure 3 -- Links estimate from MAP and normal approximation

# 3A BetaBinomial predictions from map values

bb_rand = rand.(betabin_map, 5000)

beta_89 = quantile.(bb_rand, 0.89) .+ S .- 1
beta_11 = quantile.(bb_rand, 0.11) .+ S .- 1
beta_98 = quantile.(bb_rand, 0.985) .+ S .- 1
beta_02 = quantile.(bb_rand, 0.015) .+ S .- 1
beta_50 = quantile.(bb_rand, 0.5)  .+ S .- 1

links_beta_map = plot(S, beta_98, fill=beta_02,label="", colour=:darkgrey,
    title="A. Flexible links (MAP)", title_location=:left, titlefontsize=11, framestyle=:box)
plot!(S, beta_89, fill=beta_11,label="", colour=:lightgrey)
plot!(S, beta_50, color=betab_color, label="", linewidth=4)
plot!(S, mms, linecolor=:black, label="") # Minimum number of links
plot!(S, Ms, linecolor=:black, label="") # Maximum number of links
scatter!(d[:nodes], d[:links], label="", color=:grey) # Empirical links
xaxis!(:log, "Species richness", label="")
yaxis!(:log, "Number of links")

# 3B Normal approximation of BetaBinomial
means = (Ms .- mms) .* mu_map .+ S .- 1
vars = (Ms .- mms) .* mu_map .* (1 .- mu_map) .* (1 .+ S .* (S .- 1) .* (1 / (1 + phi_map)))

approxs = Normal.(means, sqrt.(vars))
tnormal = truncated.(approxs, 0.01, Inf)

tnormal_89 = quantile.(tnormal, 0.89)
tnormal_11 = quantile.(tnormal, 0.11)
tnormal_98 = quantile.(tnormal, 0.985)
tnormal_02 = quantile.(tnormal, 0.015)
tnormal_50 = quantile.(tnormal, 0.5)

links_normal = plot(S, tnormal_98, fill=tnormal_02,label="", colour=:darkgrey,
    title="B. Normal approximation", title_location=:left, titlefontsize=11, framestyle=:box)
plot!(S, tnormal_89, fill=tnormal_11,label="", colour=:lightgrey)
plot!(S, tnormal_50, color=betab_color, label="", linewidth=4)
plot!(S, mms, linecolor=:black, label="") # Minimum number of links
plot!(S, Ms, linecolor=:black, label="") # Maximum number of links
scatter!(d[:nodes], d[:links], label="", color=:grey) # Empirical links
xaxis!(:log, "Species richness")
yaxis!(:log, "")

plot(links_beta_map, links_normal, layout=(1,2))
savefig(joinpath("figures", "betabinmap_normal_links"))




# Figure 4 - Connectance and link density

# 4A - connectance - species

## scale the "expected" distribution according to the minimum value:
bquant = LocationScale.(ms, 1 .- ms, beta_map)

# Quantiles to plot
beta015 = quantile.(bquant, 0.015)
beta985 = quantile.(bquant, 0.985)
beta11 = quantile.(bquant, 0.110)
beta89 = quantile.(bquant, 0.890)

beta500 = quantile.(bquant, 0.5)

# Empirical connectance
co_emp = d[:links] ./ (d[:nodes] .^2)


# Connectance vs species
connectance_beta = plot(S, range(beta015, stop=beta985, length=1000), color=:darkgrey,  fill=:darkgrey,
    label="", title="A", title_location=:left, titlefontsize=11, framestyle=:box) # 97% PI
plot!(S, range(beta11, stop=beta89, length=1000), color=:lightgrey, fill=:lightgrey, label="") # 89% PI
plot!(S, beta500, linecolor=betab_color, linewidth=4, label="") # Median connectance
scatter!(d[:nodes], co_emp, color=:grey, label="") # Empirical connectance
plot!(S, ms, label="", linecolor=:black) # Minimum connectance
xaxis!(:log, "Species richness")
yaxis!("Connectance")


# 4B distribution of linkage density

## scale the "expected" distribution according to the mimum value:
bquant_LS = LocationScale.(msl, S .- msl, beta_map)

beta015_LS = quantile.(bquant_LS, 0.015)
beta985_LS = quantile.(bquant_LS, 0.985)
beta11_LS = quantile.(bquant_LS, 0.11)
beta89_LS = quantile.(bquant_LS, 0.89)
beta50_LS = quantile.(bquant_LS, 0.50)

avg_degree_beta = plot(S, beta985_LS, fill=beta015_LS, color=:darkgrey, lab="", title="B",
    title_location=:left, titlefontsize=11, framestyle=:box)
plot!(S, beta11_LS, fill=beta89_LS, colour=:lightgrey, lab="",)
plot!(S, beta50_LS, linecolor=betab_color, linewidth=4, lab="")
scatter!(d.nodes, d.links ./ d.nodes, colour=:grey, lab="")
plot!(S, msl, linecolor=:black, lab="")
plot!(S, S, linecolor=:black, lab="")
xaxis!(:log, "Species richness")
yaxis!(:log, "Linkage density")


plot(connectance_beta, avg_degree_beta, label=(1,2), lab="")
savefig(joinpath("figures", "connectance_linkdens"))


# Figure 5 - Extent to which the relationship gets closer to a power law (k)

# Values of k function of s and posterior p

k_predict = zeros(Float64, (length(S), 10000))

for (i,s) in enumerate(S), j in 1:size(k_predict, 2)
    p = rand(beta_map)
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

plot(S, get_quantile(0.015), fill=get_quantile(0.985), color=:darkgrey, label="", framestyle=:box) # 97% PI
plot!(S, get_quantile(0.11), fill=get_quantile(0.89), color=:lightgrey, label="") # 89% PI
plot!(S, get_quantile(0.5), linecolor=betab_color, linewidth=4, lab = "")
xaxis!(:log, "Species richness")
yaxis!("k")
savefig(joinpath("figures", "k_powerlaw"))



# Figure 6 -- % of model prediction above or below minimum

function realistic_links(model_cf)
    realistic = zeros(Float64, (1, length(S)))
    for (i,s) in enumerate(S)
        belowmin = length(findall(model_cf[:,i] .< (s - 1)))
        abovemax = length(findall(model_cf[:,i] .> (s^2)))
        realistic[i] = 1 - (belowmin + abovemax) / size(model_cf, 1)
    end
    return(vec(realistic))
end

realistic_betab =  realistic_links(betab_cf_links)
realistic_lssl = realistic_links(lssl_cf_links)
realistic_const = realistic_links(const_cf_links)
realistic_powerlaw = realistic_links(powerlaw_cf_links)


medianspecies = median(d[:nodes], 0.5)
species05 = quantile(d[:nodes], 0.05)
species95 = quantile(d[:nodes], 0.95)

plot([medianspecies], seriestype=:vline, color=:black, lab="")
plot!([species05], seriestype=:vline, color=:grey, lab="")
plot!([species95], seriestype=:vline, color=:grey, lab="")
plot!(S, realistic_lssl, color=lssl_color, linewidth=3, label="LSSL",
    legend=:bottomright, foreground_color_legend=nothing, background_color_legend=:white)
plot!(S, realistic_const, color=const_color, linewidth=2.5, label="Constant connectance")
plot!(S, realistic_powerlaw, color=powerlaw_color, linewidth=2.5, label="Power law")
plot!(S, realistic_betab, color=betab_color, linewidth=2.5, label="Shifted beta-binomial")
xaxis!(:log, "Species richness")
yaxis!("Proportion of species or realistic link numbers", (0.3, 1.05))
savefig(joinpath("figures", "real_predict"))
