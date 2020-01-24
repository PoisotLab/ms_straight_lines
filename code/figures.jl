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

# A - connectance - species

# First try
bb_posterior
links_predict = zeros(Int64, (length(S), size(bb_posterior)[1]))

for (i,s) in enumerate(S), j in 1:size(bb_posterior)[1]
    n = s^2-(s-1)
    α = bb_posterior[j,:p]*bb_posterior[j,:theta]
    β = (1.0-bb_posterior[j,:p])*bb_posterior[j,:theta]
    links_predict[i,j] = rand(BetaBinomial(n, α, β))
end

co_predict = links_predict ./ (S.^2)

plot(S, mean(co_predict, dims = 2), linecolor = :black, lab = "Mean")


# Second try

# Parameter p
pco = zeros(Float64, (length(S), size(bb_posterior)[1]))
ms = (S .- 1) ./ (S .^2)

for (i,m) in enumerate(ms), j in 1:size(bb_posterior)[1]
    p_post = bb_posterior[j,:p]
    pco[i,j] = (1 - m) * p_post + m
end

# Parameter theta (phi)
phico = zeros(Float64, (length(S), size(bb_posterior)[1]))
for (i,m) in enumerate(ms), j in 1:size(bb_posterior)[1]
    phi_post = bb_posterior[j,:theta]
    phico[i,j] = (phi_post + m) / (1 - m)
end

# Regularized value of connectance
links = d[:links]
species = d[:nodes]
p_mean = mean(bb_posterior[:p])
phi_mean = mean(bb_posterior[:theta])

co_reg = (links .+ phi_mean * p_mean) ./ (species .^2 .+ phi_mean)

# Empirical connectance
co_emp = links ./ (species .^2)

# Connectance vs species
plot(S, mean(pco, dims=2), linecolor=:black, linewidth=4, label="")
plot!(S, ms, linecolor=:black, label="") # minimum value of connectance
scatter!(species, co_emp, label="")
xaxis!(:log, "Species richness")
yaxis!("Connectance")
savefig(joinpath("figures", "fig_04a_connectance_species"))





## B - Extent to which the relationship gets closer to a power law (k)

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
