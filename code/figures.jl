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

bb_posterior = CSV.read(joinpath("data", "beta_binomial_posterior.csv"))
@df bb_posterior density(:a)
plot!(Beta(3,7))
xaxis!((0.06, 0.13))

## compare density of posterior to pencil trick
## add value of mean to text

# Figure 3
function draw_posterior(row)
    α = row[:a]*row[:theta]
    β = (1.0-row[:a])*row[:theta]
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


########## Some tests (to be changed) ##########
p_post = bb_posterior[:a]
mean_p = mean(p_post)
median_p = median(p_post)

density(bb_posterior[:a])
density(bb_posterior[:theta])
density(bb_posterior[:a] .* bb_posterior[:theta])
density((1 .- bb_posterior[:a]) .* bb_posterior[:theta])

number_of_trials = S.*S.-(S.-1)

L_predict = zeros(Int64, (length(number_of_trials), size(bb_posterior)[1]))

for i in 1:length(number_of_trials)
    for j in 1:size(bb_posterior)[1]
        p = bb_posterior[j, :a]
        theta = bb_posterior[j, :theta]
        n = number_of_trials[i]

        α = p*theta
        β = (1.0-p)*theta
        L_predict[i,j] = rand(BetaBinomial(n, α, β))
    end
end

L_predict
co = L_predict ./ (S.^2)

density(co[, :])
maximum(co)


########## End of tests ##########



# B - Extent to which the relationship gets closer to a power law (k)

k_predict = zeros(Float64, (length(S), size(bb_posterior)[1]))

for (i,s) in enumerate(S), (j,p) in enumerate(bb_posterior[:a])
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


# quantiles: 0.015, 0.055, 0.165, 0.5, 0.835, 0.945, 0.985

get_quantile(0.015)
