import Pkg; Pkg.activate(".")

using Distributions
using Statistics
using StatsBase


# z-score of the oceanic network

# map values of posterior distribution
betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))
mu_map = median(betab_posterior[:mu])
phi_map = exp(median(betab_posterior[:phi]))

# Empirical number of species and links
Sx = 11367
Lx = 7062647

means = (Sx^2 - Sx + 1) * mu_map + Sx - 1
sds = sqrt((Sx^2 - Sx + 1) * mu_map * (1 - mu_map) * (1 + Sx* (Sx - 1) / (1 + phi_map)))

z_score = (Lx - means) / sds
