import Pkg; Pkg.activate(".")
import CSV
using DataFrames
using StatsBase
using Statistics

##
# posterior samples for beta binomial model
betab_posterior = CSV.read(joinpath("data", "posterior_distributions", "beta_binomial_posterior.csv"))
# posterior samples from previous models
lssl_posterior = CSV.read(joinpath("data", "posterior_distributions", "lssl.csv"))
const_posterior = CSV.read(joinpath("data", "posterior_distributions", "const_posterior.csv"))
powerlaw_posterior = CSV.read(joinpath("data", "posterior_distributions", "powerlaw_posterior.csv"))

## BetaBinomial
betab_posterior[:mu] |> mean
betab_posterior[:mu] |> std

betab_posterior[:phi] .|> exp |> mean
betab_posterior[:mu] .|> exp |> std

## Power Law
powerlaw_posterior[:a] .|> exp |> mean
powerlaw_posterior[:a] .|> exp |> std

powerlaw_posterior[:b] |> mean
powerlaw_posterior[:b] |> std

## Constant Connectance
const_posterior[:a] |> mean
const_posterior[:a] |> std

## LSSL
lssl_posterior[:a] .|> exp |> mean
lssl_posterior[:a] .|> exp |> std
