import Pkg; Pkg.activate(".")
import CSV
using DataFrames
using StatsPlots
##
beta_bin_df = CSV.read(joinpath(pwd(), "data", "beta_binomial_posterior.csv"))
powerlaw_df = CSV.read(joinpath(pwd(), "data", "pwrlaw_posterior.csv"))
constant_df = CSV.read(joinpath(pwd(), "data", "const_posterior.csv"))


# get the data and filter for predation only
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]


pp = beta_bin_df[:, r"y_hat"]

pp_plot = scatter(d.nodes, [pp[15,i] for i in 1:ncol(pp)] .+ (d.nodes .- 1), alpha = 0.2)
xaxis!(pp_plot, :log)
yaxis!(pp_plot, :log, "Links")
for k in rand(1:2000, 4)
    scatter!(d.nodes, [pp[k,i] for i in 1:ncol(pp)] .+ (d.nodes .- 1), alpha = 0.2)
end
