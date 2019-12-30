import Pkg; Pkg.activate(".")
import CSV
using DataFrames
using StatsPlots


# include functions
include("common_functions.jl")
##
beta_bin_df = CSV.read(joinpath(pwd(), "data", "beta_binomial_posterior.csv"))
powerlaw_df = CSV.read(joinpath(pwd(), "data", "pwrlaw_posterior.csv"))
constant_df = CSV.read(joinpath(pwd(), "data", "const_posterior.csv"))


# get the data and filter for predation only
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]



plot_posterior(constant_df)
savefig(joinpath("figures", "constant_connectance.png"))


plot_posterior(powerlaw_df)
savefig(joinpath("figures", "powerlaw_connectance.png"))

bb_plot = plot_posterior_betabin(beta_bin_df)
xaxis!(bb_plot,:identity,"Species", xlims = (3,60))
yaxis!(bb_plot,:identity,"Links", ylims = (3,100))
savefig(joinpath("figures", "beta_binomial_connectance.png"))

## make one with all three
##plot(S_L, S_LD, S_CO, layout=(1,3), size=(800,300), dpi=300, leg=false)
#xaxis!("Species")
#savefig(joinpath("figures", "relationships.png"))

pp = constant_df[:, r"y_hat"]

pp_plot = scatter(d.nodes, [pp[15,i] for i in 1:ncol(pp)],
                  alpha = 0.01, color = "black", legend = false)

for k in rand(1:2000, 45)
    scatter!(pp_plot, d.nodes, [pp[k,i] for i in 1:ncol(pp)],
            alpha = 0.01, color = "black", legend = false)
end

xaxis!(pp_plot, "Species", xlims = (0,60))
yaxis!(pp_plot, "Links", ylims = (0,100))

scatter!(d.nodes, d.links, color = "darkorange")
