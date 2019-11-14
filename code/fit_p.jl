import Pkg
Pkg.activate(".")

import CSV
using Distributions
using DataFrames
using StatsPlots
using Statistics
using StatsBase

# Get the data
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]

# Add some informations
d.p = d.links ./ (d.nodes.^2 .- (d.nodes .- 1))
d.ld = d.links ./ d.nodes
d.co = d.links ./ (d.nodes .* d.nodes)

p = fit(LogitNormal, d.p)

# General figures to show the relationships
S_L = @df d scatter(:nodes, :links, frame=:semi, c=:black)
xaxis!(S_L, :log)
yaxis!(S_L, :log, "Links")

S_LD = @df d scatter(:nodes, :ld, frame=:semi, c=:black)
xaxis!(S_LD, :log)
yaxis!(S_LD, :log, "Average degree")

S_CO = @df d scatter(:nodes, :co, frame=:semi, c=:black)
xaxis!(S_CO, :log)
yaxis!(S_CO, (0,1), "Connectance")

plot(S_L, S_LD, S_CO, layout=(1,3), size=(800,300), dpi=300, leg=false)
xaxis!("Species")
savefig(joinpath("figures", "relationships.png"))

# Fit the data
density(rand(p, 100_000), c=:lightgrey, fill=(:lightgrey, 0), frame=:semi, lab="Fit", dpi=300, size=(400,400))
density!(d.p, c=:black, ls=:dash, lab="Empirical data")
xaxis!((0, 0.5), "p")
yaxis!((0, 7), "Density")
savefig(joinpath("figures", "penciltrick.png"))
