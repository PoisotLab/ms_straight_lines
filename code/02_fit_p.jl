import Pkg
Pkg.activate(".")

import CSV
using Distributions
using DataFrames
using StatsPlots
using Statistics
using StatsBase
using DelimitedFiles

# Get the data
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]

# Add some information
d.p = d.links ./ (d.nodes.^2 .- (d.nodes .- 1))
d.ld = d.links ./ d.nodes
d.co = d.links ./ (d.nodes .* d.nodes)
d.exr = d.links .- (d.nodes .- 1)
d.exp = d.nodes.^2 .- (d.nodes .- 1)
d.pex = d.exr ./ d.exp

p = fit(Beta, d.pex)
writedlm("params.dat", params(p))
