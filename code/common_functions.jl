
## automate the extraction of a Stan model into an infdata object from Arviz
function foodweb_model_output(model_chains)
    const_stan_infdata = from_cmdstan(model_chains,
        posterior_predictive = "y_hat",
        log_likelihood = "log_lik",
        coords = Dict("id" => d.id),
        dims = Dict(
            "L" => ["id"],
            "log_lik" => ["id"],
            "y_hat" => ["id"]
    ))
end

function write_posterior(chns, path)
    bb_df = DataFrame(chns)
    CSV.write(path, bb_df)
end


function plot_posterior_betabin(post, d = d)
    pp = post[:, r"y_hat"]

    pp_plot = scatter(d.nodes, [pp[15,i] for i in 1:ncol(pp)] .+ (d.nodes .- 1),
                      alpha = 0.01, color = "black", legend = false)

    for k in rand(1:2000, 45)
        scatter!(pp_plot, d.nodes, [pp[k,i] for i in 1:ncol(pp)] .+ (d.nodes .- 1),
                alpha = 0.01, color = "black", legend = false)
    end
    xaxis!(pp_plot, :log, "Species")
    yaxis!(pp_plot, :log, "Links")

    scatter!(d.nodes, d.links, color = "darkorange")
end

function plot_posterior(post, d = d)
    pp = post[:, r"y_hat"]

    pp_plot = scatter(d.nodes, [pp[15,i] for i in 1:ncol(pp)],
                      alpha = 0.01, color = "black", legend = false)

    for k in rand(1:2000, 45)
        scatter!(pp_plot, d.nodes, [pp[k,i] for i in 1:ncol(pp)],
                alpha = 0.01, color = "black", legend = false)
    end

    xaxis!(pp_plot, "Species", xlims = (3,60))
    yaxis!(pp_plot, "Links", ylims = (3,100))

    scatter!(d.nodes, d.links, color = "darkorange")
end
