
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


neg_to_zeros(quant) = quant > 0 ? quant : 0.01 # to deal with negative quantiles
