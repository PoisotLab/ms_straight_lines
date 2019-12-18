import Pkg; Pkg.activate(".")
import CSV
using CmdStan
using ArviZ

# include functions
include("common_functions.jl")

# define cmdstan location
set_cmdstan_home!(homedir() * "/cmdstan/")

# get the data and filter for predation only
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]

### constant connectance

const constantconnectance = """
data{
    int W;
    int L[W];
    int S[W];
}
parameters{
    real a;
    real phi;
}
model{
    vector[W] mu;
    phi ~ normal( 2 , 1 );
    a ~ normal( -3 , 1 );
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i];
    }
    L ~ neg_binomial_2( mu , exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] mu;
    vector[W] y_hat;
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i];
        log_lik[i] = neg_binomial_2_lpmf( L[i] | mu[i] , exp(phi) );
        y_hat[i] = neg_binomial_2_rng(mu[i], exp(phi));
    }
}
"""

data_dict = Dict("W" => length(d.id),
    "L" => d.links,
    "S" => d.nodes)

const_conn_stan_model = Stanmodel(
    model = constantconnectance,
    nchains = 4,
    num_warmup = 1000,
    num_samples = 1000,
    name = "constant_connectance"
)
_, const_stan_chns, _ = stan(const_conn_stan_model, data_dict, summary = false);

const_stan_chns

const_stan_infdata = foodweb_model_output(const_stan_chns)

const_stan_infdata

loo(const_stan_infdata)

summary(const_stan_infdata)

plot_density(stan_infdata, var_names=["mu"])


##### power law connectance

const pwrlaw_connectance = """
data{
    int W;
    int L[W];
    int S[W];
}
parameters{
    real a;
    real b;
    real phi;
}
model{
    vector[W] mu;
    phi ~ normal( 2 , 1 );
    a ~ normal( -3 , 1 );
    b ~ normal(1, 0.2);
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i] ^ b;
    }
    L ~ neg_binomial_2( mu , exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] mu;
    vector[W] y_hat;
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i] ^ b;
        log_lik[i] = neg_binomial_2_lpmf( L[i] | mu[i] , exp(phi) );
        y_hat[i] = neg_binomial_2_rng(mu[i], exp(phi));
    }
}
"""

pwrlaw_conn_stan_model = Stanmodel(
    model = pwrlaw_connectance,
    nchains = 4,
    num_warmup = 1000,
    num_samples = 1000,
    name = "powerlaw_connectance"
)

_, pwrlaw_stan_chns, _ = stan(pwrlaw_conn_stan_model, data_dict, summary = false);

pwrlaw_stan_infdata = foodweb_model_output(pwrlaw_stan_chns)


##### beta binomial model

const betabin_connectance = """
data{
    int W;
    int L[W];
    int S[W];
}
transformed data{
    int R[W];
    int F[W];
    for ( i in 1:W ) {
        F[i] = S[i] * S[i] - (S[i] - 1);
        R[i] = L[i] - (S[i] - 1);
    }
}
parameters{
    real<lower=0,upper=1> p;
    real<lower=1> phi;
}
model{
    phi ~ normal( 3 , 0.5 );
    p ~ beta( 1.54 , 9.49 );
    R ~ beta_binomial(F ,  p, exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] y_hat;
    for ( i in 1:W ) {
        log_lik[i] = beta_binomial_lpmf( R[i] | F[i] , p * exp(phi) , (1 - p) * exp(phi) );
        y_hat[i] = beta_binomial_rng(F[i] , p * exp(phi) , (1 - p) * exp(phi) );
    }
}
"""

betabin_conn_stan_model = Stanmodel(
    model = betabin_connectance,
    nchains = 4,
    num_warmup = 1000,
    num_samples = 1000,
    name = "betabin_connectance"
)

_, betabin_stan_chns, _ = stan(betabin_conn_stan_model, data_dict, summary = false);

betabin_stan_chns

betabin_stan_infdata = foodweb_model_output(betabin_stan_chns)

summary(betabin_stan_infdata)

loo(pwrlaw_stan_infdata)
loo(const_stan_infdata)


### betabin connectance simpler

const betabin_connectance_simpler = """
data{
    int W;
    int F[W];
    int R[W];
}
parameters{
    real<lower=0,upper=1> a;
    real<lower=0> theta;
}
model{
    vector[W] pbar;
    theta ~ exponential( 3 );
    a ~ beta( 1.54 , 9.49 );
    for(i in 1:W) pbar[i] = a;
    R ~ beta_binomial(F ,  pbar * theta, (1 - pbar) * theta );
}
generated quantities{
    vector[W] log_lik;
    vector[W] y_hat;
    for ( i in 1:W ) {
        log_lik[i] = beta_binomial_lpmf( R[i] | F[i] , a * theta, (1 - a) * theta  );
        y_hat[i] = beta_binomial_rng(F[i] , a * theta, (1 - a) * theta );
    }
}
"""


data_dict_simpler = Dict(
    "W" => length(d.id),
    "F" => d.nodes .^ 2 .- (d.nodes .- 1),
    "R" => d.links      .- (d.nodes .- 1))


betabin_conn_simpler_stan_model = Stanmodel(
    model = betabin_connectance_simpler,
    nchains = 2,
    num_warmup = 1000,
    num_samples = 1000,
    name = "betabin_connectance_simpler"
)


_, betabin_stan_chns, _ = stan(betabin_conn_simpler_stan_model, data_dict_simpler, summary = false);

betabin_stan_chns

betabin_stan_chns

betabin_stan_infdata = foodweb_model_output(betabin_stan_chns)

summary(betabin_stan_infdata)
