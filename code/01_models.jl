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
    num_samples = 3000,
    name = "constant_connectance"
)
_, const_stan_chns, _ = stan(const_conn_stan_model, data_dict, summary = false);

const_stan_infdata = foodweb_model_output(const_stan_chns)



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
    num_samples = 3000,
    name = "powerlaw_connectance"
)

_, pwrlaw_stan_chns, _ = stan(pwrlaw_conn_stan_model, data_dict, summary = false);

pwrlaw_stan_infdata = foodweb_model_output(pwrlaw_stan_chns)


### beta-binomial connectance

const betabin_connectance = """
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
    theta ~ normal( 3,0.5 );
    a ~ beta( 3 , 7 );
    R ~ beta_binomial(F ,  a * exp(theta), (1 - a) * exp(theta) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] y_hat;
    for ( i in 1:W ) {
        log_lik[i] = beta_binomial_lpmf( R[i] | F[i] , a * exp(theta), (1 - a) * exp(theta)  );
        y_hat[i] = beta_binomial_rng(F[i] , a * exp(theta), (1 - a) * exp(theta) );
    }
}
"""


data_dict_simpler = Dict(
    "W" => length(d.id),
    "F" => d.nodes .^ 2 .- (d.nodes .- 1),
    "R" => d.links      .- (d.nodes .- 1))

bb_model = Stanmodel(
    model = betabin_connectance,
    nchains = 2,
    num_warmup = 1000,
    num_samples = 3000,
    name = "simple_betabin"
)


_, bb_chains , _ = stan(bb_model, data_dict_simpler, summary = false);


bb_chains_infdata = foodweb_model_output(bb_chains)

summary(bb_chains_infdata)

##### write out posterior samples as csvs

write_posterior(bb_chains, "data/beta_binomial_posterior.csv")

write_posterior(pwrlaw_stan_chns, "data/pwrlaw_posterior.csv")

write_posterior(const_stan_chns, "data/const_posterior.csv")


## could be used to plot a "ribbon"
bb_hpd = MCMCChains.hpd(bb_chains, alpha = 0.89)

bb_hpd.df

### get parameter estimates
summary(const_stan_infdata)
summary(pwrlaw_stan_infdata)
summary(bb_chains_infdata)

## calculate loo
loo(const_stan_infdata) #2798 +- 104

loo(pwrlaw_stan_infdata) # 2595 +- 49

loo(bb_chains_infdata) # 2543 +- 46


### this should produce pointwise loo calculations but I don't know how to get the numbers out
loo_pw = loo(pwrlaw_stan_infdata, pointwise = true)
