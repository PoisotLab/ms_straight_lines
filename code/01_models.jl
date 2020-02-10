import Pkg; Pkg.activate(".")
using CSV
using CmdStan
using ArviZ
using DelimitedFiles

# include functions
include("common_functions.jl")

# define cmdstan location
set_cmdstan_home!(homedir() * "/Desktop/cmdstan/")

# get the data and filter for predation only
d = CSV.read(joinpath(pwd(), "data", "network_data.dat"))
d = d[d.predation .> 0 , :]

### LSSL

const lssl = """
data{
    int W;
    int L[W];
    int S[W];
    int cf;
}
parameters{
    real a;
    real phi;
}
model{
    vector[W] mu;
    phi ~ normal( 2 , 1 );
    a ~ normal( 0.7 , 0.02 );
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i];
    }
    L ~ neg_binomial_2( mu , exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] mu;
    vector[W] y_hat;
    vector[cf] counterfactual_links;
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i];
        log_lik[i] = neg_binomial_2_lpmf( L[i] | mu[i] , exp(phi) );
        y_hat[i] = neg_binomial_2_rng(mu[i], exp(phi));
    }
    for (j in 1:cf){
     counterfactual_links[j] = neg_binomial_2_rng(exp(a) * j, exp(phi));
    }
}
"""

data_dict = Dict("W" => length(d.id),
    "L" => d.links,
    "S" => d.nodes,
    "cf" => 750)

lssl_stan_model = Stanmodel(
    model =lssl,
    nchains = 4,
    num_warmup = 1000,
    num_samples = 1000,
    name = "lssl"
)
_, lssl_stan_chns, _ = stan(lssl_stan_model, data_dict, summary = true);

### Constant connectance

const constant_connect = """
data{
    int W;
    int L[W];
    int S[W];
    int cf;
}
parameters{
    real<lower=0,upper=1> a;
    real phi;
}
model{
    vector[W] mu;
    phi ~ normal( 2 , 1 );
    a ~ beta( 3 , 7 );
    for ( i in 1:W ) {
        mu[i] = a * S[i]^2;
    }
    L ~ neg_binomial_2( mu , exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] mu;
    vector[W] y_hat;
    vector[cf] counterfactual_links;
    for ( i in 1:W ) {
        mu[i] = a * S[i] ^ 2;
        log_lik[i] = neg_binomial_2_lpmf( L[i] | mu[i] , exp(phi) );
        y_hat[i] = neg_binomial_2_rng(mu[i], exp(phi));
    }
    for (j in 1:cf){
     counterfactual_links[j] = neg_binomial_2_rng(a * j^2, exp(phi));
    }
}
"""

data_dict = Dict("W" => length(d.id),
    "L" => d.links,
    "S" => d.nodes,
    "cf" => 750)

const_connect_stan_model = Stanmodel(
    model =constant_connect,
    nchains = 4,
    num_warmup = 1000,
    num_samples = 1000,
    name = "constant_connect"
)
_, constant_connect_stan_chns, _ = stan(const_connect_stan_model, data_dict, summary = true);


##### power law connectance

const pwrlaw_connectance = """
data{
    int W;
    int L[W];
    int S[W];
    int cf;
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
    b ~ normal(2, 0.6);
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i] ^ b;
    }
    L ~ neg_binomial_2( mu , exp(phi) );
}
generated quantities{
    vector[W] log_lik;
    vector[W] mu;
    vector[W] y_hat;
    vector[cf] counterfactual_links;
    for ( i in 1:W ) {
        mu[i] = exp(a) * S[i] ^ b;
        log_lik[i] = neg_binomial_2_lpmf( L[i] | mu[i] , exp(phi) );
        y_hat[i] = neg_binomial_2_rng(mu[i], exp(phi));
    }
    for (j in 1:cf){
     counterfactual_links[j] = neg_binomial_2_rng(exp(a) * j^b, exp(phi));
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


### beta-binomial connectance
const betabin_connectance = """
data{
    int W;
    int L[W];
    int S[W];
    int cf;
}
transformed data{
    int F[W];
    int R[W];
    int M[W];
    for ( i in 1:W ) {
        M[i] = S[i] - 1;
        F[i] = S[i] * S[i] - M[i];
        R[i] = L[i]        - M[i];
    }
}
parameters{
    real<lower=0,upper=1> mu;
    real phi;
}
model{
    phi ~ normal( 3,0.5 );
    mu ~ beta( 3 , 7 );
    for (i in 1:W){
       target += beta_binomial_lpmf(  R[i] | F[i] ,  mu * exp(phi) , (1 - mu) * exp(phi));
    }
}
generated quantities{
    vector[W] log_lik;
    vector[W] y_hat;
    vector[cf] counterfactual_links;
    for ( i in 1:W ) {
        log_lik[i] = beta_binomial_lpmf( R[i] | F[i] , mu * exp(phi), (1 - mu) * exp(phi)  );
        y_hat[i] = beta_binomial_rng(F[i] , mu * exp(phi), (1 - mu) * exp(phi) ) + M[i];
    }
    for (j in 1:cf){
        counterfactual_links[j] = beta_binomial_rng(j * j - j + 1 , mu * exp(phi), (1 - mu) * exp(phi) ) + j - 1;
    }
}
"""


data_dict = Dict(
    "W" => length(d.id),
    "L" => d.links,
    "S" => d.nodes,
    "cf" => 750)

bb_model = Stanmodel(
    model = betabin_connectance,
    nchains = 2,
    num_warmup = 1000,
    num_samples = 1000,
    name = "betabin_connectance"
)


_, bb_chains , _ = stan(bb_model, data_dict, summary = true);

##### write out posterior samples as csvs

size(lssl_stan_chns)

lssl_df = DataFrame(lssl_stan_chns)
CSV.write("data/posterior_distributions/lssl.csv", lssl_df, delim=',')

constant_connect_array = DataFrame(constant_connect_stan_chns)
CSV.write("data/posterior_distributions/const_posterior.csv", constant_connect_array, delim=',')

pwrlaw_df = DataFrame(pwrlaw_stan_chns)
CSV.write("data/posterior_distributions/powerlaw_posterior.csv", pwrlaw_df, delim=',')

bb_df = DataFrame(bb_chains)
CSV.write("data/posterior_distributions/beta_binomial_posterior.csv", bb_df, delim=',')



### more posterior samples for Flexible links model

data_dict_2 = Dict(
    "W" => length(d.id),
    "L" => d.links,
    "S" => d.nodes,
    "cf" => 1500)


bb_model = Stanmodel(
        model = betabin_connectance,
        nchains = 2,
        num_warmup = 1000,
        num_samples = 1000,
        name = "betabin_connectance_bigger"
    )


_, bb_chains , _ = stan(bb_model, data_dict_2, summary = true);


## could be used to plot a "ribbon"
bb_hpd = MCMCChains.hpd(bb_chains, alpha = 0.89)

bb_hpd.df



###  calculation of LOO

bb_chains_infdata = foodweb_model_output(bb_chains)

lssl_stan_infdata = foodweb_model_output(lssl_stan_chns)

const_stan_infdata = foodweb_model_output(constant_connect_stan_chns)

pwrlaw_stan_infdata = foodweb_model_output(pwrlaw_stan_chns)


### get parameter estimates
summary(const_stan_infdata)
summary(pwrlaw_stan_infdata)
summary(bb_chains_infdata)
summary(lssl_stan_infdata)


loo(lssl_stan_infdata)
# 2940.68 +- 118.63
## calculate loo
loo(const_stan_infdata) #2798 +- 104
# 2940 +- 118 with a narrower prior

loo(pwrlaw_stan_infdata) # 2595 +- 49

loo(bb_chains_infdata) # 2543 +- 46


### this should produce pointwise loo calculations but I don't know how to get the numbers out
loo_pw = loo(pwrlaw_stan_infdata, pointwise = true)
