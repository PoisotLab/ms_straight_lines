\pagenumbering{gobble}

## Figure S1. Parameters can be estimated by Maximum Likelihood

While the full posterior distribution can be sampled using various bayesian
machinery, this is not necessary for obtaining point estimates of $\mu$
and $\phi$. A maximum likelihood estimate of each can be calculated by
rearranging eq. 4 and fitting a Beta distribution to the result:

![**Parameters can be estimated by Maximum Likelihood.** The maximum likelihood estimate of $p$ is compared to 20 samples from the posterior distribution of the flexible links model. The empirical distribution of $p$, obtained from all food webs archived on the `mangal.io` database, is also included.](figures/beta_fit.png){#fig:penciltrick}

We include this result because ecologists may wish to apply our methods for estimating $L$, $Co$ or $L/S$ without fitting a Bayesian posterior of their own. This approach loses information about the sample size of webs, but nevertheless provides a close match to both the empirical data and the bayesian posterior.

parameter  | MLE estimate  | MAP estimate
--|---|--
$\mu$  | 0.087  | 0.086 ± 0.0037
$\phi$  | 21.0  |  24.3 ± 2.4
