# 1 - Beta fit


## Parameters can be estimated by Maximum Likelihood

While the full posterior distribution can be sampled using various bayesian
machinery, this is not necessary for obtaining point estimates of $p$
and $\phi$. A maximum likelihood estimate of each can be calculated by
rearranging equation {#eq:lhat} and fitting a Beta distribution to the result:

![Parameters can be estimated by Maximum Likelihood](figures/beta_fit.png){#fig:penciltrick}

We include this result because ecologists may wish to apply our methods for estimating $L$, $Co$ or $L/S$ without fitting a Bayesian posterior of their own. This approach loses information about the sample size of webs, but nevertheless provides a close match to both the empirical data and the bayesian posterior.

parameter  | MLE estimate  | MAP estimate
--|---|--
μ  | 0.087  | 0.086 ± 0.0037
ϕ  | 21.0  |  24.3 ± 2.4




# 2 - Bigger picture box

Understanding the functions, development, and dynamics of ecological communities requires the investigation of the structure of their networks of interactions. The most challenging issue in the study of such networks is that the sampling of ecological interactions is a strenuous task, and as a result the availability of empirical data is limited. In this contribution we derive a realistic and performant statistical model to predict the number of interactions in a food web from its species richness.

Our model could be used as a first order approximation of network structure. As such, it makes the large-scale study and comparison of ecological networks more accurate and assessible. For instance, the stability and resilience of ecological communities to perturbations could be explored at the regional, continental or global scales. This is particularly relevant in the actual context of climate change and ecosystems' collapse.
