# Introduction

The L/S relationship is a fundamental one for food webs ecology

- Early predictions differ in whether this is a linear of exponential relationship
- General case is $L = aS^b$, which is easy to fit
- Has consequences for spatial scaling (Brose)
- Parameters are difficult to reason about ecologically

We know that there are hard boundaries to this system

- Cannot have more than $S^2$ interactions
- Cannot have fewer than $S-1$
- This was not used in the previous attempts to fit the relationship

We know that interactions can be viewed as stochastic events

- Interaction is the outcome of a Bernoulli process
- This lends itself to powerful approaches in probabilistic programming

In this paper

- New relationship between L and S
- Discuss how it changes network prediction
- Use it as a story telling device for a broader point about using our knowledge of the system in creative ways

# Deriving the model

@CoheBria84 hypothesized that all networks would have the same average degree,
resulting in link-species scaling expressed as

$$\hat L_\text{LSSL} = b\times S\,,$${#eq:lssl}

where $S$ is the species richness, and $b \approx 2$. @Mart92 instead suggested
that most networks should have constant connectance, expressed as $L/S^2$, and
therefore one can predict the number of interactions as

$$\hat L_\text{CC} = b\times S^2\,,$${#eq:cc}

where $b$ is a constant in $]0,1[$. Finally, @BrosOstl04 note that these two
models are instead parameterizations of the same general model, in which

$$\hat L_\text{reg} = b\times S^a\,, $${#eq:reg}

where $a$ and $b$ are constants. When the number of links and number of
interactions are transformed by their natural log, $a$ and $b$ can be estimated
with a linear regression, as done by @Mart92.

Although all of these models fit the data well enough, they neglect a
fundamental piece of ecological knowledge about food webs: as identified by
@Mart92, the number of links $L$ in a food web with $S$ nodes can be no lower
than $S-1$, and no higher than $S^2$. Another way of expressing this idea is
that because we observe a food web with $S$ species, we are guaranteed to
observed at least $S-1$ interactions. From a predictive standpoint, this means
that we need to figure out how much of the remaining interactions, out of
$S^2-(S-1)$, will be realized. Following @PoisCirt16, we suggest that each of
these interactions are instances of independent Bernoulli trials with a set
probability of success $p$, which much like $a$ and $b$ in +@eq:reg is assumed
to be a constant across all food webs.

This means that the number of predicted links can be expressed as:

$$
 \hat L = p\times\left[S^2-(S-1)\right]+(S-1)\,.
$$

This can be re-expressed as a second order polynomial:

$$\hat L = p\times S^2 + (1-p)\times S + (p-1)\,. $${#eq:L}

# Fitting the model

We have rephrased the question of connectance in food webs as the proportion of
links realized above the minimum. We use a Binomial likelihood where we consider
the number of links above the minimum as 'successes' and the number of links
between the minimum and maximum as the number of 'trials':

$$\begin{aligned}
L_i & \sim \text{Binomial}(\left[S_i^2-(S_i-1)\right], p_i)\\
\text{logit}(p_i) &\sim \text{Normal}(\mu_p, \sigma_p)\\
\mu_p & \sim \text{Normal}()\\
\sigma_p & \sim \text{Exponential}()
\end{aligned}$$

Note that this model has no deterministic component for $p$, since it is modeled
as a constant. We assume that $p$ may be described by a normal distribution on
the logit scale. The parameter $\mu_p$ replaces previous estimates of the
average connectance across all food webs. However, the variation among food webs
is not completely captured by $S$ alone, and the variation in link number is
greater than expected in a binomial distribution. This overdispersion is
captured in the hyperparameter $\sigma_p$, which also partially pools estimates
of $p_i$ towards the average value. This increases accuracy, both within sample
and when making predictions for new webs.

We selected our prior on $\mu_p$ to reflect previous estimates for the constant of connectance: we calculated the logit of @Mart92 's value. However, as no information is available for either the standard deviation of this distribution nor for $\sigma_p$, we followed the advice of (tk Simpson et al), and performed prior predictive checks. The  _tk actual values_

Two quantities are interesting to calculate from this model. First we may calculate the MAP (maximum a posteriori) estimate of average $p$ across all webs, by using the inverse logit of $\mu_p$:

$$p = \frac{e^{\mu_p}}{1 - e^{\mu_p}}$$

Secondly, we may wish to provide the predicted value of $p$ for a new web. Because

We use Stan (**tk version, ref**) which implements Bayesian inference using Hamiltonian Monte Carlo.

# Results and discussion

## The Bernoulli-based model outperforms previous solutions

**Andrew** this one is for you to write

## Connectance is constant (for large enough food webs)

Because we have an expression for the number of interactions (+@eq:L), we can
get an expression for the expected connectance, which is

$$
  \frac{\hat L}{S^2} = p\frac{S^2}{S^2} + (1-p)\frac{S}{S^2}+(p-1)\frac{1}{S^2} \,.
$$

This results in the connectance being expressed as

$$ \frac{\hat L}{S^2} = (p-1)\times S^{-2} + (1-p)\times S^{-1} + p \, .$${#eq:co}

Note that the expression of connectance is no longer a polynomial; at large
values of $S$, the terms in $S^{-1}$ and $S^{-2}$ will tend towards 0, and so
the connectance will converge towards $p$. Therefore, for large enough
ecological networks, we should expect to observe a connectance that behaves more
or less in a constant way. This result provides an interesting ecological
interpretation of $p$, namely that it represents the connectance which we expect
a network large enough for the effect of $(S-1)$ minimum links to be negligible.

Interestingly, this model still results in an expected average degree ($\hat
L/S$, the *linkage density*) for a large number of species that scales with $S$:

$$\frac{\hat L}{S} = p\times S + (1-p) + (p - 1)\times S^{-1}\,.$$

This means that the addition of $p^{-1}$ new species should increase the average
degree in the food web by 1; of course, for increasingly large values of $S$,
this may result in unrealistic average degree, as species are limited by
biological mechanisms such as handling time, capture efficiency, etc, in the
number of interactions they can establish. Seeing how most ecological networks
are reasonably small, this does not look like an unreasonable assumption.

## Only very-large food webs obey a power law

As noted by @BrosOstl04, the models of @CoheBria84 and @Mart92 results in
networks in which the relationship between $L$ and $S$ obeys a power-law, albeit
with different parameters. Our model does not make this prediction, due to the
fact that we explicitly account for the lower bound of $(S-1)$ interactions. In
+@eq:L, the term $p\times S^2$ will become increasingly important when $S$
increases, and so we can quantify the extent to which the relationship gets
closer to a power law when $S$ increases.

We do so by dividing the terms with exponents lower than 2 by the term with
exponent 2, which gives

$$k = \frac{(1-p)\times S + (p-1)}{p\times S^2}\,.$$

This will peak for small values of $S$, and then slowly decrease towards 0. We
illustrate these results in **FIGURE**, which reveals that for networks under
**500?** species, the relationship between $S$ and $L$ strongly deviates from a
power-law. Specifically, **complete -- is it under or over-estimating the number
of interactions?**.

## We can derive a measure of departure from expected number of links

Because $p$ is the probability of a number $n = S^2 - (S-1)$ of independent
Bernoulli events, we can express the variance of the number of interactions in
excess of $(S-1)$ as $n\times p\times (1-p)$. This means that given a network
with observed species richness $S$ and observed links $L$, we can calculate its
$z$-score as

$$z = \frac{(L -(S-1)) - p\times \left[S^2-(S-1)\right]}{\sqrt{p\times (1-p)\times \left[S^2-(S-1)\right]}} \,.$${#eq:z}

A network where $L = \hat L$ will have a $z$-score of 0, and any network with
more (fewer) interactions will have a positive (negative) $z$-score. This has
important practical consequences - the structure of ecological networks is often
probed for deviation from the random distribution of some measure of interest
(**ref Bascompte, Flores**), and most of these measures are in turn related to
connectance **ref P&G**; therefore, *to be continued*.

# Conclusions
