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
interactions are transformed by their natural log, $a$ and $b$ can be identified
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
probability of success $p$, which much like $a$ and $b$ in equation @{eq:reg} is
assumed to be a constant across all food webs.

This means that the number of predicted links can be expressed as:

$$
 \hat L = p\times\left[S^2-(S-1)\right]+(S-1)\,.
$$

This can be re-expressed as a second order polynomial:

$$\hat L = p\times S^2 + (1-p)\times S + (p-1)\,. $${#eq:L}

# Fitting the model

# Practical consequences

## Connectance is constant (for large enough food webs)

Because we have an expression for the number of interactions, we can get an
expression for the expected connectance, which is

$$
  \frac{\hat L}{S^2} = p\frac{S^2}{S^2} + (1-p)\frac{S}{S^2}+(p-1)\frac{1}{S^2}
$$

This results in the connectance being expressed as

$$ \frac{\hat L}{S^2} = (p-1)\times S^{-2} + (1-p)\times S^{-1} + p \, .$$ {#eq:co}

**TP** Need to re-write this section, the table is wrong

Predictions from @eq:co

1. Connectance at large $S$ is $p$, as opposed to $0^+$ for the power law
1. Linkage density at large $S$ is $1$ (because it is $1+(p1)\times S^{-1}$), whereas it is $0^+$ for the power law

the fact that LD reaches 0 with the power law is concerning, since it means that
species are expected to stop interacting in large networks - our formulation
ensures that species will establish at least one interaction. This is lower than
the 2 suggested by Coheh **Ref**, but more realistic than the $0$ of Martinez

Note TP
> I think the table below is interesting, in that it shows that the
> predictions for the previous theories where making predictions that lack sense
> either for connectance or average degree -- our does not.

| Model     | interactions | connectance         | linkage density        | limit connectance | limit linkage |
|:----------|:-------------|:--------------------|:-----------------------|:------------------|:--------------|
| LSSL      | $2S$         | $2\times S^{-1}$    | $2$                    | $0$               | $2$           |
| CC        | $aS^2$       | $a$                 | $a\times S$            | $a$               | $\infty$      |
| power law | $aS^b$       | $a\times S^{(b-2)}$ | $a\times S^{(b-1)}$    | $0$               | $0$           |
| Bernoulli | Eqn. @eq:L   | Eqn. @eq:co         | $1+(p-1)\times S^{-1}$ | $p$               | $1$           |

## Only very-large food webs obey a power law

...

## We can derive a measure of departure from expected number of links

Because $p$ is the probability of a number $n = S^2 - (S-1)$ of independent
Bernoulli events, we can express the variance of the number of interactions in
excess of $(S-1)$ as $n\times p\times (1-p)$. This means that given a network
with observed species richness $S$ and observed links $L$, we can calculate its
$z$-score as

$$z = \frac{L - p\times \left[S^2-(S-1)\right]}{\sqrt{p\times (1-p)\times \left[S^2-(S-1)\right]}} \,.$$#{eq:z}

A network where $L = \hat L$ will have a $z$-score of 0, and any network with
more (fewer) interactions will have a positive (negative) $z$-score. This has
important practical consequences - the structure of ecological networks is often
probed for deviation from the random distribution of some measure of interest
(**ref Bascompte, Flores**), and most of these measures are in turn related to
connectance **ref P&G**; therefore, *to be continued*.

# Conclusions
