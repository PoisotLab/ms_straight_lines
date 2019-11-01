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

As identified by @Mart92, the number of links $L$ in a food web with $S$ nodes
can be no lower than $S-1$, and no higher than $S^2$. Another way of expressing
this idea is that because we observe a food web with $S$ species, we are
guaranteed to observed at least $S-1$ interactions. From a predictive
standpoint, this means that we need to figure out how much of the remaining
interactions, out of $S^2-(S-1)$, will be realized. Following @PoisCirt16, we
suggest that each of these interactions are instances of independent Bernoulli
trials with a set probability of success $p$. This means that the number of
predicted links can be expressed as:

$$
 \hat L = p\times\left[S^2-(S-1)\right]+(S-1)\,.
$$

This can be re-expressed as a second order polynomial:

$$
  \hat L = p\times S^2 + (1-p)\times S + (p-1)\,.
$$

Because we have an expression for the number of interactions, we can get an
expression for the expected connectance, which is

$$
  \frac{\hat L}{S^2} = p\frac{S^2}{S^2} + (1-p)\frac{S}{S^2}+(p-1)\frac{1}{S^2}
$$

This results in the connectance being expressed as

$$ \frac{\hat L}{S^2} = (p-1)\times S^{-2} + (1-p)\times S^{-1} + p \, .$$ {#eq:co}

# Fitting the model

# Practical consequences

## Connectance is constant (for large enough food webs)

Predictions from @eq:co

1. Connectance at large $S$ is $p$, as opposed to $0^+$ for the power law
1. Linkage density at large $S$ is $1$ (because it is $1+(p1)\times S^{-1}$), whereas it is $0^+$ for the power law

the fact that LD reaches 0 with the power law is concerning, since it means that
species are expected to stop interacting in large networks - our formulation
ensures that species will establish at least one interaction. This is lower than
the 2 suggested by Coheh **Ref**, but more realistic than the $0$ of Martinez

| Model     | $L$    | $Co$ | $LD$ | $\text{lim} Co$ | $\text{lim} LD$ |
|:----------|:-------|:-----|:-----|:----------------|:----------------|
| LSSL      | $2S$   |      |      |                 |                 |
| CC        | $aS^2$ |      |      |                 |                 |
| power law | $aS^b  |      |      |                 |                 |
| Bernoulli |        |      |      |                 |                 |

## 2 ...

## 3 ...
