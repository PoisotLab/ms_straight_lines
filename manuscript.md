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
trials with a set probability of success $p$. This means that the number of predicted links can be expressed as

$$
 \hat L = p\times\left[S^2-(S-1)\right]+(S-1)
$$

# Fitting the model

# Results and discussion
