# Introduction

Community ecologists are fascinated by counting things. It is therefore no
surprise that the early food web literature paid so much attention to counting
species, counting trophic interactions, and uncovering the shape of the
relationship that binds them -- and it is undeniable that these inquiries
kickstarted what is now one of the most rapidly growing fields of ecology
[@BorrMood14]. More species always means more interactions; this scaling between
species richness $S$ and number of interactions $L$ is universal and appears
both in observed food webs and under purely neutral models of food web structure
[@CanaMouq12]. In fact, these numbers underlie most measures used to describe
the structure of a food web [@DelmBess18]. The structure of a food web, in turn,
is almost always required to understand how the community functions, develops,
and responds to changes [@McCa12; @ThomBros12], to the point where some authors
suggested that describing food webs was a necessity [@SeibCado18; @McCa07]. To
this end, a first step is to come up with an estimate for the number of existing
trophic interactions $L$, through sampling or otherwise. Although both $L$ and
$S$ can be counted in nature, the sampling of interactions is orders of
magnitude more difficult than the sampling of species [@Jord16a; @Jord16]. As a
result, we have far more information about $S$. In fact, the distribution of
species richness across the world is probably the most frequently observed and
modelled ecological phenomena. Therefore, if we can predict $L$ from $S$ in an
ecologically realistic way, we will be in a position to make first order
approximations of food web structure at large scales.

Measures of food web structure are based on three specific quantities. The first
and most straightforward is $L$, the number of trophic interactions among
species. This quantity can be large, especially in species-rich habitats, but it
cannot be arbitrarily large. It is clear to any observer of nature that of all
imaginable trophic interactions, only a fraction actually occur. If an
ecological community contains $S$ species, then the maximum number of links in
its foodweb is $S^2$: a community of omnivorous cannibals. This leads to the
second quantity: a ratio called *connectance* and defined by ecologists as $Co =
L/S^2$. Connectance has become a fundamental quantity for nearly all other
measures of food web structure and dynamics [@PascDunn06]. The third important
quantity is another ratio: *linkage density*, $L_D = L/S$. This value represents
the number of links added to the network for every additional species in the
ecological system. A closely related quantity is $L_D \times 2$, which is the
_average degree_: the average number of species with which any taxa is expected
to interact, either as predators or prey. All of these quantities capture
ecologically important aspects of a network, and all capture the information
present in a prediction of $L$ links among $S$ species. Accurate predictions of
ecological networks are extremely useful in many ecological contexts; thus it is
important to have an ecologically accurate predictive model for the underlying
value, $L$.

Because $L$ represents such a fundamental quantity, many predictive models have
been considered over the years. Here we describe three popular approaches before
describing our own proposed model. The *link-species scaling (LSSL)* model
introduced by @CoheBria84 hypothesized that all networks have the same average
degree, thats is most species should have the same number of interactions. Links
are modeled as the number of species times a constant: $\hat L_\text{LSSL} =
b\times S$, with $b \approx 2$. This model imagines that every species added to
a community increases the number of interactions by two -- for example, an
animal which consumes one resource and is consumed by one predator. Yet this
model started to show its deficiencies when data on larger food webs became
available, which revealed that $L$ increases faster than a linear function of
$S$ would. In response to the idea that the *average degree* is constant,
Martinez [@Mart92] suggested instead that *connectance* is unchanged in response
to $S$; in other words, a food web is always equally filled, regardless of
whether it has 5 or 5000 species. Under the so-called "constant connectance"
model, the number of links is proportional to the richness squared, $\hat
L_\text{CC} = b\times S^2$, where $c$ is a constant in $]0,1[$ representing the
expected value of connectance. This model can be relaxed by assuming that the
scaling of $L$ with $S$ does not necessarily follows the maximum number of
interactions, and the best fit was with a model of the form $\hat L_\text{reg} =
b\times S^a$, which is trivially a linear relationship between $\text{log}(L)$
and $\text{log}(S)$. This power law model can be parameterized in arbitrarily
complex ways, including spatial scaling and species area relationships
[@BrosOstl04]; it should further be noted that this model is a synthesis of
preview hypotheses, encompassing both the link-species scaling ($a=1, b\approx
2$) and the strict constant connectance ($a, b=2$) depending on which parameters
are fixed. Power laws are very flexible, and indeed this function matches
empirical data well -- so well that it is often treated as a "true" model which
captures the scaling of link number with species richness [@WinePian01;
@RiedRall10; @GarlCald03], and from which we should draw ecological inferences
about what shapes food webs. However, this approach is limited, because the
parameters of a power law relationship can arise from many mechanisms, and are
difficult to reason about ecologically.

But the question of how informative parameters of a power law can be is moot.
Indeed, both the general model and its variants share an important shortcoming:
they cannot be used for prediction while remaining within the bounds set by
ecological principles. In short, while they can describe the *data* adequately,
they are fundamentally unable to represent the mechanisms through which these
data emerged. This has two causes. First, models that are variations of $\hat L
\approx b\times S^a$ have no constraints --  with the exception of the "constant
connectance" model, in which $L_{CC}$ has a maximum value of $S^2$. However, we
know that the number of interactions within a food web is both lower and upper
bounded [@Mart92; @PoisGrav14]: there can be no more than $S^2$ links, and there
can be no fewer than $S-1$ links. This minimum of $S-1$ represents communities
where all species interact and at least some of the organisms are heterotrophs
[@Mart92]. Numerous simple foodwebs could have this minimal number of links --
for example, a linear food chain wherein each trophic level consists of a single
species, each of which consumes only the species below it; or a grazing
herbivore which feeds on every plant in a field. Thus interaction number is
constrained by ecological principles to be between $S^2$ and $S-1$, something
which no present model includes. Secondly, accurate predictions of $L$ from $S$
are often difficult because of how parameters are estimated. This is usually
done using a Gaussian likelihood for $L$, often after log transformation or both
$L$ and $S$. While this approach ensures that predicted values of $L$ are always
positive, it does nothing to ensure that they stay below $S^2$ and above $S-1$.
Thus a good model for $L$ should meet these two needs: a bounded expression for
the mean $\hat{L}$, as well as a bounded distribution for the likelihood.

Here we suggest a new perspective for a model of $L$ as a function of $S$ which
respects ecological bounds, and has a bounded distribution of the likelihood. We
include the minimum constraint by modelling not the total number of links, but
the number in excess of the minimum. We include the maximum constraint in a
similar fashion to the constant connectance model described above, by modelling
the proportion of flexible links which are realized in a community. In this
contribution, we show how this model not only outperforms existing efforts at
predicting the number of interactions, but also has numerous desirable
properties from which novel insights about the structure of food webs can be
derived.

# Interlude - deriving a process-based model for the number of links

Based on the ecological constraints discussed earlier, we know that $L$ is a

$$
 \hat L_{\textsc{fl}} = p\times\left[S^2-(S-1)\right]+(S-1)\,,
$$ {#eq:lhat}

where $p \in [0,1]$. When $p = 1$, $L$ is at its maximum ($S^2$), and when $p =
0$ it is at the minimum value ($S - 1$). We use the notation $L_{FL}$ to
represent that our model considers the number of "flexible" links in a food web;
that is, the number of links in excess of the minimum but below the maximum.

Our second contribution is to propose an improved means of fitting this model.
Following @PoisCirt16, we suggest that each flexible link represents an
independent Bernoulli trial with a probability  $p$ of existing, and that an
observation of $L_i$ links in community $i$ represents an aggregation of $S^2 -
(S - 1)$ such trials. If we then assume that $p$ is a constant for all links in
the same ecological community, but may vary between communities, we can model
the distribution of links directly as a shifted beta-binomial variable.

$$
[L|S,\mu, \phi] =  { S^2 - (S - 1) \choose L - (S - 1)} \frac{B(L - (S - 1) + \mu \phi, S^2 - L + (1 - \mu)\phi)}{B(\mu \phi, (1 - \mu)\phi)}
$$ {#eq:shiftBB}

Where $B$ is the Beta function, $\mu$ is the average probability of a flexible
link being realized (i.e. the average value of $p$) and $\phi$ is the
concentration around this value. The support of this distribution is limited to
only ecologically realistic values of $L$: the number of species determines the
number of possible links and it is shifted to the right by $S-1$. This means
that the problem of estimating values for $\mu$ and $\phi$ is reduced to fitting
the univariate distribution described in +@eq:shiftBB. For more detailed
explanation of the model derivation, fitting, and comparison, see Experimental
Procedures.

In this paper we will compare our flexible links model to three previous models
for $L$. We estimate parameters and compare the performance of all models using
open data from the `mangal.io` networks database. Finally, we will show how our
model for $L_{FL}$ suggests a new and more useful way of thinking about metrics
of network structure and discuss how generative models can be useful tools for
including our knowledge of a system into our predictions.

# Results and Discussion

### Flexible link model fits better and makes a plausible range of predictions

| model                | reference     | PSIS-LOO         | $\Delta \text{elpd}$ | $SE_{\Delta \text{elpd}}$ |
| -------------------- | ------------- | ---------------- | -------------------- | ------------------------- |
| Flexible links       | current paper | 2520.5 ± 44.4    | 0                    | 0                         |
| Power law            | @BrosOstl04   | 2564.3 ± 46.6    | -21.9                | 6.5                       |
| Constant             | @Mart92       | 2811.0 ± 68.3    | -145.3               | 21.1                      |
| Link-species scaling | @CoheBria84   | 39840.1 ± 2795.1 | -18659.8             | 1381.7                    |

Table: Comparison of the different modles. Pareto-smoothed important sampling values and differences relative to the maximum in the expected log predictive density for the flexible link and the three competing models. The mean and standard deviation (SD) (standard error (SE)) is given for the two metrics. {#tbl:comparison}

All models fitted well, without any problematic warnings from Stan's diagnostics
(see Experimental Procedures), but our model for flexible links outperformed
previous solutions to the problem of modelling $L$. The flexible link model,
which we fit via a beta-binomial observation model, had the most favourable
values of PSIS-LOO information criterion (Table 1) and of expected log
predictive density (ELPD), relative to the three competing models which used a
negative binomial observation model. Pareto-smoothed important sampling serves
as a guide to model selection; like other information criteria it approximates
the error in cross-validation predictions. Smaller values indicate a model which
makes better predictions. The calculation of PSIS-LOO can also furnish some
clues about potential model fits; in our case the algorithm suggested that the
constant connectance model was sensitive to extreme observations. Information
criteria are only a rough guide to model selection; as always domain expertise
should take precedence. The expected log predictive density (ELPD), on the other
hand, measures the predictive performance of the model; here, higher values
indicate more reliable predictions. This suggests that the flexible link model
will make the best predictions of $L$.

Useful predictions for $L$ must however stay within realistic boundaries
determined by ecological principles. We generated posterior predictions for all
models and visualized them against these constraints (@fig:PP_counterfactual).
The LSSL model clearly underestimated the number of links, especially in large
networks: its predictions were frequently lower than the minimum $S-1$. The
constant connectance and power law models also made many predictions below this
value, especially for small values of $S$. The flexible link model made roughly
the same predictions, but within ecologically possible values.

![**The flexible link model fits better and makes a plausible range of predictions.** The number of links is plotted as a function of species richness obtained from the posterior distributions of A) the link-species scaling, B) the constant connectance, C) the power law and D) the flexible link models. In each panel, the colored line represent the median predicted link number and the grey areas cover the 78% and 97% percentile intervals. Empirical data from the `mangal.io` database are plotted in each panel (grey dots), as well as the minimal $S-1$ and maximal $S^2$ number of links (lower and upper black lines, respectively).  ](figures/models_links.png){#fig:PP_counterfactual}

### Flexible link model makes realistic predictions for small communities

The constraints on food web structure are especially important for small
communities. This is emphasized in +@fig:real_predict, which shows that all
models other than the flexible links model fail to stay within realistic
ecological constraints. The link-species scaling model made around 29% of
unrealistic predictions of link numbers for every value of $S$ ($3 \leq S \leq
750$). The constant connectance and power law models, on the other hand, also
produced unrealistic results but for small networks only: more than 20% were
unrealistic for networks comprising less than 12 and 7 species, respectively.
Only the flexible link model never failed to predict numbers of links between
$S-1$ and $S^2$.  

![**Only the flexible link model makes realistic predictions for small communities.** Here we show the proportion of posterior predictions from each of our 4 models which fall outside ecologically realistic values. The proportion of predictions in the correct range increases with species richness for the constant connectance and LSSL models. Vertical lines show the 5%, 50% and 95% quantiles of the distribution of S, demonstrating that many communities have potentially incorrect predictions under previous models](figures/real_predict.png){#fig:real_predict}

### Parameter estimates for all models

| model                | Equation for $\hat{L}$     | parameter | interpretation                       | value | SD     |
|----------------------|----------------------------|-----------|--------------------------------------|-------|--------|
| Link-species scaling | $bS$                       | $b$       | number of links per species          | 2.2   | 0.047  |
|                      |                            | $\kappa$  | concentration of $L$ around mean     | 1.4   | 0.12   |
| Constant connectance | $bS^2$                     | $b$       | proportion of maximum links realized | 0.12  | 0.0041 |
|                      |                            | $\kappa$  | concentration of $L$ around mean     | 4.0   | 0.37   |
| Power law            | $bS^a$                     | $b$       | controls proportion of relationship  | 0.37  | 0.054  |
|                      |                            | $a$       | controls scaling of relationship     | 1.7   | 0.043  |
|                      |                            | $\kappa$  | concentration of $L$ around mean     | 4.8   | 0.41   |
| Flexible links       | $(S^2 - (S - 1))p + S - 1$ | $\mu$     | average value of $p$                 | 0.086 | 0.0037 |
|                      |                            | $\phi$    | concentration around value of $\mu$  | 24.3  | 2.4    |

Table: Parameter estimates for all models. Mean and standard deviation (SD) is given for each parameter. {#tbl:parameters}

Although we did not use the same approach to parameter estimation as previous
authors, our approach to fitting these models recovered parameter estimates that
are broadly congruent with previous work. We found very consistent values of 2.2
for $b$ of the LSSL model, which is close to the value of approximately 2 used
by @CoheBria84. Similarly, we found a value of 0.12 for $b$ of the constant
connectance model, which was consistent with the 0.14 found by @Mart92. Finally,
the parameters values we found for the power law were also comparable to the
ones found by @BrosOstl04.

For large communities, our model should behave similarly to the constant
connectance model; it is no surprise then than $\mu$ was about 0.09, which is
close to our value of 0.12 for constant connectance. In addition, we obtained a
rather large value of 24.3 for $\phi$, which shrinks the variance around the
mean of $p$ to approximately 0.003 ($var(p)=\mu(1-\mu)/(1+\phi)$). This
indicates that food webs are largely in their probabiliy of flexible links being
realized.

<!-- tk add ecological interpretation where possible -->

### Connectance and linkage density can be derived from a model for links

Of the three important quantities which describe networks ($L$, $Co$ and $L_D$)
we have directly modelled $L$ only. However, we can reuse the posterior for our
model of $L$ to parameterize a distribution for connectance ($L/S^2$) and
linkage density ($L/S$). We can derive this by noticing that +@eq:lhat can be
rearranged to show how $Co$ and $L_D$ are linear transformations of $p$:

$$ Co = \frac{\hat L}{S^2} = p\left(1 - \frac{S-1}{S^2}\right) + \frac{S-1}{S^2} ,$$ {#eq:co2}

and

$$ L_D = \frac{\hat L}{S} = p \left(S - \frac{S-1}{S} \right) +  \frac{S-1}{S},$$ {#eq:ld}

In a Beta-Binomial distribution, it is assumed that the probability of success
$p$ varies among groups of trials according to a Beta distribution. Since $p$
has a Beta distribution the linear transformations described by +@eq:co2 and
+@eq:ld also describe Beta distributions which have been shifted and scaled
according to the number of species $S$ in a community.

#### Discussion of these equations

For large ecological systems, where $S$ has a high value, +@eq:co2 and +@eq:ld
respectively approach

$$Co \approx p $$

$$ L_D \approx pS $$

This implies that the addition of $n$ species should increase the linkage
density by approximately $p\times n$. For example, the addition of 11 new
species ($p^{-1}$) should increase the linkage density in the food web by
roughly 1, meaning that each species in the original network would be expected
to develop 2 additional interactions. Of course, for increasingly large values
of $S$, this may result in an unrealistic average degree, as species are limited
by biological mechanisms such as handling time, capture efficiency, _etc_, in
the number of interactions they can establish. Most ecological networks are
however reasonably small and so this does not look like an unreasonable
assumption.

Note also that the expression of connectance is no longer a polynomial; at large
values of $S$, the value of $(S-1)/S^2$  will tend towards 0, and so the
connectance will converge towards $p$. Therefore, for large enough ecological
networks, we should expect a connectance which is independent of $S$. Thus $p$
has an interesting ecological interpretation: it represents the average
connectance of networks large enough that the proportion $(S-1)/S^{2}$ is
negligible.

In previous work it has been debated whether $Co$ is a constant across all
communities. However, some authors have found that this does not hold for small
communities. This result makes it clear that this confusion comes from
neglecting the minimum value.

#### Discussion of the shifted beta distributions.

Just as $L$ must be within ecologically meaningful bounds, $Co$ and $L_D$ must
be as well. The connectance of a food web is bounded by 0 and 1. However, the
minimum bound on links similarly imposes a lower value on connectance. This
means that the distribution for $Co$ will be a shifted beta distribution, a
transformed version of the distribution for $p$.

We can convert the distribution for $p$ into one for $Co$ by replacing $p$ with
the transformation of $Co$ as described above (+@eq:co2), and rescaling by the
new range:

$$
[Co | S, \mu, \phi] = \frac{\left(Co - \frac{S-1}{S^2}\right)^{\mu \phi - 1}\left(1 - Co\right)^{(1 - \mu)\phi - 1} }{(1 - \frac{S-1}{S^2})^{\phi - 1} \times B(\mu \phi, (1 - \mu)\phi)}
$$ {#eq:shiftBetaCo}


Similarly, we can convert the distribution for $p$ into one for $L_D$ by
replacing $p$ with the transformation of $L_D$ (+@eq:ld)


$$
[L_{D} | S, \mu, \phi] = \frac{\left(L_D - \frac{S-1} {S}\right)^{\mu \phi - 1}\left(1 - L_D\right)^{(1 - \mu)\phi - 1} }{(S - \frac{S-1}{S})^{\phi - 1} \times B(\mu \phi, (1 - \mu)\phi)}
$$ {#eq:shiftBetaLD}

In +@fig:beta_distributions, we show that the connectance and linkage density
obtained from the equations above fitted the empirical data well. Their
predictions did not exceed ecological boundaries (between $(S-1)/S^2$ and 1 for
connectance, and between $(S-1)/S$ and $S$ for the linkage density).

![**Connectance and linkage density can be derived from a model for links.** A) Connectance and B) linkage density are plotted as a function of species richness, for the maximum _a posteriori_ estimates of the flexible link model. In each panel, the colored line represent the median predicted quantity and the grey areas cover the 78% and 97% percentile intervals. Empirical data from the `mangal.io` database are plotted in each panel (grey dots). In A), the minimal $(S-1)/S^2$ connectance and in B) the minimal $(S-1)/S$ and maximum $S$ linkage density are plotted (black lines).](figures/connectance_linkdens.png){#fig:beta_distributions}

Connectance is more than the proportion of realized interactions. It has been
associated with some of the most commonly used network metrics (@PoisGrav14,
@Chag15), and contains meaningful information on the stability (@Dunn02,
@Mont06) and dynamics (@VierAlme15) of ecological communities. A probability
distribution for connectance non only accounts for the variability between
networks, but can be used to describe fundamental properties of food webs and to
identify ecological and evolutionary mechanisms shaping communities.

## Only very large food webs obey a power law

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
illustrate these results in +@fig:powerlawk, which reveals that for networks
under approximatively 120 species, the relationship between $S$ and $L$ strongly
deviates from a power law ($k > 0.1$). In large networks, the terms with
exponents lower than 2 in +@eq:L are negligible when compared to the power law
$p\times S^2$. The power law model therefore considerably under-estimates the
number of interactions, especially for small networks.

This model sheds some light on a classical result by @DunnWill02a: ecological
networks deviate most strongly from the expectations under "small world" or
"scale free" regimes when they exceed a certain connectance; this is because for
small networks, connectance is higher, and only decreases towards $p$ when the
term in $S^{-2}$ in +@eq:co vanishes.

![**Only very large food webs obey a power law.** The extent $k$ to which the relationship between $L$ and $S$ deviates from a power law is plotted as a function of species richness, obtained from the maximum _a posteriori_ estimates of the flexible link model. The colored line represent the median predicted $k$ and the grey areas cover the 78% and 97% percentile intervals.](figures/k_powerlaw.png){#fig:powerlawk}

## Normal approximation provides an analytic z-score

Ecologists are often faced with the issue of comparing several networks. Often,
they wish to know if the network they have is "unusual" relative to some
expectation. Traditionally these comparisons have been done by constructing a
Null distribution . But here we propose a means of doing so with math.

The shifted beta-binomial can be approximated by a normal distribution:

$$ L \sim Normal(\bar{L}, \sigma_L^2) $$

$$ \bar{L} = (S^2 - S + 1) \mu + S - 1$$

$$ \sigma_L^2 = (S^2 - S + 1) \mu (1 - \mu)(1 + \frac{S(S-1)}{\phi + 1})$$

This means that given a network with observed species richness $S_{obs}$ and
observed links $L_{obs}$, we can calculate its $z$-score as

$$z = \frac{L_{obs} - \bar{L}}{\sqrt{\sigma_L^2}} \,.$$ {#eq:z}

A network where $L = \hat L$ will have a $z$-score of 0, and any network with
more (fewer) interactions will have a positive (negative) $z$-score. This has
important practical consequences - the structure of ecological networks is often
probed for deviation from the random distribution of some measure of interest
(**ref Bascompte, Flores**), and most of these measures are in turn related to
connectance **ref P&G**. We suggest that the use of a $z$-score could help
identify significantly under (over) sampled networks and estimate their number
of missing (extra) links.

In (+@fig:MAPnormal), we show that the predictions made by the normal
approximation (panel B) are similar to those made by the beta distribution
parameterized with the maximum _a posteriori_ values of $\mu$ and $\phi$ (panel
A).

![**The shifted beta-binomial distribution can be approximated by a normal distribution.** The number of links is plotted as a function of species richness obtained from A) the maximum _a posteriori_ estimates of the flexible link model and B) its normal approximation. In each panel, the colored line represent the median predicted link number and the grey areas cover the 78% and 97% percentile intervals (only the 78% percentile interval is depicted in B). Empirical data from the `mangal.io` database are plotted in each panel (grey dots), as well as the minimal $S-1$ and maximal $S^2$ number of links (lower and upper black lines, respectively).](figures/betabinmap_normal_links.png){#fig:MAPnormal}

## Conclusions

Here we derived a flexible link model for the prediction of the number of links
in ecological networks using a beta-binomial distribution for $L$ and which
outperformed the three previous and more commonly used models describing the
relationship between species richness and the number of links (the link-species
scaling, the constant connectance and the power law). More importantly, we
showed how our model's parameters non only had a clear ecological
interpretation, but how they made predictions which remained within biological
boundaries.

We believe that the appropriate modeling of the number of interactions can allow
scientists to tackle a wide variety of ecological questions, which otherwise
could have been left unexplored. For instance, the functions (productivity,
resilience, _etc._) and dynamics of ecological networks at large spatial or
temporal scales could be more easily explored. It also facilitates the
conduction of network studies in regions where interaction data is lacking,
notably due to geographical and/or financial reasons.

Our ability to model the number of links in an ecological network does not
diminish the value of data collection. Among others, data on interspecific
interactions helps understanding more deeply an ecological community and the
interactions between two or more given species, as well as making better
predictions in the statistical modelling of networks.

@GarlCald03 - small food webs and large food webs behave differently, and we
didn't really knew why before

<!-- moving this to end because I don't really know where it fits in the narrative -->

##### consequences for network topology

The constraints that we discuss in this paper have important consequences for
food web topology. For example, previous work has identified that networks with
more omnivores have longer, more linear food web structure. Our approach
recognizes that link number is also constrained, independently of structure.
This means that previous work on how network topology and link number interact
may need to be updated. @WillMart04 identified that most food webs appear to be
limited in their height, as increasingly apical species require more energy
flowing in to be sustained. In addition, constraints on omnivory appear to
"linearize" food-webs; there should therefore be a limitation on the richness of
a food web, and so for small values of $S$, the difference between assuming that
there can be between $0$, or between $(S-1)$, and $S^2$ interactions is likely
biologically relevant.


# Experimental Procedures

## Bayesian model definitions

Generative models are flexible and powerful tools for understanding and
predicting natural phenomena. These models aim to create simulated data with the
same properties as observations. Creating such a model involves two key
components: a mathematical expression which represents the ecological process
being studied, and a distribution which represents our observations of this
process. Both of these components can capture our ecological understanding of a
system, including any constraints on the quantities studied.

$$
[L | a, \phi] = \text{NegBin}(a, S, \phi)
$$



## explanation of shifted beta-binomial distribution

<!-- tk: move this 2nd order polynomail down to where it actually features in an argument

This can be re-expressed as a second order polynomial:

$$\hat L = p\times S^2 + (1-p)\times S + (p-1)\,. $${#eq:L} -->


Equation [tk eq:L] implies that $\hat L_{BB}$ has a binomial distribution, with
$S^2 - S + 1$ trials and a probability $p$ of any flexible link being realized:

$$
[L | S, p] = { S^2 - (S - 1) \choose L - (S - 1)} p^{L-(S-1)}(1-p)^{S^2 - L} ,
$$

This is often termed a _shifted Binomial distribution_.

We also note that ecological communities are different in many ways besides
their number of species ($S$). Although we assume $p$ to be fixed within one
community, the precise value of $p$ will change from one community to another.
With this assumption, our likelihood becomes a shifted beta-binomial
distribution:

$$
[L|S,\mu, \phi] =  { S^2 - (S - 1) \choose L - (S - 1)} \frac{B(L - (S - 1) + \mu \phi, S^2 - L + (1 - \mu)\phi)}{B(\mu \phi, (1 - \mu)\phi)}
$${#eq:shiftBB}

Where $B$ is the beta function.

When the number of links and number of interactions are transformed by their
natural log, $a$ and $b$ can be estimated with a linear regression, as done by
@Mart92. Here, because we want to compare all our models WAIC, we were required
to use a discrete likelihood in all cases. We fit the three models above with a
negative binomial distribution for observations. This distribution is limited to
positive integers, and can vary on both sides of the mean relationship; this
gives it a similar spirit to previous work which used a normal distribution on
log-transformed variables.

In all models we use a discrete probability distribution as the likelihood, with
the number of observed links $L_i$ above the minimum as 'successes' and the
number of possible links as 'trials'. Each model tries to capture variation in
link number greater than would be predicted by $p$ alone.

Our first model uses the beta-binomial distribution for observations of $L$;
this distribution can be parameterized in terms of its mean $\mu$ and
concentration parameter, $\phi$ :

<!-- tk a better notation here; possibly imitating H&H's style?? -->

$$\begin{aligned}
\mu &\sim  \text{Beta}(3, 7)\\
\log(\phi) & \sim \text{Normal}(3, 0.5)
\end{aligned}$${#eq:betab}

We chose our prior distribution for $p$ based on @Mart92 , who gave a value of
constant connectance equal to 0.14. While this prior is "informative", it is
weakly so; as @Mart92 did not provide an estimate of the variance for his value
we chose a relatively large variation around that mean.  However, as no
information is available to inform a prior on $\phi$, we followed the advice of
(tk Simpson et al), and performed prior predictive checks. We chose prior
parameters that generated a wide range of values for $L_i$, but which did not
frequently predict webs of either maximum or minimum connectance, neither of
which are observed in nature.

### Data used to fit all models

<!-- tk describe Mangal and its awesomeness in more wholeness -->
We evaluated our model against 255 empirical foodwebs, available in the online database `mangal.io`

We use Stan (tk doi 10.18637/jss.v076.i01) which implements Bayesian inference
using Hamiltonian Monte Carlo. We ran all models using four chains and 2000
iterations per chain. All models converged with no divergent iterations.

### Model fitting and diagnostics


All models fit without any divergent iterations, which indicates that is it safe
to make inferences about the parameter estimates and to compare the models.
However, the calculation of PSIS-LOO for the LSSL model warned of problematic
values of the Pareto-k diagnostic statistic. This indicates that the model is
heavily influenced by large values. Additionally, we had to drop the largest
observation (> 50 000 links) from all datasets in order to calculate PSIS-LOO
for the LSSL model. Taken together, this suggests that the LSSL model is
insufficiently flexible to accurately reproduce the data.  



# References
