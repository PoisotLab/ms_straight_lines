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
suggested that describing food webs was a necessity [@McCa07; @SeibCado18]. To
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
are modeled as the number of species times a constant:

$$ L_{\textsc{lssl}} = b\times S \,$$ {#eq:lssl}

with $b \approx 2$. This model imagines that every species added to a community
increases the number of interactions by two -- for example, an animal which
consumes one resource and is consumed by one predator. Yet this model started to
show its deficiencies when data on larger food webs became available, which
revealed that $L$ increases faster than a linear function of $S$ would. In
response to the idea that the *average degree* is constant, @Mart92
suggested instead that *connectance* is unchanged in response to $S$; in other
words, a food web is always equally filled, regardless of whether it has 5 or
5000 species. Under the so-called "constant connectance" model, the number of
links is proportional to the richness squared,

$$ L_{\textsc{cc}} = b\times S^2\,, $$ {#eq:cc}

where $b$ is a constant in $]0,1[$ representing the expected value of
connectance. This model can be relaxed by assuming that the scaling of $L$ with
$S$ does not necessarily follows the maximum number of interactions, and the
best fit was with a model of the form

$$ L_{\textsc{pl}} = b\times S^a\,, $$ {#eq:pl}

which is trivially a linear relationship between $\text{log}(L)$ and
$\text{log}(S)$. This power law model can be parameterized in arbitrarily
complex ways, including spatial scaling and species area relationships
[@BrosOstl04]; it should further be noted that this model is a synthesis of
preview hypotheses, encompassing both the link-species scaling ($a=1, b\approx
2$) and the strict constant connectance ($a=2, 0<b<1$) depending on which parameters
are fixed. Power laws are very flexible, and indeed this function matches
empirical data well -- so well that it is often treated as a "true" model which
captures the scaling of link number with species richness [@WinePian01;
@GarlCald03; @RiedRall10], and from which we should draw ecological inferences
about what shapes food webs. However, this approach is limited, because the
parameters of a power law relationship can arise from many mechanisms, and are
difficult to reason about ecologically.

But the question of how informative parameters of a power law can be is moot.
Indeed, both the general model and its variants share an important shortcoming:
they cannot be used for prediction while remaining within the bounds set by
ecological principles. In short, while they can describe the *data* adequately,
they are fundamentally unable to represent the mechanisms through which these
data emerged. This has two causes. First, models that are variations of $L
\approx b\times S^a$ have no constraints --  with the exception of the "constant
connectance" model, in which $L_{\textsc{cc}}$ has a maximum value of $S^2$.
However, we know that the number of interactions within a food web is both lower
and upper bounded [@Mart92; @PoisGrav14]: there can be no more than $S^2$ links,
and there can be no fewer than $S-1$ links. This minimum of $S-1$ represents
communities where all species interact and at least some of the organisms are
heterotrophs [@Mart92]. Numerous simple foodwebs could have this minimal number
of links -- for example, a linear food chain wherein each trophic level consists
of a single species, each of which consumes only the species below it; or a
grazing herbivore which feeds on every plant in a field. Thus interaction number
is constrained by ecological principles to be between $S^2$ and $S-1$, something
which no present model includes. Secondly, accurate predictions of $L$ from $S$
are often difficult because of how parameters are estimated. This is usually
done using a Gaussian likelihood for $L$, often after log transformation of both
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

Based on the ecological constraints discussed earlier, we know that the number
of links $L$ is an integer such that $S-1 \le L \le S^2$. Because we know that
there are at least $S-1$ interactions, there can be at most $S^2-(S-1)$ links
*in excess* of this quantity. The $S-1$ minimum links do not need to be
modelled, because their existence is guaranteed as a pre-condition of observing
the network. The question our model should address is therefore,how many of
these $S^2-(S-1)$ "flexible" links are actually present? A second key piece of
information is that the presence of an interaction can be viewed as a discrete
stochastic event, with two outcomes (there is, or there is not, an interaction).
Assuming that all of these flexible links have the same chance of being
realized, which we call $p$, then we can write the expected number of links as

$$ L_{\textsc{fl}} = p\times\left[S^2-(S-1)\right]+(S-1)\,, $$ {#eq:lhat}

where $p \in [0,1]$. When $p = 1$, $L$ is at its maximum ($S^2$), and when $p =
0$ it is at the minimum value ($S - 1$). We use the notation $L_{\textsc{fl}}$ to
represent that our model considers the number of "flexible" links in a food web;
that is, the number of links in excess of the minimum but below the maximum.

Because we assume that every one of these flexible links is an independent
stochastic event with only two outcomes, we can follow recent literature on
probabilistic ecological networks [@PoisCirt16] and represent every flexible
link as an independent Bernoulli trial with a probability $p$ of existing.
Furthermore, the observation of $L_i$ links in community $i$ represents an
aggregation of $S^2 - (S - 1)$ such trials. If we then assume that $p$ is a
constant for all links in the same ecological community, but may vary between
communities, we can model the distribution of links directly as a shifted
beta-binomial variable.

$$ [L|S,\mu, \phi] =  { S^2 - (S - 1) \choose L - (S - 1)} \frac{B(L - (S - 1) + \mu \phi, S^2 - L + (1 - \mu)\phi)}{B(\mu \phi, (1 - \mu)\phi)} $$ {#eq:shiftBB}

Where $\mathrm{B}$ is the Beta function, $\mu$ is the average probability of a
flexible link being realized (*i.e.* the average value of $p$ across communities) and $\phi$ is the
concentration around this value. The support of this distribution is limited to
only ecologically realistic values of $L$: the number of species determines the
number of possible links and it is shifted to the right by $S-1$. This means
that the problem of estimating values for $\mu$ and $\phi$ is reduced to fitting
the univariate distribution described in +@eq:shiftBB. For more detailed
explanation of the model derivation, fitting, and comparison, see Experimental
Procedures.

In this paper we will compare our flexible links model to three previous models
for $L$. We estimate parameters and compare the performance of all models using
open data from the `mangal.io` networks database [@PoisBais16]. Finally, we will
show how our model for $L_{\textsc{fl}}$ suggests a new and more useful way of
thinking about metrics of network structure and discuss how generative models
can be useful tools for including our knowledge of a system into our
predictions.

# Results and Discussion

## Flexible links model fits better and makes a plausible range of predictions

| Model                              | eq.      | PSIS-LOO         | $\Delta \text{ELPD}$ | $SE_{\Delta \text{ELPD}}$ |
| ---------------------------------- | -------- | ---------------- | -------------------- | ------------------------- |
| Flexible links                     | @eq:lhat | 2520.5 ± 44.4    | 0                    | 0                         |
| Power law [@BrosOstl04]            | @eq:pl   | 2564.3 ± 46.6    | -21.9                | 6.5                       |
| Constant [@Mart92]                 | @eq:cc   | 2811.0 ± 68.3    | -145.3               | 21.1                      |
| Link-species scaling [@CoheBria84] | @eq:lssl | 39840.1 ± 2795.1 | -18659.8             | 1381.7                    |

Table: Comparison of the different models. Pareto-smoothed important sampling values and differences relative to the maximum in the expected log predictive density for the flexible links and the three competing models. The mean and standard deviation (SD) (standard error (SE)) is given for the two metrics. {#tbl:comparison}

All models fitted well, without any problematic warnings from Stan's diagnostics
(see Experimental Procedures), but our model for flexible links outperformed
previous solutions to the problem of modelling $L$. The flexible link model,
which we fit via a beta-binomial observation model, had the most favourable
values of PSIS-LOO information criterion (+@tbl:comparison) and of expected log
predictive density (ELPD), relative to the three competing models which used a
negative binomial observation model. Pareto-smoothed important sampling serves
as a guide to model selection [@VehtGelm17]; like other information criteria it approximates
the error in cross-validation predictions. Smaller values indicate a model which
makes better predictions. The calculation of PSIS-LOO can also provide some
clues about potential model fits; in our case the algorithm suggested that the
constant connectance model was sensitive to extreme observations. Information
criteria are only a rough guide to model selection; as always domain expertise
should take precedence. The expected log predictive density (ELPD), on the other
hand, measures the predictive performance of the model; here, higher values
indicate more reliable predictions [@VehtGelm17]. This suggests that the flexible link model
will make the best predictions of $L$.

Useful predictions for $L$ must however stay within realistic boundaries
determined by ecological principles. We generated posterior predictions for all
models and visualized them against these constraints (+@fig:PP_counterfactual).
The LSSL model clearly underestimated the number of links, especially in large
networks: its predictions were frequently lower than the minimum $S-1$. The
constant connectance and power law models also made many predictions below this
value, especially for small values of $S$. The flexible link model made roughly
the same predictions, but within ecologically possible values.

![**The flexible link model fits better and makes a plausible range of predictions.** The number of links is plotted as a function of species richness obtained from the posterior distributions of A) the link-species scaling, B) the constant connectance, C) the power law and D) the flexible link models. In each panel, the colored line represent the median predicted link number and the grey areas cover the 78% and 97% percentile intervals. Empirical data from the `mangal.io` database are plotted in each panel (grey dots), as well as the minimal $S-1$ and maximal $S^2$ number of links (lower and upper black lines, respectively).  ](figures/models_links.png){#fig:PP_counterfactual}

## The flexible links model makes realistic predictions for small communities

Constraints on food web structure are especially important for small
communities. This is emphasized in +@fig:real_predict, which shows that all
models other than the flexible links model fail to stay within realistic
ecological constraints. The link-species scaling model made around 29% of
unrealistic predictions of link numbers for every value of $S$ ($3 \leq S \leq
750$). The constant connectance and power law models, on the other hand, also
produced unrealistic results but for small networks only: more than 20% were
unrealistic for networks comprising less than 12 and 7 species, respectively.
Only the flexible links model, by design, never failed to predict numbers of
links between $S-1$ and $S^2$.  

![**Only the flexible link model makes realistic predictions for small communities.** Here we show the proportion of posterior predictions from each of our 4 models which fall outside ecologically realistic values. The proportion of predictions in the correct range increases with species richness for the constant connectance and power law models. Shaded area shows the 5%, 50% and 95% quantiles of the distribution of S, demonstrating that many communities have potentially incorrect predictions under previous models](figures/real_predict.png){#fig:real_predict}

## Parameter estimates for all models

| Model                      | parameter | interpretation                      | value | SD     |
| -------------------------- | --------- | ----------------------------------- | ----- | ------ |
| $bS$ [@CoheBria84]         | $b$       | links per species                   | 2.2   | 0.047  |
|                            | $\kappa$  | concentration of $L$ around mean    | 1.4   | 0.12   |
| $bS^2$ [@Mart92]           | $b$       | proportion of links realized        | 0.12  | 0.0041 |
|                            | $\kappa$  | concentration of $L$ around mean    | 4.0   | 0.37   |
| $bS^a$ [@BrosOstl04]       | $b$       | proportion of relationship          | 0.37  | 0.054  |
|                            | $a$       | scaling of relationship             | 1.7   | 0.043  |
|                            | $\kappa$  | concentration of $L$ around mean    | 4.8   | 0.41   |
| $(S^2 - (S - 1))p + S - 1$ | $\mu$     | average value of $p$                | 0.086 | 0.0037 |
|                            | $\phi$    | concentration around value of $\mu$ | 24.3  | 2.4    |

Table: Parameter estimates for all models. Mean and standard deviation (SD) is given for each parameter. {#tbl:parameters}

Although we did not use the same approach to parameter estimation as previous
authors, our approach to fitting these models recovered parameter estimates that
are broadly congruent with previous work. We found very consistent values of 2.2
for $b$ of the LSSL model, which is close to the original value of approximately
2 [@CoheBria84]. Similarly, we found a value of 0.12 for $b$ of the constant
connectance model, which was consistent with original estimate of 0.14
[@Mart92]. Finally, the parameters values we found for the power law were also
comparable to earlier estimates [@BrosOstl04].

For large communities, our model should behave similarly to the constant
connectance model; it is no surprise then than $\mu$ was about 0.09, which is
close to our value of 0.12 for constant connectance. In addition, we obtained a
rather large value of 24.3 for $\phi$, which shrinks the variance around the
mean of $p$ to approximately 0.003 ($var(p)=\mu(1-\mu)/(1+\phi)$). This
indicates that food webs are largely constant in their probability of flexible links being realized.

## Connectance and linkage density can be derived from a model for links

Of the three important quantities which describe networks ($L$, $Co$ and $L_D$)
we have directly modelled $L$ only. However, we can reuse the posterior for our
model of $L$ to parameterize a distribution for connectance ($L/S^2$) and
linkage density ($L/S$). We can derive this by noticing that +@eq:lhat can be
rearranged to show how $Co$ and $L_D$ are linear transformations of $p$:

$$ Co = \frac{L}{S^2} = p\left(1 - \frac{S-1}{S^2}\right) + \frac{S-1}{S^2} ,$$ {#eq:co}

and

$$ L_D = \frac{L}{S} = p \left(S - \frac{S-1}{S} \right) +  \frac{S-1}{S},$$ {#eq:ld}

In a beta-binomial distribution, it is assumed that the probability of success
$p$ varies among groups of trials according to a Beta distribution. Since $p$
has a Beta distribution the linear transformations described by +@eq:co and
+@eq:ld also describe Beta distributions which have been shifted and scaled
according to the number of species $S$ in a community.

It is worth noting that +@eq:lhat can be expressed as a second degree
polynomial, ($p\times S^2  + (1-p)\times S + (p-1)$), whose leading term is
$p\times S^2$. Therefore, for large ecological systems, where $S$ has a high
value, +@eq:co and +@eq:ld respectively approach $Co = L/S^2 \approx p$ and
$L_D = L/S \approx pS$. These are notable properties, as they imply that our
model captures both the behavior of +@eq:lsll and of +eq:cc, while having a
markedly better fit (+@tab:comparison).

This implies that the addition of $n$ species should increase the linkage
density by approximately $p\times n$. For example, the addition of 11 new
species ($p^{-1}$ according to +@tab:parameters) should increase the linkage
density in the food web by roughly 1, meaning that each species in the original
network would be expected to develop 2 additional interactions. For large enough
ecological networks, we should expect a connectance which is independent of $S$.
Thus $p$ has an interesting ecological interpretation: it represents the average
connectance of networks large enough that the proportion $(S-1)/S^{2}$ is
negligible.

Just as $L$ must be within ecologically meaningful bounds, $Co$ (+@eq:co) and
$L_D$ (+@eq:ld) must be as well. The connectance of a food web is bounded by 0
and 1. However, the minimum bound on links similarly imposes a lower value on
connectance. This means that the distribution for $Co$ will be a shifted beta
distribution, a transformed version of the distribution for $p$.

We can convert the distribution for $p$ into one for $Co$ by replacing $p$ with
the transformation of $Co$ as described above (+@eq:co), and rescaling by the
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
associated with some of the most commonly used network metrics [@PoisGrav14;
@Chag15], and contains meaningful information on the stability [@DunnWill02a;
@MontPimm06] and dynamics [@VieiAlme15] of ecological communities. A probability
distribution for connectance non only accounts for the variability between
networks, but can be used to describe fundamental properties of food webs and to
identify ecological and evolutionary mechanisms shaping communities.

## Normal approximation provides an analytic z-score

Ecologists are often faced with the issue of comparing several networks. Often,
they wish to know if the network they have is "unusual" relative to some
expectation. Traditionally these comparisons have been done by constructing a
null distribution, which involves simulating random matrices [@FortBasc06;
@BascJord03]; importantly, this approach assumes that (i) connectance is a fixed
property of the network, which does not involve any stochasticity, and (ii) the
simulated network distribution is an accurate and unbiased description of the
null distribution. Yet recent advances in the study of probabilistic ecological
networks show that connectance can also be modelled as a variable quantity
[@PoisCirt16]; given that connectance drives most of the values of food web
structure [@PoisGrav14], we provide a way to assess whether the number of links
in a network (and therefore its connectance) is surprising. We do so using maths
rather than simulations.

The shifted beta-binomial can be approximated by a normal distribution:

$$L \sim Normal(\bar{L}, \sigma_L^2)$$

$$\bar{L} = (S^2 - S + 1) \mu + S - 1$$

$$\sigma_L^2 = (S^2 - S + 1) \mu (1 - \mu)(1 + \frac{S(S-1)}{\phi + 1})$$

This means that given a network with observed species richness $S_{obs}$ and
observed links $L_{obs}$, we can calculate its $z$-score as

$$z = \frac{L_{obs} - \bar{L}}{\sqrt{\sigma_L^2}} \,.$$ {#eq:z}

A network where $L = \bar{L}$ will have a $z$-score of 0, and any network with
more (fewer) interactions will have a positive (negative) $z$-score. We suggest
that the use of a $z$-score could help identify significantly under (over)
sampled networks and estimate their number of missing (extra) links.

In (+@fig:MAPnormal), we show that the predictions made by the normal
approximation (panel B) are similar to those made by the beta distribution
parameterized with the maximum _a posteriori_ values of $\mu$ and $\phi$ (panel
A), although the later can undershoot the constraint on the minimum number of
links.

![**The shifted beta-binomial distribution can be approximated by a normal distribution.** The number of links is plotted as a function of species richness obtained from A) the maximum _a posteriori_ estimates of the flexible link model and B) its normal approximation. In each panel, the colored line represent the median predicted link number and the grey areas cover the 78% and 97% percentile intervals.](figures/betabinmap_normal_links.png){#fig:MAPnormal}

## Conclusions

Here we derived a flexible link model for the prediction of the number of links
in ecological networks using a beta-binomial distribution for $L$, and show how
it outperforms previous and more commonly used models describing this
relationship. More importantly, we showed how our model's parameters not only
have a clear ecological interpretation (specifically, the value of $p$ in
+@eq:lhat is the expected value of the connectance when $S$ is large), but how
they made predictions which remained within biological boundaries.

This model also casts new light on previous results on the structure of food
webs: small and large food webs behave differently [@GarlCald03]. Specifically,
ecological network most strongly deviate from both scale free and small world
expectations when connectance is high [@DunnWill02a]. In our model, this
behaviour is a natural prediction, as the connectance increases sharply for low
species richness (+@fig:beta_distributions), as the additive term $(S-1)S^{-2}$
in +@eq:co becomes progressively larger. In a sense, small ecological networks
are different only because due to the low values of $S$, there are only a very
limited number of flexible links, and this drives connectance to be larger.
Connectance in turn has inmplications for many ecological properties. A recent
research direction has been to reveal its impact on resistance to invasion:
denser networks with a higher connectance are comparatively more difficult to
invade [@SmitMoor16]; different levels of connectance are also associated with
different combination of primary producers, consumers, and apex predators
[@WillMart00; @WillMart04], which in turns means that different species will
have an easier success invading differently connected networks [@BaisRuss10].
Because we can infer connectance from the richness of a community, our model
also ties the invasion resistance of a network to its species richness.

Yet our model introduces a puzzling question. According to +@eq:ld, at large
values of $S$, the linkage density scales according to $p\times S$ (which is
supported by empirical data), and so species are expected to have on average
$2\times p\times S$ interactions. A useful concept in evolutionary biology is
the "Darwinian demon" [@Law79], *i.e.* an organism that would have infinite
fitness in infinite environments. Our model seems to predict the emergence of
what we call Eltonian demons, which can have arbitrarily large number of
interactions. Yet we know that constraints on handling time of preys, for
example, imposes hard limits on diet breadth [@ForiJenk17]. This result suggests
that there are other limitations to the size of food webs; indeed, the fact that
$L/S$ increases to worryingly large values only matters if ecological processes
allow $S$ to be large enough. It is known that food webs can get as high as
energy transfer allows [@ThomBros12], and as wide as competition allows
[@KefiBerl12]. In short, and as +@fig:real_predict suggests, since food webs are
likely to be constrained to remain within an acceptable richness, we have no
reason to anticipate that $p\times S$ will keep growing infinitely.

![**Stability imposes a limit on network growth**. Using +@eq:ld, we can calculate the maximum standard deviation in the strength of interactions which should ensure food web stability, $\sigma^\star = \sqrt{L_D}^{-1}$. This value falls sharply when the number of species increases, which will limit the stability of large food webs, and therefore explain why Eltonian demons should not emerge.](figures/may.png){#fig:stability}

In fact, May [@May72] suggested that a network of richness $S$ and connectance
$Co$ is stable as long as the criteria $\sigma \sqrt{S/Co} < 1$ is satisfied,
with $\sigma$ being the standard deviation of the strengths of interactions.
Under our model, $Co$ is derived from $S$, and $S/Co$ is the linkage density as
per +@eq:ld. Although this criteria is not necessarily stringent enough for the
stability of food webs [@AlleTang12; @AlleTang15], it allows deriving an
approximate value $\sigma^\star$ which is the value of $\sigma$ above which the
previous criteria is not satisfied, and the system is expected to be unstable.
This threshold is the solution to $\sigma^\star = \sqrt{L_d}^{-1}$, where $L_D$
is defined as in +@eq:ld. We illustrate this result in +@fig:stability, which
reveals that $\sigma^\star$ falls to 0 for larger species richness. These
results explain how ecological limitations (here on stability) can limit the
size of food webs, are in agreement with previous simulations, placing the
threshold for stability at about 1200 species in food webs [@AlleTang12].

Finally, our results bear important consequences for the nascent field of
studying network-areas relationships (NAR). As it has long been observed that
not all species in a food web diffuse equally through space [@HoltLawt99],
understanding how the shape of networks varies when the area increases is an
important goal, and in fact underpins the development of a macroecological
theory of food webs [@BaisGrav19]. Using a power-law as the acceptable
relationship between species and area [@Deng09; @WillGast01], the core idea of
studying NAR is to predict network structure as a consequence of the effect of
spatial scale on species richness [@GaliLurg18]. Drawing on these results, we
provide in +@fig:nar a simple illustration of the fact that, due to the
dispersal of values of $L$, the relationship between $L/S$ and area can have a
really wide confidence interval. While our simulations generally match the
empirical results on this topic [*e.g.* @WoodRuss15], they suggest that we will
observe many relationships between network structure and space ,and that picking
out the signal of network area relationships might be difficult.

![**Networks in space**. Representing the species richness as $S = k\times A^z$, with $A = 200$ and $k = 0.27$ [@GaliLurg18], we can compare the predictions of our model to that of the generally accepted power law (+@eq:pl). While our model predicts a larger linkage density in larger areas, the confidence intervals around this prediction are extremely large. - scales faster than the power law, but the confidence interval is extremely high, suggesting that we may observe either very weak, or very strong, effects of area increase.](figures/nar.png){#fig:nar}

As a conclusion, we would like to note that the relationship between $L$ and $S$
has been underpinning most of the literature on food web structure since the
1980s. Additional generations of data allows to switch from the link-species
scaling law, to constant connectance, to more general formulations based on a
power law. Our model breaks with this tradition of iterating over the same
family of relationship, and instead draws from our knowledge of ecological
processes, and from novel tools in probabilistic programming. As a result, we
provide predictions of the number of interactions which are closer to empirical
data, allow to derive new ecological insights, and can be safely assumed to
always fall within realistic values. The results presented in +@fig:stability
and +@fig:nar seem largely confirmatory, and the ability of our model to reach
the same conclusions is also a confirmation of its validity; we would like to
point out that these approaches would usually require to make inferences on the
parameters of interests, but also on the properties of a network for a given
species richness. Therefore, our model allows a real economy of parameters, as
it offers the ability to get several key elements of network structure for free
if only the species richness is known.

# Experimental Procedures

## Availability of code and data

All code and data to reproduce this article is available at `ZENODO REPO TK`.

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

This is often termed a _shifted binomial distribution_.

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

## Model fitting - data and software

<!-- tk describe Mangal and its awesomeness in more wholeness -->
We evaluated our model against 255 empirical foodwebs, available in the online
database `mangal.io` [@PoisBais16]. We queried metadata (number of nodes and
number of links) for all networks, and considered as food webs all networks
having interactions of predation and herbivory. We use Stan [@CarpGelm17a] which
implements Bayesian inference using Hamiltonian Monte Carlo. We ran all models
using four chains and 2000 iterations per chain. All models converged with no
divergent iterations.

## Model fitting - diagnostics

All models fit without any divergent iterations, which indicates that is it safe
to make inferences about the parameter estimates and to compare the models.
However, the calculation of PSIS-LOO for the LSSL model warned of problematic
values of the Pareto-k diagnostic statistic. This indicates that the model is
heavily influenced by large values. Additionally, we had to drop the largest
observation (> 50 000 links) from all datasets in order to calculate PSIS-LOO
for the LSSL model. Taken together, this suggests that the LSSL model is
insufficiently flexible to accurately reproduce the data.  



# References
