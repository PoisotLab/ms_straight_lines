### Fig. 1

**The flexible links model fits better and makes a plausible range of
predictions.** The number of links is plotted as a function of species richness
obtained from the posterior distributions of A) the link-species scaling, B) the
constant connectance, C) the power law and D) the flexible links models. In each
panel, the colored line represents the median predicted link number and the grey
areas cover the 78% and 97% percentile intervals. Empirical data from the
`mangal.io` database are plotted in each panel (grey dots), as well as the
minimal $S-1$ and maximal $S^2$ number of links (thinner and bolder black lines,
respectively). Predictions from the flexible links model are always plausible:
they stay within these biological
boundaries.

### Fig. 2

**Only the flexible links model makes realistic predictions for small
communities.** Here we show the proportion of posterior predictions from each of
our 4 models which fall outside ecologically realistic values. The proportion of
predictions in the correct range increases with species richness for the
constant connectance and power law models. Shaded area shows the 5%, 50% and 95%
quantiles of the distribution of $S$, demonstrating that many communities have
potentially incorrect predictions under previous models.

### Fig. 3

**Connectance and linkage density can be derived from a model for links.** A)
Connectance and B) linkage density are plotted as a function of species
richness, for the maximum _a posteriori_ estimates of the flexible links model.
In each panel, the colored line represents the median predicted quantity and the
grey areas cover the 78% and 97% percentile intervals. Empirical data from the
`mangal.io` database are plotted in each panel (grey dots). In A), the minimal
$(S-1)/S^2$ connectance and in B) the minimal $(S-1)/S$ and maximum $S$ linkage
density are plotted (black lines).

### Fig. 4

**Empirical distribution of food web $z$-scores** The $z$-scores of all food webs
archived on `mangal.io` have been computed using +@eq:z. Food webs with an
absolute $z$-score above 1.96 are in pink. The shaded region comprises all food
webs with an absolute $z$-score below 1.96.

### Fig. 5

**The shifted beta-binomial distribution can be approximated by a normal
distribution.** The number of links is plotted as a function of species richness
obtained from A) the maximum _a posteriori_ estimates of the flexible links
model and B) its normal approximation. In each panel, the colored line represents
the median predicted link number and the grey areas cover the 78% and 97%
percentile intervals. The minimal $S-1$ and maximal $S^2$ numbers of links are
plotted in each panel (thinner and bolder black lines,
respectively).

### Fig. 6

**Many different network-area relationships are supported by the data**.
Representing the species richness as $S = k\times A^z$ (panel A), with $A$ being
the relative area size, $k = 200$ being the maximal species richness, and $z =
0.27$ a scaling exponent [@GaliLurg18]. We then use the posterior distribution
of $L$ to predict how $L_D$ should scale with $A$. We compare the predictions of
our model to that of the generally accepted power law (+@eq:pl). While our model
predicts a larger linkage density in larger areas (panel B), the confidence
intervals around this prediction (grey areas covering the 78% and 97% percentile
intervals) are extremely large. In particular, our model scales faster than the
power law, but the confidence interval is high (due to the scaling of variance
with $S$, +@eq:bb_sigma). This suggests that we may observe either very weak, or
very strong, effects of area on networks.

### Fig. 7

**Stability imposes a limit on network size**. Using +@eq:ld, we can calculate
the maximum standard deviation in the strength of interactions which should
ensure food web stability, $\sigma^\star = 1/\sqrt{L_D}$ (panel A). The colored
line represent the median value of maximum standard deviation, based on the
posterior distribution of the flexible links model, and the grey areas cover the
78% and 97% percentile intervals. The fine and dark lines indicate the maximum
and minimum values of maximum standard deviation, respectively. The dotted line
shows the maximum for the average $L_D$, as given by +@eq:ld. The maximum
standard deviation falls sharply when the number of species increases, which
will limit the stability of large food webs, and therefore explain why Eltonian
demons should not emerge. In panel B, we show the probability of a network with
$S$ species being stable, based on draws from the posterior distribution, for
$10 \le S \le 1000$ - larger networks (thicker lines) are increasingly unlikely
to be stable.
