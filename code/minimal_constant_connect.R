library(brms)
library(tidyverse)
library(tidybayes)
library(ids)

webs <- read_csv("data/ls.csv")

foodwebs <- webs %>%
  filter(P > 0) %>%
  mutate(R = L - (S - 1),
         n = S^2 - (S - 1),
         id = as.character(id))

## alternative way to draw the bounds onto ggplot2 figures; such that the lines look good for small x

boundpoints <- tibble(S = seq(5, 700, by = 1),
                      b1 = S - 1,
                      b2 = S ^ 2) %>%
  pivot_longer(starts_with("b"), values_to = "y")


# constant connectance model ----------------------------------------------

# Martinez 92 has connectance = 0.14
logit_scaled(0.14, 0, 1) # -1.815

R_bf <- bf(R|trials(n) ~ 1 + (1|id), family = binomial(link = "logit"))

get_prior(R_bf, data = foodwebs, family = binomial(link = "logit"))

R_prior <- c(prior(normal(-1.815, 0.3), class = "Intercept"),
  prior(exponential(3),      class = "sd"))

prior_predictive <- brm(
  R_bf,
  prior = R_prior,
  data = foodwebs,
  cores = 1, chains = 1,
  sample_prior = "only")

foodwebs %>%
  add_predicted_draws(prior_predictive, n = 12) %>%
  ggplot(aes(x = S, y = .prediction + (S - 1))) +
  geom_point() +
  scale_x_log10() +
  geom_line(aes(y = y, group = name), data = boundpoints) +
  coord_trans(x = "log10", y = "log10") +
  facet_wrap(~.draw)

constant_binomial <- brm(
  R_bf,
  prior = R_prior,
  data = foodwebs,
  file = "models/constant_binomial",
  cores = 4, chains = 4,
  sample_prior = "yes")
beepr::beep(7)
constant_binomial


foodwebs %>%
  add_predicted_draws(constant_binomial, n = 100) %>%
  ggplot(aes(x = S, y = .prediction + (S - 1))) +
  geom_point(alpha = 0.12) +
  geom_point(aes(y = L), colour = "orange", data = foodwebs) +
  geom_line(aes(y = y, group = name), data = boundpoints) +
  coord_trans(x = "log10", y = "log10") +
  theme_bw()
ggsave("figures/posterior_dots.png")

# and predictions for new observations:
constant_binomial_predictions <- tibble(
  S = seq(from = min(foodwebs$S), to = max(foodwebs$S), by = 1),
  n = S^2 - (S - 1),
  id = ids::proquint(length(S))
) %>%
  add_predicted_draws(constant_binomial, n = 500, allow_new_levels = TRUE)

# plot of line marginalized over intercepts
constant_binomial_predictions %>%
  mutate(.prediction = .prediction + (S - 1)) %>%
  tidybayes::median_qi(.width = c(0.97, 0.73, 0.67)) %>%
  ggplot(aes(x = S, y = .prediction)) +
  geom_lineribbon(color = "grey", lwd = .3) +
  scale_fill_brewer(palette = "Oranges", direction = -1) +
  geom_point(aes(y = L), colour = "black", data = foodwebs, alpha = 0.3) +
  geom_line(aes(y = y, group = name), data = boundpoints) +
  coord_trans(x = "log10", y = "log10") +
  guides(fill = FALSE) +
  theme_bw()
ggsave("figures/posterior_density_interval_line.png")


# ## as a cloud of points
# constant_binomial_predictions %>%
#   mutate(.prediction = .prediction + (S - 1)) %>%
#   ggplot(aes(x = S, y = .prediction)) +
#   geom_point(color = "#E6550D", alpha = .1) +
#   geom_point(aes(y = L), colour = "black", data = foodwebs, alpha = 0.3, size = 2) +
#   geom_line(aes(y = y, group = name), data = boundpoints) +
#   coord_trans(x = "log10", y = "log10") +
#   theme_bw()
#


# connectance figure --------------------------------------------------

intercepts <- spread_draws(constant_binomial, b_Intercept)

S <- c(seq(1, 1.9, by = 0.1), 2:800)

S_and_b <- intercepts %>%
  select(.draw, b_Intercept) %>%
  sample_n(50) %>%
  expand_grid(S)

avg_network <- S_and_b %>%
  mutate(p = map_dbl(b_Intercept, ~ exp(.)/(1 + exp(.))),
         L = p * S^2 + S* (1 - p) + p - 1,
         Co = L / S ^2)

## fitted average values of connectivity
avg_network %>%
  ggplot(aes(x = S, y = Co, group = .draw)) +
  geom_line(alpha = 0.2, colour = "darkgreen") +
  scale_x_log10() +
  geom_point(aes(x = S, y = Co), inherit.aes = FALSE,
             data = foodwebs %>% mutate(Co = L / S^2)) +
  theme_bw() +
  stat_function(fun = function(x) (x - 1)/x^2, lty = 2, colour = "grey")

ggsave("figures/average_connectivity.png")


### marginalize over random effects

intercepts <- spread_draws(constant_binomial, b_Intercept, sd_id__Intercept)

b_and_sd <- intercepts %>%
  select(.draw, b_Intercept, sd_id__Intercept) %>%
  sample_n(100) %>%
  mutate(offset = rnorm(length(sd_id__Intercept), 0, sd_id__Intercept),
         p = map2_dbl(b_Intercept,offset, ~ exp(.x + .y)/(1 + exp(.x + .y))))

variation_p <-  b_and_sd %>%
  expand_grid(S)

network_variation <- variation_p %>%
  mutate(
    L = p * S^2 + S* (1 - p) + p - 1,
    Co = L / S ^2)

## as points -- not that interesting actually
network_variation  %>%
  ggplot(aes(x = S, y = Co, group = .draw)) +
  geom_point(alpha = 0.1, color = "#E6550D") +
  scale_x_log10() +
  geom_point(aes(x = S, y = Co),
             inherit.aes = FALSE,
             data = foodwebs %>% mutate(Co = L / S^2),
             col = "black", alpha = 0.7) +
  theme_bw() +
  stat_function(fun = function(x) (x - 1)/x^2)

network_var_summary <- network_variation %>%
  select(S, Co) %>%
  group_by(S) %>%
  point_interval(.width = c(.73, .89, .97))

network_var_summary %>%
  # select(-Co, -.point) %>%
  ggplot(aes(x = S, y = Co, ymin = .lower,
             ymax = .upper,
             group = .width)) +
  geom_lineribbon(colour = "darkgrey", lwd = 0.2) +
  scale_x_log10() +
  geom_line(aes(x = S, y = Co, group = .draw),
            alpha = 0.3,
            inherit.aes = FALSE,
            data = avg_network) +
  geom_point(aes(x = S, y = Co),
             inherit.aes = FALSE,
             data = foodwebs %>% mutate(Co = L / S^2),
             pch = 21, fill = "darkorange", size = 2) +
  theme_bw() +
  stat_function(fun = function(x) (x - 1)/x^2) +
  scale_fill_brewer(palette = "Greens")

ggsave("figures/marginal_connectance.png")



# posterior predictions of Co, L, Ld --------------------------------------

posterior_links_constant <- foodwebs %>%
  add_predicted_draws(constant_binomial, n = 100) %>%
  mutate(Lhat = .prediction + (S - 1))

posterior_links_constant %>%
  ggplot(aes(x = S, y = Lhat - L)) + geom_point()



posterior_links_constant %>%
  ggplot(aes(x = S, y = Lhat/S^2 - L/S^2)) + geom_point()

posterior_links_constant %>%
  mutate(Co_difference = Lhat/S^2 - L/S^2) %>%
  select(-.prediction, -Lhat) %>%
  point_interval() %>%
  ggplot(aes(x = S, y = Co_difference)) +
  geom_pointinterval() +
  scale_x_log10()



# z-scores ----------------------------------------------------------------

z_score <- function(L, s, p){
  ((L - (s - 1)) - p * (s^2 - (s - 1))) / sqrt((1 - p) * p * (s^2 - s + 1))
}

z_low <- function(p, S) {
  - sqrt(
    (p * (S^2 - S + 1) )/ (1 - p))
}

z_up <- function(p, S) {
  sqrt((1 - p) * (S^2 - S + 1) / p)
}


z_low(0.07, seq(1,300, length.out = 5))

foodwebs %>%
  mutate(z = map2_dbl(L, S, ~ z_score(L = .x, s = .y, p = 0.07))) %>%
  ggplot(aes(x = S, y = z)) +
  geom_point(alpha = 0.4) +
  # geom_line(y = -1/(1-0.07)) +
  # geom_line(y = 1/0.07) +
  scale_x_log10() +
  coord_cartesian(ylim = c(-50,150)) +
  theme_bw() +
  stat_function(fun = function(x) z_up(0.07, x)) +
  stat_function(fun = function(x) z_low(0.07, x))

ggsave("figures/zscores.png")

