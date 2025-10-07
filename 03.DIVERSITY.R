

### load cleaned data ###
fish_long <- readRDS("fish_long_cleaned.rds")

### Shannon diversity (hierarchical Bayesian model) ###
library(vegan)
library(brms)
library(ggplot2)
library(dplyr)
library(tidybayes)
library(loo)

### 1. Calculate Shannon diversity -------------------------------------------
shannon <- fish_long %>%
  group_by(survey_id, Species) %>%
  summarise(count = sum(Count), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Species, values_from = count, values_fill = list(count = 0)) %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(where(is.numeric)))) %>%
  ungroup()

meta <- fish_long %>%
  distinct(survey_id, Site, pair, Type, deployment_period, Researcher, Date)

shannon_data <- left_join(shannon, meta, by = "survey_id") %>%
  mutate(
    date_num = as.numeric(Date),
    deployment_period = factor(deployment_period, levels = c("Pre", "Post"))
  )

deployment_date <- as.Date("2023-09-07")

### 2. EXPLORATORY PLOTS ----------------------------------------------------

## a) Shannon by Reef Type
ggplot(shannon_data, aes(x = Type, y = shannon, color = Type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("Natural" = "#66BFA6", "Artificial" = "#007A87")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Reef Type", y = "Shannon Diversity Index",
       title = "Shannon Diversity by Reef Type") +
  theme_minimal()

## b) Shannon by Site
site_order <- c(
  "Aow Mao", "Aow Mao Wreck",
  "No Name Pinnacle", "No Name Wreck",
  "Hin Pee Wee", "Sattakut"
)
shannon_data$Site <- factor(shannon_data$Site, levels = site_order)

ggplot(shannon_data, aes(x = Site, y = shannon, color = Type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("Natural" = "#66BFA6", "Artificial" = "#007A87")) +
  labs(x = "Site", y = "Shannon Diversity Index",
       title = "Shannon Diversity by Site and Reef Type") +
  theme_minimal()

## c) Shannon by Site Pair (faceted)
ggplot(shannon_data, aes(x = Site, y = shannon, color = Type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Natural" = "#66BFA6", "Artificial" = "#007A87")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ pair, scales = "free_x") +
  labs(x = "Site", y = "Shannon Diversity Index",
       title = "Shannon Diversity by Site Within Pairs") +
  theme_minimal()

### 3. HIERARCHICAL BAYESIAN MODEL ------------------------------------------
priors_shannon <- c(
  set_prior("normal(1.8, 0.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("normal(0.3, 0.3)", class = "b", coef = "deployment_periodPost"),
  set_prior("normal(-0.3, 0.3)", class = "b", coef = "TypeNatural:deployment_periodPost"),
  set_prior("student_t(3, 0, 0.3)", class = "sds"),
  set_prior("student_t(3, 0, 0.3)", class = "sigma"),
  set_prior("student_t(3, 0, 0.5)", class = "sd")
)

fit_shannon_pair_spline <- brm(
  shannon ~ s(date_num, by = interaction(pair, Type), k = 6) +
    Type * deployment_period + (1 | Site),
  data = shannon_data,
  family = gaussian(),
  prior = priors_shannon,
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

### 4. MODEL OUTPUTS --------------------------------------------------------
summary(fit_shannon_pair_spline)
bayes_R2(fit_shannon_pair_spline)
loo(fit_shannon_pair_spline)
pp_check(fit_shannon_pair_spline)

### 5. POSTERIOR PREDICTIONS (HIERARCHICAL MODEL) --------------------------
pred_grid <- expand.grid(
  Site = unique(shannon_data$Site),
  Type = levels(shannon_data$Type),
  Date = seq(min(shannon_data$Date), max(shannon_data$Date), length.out = 200)
) %>%
  mutate(
    date_num = as.numeric(Date),
    deployment_period = factor(if_else(Date < deployment_date, "Pre", "Post"),
                               levels = c("Pre", "Post"))
  ) %>%
  left_join(shannon_data %>% distinct(Site, pair), by = "Site")

preds <- add_epred_draws(fit_shannon_pair_spline, newdata = pred_grid, re_formula = NULL)

### 6. PLOT — Bayesian predicted diversity trends --------------------------
ggplot(preds, aes(x = Date, y = .epred, color = Type, fill = Type)) +
  stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.2) +
  geom_vline(xintercept = deployment_date, linetype = "dashed", color = "gray40") +
  facet_wrap(~ pair, scales = "free_x") +
  scale_color_manual(values = c("Natural" = "#66BFA6", "Artificial" = "#007A87")) +
  scale_fill_manual(values = c("Natural" = "#66BFA6", "Artificial" = "#007A87")) +
  labs(
    x = "Date", y = "Predicted Shannon Diversity",
    color = "Reef Type", fill = "Reef Type",
    title = "Hierarchical Bayesian Shannon Diversity Trends by Site Pair"
  ) +
  theme_minimal(base_size = 12)

### 7. POSTERIOR SUMMARY — DiD CONTRAST ------------------------------------
post_did <- as_draws_df(fit_shannon_pair_spline) %>%
  select(`b_TypeNatural:deployment_periodPost`)

summary_stats <- post_did %>%
  summarise(
    mean = mean(`b_TypeNatural:deployment_periodPost`),
    median = median(`b_TypeNatural:deployment_periodPost`),
    lower95 = quantile(`b_TypeNatural:deployment_periodPost`, 0.025),
    upper95 = quantile(`b_TypeNatural:deployment_periodPost`, 0.975),
    prob_less0 = mean(`b_TypeNatural:deployment_periodPost` < 0),
    prob_greater0 = mean(`b_TypeNatural:deployment_periodPost` > 0)
  )
print(summary_stats)

ggplot(post_did, aes(x = `b_TypeNatural:deployment_periodPost`)) +
  geom_density(fill = "#66BFA6", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(x = "DiD Contrast (Natural vs Artificial, Post–Pre)",
       y = "Posterior density",
       title = "Posterior distribution of the DiD effect") +
  theme_minimal(base_size = 12)
