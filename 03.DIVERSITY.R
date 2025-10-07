
library(vegan)
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidybayes)
library(loo)

# assumes 00.SET defined output_dir and brm_m2()
fish_long <- readRDS("fish_long_cleaned.rds")
deployment_date <- as.Date("2023-09-07")
reef_cols <- c("Natural" = "#66BFA6", "Artificial" = "#007A87")

# 1) Shannon diversity per survey -------------------------------------------
shannon <- fish_long %>%
  group_by(survey_id, Species) %>%
  summarise(count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = count, values_fill = 0) %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(where(is.numeric)))) %>%
  ungroup()

meta <- fish_long %>%
  distinct(survey_id, site, pair, type, period, Researcher, Date)

shannon_data <- left_join(shannon, meta, by = "survey_id") %>%
  mutate(
    date_num = as.numeric(Date),                 # continuous time for spline
    period   = factor(period, levels = c("Pre","Post")),
    type     = factor(type,   levels = c("Artificial","Natural"))
  )

# 2) Exploratory plots --------------------------------------------------------
# a) by reef type
s_type <- ggplot(shannon_data, aes(x = type, y = shannon, color = type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
  scale_color_manual(values = reef_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Reef Type", y = "Shannon diversity", title = "Shannon by reef type") +
  theme_minimal()

# b) by site (keep your preferred order if present)
site_order <- c("Aow Mao","Aow Mao Wreck","No Name Pinnacle","No Name Wreck","Hin Pee Wee","Sattakut")
shannon_data$site <- factor(shannon_data$site, levels = site_order)

s_site <- ggplot(shannon_data, aes(x = site, y = shannon, color = type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
  scale_color_manual(values = reef_cols) +
  labs(x = "Site", y = "Shannon diversity", title = "Shannon by site and reef type") +
  theme_minimal()

# c) by pair
s_pair <- ggplot(shannon_data, aes(x = site, y = shannon, color = type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_color_manual(values = reef_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ pair, scales = "free_x") +
  labs(x = "Site", y = "Shannon diversity", title = "Shannon by site within pairs") +
  theme_minimal()

# 3) Hierarchical Bayesian model ---------------------------------------------
priors_shannon <- c(
  set_prior("normal(1.8, 0.5)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("student_t(3, 0, 0.3)", class = "sds"),
  set_prior("student_t(3, 0, 0.3)", class = "sigma"),
  set_prior("student_t(3, 0, 0.5)", class = "sd")
)

fit_shannon_pair_spline <- brm_m2(
  shannon ~ s(date_num, by = interaction(pair, type), k = 6) +
    type * period + (1 | site),
  data   = shannon_data,
  family = gaussian(),
  prior  = priors_shannon
)

# 4) Outputs ------------------------------------------------------------------
summary(fit_shannon_pair_spline)
bayes_R2(fit_shannon_pair_spline)
s_check <- pp_check(fit_shannon_pair_spline)

# 5) Posterior predictions ----------------------------------------------------
pred_grid <- expand.grid(
  site = unique(shannon_data$site),
  type = levels(shannon_data$type),
  Date = seq(min(shannon_data$Date), max(shannon_data$Date), length.out = 200)
) %>%
  mutate(
    date_num = as.numeric(Date),
    period   = factor(if_else(Date < deployment_date, "Pre", "Post"), levels = c("Pre","Post"))
  ) %>%
  left_join(shannon_data %>% distinct(site, pair), by = "site")

preds <- add_epred_draws(fit_shannon_pair_spline, newdata = pred_grid, re_formula = NULL)

s_bayes <- ggplot(preds, aes(x = Date, y = .epred, color = type, fill = type)) +
  stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.2) +
  geom_vline(xintercept = deployment_date, linetype = "dashed", color = "gray40") +
  facet_wrap(~ pair, scales = "free_x") +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(x = "Date", y = "Predicted Shannon diversity",
       color = "Reef type", fill = "Reef type",
       title = "Bayesian Shannon trends by site pair") +
  theme_minimal(base_size = 12)

# 6) Posterior DiD contrast ---------------------------------------------------
post_did <- as_draws_df(fit_shannon_pair_spline) %>%
  select(`b_typeNatural:periodPost`)

summary_stats <- post_did %>%
  summarise(
    mean = mean(`b_typeNatural:periodPost`),
    median = median(`b_typeNatural:periodPost`),
    lower95 = quantile(`b_typeNatural:periodPost`, 0.025),
    upper95 = quantile(`b_typeNatural:periodPost`, 0.975),
    prob_less0 = mean(`b_typeNatural:periodPost` < 0),
    prob_greater0 = mean(`b_typeNatural:periodPost` > 0)
  )
print(summary_stats)

s_did <- ggplot(post_did, aes(x = `b_typeNatural:periodPost`)) +
  geom_density(fill = reef_cols["Natural"], alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  labs(x = "DiD contrast (Natural vs Artificial, Post–Pre)",
       y = "Posterior density",
       title = "Posterior distribution of the DiD effect") +
  theme_minimal(base_size = 12)

# 7) Save plots and stats -----------------------------------------------------
ggsave(file.path(output_dir, "Fig1_Shannon_by_ReefType.png"), s_type,  width = 6,  height = 4, dpi = 600)
ggsave(file.path(output_dir, "Fig2_Shannon_by_Site.png"),     s_site,  width = 8,  height = 4, dpi = 600)
ggsave(file.path(output_dir, "Fig3_Shannon_by_Pair.png"),     s_pair,  width = 8,  height = 4, dpi = 600)
ggsave(file.path(output_dir, "Fig4_Bayesian_Shannon_Trends.png"), s_bayes, width = 8, height = 5, dpi = 600)
ggsave(file.path(output_dir, "Fig5_Posterior_DiD.png"),       s_did,   width = 6,  height = 4, dpi = 600)
ggsave(file.path(output_dir, "FigS1_Shannon_checks.png"),     s_check, width = 6,  height = 4, dpi = 600)

stats_dir <- file.path(output_dir, "stats", "diversity")
if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

sink(file.path(stats_dir, "fit_shannon_pair_spline_summary.txt"))
print(summary(fit_shannon_pair_spline))
sink()

post_summary <- as_draws_df(fit_shannon_pair_spline) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "estimate") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(estimate),
    sd = sd(estimate),
    lower95 = quantile(estimate, 0.025),
    upper95 = quantile(estimate, 0.975),
    prob_less0 = mean(estimate < 0),
    prob_greater0 = mean(estimate > 0),
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = case_when(
      prob_greater0 > 0.95 ~ "Strong positive effect",
      prob_less0 > 0.95 ~ "Strong negative effect",
      prob_greater0 > 0.75 ~ "Likely positive effect",
      prob_less0 > 0.75 ~ "Likely negative effect",
      TRUE ~ "Uncertain / mixed"
    )
  ) %>% arrange(parameter)

write.csv(post_summary, file.path(stats_dir, "posterior_summary_table.csv"), row.names = FALSE)
write.csv(bayes_R2(fit_shannon_pair_spline), file.path(stats_dir, "bayes_R2_summary.csv"), row.names = FALSE)

loo_results <- loo(fit_shannon_pair_spline)
sink(file.path(stats_dir, "loo_summary.txt")); print(loo_results); sink()
message("✅ Diversity analysis saved to: ", output_dir)

s_bayes
