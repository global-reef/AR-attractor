## Functional groups between reef types
suppressPackageStartupMessages({
  library(tidyr)
  library(lme4)
  library(glmmTMB)
  library(emmeans)
  library(ggeffects)
  library(ggplot2)
  library(tidyverse)
})

# Output folders -------------------------------------------------------------
if (!exists("output_dir")) stop("Please define output_dir before running saves.")
plots_dir <- file.path(output_dir, "figures")
rds_dir   <- file.path(output_dir, "rds")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,   showWarnings = FALSE, recursive = TRUE)

fish_long$fgroup <- factor(
  fish_long$fgroup,
  levels = c("Grazer", "Invertivore", "Mesopredator", "HTLP"),
  ordered = FALSE
)

#### Model A: NB2 type × fgroup ##############################################
# Fit negative binomial GLMM
fun_groups <- glmmTMB(
  Count ~ type * fgroup + (1 | site),
  data = fish_long,
  family = nbinom2()
)

# Save model
saveRDS(fun_groups, file.path(rds_dir, "modelA_fun_groups_nb2.rds"))

# Model summary
summary(fun_groups)
# Type II ANOVA (optional)
car::Anova(fun_groups)
# Estimated marginal means (optional)
emmeans(fun_groups, ~ type | fgroup)

# Get predicted marginal means (on response scale)
emm <- emmeans(fun_groups, ~ type | fgroup, type = "response")
saveRDS(emm,     file.path(rds_dir, "modelA_emm_object.rds"))

# Convert to a dataframe for ggplot
emm_df <- as.data.frame(emm)
saveRDS(emm_df,  file.path(rds_dir, "modelA_emm_df.rds"))

# Plot (grouped by functional group)
ggplot(emm_df, aes(x = fgroup, y = response, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = c("Natural" = "#007A87", "Artificial" = "#66BFA6")) +
  labs(
    x = "Functional Group",
    y = "Estimated Fish Count (mean ± 95% CI)",
    fill = "Reef Type"
  ) + theme_clean
ggsave(file.path(plots_dir, "modelA_emm_by_fgroup_bar.png"),
       width = 7, height = 5, dpi = 300)

# Plot (faceted by functional group)
ggplot(emm_df, aes(x = type, y = response, fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(0.7), width = 0.2) +
  facet_wrap(~ fgroup, scales = "free_y") +
  scale_fill_manual(values = c("Natural" = "#007A87", "Artificial" = "#66BFA6")) +
  labs(
    x = "Reef Type",
    y = "Estimated Fish Count (mean ± 95% CI)",
    fill = "Reef Type"
  ) + theme_clean
ggsave(file.path(plots_dir, "modelA_emm_by_type_facets.png"),
       width = 8, height = 5.5, dpi = 300)

#### Model B: NB2 with time and period #######################################
# add time in
fish_long <- fish_long %>%
  mutate(date_num = as.numeric(Date - min(Date)))

fun_groups_time <- glmmTMB(
  Count ~ type * fgroup + 
    type * (date_num + period) + 
    (1 | site),
  data = fish_long,
  family = nbinom2()
)
# Save model B
saveRDS(fun_groups_time, file.path(rds_dir, "modelB_fun_groups_time_nb2.rds"))

summary(fun_groups_time)
car::Anova(fun_groups_time)


# Get predicted values across time, for each type × fgroup combo
preds <- ggpredict(fun_groups_time, terms = c("date_num", "type", "fgroup"))
saveRDS(preds, file.path(rds_dir, "modelB_ggeffects_preds.rds"))

# mutate back to real dates
origin_date <- min(fish_long$Date)
preds$x_date <- origin_date + preds$x

ggplot(preds, aes(x = x_date, y = predicted, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  facet_wrap(~ facet, scales = "free_y") +
  scale_color_manual(values = reef_cols)  +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = "Date",
    y = "Predicted Fish Count",
    color = "Reef Type",
    fill = "Reef Type",
    title = "Modelled Fish Abundance Over Time by Functional Group"
  ) + theme_clean
ggsave(file.path(plots_dir, "modelB_time_by_fgroup.png"),
       width = 9, height = 6, dpi = 300)

# Pairwise contrasts
# Estimate marginal means for Type at each date_num and fgroup
emm_time <- emmeans(fun_groups_time, ~ type | date_num * fgroup, at = list(date_num = unique(fish_long$date_num)))
saveRDS(emm_time, file.path(rds_dir, "modelB_emm_time_object.rds"))

# Get pairwise contrasts (Artificial - Natural) at each time point
contrasts <- contrast(emm_time, method = "revpairwise")  # gives AR - NR
contrasts_df <- as.data.frame(contrasts)

# Add actual dates
origin_date <- min(fish_long$Date)
contrasts_df$date <- origin_date + contrasts_df$date_num
contrasts_df <- contrasts_df %>%
  mutate(
    ratio = exp(estimate),
    lower = exp(estimate - 1.96 * SE),
    upper = exp(estimate + 1.96 * SE)
  )
saveRDS(contrasts_df, file.path(rds_dir, "modelB_contrasts_df.rds"))

# shows ratio of AR to NR
rel_abun <- ggplot(contrasts_df, aes(x = date, y = ratio)) +
  geom_line(color = "#66BFA6", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#66BFA6", alpha = 0.3) +
  facet_wrap(~ fgroup, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title = "Relative Abundance Ratio (Artificial / Natural)",
    x = "Date",
    y = "Abundance Ratio",
    subtitle = "Shaded ribbon shows 95% CI; dashed line = equal abundance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plots_dir, "modelB_relative_abundance_ratio.png"),
       plot = rel_abun, width = 9, height = 6, dpi = 300)

# raw data over time
fish_long %>%
  group_by(Date, type, fgroup) %>%
  summarise(mean_count = mean(Count), .groups = "drop") %>%
  ggplot(aes(x = Date, y = mean_count, color = type)) +
  geom_line() +
  facet_wrap(~ fgroup, scales = "free_y") +
  scale_color_manual(values = c("Natural" = "#007A87", "Artificial" = "#66BFA6")) +
  theme_minimal()
ggsave(file.path(plots_dir, "modelB_raw_means_over_time.png"),
       width = 9, height = 6, dpi = 300)

### Model C: testing and alternatives ########################################
fish_long <- fish_long %>%
  mutate(type = factor(type), fgroup = factor(fgroup), 
         pair = factor(pair), site = factor(site), period = factor(period))

valid_fun_groups <- fish_long %>%
  filter(fgroup != "Invertivore")

fun_groups_time <- glmmTMB(
  Count ~ type * fgroup +
    type * (date_num + period) +
    pair * type * date_num +
    (1 | site), 
  data = valid_fun_groups,
  family = nbinom2()) # can change to nbinom1() “overdispersed Poisson with constant proportional overdispersion.”
saveRDS(fun_groups_time, file.path(rds_dir, "modelC_fun_groups_time_nb2_pair_interact.rds"))

m_nb  <- glmmTMB(Count ~ type*fgroup + type*(ns(date_s,3)+period)
                 + (1+date_s|site),
                 family = nbinom1(), data = valid_fun_groups)
saveRDS(m_nb, file.path(rds_dir, "modelC_m_nb_nbinom1.rds"))

summary(fun_groups_time) # frequentist model
car::Anova(fun_groups_time)

fit <- fun_groups_time # final fit
# Predictions for type × Functional Group over time
preds <- ggpredict(fit, terms = c("date_num", "type", "fgroup", "pair")) # old using marginal means
saveRDS(preds, file.path(rds_dir, "modelC_ggeffects_preds.rds"))

# Convert numeric x to Date
origin_date <- min(fish_long$Date)
preds$x_date <- origin_date + preds$x

preds_df <- as.data.frame(preds) %>%
  rename(fgroup = facet, pair = panel) %>%
  mutate(fgroup = factor(fgroup), pair = factor(pair)) %>%
  droplevels()
saveRDS(preds_df, file.path(rds_dir, "modelC_preds_df.rds"))

p_time <- ggplot(preds_df, aes(x = x_date, y = predicted, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  facet_grid(pair ~ fgroup, scales = "free_y") +  # rows = pair, cols = functional group
  scale_x_date(expand = c(0, 0)) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = "Date",
    y = "Predicted Fish Count",
    color = "Reef Type",
    fill = "Reef Type") +
  theme_clean + theme(legend.position = "bottom")
p_time
ggsave(file.path(plots_dir, "modelC_time_by_pair_and_fgroup.png"),
       plot = p_time, width = 10, height = 8, dpi = 300)
