

######### Set Up  #########
fish_long <- readRDS("fish_long_cleaned.rds")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(brms)
  library(posterior)
  library(loo)
})

priors_nb <- c(
  prior(exponential(1), class = sd),
  prior(exponential(1), class = shape),            # NB overdispersion
  prior(normal(0, 2),   class = b),                # fixed effects
  prior(normal(0, 5),   class = Intercept)         # log-scale
)



# ---- Model 1: Functional Group × Type ----

fit_fun_group <- brm_m2(
  Count ~ type * fgroup + (1 | site),
  data = fish_long,
  family = negbinomial(), adapt_delta = 0.99
)
summary(fit_fun_group)
# ---- Model 2: Attraction (Type × Period × Pair) ----
fit_fun_pairs <- brm_m2(
  Count ~ type * period * pair + fgroup + (1 | site),
  data = fish_long,
  family = negbinomial()
)

# ---- Model 3: Time-varying Nonlinear ----
fit_fun_groups_time <- brm_m2(
  Count ~ type * fgroup +
    s(date_num, by = interaction(type, fgroup, pair), k = 4) +
    type * period + (1 | pair/site),
  data = fish_long,
  family = negbinomial()
)
summary(fit_fun_group_time)
# ---- Posterior contrasts ----
hypothesis(fit_fun_pairs, "typeNatural:periodPost < 0") # attraction
bayes_R2(fit_fun_groups_time)

# ---- Predictions and Plot ----
pred_grid <- expand.grid(
  type = levels(fish_long$type),
  fgroup = levels(fish_long$fgroup),
  pair = unique(fish_long$pair),
  Date = seq(min(fish_long$Date), max(fish_long$Date), length.out = 200)
) %>%
  mutate(date_num = as.numeric(Date - min(fish_long$Date)),
         period = factor(if_else(Date < as.Date("2023-09-07"), "Pre","Post"), levels=c("Pre","Post")))

preds <- add_epred_draws(fit_fun_groups_time, newdata = pred_grid, re_formula = NA)

ggplot(preds, aes(Date, .epred, color = type, fill = type)) +
  stat_lineribbon(.width=c(0.5,0.8,0.95), alpha=0.2) +
  facet_grid(pair ~ fgroup, scales="free_y") +
  geom_vline(xintercept=as.Date("2023-09-07"), linetype="dashed", color="gray40") +
  scale_color_manual(values=c("Natural"="#66BFA6","Artificial"="#007A87")) +
  scale_fill_manual(values=c("Natural"="#66BFA6","Artificial"="#007A87")) +
  labs(x="Date", y="Predicted fish count",
       title="Bayesian modelled fish abundance by functional group and site pair") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), strip.text=element_text(face="bold"))




# comparing old and new data 


oldlongfish <- readRDS("~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/Abund_Millie/AR_Attraction_Abund/Analysis_2025.07.25/oldlongfish.rds")
old_df <- oldlongfish
new_df <- fish_long
# rename some stuff in old df 
old_df <- old_df %>%
  rename(
    site   = Site,
    type   = Type,
    fgroup = Functional_Group,
    period = deployment_period
  ) %>%
  mutate(
    site   = factor(site,   levels = levels(new_df$site)),
    pair   = factor(pair,   levels = levels(new_df$pair)),
    type   = factor(type,   levels = levels(new_df$type)),
    period = factor(period, levels = levels(new_df$period)),
    fgroup = factor(fgroup, levels = levels(new_df$fgroup))
  )

# Both: drop Invertivores before the fit
old_df <- old_df |> dplyr::filter(fgroup != "Invertivore") |> droplevels()
new_df <- fish_long |> dplyr::filter(fgroup != "Invertivore") |> droplevels()

# IMPORTANT: define date_num from the SAME origin in BOTH sets
origin <- min(dplyr::bind_rows(old_df, new_df)$Date)
old_df <- old_df |> dplyr::mutate(date_num = as.numeric(Date - origin))
new_df <- new_df |> dplyr::mutate(date_num = as.numeric(Date - origin))

# Recode period deterministically from Date (avoid accidental differences)
deploy_date <- as.Date("2023-09-07")
old_df <- old_df |> dplyr::mutate(period = factor(ifelse(Date < deploy_date,"Pre","Post"),
                                                  levels = c("Pre","Post")))
new_df <- new_df |> dplyr::mutate(period = factor(ifelse(Date < deploy_date,"Pre","Post"),
                                                  levels = c("Pre","Post")))
### what changesd ####
cmp_effort <- bind_rows(
  old_df |> mutate(which="old"),
  new_df |> mutate(which="new")
) |>
  group_by(which, pair, type, fgroup, period) |>
  summarise(n_surveys = n(),
            date_min = min(Date), date_max = max(Date),
            mean_count = mean(Count), .groups="drop")

print(cmp_effort, n = Inf)
# A lightweight key to detect new rows (tweak if you have true IDs)
key_cols <- c("site","pair","type","fgroup","Date","survey_id","Count")
only_in_new <- dplyr::anti_join(new_df |> select(all_of(key_cols)) |> distinct(),
                                old_df |> select(all_of(key_cols)) |> distinct(),
                                by = key_cols)
nrow(only_in_new)  # how many new unique rows?

library(glmmTMB)

form <- Count ~ type * fgroup +
  type * (date_num + period) +
  pair * type * date_num +
  (1 | site)

fit_old <- glmmTMB(form, data = old_df, family = nbinom2())
fit_new <- glmmTMB(form, data = new_df, family = nbinom2())

co_old <- broom.mixed::tidy(fit_old, effects="fixed")[,c("term","estimate","std.error")]
co_new <- broom.mixed::tidy(fit_new, effects="fixed")[,c("term","estimate","std.error")]

coef_cmp <- dplyr::full_join(co_old |> rename(estimate_old=estimate, se_old=std.error),
                             co_new |> rename(estimate_new=estimate, se_new=std.error),
                             by="term") |>
  dplyr::arrange(term)

print(coef_cmp, n = Inf)

make_grid <- function(df) {
  expand.grid(
    type  = levels(df$type),
    fgroup= levels(df$fgroup),
    pair  = levels(df$pair),
    Date  = seq(min(df$Date), max(df$Date), length.out = 200)
  ) |>
    dplyr::mutate(
      date_num = as.numeric(Date - origin),
      period   = factor(ifelse(Date < deploy_date, "Pre","Post"),
                        levels = c("Pre","Post"))
    )
}

grid_old <- make_grid(old_df)
grid_new <- make_grid(new_df)

pred_old <- predict(fit_old, newdata = grid_old, type="response", re.form = NA)
pred_new <- predict(fit_new, newdata = grid_new, type="response", re.form = NA)

grid_old$pred <- pred_old
grid_new$pred <- pred_new

# testing between old and new data 
library(ggplot2)

plot_grazer <- function(gr) {
  gr |> dplyr::filter(fgroup=="Grazer") |>
    ggplot(aes(Date, pred, color = type, fill = type)) +
    geom_line(linewidth=1) +
    facet_wrap(~ pair, ncol = 1, scales = "free_y") +
    geom_vline(xintercept = deploy_date, linetype="dashed", color="gray40") +
    scale_color_manual(values=c(Artificial="#555D50", Natural="#0072B2")) +
    scale_fill_manual(values =c(Artificial="#555D50", Natural="#0072B2")) +
    labs(x="Date", y="Predicted fish count", color="Reef Type", fill="Reef Type") +
    theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(), strip.text=element_text(face="bold"))
}

p_old <- plot_grazer(grid_old) + ggtitle("Grazer by pair — OLD data")
p_new <- plot_grazer(grid_new) + ggtitle("Grazer by pair — NEW data")
library(patchwork) 
p_old +  p_new

library(dplyr)
library(tidyr)
library(tibble)

deploy_date <- as.Date("2023-09-07")  # just for reference if you want it

# Base summaries
summ_old <- old_df %>%
  filter(fgroup == "Grazer") %>%
  group_by(pair, type, period) %>%
  summarise(n = n(),
            date_min = min(Date), date_max = max(Date),
            span_days = as.integer(max(Date) - min(Date)),
            mean_count = mean(Count),
            .groups = "drop") %>%
  mutate(which = "old")

summ_new <- new_df %>%
  filter(fgroup == "Grazer") %>%
  group_by(pair, type, period) %>%
  summarise(n = n(),
            date_min = min(Date), date_max = max(Date),
            span_days = as.integer(max(Date) - min(Date)),
            mean_count = mean(Count),
            .groups = "drop") %>%
  mutate(which = "new")

# Side-by-side, with deltas
cmp <- bind_rows(summ_old, summ_new) %>%
  pivot_wider(
    names_from = which,
    values_from = c(n, date_min, date_max, span_days, mean_count),
    names_sep = "."
  ) %>%
  mutate(
    dn            = n.new - n.old,
    dspan_days    = span_days.new - span_days.old,
    dmean_count   = mean_count.new - mean_count.old,
    # percent changes (safe)
    pn            = 100 * (n.new - n.old) / pmax(n.old, 1),
    pmean_count   = 100 * (mean_count.new - mean_count.old) / pmax(mean_count.old, 1e-9)
  ) %>%
  arrange(pair, type, period) %>%
  as_tibble()

# Print everything
print(cmp, n = Inf)
# Biggest absolute changes in effort
cmp %>% arrange(desc(abs(dn))) %>% select(pair, type, period, n.old, n.new, dn) %>% print(n = Inf)

# Biggest shifts in mean counts
cmp %>% arrange(desc(abs(dmean_count))) %>% 
  select(pair, type, period, mean_count.old, mean_count.new, dmean_count, pmean_count) %>%
  print(n = Inf)

# New coverage windows
cmp %>%
  mutate(date_window.old = paste(format(date_min.old), format(date_max.old), sep = " → "),
         date_window.new = paste(format(date_min.new), format(date_max.new), sep = " → ")) %>%
  select(pair, type, period, date_window.old, date_window.new, dspan_days) %>%
  arrange(desc(abs(dspan_days))) %>%
  print(n = Inf)
library(ggplot2)

# mean count change
ggplot(cmp, aes(type, interaction(pair, period), fill = dmean_count)) +
  geom_tile() +
  scale_fill_gradient2(name = "Δ mean count", low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(x = "Reef type", y = "Pair × Period", title = "Change (NEW − OLD) in mean Grazer count") +
  theme_minimal()

# effort change
ggplot(cmp, aes(type, interaction(pair, period), fill = dn)) +
  geom_tile() +
  scale_fill_gradient2(name = "Δ n surveys", low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(x = "Reef type", y = "Pair × Period", title = "Change (NEW − OLD) in survey effort") +
  theme_minimal()






#### fungroups in bayesian #### 
# --- Packages ---
library(brms)
library(ggeffects)
library(dplyr)
library(ggplot2)

# Optional: ensure factors (period should be a factor if it was in glmmTMB)
valid_fun_groups <- valid_fun_groups %>%
  mutate(
    type   = factor(type),
    fgroup = factor(fgroup),
    pair   = factor(pair),
    site   = factor(site),
    period = factor(period)
  )

# --- brms fit (NB2 analogue) ---
# Matches: Count ~ type * fgroup + type*(date_num + period) + pair*type*date_num + (1|site)
bform <- bf(
  Count ~ type * fgroup +
    type * (date_num + period) +
    pair * type * date_num +
    (1 | site)
)

priors <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd"),
  set_prior("exponential(1)", class = "shape")  # NB dispersion
)

fun_groups_time_b <- brm(
  formula = bform,
  data = valid_fun_groups,
  family = negbinomial(link = "log"),
  prior = priors,
  chains = 4, iter = 4000, warmup = 1000, cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  threads = threading(2), backend= "cmdstanr"
  seed = 42
)

# --- Predictions (population-level / marginal means)  ---
fit_b <- fun_groups_time_b

preds_b <- ggpredict(
  fit_b,
  terms = c("date_num", "type", "fgroup", "pair"),
  type = "fixed"          # population-level (re_formula = NA under the hood)
)

# --- Plot formatting identical to glmmTMB version ---
origin_date <- min(fish_long$Date, na.rm = TRUE)
preds_b$x_date <- origin_date + preds_b$x

preds_df_b <- as.data.frame(preds_b) %>%
  rename(fgroup = facet, pair = panel) %>%
  mutate(fgroup = factor(fgroup), pair = factor(pair)) %>%
  droplevels()

p_time_b <- ggplot(preds_df_b, aes(x = x_date, y = predicted, color = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  facet_grid(pair ~ fgroup, scales = "free_y") +
  scale_x_date(expand = c(0, 0)) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = "Date",
    y = "Predicted fish count",
    color = "Reef type",
    fill = "Reef type"
  ) +
  theme_clean +
  theme(legend.position = "bottom")

p_time_b
