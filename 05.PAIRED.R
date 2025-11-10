### PAIRWISE CONTRASTS AND ATTRACTION EFFECTS ###################################

suppressPackageStartupMessages({
  library(glmmTMB)
  library(broom.mixed)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

fish_long <- fish_long %>%
  mutate(type = relevel(type, ref = "Artificial"))
#### Ensure dirs exist ##########################################################
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

ensure_dir(output_dir)
ensure_dir(fits_dir)
ensure_dir(plots_dir)
### Helper functions ############################################################

write_txt <- function(x, path) cat(paste0(x, collapse = "\n"), file = path)

capture_summary <- function(model, path_txt, extra = NULL) {
  sum_text <- capture.output(summary(model))
  if (!is.null(extra))
    sum_text <- c(sum_text, "", "---- extra ----", capture.output(print(extra)))
  write_txt(sum_text, path_txt)
}

### 1) Global model #############################################################

model_pairwise <- glmmTMB(
  Count ~ type * period * pair + fgroup + (1 | site),
  data = fish_long,
  family = nbinom2()
)
capture_summary(model_pairwise,
                file.path(fits_dir, paste0("model_pairwise_", analysis_date, ".txt")))
saveRDS(model_pairwise,
        file.path(fits_dir, paste0("model_pairwise_", analysis_date, ".rds")))

### 2) Pairwise attraction effect per site pair #################################

valid_pairs <- setdiff(unique(fish_long$pair), "Sattakut")

pairwise_results <- map_dfr(valid_pairs, function(p) {
  df_sub <- filter(fish_long, pair == p)
  model <- glmmTMB(Count ~ type * period + fgroup + (1 | site),
                   data = df_sub, family = nbinom2())
  saveRDS(model, file.path(fits_dir, paste0("pair_model_", gsub("\\s+","_", p), "_", analysis_date, ".rds")))
  capture_summary(model, file.path(fits_dir, paste0("pair_model_", gsub("\\s+","_", p), "_", analysis_date, ".txt")))
  
  broom.mixed::tidy(model) %>%
    filter(term == "typeNatural:periodPost") %>%
    mutate(pair = p)
})

pairwise_results <- pairwise_results %>%
  mutate(
    estimate = -estimate,
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )

write_csv(pairwise_results %>%
            select(pair, estimate, std.error, conf.low, conf.high, p.value),
          file.path(output_dir, paste0("pairwise_attraction_", analysis_date, ".csv")))

p_pair <- ggplot(pairwise_results,
                 aes(x = estimate, y = reorder(pair, estimate))) +
  geom_point(color = "#66BFA6", size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#66BFA6") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Attraction Effect per Site Pair",
    x = "Interaction: Type × Deployment Period (log count)",
    y = "Site Pair",
    subtitle = "Positive values suggest stronger post-deployment increases on Artificial reefs"
  ) +
  theme_clean

ggsave(file.path(plots_dir, paste0("pairwise_attraction_", analysis_date, ".png")),
       p_pair, width = 8, height = 5, dpi = 300)

### 3) Attraction effect by functional group ####################################

groups <- unique(fish_long$fgroup)
combos  <- tidyr::expand_grid(pair_name = valid_pairs, group_name = groups)

pair_group_results <- pmap_dfr(combos, function(pair_name, group_name) {
  df_sub <- fish_long %>%
    filter(pair == pair_name, fgroup == group_name, !is.na(Count))
  if (nrow(df_sub) < 5) return(NULL)
  
  model <- tryCatch(
    glmmTMB(Count ~ type * period + (1 | site), data = df_sub, family = nbinom2()),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  
  saveRDS(model, file.path(fits_dir, paste0("pair_fg_",
                                            gsub("\\s+","_", pair_name), "_",
                                            gsub("\\s+","_", group_name), "_",
                                            analysis_date, ".rds")))
  
  broom.mixed::tidy(model) %>%
    filter(term == "typeNatural:periodPost") %>%
    mutate(pair = pair_name, fgroup = group_name)
})

pair_group_results <- pair_group_results %>%
  mutate(
    estimate = -estimate,
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )

write_csv(pair_group_results,
          file.path(output_dir, paste0("pair_group_attraction_", analysis_date, ".csv")))

p_pair_fg <- pair_group_results %>%
  filter(fgroup != "Invertivore") %>%
  ggplot(aes(x = estimate, y = reorder(pair, estimate), color = fgroup)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ fgroup, scales = "free_x") +
  scale_color_manual(values = c("Grazer" = "#007A87",
                                "Mesopredator" = "#66BFA6",
                                "HTLP" = "#5DA9E9")) +
  labs(
    title = "Attraction Effect (AR − NR Post) by Site Pair and Functional Group",
    subtitle = "Interaction: Type × Deployment Period (log count)",
    x = "Interaction estimate",
    y = "Site Pair",
    color = "Group"
  ) +
  theme_clean

ggsave(file.path(plots_dir, paste0("pair_group_attraction_", analysis_date, ".png")),
       p_pair_fg, width = 9, height = 5.5, dpi = 300)

### 4) Predicted total abundance over time ######################################

fish_totals <- fish_long %>%
  group_by(site, Date, type, pair, date_num, Species) %>%
  summarise(mean_count = mean(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(site, Date, type, pair, date_num) %>%
  summarise(Total_Count = ceiling(sum(mean_count)), .groups = "drop")

origin_date <- min(fish_totals$Date, na.rm = TRUE)

get_preds_site <- function(pair_name, data = fish_totals, origin = origin_date) {
  df <- data %>% filter(pair == pair_name) %>% droplevels()
  if (nrow(df) < 10) return(NULL)
  
  site_type <- df %>% distinct(site, type)
  
  m <- glmmTMB(Total_Count ~ type * date_num + (1 | site),
               data = df, family = nbinom2())
  saveRDS(m, file.path(fits_dir, paste0("time_total_",
                                        gsub("\\s+","_", pair_name), "_", analysis_date, ".rds")))
  capture_summary(m, file.path(fits_dir, paste0("time_total_",
                                                gsub("\\s+","_", pair_name), "_", analysis_date, ".txt")))
  
  newdat <- site_type %>%
    tidyr::crossing(
      date_num = seq(min(df$date_num, na.rm = TRUE),
                     max(df$date_num, na.rm = TRUE),
                     length.out = 200)
    )
  
  pr  <- predict(m, newdata = newdat, type = "link", se.fit = TRUE, re.form = NULL)
  fam <- family(m)
  eta <- pr$fit
  se  <- pr$se.fit
  
  newdat %>%
    mutate(
      predicted = fam$linkinv(eta),
      conf.low  = fam$linkinv(eta - 1.96 * se),
      conf.high = fam$linkinv(eta + 1.96 * se),
      Date      = as.Date(date_num, origin = origin),
      pair      = pair_name
    )
}

pairs <- c("Aow Mao", "No Name", "Sattakut")
all_preds <- purrr::map_dfr(pairs, get_preds_site)

write_csv(all_preds,
          file.path(output_dir, paste0("time_preds_totals_", analysis_date, ".csv")))

p_time <- ggplot(all_preds,
                 aes(Date, predicted, color = type, fill = type,
                     group = interaction(site, type))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ pair, ncol = 3, scales = "free_x") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values  = reef_cols) +
  labs(x = "Date", y = "Predicted total abundance",
       color = "Reef type", fill = "Reef type") +
  theme_clean

ggsave(file.path(plots_dir, paste0("time_preds_totals_", analysis_date, ".png")),
       p_time, width = 10, height = 4.5, dpi = 300)
