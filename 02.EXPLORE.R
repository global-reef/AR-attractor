# data exploration based on Zuur et al. 
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lubridate)
  library(glmmTMB); library(lme4); library(DHARMa); library(readr)
})

zuur_eda <- function(dat, output_dir = NULL,
                     save_plots = TRUE, save_csv = TRUE, show_plots = FALSE,
                     species_min_n = 50) {
  
  d <- dat %>%
    mutate(
      Date   = if (inherits(Date,"Date")) Date else as.Date(as.character(Date)),
      day_id = interaction(site, Date, drop = TRUE)
    )
  
  ## ---- effort (unique surveys) ----
  surv <- d %>% distinct(site, type, period, Date, researcher = Researcher, survey_id)
  effort_by_site  <- surv %>% count(site, name = "n_surveys")
  effort_hist     <- ggplot(effort_by_site, aes(n_surveys)) +
    geom_histogram(bins = 30) + labs(title = "Surveys per site")
  
  ## ---- species-level zero inflation ----
  species_zero <- d %>%
    group_by(Species) %>%
    summarise(
      n = n(),
      detections = sum(Count > 0),
      zero_rate = mean(Count == 0),
      mean_count_pos = ifelse(detections > 0, mean(Count[Count > 0]), NA_real_),
      .groups = "drop"
    ) %>%
    arrange(desc(zero_rate), desc(n))
  
  # flag likely problem species (enough data + very zero-inflated)
  species_flags <- species_zero %>%
    filter(n >= species_min_n, zero_rate >= 0.9)
  
  # per-species detection histogram (number of nonzero observations per species)
  p_species_detect_hist <- ggplot(species_zero, aes(detections)) +
    geom_histogram(bins = 40) +
    labs(title = "Detections per species (nonzero records)")
  
  ## ---- per-survey richness (observed taxa per survey_id) ----
  richness <- d %>%
    group_by(survey_id) %>%
    summarise(richness = n_distinct(Species[Count > 0]), .groups = "drop")
  p_richness <- ggplot(richness, aes(richness)) +
    geom_histogram(bins = 40) + labs(title = "Species richness per survey")
  
  ## ---- core metrics (as before, minimal) ----
  zero_rate_all <- mean(d$Count == 0, na.rm = TRUE)
  m_pois <- glm(Count ~ 1, data = d, family = poisson())
  disp   <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
  
  d$log1pC <- log1p(d$Count)
  m_day  <- lmer(log1pC ~ 1 + (1|day_id), data = d)
  sd_day <- as.data.frame(VarCorr(m_day))$sdcor[1]
  sd_res <- attr(VarCorr(m_day), "sc")
  icc_day <- sd_day^2 / (sd_day^2 + sd_res^2)
  
  # quick NB residual scan
  form_nb <- Count ~ type * fgroup + period + t_since + date_num + (1|site) + (1|day_id)
  m_nb <- glmmTMB(form_nb, data = d, family = nbinom2())
  sim  <- simulateResiduals(m_nb, n = 500)
  
  # minimal plots reused from earlier EDA
  p_hist <- ggplot(d, aes(Count)) + geom_histogram(bins = 50) +
    coord_cartesian(xlim = c(0, quantile(d$Count, 0.99, na.rm = TRUE))) +
    labs(title = "Count distribution")
  by_day <- d %>% group_by(day_id) %>% summarise(m = mean(Count), v = var(Count), .groups="drop")
  p_mv <- ggplot(by_day, aes(m, v)) + geom_point(alpha=.6) + geom_smooth(se=FALSE) +
    labs(title = "Mean vs variance by day")
  
  # save / show
  if (!is.null(output_dir) && (save_plots || save_csv)) dir.create(output_dir, TRUE, FALSE)
  if (!is.null(output_dir) && save_plots) {
    ggsave(file.path(output_dir, "EDA_hist.png"), p_hist, width=6, height=4, dpi=300)
    ggsave(file.path(output_dir, "EDA_mean_variance.png"), p_mv, width=6, height=4, dpi=300)
    ggsave(file.path(output_dir, "EDA_surveys_per_site.png"), effort_hist, width=6, height=4, dpi=300)
    ggsave(file.path(output_dir, "EDA_species_detections_hist.png"), p_species_detect_hist, width=6, height=4, dpi=300)
    ggsave(file.path(output_dir, "EDA_richness_per_survey.png"), p_richness, width=6, height=4, dpi=300)
    png(file.path(output_dir, "EDA_DHARMa_nb.png"), width=1200, height=900, res=150); plot(sim); dev.off()
  }
  if (!is.null(output_dir) && save_csv) {
    write_csv(effort_by_site, file.path(output_dir, "effort_by_site.csv"))
    write_csv(species_zero,  file.path(output_dir, "species_zero_summary.csv"))
    write_csv(species_flags, file.path(output_dir, "species_zero_flags.csv"))
    write_csv(richness,      file.path(output_dir, "survey_richness.csv"))
  }
  if (show_plots) { print(effort_hist); print(p_species_detect_hist); print(p_richness) }
  
  list(
    metrics = tibble(
      n_rows = nrow(d),
      n_surveys = nrow(surv),
      n_sites = n_distinct(surv$site),
      n_days = n_distinct(surv$Date),
      zero_rate = zero_rate_all,
      poisson_dispersion = disp,
      icc_day = icc_day
    ),
    effort = list(by_site = effort_by_site),
    species = list(summary = species_zero, flags = species_flags),
    survey = list(richness = richness)
  )
}

# run
eda <- zuur_eda(fish_long, output_dir = if (exists("output_dir")) file.path(output_dir, "EDA") else NULL)
eda$metrics
eda$species$flags   # species with ≥ species_min_n rows and zero_rate ≥ 0.9 = eels and rays, consider removing 
eda$effort # surveys per site 



####  results of eda
# Zero rate ~ 0.48 → many zeros. Expect hurdle/two-stage to help.
# Poisson dispersion ~ 86 → Poisson out. Use NB family.
# ICC_day ~ 0.037 → small but real day clustering. Keep (1|day_id).
# nb_overdispersion_ok = TRUE just says NB is handling dispersion.



# checking if we need day_id as a random effect  ##### 

df <- fish_long %>% mutate(
  Date   = if (inherits(Date,"Date")) Date else as.Date(as.character(Date)),
  day_id = interaction(site, Date, drop = TRUE)
)

m_no_day <- glmmTMB(Count ~ type*fgroup + period + t_since + date_num + (1|site),
                    data = df, family = nbinom2())
m_day    <- glmmTMB(Count ~ type*fgroup + period + t_since + date_num + (1|site) + (1|day_id),
                    data = df, family = nbinom2())
AIC(m_no_day, m_day)   # keep day

###### basic summary stats ###### 

# number of surveys per site 
fish_long %>%
  distinct(survey_id, site) %>%
  count(site) %>%
  arrange(factor(site, levels = c("Aow Mao", "Aow Mao Wreck",
                                  "No Name Pinnacle", "No Name Wreck",
                                  "Hin Pee Wee", "Sattakut")))
# total abund between reef types 
fish_long %>%
  group_by(type) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE))

# total abund by site type 
fish_long %>%
  group_by(pair, site) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE)) %>%
  arrange(factor(pair, levels = c("Aow Mao", "No Name", "Hin Pee Wee")), site)
