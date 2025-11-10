#### 03_DIVERSITY #######################
suppressPackageStartupMessages({
  library(vegan)
  library(lmerTest)
  library(car)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

#### diversity_run() ####
diversity_run <- function(fish_long,
                          reef_cols = c("Natural" = "#66BFA6", "Artificial" = "#007A87"),
                          output_dir = NULL,
                          save_plots = TRUE,
                          show_plots = FALSE) {
  #### 1) Data prep ####
  d <- fish_long %>%
    dplyr::mutate(
      Date   = if (inherits(Date, "Date")) Date else as.Date(as.character(Date)),
      type   = factor(type,   levels = c("Artificial","Natural")),
      period = factor(period, levels = c("Pre","Post")),
      site   = factor(site)
    )
  
  # per-survey species matrix
  sh_wide <- d %>%
    dplyr::group_by(survey_id, Species) %>%
    dplyr::summarise(count = sum(Count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Species, values_from = count,
                       values_fill = list(count = 0))
  
  # numeric matrix for vegan::diversity
  mat <- sh_wide %>%
    dplyr::select(-survey_id) %>%
    data.matrix() %>%
    as.matrix()
  
  shannon_data <- tibble::tibble(
    survey_id = sh_wide$survey_id,
    shannon   = vegan::diversity(mat)
  ) %>%
    dplyr::left_join(
      d %>% dplyr::distinct(survey_id, site, pair, type, period, Researcher, Date),
      by = "survey_id"
    ) %>%
    dplyr::mutate(
      date_num = as.numeric(Date),
      site     = droplevels(site),
      pair     = droplevels(pair)
    ) %>%
    dplyr::filter(!is.na(shannon), !is.na(site), !is.na(type), !is.na(period), !is.na(Date))
  
  #### 2) Plots: type, site within pairs ####
  p_type <- shannon_data %>%
    ggplot(aes(type, shannon, color = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 0.8) +
    scale_color_manual(values = reef_cols) +
    labs(x = "Reef type", y = "Shannon", title = "Shannon by reef type", color="Reef Type") +
    theme_clean
  
  p_pair <- shannon_data %>%
    ggplot(aes(site, shannon, color = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 0.9) +
    scale_color_manual(values = reef_cols) +
    facet_wrap(~ pair, scales = "free_x") +
    labs(x = "Site", y = "Shannon", color = "Reef Type",
         title = "Shannon diversity by site within pairs") +
    theme_clean +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  #### 3) Main model (time-varying mixed model) ####
  m_lmm <- lmerTest::lmer(
    shannon ~ type * (date_num + period) + (1 + date_num || site),
    data = shannon_data, REML = TRUE
  )
  
  # DiD-style contrast: Natural × Post
  coef_tab <- as.data.frame(summary(m_lmm)$coefficients)
  coef_tab$term <- rownames(coef_tab)
  did <- coef_tab %>%
    dplyr::filter(term == "typeNatural:periodPost") %>%
    dplyr::transmute(
      term,
      estimate  = Estimate,
      std.error = `Std. Error`,
      df        = df,
      t.value   = `t value`,
      p.value   = `Pr(>|t|)`,
      conf.low  = estimate + qt(0.025, df) * std.error,
      conf.high = estimate + qt(0.975, df) * std.error
    )
  
  #### 4) Site-specific cutoffs and predictions for visualization ####
  site_info <- shannon_data %>%
    dplyr::distinct(site, pair, type)
  
  cutoffs <- shannon_data %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(
      cutoff = {
        post_dates <- Date[period == "Post"]
        if (length(post_dates) == 0) as.Date(NA) else min(post_dates)
      },
      .groups = "drop"
    )
  
  site_info <- dplyr::left_join(site_info, cutoffs, by = "site")
  
  date_seq <- seq(min(shannon_data$Date), max(shannon_data$Date), length.out = 200)
  
  pred_grid <- site_info %>%
    dplyr::rowwise() %>%
    dplyr::mutate(data = list(tibble::tibble(Date = date_seq))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(
      period   = factor(dplyr::if_else(!is.na(cutoff) & Date < cutoff, "Pre", "Post"),
                        levels = c("Pre","Post")),
      date_num = as.numeric(Date)
    )
  
  pred_grid$fit_lmm <- predict(m_lmm, newdata = pred_grid,
                               re.form = NULL, allow.new.levels = TRUE)
  
  # pair-level cutoffs; keep only pairs with a genuine Pre->Post change
  pair_lines <- shannon_data %>%
    group_by(pair) %>%
    summarise(
      has_pre  = any(period == "Pre"),
      has_post = any(period == "Post"),
      pair_cut = if (has_post) min(Date[period == "Post"]) else as.Date(NA),
      .groups = "drop"
    ) %>%
    filter(has_pre & has_post & !is.na(pair_cut))
  
  p_trend <- ggplot(pred_grid, aes(Date, fit_lmm, color = type)) +
    geom_line(aes(group = interaction(site, type)), linewidth = 0.9) +
    # only draw where a cutoff exists (Aow Mao, No Name), none for Sattakut
    geom_vline(data = pair_lines, aes(xintercept = pair_cut),
               inherit.aes = FALSE, linetype = "dashed", color = "gray40") +
    facet_wrap(~ pair, scales = "free_x") +
    scale_color_manual(values = reef_cols) +
    labs(title = "Shannon trends by site pair (LMM fits)",
         x = "Date", y = "Shannon", color = "Reef Type") +
    theme_clean +
    geom_point(data = shannon_data, inherit.aes = FALSE,
               aes(Date, shannon, color = type), alpha = 0.45, size = 1)
  print(p_trend)
  
  #### 5) Save / show ####
  if (!is.null(output_dir) && save_plots) {
    dir.create(file.path(output_dir, "diversity"), recursive = TRUE, showWarnings = FALSE)
    od <- file.path(output_dir, "diversity")
    ggsave(file.path(od, "D1_Shannon_by_type.png"), p_type,  width = 6, height = 4,   dpi = 600)
    ggsave(file.path(od, "D2_Shannon_by_site_within_pairs.png"), p_pair, width = 9, height = 5, dpi = 600)
    ggsave(file.path(od, "D3_Shannon_trends_by_pair.png"), p_trend, width = 9, height = 5, dpi = 600)
    capture.output({
      cat("\nLMM summary\n"); print(summary(m_lmm))
      cat("\nType II/III (car::Anova) for LMM\n"); print(car::Anova(m_lmm))
      cat("\nDiD contrast (typeNatural:periodPost)\n"); print(did)
    }, file = file.path(od, "diversity_model_summary.txt"))
  }
  if (show_plots) { print(p_type); print(p_pair); print(p_trend) }
  
  #### 6) Return ####
  list(
    data   = shannon_data %>% dplyr::select(survey_id, site, pair, type, period, Date, shannon, date_num),
    model  = m_lmm,
    did    = did,
    plots  = list(type = p_type, pair = p_pair, trends = p_trend)
  )
}

# Run
div <- diversity_run(fish_long, output_dir = output_dir, show_plots = TRUE)
div$did
