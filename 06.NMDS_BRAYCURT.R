### BRAY–CURTIS / NMDS / PERMANOVA / SIMPER ####################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(vegan); library(ggplot2)
  library(readr); library(forcats)
})


### Prepare Bray–Curtis input ###################################################
bray_df <- fish_long %>%
  group_by(survey_id, type, Species) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0) %>%
  arrange(survey_id)

meta_df <- bray_df %>% select(survey_id, type)
sp_mat  <- bray_df %>% select(-survey_id, -type)

### Bray–Curtis and NMDS ########################################################
set.seed(42)
bray_dist <- vegdist(sp_mat, method = "bray")
nmds      <- metaMDS(sp_mat, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)

nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  mutate(survey_id = bray_df$survey_id) %>%
  left_join(meta_df, by = "survey_id")

p_nmds <- ggplot(nmds_scores, aes(NMDS1, NMDS2, fill = type)) +
  stat_ellipse(geom = "polygon", alpha = 0.3) +
  geom_point(aes(color = type), size = 2) +
  scale_fill_manual(values = reef_cols) +
  scale_color_manual(values = reef_cols) +
  labs(title = "NMDS (Bray–Curtis)", subtitle = "Community composition by reef type") +
  theme_clean

ggsave(file.path(plots_dir, paste0("nmds_bray_", analysis_date, ".png")),
       p_nmds, width = 7.2, height = 5, dpi = 300)

### PERMANOVA and dispersion check #############################################
adonis_result <- adonis2(bray_dist ~ type + pair, data = meta_df,
                         strata = meta_df$pair, permutations = 999) # changed to constrain 

disp <- betadisper(bray_dist, meta_df$type)  # check homogeneous dispersion

sink(file.path(fits_dir, paste0("permanova_", analysis_date, ".txt"))); print(adonis_result); cat("\n\nBETADISPER:\n"); print(anova(disp)); sink()

### ANOSIM (optional corroboration) ############################################
anosim_result <- anosim(bray_dist, grouping = meta_df$type, permutations = 999)
sink(file.path(fits_dir, paste0("anosim_", analysis_date, ".txt"))); print(anosim_result); sink()

### SIMPER ######################################################################
simper_result <- simper(sp_mat, group = meta_df$type, permutations = 999)
simper_sum    <- summary(simper_result)

# Species lookup and image paths
spp_lookup <- tibble::tribble(
  ~Species,            ~Functional_Group, ~Genus,                ~Species_epithet, ~sci_name,
  "Parrotfish",        "Grazer",          "Scarus",              "spp.",           "Scarus spp.",
  "Rabbitfish",        "Grazer",          "Siganus",             "spp.",           "Siganus spp.",
  "Butterflyfish",     "Grazer",          "Chaetodon",           "spp.",           "Chaetodon spp.",
  "Angelfish",         "Invertivore",     "Pomacanthus",         "spp.",           "Pomacanthus spp.",
  "Cleaner_Wrasse",    "Invertivore",     "Labroides",           "dimidiatus",     "Labroides dimidiatus",
  "Batfish",           "Invertivore",     "Ephippidae",          "spp.",           "Ephippidae spp.",
  "Thicklip",          "Invertivore",     "Hemigymnus",          "melapterus",     "Hemigymnus melapterus",
  "Red_Breast",        "Invertivore",     "Cheilinus",           "fasciatus",      "Cheilinus fasciatus",
  "Slingjaw",          "Invertivore",     "Epibulus",            "insidiator",     "Epibulus insidiator",
  "Sweetlips",         "Invertivore",     "Diagramma/Plectorhinchus", "spp.",     "Diagramma/ Plectorhinchus spp.",
  "Squirrel.Soldier",  "Invertivore",     "Holocentridae",       "spp.",           "Holocentridae spp.",
  "Triggerfish",       "Invertivore",     "Balistidae",          "spp.",           "Balistidae spp.",
  "Porcupine.Puffer",  "Invertivore",     "Diodon/Tetraodon",    "spp.",           "Diodon/ Tetraodon spp.",
 "Ray",               "Mesopredator",    "Taeniura/Neotrygon",  "spp.",           "Taeniura/ Neotrygon spp.",
 "Russels_Snapper",       "Mesopredator",    "Lutjanus",            "russellii",           "Lutjanus russellii",
 "Brown_Stripe_Snapper",       "Mesopredator",    "Lutjanus",            "vitti",           "Lutjanus vitti",
  #  "sml_snapper",       "Mesopredator",    "Lutjanus",            "spp.",           "Lutjanus (<30cm) spp.",
  "lrg_Snapper",       "HTLP",            "Lutjanus",            "spp.",           "Lutjanus (>30cm) spp.",
  "Eel",               "Mesopredator",    "Gymnothorax",         "spp.",           "Gymnothorax spp.",
  "Trevally",          "HTLP",            "Caranx",              "spp.",           "Caranx spp.",
  "Emperorfish",       "Mesopredator",    "Lethrinus",           "spp.",           "Lethrinus spp.",
  "sml_Grouper",       "Mesopredator",    "Cephalopholis/Epinephelus", "spp.",     "Cephalopholis/ Epinephelus spp.",
  "lrg_Grouper",       "HTLP",            "Epinephelus",         "spp.",           "Epinephelus (>30cm)/ Plectropomus spp.",
  "Barracuda",         "HTLP",            "Sphyraena",           "spp.",           "Sphyraena spp."
) 

# helper: attach sci_name labels, fall back to original Species
with_sci_labels <- function(df, species_col = "Species") {
  df %>%
    left_join(select(spp_lookup, Species, sci_name), by = join_by(!!species_col == Species)) %>%
    mutate(Species_label = dplyr::coalesce(sci_name, as.character(.data[[species_col]])))
}

# save a copy for provenance
readr::write_csv(spp_lookup, file.path(output_dir, paste0("species_lookup_", analysis_date, ".csv")))

# Two-group case: first element holds Artificial vs Natural
sim_df <- simper_sum[[1]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  arrange(desc(average)) %>%
  slice(1:10) %>%
  mutate(
    Direction   = ifelse(ava > avb, "Artificial", "Natural"),
    cum_percent = round(cumsum(average) / sum(average) * 100, 1),
    Significant = ifelse(p <= 0.05, "p ≤ 0.05", "n.s.")
  ) %>%
  with_sci_labels("Species") %>%
  # reorder by contribution *before* converting to italics
  mutate(Species_label = forcats::fct_reorder(Species_label, average)) %>%
  mutate(Species_label_ital = paste0("italic('", as.character(Species_label), "')"))

write_csv(sim_df, file.path(output_dir, paste0("simper_top10_", analysis_date, ".csv")))

p_simper <- ggplot(sim_df,
                   aes(x = average, y = Species_label,
                       fill = Direction, alpha = Significant)) +
  geom_col(color = "black") +
  scale_fill_manual(values = reef_cols) +
  scale_alpha_manual(values = c("p ≤ 0.05" = 1, "n.s." = 0.5)) +
  # italicize labels without affecting order
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    title = "Top 10 species driving reef-type dissimilarity",
    subtitle = "Bray–Curtis SIMPER",
    x = "Average contribution",
    y = NULL,
    fill = "Higher on",
    alpha = "Significance"
  ) +
  theme_clean
p_simper

ggsave(file.path(plots_dir, paste0("simper_top10_", analysis_date, ".png")),
       p_simper, width = 7.2, height = 5, dpi = 300)

