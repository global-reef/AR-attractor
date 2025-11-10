
######### 0. Set Analysis Date & Create Output Folder #######

# Enter the date for this analysis
analysis_date <- "2025.11.10"  # update each run
# File path (automatically join folder + date + filename)
file_path <- paste0(
  "~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/AR_Attractor/DATA/",
  analysis_date,
  "_ArtificialReefs_MASTER - data.csv"
)
raw_fish <- read.csv(file_path, stringsAsFactors=TRUE, strip.white=TRUE) 


# Create a folder named with the date inside the working directory
output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


fits_dir   <- file.path(output_dir, "fits")
plots_dir  <- file.path(output_dir, "plots")

stem_hu <- file.path(fits_dir, paste0("hurdle_", analysis_date))


####### 2. Stan / brms optimization for Apple Silicon #######

# 1. Use the modern, faster backend
options(brms.backend = "cmdstanr")

# 2. Detect available cores, leave one free for macOS
library(parallel)
n_cores <- max(1, detectCores())
options(mc.cores = n_cores)

# 3. Compile Stan models with native M1/M2 support
# (run this once if not yet done)
# cmdstanr::set_cmdstan_path(cmdstanr::install_cmdstan())

# 4. Recommended defaults for speed & stability
brms_default <- list(
  chains = 4,
  cores = 4,                # or n_cores if you prefer auto-detect
  threads = threading(2),   # use both performance + efficiency cores
  iter = 4000,
  warmup = 2000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

# Optional helper function to apply these defaults automatically:
brm_m2 <- function(formula, data, family, ...) {
  do.call(brms::brm, c(list(formula = formula, data = data, family = family),
                       brms_default, list(...)))
}
### custom theme and colour palettes ###############################################################
theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank()
  )

# colour palettes
fg_cols <- c(
  "Grazer" = "#66c2a4",
  "Invertivore" = "#41b6c4",
  "Mesopredator" = "#2c7fb8",
  "HTLP" = "#253494"
) 
reef_cols <- c("Natural" = "#66BFA6", "Artificial" = "#007A87")

