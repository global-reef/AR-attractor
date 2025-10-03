
######### 0. Set Analysis Date & Create Output Folder #######

# Enter the date for this analysis
analysis_date <- "2025.10.03"  # Update these 3 for each analysis run
# file path (adjust date for correct date)
file_path <- "~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/AR_Attractor/DATA/2025.10.03_ArtificialReefs_MASTER - data.csv"
# file_path_timed <- "~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Pelagic_Pinnacles/2_DATA/2025.06.19_TimedFishSurveys_Shallow_MASTER.csv"
# --------------------------
raw_fish <- read.csv(file_path, stringsAsFactors=TRUE, strip.white=TRUE) 


# Create a folder named with the date inside the working directory
output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
