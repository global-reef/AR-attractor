
library(dplyr)
library(tidyr)
library(forcats)
library(lubridate)

clean_data <- function(file_path) {
  # ---- Load raw data ----
  raw_fish <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  
  # ---- Remove empty rows and columns ----
  raw_fish[raw_fish == ""] <- NA
  raw_fish <- raw_fish[, colSums(!is.na(raw_fish)) > 0]
  raw_fish <- raw_fish[rowSums(!is.na(raw_fish)) > 0, ]
  
  # ---- Fix column names ----
  raw_fish <- raw_fish %>%
    rename(
      sml_Grouper = Grouper.30,
      lrg_Grouper = Grouper.30.1,
      lrg_Snapper = Snapper.30
    )
  
  # ---- Identify species columns ----
  species_cols <- c("Parrotfish","Rabbitfish","Butterflyfish","Angelfish","Cleaner_Wrasse",
                    "Batfish","Thicklip","Red_Breast","Slingjaw","Sweetlips","Squirrel.Soldier",
                    "Triggerfish","Porcupine.Puffer","Ray","Brown_Stripe_Snapper","Russels_Snapper",
                    "lrg_Snapper","Eel","Trevally","Emperorfish","sml_Grouper","lrg_Grouper","Barracuda")
  species_cols <- species_cols[species_cols %in% colnames(raw_fish)]
  
  # ---- Fix data types ----
  raw_fish[species_cols] <- lapply(raw_fish[species_cols], as.numeric)
  raw_fish$Date <- as.Date(as.character(raw_fish$Date), format = "%m/%d/%Y")
  
  # ---- Pivot to long format ----
  fish_long <- raw_fish %>%
    pivot_longer(cols = all_of(species_cols),
                 names_to = "Species",
                 values_to = "Count")
  
  # ---- Assign functional groups ----
  functional_groups <- tibble::tribble(
    ~Species, ~Functional_Group,
    "Parrotfish","Grazer","Rabbitfish","Grazer","Butterflyfish","Grazer",
    "Angelfish","Invertivore","Cleaner_Wrasse","Invertivore","Batfish","Invertivore",
    "Thicklip","Invertivore","Red_Breast","Invertivore","Slingjaw","Invertivore",
    "Sweetlips","Invertivore","Squirrel.Soldier","Invertivore","Triggerfish","Invertivore",
    "Porcupine.Puffer","Invertivore","Ray","Mesopredator","Brown_Stripe_Snapper","Mesopredator",
    "Russels_Snapper","Mesopredator","lrg_Snapper","HTLP","Eel","Mesopredator",
    "Trevally","HTLP","Emperorfish","Mesopredator","sml_Grouper","Mesopredator",
    "lrg_Grouper","HTLP","Barracuda","HTLP"
  )
  
  fish_long <- fish_long %>%
    left_join(functional_groups, by = "Species") %>%
    mutate(Count = ceiling(replace_na(Count, 0)),
           Functional_Group = factor(Functional_Group,
                                     levels = c("Grazer","Invertivore","Mesopredator","HTLP")))
  
  # ---- Remove unwanted columns ----
  fish_long <- fish_long %>%
    select(-Time, -Duration, -Depth, -Visibility, -Weather, -Current, -Boats, -total_N) %>%
    filter(Researcher != "Keisha")
  
  # ---- Fix site naming ----
  fish_long$Site <- fct_recode(fish_long$Site, "No Name Pinnacle" = "No Name")
  
  # ---- Assign site pairs ----
  fish_long <- fish_long %>%
    mutate(pair = case_when(
      Site %in% c("Aow Mao","Aow Mao Wreck") ~ "Aow Mao",
      Site %in% c("No Name Pinnacle","No Name Wreck") ~ "No Name",
      Site %in% c("Hin Pee Wee","Sattakut") ~ "Sattakut"
    )) %>%
    filter(!is.na(pair)) %>%
    mutate(Type = recode(Type, "Artifical" = "Artificial"))
  
  # ---- Check structure ----
  message("✅ Site–pair structure:")
  print(fish_long %>% count(pair, Site))
  
  # ---- Deployment period flag ----
  deployment_date <- as.Date("2023-09-07")
  fish_long <- fish_long %>%
    mutate(deployment_period = case_when(
      pair == "Sattakut" ~ "Post",
      Date < deployment_date ~ "Pre",
      TRUE ~ "Post"
    ))
  
  # ---- Handle missing counts ----
  fish_long <- fish_long %>% mutate(Count = replace_na(Count, 0))
  
  # ---- Create unique survey ID ----
  fish_long <- fish_long %>%
    mutate(survey_id = paste(Site, Date, Researcher, sep = "_"))
  
  # ---- Compute months since deployment ----
  fish_long <- fish_long %>%
    mutate(months_since_deployment =
             if_else(Date < deployment_date, 0,
                     interval(deployment_date, Date) / months(1)))
  
  # ---- Apply standardized naming (for 00.SET.R) ----
  fish_long <- fish_long %>%
    mutate(
      period   = factor(deployment_period, levels = c("Pre","Post")),
      type     = factor(Type, levels = c("Artificial","Natural")),
      t_since  = pmax(0, months_since_deployment),
      pair     = factor(pair),
      site     = factor(Site),
      fgroup   = factor(Functional_Group,
                        levels = c("Grazer","Invertivore","Mesopredator","HTLP")),
      date_num = as.numeric(Date - min(Date))
    )
  
  # ---- Final selection and ordering ----
  fish_long <- fish_long %>%
    select(site, pair, type, period, t_since, date_num,
           Researcher, survey_id, Species, fgroup, Count, Date)
  
  # ---- Save cleaned file ----
  saveRDS(fish_long, "fish_long_cleaned.rds")
  message("✅ Saved 'fish_long_cleaned.rds' with standardized names and cleaned structure.")
  return(fish_long)
}

# ---- Run cleaning ----
fish_long <- clean_data(file_path)
