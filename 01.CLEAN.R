#Load Packages
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(forcats)


clean_data <- function(file_path) {
  # Load the raw data with strings as factors
  # browser() # for troubleshooting 
  df <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  # remove blank rows and columns 
  df[df == ""] <- NA
  df <- df[, colSums(!is.na(df)) > 0]
  df <- df[rowSums(!is.na(df)) > 0, ]
  
  # Rename columns
  df <- df %>%
    rename(sml_Grouper = Grouper.30,
           lrg_Grouper = Grouper.30.1,
           lrg_Snapper = Snapper.30)
  
  # Define the species columns of interest and keep only those that exist in the data
  species_cols <- c("Parrotfish", "Rabbitfish", "Butterflyfish", "Angelfish", "Cleaner_Wrasse",
                    "Batfish", "Thicklip", "Red_Breast", "Slingjaw", "Sweetlips", "Squirrel.Soldier",
                    "Triggerfish", "Porcupine.Puffer", "Ray", "Brown_Stripe_Snapper", 
                    "Russels_Snapper", "lrg_Snapper", "Eel", "Trevally", "Emperorfish",
                    "sml_Grouper", "lrg_Grouper", "Barracuda")
  species_cols <- species_cols[species_cols %in% colnames(df)]
  
  # Fix Formats
  df[species_cols] <- lapply(df[species_cols], as.numeric)
  df$Date <- as.Date(as.character(df$Date), format = "%m/%d/%Y")
  df$Time <- format(as.POSIXct(df$Time, format = "%H:%M"), "%H:%M")
  df$Weather <- as.factor(df$Weather)
  
  # Pivot longer: create a long-format data frame with one row per survey-species
  fish_long <- df %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "Species",
      values_to = "Count"
    )
  
  # Define functional groups
  functional_groups <- data.frame(
    Species = c("Parrotfish", "Rabbitfish", "Butterflyfish", "Angelfish", "Cleaner_Wrasse",
                "Batfish", "Thicklip", "Red_Breast", "Slingjaw", "Sweetlips", "Squirrel.Soldier",
                "Triggerfish", "Porcupine.Puffer", "Ray", "Brown_Stripe_Snapper", 
                "Russels_Snapper", "lrg_Snapper", "Eel", "Trevally", "Emperorfish",
                "sml_Grouper", "lrg_Grouper", "Barracuda"),
    Functional_Group = c("Grazer", "Grazer", "Grazer", "Invertivore", "Invertivore",
                         "Invertivore", "Invertivore", "Invertivore", "Invertivore", "Invertivore",
                         "Invertivore", "Invertivore", "Invertivore", "Mesopredator", "Mesopredator",
                         "Mesopredator", "HTLP", "Mesopredator", "HTLP", "Mesopredator",
                         "Mesopredator", "HTLP", "HTLP")
  )
  
  # Join the functional groups and set order
  fish_long <- fish_long %>%
    left_join(functional_groups, by = "Species") %>%
    mutate(Functional_Group = factor(Functional_Group,
                                     levels = c("Grazer", "Invertivore", "Mesopredator", "HTLP"),
                                     ordered = TRUE))
  # round to the nearest ineger 
  fish_long <- fish_long %>% 
    mutate(Count = ceiling(Count))
  # Create a unique survey ID if not already present
  if(!"survey_id" %in% colnames(fish_long)) {
    fish_long <- fish_long %>%
      mutate(survey_id = paste(Site, Date, sep = "_"))
  }
  return(fish_long)
}


# Use the function to clean the data 
fish_long <- clean_data(file_path)
fish_long$Duration <- as.integer(fish_long$Duration)

# get rid of columns we don't need for this analyiss 
fish_long <- fish_long %>%
  select(-Time, -Duration, -Depth, -Visibility, -Weather, -Current, -Boats, -total_N)


# remove outliers (barracudas, rays, porcupine/puffers, and eels # NOT DONE RN 
fish_long <- fish_long %>%
  # filter(!Species %in% c("Barracuda", "Eel", "Porcupine.Puffer", "Ray")) %>% # for now taking this out, lets see how it goes with barracudas
  filter(Researcher != "Keisha")

# check sites 
table(fish_long$Site)
# rename no name 
fish_long$Site <- fct_recode(fish_long$Site, "No Name Pinnacle" = "No Name")

# retain sites and apply paired column 
fish_long <- fish_long %>%
  mutate(pair = case_when(
    Site %in% c("Aow Mao", "Aow Mao Wreck") ~ "Aow Mao",
    Site %in% c("No Name Pinnacle", "No Name Wreck") ~ "No Name",
    Site %in% c("Hin Pee Wee", "Sattakut") ~ "Sattakut"
  )) %>%
  filter(!is.na(pair))
fish_long <- fish_long %>%
  mutate(Type = recode(Type,
                       "Artifical" = "Artificial")) %>%
  mutate(Type = factor(Type))   # make sure it's a proper factor

# check sites and pairs 
fish_long %>%
  count(pair, Site) %>%
  arrange(pair, Site)

# add column to indicate pre-vs post deployment 
fish_long <- fish_long %>%
  mutate(deployment_period = case_when(
    pair == "Sattakut" ~ "Post",  # override
    Date < as.Date("2023-09-07") ~ "Pre",
    TRUE ~ "Post"
  ))


# fix a weird NA 
fish_long <- fish_long %>%
  mutate(Count = replace_na(Count, 0))

# Save cleaned fish_long object
saveRDS(fish_long, file = "fish_long_cleaned.rds")

