### exploratory analysis 
library(dplyr)
library(ggplot2)


# total abund between reef types 
fish_long %>%
  group_by(Type) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE))

# total abund by site type 
fish_long %>%
  group_by(pair, Site) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE)) %>%
  arrange(factor(pair, levels = c("Aow Mao", "No Name", "Hin Pee Wee")), Site)

# number of surveys per site 
fish_long %>%
  distinct(survey_id, Site) %>%
  count(Site) %>%
  arrange(factor(Site, levels = c("Aow Mao", "Aow Mao Wreck",
                                  "No Name Pinnacle", "No Name Wreck",
                                  "Hin Pee Wee", "Sattakut")))

# Visualize raw counts to check variability
ggplot(fish_long, aes(x = Type, y = Count)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 2) +
  facet_wrap(~ Functional_Group, scales = "free_y") +
  labs(
    title = "Fish Count Distributions by Functional Group and Reef Type",
    x = "Reef Type",
    y = "Raw Count"
  ) +
  theme_minimal(base_size = 12)


# check where all of the high mesopredator counts on the natural reefs are coming from 
high_mesos <- fish_long %>%
  filter(Functional_Group == "Mesopredator", Type == "Natural") %>%
  arrange(desc(Count))
high_mesos %>%
  select(Site, Date, Species, Count) %>%
  head(10)
high_mesos %>%
  group_by(Site) %>%
  summarise(
    n = n(),
    max_count = max(Count),
    mean_count = mean(Count),
    total = sum(Count)
  ) %>%
  arrange(desc(max_count))
high_mesos %>%
  group_by(Species) %>%
  summarise(
    max_count = max(Count),
    mean_count = mean(Count),
    total = sum(Count)
  ) %>%
  arrange(desc(total))
high_mesos %>%
  filter(Site == "Hin Pee Wee", Species %in% c("Russels_Snapper", "Emperorfish")) %>%
  ggplot(aes(x = Date, y = Count, color = Species)) +
  geom_point() +
  geom_line() +
  labs(title = "High Mesopredator Counts at Hin Pee Wee (Natural Reef)")

