#######################################################
##   Stroke Mortality Data Cleaning and Preparation  ##
#######################################################

# ------------------------------
# 0. Environment Setup
# ------------------------------

rm(list = ls())             # Clear all objects from environment
options(scipen = 999)       # Disable scientific notation for readability

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Set the municipality event threshold here (e.g., 100, 200, 300, ...)
EVENT_THRESHOLD <- 100
# Filenames will include this suffix, e.g., *_100.*
suffix <- paste0("_", EVENT_THRESHOLD)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Load required libraries
library(easypackages)
library(dplyr)
library(readr)
library(ggplot2)

libraries(c(
  "tidyverse", "tidylog", "vroom", "lubridate",
  "stringr", "purrr", "tidyr", "sf", "geobr",
  "readxl", "ggrepel", "weathermetrics", "zoo", "survival",
  "dlnm", "splines", "mvmeta", "metafor", "gnm", "broom", "ggh4x",
  "rlist", "mixmeta", "ckbplotr", "janitor",
  "bannerCommenter", "qs"
))

# Ensure output directory exists
dir.create("cleaning_outputs", showWarnings = FALSE, recursive = TRUE)

# ------------------------------
# 1. Load Stroke Mortality Data
# ------------------------------

# Load dataset with stroke mortality and exposure variables
df <- qread("stroke_mortality_with_exp_overall.qs")
#df <- qread("stroke_ischemic_mortality_with_exp_overall.qs")
#df <- qread("stroke_bleed_mortality_with_exp_overall.qs")
#
# Quick structure check
colnames(df)
summary(df$umidade)   # Check humidity distribution

# ------------------------------
# 2. Data Cleaning
# ------------------------------

# Cap humidity values at 100% (safeguard against sensor errors)
df$umidade[df$umidade > 100] <- 100
summary(df$umidade)

# Create key time variables and case-crossover strata
df <- df %>%
  mutate(
    year     = lubridate::year(date),
    month    = lubridate::month(date),
    dow      = weekdays(date),
    # Define strata by Year × Month × Day-of-week
    stratum = interaction(year, month, dow, drop = TRUE),
    #stratum  = factor(year):factor(month):factor(dow),
    # Calculate apparent temperature (heat index) in °C
    app_temp = weathermetrics::heat.index(
      t = temperatura,
      rh = umidade,
      temperature.metric = "celsius"
    )
  )

head(df)

# Ensure municipality codes are character (6-digit codes)
df$loc6digit <- as.character(df$loc6digit)

# ------------------------------
# 3. Join with Municipality Metadata
# ------------------------------

# Load official municipality shapefile (2020) to get names/states
states <- geobr::read_municipality(code_muni = "all", year = 2020)
glimpse(states)

# Drop geometry and keep only identifiers
states_df <- sf::st_drop_geometry(states) %>%
  dplyr::select(code_muni, abbrev_state, name_muni) %>%
  dplyr::mutate(
    code_muni6 = stringr::str_sub(code_muni, 1, 6),                # 6-digit code
    state_loc  = paste0(abbrev_state, "_", name_muni)              # "STATE_municipality"
  )

# Join mortality data with municipality metadata
df2 <- dplyr::left_join(df, states_df, by = c("loc6digit" = "code_muni6"))

# ------------------------------
# 4. Event Count Threshold
# ------------------------------

# Count total stroke deaths per location
total_events <- df2 %>%
  group_by(state_loc) %>%
  summarise(total_events = sum(events), .groups = "drop") %>%
  filter(!is.na(total_events))

# Summary: how many municipalities in total, and how many meet threshold
n_total <- n_distinct(total_events$state_loc)
n_keep  <- sum(total_events$total_events >= EVENT_THRESHOLD, na.rm = TRUE)

cat("Total municipalities:", n_total, "\n")
cat("Municipalities with ≥", EVENT_THRESHOLD, " events:", n_keep, "\n", sep = "")

# Keep only municipalities meeting the chosen threshold
df3 <- df2 %>%
  filter(state_loc %in% total_events$state_loc[total_events$total_events >= EVENT_THRESHOLD])

# ==========================================================
# Map of total stroke deaths per municipality
# ==========================================================

library(scales)
library(RColorBrewer)

# Aggregate to municipality, compute deciles
deaths_by_muni <- df3 %>%
  mutate(loc6digit = stringr::str_pad(as.character(loc6digit), width = 6, side = "left", pad = "0")) %>%
  group_by(loc6digit) %>%
  summarise(total_deaths = sum(events, na.rm = TRUE), .groups = "drop") %>%
  mutate(decile_num = dplyr::ntile(total_deaths, 10))   # 10 groups (1=lowest, 10=highest)

# Join to shapefile
muni_sf <- geobr::read_municipality(code_muni = "all", year = 2020) %>%
  mutate(
    code_muni_chr = sprintf("%07d", as.integer(code_muni)),  # zero-pad 7 digits
    code_muni6    = substr(code_muni_chr, 1, 6)              # first 6 digits
  ) %>%
  left_join(
    deaths_by_muni %>%
      mutate(loc6digit = stringr::str_pad(as.character(loc6digit), 6, side = "left", pad = "0")),
    by = c("code_muni6" = "loc6digit")
  )

# Build range labels like "{123, 456}"
decile_ranges <- deaths_by_muni %>%
  group_by(decile_num) %>%
  summarise(min_d = min(total_deaths, na.rm = TRUE),
            max_d = max(total_deaths, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(decile_num) %>%
  mutate(label = paste0("{",
                        formatC(min_d, format = "d", big.mark = ""),
                        ", ",
                        formatC(max_d, format = "d", big.mark = ""),
                        "}"))

# map labels (no "D1")
label_map <- setNames(decile_ranges$label, decile_ranges$decile_num)

# factor with 1→10 order
muni_sf$decile <- factor(muni_sf$decile_num, levels = 1:10)


decile_palette <- brewer.pal(10, "Set3")   # includes a light grey; OK per your note

p_map_deciles <- ggplot() +
  # Municipalities not included (no decile) in very light grey background
  geom_sf(data = muni_sf %>% filter(is.na(decile)), fill = "grey100", color = NA) +
  # Municipalities with stroke deaths (deciles)
  geom_sf(data = muni_sf %>% filter(!is.na(decile)),
          aes(fill = decile), color = NA) +
  # State outlines
  geom_sf(data = geobr::read_state(year = 2020), fill = NA, color = "black", size = 0.3) +
  # Discrete fill scale with labels by decile range
  scale_fill_manual(
    values = decile_palette,
    labels = label_map[levels(muni_sf$decile)],
    drop   = FALSE,
    name   = "Stroke Deaths (2003–2023)"
  ) +
  coord_sf() +
  theme_void(base_size = 14) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    plot.title        = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# File names include the threshold suffix
ggsave(paste0("cleaning_outputs/map_stroke_deaths_deciles", suffix, ".png"),
       p_map_deciles, width = 10, height = 6, dpi = 300)

# ------------------------------
# 5. Save Cleaned Datasets
# ------------------------------

# Save full dataset with all locations
#write_csv(df2, "full_data.csv")

# Save restricted dataset with ≥ threshold events
#write_csv(df3, paste0("data", suffix, ".csv"))

# ------------------------------
# 6. Split into Location-Specific Data Lists
# ------------------------------

# Create list of data frames: one per municipality (all locations)
#data_list_all  <- group_split(df2, state_loc)
#names(data_list_all) <- sapply(data_list_all, function(x) unique(x$state_loc))
#saveRDS(data_list_all, "data_list_all.rds")

# Create list of data frames: one per municipality (≥ threshold only)
data_list_thr <- group_split(df3, state_loc)
names(data_list_thr) <- sapply(data_list_thr, function(x) unique(x$state_loc))
saveRDS(data_list_thr, paste0("data_list", suffix, ".rds"))

# ------------------------------
# 7. Descriptive Summaries
# ------------------------------

# National-level 
desc_national <- df3 %>%
  summarise(
    municipalities = n_distinct(state_loc),
    observations   = n(),
    total_events   = sum(events, na.rm = TRUE),
    HI_median      = median(app_temp, na.rm = TRUE),
    HI_IQR_low     = quantile(app_temp, 0.25, na.rm = TRUE),
    HI_IQR_high    = quantile(app_temp, 0.75, na.rm = TRUE)
  )
write_csv(desc_national, paste0("cleaning_outputs/descriptive_national", suffix, ".csv"))

# --- Region-level descriptive summary

# Attach region classification to df3
muni_meta <- geobr::read_municipality(code_muni = "all", year = 2020) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(abbrev_state, name_muni, name_region)

# Merge regions into dataset
df3 <- df3 %>%
  left_join(muni_meta, by = c("abbrev_state", "name_muni"))

desc_region <- df3 %>%
  group_by(name_region) %>%
  summarise(
    municipalities = n_distinct(state_loc),
    observations   = n(),
    total_events   = sum(events, na.rm = TRUE),
    HI_median      = median(app_temp, na.rm = TRUE),
    HI_IQR_low     = quantile(app_temp, 0.25, na.rm = TRUE),
    HI_IQR_high    = quantile(app_temp, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(desc_region, paste0("cleaning_outputs/descriptive_region", suffix, ".csv"))

# --- Plot 1: Distribution of apparent temperature by region
p_exposure_region <- ggplot(df3, aes(x = app_temp, fill = name_region)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.7) +
  facet_wrap(~name_region, scales = "free_y") +
  labs(
    title = "Distribution of Heat Index by Region",
    x = "Heat Index (°C)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(paste0("cleaning_outputs/exposure_distribution_by_region", suffix, ".png"),
       plot = p_exposure_region, width = 10, height = 6, dpi = 300)

# --- Monthly seasonality: build summary
monthly_summary <- df3 %>%
  group_by(year, month) %>%
  summarise(
    mean_HI     = mean(app_temp, na.rm = TRUE),
    mean_events = mean(events, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month_date = as.Date(paste(year, month, "01", sep = "-")))

# Monthly Heat Index seasonality
p_month_HI <- ggplot(monthly_summary, aes(x = month_date, y = mean_HI)) +
  geom_line(color = "#ef6548", size = 1.2) +
  geom_point(color = "#ef6548", size = 1) +
  labs(
    title = "Monthly Mean Heat Index",
    x = "Year",
    y = "Heat Index (°C)"
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

ggsave(paste0("cleaning_outputs/monthly_heat_index", suffix, ".png"),
       plot = p_month_HI, width = 10, height = 5, dpi = 300)

# Monthly Stroke Mortality seasonality
p_month_events <- ggplot(monthly_summary, aes(x = month_date, y = mean_events)) +
  geom_line(color = "#3182bd", size = 1.2) +
  geom_point(color = "#3182bd", size = 1) +
  labs(
    title = "Monthly Mean Stroke Deaths",
    x = "Year",
    y = "Daily Stroke Deaths (mean)"
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

ggsave(paste0("cleaning_outputs/monthly_stroke_deaths", suffix, ".png"),
       plot = p_month_events, width = 10, height = 5, dpi = 300)

# ------------------------------
# 8. Housekeeping
# ------------------------------

# Remove large objects
rm(df, df2, df3, data_list_thr, states, states_df, total_events)

# Free memory
gc()
