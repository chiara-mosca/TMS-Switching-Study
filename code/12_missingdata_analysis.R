# Exact missing data counts  
library(readxl)  
library(dplyr)

raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

cat("=== TOTAL ROWS ===\n")  
cat("N rows:", nrow(raw), "\n")  
cat("N patients:", n_distinct(raw$study_id), "\n\n")

# PHQ-9 total  
cat("=== PHQ-9 TOTAL SCORE ===\n")  
cat("Missing phq9_score:", sum(is.na(raw$phq9_score)),   
    "(", round(100*mean(is.na(raw$phq9_score)), 1), "%)\n\n")

# Item-level  
cat("=== PHQ-9 ITEMS ===\n")  
items <- paste0("phq9_", 1:9)  
item_names <- c("Anhedonia", "Depressed Mood", "Sleep", "Fatigue",  
                "Appetite", "Guilt", "Concentration", "Psychomotor",   
                "Suicidal Ideation")

for (i in seq_along(items)) {  
  n_miss <- sum(is.na(raw[[items[i]]]))  
  pct <- round(100 * n_miss / nrow(raw), 1)  
  cat(sprintf("  %s (%s): %d missing (%.1f%%)\n",   
              items[i], item_names[i], n_miss, pct))  
}

total_cells <- nrow(raw) * 9  
total_miss <- sum(sapply(items, function(x) sum(is.na(raw[[x]]))))  
cat(sprintf("\n  TOTAL: %d missing out of %d item-cells (%.1f%%)\n",  
            total_miss, total_cells, 100*total_miss/total_cells))

# Visits per patient  
cat("\n=== OBSERVATIONS PER PATIENT ===\n")  
obs_per_pt <- raw %>%   
  group_by(study_id, group) %>%   
  summarise(n_obs = n(), .groups = "drop")

cat("  Overall: median =", median(obs_per_pt$n_obs),   
    ", range =", min(obs_per_pt$n_obs), "-", max(obs_per_pt$n_obs), "\n")

obs_per_pt %>%  
  group_by(group) %>%  
  summarise(  
    n_patients = n(),  
    median_obs = median(n_obs),  
    min_obs = min(n_obs),  
    max_obs = max(n_obs),  
    .groups = "drop"  
  ) %>%  
  print()  
