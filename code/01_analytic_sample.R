# =============================================================  
# 01_analytic_sample.R  
# Purpose: Create analytic sample (N=85) from N=113 TMS patients  
# Author: Chiara Mosca  
# Date: August 2026  
#  
# BEFORE FIRST USE - run these lines in the Console:  
#   install.packages("tidyverse")  
#   install.packages("readxl")  
#   install.packages("writexl")  
#  
# INPUT FILES (must be in same folder as this script):  
#   tms_clean_longitudinal.xlsx  (sheet: "in")  
#   patients_clean.xlsx  
#   tms_clean_outcomes(in).csv  
#  
# OUTPUT FILES:  
#   analytic_sample_long.xlsx  
#   analytic_sample_summary.xlsx  
#   excluded_patients.xlsx  
#  
# DECISIONS:  
#   - iTBS-only group: patients who only ever received protocol 1  
#   - Switcher group: patients who went from protocol 1 to protocol 2  
#   - Excluded TMS093 (started on protocol 2, switched to 1)  
#   - Excluded TMS084 (alternated 1 -> 2 -> 1 -> 2)  
#   - Excluded patients with protocols 4, 5, or mixed combinations  
# =============================================================

# ---- Load packages ----  
library(tidyverse)  
library(readxl)  
library(writexl)

# ---- Load data files ----  
long <- read_xlsx("tms_clean_longitudinal.xlsx", sheet = "in")  
demo <- read_xlsx("patients_clean.xlsx")  
summary_df <- read_csv("tms_clean_outcomes(in).csv", show_col_types = FALSE)

cat("=== FILES LOADED ===\n")  
cat("Longitudinal:", nrow(long), "rows,", n_distinct(long$study_id), "patients\n")  
cat("Demographics:", nrow(demo), "rows\n")  
cat("Summary:", nrow(summary_df), "rows\n")

# ---- Identify protocol combinations per patient ----  
patient_protocols <- long %>%  
  filter(!is.na(protocol_doa)) %>%  
  group_by(study_id) %>%  
  summarise(  
    protocols = paste(sort(unique(protocol_doa)), collapse = ","),  
    protocol_sequence = paste(protocol_doa, collapse = " -> "),  
    .groups = "drop"  
  )

cat("\n=== ALL PROTOCOL COMBINATIONS ===\n")  
patient_protocols %>% count(protocols) %>% print()

# ---- Define clean groups ----  
itbs_only_ids <- patient_protocols %>%  
  filter(protocols == "1") %>%  
  pull(study_id)

switcher_ids <- patient_protocols %>%  
  filter(protocols == "1,2") %>%  
  filter(!study_id %in% c("TMS093", "TMS084")) %>%  
  pull(study_id)

analytic_ids <- c(itbs_only_ids, switcher_ids)

cat("\n=== ANALYTIC SAMPLE ===\n")  
cat("iTBS-only:", length(itbs_only_ids), "\n")  
cat("Switchers:", length(switcher_ids), "\n")  
cat("Total:", length(analytic_ids), "\n")

# ---- Document all excluded patients ----  
all_ids <- unique(long$study_id)  
excluded_ids <- setdiff(all_ids, analytic_ids)

excluded_df <- patient_protocols %>%  
  filter(study_id %in% excluded_ids) %>%  
  mutate(reason = case_when(  
    study_id == "TMS093" ~ "Started on protocol 2, switched to 1 (wrong direction)",  
    study_id == "TMS084" ~ "Alternated between protocols 1 and 2 (no clean switch)",  
    protocols == "4" ~ "Protocol 4 only (10Hz, older protocol)",  
    protocols == "1,4" ~ "Mixed iTBS and 10Hz protocols",  
    protocols == "1,2,4" ~ "Mixed iTBS, bilateral, and 10Hz protocols",  
    protocols == "1,2,5" ~ "Mixed iTBS, bilateral, and 10Hz/1Hz protocols",  
    protocols == "1,4,5" ~ "Mixed iTBS, 10Hz, and 10Hz/1Hz protocols",  
    TRUE ~ "Other protocol combination"  
  ))

cat("\n=== EXCLUDED PATIENTS ===\n")  
excluded_df %>% count(reason) %>% print()

# ---- Filter longitudinal data ----  
analytic_long <- long %>%  
  filter(study_id %in% analytic_ids) %>%  
  mutate(group = if_else(study_id %in% itbs_only_ids,  
                         "iTBS_only", "switcher"))

# ---- Filter summary data ----  
analytic_summary <- summary_df %>%  
  filter(study_id %in% analytic_ids) %>%  
  mutate(group = if_else(study_id %in% itbs_only_ids,  
                         "iTBS_only", "switcher"))

# ---- Merge demographics ----  
analytic_summary <- analytic_summary %>%  
  left_join(  
    demo %>% select(study_id, age_at_first_phq, sex),  
    by = "study_id"  
  )

# ---- Final verification ----  
cat("\n=== FINAL VERIFICATION ===\n")  
cat("Longitudinal rows:", nrow(analytic_long), "\n")  
cat("Unique patients:", n_distinct(analytic_long$study_id), "\n")  
cat("Summary rows:", nrow(analytic_summary), "\n")

cat("\nGroup counts:\n")  
analytic_summary %>% count(group) %>% print()

cat("\nOutcomes by group:\n")  
analytic_summary %>% count(group, outcome) %>% print()

cat("\nBaseline PHQ-9 by group:\n")  
analytic_summary %>%  
  group_by(group) %>%  
  summarise(  
    n = n(),  
    mean_baseline = round(mean(tx_baseline, na.rm = TRUE), 1),  
    sd_baseline = round(sd(tx_baseline, na.rm = TRUE), 1),  
    mean_endpoint = round(mean(endpoint, na.rm = TRUE), 1),  
    sd_endpoint = round(sd(endpoint, na.rm = TRUE), 1),  
    .groups = "drop"  
  ) %>%  
  print()

cat("\nDemographics:\n")  
cat("Missing age:", sum(is.na(analytic_summary$age_at_first_phq)), "\n")  
cat("Missing sex:", sum(is.na(analytic_summary$sex)), "\n")

# ---- Save files as .xlsx ----  
write_xlsx(analytic_long, "analytic_sample_long.xlsx")  
write_xlsx(analytic_summary, "analytic_sample_summary.xlsx")  
write_xlsx(excluded_df, "excluded_patients.xlsx")

cat("\n=== DONE ===\n")  
cat("Saved: analytic_sample_long.xlsx\n")  
cat("Saved: analytic_sample_summary.xlsx\n")  
cat("Saved: excluded_patients.xlsx\n")  