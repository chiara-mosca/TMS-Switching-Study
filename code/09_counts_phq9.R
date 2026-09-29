library(readxl)  
library(dplyr)

raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

df <- raw %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  filter(!is.na(phq9_score) & !is.na(phq9_date))

# Only keep matched switchers  
switch_pre <- df %>% filter(group == "switcher" & phase == "pre") %>% pull(study_id) %>% unique()  
switch_post <- df %>% filter(group == "switcher" & phase == "post") %>% pull(study_id) %>% unique()  
switch_both <- intersect(switch_pre, switch_post)

df <- df %>%  
  filter(group == "iTBS_only" | study_id %in% switch_both)

# Total counts  
cat("=================================================================\n")  
cat("PHQ-9 ASSESSMENT COUNTS\n")  
cat("=================================================================\n\n")

cat(sprintf("Total patients: %d\n", n_distinct(df$study_id)))  
cat(sprintf("Total PHQ-9 assessments: %d\n\n", nrow(df)))

# By group  
itbs_df <- df %>% filter(group == "iTBS_only")  
switch_df <- df %>% filter(group == "switcher")

cat(sprintf("iTBS-Only:\n"))  
cat(sprintf("  Patients: %d\n", n_distinct(itbs_df$study_id)))  
cat(sprintf("  PHQ-9 assessments: %d\n", nrow(itbs_df)))  
cat(sprintf("  Mean assessments per patient: %.1f\n", nrow(itbs_df) / n_distinct(itbs_df$study_id)))  
cat(sprintf("  Range: %d to %d\n\n",  
            min(table(itbs_df$study_id)),  
            max(table(itbs_df$study_id))))

cat(sprintf("Switchers:\n"))  
cat(sprintf("  Patients: %d\n", n_distinct(switch_df$study_id)))  
cat(sprintf("  PHQ-9 assessments: %d\n", nrow(switch_df)))  
cat(sprintf("  Mean assessments per patient: %.1f\n", nrow(switch_df) / n_distinct(switch_df$study_id)))  
cat(sprintf("  Range: %d to %d\n\n",  
            min(table(switch_df$study_id)),  
            max(table(switch_df$study_id))))

# Switchers broken down by phase  
switch_pre_df <- switch_df %>% filter(phase == "pre")  
switch_post_df <- switch_df %>% filter(phase == "post")

cat(sprintf("Switchers Phase 1 (pre):\n"))  
cat(sprintf("  Assessments: %d\n", nrow(switch_pre_df)))  
cat(sprintf("  Mean per patient: %.1f\n\n", nrow(switch_pre_df) / n_distinct(switch_pre_df$study_id)))

cat(sprintf("Switchers Phase 2 (post):\n"))  
cat(sprintf("  Assessments: %d\n", nrow(switch_post_df)))  
cat(sprintf("  Mean per patient: %.1f\n\n", nrow(switch_post_df) / n_distinct(switch_post_df$study_id)))

cat(sprintf("TOTAL: %d + %d = %d assessments across %d patients\n",  
            nrow(itbs_df), nrow(switch_df),  
            nrow(itbs_df) + nrow(switch_df),  
            n_distinct(df$study_id)))  