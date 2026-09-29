# ============================================================================  
# Response and Remission Rates: Recalculated for 84-Patient Sample  
# TMS074 excluded | Capped at 36 sessions | Matched switchers  
# ============================================================================

if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")

library(readxl)  
library(dplyr)  
library(effsize)

# --- 1. Load and prepare data ------------------------------------------------  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

df <- raw %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  filter(!is.na(phq9_score) & !is.na(phq9_date))

# Only keep matched switchers (both phases)  
switch_pre <- df %>% filter(group == "switcher" & phase == "pre") %>% pull(study_id) %>% unique()  
switch_post <- df %>% filter(group == "switcher" & phase == "post") %>% pull(study_id) %>% unique()  
switch_both <- intersect(switch_pre, switch_post)

df <- df %>%  
  filter(group == "iTBS_only" | study_id %in% switch_both)

# Verify sample  
n_total <- n_distinct(df$study_id)  
n_itbs <- n_distinct(df$study_id[df$group == "iTBS_only"])  
n_switch <- n_distinct(df$study_id[df$group == "switcher"])

cat("=================================================================\n")  
cat("SAMPLE VERIFICATION\n")  
cat("=================================================================\n\n")  
cat(sprintf("Total patients: %d\n", n_total))  
cat(sprintf("iTBS-Only: %d\n", n_itbs))  
cat(sprintf("Switchers: %d\n", n_switch))  
cat(sprintf("TMS074 excluded: %s\n", !("TMS074" %in% df$study_id)))  
cat(sprintf("Total PHQ-9 assessments: %d\n\n", nrow(df)))

# --- 2. Extract baseline and final PHQ-9 per patient -------------------------

# Baseline = first PHQ-9 (visit 0 or earliest date)  
baseline_scores <- df %>%  
  group_by(study_id, group) %>%  
  arrange(phq9_date) %>%  
  slice(1) %>%  
  ungroup() %>%  
  select(study_id, group, baseline_phq9 = phq9_score)

# Final = last PHQ-9 (latest date within 36-session cap)  
final_scores <- df %>%  
  group_by(study_id, group) %>%  
  arrange(desc(phq9_date)) %>%  
  slice(1) %>%  
  ungroup() %>%  
  select(study_id, group, final_phq9 = phq9_score)

# Merge  
patient_df <- baseline_scores %>%  
  left_join(final_scores, by = c("study_id", "group")) %>%  
  mutate(  
    group_label = ifelse(group == "iTBS_only", "iTBS-Only", "Switchers"),  
    phq9_change = final_phq9 - baseline_phq9,  
    pct_change = (final_phq9 - baseline_phq9) / baseline_phq9 * 100,  
    response = ifelse(pct_change <= -50, 1, 0),  
    remission = ifelse(final_phq9 < 5, 1, 0)  
  )

cat("=================================================================\n")  
cat("BASELINE AND FINAL PHQ-9 SCORES\n")  
cat("=================================================================\n\n")

for (grp in c("iTBS-Only", "Switchers")) {  
  sub <- patient_df %>% filter(group_label == grp)  
  cat(sprintf("%s (n = %d):\n", grp, nrow(sub)))  
  cat(sprintf("  Baseline PHQ-9: M = %.1f (SD = %.1f)\n",  
              mean(sub$baseline_phq9), sd(sub$baseline_phq9)))  
  cat(sprintf("  Final PHQ-9:    M = %.1f (SD = %.1f)\n",  
              mean(sub$final_phq9), sd(sub$final_phq9)))  
  cat(sprintf("  Mean change:    %.2f (SD = %.2f)\n",  
              mean(sub$phq9_change), sd(sub$phq9_change)))  
  cat(sprintf("  Mean %% change:  %.1f%%\n\n",  
              mean(sub$pct_change)))  
}

# --- 3. Response rates (>=50% reduction) --------------------------------------

cat("=================================================================\n")  
cat("RESPONSE RATES (>=50%% PHQ-9 reduction from baseline)\n")  
cat("=================================================================\n\n")

itbs_resp <- patient_df %>% filter(group_label == "iTBS-Only")  
sw_resp <- patient_df %>% filter(group_label == "Switchers")

n_resp_itbs <- sum(itbs_resp$response)  
n_resp_sw <- sum(sw_resp$response)

cat(sprintf("iTBS-Only: %d/%d (%.1f%%)\n",  
            n_resp_itbs, n_itbs, n_resp_itbs / n_itbs * 100))  
cat(sprintf("Switchers: %d/%d (%.1f%%)\n\n",  
            n_resp_sw, n_switch, n_resp_sw / n_switch * 100))

# Fisher's exact test  
resp_table <- matrix(c(n_resp_itbs, n_itbs - n_resp_itbs,  
                       n_resp_sw, n_switch - n_resp_sw),  
                     nrow = 2, byrow = TRUE,  
                     dimnames = list(c("iTBS-Only", "Switchers"),  
                                     c("Response", "No Response")))

cat("Contingency table:\n")  
print(resp_table)  
cat("\n")

fisher_resp <- fisher.test(resp_table)  
cat("Fisher's exact test:\n")  
print(fisher_resp)  
cat(sprintf("\n  OR = %.2f, 95%% CI [%.2f, %.2f], p = %.3f\n\n",  
            fisher_resp$estimate,  
            fisher_resp$conf.int[1], fisher_resp$conf.int[2],  
            fisher_resp$p.value))

# Chi-square test  
chi_resp <- chisq.test(resp_table, correct = FALSE)  
cat(sprintf("  Chi-square(1) = %.2f, p = %.3f\n\n", chi_resp$statistic, chi_resp$p.value))

# --- 4. Remission rates (final PHQ-9 < 5) ------------------------------------

cat("=================================================================\n")  
cat("REMISSION RATES (final PHQ-9 < 5)\n")  
cat("=================================================================\n\n")

n_rem_itbs <- sum(itbs_resp$remission)  
n_rem_sw <- sum(sw_resp$remission)

cat(sprintf("iTBS-Only: %d/%d (%.1f%%)\n",  
            n_rem_itbs, n_itbs, n_rem_itbs / n_itbs * 100))  
cat(sprintf("Switchers: %d/%d (%.1f%%)\n\n",  
            n_rem_sw, n_switch, n_rem_sw / n_switch * 100))

# Fisher's exact test  
rem_table <- matrix(c(n_rem_itbs, n_itbs - n_rem_itbs,  
                      n_rem_sw, n_switch - n_rem_sw),  
                    nrow = 2, byrow = TRUE,  
                    dimnames = list(c("iTBS-Only", "Switchers"),  
                                    c("Remission", "No Remission")))

cat("Contingency table:\n")  
print(rem_table)  
cat("\n")

fisher_rem <- fisher.test(rem_table)  
cat("Fisher's exact test:\n")  
print(fisher_rem)  
cat(sprintf("\n  OR = %.2f, 95%% CI [%.2f, %.2f], p = %.3f\n\n",  
            fisher_rem$estimate,  
            fisher_rem$conf.int[1], fisher_rem$conf.int[2],  
            fisher_rem$p.value))

# Chi-square test  
chi_rem <- chisq.test(rem_table, correct = FALSE)  
cat(sprintf("  Chi-square(1) = %.2f, p = %.3f\n\n", chi_rem$statistic, chi_rem$p.value))

# Absolute risk difference and NNT  
ard <- (n_rem_itbs / n_itbs) - (n_rem_sw / n_switch)  
nnt <- ifelse(ard > 0, ceiling(1 / ard), NA)  
cat(sprintf("  Absolute risk difference: %.1f%%\n", ard * 100))  
cat(sprintf("  NNT: %s\n\n", ifelse(!is.na(nnt), as.character(nnt), "N/A")))

# --- 5. Within-group paired t-tests -------------------------------------------

cat("=================================================================\n")  
cat("WITHIN-GROUP PHQ-9 CHANGE (Paired t-tests)\n")  
cat("=================================================================\n\n")

for (grp in c("iTBS-Only", "Switchers")) {  
  sub <- patient_df %>% filter(group_label == grp)
  
  t_result <- t.test(sub$baseline_phq9, sub$final_phq9, paired = TRUE)  
  d_result <- cohen.d(sub$baseline_phq9, sub$final_phq9, paired = TRUE)
  
  mean_change <- mean(sub$phq9_change)  
  sd_change <- sd(sub$phq9_change)  
  n_grp <- nrow(sub)
  
  cat(sprintf("%s (n = %d):\n", grp, n_grp))  
  cat(sprintf("  Mean change: %.2f (SD = %.2f)\n", mean_change, sd_change))  
  cat(sprintf("  Paired t(%d) = %.2f, p %s\n",  
              t_result$parameter,  
              t_result$statistic,  
              ifelse(t_result$p.value < 0.001, "< .001", sprintf("= %.3f", t_result$p.value))))  
  cat(sprintf("  Cohen's d = %.2f, 95%% CI [%.2f, %.2f]\n",  
              abs(d_result$estimate),  
              d_result$conf.int[1], d_result$conf.int[2]))  
  cat(sprintf("  Effect size: %s\n\n",  
              ifelse(abs(d_result$estimate) >= 0.8, "large",  
                     ifelse(abs(d_result$estimate) >= 0.5, "medium",  
                            ifelse(abs(d_result$estimate) >= 0.2, "small", "negligible")))))  
}

# --- 6. Between-group comparisons ---------------------------------------------

cat("=================================================================\n")  
cat("BETWEEN-GROUP COMPARISONS\n")  
cat("=================================================================\n\n")

itbs_change <- patient_df$phq9_change[patient_df$group_label == "iTBS-Only"]  
sw_change <- patient_df$phq9_change[patient_df$group_label == "Switchers"]

# PHQ-9 change comparison  
t_btw <- t.test(itbs_change, sw_change, var.equal = FALSE)  
d_btw <- cohen.d(itbs_change, sw_change)

cat("PHQ-9 raw change (iTBS-Only vs Switchers):\n")  
cat(sprintf("  iTBS-Only: M = %.2f (SD = %.2f)\n",  
            mean(itbs_change), sd(itbs_change)))  
cat(sprintf("  Switchers: M = %.2f (SD = %.2f)\n",  
            mean(sw_change), sd(sw_change)))  
cat(sprintf("  Welch's t(%.1f) = %.2f, p = %.3f\n",  
            t_btw$parameter, t_btw$statistic, t_btw$p.value))  
cat(sprintf("  Cohen's d = %.2f, 95%% CI [%.2f, %.2f], %s effect\n\n",  
            d_btw$estimate,  
            d_btw$conf.int[1], d_btw$conf.int[2],  
            ifelse(abs(d_btw$estimate) >= 0.8, "large",  
                   ifelse(abs(d_btw$estimate) >= 0.5, "medium",  
                          ifelse(abs(d_btw$estimate) >= 0.2, "small", "negligible")))))

# Percent change comparison  
itbs_pct <- patient_df$pct_change[patient_df$group_label == "iTBS-Only"]  
sw_pct <- patient_df$pct_change[patient_df$group_label == "Switchers"]

t_pct <- t.test(itbs_pct, sw_pct, var.equal = FALSE)  
d_pct <- cohen.d(itbs_pct, sw_pct)

cat("PHQ-9 percent change (iTBS-Only vs Switchers):\n")  
cat(sprintf("  iTBS-Only: M = %.1f%% (SD = %.1f%%)\n",  
            mean(itbs_pct), sd(itbs_pct)))  
cat(sprintf("  Switchers: M = %.1f%% (SD = %.1f%%)\n",  
            mean(sw_pct), sd(sw_pct)))  
cat(sprintf("  Welch's t(%.1f) = %.2f, p = %.3f\n",  
            t_pct$parameter, t_pct$statistic, t_pct$p.value))  
cat(sprintf("  Cohen's d = %.2f, %s effect\n\n",  
            d_pct$estimate,  
            ifelse(abs(d_pct$estimate) >= 0.8, "large",  
                   ifelse(abs(d_pct$estimate) >= 0.5, "medium",  
                          ifelse(abs(d_pct$estimate) >= 0.2, "small", "negligible")))))

# --- 7. Summary for manuscript ------------------------------------------------

cat("=================================================================\n")  
cat("MANUSCRIPT-READY SUMMARY\n")  
cat("=================================================================\n\n")

cat(sprintf("Sample: N = %d (%d iTBS-only, %d switchers)\n", n_total, n_itbs, n_switch))  
cat(sprintf("TMS074 excluded, sessions capped at 36, matched switchers\n\n"))

cat(sprintf("RESPONSE (>=50%% PHQ-9 reduction):\n"))  
cat(sprintf("  iTBS-Only: %d/%d (%.1f%%)\n", n_resp_itbs, n_itbs, n_resp_itbs/n_itbs*100))  
cat(sprintf("  Switchers: %d/%d (%.1f%%)\n", n_resp_sw, n_switch, n_resp_sw/n_switch*100))  
cat(sprintf("  Fisher's exact p = %.3f; OR = %.2f, 95%% CI [%.2f, %.2f]\n\n",  
            fisher_resp$p.value, fisher_resp$estimate,  
            fisher_resp$conf.int[1], fisher_resp$conf.int[2]))

cat(sprintf("REMISSION (final PHQ-9 < 5):\n"))  
cat(sprintf("  iTBS-Only: %d/%d (%.1f%%)\n", n_rem_itbs, n_itbs, n_rem_itbs/n_itbs*100))  
cat(sprintf("  Switchers: %d/%d (%.1f%%)\n", n_rem_sw, n_switch, n_rem_sw/n_switch*100))  
cat(sprintf("  Fisher's exact p = %.3f; OR = %.2f, 95%% CI [%.2f, %.2f]\n",  
            fisher_rem$p.value, fisher_rem$estimate,  
            fisher_rem$conf.int[1], fisher_rem$conf.int[2]))  
cat(sprintf("  Chi-square(1) = %.2f, p = %.3f\n", chi_rem$statistic, chi_rem$p.value))  
cat(sprintf("  ARD = %.1f%%, NNT = %s\n\n",  
            ard * 100, ifelse(!is.na(nnt), as.character(nnt), "N/A")))

cat(sprintf("WITHIN-GROUP CHANGE:\n"))  
for (grp in c("iTBS-Only", "Switchers")) {  
  sub <- patient_df %>% filter(group_label == grp)  
  t_r <- t.test(sub$baseline_phq9, sub$final_phq9, paired = TRUE)  
  d_r <- cohen.d(sub$baseline_phq9, sub$final_phq9, paired = TRUE)  
  cat(sprintf("  %s: mean change = %.2f (SD = %.2f); t(%d) = %.2f, p %s; d = %.2f\n",  
              grp, mean(sub$phq9_change), sd(sub$phq9_change),  
              t_r$parameter, t_r$statistic,  
              ifelse(t_r$p.value < 0.001, "< .001", sprintf("= %.3f", t_r$p.value)),  
              abs(d_r$estimate)))  
}

cat(sprintf("\nBETWEEN-GROUP CHANGE:\n"))  
cat(sprintf("  Raw change: t(%.1f) = %.2f, p = %.3f; d = %.2f\n",  
            t_btw$parameter, t_btw$statistic, t_btw$p.value, d_btw$estimate))  
cat(sprintf("  Pct change: t(%.1f) = %.2f, p = %.3f; d = %.2f\n",  
            t_pct$parameter, t_pct$statistic, t_pct$p.value, d_pct$estimate))

cat("\nDone.\n")  