# ============================================================================  
# Item-Level HLM Models for PHQ-9 Items  
# Report: Group × Time Interaction for Each PHQ-9 Item  
# ============================================================================

# Load required libraries  
library(readxl)  
library(tidyverse)  
library(lme4)  
library(lmerTest)  
library(broom.mixed)  
library(knitr)  
library(kableExtra)

# ----------------------------------------------------------------------------  
# 1. DATA LOADING AND PREPARATION  
# ----------------------------------------------------------------------------

# Read the data  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

# Work with a copy  
df <- raw

# Inspect column names  
names(df)

# Rename PHQ-9 item columns for clarity  
item_labels <- c(  
  "phq9_1" = "Item 1: Anhedonia",  
  "phq9_2" = "Item 2: Depressed Mood",  
  "phq9_3" = "Item 3: Sleep",  
  "phq9_4" = "Item 4: Fatigue",  
  "phq9_5" = "Item 5: Appetite",  
  "phq9_6" = "Item 6: Guilt/Worthlessness",  
  "phq9_7" = "Item 7: Concentration",  
  "phq9_8" = "Item 8: Psychomotor",  
  "phq9_9" = "Item 9: Suicidal Ideation"  
)

# Define the two groups for comparison  
# "switcher" = patients who switched to bilateral iTBS/cTBS  
# "iTBS_only" = patients who stayed on left-sided iTBS only  
df <- df %>%  
  filter(group %in% c("switcher", "iTBS_only")) %>%  
  mutate(  
    group = factor(group, levels = c("iTBS_only", "switcher"),  
                   labels = c("iTBS Only", "Switcher")),  
    study_id = factor(study_id),  
    weeks = as.numeric(weeks)  
  )

# Convert PHQ-9 item columns to numeric  
item_cols <- paste0("phq9_", 1:9)  
df <- df %>%  
  mutate(across(all_of(item_cols), ~ as.numeric(.)))

# Check data structure  
str(df[, c("study_id", "group", "weeks", item_cols)])

# Summary of sample sizes  
cat("\n=== Sample Size Summary ===\n")  
df %>%  
  group_by(group) %>%  
  summarise(  
    n_patients = n_distinct(study_id),  
    n_observations = n(),  
    .groups = "drop"  
  ) %>%  
  print()

# ----------------------------------------------------------------------------  
# 2. APPROACH A: SEPARATE HLM FOR EACH PHQ-9 ITEM  
# ----------------------------------------------------------------------------

cat("\n============================================================\n")  
cat("APPROACH A: Separate HLM for Each PHQ-9 Item\n")  
cat("Model: item_score ~ weeks * group + (1 + weeks | study_id)\n")  
cat("============================================================\n\n")

# Function to run HLM for a single item  
run_item_hlm <- function(data, item_var, item_label) {  
  
  # Create a working dataset with non-missing values for this item  
  d <- data %>%  
    select(study_id, group, weeks, value = all_of(item_var)) %>%  
    filter(!is.na(value))  
  
  # Fit the model with random intercept and slope  
  # If convergence issues, fall back to random intercept only  
  model <- tryCatch({  
    lmer(value ~ weeks * group + (1 + weeks | study_id), data = d,  
         REML = TRUE,  
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  
  }, warning = function(w) {  
    message(paste0("  Warning for ", item_label, ": ", w$message))  
    message("  Attempting random intercept only model...")  
    lmer(value ~ weeks * group + (1 | study_id), data = d,  
         REML = TRUE,  
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  
  }, error = function(e) {  
    message(paste0("  Error for ", item_label, ": ", e$message))  
    message("  Fitting random intercept only model...")  
    lmer(value ~ weeks * group + (1 | study_id), data = d,  
         REML = TRUE,  
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  
  })  
  
  return(model)  
}

# Run models for all 9 items  
item_models <- list()  
item_results <- list()

for (i in 1:9) {  
  item_var <- paste0("phq9_", i)  
  item_label <- item_labels[item_var]  
  
  cat(paste0("\n--- Fitting model for ", item_label, " ---\n"))  
  
  # Fit model  
  item_models[[item_var]] <- run_item_hlm(df, item_var, item_label)  
  
  # Extract fixed effects  
  fe <- summary(item_models[[item_var]])$coefficients  
  
  # Print full model summary  
  cat("\nFull Model Summary:\n")  
  print(summary(item_models[[item_var]]))  
  
  # Extract the Group × Time interaction (weeks:groupSwitcher)  
  interaction_term <- "weeks:groupSwitcher"  
  
  if (interaction_term %in% rownames(fe)) {  
    beta <- fe[interaction_term, "Estimate"]  
    se <- fe[interaction_term, "Std. Error"]  
    t_val <- fe[interaction_term, "t value"]  
    df_val <- fe[interaction_term, "df"]  
    p_val <- fe[interaction_term, "Pr(>|t|)"]  
    ci_lower <- beta - 1.96 * se  
    ci_upper <- beta + 1.96 * se  
    
    item_results[[item_var]] <- tibble(  
      Item = item_label,  
      Beta = beta,  
      SE = se,  
      CI_Lower = ci_lower,  
      CI_Upper = ci_upper,  
      t = t_val,  
      df = df_val,  
      p = p_val  
    )  
  }  
}

# Combine results into a single table  
results_table <- bind_rows(item_results)

# Add significance indicators  
results_table <- results_table %>%  
  mutate(  
    Sig = case_when(  
      p < 0.001 ~ "***",  
      p < 0.01  ~ "**",  
      p < 0.05  ~ "*",  
      p < 0.10  ~ "†",  
      TRUE      ~ ""  
    ),  
    # Format columns for display  
    Beta_fmt = round(Beta, 4),  
    SE_fmt = round(SE, 4),  
    CI_fmt = paste0("[", round(CI_Lower, 4), ", ", round(CI_Upper, 4), "]"),  
    t_fmt = round(t, 3),  
    df_fmt = round(df, 1),  
    p_fmt = ifelse(p < 0.001, "<.001", round(p, 4))  
  )

# ----------------------------------------------------------------------------  
# 3. DISPLAY RESULTS TABLE  
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")  
cat("TABLE: Group × Time Interaction for Each PHQ-9 Item\n")  
cat("Reference group: iTBS Only; Comparison: Switcher\n")  
cat("============================================================\n\n")

# Print formatted table  
display_table <- results_table %>%  
  select(Item, `β` = Beta_fmt, `SE` = SE_fmt, `95% CI` = CI_fmt,   
         `t` = t_fmt, `df` = df_fmt, `p` = p_fmt, ` ` = Sig) %>%  
  arrange(Item)

print(kable(display_table,   
            format = "pipe",  
            align = c("l", "r", "r", "c", "r", "r", "r", "l"),  
            caption = "Group × Time (Weeks) Interaction Effects by PHQ-9 Item"))

# Also print as a clean console table  
cat("\n\nClean Results:\n")  
cat(paste(rep("-", 100), collapse = ""), "\n")  
cat(sprintf("%-35s %8s %8s %22s %8s %8s %8s %4s\n",  
            "Item", "β", "SE", "95% CI", "t", "df", "p", ""))  
cat(paste(rep("-", 100), collapse = ""), "\n")

for (j in 1:nrow(results_table)) {  
  r <- results_table[j, ]  
  cat(sprintf("%-35s %8.4f %8.4f %22s %8.3f %8.1f %8s %4s\n",  
              r$Item, r$Beta, r$SE, r$CI_fmt, r$t, r$df, r$p_fmt, r$Sig))  
}  
cat(paste(rep("-", 100), collapse = ""), "\n")  
cat("Note: † p < .10, * p < .05, ** p < .01, *** p < .001\n")  
cat("Positive β indicates greater increase (or less decrease) for Switcher group.\n")  
cat("Negative β indicates greater decrease for Switcher group relative to iTBS Only.\n\n")

# ----------------------------------------------------------------------------  
# 4. APPROACH B: SINGLE LONG-FORMAT MODEL WITH ITEM AS FACTOR  
# ----------------------------------------------------------------------------

cat("\n============================================================\n")  
cat("APPROACH B: Single Model with Item as Factor (Long Format)\n")  
cat("============================================================\n\n")

# Reshape data to long format (one row per item per observation)  
df_long <- df %>%  
  select(study_id, group, weeks, visit, all_of(item_cols)) %>%  
  pivot_longer(  
    cols = all_of(item_cols),  
    names_to = "item",  
    values_to = "score"  
  ) %>%  
  filter(!is.na(score)) %>%  
  mutate(  
    item = factor(item, levels = item_cols,  
                  labels = paste0("Item ", 1:9)),  
    obs_id = paste0(study_id, "_", visit)  
  )

# Check structure  
cat("Long-format data dimensions:", nrow(df_long), "rows\n")  
cat("Number of unique patients:", n_distinct(df_long$study_id), "\n\n")

# Fit the full model with three-way interaction: weeks × group × item  
cat("Fitting full three-way interaction model...\n")  
cat("Model: score ~ weeks * group * item + (1 + weeks | study_id)\n\n")

model_full <- tryCatch({  
  lmer(score ~ weeks * group * item + (1 + weeks | study_id),  
       data = df_long,  
       REML = TRUE,  
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e5)))  
}, warning = function(w) {  
  message("Warning: ", w$message)  
  lmer(score ~ weeks * group * item + (1 | study_id),  
       data = df_long,  
       REML = TRUE,  
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e5)))  
})

# Print model summary  
cat("\nFull Three-Way Interaction Model Summary:\n")  
print(summary(model_full))

# Extract the three-way interaction terms (weeks:group:item)  
fe_full <- as.data.frame(summary(model_full)$coefficients)  
fe_full$term <- rownames(fe_full)

three_way <- fe_full %>%  
  filter(grepl("weeks:group.*:item|weeks:item.*:group", term)) %>%  
  mutate(  
    CI_Lower = Estimate - 1.96 * `Std. Error`,  
    CI_Upper = Estimate + 1.96 * `Std. Error`,  
    Sig = case_when(  
      `Pr(>|t|)` < 0.001 ~ "***",  
      `Pr(>|t|)` < 0.01  ~ "**",  
      `Pr(>|t|)` < 0.05  ~ "*",  
      `Pr(>|t|)` < 0.10  ~ "†",  
      TRUE               ~ ""  
    )  
  )

cat("\n\nThree-Way Interaction Terms (Weeks × Group × Item):\n")  
cat("These test whether the Group × Time interaction DIFFERS across items\n")  
cat("(relative to the reference item, Item 1)\n\n")  
print(three_way %>%  
        select(term, Estimate, `Std. Error`, CI_Lower, CI_Upper,   
               `t value`, df, `Pr(>|t|)`, Sig))

# ----------------------------------------------------------------------------  
# 5. FOCUSED ANALYSIS ON ITEM 1 (ANHEDONIA)  
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")  
cat("FOCUSED ANALYSIS: Item 1 (Anhedonia)\n")  
cat("============================================================\n\n")

item1_model <- item_models[["phq9_1"]]  
cat("Full Model Summary for Item 1 (Anhedonia):\n")  
print(summary(item1_model))

# Compare Item 1 to all other items  
cat("\n\nComparison of Group × Time Interaction Across All Items:\n")  
cat("(Highlighting Item 1: Anhedonia)\n\n")

results_sorted <- results_table %>%  
  arrange(p) %>%  
  mutate(  
    Rank = row_number(),  
    Highlight = ifelse(grepl("Anhedonia", Item), ">>> ", "    ")  
  )

cat(sprintf("%4s %-4s %-35s %8s %8s %22s %8s %8s %4s\n",  
            "Rank", "", "Item", "β", "SE", "95% CI", "t", "p", ""))  
cat(paste(rep("-", 105), collapse = ""), "\n")  
for (j in 1:nrow(results_sorted)) {  
  r <- results_sorted[j, ]  
  cat(sprintf("%4d %4s %-35s %8.4f %8.4f %22s %8.3f %8s %4s\n",  
              r$Rank, r$Highlight, r$Item, r$Beta, r$SE,   
              r$CI_fmt, r$t, r$p_fmt, r$Sig))  
}  
cat(paste(rep("-", 105), collapse = ""), "\n\n")

# ----------------------------------------------------------------------------  
# 6. FDR CORRECTION FOR MULTIPLE COMPARISONS  
# ----------------------------------------------------------------------------

cat("\n============================================================\n")  
cat("MULTIPLE COMPARISON CORRECTION (FDR - Benjamini-Hochberg)\n")  
cat("============================================================\n\n")

results_fdr <- results_table %>%  
  mutate(  
    p_numeric = p,  
    p_fdr = p.adjust(p_numeric, method = "BH"),  
    Sig_FDR = case_when(  
      p_fdr < 0.001 ~ "***",  
      p_fdr < 0.01  ~ "**",  
      p_fdr < 0.05  ~ "*",  
      p_fdr < 0.10  ~ "†",  
      TRUE          ~ ""  
    )  
  ) %>%  
  arrange(p_numeric)

cat(sprintf("%-35s %10s %10s %8s %8s\n",  
            "Item", "p (raw)", "p (FDR)", "Sig_raw", "Sig_FDR"))  
cat(paste(rep("-", 80), collapse = ""), "\n")  
for (j in 1:nrow(results_fdr)) {  
  r <- results_fdr[j, ]  
  cat(sprintf("%-35s %10s %10.4f %8s %8s\n",  
              r$Item, r$p_fmt, r$p_fdr, r$Sig, r$Sig_FDR))  
}  
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# ----------------------------------------------------------------------------  
# 7. VISUALIZATION: FOREST PLOT OF GROUP × TIME INTERACTIONS  
# ----------------------------------------------------------------------------

cat("Creating forest plot of Group × Time interactions...\n\n")

forest_data <- results_table %>%  
  mutate(  
    Item_short = gsub("Item [0-9]: ", "", Item),  
    Item_num = as.numeric(gsub("Item ([0-9]):.*", "\\1", Item)),  
    significant = ifelse(p < 0.05, "Significant", "Not Significant"),  
    is_anhedonia = ifelse(grepl("Anhedonia", Item), TRUE, FALSE)  
  ) %>%  
  arrange(desc(Item_num))

p_forest <- ggplot(forest_data, aes(x = Beta, y = reorder(Item, -Item_num))) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +  
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper,   
                     color = significant),  
                 height = 0.25, linewidth = 0.8) +  
  geom_point(aes(color = significant, shape = is_anhedonia),   
             size = 3.5) +  
  scale_color_manual(values = c("Significant" = "#E41A1C",   
                                "Not Significant" = "#377EB8"),  
                     name = "Statistical\nSignificance") +  
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),  
                     guide = "none") +  
  labs(  
    title = "Group × Time Interaction by PHQ-9 Item",  
    subtitle = "Comparing Switcher vs. iTBS Only trajectories over time",  
    x = expression(paste("Interaction Coefficient (", beta, ")")),  
    y = "",  
    caption = "Note: Diamond symbol highlights Item 1 (Anhedonia).\nError bars represent 95% CI. Red = p < .05."  
  ) +  
  theme_minimal(base_size = 12) +  
  theme(  
    plot.title = element_text(face = "bold", size = 14),  
    plot.subtitle = element_text(size = 11, color = "gray40"),  
    panel.grid.minor = element_blank(),  
    legend.position = "bottom"  
  )

print(p_forest)  
ggsave("forest_plot_item_interactions.png", p_forest,   
       width = 10, height = 7, dpi = 300)

# ----------------------------------------------------------------------------  
# 8. SPAGHETTI PLOTS FOR ITEM 1 (ANHEDONIA) BY GROUP  
# ----------------------------------------------------------------------------

cat("Creating trajectory plot for Item 1 (Anhedonia)...\n\n")

d_item1 <- df %>%  
  select(study_id, group, weeks, phq9_1) %>%  
  filter(!is.na(phq9_1))

d_item1$predicted <- predict(item_models[["phq9_1"]], newdata = d_item1)

group_means <- d_item1 %>%  
  group_by(group, weeks) %>%  
  summarise(  
    mean_score = mean(phq9_1, na.rm = TRUE),  
    se_score = sd(phq9_1, na.rm = TRUE) / sqrt(n()),  
    mean_predicted = mean(predicted, na.rm = TRUE),  
    n = n(),  
    .groups = "drop"  
  )

p_item1 <- ggplot() +  
  geom_line(data = d_item1,   
            aes(x = weeks, y = phq9_1, group = study_id, color = group),  
            alpha = 0.15, linewidth = 0.3) +  
  geom_line(data = group_means,   
            aes(x = weeks, y = mean_score, color = group),  
            linewidth = 1.5) +  
  geom_point(data = group_means,   
             aes(x = weeks, y = mean_score, color = group),  
             size = 2.5) +  
  geom_smooth(data = d_item1,  
              aes(x = weeks, y = predicted, color = group),  
              method = "lm", se = TRUE, linetype = "dashed",  
              linewidth = 1, alpha = 0.2) +  
  scale_color_manual(values = c("iTBS Only" = "#2166AC",   
                                "Switcher" = "#B2182B"),  
                     name = "Group") +  
  labs(  
    title = "Item 1 (Anhedonia): Differential Trajectory by Group",  
    subtitle = "Individual patient trajectories with group means",  
    x = "Weeks from Baseline",  
    y = "PHQ-9 Item 1 Score (0-3)",  
    caption = "Solid lines = group means; Dashed lines = model-fitted trends"  
  ) +  
  scale_y_continuous(limits = c(0, 3), breaks = 0:3) +  
  theme_minimal(base_size = 12) +  
  theme(  
    plot.title = element_text(face = "bold", size = 14),  
    legend.position = "bottom"  
  )

print(p_item1)  
ggsave("item1_anhedonia_trajectories.png", p_item1,   
       width = 10, height = 7, dpi = 300)

# ----------------------------------------------------------------------------  
# 9. SUPPLEMENTARY: ALL ITEMS TRAJECTORY PANEL  
# ----------------------------------------------------------------------------

cat("Creating panel plot for all PHQ-9 items...\n\n")

df_plot_long <- df %>%  
  select(study_id, group, weeks, all_of(item_cols)) %>%  
  pivot_longer(cols = all_of(item_cols),  
               names_to = "item", values_to = "score") %>%  
  filter(!is.na(score)) %>%  
  mutate(  
    item_label = factor(item, levels = item_cols,  
                        labels = names(item_labels))  
  )

group_means_all <- df_plot_long %>%  
  group_by(group, weeks, item_label) %>%  
  summarise(mean_score = mean(score, na.rm = TRUE),  
            .groups = "drop")

p_panel <- ggplot(group_means_all,   
                  aes(x = weeks, y = mean_score, color = group)) +  
  geom_line(linewidth = 1) +  
  geom_point(size = 1.5) +  
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "dashed") +  
  facet_wrap(~ item_label, ncol = 3, scales = "free_y") +  
  scale_color_manual(values = c("iTBS Only" = "#2166AC",   
                                "Switcher" = "#B2182B"),  
                     name = "Group") +  
  labs(  
    title = "PHQ-9 Item-Level Trajectories by Treatment Group",  
    x = "Weeks from Baseline",  
    y = "Mean Item Score (0-3)"  
  ) +  
  theme_minimal(base_size = 10) +  
  theme(  
    plot.title = element_text(face = "bold", size = 14),  
    strip.text = element_text(face = "bold", size = 8),  
    legend.position = "bottom"  
  )

print(p_panel)  
ggsave("all_items_panel_trajectories.png", p_panel,   
       width = 14, height = 10, dpi = 300)

# ----------------------------------------------------------------------------  
# 10. PUBLICATION-READY TABLE (HTML/LaTeX)  
# ----------------------------------------------------------------------------

cat("\n\nGenerating publication-ready table...\n\n")

pub_table <- results_table %>%  
  arrange(as.numeric(gsub("Item ([0-9]):.*", "\\1", Item))) %>%  
  select(  
    `PHQ-9 Item` = Item,  
    `β` = Beta_fmt,  
    `SE` = SE_fmt,  
    `95% CI` = CI_fmt,  
    `t` = t_fmt,  
    `p` = p_fmt,  
    ` ` = Sig  
  )

html_table <- kable(pub_table, format = "html", align = c("l", rep("r", 5), "l"),  
                    caption = "Table X. Group × Time Interaction Effects for Each PHQ-9 Item") %>%  
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),  
                full_width = FALSE, font_size = 12) %>%  
  row_spec(1, bold = TRUE, background = "#FFFFCC") %>%  
  footnote(  
    general = c(  
      "Models: Item Score ~ Weeks × Group + (1 + Weeks | Patient ID)",  
      "Reference group: iTBS Only. Positive β indicates greater increase for Switcher group.",  
      "† p < .10, * p < .05, ** p < .01, *** p < .001"  
    ),  
    general_title = "Note: "  
  )

print(html_table)  
save_kable(html_table, file = "item_level_hlm_table.html")

# ----------------------------------------------------------------------------  
# 11. FINAL SUMMARY  
# ----------------------------------------------------------------------------

cat("\n\n============================================================\n")  
cat("SUMMARY OF FINDINGS\n")  
cat("============================================================\n\n")

sig_items <- results_table %>%  
  filter(p < 0.05)

cat("Number of PHQ-9 items with significant Group × Time interaction: ",   
    nrow(sig_items), "\n\n")

if (nrow(sig_items) > 0) {  
  cat("Significant items:\n")  
  for (j in 1:nrow(sig_items)) {  
    cat(sprintf("  %s: β = %.4f, p = %s %s\n",  
                sig_items$Item[j], sig_items$Beta[j],   
                sig_items$p_fmt[j], sig_items$Sig[j]))  
  }  
}

cat("\n\nItem 1 (Anhedonia) Results:\n")  
item1_row <- results_table %>% filter(grepl("Anhedonia", Item))  
if (nrow(item1_row) > 0) {  
  cat(sprintf("  β = %.4f\n", item1_row$Beta))  
  cat(sprintf("  SE = %.4f\n", item1_row$SE))  
  cat(sprintf("  95%% CI = %s\n", item1_row$CI_fmt))  
  cat(sprintf("  t = %.3f\n", item1_row$t))  
  cat(sprintf("  p = %s %s\n", item1_row$p_fmt, item1_row$Sig))  
  
  if (item1_row$p < 0.05) {  
    cat("\n  >> Item 1 (Anhedonia) shows a STATISTICALLY SIGNIFICANT\n")  
    cat("     differential trajectory between groups.\n")  
  } else {  
    cat("\n  >> Item 1 (Anhedonia) does NOT show a statistically significant\n")  
    cat("     differential trajectory at p < .05.\n")  
  }  
}

cat("\n\nAll output files saved:\n")  
cat("  - forest_plot_item_interactions.png\n")  
cat("  - item1_anhedonia_trajectories.png\n")  
cat("  - all_items_panel_trajectories.png\n")  
cat("  - item_level_hlm_table.html\n")  
cat("\n============================================================\n")  
cat("Analysis complete.\n")  
cat("============================================================\n")  