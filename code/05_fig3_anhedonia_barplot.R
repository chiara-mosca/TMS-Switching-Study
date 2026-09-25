# ============================================================================  
# PHQ-9 Item-Level Analysis: Group x Time Interaction  
# Pre-switch period only | Capped at 36 sessions | TMS074 excluded  
# Time measured in cumulative treatment sessions (dose-response)  
# ============================================================================

library(readxl)  
library(tidyverse)  
library(lme4)  
library(lmerTest)

# --- 1. Load & Prepare ---  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

df <- raw %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  filter(!is.na(phq9_score) & !is.na(phq9_date)) %>%  
  mutate(  
    cumulative_tx = as.numeric(cumulative_tx),  
    group = case_when(  
      group == "iTBS_only" ~ "iTBS Only",  
      group == "switcher"  ~ "Switchers",  
      TRUE ~ as.character(group)  
    )  
  ) %>%  
  filter(group %in% c("iTBS Only", "Switchers"))

# Pre-switch only  
df_pre <- df %>%  
  filter(  
    (group == "iTBS Only") |  
      (group == "Switchers" & phase == "pre")  
  )

# Only keep switchers with both phases (matched to 3-panel)  
switch_ids_pre <- df %>%  
  filter(group == "Switchers" & phase == "pre") %>%  
  pull(study_id) %>% unique()  
switch_ids_post <- df %>%  
  filter(group == "Switchers" & phase == "post") %>%  
  pull(study_id) %>% unique()  
switch_ids_both <- intersect(switch_ids_pre, switch_ids_post)

df_pre <- df_pre %>%  
  filter(group == "iTBS Only" | study_id %in% switch_ids_both)

n_itbs <- n_distinct(df_pre$study_id[df_pre$group == "iTBS Only"])  
n_switch <- n_distinct(df_pre$study_id[df_pre$group == "Switchers"])

cat(sprintf("iTBS Only: %d patients\n", n_itbs))  
cat(sprintf("Switchers: %d patients\n", n_switch))  
cat(sprintf("TMS074 excluded: %s\n", !("TMS074" %in% df_pre$study_id)))  
cat(sprintf("Time variable: cumulative_tx (sessions)\n"))  
cat(sprintf("Session range: %d to %d\n\n",  
            min(df_pre$cumulative_tx), max(df_pre$cumulative_tx)))

# --- 2. Run HLMs with cumulative_tx ---  
phq_items <- paste0("phq9_", 1:9)  
item_labels <- c("Anhedonia", "Depressed Mood", "Sleep", "Fatigue",  
                 "Appetite", "Guilt/Worthlessness", "Concentration",  
                 "Psychomotor", "Suicidal Ideation")

results_df <- tibble(  
  item = item_labels,  
  phq_col = phq_items,  
  interaction_estimate = NA_real_,  
  interaction_se = NA_real_,  
  interaction_t = NA_real_,  
  interaction_p = NA_real_  
)

for (i in seq_along(phq_items)) {  
  df_temp <- df_pre %>%  
    mutate(y = as.numeric(.data[[phq_items[i]]])) %>%  
    filter(!is.na(y), !is.na(cumulative_tx))
  
  tryCatch({  
    mod <- lmer(y ~ group * cumulative_tx + (1 | study_id), data = df_temp)  
    coefs <- summary(mod)$coefficients  
    interaction_row <- grep(":", rownames(coefs))  
    results_df$interaction_estimate[i] <- coefs[interaction_row, "Estimate"]  
    results_df$interaction_se[i] <- coefs[interaction_row, "Std. Error"]  
    results_df$interaction_t[i] <- coefs[interaction_row, "t value"]  
    results_df$interaction_p[i] <- coefs[interaction_row, "Pr(>|t|)"]  
  }, error = function(e) {  
    cat("Error for", item_labels[i], ":", e$message, "\n")  
  })  
}

# --- 3. FDR correction ---  
results_df <- results_df %>%  
  mutate(  
    p_fdr = p.adjust(interaction_p, method = "fdr"),  
    significance = factor(  
      ifelse(p_fdr < 0.05, "FDR-corrected (p < .05)", "Not significant"),  
      levels = c("FDR-corrected (p < .05)", "Not significant")  
    )  
  )

# Print results  
cat("\n=== Item-Level Results (Group x Session Interaction) ===\n\n")  
results_df %>%  
  arrange(desc(interaction_estimate)) %>%  
  mutate(across(c(interaction_estimate, interaction_se, interaction_t), ~round(., 4)),  
         p_raw = sprintf("%.4f", interaction_p),  
         p_fdr_fmt = sprintf("%.4f", p_fdr)) %>%  
  select(item, Beta = interaction_estimate, SE = interaction_se,  
         t = interaction_t, p_raw, p_fdr = p_fdr_fmt, significance) %>%  
  print(n = 9)

n_sig <- sum(results_df$p_fdr < 0.05)  
sig_items <- results_df %>% filter(p_fdr < 0.05) %>% pull(item)  
cat(sprintf("\nFDR-significant items: %d (%s)\n",  
            n_sig,  
            ifelse(n_sig > 0, paste(sig_items, collapse = ", "), "none")))

# --- 4. Plot ---  
results_df2 <- results_df %>%  
  mutate(item = fct_reorder(item, interaction_estimate, .desc = TRUE))

sig_colors <- c(  
  "FDR-corrected (p < .05)" = "#B83030",  
  "Not significant"          = "#BFBFBF"  
)

y_min <- min(results_df2$interaction_estimate, 0) - 0.008  
y_max <- max(results_df2$interaction_estimate) + 0.012

p <- ggplot(results_df2,  
            aes(x = item, y = interaction_estimate, fill = significance)) +  
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50",  
             linewidth = 0.5) +  
  geom_text(aes(label = sprintf("%.3f", interaction_estimate),  
                vjust = ifelse(interaction_estimate >= 0, -0.5, 1.5)),  
            size = 3.2, fontface = "bold", color = "gray20") +  
  scale_fill_manual(values = sig_colors, name = NULL, drop = FALSE) +  
  scale_y_continuous(  
    labels = scales::number_format(accuracy = 0.01),  
    limits = c(y_min, y_max),  
    expand = expansion(mult = c(0.02, 0.05))  
  ) +  
  labs(  
    x = NULL,  
    y = "Group x Session Interaction Coefficient",  
    title = "PHQ-9 Item-Level Differential Response to iTBS"  
  ) +  
  theme_classic(base_size = 13) +  
  theme(  
    axis.text.x   = element_text(size = 10, face = "bold", angle = 45,  
                                 hjust = 1, color = "black"),  
    axis.text.y   = element_text(size = 10, color = "black"),  
    axis.title.y  = element_text(face = "bold", size = 11,  
                                 margin = margin(r = 10)),  
    axis.line     = element_line(color = "black", linewidth = 0.5),  
    axis.ticks    = element_line(color = "black", linewidth = 0.4),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5,  
                                 margin = margin(b = 12)),  
    legend.position = "bottom",  
    legend.text     = element_text(size = 10),  
    plot.margin     = margin(12, 18, 12, 12)  
  )

print(p)

ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_item_analysis.png",  
       p, width = 10, height = 6, dpi = 300, bg = "white")  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_item_analysis.pdf",  
       p, width = 10, height = 6, bg = "white")

cat("\nSaved: phq9_item_analysis.png\n")  
cat("Saved: phq9_item_analysis.pdf\n")

# --- 5. Full results summary ---  
cat("\n=================================================================\n")  
cat("FULL RESULTS SUMMARY\n")  
cat("=================================================================\n\n")  
cat(sprintf("Sample: %d iTBS-only, %d switchers (pre-switch only)\n",  
            n_itbs, n_switch))  
cat(sprintf("Sessions capped at: 36\n"))  
cat(sprintf("TMS074 excluded: YES\n"))  
cat(sprintf("Time variable: cumulative_tx (sessions)\n"))  
cat(sprintf("Correction: FDR (Benjamini-Hochberg)\n\n"))

cat(sprintf("%-25s %8s %7s %8s %10s %10s %s\n",  
            "Item", "Beta", "SE", "t", "p (raw)", "p (FDR)", "Sig"))  
cat(paste(rep("-", 85), collapse = ""), "\n")

for (i in 1:nrow(results_df)) {  
  row <- results_df[i, ]  
  cat(sprintf("%-25s %8.4f %7.4f %8.3f %10.4f %10.4f %s\n",  
              row$item, row$interaction_estimate, row$interaction_se,  
              row$interaction_t, row$interaction_p, row$p_fdr,  
              ifelse(row$p_fdr < 0.05, "*", "")))  
}  
cat(paste(rep("-", 85), collapse = ""), "\n")

cat("\nDone.\n")  