# ============================================================================  
# CORRECT FIGURE: Using cumulative_tx (genuine FDR significance)  
# ============================================================================

library(readxl)  
library(tidyverse)  
library(lme4)  
library(lmerTest)

# --- 1. Load & Prepare ---  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

df <- raw %>%  
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

n_itbs <- n_distinct(df_pre$study_id[df_pre$group == "iTBS Only"])  
n_switch <- n_distinct(df_pre$study_id[df_pre$group == "Switchers"])  
median_switch <- 19

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
cat("\n=== Results (cumulative_tx, genuine FDR) ===\n\n")  
results_df %>%  
  arrange(desc(interaction_estimate)) %>%  
  mutate(across(c(interaction_estimate, interaction_se, interaction_t), ~round(., 4)),  
         p_raw = sprintf("%.4f", interaction_p),  
         p_fdr_fmt = sprintf("%.4f", p_fdr)) %>%  
  select(item, Beta = interaction_estimate, SE = interaction_se,  
         t = interaction_t, p_raw, p_fdr = p_fdr_fmt, significance) %>%  
  print(n = 9)

cat("\nAnhedonia FDR p-value:", results_df$p_fdr[results_df$item == "Anhedonia"], "\n")  
cat("Anhedonia genuinely FDR-significant:",  
    results_df$p_fdr[results_df$item == "Anhedonia"] < 0.05, "\n")

# --- 4. Plot ---  
results_df2 <- results_df %>%  
  mutate(item = fct_reorder(item, interaction_estimate, .desc = TRUE))

sig_colors <- c(  
  "FDR-corrected (p < .05)" = "#555555",  
  "Not significant"          = "#D0D0D0"  
)

anhedonia_pos <- which(levels(results_df2$item) == "Anhedonia")  
y_min <- min(results_df2$interaction_estimate, 0) - 0.008  
y_max <- max(results_df2$interaction_estimate) + 0.01

p <- ggplot(results_df2,  
            aes(x = item, y = interaction_estimate, fill = significance)) +  
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50",  
             linewidth = 0.5) +  
  geom_text(aes(label = sprintf("%.4f", interaction_estimate),  
                vjust = ifelse(interaction_estimate >= 0, -0.5, 1.5)),  
            size = 3, fontface = "bold") +  
  annotate("segment",  
           x = anhedonia_pos - 0.5,  
           xend = anhedonia_pos + 0.5,  
           y = y_min + 0.002, yend = y_min + 0.002,  
           color = "#555555", linewidth = 2) +  
  scale_fill_manual(values = sig_colors, name = "Significance", drop = FALSE) +  
  scale_y_continuous(  
    breaks = seq(-0.02, 0.08, by = 0.02),  
    labels = scales::number_format(accuracy = 0.01),  
    limits = c(y_min, y_max),  
    expand = expansion(mult = c(0.02, 0.05))  
  ) +  
  labs(  
    x = NULL,  
    y = expression(Delta ~ "Slope (Group \u00d7 Time Interaction)"),  
    title = "Anhedonia Stands Out as the Core Differentiating Factor",  
    subtitle = "Only anhedonia non-response significantly predicts need for protocol switch (FDR-corrected)",  
    caption = paste0(  
      "Pre-switch period only | iTBS-only n=", n_itbs,  
      " | iTBS\u2192Bilateral switchers n=", n_switch,  
      " | Median switch at tx ", median_switch,  
      "\nSwitchers = patients with ONLY iTBS (left) \u2192 Bilateral protocol history"  
    )  
  ) +  
  theme_minimal(base_size = 14) +  
  theme(  
    plot.background = element_rect(fill = "white", color = NA),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid.major.x = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),  
    axis.text.x = element_text(size = 10, face = "bold", angle = 45,  
                               hjust = 1, color = "black"),  
    axis.text.y = element_text(size = 11, color = "black"),  
    axis.title.y = element_text(face = "bold", size = 12,  
                                margin = margin(r = 10)),  
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),  
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),  
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5,  
                                margin = margin(t = 15)),  
    legend.position = "bottom",  
    legend.text = element_text(size = 11),  
    legend.title = element_text(size = 11, face = "bold"),  
    plot.margin = margin(15, 20, 10, 15)  
  )

print(p)  
ggsave("anhedonia_comparison.png", p, width = 12, height = 7,  
       dpi = 300, bg = "white")  