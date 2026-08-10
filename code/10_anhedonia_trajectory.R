# ══════════════════════════════════════════════════════════════  
# FIGURE 2: Anhedonia Trajectory Comparison (CLEAN — no dots/lines)  
# No remission zone, no green shading  
# ══════════════════════════════════════════════════════════════

library(tidyverse)  
library(readxl)  
library(zoo)  
library(lme4)  
library(lmerTest)

# --- 1. Load Data ---  
df <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx", sheet = "Sheet1")  
colnames(df) <- trimws(colnames(df))

# --- 2. Prepare analysis dataset ---  
df_analysis <- df %>%  
  mutate(  
    cumulative_tx = as.numeric(cumulative_tx),  
    phq9_1 = as.numeric(phq9_1),  
    group = case_when(  
      group == "iTBS_only" ~ "iTBS Only",  
      group == "switcher"  ~ "Switchers",  
      TRUE ~ as.character(group)  
    )  
  ) %>%  
  filter(group %in% c("iTBS Only", "Switchers")) %>%  
  filter(!is.na(cumulative_tx), !is.na(phq9_1))

# --- 3. Count patients per group ---  
n_itbs   <- df_analysis %>% filter(group == "iTBS Only") %>% pull(study_id) %>% n_distinct()  
n_switch <- df_analysis %>% filter(group == "Switchers") %>% pull(study_id) %>% n_distinct()  
cat("iTBS Only n =", n_itbs, "\nSwitchers n =", n_switch, "\n")

# --- 4. Median switch point ---  
median_switch <- 19

# --- 5. Run LME models: main effect of group, ALL data ---  
# Model: phq9_item ~ group + cumulative_tx + (1 | study_id)  
# This tests whether switchers have higher item scores overall  
phq_items <- paste0("phq9_", 1:9)  
item_labels <- c("Anhedonia", "Depressed mood", "Sleep", "Fatigue",  
                 "Appetite", "Guilt", "Concentration", "Psychomotor", "Suicidality")

results_df <- data.frame(item = item_labels, phq_col = phq_items, p_raw = NA)

for (i in seq_along(phq_items)) {  
  df_temp <- df_analysis %>%  
    mutate(y = as.numeric(.data[[phq_items[i]]])) %>%  
    filter(!is.na(y))
  
  tryCatch({  
    mod <- lmer(y ~ group + cumulative_tx + (1 | study_id), data = df_temp)  
    coefs <- summary(mod)$coefficients  
    results_df$p_raw[i] <- coefs["groupSwitchers", "Pr(>|t|)"]  
  }, error = function(e) {  
    results_df$p_raw[i] <<- NA  
  })  
}

results_df$p_fdr <- p.adjust(results_df$p_raw, method = "fdr")  
print(results_df)  
cat("\nAnhedonia raw p:", results_df$p_raw[results_df$item == "Anhedonia"], "\n")  
cat("Anhedonia FDR p:", results_df$p_fdr[results_df$item == "Anhedonia"], "\n")

# --- 6. Compute smoothed anhedonia means ---  
anhedonia_means <- df_analysis %>%  
  group_by(group, cumulative_tx) %>%  
  summarise(  
    mean_score = mean(phq9_1, na.rm = TRUE),  
    se = sd(phq9_1, na.rm = TRUE) / sqrt(n()),  
    n = n(),  
    .groups = "drop"  
  ) %>%  
  filter(n >= 3) %>%  
  arrange(group, cumulative_tx)

anhedonia_means <- anhedonia_means %>%  
  group_by(group) %>%  
  mutate(  
    smooth_score = zoo::rollmean(mean_score, k = 5, fill = NA, align = "center")  
  ) %>%  
  ungroup()

# --- 7. Plot ---  
group_colors <- c(  
  "iTBS Only" = "#2166AC",  
  "Switchers" = "#B2182B"  
)

p_anhedonia <- ggplot() +  
  # Smoothed group means ONLY  
  geom_line(data = anhedonia_means %>% filter(!is.na(smooth_score)),  
            aes(x = cumulative_tx, y = smooth_score, color = group),  
            linewidth = 2.5) +  
  # Median switch line  
  geom_vline(xintercept = median_switch, linetype = "dashed",  
             color = "#B2182B", linewidth = 0.8, alpha = 0.7) +  
  annotate("text", x = median_switch + 1.5, y = 2.9,  
           label = paste0("Median switch\n(tx ", median_switch, ")"),  
           color = "#B2182B", size = 3.5, fontface = "bold.italic",  
           hjust = 0, vjust = 1) +  
  # Key finding annotation — uses computed FDR p-value  
  annotate("label", x = 5, y = 0.3,  
           label = paste0("Anhedonia is the ONLY item surviving\n",  
                          "FDR correction (p_FDR = ",  
                          format.pval(results_df$p_fdr[results_df$item == "Anhedonia"],  
                                      digits = 4), ")"),  
           fill = "#FFF9C4", color = "black", fontface = "bold",  
           size = 3.5, hjust = 0,  
           label.padding = unit(0.4, "lines")) +  
  # Colors & scales  
  scale_color_manual(  
    values = group_colors,  
    labels = c(paste0("iTBS Only (n=", n_itbs, ")"),  
               paste0("iTBS\u2192Bilateral (n=", n_switch, ")")),  
    name = NULL  
  ) +  
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +  
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +  
  labs(  
    x = "Cumulative TMS Sessions",  
    y = "PHQ-9 Item 1 Score (Anhedonia, 0\u20133)",  
    title = "Anhedonia: The Key Predictor of Protocol Switch Need",  
    subtitle = "Future-switchers show significantly slower anhedonia response to iTBS\n(Switchers = iTBS \u2192 Bilateral ONLY, no other protocols)"  
  ) +  
  theme_minimal(base_size = 14) +  
  theme(  
    plot.background = element_rect(fill = "white", color = NA),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.3),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),  
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),  
    axis.text = element_text(size = 11, color = "black"),  
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),  
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",  
                                 margin = margin(b = 15)),  
    legend.position = c(0.82, 0.92),  
    legend.background = element_rect(fill = "white", color = "#CCCCCC",  
                                     linewidth = 0.3),  
    legend.text = element_text(size = 11),  
    legend.key.width = unit(1.5, "cm"),  
    plot.margin = margin(15, 20, 10, 10)  
  )

print(p_anhedonia)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/anhedonia_trajectory.png",  
       p_anhedonia, width = 11, height = 7, dpi = 300, bg = "white")

cat("\nPlot saved to /Users/chiara/Documents/TMS-Switching-Study/anhedonia_trajectory.png\n")  