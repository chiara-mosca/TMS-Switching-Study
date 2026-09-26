# ============================================================================  
# Anhedonia Trajectory Comparison: iTBS-Only vs Switchers  
# Capped at 36 sessions | TMS074 excluded | Time in weeks  
# HLM model-predicted lines  
# ============================================================================

library(tidyverse)  
library(readxl)  
library(lme4)  
library(lmerTest)

# --- 1. Load & Prepare Data ---  
df <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                 sheet = "Sheet1")  
colnames(df) <- trimws(colnames(df))

df_analysis <- df %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  mutate(  
    phq9_1 = as.numeric(phq9_1),  
    phq9_date = as.POSIXct(phq9_date),  
    group = case_when(  
      group == "iTBS_only" ~ "iTBS Only",  
      group == "switcher"  ~ "Switchers",  
      TRUE ~ as.character(group)  
    )  
  ) %>%  
  filter(group %in% c("iTBS Only", "Switchers")) %>%  
  filter(!is.na(phq9_1), !is.na(phq9_date))

# Only keep matched switchers (both phases)  
switch_ids_pre <- df_analysis %>%  
  filter(group == "Switchers" & phase == "pre") %>%  
  pull(study_id) %>% unique()  
switch_ids_post <- df_analysis %>%  
  filter(group == "Switchers" & phase == "post") %>%  
  pull(study_id) %>% unique()  
switch_ids_both <- intersect(switch_ids_pre, switch_ids_post)

df_analysis <- df_analysis %>%  
  filter(group == "iTBS Only" | study_id %in% switch_ids_both)

# Compute weeks from first PHQ-9 per patient  
df_analysis <- df_analysis %>%  
  group_by(study_id) %>%  
  mutate(  
    weeks = as.numeric(difftime(phq9_date, min(phq9_date), units = "weeks"))  
  ) %>%  
  ungroup()

df_analysis$group_factor <- factor(df_analysis$group, levels = c("iTBS Only", "Switchers"))

n_itbs   <- n_distinct(df_analysis$study_id[df_analysis$group == "iTBS Only"])  
n_switch <- n_distinct(df_analysis$study_id[df_analysis$group == "Switchers"])

cat(sprintf("iTBS Only: %d patients\n", n_itbs))  
cat(sprintf("Switchers: %d patients\n", n_switch))  
cat(sprintf("TMS074 excluded: %s\n", !("TMS074" %in% df_analysis$study_id)))  
cat(sprintf("Sessions capped at: 36\n"))  
cat(sprintf("Weeks range: %.1f to %.1f\n\n",  
            min(df_analysis$weeks), max(df_analysis$weeks)))

# --- 2. Run item-level LME models (main effect of group) ---  
phq_items <- paste0("phq9_", 1:9)  
item_labels <- c("Anhedonia", "Depressed Mood", "Sleep", "Fatigue",  
                 "Appetite", "Guilt/Worthlessness", "Concentration",  
                 "Psychomotor", "Suicidal Ideation")

results_df <- tibble(item = item_labels, phq_col = phq_items, p_raw = NA_real_)

for (i in seq_along(phq_items)) {  
  df_temp <- df_analysis %>%  
    mutate(y = as.numeric(.data[[phq_items[i]]])) %>%  
    filter(!is.na(y))
  
  tryCatch({  
    mod <- lmer(y ~ group_factor + weeks + (1 | study_id), data = df_temp)  
    coefs <- summary(mod)$coefficients  
    group_row <- grep("group_factor", rownames(coefs))  
    results_df$p_raw[i] <- coefs[group_row, "Pr(>|t|)"]  
  }, error = function(e) {  
    cat("Error for", item_labels[i], ":", e$message, "\n")  
  })  
}

results_df <- results_df %>%  
  mutate(p_fdr = p.adjust(p_raw, method = "fdr"))

cat("\n=== Item-Level Group Effects (FDR-corrected) ===\n\n")  
print(results_df %>% select(item, p_raw, p_fdr) %>% arrange(p_fdr))  
cat(sprintf("\nAnhedonia raw p:  %.4f\n", results_df$p_raw[results_df$item == "Anhedonia"]))  
cat(sprintf("Anhedonia FDR p:  %.4f\n\n", results_df$p_fdr[results_df$item == "Anhedonia"]))

# --- 3. Fit anhedonia HLM for trajectory lines ---  
anhedonia_model <- lmer(  
  phq9_1 ~ weeks * group_factor + (1 + weeks | study_id),  
  data = df_analysis, REML = TRUE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))  
)

if (any(grepl("failed to converge", anhedonia_model@optinfo$conv$lme4$messages))) {  
  cat("Random slope did not converge. Using random intercept only.\n\n")  
  anhedonia_model <- lmer(  
    phq9_1 ~ weeks * group_factor + (1 | study_id),  
    data = df_analysis, REML = TRUE  
  )  
} else {  
  cat("Anhedonia model converged successfully.\n\n")  
}

fe <- coef(summary(anhedonia_model))

int_itbs   <- fe[1, "Estimate"]  
slope_itbs <- fe[2, "Estimate"]  
int_sw     <- int_itbs + fe[3, "Estimate"]  
slope_sw   <- slope_itbs + fe[4, "Estimate"]

cat(sprintf("iTBS-Only: intercept = %.2f, slope = %.3f pts/wk\n", int_itbs, slope_itbs))  
cat(sprintf("Switchers: intercept = %.2f, slope = %.3f pts/wk\n\n", int_sw, slope_sw))

# --- 4. Median switch week ---  
switch_week_vals <- df_analysis %>%  
  filter(group == "Switchers" & phase == "post") %>%  
  group_by(study_id) %>%  
  summarise(switch_week = min(weeks), .groups = "drop")  
median_switch_week <- median(switch_week_vals$switch_week, na.rm = TRUE)

cat(sprintf("Median switch point: week %.1f\n\n", median_switch_week))

# --- 5. Build model-predicted lines ---  
max_wk_it <- max(df_analysis$weeks[df_analysis$group == "iTBS Only"])  
max_wk_sw <- max(df_analysis$weeks[df_analysis$group == "Switchers"])  
max_wk <- max(max_wk_it, max_wk_sw)  
xmax_plot <- ceiling(max_wk) + 0.5

legend_it <- sprintf("iTBS-Only (n = %d)", n_itbs)  
legend_sw <- sprintf("iTBS > Bilateral (n = %d)", n_switch)

fit_df <- rbind(  
  data.frame(  
    weeks = seq(0, max_wk_it, length.out = 100),  
    score = int_itbs + slope_itbs * seq(0, max_wk_it, length.out = 100),  
    Group = legend_it  
  ),  
  data.frame(  
    weeks = seq(0, max_wk_sw, length.out = 100),  
    score = int_sw + slope_sw * seq(0, max_wk_sw, length.out = 100),  
    Group = legend_sw  
  )  
)

# Clamp to 0-3 range  
fit_df$score <- pmax(fit_df$score, 0)  
fit_df$score <- pmin(fit_df$score, 3)

fit_df$Group <- factor(fit_df$Group, levels = c(legend_sw, legend_it))

# --- 6. Build plot ---  
color_itbs     <- "#2B5E8C"  
color_switcher <- "#B83030"

# FDR annotation  
fdr_p <- results_df$p_fdr[results_df$item == "Anhedonia"]  
if (fdr_p < 0.001) {  
  fdr_display <- "< .001"  
} else {  
  fdr_display <- sprintf("%.3f", fdr_p)  
}

p <- ggplot(fit_df, aes(x = weeks, y = score, color = Group)) +  
  # Model-predicted lines  
  geom_line(linewidth = 2.5) +  
  # Median switch line  
  geom_vline(xintercept = median_switch_week, linetype = "longdash",  
             color = "gray30", linewidth = 0.8) +  
  annotate("text", x = median_switch_week + 0.15, y = 2.95,  
           label = sprintf("Median switch\n(week %.1f)", median_switch_week),  
           hjust = 0, size = 3.2, fontface = "italic", color = "gray30",  
           lineheight = 1.1) +  
  # Slope annotations  
  annotate("text",  
           x = max_wk_it * 0.6,  
           y = int_itbs + slope_itbs * max_wk_it * 0.6 - 0.2,  
           label = sprintf("%+.3f pts/wk", slope_itbs),  
           color = color_itbs, size = 3.8, fontface = "bold") +  
  annotate("text",  
           x = max_wk_sw * 0.6,  
           y = int_sw + slope_sw * max_wk_sw * 0.6 + 0.2,  
           label = sprintf("%+.3f pts/wk", slope_sw),  
           color = color_switcher, size = 3.8, fontface = "bold") +  
  # FDR annotation  
  annotate("text",  
           x = xmax_plot * 0.5, y = 0.15,  
           label = sprintf("Group effect: p(FDR) = %s", fdr_display),  
           size = 3.5, fontface = "bold.italic", color = "gray30") +  
  # Colors  
  scale_color_manual(  
    values = setNames(c(color_switcher, color_itbs), c(legend_sw, legend_it)),  
    name = NULL  
  ) +  
  # Axes  
  scale_x_continuous(limits = c(0, xmax_plot),  
                     breaks = seq(0, floor(xmax_plot), 2),  
                     expand = c(0.01, 0)) +  
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5),  
                     expand = c(0, 0.05)) +  
  # Labels  
  labs(  
    title = "Anhedonia (PHQ-9 Item 1) Trajectory by Treatment Group",  
    x = "Weeks Since Treatment Start",  
    y = "PHQ-9 Item 1 Score (0-3)"  
  ) +  
  # Theme  
  theme_classic(base_size = 13) +  
  theme(  
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14,  
                                 margin = margin(b = 12)),  
    axis.title.x  = element_text(size = 11, face = "bold", color = "black",  
                                 margin = margin(t = 8)),  
    axis.title.y  = element_text(size = 11, face = "bold", color = "black",  
                                 margin = margin(r = 8)),  
    axis.text.x   = element_text(size = 10, color = "black"),  
    axis.text.y   = element_text(size = 10, color = "black"),  
    axis.line     = element_line(color = "black", linewidth = 0.5),  
    axis.ticks    = element_line(color = "black", linewidth = 0.4),  
    axis.ticks.length = unit(0.15, "cm"),  
    legend.position = c(0.75, 0.95),  
    legend.justification = c(0.5, 1),  
    legend.background = element_rect(fill = "white", color = NA),  
    legend.key.width = unit(1.5, "cm"),  
    legend.text = element_text(size = 10),  
    plot.margin = margin(12, 18, 12, 12)  
  ) +  
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

print(p)

ggsave("/Users/chiara/Documents/TMS-Switching-Study/anhedonia_trajectory.png",  
       p, width = 9, height = 6, dpi = 300, bg = "white")  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/anhedonia_trajectory.pdf",  
       p, width = 9, height = 6, bg = "white")

cat("\nSaved: anhedonia_trajectory.png\n")  
cat("Saved: anhedonia_trajectory.pdf\n")

# --- 7. Final Summary ---  
cat("\n=================================================================\n")  
cat("FINAL SUMMARY\n")  
cat("=================================================================\n\n")  
cat(sprintf("iTBS-Only: %d patients, anhedonia slope = %.3f pts/wk\n", n_itbs, slope_itbs))  
cat(sprintf("Switchers: %d patients, anhedonia slope = %.3f pts/wk\n", n_switch, slope_sw))  
cat(sprintf("Median switch: week %.1f\n", median_switch_week))  
cat(sprintf("Anhedonia group effect: p(raw) = %.4f, p(FDR) = %s\n",  
            results_df$p_raw[results_df$item == "Anhedonia"], fdr_display))  
cat(sprintf("\nTMS074 excluded: YES\n"))  
cat(sprintf("Sessions capped at: 36\n"))  
cat(sprintf("Time variable: weeks since treatment start\n"))

cat("\nDone.\n")  