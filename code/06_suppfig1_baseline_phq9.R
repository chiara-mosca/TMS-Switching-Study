# ============================================================================  
# Baseline PHQ-9 Differences Between Switchers and iTBS-Only  
# Capped at 36 sessions | TMS074 excluded  
# ============================================================================

if (!requireNamespace("effsize", quietly = TRUE)) {  
  install.packages("effsize")  
}

library(readxl)  
library(ggplot2)  
library(effsize)

# --- 1. Read and prepare data ------------------------------------------------  
file_path <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx"  
df <- as.data.frame(read_excel(file_path, sheet = "Sheet1"))

# Exclude TMS074 and cap at 36 sessions  
df <- df[df$study_id != "TMS074", ]  
df <- df[df$cumulative_tx <= 36, ]

# Filter to baseline: visit == 0  
baseline <- df[df$visit == 0, ]

# Create group labels  
baseline$group_label <- ifelse(baseline$group == "switcher", "Switchers",  
                               ifelse(baseline$group == "iTBS_only", "iTBS-Only", NA))

baseline <- baseline[!is.na(baseline$group_label), ]  
baseline <- baseline[!is.na(baseline$phq9_score), ]

# Only keep matched switchers (both phases present after cap)  
switch_ids_pre <- unique(df$study_id[df$group == "switcher" & df$phase == "pre"])  
switch_ids_post <- unique(df$study_id[df$group == "switcher" & df$phase == "post"])  
switch_ids_both <- intersect(switch_ids_pre, switch_ids_post)

baseline <- baseline[baseline$group_label == "iTBS-Only" |  
                       baseline$study_id %in% switch_ids_both, ]

cat(sprintf("TMS074 excluded: %s\n", !("TMS074" %in% baseline$study_id)))  
cat(sprintf("Patients per group:\n"))  
print(table(baseline$group_label))  
cat("\n")

# --- 2. Summary statistics ---------------------------------------------------  
groups <- split(baseline$phq9_score, baseline$group_label)

summary_stats <- data.frame(  
  group_label = c("Switchers", "iTBS-Only"),  
  n = c(length(groups[["Switchers"]]), length(groups[["iTBS-Only"]])),  
  mean_phq9 = c(mean(groups[["Switchers"]], na.rm = TRUE),  
                mean(groups[["iTBS-Only"]], na.rm = TRUE)),  
  sd_phq9 = c(sd(groups[["Switchers"]], na.rm = TRUE),  
              sd(groups[["iTBS-Only"]], na.rm = TRUE)),  
  se_phq9 = c(sd(groups[["Switchers"]], na.rm = TRUE) / sqrt(length(groups[["Switchers"]])),  
              sd(groups[["iTBS-Only"]], na.rm = TRUE) / sqrt(length(groups[["iTBS-Only"]])))  
)

cat("Summary Statistics:\n")  
print(summary_stats)  
cat("\n")

# --- 3. Welch's t-test -------------------------------------------------------  
switcher_scores <- groups[["Switchers"]]  
itbs_scores     <- groups[["iTBS-Only"]]

t_result <- t.test(switcher_scores, itbs_scores, var.equal = FALSE)  
cat("Welch's t-test:\n")  
print(t_result)  
cat("\n")

# --- 4. Cohen's d -------------------------------------------------------------  
d_result <- cohen.d(switcher_scores, itbs_scores)  
cat("Cohen's d:\n")  
print(d_result)  
cat("\n")

# --- 5. Extract values for annotation ----------------------------------------  
t_val <- round(t_result$statistic, 2)  
p_val <- t_result$p.value  
d_val <- round(abs(d_result$estimate), 2)  
df_val <- round(t_result$parameter, 1)

if (p_val < 0.001) {  
  p_display <- "< .001"  
} else {  
  p_display <- sprintf("%.3f", p_val)  
}

# Significance label  
if (p_val < 0.001) {  
  sig_label <- "***"  
} else if (p_val < 0.01) {  
  sig_label <- "**"  
} else if (p_val < 0.05) {  
  sig_label <- "*"  
} else {  
  sig_label <- "n.s."  
}

# --- 6. Set factor order -----------------------------------------------------  
summary_stats$group_label <- factor(summary_stats$group_label,  
                                    levels = c("Switchers", "iTBS-Only"))

x_labels <- c(  
  paste0("Switchers\n(n = ", summary_stats$n[summary_stats$group_label == "Switchers"], ")"),  
  paste0("iTBS-Only\n(n = ", summary_stats$n[summary_stats$group_label == "iTBS-Only"], ")")  
)

# --- 7. Build plot ------------------------------------------------------------  
bar_colors <- c("Switchers" = "#B83030", "iTBS-Only" = "#2B5E8C")

bracket_y    <- max(summary_stats$mean_phq9 + summary_stats$se_phq9) + 1.5  
star_y       <- bracket_y + 0.5  
label_offset <- 0.6

p <- ggplot(summary_stats, aes(x = group_label, y = mean_phq9, fill = group_label)) +
  
  # Bars  
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.3) +
  
  # Error bars  
  geom_errorbar(aes(ymin = mean_phq9 - se_phq9, ymax = mean_phq9 + se_phq9),  
                width = 0.10, linewidth = 0.6) +
  
  # Mean value labels  
  geom_text(aes(y = mean_phq9 + se_phq9 + label_offset,  
                label = sprintf("%.1f", mean_phq9)),  
            size = 4.5, fontface = "bold", color = "gray20") +
  
  # Significance bracket  
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,  
           linewidth = 0.4) +  
  annotate("segment", x = 1, xend = 1, y = bracket_y - 0.3, yend = bracket_y,  
           linewidth = 0.4) +  
  annotate("segment", x = 2, xend = 2, y = bracket_y - 0.3, yend = bracket_y,  
           linewidth = 0.4) +  
  annotate("text", x = 1.5, y = star_y, label = sig_label,  
           size = ifelse(sig_label == "n.s.", 4, 7), fontface = "bold") +
  
  # Colors  
  scale_fill_manual(values = bar_colors) +
  
  # X-axis labels  
  scale_x_discrete(labels = x_labels) +
  
  # Y-axis  
  scale_y_continuous(limits = c(0, 22), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  
  # Labels  
  labs(  
    title = "Baseline PHQ-9 by Treatment Group",  
    x = NULL,  
    y = "Baseline PHQ-9 Score (Mean \u00b1 SE)"  
  ) +
  
  # Theme  
  theme_classic(base_size = 13) +  
  theme(  
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14,  
                                 margin = margin(b = 12)),  
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),  
    axis.text.y   = element_text(size = 10, color = "black"),  
    axis.title.y  = element_text(size = 11, face = "bold",  
                                 margin = margin(r = 8)),  
    axis.line     = element_line(color = "black", linewidth = 0.5),  
    axis.ticks    = element_line(color = "black", linewidth = 0.4),  
    axis.ticks.length = unit(0.15, "cm"),  
    legend.position = "none",  
    plot.margin = margin(12, 18, 12, 12)  
  )

print(p)

# --- 8. Save ------------------------------------------------------------------  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/baseline_phq9_by_group.png",  
       p, width = 6, height = 6, dpi = 300, bg = "white")  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/baseline_phq9_by_group.pdf",  
       p, width = 6, height = 6, bg = "white")

cat("\nSaved: baseline_phq9_by_group.png\n")  
cat("Saved: baseline_phq9_by_group.pdf\n")

# --- 9. Full summary ---------------------------------------------------------  
cat("\n=================================================================\n")  
cat("RESULTS SUMMARY\n")  
cat("=================================================================\n\n")  
cat(sprintf("Switchers: M = %.1f (SD = %.1f), n = %d\n",  
            summary_stats$mean_phq9[1], summary_stats$sd_phq9[1], summary_stats$n[1]))  
cat(sprintf("iTBS-Only: M = %.1f (SD = %.1f), n = %d\n",  
            summary_stats$mean_phq9[2], summary_stats$sd_phq9[2], summary_stats$n[2]))  
cat(sprintf("\nWelch's t(%s) = %s, p = %s\n", df_val, t_val, p_display))  
cat(sprintf("Cohen's d = %s\n", d_val))  
cat(sprintf("TMS074 excluded: YES\n"))  
cat(sprintf("Sessions capped at: 36\n"))

cat("\nDone.\n")  