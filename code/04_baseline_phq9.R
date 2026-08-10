# ============================================================================  
# Baseline PHQ-9 Differences Between Switchers and iTBS-Only  
# ============================================================================

# Install effsize if not already installed  
if (!requireNamespace("effsize", quietly = TRUE)) {  
  install.packages("effsize")  
}

# Load required libraries  
library(readxl)  
library(ggplot2)  
library(effsize)

# --- Read in the data --------------------------------------------------------  
file_path <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx"

df <- as.data.frame(read_excel(file_path, sheet = "Sheet1"))

# --- Filter to baseline: visit == 0 ------------------------------------------  
baseline <- df[df$visit == 0, ]

cat("Total baseline (visit 0) rows:", nrow(baseline), "\n\n")

# --- Create group labels ------------------------------------------------------  
baseline$group_label <- ifelse(baseline$group == "switcher", "Switchers",  
                               ifelse(baseline$group == "iTBS_only", "iTBS-Only", NA))

# Keep only valid groups  
baseline <- baseline[!is.na(baseline$group_label), ]

# Remove rows with missing phq9_score  
baseline <- baseline[!is.na(baseline$phq9_score), ]

cat("Patients per group:\n")  
print(table(baseline$group_label))  
cat("\n")

# --- Compute summary statistics -----------------------------------------------  
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

# --- Welch's t-test -----------------------------------------------------------  
switcher_scores <- groups[["Switchers"]]  
itbs_scores     <- groups[["iTBS-Only"]]

t_result <- t.test(switcher_scores, itbs_scores, var.equal = FALSE)  
cat("Welch's t-test:\n")  
print(t_result)  
cat("\n")

# --- Cohen's d ----------------------------------------------------------------  
d_result <- cohen.d(switcher_scores, itbs_scores)  
cat("Cohen's d:\n")  
print(d_result)  
cat("\n")

# --- Extract values for annotation --------------------------------------------  
t_val <- round(t_result$statistic, 2)  
p_val <- t_result$p.value  
d_val <- round(abs(d_result$estimate), 2)

# Format p-value  
if (p_val < 0.001) {  
  p_display <- "< 0.001"  
} else {  
  p_display <- sprintf("%.2f", round(p_val, 2))  
}

subtitle_text <- paste0("Welch's t = ", t_val, ", p = ", p_display,   
                        ", Cohen's d = ", d_val)

cat("Subtitle:", subtitle_text, "\n\n")

# --- Set factor order: Switchers on left, iTBS-Only on right ------------------  
summary_stats$group_label <- factor(summary_stats$group_label,   
                                    levels = c("Switchers", "iTBS-Only"))

# --- X-axis labels with sample sizes ------------------------------------------  
x_labels <- c(  
  paste0("Switchers\n(n = ", summary_stats$n[summary_stats$group_label == "Switchers"], ")"),  
  paste0("iTBS-Only\n(n = ", summary_stats$n[summary_stats$group_label == "iTBS-Only"], ")")  
)

# --- Determine significance label ---------------------------------------------  
if (p_val < 0.001) {  
  sig_label <- "***"  
} else if (p_val < 0.01) {  
  sig_label <- "**"  
} else if (p_val < 0.05) {  
  sig_label <- "*"  
} else {  
  sig_label <- "ns"  
}

# --- Build the plot -----------------------------------------------------------  
bar_colors <- c("Switchers" = "#6B1C23", "iTBS-Only" = "#3A6FA0")

bracket_y    <- max(summary_stats$mean_phq9 + summary_stats$se_phq9) + 2.0  
star_y       <- bracket_y + 0.5  
label_offset <- 0.6

p <- ggplot(summary_stats, aes(x = group_label, y = mean_phq9, fill = group_label)) +  
  
  # Bars  
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.3) +  
  
  # Error bars (Mean +/- SE)  
  geom_errorbar(aes(ymin = mean_phq9 - se_phq9, ymax = mean_phq9 + se_phq9),  
                width = 0.10, linewidth = 0.7) +  
  
  # Mean value labels above error bars  
  geom_text(aes(y = mean_phq9 + se_phq9 + label_offset,  
                label = sprintf("%.1f", mean_phq9)),  
            size = 5.5, fontface = "bold", color = "gray30") +  
  
  # Significance bracket - horizontal line  
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,  
           linewidth = 0.5) +  
  # Left tick  
  annotate("segment", x = 1, xend = 1, y = bracket_y - 0.4, yend = bracket_y,  
           linewidth = 0.5) +  
  # Right tick  
  annotate("segment", x = 2, xend = 2, y = bracket_y - 0.4, yend = bracket_y,  
           linewidth = 0.5) +  
  # Significance label  
  annotate("text", x = 1.5, y = star_y, label = sig_label,   
           size = ifelse(sig_label == "ns", 5, 9), fontface = "bold") +  
  
  # Colors  
  scale_fill_manual(values = bar_colors) +  
  
  # X-axis labels  
  scale_x_discrete(labels = x_labels) +  
  
  # Y-axis  
  scale_y_continuous(limits = c(0, 22), breaks = seq(0, 20, 5), expand = c(0, 0)) +  
  
  # Labels  
  labs(  
    title = "Baseline PHQ-9 by Group",  
    subtitle = subtitle_text,  
    x = NULL,  
    y = "Baseline PHQ-9 Score (Mean ± SE)"  
  ) +  
  
  # Theme  
  theme_classic(base_size = 14) +  
  theme(  
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),  
    axis.text.x   = element_text(size = 13, face = "bold"),  
    axis.text.y   = element_text(size = 12),  
    axis.title.y  = element_text(size = 13),  
    legend.position = "none",  
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)  
  )

# --- Display and save ---------------------------------------------------------  
print(p)

ggsave(filename = "/Users/chiara/Documents/TMS-Switching-Study/Baseline_PHQ9_by_Group.png",  
       plot = p, width = 6, height = 6.5, dpi = 300)

cat("\nPlot saved to /Users/chiara/Documents/TMS-Switching-Study/Baseline_PHQ9_by_Group.png\n")  