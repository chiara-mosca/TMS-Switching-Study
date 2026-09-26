# ============================================================================  
# HLM + Trajectory Plot: PHQ-9 by Treatment Group  
# Capped at 36 sessions | TMS074 excluded | Time in weeks  
# Model-predicted linear trajectories (no LOESS)  
# ============================================================================

pkgs <- c("readxl", "lme4", "lmerTest", "ggplot2")  
for (p in pkgs) {  
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)  
}

library(readxl)  
library(lme4)  
library(lmerTest)  
library(ggplot2)

# ==============================================================================  
# 1. READ & PREPARE DATA  
# ==============================================================================

file_path <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx"  
df <- as.data.frame(read_excel(file_path, sheet = "Sheet1"))

# Exclude TMS074, cap at 36 sessions  
df <- df[df$study_id != "TMS074", ]  
df <- df[df$cumulative_tx <= 36, ]

# Label groups  
df$group_label <- ifelse(df$group == "switcher", "Switchers",  
                         ifelse(df$group == "iTBS_only", "iTBS-Only", NA))  
df <- df[!is.na(df$group_label) & !is.na(df$phq9_score) & !is.na(df$phq9_date), ]

# Only keep matched switchers (both phases)  
switch_ids_pre <- unique(df$study_id[df$group == "switcher" & df$phase == "pre"])  
switch_ids_post <- unique(df$study_id[df$group == "switcher" & df$phase == "post"])  
switch_ids_both <- intersect(switch_ids_pre, switch_ids_post)

df <- df[df$group_label == "iTBS-Only" | df$study_id %in% switch_ids_both, ]

# Compute weeks from first PHQ-9 per patient  
df$phq9_date <- as.POSIXct(df$phq9_date)  
df <- do.call(rbind, lapply(split(df, df$study_id), function(x) {  
  x$weeks <- as.numeric(difftime(x$phq9_date, min(x$phq9_date), units = "weeks"))  
  x  
}))  
rownames(df) <- NULL

df$group_factor <- factor(df$group_label, levels = c("iTBS-Only", "Switchers"))

n_sw <- length(unique(df$study_id[df$group_label == "Switchers"]))  
n_it <- length(unique(df$study_id[df$group_label == "iTBS-Only"]))

cat(sprintf("iTBS-Only: %d patients\n", n_it))  
cat(sprintf("Switchers: %d patients\n", n_sw))  
cat(sprintf("TMS074 excluded: %s\n", !("TMS074" %in% df$study_id)))  
cat(sprintf("Total observations: %d\n", nrow(df)))  
cat(sprintf("Weeks range: %.1f to %.1f\n\n",  
            min(df$weeks), max(df$weeks)))

# ==============================================================================  
# 2. FIT HLM  
# ==============================================================================

model_A <- lmer(phq9_score ~ weeks * group_factor + (1 + weeks | study_id),  
                data = df, REML = TRUE,  
                control = lmerControl(optimizer = "bobyqa",  
                                      optCtrl = list(maxfun = 100000)))

if (any(grepl("failed to converge", model_A@optinfo$conv$lme4$messages))) {  
  cat("Random slope did not converge. Using random intercept only.\n\n")  
  model_A <- lmer(phq9_score ~ weeks * group_factor + (1 | study_id),  
                  data = df, REML = TRUE)  
} else {  
  cat("Model converged successfully.\n\n")  
}

# ==============================================================================  
# 3. EXTRACT STATISTICS  
# ==============================================================================

fe <- as.data.frame(coef(summary(model_A)))  
colnames(fe) <- c("Estimate", "SE", "df", "t_value", "p_value")  
fe$CI_lower <- fe$Estimate - 1.96 * fe$SE  
fe$CI_upper <- fe$Estimate + 1.96 * fe$SE

vc <- as.data.frame(VarCorr(model_A))  
resid_var <- sigma(model_A)^2  
int_var <- vc$vcov[vc$var1 == "(Intercept)" & is.na(vc$var2)]  
icc_val <- int_var / (int_var + resid_var)

# Group-specific intercepts and slopes  
int_itbs   <- fe$Estimate[1]  
slope_itbs <- fe$Estimate[2]  
int_sw     <- int_itbs + fe$Estimate[3]  
slope_sw   <- slope_itbs + fe$Estimate[4]

# --- Print results ------------------------------------------------------------  
nice_names <- c("Intercept (Baseline PHQ-9, iTBS-Only)",  
                "Weekly slope (iTBS-Only)",  
                "Group difference at baseline (Switchers - iTBS)",  
                "Group x Week interaction")

cat("=================================================================\n")  
cat("MODEL SUMMARY\n")  
cat("=================================================================\n\n")  
print(summary(model_A))

cat("\n=================================================================\n")  
cat("FIXED EFFECTS\n")  
cat("=================================================================\n\n")

for (i in 1:nrow(fe)) {  
  p_str <- ifelse(fe$p_value[i] < 0.001, "< .001", sprintf("%.3f", fe$p_value[i]))  
  cat(sprintf("  %-45s\n", nice_names[i]))  
  cat(sprintf("    B = %.3f, SE = %.3f, 95%% CI [%.3f, %.3f], t(%.1f) = %.3f, p = %s\n\n",  
              fe$Estimate[i], fe$SE[i], fe$CI_lower[i], fe$CI_upper[i],  
              fe$df[i], fe$t_value[i], p_str))  
}

cat("=================================================================\n")  
cat("GROUP-SPECIFIC TRAJECTORIES\n")  
cat("=================================================================\n\n")  
cat(sprintf("  iTBS-Only: intercept = %.1f, slope = %.3f pts/wk\n", int_itbs, slope_itbs))  
cat(sprintf("  Switchers: intercept = %.1f, slope = %.3f pts/wk\n\n", int_sw, slope_sw))

cat("=================================================================\n")  
cat("RANDOM EFFECTS\n")  
cat("=================================================================\n\n")  
cat(sprintf("  Random intercept variance: %.3f (SD = %.3f)\n", int_var, sqrt(int_var)))  
if (any(vc$var1 == "weeks" & is.na(vc$var2))) {  
  sv <- vc$vcov[vc$var1 == "weeks" & is.na(vc$var2)]  
  cat(sprintf("  Random slope variance:     %.3f (SD = %.3f)\n", sv, sqrt(sv)))  
}  
cat(sprintf("  Residual variance:         %.3f (SD = %.3f)\n", resid_var, sqrt(resid_var)))  
cat(sprintf("  ICC: %.3f (%.1f%%)\n\n", icc_val, icc_val * 100))

# --- Results table ------------------------------------------------------------  
cat("=================================================================\n")  
cat("RESULTS TABLE\n")  
cat("=================================================================\n\n")  
cat(sprintf("%-45s %8s %7s %18s %8s %8s\n", "Parameter", "B", "SE", "95% CI", "t", "p"))  
cat(paste(rep("-", 92), collapse = ""), "\n")  
for (i in 1:nrow(fe)) {  
  ci <- sprintf("[%.3f, %.3f]", fe$CI_lower[i], fe$CI_upper[i])  
  ps <- ifelse(fe$p_value[i] < 0.001, "< .001", sprintf("%.3f", fe$p_value[i]))  
  cat(sprintf("%-45s %8.3f %7.3f %18s %8.3f %8s\n",  
              nice_names[i], fe$Estimate[i], fe$SE[i], ci, fe$t_value[i], ps))  
}  
cat(paste(rep("-", 92), collapse = ""), "\n")  
cat(sprintf("%-45s %8.3f\n", "iTBS-Only slope (pts/wk)", slope_itbs))  
cat(sprintf("%-45s %8.3f\n", "Switchers slope (pts/wk)", slope_sw))  
cat(sprintf("%-45s %8.3f\n", "ICC", icc_val))  
cat(paste(rep("-", 92), collapse = ""), "\n\n")

# ==============================================================================  
# 4. TRAJECTORY PLOT — HLM MODEL-PREDICTED LINES  
# ==============================================================================

color_switcher <- "#B83030"  
color_itbs     <- "#2B5E8C"

# Median switch week  
switch_week_vals <- do.call(rbind, lapply(split(df[df$group_label == "Switchers", ],  
                                                df$study_id[df$group_label == "Switchers"]), function(x) {  
                                                  post_rows <- x[x$phase == "post", ]  
                                                  if (nrow(post_rows) > 0) {  
                                                    data.frame(study_id = x$study_id[1], switch_week = min(post_rows$weeks))  
                                                  } else {  
                                                    NULL  
                                                  }  
                                                }))  
median_switch_week <- median(switch_week_vals$switch_week, na.rm = TRUE)

# Week ranges per group  
max_wk_sw <- max(df$weeks[df$group_label == "Switchers"])  
max_wk_it <- max(df$weeks[df$group_label == "iTBS-Only"])  
max_wk <- max(max_wk_sw, max_wk_it)  
xmax_plot <- ceiling(max_wk) + 0.5

# Build model-predicted lines  
legend_sw <- sprintf("iTBS > Bilateral (n = %d)", n_sw)  
legend_it <- sprintf("iTBS-Only (n = %d)", n_it)

fit_df <- rbind(  
  data.frame(  
    weeks = seq(0, max_wk_it, length.out = 100),  
    phq9  = int_itbs + slope_itbs * seq(0, max_wk_it, length.out = 100),  
    Group = legend_it  
  ),  
  data.frame(  
    weeks = seq(0, max_wk_sw, length.out = 100),  
    phq9  = int_sw + slope_sw * seq(0, max_wk_sw, length.out = 100),  
    Group = legend_sw  
  )  
)

fit_df$Group <- factor(fit_df$Group, levels = c(legend_sw, legend_it))

# --- Build plot ---------------------------------------------------------------  
p <- ggplot(fit_df, aes(x = weeks, y = phq9, color = Group)) +  
  # Remission zone  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 5,  
           fill = "#E8F5E9", alpha = 0.5) +  
  geom_hline(yintercept = 5, linetype = "dashed", color = "#66BB6A",  
             linewidth = 0.6, alpha = 0.7) +  
  annotate("text", x = xmax_plot * 0.95, y = 3,  
           label = "remission", color = "#66BB6A", size = 2.8,  
           hjust = 1, fontface = "italic", alpha = 0.7) +  
  # Model-predicted lines  
  geom_line(linewidth = 2.5) +  
  # Median switch line  
  geom_vline(xintercept = median_switch_week, linetype = "longdash",  
             color = "gray30", linewidth = 0.8) +  
  annotate("text", x = median_switch_week + 0.15, y = 26.5,  
           label = sprintf("Median switch\n(week %.1f)", median_switch_week),  
           hjust = 0, size = 3.2, fontface = "italic", color = "gray30",  
           lineheight = 1.1) +  
  # Slope annotations  
  annotate("text",  
           x = max_wk_it * 0.65, y = int_itbs + slope_itbs * max_wk_it * 0.65 - 1.5,  
           label = sprintf("%+.2f pts/wk", slope_itbs),  
           color = color_itbs, size = 3.8, fontface = "bold") +  
  annotate("text",  
           x = max_wk_sw * 0.65, y = int_sw + slope_sw * max_wk_sw * 0.65 + 1.5,  
           label = sprintf("%+.2f pts/wk", slope_sw),  
           color = color_switcher, size = 3.8, fontface = "bold") +  
  # Colors  
  scale_color_manual(  
    values = setNames(c(color_switcher, color_itbs), c(legend_sw, legend_it)),  
    name = NULL  
  ) +  
  # Axes  
  scale_x_continuous(limits = c(0, xmax_plot),  
                     breaks = seq(0, floor(xmax_plot), 2),  
                     expand = c(0.01, 0)) +  
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 25, 5),  
                     expand = c(0, 0.3)) +  
  # Labels  
  labs(  
    title = "PHQ-9 Trajectories by Treatment Group",  
    x = "Weeks Since Treatment Start",  
    y = "PHQ-9 Score"  
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

ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_trajectories_by_group.png",  
       p, width = 9, height = 6, dpi = 300, bg = "white")  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_trajectories_by_group.pdf",  
       p, width = 9, height = 6, bg = "white")

cat("\nSaved: phq9_trajectories_by_group.png\n")  
cat("Saved: phq9_trajectories_by_group.pdf\n")

# ==============================================================================  
# 5. FINAL SUMMARY  
# ==============================================================================

cat("\n=================================================================\n")  
cat("FINAL SUMMARY\n")  
cat("=================================================================\n\n")  
cat(sprintf("iTBS-Only: %d patients, intercept = %.1f, slope = %.3f pts/wk\n",  
            n_it, int_itbs, slope_itbs))  
cat(sprintf("Switchers: %d patients, intercept = %.1f, slope = %.3f pts/wk\n",  
            n_sw, int_sw, slope_sw))  
cat(sprintf("\nTMS074 excluded: YES\n"))  
cat(sprintf("Sessions capped at: 36\n"))  
cat(sprintf("Time variable: weeks since treatment start\n"))  
cat(sprintf("Median switch point: week %.1f\n", median_switch_week))  
cat(sprintf("Group x Week interaction: p = %s\n",  
            ifelse(fe$p_value[4] < 0.001, "< .001", sprintf("%.3f", fe$p_value[4]))))

cat("\nDone.\n")  