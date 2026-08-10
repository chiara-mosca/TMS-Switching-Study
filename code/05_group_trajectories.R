# ============================================================================  
# Model A: HLM + Trajectory Plot - Capped at 40 sessions  
# ============================================================================

pkgs <- c("readxl", "lme4", "lmerTest", "ggplot2", "performance")  
for (p in pkgs) {  
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)  
}

library(readxl)  
library(lme4)  
library(lmerTest)  
library(ggplot2)  
library(performance)

# ==============================================================================  
# 1. READ & PREPARE DATA  
# ==============================================================================

file_path <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx"  
df <- as.data.frame(read_excel(file_path, sheet = "Sheet1"))

df$group_label <- ifelse(df$group == "switcher", "Switchers",  
                         ifelse(df$group == "iTBS_only", "iTBS-Only", NA))  
df <- df[!is.na(df$group_label) & !is.na(df$phq9_score) & !is.na(df$cumulative_tx), ]  
df$time <- df$cumulative_tx  
df$group_factor <- factor(df$group_label, levels = c("iTBS-Only", "Switchers"))

n_sw <- length(unique(df$study_id[df$group_label == "Switchers"]))  
n_it <- length(unique(df$study_id[df$group_label == "iTBS-Only"]))

cat("Switchers n =", n_sw, "\niTBS-Only n =", n_it, "\n")  
cat("Total observations:", nrow(df), "\n\n")

# ==============================================================================  
# 2. FIT HLM  
# ==============================================================================

model_A <- lmer(phq9_score ~ time * group_factor + (1 + time | study_id),  
                data = df, REML = TRUE,  
                control = lmerControl(optimizer = "bobyqa",  
                                      optCtrl = list(maxfun = 100000)))

if (any(grepl("failed to converge", model_A@optinfo$conv$lme4$messages))) {  
  cat("Random slope did not converge. Using random intercept only.\n\n")  
  model_A <- lmer(phq9_score ~ time * group_factor + (1 | study_id),  
                  data = df, REML = TRUE)  
}

# ==============================================================================  
# 3. ALL STATISTICS  
# ==============================================================================

fe <- as.data.frame(coef(summary(model_A)))  
colnames(fe) <- c("Estimate", "SE", "df", "t_value", "p_value")  
fe$CI_lower <- fe$Estimate - 1.96 * fe$SE  
fe$CI_upper <- fe$Estimate + 1.96 * fe$SE

vc <- as.data.frame(VarCorr(model_A))  
resid_var <- sigma(model_A)^2  
int_var <- vc$vcov[vc$var1 == "(Intercept)" & is.na(vc$var2)]  
icc_val <- int_var / (int_var + resid_var)

slope_itbs <- fe$Estimate[rownames(fe) == "time"]  
slope_sw <- slope_itbs + fe$Estimate[grep("time:group_factor", rownames(fe))]

# --- Print everything ---------------------------------------------------------  
cat("=================================================================\n")  
cat("MODEL A: PHQ-9 ~ Sessions * Group + (1 + time | ID)\n")  
cat("=================================================================\n\n")  
print(summary(model_A))

cat("\n=================================================================\n")  
cat("FIXED EFFECTS\n")  
cat("=================================================================\n\n")

nice_names <- c("Intercept (Baseline PHQ-9, iTBS-Only)",  
                "Time/session slope (iTBS-Only)",  
                "Group diff at baseline (Switchers - iTBS)",  
                "Group x Time interaction")

for (i in 1:nrow(fe)) {  
  p_str <- ifelse(fe$p_value[i] < 0.001, "< .001", sprintf("%.3f", fe$p_value[i]))  
  cat(sprintf("  %-45s\n", nice_names[i]))  
  cat(sprintf("    B=%.4f, SE=%.4f, 95%%CI[%.4f, %.4f], t(%.1f)=%.3f, p=%s\n\n",  
              fe$Estimate[i], fe$SE[i], fe$CI_lower[i], fe$CI_upper[i],  
              fe$df[i], fe$t_value[i], p_str))  
}

cat("=================================================================\n")  
cat("GROUP-SPECIFIC SLOPES\n")  
cat("=================================================================\n\n")  
cat(sprintf("  iTBS-Only: B = %.4f per session\n", slope_itbs))  
cat(sprintf("  Switchers: B = %.4f per session\n\n", slope_sw))

cat("=================================================================\n")  
cat("RANDOM EFFECTS\n")  
cat("=================================================================\n\n")  
cat(sprintf("  Random intercept variance: %.4f (SD=%.4f)\n", int_var, sqrt(int_var)))  
if (any(vc$var1 == "time" & is.na(vc$var2))) {  
  sv <- vc$vcov[vc$var1 == "time" & is.na(vc$var2)]  
  cat(sprintf("  Random slope variance:     %.4f (SD=%.4f)\n", sv, sqrt(sv)))  
}  
cat(sprintf("  Residual variance:         %.4f (SD=%.4f)\n\n", resid_var, sqrt(resid_var)))

cat("=================================================================\n")  
cat("ICC\n")  
cat("=================================================================\n\n")  
cat(sprintf("  ICC = %.4f (%.1f%% between-patient variance)\n\n", icc_val, icc_val * 100))

# --- Poster table -------------------------------------------------------------  
cat("=================================================================\n")  
cat("POSTER TABLE\n")  
cat("=================================================================\n\n")  
cat(sprintf("%-45s %8s %7s %18s %8s %8s\n", "Parameter", "B", "SE", "95% CI", "t", "p"))  
cat(paste(rep("-", 92), collapse = ""), "\n")  
for (i in 1:nrow(fe)) {  
  ci <- paste0("[", sprintf("%.3f", fe$CI_lower[i]), ", ", sprintf("%.3f", fe$CI_upper[i]), "]")  
  ps <- ifelse(fe$p_value[i] < 0.001, "< .001", sprintf("%.3f", fe$p_value[i]))  
  cat(sprintf("%-45s %8.3f %7.3f %18s %8.3f %8s\n",  
              nice_names[i], fe$Estimate[i], fe$SE[i], ci, fe$t_value[i], ps))  
}  
cat(paste(rep("-", 92), collapse = ""), "\n")  
cat(sprintf("%-45s %8.4f\n", "iTBS-Only slope", slope_itbs))  
cat(sprintf("%-45s %8.4f\n", "Switchers slope", slope_sw))  
cat(sprintf("%-45s %8.4f\n", "Random intercept var", int_var))  
cat(sprintf("%-45s %8.4f\n", "Residual var", resid_var))  
cat(sprintf("%-45s %8.4f\n", "ICC", icc_val))  
cat(paste(rep("-", 92), collapse = ""), "\n\n")

# ==============================================================================  
# 4. TRAJECTORY PLOT - BOTH GROUPS CAPPED AT 40 SESSIONS  
# ==============================================================================

color_switcher <- "#6B1C23"  
color_itbs     <- "#2B547E"

# Median switch  
sw_base <- df[df$group_label == "Switchers" & !duplicated(df$study_id), ]  
median_switch <- median(sw_base$switch_tx, na.rm = TRUE)

# --- CAP BOTH GROUPS AT 40 CUMULATIVE SESSIONS -------------------------------  
cap <- 40

obs_sw <- df[df$group_label == "Switchers" & df$time <= cap, c("time", "phq9_score")]  
obs_it <- df[df$group_label == "iTBS-Only" & df$time <= cap, c("time", "phq9_score")]

loess_sw <- loess(phq9_score ~ time, data = obs_sw, span = 0.5)  
loess_it <- loess(phq9_score ~ time, data = obs_it, span = 0.5)

time_sw <- seq(min(obs_sw$time), cap, length.out = 300)  
time_it <- seq(min(obs_it$time), cap, length.out = 300)

legend_sw <- paste0("iTBS \u2192 Bilateral (n=", n_sw, ")")  
legend_it <- paste0("iTBS Only (n=", n_it, ")")

smooth_df <- rbind(  
  data.frame(time = time_sw,  
             phq9 = predict(loess_sw, data.frame(time = time_sw)),  
             Group = legend_sw),  
  data.frame(time = time_it,  
             phq9 = predict(loess_it, data.frame(time = time_it)),  
             Group = legend_it)  
)  
smooth_df <- smooth_df[!is.na(smooth_df$phq9), ]

smooth_df$Group <- factor(smooth_df$Group, levels = c(legend_sw, legend_it))

# --- Plot ---------------------------------------------------------------------  
p <- ggplot(smooth_df, aes(x = time, y = phq9, color = Group)) +  
  geom_line(linewidth = 2.5) +  
  geom_vline(xintercept = median_switch, linetype = "dashed",  
             color = color_switcher, linewidth = 1.5) +  
  annotate("text", x = median_switch - 1.5, y = 26,  
           label = paste0("Median switch\n(tx ", median_switch, ")"),  
           hjust = 1, size = 4.5, fontface = "italic", color = color_switcher) +  
  scale_color_manual(  
    values = setNames(c(color_switcher, color_itbs), c(legend_sw, legend_it)),  
    name = NULL  
  ) +  
  scale_x_continuous(limits = c(0, 45), breaks = seq(0, 40, 10),  
                     expand = c(0.01, 0)) +  
  scale_y_continuous(limits = c(0, 28), breaks = seq(0, 25, 5),  
                     expand = c(0, 0.5)) +  
  labs(title = "PHQ-9 Trajectories by Treatment Group",  
       x = "Cumulative TMS Sessions",  
       y = "PHQ-9 Score") +  
  theme_classic(base_size = 14) +  
  theme(  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12),  
    legend.position = c(0.80, 0.95),  
    legend.justification = c(0.5, 1),  
    legend.background = element_rect(fill = "white", color = NA),  
    legend.key.width = unit(1.8, "cm"),  
    legend.text = element_text(size = 12, face = "bold"),  
    plot.margin = margin(15, 25, 10, 10)  
  ) +  
  guides(color = guide_legend(override.aes = list(linewidth = 3.5)))

print(p)

ggsave(filename = "/Users/chiara/Documents/TMS-Switching-Study/Model_A_PHQ9_Trajectories.png",  
       plot = p, width = 9, height = 6.5, dpi = 300)

cat("\nPlot saved. Done.\n")  