# ============================================================================  
# Figure 2: Piecewise HLM — Pre vs Post-Switch Slope  
# CORRECTED: Random intercepts AND random slopes  
# ============================================================================

library(readxl)  
library(dplyr)  
library(lme4)  
library(lmerTest)  
library(ggplot2)

# ----------------------------------------------------------------------------  
# 1. Load and prepare data  
# ----------------------------------------------------------------------------

dat <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx")  
switch_dat <- dat %>% filter(group == "switcher")

# Create piecewise variables  
# Knot is individually determined: each patient's own switch_tx session  
switch_dat <- switch_dat %>%  
  group_by(study_id) %>%  
  mutate(  
    sessions_relative = cumulative_tx - switch_tx,  
    post = ifelse(phase == "post", 1, 0),  
    pre_sessions = ifelse(post == 0, sessions_relative, 0),  
    post_sessions = ifelse(post == 1, sessions_relative, 0)  
  ) %>%  
  ungroup()

# Also create reparameterized variables for slope change test  
switch_dat <- switch_dat %>%  
  mutate(  
    all_sessions = sessions_relative,  
    post_increment = post_sessions  
  )

n_switchers <- length(unique(switch_dat$study_id))  
cat("N switchers:", n_switchers, "\n")  
cat("Total observations:", nrow(switch_dat), "\n\n")

# ----------------------------------------------------------------------------  
# 2. Fit models — WITH RANDOM SLOPES  
# ----------------------------------------------------------------------------

cat("=================================================================\n")  
cat("FITTING MODELS WITH RANDOM INTERCEPTS AND RANDOM SLOPES\n")  
cat("=================================================================\n\n")

# --- Model 1: Piecewise (separate slopes for each phase) ---  
model_piecewise <- lmer(  
  phq9_score ~ pre_sessions + post_sessions +  
    (1 + pre_sessions + post_sessions | study_id),  
  data = switch_dat, REML = TRUE,  
  control = lmerControl(optimizer = "bobyqa",  
                        optCtrl = list(maxfun = 100000))  
)

# Check convergence  
if (!is.null(model_piecewise@optinfo$conv$lme4$messages)) {  
  cat("WARNING: Piecewise model convergence issues:\n")  
  print(model_piecewise@optinfo$conv$lme4$messages)  
  cat("\n")  
} else {  
  cat("Piecewise model converged successfully with random slopes.\n\n")  
}

# --- Model 2: Reparameterized (direct slope change test) ---  
model_change <- lmer(  
  phq9_score ~ all_sessions + post_increment +  
    (1 + all_sessions + post_increment | study_id),  
  data = switch_dat, REML = TRUE,  
  control = lmerControl(optimizer = "bobyqa",  
                        optCtrl = list(maxfun = 100000))  
)

if (!is.null(model_change@optinfo$conv$lme4$messages)) {  
  cat("WARNING: Reparameterized model convergence issues:\n")  
  print(model_change@optinfo$conv$lme4$messages)  
  cat("\n")  
} else {  
  cat("Reparameterized model converged successfully with random slopes.\n\n")  
}

# --- Model 3: Intercept-only (for model comparison) ---  
model_intercept_only <- lmer(  
  phq9_score ~ pre_sessions + post_sessions + (1 | study_id),  
  data = switch_dat, REML = TRUE  
)

# ----------------------------------------------------------------------------  
# 3. Model comparison: random slopes vs intercept-only  
# ----------------------------------------------------------------------------

cat("=================================================================\n")  
cat("MODEL COMPARISON: Random slopes vs intercept-only\n")  
cat("=================================================================\n\n")

model_comp <- anova(model_intercept_only, model_piecewise, refit = TRUE)  
print(model_comp)

cat(sprintf("\nRandom slopes model is %s (p = %s)\n\n",  
            ifelse(model_comp$`Pr(>Chisq)`[2] < .05, "SIGNIFICANTLY BETTER", "not significantly better"),  
            format.pval(model_comp$`Pr(>Chisq)`[2], digits = 4)))

# ----------------------------------------------------------------------------  
# 4. Extract all results  
# ----------------------------------------------------------------------------

fe <- summary(model_piecewise)$coefficients  
fe2 <- summary(model_change)$coefficients

pre_slope <- fe["pre_sessions", "Estimate"]  
pre_se <- fe["pre_sessions", "Std. Error"]  
pre_t <- fe["pre_sessions", "t value"]  
pre_df <- fe["pre_sessions", "df"]  
pre_p <- fe["pre_sessions", "Pr(>|t|)"]

post_slope <- fe["post_sessions", "Estimate"]  
post_se <- fe["post_sessions", "Std. Error"]  
post_t <- fe["post_sessions", "t value"]  
post_df <- fe["post_sessions", "df"]  
post_p <- fe["post_sessions", "Pr(>|t|)"]

intercept <- fe["(Intercept)", "Estimate"]

slope_change <- fe2["post_increment", "Estimate"]  
slope_change_se <- fe2["post_increment", "Std. Error"]  
slope_change_t <- fe2["post_increment", "t value"]  
slope_change_df <- fe2["post_increment", "df"]  
slope_change_p <- fe2["post_increment", "Pr(>|t|)"]

# 95% CIs  
pre_ci_lower <- pre_slope - 1.96 * pre_se  
pre_ci_upper <- pre_slope + 1.96 * pre_se  
post_ci_lower <- post_slope - 1.96 * post_se  
post_ci_upper <- post_slope + 1.96 * post_se  
change_ci_lower <- slope_change - 1.96 * slope_change_se  
change_ci_upper <- slope_change + 1.96 * slope_change_se

# Random effects  
vc <- as.data.frame(VarCorr(model_piecewise))  
resid_var <- sigma(model_piecewise)^2

# ----------------------------------------------------------------------------  
# 5. Print all results  
# ----------------------------------------------------------------------------

cat("=================================================================\n")  
cat("PIECEWISE MODEL SUMMARY (Random intercepts + random slopes)\n")  
cat("=================================================================\n\n")  
print(summary(model_piecewise))

cat("\n=================================================================\n")  
cat("REPARAMETERIZED MODEL (Direct slope change test)\n")  
cat("=================================================================\n\n")  
print(summary(model_change))

cat("\n=================================================================\n")  
cat("KEY RESULTS\n")  
cat("=================================================================\n\n")

format_p <- function(p) {  
  if (p < .001) return("< .001")  
  return(sprintf("%.4f", p))  
}

cat("PRE-SWITCH SLOPE (iTBS only phase):\n")  
cat(sprintf("  B = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",  
            pre_slope, pre_se, pre_ci_lower, pre_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n",  
            pre_df, pre_t, format_p(pre_p)))  
cat(sprintf("  Interpretation: %s improvement before switch\n\n",  
            ifelse(pre_p < .05, "SIGNIFICANT", "NO SIGNIFICANT")))

cat("POST-SWITCH SLOPE (bilateral phase):\n")  
cat(sprintf("  B = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",  
            post_slope, post_se, post_ci_lower, post_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n",  
            post_df, post_t, format_p(post_p)))  
cat(sprintf("  Interpretation: %s improvement after switch\n\n",  
            ifelse(post_p < .05, "SIGNIFICANT", "NO SIGNIFICANT")))

cat("SLOPE CHANGE (acceleration after switch):\n")  
cat(sprintf("  B = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",  
            slope_change, slope_change_se, change_ci_lower, change_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n",  
            slope_change_df, slope_change_t, format_p(slope_change_p)))  
cat(sprintf("  Interpretation: Slope %s significantly steeper after switch\n\n",  
            ifelse(slope_change_p < .05, "IS", "is NOT")))

cat("INTERCEPT (PHQ-9 at switch point):\n")  
cat(sprintf("  B = %.4f\n\n", intercept))

cat("RANDOM EFFECTS:\n")  
print(vc[, c("grp", "var1", "var2", "vcov", "sdcor")])  
cat(sprintf("\nResidual variance: %.4f (SD = %.4f)\n", resid_var, sqrt(resid_var)))

# --- Poster table ---  
cat("\n=================================================================\n")  
cat("POSTER TABLE\n")  
cat("=================================================================\n\n")  
cat(sprintf("%-40s %8s %7s %20s %10s %8s\n",  
            "Parameter", "B", "SE", "95% CI", "t(df)", "p"))  
cat(paste(rep("-", 95), collapse = ""), "\n")

params <- list(  
  list("Pre-switch slope (iTBS)", pre_slope, pre_se, pre_ci_lower, pre_ci_upper, pre_t, pre_df, pre_p),  
  list("Post-switch slope (bilateral)", post_slope, post_se, post_ci_lower, post_ci_upper, post_t, post_df, post_p),  
  list("Slope change (post - pre)", slope_change, slope_change_se, change_ci_lower, change_ci_upper, slope_change_t, slope_change_df, slope_change_p)  
)

for (param in params) {  
  ci <- sprintf("[%.3f, %.3f]", param[[4]], param[[5]])  
  tdf <- sprintf("%.2f(%.1f)", param[[6]], param[[7]])  
  ps <- format_p(param[[8]])  
  cat(sprintf("%-40s %8.4f %7.4f %20s %10s %8s\n",  
              param[[1]], param[[2]], param[[3]], ci, tdf, ps))  
}  
cat(paste(rep("-", 95), collapse = ""), "\n\n")

# ----------------------------------------------------------------------------  
# 6. Build prediction lines using CORRECTED model estimates  
# ----------------------------------------------------------------------------

pred_range_pre <- data.frame(  
  sessions_relative = seq(min(switch_dat$sessions_relative[switch_dat$post == 0]),  
                          0, length.out = 100)  
)  
pred_range_pre$predicted <- intercept +  
  pre_slope * pred_range_pre$sessions_relative

pred_range_post <- data.frame(  
  sessions_relative = seq(0, max(switch_dat$sessions_relative[switch_dat$post == 1]),  
                          length.out = 100)  
)  
pred_range_post$predicted <- intercept +  
  post_slope * pred_range_post$sessions_relative

# ----------------------------------------------------------------------------  
# 7. Create figure annotations  
# ----------------------------------------------------------------------------

# Pre-switch label  
if (pre_p < .05) {  
  pre_p_label <- sprintf("p = %.3f", pre_p)  
} else {  
  pre_p_label <- sprintf("p = %.3f, NS", pre_p)  
}

pre_annotation <- sprintf("Pre-switch slope:\n%.3f pts/session\n(%s)",  
                          pre_slope, pre_p_label)

# Post-switch / slope change label  
post_annotation <- sprintf(  
  "Post-switch slope:\n%.3f pts/session\n(slope change = %.3f, p = %s)",  
  post_slope, slope_change, format_p(slope_change_p))

# Caption  
caption_text <- sprintf(  
  "Pre-switch: %.3f (p = %s) | Post-switch: %.3f (p %s) | Slope change: %.3f (p = %s) | Random intercepts + slopes",  
  pre_slope, format_p(pre_p),  
  post_slope, format_p(post_p),  
  slope_change, format_p(slope_change_p)  
)

# ----------------------------------------------------------------------------  
# 8. Build the figure  
# ----------------------------------------------------------------------------

p <- ggplot() +  
  # Pre-switch fitted line  
  geom_line(data = pred_range_pre,  
            aes(x = sessions_relative, y = predicted, color = "Pre-switch (iTBS)"),  
            linewidth = 2.2) +  
  # Post-switch fitted line  
  geom_line(data = pred_range_post,  
            aes(x = sessions_relative, y = predicted, color = "Post-switch (Bilateral)"),  
            linewidth = 2.2) +  
  # Switch point vertical line  
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.5, color = "black") +  
  # Switch point label  
  annotate("label", x = 0, y = 27.5, label = "SWITCH POINT",  
           fontface = "bold", size = 4, fill = "lightyellow",  
           label.padding = unit(0.35, "lines")) +  
  # Pre-switch annotation  
  annotate("text", x = -15, y = 22.5, hjust = 0.5,  
           label = pre_annotation,  
           color = "#CC3333", size = 4, fontface = "bold",  
           lineheight = 1.1) +  
  # Post-switch annotation  
  annotate("text", x = 20, y = 22.5, hjust = 0.5,  
           label = post_annotation,  
           color = "#660000", size = 4, fontface = "bold",  
           lineheight = 1.1) +  
  # Color scale  
  scale_color_manual(  
    name = "Phase",  
    values = c("Pre-switch (iTBS)" = "#CC3333",  
               "Post-switch (Bilateral)" = "#660000")  
  ) +  
  # Axes  
  scale_y_continuous(limits = c(0, 28), breaks = seq(0, 25, 5)) +  
  # Labels  
  labs(  
    title = "PHQ-9 Trajectories Before and After Protocol Switch",  
    subtitle = sprintf("n = %d patients | Piecewise HLM with random intercepts and slopes",  
                       n_switchers),  
    x = "Sessions Relative to Switch Point",  
    y = "PHQ-9 Score",  
    caption = caption_text  
  ) +  
  # Theme  
  theme_minimal(base_size = 14) +  
  theme(  
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),  
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),  
    plot.caption = element_text(hjust = 0.5, color = "grey50", size = 9),  
    legend.position = "none",  
    panel.grid.minor = element_blank(),  
    axis.title = element_text(face = "bold")  
  )

# Display  
print(p)

# Save  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/figure2_piecewise_corrected.png",  
       p, width = 12, height = 7, dpi = 300)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/figure2_piecewise_corrected.pdf",  
       p, width = 12, height = 7)

cat("\n=================================================================\n")  
cat("FIGURE SAVED\n")  
cat("=================================================================\n\n")

cat(sprintf("Pre-switch slope: %.4f, p = %s — %s at alpha = .05\n",  
            pre_slope, format_p(pre_p),  
            ifelse(pre_p < .05, "SIGNIFICANT", "NOT SIGNIFICANT")))  
cat(sprintf("Post-switch slope: %.4f, p = %s\n",  
            post_slope, format_p(post_p)))  
cat(sprintf("Slope change: %.4f, p = %s\n",  
            slope_change, format_p(slope_change_p)))  
cat(sprintf("\nModel: Random intercepts + random slopes (significantly better than\n"))  
cat(sprintf("intercept-only: chi-sq(5) = %.2f, p = %s)\n",  
            model_comp$Chisq[2], format.pval(model_comp$`Pr(>Chisq)`[2], digits = 4)))

cat("\nFiles saved to /Users/chiara/Documents/TMS-Switching-Study/\n")  
cat("  - figure2_piecewise_corrected.png\n")  
cat("  - figure2_piecewise_corrected.pdf\n")  
cat("\nDone.\n")  