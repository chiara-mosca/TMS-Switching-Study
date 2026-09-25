# ============================================================================  
# Figure 2: Piecewise HLM - Pre vs Post-Switch PHQ-9 Trajectory  
# Capped at 36 sessions | Weeks relative to switch point  
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

switch_dat <- dat %>%  
  filter(group == "switcher") %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  filter(!is.na(phq9_score) & !is.na(phq9_date))

ids_with_pre <- switch_dat %>% filter(phase == "pre") %>% pull(study_id) %>% unique()  
ids_with_post <- switch_dat %>% filter(phase == "post") %>% pull(study_id) %>% unique()  
ids_both <- intersect(ids_with_pre, ids_with_post)

dropped <- setdiff(ids_with_pre, ids_with_post)  
if (length(dropped) > 0) {  
  cat(sprintf("Dropped %d additional patient(s): %s\n",  
              length(dropped), paste(dropped, collapse = ", ")))  
}

switch_dat <- switch_dat %>% filter(study_id %in% ids_both)

n_switchers <- length(unique(switch_dat$study_id))  
cat(sprintf("N switchers: %d\n", n_switchers))  
cat(sprintf("Total observations: %d\n\n", nrow(switch_dat)))

# ----------------------------------------------------------------------------  
# 2. Compute weeks relative to switch point  
# ----------------------------------------------------------------------------

switch_dates <- switch_dat %>%  
  filter(phase == "post") %>%  
  group_by(study_id) %>%  
  summarise(switch_date = min(as.POSIXct(phq9_date)), .groups = "drop")

switch_dat <- switch_dat %>%  
  left_join(switch_dates, by = "study_id") %>%  
  mutate(  
    phq9_date = as.POSIXct(phq9_date),  
    weeks_relative = as.numeric(difftime(phq9_date, switch_date, units = "weeks")),  
    post = ifelse(phase == "post", 1, 0),  
    pre_weeks = ifelse(post == 0, weeks_relative, 0),  
    post_weeks = ifelse(post == 1, weeks_relative, 0),  
    all_weeks = weeks_relative,  
    post_increment = post_weeks  
  )

cat(sprintf("Weeks relative range: %.1f to %.1f\n\n",  
            min(switch_dat$weeks_relative),  
            max(switch_dat$weeks_relative)))

# ----------------------------------------------------------------------------  
# 3. Fit models  
# ----------------------------------------------------------------------------

model_piecewise <- lmer(  
  phq9_score ~ pre_weeks + post_weeks +  
    (1 + pre_weeks + post_weeks | study_id),  
  data = switch_dat, REML = TRUE,  
  control = lmerControl(optimizer = "bobyqa",  
                        optCtrl = list(maxfun = 100000))  
)

if (!is.null(model_piecewise@optinfo$conv$lme4$messages)) {  
  cat("WARNING: Piecewise model convergence issues:\n")  
  print(model_piecewise@optinfo$conv$lme4$messages)  
  cat("\n")  
} else {  
  cat("Piecewise model converged successfully.\n\n")  
}

model_change <- lmer(  
  phq9_score ~ all_weeks + post_increment +  
    (1 + all_weeks + post_increment | study_id),  
  data = switch_dat, REML = TRUE,  
  control = lmerControl(optimizer = "bobyqa",  
                        optCtrl = list(maxfun = 100000))  
)

if (!is.null(model_change@optinfo$conv$lme4$messages)) {  
  cat("WARNING: Reparameterized model convergence issues:\n")  
  print(model_change@optinfo$conv$lme4$messages)  
  cat("\n")  
} else {  
  cat("Reparameterized model converged successfully.\n\n")  
}

model_intercept_only <- lmer(  
  phq9_score ~ pre_weeks + post_weeks + (1 | study_id),  
  data = switch_dat, REML = TRUE  
)

# ----------------------------------------------------------------------------  
# 4. Model comparison  
# ----------------------------------------------------------------------------

model_comp <- anova(model_intercept_only, model_piecewise, refit = TRUE)  
print(model_comp)  
cat(sprintf("\nRandom slopes model: %s (p = %s)\n\n",  
            ifelse(model_comp$`Pr(>Chisq)`[2] < .05, "SIGNIFICANTLY BETTER",  
                   "not significantly better"),  
            format.pval(model_comp$`Pr(>Chisq)`[2], digits = 4)))

# ----------------------------------------------------------------------------  
# 5. Extract results  
# ----------------------------------------------------------------------------

fe <- summary(model_piecewise)$coefficients  
fe2 <- summary(model_change)$coefficients

pre_slope <- fe["pre_weeks", "Estimate"]  
pre_se <- fe["pre_weeks", "Std. Error"]  
pre_t <- fe["pre_weeks", "t value"]  
pre_df <- fe["pre_weeks", "df"]  
pre_p <- fe["pre_weeks", "Pr(>|t|)"]

post_slope <- fe["post_weeks", "Estimate"]  
post_se <- fe["post_weeks", "Std. Error"]  
post_t <- fe["post_weeks", "t value"]  
post_df <- fe["post_weeks", "df"]  
post_p <- fe["post_weeks", "Pr(>|t|)"]

intercept <- fe["(Intercept)", "Estimate"]

slope_change <- fe2["post_increment", "Estimate"]  
slope_change_se <- fe2["post_increment", "Std. Error"]  
slope_change_t <- fe2["post_increment", "t value"]  
slope_change_df <- fe2["post_increment", "df"]  
slope_change_p <- fe2["post_increment", "Pr(>|t|)"]

pre_ci_lower <- pre_slope - 1.96 * pre_se  
pre_ci_upper <- pre_slope + 1.96 * pre_se  
post_ci_lower <- post_slope - 1.96 * post_se  
post_ci_upper <- post_slope + 1.96 * post_se  
change_ci_lower <- slope_change - 1.96 * slope_change_se  
change_ci_upper <- slope_change + 1.96 * slope_change_se

vc <- as.data.frame(VarCorr(model_piecewise))  
resid_var <- sigma(model_piecewise)^2

# ----------------------------------------------------------------------------  
# 6. Print results  
# ----------------------------------------------------------------------------

format_p <- function(p) {  
  if (p < .001) return("< .001")  
  return(sprintf("%.4f", p))  
}

cat("=================================================================\n")  
cat("PIECEWISE MODEL SUMMARY\n")  
cat("=================================================================\n\n")  
print(summary(model_piecewise))

cat("\n=================================================================\n")  
cat("REPARAMETERIZED MODEL\n")  
cat("=================================================================\n\n")  
print(summary(model_change))

cat("\n=================================================================\n")  
cat("KEY RESULTS\n")  
cat("=================================================================\n\n")

cat("PRE-SWITCH SLOPE (iTBS phase):\n")  
cat(sprintf("  B = %.3f pts/wk, SE = %.3f, 95%% CI [%.3f, %.3f]\n",  
            pre_slope, pre_se, pre_ci_lower, pre_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n\n",  
            pre_df, pre_t, format_p(pre_p)))

cat("POST-SWITCH SLOPE (bilateral phase):\n")  
cat(sprintf("  B = %.3f pts/wk, SE = %.3f, 95%% CI [%.3f, %.3f]\n",  
            post_slope, post_se, post_ci_lower, post_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n\n",  
            post_df, post_t, format_p(post_p)))

cat("SLOPE CHANGE (acceleration after switch):\n")  
cat(sprintf("  B = %.3f pts/wk, SE = %.3f, 95%% CI [%.3f, %.3f]\n",  
            slope_change, slope_change_se, change_ci_lower, change_ci_upper))  
cat(sprintf("  t(%.1f) = %.3f, p = %s\n\n",  
            slope_change_df, slope_change_t, format_p(slope_change_p)))

cat(sprintf("INTERCEPT (PHQ-9 at switch point): %.3f\n\n", intercept))

cat("RANDOM EFFECTS:\n")  
print(vc[, c("grp", "var1", "var2", "vcov", "sdcor")])  
cat(sprintf("\nResidual variance: %.3f (SD = %.3f)\n", resid_var, sqrt(resid_var)))

cat("\n=================================================================\n")  
cat("RESULTS TABLE\n")  
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
  cat(sprintf("%-40s %8.3f %7.3f %20s %10s %8s\n",  
              param[[1]], param[[2]], param[[3]], ci, tdf, ps))  
}  
cat(paste(rep("-", 95), collapse = ""), "\n\n")

# ----------------------------------------------------------------------------  
# 7. Build prediction lines  
# ----------------------------------------------------------------------------

x_min <- min(switch_dat$weeks_relative)  
x_max <- max(switch_dat$weeks_relative)

pred_range_pre <- data.frame(  
  weeks_relative = seq(x_min, 0, length.out = 100)  
) %>%  
  mutate(predicted = intercept + pre_slope * weeks_relative)

pred_range_post <- data.frame(  
  weeks_relative = seq(0, x_max, length.out = 100)  
) %>%  
  mutate(predicted = intercept + post_slope * weeks_relative)

# ----------------------------------------------------------------------------  
# 8. Build publication-ready figure  
# ----------------------------------------------------------------------------

pre_annotation <- sprintf("Pre-switch: %+.2f pts/wk\n(p = %s)",  
                          pre_slope, format_p(pre_p))

post_annotation <- sprintf("Post-switch: %+.2f pts/wk\n(slope change = %+.2f, p = %s)",  
                           post_slope, slope_change, format_p(slope_change_p))

col_pre_fit  <- "#CC3333"  
col_post_fit <- "#660000"

p <- ggplot() +  
  # Remission zone  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 5,  
           fill = "#E8F5E9", alpha = 0.5) +  
  geom_hline(yintercept = 5, linetype = "dashed", color = "#66BB6A",  
             linewidth = 0.6, alpha = 0.7) +  
  annotate("text", x = x_max * 0.95, y = 3,  
           label = "remission", color = "#66BB6A", size = 2.8,  
           hjust = 1, fontface = "italic", alpha = 0.7) +  
  # Pre-switch fitted line  
  geom_line(data = pred_range_pre,  
            aes(x = weeks_relative, y = predicted),  
            color = col_pre_fit, linewidth = 2.5) +  
  # Post-switch fitted line  
  geom_line(data = pred_range_post,  
            aes(x = weeks_relative, y = predicted),  
            color = col_post_fit, linewidth = 2.5) +  
  # Switch point  
  geom_vline(xintercept = 0, linetype = "longdash", linewidth = 0.8,  
             color = "gray30") +  
  annotate("text", x = 0.15, y = 26.5, label = "Switch",  
           fontface = "bold.italic", size = 3.5, hjust = 0, color = "gray30") +  
  # Pre-switch annotation  
  annotate("text",  
           x = mean(c(x_min, 0)), y = 24,  
           label = pre_annotation,  
           color = col_pre_fit, size = 3.8, fontface = "bold",  
           lineheight = 1.1, hjust = 0.5) +  
  # Post-switch annotation  
  annotate("text",  
           x = mean(c(0, x_max)), y = 24,  
           label = post_annotation,  
           color = col_post_fit, size = 3.8, fontface = "bold",  
           lineheight = 1.1, hjust = 0.5) +  
  # Axes  
  scale_x_continuous(breaks = seq(floor(x_min), ceiling(x_max), 2)) +  
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 25, 5),  
                     expand = c(0, 0.3)) +  
  labs(  
    title = sprintf("PHQ-9 Trajectory Before and After Protocol Switch (n = %d)",  
                    n_switchers),  
    x = "Weeks Relative to Switch Point",  
    y = "PHQ-9 Score"  
  ) +  
  theme_classic(base_size = 13) +  
  theme(  
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14,  
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
    plot.margin   = margin(12, 18, 12, 12)  
  )

print(p)

ggsave("/Users/chiara/Documents/TMS-Switching-Study/figure2_piecewise.png",  
       p, width = 10, height = 6, dpi = 300)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/figure2_piecewise.pdf",  
       p, width = 10, height = 6)

cat("\nSaved: figure2_piecewise.png\n")  
cat("Saved: figure2_piecewise.pdf\n")

# ----------------------------------------------------------------------------  
# 9. Final summary  
# ----------------------------------------------------------------------------

cat("\n=================================================================\n")  
cat("FINAL SUMMARY\n")  
cat("=================================================================\n\n")  
cat(sprintf("N patients: %d\n", n_switchers))  
cat(sprintf("Sessions capped at: 36\n"))  
cat(sprintf("Time variable: weeks relative to switch point\n"))  
cat(sprintf("Model: piecewise HLM with random intercepts and slopes\n\n"))  
cat(sprintf("Pre-switch slope:  %+.3f pts/wk, p = %s\n",  
            pre_slope, format_p(pre_p)))  
cat(sprintf("Post-switch slope: %+.3f pts/wk, p = %s\n",  
            post_slope, format_p(post_p)))  
cat(sprintf("Slope change:      %+.3f pts/wk, p = %s\n",  
            slope_change, format_p(slope_change_p)))

cat("\nDone.\n")  