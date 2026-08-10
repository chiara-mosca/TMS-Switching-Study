# ============================================================================  
# Figure 2: Piecewise HLM — Pre vs Post-Switch Slope (Corrected Annotations)  
# ============================================================================

library(readxl)  
library(dplyr)  
library(lme4)  
library(lmerTest)  
library(ggplot2)

# ----------------------------------------------------------------------------  
# 1. Load and prepare data  
# ----------------------------------------------------------------------------

dat <- read_excel("analytic_sample_long.xlsx")  
switch_dat <- dat %>% filter(group == "switcher")

# Create piecewise variables  
switch_dat <- switch_dat %>%  
  group_by(study_id) %>%  
  mutate(  
    sessions_relative = cumulative_tx - switch_tx,  
    post = ifelse(phase == "post", 1, 0),  
    pre_sessions = ifelse(post == 0, sessions_relative, 0),  
    post_sessions = ifelse(post == 1, sessions_relative, 0)  
  ) %>%  
  ungroup()

# ----------------------------------------------------------------------------  
# 2. Fit models  
# ----------------------------------------------------------------------------

# Piecewise model (separate slopes)  
model_piecewise <- lmer(  
  phq9_score ~ pre_sessions + post_sessions + (1 | study_id),  
  data = switch_dat, REML = TRUE  
)

# Reparameterized model (direct slope change test)  
switch_dat <- switch_dat %>%  
  mutate(  
    all_sessions = sessions_relative,  
    post_increment = post_sessions  
  )

model_change <- lmer(  
  phq9_score ~ all_sessions + post_increment + (1 | study_id),  
  data = switch_dat, REML = TRUE  
)

# Extract coefficients  
fe <- summary(model_piecewise)$coefficients  
fe2 <- summary(model_change)$coefficients

pre_slope <- fe["pre_sessions", "Estimate"]  
post_slope <- fe["post_sessions", "Estimate"]  
pre_p <- fe["pre_sessions", "Pr(>|t|)"]  
post_p <- fe["post_sessions", "Pr(>|t|)"]  
slope_change <- fe2["post_increment", "Estimate"]  
slope_change_p <- fe2["post_increment", "Pr(>|t|)"]

n_switchers <- length(unique(switch_dat$study_id))

# ----------------------------------------------------------------------------  
# 3. Build prediction lines  
# ----------------------------------------------------------------------------

pred_range_pre <- data.frame(  
  sessions_relative = seq(min(switch_dat$sessions_relative[switch_dat$post == 0]),  
                          0, length.out = 100)  
)  
pred_range_pre$predicted <- fe["(Intercept)", "Estimate"] +  
  fe["pre_sessions", "Estimate"] * pred_range_pre$sessions_relative

pred_range_post <- data.frame(  
  sessions_relative = seq(0, max(switch_dat$sessions_relative[switch_dat$post == 1]),  
                          length.out = 100)  
)  
pred_range_post$predicted <- fe["(Intercept)", "Estimate"] +  
  fe["post_sessions", "Estimate"] * pred_range_post$sessions_relative

# ----------------------------------------------------------------------------  
# 4. Create significance labels (α = .05)  
# ----------------------------------------------------------------------------

# Pre-switch label: show significance correctly  
if (pre_p < .001) {  
  pre_p_label <- "p < .001"  
} else if (pre_p < .05) {  
  pre_p_label <- sprintf("p = %.3f", pre_p)  
} else {  
  pre_p_label <- sprintf("p = %.3f, NS", pre_p)  
}

pre_annotation <- sprintf("Pre-switch slope:\n%.3f pts/session\n(%s)",  
                          pre_slope, pre_p_label)

# Post-switch / slope change label  
if (slope_change_p < .001) {  
  change_p_label <- "p < .001"  
} else {  
  change_p_label <- sprintf("p = %.4f", slope_change_p)  
}

post_annotation <- sprintf("Post-switch slope:\n%.3f pts/session\n(slope change = %.3f, %s)",  
                           post_slope, slope_change, change_p_label)

# Caption  
caption_text <- sprintf(  
  "Pre-switch slope: %.3f (p = %.3f) | Post-switch slope: %.3f (p < .001) | Slope change: %.3f (%s)",  
  pre_slope, pre_p, post_slope, slope_change, change_p_label  
)

# ----------------------------------------------------------------------------  
# 5. Build the figure  
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
    subtitle = sprintf("n = %d patients | Random intercept piecewise HLM", n_switchers),  
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
ggsave("figure2_piecewise_trajectories.png", p, width = 12, height = 7, dpi = 300)  
ggsave("figure2_piecewise_trajectories.pdf", p, width = 12, height = 7)

cat("\nFigure saved.\n")  
cat(sprintf("Pre-switch slope: %.4f, p = %.4f — %s at α = .05\n",  
            pre_slope, pre_p, ifelse(pre_p < .05, "SIGNIFICANT", "NOT SIGNIFICANT")))  
cat(sprintf("Post-switch slope: %.4f, p = %.8f\n", post_slope, post_p))  
cat(sprintf("Slope change: %.4f, p = %.4f\n", slope_change, slope_change_p))  