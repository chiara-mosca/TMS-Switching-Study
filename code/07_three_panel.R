# ══════════════════════════════════════════════════════════════  
# 3-PANEL PHQ-9 SLOPE PLOT — PUBLICATION READY  
# New dataset: analytic_sample_long.xlsx  
# ══════════════════════════════════════════════════════════════

library(dplyr)  
library(ggplot2)  
library(readxl)  
library(lme4)  
library(lmerTest)  
library(patchwork)

# ── 1. Load data ──  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

c1 <- raw %>%  
  dplyr::filter(!is.na(phq9_score) & !is.na(phq9_date)) %>%  
  arrange(study_id, phq9_date)

# ── 2. Classify patients using the group column ──  
stayer_ids <- c1 %>%  
  dplyr::filter(group == "iTBS_only") %>%  
  pull(study_id) %>%  
  unique()

switcher_ids <- c1 %>%  
  dplyr::filter(group == "switcher") %>%  
  pull(study_id) %>%  
  unique()

cat(sprintf("Stayers:   %d patients\n", length(stayer_ids)))  
cat(sprintf("Switchers: %d patients\n", length(switcher_ids)))

# ── 3. Build analysis dataset ──  
dat_stay <- c1 %>%  
  dplyr::filter(study_id %in% stayer_ids) %>%  
  mutate(panel_group = "stay_p1")

dat_switch_p1 <- c1 %>%  
  dplyr::filter(study_id %in% switcher_ids & phase == "pre") %>%  
  mutate(panel_group = "switch_p1")

dat_switch_p2 <- c1 %>%  
  dplyr::filter(study_id %in% switcher_ids & phase == "post") %>%  
  mutate(panel_group = "switch_p2")

g3 <- bind_rows(dat_stay, dat_switch_p1, dat_switch_p2)

# ── 4. Compute weeks_in_phase ──  
g3 <- g3 %>%  
  mutate(phq9_date = as.POSIXct(phq9_date)) %>%  
  group_by(study_id, panel_group) %>%  
  mutate(  
    phase_start = min(phq9_date),  
    weeks_in_phase = as.numeric(difftime(phq9_date, phase_start, units = "weeks"))  
  ) %>%  
  ungroup()

# ── 5. Set reference level ──  
g3$panel_group <- factor(g3$panel_group, levels = c("switch_p1", "stay_p1", "switch_p2"))

# ── 6. Fit HLM ──  
cat("\nFitting hierarchical linear mixed model...\n")

m3 <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3,  
  REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)

# ── 7. Extract fixed effects ──  
fe <- fixef(m3)

ref_int   <- fe["(Intercept)"]  
ref_slope <- fe["weeks_in_phase"]

int_stay_p1     <- ref_int + fe["panel_groupstay_p1"]  
slope_stay_p1   <- ref_slope + fe["weeks_in_phase:panel_groupstay_p1"]

int_switch_p1   <- ref_int  
slope_switch_p1 <- ref_slope

int_switch_p2   <- ref_int + fe["panel_groupswitch_p2"]  
slope_switch_p2 <- ref_slope + fe["weeks_in_phase:panel_groupswitch_p2"]

cat("\nPHQ-9 SLOPE PER GROUP (points/week):\n")  
cat(sprintf("  Stayers (never switched)  : %+.3f\n", slope_stay_p1))  
cat(sprintf("  Switchers phase 1         : %+.3f\n", slope_switch_p1))  
cat(sprintf("  Switchers phase 2         : %+.3f\n", slope_switch_p2))

# ── 8. Contrasts ──  
vcov_fe  <- vcov(m3)  
fe_names <- names(fe)

slope_contrast <- function(pos, neg = NULL) {  
  r <- matrix(0, nrow = 1, ncol = length(fe))  
  colnames(r) <- fe_names  
  r[1, pos] <- 1  
  if (!is.null(neg)) r[1, neg] <- -1  
  est  <- as.numeric(r %*% fe)  
  se   <- as.numeric(sqrt(r %*% vcov_fe %*% t(r)))  
  tval <- est / se  
  df_res <- df.residual(m3)  
  pval <- 2 * pt(abs(tval), df = df_res, lower.tail = FALSE)  
  ci_lo <- est - 1.96 * se  
  ci_hi <- est + 1.96 * se  
  return(c(estimate = est, ci_lo = ci_lo, ci_hi = ci_hi, p = pval))  
}

within_c <- slope_contrast("weeks_in_phase:panel_groupswitch_p2")  
btw_p1   <- slope_contrast("weeks_in_phase:panel_groupstay_p1")  
btw_p2   <- slope_contrast("weeks_in_phase:panel_groupstay_p1",  
                           "weeks_in_phase:panel_groupswitch_p2")

n_stayers   <- length(stayer_ids)  
n_switchers <- length(switcher_ids)

cat(sprintf("\nWITHIN-subject (same %d switchers, paired):\n", n_switchers))  
cat(sprintf("  phase 2 - phase 1        : %+.3f/wk  [%+.3f, %+.3f]  p=%.4f\n",  
            within_c["estimate"], within_c["ci_lo"], within_c["ci_hi"], within_c["p"]))  
cat(sprintf("\nBETWEEN-subject (%d stayers vs %d switchers):\n", n_stayers, n_switchers))  
cat(sprintf("  stayers - switch phase 1 : %+.3f/wk  [%+.3f, %+.3f]  p=%.4f\n",  
            btw_p1["estimate"], btw_p1["ci_lo"], btw_p1["ci_hi"], btw_p1["p"]))  
cat(sprintf("  stayers - switch phase 2 : %+.3f/wk  [%+.3f, %+.3f]  p=%.4f\n",  
            btw_p2["estimate"], btw_p2["ci_lo"], btw_p2["ci_hi"], btw_p2["p"]))

# ══════════════════════════════════════════════════════════════  
# 9. PUBLICATION-READY 3-PANEL PLOT  
# ══════════════════════════════════════════════════════════════

# ── Color scheme ──  
col_stay_line <- "#C6DEF1"  
col_stay_fit  <- "#2B5E8C"

col_sw1_line  <- "#FACBCB"  
col_sw1_fit   <- "#B83030"

col_sw2_line  <- "#F09090"  
col_sw2_fit   <- "#7A1515"

ymax <- max(g3$phq9_score, na.rm = TRUE) + 1

# ── Shared publication theme ──  
pub_theme <- theme_classic(base_size = 14) +  
  theme(  
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14,  
                              margin = margin(b = 6)),  
    axis.title.x = element_text(size = 12, face = "bold", color = "black",  
                                margin = margin(t = 10)),  
    axis.title.y = element_text(size = 12, face = "bold", color = "black",  
                                margin = margin(r = 10)),  
    axis.text.x  = element_text(size = 11, face = "bold", color = "black"),  
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),  
    axis.line    = element_line(color = "black", linewidth = 0.7),  
    axis.ticks   = element_line(color = "black", linewidth = 0.5),  
    axis.ticks.length = unit(0.2, "cm"),  
    plot.margin = margin(12, 18, 12, 12)  
  )

# ── Helper function to build one panel ──  
make_panel <- function(data, group_key, title_label, title_color,  
                       line_color, fit_color, intercept, slope,  
                       show_y = FALSE) {
  
  sub <- data %>%  
    dplyr::filter(panel_group == group_key) %>%  
    arrange(study_id, weeks_in_phase)
  
  n_pts <- n_distinct(sub$study_id)
  
  xmax_data <- max(sub$weeks_in_phase, na.rm = TRUE)  
  xmax_plot <- ceiling(xmax_data) + 0.5
  
  fit_df <- data.frame(x = seq(0, xmax_data, length.out = 100)) %>%  
    mutate(y = intercept + slope * x)
  
  p <- ggplot() +  
    # Spaghetti lines  
    geom_line(data = sub,  
              aes(x = weeks_in_phase, y = phq9_score, group = study_id),  
              color = line_color, linewidth = 0.65, alpha = 0.50) +  
    geom_point(data = sub,  
               aes(x = weeks_in_phase, y = phq9_score, group = study_id),  
               color = line_color, size = 2, alpha = 0.55,  
               shape = 21, fill = line_color, stroke = 0.3) +  
    # Model fit line — thick and prominent  
    geom_line(data = fit_df, aes(x = x, y = y),  
              color = fit_color, linewidth = 3) +  
    # Title  
    labs(  
      title = paste0(title_label, "\n(n = ", n_pts, " patients)"),  
      x = "Weeks Since Start of Phase"  
    ) +  
    # Manual fit-line legend  
    annotate("segment",  
             x = xmax_plot * 0.42, xend = xmax_plot * 0.56,  
             y = ymax * 0.96, yend = ymax * 0.96,  
             color = fit_color, linewidth = 2.5) +  
    annotate("text",  
             x = xmax_plot * 0.58, y = ymax * 0.96,  
             label = sprintf("model fit (%+.2f/wk)", slope),  
             hjust = 0, size = 3.5, fontface = "bold", color = fit_color) +  
    # Remission zone  
    annotate("rect", xmin = -0.3, xmax = xmax_plot, ymin = 0, ymax = 5,  
             fill = "#27AE60", alpha = 0.10) +  
    geom_hline(yintercept = 5, linetype = "dashed", color = "#27AE60",  
               linewidth = 0.9) +  
    annotate("text", x = 0, y = 2.5,  
             label = expression(italic("Remission (PHQ-9 < 5)")),  
             color = "#27AE60", size = 3.2, hjust = 0) +  
    # Scales  
    scale_x_continuous(limits = c(-0.3, xmax_plot),  
                       breaks = seq(0, floor(xmax_plot), 2)) +  
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 30, 5),  
                       expand = c(0, 0)) +  
    pub_theme +  
    theme(  
      plot.title = element_text(color = title_color)  
    )
  
  if (show_y) {  
    p <- p + labs(y = "PHQ-9 Total Score")  
  } else {  
    p <- p + labs(y = NULL)  
  }
  
  return(p)  
}

# ── Build the three panels ──  
p1 <- make_panel(g3, "stay_p1",  
                 "iTBS-Only (Stayers)", col_stay_fit,  
                 col_stay_line, col_stay_fit,  
                 int_stay_p1, slope_stay_p1,  
                 show_y = TRUE)

p2 <- make_panel(g3, "switch_p1",  
                 "Switchers \u00b7 Phase 1", col_sw1_fit,  
                 col_sw1_line, col_sw1_fit,  
                 int_switch_p1, slope_switch_p1)

p3 <- make_panel(g3, "switch_p2",  
                 "Switchers \u00b7 Phase 2", col_sw2_fit,  
                 col_sw2_line, col_sw2_fit,  
                 int_switch_p2, slope_switch_p2)

# ── Combine with patchwork ──  
combined <- p1 + p2 + p3 +  
  plot_layout(ncol = 3) +  
  plot_annotation(  
    title = "PHQ-9 Slope Trajectory by Protocol Group",  
    subtitle = sprintf(  
      paste0("Within switchers, decline accelerates after the switch ",  
             "(phase 2 \u2212 phase 1 %+.2f/wk, p = %.3f); ",  
             "stayers decline faster than switchers\u2019 phase 1 ",  
             "(between-subject %+.2f/wk, p = %.3f) and match phase 2"),  
      within_c["estimate"], within_c["p"],  
      btw_p1["estimate"], btw_p1["p"]  
    ),  
    theme = theme(  
      plot.title = element_text(face = "bold", size = 18, hjust = 0,  
                                color = "black", margin = margin(b = 4)),  
      plot.subtitle = element_text(size = 10.5, color = "gray35", hjust = 0,  
                                   margin = margin(b = 12),  
                                   lineheight = 1.2)  
    )  
  )

print(combined)

ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_slopes_3panel.png",  
       combined, width = 17, height = 6, dpi = 300)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_slopes_3panel.pdf",  
       combined, width = 17, height = 6)

cat("\n\u2713 Saved: phq9_slopes_3panel.png\n")  
cat("\u2713 Saved: phq9_slopes_3panel.pdf\n")

# ══════════════════════════════════════════════════════════════  
# 10. FULL MODEL SUMMARY OUTPUT  
# ══════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("              FULL MODEL SUMMARY\n")  
cat("══════════════════════════════════════════════════════════════\n")  
summary(m3)

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("              RESULTS TABLE\n")  
cat("══════════════════════════════════════════════════════════════\n")  
cat(sprintf(" Stayers slope:             %+.3f /week\n", slope_stay_p1))  
cat(sprintf(" Switchers Phase 1 slope:   %+.3f /week\n", slope_switch_p1))  
cat(sprintf(" Switchers Phase 2 slope:   %+.3f /week\n", slope_switch_p2))  
cat("--------------------------------------------------------------\n")  
cat(sprintf(" Within (P2-P1):   %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            within_c["estimate"], within_c["ci_lo"], within_c["ci_hi"],  
            within_c["p"], ifelse(within_c["p"] < 0.05, "*", "n.s.")))  
cat(sprintf(" Between (Stay-P1): %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            btw_p1["estimate"], btw_p1["ci_lo"], btw_p1["ci_hi"],  
            btw_p1["p"], ifelse(btw_p1["p"] < 0.05, "*", "n.s.")))  
cat(sprintf(" Between (Stay-P2): %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            btw_p2["estimate"], btw_p2["ci_lo"], btw_p2["ci_hi"],  
            btw_p2["p"], ifelse(btw_p2["p"] < 0.05, "*", "n.s.")))  
cat("══════════════════════════════════════════════════════════════\n")

cat("\nExpected poster values for comparison:\n")  
cat("  iTBS-only:       -0.85/wk\n")  
cat("  Switchers P1:    -0.25/wk\n")  
cat("  Switchers P2:    -0.82/wk\n")  
cat("  Within (P2-P1):  -0.57/wk, p = 0.012\n")  
cat("  Between (S-P1):  -0.60/wk, p = 0.022\n")

cat("\n\u2713 Analysis complete.\n")   