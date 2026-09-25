# ══════════════════════════════════════════════════════════════  
# 3-PANEL PHQ-9 SLOPE PLOT — PUBLICATION READY  
# X-axis: Weeks Since Phase Start | Capped at 36 sessions  
# Y-intercept: Actual baseline PHQ-9 (visit 0)  
# Matched switchers (same patients in both phases)  
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

# ── 2. Classify patients ──  
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

# ══════════════════════════════════════════════════════════════  
# 3. CAP AT 36 CUMULATIVE TMS SESSIONS  
#    Only keep switchers with BOTH pre and post data  
# ══════════════════════════════════════════════════════════════

dat_stay <- c1 %>%  
  dplyr::filter(study_id %in% stayer_ids) %>%  
  dplyr::filter(cumulative_tx <= 36) %>%  
  mutate(panel_group = "stay_p1")

c1_switch_capped <- c1 %>%  
  dplyr::filter(study_id %in% switcher_ids) %>%  
  dplyr::filter(cumulative_tx <= 36)

ids_with_pre <- c1_switch_capped %>%  
  dplyr::filter(phase == "pre") %>%  
  pull(study_id) %>%  
  unique()

ids_with_post <- c1_switch_capped %>%  
  dplyr::filter(phase == "post") %>%  
  pull(study_id) %>%  
  unique()

switcher_ids_both <- intersect(ids_with_pre, ids_with_post)

dropped_ids <- setdiff(ids_with_pre, ids_with_post)  
if (length(dropped_ids) > 0) {  
  cat(sprintf("\nDropped %d switcher(s) with no post-switch data after 36-session cap:\n",  
              length(dropped_ids)))  
  cat(paste("  ", dropped_ids, collapse = "\n"), "\n")  
}

c1_switch_capped <- c1_switch_capped %>%  
  dplyr::filter(study_id %in% switcher_ids_both)

dat_switch_p1 <- c1_switch_capped %>%  
  dplyr::filter(phase == "pre") %>%  
  mutate(panel_group = "switch_p1")

dat_switch_p2 <- c1_switch_capped %>%  
  dplyr::filter(phase == "post") %>%  
  mutate(panel_group = "switch_p2")

cat(sprintf("\nAfter capping and matching:\n"))  
cat(sprintf("  Stayers:      %d patients, %d observations\n",  
            n_distinct(dat_stay$study_id), nrow(dat_stay)))  
cat(sprintf("  Switchers P1: %d patients, %d observations\n",  
            n_distinct(dat_switch_p1$study_id), nrow(dat_switch_p1)))  
cat(sprintf("  Switchers P2: %d patients, %d observations\n",  
            n_distinct(dat_switch_p2$study_id), nrow(dat_switch_p2)))

g3 <- bind_rows(dat_stay, dat_switch_p1, dat_switch_p2)

# ══════════════════════════════════════════════════════════════  
# 4. COMPUTE WEEKS IN PHASE  
# ══════════════════════════════════════════════════════════════

g3 <- g3 %>%  
  mutate(phq9_date = as.POSIXct(phq9_date)) %>%  
  group_by(study_id, panel_group) %>%  
  mutate(  
    phase_start = min(phq9_date),  
    weeks_in_phase = as.numeric(difftime(phq9_date, phase_start, units = "weeks"))  
  ) %>%  
  ungroup()

# ══════════════════════════════════════════════════════════════  
# 5. BASELINE PHQ-9 (visit 0)  
# ══════════════════════════════════════════════════════════════

baseline_phq9 <- c1 %>%  
  dplyr::filter(visit == 0) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup() %>%  
  dplyr::select(study_id, baseline_phq9 = phq9_score)

g3 <- g3 %>%  
  left_join(baseline_phq9, by = "study_id")

mean_bl_stayers <- baseline_phq9 %>%  
  dplyr::filter(study_id %in% stayer_ids) %>%  
  pull(baseline_phq9) %>%  
  mean(na.rm = TRUE)

mean_bl_switchers <- baseline_phq9 %>%  
  dplyr::filter(study_id %in% switcher_ids_both) %>%  
  pull(baseline_phq9) %>%  
  mean(na.rm = TRUE)

cat(sprintf("\nMean baseline PHQ-9 (visit 0):\n"))  
cat(sprintf("  Stayers:   %.1f\n", mean_bl_stayers))  
cat(sprintf("  Switchers: %.1f\n", mean_bl_switchers))

# ── Set reference level ──  
g3$panel_group <- factor(g3$panel_group,  
                         levels = c("switch_p1", "stay_p1", "switch_p2"))

# ══════════════════════════════════════════════════════════════  
# 6. FIT HLM  
# ══════════════════════════════════════════════════════════════

cat("\nFitting hierarchical linear mixed model...\n")

m3 <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + baseline_phq9 +  
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
baseline_coef <- fe["baseline_phq9"]

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
n_switchers <- length(switcher_ids_both)

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

col_stay_line <- "#C6DEF1"  
col_stay_fit  <- "#2B5E8C"

col_sw1_line  <- "#FACBCB"  
col_sw1_fit   <- "#B83030"

col_sw2_line  <- "#F09090"  
col_sw2_fit   <- "#7A1515"

pub_theme <- theme_classic(base_size = 13) +  
  theme(  
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13,  
                                 margin = margin(b = 4)),  
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30",  
                                 margin = margin(b = 8)),  
    axis.title.x  = element_text(size = 11, face = "bold", color = "black",  
                                 margin = margin(t = 8)),  
    axis.title.y  = element_text(size = 11, face = "bold", color = "black",  
                                 margin = margin(r = 8)),  
    axis.text.x   = element_text(size = 10, color = "black"),  
    axis.text.y   = element_text(size = 10, color = "black"),  
    axis.line     = element_line(color = "black", linewidth = 0.5),  
    axis.ticks    = element_line(color = "black", linewidth = 0.4),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.margin   = margin(8, 14, 8, 8)  
  )

make_panel <- function(data, group_key, panel_label, title_label,  
                       title_color, line_color, fit_color,  
                       slope, mean_bl, show_y = FALSE) {
  
  sub <- data %>%  
    dplyr::filter(panel_group == group_key) %>%  
    arrange(study_id, weeks_in_phase)
  
  n_pts <- n_distinct(sub$study_id)
  
  xmax_data <- max(sub$weeks_in_phase, na.rm = TRUE)  
  xmax_plot <- ceiling(xmax_data) + 0.5
  
  fit_df <- data.frame(x = seq(0, xmax_data, length.out = 100)) %>%  
    mutate(y = mean_bl + slope * x)  
  fit_df$y <- pmax(fit_df$y, 0)
  
  p <- ggplot() +  
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 5,  
             fill = "#E8F5E9", alpha = 0.5) +  
    geom_hline(yintercept = 5, linetype = "dashed", color = "#66BB6A",  
               linewidth = 0.6, alpha = 0.7) +  
    geom_line(data = sub,  
              aes(x = weeks_in_phase, y = phq9_score, group = study_id),  
              color = line_color, linewidth = 0.5, alpha = 0.35) +  
    geom_point(data = sub,  
               aes(x = weeks_in_phase, y = phq9_score, group = study_id),  
               color = line_color, size = 1.5, alpha = 0.40,  
               shape = 21, fill = line_color, stroke = 0.2) +  
    geom_line(data = fit_df, aes(x = x, y = y),  
              color = fit_color, linewidth = 2.5) +  
    annotate("text", x = -0.1, y = 27, label = panel_label,  
             hjust = 0, vjust = 1, size = 5, fontface = "bold") +  
    labs(  
      title = title_label,  
      subtitle = paste0("n = ", n_pts),  
      x = "Weeks Since Phase Start"  
    ) +  
    annotate("segment",  
             x = 0.1, xend = xmax_plot * 0.18,  
             y = 24.5, yend = 24.5,  
             color = fit_color, linewidth = 2) +  
    annotate("text",  
             x = xmax_plot * 0.20, y = 24.5,  
             label = sprintf("%+.2f pts/wk", slope),  
             hjust = 0, size = 3.5, fontface = "bold", color = fit_color) +  
    annotate("text", x = xmax_plot * 0.97, y = 3,  
             label = "remission",  
             color = "#66BB6A", size = 2.8, hjust = 1,  
             fontface = "italic", alpha = 0.7) +  
    scale_x_continuous(limits = c(-0.3, xmax_plot),  
                       breaks = seq(0, floor(xmax_plot), 2),  
                       expand = c(0, 0)) +  
    scale_y_continuous(limits = c(0, 27), breaks = seq(0, 25, 5),  
                       expand = c(0, 0.3)) +  
    pub_theme +  
    theme(plot.title = element_text(color = title_color))
  
  if (show_y) {  
    p <- p + labs(y = "PHQ-9 Score")  
  } else {  
    p <- p + labs(y = NULL)  
  }
  
  return(p)  
}

p1 <- make_panel(g3, "stay_p1", "A",  
                 "iTBS-Only", col_stay_fit,  
                 col_stay_line, col_stay_fit,  
                 slope_stay_p1, mean_bl_stayers,  
                 show_y = TRUE)

p2 <- make_panel(g3, "switch_p1", "B",  
                 "Switchers - Phase 1 (iTBS)", col_sw1_fit,  
                 col_sw1_line, col_sw1_fit,  
                 slope_switch_p1, mean_bl_switchers)

p3 <- make_panel(g3, "switch_p2", "C",  
                 "Switchers - Phase 2 (Bilateral)", col_sw2_fit,  
                 col_sw2_line, col_sw2_fit,  
                 slope_switch_p2, mean_bl_switchers)

combined <- p1 + p2 + p3 +  
  plot_layout(ncol = 3) +  
  plot_annotation(  
    title = "PHQ-9 Trajectory by Treatment Protocol",  
    theme = theme(  
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5,  
                                color = "black", margin = margin(b = 8))  
    )  
  )

print(combined)

ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_slopes_3panel.png",  
       combined, width = 16, height = 5.5, dpi = 300)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_slopes_3panel.pdf",  
       combined, width = 16, height = 5.5)

cat("\n\u2713 Saved: phq9_slopes_3panel.png\n")  
cat("\u2713 Saved: phq9_slopes_3panel.pdf\n")

# ══════════════════════════════════════════════════════════════  
# 10. FULL MODEL SUMMARY  
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
cat(sprintf(" Baseline PHQ-9 coeff:      %+.3f\n", baseline_coef))  
cat("--------------------------------------------------------------\n")  
cat(sprintf(" Within (P2-P1):    %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            within_c["estimate"], within_c["ci_lo"], within_c["ci_hi"],  
            within_c["p"], ifelse(within_c["p"] < 0.05, "*", "n.s.")))  
cat(sprintf(" Between (Stay-P1): %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            btw_p1["estimate"], btw_p1["ci_lo"], btw_p1["ci_hi"],  
            btw_p1["p"], ifelse(btw_p1["p"] < 0.05, "*", "n.s.")))  
cat(sprintf(" Between (Stay-P2): %+.3f /wk [%+.3f, %+.3f] p=%.4f %s\n",  
            btw_p2["estimate"], btw_p2["ci_lo"], btw_p2["ci_hi"],  
            btw_p2["p"], ifelse(btw_p2["p"] < 0.05, "*", "n.s.")))  
cat("══════════════════════════════════════════════════════════════\n")

cat("\nSession cap verification:\n")  
cat(sprintf("  Max cumulative_tx in stayers:      %d\n",  
            max(dat_stay$cumulative_tx, na.rm = TRUE)))  
cat(sprintf("  Max cumulative_tx in switchers P1: %d\n",  
            max(dat_switch_p1$cumulative_tx, na.rm = TRUE)))  
cat(sprintf("  Max cumulative_tx in switchers P2: %d\n",  
            max(dat_switch_p2$cumulative_tx, na.rm = TRUE)))

cat(sprintf("\nMean baseline PHQ-9 (visit 0):\n"))  
cat(sprintf("  Stayers:              %.1f\n", mean_bl_stayers))  
cat(sprintf("  Switchers (P1 & P2):  %.1f\n", mean_bl_switchers))

cat(sprintf("\nPatient counts (matched):\n"))  
cat(sprintf("  Stayers:      %d\n", n_distinct(dat_stay$study_id)))  
cat(sprintf("  Switchers P1: %d\n", n_distinct(dat_switch_p1$study_id)))  
cat(sprintf("  Switchers P2: %d\n", n_distinct(dat_switch_p2$study_id)))

cat("\nWeeks in phase range per panel:\n")  
for (pg in c("stay_p1", "switch_p1", "switch_p2")) {  
  wk_range <- g3 %>% dplyr::filter(panel_group == pg) %>% pull(weeks_in_phase)  
  cat(sprintf("  %s: %.1f to %.1f weeks\n", pg,  
              min(wk_range, na.rm = TRUE), max(wk_range, na.rm = TRUE)))  
}

cat("\n\u2713 Analysis complete.\n")  