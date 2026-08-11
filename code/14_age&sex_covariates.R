# ══════════════════════════════════════════════════════════════  
# SENSITIVITY ANALYSIS: 3-PANEL PHQ-9 SLOPE PLOT  
# Controlling for Age and Sex  
# ══════════════════════════════════════════════════════════════

library(tidyverse)  
library(readxl)  
library(lme4)  
library(lmerTest)  
library(patchwork)

# ── 1. Load both files ──  
long_file <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx"  
demo_file <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_summary.xlsx"

raw <- read_excel(long_file, sheet = "Sheet1")  
colnames(raw) <- trimws(colnames(raw))

demo_raw <- read_excel(demo_file, sheet = "Sheet1")  
colnames(demo_raw) <- trimws(colnames(demo_raw))

# ── 2. Extract demographics (one row per patient) ──  
demo <- demo_raw %>%  
  select(study_id, age_at_first_phq, sex) %>%  
  mutate(  
    age = as.numeric(age_at_first_phq),  
    sex_numeric = as.numeric(sex)  
  ) %>%  
  filter(!is.na(age), !is.na(sex_numeric))

cat("Patients with demographics:", nrow(demo), "\n")  
cat("Age range:", min(demo$age), "-", max(demo$age), "\n")  
cat("Sex: Male(1) n=", sum(demo$sex_numeric == 1),  
    ", Female(2) n=", sum(demo$sex_numeric == 2), "\n\n")

# ── 3. Prepare longitudinal data ──  
c1 <- raw %>%  
  filter(!is.na(phq9_score), !is.na(phq9_date)) %>%  
  mutate(  
    phq9_score = as.numeric(phq9_score),  
    phq9_date = as.POSIXct(phq9_date)  
  ) %>%  
  arrange(study_id, phq9_date)

stayer_ids <- c1 %>% filter(group == "iTBS_only") %>% pull(study_id) %>% unique()  
switcher_ids <- c1 %>% filter(group == "switcher") %>% pull(study_id) %>% unique()

# ── 4. Build 3-panel dataset ──  
g3 <- bind_rows(  
  c1 %>% filter(study_id %in% stayer_ids) %>% mutate(panel_group = "stay_p1"),  
  c1 %>% filter(study_id %in% switcher_ids, phase == "pre") %>% mutate(panel_group = "switch_p1"),  
  c1 %>% filter(study_id %in% switcher_ids, phase == "post") %>% mutate(panel_group = "switch_p2")  
)

# ── 5. Merge demographics from analytic_sample_summary.xlsx ──  
g3 <- g3 %>%  
  left_join(demo %>% select(study_id, age, sex_numeric), by = "study_id") %>%  
  filter(!is.na(age), !is.na(sex_numeric))

cat("Final dataset:", nrow(g3), "observations\n")  
cat("  Stayers:", n_distinct(g3$study_id[g3$panel_group == "stay_p1"]), "\n")  
cat("  Switchers:", n_distinct(g3$study_id[g3$panel_group %in% c("switch_p1", "switch_p2")]), "\n\n")

# ── 6. Compute weeks_in_phase, center age ──  
g3 <- g3 %>%  
  group_by(study_id, panel_group) %>%  
  mutate(  
    phase_start = min(phq9_date),  
    weeks_in_phase = as.numeric(difftime(phq9_date, phase_start, units = "weeks"))  
  ) %>%  
  ungroup() %>%  
  mutate(age_centered = age - mean(age, na.rm = TRUE))

g3$panel_group <- factor(g3$panel_group, levels = c("switch_p1", "stay_p1", "switch_p2"))

# ══════════════════════════════════════════════════════════════  
# 7. FIT BOTH MODELS  
# ══════════════════════════════════════════════════════════════

cat("Fitting UNADJUSTED model...\n")  
m_unadj <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3, REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)

cat("Fitting ADJUSTED model (age + sex)...\n")  
m_adj <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + age_centered + sex_numeric +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3, REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)

# ══════════════════════════════════════════════════════════════  
# 8. EXTRACT RESULTS  
# ══════════════════════════════════════════════════════════════

extract_results <- function(model, label) {  
  fe <- fixef(model)  
  vcov_fe <- vcov(model)  
  fe_names <- names(fe)
  
  slope_stay <- fe["weeks_in_phase"] + fe["weeks_in_phase:panel_groupstay_p1"]  
  slope_sw1  <- fe["weeks_in_phase"]  
  slope_sw2  <- fe["weeks_in_phase"] + fe["weeks_in_phase:panel_groupswitch_p2"]
  
  int_stay <- fe["(Intercept)"] + fe["panel_groupstay_p1"]  
  int_sw1  <- fe["(Intercept)"]  
  int_sw2  <- fe["(Intercept)"] + fe["panel_groupswitch_p2"]
  
  get_contrast <- function(pos, neg = NULL) {  
    r <- matrix(0, nrow = 1, ncol = length(fe))  
    colnames(r) <- fe_names  
    r[1, pos] <- 1  
    if (!is.null(neg)) r[1, neg] <- -1  
    est  <- as.numeric(r %*% fe)  
    se   <- as.numeric(sqrt(r %*% vcov_fe %*% t(r)))  
    pval <- 2 * pt(abs(est / se), df = df.residual(model), lower.tail = FALSE)  
    c(estimate = est, ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se, p = pval)  
  }
  
  within_c <- get_contrast("weeks_in_phase:panel_groupswitch_p2")  
  btw_p1   <- get_contrast("weeks_in_phase:panel_groupstay_p1")  
  btw_p2   <- get_contrast("weeks_in_phase:panel_groupstay_p1",  
                           "weeks_in_phase:panel_groupswitch_p2")
  
  cat(sprintf("\n── %s ──\n", label))  
  cat(sprintf("  Stayers slope:           %+.3f /week\n", slope_stay))  
  cat(sprintf("  Switchers Phase 1 slope: %+.3f /week\n", slope_sw1))  
  cat(sprintf("  Switchers Phase 2 slope: %+.3f /week\n", slope_sw2))  
  cat("  ──────────────────────────────────────────────\n")  
  cat(sprintf("  Within (P2-P1):    %+.3f [%+.3f, %+.3f] p=%.4f %s\n",  
              within_c["estimate"], within_c["ci_lo"], within_c["ci_hi"],  
              within_c["p"], ifelse(within_c["p"] < 0.05, "*", "n.s.")))  
  cat(sprintf("  Between (Stay-P1): %+.3f [%+.3f, %+.3f] p=%.4f %s\n",  
              btw_p1["estimate"], btw_p1["ci_lo"], btw_p1["ci_hi"],  
              btw_p1["p"], ifelse(btw_p1["p"] < 0.05, "*", "n.s.")))  
  cat(sprintf("  Between (Stay-P2): %+.3f [%+.3f, %+.3f] p=%.4f %s\n",  
              btw_p2["estimate"], btw_p2["ci_lo"], btw_p2["ci_hi"],  
              btw_p2["p"], ifelse(btw_p2["p"] < 0.05, "*", "n.s.")))
  
  list(  
    slopes = c(stay = slope_stay, switch_p1 = slope_sw1, switch_p2 = slope_sw2),  
    intercepts = c(stay = int_stay, switch_p1 = int_sw1, switch_p2 = int_sw2),  
    contrasts = list(within = within_c, btw_p1 = btw_p1, btw_p2 = btw_p2)  
  )  
}

res_unadj <- extract_results(m_unadj, "UNADJUSTED MODEL")  
res_adj   <- extract_results(m_adj, "ADJUSTED MODEL (age + sex)")

# Covariate effects  
cat("\n── Covariate Effects ──\n")  
adj_coefs <- summary(m_adj)$coefficients  
cat(sprintf("  Age:  B = %+.3f, p = %.4f\n",  
            adj_coefs["age_centered", "Estimate"], adj_coefs["age_centered", "Pr(>|t|)"]))  
cat(sprintf("  Sex:  B = %+.3f, p = %.4f\n",  
            adj_coefs["sex_numeric", "Estimate"], adj_coefs["sex_numeric", "Pr(>|t|)"]))

# ══════════════════════════════════════════════════════════════  
# 9. COMPARISON TABLE  
# ══════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("          SENSITIVITY ANALYSIS: SIDE-BY-SIDE\n")  
cat("══════════════════════════════════════════════════════════════\n\n")  
cat(sprintf("%-25s %15s %15s\n", "", "Unadjusted", "Adjusted"))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Stayer slope",  
            res_unadj$slopes["stay"], res_adj$slopes["stay"]))  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Switcher P1 slope",  
            res_unadj$slopes["switch_p1"], res_adj$slopes["switch_p1"]))  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Switcher P2 slope",  
            res_unadj$slopes["switch_p2"], res_adj$slopes["switch_p2"]))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Within (P2-P1)",  
            res_unadj$contrasts$within["p"], res_adj$contrasts$within["p"]))  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Between (Stay-P1)",  
            res_unadj$contrasts$btw_p1["p"], res_adj$contrasts$btw_p1["p"]))  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Between (Stay-P2)",  
            res_unadj$contrasts$btw_p2["p"], res_adj$contrasts$btw_p2["p"]))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("AIC:                    %10.1f     %10.1f\n", AIC(m_unadj), AIC(m_adj)))

# ══════════════════════════════════════════════════════════════  
# 10. 3-PANEL PLOT (ADJUSTED MODEL)  
# ══════════════════════════════════════════════════════════════

col_stay_line <- "#C6DEF1"; col_stay_fit <- "#2B5E8C"  
col_sw1_line  <- "#FACBCB"; col_sw1_fit  <- "#B83030"  
col_sw2_line  <- "#F09090"; col_sw2_fit  <- "#7A1515"  
ymax <- max(g3$phq9_score, na.rm = TRUE) + 1

pub_theme <- theme_classic(base_size = 14) +  
  theme(  
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, margin = margin(b = 6)),  
    axis.title = element_text(size = 12, face = "bold", color = "black"),  
    axis.text = element_text(size = 11, face = "bold", color = "black"),  
    axis.line = element_line(color = "black", linewidth = 0.7),  
    plot.margin = margin(12, 18, 12, 12)  
  )

make_panel <- function(data, grp, title, title_col, line_col, fit_col,  
                       intercept, slope, show_y = FALSE) {  
  sub <- data %>% filter(panel_group == grp) %>% arrange(study_id, weeks_in_phase)  
  n_pts <- n_distinct(sub$study_id)  
  xmax <- ceiling(max(sub$weeks_in_phase, na.rm = TRUE)) + 0.5  
  fit_df <- tibble(x = seq(0, max(sub$weeks_in_phase), length.out = 100),  
                   y = intercept + slope * x)
  
  p <- ggplot() +  
    geom_line(data = sub, aes(x = weeks_in_phase, y = phq9_score, group = study_id),  
              color = line_col, linewidth = 0.65, alpha = 0.50) +  
    geom_point(data = sub, aes(x = weeks_in_phase, y = phq9_score),  
               color = line_col, size = 2, alpha = 0.55, shape = 21, fill = line_col) +  
    geom_line(data = fit_df, aes(x = x, y = y), color = fit_col, linewidth = 3) +  
    annotate("rect", xmin = -0.3, xmax = xmax, ymin = 0, ymax = 5,  
             fill = "#27AE60", alpha = 0.10) +  
    geom_hline(yintercept = 5, linetype = "dashed", color = "#27AE60", linewidth = 0.9) +  
    annotate("text", x = 0, y = 2.5, label = expression(italic("Remission (PHQ-9 < 5)")),  
             color = "#27AE60", size = 3.2, hjust = 0) +  
    annotate("text", x = xmax * 0.5, y = ymax * 0.96,  
             label = sprintf("slope = %+.2f/wk", slope),  
             hjust = 0, size = 3.5, fontface = "bold", color = fit_col) +  
    scale_x_continuous(limits = c(-0.3, xmax), breaks = seq(0, floor(xmax), 2)) +  
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 30, 5), expand = c(0, 0)) +  
    labs(title = paste0(title, "\n(n = ", n_pts, ")"), x = "Weeks Since Phase Start") +  
    pub_theme +  
    theme(plot.title = element_text(color = title_col))
  
  if (show_y) p <- p + labs(y = "PHQ-9 Total Score") else p <- p + labs(y = NULL)  
  p  
}

p1 <- make_panel(g3, "stay_p1", "iTBS-Only (Stayers)", col_stay_fit,  
                 col_stay_line, col_stay_fit,  
                 res_adj$intercepts["stay"], res_adj$slopes["stay"], show_y = TRUE)  
p2 <- make_panel(g3, "switch_p1", "Switchers \u00b7 Phase 1 (iTBS)", col_sw1_fit,  
                 col_sw1_line, col_sw1_fit,  
                 res_adj$intercepts["switch_p1"], res_adj$slopes["switch_p1"])  
p3 <- make_panel(g3, "switch_p2", "Switchers \u00b7 Phase 2 (Bilateral)", col_sw2_fit,  
                 col_sw2_line, col_sw2_fit,  
                 res_adj$intercepts["switch_p2"], res_adj$slopes["switch_p2"])

combined <- p1 + p2 + p3 +  
  plot_layout(ncol = 3) +  
  plot_annotation(  
    title = "PHQ-9 Trajectories by Protocol Group (Adjusted for Age & Sex)",  
    subtitle = sprintf(  
      "Within switchers (P2-P1): %+.2f/wk (p=%.3f) | Stayers vs Switchers P1: %+.2f/wk (p=%.3f)",  
      res_adj$contrasts$within["estimate"], res_adj$contrasts$within["p"],  
      res_adj$contrasts$btw_p1["estimate"], res_adj$contrasts$btw_p1["p"]  
    ),  
    theme = theme(  
      plot.title = element_text(face = "bold", size = 18, hjust = 0),  
      plot.subtitle = element_text(size = 10.5, color = "gray35", hjust = 0, margin = margin(b = 12))  
    )  
  )

print(combined)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/phq9_slopes_3panel_adjusted.png",  
       combined, width = 17, height = 6, dpi = 300)

# Full summaries  
cat("\n══ UNADJUSTED ══\n")  
print(summary(m_unadj))  
cat("\n══ ADJUSTED ══\n")  
print(summary(m_adj))  
cat("\n══ MODEL COMPARISON ══\n")  
print(anova(m_unadj, m_adj))

cat("\nIf p-values similar -> results ROBUST to age/sex.\n")  
cat("Saved: phq9_slopes_3panel_adjusted.png\n")  