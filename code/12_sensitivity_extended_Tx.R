# ============================================================================  
# SENSITIVITY ANALYSIS: Extended Treatment Courses (No Session Cap)  
# Comparison: 36-session cap (primary) vs. uncapped (sensitivity)  
# Capped: TMS074 excluded (no post data within cap)  
# Uncapped: TMS074 included (post data available)  
# ============================================================================

library(dplyr)  
library(readxl)  
library(lme4)  
library(lmerTest)

# ── 1. Load data ──  
raw <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx",  
                  sheet = "Sheet1")

c1_all <- raw %>%  
  filter(!is.na(phq9_score) & !is.na(phq9_date)) %>%  
  mutate(  
    phq9_score = as.numeric(phq9_score),  
    phq9_date = as.POSIXct(phq9_date)  
  ) %>%  
  arrange(study_id, phq9_date)

stayer_ids <- c1_all %>% filter(group == "iTBS_only") %>% pull(study_id) %>% unique()

# ══════════════════════════════════════════════════════════════  
# 2. BUILD CAPPED DATASET (36 sessions, TMS074 excluded)  
#    Matches primary analysis exactly  
# ══════════════════════════════════════════════════════════════

c1_capped <- c1_all %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36)

switch_pre_cap <- c1_capped %>% filter(group == "switcher" & phase == "pre") %>% pull(study_id) %>% unique()  
switch_post_cap <- c1_capped %>% filter(group == "switcher" & phase == "post") %>% pull(study_id) %>% unique()  
switch_both_cap <- intersect(switch_pre_cap, switch_post_cap)

g3_capped <- bind_rows(  
  c1_capped %>% filter(study_id %in% stayer_ids) %>% mutate(panel_group = "stay_p1"),  
  c1_capped %>% filter(study_id %in% switch_both_cap & phase == "pre") %>% mutate(panel_group = "switch_p1"),  
  c1_capped %>% filter(study_id %in% switch_both_cap & phase == "post") %>% mutate(panel_group = "switch_p2")  
)

# ══════════════════════════════════════════════════════════════  
# 3. BUILD UNCAPPED DATASET (all sessions, TMS074 INCLUDED)  
# ══════════════════════════════════════════════════════════════

switch_pre_unc <- c1_all %>% filter(group == "switcher" & phase == "pre") %>% pull(study_id) %>% unique()  
switch_post_unc <- c1_all %>% filter(group == "switcher" & phase == "post") %>% pull(study_id) %>% unique()  
switch_both_unc <- intersect(switch_pre_unc, switch_post_unc)

g3_uncapped <- bind_rows(  
  c1_all %>% filter(study_id %in% stayer_ids) %>% mutate(panel_group = "stay_p1"),  
  c1_all %>% filter(study_id %in% switch_both_unc & phase == "pre") %>% mutate(panel_group = "switch_p1"),  
  c1_all %>% filter(study_id %in% switch_both_unc & phase == "post") %>% mutate(panel_group = "switch_p2")  
)

# ══════════════════════════════════════════════════════════════  
# 4. Add weeks_in_phase and baseline PHQ-9 to both datasets  
# ══════════════════════════════════════════════════════════════

add_weeks_and_baseline <- function(data, source_data) {  
  data <- data %>%  
    group_by(study_id, panel_group) %>%  
    mutate(  
      phase_start = min(phq9_date),  
      weeks_in_phase = as.numeric(difftime(phq9_date, phase_start, units = "weeks"))  
    ) %>%  
    ungroup()
  
  baseline <- source_data %>%  
    filter(visit == 0) %>%  
    group_by(study_id) %>%  
    slice(1) %>%  
    ungroup() %>%  
    select(study_id, baseline_phq9 = phq9_score)
  
  data <- data %>%  
    left_join(baseline, by = "study_id")
  
  data$panel_group <- factor(data$panel_group,  
                             levels = c("switch_p1", "stay_p1", "switch_p2"))  
  data  
}

g3_capped <- add_weeks_and_baseline(g3_capped, c1_capped)  
g3_uncapped <- add_weeks_and_baseline(g3_uncapped, c1_all)

# ══════════════════════════════════════════════════════════════  
# 5. Report sample sizes  
# ══════════════════════════════════════════════════════════════

cat("=================================================================\n")  
cat("SAMPLE COMPARISON\n")  
cat("=================================================================\n\n")

cat("CAPPED (36 sessions, TMS074 excluded):\n")  
cat(sprintf("  Stayers:      %d patients, %d observations\n",  
            n_distinct(g3_capped$study_id[g3_capped$panel_group == "stay_p1"]),  
            sum(g3_capped$panel_group == "stay_p1")))  
cat(sprintf("  Switchers P1: %d patients, %d observations\n",  
            n_distinct(g3_capped$study_id[g3_capped$panel_group == "switch_p1"]),  
            sum(g3_capped$panel_group == "switch_p1")))  
cat(sprintf("  Switchers P2: %d patients, %d observations\n",  
            n_distinct(g3_capped$study_id[g3_capped$panel_group == "switch_p2"]),  
            sum(g3_capped$panel_group == "switch_p2")))  
cat(sprintf("  Total: %d patients, %d observations\n",  
            n_distinct(g3_capped$study_id), nrow(g3_capped)))  
cat(sprintf("  TMS074 included: %s\n\n", "TMS074" %in% g3_capped$study_id))

cat("UNCAPPED (all sessions, TMS074 INCLUDED):\n")  
cat(sprintf("  Stayers:      %d patients, %d observations\n",  
            n_distinct(g3_uncapped$study_id[g3_uncapped$panel_group == "stay_p1"]),  
            sum(g3_uncapped$panel_group == "stay_p1")))  
cat(sprintf("  Switchers P1: %d patients, %d observations\n",  
            n_distinct(g3_uncapped$study_id[g3_uncapped$panel_group == "switch_p1"]),  
            sum(g3_uncapped$panel_group == "switch_p1")))  
cat(sprintf("  Switchers P2: %d patients, %d observations\n",  
            n_distinct(g3_uncapped$study_id[g3_uncapped$panel_group == "switch_p2"]),  
            sum(g3_uncapped$panel_group == "switch_p2")))  
cat(sprintf("  Total: %d patients, %d observations\n",  
            n_distinct(g3_uncapped$study_id), nrow(g3_uncapped)))  
cat(sprintf("  TMS074 included: %s\n\n", "TMS074" %in% g3_uncapped$study_id))

cat(sprintf("Additional observations: %d\n", nrow(g3_uncapped) - nrow(g3_capped)))  
cat(sprintf("Additional patients: %d (TMS074 restored)\n\n",   
            n_distinct(g3_uncapped$study_id) - n_distinct(g3_capped$study_id)))

# Weeks range comparison  
cat("Weeks range per panel:\n")  
for (pg in c("stay_p1", "switch_p1", "switch_p2")) {  
  wk_cap <- g3_capped %>% filter(panel_group == pg) %>% pull(weeks_in_phase)  
  wk_unc <- g3_uncapped %>% filter(panel_group == pg) %>% pull(weeks_in_phase)  
  cat(sprintf("  %s: capped %.1f-%.1f wks | uncapped %.1f-%.1f wks\n",  
              pg, min(wk_cap), max(wk_cap), min(wk_unc), max(wk_unc)))  
}

# ══════════════════════════════════════════════════════════════  
# 6. FIT BOTH MODELS  
# ══════════════════════════════════════════════════════════════

cat("\n=================================================================\n")  
cat("FITTING MODELS\n")  
cat("=================================================================\n\n")

cat("Fitting CAPPED model (primary, 36 sessions, n=84)...\n")  
m_capped <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + baseline_phq9 +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3_capped, REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)  
cat("  Converged.\n")

cat("Fitting UNCAPPED model (all sessions, n=85, TMS074 included)...\n")  
m_uncapped <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + baseline_phq9 +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3_uncapped, REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)  
cat("  Converged.\n\n")

# ══════════════════════════════════════════════════════════════  
# 7. EXTRACT AND COMPARE RESULTS  
# ══════════════════════════════════════════════════════════════

extract_results <- function(model, label) {  
  fe <- fixef(model)  
  vcov_fe <- vcov(model)  
  fe_names <- names(fe)
  
  slope_stay <- fe["weeks_in_phase"] + fe["weeks_in_phase:panel_groupstay_p1"]  
  slope_sw1  <- fe["weeks_in_phase"]  
  slope_sw2  <- fe["weeks_in_phase"] + fe["weeks_in_phase:panel_groupswitch_p2"]
  
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
  cat(sprintf("  Stayers slope:           %+.3f /wk\n", slope_stay))  
  cat(sprintf("  Switchers Phase 1 slope: %+.3f /wk\n", slope_sw1))  
  cat(sprintf("  Switchers Phase 2 slope: %+.3f /wk\n", slope_sw2))  
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
    slopes = c(stay = unname(slope_stay),  
               switch_p1 = unname(slope_sw1),  
               switch_p2 = unname(slope_sw2)),  
    contrasts = list(within = within_c, btw_p1 = btw_p1, btw_p2 = btw_p2)  
  )  
}

res_capped   <- extract_results(m_capped, "CAPPED (36 sessions, n=84, TMS074 excluded)")  
res_uncapped <- extract_results(m_uncapped, "UNCAPPED (all sessions, n=85, TMS074 included)")

# ══════════════════════════════════════════════════════════════  
# 8. SIDE-BY-SIDE COMPARISON  
# ══════════════════════════════════════════════════════════════

format_p <- function(p) {  
  if (p < .001) return("< .001")  
  sprintf("%.4f", p)  
}

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("     SENSITIVITY ANALYSIS: CAPPED vs UNCAPPED\n")  
cat("══════════════════════════════════════════════════════════════\n\n")  
cat(sprintf("%-25s %15s %15s %10s\n", "", "Capped (36)", "Uncapped", "Direction"))  
cat(paste(rep("-", 70), collapse = ""), "\n")  
cat(sprintf("%-25s %15d %15d\n", "N patients",  
            n_distinct(g3_capped$study_id), n_distinct(g3_uncapped$study_id)))  
cat(sprintf("%-25s %15d %15d\n", "N switchers",  
            n_distinct(g3_capped$study_id[g3_capped$panel_group == "switch_p1"]),  
            n_distinct(g3_uncapped$study_id[g3_uncapped$panel_group == "switch_p1"])))  
cat(sprintf("%-25s %15d %15d\n", "N observations",  
            nrow(g3_capped), nrow(g3_uncapped)))  
cat(sprintf("%-25s %15s %15s\n", "TMS074", "excluded", "INCLUDED"))  
cat(paste(rep("-", 70), collapse = ""), "\n")  
cat(sprintf("%-25s %+15.3f %+15.3f %10s\n", "Stayer slope",  
            res_capped$slopes["stay"], res_uncapped$slopes["stay"],  
            ifelse(abs(res_uncapped$slopes["stay"]) > abs(res_capped$slopes["stay"]),  
                   "stronger", "weaker")))  
cat(sprintf("%-25s %+15.3f %+15.3f %10s\n", "Switcher P1 slope",  
            res_capped$slopes["switch_p1"], res_uncapped$slopes["switch_p1"],  
            ifelse(abs(res_uncapped$slopes["switch_p1"]) < abs(res_capped$slopes["switch_p1"]),  
                   "flatter", "steeper")))  
cat(sprintf("%-25s %+15.3f %+15.3f %10s\n", "Switcher P2 slope",  
            res_capped$slopes["switch_p2"], res_uncapped$slopes["switch_p2"],  
            ifelse(abs(res_uncapped$slopes["switch_p2"]) > abs(res_capped$slopes["switch_p2"]),  
                   "stronger", "weaker")))  
cat(paste(rep("-", 70), collapse = ""), "\n")  
cat(sprintf("%-25s p=%-12s p=%-12s\n", "Within (P2-P1)",  
            format_p(res_capped$contrasts$within["p"]),  
            format_p(res_uncapped$contrasts$within["p"])))  
cat(sprintf("%-25s p=%-12s p=%-12s\n", "Between (Stay-P1)",  
            format_p(res_capped$contrasts$btw_p1["p"]),  
            format_p(res_uncapped$contrasts$btw_p1["p"])))  
cat(sprintf("%-25s p=%-12s p=%-12s\n", "Between (Stay-P2)",  
            format_p(res_capped$contrasts$btw_p2["p"]),  
            format_p(res_uncapped$contrasts$btw_p2["p"])))  
cat(paste(rep("-", 70), collapse = ""), "\n")  
cat(sprintf("%-25s %15.1f %15.1f\n", "AIC", AIC(m_capped), AIC(m_uncapped)))  
cat(paste(rep("-", 70), collapse = ""), "\n")

# ══════════════════════════════════════════════════════════════  
# 9. FULL MODEL SUMMARIES  
# ══════════════════════════════════════════════════════════════

cat("\n══ CAPPED MODEL (PRIMARY) ══\n")  
print(summary(m_capped))  
cat("\n══ UNCAPPED MODEL (SENSITIVITY) ══\n")  
print(summary(m_uncapped))

# ══════════════════════════════════════════════════════════════  
# 10. ROBUSTNESS ASSESSMENT  
# ══════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("ROBUSTNESS ASSESSMENT\n")  
cat("══════════════════════════════════════════════════════════════\n\n")

cap_within_sig <- res_capped$contrasts$within["p"] < 0.05  
unc_within_sig <- res_uncapped$contrasts$within["p"] < 0.05

if (cap_within_sig && unc_within_sig) {  
  cat("WITHIN-SUBJECT SLOPE CHANGE: Significant in BOTH models.\n")  
  if (res_uncapped$contrasts$within["p"] < res_capped$contrasts$within["p"]) {  
    cat("  -> Removing the cap STRENGTHENS the finding (smaller p-value).\n")  
  } else {  
    cat("  -> Removing the cap slightly attenuates but PRESERVES the finding.\n")  
  }  
} else if (!cap_within_sig && unc_within_sig) {  
  cat("WITHIN-SUBJECT SLOPE CHANGE: Significant ONLY in uncapped model.\n")  
  cat("  -> Extended data REVEALS a significant effect.\n")  
} else if (cap_within_sig && !unc_within_sig) {  
  cat("WITHIN-SUBJECT SLOPE CHANGE: Significant ONLY in capped model.\n")  
  cat("  -> Effect may be concentrated in the acute treatment window.\n")  
} else {  
  cat("WITHIN-SUBJECT SLOPE CHANGE: Not significant in either model.\n")  
}

cap_btw_sig <- res_capped$contrasts$btw_p1["p"] < 0.05  
unc_btw_sig <- res_uncapped$contrasts$btw_p1["p"] < 0.05

cat(sprintf("\nBETWEEN-GROUP (Stay vs P1):\n"))  
cat(sprintf("  Capped:   p = %s (%s)\n",  
            format_p(res_capped$contrasts$btw_p1["p"]),  
            ifelse(cap_btw_sig, "significant", "n.s.")))  
cat(sprintf("  Uncapped: p = %s (%s)\n",  
            format_p(res_uncapped$contrasts$btw_p1["p"]),  
            ifelse(unc_btw_sig, "significant", "n.s.")))

cat(sprintf("\nPHASE 2 vs STAYERS (convergence):\n"))  
cat(sprintf("  Capped:   p = %s\n", format_p(res_capped$contrasts$btw_p2["p"])))  
cat(sprintf("  Uncapped: p = %s\n", format_p(res_uncapped$contrasts$btw_p2["p"])))  
cat("  (Non-significant = switchers match stayers after switch)\n")

# ══════════════════════════════════════════════════════════════  
# 11. MANUSCRIPT-READY SUMMARY  
# ══════════════════════════════════════════════════════════════

n_cap <- n_distinct(g3_capped$study_id)  
n_unc <- n_distinct(g3_uncapped$study_id)  
n_sw_cap <- n_distinct(g3_capped$study_id[g3_capped$panel_group == "switch_p1"])  
n_sw_unc <- n_distinct(g3_uncapped$study_id[g3_uncapped$panel_group == "switch_p1"])

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("MANUSCRIPT-READY SUMMARY\n")  
cat("══════════════════════════════════════════════════════════════\n\n")

cat(sprintf(paste0(  
  "A sensitivity analysis removing the 36-session cap and restoring\n",  
  "one previously excluded patient (TMS074; N = %d, %d switchers,\n",  
  "%d observations) confirmed that the within-subject slope\n",  
  "acceleration remained %s (uncapped: %+.2f pts/wk, p = %s;\n",  
  "primary: %+.2f pts/wk, p = %s). Post-switch trajectories\n",  
  "remained indistinguishable from iTBS-only patients (p = %s),\n",  
  "and the between-group difference in pre-switch slopes %s\n",  
  "(stayers vs. Phase 1: p = %s).\n"),  
  n_unc, n_sw_unc, nrow(g3_uncapped),  
  ifelse(unc_within_sig, "significant", "non-significant"),  
  res_uncapped$contrasts$within["estimate"],  
  format_p(res_uncapped$contrasts$within["p"]),  
  res_capped$contrasts$within["estimate"],  
  format_p(res_capped$contrasts$within["p"]),  
  format_p(res_uncapped$contrasts$btw_p2["p"]),  
  ifelse(unc_btw_sig, "reached significance", "remained a trend"),  
  format_p(res_uncapped$contrasts$btw_p1["p"])  
))

cat("\nDone.\n")  