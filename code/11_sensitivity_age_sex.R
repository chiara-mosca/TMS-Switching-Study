# ============================================================================  
# SENSITIVITY ANALYSIS: 3-Panel PHQ-9 Slopes  
# Controlling for Age and Sex  
# TMS074 excluded | Capped at 36 sessions | Matched switchers  
# ============================================================================

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

cat(sprintf("Patients with demographics: %d\n", nrow(demo)))  
cat(sprintf("Age range: %d - %d\n", min(demo$age), max(demo$age)))  
cat(sprintf("Sex: Male n = %d, Female n = %d\n\n",  
            sum(demo$sex_numeric == 1), sum(demo$sex_numeric == 2)))

# ── 3. Prepare longitudinal data ──  
c1 <- raw %>%  
  filter(study_id != "TMS074") %>%  
  filter(cumulative_tx <= 36) %>%  
  filter(!is.na(phq9_score), !is.na(phq9_date)) %>%  
  mutate(  
    phq9_score = as.numeric(phq9_score),  
    phq9_date = as.POSIXct(phq9_date)  
  ) %>%  
  arrange(study_id, phq9_date)

stayer_ids <- c1 %>% filter(group == "iTBS_only") %>% pull(study_id) %>% unique()  
switcher_ids <- c1 %>% filter(group == "switcher") %>% pull(study_id) %>% unique()

# Only keep matched switchers (both phases)  
switch_pre <- c1 %>% filter(group == "switcher" & phase == "pre") %>% pull(study_id) %>% unique()  
switch_post <- c1 %>% filter(group == "switcher" & phase == "post") %>% pull(study_id) %>% unique()  
switcher_ids_both <- intersect(switch_pre, switch_post)

dropped <- setdiff(switch_pre, switch_post)  
if (length(dropped) > 0) {  
  cat(sprintf("Dropped %d switcher(s): %s\n", length(dropped), paste(dropped, collapse = ", ")))  
}

cat(sprintf("TMS074 excluded: %s\n\n", !("TMS074" %in% c1$study_id)))

# ── 4. Build 3-panel dataset ──  
dat_stay <- c1 %>%  
  filter(study_id %in% stayer_ids) %>%  
  mutate(panel_group = "stay_p1")

dat_switch_p1 <- c1 %>%  
  filter(study_id %in% switcher_ids_both, phase == "pre") %>%  
  mutate(panel_group = "switch_p1")

dat_switch_p2 <- c1 %>%  
  filter(study_id %in% switcher_ids_both, phase == "post") %>%  
  mutate(panel_group = "switch_p2")

g3 <- bind_rows(dat_stay, dat_switch_p1, dat_switch_p2)

# ── 5. Merge demographics ──  
g3 <- g3 %>%  
  left_join(demo %>% select(study_id, age, sex_numeric), by = "study_id") %>%  
  filter(!is.na(age), !is.na(sex_numeric))

cat(sprintf("Final dataset: %d observations\n", nrow(g3)))  
cat(sprintf("  Stayers:      %d patients\n", n_distinct(g3$study_id[g3$panel_group == "stay_p1"])))  
cat(sprintf("  Switchers P1: %d patients\n", n_distinct(g3$study_id[g3$panel_group == "switch_p1"])))  
cat(sprintf("  Switchers P2: %d patients\n\n", n_distinct(g3$study_id[g3$panel_group == "switch_p2"])))

# ── 6. Compute weeks_in_phase, baseline PHQ-9, center age ──  
g3 <- g3 %>%  
  group_by(study_id, panel_group) %>%  
  mutate(  
    phase_start = min(phq9_date),  
    weeks_in_phase = as.numeric(difftime(phq9_date, phase_start, units = "weeks"))  
  ) %>%  
  ungroup() %>%  
  mutate(age_centered = age - mean(age, na.rm = TRUE))

# Baseline PHQ-9 (visit 0)  
baseline_phq9 <- c1 %>%  
  filter(visit == 0) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup() %>%  
  select(study_id, baseline_phq9 = phq9_score)

g3 <- g3 %>%  
  left_join(baseline_phq9, by = "study_id")

g3$panel_group <- factor(g3$panel_group, levels = c("switch_p1", "stay_p1", "switch_p2"))

# ══════════════════════════════════════════════════════════════  
# 7. FIT BOTH MODELS  
# ══════════════════════════════════════════════════════════════

cat("Fitting PRIMARY model (baseline PHQ-9 covariate)...\n")  
m_primary <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + baseline_phq9 +  
    (weeks_in_phase | study_id) +  
    (1 | study_id:panel_group),  
  data = g3, REML = FALSE,  
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))  
)

cat("Fitting ADJUSTED model (baseline PHQ-9 + age + sex)...\n\n")  
m_adj <- lmer(  
  phq9_score ~ weeks_in_phase * panel_group + baseline_phq9 + age_centered + sex_numeric +  
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
    contrasts = list(within = within_c, btw_p1 = btw_p1, btw_p2 = btw_p2)  
  )  
}

res_primary <- extract_results(m_primary, "PRIMARY MODEL (baseline PHQ-9)")  
res_adj     <- extract_results(m_adj, "ADJUSTED MODEL (baseline PHQ-9 + age + sex)")

# Covariate effects  
cat("\n── Covariate Effects (Adjusted Model) ──\n")  
adj_coefs <- summary(m_adj)$coefficients  
cat(sprintf("  Baseline PHQ-9: B = %+.3f, p %s\n",  
            adj_coefs["baseline_phq9", "Estimate"],  
            ifelse(adj_coefs["baseline_phq9", "Pr(>|t|)"] < 0.001, "< .001",  
                   sprintf("= %.3f", adj_coefs["baseline_phq9", "Pr(>|t|)"]))))  
cat(sprintf("  Age (centered): B = %+.3f, p = %.3f\n",  
            adj_coefs["age_centered", "Estimate"], adj_coefs["age_centered", "Pr(>|t|)"]))  
cat(sprintf("  Sex:            B = %+.3f, p = %.3f\n",  
            adj_coefs["sex_numeric", "Estimate"], adj_coefs["sex_numeric", "Pr(>|t|)"]))

# ══════════════════════════════════════════════════════════════  
# 9. COMPARISON TABLE  
# ══════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("          SENSITIVITY ANALYSIS: SIDE-BY-SIDE\n")  
cat("══════════════════════════════════════════════════════════════\n\n")  
cat(sprintf("%-25s %15s %15s\n", "", "Primary", "Adjusted"))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Stayer slope",  
            res_primary$slopes["stay"], res_adj$slopes["stay"]))  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Switcher P1 slope",  
            res_primary$slopes["switch_p1"], res_adj$slopes["switch_p1"]))  
cat(sprintf("%-25s %+15.3f %+15.3f\n", "Switcher P2 slope",  
            res_primary$slopes["switch_p2"], res_adj$slopes["switch_p2"]))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Within (P2-P1)",  
            res_primary$contrasts$within["p"], res_adj$contrasts$within["p"]))  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Between (Stay-P1)",  
            res_primary$contrasts$btw_p1["p"], res_adj$contrasts$btw_p1["p"]))  
cat(sprintf("%-25s p=%.4f        p=%.4f\n", "Between (Stay-P2)",  
            res_primary$contrasts$btw_p2["p"], res_adj$contrasts$btw_p2["p"]))  
cat(paste(rep("-", 57), collapse = ""), "\n")  
cat(sprintf("AIC:                    %10.1f     %10.1f\n", AIC(m_primary), AIC(m_adj)))

# ══════════════════════════════════════════════════════════════  
# 10. FULL MODEL SUMMARIES  
# ══════════════════════════════════════════════════════════════

cat("\n══ PRIMARY MODEL ══\n")  
print(summary(m_primary))  
cat("\n══ ADJUSTED MODEL ══\n")  
print(summary(m_adj))  
cat("\n══ MODEL COMPARISON ══\n")  
print(anova(m_primary, m_adj))

# ══════════════════════════════════════════════════════════════  
# 11. ROBUSTNESS ASSESSMENT  
# ══════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")  
cat("ROBUSTNESS ASSESSMENT\n")  
cat("══════════════════════════════════════════════════════════════\n\n")

# Check if conclusions change  
primary_sig <- res_primary$contrasts$within["p"] < 0.05  
adj_sig <- res_adj$contrasts$within["p"] < 0.05

if (primary_sig == adj_sig) {  
  cat("CONCLUSION: Results are ROBUST to adjustment for age and sex.\n")  
  cat("The within-subject slope change remains ")  
  cat(ifelse(adj_sig, "significant", "non-significant"))  
  cat(" after adjustment.\n")  
} else {  
  cat("WARNING: Conclusions CHANGE after adjustment for age and sex.\n")  
}

cat(sprintf("\nSlope changes:\n"))  
cat(sprintf("  Stayer:      %+.3f -> %+.3f (diff = %.3f)\n",  
            res_primary$slopes["stay"], res_adj$slopes["stay"],  
            abs(res_primary$slopes["stay"] - res_adj$slopes["stay"])))  
cat(sprintf("  Switcher P1: %+.3f -> %+.3f (diff = %.3f)\n",  
            res_primary$slopes["switch_p1"], res_adj$slopes["switch_p1"],  
            abs(res_primary$slopes["switch_p1"] - res_adj$slopes["switch_p1"])))  
cat(sprintf("  Switcher P2: %+.3f -> %+.3f (diff = %.3f)\n",  
            res_primary$slopes["switch_p2"], res_adj$slopes["switch_p2"],  
            abs(res_primary$slopes["switch_p2"] - res_adj$slopes["switch_p2"])))

cat("\nDone.\n")  