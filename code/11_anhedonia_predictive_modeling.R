# ══════════════════════════════════════════════════════════════  
# Logistic Regression & ROC Analysis:  
# Can early anhedonia non-response predict switch need?  
# ══════════════════════════════════════════════════════════════

library(tidyverse)  
library(readxl)  
library(pROC)  
library(rms)       # for Nagelkerke R²  
library(caret)     # for confusion matrix

# --- 1. Load Data ---  
df <- read_excel("/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_long.xlsx", sheet = "Sheet1")  
colnames(df) <- trimws(colnames(df))

# --- 2. Prepare data ---  
df <- df %>%  
  mutate(  
    cumulative_tx = as.numeric(cumulative_tx),  
    phq9_1 = as.numeric(phq9_1),  
    phq9_score = as.numeric(phq9_score),  
    visit = as.numeric(visit),  
    group = as.character(group)  
  ) %>%  
  filter(group %in% c("iTBS_only", "switcher"))

# --- 3. Create binary outcome: 1 = switcher, 0 = iTBS_only ---  
df <- df %>%  
  mutate(is_switcher = ifelse(group == "switcher", 1, 0))

# ══════════════════════════════════════════════════════════════  
# APPROACH 1: Anhedonia at visit 1 (~session 10, ~week 1)  
# This captures "early" response after first week of treatment  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("APPROACH 1: Anhedonia at Visit 1 (~ week 1, ~ session 10)\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

df_visit1 <- df %>%  
  filter(visit == 1) %>%  
  filter(!is.na(phq9_1)) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup()

cat("N patients at visit 1:", nrow(df_visit1), "\n")  
cat("  Switchers:", sum(df_visit1$is_switcher), "\n")  
cat("  iTBS Only:", sum(df_visit1$is_switcher == 0), "\n\n")

# Logistic regression  
model_v1 <- glm(is_switcher ~ phq9_1, data = df_visit1, family = binomial)  
summary(model_v1)

# OR and 95% CI  
or_v1 <- exp(coef(model_v1))  
ci_v1 <- exp(confint(model_v1))  
cat("\n--- Odds Ratios ---\n")  
cat("Anhedonia OR:", round(or_v1["phq9_1"], 3), "\n")  
cat("95% CI:", round(ci_v1["phq9_1", 1], 3), "-", round(ci_v1["phq9_1", 2], 3), "\n")  
cat("p-value:", round(summary(model_v1)$coefficients["phq9_1", "Pr(>|z|)"], 6), "\n")

# Nagelkerke R²  
n <- nrow(df_visit1)  
ll_null <- logLik(glm(is_switcher ~ 1, data = df_visit1, family = binomial))  
ll_full <- logLik(model_v1)  
cox_snell <- 1 - exp((2/n) * (as.numeric(ll_null) - as.numeric(ll_full)))  
nagelkerke <- cox_snell / (1 - exp((2/n) * as.numeric(ll_null)))  
cat("Nagelkerke R²:", round(nagelkerke, 4), "\n")

# Classification accuracy (threshold = 0.5)  
pred_prob_v1 <- predict(model_v1, type = "response")  
pred_class_v1 <- ifelse(pred_prob_v1 > 0.5, 1, 0)  
cm_v1 <- confusionMatrix(factor(pred_class_v1), factor(df_visit1$is_switcher), positive = "1")  
cat("Classification Accuracy:", round(cm_v1$overall["Accuracy"], 4), "\n")  
cat("Sensitivity:", round(cm_v1$byClass["Sensitivity"], 4), "\n")  
cat("Specificity:", round(cm_v1$byClass["Specificity"], 4), "\n")

# ROC analysis  
roc_v1 <- roc(df_visit1$is_switcher, pred_prob_v1, quiet = TRUE)  
cat("AUC:", round(auc(roc_v1), 4), "\n")  
cat("AUC 95% CI:", round(ci.auc(roc_v1)[1], 4), "-", round(ci.auc(roc_v1)[3], 4), "\n")

# ══════════════════════════════════════════════════════════════  
# APPROACH 2: Anhedonia CHANGE from baseline to visit 1  
# Non-response = failure to decrease anhedonia  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("APPROACH 2: Anhedonia CHANGE (baseline to visit 1)\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Get baseline (visit 0) anhedonia  
df_baseline <- df %>%  
  filter(visit == 0) %>%  
  filter(!is.na(phq9_1)) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup() %>%  
  select(study_id, phq9_1_baseline = phq9_1, is_switcher)

# Get visit 1 anhedonia  
df_v1_scores <- df %>%  
  filter(visit == 1) %>%  
  filter(!is.na(phq9_1)) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup() %>%  
  select(study_id, phq9_1_v1 = phq9_1)

# Merge and compute change  
df_change <- df_baseline %>%  
  inner_join(df_v1_scores, by = "study_id") %>%  
  mutate(  
    phq9_1_change = phq9_1_v1 - phq9_1_baseline,  # positive = worsening  
    phq9_1_pct_change = ifelse(phq9_1_baseline > 0,  
                               (phq9_1_change / phq9_1_baseline) * 100, NA),  
    anhedonia_nonresponse = ifelse(phq9_1_change >= 0, 1, 0)  # no improvement = non-response  
  )

cat("N patients with baseline + visit 1:", nrow(df_change), "\n")  
cat("  Anhedonia non-responders:", sum(df_change$anhedonia_nonresponse, na.rm = TRUE), "\n\n")

# Model with anhedonia change  
model_change <- glm(is_switcher ~ phq9_1_change, data = df_change, family = binomial)

or_change <- exp(coef(model_change))  
ci_change <- exp(confint(model_change))  
cat("--- Change Model ---\n")  
cat("OR per 1-point worsening:", round(or_change["phq9_1_change"], 3), "\n")  
cat("95% CI:", round(ci_change["phq9_1_change", 1], 3), "-",  
    round(ci_change["phq9_1_change", 2], 3), "\n")  
cat("p-value:", round(summary(model_change)$coefficients["phq9_1_change", "Pr(>|z|)"], 6), "\n")

# Nagelkerke R² for change model  
ll_null_c <- logLik(glm(is_switcher ~ 1, data = df_change, family = binomial))  
ll_full_c <- logLik(model_change)  
n_c <- nrow(df_change)  
cox_snell_c <- 1 - exp((2/n_c) * (as.numeric(ll_null_c) - as.numeric(ll_full_c)))  
nagelkerke_c <- cox_snell_c / (1 - exp((2/n_c) * as.numeric(ll_null_c)))  
cat("Nagelkerke R²:", round(nagelkerke_c, 4), "\n")

# ROC for change model  
pred_prob_change <- predict(model_change, type = "response")  
roc_change <- roc(df_change$is_switcher, pred_prob_change, quiet = TRUE)  
cat("AUC:", round(auc(roc_change), 4), "\n")  
cat("AUC 95% CI:", round(ci.auc(roc_change)[1], 4), "-", round(ci.auc(roc_change)[3], 4), "\n")

# Model with binary non-response  
model_nonresp <- glm(is_switcher ~ anhedonia_nonresponse, data = df_change, family = binomial)  
or_nr <- exp(coef(model_nonresp))  
ci_nr <- exp(confint(model_nonresp))  
cat("\n--- Binary Non-Response Model ---\n")  
cat("OR (non-response vs response):", round(or_nr["anhedonia_nonresponse"], 3), "\n")  
cat("95% CI:", round(ci_nr["anhedonia_nonresponse", 1], 3), "-",  
    round(ci_nr["anhedonia_nonresponse", 2], 3), "\n")  
cat("p-value:", round(summary(model_nonresp)$coefficients["anhedonia_nonresponse", "Pr(>|z|)"], 6), "\n")

# ══════════════════════════════════════════════════════════════  
# APPROACH 3: Anhedonia at visit 2 (~session 15, ~week 2)  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("APPROACH 3: Anhedonia at Visit 2 (~ week 2, ~ session 15)\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

df_visit2 <- df %>%  
  filter(visit == 2) %>%  
  filter(!is.na(phq9_1)) %>%  
  group_by(study_id) %>%  
  slice(1) %>%  
  ungroup()

cat("N patients at visit 2:", nrow(df_visit2), "\n")

model_v2 <- glm(is_switcher ~ phq9_1, data = df_visit2, family = binomial)  
or_v2 <- exp(coef(model_v2))  
ci_v2 <- exp(confint(model_v2))  
cat("Anhedonia OR:", round(or_v2["phq9_1"], 3), "\n")  
cat("95% CI:", round(ci_v2["phq9_1", 1], 3), "-", round(ci_v2["phq9_1", 2], 3), "\n")  
cat("p-value:", round(summary(model_v2)$coefficients["phq9_1", "Pr(>|z|)"], 6), "\n")

# Nagelkerke R²  
ll_null_v2 <- logLik(glm(is_switcher ~ 1, data = df_visit2, family = binomial))  
ll_full_v2 <- logLik(model_v2)  
n_v2 <- nrow(df_visit2)  
cox_snell_v2 <- 1 - exp((2/n_v2) * (as.numeric(ll_null_v2) - as.numeric(ll_full_v2)))  
nagelkerke_v2 <- cox_snell_v2 / (1 - exp((2/n_v2) * as.numeric(ll_null_v2)))  
cat("Nagelkerke R²:", round(nagelkerke_v2, 4), "\n")

pred_prob_v2 <- predict(model_v2, type = "response")  
pred_class_v2 <- ifelse(pred_prob_v2 > 0.5, 1, 0)  
cm_v2 <- confusionMatrix(factor(pred_class_v2), factor(df_visit2$is_switcher), positive = "1")  
cat("Classification Accuracy:", round(cm_v2$overall["Accuracy"], 4), "\n")

roc_v2 <- roc(df_visit2$is_switcher, pred_prob_v2, quiet = TRUE)  
cat("AUC:", round(auc(roc_v2), 4), "\n")  
cat("AUC 95% CI:", round(ci.auc(roc_v2)[1], 4), "-", round(ci.auc(roc_v2)[3], 4), "\n")

# ══════════════════════════════════════════════════════════════  
# APPROACH 4: Combined model — baseline + change  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("APPROACH 4: Combined Model (baseline anhedonia + change)\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

model_combined <- glm(is_switcher ~ phq9_1_baseline + phq9_1_change,  
                      data = df_change, family = binomial)  
summary(model_combined)

or_comb <- exp(coef(model_combined))  
ci_comb <- exp(confint(model_combined))  
cat("\n--- Combined Model ORs ---\n")  
for (var in c("phq9_1_baseline", "phq9_1_change")) {  
  cat(var, "OR:", round(or_comb[var], 3),  
      "  95% CI:", round(ci_comb[var, 1], 3), "-", round(ci_comb[var, 2], 3),  
      "  p:", round(summary(model_combined)$coefficients[var, "Pr(>|z|)"], 6), "\n")  
}

# Nagelkerke R²  
ll_null_comb <- logLik(glm(is_switcher ~ 1, data = df_change, family = binomial))  
ll_full_comb <- logLik(model_combined)  
cox_snell_comb <- 1 - exp((2/n_c) * (as.numeric(ll_null_comb) - as.numeric(ll_full_comb)))  
nagelkerke_comb <- cox_snell_comb / (1 - exp((2/n_c) * as.numeric(ll_null_comb)))  
cat("Nagelkerke R²:", round(nagelkerke_comb, 4), "\n")

pred_prob_comb <- predict(model_combined, type = "response")  
roc_comb <- roc(df_change$is_switcher, pred_prob_comb, quiet = TRUE)  
cat("AUC:", round(auc(roc_comb), 4), "\n")  
cat("AUC 95% CI:", round(ci.auc(roc_comb)[1], 4), "-", round(ci.auc(roc_comb)[3], 4), "\n")

# ══════════════════════════════════════════════════════════════  
# FIGURE: ROC Curves Comparison  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("Generating ROC Curve Figure...\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Build ROC data for ggplot  
roc_data <- bind_rows(  
  tibble(  
    sensitivity = roc_v1$sensitivities,  
    specificity = 1 - roc_v1$specificities,  
    model = paste0("Visit 1 Anhedonia (AUC = ", round(auc(roc_v1), 3), ")")  
  ),  
  tibble(  
    sensitivity = roc_v2$sensitivities,  
    specificity = 1 - roc_v2$specificities,  
    model = paste0("Visit 2 Anhedonia (AUC = ", round(auc(roc_v2), 3), ")")  
  ),  
  tibble(  
    sensitivity = roc_change$sensitivities,  
    specificity = 1 - roc_change$specificities,  
    model = paste0("Anhedonia Change (AUC = ", round(auc(roc_change), 3), ")")  
  ),  
  tibble(  
    sensitivity = roc_comb$sensitivities,  
    specificity = 1 - roc_comb$specificities,  
    model = paste0("Combined (AUC = ", round(auc(roc_comb), 3), ")")  
  )  
)

p_roc <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = model)) +  
  geom_line(linewidth = 1.2) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +  
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A", "#FF7F00")) +  
  labs(  
    x = "1 - Specificity (False Positive Rate)",  
    y = "Sensitivity (True Positive Rate)",  
    title = "ROC Curves: Early Anhedonia as Predictor of Protocol Switch Need",  
    subtitle = "Can anhedonia non-response at week 1-2 predict which patients will need bilateral TMS?",  
    color = NULL  
  ) +  
  theme_minimal(base_size = 14) +  
  theme(  
    plot.background = element_rect(fill = "white", color = NA),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.3),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.title = element_text(face = "bold", size = 13),  
    axis.text = element_text(size = 11, color = "black"),  
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),  
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",  
                                 margin = margin(b = 15)),  
    legend.position = c(0.65, 0.25),  
    legend.background = element_rect(fill = "white", color = "#CCCCCC", linewidth = 0.3),  
    legend.text = element_text(size = 10),  
    plot.margin = margin(15, 20, 10, 10)  
  ) +  
  coord_equal()

print(p_roc)  
ggsave("/Users/chiara/Documents/TMS-Switching-Study/roc_anhedonia_predictor.png",  
       p_roc, width = 9, height = 8, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════  
# SUMMARY TABLE  
# ══════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 60), collapse = ""), "\n")  
cat("SUMMARY TABLE\n")  
cat(paste(rep("=", 60), collapse = ""), "\n\n")

summary_table <- tibble(  
  Model = c("Visit 1 Anhedonia", "Visit 2 Anhedonia",  
            "Anhedonia Change (V0→V1)", "Combined (Baseline + Change)"),  
  OR = c(round(or_v1["phq9_1"], 3),  
         round(or_v2["phq9_1"], 3),  
         round(or_change["phq9_1_change"], 3),  
         paste0(round(or_comb["phq9_1_baseline"], 3), " / ",  
                round(or_comb["phq9_1_change"], 3))),  
  CI_95 = c(paste0(round(ci_v1["phq9_1", 1], 2), "-", round(ci_v1["phq9_1", 2], 2)),  
            paste0(round(ci_v2["phq9_1", 1], 2), "-", round(ci_v2["phq9_1", 2], 2)),  
            paste0(round(ci_change["phq9_1_change", 1], 2), "-",  
                   round(ci_change["phq9_1_change", 2], 2)),  
            "see above"),  
  p_value = c(round(summary(model_v1)$coefficients["phq9_1", "Pr(>|z|)"], 5),  
              round(summary(model_v2)$coefficients["phq9_1", "Pr(>|z|)"], 5),  
              round(summary(model_change)$coefficients["phq9_1_change", "Pr(>|z|)"], 5),  
              "see above"),  
  Nagelkerke_R2 = c(round(nagelkerke, 4), round(nagelkerke_v2, 4),  
                    round(nagelkerke_c, 4), round(nagelkerke_comb, 4)),  
  AUC = c(round(auc(roc_v1), 4), round(auc(roc_v2), 4),  
          round(auc(roc_change), 4), round(auc(roc_comb), 4))  
)

print(summary_table)

# Save summary  
write.csv(summary_table,  
          "/Users/chiara/Documents/TMS-Switching-Study/anhedonia_logistic_regression_summary.csv",  
          row.names = FALSE)

cat("\nResults saved to /Users/chiara/Documents/TMS-Switching-Study/\n")  
cat("  - roc_anhedonia_predictor.png\n")  
cat("  - anhedonia_logistic_regression_summary.csv\n")  