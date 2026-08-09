# ============================================================================  
# Table 1: Baseline Demographics  
# Journal of Affective Disorders  
# ============================================================================

required_packages <- c("readxl", "dplyr", "gtsummary", "gt", "flextable", "officer")

for (pkg in required_packages) {  
  if (!require(pkg, character.only = TRUE)) {  
    install.packages(pkg)  
    library(pkg, character.only = TRUE)  
  }  
}

# --- STEP 1: Paste the path from file.choose() here ---  
file_path <- "/Users/chiara/Documents/TMS-Switching-Study/analytic_sample_summary.xlsx"

# Example: it might look like one of these:  
# file_path <- "/Users/chiaramosca/Documents/TMS-Switching-Study/analytic_sample_summary.xlsx"  
# file_path <- "/Users/chiaramosca/Documents/TMS-Switching-Study/analytic_sample_summary.csv"  
# file_path <- "/Users/chiaramosca/Documents/TMS Switching Study/analytic_sample_summary.xlsx"

# --- STEP 2: Read file (auto-detects csv vs xlsx) ---  
if (grepl("\\.csv$", file_path)) {  
  df <- read.csv(file_path)  
} else if (grepl("\\.xlsx?$", file_path)) {  
  df <- readxl::read_excel(file_path)  
} else {  
  stop("Unrecognized file type. Must be .csv or .xlsx")  
}

cat("✓ File loaded:", nrow(df), "rows x", ncol(df), "columns\n")  
cat("Columns:", paste(names(df), collapse = ", "), "\n\n")

# ============================================================================  
# DATA PREPARATION  
# ============================================================================

df <- df %>%  
  mutate(  
    sex_label = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),  
    
    group_label = factor(  
      group,  
      levels = c("iTBS_only", "switcher"),  
      labels = c("iTBS-Only", "Switchers")  
    ),  
    
    age = age_at_first_phq,  
    baseline_phq9 = tx_baseline,  
    total_sessions = total_tx,  
    treatment_weeks = total_weeks  
  )

# ============================================================================  
# TABLE 1  
# ============================================================================

table_data <- df %>%  
  select(group_label, age, sex_label, baseline_phq9, total_sessions, treatment_weeks)

table1 <- table_data %>%  
  tbl_summary(  
    by = group_label,  
    
    label = list(  
      age             ~ "Age, M (SD)",  
      sex_label       ~ "Sex, n (%)",  
      baseline_phq9   ~ "Baseline PHQ-9, M (SD)",  
      total_sessions  ~ "Number of sessions, M (SD)",  
      treatment_weeks ~ "Treatment duration (weeks), M (SD)"  
    ),  
    
    type = list(  
      age             ~ "continuous",  
      baseline_phq9   ~ "continuous",  
      total_sessions  ~ "continuous",  
      treatment_weeks ~ "continuous",  
      sex_label       ~ "categorical"  
    ),  
    
    statistic = list(  
      all_continuous()  ~ "{mean} ({sd})",  
      all_categorical() ~ "{n} ({p}%)"  
    ),  
    
    digits = list(  
      all_continuous()  ~ 2,  
      all_categorical() ~ c(0, 1)  
    ),  
    
    missing = "no"  
    
  ) %>%  
  
  add_p(  
    test = list(  
      all_continuous()  ~ "t.test",  
      all_categorical() ~ "chisq.test"  
    ),  
    pvalue_fun = ~ style_pvalue(.x, digits = 2)  
  ) %>%  
  
  modify_header(  
    label  ~ "**Variable**",  
    stat_1 ~ "**iTBS-Only (n = {n})**",  
    stat_2 ~ "**Switchers (n = {n})**"  
  ) %>%  
  
  modify_spanning_header(  
    c("stat_1", "stat_2") ~ "**Treatment Group**"  
  ) %>%  
  
  modify_caption("**Table 1. Baseline Demographics and Clinical Characteristics**") %>%  
  
  modify_footnote(  
    all_stat_cols() ~ "Mean (SD) for continuous variables; n (%) for categorical variables."  
  ) %>%  
  
  bold_labels()

print(table1)

# ============================================================================  
# EXPORT — save to same folder as data file  
# ============================================================================

output_dir <- dirname(file_path)

# Word  
table1 %>%  
  as_flex_table() %>%  
  flextable::save_as_docx(  
    path = file.path(output_dir, "Table1_Demographics.docx"),  
    pr_section = officer::prop_section(  
      page_size = officer::page_size(orient = "portrait"),  
      page_margins = officer::page_mar(top = 1, bottom = 1, left = 0.75, right = 0.75)  
    )  
  )

cat("\n✓ Saved:", file.path(output_dir, "Table1_Demographics.docx"), "\n")

# HTML  
table1 %>%  
  as_gt() %>%  
  gt::gtsave(filename = file.path(output_dir, "Table1_Demographics.html"))

cat("✓ Saved:", file.path(output_dir, "Table1_Demographics.html"), "\n")

cat("\n✓ Done!\n")  