# Data Provenance Log

## Study: To Switch or Not to Switch?  
## Last Updated: [TODAY'S DATE]

---

## 1. Data Source

- **System:** REDCap (Zucker Hillside Hospital)  
- **Project:** TMS Clinic Database  
- **Export Date:** [DATE]  
- **Exported By:** [NAME]  
- **Original File:** `redcap_export_YYYY-MM-DD.csv`  
- **Original Dimensions:** [X] rows x [Y] columns  
- **IRB Protocol:** [NUMBER]

---

## 2. Cleaning Pipeline

### Step 1: Initial Import and Variable Selection  
- **Script:** `code/00_data_cleaning.R`  
- **Input:** `data/raw/redcap_export_YYYY-MM-DD.csv`  
- **Output:** `data/intermediate/cleaned_phq9.csv`  
- **Actions:**  
  - Selected relevant variables: [list them]  
  - Renamed variables for consistency  
  - Converted date fields to Date format  
  - Removed duplicate entries (n = ?)  
  - Filtered to MDD diagnosis only

### Step 2: PHQ-9 Data Restructuring  
- **Script:** `code/00_data_cleaning.R`  
- **Input:** `data/intermediate/cleaned_phq9.csv`  
- **Output:** `data/intermediate/item_level_long.csv`  
- **Actions:**  
  - Reshaped from wide to long format (one row per session per patient)  
  - Created session number variable  
  - Created week variable  
  - Calculated PHQ-9 total scores  
  - Extracted individual PHQ-9 items (1-9)

### Step 3: Group Assignment  
- **Script:** `code/00_data_cleaning.R`  
- **Input:** `data/intermediate/item_level_long.csv`  
- **Output:** `data/intermediate/itbs_only.csv`, `data/intermediate/switchers.csv`  
- **Actions:**  
  - Classified patients as iTBS-only (n = 41) or switchers (n = 46)  
  - **Classification criteria:**  
    - [DOCUMENT EXACTLY HOW YOU DEFINED "SWITCHER"]  
    - [How was "non-response" operationalized?]  
    - [At what session/week did switches typically occur?]  
    - [Was there a minimum number of iTBS sessions before switching?]  
  - Created switch_session variable for switchers  
  - Created pre_switch and post_switch phase indicators

### Step 4: Exclusions  
- **Script:** `code/00_data_cleaning.R`  
- **Actions:**  
  - Total charts reviewed: [N]  
  - Excluded: no MDD diagnosis (n = ?)  
  - Excluded: fewer than [X] PHQ-9 assessments (n = ?)  
  - Excluded: [other reasons] (n = ?)  
  - **Final analytic sample: N = 87**

### Step 5: Final Analysis Dataset  
- **Script:** `code/00_data_cleaning.R`  
- **Input:** intermediate files  
- **Output:** `data/final/analysis_ready.csv`  
- **Dimensions:** [X] rows x [Y] columns  
- **Verification:**  
  - Confirmed n = 41 iTBS-only, n = 46 switchers  
  - Confirmed no missing group assignments  
  - Confirmed PHQ-9 total = sum of items 1-9  
  - Spot-checked 5 random patients against original REDCap records

---

## 3. Variable Dictionary

| Variable | Type | Description | Values/Range |  
|---|---|---|---|  
| patient_id | character | De-identified patient ID | P001-P087 |  
| group | factor | Treatment group | "itbs_only", "switcher" |  
| session_num | integer | Session number | 1-[max] |  
| week | numeric | Week of treatment | 0-[max] |  
| phq9_total | integer | PHQ-9 total score | 0-27 |  
| phq9_item1 | integer | Anhedonia | 0-3 |  
| phq9_item2 | integer | Depressed mood | 0-3 |  
| phq9_item3 | integer | Sleep disturbance | 0-3 |  
| phq9_item4 | integer | Fatigue | 0-3 |  
| phq9_item5 | integer | Appetite changes | 0-3 |  
| phq9_item6 | integer | Guilt/worthlessness | 0-3 |  
| phq9_item7 | integer | Concentration | 0-3 |  
| phq9_item8 | integer | Psychomotor changes | 0-3 |  
| phq9_item9 | integer | Suicidal ideation | 0-3 |  
| switch_session | integer | Session at which protocol changed | NA for iTBS-only |  
| phase | factor | Treatment phase | "pre_switch", "post_switch", NA |  
| age | numeric | Age at treatment start | [range] |  
| sex | factor | Sex | "Male", "Female" |

---

## 4. Analysis Scripts

| Script | Purpose | Output |  
|---|---|---|  
| 00_data_cleaning.R | REDCap export to analysis-ready dataset | cleaned CSVs |  
| 01_sample_characteristics.R | Table 1, demographics, baseline comparisons | Table 1 |  
| 02_trajectory_analysis.R | HLM models, Figures 1-3 | Figures 1-3 |  
| 03_item_analysis.R | PHQ-9 item-level models, Figures 4-5 | Figures 4-5 |  
| 04_sensitivity_analyses.R | Covariate controls, robustness checks | Supplementary |  
