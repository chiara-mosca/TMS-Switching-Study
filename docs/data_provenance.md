# Data Provenance Log

## Study: To Switch or Not to Switch?: Protocol Change from iTBS to Bilateral Transcranial Magnetic Stimulation is Associated with Improved PHQ-9 Trajectory  
## Last Updated: 8/9/26 by Chiara Mosca

---

## 1. Data Source

- **System:** REDCap (Zucker Hillside Hospital)  
- **Project:** TMS Clinic Database  
- **Export Date:** 6/15/26 
- **Exported By:** Noelle Arcaro  
- **Original File:** PHQ9Project_6.15.26.csv  
- **Original Dimensions:** 1,211 rows x 132 columns  
- **IRB Protocol:** Submitted, pending approval

---

## 2. Cleaning Pipeline

### Step 1: Patient De-identification (Raw Data File -> Intermediate Data File)
- **Script:**  
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv  
- **Output:** Intermediate ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
- **Actions:**  
  - Extracted 146 unique patients from record_id
  - Assigned sequential de-identified IDs (TMS001-TMS146)
  - Mapped patient names to study_id
 
### Step 2: Initial Exclusion Criteria (Raw Data File -> Intermediate Data File)
- **Script:**
- ***Total patients in REDCap raw data file:** 146
- ***Total patients in intermediate cleaned data:** 130
- **Excluded:** 16 patients
    - 1 test/template row (TMS001)
    -  15 patients with no usable PHQ-9 data: TMS005, TMS019, TMS026, TMS028, TMS030, TMS042, TMS059, TMS068, TMS072, TMS078, TMS080, TMS089, TMS104, TMS105, TMS109
    - Reasons: no PHQ-9 forms in folder, rTMS study participants with blank course 1, discontinued courses with no assessments

### Step 3: PHQ-9 Data Extraction (Raw Data File -> Intermediate Data File)  
- **Script:**  
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv & Intermediate ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
- **Output:** phq_scores_clean.xlsx
- **Actions:**  
  - Replaced patient names with de-identified study_ids
  - Extracted course number (1-4) from redcap_event_name
  - Dropped course summary columns (demographics, protocol info, session counts)
  - Retained: study_id, course, phq9_date, protocol_doa, treatmentnumber, phq9_1, phq9_2, phq9_3, phq9_4, phq9_5, phq9_6, phq9_7, phq9_8, phq9_9, phq9_score, phq9_how difficult
  - protocol_doa coding: 0 = baseline/pre-treatment, 1=iTBS, 2=iTBS/cTBS Bilateral, 3=cTBS, 4=10 Hz, 5=10 Hz/1 Hz Bilateral
 
**Variable Dictionary - phq_scores_clean.xlsx (intermediate data file)**
| Variable | Type | Description | Values/Range |  
|---|---|---|---|  
| study_id | character | De-identified patient ID | TMS002–TMS146 |  
| course | integer | Treatment course number | 1–4 |  
| phq9_date | date | Date of PHQ-9 assessment | Various |  
| protocol_doa | integer | Protocol at date of assessment | 0=baseline/pretreatment, 1=iTBS, 2=iTBS/cTBS Bilateral, 3=cTBS, 4=10Hz, 5=10Hz/1Hz Bilateral |  
| treatmentnumber | integer | Session number within current protocol | 1–69 |  
| phq9_1 | integer | Anhedonia (little interest/pleasure) | 0–3 |  
| phq9_2 | integer | Depressed mood (feeling down) | 0–3 |  
| phq9_3 | integer | Sleep disturbance | 0–3 |  
| phq9_4 | integer | Fatigue/low energy | 0–3 |  
| phq9_5 | integer | Appetite changes | 0–3 |  
| phq9_6 | integer | Guilt/worthlessness | 0–3 |  
| phq9_7 | integer | Concentration difficulties | 0–3 |  
| phq9_8 | integer | Psychomotor changes | 0–3 |  
| phq9_9 | integer | Suicidal ideation | 0–3 |  
| phq9_score | integer | PHQ-9 total score | 0–27 |  
| phq9_how_difficult | integer | Difficulty question | 1–4 |

### Step 4: Demographics Extraction (Raw Data File -> Intermediate Data File)  
- **Script:**   
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv & Intermediate ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
- **Output:** patients_clean.xlsx (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) 
- **Actions:**  
  - Extracted sex, tmsdiagnosis, and dob from patientcourse_info rows
  - Calculated age_at_first_phq from dob and first PHQ-9 date
  - Replaced record_id with study_id
  - Retained: study_id, age_at_first_phq, sex, tmsdiagnosis
 
**Variable Dictionary - patients_clean.xlsx (intermediate data file)**
| Variable | Type | Description | Values |  
|---|---|---|---|  
| study_id | character | De-identified patient ID | TMS002–TMS146 |  
| age_at_first_phq | integer | Age at first PHQ-9 | 18–87 |  
| sex | integer | Sex | 1=Male, 2=Female |  
| tmsdiagnosis | integer | TMS diagnosis | 1=MDD, 2=Bipolar, 3=MDD/OCD |



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
