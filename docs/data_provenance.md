# Data Provenance Log

## To Switch or Not to Switch?: Protocol Change from iTBS to Bilateral Transcranial Magnetic Stimulation is Associated with Improved PHQ-9 Trajectory  
## Last Updated: 8/9/26 by CM

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

### Step 1: Patient De-identification (Raw Data File -> Initial Cleaning File)
- **Script:**  
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv  
- **Output:** ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
- **Actions:**  
  - Extracted 146 unique patients from record_id
  - Assigned sequential de-identified IDs (TMS001-TMS146)
  - Mapped patient names to study_id
 
### Step 2: Initial Exclusion Criteria (Raw Data File -> Initial Cleaning File)
- **Script:**
- ***Total patients in REDCap raw data file:** 146
- ***Total patients in initially cleaned data:** 130
- **Excluded:** 16 patients
    - 1 test/template row (TMS001)
    -  15 patients with no usable PHQ-9 data: TMS005, TMS019, TMS026, TMS028, TMS030, TMS042, TMS059, TMS068, TMS072, TMS078, TMS080, TMS089, TMS104, TMS105, TMS109
    - Reasons: no PHQ-9 forms in folder, rTMS study participants with blank course 1, discontinued courses with no assessments

### Step 3: PHQ-9 Data Extraction (Raw Data File -> Initial Cleaning File)  
- **Script:**  
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv & ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
- **Output:** phq_scores_clean.xlsx
- **Actions:**  
  - Replaced patient names with de-identified study_ids
  - Extracted course number (1-4) from redcap_event_name
  - Dropped course summary columns (demographics, protocol info, session counts)
  - Retained: study_id, course, phq9_date, protocol_doa, treatmentnumber, phq9_1, phq9_2, phq9_3, phq9_4, phq9_5, phq9_6, phq9_7, phq9_8, phq9_9, phq9_score, phq9_how difficult
  - protocol_doa coding: 0 = baseline/pre-treatment, 1=iTBS, 2=iTBS/cTBS Bilateral, 3=cTBS, 4=10 Hz, 5=10 Hz/1 Hz Bilateral
 
**Variable Dictionary - phq_scores_clean.xlsx (initial cleaning file)**
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

### Step 4: Demographics Extraction (Raw Data File -> Initial Cleaning File)  
- **Script:**   
- **Input:** Raw REDCap export (1,211 rows x 132 columns) PHQ9Project_6.15.26.csv & ID file (not included in GitHub repository due to patient data protection, contact Noelle Arcaro for file access) ID_refs_PHI.xlsx
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

### Step 5: Initial Cleaning -> Intermediate Data Files 
- **Script:**
- -**Input:** phq_scores_clean
- **Output:** tms_clean_longitudinal.xlsx & tms_clean_outcomes(in).csv 
- **Actions**
  - Filtered to first Tx course only (removed pts who had multiple courses of treatment)
  - Removed pts with missing PHQ-9 data
  - Removed pts with invalid PHQ-9 data (ex: one pt had a PHQ-9 date in year 2056?)
  - Created visit numbering
  - Calculated time variables from first PHQ-9
  - Calculated cumulative treatment
  - Assigned phase - pre (before switch) and post (after switch)
  - Assigned clinical outcomes of remitter vs. responder 

### Step 6: Intermediate Data Files -> Final Files
- **Script:** 01_analytic_sample.R (in /code)
- **Input** tms_clean_longitudinal & tms_clean_outcomes(in).csv
- **Output** analytic_sample_long.xlsx, analytic_sample_summary.xlsx (not added to repository because contains pt demographics information - contact Chiara Mosca for access), & excluded_patients.xlsx
- **Actions**
  - Strictly only included patients who were only ever on iTBS (protocol 1) for stayers and patients who only went from iTBS -> Bilateral (protocol 1 -> 2) for switchers
  - Excluded 10 Hz pts, pts who went from protocol 2 -> 1, pts who jumped between protocols 1 and 2
  - Merged demographic information from patients_clean.xlsx into analytic_sample_summary.xlsx
 
---

## 3. Final Analytic Files

### File 1: analytic_sample_long.xlsx  
- **Location:** data/final/ (in repository)  
- **Dimensions:** 567 rows × 29 columns  
- **Description:** Longitudinal PHQ-9 assessments for N=85 analytic sample patients, with one row per assessment visit

**Variable Dictionary — analytic_sample_long.xlsx**

```  
Variable              | Type      | Description                              | Values/Range  
----------------------|-----------|------------------------------------------|--------------------------------------------  
study_id              | character | De-identified patient ID                 | TMS002–TMS145  
category              | character | Original clinical category               | "Left-Sided Only", "Switched"  
outcome               | character | Clinical outcome classification          | "Remitter", "Responder", "Non-Responder"  
phq9_date             | date      | Date of PHQ-9 assessment                 | 2019-10-10 to 2026-03-26  
phq9_score            | integer   | PHQ-9 total score                        | 0–27  
protocol_doa          | integer   | Protocol at date of assessment           | 1=iTBS (left), 2=iTBS/cTBS bilateral  
protocol_name         | character | Protocol name                            | "iTBS (left)", "iTBS/cTBS bilateral"  
visit                 | integer   | Visit number (0 = baseline)              | 0–11  
weeks                 | numeric   | Weeks since first PHQ-9                  | 0–12  
days                  | integer   | Days since first PHQ-9                   | 0–84  
treatmentnumber       | integer   | Session number within current protocol   | 1–43  
cumulative_tx         | integer   | Total TMS sessions to date               | 1–53  
tx_since_last_phq     | integer   | Sessions since last PHQ-9                | 1–10  
phase                 | character | Treatment phase relative to switch       | "pre" (before switch or iTBS-only), "post" (after switch)  
switch_visit          | integer   | Visit number when protocol switched      | 1–7 (NA for iTBS-only)  
switch_week           | numeric   | Week when protocol switched              | 1–7.14 (NA for iTBS-only)  
switch_tx             | integer   | Cumulative sessions at switch            | 7–40 (NA for iTBS-only)  
consultation_phq9     | integer   | Pre-treatment consultation PHQ-9         | 5–27 (some missing)  
phq9_1                | integer   | Anhedonia (little interest/pleasure)     | 0–3  
phq9_2                | integer   | Depressed mood (feeling down)            | 0–3  
phq9_3                | integer   | Sleep disturbance                        | 0–3  
phq9_4                | integer   | Fatigue/low energy                       | 0–3  
phq9_5                | integer   | Appetite changes                         | 0–3  
phq9_6                | integer   | Guilt/worthlessness                      | 0–3  
phq9_7                | integer   | Concentration difficulties               | 0–3  
phq9_8                | integer   | Psychomotor changes                      | 0–3  
phq9_9                | integer   | Suicidal ideation                        | 0–3  
phq9_how_difficult    | integer   | Functional impairment question           | 1–4 (some missing)  
group                 | character | Group assignment                         | "iTBS_only", "switcher"  


---

### File 2: analytic_sample_summary.xlsx  
- **Location:** Not in repository (contains demographic data — contact Chiara Mosca for access)  
- **Dimensions:** 85 rows × 16 columns  
- **Description:** One row per patient with baseline/endpoint scores, outcome classification, demographics, and group assignment

**Variable Dictionary — analytic_sample_summary.xlsx**

```  
Variable              | Type      | Description                              | Values/Range  
----------------------|-----------|------------------------------------------|--------------------------------------------  
study_id              | character | De-identified patient ID                 | TMS002–TMS145  
tx_baseline           | integer   | Baseline PHQ-9 score (first assessment)  | 3–27  
endpoint              | integer   | Endpoint PHQ-9 score (last assessment)   | 0–27  
n_visits              | integer   | Total number of PHQ-9 assessments        | 1–12  
total_weeks           | numeric   | Total weeks in treatment                 | 0–12  
total_tx              | integer   | Total TMS sessions received              | 5–53  
category              | character | Original clinical category               | "Left-Sided Only", "Switched"  
consultation_phq9     | integer   | Pre-treatment consultation PHQ-9         | 5–27 (some missing)  
change                | integer   | PHQ-9 change (baseline − endpoint)       | −11 to 22  
pct_change            | numeric   | Percent change from baseline             | −100% to 100%  
remitter              | integer   | Remission achieved (endpoint ≤ 4)        | 0=No, 1=Yes  
responder             | integer   | Response achieved (≥50% reduction)       | 0=No, 1=Yes  
outcome               | character | Clinical outcome                         | "Remitter", "Responder", "Non-Responder"  
group                 | character | Group assignment                         | "iTBS_only", "switcher"  
age_at_first_phq      | integer   | Age at first PHQ-9 assessment            | 18–87  
sex                   | integer   | Sex                                      | 1=Male, 2=Female  
```

---

### File 3: excluded_patients.xlsx  
- **Location:** /data/final/ (in repository)  
- **Dimensions:** 28 rows × 4 columns  
- **Description:** Documentation of all patients excluded from analytic sample with reasons

**Variable Dictionary — excluded_patients.xlsx**

```  
Variable              | Type      | Description                              | Values  
----------------------|-----------|------------------------------------------|--------------------------------------------  
study_id              | character | De-identified patient ID                 | Various  
protocols             | character | Protocol combination codes used          | "4", "1,4", "1,2", "1,2,4", "1,2,5", "1,4,5"  
protocol_sequence     | character | Chronological protocol sequence          | e.g., "1 -> 1 -> 4 -> 4"  
reason                | character | Exclusion reason                         | See exclusion reasons below  
```

**Exclusion Reasons:**

```  
Reason                                                              | N  
--------------------------------------------------------------------|---  
Protocol 4 only (10Hz, older protocol)                              | TBD  
Mixed iTBS and 10Hz protocols                                       | TBD  
Mixed iTBS, bilateral, and 10Hz protocols                           | TBD  
Mixed iTBS, bilateral, and 10Hz/1Hz protocols                       | TBD  
Mixed iTBS, 10Hz, and 10Hz/1Hz protocols                            | TBD  
Started on protocol 2, switched to 1 (wrong direction) — TMS093     | 1  
Alternated between protocols 1 and 2 (no clean switch) — TMS084     | 1  
```

---

## 4. Final Analytic Sample Summary

```  
                   | Total | iTBS-only | Switchers  
-------------------|-------|-----------|----------  
N                  | 85    | 41        | 44  


### Sample Flow

```  
REDCap Export: 146 patients  
  ↓ −1 test row (TMS001)  
  ↓ −15 no usable PHQ-9 data  
Initial Clean: 130 patients  
  ↓ Filtered to course 1 only  
  ↓ Removed invalid/missing PHQ-9 data  
Intermediate Files: 113 patients  
  ↓ −28 excluded in 01_analytic_sample.R  
    - 26 used non-iTBS protocols (protocols 4, 5, or mixed)  
    - 1 switched bilateral → iTBS (wrong direction, TMS093)  
    - 1 alternated between protocols (no clean switch, TMS084)  
Final Sample: 85 patients (41 iTBS-only + 44 switchers)  
```

### Outcome Definitions

```  
Outcome        | Definition  
---------------|-------------------------------------------------------------  
Remitter       | Endpoint PHQ-9 ≤ 4  
Responder      | ≥50% reduction in PHQ-9 from baseline (but endpoint > 4)  
Non-Responder  | <50% reduction in PHQ-9 from baseline  
```

### Clinical Severity Cutoffs (PHQ-9)

```  
Score Range | Severity  
------------|--------------------  
0–4         | Minimal/none  
5–9         | Mild  
10–14       | Moderate  
15–19       | Moderately severe  
20–27       | Severe  
```

---

## 5. Key Decisions & Notes

1. **First course only:** Patients with multiple treatment courses (e.g., TMS021 course 1 and course 2) were filtered to course 1 only to avoid repeated-measures confounding  
2. **Protocol restriction:** Only protocol 1 (iTBS) and protocol 2 (iTBS/cTBS bilateral) were included; protocols 3 (cTBS), 4 (10Hz), and 5 (10Hz/1Hz bilateral) were excluded as they represent older or different treatment approaches  
3. **Switch direction:** Only unidirectional switches from iTBS → bilateral were included; TMS093 (bilateral → iTBS) was excluded as it represents the opposite clinical question  
4. **Clean switches only:** TMS084 was excluded for alternating between protocols (1→2→1→2) without a clear single switch point  
5. **PHQ-9 timing:** Visit 0 represents the first PHQ-9 collected during treatment (not necessarily session 1); some patients had pre-treatment consultation PHQ-9 scores stored separately  
6. **Missing item-level data:** Some individual PHQ-9 items have missing values (blank cells); total scores appear to have been calculated with available items

---

## 6. File Locations

```  
File                           | Location               | In Repository?  
-------------------------------|----------------------  |------------------  
PHQ9Project_6.15.26.csv        | /data/raw/             | No (PHI)  
ID_refs_PHI.xlsx               | Local only             | No (PHI)  
phq_scores_clean.xlsx          | /data/initial_cleaning/| Yes (initial cleaning)
patients_clean.xlsx            | Local only             | No (PHI) 
tms_clean_longitudinal.xlsx    | /data/intermediate/    | Yes (intermediate)
tms_clean_outcomes(in).csv     | /data/intermediate/    | Yes (intermediate)
analytic_sample_long.xlsx      | /data/final            | Yes (final) 
analytic_sample_summary.xlsx   | Local only.            | No (PHI) 
excluded_patients.xlsx         | /data/final/           | Yes (final) 
01_analytic_sample.R           | /code/                 | Yes (final)
```  
