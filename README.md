# TMS-Switching-Study
# To Switch or Not to Switch?: Protocol Change from iTBS to Bilateral Transcranial Magnetic Stimulation is Associated with Improved PHQ-9 Trajectory

## Authors  
Chiara Mosca, Ashkan Davani, Noelle Arcaro, Andrea Joanlanne, Anil Malhotra, Miklos Argyelan

Zucker Hillside Hospital & Feinstein Institutes for Medical Research

## Overview  
Please read data_provenance.md in docs folder for data cleaning pipeline. This repository contains analysis code and data provenance documentation for a retrospective study examining whether switching non-responsive iTBS patients to bilateral TBS is associated with improved PHQ-9 trajectories  

## Repository Structure  
    TMS-Switching-Study/  
    │  
    ├── code/   
    │  
    ├── data/  
    │   ├── raw/  
    │   │   └── README.md                          # Placeholder (raw data not shared - PHI)  
    │   │  
    │   ├── initial_cleaning/  
    │   │   ├── .gitkeep  
    │   │   └── phq_scores_clean.xlsx              # De-identified PHQ-9 item-level data  
    │   │  
    │   ├── intermediate/  
    │   │   ├── .gitkeep  
    │   │   ├── tms_clean_longitudinal.xlsx         # Longitudinal data with visit numbering  
    │   │   └── tms_clean_outcomes(in).csv          # Patient-level outcomes summary  
    │   │  
    │   └── final/  
    │       ├── .gitkeep  
    │       ├── analytic_sample_long.xlsx           # Final longitudinal analytic dataset (N=85)  
    │       └── excluded_patients.xlsx              # Excluded patients with reasons (N=28)  
    │  
    ├── docs/  
    │   └── data_provenance.md  
    │  
    ├── output/  
    │   ├── figures/                                # Publication figures 
    │   └── tables/                                 # Results tables 
    │  
    ├── .gitignore  
    ├── LICENSE  
    └── README.md  


## Reproduction  
1. Clone this repository  
2. Place raw REDCap data in `data/raw/` (requires IRB-approved access)  
3. Run scripts in numerical order:  
   - `00_data_cleaning.R` → generates cleaned datasets  
   - `01_sample_characteristics.R` → Table 1  
   - `02_trajectory_analysis.R` → Figures 1-3, HLM models  
   - `03_item_analysis.R` → Figures 4-5, item-level models  
   - `04_sensitivity_analyses.R` → supplementary analyses

## Requirements  
- R >= 4.3.0  
- Packages: tidyverse, lme4, lmerTest, here, janitor, ggplot2

Install all dependencies:  
```r  
install.packages(c("tidyverse", "lme4", "lmerTest", "here", "janitor"))  
