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
    ├── code/                       # R scripts (numbered in execution order)  
    │  
    ├── data/  
    │   ├── raw/                    # Original REDCap export (not shared)  
    │   ├── intermediate/           # Cleaned intermediate files  
    │   └── final/                  # Analysis-ready dataset  
    │  
    ├── docs/                       # Data provenance, variable dictionary  
    │  
    └── output/  
        ├── figures/                # Publication figures  
        └── tables/                 # Results tables  


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
