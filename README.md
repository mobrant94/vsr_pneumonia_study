# Maternal RSV vaccination and respiratory hospitalizations among infants in Brazil, 2022–2026

R code for the nationwide interrupted time-series analysis comparing hospital admissions and severe outcomes among infants younger than 6 months in January–August 2026, after the introduction of maternal RSV vaccination in the Brazilian Unified Health System, with those expected from the 2022–2025 trend and seasonality.

## Data

All sources are public. Monthly counts are aggregated for Brazil; no record linkage was performed.

| Source | Content | Script |
|:---|:---|:---|
| SIVEP-Gripe (OpenDataSUS) | Hospitalized SARI, RSV confirmation, ICU admission, invasive mechanical ventilation, death, breastfeeding | `01_sivep_gripe.R` |
| SIH-RD (DATASUS, via [microdatasus](https://github.com/rfsaldanha/microdatasus)) | Admissions for pneumonia (ICD-10 J12–J18) and acute bronchiolitis (J21) | `02_sih.R` |
| SINASC and SIM (DATASUS) | Live births and infant deaths by state and month (`data/`) | `03_monthly_series.R` |

`data/brazil_monthly_series.csv` contains the monthly series used in the models (infants <3 months, <6 months, and 3 to <6 months), so scripts 04 and 06 can be run without downloading the microdata.

## Scripts

| Script | Description |
|:---|:---|
| `R/00_setup.R` | Packages, study period, and helper functions |
| `R/01_sivep_gripe.R` | Hospitalized SARI cases among infants <6 months |
| `R/02_sih.R` | Pneumonia and bronchiolitis admissions among infants <6 months |
| `R/03_monthly_series.R` | Infant population, monthly series, correction of deaths for missing outcome, and SIH processing-delay correction |
| `R/04_its_models.R` | Negative binomial interrupted time-series models (spline trend as primary analysis and sensitivity analyses), expected events, and differences |
| `R/05_breastfeeding.R` | Breastfeeding among hospitalized infants (descriptive) |
| `R/06_figures.R` | Figures |

Set the working directory at the top of each script (`setwd`) to the folder of the repository and run the scripts in order. Results are written to `output/`.

R packages: `data.table`, `lubridate`, `MASS`, `sandwich`, `ggplot2`, and `microdatasus` (>= 3.0.0, only for `02_sih.R`). Analyses were performed in R 4.3.3.
