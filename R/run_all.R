# Runs the whole analysis in sequence

setwd("~/vsr_pneumonia_study")

for (f in c("R/01_sivep_gripe.R", "R/02_sih.R", "R/03_monthly_series.R",
            "R/04_its_models.R", "R/05_breastfeeding.R", "R/06_figures.R")) {
  message("---- ", f)
  source(f, encoding = "UTF-8")
}
