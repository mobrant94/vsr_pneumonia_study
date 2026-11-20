# Breastfeeding among infants hospitalized with SARI (SIVEP-Gripe field M_AMAMENTA), January to August of each year
# Descriptive only: there is no population denominator by breastfeeding status
# Outputs: output/tables/breastfeeding_completeness.csv and output/tables/breastfeeding_severity.csv

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")

s <- readRDS("data/processed/sivep_cases.rds")
s <- s[month(month) <= month(end_date) & month <= end_date]
s[, year := year(month)]
s[, period := fifelse(year == 2026, "2026", "2022-2025")]

groups <- rbind(s[lt6 == TRUE][, `:=`(age_group = "lt6", sari = "all")],
                s[lt6 == TRUE & rsv][, `:=`(age_group = "lt6", sari = "rsv")],
                s[lt3 == TRUE][, `:=`(age_group = "lt3", sari = "all")],
                s[lt3 == TRUE & rsv][, `:=`(age_group = "lt3", sari = "rsv")])

completeness <- groups[, .(cases = .N, recorded = sum(!is.na(breastfed)), breastfed = sum(breastfed %in% TRUE)),
                       by = .(age_group, sari, year)]
completeness[, `:=`(pct_recorded = 100 * recorded / cases, pct_breastfed = 100 * breastfed / recorded)]
fwrite(completeness[order(age_group, sari, year)], "output/tables/breastfeeding_completeness.csv")

# % with Wilson 95% CI (prop.test without continuity correction)
wilson <- function(x) {
  x <- x[!is.na(x)]
  ci <- prop.test(sum(x), length(x), correct = FALSE)$conf.int
  sprintf("%.1f (%.1f-%.1f)", 100 * mean(x), 100 * ci[1], 100 * ci[2])
}
severity <- groups[!is.na(breastfed), .(cases = .N, icu = wilson(icu), imv = wilson(imv), death = wilson(death)),
                   by = .(age_group, sari, period, breastfed)]
severity <- severity[order(age_group, sari, period, -breastfed)]
fwrite(severity, "output/tables/breastfeeding_severity.csv")
print(severity[age_group == "lt6"])
