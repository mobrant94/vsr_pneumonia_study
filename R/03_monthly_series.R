# Monthly series for Brazil: infant population (SINASC and SIM) and outcomes (SIH-RD and SIVEP-Gripe)
# Output: data/brazil_monthly_series.csv

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")

months <- seq(start_date, end_date, by = "month")
births <- fread("data/sinasc_births_uf_month.csv")
deaths <- fread("data/sim_infant_deaths_uf_month.csv")
births[, month := as.Date(month)]
deaths[, month := as.Date(month)]

# mean number of infants <n months in each month: births in the n previous months
# (half weight for the first and last month) minus infant deaths with the same weights
pop <- CJ(uf = ufs, month = months)
for (n in c(3, 6)) {
  w <- c(0.5, rep(1, n - 1), 0.5)
  x <- 0
  for (k in 0:n) {
    key <- pop[, .(uf, month = month %m-% months(k))]
    b <- births[key, on = .(uf, month), births]
    d <- deaths[key, on = .(uf, month), get(paste0("deaths_lt", n))]
    x <- x + w[k + 1] * (b - d)
  }
  pop[, paste0("lt", n) := x]
}
stopifnot(!anyNA(pop))
pop <- melt(pop[, .(lt3 = sum(lt3), lt6 = sum(lt6)), by = month], id.vars = "month",
            variable.name = "age_group", value.name = "infants")

# outcomes by month and age group (age groups are cumulative)
sih <- readRDS("data/processed/sih_cases.rds")
sivep <- readRDS("data/processed/sivep_cases.rds")

count_outcomes <- function(age) {
  h <- sih[get(age) == TRUE, .(bronchiolitis = sum(diagnosis == "bronchiolitis"),
                               pneumonia = sum(diagnosis == "pneumonia")), by = month]
  s <- sivep[get(age) == TRUE, .(sari_all = .N, sari_rsv = sum(rsv),
                                 icu_all = sum(icu %in% TRUE), icu_rsv = sum(icu %in% TRUE & rsv),
                                 imv_all = sum(imv %in% TRUE), imv_rsv = sum(imv %in% TRUE & rsv),
                                 death_all = sum(death %in% TRUE), death_rsv = sum(death %in% TRUE & rsv),
                                 known_all = sum(!is.na(death)), known_rsv = sum(!is.na(death) & rsv)), by = month]
  merge(merge(data.table(month = months, age_group = age), h, by = "month", all.x = TRUE), s, by = "month", all.x = TRUE)
}
br <- rbind(count_outcomes("lt3"), count_outcomes("lt6"))
br <- merge(pop, br, by = c("month", "age_group"))
setnafill(br, fill = 0, cols = names(br)[-(1:2)])

# sensitivity analysis: SIH admissions in 2026 corrected for processing delay
# p(L) = proportion of admissions in 2023-2025 processed within L months, by state of processing and diagnosis
last <- fread("data/processed/sih_last_processing.csv")
last[, last_processing := as.Date(last_processing)]
sih[, delay := months_between(month, processing)]
ref <- sih[year(month) %in% 2023:2025]
completeness <- CJ(uf_file = last$uf_file, diagnosis = c("bronchiolitis", "pneumonia"), L = 0:12)
completeness[, p := mapply(function(u, dx, l) {
  x <- ref[uf_file == u & diagnosis == dx, delay]
  if (length(x) < 200) x <- ref[diagnosis == dx, delay]   # few admissions: national distribution
  mean(x <= l)
}, uf_file, diagnosis, L)]

s26 <- merge(sih[month >= vaccine_date & month <= end_date], last, by = "uf_file")
s26[, L := pmin(12L, months_between(month, last_processing))]
s26 <- merge(s26, completeness, by = c("uf_file", "diagnosis", "L"))
nowcast <- rbind(s26[lt3 == TRUE, .(age_group = "lt3", observed = .N, corrected = sum(1 / p)), by = .(month, diagnosis)],
                 s26[lt6 == TRUE, .(age_group = "lt6", observed = .N, corrected = sum(1 / p)), by = .(month, diagnosis)])
fwrite(nowcast[order(age_group, diagnosis, month)], "output/tables/sih_delay_correction_2026.csv")

br[, `:=`(bronchiolitis_nc = as.numeric(bronchiolitis), pneumonia_nc = as.numeric(pneumonia))]
br[nowcast[diagnosis == "bronchiolitis"], on = .(month, age_group), bronchiolitis_nc := i.corrected]
br[nowcast[diagnosis == "pneumonia"], on = .(month, age_group), pneumonia_nc := i.corrected]

# infants aged 3 to <6 months = (<6 months) - (<3 months)
setorder(br, age_group, month)
b3 <- br[age_group == "lt3"]
b35 <- br[age_group == "lt6"]
for (v in setdiff(names(br), c("month", "age_group"))) b35[[v]] <- b35[[v]] - b3[[v]]
b35[, age_group := "3to5"]
br <- rbind(br, b35)

# deaths corrected for missing outcome: recorded deaths x cases / cases with known outcome
br[, death_corr_all := death_all * sari_all / pmax(1, known_all)]
br[, death_corr_rsv := death_rsv * sari_rsv / pmax(1, known_rsv)]

setorder(br, age_group, month)
fwrite(br, "data/brazil_monthly_series.csv")
