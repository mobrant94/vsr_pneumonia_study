# Hospitalized SARI among infants <6 months (SIVEP-Gripe, OpenDataSUS)
# Output: data/processed/sivep_cases.rds (one row per hospitalized case)

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")

# files used in the study (one per year; the date in the name is the extraction date)
files <- c(`2022` = "INFLUD22-23-03-2026.csv", `2023` = "INFLUD23-23-03-2026.csv",
           `2024` = "INFLUD24-23-03-2026.csv", `2025` = "INFLUD25-19-10-2026.csv",
           `2026` = "INFLUD26-09-11-2026.csv")
url <- "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG"

options(timeout = 7200)
for (y in names(files)) {
  f <- file.path("data/raw", files[[y]])
  if (!file.exists(f)) download.file(paste(url, y, files[[y]], sep = "/"), f, mode = "wb")
}

# the month of extraction is incomplete
extraction <- dmy(sub("INFLUD26-(.*)[.]csv", "\\1", files[["2026"]]))
last_month <- floor_date(extraction, "month") %m-% months(1)

vars <- c("NU_NOTIFIC", "DT_SIN_PRI", "DT_NASC", "NU_IDADE_N", "TP_IDADE", "SG_UF", "CO_MUN_RES",
          "HOSPITAL", "EVOLUCAO", "PCR_VSR", "AN_VSR", "UTI", "SUPORT_VEN", "M_AMAMENTA")

cases <- list()
for (y in names(files)) {
  s <- fread(file.path("data/raw", files[[y]]), sep = ";", select = vars, colClasses = "character")

  s[, onset := as.Date(substr(DT_SIN_PRI, 1, 10), format = "%Y-%m-%d")]
  s[, birth := as.Date(substr(DT_NASC, 1, 10), format = "%Y-%m-%d")]
  s[, month := floor_date(onset, "month")]
  s[, uf := uf_from_mun(CO_MUN_RES)]
  s[is.na(uf) & SG_UF %in% ufs, uf := SG_UF]

  # hospitalized cases with symptom onset in the year of the file
  s <- s[HOSPITAL == "1" & year(onset) == as.integer(y) & month >= start_date & month <= last_month & uf %in% ufs]

  # age at symptom onset; without date of birth, the recorded age (TP_IDADE 1 = days, 2 = months)
  s[, age := suppressWarnings(as.integer(NU_IDADE_N))]
  s[, lt3 := younger_than(birth, onset, 3)]
  s[, lt6 := younger_than(birth, onset, 6)]
  s[is.na(birth), lt3 := (TP_IDADE == "1" & age >= 0 & age < 90) | (TP_IDADE == "2" & age >= 0 & age < 3)]
  s[is.na(birth), lt6 := (TP_IDADE == "1" & age >= 0 & age < 183) | (TP_IDADE == "2" & age >= 0 & age < 6)]
  s <- s[lt6 == TRUE]
  stopifnot(!anyDuplicated(s[NU_NOTIFIC != ""], by = "NU_NOTIFIC"))

  # outcomes (unknown = NA)
  s[, rsv := PCR_VSR %in% "1" | AN_VSR %in% "1"]
  s[, icu := fifelse(UTI %in% c("1", "2"), UTI == "1", NA)]
  s[, imv := fifelse(SUPORT_VEN %in% c("1", "2", "3"), SUPORT_VEN == "1", NA)]
  s[, death := fifelse(EVOLUCAO %in% c("1", "2", "3"), EVOLUCAO %in% c("2", "3"), NA)]
  s[, breastfed := fifelse(M_AMAMENTA %in% c("1", "2"), M_AMAMENTA == "1", NA)]

  cases[[y]] <- s[, .(uf, month, onset, lt3, lt6, rsv, icu, imv, death, breastfed)]
  message(y, ": ", nrow(s), " hospitalized cases <6 months")
}

cases <- rbindlist(cases)
saveRDS(cases, "data/processed/sivep_cases.rds")
