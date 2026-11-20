# Hospitalized SARI among infants <6 months (SIVEP-Gripe, OpenDataSUS)
# Output: data/processed/sivep_cases.rds (one row per hospitalized case)

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")

raw_dir <- "data/raw"   # folder with the INFLUD files (files already downloaded are not downloaded again)

# files used in the study (one per year; the date in the name is the extraction date)
files <- c(`2022` = "INFLUD22-23-03-2026.csv", `2023` = "INFLUD23-23-03-2026.csv",
           `2024` = "INFLUD24-23-03-2026.csv", `2025` = "INFLUD25-19-10-2026.csv",
           `2026` = "INFLUD26-09-11-2026.csv")
url <- "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG"

options(timeout = 7200)
for (y in names(files)) {
  f <- file.path(raw_dir, files[[y]])
  if (!file.exists(f)) download.file(paste(url, y, files[[y]], sep = "/"), f, mode = "wb")
}

# the month of extraction is incomplete
extraction <- dmy(sub("INFLUD26-(.*)[.]csv", "\\1", files[["2026"]]))
last_month <- floor_date(extraction, "month") %m-% months(1)

# dates come as YYYY-MM-DD (or DD/MM/YYYY in some extractions)
read_date <- function(x) {
  d <- as.Date(substr(x, 1, 10), format = "%Y-%m-%d")
  d[is.na(d)] <- as.Date(substr(x[is.na(d)], 1, 10), format = "%d/%m/%Y")
  d
}

vars <- c("NU_NOTIFIC", "DT_SIN_PRI", "DT_NASC", "NU_IDADE_N", "TP_IDADE", "SG_UF", "CO_MUN_RES",
          "HOSPITAL", "EVOLUCAO", "PCR_VSR", "AN_VSR", "UTI", "SUPORT_VEN", "M_AMAMENTA",
          "DT_DIGITA", "DT_PCR", "DT_RES_AN", "DT_EVOLUCA")

cases <- list()
for (y in names(files)) {
  s <- fread(file.path(raw_dir, files[[y]]), sep = ";", select = vars, colClasses = "character")

  s[, onset := read_date(DT_SIN_PRI)]
  s[, birth := read_date(DT_NASC)]
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

  # dates used in the reporting-delay analysis: entry into the system (DT_DIGITA) and,
  # for RSV, the later of the date of entry and the date of the positive result
  s[, entered := read_date(DT_DIGITA)]
  s[, rsv_result := fifelse(PCR_VSR %in% "1", read_date(DT_PCR), read_date(DT_RES_AN))]
  s[, rsv_known := fifelse(rsv, pmax(entered, rsv_result, na.rm = TRUE), as.Date(NA))]
  s[, outcome_date := read_date(DT_EVOLUCA)]
  s[, extraction := dmy(sub("INFLUD..-(.*)[.]csv", "\\1", files[[y]]))]

  cases[[y]] <- s[, .(uf, month, onset, lt3, lt6, rsv, icu, imv, death, breastfed,
                      entered, rsv_known, outcome_date, extraction)]
  message(y, ": ", nrow(s), " hospitalized cases <6 months")
}

cases <- rbindlist(cases)
saveRDS(cases, "data/processed/sivep_cases.rds")
