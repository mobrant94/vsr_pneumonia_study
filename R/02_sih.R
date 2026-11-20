# Admissions for pneumonia (J12-J18) and acute bronchiolitis (J21) among infants <6 months (SIH-RD)
# Files are downloaded by state and processing month with the microdatasus package (>= 3.0.0)
# Outputs: data/processed/sih_cases.rds and data/processed/sih_last_processing.csv

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")
library(microdatasus)

last_processing <- as.Date("2026-09-01")   # last processing month available at the time of the analysis
vars <- c("N_AIH", "IDENT", "MUNIC_RES", "DT_INTER", "NASC", "DIAG_PRINC", "COD_IDADE", "IDADE", "ANO_CMPT", "MES_CMPT")

dir.create("data/raw/sih", showWarnings = FALSE)
for (uf in ufs) for (m in seq(start_date, last_processing, by = "month")) {
  m <- as.Date(m, origin = "1970-01-01")
  f <- sprintf("data/raw/sih/RD%s%s.rds", uf, format(m, "%y%m"))
  if (file.exists(f)) next
  d <- try(fetch_datasus(year_start = year(m), month_start = month(m), year_end = year(m), month_end = month(m),
                         uf = uf, information_system = "SIH-RD", vars = vars), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || nrow(d) == 0) {
    message("not available: ", uf, " ", format(m, "%m/%Y"))
    next
  }
  saveRDS(d, f)
}

sih <- list(); processed <- list()
for (f in list.files("data/raw/sih", full.names = TRUE)) {
  a <- as.data.table(readRDS(f))
  uf_file <- substr(basename(f), 3, 4)
  a[, processing := as.Date(sprintf("%s-%02d-01", ANO_CMPT, as.integer(MES_CMPT)))]
  processed[[f]] <- data.table(uf_file = uf_file, processing = max(a$processing))

  a[, admission := as.Date(DT_INTER, format = "%Y%m%d")]
  a[, birth := as.Date(NASC, format = "%Y%m%d")]
  a[, age := suppressWarnings(as.integer(IDADE))]
  a[, icd := toupper(gsub("[^A-Za-z0-9]", "", DIAG_PRINC))]
  a[, uf := uf_from_mun(MUNIC_RES)]
  a <- a[IDENT == "1" & !is.na(admission) & admission >= start_date & uf %in% ufs & grepl("^J1[2-8]|^J21", icd)]

  # age at admission; without date of birth, the recorded age (COD_IDADE 2 = days, 3 = months)
  a[, lt3 := younger_than(birth, admission, 3)]
  a[, lt6 := younger_than(birth, admission, 6)]
  a[is.na(birth), lt3 := (COD_IDADE == "2" & age >= 0 & age < 90) | (COD_IDADE == "3" & age >= 0 & age <= 2)]
  a[is.na(birth), lt6 := (COD_IDADE == "2" & age >= 0 & age < 183) | (COD_IDADE == "3" & age >= 0 & age <= 5)]
  a <- a[lt6 == TRUE]

  a[, diagnosis := fifelse(grepl("^J21", icd), "bronchiolitis", "pneumonia")]
  a[, month := floor_date(admission, "month")]
  sih[[f]] <- a[, .(N_AIH, uf, month, diagnosis, lt3, lt6, processing, uf_file)]
}
sih <- rbindlist(sih)
processed <- rbindlist(processed)[, .(last_processing = max(processing)), by = uf_file]

# the same AIH can be submitted in more than one processing month: keep the first
setorder(sih, N_AIH, processing)
sih <- sih[is.na(N_AIH) | N_AIH == "" | !duplicated(N_AIH)]

saveRDS(sih, "data/processed/sih_cases.rds")
fwrite(processed, "data/processed/sih_last_processing.csv")
print(processed)
