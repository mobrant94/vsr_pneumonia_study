# Packages, study parameters and a few helper functions used by all scripts

library(data.table)
library(lubridate, warn.conflicts = FALSE)

start_date   <- as.Date("2022-01-01")   # first month of the series
vaccine_date <- as.Date("2026-01-01")   # first month after the introduction of maternal RSV vaccination
end_date     <- as.Date("2026-08-01")   # last month analyzed
nw_lag       <- 12                      # lag (months) of the Newey-West variance

ufs <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT","PA",
         "PB","PE","PI","PR","RJ","RN","RS","RO","RR","SC","SP","SE","TO")
uf_code <- c(AC="12",AL="27",AM="13",AP="16",BA="29",CE="23",DF="53",ES="32",GO="52",
             MA="21",MG="31",MS="50",MT="51",PA="15",PB="25",PE="26",PI="22",PR="41",
             RJ="33",RN="24",RS="43",RO="11",RR="14",SC="42",SP="35",SE="28",TO="17")

# state of residence from the IBGE municipality code
uf_from_mun <- function(x) names(uf_code)[match(substr(x, 1, 2), uf_code)]

# event occurred before the child was n months old
younger_than <- function(birth, date, n) !is.na(birth) & birth <= date & date < birth %m+% months(n)

months_between <- function(a, b) 12L * (year(b) - year(a)) + month(b) - month(a)

for (d in c("data/raw", "data/processed", "output/tables", "output/figures"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
