# Interrupted time series: January-August 2026 vs. expected from the 2022-2025 trend and seasonality
# Negative binomial regression with log(infants) as offset and Newey-West 95% CI (12-month lag)
# Primary model: piecewise linear (spline) trend with knots at January 2023 and January 2025
# Outputs (output/tables): its_results.csv, its_monthly_expected.csv, its_diagnostics.csv

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")
library(MASS)
library(sandwich)

br <- fread("data/brazil_monthly_series.csv")
br[, month := as.Date(month)]
br[, time := months_between(start_date, month)]
br[, `:=`(month_num = month(month), post = as.integer(month >= vaccine_date),
          t = time / 12, k23 = pmax(0, time - 12) / 12, k25 = pmax(0, time - 36) / 12,
          time_c = (time - mean(time)) / 12)]

outcomes <- c("bronchiolitis", "pneumonia", "sari_rsv", "icu_rsv", "imv_rsv", "death_corr_rsv",   # primary
              "sari_all", "icu_all", "imv_all", "death_corr_all")                                  # secondary

f_spline <- y ~ t + k23 + k25 + factor(month_num) + post
f_linear <- y ~ time_c + factor(month_num) + post
models <- list(
  spline         = list(f = f_spline, start = "2022-01-01"),
  spline_delay   = list(f = f_spline, start = "2022-01-01"),   # SIH admissions corrected for processing delay
  linear         = list(f = f_linear, start = "2022-01-01"),
  linear_2023    = list(f = f_linear, start = "2023-01-01"),
  linear_2024    = list(f = f_linear, start = "2024-01-01"),
  quadratic      = list(f = y ~ time_c + I(time_c^2) + factor(month_num) + post, start = "2022-01-01"),
  no_trend       = list(f = y ~ factor(month_num) + post, start = "2022-01-01"))

results <- list(); expected <- list(); diagnostics <- list()

for (m in names(models)) for (age in c("lt6", "lt3", "3to5")) for (out in outcomes) {
  var <- if (m == "spline_delay") paste0(out, "_nc") else out
  if (!var %in% names(br)) next
  d <- br[age_group == age & month >= as.Date(models[[m]]$start)]
  d[, y := round(get(var))]

  f <- update(models[[m]]$f, . ~ . + offset(log(infants)))
  fit <- tryCatch(suppressWarnings(glm.nb(f, data = d, control = glm.control(maxit = 100))), error = function(e) NULL)
  family <- "negative binomial"
  if (is.null(fit) || !is.finite(fit$theta) || fit$theta > 1e6) {
    fit <- glm(f, family = quasipoisson(), data = d)
    family <- "quasi-Poisson"
  }
  V <- NeweyWest(fit, lag = nw_lag, prewhite = FALSE, adjust = TRUE)

  # adjusted rate ratio (2026 indicator)
  b <- coef(fit)["post"]; se <- sqrt(V["post", "post"])
  rr <- exp(b); rr_lo <- exp(b - 1.96 * se); rr_hi <- exp(b + 1.96 * se)

  # crude rate ratio: January-August 2026 vs. the same months of the reference years
  pre <- d[post == 0 & month_num <= month(end_date)]; pos <- d[post == 1]
  crude <- (sum(pos$y) / sum(pos$infants)) / (sum(pre$y) / sum(pre$infants))
  crude_se <- sqrt(1 / sum(pos$y) + 1 / sum(pre$y))

  # expected events = observed / rate ratio; differences derived from the same rate ratio
  obs <- sum(pos$y)
  results[[length(results) + 1]] <- data.table(
    model = m, age_group = age, outcome = out, family = family,
    events_ref = sum(pre$y), events_2026 = obs,
    # events from January to August per 1,000 infants (infant-months / 8 = mean number of infants)
    rate_ref = 1000 * 8 * sum(pre$y) / sum(pre$infants), rate_2026 = 1000 * 8 * obs / sum(pos$infants),
    crude_rr = crude, crude_lo = exp(log(crude) - 1.96 * crude_se), crude_hi = exp(log(crude) + 1.96 * crude_se),
    rr = rr, rr_lo = rr_lo, rr_hi = rr_hi, p = 2 * pnorm(-abs(b / se)),
    expected = obs / rr, expected_lo = obs / rr_hi, expected_hi = obs / rr_lo,
    difference = obs - obs / rr, difference_lo = obs - obs / rr_lo, difference_hi = obs - obs / rr_hi,
    pct_difference = 100 * (rr - 1), pct_lo = 100 * (rr_lo - 1), pct_hi = 100 * (rr_hi - 1))

  if (m %in% c("spline", "linear")) {
    # monthly expected rates without the vaccine (post = 0), with Newey-West 95% CI
    cf <- copy(d)[, post := 0]
    X <- model.matrix(delete.response(terms(fit)), data = cf)
    eta <- as.vector(X %*% coef(fit))
    se_eta <- sqrt(rowSums((X %*% V) * X))
    expected[[length(expected) + 1]] <- data.table(
      model = m, age_group = age, outcome = out, month = d$month, observed = d$y, infants = d$infants,
      rate_observed = 1000 * d$y / d$infants, rate_expected = 1000 * exp(eta),
      rate_expected_lo = 1000 * exp(eta - 1.96 * se_eta), rate_expected_hi = 1000 * exp(eta + 1.96 * se_eta))

    r <- residuals(fit, type = "pearson")
    diagnostics[[length(diagnostics) + 1]] <- data.table(
      model = m, age_group = age, outcome = out, family = family,
      theta = if (family == "negative binomial") fit$theta else NA,
      pearson_dispersion = sum(r^2) / df.residual(fit), acf_lag1 = acf(r, plot = FALSE, lag.max = 1)$acf[2])
  }
}

results <- rbindlist(results)
fwrite(results, "output/tables/its_results.csv")
fwrite(rbindlist(expected), "output/tables/its_monthly_expected.csv")
fwrite(rbindlist(diagnostics), "output/tables/its_diagnostics.csv")

# primary analysis, infants <6 months
print(results[model == "spline" & age_group == "lt6",
              .(outcome, observed = events_2026, expected = round(expected),
                difference = sprintf("%.0f%% (%.0f to %.0f)", pct_difference, pct_lo, pct_hi))])
