# Figures: observed and expected monthly rates, and relative differences (forest plots)
# Outputs: output/figures/*.png

setwd("~/vsr_pneumonia_study")
source("R/00_setup.R")
library(ggplot2)

res <- fread("output/tables/its_results.csv")
cf <- fread("output/tables/its_monthly_expected.csv")
cf[, month := as.Date(month)]

labels <- c(bronchiolitis = "Acute bronchiolitis admission", pneumonia = "Pneumonia admission",
            sari_rsv = "RSV-associated SARI", icu_rsv = "ICU admission, RSV-associated SARI",
            imv_rsv = "IMV, RSV-associated SARI", death_corr_rsv = "Death, RSV-associated SARI",
            sari_all = "All-cause SARI", icu_all = "ICU admission, all-cause SARI",
            imv_all = "IMV, all-cause SARI", death_corr_all = "Death, all-cause SARI")
primary <- names(labels)[1:6]
key <- c("bronchiolitis", "pneumonia", "sari_rsv", "icu_rsv")
ages <- c(lt6 = "<6 months", lt3 = "<3 months", `3to5` = "3 to <6 months")

res[, `:=`(label = factor(labels[outcome], levels = rev(labels)), age = factor(ages[age_group], levels = ages),
           source = fifelse(outcome %in% c("bronchiolitis", "pneumonia"), "Hospital admissions (SIH-SUS)",
                            "RSV-associated SARI (SIVEP-Gripe)"),
           txt = sprintf("%.0f%% (%.0f%% to %.0f%%)", pct_difference, pct_lo, pct_hi))]
cf[, `:=`(label = factor(labels[outcome], levels = labels), age = factor(ages[age_group], levels = ages))]

# observed rate (solid) and expected rate in 2026 without vaccination (dashed, 95% CI)
series_plot <- function(z) {
  ggplot(z, aes(month, rate_observed)) +
    geom_ribbon(data = z[month >= vaccine_date], aes(ymin = rate_expected_lo, ymax = rate_expected_hi), fill = "#bd5b33", alpha = .18) +
    geom_line(data = z[month >= vaccine_date], aes(y = rate_expected, linetype = "Expected without vaccination"), colour = "#bd5b33", linewidth = .7) +
    geom_line(aes(linetype = "Observed"), colour = "#126681", linewidth = .55) +
    geom_point(colour = "#126681", size = .7) +
    geom_vline(xintercept = vaccine_date, linetype = 3, colour = "grey40") +
    scale_linetype_manual(values = c("Observed" = 1, "Expected without vaccination" = 2), name = NULL,
                          guide = guide_legend(override.aes = list(colour = c("#bd5b33", "#126681")))) +
    labs(x = NULL, y = "Cases per 1,000 infants per month") +
    theme_minimal(base_size = 10.5) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", hjust = 0))
}

# relative difference between observed and expected events, January-August 2026
forest_plot <- function(z) {
  ggplot(z, aes(pct_difference, label, colour = source)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_linerange(aes(xmin = pmax(pct_lo, -100), xmax = pmin(pct_hi, 60)), linewidth = .55) +
    geom_point(shape = 15, size = 2.6) +
    geom_text(aes(x = 70, label = txt), hjust = 0, size = 3.1, colour = "grey10") +
    scale_x_continuous(limits = c(-100, 125), breaks = seq(-80, 40, 40), labels = function(x) paste0(x, "%")) +
    scale_colour_manual(values = c("#6a9fb5", "#0b3d4f"), name = NULL) +
    labs(x = "Difference relative to expected (95% CI)", y = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0))
}

for (m in c("spline", "linear")) {
  prefix <- if (m == "spline") "Figure" else "Linear_Figure"
  z <- cf[model == m & outcome %in% key]
  z[, label := droplevels(label)]

  g1 <- series_plot(z[age_group == "lt6"]) + facet_wrap(~label, ncol = 2, scales = "free_y")
  g3 <- series_plot(z[age_group != "lt6"]) + facet_grid(label ~ age, scales = "free_y", labeller = label_wrap_gen(22))
  g2 <- forest_plot(res[model == m & age_group == "lt6" & outcome %in% primary])
  g4 <- forest_plot(res[model == m & age_group != "lt6" & outcome %in% primary]) + facet_wrap(~age, ncol = 1)

  ggsave(sprintf("output/figures/%s1_series_lt6.png", prefix), g1, width = 8, height = 6.2, dpi = 300, bg = "white")
  ggsave(sprintf("output/figures/%s2_difference_lt6.png", prefix), g2, width = 8, height = 3.6, dpi = 300, bg = "white")
  ggsave(sprintf("output/figures/%s3_series_subgroups.png", prefix), g3, width = 8, height = 9.2, dpi = 300, bg = "white")
  ggsave(sprintf("output/figures/%s4_difference_subgroups.png", prefix), g4, width = 8, height = 5.6, dpi = 300, bg = "white")
}

# supplementary figures: each outcome in the three age groups (primary model)
for (i in seq_along(labels)) {
  out <- names(labels)[i]
  g <- series_plot(cf[model == "spline" & outcome == out]) + facet_wrap(~age, ncol = 1, scales = "free_y") +
    labs(title = labels[[out]]) + theme(plot.title = element_text(face = "bold", size = 11))
  ggsave(sprintf("output/figures/FigureS%02d_%s.png", i, out), g, width = 7.5, height = 7, dpi = 300, bg = "white")
}
