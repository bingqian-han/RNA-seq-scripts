# =============================================================================
# Figure S2: per-minute heart rate of the 102 sequenced individuals in each treatment, by population and shore level.
# Raw heart-rate data: Raw_heartrate_102ids.xlsx
# Input: FigureS2_per_minute_heartrate.csv
# Output: FigureS2.pdf
# =============================================================================

setwd("C:/path/to/RNAgithub/11heartrate/Formal_heartrate")
library(dplyr)
library(ggplot2)
# Read per-minute heart rate with sample information
heart <- read.csv("FigureS2_per_minute_heartrate.csv", fileEncoding = "UTF-8-BOM")
heart$population <- factor(heart$population, levels = c("North", "Central", "South"))
heart$treatment <- paste0(heart$temperature_C, " \u00B0C, ", ifelse(heart$duration_min == 15, "15 min", "3 h"))
heart$key <- paste(heart$ID, heart$temperature_C, heart$duration_min)
# Minute 0 is the calm heart rate and is not smoothed
minute0 <- heart[heart$minute == 0, ]
minute0$smoothed <- minute0$heart_rate
heart <- heart[heart$minute > 0, ]
# Smooth each individual's curve; stretches with missing values stay unsmoothed
heart$run <- ave(as.numeric(is.finite(heart$heart_rate)), heart$key, FUN = function(x) cumsum(c(1, diff(x) != 0)))
heart_15min <- heart[heart$duration_min == 15, ]
heart_3h <- heart[heart$duration_min == 180, ]
heart_15min$smoothed <- ave(heart_15min$heart_rate, heart_15min$key, heart_15min$run, FUN = function(x) if (length(x) >= 3 && all(is.finite(x))) signal::sgolayfilt(x, p = 1, n = 3) else x)
heart_3h$smoothed <- ave(heart_3h$heart_rate, heart_3h$key, heart_3h$run, FUN = function(x) if (length(x) >= 11 && all(is.finite(x))) signal::sgolayfilt(x, p = 2, n = 11) else x)
heart <- rbind(heart_15min, heart_3h)
minute0$run <- NA
heart <- rbind(minute0[, names(heart)], heart)
# Mean and SE per group and minute
heart_stats <- dplyr::summarise(dplyr::group_by(heart, treatment, population, shore_level, minute), n = sum(is.finite(smoothed)), smoothed_mean = mean(smoothed, na.rm = TRUE), smoothed_SE = sd(smoothed, na.rm = TRUE) / sqrt(n), .groups = "drop")

treatment_order <- c("17 \u00B0C, 15 min", "17 \u00B0C, 3 h", "31 \u00B0C, 15 min", "31 \u00B0C, 3 h")
heart_stats$treatment <- factor(heart_stats$treatment, levels = treatment_order)

# Plot mean +/- SE; the dashed line marks 15 min in the 31-degree, 3 h panel
fig_heart <- ggplot(heart_stats, aes(x = minute)) +
  geom_ribbon(aes(ymin = smoothed_mean - smoothed_SE, ymax = smoothed_mean + smoothed_SE, fill = population, group = interaction(population, shore_level)), alpha = 0.16) +
  geom_line(aes(y = smoothed_mean, colour = population, linetype = shore_level), linewidth = 0.5) +
  geom_segment(data = data.frame(treatment = factor("31 \u00B0C, 3 h", levels = treatment_order)), aes(x = 15, xend = 15, y = 0, yend = 108), linetype = "43", linewidth = 0.3, colour = "#444444") +
  geom_text(data = data.frame(treatment = factor("31 \u00B0C, 3 h", levels = treatment_order)), aes(x = 17.7, y = 99), label = "15 min", hjust = 0, vjust = 1, size = 6.5, size.unit = "pt", colour = "#444444") +
  facet_wrap(~ treatment, ncol = 2, scales = "free_x", axes = "all", axis.labels = "margins",
             labeller = labeller(treatment = function(x) gsub(" \u00B0C", "\u00B0C", x, fixed = TRUE))) +
  scale_x_continuous(breaks = function(limits) if (limits[2] <= 15) seq(0, 15, 3) else seq(0, 180, 30)) + scale_y_continuous(breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, 120), expand = FALSE) +
  scale_colour_manual(values = c(North = "#87aade", Central = "#C99A20", South = "#ff9999"), name = NULL) +
  scale_fill_manual(values = c(North = "#87aade", Central = "#ffe57f", South = "#ff9999"), guide = "none") +
  scale_linetype_manual(values = c(high = "solid", low = "53"), labels = c("High shore", "Low shore"), name = NULL) +
  labs(x = "Time (min)", y = "Heart rate (beats/min)") +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), plot.margin = margin(5.5, 12, 5.5, 5.5, unit = "pt"),
        strip.text.x = element_text(family = "Arial", size = 9.5, face = "bold",
                                    colour = "black", hjust = 0,
                                    margin = margin(t = 0, r = 0, b = 8, l = 0, unit = "pt")),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 6.5, colour = "black"), legend.position = "bottom", legend.text = element_text(size = 7), legend.key.width = unit(8, "mm"),
        panel.spacing.x = unit(6, "mm"), panel.spacing.y = unit(4, "mm"))
fig_heart

ggsave("FigureS2.pdf", fig_heart, width = 17, height = 12, units = "cm", device = cairo_pdf)
