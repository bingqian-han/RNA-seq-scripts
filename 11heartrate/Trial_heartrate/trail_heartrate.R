# =============================================================================
# Figures S19 and S20: per-minute heart rate in the trial experiment (North and South), first 400 min.
# Raw heart-rate data: Raw_heartrate_<population>_<temperature>.csv (12 files)
# Input: Trial_per_minute_heartrate.csv (per-minute heart rate from the raw data)
# Output: FigureS19.pdf, FigureS20.pdf
# =============================================================================

setwd("C:/path/to/RNAgithub/11heartrate/Trial_heartrate")
library(ggplot2)
# Read per-minute heart rate; 17-degree rows are the room-temperature controls
heart <- read.csv("Trial_per_minute_heartrate.csv", fileEncoding = "UTF-8-BOM")
heart$country <- unname(c(North = "NO", South = "SP")[heart$population])
heart$individual <- paste(heart$file, heart$channel, sep = "|")
heart$clean_rate <- heart$heart_rate
heart$smoothed <- NA_real_
individuals <- unique(heart$individual)

# For each individual: zeros before the last non-zero value are set to missing, then the curve is smoothed
for (i in seq_along(individuals)) {
  idx <- which(heart$individual == individuals[i])
  d <- heart[idx, ]
  nonzero <- which(is.finite(d$heart_rate) & d$heart_rate != 0)
  zeros <- which(is.finite(d$heart_rate) & d$heart_rate == 0)
  last_nonzero <- max(nonzero)
  heart$clean_rate[idx[zeros[zeros < last_nonzero]]] <- NA_real_
  # Smooth each stretch of observed values, as in Figure S2
  rate <- heart$clean_rate[idx]
  valid <- which(is.finite(rate))
  segments <- split(valid, cumsum(c(TRUE, diff(valid) != 1)))
  for (segment in segments) {
    values <- rate[segment]
    if (length(segment) >= 11) values <- pmax(0, signal::sgolayfilt(values, p = 2, n = 11))
    values[rate[segment] == 0] <- 0
    heart$smoothed[idx[segment]] <- values
  }
}

# Mean and SE of the smoothed curves for each heat treatment
countries <- c(NO = "North", SP = "South")
colors <- c(NO = "#87aade", SP = "#ff9999")
panels <- data.frame(country = c(rep("NO", 5), rep("SP", 7)), temperature_C = c(17, 29, 30, 31, 32, 17, 30, 31, 32, 33, 35, 40))
heat_panels <- panels[panels$temperature_C != 17, ]
summary_curve <- list()
for (i in seq_len(nrow(heat_panels))) {
  d <- heart[heart$country == heat_panels$country[i] & heart$temperature_C == heat_panels$temperature_C[i], ]
  ids <- unique(d$individual)
  values <- matrix(NA_real_, 401, length(ids))
  values[cbind(d$minute + 1, match(d$individual, ids))] <- d$smoothed
  average <- rowMeans(values, na.rm = TRUE)
  se <- apply(values, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
  summary_curve[[i]] <- data.frame(country = heat_panels$country[i], temperature_C = heat_panels$temperature_C[i], minute = 0:400, mean = average, lower = pmax(0, average - se), upper = average + se)
}
summary_curve <- do.call(rbind, summary_curve)
# 17-degree baseline: average of each individual's mean heart rate
room <- heart[heart$temperature_C == 17, ]
room_means <- aggregate(room$clean_rate, room[c("country", "individual")], mean, na.rm = TRUE)
baseline <- tapply(room_means$x, room_means$country, mean)
# Y-axis upper limits
raw_upper <- 180
mean_upper <- ceiling(max(summary_curve$upper, baseline, na.rm = TRUE) / 10) * 10

# Panel titles, sample sizes and baseline lines
panels$title <- paste0(countries[panels$country], " | ", panels$temperature_C, "\u00B0C")
panels$key <- paste(panels$country, panels$temperature_C)
heart$panel <- factor(panels$title[match(paste(heart$country, heart$temperature_C), panels$key)], levels = panels$title)
summary_curve$panel <- factor(panels$title[match(paste(summary_curve$country, summary_curve$temperature_C), panels$key)], levels = panels$title[panels$temperature_C != 17])
panel_n <- aggregate(heart$individual, heart[c("country", "panel")], function(x) length(unique(x)))
names(panel_n)[3] <- "n"
room_reference <- data.frame(country = heat_panels$country, reference = as.numeric(baseline[heat_panels$country]))
room_reference$panel <- factor(panels$title[match(paste(heat_panels$country, heat_panels$temperature_C), panels$key)], levels = levels(summary_curve$panel))

figure_theme <- theme_classic(base_size = 11, base_family = "Arial") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, face = "bold", colour = "black", hjust = 0, margin = margin(b = 5)),
        axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11), axis.ticks = element_line(linewidth = .3),
        axis.line = element_line(linewidth = .35), panel.spacing = grid::unit(6, "mm"), plot.margin = margin(7, 11, 7, 7), legend.position = "none")

# Figure S19: every individual's curve; grey = 17-degree
FigureS19 <- ggplot(heart, aes(minute, heart_rate, group = individual)) +
  geom_vline(xintercept = c(15, 180), linetype = "dotted", linewidth = .35, colour = "grey45") +
  geom_line(aes(colour = ifelse(temperature_C == 17, "room", country)), linewidth = .28, na.rm = TRUE) +
  geom_text(data = panel_n, aes(x = 390, y = raw_upper * .97, label = paste0("n = ", n)), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 10, size.unit = "pt", family = "Arial") +
  facet_wrap(~ panel, ncol = 3, axes = "all", axis.labels = "all") +
  scale_colour_manual(values = c(colors, room = "#555555"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 400, 100)) + scale_y_continuous(breaks = seq(0, raw_upper, 50)) +
  coord_cartesian(xlim = c(0, 400), ylim = c(0, raw_upper), expand = FALSE) +
  labs(x = "Time (min)", y = "Heart rate (beats/min)") + figure_theme

# Figure S20: mean +/- SE; dashed line = 17-degree baseline
FigureS20 <- ggplot(summary_curve, aes(minute, mean)) +
  geom_vline(xintercept = c(15, 180), linetype = "dotted", linewidth = .35, colour = "grey45") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = country), alpha = .23, colour = NA, na.rm = TRUE) +
  geom_line(aes(colour = country), linewidth = .65, na.rm = TRUE) +
  geom_hline(data = room_reference, aes(yintercept = reference), linetype = "longdash", linewidth = .5, colour = "black") +
  geom_text(data = panel_n[panel_n$panel %in% levels(summary_curve$panel), ], aes(x = 390, y = mean_upper * .97, label = paste0("n = ", n)), inherit.aes = FALSE, hjust = 1, vjust = 1, size = 10, size.unit = "pt", family = "Arial") +
  facet_wrap(~ panel, ncol = 2, axes = "all", axis.labels = "all") +
  scale_colour_manual(values = colors, guide = "none") + scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 400, 100)) + scale_y_continuous(breaks = seq(0, mean_upper, 20)) +
  coord_cartesian(xlim = c(0, 400), ylim = c(0, mean_upper), expand = FALSE) +
  labs(x = "Time (min)", y = "Heart rate (beats/min)") + figure_theme

FigureS19
ggsave("FigureS19.pdf", FigureS19, width = 18, height = 23, units = "cm", device = cairo_pdf)

FigureS20
ggsave("FigureS20.pdf", FigureS20, width = 18, height = 22.5, units = "cm", device = cairo_pdf)
