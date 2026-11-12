setwd("C:/path/to/RNAgithub/12temperature_loggers/")

# =============================================================================
# Rock temperature at the North and South sites, summer 2023 (temperature loggers, high and low shore)
#
# A: temperature every 30 mins from 1 June to 31 August (mean of the two loggers at each shore level).
# B: mean temperature between 00:00 and 05:00 of each day.
# C: number of days with at least 3 hours above 31°C (at least 7 consecutive readings above 31°C within one day).
#
# Input:  North_logger.txt (North), South_logger.txt (South); columns date, time, lowtemp1, lowtemp2, hightemp1, hightemp2
# Output: FigureS1.pdf
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)


############################################################
#Part 1: Temperature of each shore level (mean of the two loggers), 1 June to 31 August 2023

# Logger times are local time (CEST); tz = "UTC" only prevents any time zone conversion.
north <- read.delim("North_logger.txt", check.names = FALSE, stringsAsFactors = FALSE)
south <- read.delim("South_logger.txt", check.names = FALSE, stringsAsFactors = FALSE)
north$timestamp <- as.POSIXct(paste(north$date, north$time), format = "%Y/%m/%d %H:%M", tz = "UTC")
south$timestamp <- as.POSIXct(paste(south$date, south$time), format = "%Y/%m/%d %H:%M", tz = "UTC")
start_time <- as.POSIXct("2023-06-01 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-09-01 00:00:00", tz = "UTC")
north <- north[north$timestamp >= start_time & north$timestamp < end_time, ]
south <- south[south$timestamp >= start_time & south$timestamp < end_time, ]

paired <- rbind(data.frame(site = "North", shore = "Low shore", timestamp = north$timestamp, temperature = (north$lowtemp1 + north$lowtemp2) / 2),
                data.frame(site = "North", shore = "High shore", timestamp = north$timestamp, temperature = (north$hightemp1 + north$hightemp2) / 2),
                data.frame(site = "South", shore = "Low shore", timestamp = south$timestamp, temperature = (south$lowtemp1 + south$lowtemp2) / 2),
                data.frame(site = "South", shore = "High shore", timestamp = south$timestamp, temperature = (south$hightemp1 + south$hightemp2) / 2))
paired$site <- factor(paired$site, levels = c("North", "South"))
paired$shore <- factor(paired$shore, levels = c("Low shore", "High shore"))
paired$date <- as.Date(paired$timestamp, tz = "UTC")
paired$hour <- as.numeric(format(paired$timestamp, "%H", tz = "UTC")) + as.numeric(format(paired$timestamp, "%M", tz = "UTC")) / 60

# B: mean temperature between 00:00 and 05:00 (11 readings) of each day, drawn at 02:30.
nightly <- paired %>% dplyr::filter(hour <= 5) %>% dplyr::group_by(site, shore, date) %>%
  dplyr::summarise(temperature = mean(temperature), .groups = "drop") %>%
  dplyr::mutate(timestamp = as.POSIXct(date, tz = "UTC") + 2.5 * 3600)

# C: longest run of consecutive readings above 31°C in each day; a day counts when the run has at least 7 readings
# (7 readings at 30-min intervals span 3 hours).
daily <- paired %>% dplyr::arrange(site, shore, timestamp) %>% dplyr::group_by(site, shore, date) %>%
  dplyr::summarise(longest_run = max(c(0, rle(temperature > 31)$lengths[rle(temperature > 31)$values])), .groups = "drop")
counts <- daily %>% dplyr::group_by(site, shore) %>%
  dplyr::summarise(n_days = dplyr::n(), days = sum(longest_run >= 7), percent = 100 * mean(longest_run >= 7), .groups = "drop")
# Number and percentage of such days (Figure S1C).
counts
counts$label <- sprintf("%d\n(%.1f%%)", counts$days, counts$percent)


############################################################
#Part 2: Figure S1

colours <- c("High shore" = "#DC5E83", "Low shore" = "#00A1D5")
date_breaks <- as.POSIXct(c("2023-06-01", "2023-07-01", "2023-08-01", "2023-08-31"), tz = "UTC")
paper_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.25, colour = "black"),
        axis.ticks.length = grid::unit(1.1, "mm"), axis.text = element_text(size = 7, colour = "black"), axis.title = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_text(margin = margin(r = 4)), strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold", hjust = 0), panel.spacing.y = grid::unit(2.5, "mm"),
        legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.width = grid::unit(7, "mm"),
        legend.key.height = grid::unit(3, "mm"), legend.margin = margin(0, 0, 1, 0, unit = "mm"), plot.tag = element_text(size = 10, face = "bold"),
        plot.margin = margin(2, 10, 2, 2, unit = "mm"))

pA <- ggplot(paired, aes(timestamp, temperature, colour = shore)) +
  geom_line(linewidth = 0.22) +
  geom_hline(yintercept = c(17, 31), colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = end_time, y = 17, label = "17\u00b0C", hjust = -0.15, size = 2.3, family = "Arial", colour = "grey25") +
  annotate("text", x = end_time, y = 31, label = "31\u00b0C", hjust = -0.15, size = 2.3, family = "Arial", colour = "grey25") +
  facet_wrap(~site, ncol = 1) +
  scale_x_datetime(limits = c(start_time, end_time), breaks = date_breaks, labels = c("1 Jun", "1 Jul", "1 Aug", "31 Aug"), expand = expansion(mult = 0)) +
  scale_colour_manual(values = colours, breaks = c("High shore", "Low shore")) +
  scale_y_continuous(breaks = seq(0, 50, 10), expand = expansion(mult = 0)) + coord_cartesian(ylim = c(0, 55), clip = "off") +
  labs(y = "Temperature (\u00b0C)") + paper_theme + guides(colour = guide_legend(override.aes = list(linewidth = 0.6)))

pB <- ggplot(nightly, aes(timestamp, temperature, colour = shore)) +
  geom_line(linewidth = 0.3) +
  geom_hline(yintercept = 17, colour = "grey35", linetype = "dashed", linewidth = 0.25) +
  annotate("text", x = end_time, y = 17, label = "17\u00b0C", hjust = -0.15, size = 2.3, family = "Arial", colour = "grey25") +
  facet_wrap(~site, ncol = 1) +
  scale_x_datetime(limits = c(start_time, end_time), breaks = date_breaks, labels = c("1 Jun", "1 Jul", "1 Aug", "31 Aug"), expand = expansion(mult = 0)) +
  scale_colour_manual(values = colours, breaks = c("High shore", "Low shore"), guide = "none") +
  scale_y_continuous(breaks = seq(0, 25, 5), expand = expansion(mult = 0)) + coord_cartesian(ylim = c(0, 25), clip = "off") +
  labs(y = "Mean temperature, 00:00\u201305:00 (\u00b0C)") + paper_theme + theme(legend.position = "none")

pC <- ggplot(counts, aes(shore, days, fill = shore)) +
  geom_col(width = 0.62) +
  geom_text(aes(y = days + 3, label = label), size = 2.5, lineheight = 0.95, vjust = 0, family = "Arial") +
  facet_wrap(~site, nrow = 1) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_x_discrete(labels = c("Low shore" = "Low", "High shore" = "High"), expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, 92), breaks = c(0, 20, 40, 60, 80, 92), expand = expansion(mult = 0)) +
  labs(y = "Days with \u22653 h above 31 \u00b0C") +
  paper_theme + theme(legend.position = "none", plot.margin = margin(2, 1, 2, 2, unit = "mm"), panel.spacing.x = grid::unit(3, "mm"))

figure <- (pA / (pB | pC)) + plot_layout(heights = c(1.1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(1, 1, 1, 1, unit = "mm"))) &
  theme(legend.position = "top", legend.justification = "center")
figure
ggsave("FigureS1.pdf", figure, width = 180, height = 150, units = "mm", device = cairo_pdf, bg = "white")
