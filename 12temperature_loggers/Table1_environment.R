setwd("C:/path/to/RNAgithub/12temperature_loggers")

# =============================================================================
# Environmental characteristics of the three sampling sites (Table S1 and Supplementary Text 1):
# air temperature (CHELSA), sea surface temperature (Bio-ORACLE) and tidal range (nearest tide gauge).
# Input: downloaded directly from the data providers (internet connection needed)
# Output: TableS1.txt
# =============================================================================

library(terra)
library(oce)
library(xml2)
library(dplyr)

# Allow enough time for the downloads.
options(timeout = 600)

# Sampling coordinates
lat_north <- 67.19264; lon_north <- 14.35769
lat_central <- 58.83140; lon_central <- 10.99100
lat_south <- 42.26244; lon_south <- -8.75789


#### 1. Air temperature: CHELSA v2.1, 1981-2010 ####

# BIO1 (mean annual air temperature) and BIO5 (mean daily maximum air temperature of the warmest month), read from the CHELSA server
Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR")
bio1 <- terra::rast("/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio01/1981-2010/CHELSA_bio01_1981-2010_V.2.1.tif")
bio5 <- terra::rast("/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio05/1981-2010/CHELSA_bio05_1981-2010_V.2.1.tif")
names(bio1) <- "bio1"; names(bio5) <- "bio5"

# Extract air temperature from the grid cell containing each site, without interpolation
air_north <- terra::extract(c(bio1, bio5), cbind(lon_north, lat_north))
air_central <- terra::extract(c(bio1, bio5), cbind(lon_central, lat_central))
air_south <- terra::extract(c(bio1, bio5), cbind(lon_south, lat_south))

# Record the centre of each selected cell and its distance to the site (km)
cell_air_north <- terra::xyFromCell(bio1, terra::cellFromXY(bio1, cbind(lon_north, lat_north)))
cell_air_central <- terra::xyFromCell(bio1, terra::cellFromXY(bio1, cbind(lon_central, lat_central)))
cell_air_south <- terra::xyFromCell(bio1, terra::cellFromXY(bio1, cbind(lon_south, lat_south)))
air_north$latitude <- cell_air_north[, "y"]; air_north$longitude <- cell_air_north[, "x"]
air_central$latitude <- cell_air_central[, "y"]; air_central$longitude <- cell_air_central[, "x"]
air_south$latitude <- cell_air_south[, "y"]; air_south$longitude <- cell_air_south[, "x"]
air_north$km <- 111.2 * sqrt((air_north$latitude - lat_north)^2 + ((air_north$longitude - lon_north) * cos(lat_north * pi / 180))^2)
air_central$km <- 111.2 * sqrt((air_central$latitude - lat_central)^2 + ((air_central$longitude - lon_central) * cos(lat_central * pi / 180))^2)
air_south$km <- 111.2 * sqrt((air_south$latitude - lat_south)^2 + ((air_south$longitude - lon_south) * cos(lat_south * pi / 180))^2)

# Check air temperature, cell centres and distances
rbind(air_north, air_central, air_south)


#### 2. Sea surface temperature: Bio-ORACLE v3.0, 2010-2019 ####

# Download the decadal mean and the mean warmest-month surface temperature in a small box around each site (time index 1 = 2010-2019)
url_sea <- "https://erddap.bio-oracle.org/erddap/griddap/thetao_baseline_2000_2019_depthsurf.csv?"
box_north <- sprintf("[1:1:1][(%.3f):1:(%.3f)][(%.3f):1:(%.3f)]", lat_north - 0.15, lat_north + 0.15, lon_north - 0.15, lon_north + 0.15)
box_central <- sprintf("[1:1:1][(%.3f):1:(%.3f)][(%.3f):1:(%.3f)]", lat_central - 0.15, lat_central + 0.15, lon_central - 0.15, lon_central + 0.15)
box_south <- sprintf("[1:1:1][(%.3f):1:(%.3f)][(%.3f):1:(%.3f)]", lat_south - 0.15, lat_south + 0.15, lon_south - 0.15, lon_south + 0.15)
sea_north <- read.csv(paste0(url_sea, URLencode(paste0("thetao_mean", box_north, ",thetao_ltmax", box_north), reserved = TRUE)))[-1, ]
sea_central <- read.csv(paste0(url_sea, URLencode(paste0("thetao_mean", box_central, ",thetao_ltmax", box_central), reserved = TRUE)))[-1, ]
sea_south <- read.csv(paste0(url_sea, URLencode(paste0("thetao_mean", box_south, ",thetao_ltmax", box_south), reserved = TRUE)))[-1, ]
sea_north <- dplyr::mutate(sea_north, dplyr::across(c(latitude, longitude, thetao_mean, thetao_ltmax), as.numeric))
sea_central <- dplyr::mutate(sea_central, dplyr::across(c(latitude, longitude, thetao_mean, thetao_ltmax), as.numeric))
sea_south <- dplyr::mutate(sea_south, dplyr::across(c(latitude, longitude, thetao_mean, thetao_ltmax), as.numeric))

# Keep marine cells with valid values for both variables
sea_north <- dplyr::filter(sea_north, !is.na(thetao_mean), !is.na(thetao_ltmax))
sea_central <- dplyr::filter(sea_central, !is.na(thetao_mean), !is.na(thetao_ltmax))
sea_south <- dplyr::filter(sea_south, !is.na(thetao_mean), !is.na(thetao_ltmax))

# Select the marine cell nearest to each site (distance in km)
sea_north <- dplyr::mutate(sea_north, km = 111.2 * sqrt((latitude - lat_north)^2 + ((longitude - lon_north) * cos(lat_north * pi / 180))^2))
sea_central <- dplyr::mutate(sea_central, km = 111.2 * sqrt((latitude - lat_central)^2 + ((longitude - lon_central) * cos(lat_central * pi / 180))^2))
sea_south <- dplyr::mutate(sea_south, km = 111.2 * sqrt((latitude - lat_south)^2 + ((longitude - lon_south) * cos(lat_south * pi / 180))^2))
sea_north <- dplyr::slice_min(sea_north, km, n = 1, with_ties = FALSE)
sea_central <- dplyr::slice_min(sea_central, km, n = 1, with_ties = FALSE)
sea_south <- dplyr::slice_min(sea_south, km, n = 1, with_ties = FALSE)

# Check the decade (time should start with 2010), cell centres and distances
rbind(sea_north, sea_central, sea_south)


#### 3. Tidal range: nearest tide gauge, 2023 ####

# Hourly time steps for aligning observations, and 10-min time steps for locating high and low water
hours_2023 <- seq(as.POSIXct("2023-01-01 00:00", tz = "UTC"), as.POSIXct("2023-12-31 23:00", tz = "UTC"), by = "hour")
times_2023 <- seq(as.POSIXct("2023-01-01 00:00", tz = "UTC"), as.POSIXct("2023-12-31 23:50", tz = "UTC"), by = "10 min")

# North: Bodø gauge (Norwegian Mapping Authority), water level in cm
xml_bodo <- xml2::read_xml("https://vannstand.kartverket.no/tideapi.php?tide_request=stationdata&stationcode=BOO&refcode=msl&fromtime=2023-01-01T00:00&totime=2023-12-31T23:00&datatype=obs&interval=60&dst=0&tzone=0&lang=en")
nodes_bodo <- xml2::xml_find_all(xml_bodo, ".//data[@type='observation']/waterlevel")
tide_north <- data.frame(time = xml2::xml_attr(nodes_bodo, "time"), height = as.numeric(xml2::xml_attr(nodes_bodo, "value")) / 100, flag = xml2::xml_attr(nodes_bodo, "flag"))
lat_bodo <- as.numeric(xml2::xml_attr(xml2::xml_find_first(xml_bodo, ".//location"), "latitude"))
lon_bodo <- as.numeric(xml2::xml_attr(xml2::xml_find_first(xml_bodo, ".//location"), "longitude"))

# Central: Kungsvik gauge (Swedish Meteorological and Hydrological Institute), water level in cm
tide_central <- read.csv2("https://opendata-download-ocobs.smhi.se/api/version/1.0/parameter/6/station/2130/period/corrected-archive/data.csv", skip = 6, dec = ".", fill = TRUE)
tide_central <- data.frame(time = tide_central[[1]], height = as.numeric(tide_central[[2]]) / 100, flag = tide_central[[3]])
lat_kungsvik <- 58.9967; lon_kungsvik <- 11.1272

# South: Vigo gauge (GESLA, through the University of Hawaii Sea Level Center), water level in m
tide_south <- read.csv("https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_gesla.csv?time,sea_level,flag1,flag2&record_id=%22vigotg-vig-esp-cmems%22&time%3E=2023-01-01T00:00:00Z&time%3C2024-01-01T00:00:00Z", skip = 2, header = FALSE)
tide_south <- data.frame(time = tide_south[[1]], height = tide_south[[2]], flag1 = tide_south[[3]], flag2 = tide_south[[4]])
lat_vigo <- 42.243; lon_vigo <- -8.726

# Keep good-quality observations
tide_north <- dplyr::filter(tide_north, flag == "obs")
tide_central <- dplyr::filter(tide_central, flag == "G")
tide_south <- dplyr::filter(tide_south, flag1 == 1, flag2 == 1)

# Align observations to the 2023 hourly time steps
tide_north$time <- as.POSIXct(substr(sub("T", " ", tide_north$time), 1, 19), tz = "UTC")
tide_central$time <- as.POSIXct(substr(sub("T", " ", tide_central$time), 1, 19), tz = "UTC")
tide_south$time <- as.POSIXct(substr(sub("T", " ", tide_south$time), 1, 19), tz = "UTC")
height_north <- tide_north$height[match(hours_2023, tide_north$time)]
height_central <- tide_central$height[match(hours_2023, tide_central$time)]
height_south <- tide_south$height[match(hours_2023, tide_south$time)]

# Check the proportion of 2023 hours with data (should be above 0.9)
c(north = mean(!is.na(height_north)), central = mean(!is.na(height_central)), south = mean(!is.na(height_south)))

# Fit a harmonic tidal model to each gauge and predict the astronomical tide
fit_north <- oce::tidem(oce::as.sealevel(height_north, hours_2023, latitude = lat_bodo, longitude = lon_bodo, deltat = 1))
fit_central <- oce::tidem(oce::as.sealevel(height_central, hours_2023, latitude = lat_kungsvik, longitude = lon_kungsvik, deltat = 1))
fit_south <- oce::tidem(oce::as.sealevel(height_south, hours_2023, latitude = lat_vigo, longitude = lon_vigo, deltat = 1))
pred_north <- as.numeric(predict(fit_north, newdata = times_2023))
pred_central <- as.numeric(predict(fit_central, newdata = times_2023))
pred_south <- as.numeric(predict(fit_south, newdata = times_2023))

# Find predicted high and low waters
turn_north <- diff(sign(diff(pred_north)))
turn_central <- diff(sign(diff(pred_central)))
turn_south <- diff(sign(diff(pred_south)))
high_north <- which(turn_north == -2) + 1; low_north <- which(turn_north == 2) + 1
high_central <- which(turn_central == -2) + 1; low_central <- which(turn_central == 2) + 1
high_south <- which(turn_south == -2) + 1; low_south <- which(turn_south == 2) + 1

# Tidal range = each high water minus the following low water, keeping only complete pairs within 2023
range_north <- pred_north[high_north] - pred_north[low_north[findInterval(high_north, low_north) + 1]]
range_central <- pred_central[high_central] - pred_central[low_central[findInterval(high_central, low_central) + 1]]
range_south <- pred_south[high_south] - pred_south[low_south[findInterval(high_south, low_south) + 1]]
range_north <- range_north[!is.na(range_north)]
range_central <- range_central[!is.na(range_central)]
range_south <- range_south[!is.na(range_south)]

# Check the number of tidal cycles per gauge
c(north = length(range_north), central = length(range_central), south = length(range_south))

# Distance (km) from each site to its tide gauge
km_bodo <- 111.2 * sqrt((lat_bodo - lat_north)^2 + ((lon_bodo - lon_north) * cos(lat_north * pi / 180))^2)
km_kungsvik <- 111.2 * sqrt((lat_kungsvik - lat_central)^2 + ((lon_kungsvik - lon_central) * cos(lat_central * pi / 180))^2)
km_vigo <- 111.2 * sqrt((lat_vigo - lat_south)^2 + ((lon_vigo - lon_south) * cos(lat_south * pi / 180))^2)
c(north = km_bodo, central = km_kungsvik, south = km_vigo)


#### 4. Table S1 and the details reported in Supplementary Text 1 ####

# Table S1
env_table <- data.frame(
  Population = c("North", "Central", "South"),
  Site = c(sprintf("Ersvika, Norway (%.5f° N, %.5f° E)", lat_north, lon_north), sprintf("Ursholmen, Sweden (%.5f° N, %.5f° E)", lat_central, lon_central), sprintf("Canabal, Spain (%.5f° N, %.5f° W)", lat_south, abs(lon_south))),
  air_mean = sprintf("%.1f", c(air_north$bio1, air_central$bio1, air_south$bio1)),
  air_max = sprintf("%.1f", c(air_north$bio5, air_central$bio5, air_south$bio5)),
  sea_mean = sprintf("%.1f", c(sea_north$thetao_mean, sea_central$thetao_mean, sea_south$thetao_mean)),
  sea_warmest = sprintf("%.1f", c(sea_north$thetao_ltmax, sea_central$thetao_ltmax, sea_south$thetao_ltmax)),
  tide_mean = sprintf("%.2f", c(mean(range_north), mean(range_central), mean(range_south))),
  tide_max = sprintf("%.2f", c(max(range_north), max(range_central), max(range_south)))
)
names(env_table) <- c("Population", "Site",
                      "Air temperature (°C): annual mean", "Air temperature (°C): warmest-month daily maximum",
                      "Sea surface temperature (°C): annual mean", "Sea surface temperature (°C): warmest month",
                      "Tidal range (m): mean", "Tidal range (m): maximum")
env_table
write.table(env_table, "TableS1.txt", sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

# Grid cells, tide gauges and distances reported in Supplementary Text 1
env_details <- data.frame(
  Population = c("North", "Central", "South"),
  air_cell_latitude = sprintf("%.5f", c(air_north$latitude, air_central$latitude, air_south$latitude)),
  air_cell_longitude = sprintf("%.5f", c(air_north$longitude, air_central$longitude, air_south$longitude)),
  air_cell_km = sprintf("%.2f", c(air_north$km, air_central$km, air_south$km)),
  sea_cell_latitude = sprintf("%.5f", c(sea_north$latitude, sea_central$latitude, sea_south$latitude)),
  sea_cell_longitude = sprintf("%.5f", c(sea_north$longitude, sea_central$longitude, sea_south$longitude)),
  sea_cell_km = sprintf("%.2f", c(sea_north$km, sea_central$km, sea_south$km)),
  tide_gauge = c("Bodø", "Kungsvik", "Vigo"),
  tide_gauge_km = sprintf("%.1f", c(km_bodo, km_kungsvik, km_vigo)),
  tide_hours_with_data = sprintf("%.3f", c(mean(!is.na(height_north)), mean(!is.na(height_central)), mean(!is.na(height_south)))),
  tide_cycles = c(length(range_north), length(range_central), length(range_south))
)
env_details

