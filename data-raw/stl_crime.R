## Builds the package datasets `stl_crime` and `stl_boundary`.
## Run once from the package root:  source("data-raw/stl_crime.R")
## Inputs (kept in data-raw/, excluded from the built package):
##   data-raw/stl_high_dimensional_crime_features.csv
##   data-raw/stl_boundary.shp (+ .shx/.dbf/.prj), NAD83 Missouri East State Plane (US ft)
##
## Raw CSV: 652 St. Louis Police crime records (1 Mar - 30 Apr 2026), 168 columns.
## Result:  405 events at unique locations; 31 OSM features x 5 radii = 155 covariates.

library(sf)

## ---- Crime events -------------------------------------------------------
stl_crime <- read.csv("data-raw/stl_high_dimensional_crime_features.csv", check.names = FALSE)
names(stl_crime) <- gsub("\\s+", "_", trimws(names(stl_crime)))
stl_crime <- stl_crime[, setdiff(names(stl_crime),                              # temporal columns (164:168),
                                  c("hour", "weekday", "month", "weekend", "night_indicator"))]  # not used
stl_crime <- stl_crime[!duplicated(stl_crime[c("Longitude", "Latitude")]), ]   # one event per location
names(stl_crime)[names(stl_crime) == "Longitude"] <- "longitude"
names(stl_crime)[names(stl_crime) == "Latitude"]  <- "latitude"
stl_crime <- stl_crime[is.finite(stl_crime$longitude) & is.finite(stl_crime$latitude), ]
rownames(stl_crime) <- NULL

## ---- City boundary ------------------------------------------------------
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
stl_boundary <- st_read("data-raw/stl_boundary.shp", quiet = TRUE)
if (is.na(st_crs(stl_boundary))) st_crs(stl_boundary) <- 4326
stl_boundary <- st_transform(st_make_valid(stl_boundary), 4326)

usethis::use_data(stl_crime, stl_boundary, overwrite = TRUE, compress = "xz")
