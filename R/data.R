#' St. Louis crime events with built-environment features
#'
#' Crime events reported by the St. Louis Police between 1 March and
#' 30 April 2026 (mostly assaults, harassment and peace disturbances), one
#' row per unique location, with high-dimensional built-environment
#' covariates: counts of 31 OpenStreetMap feature types within 50, 100, 200,
#' 500 and 1000 m of each event. Built by `data-raw/stl_crime.R` from 652
#' records (247 repeated locations removed).
#'
#' @format A data frame with 405 rows and 163 columns:
#' * `Type`: crime type.
#' * `Reference`: police report number.
#' * `Road_Name`: block-level address.
#' * `Source`: reporting agency.
#' * `Occurance_Time`: date and time of the event (character).
#' * `latitude`, `longitude`: event location (WGS84).
#' * `crime_id`: row identifier in the raw file.
#' * 155 built-environment counts named `<feature>_<radius>m`, with radius
#'   50, 100, 200, 500 or 1000 m, for the features `atm`, `bank`, `bars`,
#'   `buildings`, `bus_stops`, `cafes`, `cinema`, `commercial`, `convenience`,
#'   `fast_food`, `fuel`, `hospital`, `industrial`, `intersection_nodes`,
#'   `motorway`, `nightclubs`, `parking`, `parks`, `pharmacy`, `playgrounds`,
#'   `primary_road`, `pubs`, `residential_road`, `restaurants`, `retail`,
#'   `schools`, `secondary_road`, `service_road`, `supermarket`,
#'   `tertiary_road`, `universities`.
#'
#' [crime_covariates()] returns the 140 non-constant covariates used as
#' candidate predictors.
#' @source St. Louis Metropolitan Police Department crime records; features
#'   from OpenStreetMap.
#' @examples
#' dim(stl_crime)
#' head(crime_covariates(stl_crime)[, 1:5])
"stl_crime"

#' St. Louis City boundary
#'
#' @format An `sf` object with one polygon (city limit, 2002; about 62 sq mi),
#'   transformed from NAD83 Missouri East State Plane (US ft) to EPSG:4326.
#' @source City of St. Louis city-limit boundary (`CITY_LIMIT_BOUND_2002`);
#'   built by `data-raw/stl_crime.R`.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) plot(sf::st_geometry(stl_boundary))
"stl_boundary"