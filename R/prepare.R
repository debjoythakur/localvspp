#' Standardize coordinate column names
#'
#' Renames the latitude/longitude columns (`latitude`/`lat` and
#' `longitude`/`lon`/`lng`/`long`, any case) to `latitude` and `longitude` and
#' drops rows with non-finite coordinates.
#'
#' @param data A data frame of events.
#' @return `data` with standardized coordinate columns.
#' @export
standardize_coords <- function(data) {
  nm  <- names(data)
  lat <- nm[tolower(nm) %in% c("latitude", "lat")][1]
  lon <- nm[tolower(nm) %in% c("longitude", "lon", "lng", "long")][1]
  if (is.na(lat) || is.na(lon)) stop("Latitude and Longitude columns not found.", call. = FALSE)
  names(data)[match(c(lat, lon), nm)] <- c("latitude", "longitude")
  data[is.finite(data$latitude) & is.finite(data$longitude), , drop = FALSE]
}

#' Which events fall inside a boundary
#'
#' @param data A data frame with `longitude`/`latitude` columns (WGS84).
#' @param boundary An `sf`/`sfc` polygon.
#' @return A logical vector, one element per row of `data`.
#' @export
inside_boundary <- function(data, boundary) {
  data <- standardize_coords(data)
  pts  <- sf::st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)
  b    <- sf::st_transform(sf::st_make_valid(boundary), 4326)
  lengths(sf::st_intersects(pts, sf::st_union(b))) > 0
}

#' Default non-covariate columns
#'
#' Identifier, text and coordinate columns that [crime_covariates()] drops.
#' @return A character vector.
#' @export
default_drop_cols <- function() {
  c("Type", "Reference", "Road_Name", "Source", "Occurance_Time", "latitude",
    "longitude", "Latitude", "Longitude", "crime_id", "geometry")
}

#' Extract the candidate covariates
#'
#' Keeps numeric, non-constant columns with more than `min_finite` finite
#' values and replaces non-finite values by the column median.
#'
#' @param data A data frame of events.
#' @param drop Columns to exclude (see [default_drop_cols()]).
#' @param min_finite Minimum number of finite values a column needs.
#' @return A data frame of covariates.
#' @export
crime_covariates <- function(data, drop = default_drop_cols(), min_finite = 5) {
  X <- as.data.frame(data)[setdiff(names(data), drop)]
  X <- X[vapply(X, is.numeric, logical(1))]
  ok <- vapply(X, function(z) sum(is.finite(z)) > min_finite &&
                 isTRUE(stats::sd(z, na.rm = TRUE) > 0), logical(1))
  X <- X[ok]
  X[] <- lapply(X, function(z) replace(z, !is.finite(z), stats::median(z[is.finite(z)])))
  X
}
