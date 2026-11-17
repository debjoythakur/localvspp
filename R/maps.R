#' Representative event in each quadrant of the study area
#'
#' Splits events at the median longitude/latitude and returns, for each
#' quadrant, the event closest to the quadrant's median location.
#'
#' @param data Data frame of events.
#' @param labels Names for the NW, NE, SW and SE quadrants.
#' @return A data frame with one row per quadrant and a `city_part` column.
#' @export
representative_locations <- function(data, labels = c("Nearby University City", "Northeast St. Louis",
                                                      "South St. Louis", "East St. Louis")) {
  data <- standardize_coords(data)
  east  <- data$longitude > stats::median(data$longitude)
  north <- data$latitude >= stats::median(data$latitude)
  data$city_part <- ifelse(north, ifelse(east, labels[2], labels[1]), ifelse(east, labels[4], labels[3]))
  out <- lapply(split(data, data$city_part), function(d) {
    d[which.min((d$longitude - stats::median(d$longitude))^2 + (d$latitude - stats::median(d$latitude))^2), ]
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Interactive map of one location with its locally selected predictors
#'
#' Requires the \pkg{leaflet} and \pkg{htmltools} packages.
#'
#' @param object A `localvs` fit with `LLI` and/or `LSI`.
#' @param row One row of [representative_locations()] (or any data-frame row
#'   with `longitude`, `latitude` and optionally `city_part`, `Road_Name`,
#'   `Type`, `Occurance_Time`).
#' @param boundary `sf` polygon drawn on the map.
#' @param large Larger layout with predictor names in two columns.
#' @return An `htmltools` tag (print it or wrap it in [htmltools::browsable()]).
#' @export
location_panel <- function(object, row, boundary, large = FALSE) {
  for (p in c("leaflet", "htmltools"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.", call. = FALSE)
  tags <- htmltools::tags
  fld  <- function(v) if (v %in% names(row)) as.character(row[[v]]) else "NA"
  s <- if (large) list(map = c("600px", "520px"), plus = 44, wt = 4,
                       wrap = "width:1200px;min-height:550px;gap:22px;padding:15px;background-color:white;",
                       box = "width:560px;min-height:520px;padding:18px;border:2px solid #333;border-radius:8px;",
                       h = "margin:0 0 12px 0;font-size:27px;", info = "font-size:19px;line-height:1.35;margin-bottom:10px;",
                       hr = "margin:12px 0;border-top:1px solid #888;", head = "font-size:21px;margin-bottom:8px;")
       else      list(map = c("330px", "360px"), plus = 30, wt = 3,
                       wrap = "width:650px;height:370px;gap:8px;margin:4px;",
                       box = "width:310px;height:360px;padding:8px;border:1.5px solid #333;border-radius:6px;font-size:12px;overflow:hidden;",
                       h = "margin:0 0 6px 0;font-size:15px;", info = "margin:0;line-height:1.15;",
                       hr = "margin:6px 0;", head = "font-size:11px;margin-bottom:3px;")
  pretty_name <- function(x) tools::toTitleCase(gsub("([0-9]+)m$", "(\\1m)", gsub("_", " ", x)))
  preds <- function(v) {
    if (!length(v)) v <- "No selected predictors"
    if (!large) return(tags$div(style = "font-size:9px;line-height:1.25;text-align:justify;word-wrap:break-word;",
                                paste(v, collapse = ", ")))
    if (identical(v, "No selected predictors")) return(tags$div(style = "font-size:18px;line-height:1.35;", v))
    tags$div(style = "display:grid;grid-template-columns:1fr 1fr;column-gap:24px;row-gap:6px;font-size:17px;line-height:1.25;",
             lapply(pretty_name(unique(v)), tags$div))
  }
  section <- function(m, label, col) {
    if (is.null(object$selection[[m]])) return(NULL)
    htmltools::tagList(tags$hr(style = s$hr),
                       tags$div(style = paste0("color:", col, ";font-weight:bold;", s$head), label),
                       preds(selected_at(object, row$longitude, row$latitude, m)))
  }
  lon <- round(row$longitude, 5); lat <- round(row$latitude, 5)
  coords <- if (large) list(tags$b("Longitude: "), lon, tags$br(), tags$b("Latitude: "), lat)
            else list(tags$b("Lon: "), lon, "   ", tags$b("Lat: "), lat)
  box <- tags$div(style = paste0("background-color:#f8f9fa;font-family:Arial;box-sizing:border-box;", s$box),
    tags$h3(fld("city_part"), style = paste0("color:#2c3e50;font-weight:bold;", s$h)),
    tags$div(style = s$info, tags$b("Road: "), fld("Road_Name"), tags$br(), tags$b("Crime: "), fld("Type"),
             tags$br(), tags$b("Time: "), fld("Occurance_Time"), tags$br(), coords),
    section("LLI", "Local LASSO selected predictors", "#008000"),
    section("LSI", "Local SCAD selected predictors", "#0055cc"))

  b  <- sf::st_transform(boundary, 4326); bb <- sf::st_bbox(b)
  m <- leaflet::leaflet(width = s$map[1], height = s$map[2])
  m <- leaflet::addProviderTiles(m, leaflet::providers$CartoDB.Positron)
  m <- leaflet::addPolygons(m, data = b, color = "black", weight = s$wt, opacity = 1,
                            fill = TRUE, fillColor = "transparent", fillOpacity = 0)
  m <- leaflet::addLabelOnlyMarkers(m, lng = row$longitude, lat = row$latitude, label = "+",
         labelOptions = leaflet::labelOptions(noHide = TRUE, direction = "center", textOnly = TRUE,
           style = list(color = "black", `font-size` = paste0(s$plus, "px"), `font-weight` = "bold")))
  m <- leaflet::fitBounds(m, bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
  tags$div(style = paste0("display:flex;align-items:flex-start;box-sizing:border-box;", s$wrap), m, box)
}

#' Map of one event with its observed features (Figure 1)
#'
#' Shows a single event on the city map together with the values of a few
#' built-environment features at that event. Two panels side by side give
#' Figure 1 of the paper. Requires the \pkg{leaflet} and \pkg{htmltools}
#' packages.
#'
#' @param data Data frame of events, e.g. [stl_crime].
#' @param i Row of `data` to show.
#' @param boundary `sf` polygon drawn on the map.
#' @param features Feature columns listed under the map.
#' @param marker_color Marker colour (an `awesomeIcons` colour name).
#' @param tag Panel tag, e.g. `"(a)"`.
#' @return An `htmltools` tag.
#' @examples
#' \dontrun{
#' htmltools::browsable(htmltools::tags$div(
#'   style = "display:grid;grid-template-columns:1fr 1fr;gap:20px;",
#'   occurrence_panel(stl_crime, 1, stl_boundary, marker_color = "blue", tag = "(a)"),
#'   occurrence_panel(stl_crime, 10, stl_boundary, marker_color = "purple", tag = "(b)")))
#' }
#' @export
occurrence_panel <- function(data, i, boundary,
                             features = c("primary_road_500m", "residential_road_500m", "service_road_500m",
                                          "buildings_50m", "parking_500m"),
                             marker_color = "blue", tag = "(a)") {
  for (p in c("leaflet", "htmltools"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.", call. = FALSE)
  tags <- htmltools::tags
  data <- standardize_coords(data)
  if (length(i) != 1 || i < 1 || i > nrow(data)) stop("i must be a single row of data.", call. = FALSE)
  pt <- data[i, , drop = FALSE]
  fld <- function(v, alt = "NA") if (v %in% names(pt)) as.character(pt[[v]]) else alt
  features <- intersect(features, names(pt))
  desc  <- paste0(fld("Type", "Unknown event"), " at ", fld("Road_Name", "Unknown road"))
  popup <- paste0("<b>Type:</b> ", fld("Type"), "<br><b>Road:</b> ", fld("Road_Name"),
                  "<br><b>Time:</b> ", fld("Occurance_Time"))

  b  <- sf::st_transform(boundary, 4326); bb <- sf::st_bbox(b)
  m <- leaflet::leaflet(width = "100%", height = 460)
  m <- leaflet::addProviderTiles(m, leaflet::providers$CartoDB.Positron)
  m <- leaflet::addPolygons(m, data = b, color = "black", weight = 3, opacity = 1,
                            fill = TRUE, fillColor = "transparent", fillOpacity = 0)
  m <- leaflet::addAwesomeMarkers(m, pt$longitude, pt$latitude, popup = popup, label = desc,
         icon = leaflet::makeAwesomeIcon(icon = "star", markerColor = marker_color,
                                         iconColor = "white", library = "fa"))
  m <- leaflet::fitBounds(m, bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])

  feat <- paste0("<p><b>Observed features:</b></p><ul style='font-size:14px;'>",
                 paste0("<li><b>", features, "</b>: ", unlist(pt[1, features]), "</li>", collapse = ""),
                 "</ul>")
  tags$div(style = "border:1px solid #bdbdbd;border-radius:6px;padding:8px;",
           tags$h3(paste0(tag, " Occurrence ", i, ": ", desc),
                   style = "text-align:center;font-family:Arial;font-size:16px;"),
           m,
           tags$div(htmltools::HTML(feat),
                    style = "background:#f7f7f7;padding:10px;margin-top:8px;border-radius:6px;border:1px solid #ddd;"))
}
