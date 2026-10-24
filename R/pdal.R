## read_points(engine = "pdal"): one PDAL pipeline per source, the view
## mapped to the point-cloud form

## PDAL dimension names to spec names, in spec order
pdal_names <- c(X = "x", Y = "y", Z = "z", Intensity = "intensity",
                ReturnNumber = "return_number", NumberOfReturns = "number_of_returns",
                Classification = "classification", ScanAngleRank = "scan_angle",
                GpsTime = "gps_time", Red = "red", Green = "green", Blue = "blue")

## the same other attributes as las_extra, under PDAL's names
pdal_extra <- c(ScanDirectionFlag = "scan_direction_flag",
                EdgeOfFlightLine = "edge_of_flight_line",
                Synthetic = "synthetic", KeyPoint = "key_point",
                Withheld = "withheld", Overlap = "overlap",
                UserData = "user_data", PointSourceId = "point_source_id",
                ScanChannel = "scanner_channel", Infrared = "nir")

## the classification flags, logical as rlas gives them
pdal_flags <- c("synthetic", "key_point", "withheld", "overlap")

read_points_pdal <- function(dsn, stages = NULL, bbox = NULL) {
  if (!requireNamespace("pdal", quietly = TRUE)) {
    stop("engine = \"pdal\" needs the pdal package: ",
         "remotes::install_github(\"hypertidy/pdal\")")
  }
  if (!is.null(bbox)) {
    stopifnot(is.numeric(bbox), length(bbox) == 4L, bbox[1] <= bbox[3], bbox[2] <= bbox[4])
  }
  local <- !grepl("^[A-Za-z][A-Za-z0-9+.-]*://", dsn)
  dsn[local] <- normalizePath(dsn[local], mustWork = TRUE)
  stages <- pdal_stage_json(stages)
  pts <- vector("list", length(dsn))
  src <- vector("list", length(dsn))
  for (i in seq_along(dsn)) {
    rd <- pdal_reader(dsn[i], bbox)
    json <- paste0("[", paste(c(rd$json, stages), collapse = ", "), "]")
    res <- pdal::pdal_pipeline(json, ids = TRUE)
    if (length(res$views) != 1L) {
      stop(sprintf("the pipeline for '%s' gave %d point views; read_points() needs one",
                   dsn[i], length(res$views)))
    }
    d <- as.data.frame(res$views[[1L]])
    id <- if (rd$window) rep(NA_real_, nrow(d)) else d$table_id + 1
    d$table_id <- NULL
    meta <- pdal_reader_metadata(res$metadata)
    s <- pdal_source(dsn[i], meta)
    ## PDAL carries Overlap for every format; LAS has it only from format 6
    if (!is.na(s$point_format) && s$point_format < 6L) d$Overlap <- NULL
    pts[[i]] <- pdal_to_points(d, source = i, id = id)
    src[[i]] <- s
  }
  points <- do.call(rbind, pts)
  sources <- do.call(rbind, src)
  rownames(points) <- rownames(sources) <- NULL
  crs <- unique(sources$crs)
  if (length(crs) > 1L) warning("files have different CRS; the form records the first")
  structure(list(points = points, sources = sources),
            meshspec = list(version = "0.1", form = "points", crs = crs[1L],
                            vertical_crs = attr(src[[1L]], "vertical_crs"),
                            engine = "pdal",
                            pipeline = if (length(stages)) stages else NULL),
            class = "pointcore_points")
}

pdal_stage_json <- function(stages) {
  if (is.null(stages)) return(character())
  if (is.list(stages)) {
    stages <- vapply(stages, function(x) {
      as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA))
    }, "")
  }
  as.character(stages)
}

## the reader stage, and whether it reads a window (COPC and EPT read only
## the nodes that meet the bounds, so table rows are not file positions)
pdal_reader <- function(dsn, bbox) {
  driver <- if (grepl("\\.copc\\.laz$", dsn, ignore.case = TRUE)) {
    "readers.copc"
  } else if (grepl("ept\\.json$", dsn)) {
    "readers.ept"
  } else {
    ""
  }
  stage <- list(filename = dsn)
  if (nzchar(driver)) stage$type <- driver
  window <- !is.null(bbox) && nzchar(driver)
  if (window) stage$bounds <- pdal_bounds(bbox)
  json <- as.character(jsonlite::toJSON(stage, auto_unbox = TRUE, digits = NA))
  ## crop exactly in every case: node bounds are wider than the window
  if (!is.null(bbox)) {
    json <- c(json, sprintf('{"type": "filters.crop", "bounds": "%s"}', pdal_bounds(bbox)))
  }
  list(json = json, window = window)
}

pdal_bounds <- function(bbox) {
  sprintf("([%.17g, %.17g], [%.17g, %.17g])", bbox[1], bbox[3], bbox[2], bbox[4])
}

pdal_to_points <- function(d, source, id) {
  known <- intersect(names(pdal_names), names(d))
  out <- d[known]
  names(out) <- pdal_names[known]
  other <- setdiff(names(d), names(pdal_names))
  if (length(other)) {
    extra <- d[other]
    nm <- pdal_extra[other]
    names(extra) <- ifelse(is.na(nm), snake_case(other), nm)
    for (f in intersect(pdal_flags, names(extra))) extra[[f]] <- extra[[f]] != 0
    out <- cbind(out, extra)
  }
  out$source <- rep(as.integer(source), nrow(out))
  out$id <- as.double(id)
  out
}

## the first reader's metadata from the pipeline metadata JSON
pdal_reader_metadata <- function(json) {
  st <- jsonlite::fromJSON(json, simplifyVector = FALSE)$stages
  i <- grep("^readers\\.", names(st))[1L]
  if (is.na(i)) return(list())
  m <- st[[i]]
  ## a stage used more than once is a list of entries
  if (is.null(names(m))) m <- m[[1L]]
  m
}

pdal_source <- function(dsn, m) {
  num <- function(x) if (is.null(x)) NA_real_ else as.double(x)
  crs <- pdal_crs(m$srs)
  out <- data.frame(dsn = dsn,
                    point_format = if (is.null(m$dataformat_id)) NA_integer_ else as.integer(m$dataformat_id),
                    n = num(m$count),
                    scale_x = num(m$scale_x), scale_y = num(m$scale_y), scale_z = num(m$scale_z),
                    offset_x = num(m$offset_x), offset_y = num(m$offset_y),
                    offset_z = num(m$offset_z),
                    crs = crs$crs, stringsAsFactors = FALSE)
  attr(out, "vertical_crs") <- crs$vertical
  out
}

## "AUTH:code" when PDAL's PROJJSON names one, else the horizontal WKT; the
## vertical part only when it has an authority code
pdal_crs <- function(srs) {
  crs <- NA_character_
  vertical <- NA_character_
  if (is.null(srs)) return(list(crs = crs, vertical = vertical))
  code <- function(x) {
    if (is.list(x) && is.list(x$id) && !is.null(x$id$code)) paste0(x$id$authority, ":", x$id$code)
    else NA_character_
  }
  j <- srs$json
  if (is.list(j)) {
    if (identical(j$type, "CompoundCRS") && length(j$components) >= 1L) {
      crs <- code(j$components[[1L]])
      if (length(j$components) >= 2L) vertical <- code(j$components[[2L]])
    } else {
      crs <- code(j)
    }
  }
  if (is.na(crs) && !is.null(srs$horizontal) && nzchar(srs$horizontal)) crs <- srs$horizontal
  list(crs = crs, vertical = vertical)
}
