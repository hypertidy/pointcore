#' Read a window of an Entwine Point Tile (EPT) dataset
#'
#' `read_ept()` reads the returns inside a bounding box from an EPT dataset,
#' local or over HTTP, into the same point-cloud form as [read_points()].
#' EPT splits a cloud into an octree of LAZ files, each depth adding detail,
#' so the window is every node file that meets the box at depths up to
#' `max_depth`. Only the node files the window needs are fetched, and they
#' are kept in `cache`.
#'
#' Each node file is a row of `sources`, with its octree `key` and `depth`,
#' so every return leads back to one file and its position (`id`) in that
#' file. The returns are in the dataset's CRS: the public USGS 3DEP datasets
#' are in Web Mercator (EPSG:3857), where x and y are stretched by 1 /
#' cos(latitude) and heights are not; see [points_transform()].
#'
#' @param url The dataset's folder or its `ept.json`, a URL or a path.
#' @param bbox `c(xmin, ymin, xmax, ymax)` in the dataset's CRS.
#' @param max_depth Deepest octree level to read; `Inf` for full density.
#' @param filter A LASlib filter string, as in [read_points()].
#' @param cache Folder for the fetched node files.
#' @returns A `pointcore_points` form; `sources` gains `key` and `depth`, and
#'   the `meshspec` attribute gains `ept` (the dataset's address and bounds).
#' @export
read_ept <- function(url, bbox, max_depth = Inf, filter = "",
                     cache = file.path(tempdir(), "pointcore-ept")) {
  base <- sub("/?ept\\.json$", "", url)
  info <- jsonlite::fromJSON(ept_fetch(paste0(base, "/ept.json"), cache, base))
  if (!identical(info$dataType, "laszip")) stop("only laszip EPT datasets are supported")
  nodes <- ept_nodes(base, info, bbox, max_depth, cache)
  if (!nrow(nodes)) stop("no EPT nodes meet the bounding box")
  files <- vapply(nodes$key, function(k) {
    ept_fetch(sprintf("%s/ept-data/%s.laz", base, k), cache, base)
  }, "")
  pts <- vector("list", length(files))
  src <- vector("list", length(files))
  crs <- ept_crs(info$srs)
  for (i in seq_along(files)) {
    h <- rlas::read.lasheader(files[i])
    d <- as.data.frame(rlas::read.las(files[i], filter = filter))
    id <- if (nzchar(filter)) filtered_ids(files[i], filter, h) else seq_len(nrow(d))
    p <- las_to_points(d, source = i, id = id)
    inside <- p$x >= bbox[1] & p$x <= bbox[3] & p$y >= bbox[2] & p$y <= bbox[4]
    pts[[i]] <- p[inside, , drop = FALSE]
    s <- las_source(sprintf("%s/ept-data/%s.laz", base, nodes$key[i]), h)
    s$crs <- crs$crs
    s$key <- nodes$key[i]
    s$depth <- nodes$depth[i]
    src[[i]] <- s
  }
  points <- do.call(rbind, pts)
  sources <- do.call(rbind, src)
  rownames(points) <- rownames(sources) <- NULL
  structure(list(points = points, sources = sources),
            meshspec = list(version = "0.1", form = "points", crs = crs$crs,
                            vertical_crs = crs$vertical,
                            ept = list(url = base, bounds = info$bounds, bbox = bbox,
                                       max_depth = max_depth)),
            class = "pointcore_points")
}

#' Transform the returns' x and y to another CRS
#'
#' Heights are left as they are, and so are `source` and `id`, so every
#' return still leads back to its file. The form's previous CRS is kept as
#' `from_crs`.
#'
#' @param x A point-cloud form.
#' @param crs The target CRS, as `sf::sf_project()` reads it.
#' @returns The form with new `x`, `y` and CRS.
#' @export
points_transform <- function(x, crs) {
  m <- attr(x, "meshspec")
  xy <- sf::sf_project(m$crs, crs, cbind(x$points$x, x$points$y), keep = TRUE)
  x$points$x <- xy[, 1]
  x$points$y <- xy[, 2]
  m$from_crs <- m$crs
  m$crs <- crs
  attr(x, "meshspec") <- m
  x
}

## the hierarchy, walked from the root: each node with its key, depth and
## count, kept when its cube meets the box (in x and y)
ept_nodes <- function(base, info, bbox, max_depth, cache) {
  b <- info$bounds
  out <- list()
  walk <- function(page) {
    h <- jsonlite::fromJSON(ept_fetch(sprintf("%s/ept-hierarchy/%s.json", base, page), cache, base))
    for (k in names(h)) {
      n <- h[[k]]
      dxyz <- as.integer(strsplit(k, "-")[[1]])
      d <- dxyz[1]
      if (d > max_depth) next
      if (!ept_meets(b, dxyz, bbox)) next
      if (n == -1) {
        walk(k)
      } else if (n > 0) {
        out[[length(out) + 1L]] <<- data.frame(key = k, depth = d, n = n)
      }
    }
  }
  walk("0-0-0-0")
  if (!length(out)) return(data.frame(key = character(), depth = integer(), n = numeric()))
  nodes <- unique(do.call(rbind, out))
  nodes[order(nodes$depth, nodes$key), , drop = FALSE]
}

ept_meets <- function(b, dxyz, bbox) {
  size <- (b[4:5] - b[1:2]) / 2^dxyz[1]
  lo <- b[1:2] + dxyz[2:3] * size
  hi <- lo + size
  lo[1] <= bbox[3] && hi[1] >= bbox[1] && lo[2] <= bbox[4] && hi[2] >= bbox[2]
}

## a local path for a dataset file, fetched once into the cache
ept_fetch <- function(path, cache, base) {
  if (!grepl("^https?://", path)) return(path)
  rel <- substring(path, nchar(base) + 2L)
  dest <- file.path(cache, gsub("[^A-Za-z0-9._-]", "_", base), rel)
  if (!file.exists(dest)) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    ## remote stores drop the odd connection: retry before giving up
    for (try in 1:4) {
      ok <- tryCatch({
        utils::download.file(path, dest, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (ok) break
      unlink(dest)
      Sys.sleep(2^try)
    }
    if (!ok) stop("could not fetch ", path)
  }
  dest
}

ept_crs <- function(srs) {
  crs <- NA_character_
  vertical <- NA_character_
  if (!is.null(srs$horizontal)) {
    crs <- paste0(if (is.null(srs$authority)) "EPSG" else srs$authority, ":", srs$horizontal)
  } else if (!is.null(srs$wkt)) {
    crs <- srs$wkt
  }
  if (!is.null(srs$vertical)) vertical <- paste0("EPSG:", srs$vertical)
  list(crs = crs, vertical = vertical)
}
