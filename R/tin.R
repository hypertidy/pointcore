#' A TIN from the point-cloud form, and the returns in each triangle
#'
#' `points_tin()` triangulates the selected returns (ground, class 2, by
#' default) with [meshcore::tri_mesh()], carrying `z` as an attribute, and
#' adds `input_row` to the vertex table: the row of `x$points` each vertex
#' came from. Input vertices come first in the triangle form (spec Decision
#' 8), so the map costs nothing. Returns sharing an (x, y) merge into one
#' vertex, which keeps the first return's row.
#'
#' `face_points()` locates returns in the TIN's triangles (in x and y) and
#' returns the link table `face`, `point`, `u`, `v`, `dz`: the triangle, the
#' row of `x$points`, the barycentric weights of the triangle's second and
#' third vertex, and the height of the return above the TIN surface. The
#' `above` and `below` thresholds keep the table to returns near the surface
#' (for example `above = 2` for ground cover on a vegetated site); returns
#' outside the TIN are dropped.
#'
#' @param x A point-cloud form from [read_points()].
#' @param classes Classifications to triangulate, or `NULL` for every return.
#' @param engine Passed to [meshcore::tri_mesh()].
#' @param ... Passed to [meshcore::tri_mesh()] (for example `max_area`).
#' @param tin A TIN from `points_tin()`.
#' @param above,below Keep returns at most this far above or below the
#'   surface, in z units.
#' @returns `points_tin()`: a `meshcore_tri` triangle form with `z` and
#'   `input_row` on the vertices. `face_points()`: a data frame.
#' @export
#' @examples
#' if (requireNamespace("cdtr", quietly = TRUE) && requireNamespace("geometry", quietly = TRUE)) {
#'   f <- system.file("extdata", "example.laz", package = "rlas")
#'   p <- read_points(f)
#'   tin <- points_tin(p, classes = NULL)
#'   fp <- face_points(tin, p)
#'   head(fp)
#' }
points_tin <- function(x, classes = 2L, engine = "cdtr", ...) {
  p <- x$points
  idx <- if (is.null(classes)) seq_len(nrow(p)) else which(p$classification %in% classes)
  if (length(idx) < 3L) stop("fewer than 3 returns selected")
  m <- meshcore::tri_mesh(p$x[idx], p$y[idx], PA = cbind(z = p$z[idx]),
                          engine = engine, crs = attr(x, "meshspec")$crs, ...)
  vr <- attr(m, "input_map")
  first <- !duplicated(vr)
  input_row <- rep(NA_integer_, nrow(m$vertices))
  input_row[vr[first]] <- idx[first]
  m$vertices$input_row <- input_row
  attr(m, "meshspec")$points_crs <- attr(x, "meshspec")$crs
  m
}

#' @rdname points_tin
#' @export
face_points <- function(tin, x, above = Inf, below = Inf) {
  p <- x$points
  loc <- meshcore::locate(tin, p$x, p$y)
  t <- tin$triangles
  z <- tin$vertices$z
  f <- loc$face
  ok <- !is.na(f)
  zs <- rep(NA_real_, nrow(p))
  zs[ok] <- (1 - loc$u[ok] - loc$v[ok]) * z[t$v0[f[ok]]] +
    loc$u[ok] * z[t$v1[f[ok]]] + loc$v[ok] * z[t$v2[f[ok]]]
  dz <- p$z - zs
  keep <- ok & dz <= above & dz >= -below
  out <- data.frame(face = f[keep], point = which(keep),
                    u = loc$u[keep], v = loc$v[keep], dz = dz[keep])
  o <- order(out$face, out$point)
  out <- out[o, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "thresholds") <- c(above = above, below = below)
  out
}

#' Per-face summaries of the returns in a TIN
#'
#' `face_stats()` turns the `face_points()` link table into one row per
#' triangle, so questions such as "which faces hold returns more than 2 m
#' above the ground" are a filter on a table, with no render. The face's own
#' shape comes first: plan `area`, `slope` and `aspect` (degrees, aspect
#' clockwise from grid north) and the height `z` at its centroid. Then the
#' returns linked to it: `n`, the largest, mean and standard deviation of
#' their height above the surface (`dz_max`, `dz_mean`, `dz_sd`), a count
#' per classification (`n_class_<code>`) and mean `intensity` when the form
#' has it. Faces with no returns get `n = 0` and `NA` summaries.
#'
#' Area and slope are in the TIN's CRS units, so they are only true in a
#' conformal projection near its centre (for example UTM); transform the
#' returns first with [points_transform()] when the source is Web Mercator.
#'
#' @param tin A TIN from `points_tin()`.
#' @param x The point-cloud form the TIN and links came from.
#' @param links A link table from `face_points()`.
#' @returns A data frame with one row per row of `tin$triangles`.
#' @export
face_stats <- function(tin, x, links) {
  t <- tin$triangles
  v <- tin$vertices
  nf <- nrow(t)
  x0 <- v$x[t$v0]; y0 <- v$y[t$v0]; z0 <- v$z[t$v0]
  e1 <- cbind(v$x[t$v1] - x0, v$y[t$v1] - y0, v$z[t$v1] - z0)
  e2 <- cbind(v$x[t$v2] - x0, v$y[t$v2] - y0, v$z[t$v2] - z0)
  nx <- e1[, 2] * e2[, 3] - e1[, 3] * e2[, 2]
  ny <- e1[, 3] * e2[, 1] - e1[, 1] * e2[, 3]
  nz <- e1[, 1] * e2[, 2] - e1[, 2] * e2[, 1]
  out <- data.frame(face = seq_len(nf),
                    area = abs(nz) / 2,
                    slope = atan2(sqrt(nx^2 + ny^2), abs(nz)) * 180 / pi,
                    ## the downhill direction is the normal's horizontal part
                    aspect = (atan2(nx * sign(nz), ny * sign(nz)) * 180 / pi) %% 360,
                    z = (z0 + v$z[t$v1] + v$z[t$v2]) / 3)
  f <- factor(links$face, levels = seq_len(nf))
  out$n <- as.vector(table(f))
  agg <- function(val, fun) {
    r <- rep(NA_real_, nf)
    s <- tapply(val, f, fun)
    r[as.integer(names(s))] <- s
    r[out$n == 0L] <- NA_real_
    r
  }
  out$dz_max <- agg(links$dz, max)
  out$dz_mean <- agg(links$dz, mean)
  out$dz_sd <- agg(links$dz, stats::sd)
  p <- x$points
  if (!is.null(p$classification)) {
    cl <- p$classification[links$point]
    for (k in sort(unique(cl))) {
      out[[paste0("n_class_", k)]] <- as.vector(table(f[cl == k]))
    }
  }
  if (!is.null(p$intensity)) out$intensity <- agg(as.double(p$intensity[links$point]), mean)
  out
}
