## a one-node EPT dataset built from rlas's example file
fake_ept <- function() {
  f <- system.file("extdata", "example.laz", package = "rlas")
  h <- rlas::read.lasheader(f)
  d <- tempfile("ept")
  dir.create(file.path(d, "ept-data"), recursive = TRUE)
  dir.create(file.path(d, "ept-hierarchy"))
  file.copy(f, file.path(d, "ept-data", "0-0-0-0.laz"))
  lo <- c(h[["Min X"]], h[["Min Y"]], h[["Min Z"]])
  span <- max(h[["Max X"]] - lo[1], h[["Max Y"]] - lo[2], h[["Max Z"]] - lo[3]) + 1
  jsonlite::write_json(list(bounds = c(lo, lo + span), dataType = "laszip", hierarchyType = "json",
                            points = h[["Number of point records"]], span = 256L,
                            srs = list(authority = "EPSG", horizontal = "26917")),
                       file.path(d, "ept.json"), auto_unbox = TRUE, digits = NA)
  jsonlite::write_json(list("0-0-0-0" = h[["Number of point records"]]),
                       file.path(d, "ept-hierarchy", "0-0-0-0.json"), auto_unbox = TRUE)
  d
}

test_that("read_ept() reads a window with its node file and ids", {
  d <- fake_ept()
  all <- read_points(system.file("extdata", "example.laz", package = "rlas"))$points
  bb <- c(quantile(all$x, 0.25), quantile(all$y, 0.25), quantile(all$x, 0.75), quantile(all$y, 0.75))
  p <- read_ept(d, bb)
  expect_s3_class(p, "pointcore_points")
  expect_equal(attr(p, "meshspec")$crs, "EPSG:26917")
  expect_equal(p$sources$key, "0-0-0-0")
  P <- p$points
  expect_true(all(P$x >= bb[1] & P$x <= bb[3] & P$y >= bb[2] & P$y <= bb[4]))
  inside <- sum(all$x >= bb[1] & all$x <= bb[3] & all$y >= bb[2] & all$y <= bb[4])
  expect_equal(nrow(P), inside)
  ## ids are positions in the node file
  expect_equal(P$z, all$z[P$id])
  expect_error(read_ept(d, bb + 1e6), "no EPT nodes")
})

test_that("points_transform() moves x and y and keeps the ids", {
  p <- read_points(system.file("extdata", "example.laz", package = "rlas"))
  attr(p, "meshspec")$crs <- "EPSG:26917"
  q <- points_transform(p, "EPSG:4326")
  expect_equal(attr(q, "meshspec")$from_crs, "EPSG:26917")
  expect_true(all(abs(q$points$x) <= 180))
  expect_equal(q$points$z, p$points$z)
  expect_equal(q$points$id, p$points$id)
})

test_that("face_stats() gives one row per face from the links", {
  skip_if_not_installed("cdtr")
  skip_if_not_installed("geometry")
  p <- read_points(system.file("extdata", "example.laz", package = "rlas"))
  tin <- points_tin(p, classes = NULL)
  fp <- face_points(tin, p)
  fs <- face_stats(tin, p, fp)
  expect_equal(nrow(fs), nrow(tin$triangles))
  expect_equal(sum(fs$n), nrow(fp))
  expect_true(all(fs$slope >= 0 & fs$slope <= 90))
  expect_true(all(fs$aspect >= 0 & fs$aspect < 360))
  f <- fp$face[which.max(fp$dz)]
  expect_equal(fs$dz_max[f], max(fp$dz[fp$face == f]))
  expect_true(all(is.na(fs$dz_max[fs$n == 0])))
  cls <- grep("^n_class_", names(fs), value = TRUE)
  expect_equal(rowSums(fs[cls]), fs$n)
})

test_that("a tilted plane has the slope and aspect it was built with", {
  skip_if_not_installed("cdtr")
  skip_if_not_installed("geometry")
  ## z falls 1 m per m eastward: slope 45 degrees, facing east (90)
  g <- expand.grid(x = 0:5, y = 0:5)
  pts <- structure(list(points = data.frame(x = g$x, y = g$y, z = -g$x, classification = 2L),
                        sources = data.frame(dsn = "plane")),
                   meshspec = list(version = "0.1", form = "points", crs = NA_character_),
                   class = "pointcore_points")
  tin <- points_tin(pts)
  fs <- face_stats(tin, pts, face_points(tin, pts))
  expect_equal(fs$slope, rep(45, nrow(fs)))
  expect_equal(fs$aspect, rep(90, nrow(fs)))
})
