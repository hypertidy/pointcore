laz <- function() system.file("extdata", "example.laz", package = "rlas")

test_that("the pdal engine gives the same form as rlas", {
  skip_if_not_installed("pdal")
  a <- read_points(laz())
  b <- read_points(laz(), engine = "pdal")
  expect_s3_class(b, "pointcore_points")
  expect_setequal(names(b$points), names(a$points))
  for (nm in names(a$points)) {
    expect_equal(b$points[[nm]], a$points[[nm]], info = nm)
  }
  expect_equal(attr(b, "meshspec")$crs, attr(a, "meshspec")$crs)
  expect_equal(b$sources$point_format, a$sources$point_format)
  expect_equal(b$sources$n, a$sources$n)
  expect_equal(b$sources$scale_x, a$sources$scale_x)
  expect_equal(b$sources$offset_z, a$sources$offset_z)
})

test_that("PDAL stages keep the ids of a full read", {
  skip_if_not_installed("pdal")
  p <- read_points(laz(), engine = "pdal")
  k <- read_points(laz(), engine = "pdal", stages = c(
    '{"type": "filters.range", "limits": "ReturnNumber[1:1]"}',
    '{"type": "filters.sort", "dimension": "Z", "order": "DESC"}'))
  expect_lt(nrow(k$points), nrow(p$points))
  expect_true(all(diff(k$points$z) <= 0))
  expect_equal(k$points$x, p$points$x[k$points$id])
  expect_equal(k$points$gps_time, p$points$gps_time[k$points$id])
  ## the same subset as the LASlib filter
  r <- read_points(laz(), filter = "-keep_first")
  expect_setequal(k$points$id, r$points$id)
  ## stages as a list
  l <- read_points(laz(), engine = "pdal",
                   stages = list(list(type = "filters.range", limits = "Classification[1:1]")))
  expect_true(all(l$points$classification == 1L))
})

test_that("bbox crops a file and keeps file positions", {
  skip_if_not_installed("pdal")
  p <- read_points(laz(), engine = "pdal")$points
  bb <- c(quantile(p$x, 0.2), quantile(p$y, 0.2), quantile(p$x, 0.8), quantile(p$y, 0.8))
  w <- read_points(laz(), engine = "pdal", bbox = bb)$points
  inside <- p$x >= bb[1] & p$x <= bb[3] & p$y >= bb[2] & p$y <= bb[4]
  expect_equal(nrow(w), sum(inside))
  expect_equal(w$z, p$z[w$id])
})

test_that("a COPC window reads its nodes, with no file positions", {
  skip_if_not_installed("pdal")
  copc <- tempfile(fileext = ".copc.laz")
  pdal::pdal_pipeline(sprintf('[%s, {"type": "writers.copc", "forward": "all", "filename": %s}]',
                              pdal::json_string(laz()), pdal::json_string(copc)))
  p <- read_points(laz(), engine = "pdal")$points
  bb <- c(quantile(p$x, 0.2), quantile(p$y, 0.2), quantile(p$x, 0.8), quantile(p$y, 0.8))
  w <- read_points(copc, engine = "pdal", bbox = bb)
  inside <- p$x >= bb[1] & p$x <= bb[3] & p$y >= bb[2] & p$y <= bb[4]
  expect_equal(nrow(w$points), sum(inside))
  expect_true(all(is.na(w$points$id)))
  expect_equal(sort(w$points$z), sort(p$z[inside]))
  ## read in full, a COPC file gives positions as usual
  f <- read_points(copc, engine = "pdal")
  expect_equal(f$points$id, seq_len(nrow(p)))
})

test_that("engine arguments are checked", {
  expect_error(read_points(laz(), bbox = c(0, 0, 1, 1)), "engine")
  skip_if_not_installed("pdal")
  expect_error(read_points(laz(), engine = "pdal", filter = "-keep_first"), "stages")
})
