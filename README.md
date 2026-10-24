# pointcore (working name)

Prototype of the point-cloud form for mesh spec v0.1: LAS, LAZ and COPC
returns as rows with ids, a TIN that keeps each vertex's return, and a link
table from every return to its triangle. Builds on
[meshcore](https://github.com/hypertidy/meshcore).

```r
remotes::install_github("hypertidy/pointcore")
```

```r
library(pointcore)
p <- read_points(system.file("extdata", "example.laz", package = "rlas"))
tin <- points_tin(p, classes = NULL)   # classes = 2L (ground) by default
fp <- face_points(tin, p, above = 2)   # returns within 2 m above the surface
```

With `engine = "pdal"` the same form comes through a PDAL pipeline
([hypertidy/pdal](https://github.com/hypertidy/pdal)): any format or URL
PDAL reads, PDAL filters as `stages`, and a `bbox` that reads only the
nodes it needs from COPC and EPT. Ids stay file positions through filters
that drop or reorder points (and are `NA` for a COPC or EPT window).

```r
g <- read_points(f, engine = "pdal",
                 stages = '{"type": "filters.range", "limits": "Classification[2:2]"}')
w <- read_points("survey.copc.laz", engine = "pdal", bbox = c(xmin, ymin, xmax, ymax))
```

Survey-scale clouds come from Entwine Point Tiles by window, and the
links summarise to one row per face:

```r
p <- read_ept("https://s3-us-west-2.amazonaws.com/usgs-lidar-public/AK_Juneau_2012", bbox)
p <- points_transform(p, "EPSG:32608")   # 3DEP EPT is Web Mercator
tin <- points_tin(p)
fp <- face_points(tin, p)
fs <- face_stats(tin, p, fp)             # area, slope, aspect, n, dz_max, n_class_*, ...
```

Each return keeps its node file (`sources$key`) and its position in that
file (`id`). With meshcore's `write_mesh_arrow()` and `mesh_viewer()` the
TIN opens in a browser where hovering a face lists its summary and its
returns: see `inst/scripts/ept-demo.R`.

With meshcore's `drape_uv()`, `as_ray_mesh()` and `pick_render()`, a pixel
of a rayvertex render leads back to its face, the face's returns and the
ortho pixel: see `inst/scripts/points-demo.R` (figure in
`inst/scripts/points-demo.png`) and the draft form in
[meshcore docs/points-form-draft.md](https://github.com/hypertidy/meshcore/blob/main/docs/points-form-draft.md).
