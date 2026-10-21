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

With meshcore's `drape_uv()`, `as_ray_mesh()` and `pick_render()`, a pixel
of a rayvertex render leads back to its face, the face's returns and the
ortho pixel: see `inst/scripts/points-demo.R` (figure in
`inst/scripts/points-demo.png`) and the draft form in
[meshcore docs/points-form-draft.md](https://github.com/hypertidy/meshcore/blob/main/docs/points-form-draft.md).
