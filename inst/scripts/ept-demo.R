## A window of public USGS 3DEP lidar to a ground TIN, per-face summaries,
## Arrow tables and the meshview browser viewer.
##
## Needs meshcore with write_mesh_arrow() (hypertidy/meshcore, branch
## mesh-arrow until merged), nanoarrow >= 0.5.0, cdtr and geometry.
library(pointcore)
library(meshcore)

base <- "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/AK_Juneau_2012"
## 400 x 400 Web Mercator units (about 210 m on the ground at 58 N)
cx <- -14963532.6
cy <- 8031006.5
bb <- c(cx - 200, cy - 200, cx + 200, cy + 200)

p <- read_ept(base, bb)                     # 136,021 returns from 19 node files
p <- points_transform(p, "EPSG:32608")      # metres on the ground (UTM 8N)
tin <- points_tin(p, classes = 2L)          # ground returns only
fp <- face_points(tin, p, below = 1)        # every return within 1 m below or any height above
fs <- face_stats(tin, p, fp)

## faces with something more than 2 m above the ground
tall <- fs[!is.na(fs$dz_max) & fs$dz_max > 2, ]
nrow(tall)

out <- file.path(tempdir(), "juneau")
write_mesh_arrow(tin, file.path(out, "mesh"), face_points = fp, face_stats = fs)
write_mesh_arrow(p, file.path(out, "points"))
mesh_viewer(out)   # serve `out` (python3 -m http.server) or choose the .arrows files
