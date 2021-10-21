#####
# adversarial experiments
#####

library(raster)

# setwd("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
# tcap.OR <- stack("../TassledCapOR/OR_TasseledCap_summer_2011.tiff")
tcap.OR <- stack("../TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
names(tcap.OR) <-c("TCB", "TCG", "TCW", "TCA")
# plot(tcap.OR$TCB)

tcap.OR_norm <- scale(tcap.OR)
# writeRaster(tcap.OR_norm, "OR_tasscap_summer_2011_cropped_normalized.tiff")

latlong <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
tcap_reproj <- projectRaster(tcap.OR_norm, crs=latlong)
writeRaster(tcap_reproj, "OR_tasscap_summer_2011_cropped_normalized_latlong.tiff")

# box <- data.frame(
#   x=c(-2000000, -2000000, -2400000, -2400000), 
#   y=c(2400000, 2800000, 2400000, 2800000)
# ) 
pts <- data.frame(
  x=runif(5, -123.2, -123),
  y=runif(5, 44, 44.2)
)
coordinates(pts) <- c("x", "y")
crs(pts) <- latlong
pts_proj <- spTransform(pts, CRSobj=crs(tcap.OR))
points(pts_proj)
# 
# tcap.OR_cropped <- crop(tcap.OR, box)
# writeRaster(tcap.OR_cropped, "OR_tasscap_summer_2011_cropped.tiff")

# 
# or.2020data <- read.delim("../checklist data/oregon2020_observations_021220_032720FIX.csv", header=TRUE, sep = ",", strip.white = TRUE)
# or.2020.WETA <- or.2020data[or.2020data$CommonName == "Western Tanager",]
# or.2020.WETA_region <- subset(or.2020.WETA, or.2020.WETA$Latitude <= 44.5)
# or.2020.WETA_region <- subset(or.2020.WETA_region, or.2020.WETA_region$Longitude <= -123)
# 
# points(or.2020.WETA_region$Longitude, or.2020.WETA_region$Latitude, col="red")
