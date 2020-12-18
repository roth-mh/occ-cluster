## Laurel Hopkins
## 11/23/2020

library(sf)
library(rgdal)  # reads entire raster into memnory  
library(raster) # does not require rasters to be held in-memory which is beneficial for 
                # analyzing and predicting to larger files

# Partially based on this post: https://geocompr.robinlovelace.net/intro.html

#########################################################
# 1. Read in GEE exports and merge -- ONLY NEED TO DO 
#    THIS ONCE
#########################################################
# read in files as rasters
files <- list.files(path="fittedImages", pattern="*.tif$", full.names=TRUE, recursive=FALSE)
ras1 <- stack(files[1])  # use stack() to read in all bands, raster() only reads in single band; brick() may be faster, but less flixible
ras2 <- stack(files[2])
ras3 <- stack(files[3])
ras4 <- stack(files[4])

# plot with spplot to view all layers
spplot(ras5)
plot(ras2[[1]])

# merge the rasters to form a single raster
x <- list(ras1, ras2, ras3, ras4)
x$filename <- 'OR_2020_fittedImage.tif'
x$overwrite <- TRUE
OR_2020 <- do.call(merge, x)
names(OR_2020) <-c("B1", "B2", "B3", "B4", "B5", "B7")

# save as a geotiff file, may not be necessary, it seems like merge saves the tif
#rf <- writeRaster(OR_2020, filename="OR_2020_fittedLandTrendr.tif", format="GTiff", overwrite=TRUE)

# plot tif
# plot(OR_2020)
plotRGB(OR_2020, r="B7", g="B4", b ="B3", stretch="lin")

#########################################################
# 2. Read in OR_2020 tif
#########################################################
OR_2020 <- stack("OR_2020_fittedImage.tif")

names(OR_2020) <-c("B1", "B2", "B3", "B4", "B5", "B7")

plotRGB(OR_2020, r="B7", g="B4", b ="B3", stretch="lin")


#########################################################
# 3. Turn species records into spatial points df
#########################################################
library(dplyr)
library(tidyr)

# read in a random species to get lat, long data
Data.directory <-"C:\\Users\\Laurel\\Documents\\Oregon State\\Research\\ICB\\"
species <- "Western Tanager"
species_file <- paste0(Data.directory, "SpeciesCounts/", species, "_Counts_062320.csv")
bird_data <- read.csv(file=species_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)

# get lat, long from Unique_Checklist_ID
bird_data$lat <- sapply(strsplit(bird_data$Unique_Checklist_ID, "_"), `[`, 1)
bird_data$long <- sapply(strsplit(bird_data$Unique_Checklist_ID, "_"), `[`, 2)
bird_data$date <- sapply(strsplit(bird_data$Unique_Checklist_ID, "_"), `[`, 3)
bird_data$year <- sapply(strsplit(bird_data$date, "-"), `[`, 1)
bird_data$lat <- as.double(bird_data$lat)
bird_data$long <- as.double(bird_data$long)
bird_data$year <- as.double(bird_data$year)

# years that have bird records
sort(unique(bird_data$year))

# save df as a csv 
write.csv(bird_data, file="OR2020_recordLocations_withYear.csv", col.names=TRUE, row.names=FALSE)

# create a SpatialPointsDataFrame with species data w/ longlat projection, then reproject to match OR_2020 CRS
prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
select <- dplyr::select #this tells R to use the select() function from the dplyr package, not the select function stored in some of the spatial packages
coords <- bird_data %>% select(long, lat)
bird_data_pt <- SpatialPointsDataFrame(coords = coords, data = bird_data, proj4string = prj)
# reproject to match OR_2020 crs
bird_data_pt <- spTransform(bird_data_pt, crs(OR_2020))
bird_data_pt@proj4string
plot(bird_data_pt)

# Check that projections match
identical(crs(bird_data_pt), crs(OR_2020)) 


#########################################################
# 4. With Oregon tif & record locations as as spatial 
#    points df, clip buffered tifs centered at species 
#    records
#########################################################

bird_data_pt_mini <- bird_data_pt[1:3,]
# define how large the buffer should be (in meters)
buffer.radius <- 2400

for (i in 1:nrow(bird_data_pt)) {
  point = st_sfc(st_point(c(bird_data_pt@coords[i,"long"], bird_data_pt@coords[i,"lat"])), crs=st_crs(OR_2020))
  buffer = st_buffer(point, dist = buffer.radius) # convert point to buffered circles
  bb = st_bbox(buffer) # create bounging box
  box = st_as_sfc(bb)
  spdf <- as_Spatial(box) # convert to SpatialPolygonsDataFrame
  test <- crop(OR_2020, spdf)
  rf <- writeRaster(test, filename=paste0("test_", i,".tif"), format="GTiff", overwrite=TRUE)
  #rf <- writeRaster(tess, filename=file.path(tmp, "test.tif"), format="GTiff", overwrite=TRUE)
  
}



files <- list.files(path="test", pattern="*.tif$", full.names=TRUE, recursive=FALSE)
input <- stack(files[3])
plot(input) # plots each layer individually

names(input) <-c("B1", "B2", "B3", "B4", "B5", "B7")
plotRGB(input, r="B7", g="B4", b ="B3", stretch="lin")

# Clipping: https://geocompr.robinlovelace.net/geometric-operations.html
# may need to first specify a coordinate reference system 
b = st_sfc(st_point(c(0, 1))) # create a point
b = st_buffer(b, dist = 1) # convert point to buffered circles
plot(b)
bb = st_bbox(b) # create bounging box
box = st_as_sfc(bb)
plot(box, add=TRUE)
