#################################################################################################
# Recreate LandTrendr summaries, given raster files of the raw bands, Tasssled Cap components, 
# and elevation data
#
# Laurel Hopkins
# 03/03/2020
#################################################################################################
remove(list=ls())
library(sf)
library(raster)
library(dplyr)
# library(velox)

get_coordinates <- function(n) {  
  us <- raster::getData('GADM', country = 'US', level = 1)
  oregon <- us[us$NAME_1 == "Oregon",]
  grid <- makegrid(oregon, n = n) # n - approximate number of cells in grid
  grid <- SpatialPoints(grid, proj4string = CRS(proj4string(oregon)))
  grid <- grid[oregon, ]
  coords <- as.data.frame(grid) 
  return(coords)
}

get_feature <- function(id) {
  features <- read.csv("features.csv", header=FALSE) # TODO
  name = features[id,1]
  metric = features[id,2]
  buffer.radius = features[id,3]
  return(list(name=name, metric=metric, buffer.radius=buffer.radius))
}

# gen_feature <- function(id, n) {  
gen_feature <- function() {
  ptm <- proc.time()
  year = "2011"
  season = "summer"
  # feature <- get_feature(id)
  # name <- feature$name
  # metric <- feature$metric
  # buffer.radius <- feature$buffer.radius
  
  # if (!metric %in% c("mean", "std", "mean_std")) {
  #   stop("Not defined metric!")
  # }
  
  #########################################################
  # 2. Read in raster file OR__tif
  #########################################################
  
  # need raster package
  # base_path = "/nfs/guille/hutchinson/data/RS_data/"
  # if (name %in% c("B1", "B2", "B3", "B4", "B5", "B7")) {    
  #   # to read in raw bands
  #   OR_tif <- stack(paste0(base_path, "OR_", year, "_fittedImage.tif"))
  #   names(OR_tif) <-c("B1", "B2", "B3", "B4", "B5", "B7")  # rename the layers, if necessary
  # } else if (name %in% c("elevation", "slope", "aspect")) {
  #   # to read in elevation, slope, aspect
  #   OR_tif <- stack(paste0(base_path, "OR_Elevation_Slope_Aspect.tif"))
  # } else if (name %in% c("TCA", "TCB", "TCG", "TCW")) {
  #   # to read in Tasseled Cap components
  #   OR_tif <- stack(paste0(base_path, "OR_TasseledCap_", season, "_", year, ".tif"))
  #   ##### here
  #   OR_tif <- stack("../TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
  #   n <- 1000000
  #   name <- "TCB"
  #   metric <- "mean"
  #   buffer.radius <- 300
  #   #####
  #   names(OR_tif) <-c("TCB", "TCG", "TCW", "TCA")
  # } else {
  #   stop("Not defined feature name!")
  # }
  # OR_tif <- stack("OR_tasscap_summer_2011_cropped.tif")
  OR_tif <- stack("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
  n <- 1000000
  name <- "TCB"
  metric <- "mean"
  buffer.radius <- 250
  #####
  names(OR_tif) <-c("TCB", "TCG", "TCW", "TCA")
  
  #########################################################
  # 3. Create lat,long spatial points df
  #########################################################
  
  # create example dataframe
  coords <- get_coordinates(n) # n - approximate number of cells in grid
  #print(dim(coords))
  long <- coords[,1]
  lat <- coords[,2]
  years <- rep(year, dim(coords)[1])
  locations <- data.frame(lat, long, years)
  coords <- locations %>% dplyr::select(long, lat)
  
  # create a SpatialPointsDataFrame w/ longlat projection
  prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")  # keep this for lat/long data
  locations_pt <- SpatialPointsDataFrame(coords = coords, data = locations, proj4string = prj)
  # print("here")
  # reproject longlat points to match projection of OR_tif
  locations_pt <- spTransform(locations_pt, crs(OR_tif))
  # print("finished")
  # Check that projections match
  if (! identical(crs(locations_pt), crs(OR_tif))) { # should print true
    stop("Projection does not match!")
  }
  
  #########################################################
  # 4. With raster &  locations as as spatial points df,
  #    extract buffered summaries
  #
  #    Clipping: https://geocompr.robinlovelace.net/geometric-operations.html
  #########################################################
  
  summary_df <- locations
  
  print("entering apply")
  ptm <- proc.time()[[3]]
  sub <- locations_pt@coords[10000:11000,]
  res <- apply(sub, MARGIN=1, FUN=construct_mean, ras=OR_tif, b.rad=buffer.radius)
  print(paste0("this took:", proc.time()[[3]]-ptm, "seconds"))
  print(sum(res))
  # TCA here
  summary_df$TCA_buffered_250 <- res
  # 
  # point = st_sfc(st_point(c(locations_pt@coords[i,"long"], locations_pt@coords[i,"lat"])))
  # buffer = st_buffer(point, dist = buffer.radius) # create buffered region around point
  # buffer_sp <- as_Spatial(buffer) # convert buffered region to SpatialPolygonsDataFrame
  # 
  # buffered_values <- unlist(raster::extract(OR_tif$TCA, buffer_sp))
  # 
  # summary_df[i, "mean"] <- mean(buffered_values, na.rm=TRUE)
  
  # old (& slow?)
  #########################################################
  # for (i in 1:nrow(locations_pt)) {  # there may be a faster way to do this with lapply
  #   if (i %% 1000 == 0){
  #     print(proc.time() - ptm)
  #     cat(i, "/", nrow(locations_pt), "\n")
  #     ptm <- proc.time()
  #   }
  #   # create buffered region
  #   
  #   point = st_sfc(st_point(c(locations_pt@coords[i,"long"], locations_pt@coords[i,"lat"])))
  #   buffer = st_buffer(point, dist = buffer.radius) # create buffered region around point
  #   buffer_sp <- as_Spatial(buffer) # convert buffered region to SpatialPolygonsDataFrame
  #   
  #   # extract values from raster in the buffered region
  #   if (name == "elevation") {
  #     buffered_values <- unlist(raster::extract(OR_tif$elevation, buffer_sp))
  #   } else if (name == "aspect") {
  #     buffered_values <- unlist(raster::extract(OR_tif$aspect, buffer_sp))
  #   } else if (name == "slope") {
  #     buffered_values <- unlist(raster::extract(OR_tif$slope, buffer_sp))
  #   } else if (name == "TCB") {
  #     buffered_values <- unlist(raster::extract(OR_tif$TCB, buffer_sp))
  #   } else if (name == "TCG") {
  #     buffered_values <- unlist(raster::extract(OR_tif$TCG, buffer_sp))
  #   } else if (name == "TCW") {
  #     buffered_values <- unlist(raster::extract(OR_tif$TCW, buffer_sp)) 
  #     buffered_values <- unlist(raster::extract(OR_tif$TCW, buffer_sp)) 
  #   } else if (name == "TCA") {
  #     buffered_values <- unlist(raster::extract(OR_tif$TCA, buffer_sp)) 
  #   } else {
  #     stop("Not defined feature name!")
  #   }    
  #   
  #   # Use list apply to calculate the mean and standard deviation for the buffered region
  #   if (metric == "mean") {
  #     summary_df[i, "mean"] <- mean(buffered_values, na.rm=TRUE)
  #   } else if (metric == "std") {
  #     summary_df[i, "std"] <- sd(buffered_values, na.rm=TRUE)
  #   } else { # both
  #     summary_df[i, "mean"] <- mean(buffered_values, na.rm=TRUE)
  #     summary_df[i, "std"] <- mean(buffered_values, na.rm=TRUE)
  #   }
  # }
  #########################################################
  
  summary_df <- dplyr::select(summary_df, -c("years"))
  if (metric %in% c("mean", "std")) {
    colnames(summary_df)[3] <- paste0(name, "_", metric, "_", buffer.radius)
  } else { # mean_std
    colnames(summary_df)[3] <- paste0(name, "_mean_", buffer.radius)
    colnames(summary_df)[4] <- paste0(name, "_std_", buffer.radius)
  }
  
  if (!dir.exists("outputs")) { 
    dir.create("outputs")
  }
  
  filename = paste0("outputs/", id, "_", name, "_", metric, "_", buffer.radius, ".csv") #TODO
  write.table(summary_df, filename, row.names=FALSE, quote = FALSE, sep=",")
  print(proc.time() - ptm)
}

construct_mean <- function(in.pt, ras, b.rad){
  
  point = st_sfc(st_point(c(in.pt[["long"]], in.pt[["lat"]])))
  buffer = st_buffer(point, dist = b.rad) # create buffered region around point
  buffer_sp <- as_Spatial(buffer) # convert buffered region to SpatialPolygonsDataFrame
  
  # here: TCA
  ptm <- proc.time()[[3]]
  buffered_values <- unlist(raster::extract(ras$TCA, buffer_sp))
  z <- proc.time()[[3]]-ptm

  # summary_df[i, "mean"] <- mean(buffered_values, na.rm=TRUE)
  # return(mean(buffered_values, na.rm=TRUE))
  
  # timing purposes
  return(z)
}



gen_feature()

# args <- commandArgs(trailingOnly = TRUE)
# gen_feature(as.integer(args[1]), as.integer(args[2]))
