# file to construct adversarial datasets
# 10.4.21
# Mark Roth

library(sf)
library(randomcoloR)

setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("run experiments/adversarial/syn_helper.R")
source("helper/helpers.R")

# function to create a 3x3 grid of points around the 
# provided checklist with a distance between each point
# specified by increment
gridAroundPoint <- function(xCoord, yCoord, site, increment=.0001){
  cids <- create_checklist_ids(9)
  xy.TL <- data.frame(x=(xCoord - increment), y=(yCoord + increment), site=site, checklist_id=cids[1])
  xy.ML <- data.frame(x=(xCoord - increment), y=(yCoord), site=site, checklist_id=cids[2])
  xy.BL <- data.frame(x=(xCoord - increment), y=(yCoord - increment), site=site, checklist_id=cids[3])
  xy.TM <- data.frame(x=xCoord, y=(yCoord + increment), site=site, checklist_id=cids[4])
  xy.MM <- data.frame(x=xCoord, y=yCoord, site=site, checklist_id=cids[5])
  xy.BM <- data.frame(x=xCoord, y=(yCoord - increment), site=site, checklist_id=cids[6])
  xy.TR <- data.frame(x=(xCoord + increment), y=(yCoord + increment), site=site, checklist_id=cids[7])
  xy.MR <- data.frame(x=(xCoord + increment), y=(yCoord), site=site, checklist_id=cids[8])
  xy.BR <- data.frame(x=(xCoord + increment), y=(yCoord - increment), site=site, checklist_id=cids[9])
  return(rbind(
    TL=xy.TL, ML=xy.ML, BL=xy.BL, TM=xy.TM, 
    MM=xy.MM,
    BM=xy.BM, TR=xy.TR, MR=xy.MR, BR=xy.BR
  ))
}

# n.centered_pts: number of centered points in a 9x9 grid
# total number of points will be 8 * n.centered_pts
construct_datapoints <- function(sp_region, n.centered_pts, OR.env, inc = .0001, toPlot=F){
  pts <- spsample(sp_region, n.centered_pts, type="random")@coords
  ptm <- proc.time()
  for(i in 1:nrow(pts)){
    if(i %% 10 == 0){
      print(i)
      print(proc.time() - ptm)
      ptm <- proc.time()
    }
    grid_df <- gridAroundPoint(pts[[i, 'x']], pts[[i, 'y']], site = i, increment = inc)
    lon <- grid_df$x
    lat <- grid_df$y
    site <- grid_df$site
    cids <- grid_df$checklist_id
    coordinates(grid_df) <- c("x", "y")
    crs(grid_df) <- latlong <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    # projecting the lat/long points into aea to extract env variables
    pts_proj <- spTransform(grid_df, CRSobj=crs(OR.env))
    
    if(i == 1){
      grid_pts.env <- data.frame(
        lon=lon, 
        lat=lat, 
        site=site,
        checklist_id=cids,
        raster::extract(OR.env, pts_proj)
      )  
    } else {
      add_df <- data.frame(
        lon=lon, 
        lat=lat, 
        site=site, 
        checklist_id=cids,
        raster::extract(OR.env, pts_proj)
      )
      grid_pts.env <- rbind(grid_pts.env, add_df)
      
      if(toPlot){
        points(add_df$lon, add_df$lat, col=randomColor())
      }
    }
    
  }
  return(grid_pts.env)
}

sampling_region <- readRDS("../TassledCapOR/sampling_region")
tcap.OR <- stack("../TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
names(tcap.OR) <-c("TCB", "TCG", "TCW", "TCA")
rounded.df <- construct_datapoints(sampling_region, 120, tcap.OR, inc = .0005, toPlot = T)

obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]

det.df <- WETA_2017[,unlist(covObj$det_cov)]

rounded_norm_det.df <- prep_syn_df(rounded.df, det.df)


########
# turning OR outline into sp
# for sampling 
########
library('maps')

define_samp_region <- function(){
  map <- map_data('state')
  OR <- subset(map, region %in% "oregon")
  coordinates(OR) <- c("long", "lat")
  crs(OR) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  OR_sf <- st_as_sf(OR, coords = coordinates(OR))
  
  # defining larger rectangle to crop the border
  box <- data.frame(
    x=c(-123, -123, -130, -130),
    y=c(44.5, 40, 44.5, 40)
  )
  coordinates(box) <- c("x", "y")
  crs(box) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  # cropping border
  cropped.border <- st_crop(OR_sf, st_bbox(box))
  separated_coord <- data.frame(
    long = unlist(map(cropped.border$geometry,1)),
    lat = unlist(map(cropped.border$geometry,2))
  )
  
  # adding a point 
  separated_coord <- rbind(separated_coord, data.frame(long=c(-123), lat=c(44.5)))
  coordinates(separated_coord) <- c("long", "lat")
  sp_region <- SpatialPolygons(list(Polygons(list(Polygon(separated_coord)), ID=1)))
  # saveRDS(sp_region, "../TassledCapOR/sampling_region")
  return(sp_region)
}

