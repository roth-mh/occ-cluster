library(sp)
library(rgdal)
library(raster)
library(gdistance)
library(sf)
library(ggplot2)

setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("site construction/NLCDvals.R")
source("helper/helpers.R")


#######################################
# NEED TO USE projection from ras.OR2 #
#######################################

# ras.OR has land cover type but ras.OR2 has the "correct" (?) projection
ras.OR <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/NLCD_2011_OR.tif")
ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/OR_2020_fittedImage.tif")

plot(ras.OR2)
plot(ras.OR)
###################################################
# DEFINING THE BOUNDARIES OF THE AREA OF INTEREST #
###################################################
coords <- matrix(nrow=2, ncol=2)
coords[1,1] <- as.double(-123)
coords[1,2] <- as.double(44.5)

coords[2,1] <- as.double(-124.5459)
coords[2,2] <- as.double(42.00583)

coords_df <- as.data.frame(coords)
colnames(coords_df) <- c("long", "lat")

######################################
# SUPPLYING LAT/LON COORDINATES WITH #
# PROPER PROJECTION & REPROJECTING ###
######################################
prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj_coords_df <- SpatialPointsDataFrame(coords = coords_df, 
                                         data = as.data.frame(c(0,00)), 
                                         proj4string = prj)

reproj_coords <- spTransform(proj_coords_df, crs(ras.OR2))

# crs(ras.OR)
# crs(ras.OR2)
crs(ras.OR) <- crs(ras.OR2)
# should be false
# identical(crs(reproj_coords), crs(ras.OR))

#############################
# MATH STUFF TO DETERMINE ###
# THE NEW COORDINATES GIVEN #
# OUR TIF IS ANGLED #########
#############################
BL <- as.data.frame(list(x=reproj_coords@bbox[1,1], y=reproj_coords@bbox[2,1]))
TL <- as.data.frame(list(x=reproj_coords@bbox[1,1], y=reproj_coords@bbox[2,2]))
TR <- as.data.frame(list(x=reproj_coords@bbox[1,2]+20000, y=reproj_coords@bbox[2,2]))
BR <- as.data.frame(list(x=reproj_coords@bbox[1,2]+20000, y=reproj_coords@bbox[2,1]))

BL$x <- BL$x - 15000
x <- rbind(BL, TL, TR, BR)

dist_mat <- as.matrix(dist(x))
THETA <- 21.5

A <- 1/(cos(THETA*pi/180)/dist_mat[1,2])
B <- sin(THETA*pi/180)*A
B.C <- dist_mat[2,3]
C <- B.C - B
D <- C*sin((90-THETA)*pi/180)
D <- D + 5000

del_x <- D*cos(THETA*pi/180)
del_y <- D*sin(THETA*pi/180)

x4 <- as.data.frame(list(x=TR$x-del_x, y=TR$y+del_y))
x3 <- as.data.frame(list(x=BL$x+del_x, y=BL$y-del_y))


# I'm not 100% sure this region captures all of the checklists we are 
# hoping to observe -- we should check how many checklists are in this
# region
coords <- rbind(x4, TR, x3, BL)

points(coords, col="white", cex=2)
#############################
#############################

#############################
# TURNING OUR COORDS INTO A #
# POLYGON; CROPPING AND #####
# MASKING ###################
#############################
c <- Polygon(coords)
l <- Polygons(list(c), 1)
sps <- SpatialPolygons(list(l), proj4string = crs(ras.OR))


test <- crop(ras.OR, sps)
cropped <- mask(test, sps)

plot(cropped)
#############################
# PLOTS & VALUE EXTRACTION ##
#############################

# plot(ras.OR)
# plot(test, add=T, col="black")
# plot(cropped, add=T, col="red")

#####################################################################################
#####################################################################################
# ##############
# # MINI EXAMPLE 1 IN SMALL BOX
# ##############
# x_val <- mean(c(cropped@extent@xmax, cropped@extent@xmin)) - 20000
# y_val <- mean(c(cropped@extent@ymax, cropped@extent@ymin)) - 20000
# 
# xmax <- x_val + 10000
# xmin <- x_val - 10000
# ymax <- y_val + 10000
# ymin <- y_val - 10000
# 
# test_ex1 <- crop(cropped, c(xmin, xmax, ymin, ymax))
# tr1 <- transition(test_ex1, transitionFunction = conductanceVals, directions = 8)
# trCorr <- geoCorrection(tr1, type = 'c')
# 
# c <- matrix(nrow=4, ncol=2)
# c[1,] <- c(xmax,ymax)
# c[2,] <- c(xmax,ymin)
# c[3,] <- c(xmin,ymin)
# c[4,] <- c(xmin,ymax)
# points(c, col="blue", cex=2)
# 
# ##############
# # MINI EXAMPLE 2 IN SMALL BOX
# ##############
# x_val2 <- WETA_reproj_coords[WETA_reproj_coords$`WETA_coords_df$checklist_id` == "S37978885",]@coords[[1]]
# y_val2 <- WETA_reproj_coords[WETA_reproj_coords$`WETA_coords_df$checklist_id` == "S37978885",]@coords[[2]]
# 
# xmax2 <- x_val2 + 10000
# xmin2 <- x_val2 - 10000
# ymax2 <- y_val2 + 10000
# ymin2 <- y_val2 - 10000
# 
# test_ex2 <- crop(cropped, c(xmin2, xmax2, ymin2, ymax2))
# tr2 <- transition(test_ex2, transitionFunction = conductanceVals, directions = 8)
# trCorr2 <- geoCorrection(tr2, type = 'c')
# 
# c2 <- matrix(nrow=4, ncol=2)
# c2[1,] <- c(xmax2,ymax2)
# c2[2,] <- c(xmax2,ymin2)
# c2[3,] <- c(xmin2,ymin2)
# c2[4,] <- c(xmin2,ymax2)
# points(c2, col="black", cex=2)
# 
# g2 <- graph_from_adjacency_matrix(trCorr2@transitionMatrix, mode = "undirected", weighted = T)
# 
# cellNum2 <- cellFromXY(trCorr, WETA_in_box@coords)
# 
# d2.A <- distances(g2, cellNum2[1], cellNum2[2])
# d2.B <- distances(g2, cellNum2[3], cellNum2[4])
# d2.C <- distances(g2, cellNum2[4], cellNum2[5])
# d2.D <- distances(g2, cellNum2[5], cellNum2[6])
# if(!(d2.A > d2.B)){
#   print("ERROR")
# }
# #####################################
# # LOAD DATA & REPROJECT WETA COORDS #
# #####################################
# o <- load.WETA()
# WETA_2017 <- o[[1]]
# covObj <- o[[2]]
# 
# WETA_coords_df <- as.data.frame(list(long=WETA_2017$longitude, lat=WETA_2017$latitude))
# WETA_coords_df$checklist_id <- WETA_2017$checklist_id
# 
# prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# WETA_proj_coords_df <- SpatialPointsDataFrame(coords = WETA_coords_df[,1:2], 
#                                               data = as.data.frame(WETA_coords_df$checklist_id),
#                                               proj4string = prj)
# 
# 
# 
# WETA_reproj_coords <- spTransform(WETA_proj_coords_df, crs(ras.OR2))
# WETA_reproj_coords
# 
# # EXAMPLE 1
# ##############################
# # SUBSET WETA FOR SMALL BOX; #
# # CALC DIST MATRIX           #
# ##############################
# WETA_in_box <- subset(WETA_reproj_coords, WETA_reproj_coords@coords[,'long'] >= xmin)
# WETA_in_box <- subset(WETA_in_box, WETA_in_box@coords[,'long'] <= xmax)
# WETA_in_box <- subset(WETA_in_box, WETA_in_box@coords[,'lat'] >= ymin)
# WETA_in_box <- subset(WETA_in_box, WETA_in_box@coords[,'lat'] <= ymax)
# 
# # EXAMPLE 2
# ##############################
# # SUBSET WETA FOR SMALL BOX; #
# # CALC DIST MATRIX           #
# ##############################
# WETA_in_box2 <- subset(WETA_reproj_coords, WETA_reproj_coords@coords[,'long'] >= xmin2)
# WETA_in_box2 <- subset(WETA_in_box2, WETA_in_box2@coords[,'long'] <= xmax2)
# WETA_in_box2 <- subset(WETA_in_box2, WETA_in_box2@coords[,'lat'] >= ymin2)
# WETA_in_box2 <- subset(WETA_in_box2, WETA_in_box2@coords[,'lat'] <= ymax2)
# 
# 
# 
# #######
# # HAVE:
# #   - graph as an igraph object
# # TODO: 
# #   - find the mapping btwn raster cell and vertex
# #   - find the vertex id for each checklist
# #   - calculate shortest path from a set of vertices 
# #       (eventually, those w/in range, or neighbors)
# #######
# tr1 <- transition(test_ex1, transitionFunction = conductanceVals, directions = 8)
# trCorr <- geoCorrection(tr1, type = 'c')
# 
# g <- graph_from_adjacency_matrix(trCorr@transitionMatrix, mode = "undirected", weighted = T)
# 
# cellNum <- cellFromXY(trCorr, WETA_in_box@coords)
# # WETA_in_box@data
# 
# 
# ############
# # EXAMPLE1; 8 CHECKLISTS
# ############
# d1.A <- distances(g, cellNum[1], cellNum[2])
# d1.B <- distances(g, cellNum[3], cellNum[4])
# d1.C <- distances(g, cellNum[4], cellNum[5])
# d1.D <- distances(g, cellNum[5], cellNum[6])
# if(!(d1.A > d1.B)){
#   print("ERROR")
# }
# 
# 
# 
# 
# WETA_in_box@data



####################
# DISP DIST MATRIX 
# FOR COMMUTE DISTANCE
####################
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   return(cormat)
# }
# 
# cdCorr <- commuteDistance(trCorr, WETA_in_box@coords)
# 
# dim <- 8
# chklsts<-WETA_in_box$`WETA_coords_df$checklist_id`
# cormat <- get_upper_tri(as.matrix(cdCorr))
# colnames(cormat) <- as.character(chklsts)
# rownames(cormat) <- as.character(chklsts)
# m <- melt(cormat)
# 
# m$value_hundT<-round(m$value/100000, 2)
# ggplot(data = m, aes(x=Var1, y=Var2, fill=value_hundT)) + 
#   geom_tile() +  geom_text(aes(label = value_hundT), color = "black", size = 3) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0))


# ######
# # Fun plotting shortest path
# ######
# xmax2 <- x_val + 5000
# xmin2 <- x_val + 3000
# ymax2 <- y_val - 3000
# ymin2 <- y_val - 4000
# 
# x <- shortestPath(trCorr, c(-2223031, 2553257), c(-2222893, 2553097), output="SpatialLines")
# test_ex2 <- crop(cropped, c(xmin2, xmax2, ymin2, ymax2))
# 
# plot(test_ex2)
# lines(x, col="red")

####################
#####################################################################################
#####################################################################################






####################
# add other raster layers
####################
hwy_50 <- shapefile("site construction/input data/Hwy_Riv_Buffers/hwy_diss_buff_50.shp")
hwy_100 <- shapefile("site construction/input data/Hwy_Riv_Buffers/hwy_diss_buff_100.shp")
hwy_200 <- shapefile("site construction/input data/Hwy_Riv_Buffers/hwy_diss_buff_200.shp")


hwy_50_rast <- shape.to.raster(hwy_50, cropped, sps)
hwy_100_rast <- shape.to.raster(hwy_100, cropped, sps)
hwy_200_rast <- shape.to.raster(hwy_200, cropped, sps)

cropped.hwy <- addLayer(cropped, 2, hwy_50_rast)
names(cropped.hwy) <- c(names(cropped.hwy)[[1]], "hwy_50")
cropped.hwy <- addLayer(cropped.hwy, 3, hwy_100_rast)
names(cropped.hwy) <- c(names(cropped.hwy)[[1]], "hwy_50", "hwy_100")
cropped.hwy <- addLayer(cropped.hwy, 4, hwy_200_rast)
names(cropped.hwy) <- c(names(cropped.hwy)[[1]], "hwy_50", "hwy_100", "hwy_200")

cropped.hwy[['hwy_50']][is.na(cropped.hwy[['hwy_50']])] <- 0
cropped.hwy[['hwy_100']][is.na(cropped.hwy[['hwy_100']])] <- 0
cropped.hwy[['hwy_200']][is.na(cropped.hwy[['hwy_200']])] <- 0

# writeValues(cropped.hwy, "site construction/cropped_hwy_rast")
outfile <- writeRaster(
  cropped.hwy, 
  filename='site construction/cropped_hwy.tif',
  format="GTiff", 
  overwrite=TRUE,
  options=c("INTERLEAVE=BAND","COMPRESS=LZW")
  )

cropped.hwy <- brick("site construction/cropped_hwy.tif")

names(cropped.hwy) <- c("NLCD_2011_OR", "hwy_50", "hwy_100", "hwy_200")

# cropped.hwy <- hwy_rast$cropped.hwy
# making NA values 0 for transition function

# 
# test_ex1 <- crop(cropped.hwy100, c(xmin, xmax, ymin, ymax))
# tr1 <- transition.BRICK(test_ex1, transitionFunction = conductanceVals, directions = 8)
# trCorr <- geoCorrection(tr1, type = 'c')


###
# adding river buffers
###
# L_riv_50 <- shapefile("site construction/input data/Hwy_Riv_Buffers/lar_riv_buff_50_diss.shp")
# L_riv_100 <- shapefile("site construction/input data/Hwy_Riv_Buffers/lar_riv_buff_100_diss.shp")
# M_riv_50 <- shapefile("site construction/input data/Hwy_Riv_Buffers/med_riv_buff_50_diss.shp")
# 
# # errors out; and L_riv_100 has values between 0 and 2 ... maybe they were combined?
# M_riv_100 <- shapefile("site construction/input data/Hwy_Riv_Buffers/med_riv_buff_100_diss.shp")
# 
# # L_50_trans <- spTransform(L_riv_50, crs(cropped))
# # crop_L_50 <- crop(L_50_trans, sps)
# # L_50_rast <- rasterize(crop_L_50, cropped)
# 
# # reproject into same coordinate system
# L_riv_50_rast <- shape.to.raster(L_riv_50, cropped, sps)
# L_riv_100_rast <- shape.to.raster(L_riv_100, cropped, sps)
# M_riv_50_rast <- shape.to.raster(M_riv_50, cropped, sps)
# # M_riv_100_rast <- shape.to.raster(M_riv_100, cropped, sps)
# 
# 
# L_riv_50_rast[is.na(L_riv_50_rast)] <- 0
# L_riv_100_rast[is.na(L_riv_100_rast)] <- 0
# M_riv_50_rast[is.na(M_riv_50_rast)] <- 0
# # M_riv_100_rast[is.na(M_riv_100_rast)] <- 0
# 
# 
# 
# writeRaster(L_riv_50_rast, 
#             "site construction/L_riv_50_rast")
# writeRaster(L_riv_100_rast, 
#             "site construction/L_riv_100_rast")
# writeRaster(M_riv_50_rast, 
#             "site construction/M_riv_50_rast")
# 
# L_riv_50_rast <- raster("site construction/L_riv_50_rast.grd")
# L_riv_100_rast <- raster("site construction/L_riv_100_rast.grd")
# M_riv_50_rast <- raster("site construction/M_riv_50_rast.grd")
# table(L_riv_100_rast)
# 
# # change names 1 by 1 because sometimes the order gets flipped around
# cropped.hwy.riv <- addLayer(cropped.hwy, 5, L_riv_50_rast)
# names(cropped.hwy.riv) <- c(names(cropped.hwy.riv)[[1]], "hwy_50", "hwy_100", "hwy_200", "L_riv_50")
# cropped.hwy.riv <- addLayer(cropped.hwy.riv, 6, L_riv_100_rast)
# names(cropped.hwy.riv) <- c(names(cropped.hwy.riv)[[1]], "hwy_50", "hwy_100", "hwy_200", "L_riv_50", "L_riv_100")
# cropped.hwy.riv <- addLayer(cropped.hwy.riv, 7, M_riv_50_rast)
# names(cropped.hwy.riv) <- c(names(cropped.hwy.riv)[[1]], "hwy_50", "hwy_100", "hwy_200", "L_riv_50", "L_riv_100", "M_riv_50")
# 
# 
# shape.to.raster <- function(shp.obj, cropped.raster, spatial.polygons){
#   # shp_obj <- shapefile(shp.file)
#   buffer_trans <- spTransform(shp.obj, crs(cropped.raster))
#   
#   crop.shp <- crop(buffer_trans, spatial.polygons)
#   buffer_rast <- rasterize(crop.shp, cropped.raster)
#   return(buffer_rast)
# }
# 
# 
# outfile <- writeRaster(
#   cropped.hwy.riv, 
#   filename='site construction/cropped_hwy_riv.tif',
#   format="GTiff", 
#   overwrite=TRUE,
#   options=c("INTERLEAVE=BAND","COMPRESS=LZW")
# )

cropped.hwy.riv <- brick("site construction/cropped_hwy_riv.tif")

names(cropped.hwy.riv) <- c("NLCD_2011_OR", "hwy_50", "hwy_100", "hwy_200", "L_riv_50", "L_riv_100", "M_riv_50")

# resistance values for highway and river buffers
# hwy 50:     3
# hwy 100:    2
# hwy 200:    1.5

# L riv 50:   3
# L riv 100:  2
# M riv 50:   2.5
# M riv 100:  1.5
cropped.hwy.riv[['hwy_50']][cropped.hwy.riv[['hwy_50']] == 1] <- 3
cropped.hwy.riv[['hwy_100']][cropped.hwy.riv[['hwy_100']] == 1] <- 2
cropped.hwy.riv[['hwy_200']][cropped.hwy.riv[['hwy_200']] == 1] <- 1.5

cropped.hwy.riv[['L_riv_50']][cropped.hwy.riv[['L_riv_50']] == 1] <- 3
cropped.hwy.riv[['L_riv_100']][cropped.hwy.riv[['L_riv_100']] == 1] <- 2
cropped.hwy.riv[['M_riv_50']][cropped.hwy.riv[['M_riv_50']] == 1] <- 2.5
# M_riv_50_rast[is.na(M_riv_50_rast)] <- 0

o <- load.WETA()
WETA_2017 <- o[[1]]
covObj <- o[[2]]

# WETA_reproj_coords <- sp::spTransform(WETA_2017[,c('latitude', 'longitude')], crs(cropped.hwy.riv))
prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
w_long_lat <- WETA_2017[,c('longitude', 'latitude')]
names(w_long_lat) <- c("long", "lat")
proj_coords_df <- SpatialPointsDataFrame(coords = w_long_lat, 
                                         data = as.data.frame(WETA_2017$checklist_id), 
                                         proj4string = prj)

WETA_reproj_coords <- spTransform(proj_coords_df, crs(cropped.hwy.riv))

##### PICKUP FROM HERE #####
############################
# find a set of checklists near a highway or river 
# and check if this conductance assignments
# group checklists appropriately

x_val2 <- WETA_reproj_coords[WETA_reproj_coords$`WETA_coords_df$checklist_id` == "",]@coords[[1]]
y_val2 <- WETA_reproj_coords[WETA_reproj_coords$`WETA_coords_df$checklist_id` == "",]@coords[[2]]

xmax2 <- x_val2 + 10000
xmin2 <- x_val2 - 10000
ymax2 <- y_val2 + 10000
ymin2 <- y_val2 - 10000

test_ex2 <- crop(cropped.hwy.riv, c(xmin2, xmax2, ymin2, ymax2))
tr2 <- transition(test_ex2, transition.BRICK = conductanceVals, directions = 8)
trCorr2 <- geoCorrection(tr2, type = 'c')

c2 <- matrix(nrow=4, ncol=2)
c2[1,] <- c(xmax2,ymax2)
c2[2,] <- c(xmax2,ymin2)
c2[3,] <- c(xmin2,ymin2)
c2[4,] <- c(xmin2,ymax2)
points(c2, col="black", cex=2)

g2 <- graph_from_adjacency_matrix(trCorr2@transitionMatrix, mode = "undirected", weighted = T)

cellNum2 <- cellFromXY(trCorr, WETA_in_box@coords)

d2.A <- distances(g2, cellNum2[1], cellNum2[2])
d2.B <- distances(g2, cellNum2[3], cellNum2[4])
d2.C <- distances(g2, cellNum2[4], cellNum2[5])
d2.D <- distances(g2, cellNum2[5], cellNum2[6])
if(!(d2.A > d2.B)){
  print("ERROR")
}

# 
# test_ex1 <- crop(cropped.hwy.riv, c(xmin, xmax, ymin, ymax))
# tr1 <- transition.BRICK(test_ex1, transitionFunction = conductanceVals, directions = 8)
# trCorr <- geoCorrection(tr1, type = 'c')



