library(sp)
library(rgdal)
library(raster)
library(gdistance)
library(sf)


ras.OR <- raster("~/Documents/Oregon State/Research/eBird/site constuction/input data/land cover data/NLCD_2011_OR.tif")
# ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/site constuction/input data/land cover data/OR_2020_fittedImage.tif")

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

identical(crs(reproj_coords), crs(ras.OR))


# points(x, col="red", cex=2)
points(coords, col="red", cex=2)

#############################
# MATH STUFF TO DETERMINE ###
# THE NEW COORDINATES GIVEN #
# OUR TIF IS ANGLED #########
#############################
BL <- as.data.frame(list(x=reproj_coords@bbox[1,1], y=reproj_coords@bbox[2,1]))
TL <- as.data.frame(list(x=reproj_coords@bbox[1,1], y=reproj_coords@bbox[2,2]))
TR <- as.data.frame(list(x=reproj_coords@bbox[1,2], y=reproj_coords@bbox[2,2]))
BR <- as.data.frame(list(x=reproj_coords@bbox[1,2], y=reproj_coords@bbox[2,1]))
x <- rbind(BL, TL, TR, BR)

dist_mat <- as.matrix(dist(x))
THETA <- 23

A <- 1/(cos(THETA*pi/180)/dist_mat[1,2])
B <- sin(THETA*pi/180)*A
B.C <- dist_mat[2,3]
C <- B.C - B
D <- C*sin((90-THETA)*pi/180)

del_x <- D*cos(THETA*pi/180)
del_y <- D*sin(THETA*pi/180)

x4 <- as.data.frame(list(x=TR$x-del_x, y=TR$y+del_y))
x3 <- as.data.frame(list(x=BL$x+del_x, y=BL$y-del_y))

coords <- rbind(x4, TR, x3, BL)
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

#############################
# PLOTS & VALUE EXTRACTION ##
#############################

pal <- colorRampPalette(c("white","black"))
pal2 <- colorRampPalette(c("green", "yellow"))
plot(ras.OR)
plot(cropped, add = T, col="white")

x_val <- mean(c(cropped@extent@xmax, cropped@extent@xmin))
y_val <- mean(c(cropped@extent@ymax, cropped@extent@ymin))

xmax <- x_val + 5000
xmin <- x_val - 5000
ymax <- y_val + 5000
ymin <- y_val - 5000

test_ex1 <- crop(cropped, c(xmin, xmax, ymin, ymax))
tr1 <- transition(test_ex1, transitionFunction = resistanceVals, directions = 8)


