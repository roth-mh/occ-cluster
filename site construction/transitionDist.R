library(sp)
library(rgdal)
library(raster)
library(gdistance)
library(sf)

source("site construction/NLCDvals.R")
source("helper/helpers.R")

#######################################
# NEED TO USE projection from ras.OR2 #
#######################################
ras.OR <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/NLCD_2011_OR.tif")
ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/OR_2020_fittedImage.tif")

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

# points(x, col="red", cex=2)
points(coords, col="red", cex=2)
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

plot(ras.OR)
plot(test, add=T, col="black")
plot(cropped, add=T, col="red")


#############################
# MINI EXAMPLE IN SMALL BOX #
#############################

x_val <- mean(c(cropped@extent@xmax, cropped@extent@xmin)) - 20000
y_val <- mean(c(cropped@extent@ymax, cropped@extent@ymin)) - 20000

xmax <- x_val + 10000
xmin <- x_val - 10000
ymax <- y_val + 10000
ymin <- y_val - 10000

test_ex1 <- crop(cropped, c(xmin, xmax, ymin, ymax))
tr1 <- transition(test_ex1, transitionFunction = resistanceVals, directions = 8)

c <- matrix(nrow =4, ncol=2)
c[1,] <- c(xmax,ymax)
c[2,] <- c(xmax,ymin)
c[3,] <- c(xmin,ymin)
c[4,] <- c(xmin,ymax)
points(c, col="blue", cex=2)
#####################################
# LOAD DATA & REPROJECT WETA COORDS #
#####################################
f_in_WETA <- "../../Class Imbalance/generate syn spec/data/linear/syn_species_1_2017.csv"
WETA_2017_all <- read.delim(f_in_WETA, header=TRUE, sep = ",")
f_in_syn_spec_form <- "../../Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
WETA_2017_region <- subset(WETA_2017_all, WETA_2017_all$latitude <= 44.5)
WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
WETA_2017 <- WETA_2017_region

WETA_coords_df <- as.data.frame(list(long=WETA_2017$longitude, lat=WETA_2017$latitude))
WETA_coords_df$checklist_id <- WETA_2017$checklist_id

prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
WETA_proj_coords_df <- SpatialPointsDataFrame(coords = WETA_coords_df[,1:2], 
                                              data = as.data.frame(WETA_coords_df$checklist_id),
                                              proj4string = prj)



WETA_reproj_coords <- spTransform(WETA_proj_coords_df, crs(ras.OR2))
WETA_reproj_coords

##############################
# SUBSET WETA FOR SMALL BOX; #
# CALC DIST MATRIX           #
##############################
xyz <- subset(WETA_reproj_coords, WETA_reproj_coords@coords[,'long'] >= xmin)
xyz <- subset(xyz, xyz@coords[,'long'] <= xmax)
xyz <- subset(xyz, xyz@coords[,'lat'] >= ymin)
xyz <- subset(xyz, xyz@coords[,'lat'] <= ymax)

cd <- commuteDistance(tr1, xyz@coords)

####################
# DISP DIST MATRIX #
####################

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
dim <- 8
chklsts<-xyz$`WETA_coords_df$checklist_id`
cormat <- get_upper_tri(as.matrix(cd))
colnames(cormat) <- as.character(chklsts)
rownames(cormat) <- as.character(chklsts)
# m <- melt(cormat)

m$value_hundT<-round(m$value/100000, 2)
ggplot(data = m, aes(x=Var1, y=Var2, fill=value_hundT)) + 
  geom_tile() +  geom_text(aes(label = value_hundT), color = "black", size = 3) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))


