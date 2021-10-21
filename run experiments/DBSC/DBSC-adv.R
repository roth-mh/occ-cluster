#######
# file to construct adversarial datasets for the DBSC alg
# 10.19.21
# Mark Roth
#######

setwd("/Users/MarkRoth/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/")
source("run experiments/adversarial/syn_helper.R")
source("helper/helpers.R")


tcap.OR <- stack("../TassledCapOR/OR_tasscap_summer_2011_cropped.tif")
names(tcap.OR) <-c("TCB", "TCG", "TCW", "TCA")

sampling_region <- readRDS("../TassledCapOR/sampling_region")


# pts <- data.frame(spsample(sampling_region, 100, type="random")@coords)
# coordinates(pts) <- c("x", "y")
# crs(pts) <- latlong <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# projecting the lat/long points into aea to extract env variables
# pts_proj <- spTransform(pts, CRSobj=crs(tcap.OR))
# 
# plot(tcap.OR$TCB)
# points(pts_proj)


construct_dbsc_dataset <- function(n_pts, min.dist, near.dist, samp_reg, OR_tif, TEST=FALSE){
  interm.df <- data.frame(x=c(), y=c())
  i <- 1
  # phase A: generate initial points 2*min dist away
  # so that when adding @ min dist, we can guarantee 
  # no pt is min_dist away
  while(nrow(interm.df) < n_pts){
    print(i)
    i <- i + 1
    n_pts.remaining <- n_pts - nrow(interm.df)
    
    # 1. sample region
    # TODO: if this takes a long time or n_pts is high,
    # maybe sample more points and filter on those that are far...
    pts <- data.frame(spsample(samp_reg, n_pts.remaining, type="random")@coords)
    
    coordinates(pts) <- c("x", "y")
    crs(pts) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    # projecting the lat/long points into aea to extract env variables
    pts_proj <- spTransform(pts, CRSobj=crs(tcap.OR))
    
    if(nrow(interm.df) > 0){
      interm.df.sp <- interm.df
      coordinates(interm.df.sp) <- c("x", "y")
      crs(interm.df.sp) <- crs(tcap.OR)
      
      pts_proj <- rbind(interm.df.sp, pts_proj)
    }
    
    # 2. calculate distance to other points
    dist <- pointDistance(pts_proj, pts_proj, lonlat = F, allpairs = T)
    
    # 3. find indices for points that are closer than min.dist
    # which(dist < min.dist, arr.ind=TRUE)
    dist <- replace(dist, !upper.tri(dist), min.dist*3 + 1)
    mat.dist <- as.data.frame(which(dist < min.dist*3, arr.ind=TRUE))
    
    if(nrow(mat.dist) != 0){
      far.pts_proj <- pts_proj[-c(mat.dist$row),]
      print(paste0(nrow(mat.dist), " close points"))
      print(paste0("number of rows in interm.df: ", nrow(interm.df)))
      print(paste0("number of rows in far.pts_proj: ", nrow(far.pts_proj@coords)))
      interm.df <- data.frame(far.pts_proj@coords)
    } else {
      print("here")
      dist <- pointDistance(pts_proj, pts_proj, lonlat = F, allpairs = T)
      dist <- replace(dist, !upper.tri(dist), min.dist*3 + 1)
      mat.dist <- as.data.frame(which(dist < min.dist*3, arr.ind=TRUE))
      print(mat.dist)
      break
    }
    

  }
  
  # phase B: add near and far points
  df <- data.frame(pts_proj@coords)
  df$site <- seq(1:nrow(df))
  
  add.df <- data.frame(x=c(),y=c())
  
  max_site <- max(df$site)
  for(row_i in 1:nrow(df)){
    near.pt <- c(x=c(),y=c())
    far.pt <- c(x=c(),y=c())
    if(runif(1) > .5){
      x.dir <- 1
      y.dir <- 0
    } else {
      x.dir <- 0
      y.dir <- 1
    }
    
    if(runif(1) > .5){
      far_pos.or.neg <- 1
      near_pos.or.neg <- -1
    } else {
      far_pos.or.neg <- -1
      near_pos.or.neg <- 1
    }
    
    #TODO: change site def
    pt <- df[row_i,]
    far.pt$x <- pt$x + x.dir * far_pos.or.neg * min.dist
    far.pt$y <- pt$y + y.dir * far_pos.or.neg * min.dist
    far.pt$site <- max_site + 1
    
    near.pt$x <- pt$x + x.dir * near_pos.or.neg * near.dist
    near.pt$y <- pt$y + y.dir * near_pos.or.neg * near.dist
    near.pt$site <- pt$site
    
    if(!TEST){
      add.df <- rbind(add.df, near.pt)  
    }
    add.df <- rbind(add.df, far.pt)
    max_site <- max_site + 1
  }
  
  comb_df <- rbind(df, add.df)
  
  return(comb_df)
}


add_env_vars <- function(df, OR_tif){
  
  ########
  # 1. pts as Spatial points for extraction
  ########
  pts_df.sp <- df
  coordinates(pts_df.sp) <- c('x', 'y')
  crs(pts_df.sp) <- crs(OR_tif)
  
  cids <- create_checklist_ids(nrow(pts_df.sp@coords))
  
  add_df <- data.frame(
    spTransform(pts_df.sp, CRSobj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")),
    checklist_id=cids,
    raster::extract(OR_tif, pts_df.sp)
  )
  
  add_df <- subset(add_df, select=-c(optional))
  
  names(add_df) <- c("site", "lon", "lat", names(add_df)[4:length(names(add_df))])
  return(add_df)
}

##################
obj <- load.WETA()
WETA_2017 <- obj[[1]]
covObj <- obj[[2]]

det.df <- WETA_2017[,unlist(covObj$det_cov)]
##################

pts_df <- construct_dbsc_dataset(200, 1000, near.dist = 5, samp_reg = sampling_region, OR_tif = tcap.OR)

pts.env_df <- add_env_vars(pts_df, OR_tif)

DBSC.adv.df <- prep_syn_df(pts.env_df, det.df)

## testing

coordinates(z$og) <- c("x", "y")
crs(z$og) <- crs(tcap.OR)
coordinates(z$add) <- c("x", "y")
crs(z$add) <- crs(tcap.OR)
comb_df <- rbind(z$og, z$add)
dist_matrix <- pointDistance(comb_df, comb_df, lonlat = F, allpairs = T)
dist <- replace(dist_matrix, !upper.tri(dist_matrix), min.dist + 1)
mat.dist <- as.data.frame(which(dist < min.dist, arr.ind=TRUE))
# 
# 


