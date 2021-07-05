###
# Helper Functions
# Mark Roth
# 6/17/20
###
library(lubridate)
library(dplyr)
library(MLmetrics)
library(auk)
library(unmarked)
library(ClustGeo)
library(data.table)
library(pracma)
library(sp)
library(SDraw)
library(Xmisc)
library(sqldf)

#######
# load covariates from a csv to an R object
#######
loadCovariates <- function(file){
  
  covars <- read.delim(file, header=FALSE, sep = ",", strip.white = TRUE)
  
  occ_cov <- lapply(covars[5,], as.character)
  det_cov <- lapply(covars[6,], as.character)
  
  occCov1 <- as.character(occ_cov[1])
  occCov2 <- as.character(occ_cov[2])
  occCov3 <- as.character(occ_cov[3])
  occCov4 <- as.character(occ_cov[4])
  occCov5 <- as.character(occ_cov[5])
  
  detCov1 <- as.character(det_cov[1])
  detCov2 <- as.character(det_cov[2])
  detCov3 <- as.character(det_cov[3])
  detCov4 <- as.character(det_cov[4])
  detCov5 <- as.character(det_cov[5])
  
  siteCovs <- as.character(occ_cov)
  obsCovs <- as.character(det_cov)
  
  return(list(occCov1=occCov1, occCov2=occCov2, occCov3=occCov3, occCov4=occCov4, occCov5=occCov5,
              detCov1=detCov1, detCov2=detCov2, detCov3=detCov3, detCov4=detCov4, detCov5=detCov5,
              siteCovs=siteCovs, obsCovs=obsCovs, det_cov=det_cov, occ_cov=occ_cov))
}

#######
# 1. selects only pertinent columns
# 2. converts date for AUK analysis
# 3. selects only 2017
#######
groomDataframe <- function(df, det_cov, occ_cov, syn_spec=FALSE){
  if(syn_spec){
    to_keep <- c(paste(det_cov, sep = " "), paste(occ_cov, sep=" "), "species_observed_syn", 
                 "observation_date", "formatted_date", "latitude", "longitude", "locality_id", 
                 "observer_id", "checklist_id", "as_date", "occupied_prob", "occupied")  
  } else {
    to_keep <- c(paste(det_cov, sep = " "), paste(occ_cov, sep=" "), "species_observed", 
                 "observation_date", "formatted_date", "latitude", "longitude", "locality_id", 
                 "observer_id", "checklist_id", "as_date")  
  }
  
  
  df$observation_date <- as.character(df$observation_date)
  df$formatted_date <- mdy(df$observation_date)
  df$as_date <- as_date(df$formatted_date)
  
  df <- subset(df, select = c(to_keep))
  
  df$occupied <- rbinom(length(df$occupied_prob), 1, df$occupied_prob)
  df$species_observed_syn <- rbinom(length(df$species_observed_syn), 1, df$species_observed_syn)
  
  df_2017 <- subset(df, df$as_date >= as_date("2017-01-01"))
  df_2017 <- subset(df_2017, df_2017$as_date < as_date("2018-01-01"))
  return(groomedDF=df_2017)
}

#########
# Rounding Lat/Long
#########
roundLatLong <- function(df, rounding_degree){
  df$rounded_lat <- round(df$latitude, digits = rounding_degree)
  df$rounded_long <- round(df$longitude, digits = rounding_degree)
  df$rounded_locality_id <- paste(as.character(df$rounded_long), as.character(df$rounded_lat), sep = "_")
  
  return(df)
}


# shows stats for the prelim table
showStats <- function(stats_dt){
  disp(sum(stats_dt$num_checklists))
  disp(nrow(stats_dt))
  disp(mean(stats_dt$AREA))
  disp(nrow(stats_dt[stats_dt$AREA > 0,]))
  disp(mean(stats_dt$num_checklists))
  disp(stats_dt[which.max(stats_dt$num_checklists)]$num_checklists)
  disp(nrow(stats_dt[stats_dt$num_checklists == 1,]))
}

#######
# get stats per site; each row is:
#######
#   1. site
#   2. num checklists
#   3. sd of occCov1
#   4. sd of occCov2
#   5. sd of occCov3
#   6. sd of occCov4
#   7. sd of occCov5
#   8. sd of lat
#   9. sd of long
#   10. site AREA calculated 1 of two ways
#     (i) if > 3 distinct checklists at a site
#       of convex hull calculated using 
#       lat/long values and Albert Equal Area projection
#       @ lat1 = 42, lat2 = 48, long1 = -115
#     (ii) if 2 distinct checklists at a site
#       distance (d) between points calculated using AEA 
#       projection then AREA = pi*(d/2)^2
#   11. vAREA = 0 (area of Voronoi polygon)
#######
siteStatsDt <- function(sites_df, covObj){
  sites <- subset(sites_df, !duplicated(site))$site
  numSites <- length(sites)
  
  coords<-data.frame(x = sites_df$longitude, y = sites_df$latitude)
  crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  deer.spdf <- SpatialPointsDataFrame(coords= coords, data = sites_df, proj4string = CRS(crs))
  proj4string(deer.spdf)
  
  Albers.crs <-CRS("+proj=aea +lat_1=42 +lat_2=47
                   +x_0=0+y_0=0 +ellps=GRS80 
                   +towgs84=0,0,0,0,0,0,0 +units=m 
                   +no_defs")
  deer.albers <-spTransform(deer.spdf, CRS=Albers.crs)
  sites_df$long_m <- deer.albers@coords[,1]
  sites_df$lat_m <- deer.albers@coords[,2]
  
  colNames <- c("site", "num_checklists", covObj$occCov1, covObj$occCov2,
                covObj$occCov3, covObj$occCov4, covObj$occCov5, "std_lat", 
                "std_long", "AREA", "vAREA")
  sites_dt <- setNames(data.table(matrix(nrow = numSites, ncol = length(colNames))), colNames)
  sites_dt <- sites_dt[, site:=as.character(site)]
  sites_dt <- sites_dt[, num_checklists:=as.character(num_checklists)]
  for(col in colNames[2:length(colNames)]){
    set(sites_dt, j = col, value = as.numeric(sites_dt[[col]]))
  }
  
  i <- 1
  for(site_i in sites){
    checklists_at_site <- subset(sites_df[sites_df$site == site_i,])
    
    avg_lat = as.numeric(sum(checklists_at_site$latitude)/nrow(checklists_at_site))
    avg_long = as.numeric(sum(checklists_at_site$longitude)/nrow(checklists_at_site))
    
    sd_cov1 = sd(unlist(checklists_at_site[covObj$occCov1]))
    sd_cov2 = sd(unlist(checklists_at_site[covObj$occCov2]))
    sd_cov3 = sd(unlist(checklists_at_site[covObj$occCov3]))
    sd_cov4 = sd(unlist(checklists_at_site[covObj$occCov4]))
    sd_cov5 = sd(unlist(checklists_at_site[covObj$occCov5]))
    
    std_lat = sd(as.numeric(checklists_at_site$latitude))
    std_long = sd(as.numeric(checklists_at_site$longitude))
    # break
    
    distinct_pts <- distinct(checklists_at_site[c("latitude", "longitude")])
    
    if(nrow(distinct_pts) > 1){
      if(nrow(distinct_pts) > 2){
        # area of convex hull
        hpts <- chull(x = checklists_at_site$long_m, y = checklists_at_site$lat_m)
        hpts <- c(hpts, hpts[1])
        xy.coords <- cbind(checklists_at_site$long_m, y = checklists_at_site$lat_m)
        chull.coords <- xy.coords[hpts,]
        area = abs(as.numeric(polyarea(x=chull.coords[,1], y = chull.coords[,2])))
      } else {      # 2 pts
        pts <- subset(checklists_at_site, select = c(lat_m, long_m))
        d <- stats::dist(pts)
        area = pi * (d[1]/2)^2
      }
    } else {
      area = 0
    }
    
    row <- c(site_i, nrow(checklists_at_site), as.double(sd_cov1), 
             as.double(sd_cov2), as.double(sd_cov3), as.double(sd_cov4), 
             as.double(sd_cov5), as.double(std_lat), as.double(std_long),
             as.double(area), vArea <- as.double(0))
    
    set(sites_dt, as.integer(i), names(sites_dt), as.list(row))
    i = i+1
  }
  sites_dt[is.na(sites_dt)] <- 0
  showStats(sites_dt)
  return(sites_dt)
}

##########
# site closure for Occ Model:
#   1. constant site covariates
#   2. no false positives (detected only if occupied)
##########
enforceClosure <- function(sites_df, occ_cov_list, sites_list){
  j<-0
  count <- 0
  # for(eBird_site in 1:numSites){  # for clust geo
  for(eBird_site in sites_list){
    j = j+1
    # if(j %% 100 == 0){
    #   print(j)
    # }
    checklists_at_site <- sites_df[sites_df$site == eBird_site,]
    
    for(occCov_i in occ_cov_list){
      checklists_at_site[occCov_i] <- mean(checklists_at_site[[occCov_i]])
    }
    
    if(j==1){
      clust_geo_df = checklists_at_site
    } else {
      clust_geo_df = rbind(clust_geo_df, checklists_at_site)
    }
    
    # if(sum(checklists_at_site$species_observed_syn) == 0){
    #   count <- count + 1
    # }
    
  }
  # disp("zero detection sites: ", as.character(count))
  # disp("prop of sites: ", as.character(count/length(sites_list)))
  # disp("occupied sites: ", as.character(length(sites_list)-count))
  return(clust_geo_df)
}

##########
# aid analysis; plot a single largest valued sites
##########
plotGroupedSite <- function(site_param, all_sites_df, title, margin, group){
  if(group == "rounded"){
    UB_lat <- site_param$rounded_lat[1] + margin
    LB_lat <- site_param$rounded_lat[1] - margin
    
    UB_long <- site_param$rounded_long[1] + margin
    LB_long <- site_param$rounded_long[1] - margin
    
    # plot the site
    plot(x = site_param$rounded_long, y = site_param$rounded_lat, col="red", xlim=c(LB_long, UB_long), ylim=c(LB_lat, UB_lat), main = title)
    # plot the all checklists
    points(x = all_sites_df$longitude, y = all_sites_df$latitude, col="green")
    # highlight the checklists @ that site
    points(x = site_param$longitude, y = site_param$latitude, col="blue")
    # add legend
    legend("bottomright", legend=c("site", "all checklists", "checklists at site"), col=c("red", "green", "blue"), pch=1)
    
  } else if(group == "clust"){
    
    UB_lat <- site_param$latitude[1] + margin
    LB_lat <- site_param$latitude[1] - margin
    
    UB_long <- site_param$longitude[1] + margin
    LB_long <- site_param$longitude[1] - margin
    
    # highlight the checklists @ that site
    plot(x = site_param$longitude, y = site_param$latitude, col="blue", main = title)#, xlim=c(LB_long, UB_long), ylim=c(LB_lat, UB_lat), )
    # plot the all checklists
    points(x = all_sites_df$longitude, y = all_sites_df$latitude, col="green")
    points(x = site_param$longitude, y = site_param$latitude, col="blue")
    # add legend
    legend("bottomright", legend=c("all checklists", "checklists at site"), col=c("green", "blue"), pch=1)
  }
}

##########
# aid analysis; plot the largest valued sites
##########
plotAllGroupedSites <- function(sites_df, sites_stats_df, title, margin, group){
  max_std_lat_site <- sites_stats_df[which.max(sites_stats_df$std_lat)]
  max_std_long_site <- sites_stats_df[which.max(sites_stats_df$std_long)]
  max_area_site <- sites_stats_df[which.max(sites_stats_df$AREA)]
  
  max_std_occCov1 <- sites_stats_df[which.max(sites_stats_df[[covObj$occCov1]])]
  max_std_occCov2 <- sites_stats_df[which.max(sites_stats_df[[covObj$occCov2]])]
  max_std_occCov3 <- sites_stats_df[which.max(sites_stats_df[[covObj$occCov3]])]
  max_std_occCov4 <- sites_stats_df[which.max(sites_stats_df[[covObj$occCov4]])]
  max_std_occCov5 <- sites_stats_df[which.max(sites_stats_df[[covObj$occCov5]])]
  
  site1 <- sites_df[sites_df$site == max_std_lat_site$site,]
  site2 <- sites_df[sites_df$site == max_std_long_site$site,]
  site3 <- sites_df[sites_df$site == max_area_site$site,]
  
  site4 <- sites_df[sites_df$site == max_std_occCov1$site,]
  site5 <- sites_df[sites_df$site == max_std_occCov2$site,]
  site6 <- sites_df[sites_df$site == max_std_occCov3$site,]
  site7 <- sites_df[sites_df$site == max_std_occCov4$site,]
  site8 <- sites_df[sites_df$site == max_std_occCov5$site,]
  
  plotGroupedSite(site1, sites_df, "site with largest lat std", margin = margin, group = group)
  plotGroupedSite(site2, sites_df, "site with largest long std", margin = margin, group = group)
  plotGroupedSite(site3, sites_df, "site with largest area", margin = margin, group = group)
  
  plotGroupedSite(site4, sites_df, "site with largest occCov1 std", margin = margin, group = group)
  plotGroupedSite(site5, sites_df, "site with largest occCov2 std", margin = margin, group = group)
  plotGroupedSite(site6, sites_df, "site with largest occCov3 std", margin = margin, group = group)
  plotGroupedSite(site7, sites_df, "site with largest occCov4 std", margin = margin, group = group)
  plotGroupedSite(site8, sites_df, "site with largest occCov5 std", margin = margin, group = group)
}

##################
# function to make 
# adding to summary 
# object easier
##################
add_to_summary <- function(res_obj, final_s, TRUE_OCC_COEFF, TRUE_DET_COEFF){
  final_s$det_coeff <- append(final_s$det_coeff, TRUE_DET_COEFF)
  final_s$occ_coeff <- append(final_s$occ_coeff, TRUE_OCC_COEFF)
  final_s$OCCprob <- append(final_s$OCCprob, mean(res_obj$occupied_prob))
  final_s$DETprob <- append(final_s$DETprob, mean(res_obj$det_prob))
  final_s$avg.MSE <- append(final_s$avg.MSE, mean(unlist(res$clustGeo9_occ_mse)))
  final_s$avg.pDet <- append(final_s$avg.pDet, mean(unlist(res_obj$p_det)))
  final_s$avg.zdOccS <- append(final_s$avg.zdOccS, mean(unlist(res_obj$zdOccS)))
  final_s$avg.zds <- append(final_s$avg.zds, mean(unlist(res_obj$zds)))
  return(final_s)
}

##########
# calculates the MSE of Occupancy given true coefficients
#   1. enforces closure assumptions
#   2. formats unmarked occu, given auk formatted object
#   3. determines formula given covariate object
#   4. runs dataset through the occu model in unmarked
#   5. calculates difference in estimated generative model
#       and true generative model
#
# if truth_df is not empty, this means MSE of probabilities will be calculated, rather than
# the MSE of coefficients
##########
calcOccMSE <- function(sites_df_occ, covariate_object, true_occ_coefficients, true_det_coefficients, syn_spec=FALSE, skip_closure=FALSE, truth_df=data.frame()){
  sites_occ <- subset(sites_df_occ, !duplicated(site))$site
  # this (v v) function is synthetic species specific
  if(skip_closure){
    closed_df <- sites_df_occ
  } else {
    closed_df <- enforceClosure(sites_df_occ, covariate_object$siteCovs, sites_occ)
  }
  
  if(syn_spec){
    spec_obs <- "species_observed_syn"
  } else {
    spec_obs <- "species_observed"
  }

  umf_AUK <- auk::format_unmarked_occu(
    closed_df,
    site_id = "site",
    response = spec_obs,
    site_covs = covariate_object$siteCovs,
    obs_covs = covariate_object$obsCovs
  )
  
  det_cov_str <- paste("", paste(covariate_object$obsCovs, collapse="+"), sep=" ~ ")
  occ_cov_str <- paste("", paste(covariate_object$siteCovs, collapse="+"), sep=" ~ ")
  
  species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
  species_formula <- as.formula(species_formula)
  
  occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")
  
  # og_syn_gen_form <- unmarked::occuPEN(formula = species_formula, occ_um, lambda = .1, pen.type = "Bayes")
  og_syn_gen_form <- unmarked::occu(formula = species_formula, occ_um)
  
  
  occ_ex_intercept <- og_syn_gen_form@estimates['state']@estimates
  det_ex_intercept <- og_syn_gen_form@estimates['det']@estimates
  
  occ_MSE <- MSE(occ_ex_intercept, true_occ_coefficients)
  det_MSE <- MSE(det_ex_intercept, true_det_coefficients)
  
  # if(nrow(truth_df) == 0){
  #   prob_occ_MSE <- 0
  #   prob_det_MSE <- 0
  #   
  # } else {
  # calculate the average checklist prob_MSE
  
  closed_df$det_int <- det_ex_intercept[[1]]
  closed_df$det_1 <- closed_df[[names(det_ex_intercept[2])]] * det_ex_intercept[[2]]
  closed_df$det_2 <- closed_df[[names(det_ex_intercept[3])]] * det_ex_intercept[[3]]
  closed_df$det_3 <- closed_df[[names(det_ex_intercept[4])]] * det_ex_intercept[[4]]
  closed_df$det_4 <- closed_df[[names(det_ex_intercept[5])]] * det_ex_intercept[[5]]
  closed_df$det_5 <- closed_df[[names(det_ex_intercept[6])]] * det_ex_intercept[[6]]
  
  closed_df$pred_det_prob <- apply(closed_df, 1, s)
  
  closed_df$occ_int <- occ_ex_intercept[[1]]
  closed_df$occ_1 <- closed_df[[names(occ_ex_intercept[2])]] * occ_ex_intercept[[2]]
  closed_df$occ_2 <- closed_df[[names(occ_ex_intercept[3])]] * occ_ex_intercept[[3]]
  closed_df$occ_3 <- closed_df[[names(occ_ex_intercept[4])]] * occ_ex_intercept[[4]]
  closed_df$occ_4 <- closed_df[[names(occ_ex_intercept[5])]] * occ_ex_intercept[[5]]
  closed_df$occ_5 <- closed_df[[names(occ_ex_intercept[6])]] * occ_ex_intercept[[6]]
  
  closed_df$pred_occupied_prob <- apply(closed_df, 1, s1)
  
  # just split this into 10 separate columns? ... :(
  t_df <- sqldf("SELECT * FROM truth_df WHERE checklist_id IN (SELECT checklist_id FROM closed_df)")
  
  closed_df <- closed_df[order(closed_df$checklist_id),]
  t_df <- t_df[order(t_df$checklist_id),]
  
  t_df$det_diff <- abs(closed_df$pred_det_prob - t_df$det_prob)
  t_df$occ_diff <- abs(closed_df$pred_occupied_prob - t_df$occupied_prob)
  
  prob_occ_MSE <- sum(t_df$occ_diff)/nrow(t_df)
  prob_det_MSE <- sum(t_df$det_diff)/nrow(t_df)
    
  # }
  
  return(list(mse=list(occ=occ_MSE, det=det_MSE), prob_mse=list(occ=prob_occ_MSE, det=prob_det_MSE), checklists=sites_df_occ, pred_form=og_syn_gen_form))
}


s <- function(x, output){
  val <- as.numeric(x[["det_int"]]) + as.numeric(x[["det_1"]]) + as.numeric(x[["det_2"]]) +
    as.numeric(x[["det_3"]]) + as.numeric(x[["det_4"]]) + as.numeric(x[["det_5"]])
  return(expit(val))
}

s1 <- function(x, output){
  val <- as.numeric(x[["occ_int"]]) + as.numeric(x[["occ_1"]]) + as.numeric(x[["occ_2"]]) +
    as.numeric(x[["occ_3"]]) + as.numeric(x[["occ_4"]]) + as.numeric(x[["occ_5"]])
  return(expit(val))
}

##########
# combines a list of list of values and their 
# corresponding named values into a single DF
##########
combineDF <- function(list_obj, list_names){
  if(length(list_obj) != length(list_names)){
    "object and names have different lengths"
    stopifnot(1 == 0)
  }
  res <- data.frame(list_obj[1], row.names = list_names[1])
  for(i in 2:length(list_obj)){
    res <- rbind(res, data.frame(list_obj[[i]], row.names = list_names[i]))  
  }
  return(res)
  
}



calc_zero_det_sites <- function(list_of_checklists){
  list_of_z_d_sites <- list()
  for(chcklsts in list_of_checklists){  
    num_z_d_sites <- 0
    for(zero_d_site in unique(chcklsts$site)){
      detections <- chcklsts[chcklsts$site == zero_d_site,]$species_observed_syn
      if(sum(detections) == 0){
        num_z_d_sites <- num_z_d_sites + 1
      }
    }
    list_of_z_d_sites <- append(list_of_z_d_sites, num_z_d_sites/length(unique(chcklsts$site)))
  }
  return(list_of_z_d_sites)
}


load.WETA <- function(all=FALSE){
  f_in_WETA <- "../../../ICB General/data generation/2017_UPDATED_COVS_df.csv"
  WETA_2017_all <- read.delim(f_in_WETA, header=TRUE, sep = ",")
  # WETA_2017 <- groomDataframe(WETA_2017, covObj$det_cov, covObj$occ_cov, syn_spec = T)
  f_in_syn_spec_form <- "../../Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
  covObj <- loadCovariates(f_in_syn_spec_form)
  covObj$siteCovs <- as.character(c("fall_nbr_TCA_mean_75", 
                                    "fall_nbr_B4_stdDev_150", 
                                    "elevation_stdDev_150", 
                                    "spring_nbr_B7_stdDev_300", 
                                    "aspect_mean_300"))
  WETA_2017_region <- subset(WETA_2017_all, WETA_2017_all$latitude <= 44.5)
  WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
  WETA_2017 <- WETA_2017_region
  if(all){
    return(list(WETA_2017_all, covObj))
  }
  return(list(WETA_2017, covObj))
}


load.WETA_filtered <- function(WETA_2017){
  f_name <- "../../../clusteredSites_2020-12-26_.csv"
  clusted_sites <- read.delim(f_name, sep=",")
  WETA_filtered <- WETA_2017
  WETA_filtered$site <- -1
  WETA_filtered$det_prob <- -1
  
  # link sites with their respective checklists
  for(row in 1:nrow(WETA_filtered)){
    # disp(row)
    long <- WETA_filtered[row, ]$longitude
    lat <- WETA_filtered[row, ]$latitude
    site_obj <- clusted_sites[clusted_sites$longitude == long & clusted_sites$latitude == lat,]
    WETA_filtered[row,]$site <- as.character(site_obj$site)
  }
  return(WETA_filtered)
}

