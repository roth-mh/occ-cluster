##############
# SKATER Helper:
#   breaking checklists into 
#   Voronoi polygons and finding
#   the neighborhood structure
##############

calcSKATER.MSE <- function(WETA_df, cov_obj, numSites, graph = FALSE, filter = TRUE, enforce_false_p = TRUE){
  if(filter){
    checklists_filtered <- filter_repeat_visits(
      WETA_df,
      min_obs = MIN_OBS, 
      max_obs = MAX_OBS,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id")
    )
  } else {
    checklists_filtered <- WETA_df
  }
  
  checklists_no_dups <- dplyr::distinct(checklists_filtered, latitude, longitude, .keep_all = TRUE)
  spatial_coords <- do.call(rbind, Map(data.frame, long=checklists_no_dups$longitude, lat=checklists_no_dups$latitude))
  sp <- SpatialPoints(spatial_coords)
  vo_poly <- voronoi.polygons(sp)
  
  if(graph){
    plot(vo_poly)
    points(x=vo_poly$x, y=vo_poly$y, col=rainbow(5))
    plot(vo_poly, col=rainbow(5))
    points(x=vo_poly$x, y=vo_poly$y, col="white", pch=20)  
  }
  sites <- findSKATERsites(vo_poly, checklists_no_dups, cov_obj, num_sites = numSites, graph = graph)
  vo_poly$site <- sites$groups
  d<-vo_poly@data
  checklists_filtered$site <- as.character(0)
  for(row in 1:nrow(checklists_filtered)){
    # disp(row)
    long <- checklists_filtered[row, ]$longitude
    lat <- checklists_filtered[row, ]$latitude
    site_obj <- d[d$x == long & d$y == lat,]
    checklists_filtered[row,]$site <- as.character(site_obj$site)
    # if(nrow(site_obj) == 0){
    #   # disp("here")
    #   disp(checklists_filtered[row,]$checklist_id)
    #   disp(checklists_filtered[row,]$latitude)
    #   disp(checklists_filtered[row,]$longitude)
    # }
    
  }
  
  SKATER_MSE <- calcOccMSE(checklists_filtered, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = TRUE, enforce_false_positives = enforce_false_p)
  # TODO: calc stats
  return(SKATER_MSE)
}
