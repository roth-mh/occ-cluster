########
# helpers to
# create syn datasets
########
attach_det_vars <- function(df, det_var.df, det_covs=c()){
  if(length(det_covs) > 0){
    col_names <- det_covs
  } else {
    col_names <- names(det_var.df)  
  }
  
  for(i in 1:length(col_names)){
    df[col_names[i]] <- sample(det_var.df[[col_names[i]]], nrow(df), replace = TRUE)
  }
  return(df)
}


############
# This is to be used to prepare
# a dataframe with points and env
# features for use in an experiment.
# it does 2 things:
############
# 1. attaches detection variables
# 2. normalizes vars
prep_syn_df <- function(df, det.df){
  # attaching detection variables
  df.det <- attach_det_vars(df, det.df)
  
  # normalizing all variables besides lat/long 
  # (environmental & detection)
  df_norm.det <- data.frame(
    longitude=df.det$lon, 
    latitude=df.det$lat, 
    site=df.det$site,
    checklist_id=df.det$checklist_id,
    locality_id=paste(as.character(df$lon), as.character(df$lat), sep = "_"),
    base::scale(subset(df.det, select=-c(checklist_id, site, lon, lat)))
  )
  
  return(df_norm.det)
}


create_checklist_ids <- function(n) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
