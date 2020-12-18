##### lambda investigation


calcOccMSELambda <- function(sites_df_occ, covariate_object, true_occ_coefficients, true_det_coefficients, syn_spec=FALSE, enforce_false_positives=TRUE, skip_closure = FALSE, lambda=NULL){
  sites_occ <- subset(sites_df_occ, !duplicated(site))$site
  # this (v v) function is synthetic species specific
  if(skip_closure){
    closed_df <- sites_df_occ
  } else {
    closed_df <- enforceClosure(sites_df_occ, covariate_object$occ_cov, sites_occ, enforce_false_positives)
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
  
  det_cov_str <- paste("", paste(covariate_object$det_cov, collapse="+"), sep=" ~ ")
  occ_cov_str <- paste("", paste(covariate_object$occ_cov, collapse="+"), sep=" ~ ")
  
  species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
  species_formula <- as.formula(species_formula)
  
  occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")
  
  # if(lambda == 0){
  #   og_syn_gen_form <- unmarked::occu(formula = species_formula, occ_um)
  # } else {
  og_syn_gen_form <- unmarked::occuPEN(formula = species_formula, occ_um, pen.type = "Bayes", lambda = lambda)
  # }
  
  
  occ_ex_intercept <- og_syn_gen_form@estimates['state']@estimates
  det_ex_intercept <- og_syn_gen_form@estimates['det']@estimates
  
  
  occ_MSE <- MSE(occ_ex_intercept, true_occ_coefficients)
  det_MSE <- MSE(det_ex_intercept, true_det_coefficients)  
  
  
  
  
  return(list(mse=list(occ=occ_MSE, det=det_MSE), checklists=sites_df_occ, pred_form=og_syn_gen_form))
}


og_data <- WETA_2017

lam_vals = c(0, .0001, .001, .01, .1, 1, 10, 100, 1000)
# occ_mse = list(e_4=0, e_3=0, e_2=0, e_1=0, e0=0, e4=0, e3=0, e2=0, e1=0)
# det_mse = list(e_4=0, e_3=0, e_2=0, e_1=0, e0=0, e4=0, e3=0, e2=0, e1=0)
occ_mse = c(0,0,0,0,0,0,0,0,0)
det_mse = c(0,0,0,0,0,0,0,0,0)
for(x in seq(100)){
  disp("iteration #", as.character(x))
  p_det <- 0
  while(p_det < .15){
    truth_df <- populateDF(kmSq, covObj$siteCovs, covObj$obsCovs, unique(kmSq$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
    og_data$species_observed_syn <- truth_df$species_observed_syn
    n_dets <- sum(truth_df$species_observed_syn)
    n_occ_sites <- sum(truth_df$occupied)
    
    p_det <- n_dets/nrow(truth_df)
    disp("number of detections: ", as.character(n_dets))
    disp("% detections: ", as.character(n_dets/nrow(truth_df)))
    # disp("% detected given its occupied: ", as.character(n_dets/n_occ_sites))
    disp("occupied probability: ", as.character(mean(truth_df$occupied_prob)))
  }
  
  i <- 0
  for(l in lam_vals){
    i <- i + 1
    clust.MSE900 <- calcClustGeoMSE(
      .8, 
      og_data, 
      covObj, 
      num_sites = 900, 
      enforce_false_p = FALSE,
      lam = l
    )
    occ_mse[i] <- occ_mse[i] + clust.MSE900$mse$occ
    det_mse[i] <- det_mse[i] + clust.MSE900$mse$det
  }
}



# lam_vals_char = as.factor(c(".0000", ".0001", ".001", ".01", ".1", "1", "10", "100", "1000"))
lam_vals_char = as.factor(c(".1", "1", "10", "100", "1000"))
plot(lam_vals_char,occ_mse[5:length(occ_mse)]/19, type="l", col="green", lwd=5, xlab="lambda", ylab="mse", main="Occ MSE x Lambda for OccuPEN")
plot(lam_vals_char,det_mse[5:length(det_mse)]/19, type="l", col="red", lwd=5, xlab="lambda", ylab="mse", main="Det MSE x Lambda for OccuPEN")

