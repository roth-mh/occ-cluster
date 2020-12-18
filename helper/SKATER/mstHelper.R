library(sf)
library(spdep)
library(dplyr)

findSKATERsites <- function(voronoi_p, uniq_checklists, cov_obj, num_sites, graph = FALSE){
  ngbhrs <- poly2nb(voronoi_p)
  env_covs_df <- subset(uniq_checklists, select = cov_obj$siteCovs)
  
  ### calculating costs
  edge_costs <- nbcosts(ngbhrs, env_covs_df)
  mst_in <- nb2listw(ngbhrs, edge_costs, style="B")
  mst <- mstree(mst_in)
  
  if(graph){
    ### the mstree plot
    par(mar=c(0,0,0,0))
    plot(x = voronoi_p$x, y = voronoi_p$y)
    plot(mst, coordinates(voronoi_p), col=2, 
         cex.lab=.6, cex.circles=0.035, fg="blue", add=TRUE)
  }
  
  # TODO: better num_sites than this!
  sites <- skater(mst[,1:2], env_covs_df, ncuts = num_sites)
  
  if(graph){
    plot(sites, coordinates(voronoi_p), cex.circles=0.035, cex.lab=0.00001)
  }
  return(sites)
}


# 
# ### the skater plot
# opar <- par(mar=c(0,0,0,0))
# # attempting to wrangle the sites
# s <- data.frame(table(sites$groups))
# s <- subset(s, s$Freq > 1)
# x <- subset(sites, sites$groups %in% s$Var1)
# plot(x, coordinates(vo_poly), cex.circles=0.035, cex.lab=.7)
# plot(x=vo_poly$x, y=vo_poly$y, col=sites$groups)

####################################
####################################
####################################